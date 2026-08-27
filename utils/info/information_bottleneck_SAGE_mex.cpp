/*
 * information_bottleneck_SAGE_mex.cpp
 *
 * Native computational kernel for the SAGE Information diagnostic.
 *
 * MATLAB:
 *   info = information_bottleneck_SAGE_mex(A,H,Y,opts)
 *
 * This MEX performs the numerical work only. It creates no figures and
 * touches no GUI objects. MATLAB/GUI code remains responsible for plotting.
 *
 * Thread controls (all optional):
 *   opts.info_threads  = 1   -> serial (default, safest with parfeval)
 *   opts.info_threads  = N>1 -> explicit OpenMP thread count
 *   opts.info_threads  = 0   -> automatic OpenMP thread count
 *   opts.info_parallel = true/false -> master on/off switch
 *
 * Production outputs include:
 *   - finite-row screening
 *   - quantized neuron/target entropy
 *   - binned neuron-target MI and normalized MI
 *   - entropy-active/effective network width
 *   - bulk KSG neuron-target MI
 *   - permutation-corrected KSG Excess MI and permutation p-values
 *
 * The returned struct intentionally matches the fields consumed by the
 * SAGE Information dashboard:
 *   info.K
 *   info.nLayers
 *   info.nTargets
 *   info.nAttributes
 *   info.keep
 *   info.H_neuron
 *   info.H_neuron_norm
 *   info.active_neuron_mask
 *   info.effective_width
 *   info.I_neuron_theta
 *   info.NMI_neuron_theta
 *   info.H_theta
 *   info.ksg.*
 *   info.uncertainty.ksg.neuron{ell}.*
 *
 * Notes:
 * - KSG values are returned in bits.
 * - Excess MI = observed KSG MI - mean(permuted-target KSG MI).
 * - All requested permutations are executed inside this single MEX call.
 * - No MATLAB<->MEX call is made inside the permutation loop.
 */

#include "mex.h"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

inline bool finite_number(double x) { return mxIsFinite(x); }
inline double qnan() { return mxGetNaN(); }

double digamma_positive(double x)
{
    double result = 0.0;
    while (x < 8.0) {
        result -= 1.0/x;
        x += 1.0;
    }
    const double inv = 1.0/x;
    const double inv2 = inv*inv;
    result += std::log(x) - 0.5*inv
        - inv2*(1.0/12.0
        - inv2*(1.0/120.0
        - inv2*(1.0/252.0
        - inv2*(1.0/240.0
        - inv2*(5.0/660.0)))));
    return result;
}

double get_opt_scalar(const mxArray* opts, const char* field, double def)
{
    if (!opts || !mxIsStruct(opts)) return def;
    const mxArray* a = mxGetField(opts,0,field);
    if (!a || mxIsEmpty(a)) return def;
    if (mxGetNumberOfElements(a)!=1 || mxIsComplex(a) ||
        !(mxIsDouble(a)||mxIsLogical(a))) return def;
    return mxGetScalar(a);
}

bool get_opt_bool(const mxArray* opts, const char* field, bool def)
{
    return get_opt_scalar(opts,field,def?1.0:0.0) != 0.0;
}

std::string get_opt_string(const mxArray* opts, const char* field,
                           const std::string& def)
{
    if (!opts || !mxIsStruct(opts)) return def;
    const mxArray* a = mxGetField(opts,0,field);
    if (!a || mxIsEmpty(a)) return def;
    char* s = mxArrayToString(a);
    if (!s) return def;
    std::string out(s);
    mxFree(s);
    std::transform(out.begin(),out.end(),out.begin(),
                   [](unsigned char c){return static_cast<char>(std::tolower(c));});
    return out;
}

std::vector<mwIndex> get_target_indices(const mxArray* opts,
                                        const char* field,
                                        mwSize P,
                                        bool allDefault)
{
    std::vector<mwIndex> idx;
    const mxArray* a = (opts && mxIsStruct(opts)) ?
        mxGetField(opts,0,field) : nullptr;

    if (!a || mxIsEmpty(a)) {
        if (allDefault) {
            for (mwIndex j=0;j<P;++j) idx.push_back(j);
        } else if (P>0) {
            idx.push_back(0);
        }
        return idx;
    }

    if (!mxIsDouble(a) || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("information_bottleneck_SAGE_mex:indices",
                          "%s must be numeric.",field);
    }

    const double* x = mxGetPr(a);
    const mwSize n = mxGetNumberOfElements(a);
    for (mwSize i=0;i<n;++i) {
        long long v = std::llround(x[i]);
        if (v>=1 && v<=static_cast<long long>(P)) {
            mwIndex z = static_cast<mwIndex>(v-1);
            if (std::find(idx.begin(),idx.end(),z)==idx.end())
                idx.push_back(z);
        }
    }
    return idx;
}

void standardize_columns(std::vector<double>& X,mwSize K,mwSize ncol)
{
    for (mwSize j=0;j<ncol;++j) {
        double mean=0.0;
        for (mwSize i=0;i<K;++i) mean += X[i+K*j];
        mean /= static_cast<double>(K);

        double ss=0.0;
        for (mwSize i=0;i<K;++i) {
            const double d=X[i+K*j]-mean;
            ss += d*d;
        }
        double sd=1.0;
        if (K>1) sd=std::sqrt(ss/static_cast<double>(K-1));
        if (!finite_number(sd)||sd<=0.0) sd=1.0;

        for (mwSize i=0;i<K;++i)
            X[i+K*j]=(X[i+K*j]-mean)/sd;
    }
}

void deterministic_jitter(std::vector<double>& X,mwSize K,mwSize ncol,
                          double amp,int offset)
{
    if (!(amp>0.0)) return;
    const double root2=std::sqrt(2.0);

    for (mwSize j=0;j<ncol;++j) {
        double mean=0.0;
        for (mwSize i=0;i<K;++i) mean += X[i+K*j];
        mean /= static_cast<double>(K);

        double ss=0.0;
        for (mwSize i=0;i<K;++i) {
            const double d=X[i+K*j]-mean;
            ss += d*d;
        }
        double sd=1.0;
        if (K>1) sd=std::sqrt(ss/static_cast<double>(K-1));
        if (!finite_number(sd)||sd<=0.0) sd=1.0;

        const double scale=amp*sd;
        for (mwSize i=0;i<K;++i) {
            const double ii=static_cast<double>(i+1);
            const double jj=static_cast<double>(j+1+offset);
            X[i+K*j] += scale*std::sin(ii*(jj+root2));
        }
    }
}

double ksg_scalar_pair(const double* x,const double* y,mwSize K,int k,
                       std::vector<double>& jointDist)
{
    if (K<=static_cast<mwSize>(k+1)) return qnan();

    double sumPsi=0.0;
    for (mwSize i=0;i<K;++i) {
        const double xi=x[i], yi=y[i];

        mwSize q=0;
        for (mwSize r=0;r<K;++r) {
            if (r==i) continue;
            const double dx=std::fabs(x[r]-xi);
            const double dy=std::fabs(y[r]-yi);
            jointDist[q++]=(dx>dy)?dx:dy;
        }

        auto begin=jointDist.begin();
        auto kth=begin+(k-1);
        auto end=begin+static_cast<std::ptrdiff_t>(K-1);
        std::nth_element(begin,kth,end);
        const double epsilon=*kth;
        if (!finite_number(epsilon)) return qnan();

        const double shrink=std::max(
            1e-12,
            32.0*std::numeric_limits<double>::epsilon()
            *std::max(1.0,epsilon));
        const double radius=std::max(0.0,epsilon-shrink);

        mwSize nx=0,ny=0;
        for (mwSize r=0;r<K;++r) {
            if (r==i) continue;
            if (std::fabs(x[r]-xi)<radius) ++nx;
            if (std::fabs(y[r]-yi)<radius) ++ny;
        }

        sumPsi += digamma_positive(static_cast<double>(nx+1))
                + digamma_positive(static_cast<double>(ny+1));
    }

    const double I_nats =
        digamma_positive(static_cast<double>(k))
        + digamma_positive(static_cast<double>(K))
        - sumPsi/static_cast<double>(K);

    double I_bits=I_nats/std::log(2.0);
    if (I_bits<0.0 && I_bits>-1e-6) I_bits=0.0;
    return I_bits;
}


struct DistanceCache
{
    mwSize K = 0;
    std::vector<double> d;       // row-major K x K absolute distances
    std::vector<double> sorted;  // row-major sorted distances, including self=0
};

DistanceCache build_distance_cache(const double* x,mwSize K)
{
    DistanceCache C;
    C.K = K;
    C.d.resize(K*K);
    C.sorted.resize(K*K);

    for (mwSize i=0;i<K;++i) {
        const double xi=x[i];
        double* row=C.d.data()+i*K;
        double* srow=C.sorted.data()+i*K;

        for (mwSize r=0;r<K;++r) {
            const double v=std::fabs(x[r]-xi);
            row[r]=v;
            srow[r]=v;
        }

        std::sort(srow,srow+K);
    }

    return C;
}

mwSize count_strict_radius(const double* sortedRow,mwSize K,double radius)
{
    if (!(radius>0.0))
        return 0;

    const double* it=std::lower_bound(sortedRow,sortedRow+K,radius);
    mwSize count=static_cast<mwSize>(it-sortedRow);

    /* The sorted row contains the self-distance 0. Remove exactly that
       observation, matching the explicit r ~= i loop in the reference code. */
    if (count>0)
        --count;

    return count;
}

double ksg_from_distance_cache(const DistanceCache& X,
                               const DistanceCache& Y,
                               const mwIndex* perm,
                               int k,
                               std::vector<double>& jointDist)
{
    const mwSize K=X.K;

    if (Y.K!=K || K<=static_cast<mwSize>(k+1))
        return qnan();

    double sumPsi=0.0;

    for (mwSize i=0;i<K;++i) {

        const mwIndex yi = perm ? perm[i] : static_cast<mwIndex>(i);

        const double* dxRow = X.d.data()+i*K;
        const double* dyRow = Y.d.data()+yi*K;

        mwSize q=0;
        for (mwSize r=0;r<K;++r) {
            if (r==i)
                continue;

            const mwIndex yr = perm ? perm[r] : static_cast<mwIndex>(r);
            const double dx=dxRow[r];
            const double dy=dyRow[yr];
            jointDist[q++]=(dx>dy)?dx:dy;
        }

        auto begin=jointDist.begin();
        auto kth=begin+(k-1);
        auto finish=begin+static_cast<std::ptrdiff_t>(K-1);
        std::nth_element(begin,kth,finish);

        const double epsilon=*kth;
        if (!finite_number(epsilon))
            return qnan();

        const double shrink=std::max(
            1e-12,
            32.0*std::numeric_limits<double>::epsilon()
            *std::max(1.0,epsilon));

        const double radius=std::max(0.0,epsilon-shrink);

        const mwSize nx=count_strict_radius(
            X.sorted.data()+i*K,K,radius);

        /* Under a permutation, {perm[r]} is still the complete set of target
           samples. Therefore the marginal y-neighbor count is obtained from
           the pre-sorted distance row of the permuted center yi. */
        const mwSize ny=count_strict_radius(
            Y.sorted.data()+yi*K,K,radius);

        sumPsi += digamma_positive(static_cast<double>(nx+1))
                + digamma_positive(static_cast<double>(ny+1));
    }

    const double I_nats =
        digamma_positive(static_cast<double>(k))
        + digamma_positive(static_cast<double>(K))
        - sumPsi/static_cast<double>(K);

    double I_bits=I_nats/std::log(2.0);

    if (I_bits<0.0 && I_bits>-1e-6)
        I_bits=0.0;

    return I_bits;
}


std::vector<uint16_t> quantize_column(const std::vector<double>& x,
                                      int nbins,
                                      const std::string& method)
{
    const mwSize K=x.size();
    std::vector<uint16_t> q(K,1);

    std::vector<double> uniq=x;
    std::sort(uniq.begin(),uniq.end());
    uniq.erase(std::unique(uniq.begin(),uniq.end()),uniq.end());

    bool integerSmall = uniq.size()<=static_cast<size_t>(nbins);
    if (integerSmall) {
        for (double v:uniq) {
            if (std::fabs(v-std::round(v))>=1e-12) {
                integerSmall=false; break;
            }
        }
    }
    if (integerSmall) {
        for (mwSize i=0;i<K;++i) {
            auto it=std::lower_bound(uniq.begin(),uniq.end(),x[i]);
            q[i]=static_cast<uint16_t>(1+std::distance(uniq.begin(),it));
        }
        return q;
    }

    if (method=="equalwidth") {
        const auto mm=std::minmax_element(x.begin(),x.end());
        const double lo=*mm.first, hi=*mm.second;
        if (!(finite_number(lo)&&finite_number(hi)&&hi>lo)) return q;
        const double span=hi-lo;
        for (mwSize i=0;i<K;++i) {
            int b=static_cast<int>(std::floor((x[i]-lo)/span*nbins));
            if (b<0) b=0;
            if (b>=nbins) b=nbins-1;
            q[i]=static_cast<uint16_t>(b+1);
        }
        return q;
    }

    /* Quantile binning matched to MATLAB's default midpoint convention.
       For sorted sample x_(k), k=1,...,K, the probability locations are

           p_k = (2*k - 1) / (2*K).

       Quantiles are linearly interpolated between these midpoint
       probabilities, with clipping to min/max at the tails.

       Bin assignment then matches MATLAB discretize default semantics:
           [edge_j, edge_{j+1}) for interior bins,
       with the last bin including the right endpoint.

       MATLAB reference:

           edges = quantile(x,linspace(0,1,nbins+1));
           edges(1)   = -inf;
           edges(end) =  inf;
           edges = unique(edges,'stable');
           q = discretize(x,edges);
    */
    std::vector<double> sx=x;
    std::sort(sx.begin(),sx.end());

    std::vector<double> edges(nbins+1);

    const double Kd = static_cast<double>(K);
    const double pFirst = 1.0/(2.0*Kd);
    const double pLast  = (2.0*Kd - 1.0)/(2.0*Kd);

    for (int b=0;b<=nbins;++b) {

        const double p = static_cast<double>(b)
                       / static_cast<double>(nbins);

        if (p <= pFirst) {
            edges[b] = sx.front();
            continue;
        }

        if (p >= pLast) {
            edges[b] = sx.back();
            continue;
        }

        /* Find k such that p_k <= p < p_(k+1).
           With zero-based index j=k-1:

               p_j = (2*j + 1)/(2*K)

           and consecutive midpoint probabilities differ by 1/K.
        */
        const double t = Kd*p - 0.5;
        mwSize j = static_cast<mwSize>(std::floor(t));

        if (j >= K-1)
            j = K-2;

        const double pLo = (2.0*static_cast<double>(j) + 1.0)
                         /(2.0*Kd);
        const double w = (p - pLo)*Kd;

        edges[b] = (1.0-w)*sx[j] + w*sx[j+1];
    }

    /* MATLAB replaces the first/last quantile edge by +/-Inf. We do not
       need to store infinity explicitly for assignment; only the unique
       interior edges determine the bins. */
    std::vector<double> uedges;
    uedges.reserve(edges.size());

    for (double e : edges) {
        if (uedges.empty() || e != uedges.back())
            uedges.push_back(e);
    }

    if (uedges.size() <= 2)
        return q;

    /* MATLAB discretize uses left-closed/right-open bins. Therefore a
       sample exactly equal to an interior edge belongs to the bin on its
       right. std::upper_bound reproduces that behavior. */
    for (mwSize i=0;i<K;++i) {
        auto it = std::upper_bound(
            uedges.begin()+1,uedges.end()-1,x[i]);

        const int b = static_cast<int>(
            std::distance(uedges.begin()+1,it));

        q[i] = static_cast<uint16_t>(b+1);
    }

    return q;
}

double entropy_discrete(const std::vector<uint16_t>& x)
{
    if (x.empty()) return qnan();
    std::unordered_map<uint16_t,mwSize> counts;
    for (uint16_t v:x) ++counts[v];

    const double n=static_cast<double>(x.size());
    double H=0.0;
    for (const auto& kv:counts) {
        const double p=kv.second/n;
        if (p>0.0) H -= p*(std::log(p)/std::log(2.0));
    }
    return H;
}

double mi_discrete(const std::vector<uint16_t>& x,
                   const std::vector<uint16_t>& y)
{
    if (x.empty() || x.size()!=y.size()) return qnan();

    std::unordered_map<uint16_t,mwSize> cx,cy;
    std::unordered_map<uint32_t,mwSize> cxy;
    for (mwSize i=0;i<x.size();++i) {
        ++cx[x[i]]; ++cy[y[i]];
        uint32_t key=(static_cast<uint32_t>(x[i])<<16)
                    |static_cast<uint32_t>(y[i]);
        ++cxy[key];
    }

    const double n=static_cast<double>(x.size());
    double I=0.0;
    for (const auto& kv:cxy) {
        const uint16_t xv=static_cast<uint16_t>(kv.first>>16);
        const uint16_t yv=static_cast<uint16_t>(kv.first&0xffffu);
        const double pxy=kv.second/n;
        const double px=cx[xv]/n;
        const double py=cy[yv]/n;
        I += pxy*(std::log(pxy/(px*py))/std::log(2.0));
    }
    if (I<0.0 && I>-1e-12) I=0.0;
    return I;
}

mxArray* make_double_row(const std::vector<double>& v)
{
    mxArray* a=mxCreateDoubleMatrix(1,v.size(),mxREAL);
    std::copy(v.begin(),v.end(),mxGetPr(a));
    return a;
}

mxArray* make_double_matrix(mwSize m,mwSize n,const std::vector<double>& v)
{
    mxArray* a=mxCreateDoubleMatrix(m,n,mxREAL);
    std::copy(v.begin(),v.end(),mxGetPr(a));
    return a;
}

mxArray* make_logical_row(const std::vector<uint8_t>& v)
{
    mxArray* a=mxCreateLogicalMatrix(1,v.size());
    mxLogical* p=mxGetLogicals(a);
    for (mwSize i=0;i<v.size();++i) p[i]=v[i]?1:0;
    return a;
}

mxArray* make_index_row(const std::vector<mwIndex>& idx)
{
    mxArray* a=mxCreateDoubleMatrix(1,idx.size(),mxREAL);
    double* p=mxGetPr(a);
    for (mwSize i=0;i<idx.size();++i) p[i]=static_cast<double>(idx[i]+1);
    return a;
}

void require_real_double_matrix(const mxArray* a,const char* name)
{
    if (!a || !mxIsDouble(a) || mxIsComplex(a) || mxIsSparse(a)) {
        mexErrMsgIdAndTxt("information_bottleneck_SAGE_mex:type",
                          "%s must be a full real double matrix.",name);
    }
}

} // namespace

void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[])
{
    if (nrhs<3 || nrhs>4)
        mexErrMsgIdAndTxt("information_bottleneck_SAGE_mex:nrhs",
            "Usage: info = information_bottleneck_SAGE_mex(A,H,Y[,opts])");
    if (nlhs>1)
        mexErrMsgIdAndTxt("information_bottleneck_SAGE_mex:nlhs",
            "One output is supported.");

    const mxArray* Amx=prhs[0];
    const mxArray* Hmx=prhs[1];
    const mxArray* Ymx=prhs[2];
    const mxArray* opts=(nrhs>=4)?prhs[3]:nullptr;

    require_real_double_matrix(Amx,"A");
    require_real_double_matrix(Ymx,"Y");
    if (!mxIsCell(Hmx) || mxGetNumberOfElements(Hmx)<1)
        mexErrMsgIdAndTxt("information_bottleneck_SAGE_mex:H",
                          "H must be a nonempty cell array.");

    const mwSize K0=mxGetM(Amx);
    const mwSize D=mxGetN(Amx);
    const mwSize P=mxGetN(Ymx);
    const mwSize L=mxGetNumberOfElements(Hmx);

    if (mxGetM(Ymx)!=K0)
        mexErrMsgIdAndTxt("information_bottleneck_SAGE_mex:size",
                          "A and Y must have the same number of rows.");

    std::vector<const mxArray*> Hsrc(L);
    std::vector<mwSize> N(L);
    for (mwIndex ell=0;ell<L;++ell) {
        Hsrc[ell]=mxGetCell(Hmx,ell);
        require_real_double_matrix(Hsrc[ell],"H{ell}");
        if (mxGetM(Hsrc[ell])!=K0)
            mexErrMsgIdAndTxt("information_bottleneck_SAGE_mex:size",
                              "All H{ell} must have K rows.");
        N[ell]=mxGetN(Hsrc[ell]);
    }

    const double* Ain=mxGetPr(Amx);
    const double* Yin=mxGetPr(Ymx);

    std::vector<uint8_t> keep0(K0,1);
    for (mwSize i=0;i<K0;++i) {
        for (mwSize j=0;j<D && keep0[i];++j)
            if (!finite_number(Ain[i+K0*j])) keep0[i]=0;
        for (mwSize j=0;j<P && keep0[i];++j)
            if (!finite_number(Yin[i+K0*j])) keep0[i]=0;
        for (mwIndex ell=0;ell<L && keep0[i];++ell) {
            const double* h=mxGetPr(Hsrc[ell]);
            for (mwSize j=0;j<N[ell];++j)
                if (!finite_number(h[i+K0*j])) {keep0[i]=0;break;}
        }
    }

    std::vector<mwIndex> rows;
    rows.reserve(K0);
    for (mwIndex i=0;i<K0;++i) if (keep0[i]) rows.push_back(i);
    const mwSize K=rows.size();
    if (K<5)
        mexErrMsgIdAndTxt("information_bottleneck_SAGE_mex:samples",
                          "Fewer than five finite samples remain.");

    std::vector<double> Y(K*P);
    for (mwSize j=0;j<P;++j)
        for (mwSize ii=0;ii<K;++ii)
            Y[ii+K*j]=Yin[rows[ii]+K0*j];

    std::vector<std::vector<double>> H(L);
    for (mwIndex ell=0;ell<L;++ell) {
        H[ell].resize(K*N[ell]);
        const double* hin=mxGetPr(Hsrc[ell]);
        for (mwSize j=0;j<N[ell];++j)
            for (mwSize ii=0;ii<K;++ii)
                H[ell][ii+K*j]=hin[rows[ii]+K0*j];
    }

    const int nbins=std::max(2,(int)std::llround(
        get_opt_scalar(opts,"nbins",10)));
    const int targetBins=std::max(2,(int)std::llround(
        get_opt_scalar(opts,"target_nbins",8)));
    const std::string quant=get_opt_string(opts,"quantization","quantile");
    const double activeThr=std::min(1.0,std::max(0.0,
        get_opt_scalar(opts,"active_entropy_threshold",0.10)));
    int k=std::max(1,(int)std::llround(get_opt_scalar(opts,"ksg_k",5)));
    if (K<=static_cast<mwSize>(k+1)) k=std::max(1,(int)K-2);
    const bool standardize=get_opt_bool(opts,"ksg_standardize",true);
    const double jitter=std::max(0.0,get_opt_scalar(opts,"ksg_jitter",1e-10));
    const bool doKsg=get_opt_bool(opts,"ksg",true);
    const bool uncertaintyNeurons=
        get_opt_bool(opts,"ksg_uncertainty_neurons",false);
    const int nPerm=std::max(0,(int)std::llround(
        get_opt_scalar(opts,"ksg_permutation_n",0)));
    const unsigned int seed=(unsigned int)std::llround(
        get_opt_scalar(opts,"random_seed",1729));

    /* Optional native threading controls.
       Default = serial (1 thread), which is safest when SAGE already uses
       MATLAB parallelism (for example parfeval workers).
         info_threads = 1  -> serial
         info_threads = 0  -> automatic maximum available OpenMP threads
         info_threads = N  -> explicit N OpenMP threads
       info_parallel = false forces serial mode regardless of info_threads.
    */
    int infoThreads =
        std::max(0,(int)std::llround(get_opt_scalar(opts,"info_threads",1)));
    const bool infoParallel =
        get_opt_bool(opts,"info_parallel",infoThreads != 1);

#ifdef _OPENMP
    if (infoThreads == 0)
        infoThreads = omp_get_max_threads();
#else
    if (infoThreads == 0)
        infoThreads = 1;
#endif

    if (!infoParallel || infoThreads < 1)
        infoThreads = 1;

    /* Match the MATLAB reference defaults: target 1 unless explicitly set. */
    const std::vector<mwIndex> ksgTargets=
        get_target_indices(opts,"ksg_target_indices",P,false);
    const std::vector<mwIndex> uncTargets=
        get_target_indices(opts,"ksg_uncertainty_target_indices",P,false);

    /* Quantize targets. */
    std::vector<std::vector<uint16_t>> Yq(P);
    std::vector<double> Htheta(P,qnan());
    for (mwSize j=0;j<P;++j) {
        std::vector<double> col(K);
        for (mwSize i=0;i<K;++i) col[i]=Y[i+K*j];
        Yq[j]=quantize_column(col,targetBins,quant);
        Htheta[j]=entropy_discrete(Yq[j]);
    }

    /* Allocate layer-level containers. */
    std::vector<std::vector<double>> Hneuron(L), HneuronNorm(L);
    std::vector<std::vector<uint8_t>> active(L);
    std::vector<double> effectiveWidth(L,qnan());
    std::vector<std::vector<double>> Ibinned(L), INMIBinned(L);
    std::vector<std::vector<double>> IKSG(L);
    std::vector<std::vector<double>> permMean(L), permP(L), excess(L);

    const double HmaxNeuron=std::log2(static_cast<double>(
        std::min<mwSize>(static_cast<mwSize>(nbins),K)));

    for (mwIndex ell=0;ell<L;++ell) {
        Hneuron[ell].assign(N[ell],qnan());
        HneuronNorm[ell].assign(N[ell],qnan());
        active[ell].assign(N[ell],0);
        Ibinned[ell].assign(N[ell]*P,qnan());
        INMIBinned[ell].assign(N[ell]*P,qnan());
        IKSG[ell].assign(N[ell]*P,qnan());
        permMean[ell].assign(N[ell]*P,qnan());
        permP[ell].assign(N[ell]*P,qnan());
        excess[ell].assign(N[ell]*P,qnan());

        for (mwSize n=0;n<N[ell];++n) {
            std::vector<double> x(K);
            for (mwSize i=0;i<K;++i) x[i]=H[ell][i+K*n];
            const auto qx=quantize_column(x,nbins,quant);

            const double hn=entropy_discrete(qx);
            Hneuron[ell][n]=hn;
            double hnorm=qnan();
            if (finite_number(HmaxNeuron) && HmaxNeuron>0.0)
                hnorm=std::min(1.0,std::max(0.0,hn/HmaxNeuron));
            HneuronNorm[ell][n]=hnorm;
            if (finite_number(hnorm) && hnorm>=activeThr) active[ell][n]=1;

            for (mwSize j=0;j<P;++j) {
                const double I=mi_discrete(qx,Yq[j]);
                Ibinned[ell][n+N[ell]*j]=I;
                if (finite_number(Htheta[j]) && Htheta[j]>0.0)
                    INMIBinned[ell][n+N[ell]*j]=I/Htheta[j];
            }
        }
        effectiveWidth[ell]=std::accumulate(active[ell].begin(),
            active[ell].end(),0.0);
    }

    /* Pre-standardize/jitter copies once for observed KSG. */
    std::vector<std::vector<double>> Hk=H;
    std::vector<double> Yk=Y;
    if (standardize) {
        for (mwIndex ell=0;ell<L;++ell) standardize_columns(Hk[ell],K,N[ell]);
        standardize_columns(Yk,K,P);
    }
    if (jitter>0.0) {
        int offset=1;
        for (mwIndex ell=0;ell<L;++ell) {
            deterministic_jitter(Hk[ell],K,N[ell],jitter,offset);
            offset += static_cast<int>(N[ell])+17;
        }
        deterministic_jitter(Yk,K,P,jitter,101);
    }

    std::vector<double> dist(K-1);

    if (doKsg) {
        for (mwIndex ell=0;ell<L;++ell) {
            for (mwIndex jj=0;jj<ksgTargets.size();++jj) {
                const mwIndex j=ksgTargets[jj];
                const double* y=Yk.data()+K*j;

                for (mwSize n=0;n<N[ell];++n) {
                    const double* x=Hk[ell].data()+K*n;
                    IKSG[ell][n+N[ell]*j]=
                        ksg_scalar_pair(x,y,K,k,dist);
                }
            }
        }
    }

    /* ---------------------------------------------------------------
       Optimized permutation-corrected neuron KSG.

       The expensive marginal distances do not change across target
       permutations. For each target/layer we therefore:

         1. build the target absolute-distance matrix once;
         2. generate/store the requested permutations once;
         3. build one neuron's distance matrix once;
         4. reuse both distance matrices for every permutation.

       Only the joint max-distance vector and kth-neighbor selection must
       still be evaluated for each permutation. This preserves the KSG-1
       statistic while removing repeated fabs(), marginal distance
       construction, and marginal neighbor scans.

       This block follows the MATLAB semantics: neuron-level uncertainty is
       evaluated only when opts.ksg_uncertainty_neurons is true.
       --------------------------------------------------------------- */
    if (doKsg && uncertaintyNeurons
            && nPerm>0 && !uncTargets.empty()) {

        std::mt19937 rng(seed+104729u);
        std::vector<mwIndex> order(K);
        std::iota(order.begin(),order.end(),0);

        for (mwIndex ell=0;ell<L;++ell) {

            std::vector<double> sum(N[ell]*P,0.0);
            std::vector<double> ge(N[ell]*P,0.0);

            for (mwIndex jj=0;jj<uncTargets.size();++jj) {

                const mwIndex j=uncTargets[jj];
                const double* y=Yk.data()+K*j;

                /* Cache target distances once for all neurons and all
                   permutations of this target. */
                const DistanceCache yCache=
                    build_distance_cache(y,K);

                /* Generate exactly nPerm permutations once. Reusing the same
                   permutation set for all neurons matches the previous native
                   implementation while allowing neuron-by-neuron caching. */
                std::vector<mwIndex> permutations(
                    static_cast<size_t>(nPerm)*K);

                for (int b=0;b<nPerm;++b) {
                    std::shuffle(order.begin(),order.end(),rng);
                    std::copy(order.begin(),order.end(),
                        permutations.begin()+static_cast<size_t>(b)*K);
                }

#ifdef _OPENMP
#pragma omp parallel for if(infoThreads > 1) num_threads(infoThreads) schedule(dynamic)
#endif
                for (mwIndex n=0;n<static_cast<mwIndex>(N[ell]);++n) {

                    std::vector<double> distLocal(K-1);

                    const mwSize nn = static_cast<mwSize>(n);
                    const mwSize ij=nn+N[ell]*j;
                    const double* x=Hk[ell].data()+K*nn;

                    /* Cache this neuron's distances once and reuse them for
                       observed and all permuted KSG evaluations. */
                    const DistanceCache xCache=
                        build_distance_cache(x,K);

                    if (!finite_number(IKSG[ell][ij])) {
                        IKSG[ell][ij]=ksg_from_distance_cache(
                            xCache,yCache,nullptr,k,distLocal);
                    }

                    double s=0.0;
                    double nge=0.0;
                    mwSize nFinite=0;

                    for (int b=0;b<nPerm;++b) {

                        const mwIndex* perm=
                            permutations.data()
                            +static_cast<size_t>(b)*K;

                        const double Ip=ksg_from_distance_cache(
                            xCache,yCache,perm,k,distLocal);

                        if (finite_number(Ip)) {
                            s += Ip;
                            ++nFinite;

                            if (Ip>=IKSG[ell][ij])
                                nge += 1.0;
                        }
                    }

                    if (nFinite>0) {
                        permMean[ell][ij]=
                            s/static_cast<double>(nFinite);

                        permP[ell][ij]=
                            (nge+1.0)
                            /static_cast<double>(nFinite+1);

                        excess[ell][ij]=
                            IKSG[ell][ij]-permMean[ell][ij];
                    }
                }
            }
        }
    }

    /* Build MATLAB output struct. */
    const char* infoFields[]={
        "K","nLayers","nTargets","nAttributes","keep",
        "H_theta","H_neuron","H_neuron_norm","active_neuron_mask",
        "effective_width","I_neuron_theta","NMI_neuron_theta",
        "ksg","uncertainty","sample","options"
    };
    mxArray* info=mxCreateStructMatrix(1,1,
        sizeof(infoFields)/sizeof(infoFields[0]),infoFields);

    mxSetField(info,0,"K",mxCreateDoubleScalar(static_cast<double>(K)));
    mxSetField(info,0,"nLayers",mxCreateDoubleScalar(static_cast<double>(L)));
    mxSetField(info,0,"nTargets",mxCreateDoubleScalar(static_cast<double>(P)));
    mxSetField(info,0,"nAttributes",mxCreateDoubleScalar(static_cast<double>(D)));

    mxArray* keepMx=mxCreateLogicalMatrix(K0,1);
    mxLogical* kp=mxGetLogicals(keepMx);
    for (mwSize i=0;i<K0;++i) kp[i]=keep0[i]?1:0;
    mxSetField(info,0,"keep",keepMx);
    mxSetField(info,0,"H_theta",make_double_row(Htheta));
    mxSetField(info,0,"effective_width",make_double_row(effectiveWidth));

    mxArray* cH=mxCreateCellMatrix(1,L);
    mxArray* cHN=mxCreateCellMatrix(1,L);
    mxArray* cAct=mxCreateCellMatrix(1,L);
    mxArray* cIB=mxCreateCellMatrix(1,L);
    mxArray* cIN=mxCreateCellMatrix(1,L);
    for (mwIndex ell=0;ell<L;++ell) {
        mxSetCell(cH,ell,make_double_row(Hneuron[ell]));
        mxSetCell(cHN,ell,make_double_row(HneuronNorm[ell]));
        mxSetCell(cAct,ell,make_logical_row(active[ell]));
        mxSetCell(cIB,ell,make_double_matrix(N[ell],P,Ibinned[ell]));
        mxSetCell(cIN,ell,make_double_matrix(N[ell],P,INMIBinned[ell]));
    }
    mxSetField(info,0,"H_neuron",cH);
    mxSetField(info,0,"H_neuron_norm",cHN);
    mxSetField(info,0,"active_neuron_mask",cAct);
    mxSetField(info,0,"I_neuron_theta",cIB);
    mxSetField(info,0,"NMI_neuron_theta",cIN);

    const char* kFields[]={
        "k","method","target_indices","uncertainty_target_indices",
        "I_neuron_theta","threads_requested","threads_used","openmp_enabled"
    };
    mxArray* ksg=mxCreateStructMatrix(1,1,8,kFields);
    mxSetField(ksg,0,"k",mxCreateDoubleScalar(k));
    mxSetField(ksg,0,"method",mxCreateString("mex-core"));
    mxSetField(ksg,0,"target_indices",make_index_row(ksgTargets));
    mxSetField(ksg,0,"uncertainty_target_indices",make_index_row(uncTargets));
    mxSetField(ksg,0,"threads_requested",
               mxCreateDoubleScalar(get_opt_scalar(opts,"info_threads",1)));
    mxSetField(ksg,0,"threads_used",
               mxCreateDoubleScalar(infoThreads));
#ifdef _OPENMP
    mxSetField(ksg,0,"openmp_enabled",mxCreateLogicalScalar(true));
#else
    mxSetField(ksg,0,"openmp_enabled",mxCreateLogicalScalar(false));
#endif

    mxArray* cK=mxCreateCellMatrix(1,L);
    for (mwIndex ell=0;ell<L;++ell)
        mxSetCell(cK,ell,make_double_matrix(N[ell],P,IKSG[ell]));
    mxSetField(ksg,0,"I_neuron_theta",cK);
    mxSetField(info,0,"ksg",ksg);

    const char* uFields[]={"ksg"};
    mxArray* uncertainty=mxCreateStructMatrix(1,1,1,uFields);
    const char* ukFields[]={"neuron"};
    mxArray* uksg=mxCreateStructMatrix(1,1,1,ukFields);
    mxArray* cUN=mxCreateCellMatrix(1,L);

    const char* nFields[]={
        "bootstrap_mean","bootstrap_lo","bootstrap_hi",
        "permutation_mean","permutation_p","excess"
    };
    for (mwIndex ell=0;ell<L;++ell) {
        mxArray* UN=mxCreateStructMatrix(1,1,6,nFields);
        std::vector<double> nanv(N[ell]*P,qnan());
        mxSetField(UN,0,"bootstrap_mean",make_double_matrix(N[ell],P,nanv));
        mxSetField(UN,0,"bootstrap_lo",make_double_matrix(N[ell],P,nanv));
        mxSetField(UN,0,"bootstrap_hi",make_double_matrix(N[ell],P,nanv));
        mxSetField(UN,0,"permutation_mean",
                   make_double_matrix(N[ell],P,permMean[ell]));
        mxSetField(UN,0,"permutation_p",
                   make_double_matrix(N[ell],P,permP[ell]));
        mxSetField(UN,0,"excess",
                   make_double_matrix(N[ell],P,excess[ell]));
        mxSetCell(cUN,ell,UN);
    }
    mxSetField(uksg,0,"neuron",cUN);
    mxSetField(uncertainty,0,"ksg",uksg);
    mxSetField(info,0,"uncertainty",uncertainty);

    const char* sFields[]={
        "K","max_empirical_entropy","target_bin_occupancy",
        "hidden_bin_occupancy"
    };
    mxArray* sample=mxCreateStructMatrix(1,1,4,sFields);
    mxSetField(sample,0,"K",mxCreateDoubleScalar(K));
    mxSetField(sample,0,"max_empirical_entropy",
               mxCreateDoubleScalar(std::log2(static_cast<double>(K))));
    mxSetField(sample,0,"target_bin_occupancy",
               mxCreateDoubleScalar(static_cast<double>(K)/targetBins));
    mxSetField(sample,0,"hidden_bin_occupancy",
               mxCreateDoubleScalar(static_cast<double>(K)/nbins));
    mxSetField(info,0,"sample",sample);
    if (opts && mxIsStruct(opts))
        mxSetField(info,0,"options",mxDuplicateArray(opts));
    else
        mxSetField(info,0,"options",mxCreateStructMatrix(1,1,0,nullptr));

    plhs[0]=info;
}
