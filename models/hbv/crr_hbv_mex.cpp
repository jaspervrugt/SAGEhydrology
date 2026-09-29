/*
 * crr_hbv_mex.cpp
 * MATLAB MEX gateway for the native HBV rainfall-runoff model.
 *
 * This gateway preserves the established crr_hbv interface while keeping
 * mxArray parsing outside the numerical kernel. Double-valued forcing
 * and initial-state arrays are passed to the core without copies; single
 * inputs are converted once at the gateway, matching the old MEX behavior.
 */

#include "mex.h"
#include "hbv.hpp"
#include "../standalone_result.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

const mxArray* field(const mxArray* s, const char* name, bool required = true)
{
    if (!s || !mxIsStruct(s) || mxGetNumberOfElements(s) != 1) {
        mexErrMsgIdAndTxt(
            "crr_hbv:Struct", "Expected a scalar structure while reading '%s'.", name);
    }

    const mxArray* a = mxGetField(s, 0, name);
    if (!a && required) {
        mexErrMsgIdAndTxt("crr_hbv:Field", "Missing field '%s'.", name);
    }
    return a;
}

double scalar(const mxArray* s, const char* name)
{
    const mxArray* a = field(s, name);
    if (mxGetNumberOfElements(a) != 1 || mxIsComplex(a) ||
        !(mxIsDouble(a) || mxIsSingle(a))) {
        mexErrMsgIdAndTxt(
            "crr_hbv:Scalar", "Field '%s' must be a real scalar single or double.", name);
    }
    return mxGetScalar(a);
}

struct InputVector {
    const double* ptr = nullptr;
    std::size_t n = 0;
    std::vector<double> converted;
};

InputVector vector_view(const mxArray* a, const char* name)
{
    if (!a || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("crr_hbv:Vector", "'%s' must be real.", name);
    }

    InputVector v;
    v.n = static_cast<std::size_t>(mxGetNumberOfElements(a));

    if (mxIsDouble(a)) {
        v.ptr = mxGetPr(a);
        return v;
    }

    if (mxIsSingle(a)) {
        const float* p = static_cast<const float*>(mxGetData(a));
        v.converted.resize(v.n);
        for (std::size_t i = 0; i < v.n; ++i) {
            v.converted[i] = static_cast<double>(p[i]);
        }
        v.ptr = v.converted.data();
        return v;
    }

    mexErrMsgIdAndTxt("crr_hbv:Vector", "'%s' must be single or double.", name);
    return v;
}

} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    const bool structured = nrhs == 5;
    if (nrhs != 4 && nrhs != 5) {
        mexErrMsgIdAndTxt("crr_hbv:nrhs", "Expected 4 inputs: t_last, z0, data, options.");
    }

    if ((structured && nlhs != 1) || (!structured && nlhs > 4)) {
        mexErrMsgIdAndTxt("crr_hbv:nlhs", "At most four outputs are supported.");
    }

    if (mxGetNumberOfElements(prhs[0]) != 1 || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("crr_hbv:tlast",
                          "First input must be a real scalar final print time.");
    }

    const double tlast = mxGetScalar(prhs[0]);
    const int ns = static_cast<int>(std::llround(tlast));
    if (!std::isfinite(tlast) || tlast < 1.0 ||
        std::abs(tlast - static_cast<double>(ns)) > 1e-12) {
        mexErrMsgIdAndTxt("crr_hbv:tlast",
                          "Final print time must be a finite integer >= 1.");
    }

    const mxArray* data = prhs[2];
    const mxArray* options = prhs[3];

    InputVector z0 = vector_view(prhs[1], "z0");
    InputVector P = vector_view(field(data, "P"), "data.P");
    InputVector Ep = vector_view(field(data, "Ep"), "data.Ep");
    InputVector T = vector_view(field(data, "T"), "data.T");

    if (P.n < static_cast<std::size_t>(ns) || Ep.n < static_cast<std::size_t>(ns) ||
        T.n < static_cast<std::size_t>(ns)) {
        mexErrMsgIdAndTxt("crr_hbv:dataSize", "P/Ep/T must contain at least ns elements.");
    }

    sage_hbv::Params p;
    p.f_c = scalar(data, "f_c");
    p.beta = scalar(data, "beta");
    p.lp = scalar(data, "lp");
    p.k_0 = scalar(data, "k_0");
    p.uzl = scalar(data, "uzl");
    p.k_1 = scalar(data, "k_1");
    p.k_2 = scalar(data, "k_2");
    p.perc = scalar(data, "perc");
    p.T_tr = scalar(data, "T_tr");
    p.f_dd = scalar(data, "f_dd");
    p.sfcf = scalar(data, "sfcf");
    p.cfr = scalar(data, "cfr");
    p.eps_s = scalar(data, "eps_s");
    p.eps_t = scalar(data, "eps_t");
    p.eps_x = scalar(data, "eps_x");

    sage_hbv::Options opt;
    opt.InitStep = scalar(options, "InitStep");
    opt.MaxStep = scalar(options, "MaxStep");
    opt.MinStep = scalar(options, "MinStep");
    opt.RelTol = scalar(options, "RelTol");
    opt.AbsTol = scalar(options, "AbsTol");
    opt.Order = scalar(options, "Order");
    opt.maxiter = static_cast<int>(std::llround(scalar(options, "maxiter")));

    int mem = 1;
    const mxArray* memA = field(options, "mem", false);
    if (memA && !mxIsEmpty(memA)) {
        mem = static_cast<int>(std::llround(mxGetScalar(memA)));
    }
    std::vector<std::string> on,jn,sn;bool wq=0,wj=0,ws=0,wjs=0,wm=0,wjm=0;if(structured){on=sage_standalone::names(prhs[4],"obs");jn=sage_standalone::names(prhs[4],"jac");sn=sage_standalone::names(prhs[4],"states");wq=sage_standalone::has(on,"Q")||sage_standalone::flag(prhs[4],"q");wj=sage_standalone::has(jn,"Q")||sage_standalone::flag(prhs[4],"jacobian");ws=sage_standalone::has(on,"SWE");wjs=sage_standalone::has(jn,"SWE");wm=sage_standalone::has(on,"SM");wjm=sage_standalone::has(jn,"SM");mem=sn.empty()?0:1;}

    int ipr = 1;
    const mxArray* iprA = field(data, "ipr", false);
    if (iprA && !mxIsEmpty(iprA)) {
        ipr = static_cast<int>(std::llround(mxGetScalar(iprA)));
    }
    if (ipr < 1) {
        ipr = 1;
    }
    if (ipr > ns + 1) {
        ipr = ns + 1;
    }

    constexpr mwSize nvar = 65;
    constexpr mwSize d = 12;
    if (z0.n != nvar) {
        mexErrMsgIdAndTxt("crr_hbv:z0",
                          "HBV augmented initial state must contain 65 elements.");
    }

    const mwSize zrows = (mem != 0) ? static_cast<mwSize>(ns + 1) : 1;
    const mwSize nq = (mem == 0 && ipr <= ns) ? static_cast<mwSize>(ns - ipr + 1) : 0;

    // Allocate MATLAB outputs once and let the core write directly into them.
    mxArray* Zmx = mxCreateDoubleMatrix(zrows, nvar, mxREAL);
    if (!structured) plhs[0] = Zmx;
    double* Z = mxGetPr(Zmx);

    double* q = nullptr;
    double* J = nullptr;
    double *swe=nullptr,*Jswe=nullptr,*sm=nullptr,*Jsm=nullptr;
    std::vector<double> qb,Jb,sb,Jsb,mb,Jmb;

    if (nlhs >= 2) {
        if (mem == 0 && nq > 0) {
            plhs[1] = mxCreateDoubleMatrix(nq, 1, mxREAL);
            q = mxGetPr(plhs[1]);
        } else {
            plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
        }
    }

    if (nlhs >= 3) {
        if (mem == 0 && nq > 0) {
            plhs[2] = mxCreateDoubleMatrix(nq, d, mxREAL);
            J = mxGetPr(plhs[2]);
        } else {
            plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
        }
    }
    if(structured&&!mem){if(wq){qb.assign(nq,0);q=qb.data();}if(wj){Jb.assign(nq*d,0);J=Jb.data();}if(ws){sb.assign(nq,0);swe=sb.data();}if(wjs){Jsb.assign(nq*d,0);Jswe=Jsb.data();}if(wm){mb.assign(nq,0);sm=mb.data();}if(wjm){Jmb.assign(nq*d,0);Jsm=Jmb.data();}}

    sage_hbv::Forcing forcing;
    forcing.P = P.ptr;
    forcing.Ep = Ep.ptr;
    forcing.T = T.ptr;
    forcing.n = std::min(P.n, std::min(Ep.n, T.n));

    sage_hbv::OutputView out;
    out.Z = Z;
    out.q = q;
    out.J = J;
    out.zrows = zrows;
    out.zcols = nvar;
    out.nq = nq;
    out.nj = (J != nullptr) ? d : 0;
    if(structured) out.nj=d;
    out.swe=swe;out.Jswe=Jswe;out.sm=sm;out.Jsm=Jsm;

    const bool fail = sage_hbv::run_into(
        ns,z0.ptr,z0.n,forcing,p,opt,mem!=0,ipr,
        structured?(wj||wjs||wjm):J!=nullptr,out);

    if (structured) {
        if(mem)plhs[0]=sage_standalone::result(prhs[4],Z,zrows,ns,ipr,5,(int)d,fail,{1});else{std::vector<sage_standalone::OutputChannel> c{{"Q",q,J},{"SWE",swe,Jswe},{"SM",sm,Jsm}};plhs[0]=sage_standalone::result(prhs[4],c,nq,(int)d,Z,zrows,ns,ipr,5,fail);}
        mxDestroyArray(Zmx);
    } else if (nlhs >= 4) {
        plhs[3] = mxCreateLogicalScalar(fail);
    }
}
