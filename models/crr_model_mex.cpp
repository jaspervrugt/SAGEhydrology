
/*
 * crr_model_mex.cpp
 * Native SAGE CRR_MODEL branch -- direct validated model cores.
 *
 * Supported built-in models:
 *   1 HYMOD, 2 HMODEL, 3 SAC-SMA, 4 Xinanjiang,
 *   5 GR4J-A, 6 HBV, 7 CFE-NWM, 21 GR4J-B.
 *
 * Native path:
 *   model preparation -> direct extracted C++ core -> q/J ->
 *   training mask -> loss/delta -> J' delta.
 *
 * V3 native objectives: 1 SAR, 2 RSS/GLS(identity), 3 NSE,
 *                       4 KGE, 5 Huber, 6 FDC, 7 JKGE.
 * Native requested outputs: q, gradient, jacobian, metrics, attribution, states.
 *
 * The MATLAB reference crr_model.m remains untouched.
 */

#include "mex.h"
#include "matrix.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "hymod/hymod.hpp"
#include "hmodel/hmodel.hpp"
#include "sacsma/sacsma.hpp"
#include "Xinanjiang/xinanjiang.hpp"
#include "gr4jA/gr4jA.hpp"
#include "gr4jB/gr4jB.hpp"
#include "hbv/hbv.hpp"
#include "cfe_nwm/cfe_nwm.hpp"

// Built-in models call the validated native cores directly.

namespace {
constexpr double QNAN = std::numeric_limits<double>::quiet_NaN();

const mxArray* fld(const mxArray* s, const char* name, bool required = true)
{
    if (!s || !mxIsStruct(s) || mxGetNumberOfElements(s) != 1) {
        mexErrMsgIdAndTxt(
            "crr_model_mex:Struct", "Expected scalar struct while reading '%s'.", name);
    }
    const mxArray* a = mxGetField(s, 0, name);
    if (!a && required) {
        mexErrMsgIdAndTxt("crr_model_mex:Field", "Missing field '%s'.", name);
    }
    return a;
}

double scl(const mxArray* a, const char* name)
{
    if (!a || mxIsComplex(a) || mxGetNumberOfElements(a) != 1 ||
        !(mxIsNumeric(a) || mxIsLogical(a))) {
        mexErrMsgIdAndTxt("crr_model_mex:Scalar", "'%s' must be a real scalar.", name);
    }
    return mxGetScalar(a);
}

bool lflag(const mxArray* s, const char* name, bool def = false)
{
    const mxArray* a = fld(s, name, false);
    if (!a || mxIsEmpty(a)) {
        return def;
    }
    if (!mxIsLogicalScalar(a)) {
        mexErrMsgIdAndTxt(
            "crr_model_mex:Request", "request.%s must be scalar logical.", name);
    }
    return mxIsLogicalScalarTrue(a);
}

bool has_requested_states(const mxArray* a)
{
    if (!a) {
        return false;
    }

    /*
     * request.states is normalized by crr_request.m to a MATLAB string
     * array. For an empty string array (strings(0,1)), the legacy C MEX API
     * can report one object element even though MATLAB's isempty() is true.
     * Therefore mxGetNumberOfElements(a)>0 is not a reliable test here.
     */
    mxArray* rhs = const_cast<mxArray*>(a);
    mxArray* lhs = nullptr;

    if (mexCallMATLAB(1, &lhs, 1, &rhs, "isempty") != 0 || !lhs) {
        mexErrMsgIdAndTxt("crr_model_mex:States",
                          "Could not evaluate isempty(request.states).");
    }

    const bool isEmpty = mxIsLogicalScalarTrue(lhs);
    mxDestroyArray(lhs);

    return !isEmpty;
}

std::vector<double> vec(const mxArray* a, const char* name)
{
    if (!a || mxIsComplex(a) || !(mxIsNumeric(a) || mxIsLogical(a))) {
        mexErrMsgIdAndTxt(
            "crr_model_mex:Vector", "'%s' must be a real numeric/logical array.", name);
    }

    const mwSize n = mxGetNumberOfElements(a);
    std::vector<double> v((size_t)n);

    switch (mxGetClassID(a)) {
    case mxDOUBLE_CLASS: {
        const double* p = mxGetPr(a);
        std::copy(p, p + n, v.begin());
        break;
    }
    case mxSINGLE_CLASS: {
        const float* p = static_cast<const float*>(mxGetData(a));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[i]);
        }
        break;
    }
    case mxLOGICAL_CLASS: {
        const mxLogical* p = mxGetLogicals(a);
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = p[i] ? 1.0 : 0.0;
        }
        break;
    }
    case mxINT8_CLASS: {
        const int8_T* p = static_cast<const int8_T*>(mxGetData(a));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[i]);
        }
        break;
    }
    case mxUINT8_CLASS: {
        const uint8_T* p = static_cast<const uint8_T*>(mxGetData(a));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[i]);
        }
        break;
    }
    case mxINT16_CLASS: {
        const int16_T* p = static_cast<const int16_T*>(mxGetData(a));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[i]);
        }
        break;
    }
    case mxUINT16_CLASS: {
        const uint16_T* p = static_cast<const uint16_T*>(mxGetData(a));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[i]);
        }
        break;
    }
    case mxINT32_CLASS: {
        const int32_T* p = static_cast<const int32_T*>(mxGetData(a));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[i]);
        }
        break;
    }
    case mxUINT32_CLASS: {
        const uint32_T* p = static_cast<const uint32_T*>(mxGetData(a));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[i]);
        }
        break;
    }
    case mxINT64_CLASS: {
        const int64_T* p = static_cast<const int64_T*>(mxGetData(a));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[i]);
        }
        break;
    }
    case mxUINT64_CLASS: {
        const uint64_T* p = static_cast<const uint64_T*>(mxGetData(a));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[i]);
        }
        break;
    }
    default:
        mexErrMsgIdAndTxt(
            "crr_model_mex:Vector", "Unsupported numeric class for '%s'.", name);
    }

    return v;
}

std::vector<mwIndex> expand_idx(const mxArray* a, const char* name)
{
    auto x = vec(a, name);
    if (x.empty()) {
        return {};
    }
    for (double v : x) {
        if (!std::isfinite(v) || std::floor(v) != v || v < 1) {
            mexErrMsgIdAndTxt("crr_model_mex:Index", "Invalid indices in %s.", name);
        }
    }
    std::vector<mwIndex> out;
    if (x.size() == 2) {
        if (x[1] < x[0]) {
            mexErrMsgIdAndTxt("crr_model_mex:Index", "Invalid range %s.", name);
        }
        for (long long i = (long long)x[0]; i <= (long long)x[1]; ++i) {
            out.push_back((mwIndex)(i - 1));
        }
    } else {
        for (size_t k = 0; k < x.size(); ++k) {
            if (k && x[k] <= x[k - 1]) {
                mexErrMsgIdAndTxt("crr_model_mex:Index", "%s must increase.", name);
            }
            out.push_back((mwIndex)((long long)x[k] - 1));
        }
    }
    return out;
}

void setf(mxArray* s, const char* name, mxArray* a)
{
    if (mxGetFieldNumber(s, name) < 0 && mxAddField(s, name) < 0) {
        mexErrMsgIdAndTxt("crr_model_mex:AddField", "Cannot add %s.", name);
    }
    mxSetField(s, 0, name, a);
}

void sets(mxArray* s, const char* name, double x)
{
    setf(s, name, mxCreateDoubleScalar(x));
}

void setv(mxArray* s, const char* name, const std::vector<double>& x)
{
    mxArray* a = mxCreateDoubleMatrix((mwSize)x.size(), 1, mxREAL);
    std::copy(x.begin(), x.end(), mxGetPr(a));
    setf(s, name, a);
}

mxArray* col(const std::vector<double>& x)
{
    mxArray* a = mxCreateDoubleMatrix((mwSize)x.size(), 1, mxREAL);
    std::copy(x.begin(), x.end(), mxGetPr(a));
    return a;
}

mxArray* mat(const std::vector<double>& x, mwSize nr, mwSize nc)
{
    mxArray* a = mxCreateDoubleMatrix(nr, nc, mxREAL);
    std::copy(x.begin(), x.end(), mxGetPr(a));
    return a;
}

struct Transform {
    std::vector<double> th, Jth;
    bool valid = true;
};

Transform transform(const std::vector<double>& x, const mxArray* mdl)
{
    Transform r;
    auto lo = vec(fld(mdl, "th_min"), "mdl.th_min");
    auto hi = vec(fld(mdl, "th_max"), "mdl.th_max");
    if (lo.size() != x.size() || hi.size() != x.size()) {
        mexErrMsgIdAndTxt("crr_model_mex:Bounds", "Bounds mismatch.");
    }
    r.th.resize(x.size());
    r.Jth.resize(x.size(), 1.0);
    int ps = (int)std::llround(scl(fld(mdl, "pspace"), "mdl.pspace"));
    for (size_t j = 0; j < x.size(); ++j) {
        double span = hi[j] - lo[j];
        if (ps == 0) {
            if (x[j] < lo[j] || x[j] > hi[j]) {
                r.valid = false;
                return r;
            }
            r.th[j] = x[j];
        } else if (ps == 1) {
            if (x[j] < 0 || x[j] > 1) {
                r.valid = false;
                return r;
            }
            r.th[j] = lo[j] + x[j] * span;
            r.Jth[j] = span;
        } else if (ps == 2) {
            double z = x[j];
            double n =
                z >= 0 ? 1.0 / (1.0 + std::exp(-z)) : std::exp(z) / (1.0 + std::exp(z));
            r.th[j] = lo[j] + n * span;
            r.Jth[j] = span * n * (1 - n);
        } else {
            mexErrMsgIdAndTxt("crr_model_mex:Pspace", "Unknown pspace.");
        }
    }
    return r;
}

void maxbas(double b,
            std::vector<double>& w,
            std::vector<double>& dw,
            int Lmax = 60,
            double ep = 1e-8)
{
    b = std::max(b, ep);
    double half = .5 * b;
    auto sm = [&](double A, double B) {
        return .5 * (A + B + std::sqrt((A - B) * (A - B) + ep * ep));
    };
    auto dsm = [&](double A, double B) {
        return .5 * (1 + (A - B) / std::sqrt((A - B) * (A - B) + ep * ep));
    };
    double den = sm(half, ep), dden = .5 * dsm(half, ep);
    std::vector<double> a(Lmax), da(Lmax);
    double sa = 0, sda = 0;
    for (int k = 0; k < Lmax; ++k) {
        double z = k - .5 * b, az = std::sqrt(z * z + ep * ep), daz = z / az * (-.5);
        double tri = 1 - az / den, ap = .5 * (tri + std::sqrt(tri * tri + ep * ep));
        double dap = .5 * (1 + tri / std::sqrt(tri * tri + ep * ep));
        double dtri = -daz / den + az / (den * den) * dden;
        a[k] = ap;
        da[k] = dap * dtri;
        sa += a[k];
        sda += da[k];
    }
    w.resize(Lmax);
    dw.resize(Lmax);
    for (int k = 0; k < Lmax; ++k) {
        w[k] = a[k] / sa;
        dw[k] = (da[k] * sa - a[k] * sda) / (sa * sa);
    }
}

std::vector<double> filt(const std::vector<double>& b, const std::vector<double>& x)
{
    std::vector<double> y(x.size(), 0);
    for (size_t i = 0; i < x.size(); ++i) {
        size_t km = std::min(i, b.size() - 1);
        double s = 0;
        for (size_t k = 0; k <= km; ++k) {
            s += b[k] * x[i - k];
        }
        y[i] = s;
    }
    return y;
}

void gr4j_uh(double x4,
             std::vector<double>& U,
             std::vector<double>& dU,
             double mts = 1.0,
             double eta = 2.5)
{
    double p = eta + 1.0;
    int N = std::max(1, (int)std::ceil((2 * x4) / mts));
    U.resize(N);
    dU.resize(N);
    auto G = [&](double t) {
        double u = t / x4;
        if (t <= 0) {
            return 0.0;
        }
        if (t < x4) {
            return std::pow(u, eta + 1.0);
        }
        if (t < 2 * x4) {
            return 1.0 - std::pow(2 - u, eta + 1.0);
        }
        return 1.0;
    };
    auto dG = [&](double t) {
        double u = t / x4;
        if (t > 0 && t < x4) {
            return -(p / x4) * std::pow(u, p);
        }
        if (t >= x4 && t < 2 * x4) {
            return -p * std::pow(2 - u, p - 1.0) * (t / (x4 * x4));
        }
        return 0.0;
    };
    double s = 0, ds = 0;
    for (int k = 0; k < N; ++k) {
        double t = (k + 1) * mts, tm = t - mts;
        double ur = G(t) - G(tm), dur = dG(t) - dG(tm);
        double epsU =
            1e-12 + 1e-6 * std::max(1.0, std::abs(ur)); // local relaxed positivity
        double rt = std::sqrt(ur * ur + epsU * epsU);
        double pos = .5 * (ur + rt), dpos = .5 * (1 + ur / rt);
        U[k] = pos;
        dU[k] = dur * dpos;
        s += U[k];
        ds += dU[k];
    }
    double ss = std::max(s, 1e-30);
    for (int k = 0; k < N; ++k) {
        double raw = U[k], dr = dU[k];
        U[k] = raw / ss;
        dU[k] = (dr * ss - raw * ds) / (ss * ss);
    }
}

typedef void (*Entry)(int, mxArray**, int, const mxArray**);

struct Sim {
    std::vector<double> q, J, Z;
    mwSize nq = 0, nj = 0, zrows = 0, zcols = 0;
    int m = 0, dode = 0;
    bool fail = false;
};

struct InputVectorView {
    const double* ptr = nullptr;
    std::size_t n = 0;
    std::vector<double> buffer;
};

InputVectorView vector_view(const mxArray* a, const char* name)
{
    InputVectorView v;
    if (!a || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("crr_model_mex:Vector", "'%s' must be real.", name);
    }
    v.n = static_cast<std::size_t>(mxGetNumberOfElements(a));
    if (mxIsDouble(a)) {
        v.ptr = mxGetPr(a);
        return v;
    }
    if (mxIsSingle(a)) {
        const float* p = static_cast<const float*>(mxGetData(a));
        v.buffer.resize(v.n);
        for (std::size_t i = 0; i < v.n; ++i) {
            v.buffer[i] = static_cast<double>(p[i]);
        }
        v.ptr = v.buffer.data();
        return v;
    }
    mexErrMsgIdAndTxt("crr_model_mex:Vector", "'%s' must be single/double.", name);
    return v;
}

template <typename Opt> Opt make_core_options(const mxArray* ode)
{
    Opt o;
    o.InitStep = scl(fld(ode, "InitStep"), "ode.InitStep");
    o.MaxStep = scl(fld(ode, "MaxStep"), "ode.MaxStep");
    o.MinStep = scl(fld(ode, "MinStep"), "ode.MinStep");
    o.RelTol = scl(fld(ode, "RelTol"), "ode.RelTol");
    o.AbsTol = scl(fld(ode, "AbsTol"), "ode.AbsTol");
    o.Order = scl(fld(ode, "Order"), "ode.Order");
    o.maxiter = static_cast<int>(std::llround(scl(fld(ode, "maxiter"), "ode.maxiter")));
    return o;
}

std::vector<double>
make_z0(const mxArray* mdl, int m, int dode, bool replicateFirst = false)
{
    auto y0 = vec(fld(mdl, "y0"), "mdl.y0");
    std::vector<double> z0(static_cast<std::size_t>(m) * (dode + 1), 0.0);
    if (replicateFirst) {
        if (y0.empty()) {
            mexErrMsgIdAndTxt("crr_model_mex:Y0", "mdl.y0 is empty.");
        }
        for (int j = 0; j < m; ++j) {
            z0[static_cast<std::size_t>(j)] = y0[0];
        }
    } else {
        if (y0.size() < static_cast<std::size_t>(m)) {
            mexErrMsgIdAndTxt("crr_model_mex:Y0", "mdl.y0 shorter than state count.");
        }
        for (int j = 0; j < m; ++j) {
            z0[static_cast<std::size_t>(j)] = y0[static_cast<std::size_t>(j)];
        }
    }
    return z0;
}

void allocate_core_outputs(Sim& R, int ns, int ipr, bool needJ, bool needStates)
{
    R.zrows = needStates ? static_cast<mwSize>(ns + 1) : 1;
    R.zcols = static_cast<mwSize>(R.m * (R.dode + 1));
    R.Z.assign(static_cast<std::size_t>(R.zrows) * R.zcols, 0.0);
    R.nq = (!needStates && ipr <= ns) ? static_cast<mwSize>(ns - ipr + 1) : 0;
    R.nj = needJ ? static_cast<mwSize>(R.dode) : 0;
    R.q.assign(static_cast<std::size_t>(R.nq), 0.0);
    if (needJ) {
        R.J.assign(static_cast<std::size_t>(R.nq) * R.dode, 0.0);
    }
}

void reconstruct_states_qj(Sim& R, const mxArray* mdl, bool needJ)
{
    auto midx = vec(fld(mdl, "idx"), "mdl.idx");
    const long long i0 = static_cast<long long>(std::llround(midx[0]));
    const long long i1 = static_cast<long long>(std::llround(midx[1]));
    if (i0 < 1 || i1 < i0 || static_cast<mwSize>(i1) > R.zrows) {
        mexErrMsgIdAndTxt("crr_model_mex:States",
                          "mdl.idx is incompatible with state history.");
    }

    R.nq = static_cast<mwSize>(i1 - i0);
    R.q.assign(static_cast<std::size_t>(R.nq), 0.0);
    const std::size_t qcol = static_cast<std::size_t>(R.m - 1);
    for (mwSize i = 0; i < R.nq; ++i) {
        const std::size_t ra = static_cast<std::size_t>(i0 - 1) + i;
        const std::size_t rb = ra + 1;
        R.q[i] = R.Z[rb + static_cast<std::size_t>(R.zrows) * qcol] -
                 R.Z[ra + static_cast<std::size_t>(R.zrows) * qcol];
    }

    if (needJ) {
        R.J.assign(static_cast<std::size_t>(R.nq) * R.dode, 0.0);
        for (int j = 0; j < R.dode; ++j) {
            const std::size_t scol = (static_cast<std::size_t>(j) + 2) * R.m - 1;
            if (scol >= R.zcols) {
                mexErrMsgIdAndTxt("crr_model_mex:J", "Sensitivity column exceeds Z.");
            }
            for (mwSize i = 0; i < R.nq; ++i) {
                const std::size_t ra = static_cast<std::size_t>(i0 - 1) + i;
                const std::size_t rb = ra + 1;
                R.J[i + static_cast<std::size_t>(R.nq) * j] =
                    R.Z[rb + static_cast<std::size_t>(R.zrows) * scol] -
                    R.Z[ra + static_cast<std::size_t>(R.zrows) * scol];
            }
        }
        R.nj = static_cast<mwSize>(R.dode);
    }
}

void scale_jacobian(Sim& R, const std::vector<double>& Jth, bool needJ)
{
    if (!needJ) {
        return;
    }
    if (R.nj != Jth.size()) {
        mexErrMsgIdAndTxt("crr_model_mex:J", "Jacobian columns mismatch.");
    }
    for (std::size_t j = 0; j < Jth.size(); ++j) {
        for (mwSize i = 0; i < R.nq; ++i) {
            R.J[i + static_cast<std::size_t>(R.nq) * j] *= Jth[j];
        }
    }
}

Sim run_builtin(int model,
                const std::vector<double>& th,
                const std::vector<double>& Jth,
                const mxArray* mdl,
                const mxArray* meteo,
                const mxArray* ode,
                bool needJ,
                bool needStates)
{
    auto midx = vec(fld(mdl, "idx"), "mdl.idx");
    if (midx.size() < 2) {
        mexErrMsgIdAndTxt("crr_model_mex:Idx", "mdl.idx needs 2 entries.");
    }
    const int ns = static_cast<int>(std::llround(scl(fld(mdl, "tout"), "mdl.tout")));
    int ipr = needStates ? 1 : static_cast<int>(std::llround(midx[0]));
    ipr = std::max(1, std::min(ipr, ns + 1));

    InputVectorView P = vector_view(fld(meteo, "P"), "meteo.P");
    InputVectorView Ep = vector_view(fld(meteo, "Ep"), "meteo.Ep");
    InputVectorView T = vector_view(fld(meteo, "T"), "meteo.T");
    const std::size_t nf = std::min(P.n, std::min(Ep.n, T.n));
    Sim R;

    if (model == 1) {
        if (th.size() != 7) {
            mexErrMsgIdAndTxt("crr_model_mex:Pars", "HYMOD needs 7 parameters.");
        }
        R.m = 7;
        R.dode = 7;
        auto z0 = make_z0(mdl, R.m, R.dode);
        allocate_core_outputs(R, ns, ipr, needJ, needStates);
        sage_hymod::Params p{};
        p.S_umax = th[0];
        p.beta = th[1];
        p.alfa = th[2];
        p.K_s = th[3];
        p.K_f = th[4];
        p.T_tr = th[5];
        p.f_dd = th[6];
        p.T_sm = 1;
        p.eps_m = 1e-6;
        p.rho = .01;
        auto opt = make_core_options<sage_hymod::Options>(ode);
        sage_hymod::Forcing F{P.ptr, Ep.ptr, T.ptr, nf};
        sage_hymod::OutputView O{R.Z.data(),
                                 R.q.empty() ? nullptr : R.q.data(),
                                 R.J.empty() ? nullptr : R.J.data(),
                                 R.zrows,
                                 R.zcols,
                                 R.nq,
                                 R.nj};
        R.fail = sage_hymod::run_into(
            ns, z0.data(), z0.size(), F, p, opt, needStates, ipr, needJ, O);
    } else if (model == 2) {
        if (th.size() != 9) {
            mexErrMsgIdAndTxt("crr_model_mex:Pars", "HMODEL needs 9 parameters.");
        }
        R.m = 6;
        R.dode = 9;
        auto z0 = make_z0(mdl, R.m, R.dode);
        allocate_core_outputs(R, ns, ipr, needJ, needStates);
        sage_hmodel::Params p{};
        p.I_max = th[0];
        p.Su_max = th[1];
        p.Q_max = th[2];
        p.a_E = th[3];
        p.a_F = th[4];
        p.a_S = 1e-6;
        p.r_f = th[5];
        p.r_s = th[6];
        p.T_tr = th[7];
        p.f_dd = th[8];
        p.T_sm = 1;
        p.eps_m = 1e-6;
        p.rho = .01;
        auto opt = make_core_options<sage_hmodel::Options>(ode);
        sage_hmodel::Forcing F{P.ptr, Ep.ptr, T.ptr, nf};
        sage_hmodel::OutputView O{R.Z.data(),
                                  R.q.empty() ? nullptr : R.q.data(),
                                  R.J.empty() ? nullptr : R.J.data(),
                                  R.zrows,
                                  R.zcols,
                                  R.nq,
                                  R.nj};
        R.fail = sage_hmodel::run_into(
            ns, z0.data(), z0.size(), F, p, opt, needStates, ipr, needJ, O);
    } else if (model == 3) {
        if (th.size() != 15) {
            mexErrMsgIdAndTxt("crr_model_mex:Pars", "SAC-SMA needs 15 parameters.");
        }
        R.m = 10;
        R.dode = 15;
        auto z0 = make_z0(mdl, R.m, R.dode);
        allocate_core_outputs(R, ns, ipr, needJ, needStates);
        sage_sacsma::Params p{};
        p.uzfwm = th[0];
        p.uztwm = th[1];
        p.lzfpm = th[2];
        p.lzfsm = th[3];
        p.lztwm = th[4];
        p.zperc = th[5];
        p.rexp = th[6];
        p.uzk = th[7];
        p.pfree = th[8];
        p.lzpk = th[9];
        p.lzsk = th[10];
        p.acm = th[11];
        p.kf = th[12];
        p.T_tr = th[13];
        p.f_dd = th[14];
        p.T_sm = 1;
        p.eps_m = 1e-6;
        p.eps_s = 1e-12;
        p.eps = 5;
        p.rho = .01;
        auto opt = make_core_options<sage_sacsma::Options>(ode);
        sage_sacsma::Forcing F{P.ptr, Ep.ptr, T.ptr, nf};
        sage_sacsma::OutputView O{R.Z.data(),
                                  R.q.empty() ? nullptr : R.q.data(),
                                  R.J.empty() ? nullptr : R.J.data(),
                                  R.zrows,
                                  R.zcols,
                                  R.nq,
                                  R.nj};
        R.fail = sage_sacsma::run_into(
            ns, z0.data(), z0.size(), F, p, opt, needStates, ipr, needJ, O);
    } else if (model == 4) {
        if (th.size() != 16) {
            mexErrMsgIdAndTxt("crr_model_mex:Pars", "Xinanjiang needs 16 parameters.");
        }
        R.m = 9;
        R.dode = 16;
        auto z0 = make_z0(mdl, R.m, R.dode);
        allocate_core_outputs(R, ns, ipr, needJ, needStates);
        double fwm = th[4], flm = th[5], Stot = th[7], W = fwm * Stot, S = (1 - fwm) * Stot;
        sage_xinanjiang::Params p{};
        p.f_p = th[0];
        p.A_im = th[1];
        p.a = th[2];
        p.b = th[3];
        p.W_max = W;
        p.LM = flm * W;
        p.c = th[6];
        p.S_max = S;
        p.S_tot = Stot;
        p.f_wm = fwm;
        p.f_lm = flm;
        p.Ex = th[8];
        p.k_i = th[9];
        p.k_g = th[10];
        p.c_i = th[11];
        p.c_g = th[12];
        p.k_f = th[13];
        p.T_tr = th[14];
        p.f_dd = th[15];
        p.T_sm = 1;
        p.eps_m = 1e-6;
        p.eps_r = 1e-12;
        p.eps = 5;
        p.rho = .01;
        auto opt = make_core_options<sage_xinanjiang::Options>(ode);
        sage_xinanjiang::Forcing F{P.ptr, Ep.ptr, T.ptr, nf};
        sage_xinanjiang::OutputView O{R.Z.data(),
                                      R.q.empty() ? nullptr : R.q.data(),
                                      R.J.empty() ? nullptr : R.J.data(),
                                      R.zrows,
                                      R.zcols,
                                      R.nq,
                                      R.nj};
        R.fail = sage_xinanjiang::run_into(
            ns, z0.data(), z0.size(), F, p, opt, needStates, ipr, needJ, O);
    } else if (model == 5) {
        if (th.size() != 8) {
            mexErrMsgIdAndTxt("crr_model_mex:Pars", "GR4J-A needs 8 parameters.");
        }
        int n1 = static_cast<int>(std::llround(scl(fld(mdl, "n1"), "mdl.n1"))),
            n2 = static_cast<int>(std::llround(scl(fld(mdl, "n2"), "mdl.n2")));
        R.m = n1 + n2 + 4;
        R.dode = 8;
        auto z0 = make_z0(mdl, R.m, R.dode);
        allocate_core_outputs(R, ns, ipr, needJ, needStates);
        sage_gr4ja::Params p{};
        p.x1 = th[0];
        p.x2 = th[1];
        p.x3 = th[2];
        p.x4 = th[3];
        p.x5 = th[4];
        p.f_p = th[5];
        p.T_tr = th[6];
        p.f_dd = th[7];
        p.kappa = 4.0 / 9.0;
        p.bg = 3.5;
        p.bR = 5;
        p.mts = 1;
        p.T_sm = 1;
        p.eps_m = 1e-6;
        p.eps_s = 1e-12;
        p.rho = .01;
        p.n1 = n1;
        p.n2 = n2;
        auto opt = make_core_options<sage_gr4ja::Options>(ode);
        sage_gr4ja::Forcing F{P.ptr, Ep.ptr, T.ptr, nf};
        sage_gr4ja::OutputView O{R.Z.data(),
                                 R.q.empty() ? nullptr : R.q.data(),
                                 R.J.empty() ? nullptr : R.J.data(),
                                 R.zrows,
                                 R.zcols,
                                 R.nq,
                                 R.nj};
        R.fail = sage_gr4ja::run_into(
            ns, z0.data(), z0.size(), F, p, opt, needStates, ipr, needJ, O);
    } else if (model == 21) {
        if (th.size() != 8) {
            mexErrMsgIdAndTxt("crr_model_mex:Pars", "GR4J-B needs 8 parameters.");
        }
        const double eta = 2.5, mts = 1.0, x4 = th[3], pexp = eta + 1.0;
        const int L = std::max(1, static_cast<int>(std::ceil(2.0 * x4 / mts)));
        std::vector<double> U(L), dU(L);
        auto cif = [=](double t) {
            if (t <= 0.0) return 0.0;
            if (t < x4) return std::pow(t / x4, pexp);
            if (t < 2.0 * x4) return 1.0 - std::pow(2.0 - t / x4, pexp);
            return 1.0;
        };
        auto dcif = [=](double t) {
            if (t > 0.0 && t < x4)
                return -(pexp / x4) * std::pow(t / x4, pexp);
            if (t >= x4 && t < 2.0 * x4)
                return -pexp * std::pow(2.0 - t / x4, pexp - 1.0) * t / (x4 * x4);
            return 0.0;
        };
        double maxRaw = 0.0;
        for (int k = 0; k < L; ++k) {
            U[k] = cif((k + 1) * mts) - cif(k * mts);
            dU[k] = dcif((k + 1) * mts) - dcif(k * mts);
            maxRaw = std::max(maxRaw, std::abs(U[k]));
        }
        const double epsU = 1e-12 + 1e-6 * std::max(1.0, maxRaw);
        double sumU = 0.0, sumDU = 0.0;
        for (int k = 0; k < L; ++k) {
            const double raw = U[k], den = std::sqrt(raw * raw + epsU * epsU);
            U[k] = 0.5 * (raw + den);
            dU[k] *= 0.5 * (1.0 + raw / den);
            sumU += U[k]; sumDU += dU[k];
        }
        const double safe = std::max(sumU, 1e-30);
        for (int k = 0; k < L; ++k) {
            dU[k] = (dU[k] * safe - U[k] * sumDU) / (safe * safe);
            U[k] /= safe;
        }
        R.m = 4 + 2 * L;
        R.dode = 8;
        auto z0 = make_z0(mdl, R.m, R.dode);
        allocate_core_outputs(R, ns, ipr, needJ, needStates);
        sage_gr4jb::Params p{};
        p.x1=th[0]; p.x2=th[1]; p.x3=th[2]; p.x4=th[3]; p.x5=th[4];
        p.f_p=th[5]; p.T_tr=th[6]; p.f_dd=th[7]; p.eta=eta; p.tau=eta;
        p.kappa=4.0/9.0; p.bg=3.5; p.bR=5.0; p.mts=mts; p.T_sm=1.0;
        p.eps_m=1e-6; p.eps_s=1e-12; p.rho=.01; p.L=L;
        p.U=U.data(); p.dUdx4=dU.data();
        auto opt = make_core_options<sage_gr4jb::Options>(ode);
        sage_gr4jb::Forcing F{P.ptr, Ep.ptr, T.ptr, nf};
        sage_gr4jb::OutputView O{R.Z.data(), R.q.empty()?nullptr:R.q.data(),
            R.J.empty()?nullptr:R.J.data(), R.zrows, R.zcols, R.nq, R.nj};
        R.fail = sage_gr4jb::run_into(ns,z0.data(),z0.size(),F,p,opt,
            needStates,ipr,needJ,O);
    } else if (model == 6) {
        if (!(th.size() == 12 || th.size() == 13)) {
            mexErrMsgIdAndTxt("crr_model_mex:Pars", "HBV needs 12/13 parameters.");
        }
        R.m = 5;
        R.dode = 12;
        auto z0 = make_z0(mdl, R.m, R.dode);
        allocate_core_outputs(R, ns, ipr, needJ, needStates);
        sage_hbv::Params p{};
        p.f_c = th[0];
        p.beta = th[1];
        p.lp = th[2];
        p.k_0 = th[3];
        p.uzl = th[4];
        p.k_1 = th[5];
        p.k_2 = th[6];
        p.perc = th[7];
        p.T_tr = th[8];
        p.f_dd = th[9];
        p.sfcf = th[10];
        p.cfr = th[11];
        p.eps_s = 1e-3;
        p.eps_t = .1;
        p.eps_x = 1e-12;
        auto opt = make_core_options<sage_hbv::Options>(ode);
        sage_hbv::Forcing F{P.ptr, Ep.ptr, T.ptr, nf};
        sage_hbv::OutputView O{R.Z.data(),
                               R.q.empty() ? nullptr : R.q.data(),
                               R.J.empty() ? nullptr : R.J.data(),
                               R.zrows,
                               R.zcols,
                               R.nq,
                               R.nj};
        R.fail = sage_hbv::run_into(
            ns, z0.data(), z0.size(), F, p, opt, needStates, ipr, needJ, O);
    } else if (model == 7) {
        if (th.size() != 15) {
            mexErrMsgIdAndTxt("crr_model_mex:Pars", "CFE-NWM needs 15 parameters.");
        }
        auto giuh = vec(fld(mdl, "giuh_ordnts"), "mdl.giuh_ordnts");
        int L = static_cast<int>(giuh.size()), K = 3;
        R.m = 3 + L + K + 1;
        R.dode = 15;
        auto z0 = make_z0(mdl, R.m, R.dode);
        allocate_core_outputs(R, ns, ipr, needJ, needStates);
        double smax = th[0], sfc = th[1] * smax, swp = th[2] * smax;
        sage_cfe_nwm::Params p{};
        p.s_max = smax;
        p.s_fc = sfc;
        p.s_wp = swp;
        p.k_sch = th[3];
        p.a1 = th[4];
        p.k_perc = th[5];
        p.lf_thr = swp + th[6] * (sfc - swp);
        p.a2 = th[7];
        p.k_lf = th[8];
        p.g_max = th[9];
        p.c_gw = th[10];
        p.mm = th[11];
        p.k_nsh = th[12];
        p.T_tr = th[13];
        p.f_dd = th[14];
        p.dT = 1;
        p.T_sm = 1;
        p.eps_m = 1e-6;
        p.eps = 5;
        p.rho = .01;
        p.K = K;
        p.L = L;
        p.giuh = giuh.data();
        auto opt = make_core_options<sage_cfe_nwm::Options>(ode);
        sage_cfe_nwm::Forcing F{P.ptr, Ep.ptr, T.ptr, nf};
        sage_cfe_nwm::OutputView O{R.Z.data(),
                                   R.q.empty() ? nullptr : R.q.data(),
                                   R.J.empty() ? nullptr : R.J.data(),
                                   R.zrows,
                                   R.zcols,
                                   R.nq,
                                   R.nj};
        R.fail = sage_cfe_nwm::run_into(
            ns, z0.data(), z0.size(), F, p, opt, needStates, ipr, needJ, O);
    } else {
        mexErrMsgIdAndTxt("crr_model_mex:Model", "Unsupported built-in model %d.", model);
    }

    if (needStates) {
        reconstruct_states_qj(R, mdl, needJ);
    }

    /*
     * HBV's optional 13th parameter is the external MAXBAS routing
     * parameter. The HBV core itself has 12 sensitivity columns.
     *
     * Apply MAXBAS only after q/J have been produced. This ordering is
     * essential when states are requested because q/J are reconstructed
     * from the stored augmented state history first.
     */
    if (model == 6 && th.size() == 13) {
        std::vector<double> w, dw;
        maxbas(th[12], w, dw);

        const auto q0 = R.q;
        R.q = filt(w, q0);

        if (needJ) {
            if (R.nj != 12) {
                mexErrMsgIdAndTxt(
                    "crr_model_mex:J",
                    "HBV core must provide 12 Jacobian columns before MAXBAS.");
            }

            std::vector<double> Jr(static_cast<std::size_t>(R.nq) * 13, 0.0);

            for (int j = 0; j < 12; ++j) {
                std::vector<double> c(R.nq);

                for (mwSize i = 0; i < R.nq; ++i) {
                    c[i] = R.J[i + static_cast<std::size_t>(R.nq) * j];
                }

                c = filt(w, c);

                for (mwSize i = 0; i < R.nq; ++i) {
                    Jr[i + static_cast<std::size_t>(R.nq) * j] = c[i];
                }
            }

            const auto dc = filt(dw, q0);

            for (mwSize i = 0; i < R.nq; ++i) {
                Jr[i + static_cast<std::size_t>(R.nq) * 12] = dc[i];
            }

            R.J.swap(Jr);
            R.nj = 13;
        }
    }

    scale_jacobian(R, Jth, needJ);
    return R;
}

struct Kge {
    double K = QNAN, r = QNAN, a = QNAN, b = QNAN;
};

Kge kge(const std::vector<double>& y, const std::vector<double>& q, double muy, double sy)
{
    Kge z;
    std::vector<double> yy, qq;
    for (size_t i = 0; i < y.size(); ++i) {
        if (std::isfinite(y[i]) && std::isfinite(q[i])) {
            yy.push_back(y[i]);
            qq.push_back(q[i]);
        }
    }
    size_t n = yy.size();
    if (n < 2) {
        return z;
    }
    double mq = std::accumulate(qq.begin(), qq.end(), 0.0) / n, ss = 0;
    for (double v : qq) {
        ss += (v - mq) * (v - mq);
    }
    double sq = std::sqrt(ss / (n - 1));
    if (std::isfinite(sy) && sy > 0) {
        z.a = sq / sy;
    }
    if (std::isfinite(muy) && muy != 0) {
        z.b = mq / muy;
    }
    if (std::isfinite(sq) && sq > 0 && std::isfinite(sy) && sy > 0) {
        double c = 0;
        for (size_t i = 0; i < n; ++i) {
            c += (yy[i] - muy) * (qq[i] - mq);
        }
        z.r = c / ((n - 1) * sy * sq);
        z.r = std::max(-1.0, std::min(1.0, z.r));
    }
    if (std::isfinite(z.r) && std::isfinite(z.a) && std::isfinite(z.b)) {
        z.K = 1 - std::sqrt((z.r - 1) * (z.r - 1) + (z.a - 1) * (z.a - 1) +
                            (z.b - 1) * (z.b - 1));
    }
    return z;
}

double huber(const std::vector<double>& r, double S, std::vector<double>* d = nullptr)
{
    if (d) {
        d->assign(r.size(), 0);
    }
    if (!(S > 0) || !std::isfinite(S)) {
        return QNAN;
    }
    double L = 0;
    bool any = false;
    for (size_t i = 0; i < r.size(); ++i) {
        if (!std::isfinite(r[i])) {
            continue;
        }
        any = true;
        double u = r[i] / S, a = std::abs(u), v = std::min(a, 1.345);
        L += v * (a - .5 * v);
        if (d) {
            (*d)[i] = -(u > 0 ? 1 : (u < 0 ? -1 : 0)) * v / S;
        }
    }
    return any ? L : QNAN;
}

void delta(int lf,
           const std::vector<double>& y,
           const std::vector<double>& q,
           std::vector<double>& d)
{
    size_t n = y.size();
    d.assign(n, 0);
    std::vector<double> e(n);
    for (size_t i = 0; i < n; ++i) {
        e[i] = y[i] - q[i];
    }
    if (lf == 1) {
        for (size_t i = 0; i < n; ++i) {
            d[i] = -(e[i] > 0 ? 1 : (e[i] < 0 ? -1 : 0));
        }
        return;
    }
    if (lf == 2) {
        for (size_t i = 0; i < n; ++i) {
            d[i] = -2 * e[i];
        }
        return;
    }
    if (lf == 3) {
        double my = std::accumulate(y.begin(), y.end(), 0.0) / n, ss = 0;
        for (double v : y) {
            ss += (v - my) * (v - my);
        }
        for (size_t i = 0; i < n; ++i) {
            d[i] = -2 * e[i] / ss;
        }
        return;
    }
    if (lf == 4) {
        double my = std::accumulate(y.begin(), y.end(), 0.0) / n,
               mq = std::accumulate(q.begin(), q.end(), 0.0) / n, sy2 = 0, sq2 = 0, c = 0;
        for (size_t i = 0; i < n; ++i) {
            sy2 += (y[i] - my) * (y[i] - my);
            sq2 += (q[i] - mq) * (q[i] - mq);
            c += (y[i] - my) * (q[i] - mq);
        }
        double sy = std::sqrt(sy2 / (n - 1)), sq = std::sqrt(sq2 / (n - 1)),
               r = c / ((n - 1) * sy * sq), v = sq / sy, z = mq / my,
               L = std::sqrt((r - 1) * (r - 1) + (v - 1) * (v - 1) + (z - 1) * (z - 1));
        for (size_t i = 0; i < n; ++i) {
            double dm = 1.0 / n, ds = (q[i] - mq) / ((n - 1) * sq),
                   dr = (y[i] - my) / ((n - 1) * sq * sy) -
                        r * (q[i] - mq) / ((n - 1) * sq * sq);
            d[i] = ((r - 1) * dr + (v - 1) * ds / sy + (z - 1) * dm / my) / L;
        }
        return;
    }
    mexErrMsgIdAndTxt("crr_model_mex:Delta", "Unsupported native loss delta.");
}

// ----------------------------
// FDC objective and derivative
// ----------------------------
double fdc_loss_native(const std::vector<double>& q, const mxArray* fdc)
{
    const size_t n = (size_t)std::llround(scl(fld(fdc, "n"), "fdc.n"));
    if (n == 0 || q.empty()) {
        return QNAN;
    }
    if (q.size() != n) {
        mexErrMsgIdAndTxt("crr_model_mex:FDC", "FDC cache length mismatch.");
    }
    auto ys = vec(fld(fdc, "ys"), "fdc.ys");
    auto Py = vec(fld(fdc, "Py"), "fdc.Py");
    double Syy = scl(fld(fdc, "S_yy"), "fdc.S_yy");
    if (ys.size() != n || Py.size() != n + 1) {
        mexErrMsgIdAndTxt("crr_model_mex:FDC", "Malformed FDC cache.");
    }
    std::vector<double> qs = q;
    std::sort(qs.begin(), qs.end());
    double Sqy = 0.0, Sqq = 0.0;
    size_t j = 0;
    const double sumY = Py.back();
    for (size_t i = 0; i < n; ++i) {
        const double qi = qs[i];
        while (j < n && ys[j] <= qi) {
            ++j;
        }
        const double left = (double)j * qi - Py[j];
        const double right = (sumY - Py[j]) - (double)(n - j) * qi;
        Sqy += left + right;
        const double w = 2.0 * (double)(i + 1) - (double)n - 1.0;
        Sqq += 2.0 * w * qi;
    }
    const double nn = (double)n * (double)n;
    return (Sqy - 0.5 * (Sqq + Syy)) / nn;
}

void fdc_delta_native(const std::vector<double>& y,
                      const std::vector<double>& q,
                      std::vector<double>& d)
{
    const size_t n = q.size();
    d.assign(n, 0.0);
    if (n == 0) {
        return;
    }
    std::vector<double> ys = y, qs = q;
    std::sort(ys.begin(), ys.end());
    std::sort(qs.begin(), qs.end());
    const double nn = (double)n * (double)n;
    for (size_t i = 0; i < n; ++i) {
        const double x = q[i];
        auto yl = std::lower_bound(ys.begin(), ys.end(), x);
        auto yu = std::upper_bound(ys.begin(), ys.end(), x);
        auto ql = std::lower_bound(qs.begin(), qs.end(), x);
        auto qu = std::upper_bound(qs.begin(), qs.end(), x);
        const long long ltY = yl - ys.begin(), gtY = ys.end() - yu;
        const long long ltQ = ql - qs.begin(), gtQ = qs.end() - qu;
        d[i] = ((double)(ltY - gtY) - (double)(ltQ - gtQ)) / nn;
    }
}

// ----------------------------
// JKGE matrix-free machinery
// ----------------------------
struct JCache {
    int method = 0;
    size_t n = 0, G = 0;
    std::vector<double> L, R, den, gid, cnt;
    std::vector<unsigned char> good;
};

JCache jkcache(const mxArray* c, int method, size_t n)
{
    JCache z;
    z.method = method;
    z.n = n;
    if (!c || !mxIsStruct(c)) {
        mexErrMsgIdAndTxt("crr_model_mex:JKGE", "dat.jkge.cache missing.");
    }
    const mxArray* ga = fld(c, "good", false);
    if (ga) {
        auto g = vec(ga, "jkge.cache.good");
        z.good.resize(g.size());
        for (size_t i = 0; i < g.size(); ++i) {
            z.good[i] = (unsigned char)(g[i] != 0);
        }
    } else {
        z.good.assign(n, 1);
    }
    if (method == 1) {
        z.L = vec(fld(c, "L"), "jkge.cache.L");
        z.R = vec(fld(c, "R"), "jkge.cache.R");
        z.den = vec(fld(c, "den"), "jkge.cache.den");
    }
    if (method == 2 || method == 4) {
        z.gid = vec(fld(c, "gid"), "jkge.cache.gid");
        z.G = (size_t)std::llround(scl(fld(c, "G"), "jkge.cache.G"));
        z.cnt = vec(fld(c, "cnt_g"), "jkge.cache.cnt_g");
    }
    return z;
}

std::vector<double> jk_benchmark(const std::vector<double>& q,
                                 int method,
                                 int nwin,
                                 const std::vector<double>& mo,
                                 const JCache& c)
{
    const size_t n = q.size();
    std::vector<double> m(n, QNAN);
    if (method == 1) {
        std::vector<double> cs(n + 1, 0.0);
        for (size_t i = 0; i < n; ++i) {
            cs[i + 1] = cs[i] + (std::isfinite(q[i]) ? q[i] : 0.0);
        }
        for (size_t i = 0; i < n; ++i) {
            size_t L = (size_t)c.L[i], R = (size_t)c.R[i];
            double den = c.den[i];
            if (den > 0) {
                m[i] = (cs[R] - cs[L - 1]) / den;
            }
        }
    } else if (method == 2) {
        std::vector<double> s(c.G, 0.0);
        for (size_t i = 0; i < n; ++i) {
            if (std::isfinite(q[i])) {
                size_t g = (size_t)c.gid[i] - 1;
                if (g < c.G) {
                    s[g] += q[i];
                }
            }
        }
        for (size_t i = 0; i < n; ++i) {
            size_t g = (size_t)c.gid[i] - 1;
            if (g < c.G && c.cnt[g] > 0) {
                m[i] = s[g] / c.cnt[g];
            }
        }
    } else if (method == 3) {
        double s = 0;
        size_t ng = 0;
        for (double v : q) {
            if (std::isfinite(v)) {
                s += v;
                ++ng;
            }
        }
        if (ng) {
            double mu = s / ng;
            std::fill(m.begin(), m.end(), mu);
        }
    } else if (method == 4) {
        std::vector<double> s(c.G, 0.0);
        for (size_t i = 0; i < n; ++i) {
            if (i < c.good.size() && c.good[i] && std::isfinite(q[i])) {
                size_t g = (size_t)c.gid[i] - 1;
                if (g < c.G) {
                    s[g] += q[i];
                }
            }
        }
        for (size_t i = 0; i < n; ++i) {
            if (i < c.good.size() && c.good[i]) {
                size_t g = (size_t)c.gid[i] - 1;
                if (g < c.G && c.cnt[g] > 0) {
                    m[i] = s[g] / c.cnt[g];
                }
            }
        }
    } else {
        mexErrMsgIdAndTxt("crr_model_mex:JKGE", "Unknown JKGE method.");
    }
    return m;
}

std::vector<double> jk_transpose(const std::vector<double>& v,
                                 const std::vector<double>& q,
                                 int method,
                                 const JCache& c)
{
    size_t n = q.size();
    std::vector<double> z(n, 0.0);
    if (method == 1) {
        std::vector<double> acc(n + 1, 0.0);
        for (size_t i = 0; i < n; ++i) {
            double den = c.den[i];
            if (std::isfinite(v[i]) && den > 0) {
                double w = v[i] / den;
                size_t L = (size_t)c.L[i], R = (size_t)c.R[i];
                acc[L - 1] += w;
                if (R < n) {
                    acc[R] -= w;
                }
            }
        }
        double s = 0;
        for (size_t i = 0; i < n; ++i) {
            s += acc[i];
            z[i] = (i < c.good.size() && c.good[i]) ? s : 0.0;
        }
    } else if (method == 2 || method == 4) {
        std::vector<double> s(c.G, 0.0);
        for (size_t i = 0; i < n; ++i) {
            if (i < c.good.size() && c.good[i] && std::isfinite(v[i])) {
                size_t g = (size_t)c.gid[i] - 1;
                if (g < c.G) {
                    s[g] += v[i];
                }
            }
        }
        for (size_t i = 0; i < n; ++i) {
            if (i < c.good.size() && c.good[i]) {
                size_t g = (size_t)c.gid[i] - 1;
                if (g < c.G && c.cnt[g] > 0) {
                    z[i] = s[g] / c.cnt[g];
                }
            }
        }
    } else if (method == 3) {
        size_t ng = 0;
        double s = 0;
        for (size_t i = 0; i < n; ++i) {
            if (i < c.good.size() && c.good[i]) {
                ++ng;
                if (std::isfinite(v[i])) {
                    s += v[i];
                }
            }
        }
        if (ng) {
            for (size_t i = 0; i < n; ++i) {
                if (i < c.good.size() && c.good[i]) {
                    z[i] = s / ng;
                }
            }
        }
    }
    return z;
}

struct JKResult {
    double JK = QNAN, M = QNAN, V = QNAN, C = QNAN;
    std::vector<double> mq, delta;
};

JKResult jk_train(const std::vector<double>& y,
                  const std::vector<double>& q,
                  const std::vector<double>& my,
                  const std::vector<mwIndex>& idx,
                  int method,
                  int nwin,
                  const std::vector<double>& mo,
                  const JCache& cache,
                  int Mdef,
                  bool needDelta)
{
    JKResult R;
    size_t N = q.size();
    R.mq = jk_benchmark(q, method, nwin, mo, cache);
    if (needDelta) {
        R.delta.assign(N, 0.0);
    }
    std::vector<unsigned char> good(N, 0);
    for (auto i : idx) {
        if (i < N) {
            good[i] = 1;
        }
    }
    std::vector<size_t> ii;
    for (size_t i = 0; i < N; ++i) {
        if (good[i] && std::isfinite(y[i]) && std::isfinite(q[i]) && std::isfinite(my[i]) &&
            std::isfinite(R.mq[i])) {
            ii.push_back(i);
        }
    }
    size_t n = ii.size();
    if (n < 2) {
        return R;
    }
    std::vector<double> ay(n), aq(n), myr(n), mqr(n), yr(n), qr(n);
    for (size_t k = 0; k < n; ++k) {
        size_t i = ii[k];
        yr[k] = y[i];
        qr[k] = q[i];
        myr[k] = my[i];
        mqr[k] = R.mq[i];
        ay[k] = yr[k] - myr[k];
        aq[k] = qr[k] - mqr[k];
    }
    double sy = 0, sq = 0;
    for (size_t k = 0; k < n; ++k) {
        sy += ay[k] * ay[k];
        sq += aq[k] * aq[k];
    }
    sy = std::sqrt(sy / n);
    sq = std::sqrt(sq / n);
    if (!(sy > 0 && sq > 0 && std::isfinite(sy) && std::isfinite(sq))) {
        return R;
    }
    double rho = 0;
    for (size_t k = 0; k < n; ++k) {
        rho += (ay[k] / sy) * (aq[k] / sq);
    }
    rho /= n;
    rho = std::max(-1.0, std::min(1.0, rho));
    double alpha = sq / sy;
    R.V = (1 - alpha) * (1 - alpha);
    R.C = (1 - rho) * (1 - rho);
    std::vector<unsigned char> gm(n, 1);
    double denM = QNAN;
    if (Mdef == 1) {
        double ma = 0;
        for (double v : yr) {
            ma += std::abs(v);
        }
        ma /= n;
        double tol = std::max(1e-8, 1e-6 * ma);
        double s = 0;
        size_t ng = 0;
        for (size_t k = 0; k < n; ++k) {
            gm[k] = (std::abs(myr[k]) > tol);
            if (gm[k]) {
                double e = 1 - mqr[k] / myr[k];
                s += e * e;
                ++ng;
            }
        }
        if (!ng) {
            return R;
        }
        R.M = s / ng;
    } else {
        double ybar = std::accumulate(yr.begin(), yr.end(), 0.0) / n;
        double den2 = 0, num2 = 0;
        for (size_t k = 0; k < n; ++k) {
            den2 += (myr[k] - ybar) * (myr[k] - ybar);
            num2 += (mqr[k] - myr[k]) * (mqr[k] - myr[k]);
        }
        denM = std::sqrt(den2);
        if (!(denM > 0 && std::isfinite(denM))) {
            return R;
        }
        R.M = num2 / (denM * denM);
    }
    double L = std::sqrt(R.M + R.V + R.C);
    if (!(L > 0 && std::isfinite(L))) {
        return R;
    }
    R.JK = 1 - L;
    if (!needDelta) {
        return R;
    }
    std::vector<double> dMmq(n, 0), dV(n, 0), dC(n, 0), dmq(n, 0), dq(n, 0);
    if (Mdef == 1) {
        size_t ng = 0;
        for (auto b : gm) {
            if (b) {
                ++ng;
            }
        }
        for (size_t k = 0; k < n; ++k) {
            if (gm[k]) {
                dMmq[k] = -2 * (1 - mqr[k] / myr[k]) / (ng * myr[k]);
            }
        }
    } else {
        for (size_t k = 0; k < n; ++k) {
            dMmq[k] = 2 * (mqr[k] - myr[k]) / (denM * denM);
        }
    }
    double S = 0;
    for (size_t k = 0; k < n; ++k) {
        S += (ay[k] / sy) * aq[k];
    }
    for (size_t k = 0; k < n; ++k) {
        dV[k] = -2 * (1 - alpha) * aq[k] / (n * sq * sy);
        double dr = (ay[k] / sy) / (n * sq) - S * aq[k] / ((double)n * n * sq * sq * sq);
        dC[k] = -2 * (1 - rho) * dr;
        dq[k] = (dV[k] + dC[k]) / (2 * L);
        dmq[k] = (dMmq[k] - dV[k] - dC[k]) / (2 * L);
        R.delta[ii[k]] += dq[k];
    }
    std::vector<double> vf(N, 0);
    for (size_t k = 0; k < n; ++k) {
        vf[ii[k]] = dmq[k];
    }
    auto bt = jk_transpose(vf, q, method, cache);
    for (size_t i = 0; i < N; ++i) {
        R.delta[i] += bt[i];
    }
    return R;
}

JKResult jk_score_given_mq(const std::vector<double>& y,
                           const std::vector<double>& q,
                           const std::vector<double>& my,
                           const std::vector<double>& mq,
                           const std::vector<mwIndex>& idx,
                           int Mdef)
{
    JKResult R;
    R.mq = mq;
    std::vector<size_t> ii;
    for (auto i : idx) {
        if (i < q.size() && std::isfinite(y[i]) && std::isfinite(q[i]) &&
            std::isfinite(my[i]) && std::isfinite(mq[i])) {
            ii.push_back(i);
        }
    }
    size_t n = ii.size();
    if (n < 2) {
        return R;
    }
    std::vector<double> ay(n), aq(n), myr(n), mqr(n), yr(n);
    for (size_t k = 0; k < n; ++k) {
        size_t i = ii[k];
        yr[k] = y[i];
        myr[k] = my[i];
        mqr[k] = mq[i];
        ay[k] = y[i] - my[i];
        aq[k] = q[i] - mq[i];
    }
    double sy = 0, sq = 0;
    for (size_t k = 0; k < n; ++k) {
        sy += ay[k] * ay[k];
        sq += aq[k] * aq[k];
    }
    sy = std::sqrt(sy / n);
    sq = std::sqrt(sq / n);
    if (!(sy > 0 && sq > 0)) {
        return R;
    }
    double rho = 0;
    for (size_t k = 0; k < n; ++k) {
        rho += (ay[k] / sy) * (aq[k] / sq);
    }
    rho = std::max(-1.0, std::min(1.0, rho / n));
    double a = sq / sy;
    R.V = (1 - a) * (1 - a);
    R.C = (1 - rho) * (1 - rho);
    if (Mdef == 1) {
        double mmq = std::accumulate(mqr.begin(), mqr.end(), 0.0) / n,
               mmy = std::accumulate(myr.begin(), myr.end(), 0.0) / n;
        if (mmy == 0) {
            return R;
        }
        double x = mmq / mmy - 1;
        R.M = x * x;
    } else {
        double ybar = std::accumulate(yr.begin(), yr.end(), 0.0) / n, den = 0, num = 0;
        for (size_t k = 0; k < n; ++k) {
            den += (myr[k] - ybar) * (myr[k] - ybar);
            num += (mqr[k] - myr[k]) * (mqr[k] - myr[k]);
        }
        if (!(den > 0)) {
            return R;
        }
        R.M = num / den;
    }
    R.JK = 1 - std::sqrt(R.M + R.V + R.C);
    return R;
}

// ----------------------------
// Metrics and attribution
// ----------------------------
struct MB {
    double SAR = QNAN, GLS = QNAN, NSE = QNAN, KGE = QNAN, Kr = QNAN, Ka = QNAN, Kb = QNAN,
           Huber = QNAN, RSS = QNAN;
};

MB metric_block(const std::vector<double>& y,
                const std::vector<double>& q,
                double muy,
                double sy,
                double TSS,
                double Sy)
{
    MB m;
    std::vector<double> yy, qq, rr;
    for (size_t i = 0; i < y.size(); ++i) {
        double r = y[i] - q[i];
        if (std::isfinite(y[i]) && std::isfinite(q[i]) && std::isfinite(r)) {
            yy.push_back(y[i]);
            qq.push_back(q[i]);
            rr.push_back(r);
        }
    }
    if (rr.empty()) {
        return m;
    }
    m.RSS = 0;
    m.SAR = 0;
    for (double r : rr) {
        m.RSS += r * r;
        m.SAR += std::abs(r);
    }
    m.GLS = m.RSS;
    if (!(std::isfinite(TSS) && TSS > 0)) {
        double mu = std::accumulate(yy.begin(), yy.end(), 0.0) / yy.size();
        TSS = 0;
        for (double v : yy) {
            TSS += (v - mu) * (v - mu);
        }
    }
    if (TSS > 0) {
        m.NSE = 1 - m.RSS / TSS;
    }
    auto kg = kge(yy, qq, muy, sy);
    m.KGE = kg.K;
    m.Kr = kg.r;
    m.Ka = kg.a;
    m.Kb = kg.b;
    m.Huber = huber(rr, Sy, nullptr);
    return m;
}

mxArray* metrics_struct(const MB& t,
                        const MB& e,
                        double Dft,
                        double Dfe,
                        const JKResult* jt,
                        const JKResult* je)
{
    const char* fn[] = {
        "SARt",   "GLSt",  "NSEt",    "KGEt",    "KGE_rt",     "KGE_alphat", "KGE_betat",
        "Hubert", "RSSt",  "JKGEt",   "JKGE_Mt", "JKGE_Vt",    "JKGE_Ct",    "SARe",
        "GLSe",   "NSEe",  "KGEe",    "KGE_re",  "KGE_alphae", "KGE_betae",  "Hubere",
        "RSSe",   "JKGEe", "JKGE_Me", "JKGE_Ve", "JKGE_Ce",    "Dfdct",      "Dfdce"};
    mxArray* s = mxCreateStructMatrix(1, 1, sizeof(fn) / sizeof(fn[0]), fn);
    double v[] = {t.SAR,
                  t.GLS,
                  t.NSE,
                  t.KGE,
                  t.Kr,
                  t.Ka,
                  t.Kb,
                  t.Huber,
                  t.RSS,
                  jt ? jt->JK : QNAN,
                  jt ? jt->M : QNAN,
                  jt ? jt->V : QNAN,
                  jt ? jt->C : QNAN,
                  e.SAR,
                  e.GLS,
                  e.NSE,
                  e.KGE,
                  e.Kr,
                  e.Ka,
                  e.Kb,
                  e.Huber,
                  e.RSS,
                  je ? je->JK : QNAN,
                  je ? je->M : QNAN,
                  je ? je->V : QNAN,
                  je ? je->C : QNAN,
                  Dft,
                  Dfe};
    for (size_t i = 0; i < sizeof(v) / sizeof(v[0]); ++i) {
        mxSetFieldByNumber(s, 0, (int)i, mxCreateDoubleScalar(v[i]));
    }
    return s;
}

mxArray* attribution_struct(const std::vector<double>& J,
                            const std::vector<double>& delta,
                            const std::vector<double>& grad,
                            const mxArray* mdl,
                            size_t n,
                            size_t d)
{
    auto lo = vec(fld(mdl, "th_min"), "mdl.th_min"),
         hi = vec(fld(mdl, "th_max"), "mdl.th_max");
    std::vector<double> At(d, 0), An(d, 0);
    for (size_t j = 0; j < d; ++j) {
        double span = hi[j] - lo[j];
        for (size_t i = 0; i < n; ++i) {
            At[j] += std::abs(J[i + n * j] * delta[i] * span);
        }
        An[j] = std::abs(span * grad[j]);
    }
    const char* f[] = {"total", "net"};
    mxArray* s = mxCreateStructMatrix(1, 1, 2, f);
    mxSetField(s, 0, "total", col(At));
    mxSetField(s, 0, "net", col(An));
    return s;
}

// ----------------------------
// State packaging
// ----------------------------
std::vector<std::string> requested_states(const mxArray* a)
{
    std::vector<std::string> out;
    if (!a || mxGetNumberOfElements(a) == 0) {
        return out;
    }
    mxArray* rhs = const_cast<mxArray*>(a);
    mxArray* c = nullptr;
    if (mexCallMATLAB(1, &c, 1, &rhs, "cellstr") != 0 || !c) {
        mexErrMsgIdAndTxt("crr_model_mex:States", "Could not parse request.states.");
    }
    size_t n = mxGetNumberOfElements(c);
    for (size_t i = 0; i < n; ++i) {
        const mxArray* x = mxGetCell(c, (mwIndex)i);
        char* p = mxArrayToString(x);
        if (p) {
            out.emplace_back(p);
            mxFree(p);
        }
    }
    mxDestroyArray(c);
    return out;
}

std::vector<std::string> state_names(int model, int nstate)
{
    std::vector<std::string> n;
    if (model == 1) {
        n = {"snow_water_equivalent",
             "soil_storage",
             "slow_reservoir",
             "fast_reservoir_1",
             "fast_reservoir_2",
             "fast_reservoir_3"};
    } else if (model == 2) {
        n = {"snow_water_equivalent",
             "interception_storage",
             "soil_storage",
             "fast_reservoir",
             "slow_reservoir"};
    } else if (model == 3) {
        n = {"snow_water_equivalent",
             "upper_tension_water",
             "upper_free_water",
             "lower_tension_water",
             "lower_primary_free_water",
             "lower_supplemental_free_water",
             "routing_reservoir_1",
             "routing_reservoir_2",
             "routing_reservoir_3"};
    } else if (model == 4) {
        n = {"snow_water_equivalent",
             "tension_water",
             "free_water",
             "interflow_reservoir",
             "baseflow_reservoir",
             "routing_reservoir_1",
             "routing_reservoir_2",
             "routing_reservoir_3"};
    } else if (model == 5) {
        n = {"snow_water_equivalent", "production_store", "routing_store"};
        for (int i = 1; i <= nstate - 3; ++i) {
            n.push_back("unit_hydrograph_storage_" + std::to_string(i));
        }
    } else if (model == 6) {
        n = {"snow_water_equivalent", "soil_moisture", "upper_zone", "lower_zone"};
    } else if (model == 7) {
        n = {"snow_water_equivalent", "soil_storage", "groundwater_storage"};
        for (int i = 1; i <= nstate - 3; ++i) {
            n.push_back("routing_storage_" + std::to_string(i));
        }
    } else if (model == 21) {
        n = {"snow_water_equivalent", "production_store", "routing_store"};
        for (int i = 1; i <= nstate - 3; ++i) {
            n.push_back("routing_memory_" + std::to_string(i));
        }
    }
    if ((int)n.size() != nstate) {
        mexErrMsgIdAndTxt(
            "crr_model_mex:States", "State catalog mismatch for model %d.", model);
    }
    return n;
}

mxArray* cellstr_array(const std::vector<std::string>& a)
{
    mxArray* c = mxCreateCellMatrix((mwSize)a.size(), 1);
    for (size_t i = 0; i < a.size(); ++i) {
        mxSetCell(c, (mwIndex)i, mxCreateString(a[i].c_str()));
    }
    return c;
}

mxArray* string_array(const std::vector<std::string>& a)
{
    mxArray *c = cellstr_array(a), *s = nullptr;
    if (mexCallMATLAB(1, &s, 1, &c, "string") != 0 || !s) {
        mxDestroyArray(c);
        mexErrMsgIdAndTxt("crr_model_mex:States", "Could not create string array.");
    }
    mxDestroyArray(c);
    return s;
}

mxArray*
states_struct(const Sim& sim, int model, const mxArray* mdl, const mxArray* reqStates)
{
    int nstate = sim.m - 1;
    if (nstate < 1 || sim.Z.empty()) {
        mexErrMsgIdAndTxt("crr_model_mex:States", "State history unavailable.");
    }
    auto names = state_names(model, nstate);
    auto req = requested_states(reqStates);
    std::vector<int> keep;
    if (req.size() == 1 && (req[0] == "all" || req[0] == "ALL")) {
        for (int j = 0; j < nstate; ++j) {
            keep.push_back(j);
        }
    } else if (req.empty()) {
        for (int j = 0; j < nstate; ++j) {
            keep.push_back(j);
        }
    } else {
        for (const auto& r : req) {
            std::string lr = r;
            std::transform(lr.begin(), lr.end(), lr.begin(), ::tolower);
            bool found = false;
            for (int j = 0; j < nstate; ++j) {
                std::string ln = names[j];
                std::transform(ln.begin(), ln.end(), ln.begin(), ::tolower);
                if (lr == ln) {
                    keep.push_back(j);
                    found = true;
                    break;
                }
            }
            if (!found) {
                mexErrMsgIdAndTxt(
                    "crr_model_mex:States", "Unknown requested state '%s'.", r.c_str());
            }
        }
    }
    auto id = vec(fld(mdl, "idx"), "mdl.idx");
    long long i0 = (long long)std::llround(id[0]), i1 = (long long)std::llround(id[1]);
    size_t nr = (size_t)(i1 - i0), nc = keep.size();
    std::vector<double> V(nr * nc);
    for (size_t j = 0; j < nc; ++j) {
        for (size_t i = 0; i < nr; ++i) {
            size_t row = (size_t)i0 + i;
            V[i + nr * j] = sim.Z[row + (size_t)sim.zrows * (size_t)keep[j]];
        }
    }
    std::vector<std::string> sel, units;
    for (int j : keep) {
        sel.push_back(names[j]);
        units.push_back("mm");
    }
    const char* f[] = {"values", "names", "units"};
    mxArray* s = mxCreateStructMatrix(1, 1, 3, f);
    mxSetField(s, 0, "values", mat(V, (mwSize)nr, (mwSize)nc));
    mxSetField(s, 0, "names", string_array(sel));
    mxSetField(s, 0, "units", string_array(units));
    return s;
}

} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("crr_model_mex:nrhs", "Need x,mdl,dat,ode,loss,request.");
    }
    if (nlhs > 2) {
        mexErrMsgIdAndTxt("crr_model_mex:nlhs", "At most two outputs.");
    }
    const mxArray *mdl = prhs[1], *dat = prhs[2], *ode = prhs[3], *loss = prhs[4],
                  *req = prhs[5];
    auto x = vec(prhs[0], "x");
    size_t d = x.size();
    int model = (int)std::llround(scl(fld(mdl, "model"), "mdl.model"));
    int mcode = (int)std::llround(scl(fld(mdl, "mcode"), "mdl.mcode"));
    int lf = (int)std::llround(scl(fld(loss, "fnc"), "loss.fnc"));
    if (model == 8) {
        mexErrMsgIdAndTxt("crr_model_mex:UserModel",
                          "User model remains external; use crr_model_cpp.m.");
    }
    if (mcode != 4) {
        mexErrMsgIdAndTxt("crr_model_mex:Mcode", "Native branch requires mdl.mcode=4.");
    }
    if (lf < 1 || lf > 7) {
        mexErrMsgIdAndTxt("crr_model_mex:Loss", "loss.fnc must be 1..7.");
    }

    bool rqQ = lflag(req, "q"), rqG = lflag(req, "gradient"), rqJ = lflag(req, "jacobian"),
         rqM = lflag(req, "metrics"), rqA = lflag(req, "attribution");

    const mxArray* st = fld(req, "states", false);
    const bool rqS = has_requested_states(st);

    bool needG = rqG || rqA, needJ = rqJ || needG;

    auto tr = transform(x, mdl);
    mxArray* out = mxCreateStructMatrix(1, 1, 0, nullptr);
    if (!tr.valid) {
        plhs[0] = mxCreateDoubleScalar(QNAN);
        if (nlhs > 1) {
            plhs[1] = out;
        } else {
            mxDestroyArray(out);
        }
        return;
    }

    auto sim = run_builtin(model, tr.th, tr.Jth, mdl, fld(dat, "meteo"), ode, needJ, rqS);

    /* A failed native integration is not a zero-flow simulation. Core
       output buffers are preallocated before integration and may contain
       valid prefixes plus untouched zeros when a step fails. Invalidate
       the complete result so no downstream loss or metric can score an
       uncomputed tail as discharge equal to zero. */
    if (sim.fail) {
        std::fill(sim.q.begin(), sim.q.end(), QNAN);
        std::fill(sim.J.begin(), sim.J.end(), QNAN);
        std::fill(sim.Z.begin(), sim.Z.end(), QNAN);
    }

    auto bad = vec(fld(dat, "bad"), "dat.bad");
    bool local = false;
    const mxArray* la = fld(mdl, "local", false);
    if (la && !mxIsEmpty(la)) {
        local = scl(la, "mdl.local") == 1;
    }
    auto tr0 = expand_idx(local ? fld(dat, "id_train") : fld(mdl, "id_train"),
                          local ? "dat.id_train" : "mdl.id_train");
    auto ev0 = expand_idx(local ? fld(dat, "id_eval") : fld(mdl, "id_eval"),
                          local ? "dat.id_eval" : "mdl.id_eval");
    auto filterIds = [&](const std::vector<mwIndex>& a) {
        std::vector<mwIndex> z;
        for (auto i : a) {
            if ((size_t)i >= bad.size()) {
                mexErrMsgIdAndTxt("crr_model_mex:Index", "Index outside dat.bad.");
            }
            if (bad[i] == 0) {
                z.push_back(i);
            }
        }
        return z;
    };
    auto idtr = filterIds(tr0), idev = filterIds(ev0);
    auto ya = vec(fld(dat, "y_n"), "dat.y_n");
    if (ya.size() != sim.q.size()) {
        mexErrMsgIdAndTxt("crr_model_mex:Length", "q and y_n differ.");
    }
    auto sel = [&](const std::vector<double>& v, const std::vector<mwIndex>& id) {
        std::vector<double> z;
        z.reserve(id.size());
        for (auto i : id) {
            z.push_back(v[i]);
        }
        return z;
    };
    auto yt = sel(ya, idtr), qt = sel(sim.q, idtr), ye = sel(ya, idev),
         qe = sel(sim.q, idev);
    std::vector<double> rt(yt.size()), re(ye.size());
    for (size_t i = 0; i < rt.size(); ++i) {
        rt[i] = yt[i] - qt[i];
    }
    for (size_t i = 0; i < re.size(); ++i) {
        re[i] = ye[i] - qe[i];
    }

    double L = QNAN, Dft = QNAN, Dfe = QNAN;
    std::vector<double> del;
    JKResult jkt, jke;
    bool haveJK = false;
    const mxArray* stats = fld(dat, "stats");
    if (rt.empty()) {
        del.clear();
    } else if (lf == 1) {
        L = 0;
        for (double e : rt) {
            L += std::abs(e);
        }
        if (needG) {
            delta(1, yt, qt, del);
        }
    } else if (lf == 2) {
        L = 0;
        for (double e : rt) {
            L += e * e;
        }
        if (needG) {
            delta(2, yt, qt, del);
        }
    } else if (lf == 3) {
        double T = scl(fld(stats, "TSSt"), "dat.stats.TSSt");
        if (std::isfinite(T) && T > 0) {
            double R = 0;
            for (double e : rt) {
                R += e * e;
            }
            L = R / T;
        }
        if (needG) {
            delta(3, yt, qt, del);
        }
    } else if (lf == 4) {
        auto z = kge(yt,
                     qt,
                     scl(fld(stats, "mut"), "dat.stats.mut"),
                     scl(fld(stats, "stdt"), "dat.stats.stdt"));
        L = 1 - z.K;
        if (needG) {
            delta(4, yt, qt, del);
        }
    } else if (lf == 5) {
        L = huber(rt, scl(fld(stats, "Syt"), "dat.stats.Syt"), needG ? &del : nullptr);
    } else if (lf == 6) {
        const mxArray* fdc = fld(dat, "fdc");
        const mxArray* ft = fld(fdc, "t");
        Dft = fdc_loss_native(qt, ft);
        L = Dft;
        if (needG) {
            fdc_delta_native(yt, qt, del);
        }
    } else if (lf == 7) {
        int method = (int)std::llround(scl(fld(loss, "method"), "loss.method"));
        int nwin = 0;
        if (method == 1 || method == 2) {
            nwin = (int)std::llround(scl(fld(loss, "n_win"), "loss.n_win"));
        }
        int Mdef = 2;
        const mxArray* ma = fld(loss, "M", false);
        if (ma && !mxIsEmpty(ma)) {
            Mdef = (int)std::llround(scl(ma, "loss.M"));
        }
        std::vector<double> mo;
        if (method == 4) {
            mo = vec(fld(fld(loss, "meta"), "mo_all"), "loss.meta.mo_all");
        }
        const mxArray* jks = fld(dat, "jkge");
        auto my = vec(fld(jks, "m_y"), "dat.jkge.m_y");
        auto cache = jkcache(fld(jks, "cache"), method, sim.q.size());
        jkt = jk_train(ya, sim.q, my, idtr, method, nwin, mo, cache, Mdef, needG);
        L = 1 - jkt.JK;
        haveJK = true;
        if (needG) {
            del.resize(idtr.size());
            for (size_t i = 0; i < idtr.size(); ++i) {
                del[i] = jkt.delta[idtr[i]];
            }
        }
    }

    std::vector<double> Jt, grad;
    if (needJ) {
        Jt.assign(idtr.size() * d, 0);
        for (size_t j = 0; j < d; ++j) {
            for (size_t i = 0; i < idtr.size(); ++i) {
                Jt[i + idtr.size() * j] = sim.J[idtr[i] + (size_t)sim.nq * j];
            }
        }
    }
    if (needG) {
        grad.assign(d, QNAN);
        if (!rt.empty()) {
            if (del.size() != idtr.size()) {
                mexErrMsgIdAndTxt("crr_model_mex:Delta", "Delta length mismatch.");
            }
            std::fill(grad.begin(), grad.end(), 0.0);
            for (size_t j = 0; j < d; ++j) {
                for (size_t i = 0; i < idtr.size(); ++i) {
                    grad[j] += Jt[i + idtr.size() * j] * del[i];
                }
            }
        }
    }

    if (rqM) {
        const mxArray* fdc = fld(dat, "fdc");
        if (!std::isfinite(Dft)) {
            Dft = fdc_loss_native(qt, fld(fdc, "t"));
        }
        const mxArray* fe = fld(fdc, "e", false);
        if (!qe.empty() && fe && scl(fld(fe, "n"), "fdc.e.n") > 0) {
            Dfe = fdc_loss_native(qe, fe);
        }
        MB mt = metric_block(yt,
                             qt,
                             scl(fld(stats, "mut"), "stats.mut"),
                             scl(fld(stats, "stdt"), "stats.stdt"),
                             scl(fld(stats, "TSSt"), "stats.TSSt"),
                             scl(fld(stats, "Syt"), "stats.Syt"));
        MB me = metric_block(ye,
                             qe,
                             scl(fld(stats, "mue"), "stats.mue"),
                             scl(fld(stats, "stde"), "stats.stde"),
                             scl(fld(stats, "TSSe"), "stats.TSSe"),
                             scl(fld(stats, "Sye"), "stats.Sye"));
        if (lf == 7 && haveJK) {
            int Mdef = 2;
            const mxArray* ma = fld(loss, "M", false);
            if (ma && !mxIsEmpty(ma)) {
                Mdef = (int)std::llround(scl(ma, "loss.M"));
            }
            auto my = vec(fld(fld(dat, "jkge"), "m_y"), "dat.jkge.m_y");
            bool ok = true;
            for (auto i : idev) {
                if (i >= jkt.mq.size() || !std::isfinite(jkt.mq[i])) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                jke = jk_score_given_mq(ya, sim.q, my, jkt.mq, idev, Mdef);
            }
        }
        mxAddField(out, "metrics");
        mxSetField(
            out,
            0,
            "metrics",
            metrics_struct(
                mt, me, Dft, Dfe, (lf == 7 ? &jkt : nullptr), (lf == 7 ? &jke : nullptr)));
    }
    if (rqQ) {
        mxAddField(out, "q");
        mxSetField(out, 0, "q", col(sim.q));
    }
    if (rqG) {
        mxAddField(out, "gradient");
        mxSetField(out, 0, "gradient", col(grad));
    }
    if (rqJ) {
        mxAddField(out, "jacobian");
        mxSetField(out, 0, "jacobian", mat(Jt, (mwSize)idtr.size(), (mwSize)d));
    }
    if (rqA) {
        mxAddField(out, "attribution");
        mxSetField(
            out, 0, "attribution", attribution_struct(Jt, del, grad, mdl, idtr.size(), d));
    }
    if (rqS) {
        mxAddField(out, "states");
        mxSetField(out, 0, "states", states_struct(sim, model, mdl, st));
    }
    plhs[0] = mxCreateDoubleScalar(L);
    if (nlhs > 1) {
        plhs[1] = out;
    } else {
        mxDestroyArray(out);
    }
}
