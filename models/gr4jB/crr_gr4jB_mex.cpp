/* Standalone MATLAB gateway for the modular GR4J-B native core. */

#include "mex.h"
#include "gr4jB.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace {

const mxArray* get_field(
    const mxArray* value,
    const char* name)
{
    if (!value ||
        !mxIsStruct(value) ||
        mxGetNumberOfElements(value) != 1) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:Struct",
            "Expected a scalar structure.");
    }
    const mxArray* result = mxGetField(value, 0, name);
    if (!result) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:Field",
            "Missing field '%s'.",
            name);
    }
    return result;
}

double get_scalar(
    const mxArray* value,
    const char* name)
{
    const mxArray* result = get_field(value, name);
    if (mxIsComplex(result) ||
        mxGetNumberOfElements(result) != 1 ||
        !mxIsNumeric(result)) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:Scalar",
            "Field '%s' must be a real scalar.",
            name);
    }
    return mxGetScalar(result);
}

struct VectorView {
    std::vector<double> owned;
    const double* data = nullptr;
    std::size_t size = 0;
};

VectorView get_vector(
    const mxArray* value,
    const char* name)
{
    if (!value ||
        mxIsComplex(value) ||
        !(mxIsDouble(value) || mxIsSingle(value))) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:Vector",
            "%s must be a real single or double vector.",
            name);
    }

    VectorView result;
    result.size = mxGetNumberOfElements(value);
    if (mxIsDouble(value)) {
        result.data = mxGetPr(value);
    } else {
        const float* source =
            static_cast<const float*>(mxGetData(value));
        result.owned.assign(source, source + result.size);
        result.data = result.owned.data();
    }
    return result;
}

} // namespace

void mexFunction(
    int nlhs,
    mxArray* plhs[],
    int nrhs,
    const mxArray* prhs[])
{
    if (nrhs != 4 || nlhs > 4) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:Usage",
            "Use [Z,q,J,fail] = crr_gr4jB("
            "t_last,z0,data,options).");
    }

    const int ns = static_cast<int>(
        std::llround(mxGetScalar(prhs[0])));
    if (ns < 1) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:Time",
            "t_last must be a positive integer.");
    }

    VectorView z0 = get_vector(prhs[1], "z0");
    VectorView P = get_vector(
        get_field(prhs[2], "P"), "P");
    VectorView Ep = get_vector(
        get_field(prhs[2], "Ep"), "Ep");
    VectorView T = get_vector(
        get_field(prhs[2], "T"), "T");
    VectorView U = get_vector(
        get_field(prhs[2], "U"), "U");
    VectorView dU = get_vector(
        get_field(prhs[2], "dUdx4"), "dUdx4");

    const int L = static_cast<int>(
        std::llround(get_scalar(prhs[2], "L")));
    constexpr int d = 8;
    const int m = 4 + 2 * L;
    const int nvar = m * (d + 1);
    if (L < 1 ||
        U.size != static_cast<std::size_t>(L) ||
        dU.size != static_cast<std::size_t>(L) ||
        z0.size != static_cast<std::size_t>(nvar)) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:Dimensions",
            "Inconsistent L, U, dUdx4, or z0 length.");
    }

    const bool mem = get_scalar(prhs[3], "mem") != 0;
    int ipr = static_cast<int>(
        std::llround(get_scalar(prhs[2], "ipr")));
    ipr = std::max(1, std::min(ipr, ns + 1));
    const std::size_t zrows =
        mem ? static_cast<std::size_t>(ns + 1) : 1u;
    const std::size_t nq =
        (!mem && ipr <= ns)
            ? static_cast<std::size_t>(ns - ipr + 1)
            : 0u;

    plhs[0] = mxCreateDoubleMatrix(
        static_cast<mwSize>(zrows),
        static_cast<mwSize>(nvar),
        mxREAL);
    double* q = nullptr;
    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(
            static_cast<mwSize>(nq), 1, mxREAL);
        q = mxGetPr(plhs[1]);
    }
    double* J = nullptr;
    if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleMatrix(
            static_cast<mwSize>(nq), d, mxREAL);
        J = mxGetPr(plhs[2]);
    }

    sage_gr4jb::Params p{};
    p.x1 = get_scalar(prhs[2], "x1");
    p.x2 = get_scalar(prhs[2], "x2");
    p.x3 = get_scalar(prhs[2], "x3");
    p.x4 = get_scalar(prhs[2], "x4");
    p.x5 = get_scalar(prhs[2], "x5");
    p.f_p = get_scalar(prhs[2], "f_p");
    p.T_tr = get_scalar(prhs[2], "T_tr");
    p.f_dd = get_scalar(prhs[2], "f_dd");
    p.eta = get_scalar(prhs[2], "eta");
    p.tau = get_scalar(prhs[2], "tau");
    p.kappa = get_scalar(prhs[2], "kappa");
    p.bg = get_scalar(prhs[2], "bg");
    p.bR = get_scalar(prhs[2], "bR");
    p.mts = get_scalar(prhs[2], "mts");
    p.T_sm = get_scalar(prhs[2], "T_sm");
    p.eps_m = get_scalar(prhs[2], "eps_m");
    p.eps_s = get_scalar(prhs[2], "eps_s");
    p.rho = get_scalar(prhs[2], "rho");
    p.L = L;
    p.U = U.data;
    p.dUdx4 = dU.data;

    sage_gr4jb::Options options{
        get_scalar(prhs[3], "InitStep"),
        get_scalar(prhs[3], "MaxStep"),
        get_scalar(prhs[3], "MinStep"),
        get_scalar(prhs[3], "RelTol"),
        get_scalar(prhs[3], "AbsTol"),
        get_scalar(prhs[3], "Order"),
        static_cast<int>(std::llround(
            get_scalar(prhs[3], "maxiter")))};
    sage_gr4jb::Forcing forcing{
        P.data,
        Ep.data,
        T.data,
        std::min(P.size, std::min(Ep.size, T.size))};
    sage_gr4jb::OutputView output{
        mxGetPr(plhs[0]),
        q,
        J,
        zrows,
        static_cast<std::size_t>(nvar),
        nq,
        static_cast<std::size_t>(d)};

    const bool fail = sage_gr4jb::run_into(
        ns,
        z0.data,
        z0.size,
        forcing,
        p,
        options,
        mem,
        ipr,
        nlhs >= 3,
        output);
    if (nlhs >= 4) {
        plhs[3] = mxCreateLogicalScalar(fail);
    }
}
