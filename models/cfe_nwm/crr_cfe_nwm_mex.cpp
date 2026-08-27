/*
 * crr_cfe_nwm_mex.cpp
 * MATLAB MEX gateway for the native CFE-NWM rainfall-runoff model.
 *
 * This gateway preserves the established crr_cfe_nwm interface while
 * keeping mxArray parsing outside the numerical kernel.
 */

#include "mex.h"
#include "cfe_nwm.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
const mxArray* field(const mxArray* s, const char* n, bool req = true)
{
    if (!s || !mxIsStruct(s) || mxGetNumberOfElements(s) != 1) {
        mexErrMsgIdAndTxt("crr_cfe_nwm:Struct", "Expected scalar struct.");
    }
    const mxArray* a = mxGetField(s, 0, n);
    if (!a && req) {
        mexErrMsgIdAndTxt("crr_cfe_nwm:Field", "Missing field '%s'.", n);
    }
    return a;
}

double scalar(const mxArray* s, const char* n)
{
    const mxArray* a = field(s, n);
    if (mxGetNumberOfElements(a) != 1 || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("crr_cfe_nwm:Scalar", "Field '%s' must be scalar.", n);
    }
    return mxGetScalar(a);
}

struct InputVector {
    const double* ptr = nullptr;
    std::size_t n = 0;
    std::vector<double> buf;
};

InputVector vector_view(const mxArray* a, const char* n)
{
    InputVector v;
    v.n = mxGetNumberOfElements(a);
    if (mxIsDouble(a)) {
        v.ptr = mxGetPr(a);
        return v;
    }
    if (mxIsSingle(a)) {
        const float* p = (const float*)mxGetData(a);
        v.buf.resize(v.n);
        for (std::size_t i = 0; i < v.n; ++i) {
            v.buf[i] = (double)p[i];
        }
        v.ptr = v.buf.data();
        return v;
    }
    mexErrMsgIdAndTxt("crr_cfe_nwm:Vector", "%s must be single/double.", n);
    return v;
}
} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 4 || nlhs > 4) {
        mexErrMsgIdAndTxt("crr_cfe_nwm:Usage", "Need t_last,z0,data,options; <=4 outputs.");
    }
    const int ns = (int)std::llround(mxGetScalar(prhs[0]));
    const mxArray* data = prhs[2];
    const mxArray* options = prhs[3];
    InputVector z0 = vector_view(prhs[1], "z0"),
                P = vector_view(field(data, "P"), "data.P"),
                Ep = vector_view(field(data, "Ep"), "data.Ep"),
                T = vector_view(field(data, "T"), "data.T");
    sage_cfe_nwm::Params p;
    p.s_max = scalar(data, "s_max");
    p.s_fc = scalar(data, "s_fc");
    p.s_wp = scalar(data, "s_wp");
    p.k_sch = scalar(data, "k_sch");
    p.a1 = scalar(data, "a1");
    p.k_perc = scalar(data, "k_perc");
    p.lf_thr = scalar(data, "lf_thr");
    p.a2 = scalar(data, "a2");
    p.k_lf = scalar(data, "k_lf");
    p.g_max = scalar(data, "g_max");
    p.c_gw = scalar(data, "c_gw");
    p.mm = scalar(data, "mm");
    p.k_nsh = scalar(data, "k_nsh");
    p.T_tr = scalar(data, "T_tr");
    p.f_dd = scalar(data, "f_dd");
    p.dT = scalar(data, "dT");
    p.T_sm = scalar(data, "T_sm");
    p.eps_m = scalar(data, "eps_m");
    p.eps = scalar(data, "eps");
    p.rho = scalar(data, "rho");
    p.K = static_cast<int>(std::llround(scalar(data, "K")));
    p.L = static_cast<int>(std::llround(scalar(data, "L")));
    InputVector giuh = vector_view(field(data, "giuh"), "data.giuh");
    p.giuh = giuh.ptr;
    sage_cfe_nwm::Options opt;
    opt.InitStep = scalar(options, "InitStep");
    opt.MaxStep = scalar(options, "MaxStep");
    opt.MinStep = scalar(options, "MinStep");
    opt.RelTol = scalar(options, "RelTol");
    opt.AbsTol = scalar(options, "AbsTol");
    opt.Order = scalar(options, "Order");
    opt.maxiter = (int)std::llround(scalar(options, "maxiter"));
    int mem = 1;
    const mxArray* ma = field(options, "mem", false);
    if (ma && !mxIsEmpty(ma)) {
        mem = (int)std::llround(mxGetScalar(ma));
    }
    int ipr = 1;
    const mxArray* ia = field(data, "ipr", false);
    if (ia && !mxIsEmpty(ia)) {
        ipr = (int)std::llround(mxGetScalar(ia));
    }
    if (ipr < 1) {
        ipr = 1;
    }
    if (ipr > ns + 1) {
        ipr = ns + 1;
    }
    const int d = 15;
    const int m = 3 + p.L + p.K + 1;
    const std::size_t nvar = (std::size_t)m * (d + 1);
    if (z0.n != nvar) {
        mexErrMsgIdAndTxt("crr_cfe_nwm:z0", "Unexpected augmented state length.");
    }
    const std::size_t zr = mem ? (std::size_t)(ns + 1) : 1u,
                      nq = (!mem && ipr <= ns) ? (std::size_t)(ns - ipr + 1) : 0u;
    plhs[0] = mxCreateDoubleMatrix((mwSize)zr, (mwSize)nvar, mxREAL);
    double *q = nullptr, *J = nullptr;
    if (nlhs >= 2) {
        plhs[1] = (!mem && nq) ? mxCreateDoubleMatrix((mwSize)nq, 1, mxREAL)
                               : mxCreateDoubleMatrix(0, 0, mxREAL);
        if (nq) {
            q = mxGetPr(plhs[1]);
        }
    }
    if (nlhs >= 3) {
        plhs[2] = (!mem && nq) ? mxCreateDoubleMatrix((mwSize)nq, d, mxREAL)
                               : mxCreateDoubleMatrix(0, 0, mxREAL);
        if (nq) {
            J = mxGetPr(plhs[2]);
        }
    }
    sage_cfe_nwm::Forcing F{P.ptr, Ep.ptr, T.ptr, std::min(P.n, std::min(Ep.n, T.n))};
    sage_cfe_nwm::OutputView O{mxGetPr(plhs[0]), q, J, zr, nvar, nq, (std::size_t)d};
    bool fail =
        sage_cfe_nwm::run_into(ns, z0.ptr, z0.n, F, p, opt, mem != 0, ipr, nlhs >= 3, O);
    if (nlhs >= 4) {
        plhs[3] = mxCreateLogicalScalar(fail);
    }
}
