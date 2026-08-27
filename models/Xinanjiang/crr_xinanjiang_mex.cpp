/*
 * crr_xinanjiang_mex.cpp
 * MATLAB MEX gateway for the native Xinanjiang rainfall-runoff model.
 *
 * This gateway preserves the established crr_xinanjiang interface while
 * keeping mxArray parsing outside the numerical kernel.
 */

#include "mex.h"
#include "xinanjiang.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
const mxArray* field(const mxArray* s, const char* n, bool req = true)
{
    if (!s || !mxIsStruct(s) || mxGetNumberOfElements(s) != 1) {
        mexErrMsgIdAndTxt("crr_xinanjiang:Struct", "Expected scalar struct.");
    }
    const mxArray* a = mxGetField(s, 0, n);
    if (!a && req) {
        mexErrMsgIdAndTxt("crr_xinanjiang:Field", "Missing field '%s'.", n);
    }
    return a;
}

double scalar(const mxArray* s, const char* n)
{
    const mxArray* a = field(s, n);
    if (mxGetNumberOfElements(a) != 1 || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("crr_xinanjiang:Scalar", "Field '%s' must be scalar.", n);
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
    mexErrMsgIdAndTxt("crr_xinanjiang:Vector", "%s must be single/double.", n);
    return v;
}
} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 4 || nlhs > 4) {
        mexErrMsgIdAndTxt("crr_xinanjiang:Usage",
                          "Need t_last,z0,data,options; <=4 outputs.");
    }
    const int ns = (int)std::llround(mxGetScalar(prhs[0]));
    const mxArray* data = prhs[2];
    const mxArray* options = prhs[3];
    InputVector z0 = vector_view(prhs[1], "z0"),
                P = vector_view(field(data, "P"), "data.P"),
                Ep = vector_view(field(data, "Ep"), "data.Ep"),
                T = vector_view(field(data, "T"), "data.T");
    sage_xinanjiang::Params p;
    p.f_p = scalar(data, "f_p");
    p.A_im = scalar(data, "A_im");
    p.a = scalar(data, "a");
    p.b = scalar(data, "b");
    p.W_max = scalar(data, "W_max");
    p.LM = scalar(data, "LM");
    p.c = scalar(data, "c");
    p.S_max = scalar(data, "S_max");
    p.S_tot = scalar(data, "S_tot");
    p.f_wm = scalar(data, "f_wm");
    p.f_lm = scalar(data, "f_lm");
    p.Ex = scalar(data, "Ex");
    p.k_i = scalar(data, "k_i");
    p.k_g = scalar(data, "k_g");
    p.c_i = scalar(data, "c_i");
    p.c_g = scalar(data, "c_g");
    p.k_f = scalar(data, "k_f");
    p.T_tr = scalar(data, "T_tr");
    p.f_dd = scalar(data, "f_dd");
    p.T_sm = scalar(data, "T_sm");
    p.eps_m = scalar(data, "eps_m");
    p.eps_r = scalar(data, "eps_r");
    p.eps = scalar(data, "eps");
    p.rho = scalar(data, "rho");

    sage_xinanjiang::Options opt;
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
    const int d = 16;
    const int m = 9;
    const std::size_t nvar = (std::size_t)m * (d + 1);
    if (z0.n != nvar) {
        mexErrMsgIdAndTxt("crr_xinanjiang:z0", "Unexpected augmented state length.");
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
    sage_xinanjiang::Forcing F{P.ptr, Ep.ptr, T.ptr, std::min(P.n, std::min(Ep.n, T.n))};
    sage_xinanjiang::OutputView O{mxGetPr(plhs[0]), q, J, zr, nvar, nq, (std::size_t)d};
    bool fail =
        sage_xinanjiang::run_into(ns, z0.ptr, z0.n, F, p, opt, mem != 0, ipr, nlhs >= 3, O);
    if (nlhs >= 4) {
        plhs[3] = mxCreateLogicalScalar(fail);
    }
}
