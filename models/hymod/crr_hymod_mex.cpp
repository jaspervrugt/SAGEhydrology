/*
 * crr_hymod_mex.cpp
 * MATLAB MEX gateway for the native HYMOD rainfall-runoff model.
 *
 * This gateway preserves the established crr_hymod interface while
 * keeping mxArray parsing outside the numerical kernel.
 */

#include "mex.h"
#include "hymod.hpp"
#include "../standalone_result.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
const mxArray* field(const mxArray* s, const char* n, bool req = true)
{
    if (!s || !mxIsStruct(s) || mxGetNumberOfElements(s) != 1) {
        mexErrMsgIdAndTxt("crr_hymod:Struct", "Expected scalar struct.");
    }
    const mxArray* a = mxGetField(s, 0, n);
    if (!a && req) {
        mexErrMsgIdAndTxt("crr_hymod:Field", "Missing field '%s'.", n);
    }
    return a;
}

double scalar(const mxArray* s, const char* n)
{
    const mxArray* a = field(s, n);
    if (mxGetNumberOfElements(a) != 1 || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("crr_hymod:Scalar", "Field '%s' must be scalar.", n);
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
    mexErrMsgIdAndTxt("crr_hymod:Vector", "%s must be single/double.", n);
    return v;
}
} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    const bool structured = nrhs == 5;
    if ((nrhs != 4 && nrhs != 5) || (structured ? nlhs != 1 : nlhs > 4)) {
        mexErrMsgIdAndTxt("crr_hymod:Usage", "Need t_last,z0,data,options; <=4 outputs.");
    }
    const int ns = (int)std::llround(mxGetScalar(prhs[0]));
    const mxArray* data = prhs[2];
    const mxArray* options = prhs[3];
    InputVector z0 = vector_view(prhs[1], "z0"),
                P = vector_view(field(data, "P"), "data.P"),
                Ep = vector_view(field(data, "Ep"), "data.Ep"),
                T = vector_view(field(data, "T"), "data.T");
    sage_hymod::Params p;
    p.S_umax = scalar(data, "S_umax");
    p.beta = scalar(data, "beta");
    p.alfa = scalar(data, "alfa");
    p.K_s = scalar(data, "K_s");
    p.K_f = scalar(data, "K_f");
    p.T_tr = scalar(data, "T_tr");
    p.f_dd = scalar(data, "f_dd");
    p.T_sm = scalar(data, "T_sm");
    p.eps_m = scalar(data, "eps_m");
    p.rho = scalar(data, "rho");

    sage_hymod::Options opt;
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
    std::vector<std::string> obsNames, jacNames, stateNames;
    bool wantQ = false, wantJ = false, wantSwe = false, wantJswe = false;
    if (structured) {
        obsNames = sage_standalone::names(prhs[4],"obs");
        jacNames = sage_standalone::names(prhs[4],"jac");
        stateNames = sage_standalone::names(prhs[4],"states");
        wantQ = sage_standalone::has(obsNames,"Q")
            || sage_standalone::flag(prhs[4],"q");
        wantJ = sage_standalone::has(jacNames,"Q")
            || sage_standalone::flag(prhs[4],"jacobian");
        wantSwe = sage_standalone::has(obsNames,"SWE");
        wantJswe = sage_standalone::has(jacNames,"SWE");
        mem = stateNames.empty() ? 0 : 1;
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
    const int d = 7;
    const int m = 7;
    const std::size_t nvar = (std::size_t)m * (d + 1);
    if (z0.n != nvar) {
        mexErrMsgIdAndTxt("crr_hymod:z0", "Unexpected augmented state length.");
    }
    const std::size_t zr = mem ? (std::size_t)(ns + 1) : 1u,
                      nq = (!mem && ipr <= ns) ? (std::size_t)(ns - ipr + 1) : 0u;
    mxArray* Zmx = mxCreateDoubleMatrix((mwSize)zr, (mwSize)nvar, mxREAL);
    if (!structured) plhs[0] = Zmx;
    double *q = nullptr, *J = nullptr, *swe = nullptr, *Jswe = nullptr;
    std::vector<double> qbuf, Jbuf, sweBuf, JsweBuf;
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
    if (structured && !mem) {
        if (wantQ) {
            qbuf.assign(nq,0.0);
            q = qbuf.data();
        }
        if (wantJ) {
            Jbuf.assign(nq*(std::size_t)d,0.0);
            J = Jbuf.data();
        }
        if (wantSwe) {
            sweBuf.assign(nq,0.0);
            swe = sweBuf.data();
        }
        if (wantJswe) {
            JsweBuf.assign(nq*(std::size_t)d,0.0);
            Jswe = JsweBuf.data();
        }
    }
    sage_hymod::Forcing F{P.ptr, Ep.ptr, T.ptr, std::min(P.n, std::min(Ep.n, T.n))};
    sage_hymod::OutputView O{mxGetPr(Zmx), q, J, swe, Jswe,
                             zr, nvar, nq, (std::size_t)d};
    bool fail =
        sage_hymod::run_into(ns, z0.ptr, z0.n, F, p, opt, mem != 0, ipr,
                             structured ? wantJ : nlhs >= 3, O);
    if (structured) {
        if (mem) {
            plhs[0] = sage_standalone::result(
                prhs[4],mxGetPr(Zmx),zr,ns,ipr,m,d,fail);
        } else {
            std::vector<sage_standalone::OutputChannel> channels{
                {"Q",q,J},{"SWE",swe,Jswe}};
            plhs[0] = sage_standalone::result(
                prhs[4],channels,nq,d,mxGetPr(Zmx),zr,ns,ipr,m,fail);
        }
        mxDestroyArray(Zmx);
    } else if (nlhs >= 4) {
        plhs[3] = mxCreateLogicalScalar(fail);
    }
}
