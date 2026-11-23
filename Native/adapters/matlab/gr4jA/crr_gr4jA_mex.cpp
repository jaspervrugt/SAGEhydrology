/*
 * crr_gr4jA_mex.cpp
 * MATLAB MEX gateway for the native GR4J-A rainfall-runoff model.
 *
 * This gateway preserves the established crr_gr4jA interface while
 * keeping mxArray parsing outside the numerical kernel.
 */

#include "mex.h"
#include "gr4jA.hpp"
#include "../standalone_result.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
const mxArray* field(const mxArray* s, const char* n, bool req = true)
{
    if (!s || !mxIsStruct(s) || mxGetNumberOfElements(s) != 1) {
        mexErrMsgIdAndTxt("crr_gr4jA:Struct", "Expected scalar struct.");
    }
    const mxArray* a = mxGetField(s, 0, n);
    if (!a && req) {
        mexErrMsgIdAndTxt("crr_gr4jA:Field", "Missing field '%s'.", n);
    }
    return a;
}

double scalar(const mxArray* s, const char* n)
{
    const mxArray* a = field(s, n);
    if (mxGetNumberOfElements(a) != 1 || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("crr_gr4jA:Scalar", "Field '%s' must be scalar.", n);
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
    mexErrMsgIdAndTxt("crr_gr4jA:Vector", "%s must be single/double.", n);
    return v;
}
} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    const bool structured = nrhs == 5;
    if ((nrhs != 4 && nrhs != 5) || (structured ? nlhs != 1 : nlhs > 4)) {
        mexErrMsgIdAndTxt("crr_gr4jA:Usage", "Need t_last,z0,data,options; <=4 outputs.");
    }
    const int ns = (int)std::llround(mxGetScalar(prhs[0]));
    const mxArray* data = prhs[2];
    const mxArray* options = prhs[3];
    InputVector z0 = vector_view(prhs[1], "z0"),
                P = vector_view(field(data, "P"), "data.P"),
                Ep = vector_view(field(data, "Ep"), "data.Ep"),
                T = vector_view(field(data, "T"), "data.T");
    sage_gr4ja::Params p;
    p.x1 = scalar(data, "x1");
    p.x2 = scalar(data, "x2");
    p.x3 = scalar(data, "x3");
    p.x4 = scalar(data, "x4");
    p.x5 = scalar(data, "x5");
    p.f_p = scalar(data, "f_p");
    p.T_tr = scalar(data, "T_tr");
    p.f_dd = scalar(data, "f_dd");
    p.kappa = scalar(data, "kappa");
    p.bg = scalar(data, "bg");
    p.bR = scalar(data, "bR");
    p.mts = scalar(data, "mts");
    p.T_sm = scalar(data, "T_sm");
    p.eps_m = scalar(data, "eps_m");
    p.eps_s = scalar(data, "eps_s");
    p.rho = scalar(data, "rho");
    p.n1 = static_cast<int>(std::llround(scalar(data, "n1")));
    p.n2 = static_cast<int>(std::llround(scalar(data, "n2")));

    sage_gr4ja::Options opt;
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
    std::vector<std::string> on,jn,sn;bool wq=0,wj=0,ws=0,wjs=0,wm=0,wjm=0;if(structured){on=sage_standalone::names(prhs[4],"obs");jn=sage_standalone::names(prhs[4],"jac");sn=sage_standalone::names(prhs[4],"states");wq=sage_standalone::has(on,"Q")||sage_standalone::flag(prhs[4],"q");wj=sage_standalone::has(jn,"Q")||sage_standalone::flag(prhs[4],"jacobian");ws=sage_standalone::has(on,"SWE");wjs=sage_standalone::has(jn,"SWE");wm=sage_standalone::has(on,"SM");wjm=sage_standalone::has(jn,"SM");mem=sn.empty()?0:1;}
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
    const int d = 8;
    const int m = p.n1 + p.n2 + 4;
    const std::size_t nvar = (std::size_t)m * (d + 1);
    if (z0.n != nvar) {
        mexErrMsgIdAndTxt("crr_gr4jA:z0", "Unexpected augmented state length.");
    }
    const std::size_t zr = mem ? (std::size_t)(ns + 1) : 1u,
                      nq = (!mem && ipr <= ns) ? (std::size_t)(ns - ipr + 1) : 0u;
    mxArray* Zmx = mxCreateDoubleMatrix((mwSize)zr, (mwSize)nvar, mxREAL);
    if (!structured) plhs[0] = Zmx;
    double *q=nullptr,*J=nullptr,*swe=nullptr,*Jswe=nullptr,*sm=nullptr,*Jsm=nullptr;std::vector<double> qb,Jb,sb,Jsb,mb,Jmb;
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
    if(structured&&!mem){if(wq){qb.assign(nq,0);q=qb.data();}if(wj){Jb.assign(nq*d,0);J=Jb.data();}if(ws){sb.assign(nq,0);swe=sb.data();}if(wjs){Jsb.assign(nq*d,0);Jswe=Jsb.data();}if(wm){mb.assign(nq,0);sm=mb.data();}if(wjm){Jmb.assign(nq*d,0);Jsm=Jmb.data();}}
    sage_gr4ja::Forcing F{P.ptr, Ep.ptr, T.ptr, std::min(P.n, std::min(Ep.n, T.n))};
    sage_gr4ja::OutputView O{mxGetPr(Zmx), q, J, zr, nvar, nq, (std::size_t)d};
    O.swe=swe;O.Jswe=Jswe;O.sm=sm;O.Jsm=Jsm;
    bool fail =
        sage_gr4ja::run_into(ns,z0.ptr,z0.n,F,p,opt,mem!=0,ipr,structured?(wj||wjs||wjm):nlhs>=3,O);
    if (structured) {
        if(mem)plhs[0]=sage_standalone::result(prhs[4],mxGetPr(Zmx),zr,ns,ipr,m,d,fail,{1});else{std::vector<sage_standalone::OutputChannel> c{{"Q",q,J},{"SWE",swe,Jswe},{"SM",sm,Jsm}};plhs[0]=sage_standalone::result(prhs[4],c,nq,d,mxGetPr(Zmx),zr,ns,ipr,m,fail);}
        mxDestroyArray(Zmx);
    } else if (nlhs >= 4) {
        plhs[3] = mxCreateLogicalScalar(fail);
    }
}
