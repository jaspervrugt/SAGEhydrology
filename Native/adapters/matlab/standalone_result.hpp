#pragma once
#include "mex.h"
#include <algorithm>
#include <cctype>
#include <cstddef>
#include <string>
#include <vector>

namespace sage_standalone {
inline std::vector<std::string> names(const mxArray* req, const char* field)
{
    std::vector<std::string> out;
    if (!req || !mxIsStruct(req)) return out;
    const mxArray* a = mxGetField(req,0,field);
    if (!a) return out;
    mxArray* rhs=const_cast<mxArray*>(a); mxArray* empty=nullptr;
    if (mexCallMATLAB(1,&empty,1,&rhs,"isempty") || !empty) return out;
    bool isEmpty=mxIsLogicalScalarTrue(empty); mxDestroyArray(empty);
    if (isEmpty) return out;
    mxArray* c=nullptr;
    if (mexCallMATLAB(1,&c,1,&rhs,"cellstr") || !c)
        mexErrMsgIdAndTxt("sage:standalone:Request","Could not parse request.%s.",field);
    for (mwIndex i=0;i<mxGetNumberOfElements(c);++i) {
        char* p=mxArrayToString(mxGetCell(c,i));
        if (p) { std::string s(p); mxFree(p); std::transform(s.begin(),s.end(),s.begin(),::toupper); out.push_back(s); }
    }
    mxDestroyArray(c); return out;
}
inline bool has(const std::vector<std::string>& a,const char* x)
{ return std::find(a.begin(),a.end(),std::string(x))!=a.end(); }
inline bool flag(const mxArray* req,const char* field)
{
    if (!req || !mxIsStruct(req)) return false;
    const mxArray* a=mxGetField(req,0,field);
    return a && mxIsLogicalScalarTrue(a);
}
inline mxArray* matrix(const std::vector<double>& v,mwSize nr,mwSize nc)
{ mxArray* a=mxCreateDoubleMatrix(nr,nc,mxREAL); std::copy(v.begin(),v.end(),mxGetPr(a)); return a; }
inline void set(mxArray* s,const char* name,mxArray* value)
{ mxAddField(s,name); mxSetField(s,0,name,value); }

struct OutputChannel {
    const char* name = nullptr;
    const double* value = nullptr;
    const double* jacobian = nullptr;
};

inline const OutputChannel* channel(
    const std::vector<OutputChannel>& channels,const std::string& name)
{
    for (const auto& item : channels) {
        if (item.name && name == item.name) return &item;
    }
    return nullptr;
}

inline mxArray* result(const mxArray* req,
                       const std::vector<OutputChannel>& channels,
                       std::size_t n,int d,const double* Z,std::size_t zr,
                       int ns,int ipr,int m,bool failed)
{
    if (!req || !mxIsStruct(req) || mxGetNumberOfElements(req)!=1)
        mexErrMsgIdAndTxt("sage:standalone:Request",
            "Fifth input must be a scalar request structure.");
    auto obs=names(req,"obs"), jac=names(req,"jac");
    if (flag(req,"q") && !has(obs,"Q")) obs.push_back("Q");
    if (flag(req,"jacobian") && !has(jac,"Q")) jac.push_back("Q");
    const bool wantStates=!names(req,"states").empty();
    mxArray* out=mxCreateStructMatrix(1,1,0,nullptr);
    set(out,"failed",mxCreateLogicalScalar(failed));
    if (!obs.empty()) {
        mxArray* s=mxCreateStructMatrix(1,1,0,nullptr);
        for (const auto& name : obs) {
            const auto* item=channel(channels,name);
            if (!item || !item->value)
                mexErrMsgIdAndTxt("sage:standalone:Observable",
                    "Requested observable '%s' is unavailable.",name.c_str());
            std::vector<double> value(item->value,item->value+n);
            set(s,name.c_str(),matrix(value,(mwSize)n,1));
        }
        set(out,"obs",s);
    }
    if (!jac.empty()) {
        mxArray* s=mxCreateStructMatrix(1,1,0,nullptr);
        for (const auto& name : jac) {
            const auto* item=channel(channels,name);
            if (!item || !item->jacobian)
                mexErrMsgIdAndTxt("sage:standalone:Jacobian",
                    "Requested Jacobian '%s' is unavailable.",name.c_str());
            std::vector<double> J(item->jacobian,item->jacobian+n*(std::size_t)d);
            set(s,name.c_str(),matrix(J,(mwSize)n,(mwSize)d));
        }
        set(out,"jac",s);
    }
    if (wantStates) {
        if (!Z || zr!=(std::size_t)(ns+1))
            mexErrMsgIdAndTxt("sage:standalone:States",
                "Full state history was not retained.");
        std::vector<double> V(n*(std::size_t)(m-1));
        for (int j=0;j<m-1;++j)
            for (std::size_t i=0;i<n;++i)
                V[i+n*j]=Z[(std::size_t)ipr+i+zr*j];
        mxArray* s=mxCreateStructMatrix(1,1,0,nullptr);
        set(s,"values",matrix(V,(mwSize)n,(mwSize)(m-1)));
        set(out,"states",s);
    }
    return out;
}

inline mxArray* result(const mxArray* req,const double* Z,std::size_t zr,
                       int ns,int ipr,int m,int d,bool failed,
                       const std::vector<int>& smIndices = {},
                       const std::vector<int>& sweIndices = {0})
{
    if (!req || !mxIsStruct(req) || mxGetNumberOfElements(req)!=1)
        mexErrMsgIdAndTxt("sage:standalone:Request","Fifth input must be a scalar request structure.");
    const auto obs=names(req,"obs"), jac=names(req,"jac");
    const bool oq=has(obs,"Q") || flag(req,"q"), os=has(obs,"SWE"),
               om=has(obs,"SM");
    const bool jq=has(jac,"Q") || flag(req,"jacobian"), js=has(jac,"SWE"),
               jm=has(jac,"SM");
    if ((om || jm) && smIndices.empty())
        mexErrMsgIdAndTxt("sage:standalone:Observable",
            "SM is not available for this model.");
    const bool wantStates=!names(req,"states").empty();
    const std::size_t n=(ipr<=ns)?(std::size_t)(ns-ipr+1):0u;
    mxArray* out=mxCreateStructMatrix(1,1,0,nullptr);
    set(out,"failed",mxCreateLogicalScalar(failed));
    if (oq || os || om) {
        mxArray* s=mxCreateStructMatrix(1,1,0,nullptr);
        if (oq) {
            std::vector<double> q(n);
            for (std::size_t i=0;i<n;++i) { std::size_t rb=(std::size_t)ipr+i,ra=rb-1; q[i]=Z[rb+zr*(m-1)]-Z[ra+zr*(m-1)]; }
            set(s,"Q",matrix(q,(mwSize)n,1));
        }
        if (os) {
            std::vector<double> swe(n,0.0);
            for (int index : sweIndices)
                for (std::size_t i=0;i<n;++i)
                    swe[i]+=Z[(std::size_t)ipr+i+
                        zr*(std::size_t)index];
            set(s,"SWE",matrix(swe,(mwSize)n,1));
        }
        if (om) {
            std::vector<double> sm(n,0.0);
            for (int index : smIndices)
                for (std::size_t i=0;i<n;++i)
                    sm[i]+=Z[(std::size_t)ipr+i+zr*(std::size_t)index];
            set(s,"SM",matrix(sm,(mwSize)n,1));
        }
        set(out,"obs",s);
    }
    if (jq || js || jm) {
        mxArray* s=mxCreateStructMatrix(1,1,0,nullptr);
        if (jq) {
            std::vector<double> J(n*(std::size_t)d);
            for (int j=0;j<d;++j) { std::size_t col=(std::size_t)(j+2)*m-1; for (std::size_t i=0;i<n;++i) { std::size_t rb=(std::size_t)ipr+i,ra=rb-1; J[i+n*j]=Z[rb+zr*col]-Z[ra+zr*col]; } }
            set(s,"Q",matrix(J,(mwSize)n,(mwSize)d));
        }
        if (js) {
            std::vector<double> J(n*(std::size_t)d,0.0);
            for (int j=0;j<d;++j)
                for (int index : sweIndices)
                    for (std::size_t i=0;i<n;++i)
                        J[i+n*(std::size_t)j]+=
                            Z[(std::size_t)ipr+i+
                              zr*((std::size_t)(j+1)*m+
                                  (std::size_t)index)];
            set(s,"SWE",matrix(J,(mwSize)n,(mwSize)d));
        }
        if (jm) {
            std::vector<double> J(n*(std::size_t)d,0.0);
            for (int j=0;j<d;++j)
                for (int index : smIndices)
                    for (std::size_t i=0;i<n;++i)
                        J[i+n*(std::size_t)j]+=
                            Z[(std::size_t)ipr+i+
                              zr*((std::size_t)(j+1)*m+
                                  (std::size_t)index)];
            set(s,"SM",matrix(J,(mwSize)n,(mwSize)d));
        }
        set(out,"jac",s);
    }
    if (wantStates) {
        std::vector<double> V(n*(std::size_t)(m-1));
        for (int j=0;j<m-1;++j) for (std::size_t i=0;i<n;++i) V[i+n*j]=Z[(std::size_t)ipr+i+zr*j];
        mxArray* s=mxCreateStructMatrix(1,1,0,nullptr);
        set(s,"values",matrix(V,(mwSize)n,(mwSize)(m-1)));
        set(out,"states",s);
    }
    return out;
}
} // namespace sage_standalone
