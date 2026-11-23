#pragma once
#include <vector>
#include <string>
#include <map>
#include <stdexcept>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <algorithm>
using mwSize=std::size_t;using mwIndex=std::size_t;
const int mxREAL=0;
struct mxArray {size_t m=0,n=0;bool object=false,logical=false;std::vector<double> data;std::vector<std::string> names;std::map<std::string,mxArray*> fields;};
inline mxArray* mxCreateDoubleMatrix(size_t m,size_t n,int){auto a=new mxArray;a->m=m;a->n=n;a->data.resize(m*n);return a;}
inline mxArray* mxCreateDoubleScalar(double v){auto a=mxCreateDoubleMatrix(1,1,0);a->data[0]=v;return a;}
inline mxArray* mxCreateLogicalScalar(bool v){auto a=mxCreateDoubleScalar(v);a->logical=true;return a;}
inline double* mxGetPr(const mxArray* a){return a?const_cast<double*>(a->data.data()):nullptr;}
inline void* mxGetData(const mxArray* a){return mxGetPr(a);}
inline double mxGetScalar(const mxArray* a){return a&&!a->data.empty()?a->data[0]:0;}
inline size_t mxGetM(const mxArray* a){return a?a->m:0;}inline size_t mxGetN(const mxArray* a){return a?a->n:0;}
inline size_t mxGetNumberOfElements(const mxArray*a){return a?(a->object?1:a->data.size()):0;}
inline bool mxIsEmpty(const mxArray*a){return !mxGetNumberOfElements(a);}inline bool mxIsComplex(const mxArray*){return false;}
inline bool mxIsStruct(const mxArray*a){return a&&a->object;}inline bool mxIsDouble(const mxArray*a){return a&&!a->object&&!a->logical;}
inline bool mxIsSingle(const mxArray*){return false;}inline bool mxIsLogical(const mxArray*a){return a&&a->logical;}
inline bool mxIsLogicalScalarTrue(const mxArray*a){return mxGetScalar(a)!=0;}
inline bool mxIsFinite(double v){return std::isfinite(v);}inline double mxGetNaN(){return std::numeric_limits<double>::quiet_NaN();}inline double mxGetEps(){return std::numeric_limits<double>::epsilon();}
inline mxArray* mxGetField(const mxArray*a,size_t,const char*k){if(!a)return nullptr;auto it=a->fields.find(k);return it==a->fields.end()?nullptr:it->second;}
inline int mxGetFieldNumber(const mxArray*a,const char*k){if(!a)return -1;auto it=std::find(a->names.begin(),a->names.end(),k);return it==a->names.end()?-1:int(it-a->names.begin());}
inline int mxAddField(mxArray*a,const char*k){auto i=mxGetFieldNumber(a,k);if(i>=0)return i;a->names.push_back(k);a->object=true;return int(a->names.size()-1);}
inline void mxDestroyArray(mxArray*a){if(a){for(auto&v:a->fields)mxDestroyArray(v.second);delete a;}}
inline void mxSetField(mxArray*a,size_t,const char*k,mxArray*v){mxAddField(a,k);auto old=mxGetField(a,0,k);if(old&&old!=v)mxDestroyArray(old);a->fields[k]=v;}
inline void mxSetFieldByNumber(mxArray*a,size_t i,int k,mxArray*v){mxSetField(a,i,a->names.at(k).c_str(),v);}
inline mxArray* mxDuplicateArray(const mxArray*a){if(!a)return nullptr;auto b=new mxArray(*a);b->fields.clear();for(auto&v:a->fields)b->fields[v.first]=mxDuplicateArray(v.second);return b;}
inline void* mxMalloc(size_t n){return std::malloc(n);}inline void mxFree(void*p){std::free(p);}
inline void mexErrMsgIdAndTxt(const char*id,const char*fmt,...){char out[2048];va_list args;va_start(args,fmt);vsnprintf(out,sizeof(out),fmt,args);va_end(args);throw std::runtime_error(std::string(id)+": "+out);}
inline void mexErrMsgTxt(const char*s){throw std::runtime_error(s);}inline int mexPrintf(const char*,...){return 0;}
