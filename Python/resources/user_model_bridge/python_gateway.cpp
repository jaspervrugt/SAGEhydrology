#include "mex.h"
void mexFunction(int,mxArray**,int,const mxArray**);
#ifdef _WIN32
#define API extern "C" __declspec(dllexport)
#else
#define API extern "C"
#endif
static thread_local std::string error;
API mxArray* sage_array(int m,int n,const double*v){auto a=mxCreateDoubleMatrix(m,n,0);if(v)std::copy(v,v+m*n,a->data.begin());return a;}
API mxArray* sage_object(){auto a=new mxArray;a->object=true;a->m=a->n=1;return a;}
API void sage_field(mxArray*a,const char*k,mxArray*v){mxSetField(a,0,k,v);}
API void sage_free(mxArray*a){mxDestroyArray(a);}
API int sage_rows(mxArray*a){return int(mxGetM(a));}API int sage_cols(mxArray*a){return int(mxGetN(a));}
API const double* sage_values(mxArray*a){return mxGetPr(a);}
API const char* sage_error(){return error.c_str();}
API int sage_call(const mxArray**in,mxArray**out){try{mexFunction(9,out,4,in);return 0;}catch(const std::exception&e){error=e.what();return 1;}}
