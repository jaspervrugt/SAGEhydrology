// MATLAB adapter for the shared HBV routing implementation.
// Copyright (c) 2026 Jasper A. Vrugt, University of California, Irvine.
#include "mex.h"
#include "hbv_routing.hpp"
#include <limits>
void mexFunction(int nlhs, mxArray* out[], int nrhs, const mxArray* in[])
{
    if (nrhs != 4 || nlhs < 1 || nlhs > 2)
        mexErrMsgIdAndTxt("hbv_route:Arguments", "Use [q,J]=hbv_route_mex(q,J,B,L).");
    for (int k=0;k<4;++k)
        if (!mxIsDouble(in[k]) || mxIsComplex(in[k]) || mxIsSparse(in[k]))
            mexErrMsgIdAndTxt("hbv_route:Type", "Real full double arrays required.");
    const std::size_t n=mxGetNumberOfElements(in[0]);
    const bool jac=mxGetNumberOfElements(in[1])!=0;
    if (n==0 || (mxGetM(in[0])!=1 && mxGetN(in[0])!=1)
        || (jac && (mxGetM(in[1])!=n || mxGetN(in[1])!=12))
        || mxGetNumberOfElements(in[2])!=1
        || mxGetNumberOfElements(in[3])!=1)
        mexErrMsgIdAndTxt("hbv_route:Shape", "Need q vector, n-by-12 J (or []), B, L.");
    const double length=mxGetScalar(in[3]);
    if (!std::isfinite(length) || length<1 || std::floor(length)!=length
        || length>std::numeric_limits<int>::max())
        mexErrMsgIdAndTxt("hbv_route:Length", "L must be a positive integer.");
    out[0]=mxCreateDoubleMatrix(n,1,mxREAL);
    if (nlhs>1) out[1]=mxCreateDoubleMatrix(jac?n:0,jac?13:0,mxREAL);
    try {
        sage_hbv_routing::route(mxGetDoubles(in[0]),
            jac?mxGetDoubles(in[1]):nullptr,n,0,mxGetScalar(in[2]),
            static_cast<int>(length),mxGetDoubles(out[0]),
            jac && nlhs>1?mxGetDoubles(out[1]):nullptr);
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("hbv_route:Failure", "%s", e.what());
    }
}
