// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/cart/crr_cart_mex.cpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

#include "mex.h"
#include "cart_daily.cpp"

/**
 * mexFunction: Marshal MATLAB inputs and return model discharge and
 * sensitivities.
 *
 * INPUT ARGUMENTS:
 * nlhs : Number of MATLAB output arguments.
 * lhs : MATLAB output-array handles, allocated by this gateway.
 * nrhs : Number of MATLAB input arguments.
 * rhs : Callback evaluating state derivatives on the model time base.
 *
 * OUTPUT ARGUMENTS:
 * lhs : MATLAB output arrays allocated by this gateway.
 */
void mexFunction(int nlhs, mxArray* lhs[], int nrhs,
                 const mxArray* rhs[])
{
    // Inputs have already been checked by the public MATLAB wrapper.
    if (nrhs != 7 || nlhs != 3)
        mexErrMsgIdAndTxt("cart:arguments","Call the MATLAB cart wrapper.");
    const int n = static_cast<int>(mxGetM(rhs[0]));
    const auto* input = mxGetPr(rhs[0]);
    std::vector<double> forcing(4*n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < 4; ++j)
            forcing[4*i+j] = input[j*n+i];

    const bool fixed = mxGetNumberOfElements(rhs[4]) == 3;
    const bool derivatives = mxIsLogicalScalarTrue(rhs[5]);
    lhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    lhs[1] = mxCreateDoubleMatrix(derivatives ? n : 0,
                                 derivatives ? cart::count : 0,mxREAL);
    lhs[2] = mxCreateDoubleMatrix(3,1,mxREAL);
    if (fixed) std::copy_n(mxGetPr(rhs[4]),3,mxGetPr(lhs[2]));
    std::vector<double> jacobian(derivatives ? n*cart::count : 0);
    #ifdef CART_UNIFORM_GATEWAY
    const int status = cart_uniform(
#else
    const int status = cart_daily(
#endif
        forcing.data(),n,mxGetPr(rhs[1]),mxGetPr(rhs[2]),
        mxGetPr(rhs[3]),mxGetPr(lhs[0]),
        derivatives ? jacobian.data() : nullptr,mxGetPr(lhs[2]),fixed,
        static_cast<int>(mxGetScalar(rhs[6])));
    if (status)
        mexErrMsgIdAndTxt("cart:training","A training mask is required.");
    for (int j = 0; derivatives && j < cart::count; ++j)
        for (int i = 0; i < n; ++i)
            mxGetPr(lhs[1])[j*n+i] = jacobian[i*cart::count+j];
}
