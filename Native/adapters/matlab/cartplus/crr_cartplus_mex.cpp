/*
 * CARTPLUS_MEX  Native CaRT+ gateway for MATLAB.
 *
 * Inputs: N-by-4 forcing [P, T, PET, training mask], five elevations,
 * five area fractions, 42 physical parameters, steps/day, infiltration
 * floor, and optionally seven fixed spin-up stores.
 * Outputs: N-by-4 flux/storage balance, optional N-by-42 discharge
 * Jacobian at fixed spin-up, and optional seven baseline spin-up stores.
 */
#include "mex.h"
#include "cartplus_bridge.cpp"

namespace {

thread_local const mxArray* zone_solver = nullptr;

/* Call the MATLAB integrator for storage and its optional derivatives. */
int solve_zone(double initial, double net, double recession, double dt,
               int derivative, double* output)
{
    auto* input = mxCreateDoubleMatrix(1, 5, mxREAL);
    auto* values = mxGetPr(input);
    values[0] = initial;
    values[1] = net;
    values[2] = recession;
    values[3] = dt;
    values[4] = derivative;
    mxArray* args[] = {const_cast<mxArray*>(zone_solver), input};
    mxArray* result = nullptr;
    auto* error = mexCallMATLABWithTrap(1, &result, 2, args, "feval");
    mxDestroyArray(input);
    const std::size_t count = derivative ? 8 : 2;
    const bool ok = !error && result && mxIsDouble(result)
        && mxGetNumberOfElements(result) == count;
    if (ok) std::copy_n(mxGetPr(result), count, output);
    if (result) mxDestroyArray(result);
    if (error) mxDestroyArray(error);
    return ok ? 0 : 1;
}

/* Require full, real, finite double arrays before exposing raw data. */
const double* input(const mxArray* value, std::size_t size)
{
    if (!mxIsDouble(value) || mxIsComplex(value) || mxIsSparse(value)
        || mxGetNumberOfElements(value) != size) {
        mexErrMsgIdAndTxt("cartplus:input", "Invalid input type or size.");
    }
    const double* values = mxGetPr(value);
    if (!std::all_of(values, values + size,
                     [](double x) { return std::isfinite(x); })) {
        mexErrMsgIdAndTxt("cartplus:input", "Inputs must be finite.");
    }
    return values;
}

/* Copy a native row-major matrix into a MATLAB column-major result. */
mxArray* matrix(const std::vector<double>& values, std::size_t rows,
                std::size_t columns)
{
    mxArray* result = mxCreateDoubleMatrix(rows, columns, mxREAL);
    auto* target = mxGetPr(result);
    for (std::size_t j = 0; j < columns; ++j) {
        for (std::size_t i = 0; i < rows; ++i) {
            target[j * rows + i] = values[i * columns + j];
        }
    }
    return result;
}

} // namespace

/* Validate, marshal, and run the same kernel used by SAGE Python. */
extern "C" void mexFunction(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[])
{
    if (nrhs < 6 || nrhs > 8 || nlhs < 1 || nlhs > 3) {
        mexErrMsgIdAndTxt("cartplus:arity",
                         "Expected six to eight inputs and 1–3 outputs.");
    }
    const auto n = mxGetM(prhs[0]);
    if (n == 0 || mxGetN(prhs[0]) != 4) {
        mexErrMsgIdAndTxt("cartplus:forcing", "Expected N-by-4 forcing.");
    }
    const auto* weather = input(prhs[0], n * 4);
    const auto* elevation = input(prhs[1], 5);
    const auto* area = input(prhs[2], 5);
    const auto* theta = input(prhs[3], 42);
    const double resolution = *input(prhs[4], 1);
    const double floor = *input(prhs[5], 1);
    if ((resolution != 1 && resolution != 24 && resolution != 96)
        || floor < 0) {
        mexErrMsgIdAndTxt("cartplus:interval", "Invalid interval or floor.");
    }
    cartplus::SpinupStates fixed{}, captured{};
    const bool has_fixed = nrhs >= 7 && !mxIsEmpty(prhs[6]);
    if (has_fixed) {
        const auto* v = input(prhs[6], 7);
        fixed = {v[0], v[1], v[2], v[3], v[4], v[5], v[6]};
    }
    std::vector<double> forcing(n * 4), output(n * 4);
    for (std::size_t j = 0; j < 4; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            forcing[i * 4 + j] = weather[j * n + i];
        }
    }
    const auto* initial = has_fixed ? &fixed : nullptr;
    const auto dt = 1.0 / resolution;
    zone_solver = nrhs == 8 && !mxIsEmpty(prhs[7]) ? prhs[7] : nullptr;
    auto callback = zone_solver ? solve_zone : nullptr;
    if (cartplus_solver_forward(n, forcing.data(), elevation, area,
                                theta, output.data(), initial,
                                &captured, floor, dt, callback)) {
        zone_solver = nullptr;
        mexErrMsgIdAndTxt("cartplus:forward", "CaRT+ simulation failed.");
    }
    plhs[0] = matrix(output, n, 4);
    if (nlhs >= 2) {
        std::vector<double> jacobian(n * 42);
        if (cartplus_solver_jacobian(n, forcing.data(), elevation, area,
                                       theta, jacobian.data(),
                                       has_fixed ? &fixed : &captured,
                                       floor, dt, callback)) {
            zone_solver = nullptr;
            mexErrMsgIdAndTxt("cartplus:jacobian", "CaRT+ Jacobian failed.");
        }
        plhs[1] = matrix(jacobian, n, 42);
    }
    if (nlhs == 3) {
        plhs[2] = matrix({captured.surface, captured.soil, captured.zone,
                          captured.wetness, captured.depressions,
                          captured.groundwater, captured.quick}, 1, 7);
    }
    zone_solver = nullptr;
}
