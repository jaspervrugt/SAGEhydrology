// MATLAB gateway to the same native ODE kernels used by SITE/SAGE Python.
#include "mex.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#ifndef _WIN32
#define __declspec(attribute)
#endif

#ifdef MATLAB_CARTPLUS
#include "../../Native/models/CaRTPlus/cartplus_continuous.cpp"
#else
#include "../../Native/models/cart/cart_ode.cpp"
#endif

namespace {
const double* data(const mxArray* value, std::size_t count) {
    if (!mxIsDouble(value) || mxIsComplex(value) || mxIsSparse(value)
        || mxGetNumberOfElements(value) != count)
        mexErrMsgIdAndTxt("cart:ode:input", "Invalid CaRT ODE input dimensions.");
    const double* values = mxGetPr(value);
    if (!std::all_of(values, values + count, [](double x) { return std::isfinite(x); }))
        mexErrMsgIdAndTxt("cart:ode:input", "CaRT ODE inputs must be finite.");
    return values;
}
}

extern "C" __declspec(dllexport) void mexFunction(int nlhs, mxArray* out[], int nrhs, const mxArray* in[]) {
    if (nrhs != 8 || nlhs < 1 || nlhs > 2 || mxGetN(in[0]) != 4 || mxGetM(in[0]) == 0)
        mexErrMsgIdAndTxt("cart:ode:input", "Expected forcing, bands, parameters and ODE settings.");
#ifdef MATLAB_CARTPLUS
    constexpr int width = 42;
#else
    constexpr int width = 11;
#endif
    const auto n = mxGetM(in[0]);
    const auto* forcing = data(in[0], 4*n);
    const auto* elevation = data(in[1], 5);
    const auto* area = data(in[2], 5);
    const auto* parameters = data(in[3], width);
    const auto* durations = data(in[4], n);
    const double m9 = *data(in[5], 1);
    const auto* options = data(in[6], 7);
    const double method = *data(in[7], 1);
    if (m9 < 1 || m9 != std::floor(m9) || m9 > (std::numeric_limits<int>::max()-19)/3
        || method < 0 || method > 11 || (method > 2 && method < 9) || method != std::floor(method)
        || options[0] <= 0 || options[1] <= 0
        || !std::all_of(durations, durations+n, [](double dt) { return dt > 0; }))
        mexErrMsgIdAndTxt("cart:ode:input", "Invalid CaRT ODE integration settings.");
    std::vector<double> weather(4*n);
    for (std::size_t i=0; i<n; ++i)
        for (std::size_t j=0; j<4; ++j)
            weather[4*i+j] = forcing[j*n+i];

    out[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    std::vector<double> jacobian(nlhs == 2 ? n*width : 0);
    double* jac = nlhs == 2 ? jacobian.data() : nullptr;
    const int solver = static_cast<int>(method);
    const int stages = static_cast<int>(m9);
#ifdef MATLAB_CARTPLUS
    std::array<double, width> seed;
    seed.fill(1.0);
    const int status = cartplus_configured_run(n, weather.data(), elevation, area, parameters,
        durations, nullptr, stages, options, seed.data(), mxGetPr(out[0]), jac, nullptr, solver, options[5], options[6]);
#else
    const int status = cart_configured_run(n, weather.data(), elevation, area, parameters,
        durations, nullptr, stages, options[0], options[1], solver, mxGetPr(out[0]), jac, nullptr, options[5], options[6]);
#endif
    if (status)
        mexErrMsgIdAndTxt("cart:ode:solve", "Native CaRT ODE integration failed.");
    if (nlhs == 2) {
        out[1] = mxCreateDoubleMatrix(n, width, mxREAL);
        auto* J = mxGetPr(out[1]);
        for (std::size_t j=0; j<width; ++j)
            for (std::size_t i=0; i<n; ++i)
                J[j*n+i] = jacobian[i*width+j];
    }
}
