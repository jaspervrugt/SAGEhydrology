// C ABI for the shared HBV routing implementation; no model equations here.
// Copyright (c) 2026 Jasper A. Vrugt, University of California, Irvine.
#include "hbv_routing.hpp"
#ifdef _WIN32
#define API __declspec(dllexport)
#else
#define API __attribute__((visibility("default")))
#endif
extern "C" API int sage_hbv_route(const double* q, const double* jac,
    int n, int start, double width, int length,
    double* routed_q, double* routed_jac) noexcept
{
    if (n < 1 || start < 0 || start >= n) return 2;
    try {
        sage_hbv_routing::route(q, jac, n, start, width, length,
                               routed_q, routed_jac);
        return 0;
    } catch (...) { return 3; }
}
extern "C" API int sage_hbv_weights(double width, int length,
    double* w, double* dw) noexcept
{
    if (!w || !dw) return 2;
    try {
        std::vector<double> weights, derivatives;
        sage_hbv_routing::weights(width, length, weights, derivatives);
        std::copy(weights.begin(), weights.end(), w);
        std::copy(derivatives.begin(), derivatives.end(), dw);
        return 0;
    } catch (...) { return 3; }
}
