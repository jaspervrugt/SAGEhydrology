// Shared differentiable HBV MAXBAS routing.
// Copyright (c) 2026 Jasper A. Vrugt, University of California, Irvine.
// Used unchanged by the MATLAB standalone/combined and Python interfaces.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#if defined(__GNUC__) && defined(__x86_64__)
#include <immintrin.h>
#endif

// Optimize only this shared FIR implementation, without fast-math or
// contraction. The surrounding ODE code retains its existing build flags.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("O3", "fp-contract=off")
#endif

namespace sage_hbv_routing {

inline void weights(double b, int length, std::vector<double>& w,
                    std::vector<double>& dw)
{
    if (!std::isfinite(b) || length < 1)
        throw std::invalid_argument("Invalid HBV routing width/window.");
    constexpr double ep = 1e-8;
    b = std::max(b, ep);
    const double half = .5 * b;
    const double root = std::sqrt((half-ep)*(half-ep) + ep*ep);
    const double den = .5 * (half + ep + root);
    const double dden = .25 * (1 + (half-ep)/root);
    std::vector<double> a(length), da(length);
    double sa = 0, sda = 0;
    for (int k = 0; k < length; ++k) {
        const double z = k-half;
        const double az = std::sqrt(z*z + ep*ep);
        const double tri = 1-az/den;
        const double r = std::sqrt(tri*tri + ep*ep);
        a[k] = .5 * (tri+r);
        da[k] = .5 * (1+tri/r)
            * (.5*z/az/den + az/(den*den)*dden);
        sa += a[k];
        sda += da[k];
    }
    w.resize(length);
    dw.resize(length);
    for (int k = 0; k < length; ++k) {
        w[k] = a[k]/sa;
        dw[k] = (da[k]*sa-a[k]*sda)/(sa*sa);
    }
}

#if defined(__GNUC__) && defined(__x86_64__)
// Runtime dispatch keeps the library portable to CPUs without AVX2. Each
// lane performs the same separate multiply and add, in the same tap order.
__attribute__((target("avx2")))
inline void filter_avx2(const std::vector<double>& w, const double* x,
                        std::size_t n, double* y)
{
    constexpr std::size_t block_size = 256;
    for (std::size_t block = 0; block < n; block += block_size) {
        const std::size_t end = std::min(n, block+block_size);
        const std::size_t last_tap = std::min(w.size(), end);
        for (std::size_t k = 0; k < last_tap; ++k) {
            const double coefficient = w[k];
            if (coefficient == 0.0) continue;
            const __m256d c = _mm256_set1_pd(coefficient);
            std::size_t i = std::max(k, block);
            for (; i+4 <= end; i+=4) {
                const __m256d a = _mm256_loadu_pd(x+i-k);
                const __m256d b = _mm256_loadu_pd(y+i);
                _mm256_storeu_pd(y+i, _mm256_add_pd(b,_mm256_mul_pd(c,a)));
            }
            for (; i < end; ++i) y[i] += coefficient*x[i-k];
        }
    }
}
#endif

inline void filter(const std::vector<double>& w, const double* x,
                   std::size_t n, double* y)
{
    std::fill_n(y, n, 0.0);
#if defined(__GNUC__) && defined(__x86_64__)
    static const bool avx2 = [] {
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx2") != 0;
    }();
    if (avx2) {
        filter_avx2(w,x,n,y);
        return;
    }
#endif
    // Portable fallback. Only exact zeros are skipped; never trim tails.
    constexpr std::size_t block_size = 256;
    for (std::size_t block = 0; block < n; block += block_size) {
        const std::size_t end = std::min(n, block+block_size);
        const std::size_t last_tap = std::min(w.size(), end);
        for (std::size_t k = 0; k < last_tap; ++k) {
            const double coefficient = w[k];
            if (coefficient == 0.0) continue;
            for (std::size_t i = std::max(k, block); i < end; ++i)
                y[i] += coefficient*x[i-k];
        }
    }
}

inline void route(const double* q, const double* jac, std::size_t n,
                  std::size_t start, double width, int length,
                  double* routed_q, double* routed_jac)
{
    if (!q || !routed_q || n == 0 || start >= n
        || (routed_jac && !jac))
        throw std::invalid_argument("Invalid HBV routing arrays/start.");
    std::vector<double> w, dw;
    weights(width, length, w, dw);
    std::copy_n(q, start, routed_q);
    filter(w, q+start, n-start, routed_q+start);
    if (routed_jac) {
        for (int j = 0; j < 12; ++j) {
            const double* input = jac+n*j;
            double* output = routed_jac+n*j;
            std::copy_n(input, start, output);
            filter(w, input+start, n-start, output+start);
        }
        std::fill_n(routed_jac+12*n, start, 0.0);
        filter(dw, q+start, n-start, routed_jac+12*n+start);
    }
}

} // namespace sage_hbv_routing

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif
