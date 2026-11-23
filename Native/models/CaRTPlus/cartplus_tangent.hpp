// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_tangent.hpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

#pragma once

#include <array>
#include <cmath>
#include <cstddef>

namespace cartplus {

inline constexpr std::size_t physical_parameter_count = 42;
#if defined(CART_SPLIT_MODEL)
inline constexpr std::size_t parameter_count = 12;
#elif defined(CART_ODE_TANGENTS) || defined(CART_BASELINE_TANGENTS)
inline constexpr std::size_t parameter_count = 11;
#elif defined(CARTPLUS_CALIBRATED_ONLY)
inline constexpr std::size_t parameter_count = 35;
#else
inline constexpr std::size_t parameter_count = physical_parameter_count;
#endif

struct Tangent {
    double value;
    std::array<double, parameter_count> slope{};

    /**
     * Tangent: Initialize the numerical scalar, routing object, or solver
     * workspace.
     *
     * INPUT ARGUMENTS:
     * v : Primal numerical value, optionally carrying tangent lanes.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    Tangent(double v = 0.0) : value(v) {}

    /**
     * chain: Combine scalar values with chain-rule derivative lanes.
     *
     * INPUT ARGUMENTS:
     * value : Primal scalar or numerical operand.
     * x : Numerical operand; scalar, vector, or forward-mode scalar.
     * factor : Multiplicative derivative or flux scaling factor.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    static Tangent chain(double value, const Tangent& x, double factor) {
        Tangent result(value);
        for (std::size_t j = 0; j < parameter_count; ++j)
            result.slope[j] = factor * x.slope[j];
        return result;
    }

    /**
     * chain: Combine scalar values with chain-rule derivative lanes.
     *
     * INPUT ARGUMENTS:
     * value : Primal scalar or numerical operand.
     * x : Numerical operand; scalar, vector, or forward-mode scalar.
     * dx : Local partial derivative with respect to x.
     * y : Numerical state/operand; augmented with sensitivities when
     * enabled.
     * dy : Local partial derivative with respect to y.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    static Tangent chain(double value, const Tangent& x, double dx,
                         const Tangent& y, double dy) {
        Tangent result(value);
        for (std::size_t j = 0; j < parameter_count; ++j)
            result.slope[j] = dx * x.slope[j] + dy * y.slope[j];
        return result;
    }

    /**
     * operator+: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator+(const Tangent& a, const Tangent& b) {
        return chain(a.value + b.value, a, 1.0, b, 1.0);
    }
    /**
     * operator-: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator-(const Tangent& a, const Tangent& b) {
        return chain(a.value - b.value, a, 1.0, b, -1.0);
    }
    /**
     * operator-: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator-(const Tangent& a) { return chain(-a.value, a, -
        1.0); }
    /**
     * operator*: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator*(const Tangent& a, const Tangent& b) {
        return chain(a.value * b.value, a, b.value, b, a.value);
    }
    /**
     * operator/: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator/(const Tangent& a, const Tangent& b) {
        const double q = a.value / b.value;
        const double inv = 1.0 / b.value;
        return chain(q, a, inv, b, -q * inv);
    }
#ifdef CART_FAST_SCALAR
    // Constants have no parameter derivatives. Avoid constructing and
    // propagating their zero tangent arrays.
    /**
     * operator+: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator+(const Tangent& a, double b) {
        Tangent result = a; result.value += b; return result;
    }
    /**
     * operator+: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator+(double a, const Tangent& b) { return b+a; }
    /**
     * operator-: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator-(const Tangent& a, double b) { return a+(-b); }
    /**
     * operator-: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator-(double a, const Tangent& b) {
        Tangent result = -b; result.value += a; return result;
    }
    /**
     * operator*: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator*(const Tangent& a, double b) {
        return chain(a.value*b,a,b);
    }
    /**
     * operator*: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator*(double a, const Tangent& b) { return b*a; }
    /**
     * operator/: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator/(const Tangent& a, double b) {
        return chain(a.value/b,a,1./b);
    }
    /**
     * operator/: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend Tangent operator/(double a, const Tangent& b) {
        const double q = a/b.value;
        return chain(q,b,-q/b.value);
    }
#endif
    /**
     * operator+=: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    Tangent& operator+=(const Tangent& b) { return *this = *this + b; }
    /**
     * operator-=: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    Tangent& operator-=(const Tangent& b) { return *this = *this - b; }
    /**
     * operator*=: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    Tangent& operator*=(const Tangent& b) { return *this = *this * b; }
    /**
     * operator<: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend bool operator<(const Tangent& a, const Tangent& b) { return a.
        value < b.value; }
    /**
     * operator>: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend bool operator>(const Tangent& a, const Tangent& b) { return b < a
        ; }
    /**
     * operator<=: Apply the scalar operation and propagate forward
     * derivative lanes.
     *
     * INPUT ARGUMENTS:
     * a : First numerical operand or elevation-band fractions.
     * b : Second numerical operand or rate coefficient.
     *
     * OUTPUT ARGUMENTS:
     * Scalar result with forward derivative lanes; comparisons return bool.
     */
    friend bool operator<=(const Tangent& a, const Tangent& b) { return !(b
        < a); }
};

using std::abs;
using std::exp;
using std::expm1;
using std::hypot;
using std::log;
using std::log1p;
using std::pow;
using std::sqrt;
using std::tanh;

/**
 * primal: Extract the real/primal value without derivative lanes.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
inline double primal(double x) { return x; }
/**
 * primal: Extract the real/primal value without derivative lanes.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
inline double primal(const Tangent& x) { return x.value; }
/**
 * abs: Take an absolute value and propagate its selected derivative.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Scalar result with forward derivative lanes; comparisons return bool.
 */
inline Tangent abs(const Tangent& x) { return x.value < 0 ? -x : x; }
/**
 * exp: Evaluate exp and propagate forward derivatives.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Scalar result with forward derivative lanes; comparisons return bool.
 */
inline Tangent exp(const Tangent& x) {
    const double y = std::exp(x.value);
    return Tangent::chain(y, x, y);
}
/**
 * expm1: Evaluate expm1 and propagate forward derivatives.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Scalar result with forward derivative lanes; comparisons return bool.
 */
inline Tangent expm1(const Tangent& x) {
    return Tangent::chain(std::expm1(x.value), x, std::exp(x.value));
}
/**
 * log: Evaluate log and propagate forward derivatives.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Scalar result with forward derivative lanes; comparisons return bool.
 */
inline Tangent log(const Tangent& x) {
    return Tangent::chain(std::log(x.value), x, 1.0 / x.value);
}
/**
 * log1p: Evaluate log1p and propagate forward derivatives.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Scalar result with forward derivative lanes; comparisons return bool.
 */
inline Tangent log1p(const Tangent& x) {
    return Tangent::chain(std::log1p(x.value), x, 1.0 / (1.0 + x.value));
}
/**
 * sqrt: Evaluate sqrt and propagate forward derivatives.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Scalar result with forward derivative lanes; comparisons return bool.
 */
inline Tangent sqrt(const Tangent& x) {
    const double y = std::sqrt(x.value);
    return y == 0 ? Tangent(y) : Tangent::chain(y, x, 0.5 / y);
}
/**
 * tanh: Evaluate tanh and propagate forward derivatives.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Scalar result with forward derivative lanes; comparisons return bool.
 */
inline Tangent tanh(const Tangent& x) {
    const double y = std::tanh(x.value);
    return Tangent::chain(y, x, 1.0 - y * y);
}
/**
 * pow: Evaluate a power and propagate forward derivative lanes.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 * exponent : Dimensionless power-law exponent.
 *
 * OUTPUT ARGUMENTS:
 * Scalar result with forward derivative lanes; comparisons return bool.
 */
inline Tangent pow(const Tangent& x, double exponent) {
#ifdef CART_FAST_POWERS
    if (exponent == 4.) {
        const double square = x.value*x.value;
        return Tangent::chain(square*square,x,4.*square*x.value);
    }
    if (exponent == 3.5 && x.value >= 0.) {
        const double root = std::sqrt(x.value);
        const double square = x.value*x.value;
        return Tangent::chain(square*x.value*root,x,3.5*square*root);
    }
#endif
    return Tangent::chain(std::pow(x.value, exponent), x,
                          exponent * std::pow(x.value, exponent - 1.0));
}
/**
 * hypot: Evaluate hypot and propagate forward derivatives.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 *
 * OUTPUT ARGUMENTS:
 * Scalar result with forward derivative lanes; comparisons return bool.
 */
inline Tangent hypot(const Tangent& x, const Tangent& y) {
    const double z = std::hypot(x.value, y.value);
    return Tangent::chain(z, x, x.value / z, y, y.value / z);
}
template<class A, class B>
auto minimum(A a,B b) -> decltype(a+b) { return b<a ? b : a; }
template<class A, class B>
auto maximum(A a,B b) -> decltype(a+b) { return a<b ? b : a; }
template<class A, class B, class C>
auto bounded(A a,B lo,C hi) -> decltype(a+lo+hi) {
    return minimum(maximum(a,lo),hi);
}

}  // namespace cartplus
