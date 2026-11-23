// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_zone.hpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

// CaRT+ near-stream zone and wetness state adapted from Xia's reference.
// Upstream revision 83682732728646669a190835d95ab3e60a8b5cb5.
// See LICENSE in this directory.
#pragma once

#include <algorithm>
#include <cmath>

#include "cartplus_tangent.hpp"

namespace cartplus {

template <typename Scalar>
struct ZoneFluxes {
    Scalar runoff;
    Scalar discharge;
    Scalar evaporation;
};

/**
 * daily_decay_fraction: Convert a daily recession rate into its release
 * fraction.
 *
 * INPUT ARGUMENTS:
 * rate : Recession or cascade release rate (1/day).
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template <typename Scalar>
inline Scalar daily_decay_fraction(Scalar rate) {
    return rate < 1e-3
        ? rate * (1.0 - rate / 2.0 + rate * rate / 6.0)
        : -expm1(-rate);
}

/**
 * advance_zone: Advance near-stream storage, runoff, and drainage.
 *
 * INPUT ARGUMENTS:
 * storage : Water store (mm); reference arguments are updated in place.
 * rain : Liquid input depth/rate on this operator time base.
 * demand : Evaporation demand on this operator time base.
 * capacity : Store or infiltration capacity in this operator.
 * residence_time : Store residence time (days).
 * dt : Integration or observation interval duration (days).
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template <typename Scalar>
inline ZoneFluxes<Scalar> advance_zone(
    Scalar& storage, Scalar rain, Scalar demand,
    Scalar capacity, Scalar residence_time, double dt = 1.0
) {
    const Scalar initial = storage / capacity;
    const Scalar recession = dt / residence_time;
    const Scalar input = maximum(rain - demand, 0.0);
    Scalar next, integrated, runoff;
    if (demand > rain) {
        const Scalar drying = (demand - rain) / capacity;
        const Scalar rate = 2.0 * drying + recession;
        const Scalar decay = daily_decay_fraction(rate);
        const Scalar x = initial * (drying / rate) * decay;
        next = initial * (1.0 - decay) / (1.0 - x);
        const Scalar integral_factor = abs(x) < 1e-3
            ? 1.0 + x / 2.0 + x * x / 3.0 : -log1p(-x) / x;
        integrated = initial * (decay / rate) * integral_factor;
        runoff = 0.0;
    } else {
        const Scalar wetting = input / capacity;
        const Scalar rate = hypot(recession, 2.0 * wetting);
        const Scalar equilibrium = 2.0 * wetting / (recession + rate);
        const Scalar denominator = 2.0 * wetting * initial + recession +
            rate;
        const Scalar ratio = 2.0 * wetting * (initial - equilibrium)
            / denominator;
        const Scalar decay = daily_decay_fraction(rate);
        const Scalar surviving = 1.0 - decay;
        next = equilibrium + (initial - equilibrium) * (1.0 - ratio)
            * surviving / (1.0 - ratio * surviving);
        const Scalar y = ratio * decay / (1.0 - ratio);
        const Scalar integral_factor = abs(y) < 1e-3
            ? 1.0 - y / 2.0 + y * y / 3.0 : log1p(y) / y;
        integrated = equilibrium + 2.0 * (initial - equilibrium) * decay
            / (denominator * (1.0 - ratio)) * integral_factor;
        runoff = input - capacity * (next - initial)
            - capacity * recession * integrated;
    }
    const Scalar discharge = capacity * recession * integrated;
    const Scalar evaporation = rain - runoff - discharge
        - (capacity * next - storage);
    storage = capacity * next;
    return {runoff, discharge, evaporation};
}

}  // namespace cartplus
