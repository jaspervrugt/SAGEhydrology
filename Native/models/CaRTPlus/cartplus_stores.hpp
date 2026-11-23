// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_stores.hpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

// CaRT+ store updates adapted from Xia's MIT-licensed reference.
// Upstream revision 83682732728646669a190835d95ab3e60a8b5cb5.
// See LICENSE in this directory.
#pragma once

#include <algorithm>
#include <cmath>

#include "cartplus_tangent.hpp"

namespace cartplus {

template <typename Scalar>
struct SurfaceFluxes {
    Scalar overflow;
    Scalar drainage;
    Scalar evaporation;
};

/**
 * advance_surface: Advance surface interception, overflow, and drainage.
 *
 * INPUT ARGUMENTS:
 * storage : Water store (mm); reference arguments are updated in place.
 * rain : Liquid input depth/rate on this operator time base.
 * demand : Evaporation demand on this operator time base.
 * capacity : Store or infiltration capacity in this operator.
 * drainage_fraction : Fraction of retained surface storage draining this
 * step.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template <typename Scalar>
inline SurfaceFluxes<Scalar> advance_surface(
    Scalar& storage, Scalar rain, Scalar demand,
    Scalar capacity, Scalar drainage_fraction
) {
    storage += rain;
    const Scalar evaporation = minimum(demand, storage);
    storage -= evaporation;
    const Scalar overflow = maximum(storage - capacity, 0.0);
    storage -= overflow;
    const Scalar drainage = drainage_fraction * storage;
    storage -= drainage;
    return {overflow, drainage, evaporation};
}

template <typename Scalar>
struct ProductionFluxes {
    Scalar runoff;
    Scalar recharge;
    Scalar evaporation;
};

/**
 * advance_pareto: Advance soil production, evaporation, and recharge.
 *
 * INPUT ARGUMENTS:
 * storage : Water store (mm); reference arguments are updated in place.
 * rain : Liquid input depth/rate on this operator time base.
 * demand : Evaporation demand on this operator time base.
 * capacity : Store or infiltration capacity in this operator.
 * percolation : Percolation coefficient in the model parameterization.
 * shape : Dimensionless storage/infiltration shape parameter.
 * evap_shape : Dimensionless evaporation response exponent.
 * net_demand : Whether rain and evaporation are netted before production.
 * dt : Integration or observation interval duration (days).
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template <typename Scalar>
inline ProductionFluxes<Scalar> advance_pareto(
    Scalar& storage, Scalar rain, Scalar demand, Scalar capacity,
    Scalar percolation, Scalar shape, Scalar evap_shape, bool net_demand,
        double dt = 1.0
) {
    const Scalar input = net_demand ? maximum(rain - demand, 0.0) : rain;
    const Scalar potential = net_demand ? maximum(demand - rain, 0.0)
                                        : demand;
    const Scalar interception = net_demand ? minimum(rain, demand) : 0.0;
    const Scalar saturation = bounded(storage / capacity, 0.0, 1.0);
    const Scalar deficit = maximum(1.0 - saturation, 1e-7);
    const Scalar exponent = shape + 1.0;
    const Scalar depth = exp(log(deficit) / exponent);
    const Scalar remaining = depth - input / (exponent * capacity);
    const Scalar volume = remaining > 0.0
        ? exp(exponent * log(maximum(remaining, 1e-30))) : 0.0;
    const Scalar retained =
        bounded(capacity * (deficit - volume), 0.0, input);
    const Scalar evap_response = tanh(potential / capacity);
    const Scalar dry_fraction = exp(evap_shape * log(deficit));
    const Scalar evaporation = minimum(
        capacity * (1.0 - dry_fraction) * evap_response
            / (1.0 + dry_fraction / deficit * evap_response),
        storage
    );
    storage = storage - evaporation + retained;
    const Scalar relative = storage / capacity;
    const Scalar squared = relative * relative;
    const Scalar recharge = -storage
        * expm1(-0.25 * dt * log1p(percolation * squared * squared));
    storage -= recharge;
    return {input - retained, recharge, evaporation + interception};
}

template <typename Scalar>
struct GroundwaterFluxes {
    Scalar discharge;
    Scalar uptake;
};

/**
 * advance_groundwater: Advance groundwater release and evaporation uptake.
 *
 * INPUT ARGUMENTS:
 * storage : Water store (mm); reference arguments are updated in place.
 * recharge : Recharge supplied to the groundwater store.
 * residual_demand : Evaporation demand remaining after upstream stores.
 * recession : Groundwater recession coefficient in the model time base.
 * uptake_factor : Dimensionless groundwater evaporation-uptake factor.
 * reference_storage : Groundwater reference storage (mm).
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template <typename Scalar>
inline GroundwaterFluxes<Scalar> advance_groundwater(
    Scalar& storage, Scalar recharge, Scalar residual_demand,
    Scalar recession, Scalar uptake_factor, Scalar reference_storage
) {
    storage += recharge;
    const Scalar x = uptake_factor * residual_demand
        / (storage + reference_storage);
    const Scalar uptake_fraction = x < 1e-3
        ? x - 0.5 * x * x + x * x * x / 6.0 : -expm1(-x);
    const Scalar uptake = storage * uptake_fraction;
    storage -= uptake;
    const Scalar next = storage / (1.0 + recession * storage);
    const Scalar discharge = storage - next;
    storage = next;
    return {discharge, uptake};
}

}  // namespace cartplus
