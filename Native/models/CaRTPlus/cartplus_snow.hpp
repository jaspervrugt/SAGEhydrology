// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_snow.hpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

// Daily CaRT+ snow step adapted from Xilin Xia's MIT-licensed bra_step.
// Upstream revision 83682732728646669a190835d95ab3e60a8b5cb5.
// See LICENSE in this directory.
#pragma once

#include <array>
#include <cstddef>

namespace cartplus {

inline constexpr std::size_t band_count = 5;

template <typename Scalar>
struct SnowParameters {
    Scalar threshold_c;
    Scalar degree_day_factor;
    Scalar snowfall_correction;
    Scalar sublimation_fraction;
    Scalar pet_melt_factor;
    Scalar lapse_rate_c_per_km;
};

template <typename Scalar>
struct SnowFluxes {
    Scalar liquid{};
    Scalar sublimation{};
    Scalar snowfall_adjustment{};
};

/**
 * advance_snow: Advance five snow bands and area-weighted liquid fluxes.
 *
 * INPUT ARGUMENTS:
 * snow : Snow-store depths for the five elevation bands (mm).
 * elevation_m : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 42 physical CaRT+ parameters in native order.
 * precipitation : Precipitation depth/rate on the selected time base.
 * temperature_c : Temperature c used by this numerical operator.
 * pet : Potential evaporation depth/rate on this operator time base.
 * dt : Integration or observation interval duration (days).
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template <typename Scalar>
SnowFluxes<Scalar> advance_snow(
    std::array<Scalar, band_count>& snow,
    const std::array<double, band_count>& elevation_m,
    const std::array<double, band_count>& area,
    const SnowParameters<Scalar>& parameters,
    double precipitation,
    double temperature_c,
    double pet, double dt = 1.0
) {
    double mean_elevation = 0.0;
    for (std::size_t band = 0; band < band_count; ++band) {
        mean_elevation += area[band] * elevation_m[band];
    }

    SnowFluxes<Scalar> fluxes;
    for (std::size_t band = 0; band < band_count; ++band) {
        const Scalar temperature = Scalar(temperature_c)
            - parameters.lapse_rate_c_per_km
                * Scalar((elevation_m[band] - mean_elevation) / 1000.0);
        const bool cold = temperature < parameters.threshold_c;
        if (cold) {
            snow[band] += parameters.snowfall_correction
                * Scalar(precipitation);
            fluxes.snowfall_adjustment += Scalar(area[band] * precipitation)
                * (parameters.snowfall_correction - Scalar(1.0));
        }

        const Scalar sublimation_demand =
            parameters.sublimation_fraction * Scalar(pet);
        const Scalar sublimation =
            sublimation_demand < snow[band] ? sublimation_demand : snow[band
                ];
        snow[band] -= sublimation;

        const Scalar melt_demand = cold ? Scalar(0.0)
            : parameters.degree_day_factor * dt
                * (temperature - parameters.threshold_c)
                + parameters.pet_melt_factor * Scalar(pet);
        const Scalar melt =
            melt_demand < snow[band] ? melt_demand : snow[band];
        snow[band] -= melt;

        fluxes.liquid += Scalar(area[band])
            * (melt + Scalar(cold ? 0.0 : precipitation));
        fluxes.sublimation += Scalar(area[band]) * sublimation;
    }
    return fluxes;
}

}  // namespace cartplus
