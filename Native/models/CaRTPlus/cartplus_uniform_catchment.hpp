// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_uniform_catchment.hpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

// Daily CAMELS-US CaRT+ forward pass based on Xia's MIT reference.
// Upstream revision f986af265633a7a980836c2e0f2124077d632420 (hourly).
// See LICENSE in this directory.
#pragma once

#include "cartplus_uniform.hpp"
#include "cartplus_snow.hpp"
#include "cartplus_stores.hpp"
#include "cartplus_zone.hpp"
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include "../diagnostic_hooks.hpp"

namespace cartplus {

/**
 * simulate_uniform: Advance uniform-input CaRT+ water balance and routing.
 *
 * INPUT ARGUMENTS:
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * p : 42 physical CaRT+ parameters in native order.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * fixed_spinup : Optional detached baseline stores for derivative replay.
 * captured_spinup : Optional buffer receiving detached spin-up stores.
 * infiltration_floor : Minimum infiltration capacity (mm/day).
 * dt : Integration or observation interval duration (days).
 * integrate_zone : Near-stream forward/sensitivity integration callback.
 * discharge_jacobian : Caller-owned buffer for discharge sensitivities.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template <typename Scalar, bool DischargeOnly = false>
inline std::vector<CatchmentDayT<Scalar>> simulate_uniform(
    const std::vector<ForcingDay>& forcing, const CatchmentParametersT<
        Scalar>& p,
    const std::array<double, band_count>& elevation,
    const std::array<double, band_count>& area,
    const SpinupStates* fixed_spinup = nullptr,
    SpinupStates* captured_spinup = nullptr,
    double infiltration_floor = 0.0,
    double dt = 1.0,
    ZoneIntegrator integrate_zone = nullptr,
    double* discharge_jacobian = nullptr
) {
    const auto n = forcing.size();
    if (!n || !(dt == 1.0 || dt == 1.0 / 24 || dt == 1.0 / 96))
        throw std::invalid_argument(
            "CaRT+ requires daily, hourly or 15-minute forcing");
    const auto steps_per_day = static_cast<std::size_t>(std::round(1.0 / dt)
        );
    const auto veg_window = 2048 * steps_per_day;
    const auto infiltration_window = 36 * steps_per_day;
    const auto zone_step = [&](Scalar& storage, Scalar rain, Scalar demand,
                               Scalar capacity, Scalar residence) {
        if (integrate_zone) {
            const Scalar initial = storage / capacity;
            const Scalar net = (rain - demand) / capacity;
            const Scalar recession = dt / residence;
            double result[8]{};
            const int derivative = std::is_same_v<Scalar, Tangent>;
            if (integrate_zone(primal(initial), primal(net), primal(
                recession),
                               dt, derivative, result)
                || !std::all_of(result, result + 8,
                                [](double x) { return std::isfinite(x); }))
                                {
                throw std::runtime_error("CaRT+ zone integration failed");
            }
            const Scalar next = zone_tangent(
                result[0], initial, net, recession, result + 2);
            const Scalar integral = zone_tangent(
                result[1], initial, net, recession, result + 3);
            const Scalar discharge = capacity * recession * integral;
            const Scalar change = capacity * next - storage;
            const Scalar runoff = net > 0
                ? rain - demand - discharge - change : Scalar(0);
            storage = capacity * next;
            return ZoneFluxes<Scalar>{runoff, discharge,
                                      rain - runoff - discharge - change};
        }
        return advance_zone(storage, rain, demand, capacity, residence, dt);
    };
    const auto spin = fixed_spinup
        ? 0 : std::min<std::size_t>(1460 * steps_per_day, n);
    const auto quick_spin = fixed_spinup
        ? 0 : std::min<std::size_t>(730 * steps_per_day, n);
    const auto baseline = std::min<std::size_t>(365 * steps_per_day, n);
    const auto diagnostics = DischargeOnly ? 0 : n;
    std::vector<Scalar> veg(n), rain(n), demand(n), sublimation(diagnostics)
        ,
        snow_change(diagnostics), snow_adjustment(diagnostics),
        surface_evap(diagnostics), surface_change(diagnostics),
        excess(n), fast(n), recharge(n), evaporation(n), zone_flow(n),
        soil_fraction(n), wetness(n),
        storage_change(diagnostics), groundwater_flow(n), routed(n), bypass(
            n),
        routing_change(diagnostics), connectivity(n), pre_spin_rain(spin),
        pre_spin_demand(spin);
    Scalar mean_veg = 0.0, mean_rain = 0.0;
    std::size_t training_days = 0;
    for (std::size_t day = 0; day < n; ++day) {
        veg[day] = logistic((forcing[day].temperature_c
            - p.vegetation_temperature) / p.vegetation_width);
        if (day < baseline) mean_veg += veg[day] / baseline;
        if (forcing[day].training) {
            mean_rain += forcing[day].precipitation;
            ++training_days;
        }
    }
    if (!training_days) throw std::invalid_argument("No CaRT+ training days"
        );
    mean_rain = maximum(mean_rain / training_days, 0.1 * dt);
    const Scalar veg_decay = exp(-dt / p.vegetation_days);
    const Scalar infiltration_decay = exp(-dt / 3.0);
    const Scalar veg_tail = (1.0 - veg_decay) * pow(veg_decay, veg_window);
    const Scalar infiltration_tail = (1.0 - infiltration_decay)
        * pow(infiltration_decay, infiltration_window);
    const Scalar surface_drainage = dt / p.surface_days;
    const Scalar infiltration_mean = (1.0-infiltration_decay)/(dt/3.0);
    Scalar veg_memory = 0.0, previous_snow = 0.0;
    Scalar surface = 0.5 * p.surface_capacity;
    std::array<Scalar, band_count> snow{};
    for (std::size_t day = 0; day < n; ++day) {
        const auto& f = forcing[day];
        veg_memory = veg_decay * veg_memory + (1.0 - veg_decay)
            * (veg[day] - mean_veg);
        if (day >= veg_window) veg_memory -= veg_tail * (veg[day-veg_window]
            - mean_veg);
        demand[day] = f.pet * (1.0 - p.vegetation_reduction
            * (1.0 - bounded(mean_veg + veg_memory, 0.0, 1.0)));
        const auto flux = uniform_snow(snow, elevation, area, p.snow,
            f.precipitation, f.temperature_c, f.pet, dt);
        rain[day] = flux.liquid;
        if constexpr (!DischargeOnly) {
            sublimation[day] = flux.sublimation;
            snow_adjustment[day] = flux.snowfall_adjustment;
            const Scalar total = std::inner_product(
                snow.begin(), snow.end(), area.begin(), Scalar(0.0));
            snow_change[day] = total - previous_snow;
            previous_snow = total;
        }
    }

    for (std::size_t day = 0; day < spin; ++day) {
        auto flux = uniform_surface(surface, rain[day], demand[day],
            p.surface_capacity, surface_drainage);
        pre_spin_rain[day] = flux.overflow + flux.drainage;
        pre_spin_demand[day] = maximum(demand[day] - flux.evaporation, 0.0);
    }
    if (fixed_spinup) surface = fixed_spinup->surface;
    if (captured_spinup) captured_spinup->surface = primal(surface);
    CART_INITIAL(0,surface);
    for (std::size_t day = 0; day < n; ++day) {
        const Scalar before = surface;
        auto flux = uniform_surface(surface, rain[day], demand[day],
            p.surface_capacity, surface_drainage);
        rain[day] = flux.overflow + flux.drainage;
        CART_TRACE(day, SURFACE_END, surface); CART_TRACE(day, SURFACE_OUT,
            rain[day]);
        CART_TRACE(day, SURFACE_ET, flux.evaporation);
        demand[day] = maximum(demand[day] - flux.evaporation, 0.0);
        if constexpr (!DischargeOnly) {
            surface_evap[day] = flux.evaporation;
            surface_change[day] = surface - before;
        }
    }
    const auto infiltration = [&](auto& input, const auto& pet, auto& runoff
        ) {
        Scalar memory = 0.0;
        Scalar initial = mean_rain;
        // Keep the unpartitioned input for the finite causal memory.
        const auto raw = input;
        for (std::size_t day = 0; day < input.size(); ++day) {
            const Scalar capacity = dt * maximum(p.infiltration_capacity *
                exp(
                -p.infiltration_wetness * (raw[day]+(memory+initial-raw[day]
                    )*infiltration_mean) / mean_rain),
                infiltration_floor);
            const Scalar available = maximum(raw[day] - pet[day], 0.0);
            const Scalar x = available / capacity;
            runoff[day] = available - capacity * (x < 1e-3
                ? x - 0.5 * x * x : -expm1(-x));
            input[day] -= runoff[day];
            memory = infiltration_decay * memory
                + (1.0 - infiltration_decay) * raw[day];
            if (day >= infiltration_window) memory -= infiltration_tail *
                raw[day-infiltration_window];
            initial *= infiltration_decay;
        }
    };
    std::vector<Scalar> spin_excess(spin);
    infiltration(pre_spin_rain, pre_spin_demand, spin_excess);
    infiltration(rain, demand, excess);
    for (std::size_t i=0;i<n;++i) {
        CART_TRACE(i, INFILTRATION_EXCESS, excess[i]); CART_TRACE(i,
            SOIL_INPUT, rain[i]);
    }

    Scalar soil = 0.7 * p.soil_capacity, zone = 0.3 * p.zone_capacity;
    Scalar wet = 0.3 * p.wetness_capacity;
    for (std::size_t day = 0; day < spin; ++day) {
        uniform_pareto<Scalar>(soil, pre_spin_rain[day], pre_spin_demand[day
            ],
            p.soil_capacity, p.soil_percolation, p.soil_shape, 2.0, true, dt
                , true);
        zone_step(zone, pre_spin_rain[day], pre_spin_demand[day],
            p.zone_capacity, p.zone_days);
        zone_step(wet, pre_spin_rain[day], pre_spin_demand[day],
            p.wetness_capacity, p.wetness_days);
    }
    if (fixed_spinup) {
        soil = fixed_spinup->soil;
        zone = fixed_spinup->zone;
        wet = fixed_spinup->wetness;
    }
    if (captured_spinup) {
        captured_spinup->soil = primal(soil);
        captured_spinup->zone = primal(zone);
        captured_spinup->wetness = primal(wet);
    }
    CART_INITIAL(1,soil); CART_INITIAL(2,zone);
    const Scalar upland = 1.0 - p.zone_fraction;
    for (std::size_t day = 0; day < n; ++day) {
        const Scalar before = upland * soil + p.zone_fraction * zone;
        const Scalar soil_before=soil, wet_before=wet;
        auto soil_flux = uniform_pareto<Scalar>(soil, rain[day], demand[day]
            ,
            p.soil_capacity, p.soil_percolation, p.soil_shape, 2.0, true, dt
                );
        auto zone_flux = zone_step(zone, rain[day], demand[day],
            p.zone_capacity, p.zone_days);
        zone_step(wet, rain[day], demand[day],
            p.wetness_capacity, p.wetness_days);
        soil_fraction[day] = bounded(0.5*(soil_before+soil) / p.
            soil_capacity, 0.0, 1.0);
        wetness[day] = 0.5*(wet_before+wet) / p.wetness_capacity;
        fast[day] = upland * soil_flux.runoff + p.zone_fraction
            * zone_flux.runoff + excess[day];
        recharge[day] = upland * soil_flux.recharge;
        evaporation[day] = upland * soil_flux.evaporation + p.zone_fraction
            * zone_flux.evaporation;
        zone_flow[day] = p.zone_fraction * zone_flux.discharge;
        CART_TRACE(day, SOIL_INTERCEPTION, minimum(rain[day],demand[day]));
        CART_TRACE(day, NEAR_INTERCEPTION, minimum(rain[day],demand[day]));
        CART_TRACE(day, SOIL_END, soil); CART_TRACE(day, SOIL_RUNOFF,
            soil_flux.runoff);
        CART_TRACE(day, SOIL_RECHARGE, soil_flux.recharge); CART_TRACE(day,
            SOIL_ET, soil_flux.evaporation);
        CART_TRACE(day, NEAR_END, zone); CART_TRACE(day, NEAR_RUNOFF,
            zone_flux.runoff);
        CART_TRACE(day, NEAR_Q, zone_flux.discharge); CART_TRACE(day,
            NEAR_ET, zone_flux.evaporation);
        if constexpr (!DischargeOnly) {
            storage_change[day] = upland * soil + p.zone_fraction * zone -
                before;
        }
    }
    Scalar depressions = 0.5 * p.depression_capacity;
    for (std::size_t day = 0; day < spin; ++day) {
        uniform_pareto<Scalar>(depressions, fast[day],
            maximum(demand[day] - evaporation[day], 0.0),
            p.depression_capacity, p.depression_percolation,
            p.depression_shape, 2.0, false, dt, true);
    }
    if (fixed_spinup) depressions = fixed_spinup->depressions;
    if (captured_spinup) captured_spinup->depressions = primal(depressions);
    CART_INITIAL(4,depressions);
    Scalar training_recharge = 0.0;
    for (std::size_t day = 0; day < n; ++day) {
        const Scalar before = depressions;
        auto flux = uniform_pareto<Scalar>(depressions, fast[day],
            maximum(demand[day] - evaporation[day], 0.0),
            p.depression_capacity, p.depression_percolation,
            p.depression_shape, 2.0, false, dt);
        fast[day] = flux.runoff;
        CART_TRACE(day, DEP_END, depressions); CART_TRACE(day, DEP_RUNOFF,
            flux.runoff);
        CART_TRACE(day, DEP_RECHARGE, flux.recharge); CART_TRACE(day, DEP_ET
            , flux.evaporation);
        recharge[day] += flux.recharge;
        evaporation[day] += flux.evaporation;
        if constexpr (!DischargeOnly) {
            storage_change[day] += depressions - before;
        }
        if (forcing[day].training) training_recharge += recharge[day];
    }
    training_recharge = maximum(training_recharge / (training_days * dt),
        1e-6);
    const Scalar equilibrium_storage = training_recharge * p.
        groundwater_days;
    const Scalar groundwater_start = fixed_spinup
        ? Scalar(fixed_spinup->groundwater) : equilibrium_storage;
    if (captured_spinup) captured_spinup->groundwater = primal(
        groundwater_start);
    CART_INITIAL(5,groundwater_start);
    const Scalar groundwater_rate =
        dt / (equilibrium_storage * p.groundwater_days);
    Scalar groundwater = groundwater_start;
    const bool skip_zero_inputs = DischargeOnly && dt != 1.0;
    UniformHydrograph<6, Scalar> uh1(p.routing_days, dt, skip_zero_inputs);
    UniformHydrograph<13, Scalar> uh2(p.routing_days, dt, skip_zero_inputs);
    for (std::size_t day = 0; day < n; ++day) {
        const Scalar before = groundwater;
        const auto flux = uniform_groundwater(groundwater, recharge[day],
            maximum(demand[day] - evaporation[day], 0.0), groundwater_rate,
            p.groundwater_uptake, groundwater_start);
        groundwater_flow[day] = flux.discharge;
        CART_TRACE(day, GW_END, groundwater); CART_TRACE(day, GW_Q, flux.
            discharge);
        CART_TRACE(day, GW_ET, flux.uptake);
        if constexpr (!DischargeOnly) {
            evaporation[day] += flux.uptake;
            storage_change[day] += groundwater - before;
        }
        const Scalar connected = p.connected_minimum
            + (1.0 - p.connected_minimum) * logistic((wetness[day]
                - p.connectivity_center + p.soil_connectivity *
                    soil_fraction[day]
                + p.groundwater_connectivity * (0.5*(before+groundwater) /
                    groundwater_start - 1.0))
                / p.connectivity_width);
        const Scalar fast_bypass = p.bypass_fraction * connected * fast[day]
            ;
        const Scalar routed_in = fast[day] - fast_bypass;
        const Scalar routed_out = uh1.advance(routed_in);
        const Scalar bypass_out = uh2.advance(fast_bypass);
        routed[day] = routed_out + 0.9 * (dt * p.exterior_inflow);
        bypass[day] = bypass_out + 0.1 * (dt * p.exterior_inflow);
        CART_TRACE(day, ROUTE_IN, routed_in); CART_TRACE(day, BYPASS_IN,
            fast_bypass);
        CART_TRACE(day, ROUTE_OUT, routed_out); CART_TRACE(day, BYPASS_OUT,
            bypass_out);
        CART_TRACE(day, EXTERIOR, dt*p.exterior_inflow);
        if constexpr (!DischargeOnly) {
            routing_change[day] = uh1.storage_change(routed_in, routed_out)
                + uh2.storage_change(fast_bypass, bypass_out);
        }
        const Scalar temperature = forcing[day].temperature_c;
        const Scalar kelvin = maximum(temperature, 0.0) + 273.15;
        const Scalar viscosity = -log(10.0) * 247.8
            * (1.0 / (kelvin - 140.0) - 1.0 / 153.15);
        connectivity[day] = log(connected) + viscosity - p.frost_reduction
            * logistic((p.frost_temperature - temperature) / p.frost_width);
    }
    Scalar quick = 0.5 * p.quick_capacity;
    const Scalar release = 0.25 * dt * pow(p.quick_capacity, -4.0);
    const Scalar exchange = dt * p.exchange * pow(30.0 / p.quick_capacity,
        3.5);
    const auto quick_step = [&](std::size_t day, Scalar* direct=nullptr) {
        return uniform_quick(quick, routed[day], bypass[day], release,
            exchange, p.quick_capacity, connectivity[day], p.quick_exponent,
                direct);
    };
    for (std::size_t day = 0; day < quick_spin; ++day) quick_step(day);
    if (fixed_spinup) quick = fixed_spinup->quick;
    if (captured_spinup) captured_spinup->quick = primal(quick);
    CART_INITIAL(6,quick);
    std::vector<CatchmentDayT<Scalar>> output;
    output.reserve(diagnostics);
    for (std::size_t day = 0; day < n; ++day) {
        const Scalar before = quick;
        Scalar direct=0;
        const auto flux = quick_step(day,&direct);
        const Scalar discharge = flux.discharge + groundwater_flow[day] +
            zone_flow[day];
        const Scalar x = maximum(discharge, 0.0) / (dt * p.
            channel_loss_capacity);
        const Scalar loss = dt * p.channel_loss_capacity * (x < 1e-3
            ? x - 0.5 * x * x : -expm1(-x));
        const Scalar flow = discharge - loss;
        CART_TRACE(day, QUICK_END, quick); CART_TRACE(day, QUICK_Q, flux.
            discharge-direct);
        CART_TRACE(day, DIRECT_Q, direct); CART_TRACE(day, CHANNEL_LOSS,
            loss);
        CART_TRACE(day, NET_EXCHANGE, flux.exchange);
        if constexpr (DischargeOnly) {
            std::copy(flow.slope.begin(), flow.slope.end(),
                      discharge_jacobian + day * parameter_count);
        } else {
            output.push_back({
                flow,
                evaporation[day] + surface_evap[day] + sublimation[day],
                snow_adjustment[day] + (dt * p.exterior_inflow) + flux.
                    exchange - loss,
                snow_change[day] + surface_change[day] + storage_change[day]
                    + routing_change[day] + quick - before
            });
        }
    }
    return output;
}

}  // namespace cartplus

