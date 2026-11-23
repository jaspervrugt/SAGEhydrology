# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/cartplus_ode.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Experimental continuous-store CaRT+ODE model.

Forcing rows contain precipitation depth, temperature, PET depth and a
training flag. Model rates and solver steps use days. The native CaRT+
parameter order is unchanged. Unit hydrographs remain discrete at forcing
boundaries, between the catchment and quick-store ODEs.
"""

import ctypes
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

from cartplus_python import hydrograph


@dataclass(frozen=True)
class Solver:
    'Integration controls in days, ' 'with depth tolerances in mm.'

    initial_step: float = .01
    maximum_step: float = .25
    relative_tolerance: float = 1e-6
    absolute_tolerance: float = 1e-8


@lru_cache(maxsize=1)
def native_entrypoint():
    """
    Bind the compiled CaRT+ODE forward and tangent entry point.

    SYNOPSIS:
        native_entrypoint(...)

    INPUT ARGUMENTS:
    None.

    OUTPUT ARGUMENTS:
    Native forward/tangent callable with the validated ABI.
    """
    library = Path(__file__).with_name("cartplus_ode_cpu.dll")
    function = ctypes.CDLL(str(library)).cartplus_ode_run
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    pointer = ctypes.POINTER(ctypes.c_double)
    function.argtypes = [
        ctypes.c_size_t, array, array, array, array, ctypes.c_double,
        ctypes.c_double, array, array, pointer,
    ]
    function.restype = ctypes.c_int
    return function


def simulate_native(
    forcing, elevation, area, parameters, resolution=1, *,
    solver=Solver(), infiltration_floor=0., jacobian=False,
):
    """
    Evaluate the native model and its optional Jacobian.

    SYNOPSIS:
        simulate_native(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    parameters : 42 physical CaRT+ parameters in native order.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    solver : Integration method or solver-control object.
    infiltration_floor : Minimum infiltration capacity (mm/day).
    jacobian : Whether to evaluate discharge sensitivities.

    OUTPUT ARGUMENTS:
    output : N-by-4 interval water-balance fluxes (mm); J : optional
    discharge Jacobian in the selected coordinates.

    NOTES:
    Forward sensitivities include parameter-dependent initial stores and the
    pilot recharge pass. The solver controls errors in both states and
    tangents. Step-size decisions are excluded from differentiation.
    """
    f, z, a, p = (
        np.ascontiguousarray(x, dtype=np.float64)
        for x in (forcing, elevation, area, parameters)
    )
    if resolution not in (1, 24, 96) or f.ndim != 2 or f.shape[1] != 4:
        raise ValueError('Use daily, hourly or 15-minute' ' CaRT+ forcing.')
    if p.shape != (42,) or z.shape != (5,) or a.shape != (5,):
        raise ValueError('Expected 42 parameters and fiv'
            'e elevation bands.')
    if not all(np.isfinite(x).all() for x in (f, z, a, p)):
        raise ValueError("CaRT+ODE inputs must be finite.")
    if not f[:, 3].any() or np.any(f[:, [0, 2]] < 0):
        raise ValueError('Training intervals and nonnega'
            'tive P/PET needed.')
    settings = np.array([
        solver.initial_step, solver.maximum_step,
        solver.relative_tolerance, solver.absolute_tolerance,
    ])
    if not np.isfinite(settings).all() or np.any(settings <= 0):
        raise ValueError("CaRT+ODE solver controls must be positive.")
    result = np.empty((len(f), 4))
    derivative = np.empty((len(f), 42)) if jacobian else None
    pointer = (
        derivative.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        if jacobian else None
    )
    status = native_entrypoint()(
        len(f), f, z, a, p, 1/resolution, infiltration_floor,
        settings, result, pointer,
    )
    if status or not np.isfinite(result).all() or (
        jacobian and not np.isfinite(derivative).all()
    ):
        raise ValueError("CaRT+ODE native integration failed.")
    return result, derivative


def interval(rhs, state, dt, solver):
    """
    Integrate to a forcing boundary and return its terminal state.

    SYNOPSIS:
        interval(...)

    INPUT ARGUMENTS:
    rhs : Callback evaluating state derivatives on the model time base.
    state : State vector (including tangent lanes when enabled).
    dt : Integration or observation interval duration (days).
    solver : Integration method or solver-control object.

    OUTPUT ARGUMENTS:
    State at the end of the requested integration duration.
    """
    result = solve_ivp(
        rhs, (0, dt), state, method="RK45", t_eval=[dt],
        first_step=min(dt, solver.initial_step, solver.maximum_step),
        max_step=min(dt, solver.maximum_step),
        rtol=solver.relative_tolerance, atol=solver.absolute_tolerance,
    )
    if not result.success or not np.isfinite(result.y).all():
        raise ValueError('CaRT+ODE failed to integrate a'
            ' forcing interval.')
    return result.y[:, -1]


def production(store, rain, pet, capacity, shape, percolation, net):
    """
    Evaluate production-store runoff, evaporation, and recharge.

    SYNOPSIS:
        production(...)

    INPUT ARGUMENTS:
    store : Current water store (mm).
    rain : Liquid input depth/rate on this operator time base.
    pet : Potential evaporation depth/rate on this operator time base.
    capacity : Store or infiltration capacity in this operator.
    shape : Dimensionless storage/infiltration shape parameter.
    percolation : Percolation coefficient in the model parameterization.
    net : Whether to remove interception before the production update.

    OUTPUT ARGUMENTS:
    Runoff, evaporation, and drainage rates on the input time base.
    """
    wetness = np.clip(store / capacity, 0, 1)
    deficit = 1 - wetness
    interception = min(rain, pet) if net else 0.
    incoming = rain - interception
    exponent = shape/(1+shape)
    threshold = .01
    retained_fraction = (
        deficit/threshold*threshold**exponent if deficit < threshold
        else deficit**exponent
    )
    runoff = incoming * (1-retained_fraction)
    evap = interception + (pet-interception) * (1-deficit**2)
    drain = max(store, 0) * np.log1p(percolation*wetness**4) / 4
    return runoff, evap, drain


def zone(store, rain, pet, capacity, days):
    """
    Evaluate near-stream runoff, evaporation, and drainage.

    SYNOPSIS:
        zone(...)

    INPUT ARGUMENTS:
    store : Current water store (mm).
    rain : Liquid input depth/rate on this operator time base.
    pet : Potential evaporation depth/rate on this operator time base.
    capacity : Store or infiltration capacity in this operator.
    days : Residence time (days).

    OUTPUT ARGUMENTS:
    Near-stream runoff, evaporation, and drainage rates.
    """
    s = np.clip(store / capacity, 0, 1)
    intercepted = min(rain, pet)
    runoff = (rain-intercepted)*s*s
    evap = intercepted+(pet-intercepted)*s*(2-s)
    flow = max(store, 0)/days
    return runoff, evap, flow


def simulate(forcing, elevation, area, parameters, resolution=1, *,
             solver=Solver(), infiltration_floor=0.):
    """
    Simulate model fluxes and optional discharge sensitivities.

    SYNOPSIS:
        simulate(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    parameters : 42 physical CaRT+ parameters in native order.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    solver : Integration method or solver-control object.
    infiltration_floor : Minimum infiltration capacity (mm/day).

    OUTPUT ARGUMENTS:
    output : N-by-4 [Q, ET, external flux, storage change] (mm/interval).

    NOTES:
    The prototype uses a short availability timescale for snow and
    interception and a 1% surface-capacity transition for overflow. Initial
    stores are specified as capacity fractions. A pilot pass estimates the
    groundwater equilibrium from training recharge.
    """
    f, elevation, area, p = (
        np.asarray(x, dtype=float)
        for x in (forcing, elevation, area, parameters)
    )
    if resolution not in (1, 24, 96) or f.ndim != 2 or f.shape[1] != 4:
        raise ValueError('Use daily, hourly or 15-minute' ' CaRT+ forcing.')
    if len(p) != 42 or elevation.shape != (5,) or area.shape != (5,):
        raise ValueError('Expected 42 parameters and fiv'
            'e elevation bands.')
    if not all(np.isfinite(x).all() for x in (f, elevation, area, p)):
        raise ValueError("CaRT+ODE inputs must be finite.")
    train = f[:, 3].astype(bool)
    if not train.any() or np.any(f[:, [0, 2]] < 0):
        raise ValueError('Training intervals and nonnega'
            'tive P/PET needed.')
    dt = 1/resolution
    (threshold, melt, snow_factor, sublimation, pet_melt, lapse,
     surface_cap, surface_days, soil_cap, soil_perc, soil_shape,
     infiltration_cap, infiltration_shape, zone_fraction, zone_cap,
     zone_days, depression_cap, depression_shape, depression_perc,
     groundwater_days, uptake, vegetation_reduction, vegetation_temp,
     vegetation_width, vegetation_days, wet_cap, wet_days, connected_min,
     connected_center, connected_width, soil_connection,
     groundwater_connection, bypass_fraction, routing_days, quick_cap,
     exchange, exterior, frost_temp, frost_width, frost_reduction,
     quick_exponent, channel_cap) = p
    upland = 1-zone_fraction
    band_offset = lapse*(elevation-area@elevation)/1000
    rain_mean = max(f[train, 0].mean()/dt, .1)
    veg_forcing = expit((f[:, 1]-vegetation_temp)/vegetation_width)
    veg_start = veg_forcing[:365*resolution].mean()
    initial = np.r_[np.zeros(5), .5*surface_cap, .7*soil_cap,
                    .3*zone_cap, .3*wet_cap, .5*depression_cap,
                    rain_mean, veg_start]
    availability_days = .02

    def rates(state, precip, warmth, pet, vegetation):
        """
        Evaluate state and cumulative-flux derivatives.

        SYNOPSIS:
            rates(...)

        INPUT ARGUMENTS:
        state : State vector (including tangent lanes when enabled).
        precip : Precipitation depth/rate on the selected time base.
        warmth : Band temperature above the snow threshold (deg C).
        pet : Potential evaporation depth/rate on this operator time base.
        vegetation : Vegetation used by this numerical operator.

        OUTPUT ARGUMENTS:
        Computed fluxes, stores, or numerical values described in the notes.
        """
        snow = np.maximum(state[:5], 0)
        surface, soil, near, wet, depression, memory, veg = state[5:12]
        cold = warmth < 0
        snowfall = cold*snow_factor*precip
        demand = pet*(1-vegetation_reduction*(1-np.clip(veg, 0, 1)))
        snow_et = np.minimum(sublimation*pet, snow/availability_days)
        available = np.maximum(snow/availability_days-snow_et, 0)
        snowmelt = np.minimum(
            available, np.where(cold, 0, melt*warmth+pet_melt*pet),
        )
        liquid = area @ (snowmelt+(~cold)*precip)
        surface_et = min(demand, max(surface, 0)/availability_days)
        surface_q = max(surface, 0)/surface_days
        excess_fraction = np.clip(
            (surface/surface_cap-.99)/.01, 0, 1,
        )
        overflow = max(liquid-surface_et-surface_q, 0)*excess_fraction
        input_rate = overflow+surface_q
        potential = max(demand-surface_et, 0)
        cap = max(infiltration_cap*np.exp(
            -infiltration_shape*max(memory, 0)/rain_mean,
        ), infiltration_floor, np.finfo(float).tiny)
        available_rain = max(input_rate-potential, 0)
        excess = available_rain+cap*np.expm1(-available_rain/cap)
        rain = input_rate-excess
        soil_q, soil_et, soil_drain = production(
            soil, rain, potential, soil_cap, soil_shape, soil_perc, True,
        )
        near_q, near_et, near_flow = zone(
            near, rain, potential, zone_cap, zone_days,
        )
        wet_q, wet_et, wet_flow = zone(
            wet, rain, potential, wet_cap, wet_days,
        )
        et = upland*soil_et+zone_fraction*near_et
        fast = upland*soil_q+zone_fraction*near_q+excess
        dep_q, dep_et, dep_drain = production(
            depression, fast, max(potential-et, 0), depression_cap,
            depression_shape, depression_perc, False,
        )
        recharge = upland*soil_drain+dep_drain
        remaining_pet = max(potential-et-dep_et, 0)
        return np.r_[
            snowfall-snow_et-snowmelt,
            liquid-surface_et-input_rate,
            rain-soil_q-soil_et-soil_drain,
            rain-near_q-near_et-near_flow,
            rain-wet_q-wet_et-wet_flow,
            fast-dep_q-dep_et-dep_drain,
            (input_rate-memory)/3,
            (vegetation-veg)/vegetation_days,
        ], (dep_q, recharge, zone_fraction*near_flow, remaining_pet,
            et+dep_et+surface_et+area@snow_et,
            (area@cold)*precip*(snow_factor-1))

    def run(reference=None):
        """
        Advance the pilot or full model over the forcing sequence.

        SYNOPSIS:
            run(...)

        INPUT ARGUMENTS:
        reference : Pilot mean recharge or reference storage used by the
            caller.

        OUTPUT ARGUMENTS:
        (fluxes, storage_change) : computed result in the operand units.
        """
        state = initial.copy()
        if reference is not None:
            state = np.r_[state, reference]
        fluxes = np.zeros((len(f), 8))
        storage_change = np.zeros(len(f))
        storage_weights = np.r_[area, 1, upland, zone_fraction, 0,
                                1, 0, 0, 1]
        for i, (precip, temp, pet, _) in enumerate(f):
            warmth = temp-band_offset-threshold
            active = reference is not None

            def rhs(_, values):
                """
                Evaluate augmented water-state and sensitivity derivatives.

                SYNOPSIS:
                    rhs(...)

                INPUT ARGUMENTS:
                _ : Reserved/unused callback argument.
                values : Numerical values or vector of state entries.

                OUTPUT ARGUMENTS:
                Computed fluxes, stores, or numerical values described in
                the notes.
                """
                change, terms = rates(
                    values, precip/dt, warmth, pet/dt, veg_forcing[i],
                )
                fast, recharge, near_flow, demand, evap, adjustment = terms
                if not active:
                    return np.r_[change, recharge]
                groundwater = max(values[12], 0)
                et_gw = groundwater*uptake*demand/(groundwater+reference)
                base = groundwater**2/(reference*groundwater_days)
                signal = values[8]/wet_cap-connected_center
                signal += soil_connection*np.clip(values[6]/soil_cap, 0, 1)
                signal += groundwater_connection*(groundwater/reference-1)
                connected = connected_min+(1-connected_min)*expit(
                    signal/connected_width,
                )
                bypass = bypass_fraction*connected*fast
                return np.r_[change, recharge-base-et_gw, fast-bypass,
                             bypass, near_flow, base, evap+et_gw,
                             adjustment, connected, recharge]

            end = interval(rhs, np.r_[state, np.zeros(8 if active else 1)],
                           dt, solver)
            if active:
                fluxes[i] = end[len(state):]
                storage_change[i] = storage_weights @ (
                    end[:len(state)] - state
                )
            else:
                fluxes[i, 7] = end[-1]
            state = end[:len(state)]
        return fluxes, storage_change

    pilot, _ = run()
    reference = max(pilot[train, 7].mean()/dt, 1e-6)*groundwater_days
    catchment, storage_change = run(reference)
    routed = hydrograph(catchment[:, 0], routing_days, dt, False)
    bypass = hydrograph(catchment[:, 1], routing_days, dt, True)
    quick = .5*quick_cap
    output = np.empty((len(f), 4))
    for i, (_, temp, _, _) in enumerate(f):
        in_quick = routed[i]/dt+.9*exterior
        in_bypass = bypass[i]/dt+.1*exterior
        kelvin = max(temp, 0)+273.15
        viscosity = -np.log(10)*247.8*(1/(kelvin-140)-1/153.15)
        signal = np.log(max(catchment[i, 6]/dt, 1e-12))+viscosity
        signal -= frost_reduction*expit((frost_temp-temp)/frost_width)

        def quick_rhs(_, values):
            """
            Advance quick storage, runoff, and net external exchange.

            SYNOPSIS:
                quick_rhs(...)

            INPUT ARGUMENTS:
            _ : Reserved/unused callback argument.
            values : Numerical values or vector of state entries.

            OUTPUT ARGUMENTS:
            (incoming - release, release + diverted, incoming - in_quick +
            diverted - in_bypass) : computed result in the operand units.
            """
            store = max(values[0], 0)
            flux = exchange*(store/quick_cap)**3.5
            diverted = max(in_bypass+.1*flux, 0)
            incoming = in_quick+.9*flux
            if store <= 0:
                incoming = max(incoming, 0)
            multiplier = signal+quick_exponent*np.log(
                max(store, 1e-3)/quick_cap,
            )
            release = .25*store*(store/quick_cap)**4*np.exp(
                6*np.tanh(multiplier/6),
            )
            return incoming-release, release+diverted, (
                incoming-in_quick+diverted-in_bypass
            )

        previous_quick = quick
        quick, quick_q, exchange_depth = interval(
            quick_rhs, [quick, 0, 0], dt, solver,
        )
        q = quick_q+catchment[i, 2]+catchment[i, 3]
        loss = -dt*channel_cap*np.expm1(-q/(dt*channel_cap))
        evap = catchment[i, 4]
        extra = catchment[i, 5]+dt*exterior+exchange_depth-loss
        routing_change = catchment[i, :2].sum()-routed[i]-bypass[i]
        change = storage_change[i]+routing_change+quick-previous_quick
        output[i] = q-loss, evap, extra, change
    return output

