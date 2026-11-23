# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/cart/cart_python.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Daily Capacityâ€“Residenceâ€“Time baseline after Xia (2026).

Adapted from the MIT-licensed xiaxilin/cartplus reference implementation.
The input parameter order follows its CaRT static parameter table.
"""

from ctypes import CDLL, POINTER, c_double, c_int
from functools import lru_cache
from pathlib import Path
import sys

import numpy as np

# Reuse the same forward-mode scalar as the CaRT+ implementation.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "CaRTPlus"))
try:
    from cartplus_python import Tangent, primal
finally:
    sys.path.pop(0)


def simulate(forcing, elevation, area, parameters, *, jacobian=False,
             initial=None, return_initial=False, mcode=1, resolution=1,
             uniform_input=False):
    """
    Simulate model fluxes and optional discharge sensitivities.

    SYNOPSIS:
        simulate(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    parameters : 11 CaRT calibration coordinates, including logarithmic
        entries.
    jacobian : Whether to evaluate discharge sensitivities.
    initial : Initial stores (mm), or None/empty to use the spin-up policy.
    return_initial : Whether to append detached spin-up stores to the return
        value.
    mcode : Backend/solver case selected by the model wrapper.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    uniform_input : Use uniform-within-interval inputs instead of step-start
        inputs.

    OUTPUT ARGUMENTS:
    q : N interval discharge depths (mm); J : N-by-11 Jacobian or None.
    Optional stores/diagnostics are appended only when requested.

    NOTES:
    P and PET are depths per interval and T is in degrees Celsius. The
    fourth column marks the training period. Parameters are in Xia's
    calibration space, including the log-transformed entries. Supply the
    returned initial stores when checking derivatives by FD.
    uniform_input=True selects Xia hourly revision f986af2, including at
    daily resolution; False retains the published step-start formulation.
    """
    if resolution not in (1, 24, 96):
        raise ValueError('CaRT supports daily, hourly an'
            'd 15-minute data.')
    dt = 1 / resolution
    forcing = np.asarray(forcing, dtype=float)
    elevation = np.asarray(elevation, dtype=float)
    area = np.asarray(area, dtype=float)
    p = np.asarray(parameters, dtype=float)
    if (forcing.ndim != 2 or forcing.shape[1] != 4
            or p.shape != (11,) or elevation.shape != (5,)
            or area.shape != (5,) or not np.isfinite(forcing).all()):
        raise ValueError(
            'Expected interval forcing, 11 ' 'parameters and 5 bands.'
        )
    if initial is not None and np.shape(initial) != (3,):
        raise ValueError("CaRT requires three initial store values.")
    train = forcing[:, 3] != 0
    if not train.any():
        raise ValueError('CaRT requires a training-perio'
            'd recharge mask.')
    if uniform_input and mcode != 4:
        raise ValueError('Uniform-input CaRT requires th' 'e C++ backend.')
    if mcode == 4:
        q, derivative, stores = _native(
            forcing, elevation, area, p, jacobian, initial, resolution,
            uniform_input,
        )
        return (q, derivative, stores) if return_initial else (
            q, derivative
        )
    if mcode not in (1, 2, 3):
        raise ValueError("Unknown CaRT solver case.")
    if initial is None:
        _, initial = _run(
            forcing, elevation, area, p, train, None, dt, resolution, mcode
        )
    if jacobian:
        p = [Tangent(value, seed) for value, seed in zip(p, np.eye(11))]
    q, _ = _run(
        forcing, elevation, area, p, train, initial, dt, resolution, mcode
    )
    discharge = np.array([primal(value) for value in q])
    derivative = (
        np.array([value.slope for value in q]) if jacobian else None
    )
    return (discharge, derivative, initial) if return_initial else (
        discharge, derivative
    )


def _run(forcing, elevation, area, p, train, initial, dt,
         resolution, mcode):
    """
    Advance the model stores, hydrographs, and discharge recursion.

    SYNOPSIS:
        _run(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    p : 11 CaRT calibration coordinates, including logarithmic entries.
    train : Boolean training-interval mask.
    initial : Initial stores (mm), or None/empty to use the spin-up policy.
    dt : Integration or observation interval duration (days).
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    mcode : Backend/solver case selected by the model wrapper.

    OUTPUT ARGUMENTS:
    (q, (s0, r0, w0)) : computed result in the operand units.
    """
    lc, x4, lx3, ltau, x2, rkp, tt, lcf, sf, lgam, lks = p
    c, x3, tau = np.exp(lc), np.exp(lx3), np.exp(ltau)
    kp = (4/9)**4 * np.exp(rkp)
    cf, gamma, ks = np.exp(lcf), np.exp(lgam), np.exp(lks)
    snow = [0.] * 5
    offsets = -gamma * (elevation - area @ elevation) / 1000
    liquid = []
    for rain, temp, pet, _ in forcing:
        total = 0.
        for k, offset in enumerate(offsets):
            warm = temp + offset - tt
            cold = warm < 0
            if cold:
                snow[k] += sf * rain
            snow[k] -= min(ks * pet, snow[k])
            melt = min(max(cf * dt * warm, 0.), snow[k])
            snow[k] -= melt
            total += area[k] * (melt if cold else rain + melt)
        liquid.append(total)

    def production(store, count):
        """
        Evaluate production-store runoff, evaporation, and recharge.

        SYNOPSIS:
            production(...)

        INPUT ARGUMENTS:
        store : Current water store (mm).
        count : Number of intervals or state entries in this helper.

        OUTPUT ARGUMENTS:
        Runoff, evaporation, and drainage rates on the input time base.
        """
        direct, recharge = [], []
        for rain, pet in zip(liquid[:count], forcing[:count, 2]):
            pnet, enet = max(rain-pet, 0.), max(pet-rain, 0.)
            s = store / c
            tp, te = np.tanh(pnet/c), np.tanh(enet/c)
            add = c*(1-s*s)*tp/(1+s*tp)
            evap = store*(2-s)*te/(1+(1-s)*te)
            available = store - evap + add
            drain = available * (1-(1+dt*kp*(available/c)**4)**-.25)
            store = available - drain
            direct.append(pnet-add)
            recharge.append(drain)
        return direct, recharge, store

    s0 = (
        production(.7*c, 1460*resolution)[2]
        if initial is None else initial[0]
    )
    direct, recharge, _ = production(s0, len(forcing))
    x9 = .9 * _route(direct, x4, False, resolution)
    x1 = .1 * _route(direct, x4, True, resolution)

    def quick(store, count):
        """
        Advance the quick store and direct routing branch.

        SYNOPSIS:
            quick(...)

        INPUT ARGUMENTS:
        store : Current water store (mm).
        count : Number of intervals or state entries in this helper.

        OUTPUT ARGUMENTS:
        (result, store) : computed result in the operand units.
        """
        result = []
        for main, bypass in zip(x9[:count], x1[:count]):
            loss = dt * x2 * (store/x3)**3.5
            available = max(store + main + .9*loss, 0.)
            store = available * (1+dt*(available/x3)**4)**-.25
            result.append(available-store+max(bypass+.1*loss, 0.))
        return result, store

    r0 = quick(.5*x3, 730*resolution)[1] if initial is None else initial[1]
    routed = quick(r0, len(forcing))[0]
    gbar = max(sum(g for g, keep in zip(recharge, train) if keep)
               / (np.count_nonzero(train)*dt), 1e-6)
    coefficient = 1/(gbar*tau*tau)
    w0 = gbar*tau if initial is None else initial[2]
    water = w0
    q = []
    for fast, g in zip(routed, recharge):
        available = water + g
        water = _recession(available, coefficient, dt, mcode)
        q.append(fast + available - water)
    return q, (s0, r0, w0)


def _recession(water, coefficient, dt, mcode):
    """
    Integrate nonlinear groundwater recession over one step.

    SYNOPSIS:
        _recession(...)

    INPUT ARGUMENTS:
    water : Available water or store depth (mm).
    coefficient : Release-law coefficient in the current interval time base.
    dt : Integration or observation interval duration (days).
    mcode : Backend/solver case selected by the model wrapper.

    OUTPUT ARGUMENTS:
    End-of-step groundwater storage (mm), preserving derivatives.
    """
    if mcode in (1, 3):
        steps = max(8, int(np.ceil(8*dt*primal(coefficient)*primal(water))))
        step = dt/steps
        for _ in range(steps):
            stage = water
            if mcode == 1:
                stage = water - .5*step*coefficient*water*water
            water -= step*coefficient*stage*stage
        return water
    from scipy.integrate import solve_ivp

    k = primal(coefficient)
    if isinstance(water, Tangent):
        dk = coefficient.slope
        def rhs(_, state):
            """
            Evaluate augmented water-state and sensitivity derivatives.

            SYNOPSIS:
                rhs(...)

            INPUT ARGUMENTS:
            _ : Reserved/unused callback argument.
            state : State vector (including tangent lanes when enabled).

            OUTPUT ARGUMENTS:
            np.r_[-k * w * w, -2 * k * w * state[1:] - dk * w * w] :
            computed result in the operand units.
            """
            w = state[0]
            return np.r_[-k*w*w, -2*k*w*state[1:] - dk*w*w]

        solution = solve_ivp(rhs, (0, dt), np.r_[water.value, water.slope],
                             rtol=1e-9, atol=1e-11)
    else:
        solution = solve_ivp(lambda _, w: -k*w*w, (0, dt), [water],
                             rtol=1e-9, atol=1e-11)
    if not solution.success:
        raise ValueError(
            'CaRT groundwater integration failed:' ' ' f'{solution.message}'
        )
    end = solution.y[:, -1]
    if isinstance(water, Tangent):
        return Tangent(end[0], end[1:])
    return end[0]


def _route(inflow, time_base, symmetric, resolution):
    """
    Convolve runoff with the selected unit-hydrograph ordinates.

    SYNOPSIS:
        _route(...)

    INPUT ARGUMENTS:
    inflow : Input series to be routed, in interval depths (mm).
    time_base : Hydrograph time base (days).
    symmetric : Whether to use the symmetric two-time-base hydrograph.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.

    OUTPUT ARGUMENTS:
    Routed interval depths (mm), with the input series length.
    """
    length = (13 if symmetric else 6)*resolution
    lag = np.arange(length+1) / (resolution*time_base)
    if symmetric:
        s = np.where(lag <= 1, .5*lag**2.5,
                     1-.5*np.maximum(2-lag, 0.)**2.5)
    else:
        s = np.minimum(lag, 1.)**2.5
    return np.convolve(inflow, np.diff(s))[:len(inflow)]


@lru_cache(maxsize=2)
def _kernel(uniform_input=False):
    """
    Load and configure the selected native model kernel.

    SYNOPSIS:
        _kernel(...)

    INPUT ARGUMENTS:
    uniform_input : Use uniform-within-interval inputs instead of step-start
        inputs.

    OUTPUT ARGUMENTS:
    Native model callable with the selected formulation ABI.
    """
    path = Path(__file__).with_name("cart_daily_v2.dll")
    if not path.is_file():
        path = path.with_name("cart_daily.dll")
    library = CDLL(str(path))
    call = library.cart_uniform if uniform_input else library.cart_daily
    pointer = POINTER(c_double)
    call.argtypes = [
        pointer, c_int, pointer, pointer, pointer, pointer, pointer,
        pointer, c_int, c_int,
    ]
    call.restype = c_int
    return call


def _native(forcing, elevation, area, p, jacobian, initial, resolution,
            uniform_input=False):
    """
    Marshal arrays to the native discharge and tangent kernel.

    SYNOPSIS:
        _native(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    p : 11 CaRT calibration coordinates, including logarithmic entries.
    jacobian : Whether to evaluate discharge sensitivities.
    initial : Initial stores (mm), or None/empty to use the spin-up policy.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    uniform_input : Use uniform-within-interval inputs instead of step-start
        inputs.

    OUTPUT ARGUMENTS:
    (q, j, stores) : computed result in the operand units.
    """
    def pointer(x):
        """
        Return a native double pointer for an optional array.

        SYNOPSIS:
            pointer(...)

        INPUT ARGUMENTS:
        x : Numerical operand; scalar, vector, or forward-mode scalar.

        OUTPUT ARGUMENTS:
        Pointer to contiguous doubles, or a null pointer for None.
        """
        return x.ctypes.data_as(POINTER(c_double))

    f, z, a, p = map(np.ascontiguousarray, (forcing, elevation, area, p))
    q = np.empty(len(f))
    j = np.empty((len(f), 11)) if jacobian else None
    stores = np.empty(3) if initial is None else np.array(initial, float)
    status = _kernel(uniform_input)(
        pointer(f), len(f), pointer(z), pointer(a), pointer(p),
        pointer(q), pointer(j) if j is not None else None,
        pointer(stores), int(initial is not None), resolution,
    )
    if status:
        raise ValueError('CaRT could not initialize the ' 'routing stores.')
    return q, j, stores
