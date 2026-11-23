# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/cartplus_cpu.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
CaRT+ forward simulation and discharge sensitivities.
"""
import ctypes
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from math import ceil
from pathlib import Path

import numpy as np


@lru_cache(maxsize=2)
def _entrypoint(library, regional=False, subdaily=False):
    """
    Bind the selected native entry point and its C ABI.

    SYNOPSIS:
        _entrypoint(...)

    INPUT ARGUMENTS:
    library : Native library path or loaded backend handle.
    regional : Use regional elevation-band inputs.
    subdaily : Whether the forcing uses more than one interval per day.

    OUTPUT ARGUMENTS:
    Native callable with argument and return types configured.
    """
    dll = ctypes.CDLL(str(library))
    forward = (dll.cartplus_interval_forward if subdaily else
               dll.cartplus_regional_forward if regional else
               dll.cartplus_us_forward)
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    pointer = ctypes.POINTER(ctypes.c_double)
    forward.argtypes = [ctypes.c_size_t, *([array] * 5), pointer, pointer]
    if regional or subdaily:
        forward.argtypes += [ctypes.c_double]
    if subdaily:
        forward.argtypes += [ctypes.c_double]
    forward.restype = ctypes.c_int
    return forward


@lru_cache(maxsize=2)
def _jacobian_entrypoint(
    library, regional=False, subdaily=False, calibrated_only=False,
):
    """
    Bind the native forward-sensitivity entry point.

    SYNOPSIS:
        _jacobian_entrypoint(...)

    INPUT ARGUMENTS:
    library : Native library path or loaded backend handle.
    regional : Use regional elevation-band inputs.
    subdaily : Whether the forcing uses more than one interval per day.
    calibrated_only : Whether to retain only calibrated parameter
        directions.

    OUTPUT ARGUMENTS:
    Native forward-sensitivity callable with its C ABI configured.
    """
    dll = ctypes.CDLL(str(library))
    prefix = "interval" if subdaily else "regional" if regional else "us"
    suffix = "35" if calibrated_only else ""
    derivative = getattr(dll, f"cartplus_{prefix}_jacobian{suffix}")
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    pointer = ctypes.POINTER(ctypes.c_double)
    derivative.argtypes = [ctypes.c_size_t, *([array] * 5), pointer]
    if regional or subdaily:
        derivative.argtypes += [ctypes.c_double]
    if subdaily:
        derivative.argtypes += [ctypes.c_double]
    derivative.restype = ctypes.c_int
    return derivative


def simulate(
    library, forcing, elevation, area, parameters, *, jacobian=False,
    numerical=False, infiltration_floor=0.0, resolution=1, mcode=4,
    calibrated_only=False, uniform_input=False,
):
    """
    Simulate model fluxes and optional discharge sensitivities.

    SYNOPSIS:
        simulate(...)

    INPUT ARGUMENTS:
    library : Native library path or loaded backend handle.
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    parameters : 42 physical CaRT+ parameters in native order.
    jacobian : Whether to evaluate discharge sensitivities.
    numerical : Use numerical rather than analytic discharge sensitivities.
    infiltration_floor : Minimum infiltration capacity (mm/day).
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    mcode : Backend/solver case selected by the model wrapper.
    calibrated_only : Whether to retain only calibrated parameter
        directions.
    uniform_input : Use uniform-within-interval inputs instead of step-start
        inputs.

    OUTPUT ARGUMENTS:
    output : N-by-4 interval water-balance fluxes (mm); J : optional
    discharge Jacobian in the selected coordinates.

    NOTES:
    Cases 1-3 integrate the near-stream and wetness ODEs in Python. Cases
    1-3 use the standalone Python model. Case 4 loads C++. Jacobians are
    discharge derivatives at fixed spin-up stores. A calibrated-only DLL
    returns its 35 seeded physical columns.
    """
    if mcode not in (1, 2, 3, 4):
        raise ValueError("CaRT+ mcode must be 1, 2, 3 or 4.")
    if calibrated_only and (mcode != 4 or numerical):
        raise ValueError('The 35-column Jacobian needs n'
            'ative analytics.')
    if resolution not in (1, 24, 96):
        raise ValueError('Use daily, hourly or 15-minute' ' CaRT+ forcing.')
    forcing, elevation, area, parameters = (
        np.ascontiguousarray(value, dtype=np.float64)
        for value in (forcing, elevation, area, parameters)
    )
    if forcing.ndim != 2 or forcing.shape[1] != 4 or not len(forcing):
        raise ValueError(
            'Expected interval precipitatio' 'n, temperature, PET and mask.'
        )
    if (elevation.shape != (5,) or area.shape != (5,)
            or parameters.shape != (42,)):
        raise ValueError('Expected five bands and 42 phy'
            'sical parameters.')
    if not all(np.isfinite(v).all()
               for v in (forcing, elevation, area, parameters)):
        raise ValueError("CaRT+ inputs must be finite.")
    if mcode != 4:
        return standalone(
            forcing, elevation, area, parameters, resolution,
            infiltration_floor, mcode, jacobian, numerical,
            uniform_input=uniform_input,
        )
    subdaily = resolution != 1
    extra = (infiltration_floor,) if infiltration_floor or subdaily else ()
    if subdaily:
        extra += (1.0 / resolution,)
    forward = _entrypoint(library, bool(infiltration_floor), subdaily)
    spinup = (ctypes.c_double * 7)()

    def run(p, fixed=None, captured=None):
        """
        Advance the pilot or full model over the forcing sequence.

        SYNOPSIS:
            run(...)

        INPUT ARGUMENTS:
        p : Model parameter vector; CaRT calibration or CaRT+ physical
            units.
        fixed : Detached spin-up stores, or None/empty to determine them.
        captured : Optional collector receiving baseline spin-up stores.

        OUTPUT ARGUMENTS:
        output : computed result in the operand units.
        """
        output = np.empty_like(forcing)
        status = forward(len(forcing), forcing, elevation, area, p, output,
                         fixed, captured, *extra)
        if status or not np.isfinite(output).all():
            raise ValueError("CaRT+ native simulation failed.")
        return output

    baseline = run(parameters, captured=spinup)
    if not jacobian:
        return baseline, None
    result = np.empty((len(forcing), 35 if calibrated_only else 42))
    if not numerical:
        derivative = _jacobian_entrypoint(
            library, bool(infiltration_floor), subdaily, calibrated_only,
        )
        status = derivative(
            len(forcing), forcing, elevation, area, parameters,
            result, spinup, *extra,
        )
        if status or not np.isfinite(result).all():
            raise ValueError("CaRT+ native discharge Jacobian failed.")
        return baseline, result
    for j, value in enumerate(parameters):
        step = 1e-5 * max(abs(value), 1e-3)
        plus, minus = parameters.copy(), parameters.copy()
        plus[j] += step
        minus[j] -= step
        result[:, j] = (
            run(plus, spinup)[:, 0] - run(minus, spinup)[:, 0]
        ) / (2 * step)
    return baseline, result


@lru_cache(maxsize=1)
def python_model(uniform_input=False):
    """
    Load the selected sequential CaRT+ Python implementation.

    SYNOPSIS:
        python_model(...)

    INPUT ARGUMENTS:
    uniform_input : Use uniform-within-interval inputs instead of step-start
        inputs.

    OUTPUT ARGUMENTS:
    Loaded Python module implementing the selected formulation.
    """
    path = Path(__file__).with_name(
        "cartplus_uniform_python.py" if uniform_input else
            "cartplus_python.py"
    )
    spec = spec_from_file_location("_cartplus_python", path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def standalone(forcing, elevation, area, p, resolution, floor, mcode,
               jacobian, numerical, fixed=None, uniform_input=False):
    """
    Run sequential CaRT+ and optional fixed-spin-up derivatives.

    SYNOPSIS:
        standalone(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    p : 42 physical CaRT+ parameters in native order.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    floor : Minimum infiltration capacity in daily-rate units (mm/day).
    mcode : Backend/solver case selected by the model wrapper.
    jacobian : Whether to evaluate discharge sensitivities.
    numerical : Use numerical rather than analytic discharge sensitivities.
    fixed : Detached spin-up stores, or None/empty to determine them.
    uniform_input : Use uniform-within-interval inputs instead of step-start
        inputs.

    OUTPUT ARGUMENTS:
    Computed fluxes, stores, or numerical values described in the notes.
    """
    model = python_model(uniform_input)

    def zone(initial, net, recession, dt):
        """
        Evaluate near-stream runoff, evaporation, and drainage.

        SYNOPSIS:
            zone(...)

        INPUT ARGUMENTS:
        initial : Initial stores (mm), or None/empty to use the spin-up
            policy.
        net : Whether to remove interception before the production update.
        recession : Groundwater recession coefficient in the model time
            base.
        dt : Integration or observation interval duration (days).

        OUTPUT ARGUMENTS:
        Near-stream runoff, evaporation, and drainage rates.
        """
        inputs = (initial, net, recession)
        seeded = any(isinstance(x, model.Tangent) for x in inputs)
        result = integrate_zone(
            *(model.primal(x) for x in inputs), dt, seeded, mcode,
        )
        if not seeded:
            return result
        return [
            model.Tangent(result[i], sum(
                result[2+2*j+i]*model.Tangent.lift(x).slope
                for j, x in enumerate(inputs)
            ))
            for i in range(2)
        ]

    def run(values, states):
        """
        Advance the pilot or full model over the forcing sequence.

        SYNOPSIS:
            run(...)

        INPUT ARGUMENTS:
        values : Numerical values or vector of state entries.
        states : Hydrological state history or detached spin-up stores.

        OUTPUT ARGUMENTS:
        model.simulate(forcing, elevation, area, values, resolution, floor,
        zone, states) : computed result in the operand units.
        """
        return model.simulate(
            forcing, elevation, area, values, resolution, floor,
            zone, states,
        )

    baseline, spinup = run(p, fixed)
    if not np.isfinite(baseline).all():
        raise ValueError("CaRT+ standalone simulation failed.")
    if not jacobian:
        return baseline, None
    if not numerical:
        seeded = [model.Tangent(x, s) for x, s in zip(p, np.eye(len(p)))]
        result, _ = run(seeded, spinup)
        slopes = np.vstack([
            np.broadcast_to(model.Tangent.lift(x).slope, p.shape)
            for x in result[:, 0]
        ])
        if not np.isfinite(slopes).all():
            raise ValueError('CaRT+ standalone discharge Jac'
                'obian failed.')
        return baseline, slopes
    slopes = np.empty((len(forcing), len(p)))
    for j, value in enumerate(p):
        step = 1e-5 * max(abs(value), 1e-3)
        plus, minus = p.copy(), p.copy()
        plus[j] += step
        minus[j] -= step
        slopes[:, j] = (
            run(plus, spinup)[0][:, 0] - run(minus, spinup)[0][:, 0]
        ) / (2*step)
    return baseline, slopes



def integrate_zone(initial, net, recession, dt, derivative, mcode):
    """
    Integrate a near-stream store with the selected solver.

    SYNOPSIS:
        integrate_zone(...)

    INPUT ARGUMENTS:
    initial : Initial stores (mm), or None/empty to use the spin-up policy.
    net : Whether to remove interception before the production update.
    recession : Groundwater recession coefficient in the model time base.
    dt : Integration or observation interval duration (days).
    derivative : Requested partial derivative or derivative-output flag.
    mcode : Backend/solver case selected by the model wrapper.

    OUTPUT ARGUMENTS:
    Flattened endpoint stores and local forward partial derivatives.

    NOTES:
    Time is normalized to [0, 1]. The optional three sensitivity columns
    correspond to initial storage, net input, and recession.
    """
    y = np.zeros((2, 4 if derivative else 1))
    y[0, 0] = initial
    if derivative:
        y[0, 1] = 1

    def rhs(time, state):
        """
        Evaluate augmented water-state and sensitivity derivatives.

        SYNOPSIS:
            rhs(...)

        INPUT ARGUMENTS:
        time : Current model time (days).
        state : State vector (including tangent lanes when enabled).

        OUTPUT ARGUMENTS:
        rate.ravel() : computed result in the operand units.
        """
        state = state.reshape(y.shape)
        s = state[0, 0]
        if net >= 0:
            shape, slope = 1 - s*s, -2*s
        else:
            shape, slope = s*(2-s), 2-2*s
        rate = np.empty_like(state)
        rate[0] = (net*slope-recession)*state[0]
        rate[0, 0] = net*shape-recession*s
        rate[1] = state[0]
        if derivative:
            rate[0, 2] += shape
            rate[0, 3] -= s
        return rate.ravel()

    state = y.ravel()
    if mcode == 2:
        from scipy.integrate import solve_ivp

        sol = solve_ivp(rhs, (0, 1), state, method="RK45",
                        rtol=1e-10, atol=1e-12)
        if not sol.success:
            raise ValueError(sol.message)
        state = sol.y[:, -1]
    else:
        steps = max(1, ceil(16*dt), ceil(4*(abs(net)+recession)))
        step = 1 / steps
        for _ in range(steps):
            rate = rhs(0, state)
            if mcode == 1:
                rate = rhs(0, state + .5*step*rate)
            state = state + step*rate
    return state.reshape(y.shape).ravel(order="F")

