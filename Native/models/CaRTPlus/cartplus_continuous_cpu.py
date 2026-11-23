# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/cartplus_continuous_cpu.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Python binding for the native, fully coupled CaRT+ integrator.
"""

import ctypes
from functools import lru_cache
from pathlib import Path

import numpy as np

from cartplus_continuous import _inputs
from cartplus_inputs import PARAMETER_MAP

CALIBRATED_COLUMNS = tuple(k for k, (group, _, _) in enumerate(
    PARAMETER_MAP)
                           if group != 2)


_PROFILED_SOLVERS = ('rk2_replay', 'rk2_optimized', 'rk2_analytic_quick')
_EXPERIMENTAL_SOLVERS = ('rk2_workspace', 'rk2_exact_vegetation',
    'rk2_analytic_routing',
                         'rk2_replay', 'rk2_optimized',
                             'rk2_analytic_quick')
_SOLVERS = ('dop853', 'rk2', 'rk23', 'rk45', 'heun_fixed', 'euler', *
    _EXPERIMENTAL_SOLVERS)


@lru_cache(maxsize=None)
def _entrypoint(calibrated_only, solver, configured=False, keep_tangents=False):
    """
    Bind the selected native entry point and its C ABI.

    SYNOPSIS:
        _entrypoint(...)

    INPUT ARGUMENTS:
    calibrated_only : Whether to retain only calibrated parameter
        directions.
    solver : Integration method or solver-control object.
    configured : Whether to bind the entry point with explicit step limits.

    OUTPUT ARGUMENTS:
    Native callable with argument and return types configured.
    """
    suffix = '_35' if calibrated_only else ''
    if solver in _EXPERIMENTAL_SOLVERS:
        suffix += '_experimental'
    path = Path(__file__).with_name(f'cartplus_continuous_cpu{suffix}.dll')
    if solver not in _EXPERIMENTAL_SOLVERS:
        optimized = path.with_name(
            f'cartplus_continuous_cpu{suffix}_v3.dll')
        if optimized.is_file():
            path = optimized
    if configured or (solver not in _EXPERIMENTAL_SOLVERS
                      and not path.is_file()):
        path = path.with_name('cartplus_numerics_v1' + ('_35' if
            calibrated_only else '') + '.dll')
    library = ctypes.CDLL(str(path))
    entry = 'cartplus_continuous_run' if solver == 'dop853' else \
        f'cartplus_continuous_{solver}_run'
    if keep_tangents:
        entry = 'cartplus_configured_predict'
    elif configured:
        entry = 'cartplus_configured_run'
    call = getattr(library, entry)
    pointer = ctypes.POINTER(ctypes.c_double)
    call.argtypes = [ctypes.c_size_t, *([pointer]*6), ctypes.c_int, *([
        pointer]*5)]
    if solver in _PROFILED_SOLVERS:
        call.argtypes += [pointer]
    if configured:
        call.argtypes += [ctypes.c_int, ctypes.c_double, ctypes.c_double]
    call.restype = ctypes.c_int
    return call


def simulate(forcing, elevation, area, parameters, *, jacobian=False,
             initial=None, resolution=1, interval_days=None, refinement=16,
             rtol=1e-8, atol=1e-10, infiltration_floor=0.,
             availability_days=.02, transition_fraction=.01,
             sensitivity_coordinates='physical', return_diagnostics=False,
             calibrated_only=False, solver='dop853', max_step=None,
                 min_step=0., keep_tangents=False):
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
    jacobian : Whether to evaluate discharge sensitivities.
    keep_tangents : Use the augmented solver without returning a Jacobian.
    initial : Initial stores (mm), or None/empty to use the spin-up policy.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    interval_days : Positive interval duration (days), scalar or length N.
    refinement : Positive number of routing stages per X4 time base.
    rtol : Positive relative integration-error tolerance (-).
    atol : Positive absolute integration-error tolerance for augmented
        states.
    infiltration_floor : Minimum infiltration capacity (mm/day).
    availability_days : Water-availability regularization time (days).
    transition_fraction : Dimensionless transition smoothing width.
    sensitivity_coordinates : Sensitivity directions: physical or
        transformed coordinates.
    return_diagnostics : Whether to append solver/water-balance diagnostics.
    calibrated_only : Whether to retain only calibrated parameter
        directions.
    solver : Integration method or solver-control object.
    max_step : Maximum/fixed solver step (days), or None for backend
        default.
    min_step : Nonnegative lower solver-step bound (days).

    OUTPUT ARGUMENTS:
    q : N discharge depths (mm/interval); J : discharge Jacobian or None.
    Diagnostics are appended only when requested.
    """
    if solver not in _SOLVERS:
        raise ValueError(f"solver must be one of {_SOLVERS}")
    keep_tangents = keep_tangents and not jacobian
    if keep_tangents:
        try:
            _entrypoint(calibrated_only, solver, True, True)
        except (AttributeError, OSError):
            # Older installed libraries only expose tangent mode with J.
            q, *outputs = simulate(
                forcing, elevation, area, parameters, jacobian=True,
                initial=initial, resolution=resolution, interval_days=interval_days,
                refinement=refinement, rtol=rtol, atol=atol,
                infiltration_floor=infiltration_floor,
                availability_days=availability_days,
                transition_fraction=transition_fraction,
                sensitivity_coordinates=sensitivity_coordinates,
                return_diagnostics=return_diagnostics,
                calibrated_only=calibrated_only, solver=solver,
                max_step=max_step, min_step=min_step,
            )
            return (q, None, *outputs[1:])
    configured = keep_tangents or max_step is not None or min_step != 0 or solver in (
        'rk45','heun_fixed','euler')
    if configured:
        if solver in _EXPERIMENTAL_SOLVERS:
            raise ValueError('Step limits are not supported '
                'for experimental solvers.')
        max_step = (1/96 if solver in ('rk45','heun_fixed','euler') else
            1e100) if max_step is None else float(max_step)
        if not np.isfinite([max_step,min_step]).all() or not 0 <= min_step \
            <= max_step or max_step <= 0:
            raise ValueError('Require 0 <= min_step <= max_s'
                'tep, in days.')
    f, z, a, p, dt, initial, seed = _inputs(
        forcing, elevation, area, parameters, initial, resolution,
            interval_days,
        refinement, rtol, atol, infiltration_floor, availability_days,
        transition_fraction, sensitivity_coordinates)
    pointer = ctypes.POINTER(ctypes.c_double)

    def address(array):
        """
        Return a native double pointer for an optional array.

        SYNOPSIS:
            address(...)

        INPUT ARGUMENTS:
        array : Contiguous numerical array, or None for a null native
            pointer.

        OUTPUT ARGUMENTS:
        Pointer to contiguous doubles, or a null pointer for None.
        """
        return None if array is None else array.ctypes.data_as(pointer)

    f, z, a, p, dt, seed = map(np.ascontiguousarray, (f, z, a, p, dt, seed))
    if initial is not None:
        initial = np.ascontiguousarray(initial)
    opt = np.array([rtol, atol, infiltration_floor, availability_days,
        transition_fraction])
    q = np.empty(len(f))
    width = len(CALIBRATED_COLUMNS) if calibrated_only else 42
    j = np.empty((len(f), width)) if jacobian else None
    balance = np.empty(len(f)) if return_diagnostics else None
    stats = np.empty((2, 2)) if return_diagnostics and solver in \
        _PROFILED_SOLVERS else None
    extra = (address(stats),) if solver in _PROFILED_SOLVERS else ()
    if configured:
        codes = dict(dop853=0, rk2=1, rk23=2, rk45=9, heun_fixed=10, euler
            =11)
        extra = (codes[solver],max_step,min_step)
    status = _entrypoint(calibrated_only, solver, configured, keep_tangents)(len(f), *map
        (address, (f, z, a, p, dt, initial)),
                           refinement, *map(address, (opt, seed, q, j,
                               balance)), *extra)
    if status or not all(np.isfinite(v).all() for v in (q, j, balance) if
        v is not None):
        raise RuntimeError('Native coupled CaRT+ integration failed.')
    if return_diagnostics:
        diagnostics = {'balance_residual_mm': balance}
        if stats is not None:
            diagnostics['rk2_trials'] = stats.astype(np.int64)
        return q, j, diagnostics
    return q, j
