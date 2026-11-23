# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/cart/cart_ode_cpu.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Python binding for the native event-aware coupled CaRT integrator.
"""

import ctypes
from functools import lru_cache
from pathlib import Path

import numpy as np

from cart_ode import _inputs

_POINTER = ctypes.POINTER(ctypes.c_double)
_SOLVERS = {'dop853': 0, 'rk2': 1, 'rk23': 2, 'rk45': 9, 'heun_fixed': 10,
    'euler': 11}


@lru_cache(maxsize=2)
def _entrypoint(configured=False):
    """
    Bind the selected native entry point and its C ABI.

    SYNOPSIS:
        _entrypoint(...)

    INPUT ARGUMENTS:
    configured : Whether to bind the entry point with explicit step limits.

    OUTPUT ARGUMENTS:
    Native callable with argument and return types configured.
    """
    path = Path(__file__).with_name('cart_continuous_cpu_v2.dll')
    if not path.is_file():
        path = path.with_name('cart_continuous_cpu.dll')
    if configured or not path.is_file():
        path = path.with_name('cart_numerics_v1.dll')
    library = ctypes.CDLL(str(path))
    call = library.cart_configured_run if configured else library. \
        cart_continuous_run
    call.argtypes = [ctypes.c_size_t, *([_POINTER]*6), ctypes.c_int,
                     ctypes.c_double, ctypes.c_double, ctypes.c_int,
                     *([_POINTER]*3)]
    if configured:
        call.argtypes += [ctypes.c_double, ctypes.c_double]
    call.restype = ctypes.c_int
    return call


def simulate(forcing, elevation, area, parameters, *, jacobian=False,
             initial=None, resolution=1, refinement=16, rtol=1e-8,
             atol=1e-10, interval_days=None, return_diagnostics=False,
             solver='dop853', max_step=None, min_step=0.):
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
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    refinement : Positive number of routing stages per X4 time base.
    rtol : Positive relative integration-error tolerance (-).
    atol : Positive absolute integration-error tolerance for augmented
        states.
    interval_days : Positive interval duration (days), scalar or length N.
    return_diagnostics : Whether to append solver/water-balance diagnostics.
    solver : Integration method or solver-control object.
    max_step : Maximum/fixed solver step (days), or None for backend
        default.
    min_step : Nonnegative lower solver-step bound (days).

    OUTPUT ARGUMENTS:
    q : N interval discharge depths (mm); J : N-by-11 Jacobian or None.
    Optional stores/diagnostics are appended only when requested.
    """
    if solver not in _SOLVERS:
        raise ValueError("solver must be 'dop853', 'rk2', or 'rk23'")
    configured = max_step is not None or min_step != 0 or solver in (
        'rk45','heun_fixed','euler')
    if configured:
        max_step = (1/96 if solver in ('rk45','heun_fixed','euler') else
            1e100) if max_step is None else float(max_step)
        if not np.isfinite([max_step,min_step]).all() or not 0 <= min_step \
            <= max_step or max_step <= 0:
            raise ValueError('Require 0 <= min_step <= max_s'
                'tep, in days.')
    f, z, a, p, dt, initial = _inputs(
        forcing, elevation, area, parameters, initial, resolution,
        interval_days, refinement, rtol, atol)

    def pointer(array):
        """
        Return a native double pointer for an optional array.

        SYNOPSIS:
            pointer(...)

        INPUT ARGUMENTS:
        array : Contiguous numerical array, or None for a null native
            pointer.

        OUTPUT ARGUMENTS:
        Pointer to contiguous doubles, or a null pointer for None.
        """
        return None if array is None else array.ctypes.data_as(_POINTER)

    inputs = [np.ascontiguousarray(x) if x is not None else None
              for x in (f, z, a, p, dt, initial)]
    q = np.empty(len(f))
    j = np.empty((len(f), 11)) if jacobian else None
    balance = np.empty(len(f)) if return_diagnostics else None
    extra = (max_step,min_step) if configured else ()
    status = _entrypoint(configured)(len(f), *map(pointer, inputs),
        refinement, rtol, atol,
                           _SOLVERS[solver], *map(pointer, (q, j, balance)
                               ), *extra)
    if status or not all(np.isfinite(v).all() for v in (q, j, balance) if
        v is not None):
        raise RuntimeError('Native coupled CaRT ' f'{solver}'
            ' integration failed.')
    if return_diagnostics:
        return q, j, {'balance_residual_mm': balance}
    return q, j
