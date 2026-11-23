# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_models.py
#
# Public rainfall-runoff models with MATLAB-compatible case selection.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Public rainfall-runoff models with MATLAB-compatible case selection.

mcode=1: Python adaptive RK2 (Heun/Euler); 2: Python Dormand-Prince 5(4);
3: Python explicit Euler with each model's MATLAB substep count;
4: the existing native C++ implementation. Default remains 4.
"""

import importlib
import numpy as np
from numba import njit
from crr_native import parameters, STATE_COUNTS, PARAM_COUNTS
from crr_integrators import fixed_or_adaptive, dormand_prince

_NAMES = {
    1: "hymod",
    2: "hmodel",
    3: "sacsma",
    4: "xinanjiang",
    5: "gr4ja",
    6: "hbv",
    7: "cfe_nwm",
}
_STEPS = {1: 50, 2: 500, 3: 1000, 4: 500, 5: 1000, 6: 500, 7: 500}

# Import at module scope so Numba compiles the existing independent HYMOD
# RHS.
from site_hymod import _rhs as _hymod_rhs


@njit(cache=True)
def _hymod_augmented(z, forcing, values):
    """Evaluate the HYMOD state and sensitivity derivatives."""
    state = z.reshape((8, 7)).T
    derivative = _hymod_rhs(
        state, values[:7], forcing[0], forcing[1], forcing[2], values[7:]
    )
    return derivative.T.copy().reshape(56)


def simulate(
    model,
    theta,
    forcing,
    *,
    mcode=4,
    y0=None,
    model_config=None,
    settings=(0.01, 1.0, 1e-4, 1e-3, 1e-3, 2.0, 10000.0),
    routing_start=0,
    route_hbv=True
):
    """Simulate the selected model and requested sensitivity outputs."""
    model = int(model)
    if model == 8:
        if mcode != 4:
            raise ValueError(
                "The MATLAB-compatible authored user_model supports only"
                " mcode=4"
            )
        from crr_native import simulate as native

        return native(
            model,
            theta,
            forcing,
            y0=y0,
            model_config=model_config,
            settings=settings,
            routing_start=routing_start,
        )
    if model not in _NAMES or mcode not in (1, 2, 3, 4):
        raise ValueError("Choose public model 1..7 and mcode 1..4")
    if mcode == 4 and model == 6 and not route_hbv:
        return _native_hbv_core(theta, forcing, y0, model_config, settings)
    if mcode == 4:
        from crr_native import simulate as native

        return native(
            model,
            theta,
            forcing,
            y0=y0,
            model_config=model_config,
            settings=settings,
            routing_start=routing_start,
        )
    theta = np.asarray(theta, dtype=float)
    forcing = np.ascontiguousarray(forcing, dtype=float)
    settings = np.asarray(settings, dtype=float)
    if theta.shape != (
        PARAM_COUNTS[model] - (model == 6 and not route_hbv),
    ):
        raise ValueError("Incorrect parameter count")
    if forcing.ndim != 2 or forcing.shape[1] != 3 or not len(forcing):
        raise ValueError("Need nonempty N-by-3 forcing")
    if (
        settings.shape != (7,)
        or min(settings[:5]) <= 0
        or settings[2] > settings[1]
        or settings[5] <= 0
        or settings[6] < 1
    ):
        raise ValueError("Invalid solver settings")
    p = parameters(model, theta, model_config)
    p["giuh"] = np.asarray(
        (model_config or {}).get(
            "giuh_ordnts", [0.1, 0.35, 0.2, 0.14, 0.1, 0.06, 0.05]
        ),
        dtype=float,
    )
    p["L"] = len(p["giuh"])
    if model == 7 and (
        int(p["K"]) != p["K"]
        or p["K"] < 1
        or not np.isfinite(p["giuh"]).all()
        or np.any(p["giuh"] < 0)
        or not np.isclose(p["giuh"].sum(), 1)
    ):
        raise ValueError("Invalid CFE routing configuration")
    m = STATE_COUNTS[model] if model != 7 else 4 + int(p["K"]) + p["L"]
    d = 12 if model == 6 else len(theta)
    y0 = np.full(m, 1e-5) if y0 is None else np.asarray(y0, dtype=float)
    if y0.shape != (m,) or not all(
        np.all(np.isfinite(v)) for v in (theta, forcing, settings, y0)
    ):
        raise ValueError("Invalid model inputs")
    z = np.zeros(m * (d + 1))
    z[:m] = y0
    if model == 1:
        rhs = _hymod_augmented
        values = np.r_[theta, [p["T_sm"], p["eps_m"], p["rho"]]]
    else:
        equations = importlib.import_module(
            "crr_" + _NAMES[model] + "_equations"
        )
        rhs = equations.augmented
        values = tuple(p[k] for k in equations.FIELDS)
    if mcode == 2:
        history = dormand_prince(
            rhs,
            forcing,
            z,
            values,
            settings,
            float(len(forcing)),
            model == 4,
        )
    else:
        history = fixed_or_adaptive(
            rhs,
            forcing,
            z,
            values,
            settings,
            _STEPS[model] if mcode == 3 else 0,
        )
    discharge = np.diff(history[:, m - 1])
    jacobian = np.diff(history[:, m - 1 :: m], axis=0)[:, 1:]
    if model == 6 and route_hbv:
        from hbv_routing import route

        discharge, jacobian = route(
            discharge,
            jacobian,
            theta[12],
            (model_config or {}).get("routing_length", 60),
            routing_start,
        )
    return discharge, jacobian


def _native_hbv_core(theta, forcing, y0, config, settings):
    """
    Existing native entry point, without Python's external MAXBAS
    filter.
    """
    import ctypes as C
    from crr_native import library

    theta = np.asarray(theta, float)
    forcing = np.asfortranarray(forcing, dtype=float)
    if (
        theta.shape != (12,)
        or forcing.ndim != 2
        or forcing.shape[1] != 3
        or not len(forcing)
    ):
        raise ValueError("Invalid raw HBV inputs")
    lib, manifest = library()
    p = parameters(6, theta, config)
    packed = np.asarray(
        [p[k] for k in manifest["6"]["fields"]], dtype=float
    )
    initial = (
        np.full(5, 1e-5)
        if y0 is None
        else np.ascontiguousarray(y0, dtype=float)
    )
    settings = np.ascontiguousarray(settings, dtype=float)
    giuh = np.zeros(1)
    if (
        initial.shape != (5,)
        or settings.shape != (7,)
        or not all(
            np.isfinite(v).all()
            for v in (forcing, packed, initial, settings)
        )
    ):
        raise ValueError("Invalid raw HBV inputs")
    if (
        min(settings[:5]) <= 0
        or settings[2] > settings[1]
        or settings[5] <= 0
        or not 1 <= settings[6] <= np.iinfo(np.int32).max
    ):
        raise ValueError("Invalid solver settings")
    q = np.empty(len(forcing))
    J = np.empty((len(forcing), 12), order="F")
    ptr = C.POINTER(C.c_double)

    def addr(x):
        """Return a native pointer to the supplied numeric array."""
        return x.ctypes.data_as(ptr)

    status = lib.site_crr_run(
        6,
        len(forcing),
        5,
        12,
        *[addr(forcing[:, i]) for i in range(3)],
        addr(packed),
        len(packed),
        addr(settings),
        addr(initial),
        addr(giuh),
        1,
        addr(q),
        addr(J)
    )
    if status:
        raise ValueError("Native HBV integration failed")
    return q, J
