# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/test_cartplus_jacobian.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Check all CaRT+ solver Jacobians against five-point differences.
"""

import ctypes
import json
from pathlib import Path

import numpy as np

from cartplus_cpu import _entrypoint, integrate_zone, python_model, simulate
from cartplus_inputs import physical_parameters


def check(resolution, mcode, library_name="cartplus_cpu.dll", balance_atol
    =1e-10):
    """
    Check native/Python model discharge and Jacobian agreement.

    SYNOPSIS:
        check(...)

    INPUT ARGUMENTS:
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    mcode : Backend/solver case selected by the model wrapper.
    library_name : Filename of the native model library.
    balance_atol : Absolute water-balance validation tolerance (mm).

    OUTPUT ARGUMENTS:
    (np.max(abs(analytic - reference)), np.max(abs(numerical - reference)))
    : computed result in the operand units.
    """
    root = Path(__file__).resolve().parent
    spec = json.loads((root / "parameters.json").read_text())

    def defaults(group):
        """
        Read baseline parameter values in the selected group.

        SYNOPSIS:
            defaults(...)

        INPUT ARGUMENTS:
        group : Group used by this numerical operator.

        OUTPUT ARGUMENTS:
        {key: value[2] for key, value in spec[group].items()} : computed
        result in the operand units.
        """
        return {key: value[2] for key, value in spec[group].items()}

    parameters = physical_parameters(
        defaults("static"), defaults("law"), spec["shared"],
    )
    t = np.arange(10.)
    p = (2 + 1.5*np.sin(t) + 8*(np.sin(2*t) > .6)) / resolution
    temperature = 4 - 14*np.cos(t/2)
    pet = (.7 + .5*np.cos(t)) / resolution
    f = np.column_stack((p, temperature, pet, (t >= 1) & (t < 8)))
    z = np.array([100., 220., 350., 510., 740.])
    area = np.array([.1, .2, .3, .25, .15])
    library = root / (library_name if mcode == 4 else "missing.dll")
    options = dict(resolution=resolution, mcode=mcode, jacobian=True)
    out, analytic = simulate(library, f, z, area, parameters, **options)
    _, numerical = simulate(
        library, f, z, area, parameters, numerical=True, **options,
    )
    if mcode < 4:
        model = python_model()

        def zone(initial, net, recession, dt):
            """
            Evaluate near-stream runoff, evaporation, and drainage.

            SYNOPSIS:
                zone(...)

            INPUT ARGUMENTS:
            initial : Initial stores (mm), or None/empty to use the spin-up
                policy.
            net : Whether to remove interception before the production
                update.
            recession : Groundwater recession coefficient in the model time
                base.
            dt : Integration or observation interval duration (days).

            OUTPUT ARGUMENTS:
            Near-stream runoff, evaporation, and drainage rates.
            """
            return integrate_zone(initial, net, recession, dt, False, mcode)

        _, states = model.simulate(
            f, z, area, parameters, resolution, 0., zone,
        )

        def discharge(x):
            """
            Extract interval discharge from a model callback.

            SYNOPSIS:
                discharge(...)

            INPUT ARGUMENTS:
            x : Numerical operand; scalar, vector, or forward-mode scalar.

            OUTPUT ARGUMENTS:
            model.simulate(f, z, area, x, resolution, 0.0, zone,
            states)[0][:, 0] : computed result in the operand units.
            """
            return model.simulate(
                f, z, area, x, resolution, 0., zone, states,
            )[0][:, 0]
    else:
        interval = resolution != 1
        extra = (0., 1/resolution) if interval else ()
        forward = _entrypoint(library, False, interval)
        states = (ctypes.c_double*7)()
        initial = np.empty_like(f)
        assert forward(
            len(f), f, z, area, parameters, initial, None, states, *extra,
        ) == 0

        def discharge(x):
            """
            Extract interval discharge from a model callback.

            SYNOPSIS:
                discharge(...)

            INPUT ARGUMENTS:
            x : Numerical operand; scalar, vector, or forward-mode scalar.

            OUTPUT ARGUMENTS:
            result[:, 0] : computed result in the operand units.
            """
            result = np.empty_like(f)
            assert forward(
                len(f), f, z, area, x, result, states, None, *extra,
            ) == 0
            return result[:, 0]

    reference = np.empty_like(analytic)
    for j, value in enumerate(parameters):
        h = np.zeros_like(parameters)
        h[j] = 3e-6 * max(abs(value), 1e-3)
        plus2 = discharge(parameters + 2*h)
        plus1 = discharge(parameters + h)
        minus1 = discharge(parameters - h)
        minus2 = discharge(parameters - 2*h)
        reference[:, j] = (-plus2+8*plus1-8*minus1+minus2)/(12*h[j])
    for jacobian in (analytic, numerical):
        np.testing.assert_allclose(
            jacobian, reference, rtol=.005, atol=2e-5,
        )
    balance = p + out[:, 2] - out[:, 0] - out[:, 1] - out[:, 3]
    np.testing.assert_allclose(balance, 0, atol=balance_atol)
    return np.max(abs(analytic-reference)), np.max(abs(numerical-reference))


if __name__ == "__main__":
    for resolution in (1, 24, 96):
        for mcode in (1, 2, 3, 4):
            analytic, numerical = check(resolution, mcode)
            print(f"{resolution:2} steps/day, case {mcode}: "
                  f"analytic {analytic:.3g}, numerical {numerical:.3g}",
                  flush=True)
