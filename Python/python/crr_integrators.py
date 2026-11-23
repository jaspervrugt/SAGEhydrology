# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_integrators.py
#
# Integration paths corresponding to MATLAB model cases 1, 2 and 3.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Integration paths corresponding to MATLAB model cases 1, 2 and 3.

All methods evolve the full augmented state, including analytic
sensitivities.
The native C++ path (case 4) remains in crr_native.
"""

import numpy as np
from numba import njit

# Continuous-extension coefficients used by MATLAB's Dormand-Prince path.
# SciPy uses a different quartic extension, so its RK45.P is not
# interchangeable.
_DENSE = np.array(
    [
        [1.0, -183 / 64, 37 / 12, -145 / 128],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 1500 / 371, -1000 / 159, 1000 / 371],
        [0.0, -125 / 32, 125 / 12, -375 / 64],
        [0.0, 9477 / 3392, -729 / 106, 25515 / 6784],
        [0.0, -11 / 7, 11 / 3, -55 / 28],
        [0.0, 3 / 2, -4.0, 5 / 2],
    ]
)


@njit(cache=True, error_model="numpy")
def fixed_or_adaptive(rhs, forcing, z, values, settings, steps):
    """Integrate the augmented state with Euler or adaptive RK2."""
    hin, hmax, hmin, rtol, atol, order, maxreject = settings
    history = np.empty((len(forcing) + 1, len(z)))
    history[0] = z
    carry = max(hmin, min(hin, hmax))
    for day in range(len(forcing)):
        f = forcing[day]
        if steps:
            for it in range(steps):
                z = z + rhs(z, f, values) / steps
        else:
            # Each meteorological row is a forcing discontinuity. Use a
            # local
            # clock and restart the controller; carry every augmented state.
            t = 0.0
            end = 1.0
            h = min(max(hmin, min(hin, hmax)), 1.0)
            rejected = 0
            while t < end:
                e = rhs(z, f, values)
                ze = z + h * e
                zh = z + 0.5 * h * (e + rhs(ze, f, values))
                if not np.all(np.isfinite(zh)) or np.max(np.abs(zh)) > 1e12:
                    raise ValueError("Model integration failed")
                error = np.sqrt(
                    np.sum(
                        (
                            (1.0 / (rtol * np.abs(zh) + atol))
                            * np.abs(ze - zh)
                        )
                        ** 2
                    )
                    / len(z)
                )
                if error <= 1 or h <= hmin:
                    z = zh
                    t += h
                    rejected = 0
                else:
                    rejected += 1
                factor = (
                    5.0
                    if error == 0
                    else max(0.2, min(5.0, 0.9 * error ** (-1 / order)))
                )
                carry = max(hmin, min(h * factor, hmax))
                h = min(carry, end - t)
                if rejected >= maxreject:
                    raise ValueError("Model rejected-step limit reached")
        if not np.all(np.isfinite(z)) or np.max(np.abs(z)) > 1e12:
            raise ValueError("Model integration failed")
        history[day + 1] = z
    return history


@njit(cache=True, error_model="numpy")
def _dp_rhs(rhs, z, forcing, values, nonnegative_soil):
    """Evaluate the Dormand-Prince right-hand side with state guards."""
    derivative = rhs(z, forcing, values)
    if nonnegative_soil:
        for index in range(1, 3):
            if z[index] <= 0.0:
                derivative[index] = max(derivative[index], 0.0)
    return derivative


@njit(cache=True, error_model="numpy")
def _dormand_prince_interval(
    rhs, forcing, z, values, settings, final, nonnegative_soil=False
):
    """
    DP5(4) with infinity-norm error control and quartic dense
    output.

    Use MATLAB-compatible component scaling, step rejection and
    growth limits;
    SciPy RK45's RMS controller would give different accepted steps.

    Xinanjiang enables nonnegative_soil for raw W and Sw (indices 1, 2),
    corresponding to ode45 NonNegative=[2,3]. Sensitivities stay signed.
    """
    hin, hmax, unused, rtol, atol, unused_order, unused_maxreject = settings
    if nonnegative_soil and (z[1] < 0.0 or z[2] < 0.0):
        raise ValueError("Nonnegative soil states need nonnegative y0")
    history = np.empty((len(forcing) + 1, len(z)))
    history[0] = z
    t = 0.0
    absh = min(hmax, hin)
    next_output = 1
    f1 = _dp_rhs(rhs, z, forcing[0], values, nonnegative_soil)
    while t < final:
        hmin = 16 * np.spacing(t)
        absh = max(hmin, min(absh, hmax))
        h = absh
        done = 1.1 * absh >= final - t
        if done:
            h = final - t
            absh = abs(h)
        first_try = True
        while True:
            f2 = _dp_rhs(
                rhs,
                z + h * ((1 / 5) * f1),
                forcing[
                    min(int(np.floor(t + h * (1 / 5))), len(forcing) - 1)
                ],
                values,
                nonnegative_soil,
            )
            f3 = _dp_rhs(
                rhs,
                z + h * (3 / 40 * f1 + 9 / 40 * f2),
                forcing[
                    min(int(np.floor(t + h * (3 / 10))), len(forcing) - 1)
                ],
                values,
                nonnegative_soil,
            )
            f4 = _dp_rhs(
                rhs,
                z + h * (44 / 45 * f1 - 56 / 15 * f2 + 32 / 9 * f3),
                forcing[
                    min(int(np.floor(t + h * (4 / 5))), len(forcing) - 1)
                ],
                values,
                nonnegative_soil,
            )
            f5 = _dp_rhs(
                rhs,
                z
                + h
                * (
                    19372 / 6561 * f1
                    - 25360 / 2187 * f2
                    + 64448 / 6561 * f3
                    - 212 / 729 * f4
                ),
                forcing[
                    min(int(np.floor(t + h * (8 / 9))), len(forcing) - 1)
                ],
                values,
                nonnegative_soil,
            )
            f6 = _dp_rhs(
                rhs,
                z
                + h
                * (
                    9017 / 3168 * f1
                    - 355 / 33 * f2
                    + 46732 / 5247 * f3
                    + 49 / 176 * f4
                    - 5103 / 18656 * f5
                ),
                forcing[min(int(np.floor(t + h)), len(forcing) - 1)],
                values,
                nonnegative_soil,
            )
            tn = final if done else t + h
            h = tn - t
            zn = z + h * (
                35 / 384 * f1
                + 500 / 1113 * f3
                + 125 / 192 * f4
                - 2187 / 6784 * f5
                + 11 / 84 * f6
            )
            f7 = _dp_rhs(
                rhs,
                zn,
                forcing[min(int(np.floor(tn)), len(forcing) - 1)],
                values,
                nonnegative_soil,
            )
            errvec = (
                71 / 57600 * f1
                - 71 / 16695 * f3
                + 71 / 1920 * f4
                - 17253 / 339200 * f5
                + 22 / 525 * f6
                - 1 / 40 * f7
            )
            err = absh * np.max(
                np.abs(errvec)
                / np.maximum(np.maximum(np.abs(z), np.abs(zn)), atol / rtol)
            )
            positivity_reject = False
            if nonnegative_soil and err <= rtol:
                negative_error = max(0.0, -zn[1], -zn[2]) / (atol / rtol)
                if negative_error > rtol:
                    err = negative_error
                    positivity_reject = True
            if not np.all(np.isfinite(zn)) or np.max(np.abs(zn)) > 1e12:
                raise ValueError("Model integration failed")
            if err <= rtol:
                break
            if absh <= hmin:
                raise ValueError("Dormand-Prince tolerance not met")
            absh = (
                max(hmin, absh * max(0.1, 0.8 * (rtol / err) ** 0.2))
                if first_try and not positivity_reject
                else max(hmin, 0.5 * absh)
            )
            h = absh
            first_try = False
            done = False
        reset_derivative = False
        if nonnegative_soil:
            for index in range(1, 3):
                if zn[index] < 0.0:
                    zn[index] = 0.0
                    reset_derivative = True
        while next_output <= int(final) and next_output <= tn:
            if next_output == tn:
                history[next_output] = zn
            else:
                x = (next_output - t) / h
                coefficients = np.empty(7)
                for j in range(7):
                    coefficients[j] = x * _DENSE[j, 0] + x * x * (
                        _DENSE[j, 1] + x * (_DENSE[j, 2] + x * _DENSE[j, 3])
                    )
                history[next_output] = z + h * (
                    f1 * coefficients[0]
                    + f3 * coefficients[2]
                    + f4 * coefficients[3]
                    + f5 * coefficients[4]
                    + f6 * coefficients[5]
                    + f7 * coefficients[6]
                )
            if nonnegative_soil:
                for index in range(1, 3):
                    history[next_output, index] = max(
                        history[next_output, index], 0.0
                    )
            next_output += 1
        t = tn
        z = zn
        f1 = (
            _dp_rhs(
                rhs,
                z,
                forcing[min(int(np.floor(t)), len(forcing) - 1)],
                values,
                nonnegative_soil,
            )
            if reset_derivative
            else f7
        )
        if first_try:
            growth_denominator = 1.25 * (err / rtol) ** 0.2
            if growth_denominator > 0.2:
                absh = absh / growth_denominator
            else:
                absh = 5.0 * absh
    for index in range(next_output, len(history)):
        history[index] = history[next_output - 1]
    return history


@njit(cache=True, error_model="numpy")
def dormand_prince(
    rhs, forcing, z, values, settings, final, nonnegative_soil=False
):
    """Integrate piecewise-constant forcing on separate unit intervals.

    Restart the DP controller and its FSAL derivative at each forcing
    change. State and sensitivity values carry forward without resetting.
    Local solver time avoids large-calendar-time rounding in step sizes.
    """
    if final != len(forcing):
        raise ValueError("Final time must match the forcing interval count")
    history = np.empty((len(forcing) + 1, len(z)))
    history[0] = z
    for interval in range(len(forcing)):
        segment = _dormand_prince_interval(
            rhs,
            forcing[interval : interval + 1],
            z,
            values,
            settings,
            1.0,
            nonnegative_soil,
        )
        z = segment[-1].copy()
        history[interval + 1] = z
    return history
