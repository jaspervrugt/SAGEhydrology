# ==========================================================================
# SAGEhydrology | Python source
# File: python/site_hymod.py
#
# Daily, snow-enabled SITE/SAGE HYMOD and SITE algorithm 5.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Daily, snow-enabled SITE/SAGE HYMOD and SITE algorithm 5.

The equations and analytic sensitivities follow SAGEhydrology/models/hymod.
Backends: Numba compiles Python loops; cpp calls the existing SAGE C++ core.
"""

from dataclasses import dataclass
import numpy as np
from numba import njit

LOWER = np.array([50.0, 0.1, 1e-5, 1e-5, 0.1, -3.0, 0.1])
UPPER = np.array([5000.0, 10.0, 1.0, 1.0, 100.0, 3.0, 10.0])
PARAMETERS = ("s_umax", "beta", "alfa", "K_s", "K_f", "T_tr", "f_dd")


@njit(cache=True)
def _rhs(z, th, P, Ep, T, config):
    """Evaluate HYMOD state derivatives and forward sensitivities."""
    cap, beta, alpha, ks, kf, threshold, dd = th
    width, eps, rho = config
    swe0, su, ss, f1, f2, f3 = z[:6, 0]
    swe = 0.5 * (swe0 + np.sqrt(swe0 * swe0 + eps * eps))
    dswe = 0.5 * (1 + swe0 / np.sqrt(swe0 * swe0 + eps * eps))
    width = max(width, eps)
    a = T - threshold
    tan = np.tanh(a / width)
    snow = 0.5 * (1 - tan)
    pos = 0.5 * (a + np.sqrt(a * a + width * width))
    pot = dd * pos
    delta = swe - pot
    root = np.sqrt(delta * delta + eps * eps)
    melt = 0.5 * (swe + pot - root)
    ms = 0.5 * (1 - delta / root)
    mp = 0.5 * (1 + delta / root)
    liquid = P * (1 - snow) + melt
    ratio = su / cap
    r1 = np.sqrt((ratio - 1) ** 2 + 1e-12)
    r0 = np.sqrt(ratio**2 + 1e-12)
    ur = ratio - 0.5 * (ratio - 1 + r1) + 0.5 * (-ratio + r0)
    dur = 1 - 0.5 * (1 + (ratio - 1) / r1) + 0.5 * (-1 + ratio / r0)
    raw = 1 - ur
    rd = np.sqrt((raw - 1e-6) ** 2 + 1e-20)
    deficit = 0.5 * (raw + 1e-6 + rd)
    ddef = 0.5 * (1 + (raw - 1e-6) / rd)
    power = np.exp(beta * np.log(deficit))
    runoff_factor = 1 - power
    qu = liquid * runoff_factor
    ea = Ep * ur * (1 + rho) / (ur + rho)
    qu_ur = liquid * power * beta / deficit * ddef
    ea_ur = Ep * (1 + rho) * rho / (ur + rho) ** 2
    qu_su = qu_ur * dur / cap
    ea_su = ea_ur * dur / cap
    qu_swe = runoff_factor * ms
    fx = np.zeros((7, 7))
    fx[0, 0] = -ms
    fx[1, 0] = ms - qu_swe
    fx[1, 1] = -qu_su - ea_su
    fx[2, 0], fx[2, 1], fx[2, 2] = (
        (1 - alpha) * qu_swe,
        (1 - alpha) * qu_su,
        -ks,
    )
    fx[3, 0], fx[3, 1], fx[3, 3] = alpha * qu_swe, alpha * qu_su, -kf
    fx[4, 3], fx[4, 4] = kf, -kf
    fx[5, 4], fx[5, 5] = kf, -kf
    fx[6, 2], fx[6, 5] = ks, kf
    fx[:, 0] *= dswe
    ft = np.zeros((7, 7))
    qc = qu_ur * dur * (-su / cap**2)
    ec = ea_ur * dur * (-su / cap**2)
    ft[1, 0], ft[2, 0], ft[3, 0] = -qc - ec, (1 - alpha) * qc, alpha * qc
    qb = -liquid * power * np.log(deficit)
    ft[1, 1], ft[2, 1], ft[3, 1] = -qb, (1 - alpha) * qb, alpha * qb
    ft[2, 2], ft[3, 2] = -qu, qu
    ft[2, 3], ft[6, 3] = -ss, ss
    ft[3, 4], ft[4, 4], ft[5, 4], ft[6, 4] = -f1, f1 - f2, f2 - f3, f3
    snow_t = P * 0.5 * (1 - tan * tan) / width
    melt_t = -mp * dd * 0.5 * (1 + a / np.sqrt(a * a + width * width))
    liq_t = -snow_t + melt_t
    ft[0, 5], ft[1, 5] = snow_t - melt_t, liq_t * (1 - runoff_factor)
    ft[2, 5], ft[3, 5] = (
        1 - alpha
    ) * runoff_factor * liq_t, alpha * runoff_factor * liq_t
    melt_dd = mp * pos
    ft[0, 6], ft[1, 6] = -melt_dd, melt_dd * (1 - runoff_factor)
    ft[2, 6], ft[3, 6] = (
        1 - alpha
    ) * runoff_factor * melt_dd, alpha * runoff_factor * melt_dd
    out = np.empty((7, 8))
    for i in range(7):
        for j in range(7):
            value = ft[i, j]
            for k in range(7):
                value += fx[i, k] * z[k, j + 1]
            out[i, j + 1] = value
    out[0, 0], out[1, 0], out[2, 0] = (
        P * snow - melt,
        liquid - qu - ea,
        (1 - alpha) * qu - ks * ss,
    )
    out[3, 0], out[4, 0], out[5, 0], out[6, 0] = (
        alpha * qu - kf * f1,
        kf * (f1 - f2),
        kf * (f2 - f3),
        ks * ss + kf * f3,
    )
    return out


@njit(cache=True)
def _integrate(theta, forcing, y0, config, settings):
    """Integrate the HYMOD augmented state over the forcing record."""
    hin, hmax, hmin, rtol, atol, order, maxreject = settings
    z = np.zeros((7, 8))
    z[:, 0] = y0
    q = np.empty(len(forcing))
    jac = np.empty((len(forcing), 7))
    carry = max(hmin, min(hin, hmax))
    for day in range(len(forcing)):
        prev = z[6, :].copy()
        t, end = float(day), float(day + 1)
        h = min(carry, 1.0)
        rejected = 0
        P, Ep, T = forcing[day]
        while t < end:
            e = _rhs(z, theta, P, Ep, T, config)
            ze = z + h * e
            zh = z + 0.5 * h * (e + _rhs(ze, theta, P, Ep, T, config))
            if not np.all(np.isfinite(zh)) or np.max(np.abs(zh)) > 1e12:
                raise ValueError("HYMOD integration failed")
            wrms = np.sqrt(
                np.mean(((ze - zh) / (rtol * np.abs(zh) + atol)) ** 2)
            )
            if wrms <= 1 or h <= hmin:
                z = zh
                t += h
                rejected = 0
            else:
                rejected += 1
            fac = (
                5.0
                if wrms <= 0
                else min(5.0, max(0.2, 0.9 * wrms ** (-1 / order)))
            )
            carry = max(hmin, min(h * fac, hmax))
            tol = 10 * np.finfo(np.float64).eps * max(1.0, abs(end), abs(t))
            if end - t <= tol:
                break
            h = min(carry, end - t)
            if h < tol or rejected >= maxreject:
                raise ValueError("HYMOD step limit reached")
        q[day] = z[6, 0] - prev[0]
        jac[day, :] = z[6, 1:] - prev[1:]
    return q, jac


def simulate(
    theta,
    forcing,
    *,
    y0=None,
    config=(1.0, 1e-6, 0.01),
    backend="numba",
    settings=(0.01, 1.0, 1e-4, 1e-3, 1e-3, 2.0, 10000.0)
):
    """Return daily discharge (mm/day) and its native-parameter Jacobian.

    Forcing columns are precipitation, PET (mm/day), temperature (deg C).
    backend='cpp' calls the C++ core through ctypes; default is 'numba'.
    Initial stores default to 1e-5 mm, as in SAGE. All seven states and
    49 sensitivity states control the adaptive Heun error estimate.
    """
    theta = np.asarray(theta, dtype=float)
    forcing = np.ascontiguousarray(forcing, dtype=float)
    y0 = np.full(7, 1e-5) if y0 is None else np.asarray(y0, dtype=float)
    if (
        theta.shape != (7,)
        or y0.shape != (7,)
        or forcing.ndim != 2
        or forcing.shape[1] != 3
    ):
        raise ValueError(
            "Expected seven parameters/states and an N-by-3 forcing array"
        )
    if not all(
        np.all(np.isfinite(v))
        for v in (theta, forcing, y0, config, settings)
    ):
        raise ValueError("Model inputs must be finite")
    if (
        theta[0] <= 0
        or theta[1] <= 0
        or np.any(theta[3:5] <= 0)
        or theta[6] <= 0
        or not 0 <= theta[2] <= 1
    ):
        raise ValueError("Invalid HYMOD parameters")
    if (
        len(settings) != 7
        or min(settings[:5]) <= 0
        or settings[2] > settings[1]
        or settings[5] <= 0
        or settings[6] < 1
    ):
        raise ValueError("Invalid solver settings")
    if len(config) != 3 or min(config) <= 0:
        raise ValueError("Smoothing settings must be positive")
    if backend == "cpp":
        from native_hymod import integrate

        return integrate(
            theta,
            forcing,
            y0,
            np.asarray(config, float),
            np.asarray(settings, float),
        )
    if backend != "numba":
        raise ValueError("backend must be numba or cpp")
    return _integrate(
        theta,
        forcing,
        y0,
        np.asarray(config, float),
        np.asarray(settings, float),
    )


def rss_objective(
    forcing, observed, train, *, lower=LOWER, upper=UPPER, **kwargs
):
    """Objective/analytic gradient in normalized [0,1] parameter space."""
    observed = np.asarray(observed, float)
    train = (
        np.asarray(train, bool) & np.isfinite(observed) & (observed >= 0)
    )
    if (
        observed.shape != (len(forcing),)
        or train.shape != observed.shape
        or not np.any(train)
    ):
        raise ValueError(
            "Need at least one valid training observation and matching"
            " lengths"
        )
    lower, upper = np.asarray(lower, float), np.asarray(upper, float)
    if lower.shape != (7,) or upper.shape != (7,) or np.any(lower >= upper):
        raise ValueError("Expected seven ordered parameter bounds")

    def fun(x):
        """Evaluate the local objective callback."""
        try:
            q, J = simulate(lower + x * (upper - lower), forcing, **kwargs)
        except ValueError:
            return 1e20, np.zeros(7)
        residual = q[train] - observed[train]
        return float(residual @ residual), 2 * (J[train].T @ residual) * (
            upper - lower
        )

    return fun


@dataclass
class Result:
    x: np.ndarray
    f: float
    iterations: int
    exitflag: int
    evaluations: int
    trajectory: np.ndarray


def algorithm5(
    fun,
    x0,
    *,
    max_iter=1000,
    memory=5,
    tol_pg=1e-6,
    tol_rel_f=1e-6,
    buffer=1e-6,
    c1=1e-4,
    ls_rho=0.5,
    ls_max=40,
    alpha0=1.0,
    iteration_callback=None
):
    """Port of SITE lbfgsb_minimize, not SciPy's different L-BFGS-B."""
    if max_iter < 0 or memory < 1 or not 0 <= buffer < 0.5:
        raise ValueError("Invalid optimizer settings")
    lo, hi = buffer, 1 - buffer

    def project(value):
        """Project an optimizer candidate onto the admissible bounds."""
        value = np.clip(value, lo, hi)
        value = np.where(
            abs(value - lo) <= 8 * np.spacing(max(1, abs(lo))), lo, value
        )
        return np.where(
            abs(value - hi) <= 8 * np.spacing(max(1, abs(hi))), hi, value
        )

    x = project(np.asarray(x0, float))
    f, g = fun(x)
    evaluations = 1
    traj = [f]
    history = []
    stall = 0
    flag = 0

    def direction(p):
        """Compute the bound-constrained optimizer search direction."""
        p = p.copy()
        p[((x <= lo) & (p < 0)) | ((x >= hi) & (p > 0))] = 0
        return p

    for k in range(1, max_iter + 1):
        active = ((x <= lo) & (g >= 0)) | ((x >= hi) & (g <= 0))
        pg_vec = g.copy()
        pg_vec[x <= lo] = np.minimum(0, g[x <= lo])
        pg_vec[x >= hi] = np.maximum(0, g[x >= hi])
        pg = np.linalg.norm(pg_vec)
        if pg < tol_pg:
            flag = 1
            break
        if k > 5 and abs(traj[-1] - traj[-5]) <= tol_rel_f * max(1, abs(f)):
            flag = 2
            break
        r = pg_vec.copy()
        coefficients = []
        for s, y, rho in reversed(history):
            a = rho * (s @ r)
            coefficients.append(a)
            r -= a * y
        if history:
            s, y, _ = history[-1]
            r *= np.clip((y @ s) / (y @ y), 1e-4, 1e4)
        for (s, y, rho), a in zip(history, reversed(coefficients)):
            r += s * (a - rho * (y @ r))
        p = direction(-r)
        if (
            np.linalg.norm(p) < 1e-14
            or not np.isfinite(g @ p)
            or g @ p >= 0
        ):
            history = []
            p = direction(-g)
            if np.linalg.norm(p) < 1e-14:
                flag = 3
                break
        gp = g @ p
        if abs(gp) < 1e-12 and np.linalg.norm(g) < 1e-6:
            flag = 5
            break

        def search(p):
            """Perform the optimizer line search."""
            nonlocal evaluations
            alpha = min(alpha0, 1 / max(1, np.linalg.norm(p, ord=np.inf)))
            previous = None
            if not np.isfinite(g @ p) or g @ p >= 0:
                return x, f, g, alpha, False
            for _ in range(ls_max):
                xn = project(x + alpha * p)
                step = xn - x
                if not np.any(step):
                    break
                slope = g @ step
                duplicate = previous is not None and np.array_equal(
                    xn, previous
                )
                previous = xn
                if not duplicate and np.isfinite(slope) and slope < 0:
                    try:
                        fn, gn = fun(xn)
                    except __import__(
                        "simulation_failure"
                    ).IntegrationFailure:
                        evaluations += 1
                        alpha *= ls_rho
                        continue
                    evaluations += 1
                    if (
                        np.isfinite(fn)
                        and np.all(np.isfinite(gn))
                        and fn <= f + c1 * slope
                    ):
                        return xn, fn, gn, alpha, True
                alpha *= ls_rho
            return x, f, g, alpha, False

        xn, fn, gn, alpha, ok = search(p)
        if not ok and history:
            history = []
            xn, fn, gn, alpha, ok = search(direction(-g))
        if not ok:
            flag = 4
            break
        s, y = xn - x, gn - g
        new_active = ((xn <= lo) & (gn >= 0)) | ((xn >= hi) & (gn <= 0))
        changed = not np.array_equal(active, new_active)
        y[new_active] = 0
        ys = y @ s
        sny = np.linalg.norm(s) * np.linalg.norm(y)
        if changed:
            history = []
            stall = 0
        elif (
            np.isfinite(ys)
            and np.isfinite(sny)
            and ys > 1e-10 * sny
            and ys > 1e-12
            and np.linalg.norm(y) > 1e-10
        ):
            history.append((s, y, 1 / ys))
            history = history[-memory:]
            stall = 0
        else:
            history = []
            stall += 1
        x, f, g = xn, fn, gn
        traj.append(f)
        if iteration_callback is not None:
            iteration_callback(
                dict(
                    iteration=len(traj) - 1,
                    objective=float(f),
                    evaluations=evaluations,
                    alpha=float(alpha),
                )
            )
        if stall >= 10:
            flag = 6
            break
        if np.linalg.norm(g) < 10 * tol_pg:
            flag = 7
            break
    return Result(
        x, float(f), len(traj) - 1, flag, evaluations, np.asarray(traj)
    )
