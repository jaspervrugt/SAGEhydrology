# ==========================================================================
# SAGEhydrology | Python source
# File: gui_python/sage_science.py
#
# Array translations of SAGE diagnostics/metrics; no training frameworks.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Array translations of SAGE diagnostics/metrics; no training
frameworks.

References: src/observation_loss.m, src/delta_n.m,
src/build_split.m,
utils/metrics/{huber_loss,fdc_metrics_cached}.m,
utils/info/sage_attribution.m,
utils/info/sage_information_bottleneck_matlab.m (KSG-1), and
utils/diagnostics/analyze_discharge_differencing.m.
"""

from sage_objective_error import UndefinedObjective
from math import comb
import numpy as np


def observation_loss(y, q, method=2, formulation=1):
    y = np.asarray(y, float)
    q = np.asarray(q, float)
    n = len(y)
    r = y - q
    if not n:
        raise UndefinedObjective("No valid observations in the selected period.")
    if method == 1:
        return float(np.abs(r).sum()), -np.sign(r)
    if method == 2:
        return float(r @ r), -2 * r
    if method == 3:
        tss = np.sum((y - y.mean()) ** 2)
        if tss <= 0:
            raise UndefinedObjective("NSE requires variable observations.")
        return float(r @ r / tss), -2 * r / tss
    if method == 4:
        yc = (y - y[0]) - np.mean(y - y[0])
        qc = (q - q[0]) - np.mean(q - q[0])
        sy = np.linalg.norm(yc)
        sq = np.linalg.norm(qc)
        if min(sy, sq) <= 0 or y.mean() == 0:
            raise UndefinedObjective(
                "KGE requires variable, nonzero-mean observations and"
                " simulations."
            )
        corr = np.clip(yc @ qc / (sy * sq), -1, 1)
        alpha = sq / sy
        beta = q.mean() / y.mean()
        value = np.sqrt(
            (corr - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2
        )
        if value == 0:
            return 0.0, np.zeros(n)
        dr = yc / (sy * sq) - corr * qc / (sq * sq)
        return (
            float(value),
            (
                (corr - 1) * dr
                + (alpha - 1) * qc / (sy * sq)
                + (beta - 1) / (n * y.mean())
            )
            / value,
        )
    if method == 5:
        scale = 1.482602218505602 * np.median(np.abs(y - np.median(y)))
        if not np.isfinite(scale) or scale <= 0:
            scale = np.median(y[y > 0]) if np.any(y > 0) else 0
        scale = max(scale, 1e-6, 1e-6 * np.mean(np.abs(y)))
        u = r / scale
        v = np.minimum(np.abs(u), 1.345)
        return (
            float(np.sum(v * (np.abs(u) - 0.5 * v))),
            -np.sign(u) * v / scale,
        )
    if method == 6:
        ys = np.sort(y)
        order = np.argsort(q, kind="stable")
        qs = q[order]
        delta = np.empty(n)
        if formulation in (2, 3):
            if formulation == 2:
                diff = qs - ys
                derivative = 2 * diff / n
            else:
                if not np.any(ys > 0):
                    raise UndefinedObjective(
                        "Log-FDC requires positive observed discharge."
                    )
                q0 = 0.01 * np.median(ys[ys > 0])
                if np.any(qs + q0 <= 0):
                    raise UndefinedObjective("Invalid log-FDC discharge.")
                diff = np.log(qs + q0) - np.log(ys + q0)
                derivative = 2 * diff / (n * (qs + q0))
            delta[order] = derivative
            return float(np.mean(diff * diff)), delta
        # Integrate the squared empirical-CDF difference directly. The
        # equivalent pair-distance subtraction loses precision near equality.
        combined = np.r_[qs, ys]
        merge = np.argsort(combined, kind="stable")
        balance = np.cumsum(np.where(merge < n, 1, -1))[:-1] / n
        value = np.sum(np.diff(combined[merge]) * balance * balance)

        def signs(sorted_values):
            return (
                np.searchsorted(sorted_values, q, side="left")
                + np.searchsorted(sorted_values, q, side="right")
                - n
            )

        return float(value), (signs(ys) - signs(qs)) / (n * n)
    raise ValueError("Unknown observation loss.")


def benchmark(values, method, window, dates, support=None):
    """Matrix-free B and B' from jkge_benchmark.m / jkge_transpose.m."""
    v = np.asarray(values, float)
    n = len(v)
    good = np.isfinite(v)
    cache_good = good if support is None else np.asarray(support, bool)
    if method == 1:
        half = (window - 1) // 2
        left = np.maximum(0, np.arange(n) - half)
        right = np.minimum(n, np.arange(n) + half + 1)
        counts = np.r_[0, np.cumsum(cache_good.astype(int))]
        den = counts[right] - counts[left]

        def apply(x):
            cs = np.r_[0, np.cumsum(np.where(good, x, 0))]
            return np.divide(
                cs[right] - cs[left],
                den,
                out=np.full(n, np.nan),
                where=den > 0,
            )

        def transpose(x):
            w = np.divide(x, den, out=np.zeros(n), where=den > 0)
            acc = np.zeros(n + 1)
            np.add.at(acc, left, w)
            np.add.at(acc, right, -w)
            return np.where(good, np.cumsum(acc[:-1]), 0)

    else:
        if method == 2:
            group = np.arange(n) // window
        elif method == 3:
            group = np.zeros(n, int)
        elif method == 4:
            group = (
                np.asarray(dates).astype("datetime64[M]").astype(int) % 12
            )
        else:
            raise ValueError("Unknown JKGE benchmark operator.")
        count = np.bincount(
            group[cache_good], minlength=int(group.max()) + 1
        )

        def reduce(x, selected):
            sums = np.bincount(
                group[selected],
                weights=np.asarray(x)[selected],
                minlength=len(count),
            )
            return np.divide(
                sums,
                count,
                out=np.full(len(count), np.nan),
                where=count > 0,
            )[group]

        def apply(x):
            if method == 3:
                return np.full(n, np.mean(np.asarray(x)[good]))
            selected = cache_good if method == 4 else good
            sums = np.bincount(
                group[selected],
                weights=np.asarray(x)[selected],
                minlength=len(count),
            )
            result = reduce(x, selected)
            return np.where(cache_good, result, np.nan) if method == 4 else result

        def transpose(x):
            selected = cache_good if method == 4 else good
            # The forward operators 1--3 use every finite simulation value,
            # even where observations are missing. Differentiate that map.
            if method == 3:
                return np.where(good, np.sum(np.where(np.isfinite(x), x, 0)) / np.sum(good), 0)
            source = cache_good if method == 4 else np.ones(n, dtype=bool)
            answer = reduce(np.where(np.isfinite(x), x, 0), source)
            return np.where(selected & np.isfinite(answer), answer, 0)

    return apply(v), apply, transpose


def full_loss(y, q, mask, c, dates, components=False):
    """
    Return complete-series dL/dq, retaining JKGE's benchmark chain
    rule.
    """
    y = np.asarray(y, float)
    q = np.asarray(q, float)
    mask = (
        np.asarray(mask, bool) & np.isfinite(y) & (y >= 0) & np.isfinite(q)
    )
    from sage_zscore import transform, loss as zscore_loss

    response = c.get("_loss_response", "Q")
    if (
        transform(c, response) == "zscore"
        and not components
        and not c.get("_physical_diagnostics")
    ):
        value, partial = zscore_loss(y[mask], q[mask])
        delta = np.zeros(len(y))
        delta[mask] = partial
        return value, delta
    if c["loss"] != 7:
        value, d = observation_loss(
            y[mask], q[mask], c["loss"], c["fdc_formulation"]
        )
        delta = np.zeros(len(y))
        delta[mask] = d
        return value, delta
    # MATLAB forms JKGE benchmarks over the scored record, excluding
    # spin-up. Including spin-up changes the earliest period's moving mean.
    if (
        not c.get("_jkge_scored")
        and c.get("train_start")
        and c.get("eval_start")
    ):
        start = (
            min(c["train_start"], c["eval_start"])
            if c.get("split", "manual") == "manual"
            else c.get("common_start", c["train_start"])
        )
        keep = np.asarray(dates) >= np.datetime64(start)
        if not keep.all():
            answer = full_loss(
                y[keep],
                q[keep],
                mask[keep],
                dict(c, _jkge_scored=True),
                np.asarray(dates)[keep],
                components,
            )
            if components:
                return answer
            value, partial = answer
            delta = np.zeros(len(y))
            delta[keep] = partial
            return value, delta
    method = c["jkge_operator"]
    window = int(c["jkge_window"]) * int(c["resolution"])
    if method == 1 and window % 2 == 0:
        window += 1  # Same centered odd-window rule as MATLAB prep_stats.
    observed = np.where(y >= 0, y, np.nan)
    ranges = c.get("unscored_ranges", "")
    if ranges.strip():
        from input_period import unscored_mask

        observed[unscored_mask(dates, ranges)] = np.nan
    my, _, _ = benchmark(observed, method, window, dates)
    # prep_stats.m passes the observed-record cache through both paths.
    mq, _, transpose = benchmark(
        q, method, window, dates, support=np.isfinite(observed)
    )
    mask = mask & np.isfinite(my) & np.isfinite(mq)
    n = np.sum(mask)
    if n < 2:
        raise UndefinedObjective("JKGE needs at least two valid samples.")
    yr = y[mask]
    ay = yr - my[mask]
    aq = q[mask] - mq[mask]
    if method == 3:
        fy, fq = y[np.isfinite(y)], q[np.isfinite(q)]
        ay = (y[mask] - fy[0]) - np.mean(fy - fy[0])
        aq = (q[mask] - fq[0]) - np.mean(fq - fq[0])
    sy = np.sqrt(np.mean(ay * ay))
    sq = np.sqrt(np.mean(aq * aq))
    if min(sy, sq) <= 0:
        raise UndefinedObjective("JKGE requires nonzero residual variability.")
    zy = ay / sy
    rho = np.clip(np.mean(aq / sq * zy), -1, 1)
    alpha = sq / sy
    dm = np.zeros(n)
    if c["jkge_m"] == 1:
        m = my[mask]
        valid = np.abs(m) > max(1e-8, 1e-6 * np.mean(np.abs(yr)))
        if not np.any(valid):
            raise UndefinedObjective("JKGE benchmark is zero.")
        err = 1 - mq[mask][valid] / m[valid]
        M = np.mean(err * err)
        dm[valid] = -2 * err / (np.sum(valid) * m[valid])
    else:
        den = np.sum((my[mask] - yr.mean()) ** 2)
        if not np.isfinite(den) or den <= 0:
            raise UndefinedObjective(
                "Norm-based JKGE is undefined for a constant benchmark;"
                " choose another operator or M definition."
            )
        diff = mq[mask] - my[mask]
        M = diff @ diff / den
        dm = 2 * diff / den
    if components:
        return dict(
            jkge_m=float(M),
            jkge_v=float((1 - alpha) ** 2),
            jkge_c=float((1 - rho) ** 2),
        )
    value = float(np.sqrt(M + (1 - alpha) ** 2 + (1 - rho) ** 2))
    if value == 0:
        return 0.0, np.zeros(len(y))
    dv = -2 * (1 - alpha) * aq / (n * sq * sy)
    dr = zy / (n * sq) - np.sum(zy * aq) * aq / (n * n * sq**3)
    dc = -2 * (1 - rho) * dr
    direct = np.zeros(len(y))
    indirect = np.zeros(len(y))
    direct[mask] = (dv + dc) / (2 * value)
    indirect[mask] = (dm - dv - dc) / (2 * value)
    return value, direct + transpose(indirect)


def diagnostic_scores(y, q, mask, c, dates, baseline_cache=None):
    """
    Scenario metrics used by MATLAB's Metrics/ECDF/Generalization
    tabs.
    """
    mask = (
        np.asarray(mask, bool)
        & np.isfinite(y)
        & np.isfinite(q)
        & (np.asarray(y) >= 0)
    )
    result = {}
    if not mask.any():
        return result
    yt = np.asarray(y)[mask]
    qt = np.asarray(q)[mask]
    if baseline_cache is not None:
        signature = (mask.tobytes(), yt.tobytes())
        if baseline_cache.get("signature") != signature:
            baseline_cache.clear()
            baseline_cache["signature"] = signature
    for form, name in ((1, "s_fdc"), (2, "s_p"), (3, "s_logp")):
        try:
            simulated = qt
            if form == 3 and c.get("_loss_response", "Q") == "Q":
                settings = c.get("settings")
                if settings is not None and len(settings) > 4:
                    # Discharge is the difference of two cumulative states.
                    # The local solver tolerance is not a bound for the
                    # difference. Allow a small guard band for accumulated
                    # error in the log-FDC diagnostic only.
                    tolerance = 3 * settings[4]
                    minimum = qt.min()
                    if np.isfinite(tolerance) and -tolerance <= minimum < 0:
                        simulated = np.maximum(qt, 0)
            value, _ = observation_loss(yt, simulated, 6, form)
            if baseline_cache is not None and form in baseline_cache:
                baseline = baseline_cache[form]
            else:
                baseline, _ = observation_loss(
                    yt, np.full(len(yt), np.median(yt)), 6, form
                )
                if baseline_cache is not None:
                    baseline_cache[form] = baseline
            result[name] = (
                float(1 - value / baseline) if baseline > 0 else None
            )
            result[{1: "d_fdc", 2: "d_p", 3: "d_logp"}[form]] = float(value)
        except ValueError:
            result[name] = None
    for method, key in ((1, "sar"), (2, "gls"), (5, "huber")):
        try:
            result[key] = float(observation_loss(yt, qt, method)[0])
        except ValueError:
            result[key] = None
    if np.std(yt) > 0 and np.std(qt) > 0:
        result["kge_r"] = float(np.corrcoef(yt, qt)[0, 1])
    if np.std(yt) > 0:
        result["kge_alpha"] = float(np.std(qt) / np.std(yt))
    if np.mean(yt) != 0:
        result["kge_beta"] = float(np.mean(qt) / np.mean(yt))
    try:
        loss, _ = full_loss(
            y, q, mask, dict(c, loss=7, _physical_diagnostics=True), dates
        )
        result["jkge"] = float(1 - loss)
        if c["loss"] == 7:
            result.update(
                full_loss(
                    y,
                    q,
                    mask,
                    dict(c, _physical_diagnostics=True),
                    dates,
                    components=True,
                )
            )
    except ValueError:
        result["jkge"] = None
    return result


def temporal_masks(dates, c, precipitation=None):
    """Apply scoring exclusions after assigning the temporal split."""
    from input_period import unscored_mask

    train, evaluation = _temporal_split(dates, c, precipitation)
    ranges = c.get("unscored_ranges", "")
    if not ranges.strip():
        return train, evaluation
    excluded = unscored_mask(dates, ranges)
    if np.any(train) and not np.any(train & ~excluded):
        raise ValueError("Unscored ranges remove the entire training period.")
    if np.any(evaluation) and not np.any(evaluation & ~excluded):
        raise ValueError("Unscored ranges remove the entire evaluation period.")
    return train & ~excluded, evaluation & ~excluded


def _temporal_split(dates, c, precipitation=None):
    """
    Scored-index split rules from build_split.m, with MATLAB
    permutations.
    """
    from matlab_random import split_permutation
    from sage_temporal_validation import validate_design
    validate_design(c)

    method = c["split"]
    dates = np.asarray(dates)
    train = np.zeros(len(dates), bool)
    if method == "manual":

        def period(start, end):
            a = np.datetime64(c[start])
            b = np.datetime64(c[end]) + np.timedelta64(1, "D")
            if b <= a:
                raise ValueError(
                    "Period end must be on or after its start."
                )
            return (dates >= a) & (dates < b)

        train = period("train_start", "train_end")
        evaluation = np.zeros(len(dates), bool)
        if c.get("evaluation_enabled", True):
            evaluation = period("eval_start", "eval_end")
            if np.any(train & evaluation):
                raise ValueError(
                    "Training and evaluation periods overlap; use separate"
                    " dates or disable evaluation."
                )
        return train, evaluation
    pool = np.flatnonzero(
        (dates >= np.datetime64(c["common_start"]))
        & (dates < np.datetime64(c["common_end"]) + np.timedelta64(1, "D"))
    )
    n = len(pool)
    if n < 2:
        raise ValueError(
            "The common period needs at least two scored samples."
        )
    fraction = c["train_fraction"]
    round_positive = lambda v: int(np.floor(v + 0.5))
    if (
        method in ("random", "deterministic_block", "random_block")
        and not 0 < fraction < 1
    ):
        raise ValueError("Training fraction must be between 0 and 1.")
    if method == "random" or method in ("deterministic_kfold", "random_kfold"):
        # Calendar days are indivisible: every subdaily sample inherits
        # its day's assignment. The same seed yields the same day split
        # at daily, hourly and 15-minute resolution.
        days, day_index = np.unique(dates[pool].astype("datetime64[D]"), return_inverse=True)
        day_count = len(days)
        if day_count < 2:
            raise ValueError("Select at least two calendar days for this split method.")
        if method == "random":
            order = split_permutation(day_count, c)
            train_days = order[:max(1,min(day_count-1,round_positive(fraction*day_count)))]
        else:
            if not 2 <= c["folds"] <= day_count or not 1 <= c["fold"] <= c["folds"]:
                raise ValueError("Fold count must be between 2 and the number of calendar days; held-out fold must be within that count.")
            order = split_permutation(day_count,c) if method == "random_kfold" else np.arange(day_count)
            edges = np.floor(np.linspace(0,day_count,c["folds"]+1)+.5).astype(int)
            evaluation_days = order[edges[c["fold"]-1]:edges[c["fold"]]]
            train_days = np.setdiff1d(np.arange(day_count),evaluation_days)
        chosen = np.flatnonzero(np.isin(day_index,train_days))
    elif method in ("deterministic_block", "random_block"):
        # MATLAB block is a number of scored samples, including subdaily
        # data.
        block = int(c["block_days"])
        if block < 1:
            raise ValueError(
                "Block size must be a positive number of samples."
            )
        blocks = [
            np.arange(i, min(i + block, n)) for i in range(0, n, block)
        ]
        if c.get("block_unit") == "water_years":
            from sage_water_year_blocks import water_year_blocks
            blocks = water_year_blocks(dates[pool], c)
            block = min(len(part) for part in blocks)
        if len(blocks) < 2:
            raise ValueError("The common period needs at least two blocks.")
        if method == "deterministic_block":
            if fraction * n < block:
                raise ValueError("Training fraction is too small for one complete block; increase the fraction or reduce block size.")
            # Select the eligible whole-block boundary closest to the
            # requested sample count. A short final block must not carry
            # the same weight as a complete block.
            boundaries = np.cumsum([len(part) for part in blocks])[:-1]
            boundary = int(np.argmin(np.abs(boundaries - fraction * n))) + 1
            chosen = np.concatenate(blocks[:boundary])
        else:
            order = split_permutation(len(blocks), c)
            selected = []
            total = 0
            for index in order:
                if total < round_positive(fraction * n):
                    selected.append(index)
                    total += len(blocks[index])
            if len(selected) == len(blocks):
                selected.remove(order[-1])
            if not selected:
                selected = [order[0]]
            chosen = np.concatenate([blocks[i] for i in selected])
    elif method == "rainfall_block":
        if precipitation is None:
            raise ValueError(
                "Load a basin rainfall preview to rank its water years."
            )
        p = np.asarray(precipitation, float)
        if p.shape != dates.shape:
            raise ValueError("Precipitation must match the split timeline.")
        years = dates[pool].astype("datetime64[Y]").astype(int) + 1970
        months = dates[pool].astype("datetime64[M]").astype(int) % 12 + 1
        wy = years + (months >= 10)
        unique = np.unique(wy)
        if len(unique) < 2:
            raise ValueError(
                "Rainfall-block splitting needs at least two water years"
                " (October–September)."
            )
        totals = np.array(
            [np.nansum(p[pool][wy == year]) for year in unique]
        )
        order = np.argsort(-totals, kind="stable")
        chosen = np.flatnonzero(np.isin(wy, unique[order[::2]]))
    else:
        raise ValueError(f"Unknown temporal split: {method}")
    train[pool[chosen]] = True
    valid = np.zeros(len(dates), bool)
    valid[pool] = True
    return train, valid & ~train


def attribution(jacobian, delta, lower, upper, parameter_derivative=None):
    contributions = (
        np.asarray(jacobian)
        * np.asarray(delta)[:, None]
        * ((np.asarray(upper) - lower) if parameter_derivative is None else parameter_derivative)
    )
    return np.sum(np.abs(contributions), axis=0), np.abs(
        np.sum(contributions, axis=0)
    )


def ksg_mi(x, y, k=5, jitter=1e-10, max_dim=4):
    """Literal KSG-1 counting/jitter rules from SAGE's MATLAB reference."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.ndim == 1:
        x = x[:, None]
    if y.ndim == 1:
        y = y[:, None]
    good = np.isfinite(x).all(axis=1) & np.isfinite(y).all(axis=1)
    x = x[good].copy()
    y = y[good].copy()
    n = len(x)
    if n <= k + 1 or x.shape[1] + y.shape[1] > max_dim:
        return np.nan
    for values, offset in ((x, 1), (y, 101)):
        sd = values.std(axis=0, ddof=1)
        sd[~np.isfinite(sd) | (sd <= 0)] = 1
        values -= values.mean(axis=0)
        values /= sd
        for j in range(values.shape[1]):
            s = values[:, j].std(ddof=1)
            if not np.isfinite(s) or s <= 0:
                s = 1
            values[:, j] += (
                jitter
                * s
                * np.sin(
                    np.arange(1, n + 1) * (j + 1 + offset + np.sqrt(2))
                )
            )
    # psi at positive integers: H_(n-1) - Euler's constant; constants
    # cancel.
    harmonic = np.r_[0, np.cumsum(1 / np.arange(1, n + 1))]
    terms = []
    for i in range(n):
        dx = np.max(np.abs(x - x[i]), axis=1)
        dy = np.max(np.abs(y - y[i]), axis=1)
        dx[i] = np.inf
        dy[i] = np.inf
        epsilon = np.partition(np.maximum(dx, dy), k - 1)[k - 1]
        radius = max(
            0.0, epsilon - max(1e-12, 32 * np.spacing(max(1.0, epsilon)))
        )
        terms.append(
            harmonic[np.sum(dx < radius)] + harmonic[np.sum(dy < radius)]
        )
    value = (harmonic[k - 1] + harmonic[n - 1] - np.mean(terms)) / np.log(2)
    return 0.0 if -1e-6 < value < 0 else float(value)


def discharge_differencing(q, half_window=100):
    q = np.asarray(q, float).copy()
    q[~np.isfinite(q) | (q < 0)] = np.nan
    variance = np.diff(q, n=3) ** 2 / comb(6, 3)
    discharge = np.mean(
        np.lib.stride_tricks.sliding_window_view(q, 4), axis=1
    )
    good = (
        np.isfinite(discharge)
        & np.isfinite(variance)
        & (np.sqrt(variance) >= 1e-10)
    )
    discharge = discharge[good]
    variance = variance[good]
    order = np.argsort(discharge, kind="stable")
    width = 2 * half_window + 1
    if len(discharge) < width:
        raise ValueError(
            f"At least {width} valid differenced samples are required."
        )
    smooth_q = discharge[order][half_window : len(discharge) - half_window]
    smooth_sigma = np.sqrt(
        np.convolve(variance[order], np.ones(width) / width, mode="valid")
    )
    mean = np.nanmean(q)
    X = np.column_stack([smooth_q, np.full(len(smooth_q), mean)])
    # analyze_discharge_differencing/local_nnls2 checks the unconstrained
    # solution and each boundary of the two-dimensional nonnegative
    # quadrant.
    candidates = [
        np.zeros(2),
        np.array(
            [max(0.0, X[:, 0] @ smooth_sigma / (X[:, 0] @ X[:, 0])), 0.0]
        ),
        np.array(
            [0.0, max(0.0, X[:, 1] @ smooth_sigma / (X[:, 1] @ X[:, 1]))]
        ),
    ]
    fit = np.linalg.lstsq(X, smooth_sigma, rcond=None)[0]
    if np.all(fit >= 0):
        candidates.append(fit)
    fit = min(candidates, key=lambda p: np.sum((smooth_sigma - X @ p) ** 2))
    return (
        discharge,
        np.sqrt(variance),
        smooth_q,
        smooth_sigma,
        fit[0],
        fit[1] * mean,
    )


def variogram(latlon, values, max_range=1500.0, bins=30):
    ll = np.deg2rad(np.asarray(latlon, float))
    values = np.asarray(values, float)
    i, j = np.triu_indices(len(values), 1)
    a = (
        np.sin((ll[i, 0] - ll[j, 0]) / 2) ** 2
        + np.cos(ll[i, 0])
        * np.cos(ll[j, 0])
        * np.sin((ll[i, 1] - ll[j, 1]) / 2) ** 2
    )
    distance = 2 * 6371 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
    semivariance = 0.5 * (values[i] - values[j]) ** 2
    good = (
        np.isfinite(distance)
        & np.isfinite(semivariance)
        & (distance <= 1.15 * max_range)
    )
    edges = np.linspace(0, 1.15 * max_range, bins + 1)
    counts = np.histogram(distance[good], edges)[0]
    sums = np.histogram(distance[good], edges, weights=semivariance[good])[
        0
    ]
    gamma = np.divide(
        sums, counts, out=np.full(bins, np.nan), where=counts > 0
    )
    return (edges[:-1] + edges[1:]) / 2, gamma, counts
