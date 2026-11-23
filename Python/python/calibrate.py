# ==========================================================================
# SAGEhydrology | Python source
# File: python/calibrate.py
#
# Standalone Python CAMELS-US daily HYMOD calibration using SITE
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Standalone Python CAMELS-US daily HYMOD calibration using SITE
algorithm 5.
"""

import argparse
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from camels_us import read_daily
from site_hymod import LOWER, UPPER, rss_objective, algorithm5, simulate

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--gauges", nargs="+", required=True)
    parser.add_argument(
        "--product", choices=("daymet", "maurer", "nldas"), default="nldas"
    )
    parser.add_argument("--pet", type=int, choices=(1, 2, 3), default=1)
    parser.add_argument(
        "--backend", choices=("numba", "cpp"), default="numba"
    )
    parser.add_argument("--starts", type=int, default=10)
    parser.add_argument("--max-iter", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261012)
    parser.add_argument(
        "--output", type=Path, default=Path("python_results.json")
    )
    args = parser.parse_args()
    if args.starts < 1 or args.max_iter < 1:
        parser.error("starts and max-iter must be positive")
    rng = np.random.default_rng(args.seed)
    # SITE's interior Latin hypercube; no cross-language RNG equivalence
    # assumed.
    starts = np.array(
        [
            0.1
            + 0.8
            * (rng.permutation(args.starts) + rng.random(args.starts))
            / args.starts
            for _ in range(7)
        ]
    ).T
    report = dict(
        product=args.product,
        pet=args.pet,
        backend=args.backend,
        seed=args.seed,
        algorithm=5,
        basins=[],
    )
    for gauge in args.gauges:
        data = read_daily(
            args.data_root, gauge, product=args.product, pet=args.pet
        )
        fun = rss_objective(
            data["forcing"],
            data["observed"],
            data["train"],
            backend=args.backend,
        )
        tick = perf_counter()
        fun(starts[0])
        first = perf_counter() - tick
        trials = []
        for x0 in starts:
            tick = perf_counter()
            out = algorithm5(fun, x0, max_iter=args.max_iter)
            trials.append(
                dict(
                    seconds=perf_counter() - tick,
                    rss=out.f,
                    x=out.x.tolist(),
                    iterations=out.iterations,
                    exitflag=out.exitflag,
                    evaluations=out.evaluations,
                )
            )
        best = min(trials, key=lambda t: t["rss"])
        theta = LOWER + np.array(best["x"]) * (UPPER - LOWER)
        q, _ = simulate(theta, data["forcing"], backend=args.backend)

        def nse(mask):
            """Compute Nash-Sutcliffe efficiency from observed and
            simulated flow.
            """
            y = data["observed"][mask]
            tss = np.sum((y - np.mean(y)) ** 2)
            return (
                float(1 - np.sum((q[mask] - y) ** 2) / tss)
                if tss > 0
                else None
            )

        result = dict(
            gauge=data["gauge"],
            theta=theta.tolist(),
            nse_train=nse(data["train"]),
            nse_eval=nse(data["evaluation"]),
            first_call_seconds=first,
            trials=trials,
        )
        report["basins"].append(result)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(
            f"{data['gauge']}: RSS {best['rss']:.6g}, train NSE"
            f" {result['nse_train']:.4f}",
            flush=True,
        )
