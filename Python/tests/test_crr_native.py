# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_crr_native.py
#
# Basic native dispatch and analytic-gradient checks for all seven
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Basic native dispatch and analytic-gradient checks for all seven
models.
"""

import sys
from pathlib import Path
import unittest
import numpy as np
from build_crr import library_path
from crr_native import simulate, rss_objective

from site_schema import session_defaults


@unittest.skipUnless(library_path().is_file(), "Build CRR library first")
class CRRTests(unittest.TestCase):
    def test_all_models_and_rss_gradient(self):
        """Verify all models and rss gradient."""
        profiles = session_defaults()["model_parameters"]
        forcing = np.column_stack(
            (
                np.array(
                    [2.0, 0.0, 8.0, 1.0, 4.0, 0.0, 5.0, 2.0, 1.0, 0.0]
                ),
                np.full(10, 2.0),
                np.linspace(-2, 8, 10),
            )
        )
        for model in range(1, 8):
            with self.subTest(model=model):
                rows = profiles[str(model)]["parameters"]
                lo = np.array([p["lower"] for p in rows])
                hi = np.array([p["upper"] for p in rows])
                theta = lo + 0.4 * (hi - lo)
                q, J = simulate(model, theta, forcing)
                self.assertEqual(J.shape, (10, len(theta)))
                self.assertTrue(np.all(np.isfinite(q)))
                self.assertTrue(np.all(np.isfinite(J)))
                fun = rss_objective(
                    model,
                    forcing,
                    q + 0.2,
                    np.ones(10, bool),
                    lower=lo,
                    upper=hi,
                )
                f, g = fun(np.full(len(theta), 0.4))
                np.testing.assert_allclose(f, 0.4, rtol=1e-7)
                np.testing.assert_allclose(
                    g,
                    2 * J.T @ np.full(10, -0.2) * (hi - lo),
                    rtol=1e-7,
                    atol=1e-7,
                )

    def test_derivatives_tight_solver(self):
        """Verify derivatives tight solver."""
        profiles = session_defaults()["model_parameters"]
        forcing = np.tile([5.0, 2.0, 4.0], (15, 1))
        for model in range(1, 8):
            with self.subTest(model=model):
                rows = profiles[str(model)]["parameters"]
                theta = np.array(
                    [
                        p["lower"] + 0.4 * (p["upper"] - p["lower"])
                        for p in rows
                    ]
                )
                kwargs = dict(
                    settings=(0.01, 0.2, 1e-7, 1e-7, 1e-7, 2.0, 10000.0)
                )
                q, J = simulate(model, theta, forcing, **kwargs)
                for j in ([0, 12] if model == 6 else [0]):
                    h = 1e-4 * max(1, abs(theta[j]))
                    a = theta.copy()
                    b = theta.copy()
                    a[j] += h
                    b[j] -= h
                    qa, _ = simulate(model, a, forcing, **kwargs)
                    qb, _ = simulate(model, b, forcing, **kwargs)
                    np.testing.assert_allclose(
                        J[:, j], (qa - qb) / (2 * h), rtol=0.03, atol=1e-5
                    )


if __name__ == "__main__":
    unittest.main()
