# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_native_hymod.py
#
# Native wrapper checks; real CAMELS parity is enforced by
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Native wrapper checks; real CAMELS parity is enforced by
benchmark_native.
"""

import unittest
import numpy as np
from build_native import library_path
from site_hymod import simulate


@unittest.skipUnless(library_path().is_file(), "Build native HYMOD first")
class NativeTests(unittest.TestCase):
    def test_states_and_parameter_derivatives(self):
        """Verify states and parameter derivatives."""
        rng = np.random.default_rng(21)
        forcing = np.column_stack(
            (
                rng.uniform(0, 12, 30),
                rng.uniform(0, 4, 30),
                rng.uniform(-8, 12, 30),
            )
        )
        theta = np.array([500.0, 1.5, 0.6, 0.05, 1.0, 0.0, 3.0])
        settings = (0.01, 0.2, 1e-6, 1e-7, 1e-7, 2.0, 10000.0)
        q, J = simulate(theta, forcing, backend="cpp", settings=settings)
        nq, nJ = simulate(theta, forcing, settings=settings)
        np.testing.assert_allclose(q, nq, rtol=1e-6, atol=1e-7)
        np.testing.assert_allclose(J, nJ, rtol=1e-5, atol=1e-7)
        for j in range(7):
            h = 1e-4 * max(1, abs(theta[j]))
            plus = theta.copy()
            minus = theta.copy()
            plus[j] += h
            minus[j] -= h
            qp, _ = simulate(
                plus, forcing, backend="cpp", settings=settings
            )
            qm, _ = simulate(
                minus, forcing, backend="cpp", settings=settings
            )
            np.testing.assert_allclose(
                J[:, j], (qp - qm) / (2 * h), rtol=0.02, atol=1e-4
            )

    def test_noncontiguous_forcing(self):
        """Verify noncontiguous forcing."""
        forcing = np.tile([5.0, 2.0, 4.0], (20, 1))
        theta = np.array([500.0, 1.5, 0.6, 0.05, 1.0, 0.0, 3.0])
        q, J = simulate(theta, forcing[::2], backend="cpp")
        expected, eJ = simulate(theta, forcing[::2].copy(), backend="cpp")
        np.testing.assert_array_equal(q, expected)
        np.testing.assert_array_equal(J, eJ)

    def test_invalid_input(self):
        """Verify invalid input."""
        with self.assertRaises(ValueError):
            simulate(np.zeros(7), np.zeros((3, 3)), backend="cpp")
        with self.assertRaises(ValueError):
            simulate(
                np.array([500.0, 1.5, 0.6, 0.05, 1.0, 0.0, 3.0]),
                np.zeros((0, 3)),
                backend="cpp",
            )


if __name__ == "__main__":
    unittest.main()
