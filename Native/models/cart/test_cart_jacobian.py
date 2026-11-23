# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/cart/test_cart_jacobian.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Check the CaRT daily implementations and detached-store derivatives.
"""

import unittest

import numpy as np

from cart_python import simulate


class CartJacobianTests(unittest.TestCase):
    def test_daily_backends_and_centered_differences(self):
        """
        Verify daily backends and centered differences.

        SYNOPSIS:
            test_daily_backends_and_centered_differences(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        rng = np.random.default_rng(52)
        forcing = np.column_stack([
            rng.gamma(2, 3, 36), rng.uniform(-6, 15, 36),
            rng.uniform(.1, 4, 36), np.arange(36) < 24,
        ])
        elevation = np.arange(5) * 200
        area = np.full(5, .2)
        parameters = np.array([
            np.log(225), 1.26, np.log(33), np.log(78), -.35, 0., 0.,
            np.log(3), 1., np.log(6.5), np.log(.001),
        ])
        for resolution in (1, 24, 96):
            data = forcing.copy()
            data[:, [0, 2]] /= resolution
            q, jac = simulate(
                data, elevation, area, parameters, jacobian=True,
                resolution=resolution, mcode=4,
            )
            for case in (1, 2, 3, 4):
                with self.subTest(case=case, resolution=resolution):
                    output, derivative, stores = simulate(
                        data, elevation, area, parameters, jacobian=True,
                        mcode=case, resolution=resolution,
                        return_initial=True,
                    )
                    if case == 2 and resolution != 24:
                        np.testing.assert_allclose(output, q, atol=1e-8)
                        np.testing.assert_allclose(
                            derivative, jac, atol=1e-8
                        )
                    for index, delta in enumerate(np.eye(11) * 1e-5):
                        kwargs = dict(initial=stores, mcode=case,
                                      resolution=resolution)
                        plus, _ = simulate(
                            data, elevation, area, parameters+delta,
                            **kwargs
                        )
                        minus, _ = simulate(
                            data, elevation, area, parameters-delta,
                            **kwargs
                        )
                        np.testing.assert_allclose(
                            derivative[:, index], (plus-minus)/2e-5,
                            atol=1e-7, rtol=1e-6,
                        )


if __name__ == "__main__":
    unittest.main()
