# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/cart/test_cart_ode.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Independent conservation, event, tangent and time-resolution checks.
"""

import unittest

import numpy as np

from cart_ode import routing_fractions, simulate
from cart_python import simulate as split


P = np.array([np.log(225), 1.26, np.log(33), np.log(78), -.35,
              0., 0., np.log(3), 1., np.log(6.5), np.log(.001)])
Z = np.array([100., 220., 350., 510., 740.])
A = np.array([.1, .2, .3, .25, .15])


def forcing():
    """
    Partition band precipitation and snowmelt into liquid supply.

    SYNOPSIS:
        forcing(...)

    INPUT ARGUMENTS:
    None.

    OUTPUT ARGUMENTS:
    np.column_stack((2 + 8 * (np.sin(2 * t) > 0.6), 4 - 14 * np.cos(t / 2),
    0.7 + 0.5 * np.cos(t), t < 6)) : computed result in the operand units.
    """
    t = np.arange(8)
    return np.column_stack((2+8*(np.sin(2*t) > .6), 4-14*np.cos(t/2),
                            .7+.5*np.cos(t), t < 6))


class ContinuousCartTests(unittest.TestCase):
    def test_all_parameter_tangents_with_snow_events(self):
        """
        Verify all parameter tangents with snow events.

        SYNOPSIS:
            test_all_parameter_tangents_with_snow_events(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        f = forcing()
        for detached in (False, True):
            options = dict(refinement=4, rtol=2e-12, atol=2e-14)
            if detached:
                options['initial'] = simulate(f, Z, A, P,
                    return_diagnostics=True,
                                               **options)[2]['initial']
            q, jac, diag = simulate(f, Z, A, P, jacobian=True,
                                   return_diagnostics=True, **options)
            self.assertGreaterEqual(diag['snow_events'], 5)
            self.assertLess(np.max(np.abs(diag['balance_residual_mm'])),
                1e-8)
            self.assertGreaterEqual(diag['states'].min(), -1e-9)
            for k in range(11):
                h = 1e-4
                direction = np.eye(11)[k]*h
                plus = simulate(f, Z, A, P+direction, **options)[0]
                minus = simulate(f, Z, A, P-direction, **options)[0]
                fd = (plus-minus)/(2*h)
                np.testing.assert_allclose(jac[:, k], fd, atol=1e-6, rtol=
                    2e-5,
                                           err_msg=
    f'coordinate {k}, detached={detached}')

    def test_constant_forcing_subdivision_and_tolerances(self):
        """
        Verify constant forcing subdivision and tolerances.

        SYNOPSIS:
            test_constant_forcing_subdivision_and_tolerances(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        f = forcing()
        q, j = simulate(f, Z, A, P, jacobian=True, refinement=4,
                        rtol=1e-10, atol=1e-12)
        for resolution in (24, 96):
            hourly = np.repeat(f, resolution, axis=0)
            hourly[:, [0, 2]] /= resolution
            fine, tangent = simulate(hourly, Z, A, P, jacobian=True,
                                     refinement=4, resolution=resolution,
                                     rtol=1e-10, atol=1e-12)
            np.testing.assert_allclose(fine.reshape(-1, resolution).sum(1)
                , q,
                                       atol=3e-8, rtol=3e-8)
            np.testing.assert_allclose(tangent.reshape(-1, resolution, 11)
                .sum(1), j,
                                       atol=3e-7, rtol=3e-6)
        coarse, _ = simulate(f, Z, A, P, refinement=4, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(coarse, q, atol=2e-6, rtol=2e-6)
        # Irregular observation times preserve the identical continuous
        # forcing.
        durations = np.tile([.2, .3, .5], len(f))
        irregular = np.repeat(f, 3, axis=0)
        irregular[:, [0, 2]] *= durations[:, None]
        qi, ji = simulate(irregular, Z, A, P, refinement=4, jacobian=True,
                          interval_days=durations, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(qi.reshape(-1, 3).sum(1), q, atol=3e-8,
            rtol=3e-8)
        np.testing.assert_allclose(ji.reshape(-1, 3, 11).sum(1), j, atol=
            3e-7, rtol=3e-6)

    def test_cascade_impulse_moments(self):
        # Residence time in stage j is Erlang(j, m/X4). Its mixture
        # moments independently check cascade exit fractions and refinement.
        """
        Verify cascade impulse moments.

        SYNOPSIS:
            test_cascade_impulse_moments(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        errors = []
        for m in (4, 16, 64):
            for branch, eta in enumerate(routing_fractions(m)):
                weights = eta*np.r_[1., np.cumprod(1-eta[:-1])]
                self.assertAlmostEqual(weights.sum(), 1., places=13)
                self.assertEqual(eta[-1], 1.)
                lag = np.arange(1, len(eta)+1)/m
                mean = weights@lag
                target = 1 if branch else 5/7
                if branch == 0:
                    errors.append(abs(mean-target))
        self.assertTrue(errors[2] < errors[1] < errors[0])

    def test_floor_dry_and_nondifferentiable_flags(self):
        """
        Verify floor dry and nondifferentiable flags.

        SYNOPSIS:
            test_floor_dry_and_nondifferentiable_flags(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        f = np.array([[0., 15., 0., 1.], [0., 15., 0., 1.]])
        initial = np.zeros(8+3*2)
        q, j, diag = simulate(f, Z, A, P, initial=initial, refinement=2,
                              jacobian=True, return_diagnostics=True)
        np.testing.assert_array_equal(q, 0.)
        np.testing.assert_array_equal(j, 0.)
        self.assertEqual(diag['gbar'], 1e-6)
        np.testing.assert_array_equal(diag['gbar_gradient'], 0.)
        f[:, 1] = 0.
        _, _, diag = simulate(f, np.zeros(5), A, P, refinement=2,
                              return_diagnostics=True)
        self.assertEqual(diag['threshold_coincidences'], [0, 1])

    def test_simultaneous_depletion_and_boundary_flag(self):
        """
        Verify simultaneous depletion and boundary flag.

        SYNOPSIS:
            test_simultaneous_depletion_and_boundary_flag(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        f = np.array([[2., 8., 1., 1.], [1., 10., 1., 1.]])
        initial = np.zeros(8+3*4)
        initial[:5] = 4.2
        initial[5:8] = [150., 12., 5.]
        options = dict(initial=initial, refinement=4, rtol=2e-12, atol=
            2e-14)
        q, j, d = simulate(f, np.zeros(5), A, P, jacobian=True,
                           return_diagnostics=True, **options)
        self.assertEqual(d['snow_events'], 5)
        direction = np.eye(11)[7]*1e-4
        fd = (simulate(f, np.zeros(5), A, P+direction, **options)[0]
              -simulate(f, np.zeros(5), A, P-direction, **options)[0])/2e-4
        np.testing.assert_allclose(j[:, 7], fd, atol=1e-6, rtol=2e-5)
        # Exactly exhausted at the forcing boundary: report ambiguity.
        initial[:5] = 3.*8.+.001
        _, _, d = simulate(f, np.zeros(5), A, P, return_diagnostics=True,
            **options)
        self.assertEqual(len(d['boundary_event_coincidences']), 5)

    def test_split_approaches_continuous_solution(self):
        # Rain-only data isolate finite-update error from snow events.
        # Many bins preserve the author's two hydrographs as dt -> 0;
        # finite m still contributes a distinct routing approximation.
        """
        Verify split approaches continuous solution.

        SYNOPSIS:
            test_split_approaches_continuous_solution(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        f = forcing()
        f[:, 1] = 20.
        initial = np.zeros(8+3*64)
        initial[5:8] = [150., 12., 5.]
        continuous = simulate(f, Z, A, P, refinement=64, initial=initial,
                              rtol=1e-10, atol=1e-12)[0]
        errors = []
        for res in (1, 24, 96):
            data = np.repeat(f, res, axis=0)
            data[:, [0, 2]] /= res
            interval = split(data, Z, A, P, initial=initial[5:8],
                             resolution=res, mcode=4)[0]
            errors.append(np.linalg.norm(interval.reshape(-1, res).sum(1)-
                continuous))
        self.assertLess(errors[-1], errors[0])

    def test_invalid_input(self):
        """
        Verify invalid input.

        SYNOPSIS:
            test_invalid_input(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        for options in ({'refinement': 0}, {'refinement': 1.5},
                        {'initial': np.zeros(3)}, {'resolution': 2},
                        {'interval_days': 0}, {'interval_days': [1., 2.]}):
            with self.assertRaises(ValueError):
                simulate(forcing(), Z, A, P, **options)


if __name__ == '__main__':
    unittest.main()
