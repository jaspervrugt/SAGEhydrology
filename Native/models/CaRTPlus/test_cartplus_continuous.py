# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/test_cartplus_continuous.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Conservation, sensitivities and observation-interval independence.
"""

import json
from pathlib import Path
import unittest

import numpy as np

from cartplus_continuous import simulate
from cartplus_inputs import PARAMETER_MAP, physical_parameters


def parameters():
    """
    Load the reference physical parameter vector.

    SYNOPSIS:
        parameters(...)

    INPUT ARGUMENTS:
    None.

    OUTPUT ARGUMENTS:
    physical_parameters(*({k: v[2] for k, v in spec[g].items()} for g in
    ('static', 'law')), spec['shared']) : computed result in the operand
    units.
    """
    spec = json.loads(Path(__file__).with_name('parameters.json').
        read_text())
    return physical_parameters(*({k: v[2] for k, v in spec[g].items()}
                                  for g in ('static', 'law')), spec[
                                      'shared'])


def forcing():
    """
    Partition band precipitation and snowmelt into liquid supply.

    SYNOPSIS:
        forcing(...)

    INPUT ARGUMENTS:
    None.

    OUTPUT ARGUMENTS:
    np.array([[5.0, -5.0, 1.0, 1.0], [12.0, -2.0, 0.8, 1.0], [8.0, 12.0,
    1.1, 1.0], [2.0, 16.0, 1.4, 0.0]]) : computed result in the operand
    units.
    """
    return np.array([[5., -5., 1., 1.], [12., -2., .8, 1.],
                     [8., 12., 1.1, 1.], [2., 16., 1.4, 0.]])


Z = np.array([100., 220., 350., 510., 740.])
A = np.array([.1, .2, .3, .25, .15])


def perturb(p, direction, step):
    """
    Apply a parameter perturbation for a derivative check.

    SYNOPSIS:
        perturb(...)

    INPUT ARGUMENTS:
    p : Model parameter vector; CaRT calibration or CaRT+ physical units.
    direction : Parameter derivative direction or perturbation index.
    step : Step used by this numerical operator.

    OUTPUT ARGUMENTS:
    out : computed result in the operand units.
    """
    out = p.copy()
    for k, (_, _, mult) in enumerate(PARAMETER_MAP):
        out[k] = p[k]*np.exp(step*direction[k]) if mult is not None else p \
            [k]+step*direction[k]
    return out


class ContinuousPlusTests(unittest.TestCase):
    def test_tangent_direction_initialization_and_conservation(self):
        """
        Verify tangent direction initialization and conservation.

        SYNOPSIS:
            test_tangent_direction_initialization_and_conservation(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        p, f = parameters(), forcing()
        direction = np.sin(np.arange(42)+.2)/5
        opts = dict(refinement=2, rtol=2e-11, atol=2e-13,
                    sensitivity_coordinates='calibration')
        for detached in (False, True):
            if detached:
                opts['initial'] = diagnostics['initial']
            q, j, diagnostics = simulate(f, Z, A, p, jacobian=True,
                                         return_diagnostics=True, **opts)
            fd = (simulate(f, Z, A, perturb(p, direction, 1e-3), **opts)[0]
                  -simulate(f, Z, A, perturb(p, direction, -1e-3), **opts)
                      [0])/2e-3
            np.testing.assert_allclose(j@direction, fd, atol=2e-5, rtol=
                2e-4)
            self.assertEqual(j.shape, (4, 42))
            self.assertEqual(diagnostics['state_jacobian'].shape, (4, 20,
                42))
            self.assertLess(np.max(np.abs(diagnostics[
                'balance_residual_mm'])), 1e-8)
            np.testing.assert_allclose(q+diagnostics['channel_loss_mm'],
                                       diagnostics['preloss_depth_mm'],
                                        atol=1e-12)
            self.assertGreaterEqual(diagnostics['states'].min(), -1e-8)

    def test_physical_and_calibration_chain_rule(self):
        """
        Verify physical and calibration chain rule.

        SYNOPSIS:
            test_physical_and_calibration_chain_rule(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        p, f = parameters(), forcing()[:2]
        q, jp = simulate(f, Z, A, p, jacobian=True, refinement=2,
                          rtol=1e-10, atol=1e-12)
        qc, jc = simulate(f, Z, A, p, jacobian=True, refinement=2,
                          rtol=1e-10, atol=1e-12, sensitivity_coordinates=
                              'calibration')
        scale = np.array([v if mult is not None else 1. for v, (_, _, mult
            ) in zip(p, PARAMETER_MAP)])
        np.testing.assert_allclose(qc, q, atol=1e-7, rtol=1e-7)
        np.testing.assert_allclose(jc, jp*scale, atol=1e-6, rtol=1e-5)

    def test_forcing_subdivision_and_irregular_intervals(self):
        """
        Verify forcing subdivision and irregular intervals.

        SYNOPSIS:
            test_forcing_subdivision_and_irregular_intervals(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        p, f = parameters(), forcing()
        opts = dict(refinement=2, rtol=2e-11, atol=2e-13,
                    sensitivity_coordinates='calibration')
        q, j = simulate(f, Z, A, p, jacobian=True, **opts)
        for res in (24, 96):
            sub = np.repeat(f, res, axis=0)
            sub[:, [0, 2]] /= res
            qs, js = simulate(sub, Z, A, p, jacobian=True, resolution=res,
                **opts)
            np.testing.assert_allclose(qs.reshape(-1, res).sum(1), q, atol
                =2e-7, rtol=2e-7)
            np.testing.assert_allclose(js.reshape(-1, res, 42).sum(1), j,
                atol=3e-6, rtol=3e-5)
        durations = np.tile([.2, .3, .5], len(f))
        sub = np.repeat(f, 3, axis=0)
        sub[:, [0, 2]] *= durations[:, None]
        qs, js = simulate(sub, Z, A, p, jacobian=True, interval_days=
            durations, **opts)
        np.testing.assert_allclose(qs.reshape(-1, 3).sum(1), q, atol=2e-7,
            rtol=2e-7)
        np.testing.assert_allclose(js.reshape(-1, 3, 42).sum(1), j, atol=
            3e-6, rtol=3e-5)

    def test_floor_and_phase_flags(self):
        """
        Verify floor and phase flags.

        SYNOPSIS:
            test_floor_and_phase_flags(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        p = parameters()
        f = np.array([[0., 15., 0., 1.]])
        _, _, d = simulate(f, Z, A, p, initial=np.zeros(20), refinement=2,
                            jacobian=True, return_diagnostics=True)
        self.assertAlmostEqual(d['reference_storage'], 1e-6*p[19])
        np.testing.assert_allclose(d['reference_gradient'][np.arange(42)
            != 19], 0.)
        f[0, 1] = p[0]
        _, _, d = simulate(f, np.zeros(5), A, p, refinement=2,
            return_diagnostics=True)
        self.assertEqual(d['threshold_coincidences'], [0])

    def test_active_groundwater_uptake_drydown(self):
        """
        Verify active groundwater uptake drydown.

        SYNOPSIS:
            test_active_groundwater_uptake_drydown(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        p = parameters()
        f = np.column_stack((np.zeros(6), np.full(6, 15.),
                              np.full(6, 4.), np.ones(6)))
        initial = np.zeros(20)
        initial[[6, 7, 8, 9, 10, 11, 12, 13]] = [150., 4.5, 15., .1, 1.,
            .5, 5., 16.5]
        opts = dict(initial=initial, refinement=2, rtol=2e-12, atol=2e-14,
                    sensitivity_coordinates='calibration')
        q, j = simulate(f, Z, A, p, jacobian=True, **opts)
        self.assertGreater(np.linalg.norm(j[:, 20]), 1e-5)
        for k in (20, 21, 24, 31):
            direction = np.eye(42)[k]
            fd = (simulate(f, Z, A, perturb(p, direction, 1e-3), **opts)[0]
                  -simulate(f, Z, A, perturb(p, direction, -1e-3), **opts)
                      [0])/2e-3
            np.testing.assert_allclose(j[:, k], fd, atol=1e-6, rtol=2e-4)

    def test_invalid_inputs(self):
        """
        Verify invalid inputs.

        SYNOPSIS:
            test_invalid_inputs(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        for options in ({'refinement': 0}, {'availability_days': 0},
                        {'transition_fraction': 2}, {'interval_days': -1},
                        {'initial': np.zeros(3)}, {
                            'sensitivity_coordinates': 'unknown'}):
            with self.assertRaises(ValueError):
                simulate(forcing(), Z, A, parameters(), **options)


if __name__ == '__main__':
    unittest.main()
