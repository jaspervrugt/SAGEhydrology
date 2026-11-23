# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_attribute_preprocessing.py
#
# Attribute preprocessing regression checks.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

import json
import unittest
import numpy as np
from sage_attribute_preprocessing import METHODS, fit, apply, settings
from sage_network import standardize

class AttributePreprocessingTests(unittest.TestCase):

    def test_reporting_uses_effective_methods(self):
        """Test reporting uses effective methods."""
        from sage_attribute_preprocessing import configuration_summary
        c = dict(attribute_preprocessing='custom', attribute_scaling='robust', attribute_fit_cohort='training', attribute_missing='median', attributes=['area', 'slope'], attribute_methods_enabled=True, attribute_methods={'area': 'log_zscore', 'not_selected': 'none'})
        summary = dict(configuration_summary(c))
        self.assertEqual(summary['Global method'], '4 Robust scaling: median / IQR')
        self.assertTrue(summary['area'].startswith('5 '))
        self.assertNotIn('not_selected', summary)
        c['attribute_methods_enabled'] = False
        self.assertNotIn('area', dict(configuration_summary(c)))
        c['attribute_preprocessing'] = 'matlab'
        summary = dict(configuration_summary(c))
        self.assertTrue(summary['Global method'].startswith('1 '))
        self.assertEqual(summary['Std deviation'], 'Sample (N - 1)')

    def options(self, method='zscore', **kw):
        """Options."""
        return dict(dict(method=method, cohort='training', missing='median', attributes=['a', 'b'], overrides={}), **kw)

    def test_all_transforms_are_serializable_and_reusable(self):
        """Test all transforms are serializable and reusable."""
        raw = np.array([[0.0, 1.0], [1.0, 3.0], [3.0, 2.0], [5.0, 8.0], [100.0, 1000.0]])
        for method in METHODS:
            with self.subTest(method=method):
                out, state = fit(raw, [0, 1, 2, 3], self.options(method))
                self.assertTrue(np.isfinite(out).all())
                loaded = json.loads(json.dumps(state, allow_nan=False))
                np.testing.assert_allclose(out, apply(raw, loaded), rtol=0, atol=0)
                changed = raw.copy()
                changed[-1] *= 2
                _, other = fit(changed, [0, 1, 2, 3], self.options(method))
                self.assertEqual(state, other)

    def test_known_scaling_values_and_extrapolation(self):
        """Test known scaling values and extrapolation."""
        x = np.array([[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]])
        for method, target in (('minmax01', [0, 1, 2]), ('minmax11', [-1, 1, 3])):
            out, _ = fit(x, [0, 1], self.options(method))
            np.testing.assert_allclose(out[:, 0], target)
        out, _ = fit(x, [0, 1], self.options())
        np.testing.assert_allclose(out[:2].mean(0), 0, atol=1e-15)
        np.testing.assert_allclose(out[:2].std(0, ddof=1), 1)

    def test_missing_values_and_per_attribute_overrides(self):
        """Test missing values and per attribute overrides."""
        x = np.array([[0.0, 2.0], [2.0, 4.0], [np.nan, 10.0]])
        out, state = fit(x, [0, 1], self.options('minmax01', overrides={'b': 'none'}))
        self.assertEqual(out[2, 0], 0.5)
        np.testing.assert_array_equal(out[:, 1], x[:, 1])
        with self.assertRaises(ValueError):
            fit(x, [0, 1], self.options(missing='reject'))
        with self.assertRaises(ValueError):
            apply([[np.nan, 1]], dict(state, missing='reject'))
        with self.assertRaises(ValueError):
            fit([[-1.0, 2.0], [2.0, 4.0]], [0, 1], self.options('log_zscore'))

    def test_legacy_presets_unchanged(self):
        """Test legacy presets unchanged."""
        x = np.array([[1.0, 2.0], [3.0, 4.0], [100.0, 200.0]])
        for mode in ('matlab', 'training'):
            out, state = standardize(x, [0, 1], mode)
            np.testing.assert_allclose(out, apply(x, state))
        _, matlab = standardize(x, [0, 1], 'matlab')
        np.testing.assert_allclose(matlab['scale'], x.std(0, ddof=1))
        _, legacy = standardize(x, [0, 1], 'training')
        np.testing.assert_allclose(legacy['scale'], x[:2].std(0))

    def test_per_attribute_checkbox_controls_effective_methods(self):
        """Test per attribute checkbox controls effective methods."""
        c = dict(attribute_preprocessing='custom', attribute_scaling='zscore', attribute_fit_cohort='training', attribute_missing='reject', attribute_methods={'a': 'none'}, attribute_methods_enabled=True)
        self.assertEqual(settings(c)['overrides'], {'a': 'none'})
        c['attribute_methods_enabled'] = False
        self.assertEqual(settings(c)['overrides'], {})
        self.assertEqual(c['attribute_methods'], {'a': 'none'})
        c.pop('attribute_methods_enabled')
        self.assertEqual(settings(c)['overrides'], {'a': 'none'})

    def test_custom_nonlinear_checkpoint_resume(self):
        """Test custom nonlinear checkpoint resume."""
        from pathlib import Path
        import tempfile
        from synthetic_training import configuration
        from run_sage import execute
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            c = configuration(root / 'data')
            c.update(attribute_preprocessing='custom', attribute_scaling='quantile_normal', attribute_fit_cohort='training', attribute_missing='reject')
            full = execute(c, root / 'full')
            execute(dict(c, max_iter=1), root / 'short')
            resumed = execute(dict(c, resume_checkpoint=str(root / 'short/training_checkpoint.npz')), root / 'resumed')
            np.testing.assert_allclose([h['objective'] for h in full['history']], [h['objective'] for h in resumed['history']], rtol=0, atol=1e-12)
            np.testing.assert_allclose(full['parameter_field']['normalized'], resumed['parameter_field']['normalized'], rtol=0, atol=1e-12)

    def test_constant_columns(self):
        """Test constant columns."""
        for method in METHODS:
            out, state = fit(np.ones((4, 2)), [0, 1, 2], self.options(method))
            self.assertTrue(np.isfinite(out).all(), method)
            np.testing.assert_allclose(apply(np.ones((4, 2)), state), out)
if __name__ == '__main__':
    unittest.main()
