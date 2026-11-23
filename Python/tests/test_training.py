# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_training.py
#
# Portable headless SAGE entry point.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Verify standalone training, parallel results and checkpoint resume."""

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
from synthetic_training import configuration
from run_sage import execute
from sage_parameter_scaling import physical, normalized, derivative


class TrainingTests(unittest.TestCase):
    def test_serial_parallel_and_resume(self):
        """Match complete, resumed and parallel runs on artificial basins."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = configuration(root / "data")
            full = execute(config, root / "full")
            short_config = dict(config, max_iter=1)
            execute(short_config, root / "short")
            resumed_config = dict(
                config,
                resume_checkpoint=str(
                    root / "short/training_checkpoint.npz"
                ),
            )
            resumed = execute(resumed_config, root / "resumed")
            parallel = execute(
                dict(config, parallel=True), root / "parallel"
            )
            for other in (resumed, parallel):
                np.testing.assert_allclose(
                    [h["objective"] for h in full["history"]],
                    [h["objective"] for h in other["history"]],
                    rtol=0,
                    atol=1e-12,
                )
                np.testing.assert_allclose(
                    full["parameter_field"]["normalized"],
                    other["parameter_field"]["normalized"],
                    rtol=0,
                    atol=1e-12,
                )
            self.assertEqual(len(full["basins"]), 4)
            self.assertEqual(len(full["training_basins"]), 3)
            self.assertLess(
                full["history"][-1]["objective"],
                full["history"][0]["objective"],
            )
            self.assertEqual(full["attribute_scaler"]["mode"], "matlab")
            self.assertFalse(
                any(n.startswith("PySide") for n in sys.modules)
            )
            self.assertFalse(any(n == "private_crr" for n in sys.modules))

    def test_log_parameter_chain_rule(self):
        """Verify round trips and derivatives for linear/log mappings."""
        c = dict(
            lower=[0.01, -3],
            upper=[100, 3],
            parameter_scales=["log", "linear"],
        )
        x = np.array([0.4, 0.7])
        theta = physical(c, x)
        np.testing.assert_allclose(normalized(c, theta), x)
        h = 1e-6
        np.testing.assert_allclose(
            derivative(c, theta),
            (physical(c, x + h) - physical(c, x - h)) / (2 * h),
            rtol=1e-8,
        )

    def test_private_model_is_rejected_before_execution(self):
        """Prevent a public source run from invoking unavailable models."""
        with self.assertRaisesRegex(ValueError, "public models"):
            execute({"model": 16}, "unused-private-model-output")
