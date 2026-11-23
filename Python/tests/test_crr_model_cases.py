# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_crr_model_cases.py
#
# Cross-language discharge/Jacobian regressions generated with MATLAB
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Cross-language discharge/Jacobian regressions generated with MATLAB
R2026a.
"""

from pathlib import Path
import json, unittest
import numpy as np
from crr_models import simulate


class MatlabModelCases(unittest.TestCase):
    def test_seven_models_four_cases_three_resolutions(self):
        """Verify seven models four cases three resolutions."""
        path = Path(__file__).parent / "tests/data/matlab_model_cases.npz"
        with np.load(path, allow_pickle=False) as fixture:
            manifest = json.loads(str(fixture["manifest"]))
            self.assertEqual(len(manifest), 84)
            for row in manifest:
                with self.subTest(**row):
                    key = row["key"]
                    kwargs = dict(
                        mcode=row["mcode"],
                        y0=fixture[key + "_y0"],
                        model_config=json.loads(
                            str(fixture[key + "_config"])
                        ),
                    )
                    if row["model"] == 6:
                        kwargs["route_hbv"] = False
                    q, j = simulate(
                        row["model"],
                        fixture[key + "_theta"],
                        fixture[key + "_forcing"],
                        **kwargs,
                    )
                    np.testing.assert_allclose(
                        q, fixture[key + "_q"], atol=1e-9, rtol=1e-10
                    )
                    np.testing.assert_allclose(
                        j, fixture[key + "_J"], atol=1e-9, rtol=1e-10
                    )
                    if row["mcode"] == 2:
                        self.assertGreater(abs(q[-1]), 0)

    def test_single_interval_ode45_endpoints(self):
        """Verify single interval ode45 endpoints."""
        with np.load(
            Path(__file__).parent / "tests/data/matlab_model_cases.npz",
            allow_pickle=False,
        ) as fixture:
            for row in json.loads(str(fixture["endpoints"])):
                with self.subTest(model=row["model"]):
                    key = row["key"]
                    kwargs = dict(
                        mcode=2,
                        y0=fixture[key + "_y0"],
                        model_config=json.loads(
                            str(fixture[key + "_config"])
                        ),
                    )
                    if row["model"] == 6:
                        kwargs["route_hbv"] = False
                    q, j = simulate(
                        row["model"],
                        fixture[key + "_theta"],
                        fixture[key + "_forcing"],
                        **kwargs,
                    )
                    np.testing.assert_allclose(
                        q, fixture[key + "_q"], atol=1e-9, rtol=1e-10
                    )
                    np.testing.assert_allclose(
                        j, fixture[key + "_J"], atol=1e-9, rtol=1e-10
                    )

    def test_default_remains_existing_native_cpp(self):
        """Verify default remains existing native cpp."""
        from crr_native import simulate as native

        with np.load(
            Path(__file__).parent / "tests/data/matlab_model_cases.npz",
            allow_pickle=False,
        ) as fixture:
            for model in range(1, 8):
                key = f"m{model}_c4_r1"
                theta = fixture[key + "_theta"]
                if model == 6:
                    theta = np.r_[theta, 24.0]
                kwargs = dict(
                    y0=fixture[key + "_y0"],
                    model_config=json.loads(str(fixture[key + "_config"])),
                )
                expected = native(
                    model, theta, fixture[key + "_forcing"], **kwargs
                )
                actual = simulate(
                    model, theta, fixture[key + "_forcing"], **kwargs
                )
                for a, b in zip(actual, expected):
                    np.testing.assert_array_equal(a, b)


if __name__ == "__main__":
    unittest.main()
