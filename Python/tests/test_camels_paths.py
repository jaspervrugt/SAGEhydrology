# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_camels_paths.py
#
# Official CAMELS-US Daymet/CIDA filename and installer regressions.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Official CAMELS-US Daymet/CIDA filename and installer regressions."""

from pathlib import Path
import unittest
import numpy as np
import test_regional_data as helpers
from camels_us import forcing_files, read_daily
from dataset_install import _destination


class CamelsPathTests(unittest.TestCase):
    setUp = helpers.RegionalTests.setUp
    tearDown = helpers.RegionalTests.tearDown

    def files(self, label="cida"):
        """Create or resolve the test input files."""
        folder = self.root / "daily/v1p2/forcing/daymet/01"
        folder.mkdir(parents=True, exist_ok=True)
        forcing = folder / f"01022500_lump_{label}_forcing_leap.txt"
        forcing.write_text(
            "40\n100\n100000000\nYear Mnth Day Hr dayl(s) "
            "prcp(mm/day) srad(W/m2) swe(mm) tmax(C) tmin(C) vp(Pa)\n"
            + "".join(
                f"2005 10 {day} 0 40000 2 100 0 10 0 800\n"
                for day in range(1, 5)
            )
        )
        folder = self.root / "daily/v1p2/streamflow/01"
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "01022500_streamflow_qc.txt").write_text(
            "".join(f"01022500 2005 10 {day} 1 A\n" for day in range(1, 5))
        )
        return forcing

    def test_official_cida_files_load_as_daymet(self):
        """Verify official cida files load as daymet."""
        file = self.files()
        data = read_daily(
            self.root, "01022500", product="daymet", **self.period
        )
        self.assertEqual(data["forcing_file"], str(file))
        np.testing.assert_array_equal(data["forcing"][:, 0], 2)
        self.assertEqual(data["train"].sum(), 2)

    def test_noncanonical_alias_is_not_a_matlab_source(self):
        """Require the exact CIDA source pattern used by MATLAB."""
        self.files("daymet")
        with self.assertRaisesRegex(FileNotFoundError, "found 0"):
            read_daily(
                self.root, "01022500", product="daymet", **self.period
            )

    def test_exact_cida_source_ignores_noncanonical_alias(self):
        """Select the MATLAB CIDA file even if an alias is present."""
        expected = self.files("cida")
        self.files("daymet")
        data = read_daily(
            self.root, "01022500", product="daymet", **self.period
        )
        self.assertEqual(data["forcing_file"], str(expected))
        np.testing.assert_array_equal(data["forcing"][:, 0], 2)

    def test_installer_routes_cida_to_daymet_and_catalogue_agrees(self):
        """Verify installer routes cida to daymet and catalogue agrees."""
        file = self.files()
        self.assertEqual(
            _destination(file, "CAMELS_US", helpers.SAGE),
            Path("daily/v1p2/forcing/daymet/01") / file.name,
        )
        from regional_data import catalogue

        self.assertEqual(
            catalogue(
                self.root, "CAMELS_US", helpers.SAGE, product="daymet"
            ),
            [("01022500", "", "")],
        )
        self.assertEqual(forcing_files(self.root, "nldas"), [])


if __name__ == "__main__":
    unittest.main()
