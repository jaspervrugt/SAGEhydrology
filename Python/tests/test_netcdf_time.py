# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_netcdf_time.py
#
# Compare fast dates and fallback calendars against netCDF4.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Compare fast dates and fallback calendars against netCDF4."""

import unittest
import numpy as np
import pandas as pd
from netCDF4 import num2date
from netcdf_time import decode


class NetCDFTimeTests(unittest.TestCase):
    def test_integral_modern_clocks(self):
        """Verify integral modern clocks."""
        for unit, origin in (
            ("hours", "1980-01-01 09:00:00"),
            ("minutes", "1950-01-01"),
            ("days", "1980-01-01 12:00:00"),
            ("seconds", "2000-02-28 23:59:59"),
        ):
            values = np.array([-1, 0, 1, 24, 365, 100000], dtype=float)
            units = unit + " since " + origin
            expected = pd.DatetimeIndex(
                [str(v) for v in num2date(values, units)]
            )
            np.testing.assert_array_equal(decode(values, units), expected)

    def test_fractional_and_nonstandard_calendars(self):
        """Verify fractional and nonstandard calendars."""
        for values, units, calendar in (
            ([0.125, 0.5, 1.75], "seconds since 2000-01-01", "standard"),
            ([0, 1, 2, 3], "days since 2000-02-28", "noleap"),
            ([0, 1, 2], "days since 2000-01-01", "360_day"),
        ):
            expected = pd.DatetimeIndex(
                [str(v) for v in num2date(values, units, calendar=calendar)]
            )
            np.testing.assert_array_equal(
                decode(values, units, calendar), expected
            )


if __name__ == "__main__":
    unittest.main()
