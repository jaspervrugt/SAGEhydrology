# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_netcdf_buffered.py
#
# Buffered I/O must preserve masks and retain native thread safety.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Buffered I/O must preserve masks and retain native thread safety."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import numpy as np
from netCDF4 import Dataset
from netcdf_io import dataset


class BufferedNetCDFTests(unittest.TestCase):
    def test_masks_values_and_parallel_reads(self):
        """Verify masks values and parallel reads."""
        with TemporaryDirectory() as folder:
            path = Path(folder) / "fixture.nc"
            with Dataset(path, "w") as handle:
                handle.createDimension("time", 4)
                time = handle.createVariable("time", "i4", ("time",))
                time.units = "hours since 1980-01-01 09:00:00"
                time[:] = np.arange(4)
                q = handle.createVariable(
                    "q_obs", "f4", ("time",), fill_value=-9999
                )
                q[:] = [1.25, -9999, 3.5, 0]

            def read(buffered):
                """Read the configured basin data and reference time
                window.
                """
                with dataset(path, buffered=buffered) as handle:
                    return (
                        handle.variables["time"].units,
                        np.ma.filled(handle.variables["q_obs"][:], np.nan),
                    )

            units, expected = read(False)
            with ThreadPoolExecutor(max_workers=8) as pool:
                futures = [pool.submit(read, True) for _ in range(16)]
                for future in futures:
                    actual_units, actual = future.result()
                    self.assertEqual(actual_units, units)
                    np.testing.assert_array_equal(actual, expected)


if __name__ == "__main__":
    unittest.main()
