# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_timeseries_storage.py
#
# Canonical series precision policy and preserved calendar/mask types.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Canonical series precision policy and preserved calendar/mask types."""

import unittest
import numpy as np
from timeseries_storage import store_timeseries


class StorageTests(unittest.TestCase):
    def test_numeric_storage_and_metadata(self):
        """Verify numeric storage and metadata."""
        dates = np.array(["2000-01-01"], dtype="datetime64[D]")
        mask = np.array([True])
        attrs = np.array([0.123456789], dtype=np.float64)
        data = dict(
            forcing=np.array([[0.1, 0.2, 0.3]]),
            observed=np.array([np.nan]),
            observed_SWE=np.array([1.1]),
            dates=dates,
            train=mask,
            attributes=attrs,
        )
        stored = store_timeseries(data)
        for key in ["forcing", "observed", "observed_SWE"]:
            self.assertEqual(stored[key].dtype, np.float32)
            self.assertEqual(stored[key].nbytes, data[key].nbytes // 2)
        self.assertTrue(np.isnan(stored["observed"][0]))
        self.assertIs(stored["dates"], dates)
        self.assertIs(stored["train"], mask)
        self.assertIs(stored["attributes"], attrs)
        self.assertEqual(data["forcing"].dtype, np.float64)

    def test_idempotent_no_extra_allocation(self):
        """Verify idempotent no extra allocation."""
        series = np.array([[0.1, 0.2, 0.3]], dtype=np.float32)
        self.assertIs(
            store_timeseries(dict(forcing=series))["forcing"], series
        )


if __name__ == "__main__":
    unittest.main()
