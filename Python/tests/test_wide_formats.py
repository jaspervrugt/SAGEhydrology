# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_wide_formats.py
#
# Test wide formats.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

import os, sys, tempfile, unittest
from pathlib import Path
import numpy as np
import pandas as pd
import wide_data_cache as cache


class WideTests(unittest.TestCase):
    def setUp(self):
        """Prepare isolated fixtures for each test."""
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)
        self.old = os.environ.get("LOCALAPPDATA")
        os.environ["LOCALAPPDATA"] = str(self.base)

    def tearDown(self):
        """Release fixtures and restore the state changed by the test."""
        cache._mapped.cache_clear()
        if self.old is None:
            os.environ.pop("LOCALAPPDATA", None)
        else:
            os.environ["LOCALAPPDATA"] = self.old
        self.tmp.cleanup()

    def test_separators_missing_and_invalidation(self):
        """Verify separators missing and invalidation."""
        for sep in [",", ";", "\t", " "]:
            p = self.base / ("table" + str(ord(sep)) + ".csv")
            pd.DataFrame(
                {
                    "gauge_id": ["2000-01-01", "2000-01-02"],
                    "001": [1, "bad"],
                    "002": [3, 4],
                }
            ).to_csv(p, sep=sep, index=False)
            s = cache.series(p, "1", date="gauge_id")
            np.testing.assert_allclose(s, [1, np.nan], equal_nan=True)
            np.testing.assert_allclose(
                cache.series(p, "2", date="gauge_id"), [3, 4]
            )
            p.write_text(p.read_text().replace("3", "5"))
            np.testing.assert_allclose(
                cache.series(p, "2", date="gauge_id"), [5, 4]
            )

    def test_sequential_across_chunks(self):
        """Verify sequential across chunks."""
        p = self.base / "seq.csv"
        pd.DataFrame(
            {"day": ["2000-01-01"] * 2200, "A": np.arange(2200)}
        ).to_csv(p, index=False)
        s = cache.series(p, "A", date="day", sequential=True)
        self.assertEqual(
            s.index[-1],
            pd.Timestamp("2000-01-01") + pd.Timedelta(days=2199),
        )

    def test_duplicate_normalized_dates_rejected(self):
        """Verify duplicate normalized dates rejected."""
        p = self.base / "dup.csv"
        p.write_text("date,A\n2000-01-01 01:00,1\n2000-01-01 02:00,2\n")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            cache.series(p, "A", date="date")

    def test_ymd(self):
        """Verify ymd."""
        p = self.base / "ymd.csv"
        p.write_text("year,month,day,A\n2000,1,1,1\n2000,1,2,2\n")
        np.testing.assert_array_equal(
            cache.series(p, "A", ("year", "month", "day")), [1, 2]
        )

    def test_selected_columns_share_one_index(self):
        """Verify selected columns share one index."""
        p = self.base / "selected.csv"
        p.write_text("year,month,day,001,002,003\n2000,1,1,1,2,3\n")
        first = cache.series(
            p, "1", ("year", "month", "day"), gauges=["001", "002"]
        )
        second = cache.series(
            p, "002", ("year", "month", "day"), gauges=["001", "002"]
        )
        np.testing.assert_array_equal(first, [1])
        np.testing.assert_array_equal(second, [2])
        files = list(
            (self.base / "SAGEhydrology/Python/wide_tables").glob(
                "*/table.json"
            )
        )
        import json

        self.assertEqual(len(files), 1)
        self.assertEqual(
            json.loads(files[0].read_text())["columns"], ["001", "002"]
        )


if __name__ == "__main__":
    unittest.main()
