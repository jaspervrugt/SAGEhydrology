# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_wide_data_cache.py
#
# Verify one CSV parse, concurrent reuse and source invalidation.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Verify one CSV parse, concurrent reuse and source invalidation."""

import os
from pathlib import Path
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch
import numpy as np
import pandas as pd
from wide_data_cache import series, _mapped


class WideDataCacheTest(unittest.TestCase):
    def test_shared_columns_and_invalidation(self):
        """Verify shared columns and invalidation."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "wide.csv"
            source.write_text(
                "year,month,day,001,A2\n2000,1,1,1.2,3\n"
                "2000,1,2,-999,4\n2000,1,3,,5\n",
                encoding="utf-8",
            )
            with patch.dict(os.environ, LOCALAPPDATA=str(root / "cache")):
                raw = pd.read_csv
                with patch(
                    "wide_data_cache.pd.read_csv", wraps=raw
                ) as count:
                    with ThreadPoolExecutor(max_workers=16) as pool:
                        result = list(
                            pool.map(
                                lambda _: series(
                                    source, "1", ("year", "month", "day")
                                ),
                                range(32),
                            )
                        )
                    assert count.call_count == 1, count.call_count
                for value in result:
                    np.testing.assert_array_equal(
                        value.to_numpy(), [1.2, -999, np.nan]
                    )
                np.testing.assert_array_equal(
                    series(
                        source, "A2", ("year", "month", "day")
                    ).to_numpy(),
                    [3, 4, 5],
                )
                with patch(
                    "wide_data_cache.pd.read_csv",
                    side_effect=AssertionError("reparse"),
                ):
                    series(source, "001", ("year", "month", "day"))
                source.write_text(source.read_text().replace("1.2", "9.25"))
                assert (
                    series(source, "001", ("year", "month", "day")).iloc[0]
                    == 9.25
                )
                try:
                    series(source, "missing", ("year", "month", "day"))
                except ValueError:
                    pass
                else:
                    raise AssertionError("missing ID")
                del result, value
                _mapped.cache_clear()
                print(
                    "PASS: one parse, concurrent reuse, exact values, "
                    "source invalidation."
                )


if __name__ == "__main__":
    unittest.main()
