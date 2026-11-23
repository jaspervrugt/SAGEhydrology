# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_reader_arrays.py
#
# Preserve missing-value semantics in the shared numeric reader.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Preserve missing-value semantics in the shared numeric reader."""

import unittest
import numpy as np
import pandas as pd
from reader_arrays import mask_missing_values


class ReaderArraysTests(unittest.TestCase):
    def test_numeric_dtypes_and_missing_values_match_legacy(self):
        """Verify numeric dtypes and missing values match legacy."""
        for dtype in (np.float32, np.float64, np.int32, np.int64):
            with self.subTest(dtype=dtype):
                values = np.array(
                    [[1, -999], [-9999, 3], [0, 7]], dtype=dtype
                )
                frame = pd.DataFrame(values, columns=["P", "Q"])
                expected = frame.mask(frame.isin([-999, -9999]))
                pd.testing.assert_frame_equal(
                    mask_missing_values(frame), expected, check_exact=True
                )

    def test_nan_infinity_labels_and_source_ownership(self):
        """Verify nan infinity labels and source ownership."""
        frame = pd.DataFrame(
            [[np.nan, np.inf], [-np.inf, -999], [-9999, -999.1]],
            columns=["temperature", "flow"],
            index=pd.date_range("2000-01-01", periods=3, freq="h"),
        )
        original = frame.copy(deep=True)
        expected = frame.mask(frame.isin([-999, -9999]))
        result = mask_missing_values(frame)
        pd.testing.assert_frame_equal(result, expected, check_exact=True)
        result.iloc[0, 0] = 12
        pd.testing.assert_frame_equal(frame, original, check_exact=True)

    def test_mixed_strings_keep_original_membership_rules(self):
        """Verify mixed strings keep original membership rules."""
        frame = pd.DataFrame(dict(P=[1, -999], label=["-999", "station"]))
        pd.testing.assert_frame_equal(
            mask_missing_values(frame),
            frame.mask(frame.isin([-999, -9999])),
            check_exact=True,
        )


if __name__ == "__main__":
    unittest.main()
