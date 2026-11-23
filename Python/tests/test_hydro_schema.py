# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_hydro_schema.py
#
# Independent MATLAB derivations and generic source-rule regressions.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Independent MATLAB derivations and generic source-rule regressions."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from hydro_schema import resolve
from hydro_derivations import derive, variable
from hydro_io import csv_raw, native, timestamps


class HydroSchemaTests(unittest.TestCase):
    def test_identical_timestamp_columns_share_validated_immutable_clock(
        self,
    ):
        """Verify identical timestamp columns share validated immutable
        clock.
        """
        from hydro_io import fixed_clock

        a = fixed_clock([b"2000-01-01,1"], "yyyy-MM-dd", ",")
        b = fixed_clock([b"2000-01-01,900"], "yyyy-MM-dd", ",")
        self.assertIs(a, b)
        with self.assertRaises(ValueError):
            b[0] = np.datetime64("2001-01-01")
        c = fixed_clock([b"2000-01-02,1"], "yyyy-MM-dd", ",")
        self.assertEqual(c[0], np.datetime64("2000-01-02"))

    def test_named_month_clock_validation(self):
        """Verify named month clock validation."""
        from hydro_io import fixed_clock

        rows = [b"29-Feb-2000 23:59:59,1", b"01-mAr-2000 00:00:00,2"]
        actual = fixed_clock(rows, "dd-MMM-yyyy HH:mm:ss", ",")
        np.testing.assert_array_equal(
            actual,
            np.array(
                ["2000-02-29T23:59:59", "2000-03-01T00:00:00"],
                dtype="datetime64[s]",
            ),
        )
        self.assertIsNone(
            fixed_clock(
                [b"01-XYZ-2000 00:00:00,1"], "dd-MMM-yyyy HH:mm:ss", ","
            )
        )
        quoted = fixed_clock(
            [b'"01-Jan-2000 00:00:00",1'], "dd-MMM-yyyy HH:mm:ss", ","
        )
        self.assertEqual(quoted[0], np.datetime64("2000-01-01"))
        self.assertIsNone(
            fixed_clock(
                [b'"01-Jan-2000 00:00:00x",1'], "dd-MMM-yyyy HH:mm:ss", ","
            )
        )
        with self.assertRaises(ValueError):
            fixed_clock(
                [b"29-Feb-1900 00:00:00,1"], "dd-MMM-yyyy HH:mm:ss", ","
            )

    def test_calendar_components_match_original_validation(self):
        """Verify calendar components match original validation."""
        from hydro_io import calendar_clock, timestamps

        frame = pd.DataFrame(
            dict(
                y=["2000", "1900"],
                m=["2", "3"],
                d=["29", "1"],
                h=["23", "0"],
            )
        )
        spec = dict(
            mode="ymd_columns",
            year_column="y",
            month_column="m",
            day_column="d",
            hour_column="h",
        )
        actual = timestamps(frame, spec)
        expected = pd.to_datetime(
            dict(year=[2000, 1900], month=[2, 3], day=[29, 1], hour=[23, 0])
        ).to_numpy(dtype="datetime64[s]")
        np.testing.assert_array_equal(actual, expected)
        frame.loc[0, "h"] = "24"
        self.assertEqual(
            timestamps(frame, spec)[0], np.datetime64("2000-03-01T00:00:00")
        )
        frame.loc[0, "h"] = "0.5"
        self.assertEqual(
            timestamps(frame, spec)[0], np.datetime64("2000-02-29T00:30:00")
        )
        frame.loc[1, "m"] = "2"
        frame.loc[1, "d"] = "29"
        with self.assertRaises(ValueError):
            timestamps(frame, spec)

    def test_window_preserves_exact_values_and_checks_full_clock(self):
        """Verify window preserves exact values and checks full clock."""
        import hydro_io

        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "source.csv"
            path.write_text(
                "date,value\n2000-01-01,1.0000000596046448\n"
                "2000-01-02,2.0000001192092896\n2000-01-03,3\n"
            )
            schema = dict(
                file=dict(delimiter=","),
                strict_time=True,
                time=dict(
                    mode="column", column="date", input_format="yyyy-MM-dd"
                ),
                variables=dict(P=dict(source="value")),
            )
            target = np.array(["2000-01-02"], dtype="datetime64[s]")
            with patch.object(hydro_io, "_CSV_WINDOW_MIN_BYTES", 0):
                result = hydro_io.csv_raw(path, schema, target)
                # Invalid chronology outside the requested window still
                # fails, rather than being hidden by the optimization.
                path.write_text(path.read_text() + "1999-01-01,4\n")
                with self.assertRaisesRegex(ValueError, "increasing"):
                    hydro_io.csv_raw(path, schema, target)
            self.assertEqual(result["P"][0], float("2.0000001192092896"))

    def test_wide_selected_basins_share_source_batch_and_invalidate(self):
        """Verify wide selected basins share source batch and invalidate."""
        import hydro_io
        from local_data import capture_sources

        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "wide.csv"
            path.write_text("date,1,2,3\n2000-01-01,1,2,3\n")
            spec = dict(delimiter=",")
            available = {str(i): str(i) for i in (1, 2, 3)}
            with capture_sources(["1", "2", "3"]), patch.object(
                hydro_io, "table", wraps=hydro_io.table
            ) as parsed:
                for column in available:
                    frame = hydro_io.wide_table(
                        path, spec, ["date"], column, available, {}
                    )
                    self.assertEqual(frame[column].iloc[0], int(column))
                self.assertEqual(parsed.call_count, 1)
                path.write_text("date,1,2,3\n2000-01-01,10,20,30\n")
                frame = hydro_io.wide_table(
                    path, spec, ["date"], "1", available, {}
                )
                self.assertEqual(frame["1"].iloc[0], 10)
                self.assertEqual(parsed.call_count, 2)

    def test_packaged_catalogue_uses_bundled_resources(self):
        """Verify packaged catalogue uses bundled resources."""
        from hydro_schema import catalogue
        import sys

        with tempfile.TemporaryDirectory() as folder:
            resource = Path(folder) / "resources/sage"
            resource.mkdir(parents=True)
            (resource / "hydro_schemas.json").write_text('{"fixture": 1}')
            with patch.object(sys, "_MEIPASS", folder, create=True):
                self.assertEqual(catalogue(None), {"fixture": 1})

    def test_all_derivations_against_matlab_materializer(self):
        """Verify all derivations against matlab materializer."""
        path = (
            Path(__file__).parent
            / "tests/data"
            / "matlab_hydro_derivations.npz"
        )
        with np.load(path, allow_pickle=False) as saved:
            cases = json.loads(str(saved["metadata"].item()))
            raw = {
                k: saved[k]
                for k in saved.files
                if k not in ("metadata", "expected", "dates")
            }
            for index, case in enumerate(cases):
                with self.subTest(derivation=case["name"]):
                    data = derive(
                        dict(raw),
                        case["schema"],
                        saved["dates"],
                        np.array([45.0, 100.0, 1e8]),
                    )
                    values = variable(
                        data["Ep"],
                        case["schema"]["variables"]["Ep"],
                        [45.0, 100.0, 1e8],
                        86400,
                    )
                    np.testing.assert_array_equal(
                        values.astype(np.float32), saved["expected"][index]
                    )

    def test_profile_merge_and_unique_match(self):
        """Verify profile merge and unique match."""
        schema = dict(
            timeline=dict(step=86400),
            variables=dict(Q=dict(source="flow", units="m3/s")),
            profiles=dict(
                hourly=dict(
                    match=dict(dt=24),
                    schema=dict(
                        timeline=dict(step=3600),
                        variables=dict(Q=dict(scale=1000)),
                    ),
                )
            ),
        )
        result = resolve(schema, dict(dt=24))
        self.assertEqual(result["timeline"]["step"], 3600)
        self.assertEqual(result["variables"]["Q"]["source"], "flow")
        self.assertEqual(result["variables"]["Q"]["scale"], 1000)
        with self.assertRaisesRegex(ValueError, "Expected one"):
            resolve(schema, dict(dt=96))

    def test_native_missing_then_units_clip_fill_transform(self):
        """Verify native missing then units clip fill transform."""
        spec = dict(
            missing_values=[-999],
            invalid_le=-99,
            valid_min=-2,
            scale=2,
            offset=-1,
            clip_min=-3,
            fill_missing=7,
            area_normalize=True,
            transform=dict(matlab_function="@(x,context)abs(x)"),
        )
        raw = native([-999, -100, -3, -2, 4], spec)
        result = variable(raw, spec, [0, 0, 1e6], 1000)
        np.testing.assert_array_equal(result, [7, 7, 7, 3, 7])

    def test_named_csv_alias_quality_and_calendar_join(self):
        """Verify named csv alias quality and calendar join."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "flow.csv"
            path.write_text(
                "date,alternate,flag\n2000-01-01,2,0\n2000-01-02,3,1\n"
            )
            schema = dict(
                file={},
                time=dict(
                    mode="column", column="date", join="calendar_day"
                ),
                variables=dict(
                    Q=dict(
                        sources=["absent", "alternate"],
                        source_operation="aliases",
                        quality=dict(
                            sources=["flag"], valid_min=0, valid_max=0
                        ),
                    )
                ),
            )
            dates = np.array(
                [
                    "2000-01-01T09:00:00",
                    "2000-01-02T09:00:00",
                    "2000-01-03T09:00:00",
                ],
                dtype="datetime64[s]",
            )
            result = csv_raw(path, schema, dates)
            np.testing.assert_array_equal(result["Q"], [2, np.nan, np.nan])
            path.write_text("date,alternate\n2000-01-01,2\n")
            with self.assertRaises(ValueError):
                csv_raw(path, schema, dates)

    def test_partial_first_day_uses_final_slots(self):
        """Verify partial first day uses final slots."""
        frame = pd.DataFrame(dict(y=[2000] * 3, m=[1] * 3, d=[1, 1, 2]))
        spec = dict(
            mode="ymd_row_sequence",
            year_column="y",
            month_column="m",
            day_column="d",
            steps_per_day=24,
            partial_first="final_slots",
        )
        np.testing.assert_array_equal(
            timestamps(frame, spec),
            np.array(
                [
                    "2000-01-01T22:00:00",
                    "2000-01-01T23:00:00",
                    "2000-01-02T00:00:00",
                ],
                dtype="datetime64[s]",
            ),
        )

    def test_blank_delimiter_placeholder_is_missing_series(self):
        """Verify blank delimiter placeholder is missing series."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.csv"
            path.write_text(",,,\n,,,\n")
            schema = dict(variables=dict(Q=dict(source="flow")))
            result = csv_raw(path, schema, np.arange(3))
            self.assertTrue(np.isnan(result["Q"]).all())


if __name__ == "__main__":
    unittest.main()
