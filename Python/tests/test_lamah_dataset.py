# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_lamah_dataset.py
#
# Small independent LamaH fixtures; no downloads or user data required.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Small independent LamaH fixtures; no downloads or user data required."""

import csv
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from lamah_dataset import available, basin_ids, resolve_root, reuse, rows


class LamaHTests(unittest.TestCase):
    def setUp(self):
        """Prepare isolated fixtures for each test."""
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.parent = Path(self.temp.name)
        self.root = self.parent / "LAMAH_CE"
        self.root.mkdir()
        self.table(
            "Gauge_attributes.csv",
            ["ID", "country", "name", "lon", "lat", "area_gov"],
            [
                ["1", "AUT", "Austrian gauge", 4285776, 2684573, 10],
                ["2", "DEU", "German gauge", 4285776, 2684573, 20],
                ["3", "AUT", "Second Austrian gauge", 4285776, 2684573, 30],
                ["4", "AUT", "No total-upstream catchment", 0, 0, 0],
            ],
        )
        self.table(
            "Catchment_attributes.csv",
            ["ID", "area_calc"],
            [["1", 10], ["2", 20], ["3", 30]],
        )
        for folder in ("forcing", "streamflow"):
            target = self.root / "daily" / "timeseries" / folder
            target.mkdir(parents=True)
            for gauge in ("1", "2", "3"):
                (target / f"ID_{gauge}.csv").write_text(
                    "independent immutable fixture\n", encoding="utf-8"
                )

    def table(self, name, fields, values):
        """Read the selected source table."""
        with (self.root / name).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter=";")
            writer.writerow(fields)
            writer.writerows(values)

    def test_daily_hourly_archive_checksums(self):
        """Verify daily hourly archive checksums."""
        from dataset_install import download_options

        resource = Path(__file__).resolve().parents[1] / "resources/sage"
        for region in ("LAMAH_CE", "LAMAH_AT"):
            daily = download_options(resource, region, 1)[0]
            hourly = download_options(resource, region, 24)[0]
            self.assertEqual(
                daily["md5"], "69fd2733e969513403f923ecc5eaa3dc"
            )
            self.assertEqual(
                hourly["md5"], "121b2292b288ef1c2ef4d96250be77ac"
            )
            self.assertIn("daily_hourly", hourly["name"])
            self.assertEqual(download_options(resource, region, 96), [])

    def test_canonical_catalogue_filters_austria(self):
        """Verify canonical catalogue filters austria."""
        from regional_data import catalogue, region_registry

        resource = Path(__file__).resolve().parents[1] / "resources/sage"
        registry = region_registry(str(resource))
        self.assertNotIn("CAMELS" + "_AT", registry)
        self.assertIn("LAMAH_CE", registry)
        self.assertIn("LAMAH_AT", registry)
        for region, expected in (
            ("LAMAH_CE", {"1", "2", "3"}),
            ("LAMAH_AT", {"1", "3"}),
        ):
            actual = catalogue(self.root, region, str(resource))
            self.assertEqual({row[0] for row in actual}, expected)

    def test_country_and_catchment_intersection(self):
        """Verify country and catchment intersection."""
        self.assertEqual(basin_ids(self.root, "LAMAH_AT"), {"1", "3"})
        self.assertEqual(basin_ids(self.root, "LAMAH_CE"), {"1", "2", "3"})

    def test_reuse_filters_metadata_and_preserves_source(self):
        """Verify reuse filters metadata and preserves source."""
        before = (self.root / "Gauge_attributes.csv").read_bytes()
        target = reuse(self.root, self.parent, "LAMAH_AT")
        self.assertTrue(available(target, "LAMAH_AT"))
        self.assertEqual(
            {r["country"] for r in rows(target / "Gauge_attributes.csv")},
            {"AUT"},
        )
        self.assertFalse(
            (target / "daily/timeseries/forcing/ID_2.csv").exists()
        )
        self.assertEqual(
            (self.root / "Gauge_attributes.csv").read_bytes(), before
        )
        self.assertEqual(
            (target / "daily/timeseries/forcing/ID_1.csv").read_bytes(),
            (self.root / "daily/timeseries/forcing/ID_1.csv").read_bytes(),
        )

    def test_both_resolutions_share_canonical_folder(self):
        """Verify both resolutions share canonical folder."""
        daily = reuse(self.root, self.parent, "LAMAH_AT", 1)
        for folder in ("forcing", "streamflow"):
            source = self.root / "daily/timeseries" / folder
            target = self.root / "hourly/timeseries" / folder
            target.mkdir(parents=True)
            for file in source.iterdir():
                (target / file.name).write_bytes(file.read_bytes())
        hourly = reuse(self.root, self.parent, "LAMAH_AT", 24)
        self.assertEqual(daily, self.parent / "LAMAH_AT")
        self.assertEqual(hourly, daily)
        self.assertTrue(available(daily, "LAMAH_AT", 1))
        self.assertTrue(available(daily, "LAMAH_AT", 24))
        self.assertEqual(reuse(self.root, self.parent, "LAMAH_AT"), daily)
        self.assertFalse(list(self.parent.glob("LAMAH_AT_*")))

    def test_archive_install_adds_resolution_to_same_folder(self):
        """Verify archive install adds resolution to same folder."""
        from dataset_install import install_archives

        resource = Path(__file__).resolve().parents[1] / "resources/sage"
        parent = self.parent / "archive_install"
        for folder, variable in (
            ("forcing", "prec"),
            ("streamflow", "qobs"),
        ):
            for file in (self.root / "daily/timeseries" / folder).iterdir():
                file.write_text(f"YYYY;MM;DD;{variable}\n2005;10;1;1\n")
        daily = install_archives(
            [self.root], parent, "LAMAH_AT", resource, resolution=1
        )
        for folder in ("forcing", "streamflow"):
            source = self.root / "daily/timeseries" / folder
            target = self.root / "hourly/timeseries" / folder
            target.mkdir(parents=True)
            for file in source.iterdir():
                (target / file.name).write_text(
                    file.read_text()
                    .replace("DD;", "DD;HH;")
                    .replace("2005;10;1;", "2005;10;1;0;")
                )
        hourly = install_archives(
            [self.root], parent, "LAMAH_AT", resource, resolution=24
        )
        self.assertEqual(daily, parent / "LAMAH_AT")
        self.assertEqual(hourly, daily)
        self.assertTrue(available(hourly, "LAMAH_AT", 1))
        self.assertTrue(available(hourly, "LAMAH_AT", 24))
        self.assertFalse(list(parent.glob("LAMAH_AT_*")))

    def test_required_pairs_and_resolution(self):
        """Verify required pairs and resolution."""
        self.assertTrue(available(self.root, "LAMAH_AT"))
        (self.root / "daily/timeseries/forcing/ID_3.csv").unlink()
        self.assertFalse(available(self.root, "LAMAH_AT"))
        self.assertFalse(available(self.root, "LAMAH_CE", 24))
        self.assertFalse(available(self.root, "LAMAH_CE", 96))

    def test_discover_canonical_subset_and_ce(self):
        """Verify discover canonical subset and ce."""
        subset = reuse(self.root, self.parent, "LAMAH_AT")
        self.assertEqual(
            resolve_root(self.parent / "LAMAH_AT", "LAMAH_AT"), subset
        )
        self.assertEqual(
            resolve_root(self.parent / "LAMAH_CE", "LAMAH_CE"), self.root
        )


if __name__ == "__main__":
    unittest.main()
