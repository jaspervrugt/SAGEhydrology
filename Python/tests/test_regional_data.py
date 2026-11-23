# ==========================================================================
# SAGEhydrology | Python source
# File: tests/test_regional_data.py
#
# Reader and installer regression checks using small, independent
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Reader and installer regression checks using small, independent
fixtures.
"""

from pathlib import Path
import tempfile
import threading
import time
import unittest
import zipfile
import numpy as np
import pandas as pd
from regional_data import read_region, region_registry, reader_issues
from dataset_install import (
    extract,
    install_archives,
    InstallCancelled,
    download_options,
)

SAGE = __import__("sage_paths").SAGE_RESOURCES


class RegionalTests(unittest.TestCase):
    def setUp(self):
        """Prepare isolated fixtures for each test."""
        base = (
            Path(__file__).resolve().parents[2]
            / "benchmarks/Python/reader_tests"
        )
        base.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=base)
        self.root = Path(self.temp.name).resolve()
        assert self.root.is_relative_to(base.resolve())
        self.period = dict(
            train_start="2005-10-01",
            train_end="2005-10-02",
            eval_start="2005-10-03",
            eval_end="2005-10-04",
            spinup=0,
        )

    def tearDown(self):
        """Release fixtures and restore the state changed by the test."""
        for attempt in range(15):
            try:
                self.temp.cleanup()
                return
            except PermissionError:
                if attempt == 14:
                    raise
                time.sleep(0.2)

    def poland_archive(self):
        """Construct the Poland dataset archive fixture."""
        archive = self.root / "poland.zip"
        dates = pd.date_range("2005-10-01", periods=4)
        frame = pd.DataFrame(
            dict(
                date=dates,
                precipitation_mean=[1, 2, 3, 4],
                temperature_mean=[10] * 4,
                minimum_temperature_mean=[5] * 4,
                maximum_temperature_mean=[15] * 4,
                discharge_spec_obs=[1, -999, 3, 4],
            )
        )
        attributes = pd.DataFrame(
            dict(
                gauge_id=["149180020"],
                gauge_lon=[50],
                gauge_lat=[19],
                area_metadata=[100],
                area_km2=[100],
                gauge_elev=[200],
            )
        )
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr(
                "CAMELS-PL/timeseries/CAMELS_PL_hydromet_timeseries"
                "_149180020.csv",
                frame.to_csv(index=False),
            )
            handle.writestr(
                "CAMELS-PL/CAMELS_PL_topographic_attributes.csv",
                attributes.to_csv(index=False),
            )
            handle.writestr(
                "CAMELS-PL/gauge_information_PL.csv",
                attributes.to_csv(index=False),
            )
        return archive

    def test_install_poland_then_read_units_masks_and_latitude(self):
        """Verify install poland then read units masks and latitude."""
        installed = install_archives(
            [self.poland_archive()], self.root / "data", "CAMELS_PL", SAGE
        )
        data = read_region(
            installed,
            "149180020",
            region="CAMELS_PL",
            sage_root=SAGE,
            **self.period
        )
        np.testing.assert_array_equal(data["forcing"][:, 0], [1, 2, 3, 4])
        self.assertTrue(np.all(data["forcing"][:, 1] > 0))
        np.testing.assert_array_equal(
            data["train"], [True, False, False, False]
        )
        np.testing.assert_array_equal(
            data["evaluation"], [False, False, True, True]
        )
        self.assertEqual(data["observed"][0], 1)

    def test_installs_never_replace_existing_directory(self):
        """Verify installs never replace existing directory."""
        archive = self.poland_archive()
        parent = self.root / "data"
        parent.mkdir()
        existing = parent / "CAMELS_PL"
        existing.mkdir()
        sentinel = existing / "keep.txt"
        sentinel.write_text("original")
        a = install_archives([archive], parent, "CAMELS_PL", SAGE)
        b = install_archives([archive], parent, "CAMELS_PL", SAGE)
        self.assertNotEqual(a, b)
        self.assertEqual(sentinel.read_text(), "original")

    def test_archive_traversal_rejected(self):
        """Verify archive traversal rejected."""
        archive = self.root / "bad.zip"
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr("../escape.txt", "bad")
        with self.assertRaises(ValueError):
            extract(archive, self.root / "extracted")
        self.assertFalse((self.root / "escape.txt").exists())

    def test_cancel_stops_installation(self):
        """Verify cancel stops installation."""
        stop = threading.Event()
        stop.set()
        with self.assertRaises(InstallCancelled):
            install_archives(
                [self.poland_archive()],
                self.root / "data",
                "CAMELS_PL",
                SAGE,
                stop=stop,
            )

    def test_all_region_training_routes_and_download_recipe(self):
        """Verify all region training routes and download recipe."""
        registry = region_registry(str(SAGE))
        self.assertEqual(len(registry), 47)
        self.assertIn("LAMAH_CE", registry)
        self.assertIn("LAMAH_AT", registry)
        self.assertNotIn("CAMELS" + "_AT", registry)
        for code, meta in registry.items():
            resolution = 1 if "Daily" in meta["resolutions"] else 24
            self.assertEqual(
                reader_issues(SAGE, code, resolution), [], code
            )
        options = download_options(SAGE, "CAMELS_PL")
        self.assertEqual(
            options[0]["md5"], "0ff7d9635c941f4aac49c191376502c9"
        )

    def test_australia_wide_file_product_selection(self):
        """Verify australia wide file product selection."""
        base = self.root / "daily"
        dates = pd.date_range("2005-10-01", periods=4)
        for path, value in [
            ("precipitation/precipitation_AGCD.csv", 7),
            ("evaporative_demand/evap_syn_SILO.csv", 2),
            ("temperature/tmin_SILO.csv", 5),
            ("temperature/tmax_SILO.csv", 15),
            ("streamflow/streamflow_mmd.csv", 1),
        ]:
            file = base / path
            file.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                dict(
                    year=dates.year,
                    month=dates.month,
                    day=dates.day,
                    **{"001": [value] * 4}
                )
            ).to_csv(file, index=False)
        data = read_region(
            self.root,
            "001",
            region="CAMELS_AU",
            sage_root=SAGE,
            precipitation=2,
            **self.period
        )
        np.testing.assert_array_equal(
            data["forcing"], np.tile([7, 2, 10], (4, 1))
        )

    def test_brazil_multifile_join_and_volumetric_discharge(self):
        """Verify brazil multifile join and volumetric discharge."""
        dates = pd.date_range("2005-10-01", periods=4)
        base = self.root / "daily"
        sources = {
            "precipitation/123_precipitation.txt": dict(
                p_brdwgd=[1] * 4,
                p_chirps=[2] * 4,
                p_cpc=[3] * 4,
                p_era5land=[4] * 4,
                p_mswep=[5] * 4,
            ),
            "pet/123_potential_evapotransp.txt": dict(
                pet_gleam=[2] * 4, pet_era5land=[4] * 4
            ),
            "temperature/123_temperature.txt": dict(
                tmax_cpc=[15] * 4,
                tmax_era5land=[15] * 4,
                tmax_brdwgd=[15] * 4,
                tmean_era5land=[10] * 4,
                tmin_cpc=[5] * 4,
                tmin_era5land=[5] * 4,
                tmin_brdwgd=[5] * 4,
            ),
            "streamflow/123_streamflow.txt": dict(
                streamflow_m3s=[1] * 4,
                qual_control_by_ana=[1, 0, 1, 1],
                qual_flag=[1] * 4,
            ),
        }
        for path, columns in sources.items():
            file = base / path
            file.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(
                dict(
                    year=dates.year,
                    month=dates.month,
                    day=dates.day,
                    **columns
                )
            ).to_csv(file, index=False, sep=" ")
        pd.DataFrame(
            dict(
                gauge_id=[123],
                gauge_lat=[-20],
                gauge_elev=[100],
                area_km2=[100],
            )
        ).to_csv(self.root / "gauge_information.txt", index=False, sep="\t")
        data = read_region(
            self.root,
            "123",
            region="CAMELS_BR",
            sage_root=SAGE,
            precipitation=0,
            pet=0,
            **self.period
        )
        np.testing.assert_array_equal(
            data["forcing"], np.tile([3, 3, 10], (4, 1))
        )
        self.assertAlmostEqual(data["observed"][0], 0.864)
        self.assertFalse(data["train"][1])
        for precip, evap, expected in [
            (1, 1, [1, 2, 10]),
            (5, 2, [5, 4, 10]),
        ]:
            selected = read_region(
                self.root,
                "123",
                region="CAMELS_BR",
                sage_root=SAGE,
                precipitation=precip,
                pet=evap,
                **self.period
            )
            np.testing.assert_array_equal(
                selected["forcing"], np.tile(expected, (4, 1))
            )

    def test_brazil_static_metadata_reused_between_basins(self):
        """Verify brazil static metadata reused between basins."""
        from unittest.mock import patch
        import regional_data

        pd.DataFrame(
            dict(gauge_id=[123, 456], area=[100, 200], elev_mean=[10, 20])
        ).to_csv(
            self.root / "camels_br_topography.txt", index=False, sep=" "
        )
        pd.DataFrame(
            dict(gauge_id=[123, 456], gauge_lat=[-20, -21])
        ).to_csv(self.root / "camels_br_location.txt", index=False, sep=" ")
        profile = {"meta": {"acronym": "BR"}, "auxiliary_tables": []}
        regional_data._static_table.cache_clear()
        original = regional_data._parse_table
        with patch.object(
            regional_data, "_parse_table", wraps=original
        ) as parse:
            first = regional_data._metadata(profile, self.root, "123")
            second = regional_data._metadata(profile, self.root, "456")
        self.assertEqual(parse.call_count, 2)
        self.assertEqual(first, dict(area=100, elev=10, lat=-20))
        self.assertEqual(second, dict(area=200, elev=20, lat=-21))

    def test_regional_static_fallbacks_reused_between_basins(self):
        """Verify regional static fallbacks reused between basins."""
        import regional_data
        from unittest.mock import patch

        fixtures = [
            (
                "SE",
                "catchments_physical_properties.csv",
                dict(
                    ID=[123, 456],
                    Latitude_WGS84=[50, 51],
                    Area_km2=[100, 200],
                    Elevation_mabsl=[10, 20],
                ),
            ),
            (
                "PL",
                "CAMELS_PL_topographic_attributes.csv",
                dict(
                    gauge_id=[123, 456],
                    gauge_lon=[50, 51],
                    area_metadata=[100, 200],
                    gauge_elev=[10, 20],
                ),
            ),
            (
                "EU",
                "CHAMELLION_topographic_attributes.csv",
                dict(
                    gauge_id=[123, 456],
                    gauge_lat=[50, 51],
                    area=[100, 200],
                    gauge_elev=[10, 20],
                ),
            ),
            (
                "NZ",
                "1.CAMELS_NZ_Catchment_information.csv",
                {
                    "RID": [123, 456],
                    "Latitude (WGS 84)": [50, 51],
                    "uparea": [100, 200],
                    "elevation": [10, 20],
                },
            ),
            (
                "KRH",
                "attributes_general.csv",
                dict(STAID=[123, 456], Lat_gage=[50, 51], Area=[100, 200]),
            ),
        ]
        for acronym, name, values in fixtures:
            with self.subTest(region=acronym):
                pd.DataFrame(values).to_csv(self.root / name, index=False)
                regional_data._static_table.cache_clear()
                profile = dict(
                    meta=dict(acronym=acronym), auxiliary_tables=[]
                )
                with patch.object(
                    regional_data,
                    "_parse_table",
                    wraps=regional_data._parse_table,
                ) as parser:
                    first = regional_data._metadata(
                        profile, self.root, "123"
                    )
                    second = regional_data._metadata(
                        profile, self.root, "456"
                    )
                self.assertEqual(parser.call_count, 1)
                self.assertEqual(first["area"], 100)
                self.assertEqual(second["area"], 200)

    def test_canada_netcdf_units_and_quality_mask(self):
        """Verify canada netcdf units and quality mask."""
        from netCDF4 import Dataset

        for folder, name, variables in [
            (
                "forcing/daymet",
                "CAN_01A_daymet_lumped.nc",
                dict(
                    prcp=[2] * 4, pet=[1] * 4, tmin=[5] * 4, tmax=[15] * 4
                ),
            ),
            (
                "discharge",
                "CAN_01A_daily_flow_observations.nc",
                dict(
                    q_obs=[1] * 4,
                    q_obs_is_ice_affected=[0, 1, 0, 0],
                    q_obs_is_partial_day=[0] * 4,
                ),
            ),
        ]:
            directory = self.root / "daily" / folder
            directory.mkdir(parents=True, exist_ok=True)
            with Dataset(directory / name, "w") as handle:
                handle.createDimension("time", 4)
                clock = handle.createVariable("time", "f8", ("time",))
                clock.units = "days since 2005-10-01"
                clock[:] = np.arange(4)
                for key, value in variables.items():
                    handle.createVariable(key, "f8", ("time",))[:] = value
        pd.DataFrame(
            dict(
                Station_id=["01A"],
                Station_lat=[50],
                Basin_area_km2=[100],
            )
        ).to_csv(self.root / "camels-spat-metadata.csv", index=False)
        pd.DataFrame({"Attribute": ["elev_mean"], "CAN_01A": [100]}).to_csv(
            self.root / "attributes-lumped.csv", index=False
        )
        data = read_region(
            self.root,
            "01A",
            region="CAMELS_CA",
            sage_root=SAGE,
            **self.period
        )
        np.testing.assert_array_equal(
            data["forcing"], np.tile([2, 1, 10], (4, 1))
        )
        self.assertAlmostEqual(data["observed"][0], 0.864)
        self.assertFalse(data["train"][1])

    def test_korean_hourly_sign_units_and_complete_days(self):
        """Verify korean hourly sign units and complete days."""
        directory = self.root / "hourly" / "timeseries"
        directory.mkdir(parents=True)
        pd.DataFrame(
            dict(
                DateTime=pd.date_range(
                    "2005-10-01", periods=96, freq="h"
                ).strftime("%d-%b-%Y %H:%M:%S"),
                total_precipitation=[2] * 96,
                potential_evaporation=[-0.1] * 96,
                temperature_2m=[10] * 96,
                streamflow=[1] * 96,
            )
        ).to_csv(directory / "1001620.csv", index=False)
        pd.DataFrame(
            dict(gauge_id=[1001620], gauge_lat=[37], area_km2=[100])
        ).to_csv(self.root / "gauge_information.txt", sep="\t", index=False)
        data = read_region(
            self.root,
            "1001620",
            region="CAMELSH_KR",
            sage_root=SAGE,
            resolution=24,
            **self.period
        )
        self.assertEqual(len(data["dates"]), 96)
        self.assertEqual(data["train"].sum(), 48)
        self.assertEqual(data["evaluation"].sum(), 48)
        np.testing.assert_allclose(
            data["forcing"],
            np.tile(np.array([2, 0.1, 10]), (96, 1)),
            rtol=1e-7,
        )
        self.assertTrue(np.issubdtype(data["forcing"].dtype, np.floating))
        self.assertTrue(np.issubdtype(data["observed"].dtype, np.floating))
        self.assertAlmostEqual(data["observed"][0], 0.036)

    def test_austrian_raw_archive_keeps_forcing_and_flow_separate(self):
        """Verify austrian raw archive keeps forcing and flow separate."""
        archive = self.root / "lamah.zip"
        dates = pd.date_range("2005-10-01", periods=4)
        clock = dict(YYYY=dates.year, MM=dates.month, DD=dates.day)
        meteo = pd.DataFrame(
            dict(
                **clock,
                prec=[3] * 4,
                total_et=[-2] * 4,
                **{"2m_temp_mean": [10] * 4}
            )
        )
        flow = pd.DataFrame(dict(**clock, qobs=[1] * 4))
        gauges = pd.DataFrame(
            dict(
                ID=[1],
                name=["Test gauge"],
                country=["AUT"],
                elev=[420],
                lon=[4285776],
                lat=[2684573],
            )
        )
        catchments = pd.DataFrame(
            dict(ID=[1], area_calc=[100], elev_mean=[200])
        )
        with zipfile.ZipFile(archive, "w") as handle:
            for path, frame in [
                (
                    "A_basins_total_upstrm/2_timeseries/daily/ID_1.csv",
                    meteo,
                ),
                ("D_gauges/2_timeseries/daily/ID_1.csv", flow),
                ("D_gauges/1_attributes/Gauge_attributes.csv", gauges),
                (
                    (
                        "A_basins_total_upstrm/1_attributes/"
                        "Catchment_attributes.csv"
                    ),
                    catchments,
                ),
            ]:
                handle.writestr(path, frame.to_csv(index=False, sep=";"))
        installed = install_archives(
            [archive], self.root / "data", "LAMAH_CE", SAGE
        )
        self.assertTrue(
            (installed / "daily/timeseries/forcing/ID_1.csv").exists()
        )
        self.assertTrue(
            (installed / "daily/timeseries/streamflow/ID_1.csv").exists()
        )
        data = read_region(
            installed,
            "1",
            region="LAMAH_CE",
            sage_root=SAGE,
            pet=4,
            **self.period
        )
        np.testing.assert_array_equal(
            data["forcing"], np.tile([3, 2, 10], (4, 1))
        )
        self.assertAlmostEqual(data["observed"][0], 0.864)


if __name__ == "__main__":
    unittest.main()
