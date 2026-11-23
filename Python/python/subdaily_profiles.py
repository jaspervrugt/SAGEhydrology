# ==========================================================================
# SAGEhydrology | Python source
# File: python/subdaily_profiles.py
#
# Remaining SAGE hourly profiles and joins, without MATLAB execution.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Remaining SAGE hourly profiles and joins, without MATLAB execution."""

from pathlib import Path
import re
import numpy as np
import pandas as pd
from regional_data import _table, _metadata, table_profile, strings
from regional_subdaily import hourly_pet
from reader_arrays import mask_missing_values

EXTRA = {
    "LAMAH_CE",
    "LAMAH_AT",
    "LAMAH_CE",
    "CAMELS_IS",
    "CAMELS_CA",
    "CAMELS_SPAT",
    "CAMELS_CZ",
    "CAMELS_GB",
    "CAMELS_NZ",
    "CAMELS_US",
}


def one(root, pattern):
    """Resolve a single subdaily profile source."""
    files = list(Path(root).glob(pattern))
    if len(files) != 1:
        raise FileNotFoundError(
            f"Expected one {Path(root)/pattern}; found {len(files)}"
        )
    return files[0]


def locations(root, region):
    """Resolve the subdaily profile file locations."""
    from lamah_dataset import reader_region

    region = reader_region(region)
    base = Path(root) / "hourly"
    if region == "LAMAH_CE":
        return base / "timeseries" / "forcing", "ID_{gauge}.csv"
    if region == "CAMELS_IS":
        return base / "forcing", "ID_{gauge}.csv"
    if region in ("CAMELS_CA", "CAMELS_SPAT"):
        return (
            base / "forcing",
            ("CAN_" if region == "CAMELS_CA" else "")
            + "{gauge}_rdrs_lumped.nc",
        )
    if region == "CAMELS_GB":
        return (
            base / "timeseries",
            "camels_gb_v2_hydromet_hourly_timeseries_{gauge}_*.csv",
        )
    if region == "CAMELS_NZ":
        return (
            base / "timeseries" / "precipitation",
            "precipitation_station_id_{gauge}.csv",
        )
    if region == "CAMELS_CZ":
        return base / "timeseries", "camelscz_{gauge}.csv"
    if region == "CAMELS_US":
        return base / "forcing", "{gauge}_hourly_nldas.csv"
    return base / "timeseries", "{gauge}.csv"


def inventory(root, region):
    """Inspect the files available for the regional reader."""
    folder, pattern = locations(root, region)
    regex = (
        re.escape(pattern)
        .replace(re.escape("{gauge}"), r"(?P<gauge>.+?)")
        .replace(r"\*", ".*")
    )
    gauges = {
        re.fullmatch(regex, f.name).group("gauge")
        for f in folder.glob("**/" + pattern.replace("{gauge}", "*"))
        if re.fullmatch(regex, f.name)
    }
    if region == "CAMELS_US":
        gauges &= {
            f.name.removesuffix("-usgs-hourly.csv")
            for f in (Path(root) / "hourly/streamflow").glob(
                "*-usgs-hourly.csv"
            )
        }
    return [(g, "", "") for g in sorted(gauges)]


def indexed(
    file,
    date="date",
    ymd=False,
    strict=True,
    row_sequence=False,
    columns=None,
):
    """Read a gauge-indexed subdaily source."""
    frame = _table(file, columns=columns, transient=True)
    if ymd:
        clock = dict(year=frame["YYYY"], month=frame["MM"], day=frame["DD"])
        if not row_sequence:
            clock["hour"] = frame["hh"]
        if "mm" in frame:
            clock["minute"] = frame["mm"]
        dates = pd.to_datetime(clock)
        if row_sequence:
            groups = dates.ne(dates.shift()).cumsum()
            slots = dates.groupby(groups).cumcount()
            counts = dates.groupby(groups).transform("size")
            if (counts > 24).any():
                raise ValueError(
                    f"More than 24 discharge rows per day in {file}"
                )
            slots = slots + np.where(
                groups == groups.iloc[0], 24 - counts, 0
            )
            dates = dates + pd.to_timedelta(slots, unit="h")
    else:
        text = frame[date].astype(str)
        try:
            dates = pd.to_datetime(text, format="ISO8601")
        except (ValueError, TypeError):
            dates = pd.to_datetime(text, format="mixed")
    frame.index = pd.DatetimeIndex(dates)
    if frame.index.duplicated().any():
        if strict:
            raise ValueError(f"Duplicate timestamps in {file}")
        # Iceland declares strict_time=false. MATLAB ismember selects the
        # first
        # source occurrence when mapping these publisher timestamps to the
        # grid.
        frame = frame.loc[~frame.index.duplicated(keep="first")]
    return frame.sort_index()


def nc(file, origin, columns=None):
    """Read the configured subdaily NetCDF data."""
    from local_data import require_local

    require_local(file)
    from netcdf_time import decode
    from netcdf_io import dataset as netcdf_dataset

    with netcdf_dataset(file, buffered=columns is not None) as dataset:
        clock = dataset.variables.get("time")
        if clock is None:
            raise ValueError(f"NetCDF time coordinate missing: {file}")
        index = decode(
            clock[:],
            getattr(clock, "units", origin),
            getattr(clock, "calendar", "standard"),
        ).floor("h")
        columns = {
            key: np.ma.filled(v[:].astype(float), np.nan).reshape(-1)
            for key, v in dataset.variables.items()
            if (columns is None or key in columns)
            and np.issubdtype(v.dtype, np.number)
            and v.size == len(index)
        }
    return mask_missing_values(pd.DataFrame(columns, index=index))


def daily_join(frame, daily, field):
    """Join daily source values onto the subdaily reference grid."""
    values = pd.to_numeric(daily[field], errors="coerce").replace(
        [-999, -9999], np.nan
    )
    return values.reindex(frame.index.normalize()).to_numpy(float)


def load(root, gauge, region, sage, pet, precipitation, temperature):
    """Load the selected data or metadata resource."""
    from lamah_dataset import reader_region

    region = reader_region(region)
    root = Path(root)
    gauge = str(gauge)
    folder, pattern = locations(root, region)
    if region in ("CAMELS_CA", "CAMELS_SPAT"):
        from camels_spat_io import inputs

        file, flow_file = inputs(root, region, gauge, 24)
    else:
        file = one(folder, "**/" + pattern.replace("{gauge}", gauge))
    profile = table_profile(sage, region)
    meta = _metadata(profile, root, gauge, file)
    if region in ("CAMELS_CA", "CAMELS_SPAT"):
        meteo = nc(
            file,
            "hours since 1980-01-01 09:00:00",
            columns={"RDRS_v2.1_A_PR0_SFC", "pet", "RDRS_v2.1_P_TT_1.5m"},
        )
        full = "CAN_" + gauge if region == "CAMELS_CA" else gauge
        flow = nc(
            flow_file,
            "minutes since 1950-01-01",
            columns={
                "q_obs",
                "q_obs_is_ice_affected",
                "q_obs_is_below_sensor_level",
                "q_obs_is_malfunction_affected",
                "q_obs_is_backwater_affected",
            },
        )
        q = pd.to_numeric(flow["q_obs"], errors="coerce").copy()
        flags = (
            ["q_obs_is_ice_affected", "q_obs_is_below_sensor_level"]
            if full.startswith("CAN_")
            else [
                "q_obs_is_ice_affected",
                "q_obs_is_malfunction_affected",
                "q_obs_is_backwater_affected",
            ]
        )
        for flag in flags:
            if flag in flow:
                q = q.where(flow[flag] == 0)
        frame = pd.DataFrame(
            dict(
                P=meteo["RDRS_v2.1_A_PR0_SFC"] * 3600,
                Ep=meteo["pet"] * 3600,
                T=meteo["RDRS_v2.1_P_TT_1.5m"] - 273.15,
                Q=q,
            )
        )
        volume = True
    elif region in ("LAMAH_CE", "CAMELS_IS"):
        meteo = indexed(
            file,
            ymd=True,
            strict=region != "CAMELS_IS",
            columns=[
                "YYYY",
                "MM",
                "DD",
                "hh",
                "mm",
                "2m_temp",
                "prec",
                "pet",
                "total_et",
                "2m_dp_temp",
                "10m_wind_u",
                "10m_wind_v",
                "surf_net_solar_rad",
                "surf_net_therm_rad",
                "surf_press",
            ],
        )
        flowroot = (
            root
            / "hourly"
            / (
                "timeseries/streamflow"
                if region == "LAMAH_CE"
                else "discharge"
            )
        )
        flow = indexed(
            one(flowroot, "ID_" + gauge + ".csv"),
            ymd=True,
            strict=region != "CAMELS_IS",
            row_sequence=region == "CAMELS_IS",
            columns=["YYYY", "MM", "DD", "hh", "mm", "qobs", "qc_flag"],
        )
        t = meteo["2m_temp"].to_numpy(float)
        ep = np.zeros(len(t))
        if region == "CAMELS_IS" and pet in (4, 5):
            ep = meteo["pet" if pet == 4 else "total_et"].to_numpy(
                float, copy=True
            )
            if pet == 4:
                ep[ep < 0] = np.nan
        elif pet == 4:
            ep = np.abs(meteo["total_et"].to_numpy(float))
        elif pet:
            dew = meteo["2m_dp_temp"].to_numpy(float)
            u = meteo["10m_wind_u"].to_numpy(float)
            v = meteo["10m_wind_v"].to_numpy(float)
            solar = meteo["surf_net_solar_rad"].to_numpy(float)
            thermal = meteo["surf_net_therm_rad"].to_numpy(float)
            pressure = meteo["surf_press"].to_numpy(float)
            tk = t + 273.15
            gamma = 0.000665 * pressure / 1000
            wind = np.hypot(u, v) * 4.87 / np.log(67.8 * 10 - 5.42)
            ea = 0.6108 * np.exp(17.27 * dew / (dew + 237.3))
            es = 0.6108 * np.exp(17.27 * t / (t + 237.3))
            slope = 4098 * es / (t + 237.3) ** 2
            rn = (solar - thermal) * 0.0036
            g = np.where(solar > 0, 0.1, 0.5) * rn
            if pet == 1:
                ep = (
                    0.408 * slope * (rn - g)
                    + gamma * (37 / tk) * wind * np.maximum(es - ea, 0)
                ) / (slope + gamma * (1 + 0.34 * wind))
            elif pet == 2:
                ep = 1.26 * slope / (slope + gamma) * 0.408 * (rn - g)
            elif pet == 3:
                ep = (
                    0.61
                    * slope
                    / (slope + gamma)
                    * 0.408
                    * np.maximum(solar, 0)
                    * 0.0036
                    - 0.12 / 24
                )
            else:
                raise ValueError("Invalid LamaH hourly PET selection")
        q = (
            flow["qobs"].where(flow["qc_flag"] <= 100)
            if region == "CAMELS_IS" and "qc_flag" in flow
            else flow["qobs"]
        )
        frame = pd.DataFrame(
            dict(
                P=meteo["prec"],
                Ep=pd.Series(np.maximum(ep, 0), index=meteo.index),
                T=meteo["2m_temp"],
                Q=q,
            )
        )
        volume = True
    elif region == "CAMELS_GB":
        meteo = indexed(file)
        daily = indexed(
            one(
                root / "daily" / "timeseries",
                "*hydromet*_" + gauge + "_*.csv",
            )
        )
        frame = pd.DataFrame(
            dict(
                P=meteo[
                    ["precipitation_cehgear", "precipitation_gradgb"][
                        precipitation - 1
                    ]
                ],
                Ep=daily_join(meteo, daily, ["pet", "peti"][pet - 1]) / 24,
                T=daily_join(meteo, daily, "temperature"),
                Q=meteo["discharge_spec"],
            ),
            index=meteo.index,
        )
        volume = False
    elif region == "CAMELS_CZ":
        meteo = indexed(file, date="time")
        text = profile.get("text", "")

        def choices(name):
            """Return available regional input choices."""
            return (
                profile["file_options"][name]
                if "file_options" in profile
                else strings(
                    re.search(name + r"\s*=\s*\{(.*?)\}", text, re.S).group(
                        1
                    )
                )
            )

        if pet == 1:
            ep = meteo["potential_evaporation_ERA5_LAND"]
        else:
            daily = indexed(
                one(
                    root / "daily" / "timeseries",
                    "camelscz_" + gauge + ".csv",
                )
            )
            ep = daily_join(meteo, daily, choices("eDaily")[pet - 2]) / 24
        frame = pd.DataFrame(
            dict(
                P=meteo[choices("pHourly")[precipitation - 1]],
                Ep=ep,
                T=meteo[choices("tHourly")[temperature - 1]],
                Q=meteo["streamflow"],
            ),
            index=meteo.index,
        )
        volume = False
    elif region == "CAMELS_NZ":
        base = root / "hourly" / "timeseries"
        meteo = indexed(file, "time")
        t = indexed(
            one(
                base / "temperature",
                "temperature_station_id_" + gauge + ".csv",
            ),
            "time",
        )["temperature"]
        if t.median() > 100:
            t = t - 273.15
        ep = indexed(
            one(base / "pet", "PET_station_id_" + gauge + ".csv"), "time"
        )["PET"]
        q = indexed(
            one(base / "streamflow", "flow_station_id_" + gauge + ".csv"),
            "time",
        )["flow"]
        frame = pd.DataFrame(
            dict(P=meteo["precipitation"], Ep=ep, T=t, Q=q)
        )
        volume = True
    elif region == "CAMELS_US":
        wanted = [
            "date",
            "total_precipitation",
            "potential_evaporation",
            "temperature",
        ]
        if pet != 1:
            wanted += [
                "longwave_radiation",
                "pressure",
                "shortwave_radiation",
                "specific_humidity",
                "wind_u",
                "wind_v",
            ]
        meteo = indexed(file, columns=wanted)
        flow = indexed(
            one(root / "hourly" / "streamflow", gauge + "-usgs-hourly.csv"),
            columns=["date", "QObs(mm/h)", "QObs_CAMELS(mm/h)"],
        )
        q = (
            pd.to_numeric(flow["QObs(mm/h)"], errors="coerce").where(
                lambda v: v >= 0
            )
            if "QObs(mm/h)" in flow
            else pd.Series(np.nan, index=flow.index)
        )
        if "QObs_CAMELS(mm/h)" in flow:
            q = q.fillna(
                pd.to_numeric(
                    flow["QObs_CAMELS(mm/h)"], errors="coerce"
                ).where(lambda v: v >= 0)
            )
        frame = pd.DataFrame(
            dict(
                P=meteo["total_precipitation"],
                Ep=meteo["potential_evaporation"],
                T=meteo["temperature"],
                Q=q,
            )
        )
        volume = False
        if pet != 1:
            frame["Ep"] = pd.Series(
                hourly_pet(
                    *[
                        meteo[key].to_numpy(float)
                        for key in (
                            "longwave_radiation",
                            "pressure",
                            "shortwave_radiation",
                            "specific_humidity",
                            "temperature",
                            "wind_u",
                            "wind_v",
                        )
                    ],
                    pet,
                ),
                index=meteo.index,
            )
    else:
        raise ValueError(region)
    frame = frame.apply(pd.to_numeric, errors="coerce")
    frame = mask_missing_values(frame)
    frame["Ep"] = frame["Ep"].clip(lower=0)
    if volume:
        if meta.get("area", 0) <= 0:
            raise ValueError(
                "Positive basin area required for hourly discharge"
            )
        frame["Q"] = frame["Q"] * 3.6 / meta["area"]
    frame["Q"] = frame["Q"].where(frame["Q"] >= 0)
    return frame.reset_index(names="date"), file
