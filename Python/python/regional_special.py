# ==========================================================================
# SAGEhydrology | Python source
# File: python/regional_special.py
#
# Adapters for wide, multi-file and NetCDF regional datasets.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Adapters for wide, multi-file and NetCDF regional datasets."""

from pathlib import Path
from functools import lru_cache
import re
import numpy as np
import pandas as pd

SPECIAL = {
    "CAMELS_AU",
    "CAMELS_BR",
    "CAMELS_CA",
    "CAMELS_SPAT",
    "CAMELS_CH",
    "CAMELS_CL",
    "CAMELS_IND",
    "CAMELS_NZ",
    "HYDRO_CIS",
    "BULL_ES",
}


def _one(root, pattern):
    """Read one source value or column."""
    files = list(Path(root).glob(pattern))
    if len(files) != 1:
        raise FileNotFoundError(
            f"Expected one {Path(root)/pattern}; found {len(files)}"
        )
    return files[0]


def _indexed(path, date=None, ymd=None, skiprows=0):
    """Read a gauge-indexed source table."""
    from regional_data import _table

    frame = _table(path, skiprows)
    if ymd:
        dates = pd.to_datetime(
            dict(year=frame[ymd[0]], month=frame[ymd[1]], day=frame[ymd[2]])
        )
    else:
        dates = pd.to_datetime(frame[date], errors="raise")
    if dates.duplicated().any():
        raise ValueError(f"Duplicate timestamps in {path}")
    frame.index = pd.DatetimeIndex(dates).normalize()
    return frame


def _wide(path, gauge, date=None, ymd=None, sequential=False):
    """Read a wide-format basin source table."""
    path = Path(path)
    if (
        not path.exists()
        and (path.parent / "streamflow" / path.name).exists()
    ):
        path = path.parent / "streamflow" / path.name
    from local_data import require_local

    require_local(path)
    from wide_data_cache import series
    from local_data import selected_gauges

    # Match MATLAB's wide-input boundary before derived calculations.
    return (
        series(
            path,
            str(gauge),
            tuple(ymd or ()),
            date=date,
            sequential=sequential,
            gauges=selected_gauges(),
        )
        .astype(np.float32)
        .astype(float)
    )


@lru_cache(maxsize=64)
def _wide_cached(path, gauge, date, ymd, sequential, mtime, size):
    # Wide regional files contain hundreds/thousands of catchment columns.
    # Parse only the requested station and clock, rather than allocating an
    # entire regional matrix for every basin or forcing option.
    """Read a wide-format source through the shared cache."""
    with Path(path).open(encoding="utf-8-sig") as stream:
        first = stream.readline()
    sep = (
        ";"
        if first.count(";") > first.count(",")
        else "," if "," in first else "\t" if "\t" in first else r"\s+"
    )
    header = pd.read_csv(
        path, sep=sep, nrows=0, encoding="utf-8-sig"
    ).columns
    column = next(
        (
            c
            for c in header
            if c.strip() == gauge
            or c.strip().lstrip("0") == gauge.lstrip("0")
        ),
        None,
    )
    if column is None:
        raise ValueError(f"Basin column {gauge} absent in {path}")
    clock = list(ymd) if ymd else [date]
    wanted = {name.strip() for name in (column, *clock)}
    frame = pd.read_csv(
        path,
        sep=sep,
        usecols=lambda name: name.strip() in wanted,
        encoding="utf-8-sig",
    )
    frame.columns = [c.strip() for c in frame]
    if ymd:
        dates = pd.to_datetime(
            dict(year=frame[ymd[0]], month=frame[ymd[1]], day=frame[ymd[2]])
        )
    elif sequential:
        dates = pd.date_range(
            pd.to_datetime(frame[date].iloc[0]),
            periods=len(frame),
            freq="D",
        )
    else:
        dates = pd.to_datetime(frame[date], errors="raise")
    dates = pd.DatetimeIndex(dates).normalize()
    if dates.duplicated().any():
        raise ValueError(f"Duplicate timestamps in {path}")
    return pd.Series(
        pd.to_numeric(frame[column.strip()], errors="coerce").to_numpy(),
        index=dates,
    )


def _nc(path, time="time", origin=None, columns=None):
    """Read the selected NetCDF source."""
    from local_data import require_local

    require_local(path)
    from netcdf_time import decode
    from netcdf_io import dataset as netcdf_dataset

    with netcdf_dataset(path, buffered=columns is not None) as dataset:
        clock = dataset.variables.get(time)
        units = getattr(clock, "units", origin)
        if not units:
            raise ValueError(f"NetCDF time units absent in {path}")
        if clock is None:
            size = max(len(d) for d in dataset.dimensions.values())
            coordinate = np.arange(size)
        else:
            coordinate = clock[:]
        dates = decode(
            coordinate, units, getattr(clock, "calendar", "standard")
        ).normalize()
        if dates.duplicated().any():
            raise ValueError(f"Duplicate daily timestamps in {path}")
        selected = columns
        columns = {}
        for key, var in dataset.variables.items():
            if selected is not None and key not in selected:
                continue
            values = var[:]
            if values.size == len(dates) and np.issubdtype(
                values.dtype, np.number
            ):
                columns[key] = np.ma.filled(
                    values.astype(float), np.nan
                ).reshape(-1)
    return pd.DataFrame(columns, index=dates)


def inventory(root, region, profile):
    """Inspect the files available for the regional reader."""
    base = Path(root) / "daily"
    patterns = {
        "CAMELS_BR": ("precipitation", "{gauge}_precipitation.txt"),
        "CAMELS_CH": (
            "timeseries/observation_based",
            "CAMELS_CH_obs_based_{gauge}.csv",
        ),
        "CAMELS_NZ": (
            "timeseries/precipitation",
            "daily_precipitation_station_id_{gauge}.csv",
        ),
        "CAMELS_IND": ("catchment_mean_forcings", "{gauge}.csv"),
        "HYDRO_CIS": ("timeseries", "{gauge}.nc"),
        "BULL_ES": ("timeseries/AEMET", "AEMET_{gauge}.csv"),
    }
    if region in ("CAMELS_CA", "CAMELS_SPAT"):
        folder = base / "forcing"
        pattern = (
            "CAN_{gauge}_daymet_lumped.nc"
            if region == "CAMELS_CA"
            else "{gauge}_daymet_lumped.nc"
        )
        files = folder.glob("**/" + pattern.replace("{gauge}", "*"))
    elif region in ("CAMELS_AU", "CAMELS_CL"):
        path = (
            base
            / "streamflow"
            / (
                "streamflow_mmd.csv"
                if region == "CAMELS_AU"
                else "3_CAMELScl_streamflow_mm.txt"
            )
        )
        from regional_data import _table

        # Read just the header of potentially very large wide files.
        with path.open(encoding="utf-8-sig") as handle:
            line = handle.readline()
        delimiter = "," if "," in line else "\t"
        columns = pd.read_csv(path, sep=delimiter, nrows=0).columns
        return [
            (c, "", "")
            for c in columns
            if c not in ("year", "month", "day", "gauge_id", "date")
        ]
    else:
        folder, pattern = patterns[region]
        files = (base / folder).glob(pattern.replace("{gauge}", "*"))
    regex = re.escape(pattern).replace(
        re.escape("{gauge}"), r"(?P<gauge>.+)"
    )
    return sorted(
        {
            (re.fullmatch(regex, file.name).group("gauge"), "", "")
            for file in files
            if re.fullmatch(regex, file.name)
        }
    )


def load(
    root, gauge, region, profile, pet=1, temperature=1, precipitation=1
):
    """Load the selected data or metadata resource."""
    from regional_data import strings, _metadata

    base = Path(root) / "daily"
    gauge = str(gauge)
    text = profile.get("text", "")
    units = "mm/day"
    metadata = {}
    files = []

    def column(frame, name):
        """Extract a numeric column from the regional source."""
        if name not in frame:
            raise ValueError(f"Missing variable {name}")
        values = pd.to_numeric(frame[name], errors="coerce")
        return values.mask(values.isin([-999, -9999]))

    def read(pattern, date="date", ymd=None):
        """Read the configured basin data and reference time window."""
        file = _one(base, pattern.replace("{gauge}", gauge))
        files.append(file)
        return _indexed(file, date, ymd)

    if region == "CAMELS_AU":
        paths = {}
        for key, selection in [
            ("P", precipitation),
            ("Ep", pet),
            ("Tmin", temperature),
            ("Tmax", temperature),
        ]:
            options = profile.get("file_options", {}).get(key)
            if options is None:
                options = strings(
                    re.search(
                        r"S\.variables\."
                        + key
                        + r"\.files\s*=\s*\{(.*?)\}",
                        text,
                        re.S,
                    ).group(1)
                )
            paths[key] = options[int(selection) - 1]
        series = {
            k: _wide(base / path, gauge, ymd=["year", "month", "day"])
            for k, path in paths.items()
        }
        series["Q"] = _wide(
            base / "streamflow_mmd.csv", gauge, ymd=["year", "month", "day"]
        )
        frame = pd.DataFrame(series)
        frame[frame <= -99] = np.nan
        frame["T"] = 0.5 * (
            frame["Tmin"].astype(float) + frame["Tmax"].astype(float)
        )
        files = [base / paths["P"]]
    elif region == "CAMELS_CL":

        def array(name):
            """Convert regional input values to an array."""
            return (
                profile["file_options"][name]
                if "file_options" in profile
                else strings(
                    re.search(name + r"\s*=\s*\{(.*?)\}", text, re.S).group(
                        1
                    )
                )
            )

        p = array("pFiles")[int(precipitation) - 1]
        ep = array("epFiles")[int(pet) - 1]
        series = {
            "P": _wide(base / p, gauge, date="gauge_id"),
            "Ep": _wide(
                base / ep, gauge, date="gauge_id", sequential=pet == 2
            ),
            "Q": _wide(
                base / "3_CAMELScl_streamflow_mm.txt",
                gauge,
                date="gauge_id",
            ),
        }
        if temperature == 2:
            series["T"] = 0.5 * (
                _wide(
                    base / "temperature/8_CAMELScl_tmin_cr2met.txt",
                    gauge,
                    date="gauge_id",
                )
                + _wide(
                    base / "temperature/9_CAMELScl_tmax_cr2met.txt",
                    gauge,
                    date="gauge_id",
                )
            )
        else:
            series["T"] = _wide(
                base / "temperature/10_CAMELScl_tmean_cr2met.txt",
                gauge,
                date="gauge_id",
            )
        frame = pd.DataFrame(series)
        frame[frame <= -99] = np.nan
        files = [base / p]
    elif region == "CAMELS_BR":
        p = read(
            "precipitation/{gauge}_precipitation.txt",
            ymd=["year", "month", "day"],
        )
        ep = read(
            "pet/{gauge}_potential_evapotransp.txt",
            ymd=["year", "month", "day"],
        )
        t = read(
            "temperature/{gauge}_temperature.txt",
            ymd=["year", "month", "day"],
        )
        q = read(
            "streamflow/{gauge}_streamflow.txt",
            ymd=["year", "month", "day"],
        )
        pcols = ["p_brdwgd", "p_chirps", "p_cpc", "p_era5land", "p_mswep"]
        ecols = ["pet_gleam", "pet_era5land"]
        # MATLAB's product selectors reserve zero for the finite-product
        # mean.
        # Mask each product's missing values before averaging.
        pvalue = (
            pd.DataFrame({c: column(p, c) for c in pcols}).mean(axis=1)
            if precipitation == 0
            else column(p, pcols[precipitation - 1])
        )
        evalue = (
            pd.DataFrame({c: column(ep, c) for c in ecols}).mean(axis=1)
            if pet == 0
            else column(ep, ecols[pet - 1])
        )
        tcols = [
            ["tmean_era5land"],
            ["tmax_cpc", "tmin_cpc"],
            ["tmax_era5land", "tmin_era5land"],
            ["tmax_brdwgd", "tmin_brdwgd"],
        ][temperature - 1]
        qvalue = column(q, "streamflow_m3s")
        for quality in ("qual_control_by_ana", "qual_flag"):
            if quality in q:
                qvalue = qvalue.where(column(q, quality) > 0)
        frame = pd.DataFrame(
            dict(
                P=pvalue,
                Ep=evalue,
                T=t[tcols].mean(axis=1, skipna=False),
                Q=qvalue,
            )
        )
        units = "m3/s"
    elif region == "CAMELS_CH":
        obs = read(
            "timeseries/observation_based/CAMELS_CH_obs_based_{gauge}.csv"
        )
        sim = read(
            "timeseries/simulation_based/CAMELS_CH_sim_based_{gauge}.csv"
        )
        frame = pd.DataFrame(
            dict(
                P=column(obs, "precipitation(mm/d)"),
                T=column(obs, "temperature_mean(degC)"),
                Ep=column(sim, "pet_sim(mm/d)"),
                Q=column(obs, "discharge_spec(mm/d)"),
            )
        )
        if precipitation == 2:
            frame["P"] = column(sim, "precipitation_sim(mm/d)")
        if temperature == 2:
            frame["T"] = 0.5 * (
                column(obs, "temperature_min(degC)")
                + column(obs, "temperature_max(degC)")
            )
    elif region == "CAMELS_IND":
        meteo = read(
            "catchment_mean_forcings/{gauge}.csv",
            ymd=["year", "month", "day"],
        )
        q = _wide(
            base / "streamflow_timeseries" / "streamflow_observed.csv",
            gauge,
            ymd=["year", "month", "day"],
        )
        t = column(
            meteo,
            ["tavg(C)", "tavg(C)", "tmax(C)", "tmin(C)"][temperature - 1],
        )
        if temperature == 2:
            t = 0.5 * (column(meteo, "tmin(C)") + column(meteo, "tmax(C)"))
        frame = pd.DataFrame(
            dict(
                P=column(meteo, "prcp(mm/day)"),
                Ep=column(
                    meteo, ["pet(mm/day)", "pet_gleam(mm/day)"][pet - 1]
                ),
                T=t,
                Q=q,
            )
        )
        units = "m3/s"
    elif region == "CAMELS_NZ":
        p = read(
            "timeseries/precipitation/daily_precipitation_station_i"
            "d_{gauge}.csv",
            "time",
        )
        t = read(
            "timeseries/temperature/daily_temperature_station_id_{g"
            "auge}.csv",
            "time",
        )
        q = read(
            "timeseries/streamflow/daily_flow_station_id_{gauge}.csv",
            "time",
        )
        temperature_values = column(t, "temperature")
        if temperature_values.median() > 100:
            temperature_values = temperature_values - 273.15
        frame = pd.DataFrame(
            dict(
                P=column(p, "precipitation"),
                T=temperature_values,
                Q=column(q, "flow"),
            )
        )
        units = "m3/s"
        profile["variables"]["Ep"] = dict(
            derive="oudin",
            coefficient=1 / (100 * 2.45),
            input_precision="double",
        )
    elif region == "BULL_ES":
        products = ["AEMET/AEMET", "ERA5_Land/ERA5_Land", "EMO1_arc/EMO1"]
        tables = {
            choice: read(
                "timeseries/" + products[choice - 1] + "_{gauge}.csv"
            )
            for choice in set((precipitation, temperature, pet))
        }
        q = read("timeseries/streamflow/streamflow_{gauge}.csv")
        frame = pd.DataFrame(
            dict(
                P=column(tables[precipitation], "total_precipitation"),
                Ep=column(tables[pet], "potential_evapotranspiration"),
                T=column(tables[temperature], "temperature_mean"),
                Q=column(q, "streamflow"),
            )
        )
    elif region == "HYDRO_CIS":
        file = _one(base / "timeseries", gauge + ".nc")
        files = [file]
        meteo = _nc(file, "index", "days since 2008-01-01")
        suffix = ["e5l", "e5"][temperature - 1]
        frame = pd.DataFrame(
            dict(
                P=column(
                    meteo,
                    [
                        "prcp_e5l",
                        "prcp_e5",
                        "prcp_mswep",
                        "prcp_gpcp",
                        "prcp_imerg",
                    ][precipitation - 1],
                ),
                T=0.5
                * (
                    column(meteo, "t_min_" + suffix)
                    + column(meteo, "t_max_" + suffix)
                ),
                Q=column(meteo, "q_mm_day"),
            )
        )
        if pet == 2:
            frame["Ep"] = column(meteo, "Ep")
        else:
            profile["variables"]["Ep"] = dict(derive="oudin")
        frame["P"] = frame["P"].where(frame["P"] >= -1e-3).clip(lower=0)
    elif region in ("CAMELS_CA", "CAMELS_SPAT"):
        full = "CAN_" + gauge if region == "CAMELS_CA" else gauge
        from camels_spat_io import inputs

        file, flow_file = inputs(root, region, gauge, 1)
        files = [file]
        meteo = _nc(
            file,
            origin="days since 1980-01-01 12:00:00",
            columns={"prcp", "pet", "tmin", "tmax"},
        )
        flow = _nc(
            flow_file,
            origin="minutes since 1950-01-01",
            columns={
                "q_obs",
                "q_obs_is_ice_affected",
                "q_obs_is_partial_day",
                "q_obs_is_malfunction_affected",
                "q_obs_is_backwater_affected",
            },
        )
        q = column(flow, "q_obs")
        quality = (
            ["q_obs_is_ice_affected", "q_obs_is_partial_day"]
            if full.startswith("CAN_")
            else [
                "q_obs_is_ice_affected",
                "q_obs_is_malfunction_affected",
                "q_obs_is_backwater_affected",
            ]
        )
        for flag in quality:
            if flag in flow:
                q = q.where(column(flow, flag) == 0)
        frame = pd.DataFrame(
            dict(
                P=column(meteo, "prcp"),
                Ep=column(meteo, "pet"),
                T=0.5 * (column(meteo, "tmin") + column(meteo, "tmax")),
                Q=q,
            )
        )
        units = "m3/s"
    else:
        raise ValueError(region)
    for key in ("P", "T", "Q"):
        profile["variables"][key] = dict(
            source=key,
            units=(
                units if key == "Q" else "degC" if key == "T" else "mm/day"
            ),
        )
    if "Ep" in frame:
        profile["variables"]["Ep"] = dict(source="Ep", units="mm/day")
    metadata = _metadata(profile, root, gauge, files[0])
    return (
        frame.reset_index(names="date"),
        frame.index.to_numpy(dtype="datetime64[D]"),
        files[0],
        metadata,
    )
