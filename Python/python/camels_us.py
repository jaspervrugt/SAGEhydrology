# ==========================================================================
# SAGEhydrology | Python source
# File: python/camels_us.py
#
# Read raw CAMELS-US daily forcing/discharge and reproduce SAGE daily
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Read raw CAMELS-US daily forcing/discharge and reproduce SAGE daily
PET.
"""

from datetime import date
from pathlib import Path
import numpy as np


def _daily_dates(ymd):
    """Build publisher calendar dates without one Python string per row."""
    values = np.asarray(ymd).astype(np.int64)
    year, month, day = values.T
    if np.any((month < 1) | (month > 12) | (day < 1) | (day > 31)):
        raise ValueError("Invalid CAMELS calendar date")
    months = (year - 1970).astype("datetime64[Y]").astype(
        "datetime64[M]"
    ) + (month - 1).astype("timedelta64[M]")
    dates = months.astype("datetime64[D]") + (day - 1).astype(
        "timedelta64[D]"
    )
    if np.any(dates.astype("datetime64[M]") != months):
        raise ValueError("Invalid CAMELS calendar date")
    return dates


def daily_pet(
    tmax,
    tmin,
    radiation,
    vapor,
    doy,
    latitude,
    elevation,
    method=1,
    *,
    wind2=2,
):
    """SAGE et0_fao56_daily equations; wind speed defaults to 2 m/s."""
    mean = 0.5 * (tmax + tmin)
    rs = radiation * 0.0864
    ea = vapor / 1000
    es = 0.5 * (
        0.6108 * np.exp(17.27 * tmax / (tmax + 237.3))
        + 0.6108 * np.exp(17.27 * tmin / (tmin + 237.3))
    )
    slope = (
        4098
        * 0.6108
        * np.exp(17.27 * mean / (mean + 237.3))
        / (mean + 237.3) ** 2
    )
    pressure = 101.3 * ((293 - 0.0065 * elevation) / 293) ** 5.26
    gamma = 0.000665 * pressure
    phi = np.deg2rad(latitude)
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365)
    delta = 0.409 * np.sin(2 * np.pi * doy / 365 - 1.39)
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(delta), -1, 1))
    ra = (
        24
        * 60
        / np.pi
        * 0.0820
        * dr
        * (
            ws * np.sin(phi) * np.sin(delta)
            + np.cos(phi) * np.cos(delta) * np.sin(ws)
        )
    )
    rso = (0.75 + 2e-5 * elevation) * ra
    rn = 0.77 * rs - 4.903e-9 * (
        (tmax + 273.16) ** 4 + (tmin + 273.16) ** 4
    ) / 2 * (0.34 - 0.14 * np.sqrt(np.maximum(ea, 0))) * (
        1.35 * np.clip(rs / np.maximum(rso, np.finfo(float).eps), 0, 1)
        - 0.35
    )
    if method == 1:
        pet = (
            0.408 * slope * rn
            + gamma * (900 / (mean + 273)) * wind2 * (es - ea)
        ) / (slope + gamma * (1 + 0.34 * wind2))
    elif method == 2:
        pet = 1.26 * slope / (slope + gamma) * 0.408 * rn
    elif method == 3:
        pet = 0.61 * slope / (slope + gamma) * 0.408 * rs - 0.12
    else:
        raise ValueError(
            "PET method must be 1 (PM), 2 (PT), or 3 (Makkink)"
        )
    return np.maximum(pet, 0)


def _unique(root, pattern):
    """Return unique values while preserving their original order."""
    matches = list(root.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one file matching {root / pattern}; found"
            f" {len(matches)}"
        )
    return matches[0]


def forcing_files(root, product, gauge="*"):
    """Official Daymet archives label their forcing files 'cida'."""
    if product not in ("daymet", "maurer", "nldas"):
        raise ValueError("Invalid CAMELS-US forcing product")
    folder = Path(root) / "daily/v1p2/forcing" / product
    labels = ("cida", "daymet") if product == "daymet" else (product,)
    return sorted(
        {
            p
            for label in labels
            for p in folder.glob(f"*/{gauge}_lump_{label}_forcing_leap.txt")
        }
    )


def read_daily(
    root,
    gauge,
    *,
    product="nldas",
    pet=1,
    train_start="1999-10-01",
    train_end="2008-09-30",
    eval_start="1989-10-01",
    eval_end="1999-09-30",
    spinup=365,
    allow_invalid=False,
    meteo_only=False,
    _store_single=True,
):
    """Read CAMELS-US daily inputs through the shared schema engine."""
    from regional_data import read_region

    return read_region(
        root,
        gauge,
        region="CAMELS_US",
        resolution=1,
        product=product,
        pet=pet,
        train_start=train_start,
        train_end=train_end,
        eval_start=eval_start,
        eval_end=eval_end,
        spinup=spinup,
        allow_invalid=allow_invalid,
        meteo_only=meteo_only,
    )


def read_discharge(root, gauge, data):
    """Read discharge using the same schema and conversion as MATLAB."""
    from hydro_schema import configuration, identifier
    from hydro_io import read_raw
    from hydro_derivations import variable
    from local_data import capture_sources

    _, schema, paths, options = configuration(
        None, "CAMELS_US", 1, "nldas", 1, 1, 1
    )
    with capture_sources() as sources:
        raw, used = read_raw(
            Path(root) / paths["discharge"],
            identifier(gauge, schema["id"]),
            schema,
            data["dates"],
            options,
        )
    q = variable(
        raw["Q"],
        schema["variables"]["Q"],
        [np.nan, np.nan, float(data["area"])],
        86400,
    )
    valid = np.isfinite(q) & (q >= 0)
    q[~valid] = np.nan
    q = q.astype(np.float32).astype(float)
    from timeseries_storage import store_timeseries

    return store_timeseries(
        dict(
            data,
            observed=q,
            matlab_observed=q.copy(),
            train=data["train"] & valid,
            evaluation=data["evaluation"] & valid,
            streamflow_file=str(used[0]),
            _source_files=sorted(
                set(data.get("_source_files", [])) | sources
            ),
        )
    )
