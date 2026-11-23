# ==========================================================================
# SAGEhydrology | Python source
# File: python/regional_subdaily.py
#
# Subdaily regional readers. Fluxes use mm per model time step.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Subdaily regional readers. Fluxes use mm per model time step."""

from pathlib import Path
from camels_spat_io import inputs as spat_inputs
import numpy as np
import pandas as pd
from regional_data import table_profile, _table, _metadata

SUPPORTED = {
    "CAMELS_DEH",
    "CAMELSH_KR",
    "CAMELSH_US",
    "CAMELS_LUX",
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


def inventory(root, region, resolution):
    """Inspect the files available for the regional reader."""
    from subdaily_profiles import EXTRA, inventory as extra_inventory

    if region in EXTRA:
        return extra_inventory(root, region)
    base = (
        Path(root)
        / ("hourly" if resolution == 24 else "15min")
        / "timeseries"
    )
    patterns = {
        "CAMELS_DEH": "CAMELS_DE_1h_hydromet_timeseries_*.csv",
        "CAMELSH_KR": "*.csv",
        "CAMELSH_US": "*.nc",
        "CAMELS_LUX": (
            "CAMELS_LUX_hydromet_timeseries_"
            f"{('hourly' if resolution == 24 else '15min')}"
            "_ID_*.csv"
        ),
    }
    files = list(base.glob(patterns[region]))
    result = []
    for file in files:
        gauge = file.stem
        if region == "CAMELS_DEH":
            gauge = gauge.removeprefix("CAMELS_DE_1h_hydromet_timeseries_")
        if region == "CAMELS_LUX":
            gauge = gauge.rsplit("_ID_", 1)[1]
        result.append((gauge, "", ""))
    return sorted(result)


def hourly_pet(longwave, pressure, radiation, humidity, t, u, v, method):
    """Construct hourly PET using the regional source convention."""
    tk = t + 273.15
    gamma = 0.000665 * pressure / 1000
    wind = np.hypot(u, v) * 4.87 / np.log(67.8 * 10 - 5.42)
    ea = humidity * pressure / (0.622 + 0.378 * humidity) / 1000
    es = 0.6108 * np.exp(17.27 * t / (t + 237.3))
    slope = 4098 * es / (t + 237.3) ** 2
    rn = (
        0.77 * radiation + longwave - 0.98 * 5.670374419e-8 * tk**4
    ) * 0.0036
    g = np.where(radiation > 0, 0.1, 0.5) * rn
    if method == 1:
        ep = (
            0.408 * slope * (rn - g)
            + gamma * (37 / tk) * wind * np.maximum(es - ea, 0)
        ) / (slope + gamma * (1 + 0.34 * wind))
    elif method == 2:
        ep = 1.26 * slope / (slope + gamma) * 0.408 * (rn - g)
    elif method == 3:
        ep = (
            0.61 * slope / (slope + gamma) * 0.408 * radiation * 0.0036
            - 0.12 / 24
        )
    else:
        raise ValueError("Unknown hourly PET method")
    return np.maximum(ep, 0)


def read(
    root,
    gauge,
    region,
    sage_root,
    resolution,
    pet,
    precipitation,
    temperature,
    period,
    allow_invalid=False,
):
    """Read every subdaily region through the shared schema engine."""
    from regional_data import read_region

    return read_region(
        root,
        gauge,
        region=region,
        sage_root=sage_root,
        resolution=resolution,
        pet=pet,
        precipitation=precipitation,
        temperature=temperature,
        allow_invalid=allow_invalid,
        **period,
    )
