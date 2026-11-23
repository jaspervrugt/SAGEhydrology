# ==========================================================================
# SAGEhydrology | Python source
# File: python/hydro_derivations.py
#
# Schema-defined hydrologic transformations in MATLAB operation order.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Schema-defined hydrologic transformations in MATLAB operation order."""

import numpy as np
import pandas as pd


def variable(value, spec, aux, step):
    """Resolve a variable required by a derived data column."""
    value = np.asarray(value, dtype=float).copy()
    value[np.isin(value, spec.get("missing_values", [-9999, -999]))] = (
        np.nan
    )
    value[
        (value < float(spec.get("valid_min", -np.inf)))
        | (value > float(spec.get("valid_max", np.inf)))
    ] = np.nan
    value = value * spec.get("scale", 1) + spec.get("offset", 0)
    finite = np.isfinite(value)
    if "clip_min" in spec:
        value[finite & (value < float(spec["clip_min"]))] = float(
            spec["clip_min"]
        )
    if "clip_max" in spec:
        value[finite & (value > float(spec["clip_max"]))] = float(
            spec["clip_max"]
        )
    if spec.get("fill_missing") is not None and spec["fill_missing"] != []:
        value[~np.isfinite(value)] = spec["fill_missing"]
    if spec.get("area_normalize"):
        if not np.isfinite(aux[2]) or aux[2] <= 0:
            raise ValueError("Positive basin area required for discharge")
        value = value * step * 1000 / aux[2]
    transform = spec.get("transform")
    if transform:
        if (
            transform.get("matlab_function", "").replace(" ", "")
            != "@(x,context)abs(x)"
        ):
            raise ValueError(f"Unsupported schema transform: {transform}")
        value = np.abs(value)
    return value


def derive(raw, schema, dates, aux):
    """Compute the requested meteorological derived variable."""
    variables = schema["variables"]
    for field, spec in variables.items():
        if spec.get("derive", "").lower() == "constant":
            raw[field] = np.full(len(dates), float(spec["value"]))
    for field, spec in variables.items():
        if field in raw or not spec.get("derive"):
            continue
        method = spec["derive"].lower()
        if method == "mean_tmin_tmax":
            value = 0.5 * (raw["Tmin"] + raw["Tmax"])
        elif method == "weighted_sum":
            value = np.zeros(len(dates))
            for name, weight in zip(spec["inputs"], spec["weights"]):
                value = value + weight * raw[name]
        elif method == "temperature_auto_celsius":
            value = raw[spec.get("input", "TemperatureRaw")].copy()
            finite = value[np.isfinite(value)]
            if len(finite) and np.median(finite) > 100:
                value -= 273.15
        elif method in ("oudin", "hargreaves"):
            from regional_data import extraterrestrial

            lat = aux[0]
            if not np.isfinite(lat):
                lat = (
                    60
                    if method == "oudin"
                    else spec.get("latitude_fallback", np.nan)
                )
            radiation = extraterrestrial(dates, lat)
            t = raw["T"]
            if method == "oudin":
                if spec.get("input_precision") != "double":
                    t = t.astype(np.float32).astype(float)
                value = (
                    spec.get("coefficient", 0.408 / 100)
                    * radiation
                    * np.maximum(t + 5, 0)
                )
                value = value.astype(np.float32).astype(float)
                value[~np.isfinite(value) | (value < 0)] = 0
            else:
                value = (
                    0.0023
                    * (radiation / 2.45)
                    * (t + 17.8)
                    * np.sqrt(np.maximum(0, raw["Tmax"] - raw["Tmin"]))
                )
                value[value < 0] = 0
                if "fill_missing" in spec:
                    value[~np.isfinite(value)] = spec["fill_missing"]
        elif method == "vapor_pressure_dewpoint":
            t = raw[spec.get("input", "Tdew")]
            value = 1000 * 0.6108 * np.exp(17.27 * t / (t + 237.3))
        elif method == "vapor_pressure_relative_humidity":
            t = raw[spec.get("temperature_input", "T")]
            rh = np.clip(
                raw[spec.get("humidity_input", "RelativeHumidity")], 0, 100
            )
            value = (
                1000 * (rh / 100) * 0.6108 * np.exp(17.27 * t / (t + 237.3))
            )
        elif method == "mixing_ratio_to_specific_humidity":
            ratio = raw[spec.get("input", "MixingRatio")] * spec.get(
                "input_scale", 1
            )
            value = ratio / (1 + ratio)
        elif method == "wind_uv_to_2m":
            value = np.hypot(raw["WindU10"], raw["WindV10"]) * (
                4.87 / np.log(67.8 * 10 - 5.42)
            )
        elif method == "fao56_daily":
            from camels_us import daily_pet

            radiation = raw["Radiation"] * spec.get("radiation_scale", 1)
            if variables["Radiation"].get("averaging") == "daylight":
                daylength = raw["DayLength"]
                if np.any(
                    ~np.isfinite(daylength)
                    | (daylength < 0)
                    | (daylength > 86400)
                ):
                    raise ValueError("Invalid daylight duration")
                radiation = radiation * daylength / 86400
            if not np.all(np.isfinite(aux[:2])):
                raise ValueError(
                    "Daily PET requires latitude and elevation"
                )
            value = daily_pet(
                raw["Tmax"],
                raw["Tmin"],
                radiation,
                raw["VaporPressure"],
                pd.DatetimeIndex(dates).dayofyear.to_numpy(),
                aux[0],
                aux[1],
                spec["method"],
                wind2=raw["Wind2"],
            )
        elif method == "lamah_hourly":
            value = lamah_pet(raw, spec["method"])
        elif method in ("fao56_hourly", "fao56_hourly_estimated_longwave"):
            from regional_subdaily import hourly_pet

            longwave = raw.get("Longwave")
            if method.endswith("estimated_longwave"):
                tk = raw["T"] + 273.15
                ea = raw["VaporPressure"] / 1000
                clear = np.clip(
                    1.24 * np.maximum(ea / tk, 0) ** (1 / 7), 0.55, 1
                )
                cloud = np.clip(raw["CloudCover"] / 100, 0, 1)
                sky = np.clip(clear * (1 + 0.22 * cloud**2), 0.55, 1)
                longwave = sky * 5.670374419e-8 * tk**4
            value = hourly_pet(
                longwave,
                raw["Pressure"] * spec.get("pressure_scale", 1),
                raw["Radiation"] * spec.get("radiation_scale", 1),
                raw["SpecificHumidity"],
                raw["T"],
                raw["WindU10"],
                raw["WindV10"],
                spec["method"],
            )
        else:
            raise ValueError(f"Unknown hydro derivation: {method}")
        raw[field] = value
    return raw


def lamah_pet(raw, method):
    """Compute PET using the LamaH reader convention."""
    t = raw["T"]
    dew = raw["Tdew"]
    wind = (
        np.hypot(raw["WindU10"], raw["WindV10"])
        * 4.87
        / np.log(67.8 * 10 - 5.42)
    )
    rn = (raw["NetSolar"] - raw["NetThermal"]) * 0.0036
    rs = np.maximum(raw["NetSolar"], 0) * 0.0036
    gamma = 0.000665 * raw["Pressure"] / 1000
    es = 0.6108 * np.exp(17.27 * t / (t + 237.3))
    ea = 0.6108 * np.exp(17.27 * dew / (dew + 237.3))
    slope = 4098 * es / (t + 237.3) ** 2
    g = np.where(rs > 0, 0.1, 0.5) * rn
    if method == 1:
        value = (
            0.408 * slope * (rn - g)
            + gamma * (37 / (t + 273.15)) * wind * np.maximum(es - ea, 0)
        ) / (slope + gamma * (1 + 0.34 * wind))
    elif method == 2:
        value = 1.26 * slope / (slope + gamma) * 0.408 * (rn - g)
    elif method == 3:
        value = 0.61 * slope / (slope + gamma) * 0.408 * rs - 0.12 / 24
    else:
        raise ValueError("Unknown hourly PET method")
    return np.maximum(value, 0)
