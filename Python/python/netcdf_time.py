# ==========================================================================
# SAGEhydrology | Python source
# File: python/netcdf_time.py
#
# Avoid per-timestamp Python objects for integral modern CF clocks.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Avoid per-timestamp Python objects for integral modern CF clocks."""

import re
import numpy as np
import pandas as pd


def decode(values, units, calendar="standard"):
    """Match netCDF4 dates, accelerating only exact integral coordinates."""
    match = re.fullmatch(
        r"(seconds?|minutes?|hours?|days?) since (.+)", str(units).strip()
    )
    raw = np.ma.asarray(values)
    if (
        match
        and calendar in ("standard", "gregorian", "proleptic_gregorian")
        and not np.ma.getmaskarray(raw).any()
    ):
        numeric = np.asarray(raw)
        if (
            np.isfinite(numeric).all()
            and np.equal(numeric, np.trunc(numeric)).all()
        ):
            try:
                epoch = pd.Timestamp(match[2])
                if epoch.tzinfo is not None:
                    epoch = epoch.tz_convert("UTC").tz_localize(None)
                unit = {"s": "s", "m": "m", "h": "h", "d": "D"}[match[1][0]]
                if epoch.nanosecond:
                    raise ValueError(
                        "Use CF rounding for submicrosecond epochs"
                    )
                index = pd.DatetimeIndex(
                    epoch
                    + pd.to_timedelta(numeric.astype("int64"), unit=unit)
                )
                # The mixed Gregorian calendar differs before 1582.
                if epoch >= pd.Timestamp("1582-10-15") and (
                    not len(index)
                    or index.min() >= pd.Timestamp("1582-10-15")
                ):
                    return index
            except (ValueError, OverflowError, TypeError):
                pass
    from netCDF4 import num2date

    dates = num2date(values, units, calendar=calendar)
    return pd.DatetimeIndex([str(value) for value in dates])
