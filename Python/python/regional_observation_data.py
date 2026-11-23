# ==========================================================================
# SAGEhydrology | Python source
# File: python/regional_observation_data.py
#
# Extract auxiliary observations without reopening a combined source file.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Extract auxiliary observations without reopening a combined source file.
"""

import numpy as np
import pandas as pd


def arrays(frame, region):
    """Preserve the existing SWE and soil-layer storage conversions."""
    if region not in ("HYD_RESPONSES", "MACH_US"):
        return {}
    columns = {
        "SWE": "SWECLQMD_mean" if region == "HYD_RESPONSES" else "SWE"
    }
    frame = frame.copy(deep=False)
    layers = ["swvl1_mean_mean", "swvl2_mean_mean", "swvl3_mean_mean"]
    if region == "HYD_RESPONSES" and all(k in frame for k in layers):
        values = (
            frame[layers]
            .apply(pd.to_numeric, errors="coerce")
            .to_numpy(float, copy=True)
        )
        values[np.isin(values, [-999, -99, -99.9, -99.99])] = np.nan
        frame["_SM_mm"] = values @ np.array([70.0, 210.0, 720.0])
        columns["SM"] = "_SM_mm"
    result = {}
    for name in ("SWE", "SM"):
        column = columns.get(name)
        if column is None or column not in frame:
            result["_source_" + name] = np.full(len(frame), np.nan)
            continue
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(
            float, copy=True
        )
        values[np.isin(values, [-999, -9999, -99, -99.9, -99.99])] = np.nan
        values[np.isfinite(values) & (values < 0)] = np.nan
        result["_source_" + name] = values
    return result
