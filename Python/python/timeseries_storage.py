# ==========================================================================
# SAGEhydrology | Python source
# File: python/timeseries_storage.py
#
# MATLAB-compatible storage for canonical input time series.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""MATLAB-compatible storage for canonical input time series.

Derivation and unit conversion occur before this boundary. Solvers and loss
functions promote inputs to double for computation; dates, masks, attributes
and parameters retain their existing types.
"""

import numpy as np


def store_timeseries(data):
    """Store basin time series with the configured retention policy."""
    result = dict(data)
    for key in (
        "forcing",
        "observed",
        "observed_SWE",
        "observed_SM",
        "_source_SWE",
        "_source_SM",
    ):
        if key in result:
            result[key] = np.asarray(result[key], dtype=np.float32)
    return result
