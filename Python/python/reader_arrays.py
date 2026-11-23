# ==========================================================================
# SAGEhydrology | Python source
# File: python/reader_arrays.py
#
# Shared numeric missing-value handling for regional data readers.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Shared numeric missing-value handling for regional data readers."""

import numpy as np


def mask_missing_values(frame):
    """Mask publisher sentinels without hashing every numeric value.

    Preserve DataFrame labels, column dtypes and existing missing values.
    Mixed/string tables retain pandas' original membership semantics.
    """
    values = frame.to_numpy(copy=False)
    if not np.issubdtype(values.dtype, np.number):
        return frame.mask(frame.isin([-999, -9999]))
    missing = (values == -999) | (values == -9999)
    return frame.mask(missing)
