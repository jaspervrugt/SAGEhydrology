# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_water_year_blocks.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_water_year_blocks."""
import numpy as np

def water_year_blocks(dates, c):
    """Construct complete calendar-based water-year blocks."""
    months = dates.astype('datetime64[M]').astype(int) % 12 + 1
    years = dates.astype('datetime64[Y]').astype(int) + 1970
    month = int(c.get('water_year_start_month', 10))
    if not 1 <= month <= 12:
        raise ValueError('Water-year start month must be between 1 and 12.')
    width = float(c.get('block_years', 1))
    if not np.isfinite(width) or width != int(width) or width < 1:
        raise ValueError('Block size in water years must be a positive whole number.')
    labels = years - (months < month)
    first = int(labels[0])
    start = np.datetime64(f'{first:04d}-{month:02d}-01')
    final = int(labels[-1]) + 1
    end = np.datetime64(f'{final:04d}-{month:02d}-01')
    if dates[0].astype('datetime64[D]') != start or dates[-1].astype('datetime64[D]') != end - np.timedelta64(1, 'D'):
        raise ValueError(f'Water-year blocks require complete years: start on {month:02d}-01 and end on the preceding day of a later year.')
    groups = (labels - first) // int(width)
    return [np.flatnonzero(groups == group) for group in np.unique(groups)]
