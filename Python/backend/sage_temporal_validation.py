# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_temporal_validation.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_temporal_validation."""
from datetime import date
import math

def validate_design(c):
    """Validate temporal split fields and their consistency."""
    method = c['split']
    if method not in ('manual', 'random', 'deterministic_block', 'random_block', 'deterministic_kfold', 'random_kfold', 'rainfall_block'):
        raise ValueError('Choose a supported train/eval method.')

    def integer(key, minimum, maximum):
        """Validate an integer-valued temporal-design setting."""
        value = float(c[key])
        if not math.isfinite(value) or value != int(value) or (not minimum <= value <= maximum):
            raise ValueError(f"{key.replace('_', ' ').capitalize()} must be a whole number between {minimum} and {maximum}.")
        return int(value)
    resolution = integer('resolution', 1, 96)
    if resolution not in (1, 24, 96):
        raise ValueError('Choose daily, hourly or 15-minute resolution.')
    integer('spinup', 0, 3650)

    def interval(a, b):
        """Validate a temporal date interval."""
        start, end = (date.fromisoformat(c[a]), date.fromisoformat(c[b]))
        if start > end:
            raise ValueError('Start date must be on or before its end date.')
        return (start, end)
    if method == 'manual':
        a, b = interval('train_start', 'train_end')
        if c.get('evaluation_enabled', True):
            x, y = interval('eval_start', 'eval_end')
            if a <= y and x <= b:
                raise ValueError('Training and evaluation dates must not overlap.')
        return
    a, b = interval('common_start', 'common_end')
    day_count = (b - a).days + 1
    count = day_count * resolution
    if count < 2:
        raise ValueError('Select at least two scored time steps.')
    if method in ('random', 'random_kfold', 'deterministic_kfold') and day_count < 2:
        raise ValueError('Select at least two calendar days for this split method.')
    if method in ('random', 'random_block', 'random_kfold'):
        integer('split_seed', 0, 4294967295)
    if method in ('random', 'random_block', 'deterministic_block'):
        fraction = float(c['train_fraction'])
        if not math.isfinite(fraction) or not 0 < fraction < 1:
            raise ValueError('Training fraction must be greater than 0 and less than 1.')
    if method in ('deterministic_block', 'random_block') and c.get('block_unit', 'days') == 'days':
        block = integer('block_days', 1, count - 1)
        if fraction * count < block or (1 - fraction) * count < block:
            lower = math.ceil(100 * block / count) / 100
            upper = math.floor(100 * (1 - block / count)) / 100
            raise ValueError(f'This block size requires a training fraction between {lower:.2f} and {upper:.2f}. Reduce block size for a wider range.')
    if method in ('deterministic_block', 'random_block') and c.get('block_unit') == 'water_years':
        import numpy as np
        from sage_water_year_blocks import water_year_blocks
        days = np.arange(np.datetime64(a), np.datetime64(b) + np.timedelta64(1, 'D'))
        parts = water_year_blocks(days, c)
        if len(parts) < 2:
            raise ValueError('Choose at least two water-year blocks.')
        smallest = min((len(part) for part in parts))
        if min(fraction, 1 - fraction) * len(days) < smallest:
            raise ValueError('Training and evaluation must each allow at least one water-year block.')
    if method in ('deterministic_kfold', 'random_kfold'):
        folds = integer('folds', 2, day_count)
        integer('fold', 1, folds)
    if method == 'rainfall_block':
        first = a.year if (a.month, a.day) <= (10, 1) else a.year + 1
        last = b.year if (b.month, b.day) >= (9, 30) else b.year - 1
        if last - first < 2:
            raise ValueError('Rainfall-block splitting requires at least two complete October–September water years.')
