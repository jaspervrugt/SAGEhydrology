# ==========================================================================
# SAGEhydrology | Python source
# File: python/input_period.py
#
# Scientific period settings shared by input caches and GUI timelines.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Scientific period settings shared by input caches and GUI timelines."""

from datetime import date
import re

import numpy as np


def unscored_mask(dates, text):
    """Return the samples in inclusive calendar-date exclusions."""
    dates = np.asarray(dates)
    mask = np.zeros(dates.shape, dtype=bool)
    for item in re.split(r"[;,]", text) if text.strip() else ():
        try:
            first, last = (
                date.fromisoformat(value.strip()) for value in re.split(r"\||\.\.", item)
            )
            if last < first:
                raise ValueError
        except ValueError:
            raise ValueError(
                "Use YYYY-MM-DD | YYYY-MM-DD for each unscored range, "
                "with the earlier date first. Separate ranges with semicolons."
            ) from None
        mask |= (dates >= np.datetime64(first)) & (
            dates < np.datetime64(last) + np.timedelta64(1, "D")
        )
    return mask


PERIOD_KEYS = {
    "unscored_ranges",
    "split",
    "spinup",
    "resolution",
    "train_start",
    "train_end",
    "eval_start",
    "eval_end",
    "common_start",
    "common_end",
    "evaluation_enabled",
    "block_days",
    "block_unit",
    "block_years",
    "water_year_start_month",
    "train_fraction",
    "split_seed",
    "folds",
    "fold",
    "temporal_rng",
}
