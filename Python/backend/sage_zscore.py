# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_zscore.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_zscore."""
import numpy as np
from sage_objective_error import UndefinedObjective

def transform(c, name):
    """Transform."""
    settings = c.get('observation_transforms', {})
    if not isinstance(settings, dict):
        raise ValueError('Observation transforms must be a named mapping.')
    if set(settings) - {'Q', 'SWE', 'SM'}:
        raise ValueError('Unknown observation transform response.')
    for value in settings.values():
        if value not in ('none', 'zscore'):
            raise ValueError('Use none or zscore for each response.')
    return settings.get(name, 'none')

def eligible(y):
    """Eligible."""
    y = np.asarray(y, float)
    return len(y) >= 2 and np.any(y != y[0])

def loss(y, sim):
    """Return objective and dL/dsim, including mean/scale derivatives."""
    y = np.asarray(y, float)
    sim = np.asarray(sim, float)
    if y.shape != sim.shape or not np.isfinite(y).all():
        raise ValueError('Z-score loss requires finite paired vectors.')
    if not np.isfinite(sim).all() or not eligible(y):
        raise UndefinedObjective('Z-score loss needs two variable paired samples.')
    centered = sim - sim[0] - np.mean(sim - sim[0])
    scale = np.sqrt(np.mean(centered * centered))
    if not np.isfinite(scale) or scale <= 0:
        raise UndefinedObjective('Z-score loss: simulation has zero variability.')
    zs = centered / scale
    yc = y - y[0] - np.mean(y - y[0])
    zy = yc / np.sqrt(np.mean(yc * yc))
    residual = zs - zy
    delta = 2 / (len(y) * scale) * (residual - residual.mean() - zs * np.mean(residual * zs))
    return (float(np.mean(residual * residual)), delta)
