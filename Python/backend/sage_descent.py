# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_descent.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_descent."""
import numpy as np

def update(arrays, gradient, first, second, step, c, rate):
    """step is one-based, as in SAGE's optimizer update counter."""
    update2 = old2 = weights2 = 0.0
    for i, (array, g) in enumerate(zip(arrays, gradient)):
        old = array.copy()
        old2 += float(np.sum(old * old))
        if i % 2 == 0:
            weights2 += float(np.sum(old * old))
        if c['optimizer'] == 2:
            first[i] *= c['beta1']
            first[i] += (1 - c['beta1']) * g
            second[i] *= c['beta2']
            second[i] += (1 - c['beta2']) * g * g
            direction = first[i] / (1 - c['beta1'] ** step) / (np.sqrt(second[i] / (1 - c['beta2'] ** step)) + c['epsilon'])
            if i % 2 == 0:
                array *= 1 - rate * c['weight_decay']
        else:
            magnitude = np.linalg.norm(g)
            direction = g / magnitude if magnitude > 0 else g
        array -= rate * direction
        update2 += float(np.sum((array - old) ** 2))
    absolute = float(np.sqrt(update2))
    return dict(updateNorm=absolute, relativeUpdateNorm=absolute / max(float(np.sqrt(old2)), np.finfo(float).eps), moment1Norm=float(np.sqrt(sum((np.sum(v * v) for v in first)))) if c['optimizer'] == 2 else None, moment2Norm=float(np.sqrt(sum((np.sum(v * v) for v in second)))) if c['optimizer'] == 2 else None, weightDecayNorm=float(rate * c['weight_decay'] * np.sqrt(weights2)) if c['optimizer'] == 2 else 0.0)
