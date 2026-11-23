# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_solver_recovery.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_solver_recovery."""
import numpy as np
from sage_descent import update

class RecoveryExhausted(RuntimeError):
    """No acceptable proposal was found within the bounded retry budget."""

def capture(network, first, second, gradient, step, rate):
    """Capture."""
    return dict(base=[v.copy() for v in network.arrays], first=[v.copy() for v in first], second=[v.copy() for v in second], gradient=[v.copy() for v in gradient], step=int(step), rate=float(rate), factor=1.0)

def restore(proposal, network, first, second):
    """Restore."""
    for target, source in [(network.arrays, proposal['base']), (first, proposal['first']), (second, proposal['second'])]:
        for array, saved in zip(target, source):
            array[:] = saved

def retry(proposal, network, first, second, config):
    """Retry."""
    restore(proposal, network, first, second)
    proposal['factor'] *= 0.5
    return update(network.arrays, proposal['gradient'], first, second, proposal['step'], config, proposal['rate'] * proposal['factor'])
