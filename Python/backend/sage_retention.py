# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_retention.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_retention."""

def selection(config, training, evaluation):
    """Selection."""
    roles = {'train': list(training), 'eval': list(evaluation)}
    explicit = config.get('postprocess_selection')
    if explicit is not None:
        allowed = set(training) | set(evaluation)
        return {g: [p for p in periods if p in ('train', 'evaluation')] for g, periods in explicit.items() if g in allowed and periods}
    diagnostic = config.get('diagnostic_selection', {})
    requested = diagnostic.get('selection')
    counts = diagnostic.get('counts', {})
    flags = diagnostic.get('scenarios', {})
    result = {}
    for role, ids in roles.items():
        prefix = 't' if role == 'train' else 'e'
        periods = [p for p, tag in [('train', prefix + 't'), ('evaluation', prefix + 'e')] if flags.get(tag, True)]
        if not periods:
            continue
        selected = [g for g in ids if requested and requested.get(g)]
        wanted = min(len(ids), max(0, int(counts.get(role, len(selected) if requested is not None else 5))))
        chosen = (selected + [g for g in ids if g not in selected])[:wanted]
        for gauge in chosen:
            retained = [p for p in (requested.get(gauge, periods) if requested is not None else periods) if p in periods]
            if retained:
                result[gauge] = retained
    return result

def restrict_periods(series, periods):
    """
    Keep continuous model arrays; expose only requested scored
    scenarios.
    """
    import numpy as np
    result = dict(series)
    for period in ('train', 'evaluation'):
        if period not in periods:
            result[period] = np.zeros_like(series[period], dtype=bool)
    return result
