# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_quality_summary.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_quality_summary."""
from collections import Counter

def summary(rows):
    """Summary."""
    excluded = [row for row in rows if not row.get('eligible', True)]
    lines = [f'Data-quality screening: {len(excluded)} of {len(rows)} basins excluded.']
    reasons = Counter((row.get('reason') or 'Unspecified exclusion' for row in excluded))
    for reason, count in reasons.most_common(2):
        lines.append(f'Reason ({count}): {reason}.')
    if len(reasons) > 2:
        lines.append(f'Additional exclusion reasons: {len(reasons) - 2}.')
    if excluded:
        lines.append('Full gauge details: Region / Run notes and Data quality table.')
    return lines
