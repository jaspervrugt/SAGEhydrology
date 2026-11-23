# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_report_metadata.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_report_metadata."""

def snapshot(c):
    """Snapshot."""
    from regional_maps import attributes, zones
    rows = attributes(c['data_root'], c['region'], c['sage_root'])
    zone = zones(rows, c['region'], c.get('min_per_zone', 5), c.get('max_zones', 8)) if rows else {}
    selected = {g: rows[g] for g in c['gauges'] if g in rows}
    from sage_basin_names import display_names
    names = display_names(c, c['gauges'])
    return dict(basin_metadata=selected, basin_names=names, basin_zones={g: zone[g] for g in c['gauges'] if g in zone}, zone_names=list(dict.fromkeys(zone.values())))
