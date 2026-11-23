# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_basin_presets.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_basin_presets."""
from pathlib import Path
import re
from site_schema import region_catalogue
BUNDLED = Path(__file__).parent / 'assets' / 'basin_presets'

def bundled_file(region, name):
    """Bundled file."""
    path = BUNDLED / region / Path(name).name
    return path.resolve() if path.is_file() else None

def python_split_path(c):
    """Python split path."""
    value = c.get('basin_split_file', '')
    if not value:
        return value
    path = Path(value)
    bundled = bundled_file(c['region'], path.name)
    source = Path(c.get('sage_root', ''))
    standard = not path.is_absolute() and (not path.is_file()) or ('sagehydrology' in [part.lower() for part in path.parts] and 'regions' in [part.lower() for part in path.parts]) or path.is_relative_to(source)
    return str(bundled) if bundled and standard else value

def read_ids(path, region):
    """Read ids."""
    ids = [x for x in re.split('[\\s,;]+', Path(path).read_text(encoding='utf-8-sig').strip()) if x]
    if region == 'CAMELS_US':
        ids = [re.sub('\\.0+$', '', x).zfill(8) for x in ids]
    return ids

def match_ids(ids, available, region=None):
    """
    Match numeric station IDs independent of padding; retain
    catalogue spelling.

    Exact matches win. Ambiguous numeric aliases are deliberately
    not guessed.
    File order is retained because it defines prescribed train/eval
    membership.
    """
    available = {str(g) for g in available}
    numeric = {}
    for gauge in available:
        if gauge.isascii() and gauge.isdigit():
            numeric.setdefault(str(int(gauge)), []).append(gauge)
    aliases_by_region = {}
    if region:
        from regional_maps import identifier
        for gauge in available:
            aliases_by_region.setdefault(identifier(gauge, region), []).append(gauge)
    result = []
    for value in ids:
        gauge = str(value)
        regional = aliases_by_region.get(identifier(gauge, region), []) if region else []
        if gauge not in available and len(regional) == 1:
            gauge = regional[0]
        aliases = numeric.get(str(int(gauge)), []) if gauge.isascii() and gauge.isdigit() else []
        result.append(gauge if gauge in available or len(aliases) != 1 else aliases[0])
    return result

def resolve_file(c):
    """Resolve file."""
    path = Path(python_split_path(c))
    if path.is_file():
        return path
    if path.parent.name == 'P06_CaRTPlus_US_daily':
        renamed = {'xia_gb_split.txt': 'GB_669_Xia_602_67_10pct_split.txt', 'xia_aus_split.txt': 'AUS_561_Xia_505_56_10pct_split.txt'}.get(path.name)
        if renamed and path.with_name(renamed).is_file():
            return path.with_name(renamed)
    meta = region_catalogue(c['sage_root'])[c['region']]
    for base in (Path(meta['source']).parent, Path(c['data_root'])):
        candidate = base / path
        if candidate.is_file():
            return candidate
    raise ValueError('Basin split file not found: ' + str(path))
_legacy_presets = presets

def presets(c):
    """Presets."""
    from xia_camels_us import matches
    if c['region'] == 'CAMELS_US' and matches(c.get('data_root', '')):
        return [dict(total=671, training=603, evaluation=68, universe=str(bundled_file('CAMELS_US', 'US_671_basins.txt')), split=str(bundled_file('CAMELS_US', 'US_671_Xia_603_68_split.txt')))]
    choices = _legacy_presets(c)
    from runtime_metadata import load
    data = load(str(c.get('sage_root', '')))
    preferred = data['regions'][c['region']].get('basin_counts', {}).get(str(c['resolution']), {}).get('universe') if data else None
    return sorted(choices, key=lambda p: (p['total'] != preferred, p['total'])) if data else choices
