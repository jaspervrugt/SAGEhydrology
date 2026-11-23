# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_basin_names.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_basin_names."""
import csv, re, warnings
from functools import lru_cache
from pathlib import Path
from gauge_names import standardize, text
BUNDLED = Path(__file__).parent / 'assets' / 'sage_station_names'

def lookup_path(override, c):
    """Lookup path."""
    candidates = [Path(override), Path(c.get('sage_root', '')) / 'utils' / 'metadata' / override, BUNDLED / Path(override).name]
    return next((path for path in candidates if path.is_file()), None)

@lru_cache(maxsize=32)
def _lookup(file, mtime, size):
    """Lookup."""
    with Path(file).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        columns = {str(name).lower(): name for name in reader.fieldnames or []}
        identity = next((columns[key] for key in ('gauge_id', 'station_id', 'official_id', 'id') if key in columns), None)
        name = next((columns[key] for key in ('gauge_name', 'station_name', 'name') if key in columns), None)
        if identity is None or name is None:
            raise ValueError('Station-name table requires gauge_id and gauge_name columns.')
        result = {}
        for row in reader:
            key = re.sub('\\.0+$', '', text(row.get(identity)).upper())
            value = text(row.get(name))
            if key and value:
                result.setdefault(key, value)
    return result

def apply_overrides(names, schema, c):
    """Apply overrides."""
    override = schema.get('metadata', {}).get('name_override_file')
    if not override:
        return names
    path = lookup_path(override, c)
    if path is None:
        warnings.warn('Station-name lookup not found: ' + str(override), RuntimeWarning)
        return names
    try:
        stat = path.stat()
        lookup = _lookup(str(path.resolve()), stat.st_mtime_ns, stat.st_size)
        return {g: lookup.get(lookup_identity(g, schema), name) for g, name in names.items()}
    except (OSError, ValueError, csv.Error) as exc:
        warnings.warn('Could not apply station-name lookup ' + str(path) + ': ' + str(exc), RuntimeWarning)
        return names

def lookup_identity(gauge, schema):
    """
    Use MATLAB's output ID for lookup, retaining Python's source ID
    keys.
    """
    value = str(gauge).strip()
    identity = schema.get('id', {})
    if identity.get('lowercase'):
        value = value.lower()
    if identity.get('uppercase'):
        value = value.upper()
    from sage_attribute_reader import pairs
    for pattern, replacement in pairs(identity.get('output_regex', [])):
        value = re.sub(pattern, replacement, value)
    return re.sub('\\.0+$', '', value.upper())

def metadata_names(rows, schema, c):
    """Metadata names."""
    meta = schema.get('metadata', {})
    region = schema.get('region', c['region'])

    def seq(value):
        """Seq."""
        return value if isinstance(value, list) else [value] if value else []

    def column(row, candidates):
        """Column."""
        lookup = {str(k).lower(): k for k in row}
        return next((lookup[str(n).lower()] for n in seq(candidates) if str(n).lower() in lookup), None)

    def title(value, apostrophe):
        """Title."""
        result = []
        capitalize = True
        for ch in value.lower().strip():
            if ch.isalpha():
                result.append(ch.upper() if capitalize else ch)
                capitalize = False
            else:
                result.append(ch)
                if ch in ' -/(' or (ch == "'" and apostrophe):
                    capitalize = True
        return ''.join(result)
    result = {}
    separator = meta.get('name_separator', ': ')
    for gauge, row in rows.items():
        fallback = meta.get('name_fallback_prefix', '') + gauge
        components = seq(meta.get('name_components'))
        if components:
            primary = column(row, [components[0]])
            value = text(row.get(primary))
            value = separator.join([value] + [text(row.get(column(row, [key]))) for key in components[1:] if text(row.get(column(row, [key])))]) if value else fallback
            secondary = [text(row.get(column(row, [key]))) for key in seq(meta.get('parenthetical_components'))]
            secondary = [s for s in secondary if s]
            if secondary:
                value += ' (' + separator.join(secondary) + ')'
        else:
            key = column(row, meta.get('name_sources', []))
            value = text(row.get(key)) if key is not None else fallback
        if meta.get('standardize', True):
            value = standardize(value, region, gauge)
        transform = meta.get('name_transform', '')
        if transform in ('title_case', 'title_case_words'):
            value = title(value, transform == 'title_case')
        result[gauge] = value
    return apply_overrides(result, schema, c)

def display_names(c, gauges=None):
    """Display names."""
    from xia_camels_us import matches, project_region, standard_root
    root = Path(c['data_root'])
    if c['region'] == 'CAMELS_US' and matches(root):
        source = standard_root(root, c['region']) / 'camels_name.txt'
        names = {}
        if source.is_file():
            with source.open(encoding='utf-8-sig', newline='') as handle:
                for row in csv.DictReader(handle, delimiter=';'):
                    gauge = row['gauge_id'].zfill(8)
                    names[gauge] = standardize(row['gauge_name'], c['region'], gauge)
        if gauges is None:
            return names
        return {g: names.get(g, '') for g in gauges}
    xia_gb = c['region'] == 'CAMELS_GB' and project_region(root) == 'CAMELS_GB'
    if xia_gb:
        c = dict(c, data_root=str(standard_root(root, c['region'])))
    from sage_attribute_reader import read, metadata
    try:
        names = {g: r['display_name'] for g, r in read(c).items()}
    except (FileNotFoundError, ValueError, KeyError) as exc:
        if xia_gb:
            return {g: '' for g in gauges or []}
        schema = metadata().get(c['region'], {}).get('attributes', {})
        if gauges is None or not schema.get('metadata', {}).get('name_override_file'):
            raise
        warnings.warn('Full basin attributes unavailable; using stored station names: ' + str(exc), RuntimeWarning)
        names = apply_overrides({g: standardize('', c['region'], g) for g in gauges}, schema, c)
    if gauges is None:
        return names
    selected = {g: names.get(g, standardize('', c['region'], g)) for g in gauges}
    return apply_overrides(selected, metadata().get(c['region'], {}).get('attributes', {}), c)
