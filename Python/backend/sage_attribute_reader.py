# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_attribute_reader.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_attribute_reader."""
from copy import deepcopy
from functools import lru_cache
import json
import csv
from io import BytesIO, StringIO
from pathlib import Path
import re
import numpy as np
import pandas as pd
from regional_maps import identifier

@lru_cache(maxsize=1)
def metadata():
    """Metadata."""
    return json.loads(Path(__file__).with_name('sage_region_metadata.json').read_text(encoding='utf-8'))

def sequence(value):
    """Sequence."""
    return value if isinstance(value, list) else [value] if value else []

def pairs(value):
    """Pairs."""
    if not value:
        return []
    if isinstance(value[0], list):
        return list(zip(*value)) if len(value) == 2 else value
    half = len(value) // 2
    return list(zip(value[:half], value[half:]))

def valid_name(value):
    """Valid name."""
    value = re.sub('\\s+(\\w)', lambda m: m[1].upper(), str(value).strip())
    value = re.sub('[^a-zA-Z0-9_]', '_', value)
    return value if value and value[0].isalpha() else 'x' + value

def read(c):
    """Read."""
    schema = deepcopy(metadata()[c['region']]['attributes'])
    if c['resolution'] == 24 and 'hourly' in schema.get('profiles', {}):
        schema.update(schema['profiles']['hourly']['schema'])
    tables = sequence(schema['tables'])
    root = Path(c['data_root'])
    for folder in sequence(schema.get('root_candidates', [])):
        candidate = root / str(folder).replace('\\', '/')
        if all(((candidate / t.get('file', '')).is_file() for t in tables if t.get('required', True) and t.get('file'))):
            root = candidate
            break
    identity = schema['id']
    idcol = identity['column']

    def normalize(value):
        """Normalize."""
        if pd.isna(value):
            return ''
        text = str(value).strip()
        for pattern, replacement in pairs(identity.get('regex', [])):
            replacement = re.sub('\\$(\\d+)', lambda m: '\\g<' + m[1] + '>', replacement)
            text = re.sub(pattern, replacement, text)
        if identity.get('uppercase'):
            text = text.upper()
        if identity.get('lowercase'):
            text = text.lower()
        if identity.get('numeric_canonical') and re.fullmatch('\\d+', text):
            text = str(int(text))
        if identity.get('pad_width'):
            text = text.zfill(int(identity['pad_width']))
        return identifier(text, c['region'])

    def table(spec, path):
        """Table."""
        if not path.is_file():
            if spec.get('required', True):
                raise FileNotFoundError('Required SAGE attribute table: ' + str(path))
            return None
        from local_data import require_local
        require_local(path)
        from local_data import attribute_bytes
        raw = attribute_bytes(path)
        if path.suffix.lower() in ('.xlsx', '.xls'):
            try:
                frame = pd.read_excel(path, sheet_name=spec.get('sheet', 0), dtype=str)
            except ImportError:
                from sage_xlsx import read_sheet
                frame = read_sheet(path, spec.get('sheet', 0))
        else:
            delimiter = spec.get('delimiter')
            delimiter = {'\\t': '\t', '\\s+': '\\s+'}.get(delimiter, delimiter)
            if not delimiter:
                header = raw.splitlines()[0].decode(spec.get('encoding') or 'utf-8-sig') if raw else ''
                delimiter = ';' if header.count(';') > header.count(',') else ',' if ',' in header else '\t' if '\t' in header else '\\s+'
            comment = spec.get('comment_style')
            comment = comment if isinstance(comment, str) and len(comment) == 1 else None
            if 'hydrometry' in path.name.lower() and delimiter == ',':
                with StringIO(raw.decode(spec.get('encoding') or 'utf-8-sig'), newline='') as handle:
                    records = csv.reader(handle)
                    columns = next(records)
                    values = []
                    for row in records:
                        if len(row) > len(columns):
                            if 'comment' not in columns[-1].lower():
                                raise ValueError('Unexpected extra attribute columns: ' + str(path))
                            row = row[:len(columns) - 1] + [','.join(row[len(columns) - 1:])]
                        values.append(row + [None] * max(0, len(columns) - len(row)))
                frame = pd.DataFrame(values, columns=columns).replace('', np.nan)
            else:
                frame = pd.read_csv(BytesIO(raw), sep=delimiter, engine='c' if len(delimiter) == 1 else 'python', dtype=str, encoding=spec.get('encoding') or 'utf-8-sig', comment=comment, quoting=csv.QUOTE_NONE if spec.get('layout') == 'transposed_attributes' else csv.QUOTE_MINIMAL)
        if spec.get('drop_empty_rows'):
            frame = frame.dropna(axis=0, how='all')
        layout = spec.get('layout', '')
        if layout in ('transposed_attributes', 'catalog_attribute_matrix'):
            if layout == 'transposed_attributes':
                frame.columns = [str(n).strip('" ') for n in frame]
                frame = frame.map(lambda v: v.strip('" ') if isinstance(v, str) else v)
            key = spec.get('attribute_column') or frame.columns[0]
            start = int(spec.get('metadata_columns', 1)) if layout == 'catalog_attribute_matrix' else 1
            matrix = frame.set_index(key).iloc[:, max(0, start - 1):]
            matrix.index = [valid_name(n) for n in matrix.index]
            for derived in sequence(spec.get('derived', {})):
                if derived and derived.get('operation') == 'row_sum':
                    matrix.loc[derived['target']] = matrix.loc[derived['rows']].apply(pd.to_numeric, errors='coerce').sum(axis=0, min_count=1)
            frame = matrix.T.rename_axis(idcol).reset_index()
            if layout == 'catalog_attribute_matrix':
                desired = set(metadata()[spec.get('catalog_region', c['region'])]['catalog']['names']) - set(spec.get('exclude_catalog', []))
                frame = frame[[idcol] + [n for n in frame if n in desired]]
        frame = frame.rename(columns=dict(pairs(spec.get('column_renames', []))))
        filt = spec.get('row_filter', {})
        if filt:
            frame = frame[frame[filt['column']].astype(str) == str(filt['value'])]
        if spec.get('make_valid_names'):
            frame.columns = [valid_name(n) for n in frame]
        if spec.get('make_unique_names'):
            used = set()
            columns = []
            for name in frame:
                original = name
                i = 1
                while name in used:
                    name = f'{original}_{i}'
                    i += 1
                used.add(name)
                columns.append(name)
            frame.columns = columns
        key = next((k for k in sequence(spec['keys']) if k in frame), None)
        if key is None:
            raise ValueError(f"{path.name}: missing basin identifier column {spec['keys']}")
        target = spec.get('key_target') or idcol
        frame = frame.rename(columns={key: target})
        frame[target] = frame[target].map(normalize)
        if spec.get('drop_missing_key'):
            frame = frame[frame[target] != '']
        keep = sequence(spec.get('keep_columns', []))
        if keep:
            frame = frame[[target] + [n for n in keep if n in frame and n != target]]
        widen = sequence(spec.get('widen', {}).get('keys', []))
        if widen and frame[target].duplicated().any():
            rows = {}
            for row in frame.to_dict('records'):
                values = rows.setdefault(row[target], {target: row[target]})
                suffix = re.sub('[^a-z0-9]+', '_', '_'.join((str(row[k]) for k in widen)).lower()).strip('_') or 'x'
                for k, v in row.items():
                    if k not in [target] + widen:
                        values[valid_name(k + '_' + suffix)] = v
            frame = pd.DataFrame(rows.values())
        if spec.get('unique_key'):
            frame = frame.drop_duplicates(target, keep='first')
        exclude = spec.get('exclude_columns_regex')
        if exclude:
            frame = frame[[n for n in frame if n == target or not re.search(exclude, n)]]
        prefix = spec.get('prefix_columns', '')
        if prefix:
            frame = frame.rename(columns={n: prefix + n for n in frame if n != target and (not n.startswith(prefix))})
        exceptions = set(sequence(spec.get('numeric_text_exceptions', [])))
        for name in frame:
            if name == target or name in exceptions:
                continue
            numeric = pd.to_numeric(frame[name], errors='coerce')
            present = frame[name].notna() & (frame[name].astype(str).str.strip() != '')
            if present.any() and numeric[present].notna().all():
                frame[name] = numeric
            if spec.get('boolean_text') and present.any() and frame[name].astype(str).str.lower()[present].isin(['true', 'false', '0', '1']).all():
                frame[name] = frame[name].astype(str).str.lower().isin(['true', '1']).astype(float)
            if name in sequence(spec.get('absolute_columns', [])):
                frame[name] = numeric.abs()
            source_name = name.removeprefix(prefix) if prefix else name
            if any((source_name.startswith(p) for p in sequence(spec.get('fill_missing_zero_prefixes', [])))):
                frame[name] = numeric.fillna(0)
        if spec.get('drop_empty_columns'):
            frame = frame.dropna(axis=1, how='all')
        if spec.get('drop_empty_rows'):
            frame = frame.dropna(axis=0, how='all')
        fallback = spec.get('fallback_file')
        if fallback and (root / fallback).is_file():
            other = table({k: v for k, v in spec.items() if k != 'fallback_file'}, root / fallback)
            frame = frame.set_index(target).combine_first(other.set_index(target)).reset_index()
        return frame
    result = None
    for original in tables:
        files = sorted(root.glob(original['pattern'])) if original.get('pattern') else [root / original['file']]
        for path in files:
            spec = deepcopy(original)
            if spec.get('prefix_from_file'):
                prefix = re.sub('^catchments_hydrological_signatures_?', 'hs_', path.stem)
                spec['prefix_columns'] = valid_name(re.sub('[^a-zA-Z0-9]', '_', prefix).rstrip('_')) + '_'
            frame = table(spec, path)
            if frame is None:
                continue
            if result is None:
                result = frame
                continue
            key = spec.get('join_key') or idcol
            if key not in result or key not in frame:
                raise ValueError('Attribute join key is missing: ' + key)
            overlap = [n for n in frame if n in result and n != key]
            policy = spec.get('duplicate_policy', 'suffix')
            if policy == 'drop':
                frame = frame.drop(columns=overlap)
            elif policy == 'coalesce':
                merged = result.merge(frame, on=key, how=spec.get('join_type', 'left'), suffixes=('', '__fallback'))
                for n in overlap:
                    merged[n] = merged[n].combine_first(merged.pop(n + '__fallback'))
                result = merged
                continue
            elif policy == 'replace_left':
                result = result.drop(columns=overlap)
            else:
                occupied = set(result) | set(frame)
                renames = {}
                for n in overlap:
                    candidate = n + spec.get('duplicate_suffix', '_right')
                    i = 1
                    while candidate in occupied:
                        candidate = n + spec.get('duplicate_suffix', '_right') + '_' + str(i)
                        i += 1
                    renames[n] = candidate
                    occupied.add(candidate)
                frame = frame.rename(columns=renames)
            how = {'full': 'outer', 'outer': 'outer', 'inner': 'inner', 'right': 'right'}.get(spec.get('join_type'), 'left')
            result = result.merge(frame, on=key, how=how)
    if result is None:
        raise ValueError('No installed SAGE attribute tables.')
    additions = {}
    for alias in sequence(schema.get('aliases', [])):
        if alias['target'] in result:
            continue
        source = next((n for n in sequence(alias.get('sources', [])) if n in result or n in additions), None)
        if source is not None:
            additions[alias['target']] = result[source] if source in result else additions[source]
        elif alias.get('required') and alias.get('default') is None:
            raise ValueError('Required attribute missing: ' + alias['target'])
        elif alias.get('default') is not None:
            additions[alias['target']] = alias['default']
    if additions:
        result = pd.concat([result, pd.DataFrame(additions, index=result.index)], axis=1)
    if c['region'] == 'CAMELS_ET' and 'Shrub_bush' in result and ('Shrubbush' not in result):
        result = result.rename(columns={'Shrub_bush': 'Shrubbush'})
    rows = {}
    for row in result.to_dict('records'):
        gauge = normalize(row[idcol])
        if gauge:
            rows.setdefault(gauge, row)
    from regional_maps import attributes
    geographic = attributes(c['data_root'], c['region'], c['sage_root'])
    for gauge, row in rows.items():
        for key, value in geographic.get(gauge, {}).items():
            if key in ('gauge_lat', 'gauge_lon') or key not in row:
                row[key] = value
    from sage_basin_names import metadata_names
    names = metadata_names(rows, schema, c)
    for gauge, row in rows.items():
        row['display_name'] = names[gauge]
    return rows
