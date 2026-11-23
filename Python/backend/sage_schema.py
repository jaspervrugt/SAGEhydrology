# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_schema.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_schema."""
from copy import deepcopy
import csv
import re
from pathlib import Path
import numpy as np
from engine import validate
from site_schema import session_defaults as site_defaults
from regional_maps import attributes as regional_attributes
DEFAULT_ATTRIBUTES = ['p_mean', 'pet_mean', 'frac_snow', 'aridity', 'high_prec_freq', 'high_prec_dur', 'low_prec_freq', 'low_prec_dur', 'carbonate_rocks_frac', 'geol_permeability', 'soil_depth_pelletier', 'soil_depth_statsgo', 'soil_porosity', 'soil_conductivity', 'max_water_content', 'sand_frac', 'silt_frac', 'clay_frac', 'elev_mean', 'slope_mean', 'area_gages2', 'frac_forest', 'lai_max', 'lai_diff', 'gvf_max', 'gvf_diff']

def session_defaults():
    """Return the default SAGE run configuration."""
    s = site_defaults()
    from sage_paths import OUTPUT_ROOT
    s.update(application='SAGE', block_unit='days', block_years=1, water_year_start_month=10, output_root=str(OUTPUT_ROOT), optimizer=2, max_iter=250, starts=1, jkge_window=7, attributes=DEFAULT_ATTRIBUTES[:], hidden_widths=[32], activations=['tanh'], network_seed=0, basin_train_fraction=0.8, basin_split_seed=1, learning_rate=0.01, beta1=0.9, beta2=0.999, epsilon=1e-08, weight_decay=0.0001, clip_norm=1.0, lr_schedule='constant', lr_min=0.0001, lr_hold=0.2)
    from basin_pool import cpu_inventory
    cpus = cpu_inventory()
    workers = max(1, cpus['physical'] or cpus['logical'] or 1)
    s.update(parallel=True, workers=workers, max_pending=2 * workers)
    s['resume_checkpoint'] = ''
    s['postprocess_selection'] = None
    s['evaluation_enabled'] = True
    s['temporal_rng'] = 'matlab'
    s['unscored_ranges'] = ''
    s.update(attribute_preprocessing='matlab', attribute_scaling='zscore', attribute_fit_cohort='training', attribute_missing='reject', attribute_methods={}, attribute_methods_enabled=False, basin_sampling='random', basin_split_file='')
    s.update(observations=['Q'], observation_root='', observation_weights=dict(Q=1.0, SWE=0.5, SM=0.5), observation_transforms=dict(Q='none', SWE='none', SM='none'), observation_normalization='auto', observation_scales=dict(Q=1.0, SWE=1.0, SM=1.0))
    return s

def strip_legacy_user_model_definition(config):
    """Remove editor metadata saved by older SAGE Python sessions."""
    models = config.get('model_parameters', {})
    if 'definition' in models.get('99', {}):
        config['model_parameters'] = dict(models)
        config['model_parameters']['99'] = {key: value for key, value in models['99'].items() if key != 'definition'}

def load_session(path):
    """Load and validate a saved SAGE configuration."""
    if Path(path).suffix.lower() == '.mat':
        from sage_mat_session import configuration
        return configuration(path, session_defaults())
    document = __import__('sage_paths').load_result(path)
    if not isinstance(document, dict) or document.get('application') != 'SAGE':
        raise ValueError('Select a SAGE session JSON, MATLAB session MAT, or SAGE run.json.')
    if isinstance(document.get('config'), dict):
        document = dict(document['config'], application='SAGE')
    s = session_defaults()
    s.update(document)
    strip_legacy_user_model_definition(s)
    s['parameter_scales'] = document.get('parameter_scales', ['linear'] * len(s['lower']))
    if 'attribute_methods_enabled' not in document:
        s['attribute_methods_enabled'] = bool(document.get('attribute_methods', {}))
    if not isinstance(s.get('config'), (list, tuple)) or len(s['config']) != 3:
        raise ValueError('Invalid SAGE numerical configuration: expected three smoothing values.')
    template = session_defaults()
    for key, default in template.items():
        if key in s and type(default) in (int, float) and isinstance(s[key], str):
            try:
                number = float(s[key])
                if not np.isfinite(number) or (type(default) is int and number != int(number)):
                    raise ValueError()
                s[key] = type(default)(number)
            except (ValueError, TypeError):
                raise ValueError(f'Invalid numeric session setting: {key}') from None
    for key in ('lower', 'upper', 'y0', 'settings', 'config'):
        try:
            s[key] = [float(value) for value in s[key]]
        except (ValueError, TypeError):
            raise ValueError(f'Invalid numeric session array: {key}') from None
    from sage_paths import SOFTWARE_ROOT, SAGE_RESOURCES, PYTHON_ROOT
    s.update(software_root=str(SOFTWARE_ROOT), site_root=str(PYTHON_ROOT), sage_root=str(SAGE_RESOURCES))
    if 'SITEhydrology' in str(s.get('output_root', '')):
        s['output_root'] = str(__import__('sage_paths').OUTPUT_ROOT)
    if 'temporal_rng' not in document:
        s['temporal_rng'] = 'numpy'
    if 'attribute_preprocessing' not in document:
        s['attribute_preprocessing'] = 'training'
    from sage_parameter_scaling import validate as validate_scales, modes
    validate_scales(s)
    s['parameter_scales'] = modes(s)
    return s

def attribute_rows(s):
    """Read the selected basin attribute rows."""
    from sage_attribute_reader import read, metadata
    import xia_camels_us
    if xia_camels_us.project_region(s['data_root']) == s['region']:
        return xia_camels_us.attributes(s['data_root'])
    if s['region'] in metadata():
        rows = read(s)
        if s['region'] == 'CAMELS_AU' and 'npp_diff' in s.get('attributes', []):
            for row in rows.values():
                monthly = np.array([row.get(f'npp_{month}', np.nan) for month in range(1, 13)], dtype=float)
                row['npp_diff'] = float(np.ptp(monthly))
        return rows
    if s['region'] != 'CAMELS_US':
        rows = regional_attributes(s['data_root'], s['region'], s['sage_root'])
        from regional_maps import identifier
        from regional_data import _table
        import pandas as pd
        root = Path(s['data_root'])
        files = []
        for folder in (root, root / 'attributes', root / 'daily' / 'attributes', root / 'hourly' / 'attributes'):
            if folder.exists():
                files.extend((p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in ('.csv', '.txt', '.xlsx') and any((t in p.name.lower() for t in ('attribute', 'topograph', 'climate', 'landcover', 'soil', 'geology', 'simulation_benchmark')))))
        keys = ('gauge_id', 'gauge', 'id', 'station_id', 'basin_id', 'catchment_id', 'station_number', 'site_code', 'code_station', 'sit_code_h3', 'gauge_number', 'river_id', 'rid')
        for file in sorted(set(files)):
            frame = pd.read_excel(file) if file.suffix.lower() == '.xlsx' else _table(file)
            lookup = {str(col).lower().strip(): col for col in frame.columns}
            key = next((lookup[k] for k in keys if k in lookup), None)
            if key is None:
                continue
            for row in frame.to_dict('records'):
                gauge = identifier(row[key], s['region'])
                target = rows.setdefault(gauge, {})
                canonical = {'gauge_lat', 'gauge_lon', 'elev_mean', 'aridity', 'frac_snow', 'p_mean', 'pet_mean'}
                target.update({str(k): v for k, v in row.items() if k != key and (not (k in canonical and k in target))})
        return rows
    rows = {}
    for name in ('topo', 'clim', 'soil', 'geol', 'vege', 'hydro'):
        path = Path(s['data_root']) / f'camels_{name}.txt'
        if path.is_file():
            with path.open(encoding='utf-8-sig') as f:
                for r in csv.DictReader(f, delimiter=';'):
                    rows.setdefault(r['gauge_id'], {}).update(r)
    return rows

def attribute_catalog(s, *, include_extra=False):
    """Static catalog labels/default/selectable IDs; no MATLAB execution."""
    from xia_camels_us import attribute_names, project_region, ready
    root = s.get('data_root', '')
    metadata = Path(__file__).with_name('sage_region_metadata.json')
    if metadata.is_file():
        from sage_attribute_reader import metadata as cached_metadata
        region = cached_metadata().get(s['region'], {})
        catalog = region.get('catalog')
        if catalog:
            names = catalog['names']
            labels = dict(zip(names, catalog['labels']))
            if s['region'] == 'CAMELS_AU':
                labels['npp_diff'] = 'Seasonal amplitude of net primary productivity'
            if include_extra:
                labels.update(region.get('extra_attribute_labels', {}))
            if s['region'] == 'CAMELS_GB' and project_region(root) == 'CAMELS_GB' and ready(root):
                names = attribute_names(root)
                return ({name: labels.get(name, name) for name in names}, names, [])
            selectable = np.atleast_1d(catalog.get('selectable', [True] * len(names)))
            return (labels, [names[int(i) - 1] for i in np.atleast_1d(catalog['default_ids'])], {n for n, enabled in zip(names, selectable) if not enabled})
    return ({}, [], set())
from functools import lru_cache

@lru_cache(maxsize=32)
def _prescribed_partition(gauges, region, path, stamp, inline, n):
    """Read the prescribed training and evaluation basin split."""
    from sage_basin_presets import match_ids
    import re
    text = ' '.join(inline) if inline else Path(path).read_text(encoding='utf-8-sig')
    ids = [g for g in re.split('[\\s,;]+', text.strip()) if g]
    if region == 'CAMELS_US':
        ids = [g.zfill(8) for g in ids]
    ids = match_ids(ids, gauges, region)
    if len(ids) != len(gauges) or len(set(ids)) != len(gauges) or set(ids) != set(gauges):
        raise ValueError('Split file must list each selected basin once, training basins first.')
    positions = {g: i for i, g in enumerate(gauges)}
    order = tuple((positions[g] for g in ids))
    return (order[:n], order[n:])

def basin_partition(c):
    """Resolve training and evaluation basin membership."""
    count = len(c['gauges'])
    n = max(1, min(count, int(np.floor(count * c['basin_train_fraction'] + 1e-10))))
    rng = np.random.RandomState(c['basin_split_seed'])
    mode = c.get('basin_sampling', 'random')
    if mode == 'file':
        from sage_basin_presets import resolve_file
        inline = tuple(c.get('basin_split_ids') or ())
        path = resolve_file(c) if not inline else None
        info = path.stat() if path is not None else None
        stamp = (info.st_mtime_ns, info.st_size, info.st_ctime_ns) if info else None
        train, evaluation = _prescribed_partition(tuple(c['gauges']), c['region'], str(path), stamp, inline, n)
        return (np.array(train, int), np.array(evaluation, int))
    elif mode == 'zone':
        from regional_maps import zones
        groups = zones(attribute_rows(c), c['region'], c['min_per_zone'], c['max_zones'])
        if any((g not in groups for g in c['gauges'])):
            raise ValueError('Zone sampling needs zone attributes for each basin.')
        names = list(dict.fromkeys((groups[g] for g in c['gauges'])))
        selected = []
        for zone in names:
            members = np.array([i for i, g in enumerate(c['gauges']) if groups[g] == zone])
            members = rng.permutation(members)
            selected.extend(members[:int(np.floor(n * len(members) / count + 0.5))])
        if len(selected) > n:
            selected = list(rng.permutation(selected)[:n])
        if len(selected) < n:
            other = np.setdiff1d(np.arange(count), selected)
            selected.extend(rng.permutation(other)[:n - len(selected)])
        train = np.array(selected, int)
    else:
        train = rng.permutation(count)[:n]
    return (np.sort(train), np.setdiff1d(np.arange(count), train))

def attribute_matrix(s):
    """Construct the selected basin attribute matrix."""
    if not s['attributes']:
        raise ValueError('Select at least one catchment attribute.')
    rows = attribute_rows(s)
    matrix = []
    for gauge in s['gauges']:
        if gauge not in rows:
            raise ValueError(f'No attributes installed for {gauge}.')
        values = []
        for key in s['attributes']:
            try:
                value = float(rows[gauge].get(key, np.nan))
                values.append(np.nan if value in (-999, -9999, -99999) else value)
            except (ValueError, TypeError):
                values.append(np.nan)
        matrix.append(values)
    return np.array(matrix)

def execution_issues(s):
    """Collect configuration issues that prevent execution."""
    from regional_data import reader_issues
    issues = reader_issues(s['sage_root'], s['region'], s['resolution'])
    from sage_parameter_scaling import validate as validate_scales
    try:
        validate_scales(s)
    except (ValueError, TypeError, KeyError) as error:
        issues.append(str(error))
    if s.get('temporal_rng', 'matlab') not in ('matlab', 'numpy'):
        issues.append('Unknown temporal permutation policy.')
    seed = s['split_seed']
    if not np.isfinite(seed) or int(seed) != seed or (not 0 <= seed <= 2 ** 32 - 1):
        issues.append('Temporal split seed must be an integer from 0 to 4294967295.')
    if s['backend'] != 'cpp':
        issues.append('SAGE training requires the native C++ backend.')
    numeric = ('learning_rate', 'lr_min', 'lr_hold', 'beta1', 'beta2', 'epsilon', 'weight_decay', 'clip_norm', 'basin_train_fraction', 'train_fraction')
    if any((not np.isfinite(s[key]) for key in numeric)):
        issues.append('Training controls must be finite numbers.')
    if s['loss'] not in range(1, 8):
        issues.append('Unknown observation loss.')
    if s['optimizer'] not in (1, 2):
        issues.append('Choose gradient descent or Adam/AdamW.')
    if s['split'] == 'manual':
        try:
            from sage_science import temporal_masks
            first = min(s['train_start'], s['eval_start']) if s.get('evaluation_enabled', True) else s['train_start']
            last = max(s['train_end'], s['eval_end']) if s.get('evaluation_enabled', True) else s['train_end']
            temporal_masks(np.array([np.datetime64(first), np.datetime64(last)]), s)
            if s.get('evaluation_enabled', True) and s['train_start'] <= s['eval_end'] and (s['eval_start'] <= s['train_end']):
                issues.append('Training and evaluation periods overlap; use separate dates or disable evaluation.')
        except ValueError as exc:
            issues.append(str(exc))
    if s['split'] != 'manual':
        from datetime import date
        try:
            if date.fromisoformat(s['common_end']) < date.fromisoformat(s['common_start']):
                issues.append('Common period end must follow its start.')
        except ValueError:
            issues.append('Use YYYY-MM-DD common-period dates.')
        if s['split'] in ('random', 'deterministic_block', 'random_block') and (not 0 < s['train_fraction'] < 1):
            issues.append('Temporal training fraction must be in (0,1).')
        if s['split'] in ('deterministic_block', 'random_block') and int(s['block_days']) < 1:
            issues.append('Temporal block size must be positive.')
        if 'kfold' in s['split'] and (s['folds'] < 2 or not 1 <= s['fold'] <= s['folds']):
            issues.append('Invalid temporal fold count/index.')
    if s['starts'] != 1:
        issues.append('SAGE currently supports one shared-network initialization per run.')
    if not s.get('attributes'):
        issues.append('Select catchment attributes.')
    if not 0 < s['basin_train_fraction'] <= 1:
        issues.append('Basin training fraction must be in (0,1].')
    if not 0 < s['beta1'] < 1 or not 0 < s['beta2'] < 1:
        issues.append('Moment decays must be in (0,1).')
    if min(s['learning_rate'], s['epsilon'], s['lr_min']) <= 0:
        issues.append('Learning rates and epsilon must be positive.')
    if s['lr_min'] > s['learning_rate']:
        issues.append('Minimum learning rate cannot exceed the initial learning rate.')
    if min(s['weight_decay'], s['clip_norm']) < 0:
        issues.append('Weight decay and clipping norm cannot be negative.')
    if not 0 <= s['lr_hold'] < 1:
        issues.append('Hold fraction must be in [0,1).')
    if s['split'] in ('random', 'random_kfold') and s['loss'] == 7 and (s['jkge_operator'] < 3):
        issues.append('Random temporal sampling requires JKGE operator 3 or 4, as in MATLAB.')
    if s['lr_schedule'] not in ('constant', 'cosine', 'hold_cosine'):
        issues.append('Unknown learning-rate schedule.')
    if s.get('attribute_preprocessing', 'training') not in ('matlab', 'training', 'custom'):
        issues.append('Unknown attribute preprocessing policy.')
    if s.get('basin_sampling', 'random') not in ('random', 'zone', 'file'):
        issues.append('Unknown basin sampling policy.')
    if s['jkge_operator'] not in (1, 2, 3, 4) or s['jkge_m'] not in (1, 2) or s['jkge_window'] < 1:
        issues.append('Invalid JKGE benchmark controls.')
    if s.get('resume_checkpoint') and (not Path(s['resume_checkpoint']).is_file()):
        issues.append('The resume checkpoint does not exist.')
    from sage_observations import available
    if s.get('observation_root') and (not Path(s['observation_root']).exists()):
        issues.append('The auxiliary observation file or folder does not exist.')
    from sage_zscore import transform
    try:
        transform(s, 'Q')
    except ValueError as error:
        issues.append(str(error))
    selected = s.get('observations', ['Q'])
    availability = available(s)
    if not selected or any((name not in ('Q', 'SWE', 'SM') for name in selected)):
        issues.append('Choose Q, SWE and/or SM observation targets.')
    for name in selected:
        if not availability.get(name, False):
            issues.append(f'{name} is unavailable for this dataset/model. Choose a supported profile or supply observation CSVs in mm.')
    if selected != ['Q'] or s.get('observation_transforms', {}).get('Q') == 'zscore':
        from build_observables import library_path
        if not library_path().is_file():
            issues.append('Build named native observations: python python/build_observables.py')
    if s.get('observation_normalization', 'auto') not in ('auto', 'none', 'benchmark', 'initial_loss', 'manual'):
        issues.append('Unknown joint-loss normalization.')
    if any((not np.isfinite(v) or v < 0 for v in s.get('observation_weights', {}).values())):
        issues.append('Observation weights must be finite and nonnegative.')
    if any((not np.isfinite(v) or v <= 0 for v in s.get('observation_scales', {}).values())):
        issues.append('Manual observation scales must be finite and positive.')
    return issues

def worker_configuration(s):
    """Construct the validated configuration passed to workers."""
    issues = execution_issues(s)
    if issues:
        raise ValueError('\n'.join(issues))
    from hydrologic_data_actions import configure
    c = configure(deepcopy(s))
    strip_legacy_user_model_definition(c)
    if not c.get('resume_checkpoint'):
        for key in list(c):
            if key.startswith('_observation_'):
                del c[key]
    c = validate(c)
    from sage_parameter_scaling import modes
    c['parameter_scales'] = modes(c)
    if c['model'] == 99:
        from user_model_native import compile_model, directory
        folder = directory()
        c.setdefault('model_config', {}).update(_user_model_folder=str(folder), _user_model_library=str(compile_model(folder)), _resolution=c['resolution'])
    from sage_network import Network
    Network(len(c['attributes']), len(c['lower']), c['hidden_widths'], c['activations'], c['network_seed'])
    from sage_network import standardize
    from hydrologic_data_actions import exclusion_reason
    from sage_quality import active_configuration
    checked = active_configuration(c, [dict(gauge=g, eligible=not exclusion_reason(c['data_actions'], g)) for g in c['gauges']])
    train, _ = basin_partition(checked)
    standardize(attribute_matrix(checked), train, __import__('sage_attribute_preprocessing').settings(c))
    if c.get('resume_checkpoint'):
        from sage_checkpoint import load_training
        load_training(c['resume_checkpoint'], c)
    return c
