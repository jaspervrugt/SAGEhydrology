# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_matlab_project.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_matlab_project."""
from copy import deepcopy
from datetime import date
from pathlib import Path
import re
from sage_schema import attribute_catalog
from site_schema import model_catalogue
from regional_data import region_registry
from regional_maps import identifier

def translate(cfg, base, path):
    """
    Shared declarative configuration translator for JSON and MAT
    sessions.
    """
    import numpy as np
    from sage_schema import session_defaults
    path = Path(path)
    c = session_defaults()
    for key in ('software_root', 'site_root', 'sage_root', 'output_root', 'workers', 'max_pending', 'verbose_progress'):
        if key in base:
            c[key] = deepcopy(base[key])
    c['region'] = cfg.get('region', c['region'])
    registry = region_registry(c['sage_root'])
    if c['region'] not in registry:
        raise ValueError('Unknown dataset ' + c['region'])
    c['data_root'] = str(Path(c['software_root']) / 'Data' / registry[c['region']]['data_root'])
    if cfg.get('dirD') and Path(cfg['dirD']).is_dir():
        c['data_root'] = cfg['dirD']
    if cfg.get('dataset') in ('xia_camels_us', 'xia_camels_gb'):
        from xia_camels_us import directory, ready
        variant = 'gb' if cfg['dataset'] == 'xia_camels_gb' else 'us'
        root = directory(c['software_root'], variant)
        if not ready(root):
            raise ValueError('Use Data Install on the Project 6 page first.')
        c['data_root'] = str(root)
    prd = cfg.get('prd', {})
    mdl = cfg.get('mdl', {})
    misc = cfg.get('misc', {})
    bas = cfg.get('bas', {})
    c['resolution'] = int(prd.get('dt', 1))
    c['model'] = int(mdl.get('model', 1))
    if 'backend' in mdl:
        c['backend'] = mdl['backend']
    catalog = model_catalogue(c['sage_root'], c['resolution'])
    model = catalog.get(str(c['model']))
    if model is None:
        raise ValueError('Unsupported hydrological model.')
    c.update(model_parameters=catalog, lower=[p['lower'] for p in model['parameters']], upper=[p['upper'] for p in model['parameters']], parameter_scales=['linear'] * len(model['parameters']), y0=[1e-05] * model['states'])
    for old, new in [('dts', 'train_start'), ('dte', 'train_end'), ('des', 'eval_start'), ('dee', 'eval_end')]:
        if old in prd:
            day, month, year = map(int, prd[old])
            c[new] = date(year, month, day).isoformat()
    c['split'] = prd.get('method', 'manual')
    c['spinup'] = int(prd.get('spinup', 365))
    c['unscored_ranges'] = prd.get('unscored_ranges', '')
    c['evaluation_enabled'] = c['split'] != 'manual' or not (prd.get('eval_was_empty', False) or mdl.get('mode', 4) in (1, 3))
    for old, new in [('ds', 'common_start'), ('de', 'common_end')]:
        if old in prd:
            day, month, year = map(int, prd[old])
            c[new] = date(year, month, day).isoformat()
    for old, new in [('block_size', 'block_days'), ('train_frac', 'train_fraction'), ('seed', 'split_seed'), ('n_folds', 'folds'), ('fold', 'fold')]:
        if old in prd and np.size(prd[old]):
            c[new] = prd[old]
    for old, new in [('block', 'block_days'), ('tr_frac', 'train_fraction')]:
        if old in prd and np.size(prd[old]):
            c[new] = prd[old]
    c['basin_sampling'] = bas.get('sample', 'file' if cfg.get('file_split') else 'random')
    for old, new in [('seed', 'basin_split_seed'), ('min_per_zone', 'min_per_zone'), ('max_zones', 'max_zones')]:
        if old in bas and np.size(bas[old]):
            c[new] = bas[old]
    c['basin_train_fraction'] = float(bas.get('K_t', bas.get('K', 1))) / float(bas.get('K', 1))
    source = Path(registry[c['region']]['source']).parent

    def ids_file(name):
        """Ids file."""
        from sage_basin_presets import bundled_file
        packaged = bundled_file(c['region'], name)
        candidates = ([packaged] if packaged else []) + [source / name, Path(c['data_root']) / name, path.parent / name]
        found = next((p for p in candidates if p.is_file()), None)
        if found is None:
            raise ValueError('Cannot find basin list ' + name)
        return found
    if cfg.get('_import_gauges'):
        c['gauges'] = cfg['_import_gauges']
    elif cfg.get('file_univ'):
        c['gauges'] = [identifier(g, c['region']) for g in re.split('[\\s,;]+', ids_file(cfg['file_univ']).read_text().strip()) if g]
    if cfg.get('file_split') and c['basin_sampling'] == 'file' and (not cfg.get('_import_gauges')):
        c['basin_split_file'] = str(ids_file(cfg['file_split']))
    if cfg.get('dataset') == 'xia_camels_us':
        from xia_camels_us import catalogue
        gauges = [row[0] for row in catalogue(c['data_root'])]
        from sage_basin_presets import presets, read_ids
        split = presets(c)[0]['split']
        ids = read_ids(split, c['region'])
        if len(gauges) != 671 or len(ids) != 671 or set(ids) != set(gauges):
            raise ValueError('Unexpected Project 6 basin membership.')
        c['gauges'] = gauges
        c['basin_sampling'] = 'file'
        c['basin_split_file'] = split
    labels, defaults, excluded = attribute_catalog(c)
    names = list(labels)
    ids = np.atleast_1d(bas.get('id_attr', []))
    if any((int(i) != i or not 1 <= i <= len(names) for i in ids)):
        raise ValueError('Invalid MATLAB attribute catalog index.')
    if 'attribute_names' in bas:
        c['attributes'] = list(bas['attribute_names'])
        if len(set(c['attributes'])) != len(c['attributes']) or any((name not in labels for name in c['attributes'])):
            raise ValueError('Invalid project attribute names.')
    else:
        c['attributes'] = [names[int(i) - 1] for i in ids] if 'id_attr' in bas else defaults
    if any((n in excluded for n in c['attributes'])):
        raise ValueError('Project selects attributes disabled by the MATLAB catalog.')
    net = misc.get('net', {})
    alg = misc.get('alg', {})
    loss = misc.get('loss', {})
    meteo = misc.get('meteo', {})
    c['hidden_widths'] = np.atleast_1d(net.get('h', [32])).astype(int).tolist()
    c['activations'] = np.atleast_1d(net.get('tf', ['tanh'])).tolist()
    c['network_seed'] = net.get('seed', 0)
    for old, new in [('method', 'optimizer'), ('i_max', 'max_iter'), ('lr', 'learning_rate'), ('beta_1', 'beta1'), ('beta_2', 'beta2'), ('vareps', 'epsilon'), ('wdecay', 'weight_decay'), ('clipn', 'clip_norm'), ('lr_min', 'lr_min'), ('lr_scheme', 'lr_schedule'), ('lr_hold', 'lr_hold')]:
        if old in alg:
            c[new] = alg[old]
    for old, new in [('fnc', 'loss'), ('method', 'jkge_operator'), ('n_win', 'jkge_window'), ('M', 'jkge_m'), ('observed', 'observations'), ('weight', 'observation_weights'), ('transform', 'observation_transforms')]:
        if old in loss:
            c[new] = loss[old]
    normalization = loss.get('normalization', {})
    if normalization:
        c['observation_normalization'] = normalization.get('method', 'auto')
    if normalization.get('scale'):
        c['observation_scales'].update(normalization['scale'])
    c['observations'] = np.atleast_1d(c['observations']).tolist()
    c['fdc_formulation'] = loss.get('fdc', {}).get('formulation', 1)
    profile_name = {1: 'Daily', 24: 'Hourly', 96: '15-minute'}.get(c['resolution'], 'Daily')
    profiles = registry[c['region']]['profiles']
    profile = profiles.get(profile_name, next(iter(profiles.values())))
    controls = profile.get('meteo', {})
    for old, new in [('data', 'product'), ('precip', 'precipitation'), ('temp', 'temperature'), ('pet', 'pet')]:
        control = controls.get(new, {})
        items = control.get('items', [])
        index = int(meteo.get(old, 1)) - 1
        if not items:
            continue
        if len(items) == 1 and (not control.get('enabled', True)):
            index = 0
        if not 0 <= index < len(items):
            raise ValueError('Invalid imported meteorology choice: ' + old)
        label = items[index]
        c[new] = next((p for p in ('daymet', 'maurer', 'nldas') if p in label.lower()), label) if new == 'product' else int(re.match('^(\\d+)\\s', label)[1]) if re.match('^(\\d+)\\s', label) else index + 1
    c['parallel'] = mdl.get('calc', '') in ('par', 'parfeval')
    for key in ('network_seed', 'basin_split_seed', 'split_seed', 'block_days', 'folds', 'fold', 'optimizer', 'max_iter', 'jkge_operator', 'jkge_m', 'jkge_window', 'fdc_formulation', 'spinup', 'min_per_zone', 'max_zones'):
        c[key] = int(c[key])
    ode = cfg.get('ode', {})
    for i, key in enumerate(('InitStep', 'MaxStep', 'MinStep', 'RelTol', 'AbsTol', 'Order', 'maxiter')):
        if key in ode:
            c['settings'][i] = float(ode[key])
    options = (cfg.get('model_options') or {}).get('model_' + str(c['model']), {})
    if options.get('config'):
        c['model_config'] = deepcopy(options['config'])
    if mdl.get('config'):
        c.setdefault('model_config', {}).update(mdl['config'])
    for i, key in enumerate(('T_sm', 'eps_m', 'rho')):
        if key in c.get('model_config', {}):
            c['config'][i] = c['model_config'][key]
    y0 = ode.get('y0', options.get('y0'))
    if y0 is not None and np.size(y0):
        c['y0'] = np.atleast_1d(y0).astype(float).tolist()
    c['diagnostics']['attribution'] = bool(misc.get('attr', False))
    c.update(application='SAGE', attribute_preprocessing='matlab', resume_checkpoint='', postprocess_selection=None)
    return c
