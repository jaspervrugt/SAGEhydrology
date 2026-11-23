# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_observations.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_observations."""
from sage_parameter_scaling import physical
from pathlib import Path
import numpy as np
import pandas as pd
from sage_science import full_loss
from sage_objective_error import UndefinedObjective
from sage_zscore import transform, eligible

def available(c):
    """Determine which observed response types are available."""
    if c.get('observation_root'):
        return {'Q': True, 'SWE': True, 'SM': c['model'] != 1 or bool(c.get('model_config', {}).get('sm_ind'))}
    return dict(Q=True, SWE=c['region'] in ('CAMELS_CL', 'CAMELS_CH', 'HYD_RESPONSES', 'MACH_US'), SM=c['region'] == 'HYD_RESPONSES' and c['model'] != 1)

def load(c, gauge, data):
    """Load the selected data or metadata resource."""
    result = dict(Q=np.asarray(data['observed'], dtype=np.float32))
    if not c.get('observation_root') and '_source_SWE' in data:
        for name in ('SWE', 'SM'):
            result[name] = np.asarray(data['_source_' + name], dtype=np.float32)
        return result
    dates = pd.DatetimeIndex(data['dates'])
    from regional_data import _table, _dates, table_profile
    root = Path(c['data_root'])
    frame = None
    clock = None
    columns = {}
    custom = c.get('observation_root', '')
    if custom:
        path = Path(custom)
        if path.is_dir():
            choices = [p for p in (path / (gauge + '.csv'), path / (gauge + '_observations.csv')) if p.is_file()]
            if len(choices) > 1:
                raise ValueError(f'{gauge}: ambiguous observation sidecars.')
            path = choices[0] if choices else path / (gauge + '.csv')
        if path.is_file():
            frame = _table(path)
            if 'gauge_id' in frame:
                ids = frame['gauge_id'].astype(str).str.removesuffix('.0')
                target = gauge
                if c['region'] == 'CAMELS_US':
                    ids = ids.str.zfill(8)
                frame = frame.loc[ids == target].copy()
            column = next((k for k in frame.columns if k.lower() in ('date', 'time', 'datetime', 'timestamp')), None)
            if column is None:
                raise ValueError('Observation CSV requires a date/time column.')
            clock = pd.to_datetime(frame[column], errors='raise')
            for key in ('SWE', 'SM'):
                col = next((k for k in (key + '_mm', key) if k in frame), None)
                if col:
                    columns[key] = col
    elif c['region'] == 'CAMELS_CL':
        from regional_special import _wide
        path = root / 'daily' / 'swe' / '13_CAMELScl_swe.txt'
        if not path.is_file():
            path = root / 'swe' / '13_CAMELScl_swe.txt'
        if path.is_file():
            values = _wide(path, gauge, date='gauge_id').reindex(dates).to_numpy(float)
            values[values <= -99] = np.nan
            result['SWE'] = values
    elif c['region'] == 'CAMELS_CH':
        path = root / 'daily' / 'timeseries' / 'observation_based' / f'CAMELS_CH_obs_based_{gauge}.csv'
        if not path.is_file():
            path = root / 'timeseries' / 'observation_based' / f'CAMELS_CH_obs_based_{gauge}.csv'
        if path.is_file():
            frame = _table(path)
            clock = pd.to_datetime(frame['date'])
            columns = {'SWE': 'swe(mm)'}
    elif c['region'] in ('HYD_RESPONSES', 'MACH_US'):
        path = Path(data['forcing_file'])
        frame = _table(path)
        profile = table_profile(c['sage_root'], c['region'], c['pet'], c.get('temperature', 1), c.get('precipitation', 1))
        clock = _dates(frame, profile)
        columns = {'SWE': 'SWECLQMD_mean' if c['region'] == 'HYD_RESPONSES' else 'SWE'}
        if c['region'] == 'HYD_RESPONSES' and all((k in frame for k in ('swvl1_mean_mean', 'swvl2_mean_mean', 'swvl3_mean_mean'))):
            values = frame[['swvl1_mean_mean', 'swvl2_mean_mean', 'swvl3_mean_mean']].apply(pd.to_numeric, errors='coerce').to_numpy(float, copy=True)
            values[np.isin(values, [-999, -99, -99.9, -99.99])] = np.nan
            frame['_SAGE_SM_mm'] = values @ np.array([70.0, 210.0, 720.0])
            columns['SM'] = '_SAGE_SM_mm'
    if frame is not None:
        index = pd.DatetimeIndex(clock)
        if index.duplicated().any():
            raise ValueError(f'{gauge}: duplicate observation timestamps.')
        for key, column in columns.items():
            if column not in frame:
                continue
            values = pd.to_numeric(frame[column], errors='coerce').to_numpy(float, copy=True)
            values[np.isin(values, [-999, -9999, -99, -99.9, -99.99])] = np.nan
            values[np.isfinite(values) & (values < 0)] = np.nan
            result[key] = pd.Series(values, index=index).reindex(dates).to_numpy()
    for key in ('SWE', 'SM'):
        if key not in result:
            result[key] = np.full(len(dates), np.nan)
    return {name: np.asarray(values, dtype=np.float32) for name, values in result.items()}

def names(c):
    """Return the selected observed-response names."""
    return c.get('observations', ['Q'])

def use_joint(c):
    """Select named channels for multivariate or standardized objectives."""
    return names(c) != ['Q'] or transform(c, 'Q') == 'zscore'

def normalize(c, evaluator, gauges, x, stop=None):
    """Resolve the joint response-loss normalization."""
    selected = names(c)
    if selected == ['Q'] and transform(c, 'Q') == 'none':
        return
    weights = np.array([c.get('observation_weights', {}).get(k, 1.0) for k in selected], float)
    if weights.sum() <= 0:
        weights[:] = 1
    weights /= weights.sum()
    counts = {k: 0 for k in selected}
    losses = {k: [] for k in selected}
    method = c.get('observation_normalization', 'auto')
    if len(selected) == 1:
        method = 'none'
    elif method == 'auto':
        method = 'benchmark' if c['loss'] in (1, 2, 5, 6) else 'none'
    if method == 'benchmark' and c['loss'] in (4, 7) and any((transform(c, k) == 'none' for k in selected)):
        raise ValueError('SAGE benchmark normalization is undefined for KGE/JKGE; use none, initial_loss or manual.')
    if c.get('_observation_scales'):
        return
    for index, gauge in enumerate(gauges):
        if stop is not None and stop.is_set():
            from engine import Cancelled
            raise Cancelled()
        data, _, kwargs = evaluator.prepare(gauge)
        observed = data['observations']
        sim = None
        if method == 'initial_loss':
            from crr_observables import simulate
            theta = physical(c, x[index])
            sim, _ = simulate(c['model'], theta, data['forcing'], **kwargs)
        for key in selected:
            y = observed[key]
            mask = data['train'] & np.isfinite(y) & (y >= 0)
            if not mask.any():
                continue
            standardized = transform(c, key) == 'zscore'
            if standardized and (not eligible(y[mask])):
                continue
            counts[key] += 1
            if method == 'benchmark' and standardized:
                losses[key].append(2.0)
                continue
            if method in ('benchmark', 'initial_loss'):
                q = np.full(len(y), np.mean(y[mask])) if method == 'benchmark' else sim[key]
                value, _ = full_loss(y, q, mask, dict(c, _loss_response=key), data['dates'])
                losses[key].append(value)
    if any((counts[k] == 0 for k in selected)):
        raise ValueError('No training observations for: ' + ', '.join((k for k in selected if counts[k] == 0)))
    scales = {}
    for key in selected:
        value = c.get('observation_scales', {}).get(key, 1.0) if method == 'manual' else np.mean(losses[key]) if method in ('benchmark', 'initial_loss') else 1.0
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f'Cannot prepare a positive {method} scale for {key}.')
        scales[key] = max(abs(float(value)), 1e-12)
    c['_observation_scales'] = scales
    c['_observation_counts'] = counts
    c['_observation_weights'] = dict(zip(selected, weights.tolist()))
    c['_observation_n_total'] = len(gauges)

def joint_loss(c, data, sim, jac, period='train'):
    """Combine response losses and their parameter gradients."""
    selected = names(c)
    entries = []
    for key in selected:
        y = data['observations'][key]
        mask = data[period] & np.isfinite(y) & (y >= 0)
        if not mask.any():
            continue
        if transform(c, key) == 'zscore' and (not eligible(y[mask])):
            continue
        try:
            value, delta = full_loss(y, sim[key], mask, dict(c, _loss_response=key), data['dates'])
        except UndefinedObjective:
            continue
        if not np.isfinite(value) or not np.isfinite(delta).all():
            continue
        entries.append((key, value, delta, jac[key]))
    if not entries:
        if period == 'train' and '_observation_counts' in c and any((transform(c, key) == 'zscore' for key in selected)):
            return (0.0, np.zeros(len(c['lower'])), [])
        return (float('nan'), np.full(len(c['lower']), np.nan), [])
    if len(selected) == 1:
        key = selected[0]
        coefficients = {key: c.get('_observation_n_total', 1) / c.get('_observation_counts', {}).get(key, 1)}
    elif '_observation_scales' in c:
        coefficients = {k: c['_observation_weights'][k] / c['_observation_scales'][k] * c['_observation_n_total'] / c['_observation_counts'][k] for k in selected}
    else:
        weights = {k: c.get('observation_weights', {}).get(k, 1.0) for k, _, _, _ in entries}
        total = sum(weights.values())
        coefficients = {k: (weights[k] / total if total > 0 else 1 / len(entries)) / c.get('observation_scales', {}).get(k, 1.0) for k in weights}
    value = sum((coefficients[k] * v for k, v, _, _ in entries))
    gradient = sum((coefficients[k] * (J.T @ delta) for k, _, delta, J in entries))
    return (value, gradient, [(k, coefficients[k], delta, J) for k, _, delta, J in entries])
