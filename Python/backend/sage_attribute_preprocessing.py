# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_attribute_preprocessing.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_attribute_preprocessing."""
import numpy as np
METHODS = {'zscore': 'Z-score standardization', 'minmax01': 'Min–max [0, 1]', 'minmax11': 'Min–max [−1, 1]', 'robust': 'Robust scaling: median / IQR', 'log_zscore': 'Log1p + z-score (nonnegative attributes)', 'power': 'Yeo–Johnson power transform + z-score', 'quantile_uniform': 'Quantile transform: uniform [0, 1]', 'quantile_normal': 'Quantile transform: normal', 'none': 'None (original values)'}

def settings(c):
    """Resolve the effective attribute-preprocessing settings."""
    if c.get('attribute_preprocessing', 'training') != 'custom':
        return c.get('attribute_preprocessing', 'training')
    return dict(method=c.get('attribute_scaling', 'zscore'), cohort=c.get('attribute_fit_cohort', 'training'), missing=c.get('attribute_missing', 'reject'), attributes=c.get('attributes', []), overrides=c.get('attribute_methods', {}) if c.get('attribute_methods_enabled', bool(c.get('attribute_methods', {}))) else {})

def nonlinear(values, state):
    """Apply the fitted nonlinear attribute transformation."""
    method = state['method']
    if method == 'log_zscore':
        if np.any(values < 0):
            raise ValueError('Log1p preprocessing requires nonnegative attribute values.')
        return np.log1p(values)
    if method == 'power':
        from scipy.stats import yeojohnson
        return yeojohnson(values, lmbda=state['lambda'])
    if method.startswith('quantile_'):
        result = np.interp(values, state['knots'], state['probabilities'])
        if method == 'quantile_normal':
            from scipy.special import ndtri
            result = ndtri(np.clip(result, 1e-07, 1 - 1e-07))
        return result
    return values

def fit(raw, training_rows, options):
    """Fit attribute transformations and transform the basin attributes.

    Parameters
    ----------
    raw : array-like, shape (basins, attributes)
        Original catchment attributes; columns retain their selected order.
    training_rows : array-like of int
        Row indices of training basins, used for fitting and imputation.
    options : dict
        Scaling method, fitting cohort, missing-value policy, attribute
        names
        and optional per-attribute overrides.

    Returns
    -------
    transformed : ndarray
        Transformed attributes for all supplied basins.
    state : dict
        Serializable fitted statistics for reuse on evaluation or new
        basins.
    """
    x = np.asarray(raw, float)
    train = np.asarray(training_rows, int)
    if x.ndim != 2 or not len(train):
        raise ValueError('Attribute preprocessing requires a matrix and training basins.')
    if options.get('cohort') not in ('training', 'selected'):
        raise ValueError('Choose training basins or all selected basins for preprocessing.')
    if options.get('missing') not in ('reject', 'median'):
        raise ValueError('Choose a supported missing-value policy.')
    if options['missing'] == 'reject' and (not np.isfinite(x).all()):
        raise ValueError('Selected attributes contain missing/nonfinite values; choose median imputation or complete attributes.')
    if options.get('method') not in METHODS or not isinstance(options.get('overrides', {}), dict):
        raise ValueError('Choose a supported attribute scaling method and valid per-attribute settings.')
    training = x[train]
    if np.any(~np.isfinite(training).any(axis=0)):
        raise ValueError('Every attribute needs at least one finite training-basin value.')
    median = np.nanmedian(np.where(np.isfinite(training), training, np.nan), axis=0)
    filled = np.where(np.isfinite(x), x, median)
    fitted = filled[train] if options['cohort'] == 'training' else filled
    columns = []
    names = options.get('attributes', [])
    for j in range(x.shape[1]):
        name = names[j] if j < len(names) else str(j + 1)
        method = options.get('overrides', {}).get(name, options['method'])
        if method not in METHODS:
            raise ValueError(f'Unknown preprocessing method for {name}: {method}')
        values = fitted[:, j]
        state = dict(method=method, center=0.0, scale=1.0, name=name)
        if method == 'power':
            from scipy.stats import yeojohnson_normmax
            state['lambda'] = float(yeojohnson_normmax(values)) if np.ptp(values) > 0 else 1.0
        if method.startswith('quantile_'):
            knots, first, counts = np.unique(np.sort(values), return_index=True, return_counts=True)
            ranks = (first + (counts - 1) / 2) / max(1, len(values) - 1)
            if len(knots) == 1:
                ranks = np.array([0.5])
            state.update(knots=knots.tolist(), probabilities=ranks.tolist())
        transformed = nonlinear(values, state)
        if method in ('zscore', 'log_zscore', 'power'):
            state.update(center=float(transformed.mean()), scale=float(transformed.std(ddof=1)) if len(values) > 1 else 1.0)
        elif method in ('minmax01', 'minmax11'):
            state.update(center=float(values.min()), scale=float(np.ptp(values)))
        elif method == 'robust':
            q1, q3 = np.percentile(values, [25, 75])
            state.update(center=float(np.median(values)), scale=float(q3 - q1))
        if not np.isfinite(state['scale']) or state['scale'] == 0:
            state['scale'] = 1.0
        columns.append(state)
    scaler = dict(version=2, mode='custom', missing=options['missing'], cohort=options['cohort'], median=median.tolist(), mean=[v['center'] for v in columns], scale=[v['scale'] for v in columns], columns=columns)
    return (apply(x, scaler), scaler)

def apply(raw, scaler):
    """Apply saved attribute transformations to new basin attributes."""
    x = np.asarray(raw, float)
    if scaler.get('mode') == 'matlab' or scaler.get('missing') == 'reject':
        if not np.isfinite(x).all():
            raise ValueError('Saved attribute preprocessing rejects missing/nonfinite values.')
    filled = np.where(np.isfinite(x), x, scaler['median'])
    if 'columns' not in scaler:
        return (filled - np.asarray(scaler['mean'])) / np.asarray(scaler['scale'])
    if filled.shape[1] != len(scaler['columns']):
        raise ValueError('Saved attribute preprocessing has a different attribute count.')
    result = np.empty_like(filled)
    for j, state in enumerate(scaler['columns']):
        v = (nonlinear(filled[:, j], state) - state['center']) / state['scale']
        if state['method'] == 'minmax11':
            v = 2 * v - 1
        result[:, j] = v
    if not np.isfinite(result).all():
        raise ValueError('Attribute preprocessing produced nonfinite values.')
    return result

def configuration_summary(c):
    """
    Describe the effective setup, including active per-attribute
    overrides.
    """
    numbers = {key: i for i, key in enumerate(METHODS, 1)}

    def label(key):
        """Format a numbered preprocessing method label."""
        return f'{numbers[key]} {METHODS[key]}'
    effective = settings(c)
    mode = c.get('attribute_preprocessing', 'training')
    if isinstance(effective, str):
        matlab = effective == 'matlab'
        effective = dict(method='zscore', cohort='selected' if matlab else 'training', missing='reject' if matlab else 'median', overrides={})
    method = effective['method']
    rows = [('Preset', {'matlab': 'MATLAB compatibility', 'training': 'Training basins only (legacy)', 'custom': 'Custom'}.get(mode, mode)), ('Global method', label(method)), ('Fit cohort', 'All selected training + evaluation basins' if effective['cohort'] == 'selected' else 'Training basins only'), ('Missing vals', 'Reject missing/nonfinite values' if effective['missing'] == 'reject' else 'Fill with training-basin median')]
    if method in ('zscore', 'log_zscore', 'power'):
        rows.append(('Std deviation', 'Population (N)' if mode == 'training' else 'Sample (N - 1)'))
    overrides = {key: value for key, value in effective.get('overrides', {}).items() if key in c.get('attributes', [])}
    rows.append(('Per attribute', 'Enabled' if mode == 'custom' and c.get('attribute_methods_enabled', bool(c.get('attribute_methods', {}))) else 'Disabled'))
    if overrides:
        rows.append(('Other attrs', '0 Global method'))
        rows.extend(((str(name), label(overrides[name])) for name in c.get('attributes', []) if name in overrides))
    else:
        rows.append(('All attrs', 'Use global method'))
    return rows
