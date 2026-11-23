# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_quality.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_quality."""
import numpy as np

def _ratio(q, p):
    """Compute a data-coverage ratio."""
    return float(q / p) if np.isfinite(p) and p > 0 else None

def summarize(data, c, role='Training', documented_reason=''):
    """Summarize the quality of a basin input record."""
    f = np.asarray(data['forcing'])
    q = np.asarray(data['observed'])
    dates = np.asarray(data['dates'])
    p, ep, t = f.T[:3]
    train = np.asarray(data['train'], bool)
    evaluation = np.asarray(data['evaluation'], bool)
    finite = np.isfinite(q) & (q >= 0)
    bad = ~finite
    n = len(q)
    changes = np.diff(np.r_[False, bad, False].astype(int))
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    row = dict(gauge=str(data['gauge']), use=role, n_meteo=len(f), n_q=n, n_bad_q=int(bad.sum()), longest_q_gap=int(max(ends - starts, default=0)))
    for name, values in [('p', p), ('ep', ep), ('t', t), ('q', q)]:
        valid = np.isfinite(values) & (values >= 0 if name == 'q' else True)
        row['coverage_' + name] = 100 * float(valid.mean()) if len(valid) else 0.0
    complete = bool(np.isfinite(f).all())
    reasons = []
    insufficient = False
    constant = False
    operative = [('train', train)]
    if not (c['split'] == 'manual' and (not c.get('evaluation_enabled', True))):
        operative.append(('eval', evaluation))
    for name, mask in [('train', train), ('eval', evaluation)]:
        valid = mask & finite
        coverage = 100 * valid.sum() / mask.sum() if mask.any() else None
        variable = bool(valid.sum() >= 2 and np.ptp(q[valid]) > 0)
        row['coverage_q_' + name] = coverage
        row['discharge_variable_' + name] = variable
        pair = valid & np.isfinite(p)
        row['runoff_ratio_' + name] = _ratio(q[pair].sum(), p[pair].sum())
        if any((name == key for key, _ in operative)):
            insufficient |= coverage is None or coverage < 5
            constant |= not variable
    if not complete:
        reasons.append('incomplete forcing')
    row['n_negative_p'] = int(np.sum(p < 0))
    row['n_negative_ep'] = int(np.sum(ep < 0))
    if insufficient:
        reasons.append('less than 5% discharge in one or both periods')
    elif constant:
        reasons.append('constant discharge in one or both periods')
    pair = (train | evaluation) & finite & np.isfinite(p)
    mp = float(p[pair].mean()) if pair.any() else np.nan
    mq = float(q[pair].mean()) if pair.any() else np.nan
    ep_pair = pair & np.isfinite(ep)
    me = float(ep[ep_pair].mean()) if ep_pair.any() else np.nan
    years = dates.astype('datetime64[Y]').astype(int) + 1970
    months = dates.astype('datetime64[M]').astype(int) % 12 + 1
    water_year = years + (months >= 10)
    annual = []
    for year in np.unique(water_year[train | evaluation]):
        expected = (water_year == year) & (train | evaluation) & np.isfinite(p)
        valid = expected & finite
        if expected.any() and 100 * valid.sum() / expected.sum() >= 30 and (p[valid].sum() > 0):
            annual.append(float(q[valid].sum() / p[valid].sum()))
    alert = any((v is not None and (v < 0 or v > 1.25) for v in [row['runoff_ratio_train'], row['runoff_ratio_eval']]))
    row.update(forcing_complete=complete, discharge_complete=not (insufficient or constant), eligible=not reasons and (not documented_reason), reason=documented_reason or ' and '.join(reasons), documented_exclusion=bool(documented_reason), hydrologic_alert=alert, mean_p=mp, mean_ep=me, mean_q=mq, runoff_ratio=_ratio(mq, mp), aridity_ratio=_ratio(me, mp), water_balance_remainder=mp - mq, runoff_ratio_annual_mean=float(np.mean(annual)) if annual else None, runoff_ratio_annual_std=float(np.std(annual, ddof=1)) if len(annual) > 1 else 0.0 if annual else None, runoff_ratio_annual_years=len(annual), n_bad_forcing=int((~np.isfinite(f).all(axis=1)).sum()))
    return {key: None if isinstance(value, (float, np.floating)) and (not np.isfinite(value)) else value.item() if isinstance(value, np.generic) else value for key, value in row.items()}

def read(c, gauge, meteo_only=False):
    """Read the configured basin data and reference time window."""
    if c.get('_run_data_ready'):
        from sage_run_data import read as run_read
        return run_read(c, gauge)
    from sage_prepared_data import get, put
    cached = get(c, gauge)
    if cached is not None:
        cached['_prepared_cache_hit'] = True
        return cached
    from regional_data import read_region
    from sage_science import temporal_masks
    periods = {k: c[k] for k in ('train_start', 'train_end', 'eval_start', 'eval_end')}
    if c['split'] == 'manual' and (not c.get('evaluation_enabled', True)):
        periods.update(eval_start=c['train_start'], eval_end=c['train_end'])
    if c['split'] != 'manual':
        periods.update(train_start=c['common_start'], train_end=c['common_end'], eval_start=c['common_start'], eval_end=c['common_end'])
    data = read_region(c['data_root'], gauge, **{k: c[k] for k in ('region', 'sage_root', 'resolution', 'product', 'pet', 'precipitation', 'temperature', 'spinup', 'data_actions') if k in c}, **periods, allow_invalid=True, meteo_only=meteo_only, basin_selection=c.get('gauges', [gauge]))
    data['train'], data['evaluation'] = temporal_masks(data['dates'], c, precipitation=data['forcing'][:, 0])
    if not meteo_only:
        put(c, gauge, data)
    return data

def screen(c, progress=lambda *args: None, stop=None, reader_progress=None):
    """Screen basin inputs against the configured quality rules."""
    from hydrologic_data_actions import configure, exclusion_reason, normalize_gauge
    c = configure(c)
    from sage_schema import basin_partition
    from engine import Cancelled
    train, evaluation = basin_partition(c)
    training = set(train)
    rows = []
    gauges = c.get('excluded_gauges', [])
    reasons = c.get('exclusion_reason', 'documented regional data-quality exclusion')
    if isinstance(reasons, str):
        reasons = [reasons] * len(gauges)
    excluded = {normalize_gauge(g): reason for g, reason in zip(gauges, reasons)}
    for gauge in c['gauges']:
        reason = exclusion_reason(c['data_actions'], gauge)
        if reason:
            excluded[normalize_gauge(gauge)] = reason
    from xia_camels_us import matches
    separate = c['region'] == 'CAMELS_US' and c['resolution'] == 1 and (not matches(c['data_root']))
    total = len(c['gauges'])
    import threading
    from sage_input_loading import prepared_inputs

    def prepare(settings, gauge):
        """Prepare basin inputs and statistics for evaluation."""
        if normalize_gauge(gauge) in excluded:
            return {'_prepared_cache_hit': True}
        return read(settings, gauge, meteo_only=separate)
    inputs = prepared_inputs(c, stop or threading.Event(), prepare)
    for i, (gauge, data, read_error) in enumerate(inputs):
        if stop is not None and stop.is_set():
            raise Cancelled()
        role = 'Training' if i in training else 'Evaluation'
        if normalize_gauge(gauge) in excluded:
            rows.append(dict(gauge=gauge, use=role, eligible=False, reason=excluded[normalize_gauge(gauge)], documented_exclusion=True, data_actions=[a for a in c['data_actions'] if normalize_gauge(a['gauge_id']) == normalize_gauge(gauge)]))
            progress(gauge, i + 1, total)
            if reader_progress:
                reader_progress(7, gauge, i + 1, total)
            continue
        try:
            if read_error is not None:
                raise read_error
            if separate and (not data.get('_prepared_cache_hit', False)):
                from camels_us import read_discharge
                from sage_science import temporal_masks
                data = read_discharge(c['data_root'], gauge, data)
                data['train'], data['evaluation'] = temporal_masks(data['dates'], c, precipitation=data['forcing'][:, 0])
                from sage_prepared_data import put
                put(c, gauge, data)
            if '_run_data_directory' in c:
                from sage_observations import load, use_joint
                from sage_run_data import publish
                targets = load(c, gauge, data) if use_joint(c) else None
                publish(c, gauge, data, targets)
            row = summarize(data, c, role, excluded.get(normalize_gauge(gauge), ''))
            row['data_actions'] = [a for a in c['data_actions'] if normalize_gauge(a['gauge_id']) == normalize_gauge(gauge)]
        except (ValueError, FileNotFoundError, KeyError) as exc:
            row = dict(gauge=gauge, use=role, eligible=False, reason=str(exc), read_error=True)
        rows.append(row)
        progress(gauge, i + 1, total)
        if reader_progress:
            reader_progress(7, gauge, i + 1, total)
    if reader_progress and total:
        reader_progress(8, c['gauges'][-1], total, total)
    return rows

def active_configuration(c, rows):
    """Build the run configuration from eligible basins."""
    from sage_schema import basin_partition
    train, evaluation = basin_partition(c)
    eligible = {r['gauge'] for r in rows if r['eligible']}
    if len(eligible) == len(c['gauges']):
        return dict(c)
    nt = [c['gauges'][i] for i in train if c['gauges'][i] in eligible]
    ne = [c['gauges'][i] for i in evaluation if c['gauges'][i] in eligible]
    if not nt:
        raise ValueError('Data-quality screening left no eligible training basins.')
    active = dict(c, gauges=nt + ne, basin_sampling='file', basin_split_ids=nt + ne, basin_train_fraction=len(nt) / (len(nt) + len(ne)))
    return active
