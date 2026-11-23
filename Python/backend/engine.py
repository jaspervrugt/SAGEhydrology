# ==========================================================================
# SAGEhydrology | Python source
# File: backend/engine.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from engine."""
import json
import math
from pathlib import Path
import sys
from time import sleep
import numpy as np
ROOT = Path(getattr(sys, '_MEIPASS', Path(__file__).resolve().parents[1]))
from sage_paths import SOFTWARE_ROOT
from site_hymod import LOWER, UPPER

class Cancelled(Exception):
    pass

def defaults():
    """Defaults."""
    return dict(data_root=str(SOFTWARE_ROOT / 'Data' / 'CAMELS_US'), output_root=str(ROOT / 'gui_python' / 'runs'), product='nldas', pet=1, backend='cpp', gauges=['01022500'], starts=10, max_iter=1000, seed=20261012, train_start='1999-10-01', train_end='2008-09-30', eval_start='1989-10-01', eval_end='1999-09-30', spinup=365, lower=LOWER.tolist(), upper=UPPER.tolist(), y0=[1e-05] * 7, parallel=False, workers=1, max_pending=2, verbose_progress=True, settings=[0.01, 1.0, 0.0001, 0.001, 0.001, 2.0, 10000.0], config=[1.0, 1e-06, 0.01])

def validate(config):
    """Validate."""
    if int(config.get('model', 1)) not in (*range(1, 8), 99):
        raise ValueError('Supported models are public models 1-7 and user_model (99)')
    c = dict(config)
    from basin_pool import cpu_inventory
    c['parallel'] = bool(c.get('parallel', False))
    c['workers'] = int(c.get('workers', 1))
    c['max_pending'] = int(c.get('max_pending', 2))
    if not 1 <= c['workers'] <= cpu_inventory()['logical']:
        raise ValueError('Worker count must be within the available logical CPU count.')
    if c['max_pending'] < 1:
        raise ValueError('Maximum pending jobs must be positive.')
    model = int(c.get('model', 1))
    c['model'] = model
    from crr_native import PARAM_COUNTS, STATE_COUNTS
    if model == 99:
        from user_model_native import metadata
        folder = c.get('model_config', {}).get('_user_model_folder')
        user = metadata(folder, c.get('resolution', 1))
        d, m = (len(user['parameters']), user['states'])
    else:
        if model not in PARAM_COUNTS:
            raise ValueError('Unknown CRR model')
        d = PARAM_COUNTS[model]
        m = STATE_COUNTS[model]
    if c['backend'] == 'numba' and model != 1:
        raise ValueError('Numba backend supports HYMOD only; select native C++.')
    if model == 7:
        cfg = c.get('model_config', {})
        order = cfg.get('K', 3)
        ordinates = cfg.get('giuh_ordnts', [0.1, 0.35, 0.2, 0.14, 0.1, 0.06, 0.05])
        if not isinstance(order, (int, float)) or not np.isfinite(order) or order < 1 or (int(order) != order):
            raise ValueError('CFE Nash order K must be an integer >= 1.')
        m = 4 + int(order) + len(ordinates)
    for field in ('starts', 'max_iter', 'seed', 'spinup'):
        c[field] = int(c[field])
    if c['starts'] < 1 or c['max_iter'] < 1 or c['spinup'] < 0 or (c['seed'] < 0):
        raise ValueError('Starts/iterations must be positive; spin-up and seed nonnegative.')
    c['pet'] = int(c['pet'])
    region = c.get('region', 'CAMELS_US')
    if c['backend'] not in ('cpp', 'numba') or (region == 'CAMELS_US' and (c['product'] not in ('nldas', 'daymet', 'maurer') or c['pet'] not in (1, 2, 3))):
        raise ValueError('Invalid backend, forcing product, or PET selection.')
    if not Path(c['data_root']).is_dir():
        raise ValueError('Choose the installed regional data directory.')
    c['gauges'] = [str(g).zfill(8) if region == 'CAMELS_US' else str(g) for g in c['gauges']]
    if region == 'CAMELS_US' and any((len(g) != 8 or not g.isdigit() for g in c['gauges'])):
        raise ValueError('Select valid eight-digit CAMELS-US gauges.')
    if not c['gauges'] or any((not g or any((ch in g for ch in '/\\:')) or g in ('.', '..') for g in c['gauges'])):
        raise ValueError('Select at least one valid regional gauge.')
    if region != 'CAMELS_US':
        from regional_data import reader_issues
        issues = reader_issues(c['sage_root'], region, c.get('resolution', 1))
        if issues:
            raise ValueError('\n'.join(issues))
    if len(c['gauges']) != len(set(c['gauges'])):
        raise ValueError('Gauge selection contains duplicates.')
    for a, b in (('train_start', 'train_end'), ('eval_start', 'eval_end')):
        from datetime import date
        start, end = (date.fromisoformat(c[a]), date.fromisoformat(c[b]))
        if end < start:
            raise ValueError('Period end must follow period start.')
    for field, size in (('lower', d), ('upper', d), ('y0', m), ('settings', 7), ('config', 3)):
        values = np.asarray(c[field], float)
        if values.shape != (size,) or not np.all(np.isfinite(values)):
            raise ValueError(f'{field} must contain {size} finite numbers.')
        c[field] = values.tolist()
    lo, hi = (np.array(c['lower']), np.array(c['upper']))
    if np.any(lo >= hi):
        raise ValueError('Every parameter lower bound must be below its upper bound.')
    if model == 1 and (lo[0] <= 0 or lo[1] <= 0 or np.any(lo[3:5] <= 0) or (lo[6] <= 0) or (lo[2] < 0) or (hi[2] > 1)):
        raise ValueError('Capacity, beta, recessions, and melt factor must be positive; alpha must be in [0,1].')
    if min(c['config']) <= 0:
        raise ValueError('Smoothing settings must be positive.')
    s = c['settings']
    if min(s[:5]) <= 0 or s[2] > s[1] or s[5] <= 0 or (s[6] < 1) or (s[6] != int(s[6])):
        raise ValueError('Invalid solver settings or rejection limit.')
    if np.any(np.array(c['y0'])[-9 if model in (14, 17) else -6 if model in (15, 16) else 0:] < 0) or (model == 1 and c['y0'][1] > lo[0]):
        raise ValueError('Initial stores must be nonnegative; soil storage must not exceed the minimum capacity.')
    if c['backend'] == 'cpp' and model not in (18, 19, 20, 99):
        from build_crr import library_path
        if not library_path().is_file():
            raise ValueError('Build the C++ library on the Paths tab, or select Numba.')
    from hydrologic_data_actions import configure
    return configure(c)

def json_scalar(value):
    """Convert NumPy numeric scalars at strict JSON boundaries."""
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    raise TypeError(f'Object of type {type(value).__name__} is not JSON serializable')

from tempfile import NamedTemporaryFile

def checkpoint_locked(exc):
    return isinstance(exc, PermissionError) or getattr(exc, "winerror", None) in (32, 33)

def atomic_json(path, data, *, attempts=100):
    path = Path(path)
    payload = json.dumps(
        data, indent=2, allow_nan=False, default=json_scalar
    ) + "\n"
    temporary = None
    staged = False
    try:
        for attempt in range(attempts):
            try:
                if not staged:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with NamedTemporaryFile(
                        mode="w", encoding="utf-8", dir=path.parent,
                        prefix=f".{path.name}.", suffix=".tmp", delete=False,
                    ) as stream:
                        temporary = Path(stream.name)
                        stream.write(payload)
                    staged = True
                temporary.replace(path)
                return
            except OSError as exc:
                if not (checkpoint_locked(exc) or isinstance(exc, FileNotFoundError)):
                    raise
                if attempt == attempts - 1:
                    raise
                if not staged or isinstance(exc, FileNotFoundError):
                    if temporary is not None:
                        temporary.unlink(missing_ok=True)
                    temporary = None
                    staged = False
                sleep(0.15)
    finally:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except PermissionError:
                pass

def score(y, q, mask):
    """Score."""
    mask = mask & np.isfinite(y) & np.isfinite(q)
    y, q = (y[mask], q[mask])
    if len(y) < 2:
        return dict(nse=None, kge=None, rss=float(np.sum((y - q) ** 2)) if len(y) else None, count=len(y))
    rss = float(np.sum((y - q) ** 2))
    tss = float(np.sum((y - y.mean()) ** 2))
    nse = 1 - rss / tss if tss > 0 else None
    kge = None
    if y.std() > 0 and q.std() > 0 and (y.mean() != 0):
        corr = float(np.corrcoef(y, q)[0, 1])
        ratio = float(q.std() / y.std())
        bias = float(q.mean() / y.mean())
        kge = 1 - math.sqrt((corr - 1) ** 2 + (ratio - 1) ** 2 + (bias - 1) ** 2)
    return dict(nse=nse, kge=kge, rss=rss, count=len(y))
