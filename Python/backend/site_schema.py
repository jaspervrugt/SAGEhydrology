# ==========================================================================
# SAGEhydrology | Python source
# File: backend/site_schema.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from site_schema."""
from copy import deepcopy
from pathlib import Path
import re
from engine import ROOT, defaults
MODEL_NAMES = {1: 'HYMOD', 2: 'HMODEL', 3: 'SAC-SMA', 4: 'Xinanjiang', 5: 'GR4J-A', 6: 'HBV', 7: 'CFE-NWM', 99: 'user_model'}

LOSSES = {1: 'Sum of absolute residuals', 2: 'Generalized least squares / RSS', 3: 'Nash–Sutcliffe efficiency', 4: 'Kling–Gupta efficiency', 5: 'Huber robust estimator', 6: 'Flow-duration curve divergence', 7: 'Jawad–Kling–Gupta efficiency'}

def region_catalogue(sage_root):
    """Region catalogue."""
    from regional_data import region_registry
    return region_registry(str(Path(sage_root).resolve()))

def session_defaults():
    """Session defaults."""
    old = defaults()
    from sage_paths import SAGE_RESOURCES, SOFTWARE_ROOT
    sage = str(SAGE_RESOURCES)
    return dict(schema_version=2, notes='', software_root=str(SOFTWARE_ROOT), site_root=str(ROOT), sage_root=sage, data_root=old['data_root'], output_root=str(ROOT.parent / 'results' / 'Python'), region='CAMELS_US', resolution=1, split='manual', train_start=old['train_start'], train_end=old['train_end'], eval_start=old['eval_start'], eval_end=old['eval_end'], spinup=365, common_start='1989-10-01', common_end='2008-09-30', block_days=365, train_fraction=0.8, split_seed=1, folds=5, fold=1, product='nldas', pet=1, model=1, model_parameters=model_catalogue(sage), lower=old['lower'], upper=old['upper'], y0=old['y0'], settings=old['settings'], config=old['config'], loss=2, jkge_m=2, jkge_window=31, jkge_operator=1, fdc_formulation=1, optimizer=5, starts=10, max_iter=1000, seed=20261012, checkpoint_every=1, parallel=False, workers=1, max_pending=2, backend='cpp', finisher='lbfgsb', verbose_progress=False, min_per_zone=5, max_zones=8, optimizer_options=dict(memory=5, tol_pg=1e-06, tol_rel_f=1e-06, buffer=1e-06, c1=0.0001, ls_rho=0.5, ls_max=40, alpha0=1.0), gauges=['01022500'], diagnostics=dict(train=True, evaluation=True, hydrographs=True, fdc=True, convergence=True, attribution=False, selected_count=5))

def model_catalogue(sage_root, resolution=1):
    """Model catalogue."""
    from runtime_metadata import load
    data = load(str(sage_root))
    if data is None:
        return dict()
    result = {k: deepcopy(v) for k, v in data['models'].items() if str(k) not in {'11', '12', '13', '14', '15', '16', '17', '20'}}
    if resolution > 1:
        for meta in result.values():
            for parameter in meta['parameters']:
                unit = parameter['unit'].lower().replace(' ', '')
                factor = resolution if unit in ('d', 'day', 'days') else 1 / resolution if '/d' in unit or unit in ('d^-1', 'd-1', 'day^-1', 'day-1') else 1
                for key in ('lower', 'upper'):
                    parameter[key] *= factor
                token = 'h' if resolution == 24 else '15 min'
                if factor != 1:
                    parameter['unit'] = token if factor == resolution else re.sub('/days?(?=$|/)|/d(?=$|/)', '/' + token, parameter['unit'])
    result = {key: value for key, value in result.items() if str(key) in set("1234567")}
    from user_model_native import directory, metadata

    folder = directory()
    if folder is not None:
        result['99'] = {
            key: value for key, value in metadata(folder, resolution).items()
            if key != 'definition'
        }
    return result
