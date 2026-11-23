# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_print_view.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_print_view."""
import numpy as np
SCENARIOS = ('tt', 'te', 'et', 'ee')

def finite(values):
    """Finite."""
    return [float(v) for v in values if v is not None and np.isfinite(v)]

def scenario_population(records, scenario):
    """Keep a consistent requested basin cohort for every diagnostic."""
    roles = {'training'} if scenario[0] == 't' else {'held out', 'evaluation'}
    return [r for r in records if r.get('basin_role') in roles]

def failed_in_scenario(records, scenario):
    """Identify explicit integration failures, never infer from loss size."""
    return [r for r in scenario_population(records, scenario) if r.get('failed')]

def summarize(records, c):
    """Summarize."""
    from sage_observations import names
    result = {}
    for scenario in SCENARIOS:
        role = 'training' if scenario[0] == 't' else 'held out'
        period = 'train' if scenario[1] == 't' else 'evaluation'
        population = [r for r in records if r.get('basin_role') == role or (role == 'held out' and r.get('basin_role') == 'evaluation')]
        failures = failed_in_scenario(population, scenario)
        complete = bool(population) and (not failures)
        observations = {}
        for name in names(c):
            rows = [(r.get('observations', {}).get(name, {}).get(period, {}) if name != 'Q' else r.get(period, {})) if not r.get('failed') else {} for r in population]
            data = {}
            for metric in ('nse', 'kge', 'jkge', 's_fdc', 's_p', 's_logp', 'loss', 'rss'):
                values = finite((r.get(metric) for r in rows))
                data[metric] = float(np.sum(values)) if metric in ('loss', 'rss') and values else float(np.median(values)) if values else None
            for metric in ('nse', 'kge', 'jkge', 's_fdc', 's_p', 's_logp'):
                scores = finite((r.get(metric) for r in rows))
                data['sib_' + metric] = float(np.mean(1 - np.asarray(scores))) if scores else None
            nse = finite((r.get('nse') for r in rows))
            data['sib'] = float(np.mean(1 - np.asarray(nse))) if nse else None
            performance = {}
            for metric in ('loss', 'sar', 'gls', 'rss', 'huber', 'nse', 'kge', 'jkge', 's_fdc', 's_p', 's_logp', 'kge_r', 'kge_alpha', 'kge_beta', 'jkge_m', 'jkge_v', 'jkge_c'):
                sample = finite((r.get(metric) for r in rows))
                performance[metric] = dict(mean=float(np.mean(sample)) if sample else None, median=float(np.median(sample)) if sample else None)
            data['performance'] = performance
            data['n'] = len(nse)
            data['n_total'] = len(population)
            data['n_failed'] = len(failures)
            data['complete'] = complete
            if failures and (not c.get('omit_failed_basins')):
                data['loss'] = None
                data['rss'] = None
                for metric in data:
                    if metric == 'sib' or metric.startswith('sib_'):
                        data[metric] = None
            observations[name] = data
        total = []
        for r in population:
            if r.get('failed'):
                continue
            entries = []
            for name in names(c):
                value = (r.get(period, {}) if name == 'Q' else r.get('observations', {}).get(name, {}).get(period, {})).get('loss')
                if value is None or not np.isfinite(value):
                    continue
                if len(names(c)) == 1:
                    coefficient = c.get('_observation_n_total', 1) / c.get('_observation_counts', {}).get(name, 1)
                else:
                    coefficient = c.get('_observation_weights', {}).get(name, c.get('observation_weights', {}).get(name, 1) / sum((c.get('observation_weights', {}).get(k, 1) for k in names(c)))) / c.get('_observation_scales', {}).get(name, 1) * c.get('_observation_n_total', 1) / c.get('_observation_counts', {}).get(name, 1)
                entries.append(coefficient * value)
            if entries:
                total.append(sum(entries))
        result[scenario] = dict(observations=observations, loss=float(sum(total)) if total and (complete or c.get('omit_failed_basins')) else None, mean_loss=float(np.mean(total)) if total and (complete or c.get('omit_failed_basins')) else None, n_total=len(population), n_failed=len(failures), complete=complete, omit_failed_basins=bool(c.get('omit_failed_basins')))
    return result

def values(records, scenario, metric, observation='Q'):
    """Values."""
    role = 'training' if scenario[0] == 't' else 'held out'
    period = 'train' if scenario[1] == 't' else 'evaluation'
    rows = [r for r in records if r.get('basin_role') == role or (role == 'held out' and r.get('basin_role') == 'evaluation')]
    rows = [r for r in rows if not r.get('failed')]
    return finite(((r.get(period, {}) if observation == 'Q' else r.get('observations', {}).get(observation, {}).get(period, {})).get(metric) for r in rows))
