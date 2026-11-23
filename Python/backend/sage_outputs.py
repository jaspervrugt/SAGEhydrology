# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_outputs.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_outputs."""
from sage_parameter_scaling import derivative
from time import perf_counter
import numpy as np
from engine import score
from sage_science import full_loss, attribution, diagnostic_scores

def prediction(c, data, kwargs, theta):
    """Prediction."""
    model_config = kwargs.get('model_config') or {}
    import crr_native
    return crr_native.simulate(c['model'], theta, data['forcing'], **kwargs)[0]

def basin_output(c, gauge, data, kwargs, theta, role, history=(), started=None, *, compute_attribution=True):
    """Basin output."""
    import crr_native
    from sage_observations import names, joint_loss, use_joint
    named = use_joint(c)
    if named:
        from crr_observables import simulate
        simulated, jac = simulate(c['model'], theta, data['forcing'], **kwargs)
        q = simulated['Q']
        J = jac['Q']
    elif not compute_attribution:
        q, J = (prediction(c, data, kwargs, theta), None)
    else:
        q, J = crr_native.simulate(c['model'], theta, data['forcing'], **kwargs)
    if not np.isfinite(q).all():
        raise ValueError(f'{gauge}: nonfinite final discharge.')
    totals = []
    nets = []
    scores = {}
    for period in ('train', 'evaluation'):
        mask = data[period] & np.isfinite(data['observed']) & (data['observed'] >= 0)
        if J is not None and mask.any():
            try:
                _, delta = full_loss(data['observed'], q, mask, c, data['dates'])
                total, net = attribution(J, delta, c['lower'], c['upper'], parameter_derivative=derivative(c, theta))
            except ValueError:
                total = net = np.full(len(theta), np.nan)
        else:
            total = net = np.full(len(theta), np.nan)
        totals.append(total)
        nets.append(net)
        scores[period] = score(data['observed'], q, data[period])
        scores[period].update(diagnostic_scores(data['observed'], q, data[period], c, data['dates']))
    result = dict(gauge=gauge, model=c['model'], theta=np.asarray(theta).tolist(), **scores, seconds=perf_counter() - started if started is not None else 0.0, status='Shared network', basin_role=role, series_file=f'{gauge}.npz', trials=[dict(trajectory=[h['objective'] for h in history])], data_quality=dict(samples=len(q), missing_observations=int(np.sum(~np.isfinite(data['observed']) | (data['observed'] < 0))), forcing_nonfinite=int(np.sum(~np.isfinite(data['forcing'])))))
    series = dict(dates=data['dates'], forcing=data['forcing'], observed=data['observed'], simulated=q, train=data['train'], evaluation=data['evaluation'], attribution_total=totals, attribution_net=nets)
    if named:
        result['observations'] = {}
        result['joint_objective'] = {}
        joint_totals = []
        joint_nets = []
        for name in names(c):
            y = data['observations'][name]
            predicted = simulated[name]
            series['observed_' + name] = y
            series['simulated_' + name] = predicted
            result['observations'][name] = {}
            for period in ('train', 'evaluation'):
                metrics = score(y, predicted, data[period])
                metrics.update(diagnostic_scores(y, predicted, data[period], dict(c, _loss_response=name), data['dates']))
                try:
                    metrics['loss'] = float(full_loss(y, predicted, data[period], dict(c, _loss_response=name), data['dates'])[0])
                except ValueError:
                    metrics['loss'] = None
                result['observations'][name][period] = metrics
            per_total = []
            per_net = []
            for period in ('train', 'evaluation'):
                try:
                    _, delta = full_loss(y, predicted, data[period], dict(c, _loss_response=name), data['dates'])
                    total, net = attribution(jac[name], delta, c['lower'], c['upper'], parameter_derivative=derivative(c, theta))
                except ValueError:
                    total = net = np.full(len(theta), np.nan)
                per_total.append(total)
                per_net.append(net)
            series['attribution_total_' + name] = per_total
            series['attribution_net_' + name] = per_net
        for period in ('train', 'evaluation'):
            try:
                value, _, parts = joint_loss(c, data, simulated, jac, period)
                if not parts:
                    raise ValueError('No scored observations in this period.')
                result['joint_objective'][period] = float(value)
                contribution = sum((weight * jacobian * delta[:, None] for _, weight, delta, jacobian in parts)) * derivative(c, theta)
                joint_totals.append(np.sum(np.abs(contribution), axis=0))
                joint_nets.append(np.abs(np.sum(contribution, axis=0)))
            except ValueError:
                result['joint_objective'][period] = None
                joint_totals.append(np.full(len(theta), np.nan))
                joint_nets.append(np.full(len(theta), np.nan))
        series['attribution_total'] = joint_totals
        series['attribution_net'] = joint_nets
    return (result, series)
