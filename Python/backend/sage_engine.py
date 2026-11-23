# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_engine.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_engine."""
from sage_parameter_scaling import physical, derivative
from contextlib import closing
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from collections import OrderedDict
from copy import deepcopy
import multiprocessing as mp
import os
from pathlib import Path
from time import perf_counter
import traceback
import json
import numpy as np
from engine import checkpoint_locked, atomic_json, score, Cancelled
from sage_network import Network, standardize
from sage_schema import worker_configuration, attribute_matrix, basin_partition
from sage_science import full_loss, temporal_masks, attribution, diagnostic_scores
from sage_checkpoint import save_npz, save_training, load_training
from sage_descent import update as descent_update
from sage_observations import load as load_observations, joint_loss, normalize, use_joint
from sage_outputs import basin_output, prediction
_EVALUATOR = None

class Evaluator:

    def __init__(self, c, elevation_inputs=None):
        """Initialize the evaluator or network state."""
        self.c = c
        self.elevation_inputs = elevation_inputs or {}
        self.cache = OrderedDict()

    def prepare(self, gauge):
        """Prepare basin inputs and statistics for evaluation."""
        if gauge in self.cache:
            self.cache.move_to_end(gauge)
        else:
            while len(self.cache) >= 8:
                self.cache.popitem(last=False)
            c = self.c
            periods = {k: c[k] for k in ('train_start', 'train_end', 'eval_start', 'eval_end')}
            if c['split'] == 'manual' and (not c.get('evaluation_enabled', True)):
                periods.update(eval_start=c['train_start'], eval_end=c['train_end'])
            if c['split'] != 'manual':
                periods.update(train_start=c['common_start'], train_end=c['common_end'], eval_start=c['common_start'], eval_end=c['common_end'])
            from sage_quality import read as read_prepared
            data = read_prepared(c, gauge)
            from cartplus_basin_inputs import prepare_cartplus_bands
            if gauge in self.elevation_inputs:
                data.update(self.elevation_inputs[gauge])
            else:
                prepare_cartplus_bands(c, gauge, data)
            forcing = np.asarray(data['forcing'])
            bad_forcing = ~np.isfinite(forcing).all(axis=1)
            if bad_forcing.any():
                first = np.flatnonzero(bad_forcing)[0]
                raise ValueError(f"Gauge {gauge}: missing or invalid meteorological forcing on {data['dates'][first]}; P, PET, T = {forcing[first].tolist()}")
            if not np.asarray(data['train']).any():
                raise ValueError('No valid training discharge')
            data['train'], data['evaluation'] = temporal_masks(data['dates'], c, precipitation=data['forcing'][:, 0])
            targets = {key.removeprefix('_target_'): value for key, value in data.items() if key.startswith('_target_')}
            data['observations'] = targets or (load_observations(c, gauge, data) if use_joint(c) else dict(Q=data['observed']))
            cfg = dict(c.get('model_config', {}))
            cfg.update(T_sm=c['config'][0], eps_m=c['config'][1], rho=c['config'][2], routing_length=60 * c['resolution'])
            kwargs = dict(y0=c['y0'], model_config=cfg, settings=c['settings'], routing_start=c['spinup'] * c['resolution'])

            def fun(x):
                """Evaluate the local objective callback."""
                if use_joint(c):
                    from crr_observables import simulate
                    theta = physical(c, x)
                    dpdx = derivative(c, theta)
                    simulated, jac = simulate(c['model'], theta, data['forcing'], **kwargs)
                    data['_predicted'] = simulated
                    value, gradient, parts = joint_loss(c, data, simulated, jac)
                    contribution = sum((weight * jacobian * delta[:, None] for _, weight, delta, jacobian in parts), np.zeros((len(data['forcing']), len(x)))) * dpdx
                    data['_attribution'] = (np.sum(np.abs(contribution), axis=0), np.abs(np.sum(contribution, axis=0)))
                    return (value, gradient * dpdx)
                import crr_native
                theta = physical(c, x)
                dpdx = derivative(c, theta)
                q, J = crr_native.simulate(c['model'], theta, data['forcing'], **kwargs)
                data['_predicted'] = {'Q': q}
                mask = data['train'] & np.isfinite(data['observed']) & (data['observed'] >= 0)
                value, delta = full_loss(data['observed'], q, mask, c, data['dates'])
                from sage_science import attribution
                data['_attribution'] = attribution(J, delta, c['lower'], c['upper'], parameter_derivative=dpdx)
                return (value, J.T @ delta * dpdx)
            self.cache[gauge] = (data, fun, kwargs)
        return self.cache[gauge]

    def evaluate(self, gauge, x):
        """Evaluate a basin objective and its analytic gradient."""
        from simulation_failure import IntegrationFailure
        data, fun, _ = self.prepare(gauge)
        try:
            f, g = fun(x)
            if not np.isfinite(f) or not np.all(np.isfinite(g)):
                raise IntegrationFailure('Nonfinite objective or gradient')
            return (f, g)
        except IntegrationFailure as exc:
            data.pop('_predicted', None)
            data.pop('_attribution', None)
            error = IntegrationFailure(f'{gauge}: {exc}')
            error.gauge = gauge
            error.normalized_parameters = np.asarray(x).tolist()
            error.physical_parameters = physical(self.c, x).tolist()
            raise error from exc

    def failed_record(self, gauge, x, reason):
        """Failed record."""
        data, _, _ = self.prepare(gauge)
        data.pop('_predicted', None)
        data.pop('_attribution', None)
        periods = {p: dict(loss=None) for p in ('train', 'evaluation')}
        record = dict(gauge=gauge, failed=True, failure_reason=reason, normalized_parameters=np.asarray(x).tolist(), physical_parameters=physical(self.c, x).tolist(), observations={name: periods.copy() for name in data['observations']}, **periods)
        return (None, None, record)

    def evaluate_with_metrics(self, gauge, x, predict_only=False, omit_failure=False, *, postprocess=False):
        """Evaluate with metrics."""
        from simulation_failure import IntegrationFailure
        try:
            return self._evaluate_with_metrics(gauge, x, predict_only, postprocess)
        except IntegrationFailure as exc:
            if not predict_only and (not omit_failure):
                raise
            return self.failed_record(gauge, x, str(exc))

    def _evaluate_with_metrics(self, gauge, x, predict_only, postprocess):
        """Evaluate a basin and compute its diagnostic metrics."""
        if predict_only:
            data, _, kwargs = self.prepare(gauge)
            theta = physical(self.c, x)
            if use_joint(self.c):
                from crr_observables import simulate
                data['_predicted'], _ = simulate(self.c['model'], theta, data['forcing'], **kwargs)
            else:
                if postprocess:
                    predicted = prediction(self.c, data, kwargs, theta)
                else:
                    import crr_native
                    predicted, _ = crr_native.simulate(self.c['model'], theta, data['forcing'], **kwargs)
                data['_predicted'] = {'Q': predicted}
            f, g = (0.0, np.zeros(len(theta)))
        else:
            f, g = self.evaluate(gauge, x)
            data, _, _ = self.prepare(gauge)
        metrics = {}
        observations = {}
        for name, observed in data['observations'].items():
            predicted = data['_predicted'][name]
            values = {}
            for period in ('train', 'evaluation'):
                values[period] = score(observed, predicted, data[period])
                values[period].update(diagnostic_scores(observed, predicted, data[period], dict(self.c, _loss_response=name), data['dates'], data.setdefault('_diagnostic_baselines', {}).setdefault((name, period), {})))
                try:
                    values[period]['loss'] = float(full_loss(observed, predicted, data[period], dict(self.c, _loss_response=name), data['dates'])[0])
                except ValueError:
                    values[period]['loss'] = None
            if name == 'Q':
                metrics.update(values)
            observations[name] = values
        record = dict(gauge=gauge, **metrics, observations=observations)
        if not predict_only and '_attribution' in data:
            record['_attribution'] = [v.tolist() for v in data['_attribution']]
        selection = self.c.get('diagnostic_selection', {}).get('selection', {})
        if selection.get(gauge):
            record['_series'] = dict(dates=data['dates'], forcing=data['forcing'], observed=data['observed'], simulated=data['_predicted']['Q'], train=data['train'], evaluation=data['evaluation'])
            for name, y in data['observations'].items():
                record['_series']['observed_' + name] = y
                record['_series']['simulated_' + name] = data['_predicted'][name]
        return (f, g, record)

def stop_failed_training(out, report, iteration, count):
    """Record a fatal empty training cohort and block report export."""
    message = f'Iteration {iteration}: all {count} training basins failed integration. No valid loss or gradient remains. Training stopped; figures and PPTX are disabled.'
    report['postprocessing_blocked'] = message
    report['postprocessing_complete'] = False
    with (out / 'integration_failure_notes.txt').open('a', encoding='utf-8') as stream:
        stream.write(message + '\n')
    raise ValueError(message)

def document_failures(out, report, iteration, training, evaluation):
    """Append reproducible numerical failure details without fake losses."""
    events = []
    for rows, role in ((training, 'training'), (evaluation, 'held out')):
        for row in rows:
            record = row[2]
            if record.get('failed'):
                events.append(dict(iteration=iteration, basin_role=role, gauge=record['gauge'], reason=record['failure_reason'], normalized_parameters=record['normalized_parameters'], physical_parameters=record['physical_parameters']))
    if not events:
        return
    report.setdefault('numerical_failures', []).extend(events)
    with (out / 'integration_failures.jsonl').open('a', encoding='utf-8') as stream:
        for event in events:
            stream.write(json.dumps(event, allow_nan=False) + '\n')
    with (out / 'integration_failure_notes.txt').open('a', encoding='utf-8') as stream:
        for event in events:
            stream.write(f"Iteration {iteration}; basin {event['gauge']} ({event['basin_role']}): {event['reason']}. Loss/metrics missing; excluded for this iteration.\n")

def _initialize(c, elevation_inputs=None):
    """Initialize a basin worker with the shared run configuration."""
    global _EVALUATOR
    if not c.get('_run_data_ready'):
        raise RuntimeError('The coordinator must prepare basin inputs before workers start.')
    _EVALUATOR = Evaluator(c, elevation_inputs)

def _task(gauge, x, metrics=False):
    """Evaluate one basin task in a worker process."""
    return _EVALUATOR.evaluate_with_metrics(gauge, x, predict_only=metrics == 'predict', omit_failure=metrics == 'omit') if metrics else _EVALUATOR.evaluate(gauge, x)

def evaluate_batch(evaluator, executor, gauges, x, limit, stop, metrics=False, progress=None):
    """Evaluate a basin batch and collect its objectives and gradients."""
    if executor is None:
        rows = []
        for gauge, theta in zip(gauges, x):
            if stop.is_set():
                raise Cancelled()
            rows.append(evaluator.evaluate_with_metrics(gauge, theta, predict_only=metrics == 'predict', omit_failure=metrics == 'omit') if metrics else evaluator.evaluate(gauge, theta))
            if progress:
                progress(len(rows), len(gauges))
        return rows
    results = [None] * len(gauges)
    pending = {}
    index = 0
    completed = 0
    while pending or index < len(gauges):
        if stop.is_set():
            for future in pending:
                future.cancel()
            raise Cancelled()
        while index < len(gauges) and len(pending) < limit:
            future = executor.submit(_task, gauges[index], x[index], metrics)
            pending[future] = index
            index += 1
        done, _ = wait(pending, timeout=0.1, return_when=FIRST_COMPLETED)
        for future in done:
            position = pending.pop(future)
            try:
                results[position] = future.result()
            except Exception:
                for other in pending:
                    other.cancel()
                while pending:
                    done_old, _ = wait(pending, timeout=0.1)
                    for other in done_old:
                        pending.pop(other, None)
                    if stop.is_set():
                        raise Cancelled()
                raise
            completed += 1
            if progress:
                progress(completed, len(gauges))
    return results

def _postprocess_basin(evaluator, gauge, x, role, retained, clock):
    """Postprocess basin."""
    c = evaluator.c
    basin_clock = perf_counter()
    data, _, kwargs = evaluator.prepare(gauge)
    theta = physical(c, x)
    from simulation_failure import IntegrationFailure
    series = None
    try:
        if retained:
            result, series = basin_output(c, gauge, data, kwargs, theta, role, started=clock, compute_attribution=c.get('diagnostics', {}).get('attribution', False))
        else:
            _, _, result = evaluator.evaluate_with_metrics(gauge, x, predict_only=True, postprocess=True)
            result.update(model=c['model'], theta=theta.tolist(), basin_role=role, status='Numerical failure' if result.get('failed') else 'Shared network', series_file=None, seconds=perf_counter() - basin_clock, trials=[])
    except IntegrationFailure as exc:
        _, _, result = evaluator.failed_record(gauge, x, str(exc))
        result.update(model=c['model'], theta=theta.tolist(), basin_role=role, status='Numerical failure', series_file=None, seconds=perf_counter() - basin_clock, trials=[])
    return (result, series)

def _postprocess_task(*args):
    """Postprocess task."""
    return _postprocess_basin(_EVALUATOR, *args)

def _postprocess_results(executor, evaluator, tasks, limit, stop):
    """Postprocess results."""
    if executor is None:
        for task in tasks:
            if stop.is_set():
                raise Cancelled()
            yield _postprocess_basin(evaluator, *task)
        return
    pending = set()
    tasks = iter(tasks)
    try:
        while True:
            if stop.is_set():
                raise Cancelled()
            while len(pending) < limit:
                item = next(tasks, None)
                if item is None:
                    break
                pending.add(executor.submit(_postprocess_task, *item))
            if not pending:
                return
            done, _ = wait(pending, timeout=0.1, return_when=FIRST_COMPLETED)
            for future in done:
                if stop.is_set():
                    raise Cancelled()
                pending.remove(future)
                yield future.result()
    finally:
        for future in pending:
            future.cancel()

def learning_rate(c, step):
    """Return the learning rate for the selected iteration."""
    if c['lr_schedule'] == 'constant':
        return c['learning_rate']
    fraction = min(1, step / max(1, c['max_iter'] - 2))
    hold = c['lr_hold'] if c['lr_schedule'] == 'hold_cosine' else 0
    progress = max(0, (fraction - hold) / (1 - hold))
    return c['lr_min'] + 0.5 * (c['learning_rate'] - c['lr_min']) * (1 + np.cos(np.pi * progress))

def run(config, output_dir, events, stop):
    """Execute SAGE training and persist the run results.

    Parameters
    ----------
    config : dict
        Validated model, data, network, optimizer and output settings.
    output_dir : path-like
        Destination directory for this isolated run.
    events : queue-like
        Receives progress, timing, diagnostics and completion events.
    stop : event-like
        Signals a requested stop; the configured stop action controls
        saving.

    Notes
    -----
    The run writes its report, network and requested basin time series to
    output_dir. Failures are recorded in the report and completion events.
    """
    report = dict(application='SAGE', config=config, basins=[], errors=[], history=[], trace=[], status='running')
    out = Path(output_dir)
    created = False
    stop_choice = None
    executor = None
    clock = perf_counter()
    report['timing'] = {}
    phase = 'data'
    phase_started = clock
    selection_published = False
    training_started = False
    best_arrays = None

    def mark(next_phase):
        """Close the current runtime phase and start the next phase."""
        nonlocal phase, phase_started
        now = perf_counter()
        report['timing'][phase] = report['timing'].get(phase, 0.0) + now - phase_started
        phase = next_phase
        phase_started = now
        emit('timing', timing=dict(report['timing']))

    def emit(kind, **kw):
        """Publish a structured run event."""
        event = dict(kind=kind, **kw)
        events.put(event)
        from sage_log import bounded_lines, event as format_event
        message = format_event(event)
        if created and message:
            with (out / 'sage.log').open('a', encoding='utf-8') as handle:
                handle.write('\n'.join(bounded_lines(message)) + '\n')

    def finish_training(status='completed'):
        """Publish the retained iteration once, before postprocessing."""
        nonlocal selection_published
        if selection_published:
            return
        from sage_model_selection import completed, selection_text
        selected = report.get('model_selection')
        done = completed(report)
        if selected is not None:
            selected.update(completed_iterations=done, maximum_iterations=c['max_iter'])
        report['training_complete'] = True
        emit('training_complete', done=done, total=c['max_iter'], stage=13 + band_step, status=status, selection=selected)
        with (out / 'demo_SAGE_iterations.txt').open('a', encoding='utf-8') as stream:
            stream.write('\n' + selection_text(selected, done) + '\n')
        selection_published = True
    try:
        c = worker_configuration(config)
        from sage_run_data import create
        create(c)
        c['_run_data_ready'] = False
        report['config'] = c
        out.mkdir(parents=True, exist_ok=False)
        created = True
        from sage_log import header
        emit('log', message=header(c))
        if c.get('data_action_note'):
            from sage_log import automatic_actions
            emit('log', message=automatic_actions(c['data_action_note']))
        for number, label in [(1, 'CRR backend'), (2, 'read basin attributes'), (3, 'use prepared basins'), (4, 'read model'), (5, 'read numerical settings'), (6, 'build train/eval time split')]:
            emit('log', message=f'({number:02d}) {label} ...')
        from sage_quality import screen, active_configuration
        quality = screen(c, stop=stop, reader_progress=lambda stage, gauge, done, total: emit('data_progress', stage=stage, gauge=gauge, done=done, total=total))
        from cartplus_basin_inputs import prepare_cartplus_bands, uses_elevation_bands
        band_step = int(uses_elevation_bands(c))
        elevation_inputs = {}
        if band_step:
            emit('log', message='(09) read elevation band data ...')
            for row in quality:
                if stop.is_set():
                    raise Cancelled()
                if not row['eligible']:
                    continue
                gauge = row['gauge']
                inputs = {}
                try:
                    prepare_cartplus_bands(c, gauge, inputs)
                except Exception:
                    continue
                elevation_inputs[gauge] = inputs
        c['_run_data_ready'] = True
        mark('initialization')
        report['data_quality'] = quality
        emit('quality', rows=quality)
        emit('log', message=f'({9 + band_step:02d}) check data consistency ...')
        unreadable = [row for row in quality if row.get('read_error')]
        if unreadable:
            raise ValueError('Cannot read requested basin ' + unreadable[0]['gauge'] + ': ' + unreadable[0]['reason'])
        if not c.get('resume_checkpoint'):
            report['requested_config'] = dict(c)
            c = active_configuration(c, quality)
            report['config'] = c
        elif not all((row['eligible'] for row in quality)):
            raise ValueError('Checkpoint basin data no longer pass data-quality screening; start a new run to change the basin population.')
        from sage_quality_summary import summary
        for line in summary(quality):
            emit('log', message='      ' + line)
        count = len(c['gauges'])
        train, held_out = basin_partition(c)
        n = len(train)
        A, scaler = standardize(attribute_matrix(c), train, __import__('sage_attribute_preprocessing').settings(c))
        report.update(training_basins=[c['gauges'][i] for i in train], held_out_basins=[c['gauges'][i] for i in held_out], attribute_scaler=scaler)
        from sage_report_metadata import snapshot as geographic_snapshot
        report.update(geographic_snapshot(c))
        network = Network(A.shape[1], len(c['lower']), c['hidden_widths'], c['activations'], c['network_seed'])
        first = [np.zeros_like(v) for v in network.arrays]
        second = [v.copy() for v in first]
        start_step = 0
        best = np.inf
        best_arrays = None
        pending_proposal = None
        optimizer_step = 0
        rate_scale = 1.0
        if c.get('resume_checkpoint'):
            meta, arrays = load_training(c['resume_checkpoint'], c)
            for value, saved in zip(network.arrays, arrays['current']):
                value[:] = saved
            first = arrays['first']
            second = arrays['second']
            best_arrays = arrays['best']
            best = meta['best']
            start_step = meta['next_iteration']
            pending_proposal = arrays.get('recovery')
            optimizer_step = meta.get('optimizer_step', start_step)
            rate_scale = meta.get('rate_scale', 1.0)
            for key in ('history', 'trace', 'attribute_scaler', 'training_basins', 'held_out_basins'):
                report[key] = meta[key]
            for key, value in meta['config'].items():
                if key.startswith('_observation_'):
                    c[key] = value
            for key in ('rejected_updates', 'recovered_updates'):
                report[key] = meta.get(key, [])
            scaler = meta['attribute_scaler']
            from sage_model_selection import select
            selected = meta.get('model_selection') or select(report, best)
            if selected is not None:
                report['model_selection'] = selected
            for key in ('best_metrics', 'best_parameter_field'):
                if meta.get(key) is not None:
                    report[key] = meta[key]
            raw = attribute_matrix(c)
            from sage_attribute_preprocessing import apply as apply_attributes
            A = apply_attributes(raw, scaler)
        from sage_retention import selection as series_selection
        retained = series_selection(c, [c['gauges'][i] for i in train], [c['gauges'][i] for i in held_out])
        c['diagnostic_selection'] = dict(c.get('diagnostic_selection', {}), selection=retained)
        report['series_retention'] = dict(policy='selected postprocessor basins only', gauges=list(retained))
        evaluator = Evaluator(c, elevation_inputs)
        if use_joint(c):
            emit('log', message='      Preparing fixed SAGE observation scales and training-basin availability counts…')
            initial_x, _ = network.forward(A[train])
            normalize(c, evaluator, [c['gauges'][i] for i in train], initial_x, stop)
        emit('log', message=f'({10 + band_step:02d}) setup parallel pool ...')
        workers = min(c['workers'], len(train)) if c['parallel'] else 1
        if c['parallel']:
            for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
                os.environ[key] = '1'
            executor = ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'), initializer=_initialize, initargs=(c, elevation_inputs))
        emit('pool', workers=workers, max_pending=c['max_pending'], logical=os.cpu_count() or 1)
        emit('log', message=f'({11 + band_step:02d}) prepare discharge statistics ...')
        emit('log', message=f'({12 + band_step:02d}) initialize performance diagnostics ...')
        emit('log', message=f"      Shared-network training: {n} training basins; {count - n} held-out basins.\n      Loss = mean selected basin objective ({c['loss']}).")
        mark('training')
        emit('diagnostic_history', history=report['history'], config=c)
        training_started = True
        emit('iteration_progress', done=start_step, total=c['max_iter'], stage=13 + band_step)
        elapsed_offset = report['history'][-1]['seconds'] if start_step else 0.0
        from sage_iteration_format import Formatter
        iteration_formatter = Formatter()
        with (out / 'demo_SAGE_iterations.txt').open('w', encoding='utf-8', newline='') as transcript:
            for previous in report['history']:
                if 'diagnostics' in previous:
                    transcript.write(previous.get('iteration_text') or iteration_formatter.format(previous, c))
                if 'diagnostics' in previous and 'iteration_text' in previous:
                    iteration_formatter.format(previous, c)
        previous_update = report['history'][-1].get('next_update', {}) if start_step else {}
        from sage_stop import requested, save_network
        for step in range(start_step, c['max_iter']):
            stop_choice = requested(out)
            if stop_choice is not None:
                if not report['history']:
                    raise Cancelled()
                break
            iteration_clock = perf_counter()
            if stop.is_set():
                raise Cancelled()
            emit('iteration_progress', done=step, total=c['max_iter'], current=step + 1, stage=13 + band_step)
            x, cache = network.forward(A[train])
            gauges = [c['gauges'][i] for i in train]
            progress_clock = [0.0]
            rate = learning_rate(c, step) * rate_scale

            def training_progress(done, total, phase='training'):
                """Publish progress for the active training iteration."""
                now = perf_counter()
                if done in (0, 1, total) or now - progress_clock[0] >= 0.25:
                    progress_clock[0] = now
                    emit('training_progress', iteration=step + 1, done=done, total=total, phase=phase, training_total=n, evaluation_total=len(held_out), learning_rate=float(rate))
            training_progress(0, len(gauges))
            x, cache = network.forward(A[train])
            rows = evaluate_batch(evaluator, executor, gauges, x, c['max_pending'], stop, metrics='omit', progress=training_progress)
            valid = [row for row in rows if not row[2].get('failed')]
            if not valid:
                document_failures(out, report, step + 1, rows, [])
                stop_failed_training(out, report, step + 1, len(rows))
            value = float(np.mean([r[0] for r in valid]))
            if not np.isfinite(value):
                raise ValueError('Nonfinite mean training objective.')
            basin_gradients = np.zeros_like(x)
            for index, row in enumerate(rows):
                if not row[2].get('failed'):
                    basin_gradients[index] = row[1] / len(valid)
            grad = network.backward(cache, basin_gradients)
            norm = float(np.sqrt(sum((np.sum(v * v) for v in grad))))
            if not np.isfinite(norm):
                raise ValueError('Nonfinite network gradient.')
            rate = learning_rate(c, step) * rate_scale
            item = dict(iteration=step + 1, objective=value, gradient_norm=norm, learning_rate=float(rate), seconds=elapsed_offset + perf_counter() - clock)
            report['history'].append(item)
            field, information_cache = network.forward(A)
            information_layers = [np.tanh(z) if network.activations[i] == 'tanh' else np.maximum(z, 0) for i, (_, z) in enumerate(information_cache[0][:-1])]
            emit('information_snapshot', snapshot=dict(iteration=step + 1, gauges=c['gauges'], attributes=A, layers=information_layers, normalized=field, activations=list(network.activations)))
            report['parameter_field'] = dict(gauges=c['gauges'], normalized=field.tolist(), source=f'current network at iteration {step + 1}')
            emit('parameter_field', field=report['parameter_field'])
            if len(held_out):
                training_progress(0, len(held_out), 'evaluation')
            evaluation_rows = evaluate_batch(evaluator, executor, [c['gauges'][i] for i in held_out], field[held_out], c['max_pending'], stop, metrics='predict', progress=lambda done, total: training_progress(done, total, 'evaluation'))
            training_progress(count, count, 'basins')
            metrics = [dict(row[2], basin_role=role) for batch, role in [(rows, 'training'), (evaluation_rows, 'held out')] for row in batch]
            failures = [r for r in metrics if r.get('failed')]
            if failures:
                item['integration_failures'] = failures
                document_failures(out, report, step + 1, rows, evaluation_rows)
                emit('log', message=f'      Integration failures: {len(failures)} basins; details recorded in Run Notes and integration_failures.jsonl.')
            live_series = []
            for record in metrics:
                series = record.pop('_series', None)
                if series is not None:
                    from sage_retention import restrict_periods
                    live_series.append(dict(record, iteration=step + 1, series_data=restrict_periods(series, retained[record['gauge']])))
            from sage_print_view import summarize
            item['diagnostics'] = summarize(metrics, dict(c, omit_failed_basins=True))
            attribution_rows = [r.pop('_attribution') for r in metrics if '_attribution' in r]
            if attribution_rows:
                reduced = np.nanmedian(np.abs(np.asarray(attribution_rows)), axis=0)
                item['attribution'] = dict(total=reduced[0].tolist(), net=reduced[1].tolist(), statistic='median(abs(.))', basins='training')
                absolute = np.abs(np.asarray(attribution_rows))
                item['attribution']['summaries'] = {}
                for label, reducer in [('median(abs(.))', np.nanmedian), ('mean(abs(.))', np.nanmean), ('max(abs(.))', np.nanmax)]:
                    summary = reduced if label == 'median(abs(.))' else reducer(absolute, axis=0)
                    item['attribution']['summaries'][label] = dict(total=summary[0].tolist(), net=summary[1].tolist())
            has_eval = c.get('evaluation_enabled', True) if c['split'] == 'manual' else True
            item['scenarios'] = ['tt'] + (['te'] if has_eval else []) + (['et'] + (['ee'] if has_eval else []) if len(held_out) else [])
            item['network_parameters'] = sum((v.size for v in network.arrays))
            item['iteration_seconds'] = perf_counter() - iteration_clock
            item['iteration_text'] = iteration_formatter.format(item, c)
            with (out / 'demo_SAGE_iterations.txt').open('a', encoding='utf-8', newline='') as transcript:
                transcript.write(item['iteration_text'])
            report['live_metrics'] = metrics
            report['iteration_parameter_field'] = report['parameter_field']
            emit('metrics', rows=metrics, iteration=step + 1)
            quantiles = np.quantile(x, [0.05, 0.15, 0.25, 0.5, 0.75, 0.85, 0.95], axis=0, method='hazen').T.tolist()
            report['trace'].append(dict(iteration=step + 1, quantiles=quantiles))
            improved = value < best
            if improved:
                best = value
                best_arrays = [v.copy() for v in network.arrays]
                from sage_model_selection import record
                selected = record(report, item)
                selected.update(training_basins_used=len(valid), training_basins_total=len(rows))
                report['best_metrics'] = deepcopy(metrics)
                report['best_parameter_field'] = deepcopy(report['parameter_field'])
                report['best_parameter_field']['source'] = f'best training network, iteration {step + 1}'
            item['iteration_complete'] = True
            clipped = bool(c['clip_norm'] > 0 and norm > c['clip_norm'])
            if clipped:
                grad = [v * c['clip_norm'] / (norm + 1e-08) for v in grad]
            item['optimizer'] = dict(previous_update, gradNormPre=norm, gradNormPost=float(np.sqrt(sum((np.sum(v * v) for v in grad)))), gradientClipped=clipped)
            from sage_solver_recovery import capture
            optimizer_step += 1
            pending_proposal = capture(network, first, second, grad, optimizer_step, rate)
            previous_update = descent_update(network.arrays, grad, first, second, optimizer_step, c, rate)
            item['next_update'] = previous_update
            save_training(out / 'training_checkpoint.npz', c, network, first, second, best_arrays, step + 1, best, report, pending=pending_proposal, optimizer_step=optimizer_step, rate_scale=rate_scale)
            if improved:
                try:
                    save_network(out, best_arrays, c['attributes'], scaler, selection=report.get('model_selection'))
                except PermissionError as exc:
                    emit('log', message=f'Could not publish network.npz during training. The resume checkpoint was saved. {exc}')
            atomic_json(out / 'run.json', report)
            emit('live_series', rows=live_series, iteration=step + 1)
            emit('sage_iteration', **item)
            emit('iteration_progress', done=step + 1, total=c['max_iter'], stage=13 + band_step)
            stop_choice = requested(out)
            if stop_choice is not None:
                break
        finish_training('completed' if stop_choice is None else 'stopped after')
        for array, saved in zip(network.arrays, best_arrays):
            array[:] = saved
        save_network(out, best_arrays, c['attributes'], scaler, selection=report.get('model_selection'))
        x, _ = network.forward(A)
        report['parameter_field'] = dict(gauges=c['gauges'], normalized=x.tolist(), source=f"best training network, iteration {report['model_selection']['iteration']}")
        emit('parameter_field', field=report['parameter_field'])
        if report.get('best_metrics'):
            report['live_metrics'] = report['best_metrics']
            emit('metrics', rows=report['best_metrics'], iteration=report['model_selection']['iteration'])
            emit('live_series', rows=[], iteration=report['model_selection']['iteration'])
        if stop_choice is not None:
            report['stop_requested'] = True
            report['postprocess_after_stop'] = stop_choice['postprocess']
            if not stop_choice['postprocess']:
                report['status'] = 'stopped'
                emit('log', message='Stopped after the completed iteration; retained network saved; postprocessing skipped.')
                return
        mark('postprocessing')
        stage = 14 + band_step
        total = len(c['gauges'])
        emit('postprocessing_started', stage=stage, total=total, selection=report.get('model_selection'))
        from sage_retention import restrict_periods
        tasks = ((gauge, x[i], 'training' if i in train else 'held out', gauge in retained, clock) for i, gauge in enumerate(c['gauges']))
        checkpoint_delayed = False
        with closing(_postprocess_results(executor, evaluator, tasks, c['max_pending'], stop)) as results:
            for result, series in results:
                gauge = result['gauge']
                if series is not None:
                    result['trials'] = [dict(trajectory=[row['objective'] for row in report['history']])]
                    save_npz(out / result['series_file'], **restrict_periods(series, retained[gauge]))
                if result.get('failed'):
                    emit('log', message=f"      {gauge}: final simulation failed: {result['failure_reason']}; no series saved.")
                report['basins'].append(result)
                try:
                    atomic_json(out / 'run.json', report, attempts=5)
                except OSError as exc:
                    if not checkpoint_locked(exc):
                        raise
                    if not checkpoint_delayed:
                        emit('log', message='      Run checkpoint is temporarily locked. Postprocessing continues and will retry the save.')
                    checkpoint_delayed = True
                else:
                    checkpoint_delayed = False
                emit('result', result=result, directory=str(out))
                emit('postprocessing_progress', stage=stage, done=len(report['basins']), total=total)
        order = {gauge: i for i, gauge in enumerate(c['gauges'])}
        report['basins'].sort(key=lambda row: order[row['gauge']])
        report['postprocessing_complete'] = True
        report['status'] = 'stopped' if stop_choice is not None else 'complete'
    except __import__('sage_solver_recovery').RecoveryExhausted as exc:
        report['status'] = 'stopped'
        report['stop_reason'] = 'numerical recovery exhausted'
        report['recovery_exhausted'] = True
        report['errors'].append(dict(message=str(exc)))
        if best_arrays is not None:
            from sage_stop import save_network
            save_network(out, best_arrays, c['attributes'], scaler, selection=report.get('model_selection'))
            save_training(out / 'training_checkpoint.npz', c, network, first, second, best_arrays, step, best, report, pending=None, optimizer_step=optimizer_step, rate_scale=rate_scale)
        emit('log', message=str(exc))
    except Cancelled:
        report['status'] = 'stopped'
        if best_arrays is not None:
            from sage_stop import save_network
            save_network(out, best_arrays, c['attributes'], scaler, selection=report.get('model_selection'))
        else:
            emit('log', message='Stopped before a completed iteration; no retained network is available for postprocessing.')
    except Exception as exc:
        report['status'] = 'failed'
        report['errors'].append(dict(message=str(exc)))
        emit('error', message=str(exc) + '\n' + traceback.format_exc(), fatal=True)
    finally:
        if created and training_started and (not selection_published):
            if 'best_arrays' in locals():
                finish_training('stopped after' if report['status'] == 'stopped' else 'ended after')
                if report.get('best_metrics'):
                    report['live_metrics'] = report['best_metrics']
                    emit('metrics', rows=report['best_metrics'], iteration=report['model_selection']['iteration'])
                    emit('live_series', rows=[], iteration=report['model_selection']['iteration'])
                if report.get('best_parameter_field'):
                    report['parameter_field'] = report['best_parameter_field']
                    emit('parameter_field', field=report['parameter_field'])
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        if 'evaluator' in locals():
            evaluator.cache.clear()
        if 'c' in locals():
            from sage_run_data import cleanup
            cleanup(c)
        mark('finished')
        report['timing']['total'] = perf_counter() - clock
        emit('timing', timing=dict(report['timing']))
        if created:
            for key in ('config', 'requested_config'):
                if key in report:
                    report[key] = {name: value for name, value in report[key].items() if not name.startswith('_run_data_')}
            try:
                atomic_json(out / 'run.json', report)
            except Exception as exc:
                report['status'] = 'checkpoint failed'
                emit('error', message='Final checkpoint could not be replaced: ' + str(exc), fatal=True)
        emit('finished', status=report['status'], directory=str(out), count=len(report['basins']), postprocessing_complete=bool(report.get('postprocessing_complete')))
