# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_checkpoint.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_checkpoint."""
import json
from pathlib import Path
from time import sleep
import numpy as np
from sage_network import Network

def save_npz(path, *, attempts=100, **arrays):
    """Save an array archive for a training checkpoint."""
    path = Path(path)
    from tempfile import NamedTemporaryFile
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        for attempt in range(attempts):
            try:
                if temporary is None or not temporary.exists():
                    with NamedTemporaryFile(dir=path.parent, prefix=path.stem + '.', suffix='.tmp.npz', delete=False) as stream:
                        temporary = Path(stream.name)
                        try:
                            np.savez_compressed(stream, **arrays)
                        except (PermissionError, FileNotFoundError):
                            incomplete = temporary
                            temporary = None
                            stream.close()
                            incomplete.unlink(missing_ok=True)
                            raise
                temporary.replace(path)
                return
            except (PermissionError, FileNotFoundError):
                if attempt + 1 == attempts:
                    raise
                sleep(0.15)
    finally:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except PermissionError:
                pass

def save_training(path, c, network, first, second, best_arrays, step, best, report, pending=None, optimizer_step=None, rate_scale=1.0):
    """Save the network and optimizer training checkpoint."""
    payload = {f'current_{i}': v for i, v in enumerate(network.arrays)}
    for prefix, arrays in (('first', first), ('second', second), ('best', best_arrays)):
        payload.update({f'{prefix}_{i}': v for i, v in enumerate(arrays)})
    metadata = dict(version=2, config=c, next_iteration=step, best=best, history=report['history'], trace=report.get('trace', []), attribute_scaler=report['attribute_scaler'], training_basins=report['training_basins'], held_out_basins=report['held_out_basins'], model_selection=report.get('model_selection'), best_metrics=report.get('best_metrics'), best_parameter_field=report.get('best_parameter_field'))
    metadata.update(optimizer_step=step if optimizer_step is None else optimizer_step, rate_scale=rate_scale, rejected_updates=report.get('rejected_updates', []), recovered_updates=report.get('recovered_updates', []))
    if pending is not None:
        metadata['recovery'] = {k: pending[k] for k in ('step', 'rate', 'factor')}
        for name in ('base', 'first', 'second', 'gradient'):
            payload.update({f'recovery_{name}_{i}': v for i, v in enumerate(pending[name])})
    from engine import json_scalar
    payload['metadata'] = np.array(json.dumps(metadata, allow_nan=False, default=json_scalar))
    save_npz(path, **payload)

def load_training(path, c):
    """Load and validate a training-resume checkpoint."""
    path = __import__('sage_paths').relocate(str(path))
    with np.load(path, allow_pickle=False) as data:
        meta = __import__('sage_paths').relocate(json.loads(str(data['metadata'])))
        if meta.get('version') not in (1, 2):
            raise ValueError('Unsupported training checkpoint version.')
        original = meta['config']
        from sage_parameter_scaling import modes
        if modes(original) != modes(c):
            raise ValueError('Parameter scaling differs from checkpoint; start a fresh run after changing Scale.')
        from hydrologic_data_actions import effective_actions
        if effective_actions(original.get('data_actions', []), c['gauges']) != effective_actions(c.get('data_actions', []), c['gauges']):
            raise ValueError('Regional data corrections differ from the training checkpoint; start a fresh training run with corrected data.')
        if c.get('split') in ('random', 'random_block', 'random_kfold') and original.get('temporal_rng', 'numpy') != c.get('temporal_rng', 'matlab'):
            raise ValueError('Resume permutation policy differs from checkpoint. Older Python runs require NumPy legacy; start a new MATLAB-matched run to change its split.')
        from sage_zscore import transform
        if any((transform(original, key) != transform(c, key) for key in ('Q', 'SWE', 'SM'))):
            raise ValueError('Resume response objective differs from checkpoint.')
        mutable = {'max_iter', 'parallel', 'workers', 'max_pending', 'output_root', 'resume_checkpoint', 'notes', 'diagnostics', 'diagnostic_selection', 'optimizer_options', 'postprocess_selection', '_run_data_directory', '_run_data_ready'}
        for key, value in original.items():
            if key == 'model_config' and c.get('model') in (18, 19, 20):
                continue
            if key == 'basin_split_file' and value and c.get(key) and (value != c[key]):
                from sage_basin_presets import resolve_file
                if resolve_file(original).read_bytes() == resolve_file(c).read_bytes():
                    continue
            if key not in mutable and key != 'parameter_scales' and (key != 'observation_transforms') and (not key.startswith('_observation_')) and (c.get(key) != value):
                raise ValueError(f'Resume setting differs from checkpoint: {key}')
        if c['max_iter'] <= meta['next_iteration']:
            raise ValueError('Increase maximum iterations to extend this checkpoint.')
        arrays = {prefix: [data[f'{prefix}_{i}'].copy() for i in range(len([k for k in data.files if k.startswith(prefix + '_')]))] for prefix in ('current', 'first', 'second', 'best')}
        if meta.get('recovery'):
            arrays['recovery'] = dict(meta['recovery'])
            for name in ('base', 'first', 'second', 'gradient'):
                arrays['recovery'][name] = [data[f'recovery_{name}_{i}'].copy() for i in range(len(arrays['current']))]
    return (meta, arrays)

def load_network(directory):
    """Load the saved network parameter arrays."""
    directory = Path(__import__('sage_paths').relocate(str(directory)))
    report = __import__('sage_paths').load_result(directory / 'run.json')
    c = report['config']
    net = Network(len(c['attributes']), len(c['lower']), c['hidden_widths'], c['activations'], c['network_seed'])
    c.setdefault('temporal_rng', 'numpy')
    with np.load(directory / 'network.npz', allow_pickle=False) as data:
        for i, array in enumerate(net.arrays):
            array[:] = data[f'array_{i}']
    return (net, report)
