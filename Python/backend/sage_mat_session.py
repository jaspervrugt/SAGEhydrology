# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_mat_session.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_mat_session."""
from copy import deepcopy
from io import BytesIO
from pathlib import Path
import sys
import numpy as np
from scipy.io import loadmat
from scipy.io.matlab import MatlabOpaque, varmats_from_mat

def _string_storage(value):
    """
    MATLAB string's version-1 packed UTF-16 storage (MCOS data
    property).
    """
    a = np.asarray(value)
    if a.dtype != np.dtype('uint64'):
        raise ValueError('Unsupported MATLAB string storage.')
    a = a.ravel()
    ndim = int(a[1]) if len(a) > 1 else 0
    if len(a) < 4 or a[0] != 1 or (not 2 <= ndim <= 8):
        raise ValueError('Unsupported MATLAB string serialization.')
    count = int(np.prod(a[2:2 + ndim]))
    lengths = a[2 + ndim:2 + ndim + count]
    if count > 100000 or len(lengths) != count:
        raise ValueError('Invalid MATLAB string dimensions.')
    raw = a[2 + ndim + count:].astype('<u8').tobytes()
    strings = []
    start = 0
    for length in lengths:
        if length == np.iinfo(np.uint64).max:
            strings.append('')
            continue
        end = start + 2 * int(length)
        if end > len(raw):
            raise ValueError('Truncated MATLAB string data.')
        strings.append(raw[start:end].decode('utf-16-le'))
        start = end
    return strings[0] if len(strings) == 1 else strings

def read(path):
    """Read."""
    path = Path(path)
    with path.open('rb') as handle:
        header = handle.read(128)
    if header.startswith(b'MATLAB 7.3'):
        try:
            import h5py
        except ImportError:
            vendor = Path(__file__).resolve().parents[2] / 'local_dependencies/Python/vendor'
            sys.path.insert(0, str(vendor))
            try:
                import h5py
            except ImportError as exc:
                raise ValueError('v7.3 sessions require h5py. Install it with: python -m pip install h5py') from exc
        with h5py.File(path, 'r') as f:
            if 'session' not in f:
                raise ValueError('MAT-file does not contain the variable "session".')
            wrapper = f.get('#subsystem#/MCOS')
            storage = []
            if wrapper is not None:
                storage = [f[r] for r in wrapper[()].ravel() if r]

            def decode(node, depth=0):
                """Decode."""
                if depth > 40:
                    raise ValueError('MAT session nesting is too deep.')
                kind = node.attrs.get('MATLAB_class', b'')
                kind = kind.decode() if isinstance(kind, bytes) else str(kind)
                if kind in ('datetime', 'duration', 'function_handle', 'opaque'):
                    return None
                if isinstance(node, h5py.Group):
                    return {k: decode(v, depth + 1) for k, v in node.items() if k != 'createdAt'}
                if node.attrs.get('MATLAB_empty', False):
                    return []
                a = node[()]
                if kind == 'string':
                    data = np.asarray(a).ravel()
                    index = int(data[-2]) + 1
                    if len(data) != 6 or not 0 <= index < len(storage):
                        raise ValueError('Unsupported MATLAB string object reference.')
                    return _string_storage(storage[index][()])
                if h5py.check_dtype(ref=node.dtype) is not None:
                    items = [decode(f[r], depth + 1) if r else None for r in np.asarray(a).ravel(order='C')]
                    return items[0] if len(items) == 1 else items
                if kind == 'char':
                    a = np.asarray(a).T
                    if a.ndim < 2 or a.shape[0] == 1:
                        return np.asarray(a, dtype='<u2').tobytes().decode('utf-16-le').rstrip('\x00')
                    return [np.asarray(row, dtype='<u2').tobytes().decode('utf-16-le').rstrip('\x00 ') for row in a]
                a = np.asarray(a).T
                if kind == 'logical':
                    a = a.astype(bool)
                if a.size == 1:
                    return a.item()
                return a.squeeze().tolist()
            return decode(f['session'])
    document = loadmat(path, simplify_cells=True)
    if 'session' not in document:
        raise ValueError('MAT-file does not contain the variable "session".')
    storage = []
    if '__function_workspace__' in document:
        b = document['__function_workspace__'].tobytes()
        for _, stream in varmats_from_mat(BytesIO(b' ' * 124 + b[:4] + b[8:])):
            container = loadmat(stream, simplify_cells=True).get('__function_workspace__')
            if isinstance(container, np.ndarray) and container.dtype == object and container.size:
                obj = container.ravel()[0]
                if hasattr(obj, 'MCOS'):
                    storage = list(obj.MCOS['_ObjectMetadata'][0].ravel())
                    break

    def clean(value):
        """Clean."""
        if isinstance(value, MatlabOpaque):
            name = value['_Class'][0]
            data = value['_ObjectMetadata'][0]
            if name == 'string':
                index = int(np.asarray(data).ravel()[-2]) + 1
                if 0 <= index < len(storage):
                    return _string_storage(storage[index])
                raise ValueError('MATLAB string data are missing from this session.')
            return None
        if isinstance(value, dict):
            return {k: clean(v) for k, v in value.items() if k != 'createdAt'}
        if isinstance(value, np.ndarray):
            return clean(value.item()) if value.size == 1 else [clean(v) for v in value.tolist()]
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, (tuple, list)):
            return [clean(v) for v in value]
        return value
    return clean(document['session'])

def configuration(path, base):
    """Configuration."""
    from sage_matlab_project import translate
    from regional_maps import identifier
    session = read(path)
    if not isinstance(session, dict) or not isinstance(session.get('cfg'), dict):
        raise ValueError('Session does not contain a valid cfg structure.')
    cfg = deepcopy(session['cfg'])
    warnings = []
    prepared = session.get('preparedBasins') or {}
    if not isinstance(prepared, dict):
        prepared = {}
    gauges = np.atleast_1d(prepared.get('id_gauge', []))
    fingerprint = session.get('dataFingerprint') or {}
    compatible = bool(session.get('runtimePrepared')) and len(gauges) > 0 and isinstance(fingerprint, dict) and bool(fingerprint)
    for saved, actual in [('region', cfg.get('region')), ('dt', cfg.get('prd', {}).get('dt')), ('K', cfg.get('bas', {}).get('K')), ('Kt', cfg.get('bas', {}).get('K_t')), ('Ke', cfg.get('bas', {}).get('K_e')), ('sample', cfg.get('bas', {}).get('sample')), ('fileUniverse', cfg.get('file_univ'))]:
        if saved in fingerprint and fingerprint[saved] != actual:
            compatible = False
            warnings.append('Prepared basin fingerprint differs: ' + saved)
    if compatible and fingerprint.get('universeExists'):
        from regional_data import region_registry
        from datetime import datetime
        source = Path(region_registry(base['sage_root'])[cfg['region']]['source']).parent
        candidates = [Path(str(fingerprint.get('universePath', ''))), source / str(cfg.get('file_univ', ''))]
        universe = next((p for p in candidates if p.is_file()), None)
        if universe is None:
            compatible = False
            warnings.append('The saved basin-universe file is unavailable.')
        else:
            stat = universe.stat()
            stamp = datetime.fromtimestamp(stat.st_mtime)
            datenum = stamp.toordinal() + 366 + (stamp.hour * 3600 + stamp.minute * 60 + stamp.second + stamp.microsecond / 1000000.0) / 86400
            if np.isfinite(fingerprint.get('universeBytes', np.nan)) and stat.st_size != fingerprint['universeBytes']:
                compatible = False
                warnings.append('The saved basin-universe file size changed.')
            if np.isfinite(fingerprint.get('universeDateNum', np.nan)) and abs(datenum - fingerprint['universeDateNum']) > 1.1 / 86400:
                compatible = False
                warnings.append('The saved basin-universe file modification time changed.')
    if compatible:
        cfg['_import_gauges'] = [identifier(g, cfg['region']) for g in gauges]
        if len(set(cfg['_import_gauges'])) != len(gauges):
            raise ValueError('Prepared session contains duplicate basin IDs.')
        kt = int(prepared.get('K_t', cfg['bas']['K_t']))
        ke = int(prepared.get('K_e', len(gauges) - kt))
        if kt < 1 or kt + ke != len(gauges):
            raise ValueError('Prepared basin counts do not match the saved basin IDs.')
        cfg['bas'].update(K=len(gauges), K_t=kt, K_e=ke)
    elif session.get('runtimePrepared'):
        warnings.append('Saved basin preparation is stale or has no fingerprint; the basin division must be prepared again.')
    c = translate(cfg, base, Path(path))
    c['basin_split_ids'] = []
    if compatible:
        c['basin_sampling'] = 'file'
        c['basin_split_file'] = ''
        c['basin_split_ids'] = cfg['_import_gauges']
        c['basin_train_fraction'] = kt / len(gauges)
    key = f"m{c['model']:02d}_dt{c['resolution']:02d}"
    overrides = session.get('parOverride') or cfg.get('par_override') or {}
    if not isinstance(overrides, dict):
        raise ValueError('Invalid parameter override structure in the session.')
    override = overrides.get(key, {})
    if override:
        lo = np.atleast_1d(override.get('min', [])).astype(float)
        hi = np.atleast_1d(override.get('max', [])).astype(float)
        if len(lo) != len(c['lower']) or hi.shape != lo.shape or (not np.isfinite(lo).all()) or (not np.isfinite(hi).all()) or np.any(lo >= hi):
            raise ValueError('Saved parameter bounds do not match the selected model.')
        c['lower'] = lo.tolist()
        c['upper'] = hi.tolist()
    note = session.get('runNote', cfg.get('runNote', ''))
    c['notes'] = '\n'.join((str(v) for v in note)) if isinstance(note, list) else str(note or '')
    c['matlab_import'] = dict(source=str(Path(path).resolve()), schema_version=session.get('schemaVersion', 1), warnings=warnings, prepared_restored=compatible, selected_left_tab=session.get('selectedLeftTab', 'Notes'), selected_right_tab=session.get('selectedRightTab', 'Summary'))
    if session.get('schemaVersion', 1) > 3:
        warnings.append('This session uses a newer MATLAB schema; recognized settings were imported.')
    c['output_root'] = base['output_root']
    c['resume_checkpoint'] = ''
    plot = cfg.get('misc', {}).get('plot', {})
    selection = plot.get('gaugescen', {})
    if isinstance(selection, dict) and selection:
        c['postprocess_selection'] = {}
        for scenario, ids in selection.items():
            if scenario not in ('tt', 'te', 'et', 'ee'):
                continue
            period = 'train' if scenario[-1] == 't' else 'evaluation'
            for gauge in np.atleast_1d(ids):
                target = c['postprocess_selection'].setdefault(identifier(gauge, c['region']), [])
                if period not in target:
                    target.append(period)
    return c
