# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_prepared_data.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_prepared_data."""
from pathlib import Path
import hashlib, json, os, tempfile, uuid
import numpy as np
from timeseries_storage import store_timeseries
CACHE_VERSION = 10
_LEGACY_READERS = {7: (('fbfeff5ef9a0f26c2c2b96d626039f5f7b798049bf50451a3c61e53e87fce58d', False), ('1564aa1f21f48b9863ba44d5a93a25fc296890c43138fc638032feb267cfcaca', True))}

def supported(c):
    """Report whether the canonical reader supports prepared-array reuse."""
    return True

def _legacy_supported(c):
    """Legacy supported."""
    if c.get('region') in ('CAMELS_CA', 'CAMELS_SPAT'):
        return c.get('resolution', 1) in (1, 24)
    if c.get('region') == 'CAMELS_DEH' and c.get('resolution') == 24:
        return True
    if c.get('resolution', 1) != 1:
        return False
    if c.get('region') in ('CAMELS_US', 'CAMELS_AU', 'CAMELS_BR'):
        return True
    from regional_special import SPECIAL
    if c.get('region') in SPECIAL:
        return False
    from regional_data import table_profile
    try:
        profile = table_profile(c['sage_root'], c['region'], **{k: c.get(k, 1) for k in ('pet', 'temperature', 'precipitation')})
    except (ValueError, KeyError):
        return False
    return profile.get('layout') == 'one_file_per_basin' and profile.get('format') == 'csv' and (not profile.get('flow_pattern'))

def location(c, gauge, *, reader=None, legacy=False):
    """Locate one model-independent basin/profile cache, as in MATLAB."""
    if reader is not None or legacy:
        return _legacy_location(c, gauge, reader=reader, legacy=legacy)
    import re
    fields = ('region', 'data_root', 'resolution', 'product', 'pet', 'precipitation', 'temperature', 'data_actions')
    profile = {k: c.get(k) for k in fields}
    from hydrologic_data_actions import effective_actions, load
    actions = c.get('data_actions')
    if actions is None:
        actions = load(c.get('region', ''), c.get('resolution', 1), c.get('sage_root'))
    profile['data_actions'] = effective_actions(actions, [gauge])
    key = hashlib.sha256(json.dumps(profile, sort_keys=True, default=str).encode()).hexdigest()[:24]
    region = re.sub('[^A-Za-z0-9_-]', '_', str(c.get('region', 'data')))
    basin = re.sub('[^A-Za-z0-9_-]', '_', str(gauge))
    return Path(os.environ.get('LOCALAPPDATA', tempfile.gettempdir())) / 'SAGEhydrology/Python/canonical_inputs' / f'v{CACHE_VERSION}' / region / f'{basin}_{key}'

def _request(c):
    """Request."""
    from input_period import PERIOD_KEYS
    return {k: c.get(k) for k in sorted(PERIOD_KEYS)}

def _manifest(c, sources):
    """Manifest."""
    from reader_identity import metadata_digest
    return dict(cacheVersion=CACHE_VERSION, request=_request(c), schema=metadata_digest(c['sage_root'], c.get('region')), sources=sources)

def _legacy_location(c, gauge, *, reader=None, legacy=False):
    """Address an explicitly reviewed old entry during one-way migration."""
    from input_period import PERIOD_KEYS
    keys = PERIOD_KEYS | {'region', 'data_root', 'sage_root', 'product', 'pet', 'precipitation', 'temperature', 'data_actions'}
    from reader_identity import reader_digest, metadata_digest
    settings = {k: c.get(k) for k in sorted(keys) if k != 'sage_root'}
    payload = dict(version=5, reader=reader or reader_digest(), metadata=metadata_digest(c['sage_root'], None if legacy else c.get('region')), gauge=str(gauge), settings=settings)
    key = hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()
    folder = Path(os.environ.get('LOCALAPPDATA', tempfile.gettempdir())) / 'SAGEhydrology/Python/prepared_data'
    return folder / key

def stamp(path):
    """Return the path, modification time and size used for invalidation."""
    p = Path(path)
    s = p.stat()
    return [str(p.resolve()), s.st_mtime_ns, s.st_size]

def dependencies(c, data):
    """Validate exactly the generic reader's declared data dependencies."""
    paths = {str(Path(p).resolve()) for p in data.get('_source_files', [])}
    paths.update((str(data[k]) for k in ('forcing_file', 'streamflow_file') if data.get(k)))
    paths.update((str(Path(p).parent) for p in tuple(paths)))
    paths.add(str(Path(c['data_root']).resolve()))
    return [stamp(p) for p in sorted(paths)]

def get(c, gauge, *, status=None):
    """Load canonical arrays only when every recorded dependency matches.

    Return None on a missing, stale or unreadable cache entry.
    """
    status = {} if status is None else status
    current = location(c, gauge)
    result = _get(c, current, status)
    if result is not None:
        status['reason'] = 'reused'
        return result
    for base in _legacy_candidates(c, gauge):
        result = _get(c, base, status, legacy=True)
        if result is not None:
            _promote(c, current, result, status['_sources'])
            status['reason'] = 'reused'
            return result
    return None

def _legacy_candidates(c, gauge):
    """Legacy candidates."""
    return [_legacy_location(c, gauge, reader=reader, legacy=legacy) for reader, legacy in _LEGACY_READERS.get(CACHE_VERSION, ())]

def _promote(c, current, data, sources):
    """Write existing arrays into MATLAB's single atomic payload pattern."""
    try:
        _save(current, data, _manifest(c, sources))
    except OSError:
        pass

def _get(c, base, status, *, legacy=False):
    """Get."""
    from zipfile import BadZipFile
    try:
        with np.load(base.with_suffix('.npz'), allow_pickle=False) as saved:
            manifest = json.loads(base.with_suffix('.json').read_text()) if legacy else json.loads(str(saved['_cache_manifest'].item()))
            if not legacy:
                expected = _manifest(c, manifest['sources'])
                checks = (('cacheVersion', 'cache_version_changed'), ('request', 'request_changed'), ('schema', 'input_schema_changed'))
                for key, reason in checks:
                    if manifest.get(key) != expected[key]:
                        status['reason'] = reason
                        return None
            if any((stamp(row[0]) != row for row in manifest['sources'])):
                status['reason'] = 'source_changed'
                return None
            data = {k: saved[k] for k in saved.files if k != '_cache_manifest'}
            status['_sources'] = manifest['sources']
        for key in ('gauge', 'forcing_file', 'streamflow_file'):
            if key in data:
                data[key] = str(data[key].item())
        if c.get('region') == 'CAMELS_US' and c.get('resolution', 1) == 1 and ('area' not in data):
            status['reason'] = 'missing_replay_metadata'
            return None
        return store_timeseries(data)
    except FileNotFoundError:
        status.setdefault('reason', 'not_cached_for_inputs')
        return None
    except (OSError, ValueError, KeyError, EOFError, BadZipFile):
        status['reason'] = 'unreadable_cache'
        return None

def put(c, gauge, data):
    """Atomically cache numeric input arrays and their source manifest.

    Exclude model results and object arrays. Cache write failures must never
    prevent scientific input loading.
    """
    data = dict(data)
    for key in ('forcing', 'observed', '_source_SWE', '_source_SM'):
        if key in data:
            data[key] = np.asarray(data[key], dtype=np.float32).astype(float)
    if not supported(c):
        return 'unsupported'
    if len(data.get('_source_files', [])) == 0 and (not _legacy_supported(c)):
        return 'missing_source_manifest'
    base = location(c, gauge)
    try:
        sources = dependencies(c, data)
        keys = ('gauge', 'area', 'dates', 'forcing', 'matlab_SW', 'matlab_DayLength', 'matlab_observed', 'observed', 'train', 'evaluation', 'forcing_file', 'streamflow_file')
        payload = {k: data[k] for k in keys if k in data}
        payload.update({k: v for k, v in data.items() if k.startswith('_source_')})
        if any((np.asarray(v).dtype.hasobject for v in payload.values())):
            return 'unsupported_array_type'
        _save(base, payload, _manifest(c, sources))
        return 'written'
    except OSError:
        return 'write_failed'

def _save(base, data, manifest):
    """Atomically replace one file containing both arrays and signatures."""
    temporary = base.parent / uuid.uuid4().hex
    try:
        base.parent.mkdir(parents=True, exist_ok=True)
        text = json.dumps(manifest, sort_keys=True, default=str)
        payload = dict(data, _cache_manifest=np.array(text))
        with temporary.with_suffix('.npz').open('wb') as handle:
            np.savez(handle, **payload)
        os.replace(temporary.with_suffix('.npz'), base.with_suffix('.npz'))
        try:
            temporary.with_suffix('.json').write_text(text)
            os.replace(temporary.with_suffix('.json'), base.with_suffix('.json'))
        except OSError:
            pass
    finally:
        for suffix in ('.npz', '.json'):
            try:
                temporary.with_suffix(suffix).unlink(missing_ok=True)
            except OSError:
                pass
