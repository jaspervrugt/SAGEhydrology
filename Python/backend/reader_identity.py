# ==========================================================================
# SAGEhydrology | Python source
# File: backend/reader_identity.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from reader_identity."""
from functools import lru_cache
import hashlib
import importlib.util
import marshal
from pathlib import Path
import types
MODULES = ('timeseries_storage', 'hydro_reader', 'hydro_schema', 'hydro_io', 'hydro_aux', 'hydro_derivations', 'regional_data', 'reader_arrays', 'regional_special', 'regional_subdaily', 'subdaily_profiles', 'camels_us', 'local_data', 'netcdf_time', 'netcdf_io', 'wide_data_cache', 'runtime_metadata', 'sage_quality', 'sage_science', 'hydrologic_data_actions', 'regional_observation_data', 'camels_spat_io', 'sage_period')

def normalized(code):
    """Normalized."""
    import dis
    docs = {0} if code.co_flags & 2 else set()
    locations = set()
    instructions = list(dis.get_instructions(code))
    for before, after in zip(instructions, instructions[1:]):
        if before.opname == 'LOAD_CONST' and after.opname == 'STORE_NAME' and (after.argval == '__doc__'):
            docs.add(before.arg)
        if before.opname == 'LOAD_CONST' and after.opname == 'STORE_NAME' and (after.argval == '__firstlineno__'):
            uses = sum((i.opname == 'LOAD_CONST' and i.arg == before.arg for i in instructions))
            if uses == 1:
                locations.add(before.arg)
    return code.replace(co_filename='', co_firstlineno=1, co_linetable=b'', co_consts=tuple((normalized(v) if isinstance(v, types.CodeType) else None if i in docs and isinstance(v, str) else v if i not in locations else 0 for i, v in enumerate(code.co_consts))))

@lru_cache(maxsize=1)
def reader_digest():
    """Reader digest."""
    digest = hashlib.sha256(b'prepared-reader-v2')
    for name in MODULES:
        spec = importlib.util.find_spec(name)
        if spec is None:
            raise ImportError(name)
        code = spec.loader.get_code(name)
        digest.update(name.encode())
        digest.update(marshal.dumps(normalized(code)))
    return digest.hexdigest()

def metadata_file(root):
    """Metadata file."""
    path = Path(root) / 'runtime_metadata.json'
    if not path.is_file():
        import runtime_metadata
        path = Path(runtime_metadata.__file__).resolve().parents[1] / 'resources/sage/runtime_metadata.json'
    return path

@lru_cache(maxsize=16)
def _metadata_digest(path, modified, size, region=None):
    """Metadata digest."""
    raw = Path(path).read_bytes()
    if region is None:
        return hashlib.sha256(raw).hexdigest()
    import json
    data = json.loads(raw)
    for key in ('models', 'bounds', 'variogram_ranges', 'downloads'):
        data.pop(key, None)
    regions = {region}
    if region.startswith('LAMAH_'):
        regions.add('LAMAH_CE')
    if region == 'CAMELS_DEH':
        regions.add('CAMELS_DE')
    for key in ('regions', 'table_profiles', 'checks'):
        if key in data:
            data[key] = {k: v for k, v in data[key].items() if k.split(':', 1)[0] in regions}
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()

def metadata_digest(root, region=None):
    """Metadata digest."""
    generic = Path(root) / 'hydro_schemas.json'
    if region is not None and generic.is_file():
        import json
        from hydro_schema import catalogue
        record = catalogue(str(root)).get(region)
        if record is not None:

            def inputs(value):
                """Inputs."""
                if isinstance(value, dict):
                    return {k: inputs(v) for k, v in value.items() if k not in ('source_file', 'progress')}
                if isinstance(value, list):
                    return [inputs(v) for v in value]
                return value
            selected = {k: record[k] for k in ('meteo', 'discharge')}
            selected['paths'] = {k: v['paths'] for k, v in record['resolutions'].items()}
            return hashlib.sha256(json.dumps(inputs(selected), sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    path = metadata_file(root)
    info = path.stat()
    digest = _metadata_digest(str(path), info.st_mtime_ns, info.st_size, region)
    return digest
