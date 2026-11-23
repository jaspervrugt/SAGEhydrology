# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_run_data.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_run_data."""
from pathlib import Path
import hashlib
import json
import os
from tempfile import mkdtemp
import numpy as np
KEY = '_run_data_directory'
ARRAYS = ('dates', 'forcing', 'observed', 'train', 'evaluation')
TEXT = ('gauge', 'forcing_file', 'streamflow_file')

def create(c):
    """Create an isolated store before screening starts."""
    import tempfile
    base = Path(os.environ.get('LOCALAPPDATA', tempfile.gettempdir()))
    base = base / 'SAGEhydrology/Python/run_inputs'
    base.mkdir(parents=True, exist_ok=True)
    c[KEY] = mkdtemp(prefix='run_', dir=base)

def folder(c, gauge):
    """Folder."""
    return Path(c[KEY]) / hashlib.sha256(str(gauge).encode()).hexdigest()[:24]

def publish(c, gauge, data, targets=None):
    """Publish screened arrays once; metadata is the completion marker."""
    if KEY not in c:
        return
    target = folder(c, gauge)
    target.mkdir(parents=True, exist_ok=True)
    arrays = {key: data[key] for key in ARRAYS}
    for name, values in (targets or {}).items():
        arrays['_target_' + name] = values
    for key, values in arrays.items():
        values = np.asarray(values)
        if values.dtype.hasobject:
            raise ValueError(f'Non-numeric prepared basin input: {key}')
        np.save(target / (key + '.npy'), values, allow_pickle=False)
    document = dict(arrays=list(arrays), text={key: str(data[key]) for key in TEXT if key in data})
    (target / 'basin.json').write_text(json.dumps(document), encoding='utf-8')

def read(c, gauge):
    """Map prepared inputs or fail; never fall back to source readers."""
    target = folder(c, gauge)
    try:
        document = json.loads((target / 'basin.json').read_text(encoding='utf-8'))
        result = dict(document['text'])
        for key in document['arrays']:
            result[key] = np.load(target / (key + '.npy'), mmap_mode='r', allow_pickle=False)
        return result
    except (OSError, ValueError, KeyError) as exc:
        raise RuntimeError(f'Prepared inputs missing or invalid for basin {gauge}; training will not reread dataset files.') from exc

def cleanup(c):
    """Remove only this store after worker shutdown releases mappings."""
    import shutil
    value = c.get(KEY)
    if not value:
        return
    target = Path(value).resolve()
    if target.parent.name != 'run_inputs' or not target.name.startswith('run_'):
        raise ValueError('Unexpected run-input store path')
    try:
        shutil.rmtree(target)
    except OSError:
        pass
