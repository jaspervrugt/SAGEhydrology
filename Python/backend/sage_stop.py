# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_stop.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_stop."""
import json
from pathlib import Path
import numpy as np

def requested(directory):
    """Return the pending run-stop action."""
    path = Path(directory) / 'stop_request.json'
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding='utf-8'))

def save_network(directory, arrays, attributes, scaler, selection=None):
    """Save network arrays for the stopped or completed run."""
    from sage_checkpoint import save_npz
    extra = {}
    if selection is not None:
        extra['model_selection_json'] = np.array(json.dumps(selection))
    save_npz(Path(directory) / 'network.npz', **{f'array_{i}': v for i, v in enumerate(arrays)}, attributes=np.array(attributes), mean=np.array(scaler['mean']), scale=np.array(scaler['scale']), median=np.array(scaler['median']), attribute_scaler_json=np.array(json.dumps(scaler, allow_nan=False)), **extra)
