# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_observables.py
#
# Public native CRR adapter.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Q/SWE/SM and analytic Jacobians, following SAGE src/crr_result.m."""
import ctypes as C
from functools import lru_cache
import json
import numpy as np
from build_observables import library_path
from crr_native import parameters, STATE_COUNTS, PARAM_COUNTS

@lru_cache(maxsize=1)
def library():
    """Load and configure the compiled numerical library."""
    path = library_path()
    if not path.is_file():
        raise FileNotFoundError('Build named observations: python python/build_observables.py')
    lib = C.CDLL(str(path.resolve()))
    ptr = C.POINTER(C.c_double)
    lib.sage_observables_run.argtypes = [C.c_int] * 4 + [ptr] * 4 + [C.c_int] + [ptr] * 3 + [C.c_int, C.POINTER(C.c_int), C.c_int, ptr, ptr]
    lib.sage_observables_run.restype = C.c_int
    return (lib, json.loads(path.with_suffix('.json').read_text())['models'])

def simulate(model, theta, forcing, *, y0=None, model_config=None, settings=(0.01, 1.0, 0.0001, 0.001, 0.001, 2.0, 10000.0), routing_start=0):
    """Simulate the selected model and requested sensitivity outputs."""
    if int(model) == 99:
        from user_model_native import simulate as user_simulate
        q, jacobian = user_simulate(theta, forcing, y0=y0, model_config=model_config, settings=settings, routing_start=routing_start)
        return ({'Q': q}, {'Q': jacobian})
    model = int(model)
    theta = np.asarray(theta, float)
    f = np.asfortranarray(forcing, dtype=float)
    cfg = model_config or {}
    s = np.ascontiguousarray(settings, float)
    if theta.shape != (PARAM_COUNTS[model],) or f.ndim != 2 or f.shape[1] != 3:
        raise ValueError('Incorrect named-observation input dimensions.')
    if not np.isfinite(theta).all() or not np.isfinite(f).all() or (not np.isfinite(s).all()):
        raise ValueError('Nonfinite native inputs.')
    if s.shape != (7,) or min(s[:5]) <= 0 or s[2] > s[1] or (s[5] <= 0) or (not 1 <= s[6] <= np.iinfo(np.int32).max):
        raise ValueError('Invalid solver settings.')
    lib, manifest = library()
    p = parameters(model, theta, cfg)
    giuh = np.ascontiguousarray(cfg.get('giuh_ordnts', [0.1, 0.35, 0.2, 0.14, 0.1, 0.06, 0.05]), float)
    if np.any(giuh < 0) or not np.isfinite(giuh).all() or (not np.isclose(giuh.sum(), 1)):
        raise ValueError('Invalid GIUH ordinates.')
    p['L'] = len(giuh)
    m = STATE_COUNTS[model] if model != 7 else 4 + int(p['K']) + len(giuh)
    initial = np.full(m, 1e-05) if y0 is None else np.ascontiguousarray(y0, float)
    if initial.shape != (m,):
        raise ValueError('Incorrect initial-state dimensions.')
    if not np.isfinite(initial).all():
        raise ValueError('Nonfinite initial states.')
    packed = np.ascontiguousarray([p[name] for name in manifest[str(model)]['fields']], float)
    d = 12 if model == 6 else PARAM_COUNTS[model]
    n = len(f)
    obs = np.empty((n, 3), order='F')
    J = np.empty((n, d, 3), order='F')
    soil = np.ascontiguousarray(np.asarray(cfg.get('sm_ind', []), np.int32) - 1)
    nsoil = len(soil) if 'sm_ind' in cfg else -1
    ptr = C.POINTER(C.c_double)
    address = lambda a: a.ctypes.data_as(ptr)
    status = lib.sage_observables_run(model, n, m, d, *[address(f[:, i]) for i in range(3)], address(packed), len(packed), address(s), address(initial), address(giuh), len(giuh), soil.ctypes.data_as(C.POINTER(C.c_int)), nsoil, address(obs), address(J))
    if status:
        from simulation_failure import IntegrationFailure
        raise IntegrationFailure(f'Named native model failed with status {status}.')
    if model == 6:
        if not 0 <= routing_start < n:
            raise ValueError('Invalid HBV routing start.')
        jac = np.zeros((n, 13, 3), order='F')
        jac[:, :12, :] = J
        from hbv_routing import route
        obs[:, 0], jac[:, :, 0] = route(obs[:, 0], J[:, :, 0], theta[12], cfg.get('routing_length', 60), routing_start)
        J = jac
    names = ['Q', 'SWE'] + (['SM'] if model != 1 or nsoil > 0 else [])
    return ({name: obs[:, i] for i, name in enumerate(names)}, {name: J[:, :, i] for i, name in enumerate(names)})
