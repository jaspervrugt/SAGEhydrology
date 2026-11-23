# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_native.py
#
# Public native CRR adapter.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Seven-model native dispatch, analytic sensitivities, and C++ RSS
gradients.
"""
import ctypes as C
from functools import lru_cache
import json
import numpy as np
from build_crr import library_path
PARAM_COUNTS = {1: 7, 2: 9, 3: 15, 4: 16, 5: 8, 6: 13, 7: 15}
STATE_COUNTS = {1: 7, 2: 6, 3: 10, 4: 9, 5: 10, 6: 5, 7: 14}

@lru_cache(maxsize=1)
def library():
    """Load and configure the compiled numerical library."""
    path = library_path().resolve()
    if not path.is_file():
        raise FileNotFoundError('Build the CRR library: python python/build_crr.py')
    lib = C.CDLL(str(path))
    ptr = C.POINTER(C.c_double)
    lib.site_crr_run.argtypes = [C.c_int] * 4 + [ptr] * 4 + [C.c_int] + [ptr] * 3 + [C.c_int] + [ptr] * 2
    lib.site_crr_run.restype = C.c_int
    lib.site_crr_rss.argtypes = [C.c_int, C.c_int, ptr, ptr, ptr, C.POINTER(C.c_ubyte), ptr, ptr, ptr]
    lib.site_crr_rss.restype = None
    manifest = json.loads((path.parent / 'crr_build.json').read_text())['models']
    return (lib, manifest)

def parameters(model, theta, config=None):
    """Return the selected model parameter specification."""
    t = np.asarray(theta, float)
    p = {}
    cfg = dict(T_sm=1.0, eps_m=1e-06, rho=0.01, eps_s=1e-12, eps_r=1e-12, eps=5.0, a_S=1e-06, kappa=4 / 9, bg=3.5, bR=5.0, mts=1.0, n1=2, n2=4, dT=1.0, K=3)
    if model == 6:
        cfg.update(eps_s=0.001, eps_t=0.1, eps_x=1e-12)
    if config:
        cfg.update(config)
    if model == 1:
        names = ['S_umax', 'beta', 'alfa', 'K_s', 'K_f', 'T_tr', 'f_dd']
    elif model == 2:
        names = ['I_max', 'Su_max', 'Q_max', 'a_E', 'a_F', 'r_f', 'r_s', 'T_tr', 'f_dd']
    elif model == 3:
        names = ['uzfwm', 'uztwm', 'lzfpm', 'lzfsm', 'lztwm', 'zperc', 'rexp', 'uzk', 'pfree', 'lzpk', 'lzsk', 'acm', 'kf', 'T_tr', 'f_dd']
    elif model == 4:
        names = ['f_p', 'A_im', 'a', 'b', 'f_wm', 'f_lm', 'c', 'S_tot', 'Ex', 'k_i', 'k_g', 'c_i', 'c_g', 'k_f', 'T_tr', 'f_dd']
        p.update(W_max=t[4] * t[7], LM=t[5] * t[4] * t[7], S_max=(1 - t[4]) * t[7])
    elif model == 5:
        names = ['x1', 'x2', 'x3', 'x4', 'x5', 'f_p', 'T_tr', 'f_dd']
    elif model == 6:
        names = ['f_c', 'beta', 'lp', 'k_0', 'uzl', 'k_1', 'k_2', 'perc', 'T_tr', 'f_dd', 'sfcf', 'cfr']
    elif model == 7:
        names = ['s_max', 's_fc', 's_wp', 'k_sch', 'a1', 'k_perc', 'lf_thr', 'a2', 'k_lf', 'g_max', 'c_gw', 'mm', 'k_nsh', 'T_tr', 'f_dd']
    else:
        raise ValueError('Unknown model')
    p.update(cfg)
    p.update(zip(names, t))
    if model == 7:
        p['s_fc'] = t[1] * t[0]
        p['s_wp'] = t[2] * t[0]
        p['lf_thr'] = p['s_wp'] + t[6] * (p['s_fc'] - p['s_wp'])
    return p

def maxbas(b, length=60):
    """Shared native MAXBAS weights and routing-width derivatives."""
    from hbv_routing import weights
    return weights(b, length)

def simulate(model, theta, forcing, *, y0=None, model_config=None, settings=(0.01, 1.0, 0.0001, 0.001, 0.001, 2.0, 10000.0), routing_start=0):
    """Simulate the selected model and requested sensitivity outputs."""
    if int(model) == 99:
        from user_model_native import simulate as user_simulate
        return user_simulate(theta, forcing, y0=y0, model_config=model_config, settings=settings, routing_start=routing_start)
    if int(model) == 20:
        from cartplus_basin_inputs import simulate_cartsplit
        return simulate_cartsplit(theta, forcing, model_config)
    if int(model) == 19:
        from cartplus_basin_inputs import simulate_cart
        return simulate_cart(theta, forcing, model_config)
    model = int(model)
    theta = np.asarray(theta, float)
    forcing = np.asfortranarray(forcing, dtype=float)
    s = np.ascontiguousarray(settings, dtype=float)
    if model not in PARAM_COUNTS or theta.shape != (PARAM_COUNTS[model],):
        raise ValueError('Incorrect parameter count')
    if forcing.ndim != 2 or forcing.shape[1] != 3 or (not len(forcing)):
        raise ValueError('Need N-by-3 daily forcing')
    if not all((np.all(np.isfinite(v)) for v in (theta, forcing, s))):
        raise ValueError('Nonfinite model inputs')
    if s.shape != (7,) or min(s[:5]) <= 0 or s[2] > s[1] or (s[5] <= 0) or (not 1 <= s[6] <= np.iinfo(np.int32).max):
        raise ValueError('Invalid solver settings')
    lib, manifest = library()
    p = parameters(model, theta, model_config)
    giuh = np.ascontiguousarray((model_config or {}).get('giuh_ordnts', [0.1, 0.35, 0.2, 0.14, 0.1, 0.06, 0.05]), dtype=float)
    if np.any(giuh < 0) or not np.all(np.isfinite(giuh)) or (not np.isclose(giuh.sum(), 1)):
        raise ValueError('GIUH ordinates must sum to one')
    p['L'] = len(giuh)
    m = STATE_COUNTS[model] if model != 7 else 4 + int(p['K']) + len(giuh)
    initial = np.full(m, 1e-05) if y0 is None else np.ascontiguousarray(y0, dtype=float)
    if initial.shape != (m,) or not np.all(np.isfinite(initial)):
        raise ValueError(f'Model {model} requires {m} initial states')
    packed = np.ascontiguousarray([p[name] for name in manifest[str(model)]['fields']], dtype=float)
    if not np.all(np.isfinite(packed)):
        raise ValueError('Nonfinite model configuration')
    d = int(manifest[str(model)]['core_parameters'])
    q = np.empty(len(forcing))
    J = np.empty((len(forcing), d), order='F')
    ptr = C.POINTER(C.c_double)

    def address(v):
        """Return a native pointer to the supplied numeric array."""
        return v.ctypes.data_as(ptr)
    status = lib.site_crr_run(model, len(forcing), m, d, *[address(forcing[:, i]) for i in range(3)], address(packed), len(packed), address(s), address(initial), address(giuh), len(giuh), address(q), address(J))
    if status:
        from simulation_failure import IntegrationFailure
        raise IntegrationFailure(f'Native model {model} failed with status {status}')
    if model == 6:
        from hbv_routing import route
        q, J = route(q, J, theta[12], (model_config or {}).get('routing_length', 60), routing_start)
    return (q, J)

def rss_objective(model, forcing, observed, train, *, lower, upper, **kwargs):
    """Evaluate residual sum of squares and its analytic gradient."""
    y = np.ascontiguousarray(observed, float)
    mask = np.ascontiguousarray(np.asarray(train, bool) & np.isfinite(y) & (y >= 0), dtype=np.uint8)
    lo, hi = (np.asarray(lower, float), np.asarray(upper, float))
    scale = np.ascontiguousarray(hi - lo)
    if y.shape != (len(forcing),) or mask.shape != y.shape or (not mask.any()) or (lo.shape != (PARAM_COUNTS[model],)) or np.any(scale <= 0):
        raise ValueError('Invalid objective arrays/bounds')
    lib, _ = library()
    ptr = C.POINTER(C.c_double)

    def addr(v):
        """Return a native pointer to the supplied numeric array."""
        return v.ctypes.data_as(ptr)

    def fun(x):
        """Evaluate the local objective callback."""
        try:
            q, J = simulate(model, lo + np.asarray(x) * scale, forcing, **kwargs)
        except __import__('simulation_failure').IntegrationFailure:
            raise
        except ValueError:
            return (1e+20, np.zeros(len(lo)))
        f = C.c_double()
        g = np.empty(len(lo))
        J = np.asfortranarray(J)
        lib.site_crr_rss(len(q), len(lo), addr(q), addr(J), addr(y), mask.ctypes.data_as(C.POINTER(C.c_ubyte)), addr(scale), C.byref(f), addr(g))
        return (f.value, g)
    return fun
