# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_parameter_scaling.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_parameter_scaling."""
import numpy as np
_CART_LOG_INDICES = {18: (0, 2, 3, 5, 7, 9, 11, 12, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 31, 32), 19: (0, 2, 3, 5, 7, 9, 10), 20: (0, 2, 3, 5, 7, 9, 10)}

def cart_log_indices(model):
    """Cart log indices."""
    return _CART_LOG_INDICES.get(int(model), ())

def modes(c):
    """Resolve each parameter scaling mode."""
    n = len(c['lower'])
    values = c.get('parameter_scales')
    if values is None:
        return ['linear'] * n
    if not isinstance(values, (list, tuple)) or len(values) != n:
        raise ValueError(f'Parameter Scale must contain {n} Linear/Log choices.')
    result = [str(value).lower() for value in values]
    encoded = cart_log_indices(c.get('model', 0))
    if any((value not in ('linear', 'log') and (not (value == 'physical_linear' and index in encoded)) for index, value in enumerate(result))):
        raise ValueError('Parameter Scale must be Linear or Log.')
    return result

def validate(c):
    """Validate parameter scaling and physical bounds."""
    lo = np.asarray(c['lower'], float)
    hi = np.asarray(c['upper'], float)
    if lo.ndim != 1 or lo.shape != hi.shape or (not np.isfinite(lo).all()) or (not np.isfinite(hi).all()) or np.any(hi <= lo):
        raise ValueError('Parameter bounds must be finite, with Min < Max.')
    selected = modes(c)
    log = np.array([value == 'log' for value in selected], bool)
    if np.any(log & (lo <= 0)):
        rows = ', '.join((str(i + 1) for i in np.flatnonzero(log & (lo <= 0))))
        raise ValueError(f'Log scale requires Min > 0 (parameter rows {rows}).')
    physical_linear = np.asarray([value == 'physical_linear' for value in selected])
    if np.any(physical_linear):
        with np.errstate(over='ignore', under='ignore'):
            physical_lo, physical_hi = (np.exp(lo[physical_linear]), np.exp(hi[physical_linear]))
        if np.any(physical_lo <= 0) or np.any(~np.isfinite(physical_hi)) or np.any(physical_hi <= physical_lo):
            raise ValueError('CaRT physical bounds must be finite, positive and distinct.')
    return (lo, hi, log)

def physical(c, x):
    """Map normalized parameters to physical values.

    Parameters
    ----------
    c : dict
        Physical lower/upper bounds and each parameter's Linear/Log mode.
    x : array-like
        Normalized parameters with the parameter dimension on the last axis.

    Returns
    -------
    ndarray
        Physical parameters with the input shape and the configured mapping.
    """
    lo, hi, log = validate(c)
    x = np.asarray(x, float)
    if x.ndim < 1 or x.shape[-1] != len(lo):
        raise ValueError('Normalized parameter dimensions do not match bounds.')
    result = lo + x * (hi - lo)
    result[..., log] = _log_physical(lo[log], hi[log], x[..., log])
    for index, scale in enumerate(modes(c)):
        if scale == 'physical_linear':
            result[..., index] = _cart_physical_linear(lo[index], hi[index], x[..., index])
    return result

def normalized(c, p):
    """Map physical parameters to normalized values."""
    lo, hi, log = validate(c)
    p = np.asarray(p, float)
    if p.ndim < 1 or p.shape[-1] != len(lo):
        raise ValueError('Physical parameter dimensions do not match bounds.')
    if np.any(p[..., log] <= 0):
        raise ValueError('Physical parameters on Log scale must be positive.')
    result = (p - lo) / (hi - lo)
    result[..., log] = (np.log(p[..., log]) - np.log(lo[log])) / (np.log(hi[log]) - np.log(lo[log]))
    for index, scale in enumerate(modes(c)):
        if scale == 'physical_linear':
            result[..., index] = _cart_normalized(lo[index], hi[index], p[..., index])
    return result

def derivative(c, p):
    """dp/dx evaluated at physical p, for the model Jacobian chain rule."""
    lo, hi, log = validate(c)
    p = np.asarray(p, float)
    result = np.broadcast_to(hi - lo, p.shape).copy()
    result[..., log] = p[..., log] * (np.log(hi[log]) - np.log(lo[log]))
    for index, scale in enumerate(modes(c)):
        if scale == 'physical_linear':
            result[..., index] = (np.exp(hi[index]) - np.exp(lo[index])) / np.exp(p[..., index])
    return result

def physical_column(c, index, x):
    """Map one normalized parameter column to physical values."""
    lo, hi, log = validate(c)
    x = np.asarray(x, float)
    if modes(c)[index] == 'physical_linear':
        return _cart_physical_linear(lo[index], hi[index], x)
    return _log_physical(lo[index], hi[index], x) if log[index] else lo[index] + x * (hi[index] - lo[index])

def _cart_physical_linear(lo, hi, x):
    """Cart physical linear."""
    value = np.log((1 - x) * np.exp(lo) + x * np.exp(hi))
    return np.where(x == 0, lo, np.where(x == 1, hi, value))

def _cart_normalized(lo, hi, value):
    """Cart normalized."""
    lower = np.exp(lo)
    return (np.exp(value) - lower) / (np.exp(hi) - lower)

def _log_physical(lo, hi, x):
    """Evaluate the logarithmic parameter back-transformation."""
    values = lo * np.exp(x * (np.log(hi) - np.log(lo)))
    values = np.where((x >= 0) & (x <= 1), np.minimum(np.maximum(values, lo), hi), values)
    return np.where(x == 0, lo, np.where(x == 1, hi, values))
