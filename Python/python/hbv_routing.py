# ==========================================================================
# SAGEhydrology | Python source
# File: python/hbv_routing.py
#
# Route discharge and all 13 sensitivities using shared C++ code.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Route discharge and all 13 sensitivities using shared C++ code."""

import ctypes as C
from functools import lru_cache
import numpy as np
from build_hbv_routing import library_path


@lru_cache(maxsize=1)
def library():
    """Load and configure the compiled numerical library."""
    path = library_path()
    if not path.is_file():
        raise FileNotFoundError(
            "Build shared HBV routing: python python/build_hbv_routing.py"
        )
    lib = C.CDLL(str(path.resolve()))
    ptr = C.POINTER(C.c_double)
    lib.sage_hbv_route.argtypes = [
        ptr,
        ptr,
        C.c_int,
        C.c_int,
        C.c_double,
        C.c_int,
        ptr,
        ptr,
    ]
    lib.sage_hbv_route.restype = C.c_int
    lib.sage_hbv_weights.argtypes = [C.c_double, C.c_int, ptr, ptr]
    lib.sage_hbv_weights.restype = C.c_int
    return lib


def _length(value):
    """Resolve the HBV routing length."""
    if not np.isfinite(value) or value < 1 or int(value) != value:
        raise ValueError("Routing length must be a positive integer.")
    if value > np.iinfo(np.int32).max:
        raise ValueError("Routing length exceeds the native integer range.")
    return int(value)


def weights(width, length=60):
    """Return the unchanged smoothed MAXBAS weights and derivatives."""
    length = _length(length)
    if not np.isfinite(width):
        raise ValueError("Routing width must be finite.")
    w, dw = np.empty(length), np.empty(length)
    ptr = C.POINTER(C.c_double)
    status = library().sage_hbv_weights(
        width, length, w.ctypes.data_as(ptr), dw.ctypes.data_as(ptr)
    )
    if status:
        raise ValueError(f"Native HBV weights failed: {status}.")
    return w, dw


def route(q, jacobian, width, length=60, start=0):
    """Preserve the prefix and reset routing at the spin-up boundary."""
    q = np.ascontiguousarray(q, dtype=float)
    if q.ndim != 1 or not len(q) or len(q) > np.iinfo(np.int32).max:
        raise ValueError("Discharge must be a nonempty vector.")
    if (
        not np.isfinite(start)
        or int(start) != start
        or not 0 <= start < len(q)
    ):
        raise ValueError("Invalid HBV routing start.")
    if not np.isfinite(width):
        raise ValueError("Routing width must be finite.")
    jac = None
    result_jac = None
    if jacobian is not None:
        jac = np.asfortranarray(jacobian, dtype=float)
        if jac.shape != (len(q), 12):
            raise ValueError("Unrouted HBV Jacobian must have 12 columns.")
        result_jac = np.empty((len(q), 13), order="F")
    result_q = np.empty_like(q)
    ptr = C.POINTER(C.c_double)
    status = library().sage_hbv_route(
        q.ctypes.data_as(ptr),
        jac.ctypes.data_as(ptr) if jac is not None else None,
        len(q),
        int(start),
        width,
        _length(length),
        result_q.ctypes.data_as(ptr),
        result_jac.ctypes.data_as(ptr) if result_jac is not None else None,
    )
    if status:
        raise ValueError(f"Native HBV routing failed: {status}.")
    return result_q, result_jac
