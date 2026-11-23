# ==========================================================================
# SAGEhydrology | Python source
# File: python/native_hymod.py
#
# ctypes wrapper: Python calls the same HYMOD C++ source used by SAGE
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
ctypes wrapper: Python calls the same HYMOD C++ source used by SAGE
MEX.
"""

import ctypes
from functools import lru_cache
import numpy as np
from build_native import library_path


@lru_cache(maxsize=1)
def _library():
    """Load and configure the compiled numerical library."""
    path = library_path().resolve()
    if not path.is_file():
        raise FileNotFoundError(
            "Build native HYMOD first: python python/build_native.py"
            f" (missing {path})"
        )
    lib = ctypes.CDLL(str(path))
    pointer = ctypes.POINTER(ctypes.c_double)
    lib.site_hymod_run.argtypes = [ctypes.c_int] + [pointer] * 9
    lib.site_hymod_run.restype = ctypes.c_int
    return lib


def integrate(theta, forcing, y0, config, settings):
    """
    Validated input arrays enter here; outputs are discharge and
    native J.
    """
    if (
        len(forcing) < 1
        or len(forcing) > np.iinfo(np.int32).max
        or settings[6] > np.iinfo(np.int32).max
    ):
        raise ValueError(
            "C++ day count/iteration limit must fit a positive int32"
        )
    # C++ expects contiguous forcing vectors and column-major Jacobians.
    forcing = np.asfortranarray(forcing, dtype=np.float64)
    arrays = [forcing[:, i] for i in range(3)]
    arrays += [
        np.ascontiguousarray(v, dtype=np.float64)
        for v in (theta, config, settings, y0)
    ]
    q = np.empty(len(forcing), dtype=np.float64)
    jac = np.empty((len(forcing), 7), dtype=np.float64, order="F")
    pointer = ctypes.POINTER(ctypes.c_double)
    arguments = [v.ctypes.data_as(pointer) for v in (*arrays, q, jac)]
    status = _library().site_hymod_run(len(forcing), *arguments)
    if status:
        raise ValueError(
            f"Native HYMOD integration failed (status {status})"
        )
    return q, jac
