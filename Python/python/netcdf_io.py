# ==========================================================================
# SAGEhydrology | Python source
# File: python/netcdf_io.py
#
# Keep NetCDF open, reads, conversions and close under one process lock.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Keep NetCDF open, reads, conversions and close under one process lock."""

from contextlib import contextmanager
from threading import RLock
from pathlib import Path

_LOCK = RLock()


@contextmanager
def dataset(path, buffered=False):
    """Yield a dataset while excluding other native NetCDF reader calls."""
    # netCDF-C is not thread safe, even for distinct read-only files.
    # Independent worker processes retain their own lock and parallelism.
    # Ordinary file I/O is thread safe. Keep its provider waits outside
    # the native library lock, with a strict per-file memory bound.
    payload = None
    if buffered:
        source = Path(path)
        if source.stat().st_size <= 64 << 20:
            payload = source.read_bytes()
    with _LOCK:
        from netCDF4 import Dataset

        handle = (
            Dataset(Path(path).name, memory=payload)
            if payload is not None
            else Dataset(path)
        )
        with handle:
            yield handle
