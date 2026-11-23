# ==========================================================================
# SAGEhydrology | Python source
# File: python/local_data.py
#
# Fail clearly when Windows cloud placeholders would stall data readers.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Fail clearly when Windows cloud placeholders would stall data readers."""

from pathlib import Path
from contextlib import contextmanager
from contextvars import ContextVar

_SOURCES = ContextVar("reader_source_files", default=None)
_GAUGES = ContextVar("reader_basin_selection", default=())


def selected_gauges():
    """Return basin IDs selected by the current reader context."""
    return _GAUGES.get()


@contextmanager
def capture_sources(gauges=()):
    """
    Record every declared source read, including cache hits and NetCDF.
    """
    sources = set()
    token = _SOURCES.set(sources)
    gauge_token = _GAUGES.set(tuple(map(str, gauges)))
    try:
        yield sources
    finally:
        _SOURCES.reset(token)
        _GAUGES.reset(gauge_token)


class DataNotLocalError(FileNotFoundError):
    pass


def require_local(path):
    """Record a source dependency and reject online-only cloud files.

    Raise DataNotLocalError before a provider download can stall a read.
    """
    file = Path(path)
    info = file.stat()
    # FILE_ATTRIBUTE_RECALL_ON_DATA_ACCESS: provider must retrieve content.
    if getattr(info, "st_file_attributes", 0) & 0x400000:
        raise DataNotLocalError(
            f"Data file is online-only or partly downloaded: {file}. In"
            " Dropbox/File Explorer, make this file available offline"
            " before reading it."
        )
    sources = _SOURCES.get()
    if sources is not None:
        sources.add(str(file.resolve()))
    return file


from functools import lru_cache


@lru_cache(maxsize=64)
def _attribute_bytes(path, mtime_ns, size):
    """Read bytes from a local catchment-attribute file."""
    return Path(path).read_bytes()


def attribute_bytes(path):
    """
    Reuse static metadata bytes until the file changes; never cache
    forcing.
    """
    file = require_local(path)
    info = file.stat()
    return _attribute_bytes(
        str(file.resolve()), info.st_mtime_ns, info.st_size
    )
