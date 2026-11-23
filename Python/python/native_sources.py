# ==========================================================================
# SAGEhydrology | Python source
# File: python/native_sources.py
#
# Find source inputs; compiled applications do not need C++ sources.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Find source inputs; compiled applications do not need C++ sources."""

from pathlib import Path


def native_root(sage_root=None):
    """Resolve Native from a checkout or an explicitly supplied path.

    Legacy --sage-root arguments may name Matlab or resources/sage.
    Search their parents for Native; never use a duplicate model core.
    Missing source fails before invoking the compiler.
    """
    starts = [Path(__file__).resolve().parents[2]]
    if sage_root is not None:
        starts.insert(0, Path(sage_root).resolve())
    for start in starts:
        for parent in (start, *start.parents):
            for candidate in (parent / "Native", parent):
                if (candidate / "models/hbv/hbv_routing.hpp").is_file():
                    return candidate
    raise FileNotFoundError(
        "Shared Native sources are missing. Keep Native, Matlab and "
        "Python together in the SAGEhydrology checkout."
    )


def library_directory():
    """Keep checkout build products separate from packaged libraries."""
    import sys

    if getattr(sys, "frozen", False):
        return Path(sys._MEIPASS) / "native/build"
    return (
        Path(__file__).resolve().parents[2]
        / "release_artifacts/native/Python"
        / sys.platform
    )
