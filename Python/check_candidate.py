# ==========================================================================
# SAGEhydrology | Python source
# File: check_candidate.py
#
# Portable headless SAGE entry point.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Run the candidate tests with only this candidate on the module path."""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import run_sage


def main():
    """Discover public numerical tests and return their process status."""
    suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return not result.wasSuccessful()


if __name__ == "__main__":
    raise SystemExit(main())
