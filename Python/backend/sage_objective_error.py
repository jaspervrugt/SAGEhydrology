# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_objective_error.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_objective_error."""

class UndefinedObjective(ValueError):
    """The selected samples do not define this objective or its derivative."""
