# ==========================================================================
# SAGEhydrology | Python source
# File: python/sage_paths.py
#
# Sage paths.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

from pathlib import Path
import os, sys

FROZEN = bool(getattr(sys, "frozen", False))
PYTHON_ROOT = Path(
    getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[1])
)
SOFTWARE_ROOT = PYTHON_ROOT.parent.parent
if FROZEN:
    SOFTWARE_ROOT = next(
        (
            p
            for p in Path(sys.executable).resolve().parents
            if (p / "Data").is_dir()
        ),
        SOFTWARE_ROOT,
    )
if os.environ.get("SAGE_SOFTWARE_ROOT"):
    SOFTWARE_ROOT = (
        Path(os.environ["SAGE_SOFTWARE_ROOT"]).expanduser().resolve()
    )
OUTPUT_ROOT = (
    (
        Path(os.environ.get("LOCALAPPDATA", Path.home() / ".local/share"))
        / "SAGEhydrology/Python/results"
    )
    if FROZEN
    else PYTHON_ROOT.parent / "results/Python"
)
SAGE_RESOURCES = PYTHON_ROOT / "resources" / "sage"
DATA_ROOT = SOFTWARE_ROOT / "Data"


def relocate(value):
    """
    Translate legacy installation/result paths without changing
    saved files.
    """
    if isinstance(value, dict):
        result = {key: relocate(item) for key, item in value.items()}
        if "sage_root" in result:
            result["sage_root"] = str(SAGE_RESOURCES)
        if "site_root" in result:
            result["site_root"] = str(PYTHON_ROOT)
        return result
    if isinstance(value, list):
        return [relocate(item) for item in value]
    if isinstance(value, tuple):
        return tuple(relocate(item) for item in value)
    if not isinstance(value, str):
        return value
    normalized = value.replace("\\", "/")
    sage = str(PYTHON_ROOT.parent).replace("\\", "/")
    legacy = str(
        SOFTWARE_ROOT / "SITEhydrology/gui_python/sage_runs"
    ).replace("\\", "/")
    project = PYTHON_ROOT.parent
    moved = (
        (legacy, OUTPUT_ROOT / "legacy_SITE_runs"),
        (str(project / "Python/results"), OUTPUT_ROOT),
        (str(project / "Matlab/results"), project / "results/Matlab"),
    )
    for old, new in moved:
        old = old.replace("\\", "/")
        if (
            normalized.lower() == old.lower()
            or normalized.lower().startswith(old.lower() + "/")
        ):
            return str(new / normalized[len(old) :].lstrip("/"))
    if normalized.lower() == sage.lower():
        return str(PYTHON_ROOT.parent / "Matlab")
    if normalized.lower().startswith(sage.lower() + "/"):
        rest = normalized[len(sage) + 1 :]
        roots = {
            "matlab",
            "python",
            "native",
            "user_model",
            "results",
            "benchmarks",
            "release_artifacts",
            "maintenance_archive",
            "local_dependencies",
        }
        if rest.split("/")[0].lower() not in roots:
            return str(PYTHON_ROOT.parent / "Matlab" / rest)
    return value


def load_result(path):
    """Load a saved run and relocate its stored paths."""
    import json

    return relocate(
        json.loads(
            Path(relocate(str(path))).read_text(encoding="utf-8-sig")
        )
    )
