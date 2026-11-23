# ==========================================================================
# SAGEhydrology | Python source
# File: python/runtime_metadata.py
#
# Versioned, language-neutral metadata owned by the Python distribution.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Versioned, language-neutral metadata owned by the Python distribution."""

from pathlib import Path
from functools import lru_cache
from copy import deepcopy
import json


def clean_unit(value):
    """
    Display readable temperature and metric powers without changing
    data.
    """
    import re

    text = str(value)
    # Repeated mojibake keeps the final degree glyph but adds a growing
    # non-ASCII prefix. Match only that corrupt prefix before C or F;
    # preserve multiplication signs, exponents and the time unit.
    text = re.sub(
        r"[\u00c3\u00c2\u00e2][^\x00-\x7f]*\u00b0(?=[CF]\b)",
        "\u00b0",
        text,
    )
    text = text.replace("\ufffdC", "\u00b0C")
    text = re.sub(
        r"\bdeg(?:ree)?s?\s*([CF])\b",
        lambda m: "°" + m[1].upper(),
        text,
        flags=re.I,
    )
    text = re.sub(r"\bdeg\b", "°", text)
    # Handle explicit integer powers, then legacy compact metric units.
    superscripts = str.maketrans("0123456789+-−", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁻")
    text = re.sub(
        r"\^(?:\{([+−-]?\d+)\}|([+−-]?\d+))(?![\d.])",
        lambda m: (m[1] or m[2]).translate(superscripts),
        text,
    )
    text = re.sub(
        r"\b(km|cm|mm|dm|m)(-1|-?[23])(?![A-Za-z0-9.])",
        lambda m: m[1] + m[2].translate(superscripts),
        text,
    )
    return text.replace("mm/day", "mm/d").replace("mm/°C/day", "mm/d/°C")


def display_unit(value):
    """Use slash notation for rates shown in model tables."""
    import re

    text = clean_unit(value)
    text = text.replace("mm/°C/day", "mm/d/°C")
    text = re.sub(r"/day\b", "/d", text)
    text = re.sub(r"\s+(d|day|h|min)⁻¹(?!\w)", lambda match: "/" + match[1].replace("day", "d"), text)
    return re.sub(r"^(d|day|h|min)⁻¹$", lambda match: "1/" + match[1].replace("day", "d"), text)


@lru_cache(maxsize=4)
def load(sage_root):
    """Load the selected data or metadata resource."""
    path = Path(sage_root) / "runtime_metadata.json"
    if not path.is_file():
        path = (
            Path(__file__).resolve().parents[1]
            / "resources"
            / "sage"
            / "runtime_metadata.json"
        )
    data = (
        json.loads(path.read_text(encoding="utf-8"))
        if path.is_file()
        else None
    )
    for model in (data or {}).get("models", {}).values():
        for parameter in model.get("parameters", []):
            parameter["unit"] = clean_unit(parameter.get("unit", ""))
    return data


def regions(sage_root):
    """Return the runtime regional dataset metadata."""
    data = load(str(sage_root))
    if data is None:
        return None
    result = deepcopy(data["regions"])
    resource = Path(sage_root)
    if not (resource / "runtime_metadata.json").is_file():
        resource = (
            Path(__file__).resolve().parents[1] / "resources" / "sage"
        )
    for meta in result.values():
        meta["source"] = str(resource / meta["source"])
    return result


def variogram_range(sage_root, region):
    """Return the configured regional variogram range."""
    data = load(str(sage_root))
    return (
        data["variogram_ranges"].get(region, 1500.0)
        if data
        else (3750.0 if region == "CAMELS_US" else 1500.0)
    )
