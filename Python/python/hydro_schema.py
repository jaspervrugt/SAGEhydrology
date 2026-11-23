# ==========================================================================
# SAGEhydrology | Python source
# File: python/hydro_schema.py
#
# Resolve exported MATLAB hydro schemas without regional reader code.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Resolve exported MATLAB hydro schemas without regional reader code."""

from copy import deepcopy
from functools import lru_cache
from pathlib import Path
import json
import re
import sys


def catalogue(root):
    """Return the hydrologic reader catalogue."""
    path = Path(root or "") / "hydro_schemas.json"
    if not path.is_file():
        base = Path(
            getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[1])
        )
        path = base / "resources/sage/hydro_schemas.json"
    info = path.stat()
    return _catalogue(str(path), info.st_mtime_ns, info.st_size)


@lru_cache(maxsize=4)
def _catalogue(path, modified, size):
    """Load the hydrologic reader catalogue."""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def merge(base, override):
    """Merge reader-profile settings."""
    result = deepcopy(base)
    for key, value in override.items():
        result[key] = (
            merge(result[key], value)
            if isinstance(result.get(key), dict) and isinstance(value, dict)
            else deepcopy(value)
        )
    return result


def items(value):
    """Normalize a reader-profile collection to its entries."""
    return value if isinstance(value, list) else [value]


def resolve(schema, options):
    """Port resolve_hydro_schema's recursive merge and unique match rule."""
    schema = deepcopy(schema)
    if "variants" in schema:
        key = "p{pet}_t{temperature}_r{precipitation}".format(**options)
        if key not in schema["variants"]:
            raise ValueError(f"Unsupported meteorological selection: {key}")
        schema = schema["variants"][key]
    if "profiles" in schema:
        profiles = schema.pop("profiles")
        matches = [
            (name, p["schema"])
            for name, p in profiles.items()
            if all(
                options.get(k) in items(v) for k, v in p["match"].items()
            )
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one schema profile; found {len(matches)}"
            )
        name, override = matches[0]
        schema = merge(schema, override)
        schema["selected_profile"] = name
    schema.setdefault("layout", "one_file_per_basin")
    schema.setdefault("format", "csv")
    schema.setdefault("missing_values", [-9999, -999])
    schema.setdefault("strict_time", True)
    schema.setdefault("file", {})
    schema["file"].setdefault("delimiter", ",")
    schema["file"].setdefault("header_lines", 0)
    schema.setdefault("id", {})
    schema.setdefault("aux", {})
    schema.setdefault("time", {})
    schema["timeline"].setdefault("step", 86400)
    for spec in schema["variables"].values():
        spec.setdefault("missing_values", schema["missing_values"])
    return schema


def identifier(value, spec, metadata=False):
    """Normalize a reader-profile identifier."""
    value = str(value).strip()
    if spec.get("strip_decimal", True):
        value = re.sub(r"\.0+$", "", value)
    if metadata:
        value = value.upper()
        prefix = spec.get("strip_prefix", "")
        if prefix:
            value = re.sub("^" + re.escape(prefix), "", value, flags=re.I)
    if spec.get("strip_leading_zeros"):
        value = value.lstrip("0") or "0"
    value = value.zfill(int(spec.get("pad_width", 0)))
    prefix = spec.get("prefix", "")
    if prefix and not value.startswith(prefix):
        value = prefix + value
    return value


def configuration(
    root, region, resolution, product, pet, temperature, precipitation
):
    """Resolve the reader profile for the requested configuration."""
    record = catalogue(str(root or ""))[region]
    product_text = str(product).strip().lower()
    number = re.match(r"^(\d+)", product_text)
    data = (
        int(number[1])
        if number
        else {"daymet": 1, "maurer": 2, "nldas": 3}.get(product_text, 1)
    )
    options = dict(
        dt=resolution,
        pet=pet,
        temp=temperature,
        temperature=temperature,
        precip=precipitation,
        precipitation=precipitation,
        data=data,
    )
    resolutions = [
        x for x in record["resolutions"].values() if x["dt"] == resolution
    ]
    if len(resolutions) != 1:
        raise ValueError(f"Unsupported {region} resolution: {resolution}")
    return (
        resolve(record["meteo"], options),
        resolve(record["discharge"], dict(dt=resolution)),
        resolutions[0]["paths"],
        options,
    )
