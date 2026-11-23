# ==========================================================================
# SAGEhydrology | Python source
# File: python/hydro_aux.py
#
# Declarative auxiliary metadata operations from read_hydro_timeseries.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Declarative auxiliary metadata operations from read_hydro_timeseries."""

from functools import lru_cache
from pathlib import Path
import re
import numpy as np
from hydro_io import table, numeric, file
from hydro_schema import identifier, items
from local_data import require_local


@lru_cache(maxsize=64)
def _cached_table(path, stamp, size, delimiter):
    """Read and cache an auxiliary metadata table."""
    return table(
        path, {"delimiter": delimiter} if delimiter else {}, static=True
    )


def metadata_table(path, spec):
    """Load the auxiliary metadata table for a reader profile."""
    path = require_local(path)
    info = path.stat()
    return _cached_table(
        str(path), info.st_mtime_ns, info.st_size, spec.get("delimiter")
    )


def auxiliary(root, data_root, gauge, schema, existing=None):
    """Use declared tables first, then series and file-header overrides."""
    aux = np.array(
        existing if existing is not None else [np.nan] * 3, dtype=float
    )
    spec = schema.get("aux", {})
    tables = spec.get("tables", {})
    if spec.get("values"):
        aux[:] = np.asarray(spec["values"]).reshape(-1, 3)[0]
    elif tables and existing is None:
        aux[:] = np.nan
        for entry in tables.values():
            path = (Path(root) / entry["file"].replace("\\", "/")).resolve()
            # GB v1/v2 releases place the same topographic fields either
            # beside daily timeseries or at the dataset root.
            if (
                not path.is_file()
                and path.name.lower()
                == "camels_gb_topographic_attributes.csv"
            ):
                candidates = [
                    base / name
                    for base in (path.parent, Path(data_root))
                    for name in (
                        "CAMELS_GB_topographic_attributes.csv",
                        "camels_gb_v2_topographic_attributes.csv",
                    )
                ]
                path = next((p for p in candidates if p.is_file()), path)
            frame = metadata_table(path, entry)
            request = identifier(gauge, {}, metadata=True)
            if entry.get("request_prefix_length"):
                request = request[: entry["request_prefix_length"]]
            if entry.get("layout") == "wide_attributes":
                columns = [
                    c for c in frame if c != entry["attribute_column"]
                ]
                prefixes = entry.get("column_prefix")
                if prefixes:
                    columns = [
                        c
                        for c in columns
                        if any(
                            c.upper().startswith(p.upper())
                            for p in items(prefixes)
                        )
                    ]
                matches = [
                    c
                    for c in columns
                    if identifier(c, entry, metadata=True) == request
                ]
                if not matches:
                    continue
                for j, field in enumerate(("lat", "elev", "area")):
                    if not np.isfinite(aux[j]) and entry.get(
                        field + "_row"
                    ):
                        selected = (
                            frame[entry["attribute_column"]].str.lower()
                            == entry[field + "_row"].lower()
                        )
                        if not selected.any():
                            raise ValueError(
                                "Missing metadata attribute: "
                                + entry[field + "_row"]
                            )
                        aux[j] = numeric(frame.loc[selected, matches[0]])[
                            0
                        ] * entry.get(field + "_scale", 1)
                continue
            ids = frame[entry["key"]].map(
                lambda x: identifier(x, entry, metadata=True)
            )
            selected = frame[ids == request]
            if selected.empty:
                continue
            for j, field in enumerate(("lat", "elev", "area")):
                if np.isfinite(aux[j]):
                    continue
                if (
                    field == "lat"
                    and entry.get("lat_derive") == "lambert93_to_wgs84"
                ):
                    # Same ellipsoid and Lambert-93 constants as MATLAB.
                    x = numeric(selected[entry["lat_x"]])[0]
                    y = numeric(selected[entry["lat_y"]])[0]
                    n, c, xs, ys = (
                        0.7256077650532670,
                        11754255.426096,
                        700000.0,
                        12655612.049876,
                    )
                    radius = np.hypot(x - xs, y - ys)
                    lat_iso = -np.log(abs(radius / c)) / n
                    eccentricity = 0.0818191910428158
                    phi = 2 * np.arctan(np.exp(lat_iso)) - np.pi / 2
                    for _ in range(12):
                        phi = (
                            2
                            * np.arctan(
                                (
                                    (1 + eccentricity * np.sin(phi))
                                    / (1 - eccentricity * np.sin(phi))
                                )
                                ** (eccentricity / 2)
                                * np.exp(lat_iso)
                            )
                            - np.pi / 2
                        )
                    value = np.rad2deg(phi)
                elif entry.get(field):
                    value = numeric(selected[entry[field]])[0]
                else:
                    continue
                aux[j] = value * entry.get(field + "_scale", 1)
    elif existing is None:
        # MATLAB's caller supplies these bas.zone fields before the reader.
        for path in (Path(data_root) / "gauge_information.txt",):
            if not path.is_file():
                continue
            frame = metadata_table(path, {})
            if "gauge_id" not in frame:
                continue
            clean = lambda x: re.sub(
                r"^0+(?=[0-9])", "", re.sub(r"\.0+$", "", str(x))
            )
            selected = frame[frame.gauge_id.map(clean) == clean(gauge)]
            if not selected.empty:
                for j, column in enumerate(
                    ("gauge_lat", "gauge_elev", "area_km2")
                ):
                    if column in selected:
                        aux[j] = numeric(selected[column])[0]
                if np.isfinite(aux[2]) and aux[2] < 1e7:
                    aux[2] *= 1e6
                break
    for field, entry in spec.get("series", {}).items():
        j = ("lat", "elev", "area").index(field)
        if np.isfinite(aux[j]):
            continue
        path = file(root, entry["file_pattern"], gauge)
        value = numeric(
            table(path, columns=[entry["source"]])[entry["source"]]
        )
        value = np.median(value[np.isfinite(value)])
        if entry.get("reducer") != "median":
            raise ValueError(
                f'Unsupported auxiliary reducer: {entry.get("reducer")}'
            )
        if entry.get("derive") == "pressure_elevation":
            pressure = (
                value / 1000
                if value > 20000
                else value / 10 if value > 200 else value
            )
            value = (
                44330 * (1 - (pressure / 101.325) ** 0.1903)
                if 0 < pressure <= 110
                else np.nan
            )
        aux[j] = value * entry.get("scale", 1)
    header = spec.get("file_header")
    if header:
        path = file(root, schema["file"]["pattern"], gauge)
        with path.open(encoding="utf-8-sig") as stream:
            lines = [stream.readline() for _ in range(max(header["lines"]))]
        for j, line in enumerate(header["lines"]):
            match = re.search(
                r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", lines[line - 1]
            )
            if match:
                aux[j] = float(match[0]) * header.get("scale", [1, 1, 1])[j]
    return aux
