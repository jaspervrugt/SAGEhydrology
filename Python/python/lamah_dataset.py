# SAGEhydrology | Python source
# ==========================================================================
# SITEhydrology | Python source
# File: python/lamah_dataset.py
#
# Separate regional selections and installations for LamaH-CE/Austria.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

# ==========================================================================
# SITEhydrology / SAGEhydrology | Python source
# Country subsets of the published LamaH-CE dataset; source IDs stay intact.
# Copyright (c) 2026 Jasper A. Vrugt | University of California, Irvine
# ==========================================================================

"""Separate regional selections and installations for LamaH-CE/Austria."""

import csv
from functools import lru_cache
import json
import os
from pathlib import Path
import shutil

REGIONS = ("LAMAH_CE", "LAMAH_AT")


def reader_region(region):
    """Reuse the established equations, units and time-series reader."""
    return "LAMAH_CE" if region in REGIONS else region


def resolve_root(root, region, resolution=1):
    """Discover installed LamaH selections."""
    root = Path(root)
    if region not in REGIONS:
        return root
    stream = "daily" if resolution == 1 else "hourly"
    if (root / "Gauge_attributes.csv").is_file() and all(
        (root / stream / "timeseries" / folder).is_dir()
        for folder in ("forcing", "streamflow")
    ):
        return root
    candidates = [root, root.parent / region]
    stream = "daily" if resolution == 1 else "hourly"
    for candidate in candidates:
        # Readers resolve a location cheaply. Full pair validation belongs
        # to installation/readiness, never to each basin evaluation.
        if (candidate / "Gauge_attributes.csv").is_file() and all(
            (candidate / stream / "timeseries" / folder).is_dir()
            for folder in ("forcing", "streamflow")
        ):
            return candidate
    return root


def rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=";"))


def basin_ids(root, region):
    """Intersect country metadata with total-upstream catchment IDs."""
    root = Path(root)
    gauge = root / "Gauge_attributes.csv"
    catch = root / "Catchment_attributes.csv"
    return _basin_ids(
        str(root.resolve()),
        region,
        gauge.stat().st_mtime_ns,
        catch.stat().st_mtime_ns,
    )


@lru_cache(maxsize=32)
def _basin_ids(root, region, gauge_stamp, catch_stamp):
    root = Path(root)
    gauges = rows(root / "Gauge_attributes.csv")
    catchments = {r["ID"] for r in rows(root / "Catchment_attributes.csv")}
    return frozenset(
        {
            r["ID"]
            for r in gauges
            if r["ID"] in catchments
            and (region != "LAMAH_AT" or r["country"].upper() == "AUT")
        }
    )


def available(root, region, resolution=1):
    """Require the exact forcing/discharge pairs, not just file counts."""
    return _complete(
        resolve_root(root, region, resolution), region, resolution
    )


def _complete(root, region, resolution):
    root = Path(root)
    if resolution not in (1, 24):
        return False
    try:
        ids = basin_ids(root, region)
    except (OSError, KeyError, ValueError):
        return False
    if not ids:
        return False
    stream = "daily" if resolution == 1 else "hourly"
    return all(
        (root / stream / "timeseries" / folder / f"ID_{g}.csv").is_file()
        for folder in ("forcing", "streamflow")
        for g in ids
    )


def find_existing(parent, region, resolution=1):
    parent = Path(parent)
    candidates = [parent, parent / "LAMAH_CE"]
    for candidate in candidates:
        candidate = resolve_root(candidate, region, resolution)
        if _complete(candidate, region, resolution):
            return candidate
    return None


def _filter_table(source, target, allowed):
    data = rows(source)
    if not data:
        raise ValueError(f"Empty LamaH metadata: {source}")
    with Path(target).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, data[0].keys(), delimiter=";")
        writer.writeheader()
        writer.writerows(r for r in data if r["ID"] in allowed)


def finalize(root, region, resolution=1):
    """Write subset metadata and geographic coordinates for both readers."""
    root = Path(root)
    allowed = basin_ids(root, region)
    if not allowed:
        raise ValueError("No usable LamaH basins match this selection")
    for name in ("Gauge_attributes.csv", "Catchment_attributes.csv"):
        source = root / name
        staged = source.with_suffix(".subset.tmp")
        _filter_table(source, staged, allowed)
        staged.replace(source)
    from pyproj import Transformer

    projection = Transformer.from_crs(3035, 4326, always_xy=True)
    with (root / "gauge_information.txt").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter=";")
        writer.writerow(
            ["gauge_id", "gauge_name", "gauge_lat", "gauge_lon", "area"]
        )
        areas = {
            r["ID"]: r for r in rows(root / "Catchment_attributes.csv")
        }
        for r in rows(root / "Gauge_attributes.csv"):
            lon, lat = projection.transform(
                float(r["lon"]), float(r["lat"])
            )
            area = areas[r["ID"]].get("area_calc", r.get("area_gov", ""))
            writer.writerow([r["ID"], r["name"], lat, lon, area])
    if not available(root, region, resolution):
        raise ValueError("Incomplete LamaH forcing/discharge pairs")
    manifest = dict(
        dataset=region,
        source_dataset="LamaH-CE",
        source_record="https://zenodo.org/records/5153305",
        country_filter="AUT" if region == "LAMAH_AT" else None,
        basins=len(allowed),
        resolution=resolution,
        resolutions=[n for n in (1, 24) if _complete(root, region, n)],
        status="installed",
    )
    (root / "lamah_installation.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return root


def reuse(source, parent, region, resolution=1, progress=None, stop=None):
    """Install locally without a download; link immutable time series."""
    from dataset_install import _check

    source, parent = Path(source), Path(parent)
    allowed = basin_ids(source, region)
    parent.mkdir(parents=True, exist_ok=True)
    target = parent / region
    if _complete(target, region, resolution):
        return target
    target.mkdir(exist_ok=True)
    for name in ("Gauge_attributes.csv", "Catchment_attributes.csv"):
        if not (target / name).exists():
            shutil.copy2(source / name, target / name)
    stream = "daily" if resolution == 1 else "hourly"
    files = [
        Path(stream) / "timeseries" / folder / f"ID_{g}.csv"
        for folder in ("forcing", "streamflow")
        for g in sorted(allowed, key=int)
    ]
    try:
        for i, relative in enumerate(files):
            _check(stop)
            dest = target / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                if dest.read_bytes() != (source / relative).read_bytes():
                    raise ValueError(f"Conflicting LamaH file: {dest}")
                continue
            try:
                os.link(source / relative, dest)
            except OSError:
                shutil.copy2(source / relative, dest)
            if progress:
                progress("Installing LamaH basin files", i + 1, len(files))
        return finalize(target, region, resolution)
    except Exception as exc:
        (target / "lamah_installation.json").write_text(
            json.dumps(dict(status="incomplete", error=str(exc))),
            encoding="utf-8",
        )
        raise
