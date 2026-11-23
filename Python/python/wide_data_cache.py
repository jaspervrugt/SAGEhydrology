# ==========================================================================
# SAGEhydrology | Python source
# File: python/wide_data_cache.py
#
# Bounded-memory, source-stamped cache for wide daily numeric tables.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Bounded-memory, source-stamped cache for wide daily numeric tables."""

from pathlib import Path
from functools import lru_cache
from threading import RLock
import hashlib
import json
import os
import shutil
import tempfile
import time

import numpy as np
import pandas as pd

_LOCK = RLock()


def _remove_staging(path, cache):
    """Remove the temporary wide-cache staging directory."""
    target = path.resolve()
    if target.parent != cache.resolve() or not target.name.startswith(
        "tmp_"
    ):
        raise ValueError("Unexpected wide-table staging path")
    shutil.rmtree(target, ignore_errors=True)


def _build(
    path,
    clock,
    cache,
    target,
    signature,
    date=None,
    sequential=False,
    gauges=(),
):
    """Parse selected columns in bounded chunks and publish a disk index.

    Validate the daily clock and source stamp before publishing. A changed
    source cannot become a reusable cache entry.
    """
    staging = Path(tempfile.mkdtemp(prefix="tmp_", dir=cache))
    matrix = None
    try:
        with path.open("rb") as handle:
            rows = (
                sum(
                    block.count(b"\n")
                    for block in iter(
                        lambda: handle.read(4 * 1024 * 1024), b""
                    )
                )
                + 1
            )
        offset = 0
        dates = []
        columns = None
        with path.open(encoding="utf-8-sig") as handle:
            first = handle.readline()
        sep = (
            ";"
            if first.count(";") > first.count(",")
            else "," if "," in first else "\t" if "\t" in first else r"\s+"
        )
        start = None
        wanted = {str(g).lstrip("0") for g in gauges}
        with pd.read_csv(
            path,
            sep=sep,
            engine="c",
            chunksize=1024,
            encoding="utf-8-sig",
            low_memory=False,
            usecols=(
                (
                    lambda c: str(c).strip() in clock
                    or str(c).strip().lstrip("0") in wanted
                )
                if wanted
                else None
            ),
        ) as reader:
            for chunk in reader:
                chunk.columns = [str(c).strip() for c in chunk.columns]
                if columns is None:
                    columns = [c for c in chunk if c not in clock]
                    matrix = np.lib.format.open_memmap(
                        staging / "values.npy",
                        mode="w+",
                        dtype="float64",
                        shape=(rows, len(columns)),
                        fortran_order=True,
                    )
                end = offset + len(chunk)
                numeric = chunk[columns].copy()
                for name in numeric.select_dtypes(exclude=np.number):
                    numeric[name] = pd.to_numeric(
                        numeric[name], errors="coerce"
                    )
                matrix[offset:end] = numeric.to_numpy(dtype=float)
                if date is None:
                    stamp = pd.to_datetime(
                        dict(
                            year=chunk[clock[0]],
                            month=chunk[clock[1]],
                            day=chunk[clock[2]],
                        )
                    ).to_numpy()
                elif sequential:
                    if start is None:
                        start = pd.to_datetime(chunk[date].iloc[0])
                    stamp = pd.date_range(
                        start + pd.Timedelta(days=offset),
                        periods=len(chunk),
                        freq="D",
                    ).to_numpy()
                else:
                    stamp = pd.to_datetime(
                        chunk[date], errors="raise"
                    ).to_numpy()
                dates.append(pd.DatetimeIndex(stamp).normalize().to_numpy())
                offset = end
        if matrix is None:
            raise ValueError(f"Empty wide table: {path}")
        matrix.flush()
        del matrix
        matrix = None
        dates = np.concatenate(dates)
        if len(np.unique(dates)) != len(dates):
            raise ValueError(f"Duplicate timestamps in {path}")
        if (path.stat().st_mtime_ns, path.stat().st_size) != signature:
            raise ValueError(f"Source changed while reading {path}")
        np.save(staging / "dates.npy", dates, allow_pickle=False)
        (staging / "table.json").write_text(
            json.dumps(dict(columns=columns, rows=offset)), encoding="utf-8"
        )
        for attempt in range(8):
            try:
                os.replace(staging, target)
                break
            except OSError:
                if (target / "table.json").is_file():
                    break
                if attempt == 7:
                    raise
                time.sleep(0.05 * (attempt + 1))
    finally:
        if matrix is not None:
            del matrix
        if staging.exists():
            _remove_staging(staging, cache)


@lru_cache(maxsize=8)
def _mapped(folder):
    """Open a completed numeric index as a read-only memory map."""
    folder = Path(folder)
    document = json.loads(
        (folder / "table.json").read_text(encoding="utf-8")
    )
    values = np.load(
        folder / "values.npy", mmap_mode="r", allow_pickle=False
    )
    dates = np.load(folder / "dates.npy", allow_pickle=False)
    return document, values, dates


def series(
    path, gauge, clock=None, *, date=None, sequential=False, gauges=()
):
    """Return one basin without reparsing the regional CSV."""
    path = Path(path).resolve()
    info = path.stat()
    signature = (info.st_mtime_ns, info.st_size)
    clock = tuple(clock) if clock else (date,)
    if not all(clock):
        raise ValueError(
            "A date column or year/month/day columns are required"
        )
    identity = (
        [2, str(path), *signature, list(clock)]
        if date is None and not sequential
        else [3, str(path), *signature, list(clock), date, sequential]
    )
    key = hashlib.sha256(json.dumps(identity).encode()).hexdigest()[:24]
    base = Path(os.environ.get("LOCALAPPDATA", tempfile.gettempdir()))
    cache = base / "SAGEhydrology/Python/wide_tables"
    with _LOCK:
        cache.mkdir(parents=True, exist_ok=True)
    target = cache / key
    # Reuse an existing complete matrix. Otherwise index only this run's
    # selected basins, matching MATLAB's batched wide-column selection.
    gauges = (
        tuple(
            sorted({str(g).strip().lstrip("0") for g in (*gauges, gauge)})
        )
        if gauges
        else ()
    )
    if gauges and not (target / "table.json").is_file():
        key = hashlib.sha256(
            json.dumps([4, identity, gauges]).encode()
        ).hexdigest()[:24]
        target = cache / key
    with _LOCK:
        if not (target / "table.json").is_file():
            _build(
                path,
                clock,
                cache,
                target,
                signature,
                date,
                sequential,
                gauges,
            )
        document, values, dates = _mapped(str(target))
    column = next(
        (
            i
            for i, name in enumerate(document["columns"])
            if name == gauge or name.lstrip("0") == gauge.lstrip("0")
        ),
        None,
    )
    if column is None:
        raise ValueError(f"Basin column {gauge} absent in {path}")
    return pd.Series(
        values[: document["rows"], column].copy(),
        index=pd.DatetimeIndex(dates).normalize(),
    )
