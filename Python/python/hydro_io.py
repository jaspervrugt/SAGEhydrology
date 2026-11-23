# ==========================================================================
# SAGEhydrology | Python source
# File: python/hydro_io.py
#
# MATLAB generic reader storage layouts, dates, and native source rules.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""MATLAB generic reader storage layouts, dates, and native source rules."""

from functools import lru_cache
from collections import OrderedDict
from io import BytesIO
from pathlib import Path
import glob
import re
import json
import hashlib
from threading import RLock
import numpy as np
import pandas as pd
from hydro_schema import identifier, items
from local_data import require_local, attribute_bytes, selected_gauges

_WIDE_TABLES = OrderedDict()
_WIDE_LOCK = RLock()
_WIDE_BYTES = 64 * 1024 * 1024
_CSV_WINDOW_MIN_BYTES = 8_000_000
_CLOCK_TABLES = OrderedDict()
_CLOCK_BYTES = 32 * 1024 * 1024


def file(root, pattern, gauge):
    """Resolve the gauge-specific input filename."""
    padded = f"{int(gauge):02d}" if str(gauge).isdigit() else str(gauge)
    pattern = pattern.replace("{gauge02}", padded).replace(
        "{gauge}", str(gauge)
    )
    path = str(Path(root) / pattern.replace("\\", "/"))
    matches = glob.glob(path) if "*" in path or "?" in path else [path]
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one source for {path}; found {len(matches)}"
        )
    return require_local(Path(matches[0]).resolve())


def table(
    path,
    spec=None,
    columns=None,
    static=False,
    numeric_columns=(),
    nrows=None,
    row_start=0,
    source_names=None,
    integer_columns=(),
):
    """Read the selected source table."""
    spec = spec or {}
    delimiter = spec.get("delimiter")
    if delimiter is None:
        with Path(path).open(encoding="utf-8-sig") as stream:
            line = stream.readline()
        delimiter = "\t" if "\t" in line else ";" if ";" in line else ","
    delimiter = "\t" if delimiter == "\\t" else delimiter
    if isinstance(delimiter, list):
        delimiter = r"\s+"
    names = spec.get("variable_names") or None
    skip = int(spec.get("header_lines", 0))
    if row_start:
        skip += row_start + (0 if names else 1)
        names = names or list(source_names)
    kwargs = dict(
        sep=delimiter,
        skiprows=skip,
        encoding="utf-8-sig",
        usecols=columns,
        names=names,
        float_precision="round_trip",
        low_memory=False,
        nrows=nrows,
        na_values=["--", "missing", " ", "  ", "   ", "    "],
        skipinitialspace=True,
    )
    payload = (
        BytesIO(attribute_bytes(path)) if static else require_local(path)
    )
    # Keep IDs as strings; scientific columns are converted explicitly.
    dtype = (
        {name: float for name in numeric_columns}
        if numeric_columns
        else str
    )
    if integer_columns:
        dtype = dict(dtype) if isinstance(dtype, dict) else {}
        dtype.update({name: np.int64 for name in integer_columns})
    return pd.read_csv(payload, **kwargs, dtype=dtype)


def numeric(value):
    """Convert source values to a numeric array."""
    source = np.asarray(value)
    if np.issubdtype(source.dtype, np.number):
        return np.asarray(source, dtype=float).copy()
    text = pd.Series([str(x).strip() for x in source.ravel()], dtype=object)
    missing = text.str.lower().isin(
        [
            "",
            "nan",
            "na",
            "n/a",
            "null",
            "none",
            "-nan",
            "missing",
            "--",
            "<na>",
        ]
    )
    try:
        result = text.where(~missing, "nan").to_numpy().astype(float)
    except ValueError as error:
        raise ValueError(
            "Unexpected text in numeric hydro column"
        ) from error
    if np.any(np.isnan(result) & ~missing.to_numpy()):
        bad = text[np.isnan(result) & ~missing.to_numpy()].iloc[0]
        raise ValueError(
            f"Unexpected text in numeric hydro column: {bad!r}"
        )
    return result


def native(value, spec):
    """Read the source using its native column conventions."""
    value = numeric(value)
    value[
        np.isin(value, items(spec.get("missing_values", [-9999, -999])))
    ] = np.nan
    if spec.get("invalid_le") is not None and spec["invalid_le"] != []:
        value[value <= float(spec["invalid_le"])] = np.nan
    return value


def timestamps(frame, spec):
    """Construct timestamps from the source date columns."""
    mode = spec.get("mode", "column")
    if mode == "column":
        fmt = spec.get("input_format", "")
        mapping = [
            ("yyyy", "%Y"),
            ("MMM", "%b"),
            ("MM", "%m"),
            ("dd", "%d"),
            ("HH", "%H"),
            ("hh", "%H"),
            ("mm", "%M"),
            ("ss", "%S"),
        ]
        for matlab, python in mapping:
            fmt = fmt.replace(matlab, python)
        clock = pd.to_datetime(frame[spec["column"]], format=fmt or "mixed")
        return clock.to_numpy(dtype="datetime64[s]")

    def component(name):
        """Read one configured forcing or discharge component."""
        value = np.asarray(frame[spec[name + "_column"]])
        if value.dtype.kind == "b":
            return numeric(value)
        try:
            return value.astype(float)
        except (ValueError, TypeError):
            return numeric(value)

    components = {k: component(k) for k in ("year", "month", "day")}
    for k in ("hour", "minute", "second"):
        if spec.get(k + "_column"):
            components[k] = component(k)
    clock = calendar_clock(components)
    if clock is None:
        clock = pd.to_datetime(components).to_numpy(dtype="datetime64[s]")
    if mode == "ymd_row_sequence":
        steps = int(spec["steps_per_day"])
        starts = np.flatnonzero(np.r_[True, clock[1:] != clock[:-1]])
        for block, (begin, end) in enumerate(
            zip(starts, np.r_[starts[1:], len(clock)])
        ):
            count = end - begin
            if count > steps:
                raise ValueError("Too many rows per day")
            offset = (
                steps - count
                if block == 0 and spec.get("partial_first") == "final_slots"
                else 0
            )
            clock[begin:end] += (
                np.arange(offset, offset + count) * (86400 // steps)
            ).astype("timedelta64[s]")
    return clock


def calendar_clock(parts):
    """Decode ordinary integral calendar components without table assembly.

    Fractional, missing, invalid, or overflowing components retain pandas'
    original validation and normalization path.
    """
    limits = dict(
        year=(1, 9999),
        month=(1, 12),
        day=(1, 31),
        hour=(0, 23),
        minute=(0, 59),
        second=(0, 59),
    )
    for name, value in parts.items():
        low, high = limits[name]
        if np.any(~np.isfinite(value)) or np.any(
            (value < low) | (value > high) | (value != np.floor(value))
        ):
            return None
    index = (parts["year"].astype("int64") - 1970) * 12
    index += parts["month"].astype("int64") - 1
    first = index.astype("datetime64[M]").astype("datetime64[D]")
    following = (index + 1).astype("datetime64[M]").astype("datetime64[D]")
    day = parts["day"].astype("int64")
    if np.any(day > (following - first).astype("int64")):
        return None
    clock = (first + (day - 1).astype("timedelta64[D]")).astype(
        "datetime64[s]"
    )
    for name, scale in (("hour", 3600), ("minute", 60), ("second", 1)):
        if name in parts:
            clock += (parts[name].astype("int64") * scale).astype(
                "timedelta64[s]"
            )
    return clock


def align(source, target, values, schema):
    """Align source values to the requested reference timestamps."""
    if schema.get("strict_time", True) and np.any(
        np.diff(source).astype("int64") <= 0
    ):
        raise ValueError(
            "Hydro source timestamps are not strictly increasing"
        )
    if schema["time"].get("join") == "calendar_day":
        source, target = source.astype("datetime64[D]"), target.astype(
            "datetime64[D]"
        )
    unique, first = np.unique(source, return_index=True)
    position = np.searchsorted(unique, target)
    found = position < len(unique)
    found[found] &= unique[position[found]] == target[found]
    index = np.full(len(target), -1, dtype=int)
    index[found] = first[position[found]]
    result = np.full(len(target), np.nan)
    found = index >= 0
    result[found] = values[index[found]]
    return result


def sources(spec, available):
    """Resolve the configured input data sources."""
    requested = (
        [spec["source"]] if spec.get("source") else spec.get("sources", [])
    )
    present = [x for x in requested if x in available]
    if requested and not present:
        raise ValueError(f"Missing hydro source columns: {requested}")
    if spec.get("source_operation", "aliases") == "aliases":
        present = present[:1]
    return present


def combine(columns, spec):
    """Combine source columns according to the reader profile."""
    if len(columns) == 1:
        return columns[0]
    matrix = np.column_stack(columns)
    operation = spec.get("source_operation", "aliases")
    if operation == "coalesce":
        result = np.full(len(matrix), np.nan)
        for column in columns:
            use = ~np.isfinite(result) & np.isfinite(column)
            result[use] = column[use]
        return result
    if operation == "mean":
        finite = np.isfinite(matrix)
        count = finite.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(finite, matrix, 0).sum(axis=1) / count
    raise ValueError(f"Unsupported source operation: {operation}")


def quality_mask(get, spec, path, size):
    """Construct the source-specific data-quality mask."""
    quality = spec.get("quality")
    valid = np.ones(size, dtype=bool)
    if not quality:
        return valid
    names = quality["sources"]
    for prefix, selected in quality.get("file_prefix_sources", {}).items():
        if path.stem.startswith(prefix + "_"):
            names = selected
            break
    for name in items(names):
        value = numeric(get(name))
        if quality.get("require_finite", True):
            valid &= np.isfinite(value)
        valid &= (value >= float(quality.get("valid_min", -np.inf))) & (
            value <= float(quality.get("valid_max", np.inf))
        )
    return valid


def csv_raw(path, schema, target):
    """Read raw values from a delimited source file."""
    with path.open(encoding="utf-8-sig") as stream:
        if not any(any(c.isalnum() for c in line) for line in stream):
            return {
                k: np.full(len(target), np.nan)
                for k, v in schema["variables"].items()
                if v.get("source") or v.get("sources")
            }
    header = table(path, schema["file"], nrows=0)
    if path.name.startswith("CAMELS_GB_hydromet_timeseries_"):
        # GB v2 retains the v1 CEH-GEAR/CHESS products with explicit names.
        # Resolve only absent legacy names; never replace an explicit
        # choice.
        from copy import deepcopy

        aliases = {
            "precipitation": "precipitation_cehgear",
            "temperature": "temperature_chess",
            "pet": "pet_chess",
            "peti": "peti_chess",
        }
        schema = deepcopy(schema)
        for spec in schema["variables"].values():
            source = spec.get("source")
            if (
                source in aliases
                and source not in header.columns
                and aliases[source] in header.columns
            ):
                spec["source"] = aliases[source]
    numeric_columns = set()
    for spec in schema["variables"].values():
        numeric_columns.update(sources(spec, header.columns))
        quality = spec.get("quality")
        if quality:
            numeric_columns.update(items(quality["sources"]))
    ts = schema["time"]
    time_columns = (
        [ts["column"]]
        if ts.get("mode", "column") == "column"
        else [v for k, v in ts.items() if k.endswith("_column")]
    )
    window = csv_window(path, schema, target, header, time_columns)
    if window is not None:
        clock, begin, count = window
        if count == 0:
            return {
                k: np.full(len(target), np.nan)
                for k, v in schema["variables"].items()
                if sources(v, header.columns)
            }
        frame = table(
            path,
            schema["file"],
            columns=list(numeric_columns),
            numeric_columns=numeric_columns,
            row_start=begin,
            nrows=count,
            source_names=header.columns,
        )
    else:
        frame = table(
            path,
            schema["file"],
            columns=list(numeric_columns | set(time_columns)),
            numeric_columns=numeric_columns,
        )
        clock = timestamps(frame, schema["time"])
    result = {}
    for field, spec in schema["variables"].items():
        selected = sources(spec, frame.columns)
        if not selected:
            continue
        value = combine(
            [native(frame[name], spec) for name in selected], spec
        )
        value[~quality_mask(frame.__getitem__, spec, path, len(frame))] = (
            np.nan
        )
        result[field] = align(clock, target, value, schema)
    return result


def csv_window(path, schema, target, header, time_columns):
    """Parse dates once, then convert only requested numeric CSV rows.

    Retain full timestamp validation and the exact round-trip conversion.
    Physical-line counting guards against blank and multiline CSV records;
    these layouts retain the general reader. No regular spacing is assumed.
    """
    if (
        path.stat().st_size < _CSV_WINDOW_MIN_BYTES
        or not schema.get("strict_time", True)
        or not len(target)
    ):
        return None
    # skiprows counts physical lines, whereas pandas omits blank records
    # and joins quoted multiline records. Use row offsets only when these
    # physical and logical record counts agree.
    raw = path.read_bytes()
    lines = raw.count(b"\n") + bool(raw and not raw.endswith(b"\n"))
    headers = int(schema["file"].get("header_lines", 0))
    headers += not bool(schema["file"].get("variable_names"))
    clock = None
    fmt = schema["time"].get("input_format")
    delimiter = schema["file"].get("delimiter", ",")
    if (
        len(time_columns) == 1
        and time_columns[0] == header.columns[0]
        and fmt in (
            "yyyy-MM-dd",
            "yyyy-MM-dd HH:mm:ss",
            "dd-MMM-yyyy HH:mm:ss",
        )
        and isinstance(delimiter, str)
        and len(delimiter) == 1
    ):
        records = raw.splitlines()[headers:]
        clock = fixed_clock(records, fmt, delimiter)
    if clock is None:
        if schema["time"].get("mode", "column") != "column":
            try:
                frame = table(
                    path,
                    schema["file"],
                    columns=time_columns,
                    integer_columns=time_columns,
                )
            except (ValueError, TypeError):
                frame = table(path, schema["file"], columns=time_columns)
        else:
            frame = table(path, schema["file"], columns=time_columns)
        clock = timestamps(frame, schema["time"])
    if lines != len(clock) + headers:
        return None
    if np.any(np.diff(clock).astype("int64") <= 0):
        raise ValueError(
            "Hydro source timestamps are not strictly increasing"
        )
    source, request = clock, target
    if schema["time"].get("join") == "calendar_day":
        source = clock.astype("datetime64[D]")
        request = target.astype("datetime64[D]")
    begin = int(np.searchsorted(source, request.min(), side="left"))
    end = int(np.searchsorted(source, request.max(), side="right"))
    return clock[begin:end], begin, end - begin


def fixed_clock(records, fmt, delimiter):
    """Decode schema-declared fixed dates with vectorized validation."""
    if not records:
        return None
    named = fmt == "dd-MMM-yyyy HH:mm:ss"
    width = 20 if named else 19 if "HH" in fmt else 10
    quoted = records[0].startswith(b'"')
    span = width + 1 + 2 * quoted
    prefix = np.asarray([row[:span] for row in records], dtype=f"S{span}")
    key = (fmt, delimiter, hashlib.sha256(memoryview(prefix)).digest())
    with _WIDE_LOCK:
        if key in _CLOCK_TABLES:
            _CLOCK_TABLES.move_to_end(key)
            return _CLOCK_TABLES[key]
    fields = prefix.view("uint8").reshape(-1, span)
    if quoted:
        if not (
            np.all(fields[:, 0] == 34)
            and np.all(fields[:, width + 1] == 34)
            and np.all(fields[:, width + 2] == ord(delimiter))
        ):
            return None
        chars = fields[:, 1 : width + 1]
    else:
        if not np.all(fields[:, width] == ord(delimiter)):
            return None
        chars = fields[:, :width]
    separators = (
        {2: 45, 6: 45, 11: 32, 14: 58, 17: 58} if named else {4: 45, 7: 45}
    )
    if not named and width == 19:
        separators.update({10: 32, 13: 58, 16: 58})
    digits = [
        i
        for i in range(width)
        if i not in separators and (not named or i not in (3, 4, 5))
    ]
    valid = all(
        np.all(chars[:, i] == value) for i, value in separators.items()
    ) and np.all((chars[:, digits] >= 48) & (chars[:, digits] <= 57))
    if not valid:
        return None
    if named:
        tokens = np.bitwise_or(chars[:, 3:6], 32).copy().view("S3").ravel()
        months = np.zeros(len(chars), dtype="uint8")
        for number, name in enumerate(
            (
                b"jan",
                b"feb",
                b"mar",
                b"apr",
                b"may",
                b"jun",
                b"jul",
                b"aug",
                b"sep",
                b"oct",
                b"nov",
                b"dec",
            ),
            1,
        ):
            months[tokens == name] = number
        if np.any(months == 0):
            return None
        iso = np.empty((len(chars), 19), dtype="uint8")
        iso[:, :4] = chars[:, 7:11]
        iso[:, 4] = iso[:, 7] = 45
        iso[:, 5] = months // 10 + 48
        iso[:, 6] = months % 10 + 48
        iso[:, 8:10] = chars[:, :2]
        iso[:, 10] = 32
        iso[:, 11:] = chars[:, 12:20]
        clock = iso.view("S19").ravel().astype("datetime64[s]")
    else:
        clock = (
            chars[:, :width]
            .copy()
            .view(f"S{width}")
            .ravel()
            .astype("datetime64[s]")
        )
    # Separate meteo variables often repeat the same complete clock.
    # Reuse only byte-identical, already validated timestamp columns.
    clock.flags.writeable = False
    if clock.nbytes <= _CLOCK_BYTES:
        with _WIDE_LOCK:
            _CLOCK_TABLES[key] = clock
            while (
                sum(v.nbytes for v in _CLOCK_TABLES.values()) > _CLOCK_BYTES
            ):
                _CLOCK_TABLES.popitem(last=False)
    return clock


def netcdf_raw(path, schema, target):
    """Read raw values from a NetCDF source file."""
    from netcdf_io import dataset

    with dataset(path, buffered=True) as nc:
        ts = schema["time"]
        if ts["mode"] == "index":
            step = int(schema["timeline"]["step"])
            start = (
                int(
                    (
                        target[0] - np.datetime64(ts["file_start"], "s")
                    ).astype("int64")
                )
                // step
            )
            if start < 0:
                raise ValueError(
                    "Requested window starts before NetCDF source"
                )
            section = slice(start, start + len(target))
            clock = target
        else:
            names = list(nc.variables)
            preferred = [ts["column"], "time", "date", "index"]
            name = next(
                (
                    n
                    for p in preferred
                    for n in names
                    if n.lower() == p.lower()
                ),
                None,
            )
            if name is None:
                name = next(
                    (
                        n
                        for n in names
                        if " since " in getattr(nc[n], "units", "")
                    ),
                    None,
                )
            if name is None and ts.get("implicit_index"):
                values = np.arange(
                    max(len(x) for x in nc.dimensions.values())
                )
                units = ts["units"] + " since " + ts["origin"]
            elif name is None:
                raise ValueError(f"Missing NetCDF time coordinate: {path}")
            else:
                values = np.asarray(nc[name][:]).ravel()
                units = getattr(
                    nc[name],
                    "units",
                    ts.get("units", "days")
                    + " since "
                    + ts.get("origin", "1970-01-01"),
                )
            match = re.match(
                r"\s*(seconds?|minutes?|hours?|days?)\s+since\s+(.+)",
                units,
                re.I,
            )
            if not match:
                raise ValueError(f"Unsupported NetCDF time units: {units}")
            seconds = {"s": 1, "m": 60, "h": 3600, "d": 86400}[
                match[1][0].lower()
            ]
            origin = np.datetime64(
                pd.Timestamp(match[2]).to_datetime64(), "s"
            )
            clock = origin + np.rint(values * seconds).astype(
                "timedelta64[s]"
            )
            if ts.get("snap"):
                clock = clock.astype(
                    "datetime64["
                    + {"day": "D", "hour": "h"}[ts["snap"]]
                    + "]"
                ).astype("datetime64[s]")
            section = slice(None)

        def get(name):
            """Retrieve the requested reader input."""
            value = nc[name][:].squeeze()[section]
            return np.asarray(
                np.ma.filled(value.astype(float), np.nan)
            ).ravel()

        result = {}
        for field, spec in schema["variables"].items():
            selected = sources(spec, nc.variables)
            if selected:
                value = combine(
                    [native(get(name), spec) for name in selected], spec
                )
                value[~quality_mask(get, spec, path, len(value))] = np.nan
                result[field] = align(clock, target, value, schema)
        return result


def read_raw(root, gauge, schema, target, options):
    """Read source values using the configured storage format."""
    layout = schema["layout"]
    if layout == "one_file_per_basin":
        path = file(root, schema["file"]["pattern"], gauge)
        loader = csv_raw if schema["format"] == "csv" else netcdf_raw
        return loader(path, schema, target), [path]
    if layout == "multi_file_per_basin":
        result, paths = {}, []
        for name, fs in schema["files"].items():
            selected = {
                k: v
                for k, v in schema["variables"].items()
                if v.get("file") == name
            }
            if not selected:
                continue
            sub = dict(
                schema,
                layout="one_file_per_basin",
                file=fs,
                format=fs.get("format", schema["format"]),
                variables=selected,
                time=dict(schema["time"], **fs.get("time", {})),
            )
            part, used = read_raw(root, gauge, sub, target, options)
            if result.keys() & part.keys():
                raise ValueError(
                    "Hydro variable belongs to multiple named files"
                )
            result.update(part)
            paths.extend(used)
        return result, paths
    if layout != "wide_files":
        raise ValueError(f"Unknown hydro layout: {layout}")
    result, paths = {}, []
    for field, spec in schema["variables"].items():
        if not spec.get("files"):
            continue
        choices = items(spec["files"])
        choice = int(
            options.get(spec.get("selector"), spec.get("default", 1))
        )
        if not 1 <= choice <= len(choices):
            raise ValueError(
                f"Invalid wide-file selection for {field}: {choice}"
            )
        path = file(root, choices[choice - 1], gauge)
        fs = schema["file"]
        # Read only the requested basin, exactly as MATLAB's wide layout.
        with path.open(encoding="utf-8-sig") as stream:
            names = (
                stream.readline().strip().split(fs.get("delimiter", ","))
            )
        id_spec = dict(
            schema.get("id", {}),
            strip_prefix=schema.get("id", {}).get("prefix", ""),
        )
        request = identifier(gauge, id_spec, metadata=True)
        normalized = {
            identifier(n.strip().strip('"'), id_spec, metadata=True): (
                n.strip().strip('"')
            )
            for n in names
        }
        column = normalized.get(request)
        if column is None:
            raise ValueError(f"Missing wide basin column {gauge} in {path}")
        ts = dict(schema["time"])
        if ts["mode"] == "column":
            dates = [ts["column"]]
        else:
            dates = [ts[k + "_column"] for k in ("year", "month", "day")]
        frame = wide_table(path, fs, dates, column, normalized, id_spec)
        clock = timestamps(frame, ts)
        if spec.get("time_mode") == "row_sequence_from_first":
            clock = clock[0] + (
                np.arange(len(clock)) * int(schema["timeline"]["step"])
            ).astype("timedelta64[s]")
        value = native(frame[column], spec).astype(np.float32).astype(float)
        result[field] = align(clock, target, value, schema)
        paths.append(path)
    return result, paths


def wide_table(path, spec, dates, column, available, id_spec):
    """Read a bounded batch of selected basin columns once per source.

    MATLAB reads the selected wide-file columns together. Python's basin
    coordinator calls the reader individually, so retain that shared batch
    across calls. Source stamps invalidate it; canonical disk caches still
    provide restart reuse. Retain at most 64 MiB per process.
    """
    wanted = {column}
    for gauge in selected_gauges():
        name = available.get(identifier(gauge, id_spec, metadata=True))
        if name is not None:
            wanted.add(name)
    ordered = sorted(wanted)
    begin = ordered.index(column) // 32 * 32
    columns = ordered[begin : begin + 32]
    info = require_local(path).stat()
    key = (
        str(path),
        info.st_mtime_ns,
        info.st_size,
        json.dumps(spec, sort_keys=True),
        tuple(dates),
        tuple(columns),
    )
    with _WIDE_LOCK:
        if key in _WIDE_TABLES:
            _WIDE_TABLES.move_to_end(key)
            return _WIDE_TABLES[key][0]
        frame = table(
            path, spec, columns=dates + columns, numeric_columns=columns
        )
        size = int(frame.memory_usage(deep=True).sum())
        if size <= _WIDE_BYTES:
            _WIDE_TABLES[key] = (frame, size)
            while sum(v[1] for v in _WIDE_TABLES.values()) > _WIDE_BYTES:
                _WIDE_TABLES.popitem(last=False)
        return frame
