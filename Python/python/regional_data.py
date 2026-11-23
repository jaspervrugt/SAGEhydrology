# ==========================================================================
# SAGEhydrology | Python source
# File: python/regional_data.py
#
# Regional inventory and the shared schema-driven input reader.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Regional inventory and the shared schema-driven input reader.

Production input reads consume exported MATLAB schema objects. Legacy
schema-text helpers remain for inventory and installer compatibility.
"""

from datetime import datetime
from functools import lru_cache
from pathlib import Path
from camels_spat_io import inputs as spat_inputs
import re
import numpy as np
import pandas as pd

GRDC = {"AR", "BE", "EE", "IE", "JM", "NO", "PR", "NA", "ZA"}


def clean_source(path):
    """Read declarative MATLAB text without executing it.

    Strip comments while retaining percent signs inside quoted values.
    """
    text = Path(path).read_text(encoding="utf-8-sig")
    lines = []
    for line in text.splitlines():
        quoted = False
        end = len(line)
        i = 0
        while i < len(line):
            if line[i] == "'":
                if quoted and i + 1 < len(line) and line[i + 1] == "'":
                    i += 2
                    continue
                quoted = not quoted
            if line[i] == "%" and not quoted:
                end = i
                break
            i += 1
        lines.append(line[:end])
    return re.sub(r"\.\.\.\s*\n", "", "\n".join(lines))


def strings(text):
    """Extract MATLAB single-quoted strings from a schema fragment."""
    return re.findall(r"'([^']*)'", text)


@lru_cache(maxsize=4)
def region_registry(sage_root):
    """Return installed region metadata keyed by canonical region code."""
    result = {}
    for path in sorted(
        (Path(sage_root) / "regions").rglob("region_config_*.m")
    ):
        text = clean_source(path)

        def value(key, default=""):
            """Resolve a selected regional option value."""
            match = re.search(r"R\." + key + r"\s*=\s*'([^']*)'", text)
            if match:
                return match.group(1)
            block = re.search(r"R\s*=\s*struct\((.*?)\);", text, re.S)
            match = (
                re.search("'" + key + r"'\s*,\s*'([^']*)'", block.group(1))
                if block
                else None
            )
            return match.group(1) if match else default

        code = value("code")
        name = value("name")
        acronym = value("acronym")
        inherited = re.search(
            r"region_config_GRDC\(\s*'([^']*)'\s*,\s*'([^']*)'", text
        )
        if inherited:
            acronym, name = inherited.groups()
            code = "CAMELS_" + acronym
            text = (
                clean_source(
                    Path(sage_root)
                    / "regions"
                    / "shared"
                    / "region_config_GRDC.m"
                )
                + "\n"
                + text
            )
        if not code:
            continue
        resolution_match = re.search(
            r"R\.resolutions\s*=\s*\{(.*?)\}", text, re.S
        )
        resolutions = (
            strings(resolution_match.group(1))
            if resolution_match
            else ["Daily"]
        )
        profiles = {}
        # X is repeatedly reset before each resolution definition.
        blocks = [
            (block, False)
            for block in re.split(r"X\s*=\s*struct\(\)\s*;", text)[1:] or [
                text
            ]
        ]
        # CAMELS-CZ inherits its hourly controls from X and overrides them
        # in H.
        for derived in re.finditer(
            r"(\w+)\s*=\s*X\s*;(.*?)R\.by_resolution\.\w+\s*=\s*\1\s*;",
            text,
            re.S,
        ):
            blocks.append(
                (
                    re.sub(
                        r"\b" + derived.group(1) + r"\.",
                        "X.",
                        derived.group(2),
                    ),
                    True,
                )
            )
        for block, inherits in blocks:
            label = re.search(r"X\.label\s*=\s*'([^']*)'", block)
            if not label:
                continue
            controls = {}
            for key in ("product", "precipitation", "temperature", "pet"):
                match = re.search(
                    r"X\.meteo\." + key + r"\.items\s*=\s*\{(.*?)\};",
                    block,
                    re.S,
                )
                struct_control = re.search(
                    r"X\.meteo\." + key + r"\s*=\s*struct\((.*?)\);",
                    block,
                    re.S,
                )
                if not match and struct_control:
                    match = re.search(
                        r"'items'\s*,\s*\{\{(.*?)\}\}",
                        struct_control.group(1),
                        re.S,
                    )
                default = re.search(
                    r"X\.meteo\." + key + r"\.default\s*=\s*'([^']*)'",
                    block,
                )
                base = (
                    profiles.get("Daily", {}).get("meteo", {}).get(key, {})
                    if inherits
                    else {}
                )
                enabled = re.search(
                    r"X\.meteo\." + key + r"\.enabled\s*=\s*(true|false)",
                    block,
                )
                controls[key] = dict(
                    items=(
                        strings(match.group(1))
                        if match
                        else base.get("items", [])
                    ),
                    enabled=(
                        enabled.group(1) == "true"
                        if enabled
                        else base.get("enabled", False)
                    ),
                    default=(
                        default.group(1)
                        if default
                        else base.get("default", "")
                    ),
                )
            dates = (
                dict(profiles.get("Daily", {}).get("dates", {}))
                if inherits
                else {}
            )
            for key in (
                "train_start",
                "train_end",
                "eval_start",
                "eval_end",
            ):
                match = re.search(
                    r"X\.period\.manual\." + key + r"\s*=\s*'([^']*)'",
                    block,
                )
                if not match:
                    period_struct = re.search(
                        r"X\.period\.manual\s*=\s*struct\((.*?)\);",
                        block,
                        re.S,
                    )
                    if period_struct:
                        match = re.search(
                            "'" + key + r"'\s*,\s*'([^']*)'",
                            period_struct.group(1),
                        )
                if match:
                    dates[key] = (
                        datetime.strptime(match.group(1), "%d/%m/%Y")
                        .date()
                        .isoformat()
                    )
            profiles[label.group(1)] = dict(meteo=controls, dates=dates)
        if inherited:
            date_values = re.findall(
                r"'([0-9]{2}/[0-9]{2}/[0-9]{4})'", text
            )
            if len(date_values) >= 4:
                profiles["Daily"]["dates"] = dict(
                    zip(
                        (
                            "eval_start",
                            "eval_end",
                            "train_start",
                            "train_end",
                        ),
                        [
                            datetime.strptime(d, "%d/%m/%Y")
                            .date()
                            .isoformat()
                            for d in date_values[-4:]
                        ],
                    )
                )
        result[code] = dict(
            name=name,
            dataset=value("dataset", code.replace("_", "-")),
            default_resolution=value("default_resolution", resolutions[0]),
            acronym=acronym,
            data_root=value("data_root", code),
            source=str(path.resolve()),
            resolutions=resolutions,
            profiles=profiles,
            meteo=next(iter(profiles.values()), {}).get("meteo", {}),
        )
    if "CAMELS_SPAT" in result and "CAMELS_CA" in result:
        from copy import deepcopy

        result["CAMELS_SPAT"].update(
            dataset="CAMELS-SPAT",
            resolutions=result["CAMELS_CA"]["resolutions"][:],
            profiles=deepcopy(result["CAMELS_CA"]["profiles"]),
            meteo=deepcopy(result["CAMELS_CA"]["meteo"]),
        )
        text = clean_source(result["CAMELS_SPAT"]["source"])
        for key in ("train_start", "train_end", "eval_start", "eval_end"):
            match = re.search(
                r"X\.period\.manual\." + key + r"\s*=\s*'([^']*)'", text
            )
            if match:
                result["CAMELS_SPAT"]["profiles"]["Hourly"]["dates"][
                    key
                ] = (
                    datetime.strptime(match.group(1), "%d/%m/%Y")
                    .date()
                    .isoformat()
                )
    return result


def _section(source, name="meteo"):
    """Extract a named section from regional metadata."""
    text = clean_source(source)
    match = re.search(
        r"function\s+S\s*=\s*local_" + name + r"_schema\(\)", text
    )
    if not match:
        match = re.search(
            r"function\s+\[.*?\]\s*=\s*local_schemas\(\)", text
        )
    if not match and "R.code = 'CAMELS_KR'" in text:
        return text
    if not match:
        raise ValueError(f"No declarative {name} schema in {source}")
    return re.split(r"\nfunction\s", text[match.end() :])[0]


def _literal(text, key, default=None):
    """Parse a literal value in regional metadata."""
    match = re.search(re.escape(key) + r"\s*=\s*'([^']*)'", text)
    if match:
        return match.group(1)
    match = re.search(
        re.escape(key) + r"\s*=\s*fullfile\((.*?)\);", text, re.S
    )
    return str(Path(*strings(match.group(1)))) if match else default


def table_profile(sage_root, region, pet=1, temperature=1, precipitation=1):
    """Resolve the selected PET, temperature and precipitation schema.

    Prefer bundled runtime metadata; retain the legacy schema fallback.
    Reject unsupported forcing selections instead of substituting data.
    """
    from lamah_dataset import reader_region

    region = reader_region(region)
    registry = region_registry(str(Path(sage_root).resolve()))
    meta = registry[region]
    text = _section(
        registry["CAMELS_CA"]["source"]
        if region == "CAMELS_SPAT"
        else meta["source"]
    )
    variables = {}
    for match in re.finditer(
        r"S\.variables\.(\w+)\s*=\s*(?:local_variable|variable)\((.*?)\);",
        text,
        re.S,
    ):
        args = strings(match.group(2))
        if args:
            variables.setdefault(
                match.group(1),
                dict(
                    source=args[0], units=args[1] if len(args) > 1 else ""
                ),
            )
    for match in re.finditer(
        r"S\.variables\.(\w+)\.sources\s*=\s*\{(.*?)\};", text, re.S
    ):
        key = match.group(1)
        variables.setdefault(key, {})["sources"] = strings(match.group(2))
        variables[key]["operation"] = _literal(
            text, "S.variables." + key + ".source_operation", "aliases"
        )
    for match in re.finditer(
        r"S\.variables\.(\w+)\.derive\s*=\s*\'(.*?)\'", text
    ):
        variables.setdefault(match.group(1), {})["derive"] = match.group(2)
    for match in re.finditer(
        r"S\.variables\.(\w+)\.units\s*=\s*'([^']*)'", text
    ):
        variables.setdefault(match.group(1), {})["units"] = match.group(2)
    pattern = _literal(text, "S.file.pattern")
    profile = dict(
        meta=meta,
        text=text,
        variables=variables,
        pattern=pattern,
        layout=_literal(text, "S.layout"),
        format=_literal(text, "S.format"),
        date=_literal(text, "S.time.column"),
        date_format=_literal(text, "S.time.input_format", ""),
        ymd=[
            _literal(text, "S.time." + key + "_column")
            for key in ("year", "month", "day")
        ],
    )

    def source(key, column, units="mm/day"):
        """Resolve the selected regional data source."""
        variables[key] = dict(source=column, units=units)

    def choices(key, columns, index):
        """Return available regional input choices."""
        if not 1 <= int(index) <= len(columns):
            raise ValueError(f"Invalid {region} {key} selection")
        source(
            key,
            columns[int(index) - 1],
            "degC" if key == "T" else ("mm/day"),
        )

    if region == "CAMELS_GB":
        profile.update(
            layout="one_file_per_basin",
            pattern="CAMELS_GB_hydromet_timeseries_{gauge}_*.csv",
        )
        source("P", "precipitation")
        source("T", "temperature", "degC")
        source("Q", "discharge_spec")
        choices("Ep", ["pet", "peti"], pet or 1)
    elif region == "CAMELS_FR":
        profile["skiprows"] = 7
        choices("Ep", ["tsd_pet_pm", "tsd_pet_pe", "tsd_pet_ou"], pet or 1)
        choices(
            "T",
            ["tsd_temp", "tsd_temp", "tsd_temp_min", "tsd_temp_max"],
            temperature or 1,
        )
        if temperature == 2:
            variables["T"] = dict(derive="mean_tmin_tmax")
            source("Tmin", "tsd_temp_min")
            source("Tmax", "tsd_temp_max")
    elif region == "CAMELS_FI":
        variables["T"] = dict(
            sources=[
                ["temperature_mean", "temperature_gmin"],
                ["temperature_min", "temperature_mean"],
                ["temperature_max", "temperature_mean"],
                ["temperature_gmin", "temperature_mean"],
            ][(temperature or 1) - 1]
        )
        variables["Ep"] = dict(
            sources=[
                ["pet", "pet_fmi", "pet_singer"],
                ["pet_fmi", "pet", "pet_singer"],
                ["pet_singer", "pet", "pet_fmi"],
            ][(pet or 1) - 1],
            fill_missing=0,
            clip_min=0,
        )
    elif region == "CAMELS_LUX":
        profile["pattern"] = (
            "CAMELS_LUX_hydromet_timeseries__daily_ID_{gauge02}.csv"
        )
        choices(
            "P",
            ["RR_rad", "RR_stn", "tp", "RR_min_rad", "RR_max_rad"],
            precipitation or 1,
        )
        choices("T", ["t2m", "T_stn"], temperature or 1)
        choices("Ep", ["PET_PM", "PET_Oudin"], pet or 1)
    elif region == "EStreams_PT":
        match = re.search(r"pet\s*=\s*\{(.*?)\}", text, re.S)
        choices("Ep", strings(match.group(1)), pet or 1)
    elif region == "CAMELS_CZ":
        profile["pattern"] = "camelscz_{gauge}.csv"
        for key, var, selection in [
            ("P", "pVars", precipitation),
            ("T", "tVars", temperature),
            ("Ep", "eVars", pet),
        ]:
            match = re.search(var + r"\s*=\s*\{(.*?)\}", text, re.S)
            choices(key, strings(match.group(1)), selection or 1)
    elif region == "CAMELS_ES":
        source("P", "total_precipitation_sum")
        source("T", "temperature_2m_mean", "degC")
        source("Ep", "potential_evaporation_sum")
    elif region == "CAMELS_IL":
        if pet == 4:
            source("Ep", "potential_evaporation_sum_FAO_PENMAN_MONTEITH")
        elif pet == 5:
            source("Ep", "potential_evaporation_sum_ERA5_LAND")
        else:
            variables["Ep"] = dict(derive="fao56_daily")
            source("Tmin", "temperature_2m_min")
            source("Tmax", "temperature_2m_max")
    elif region == "CAMELS_KR":
        profile.update(
            pattern="CAMELS_KR_Meteorological_timeseries_{gauge}.csv",
            flow_pattern=(
                "../discharge/CAMELS_KR_Hydrological_timeseries_{gauge}.csv"
            ),
            date="date",
            format="csv",
            layout="one_file_per_basin",
            file_folder="forcing",
        )
        source("P", "prec")
        source("T", "temp_avg", "degC")
        choices("Ep", ["pet", "pet_gleam"], pet or 1)
        source("Q", "discharge_spec")
    elif region == "CAMELS_ET":
        profile["date"] = "date"
        profile["date_format"] = "yyyy-MM-dd"
        variables["T"] = dict(derive="mean_tmin_tmax")
    elif region == "CAMELS_DE":
        profile["pattern"] = "CAMELS_DE_hydromet_timeseries_DE{gauge}.csv"
        source("P", "precipitation_mean")
        source("T", "temperature_mean", "degC")
        source("Q", "discharge_spec_obs")
        source("Tmin", "temperature_min")
        source("Tmax", "temperature_max")
        variables["Ep"] = dict(derive="fao56_daily")
    elif region in ("LAMAH_CE", "CAMELS_IS"):
        profile["pattern"] = "forcing/ID_{gauge}.csv"
        profile["flow_pattern"] = (
            "streamflow/ID_{gauge}.csv"
            if region == "LAMAH_CE"
            else "discharge/ID_{gauge}.csv"
        )
        source("Q", "qobs", "m3/s")
        variables["Ep"] = dict(derive="fao56_daily")
        if region == "CAMELS_IS" and pet in (4, 5):
            variables["Ep"] = dict(
                source="pet" if pet == 4 else ("total_et")
            )
        elif pet == 4:
            variables["Ep"] = dict(source="total_et", absolute=True)
        source("Tmin", "2m_temp_min")
        source("Tmax", "2m_temp_max")
        if temperature == 2:
            variables["T"] = dict(derive="mean_tmin_tmax")
    elif region == "HYD_RESPONSES":
        variables.pop("Ep", None)
        variables["Ep"] = dict(
            source="pev_mean_mean", absolute=True, clip_min=0
        )
        source("T", "TabsD_mean", "degC")
        source("P", "RhiresD_mean")
        if precipitation == 2:
            source("P", "tp_mean_mean")
        if temperature == 2:
            source("T", "t2m_mean_mean", "degC")
        elif temperature == 3:
            variables["T"] = dict(derive="mean_tmin_tmax")
            source("Tmin", "TminD_mean")
            source("Tmax", "TmaxD_mean")
    return profile


def reader_issues(sage_root, region, resolution=1):
    """List unsupported reader settings for the selected region."""
    if resolution != 1:
        from regional_subdaily import SUPPORTED

        return (
            []
            if region in SUPPORTED
            and resolution in (24, 96)
            and (resolution == 24 or region == "CAMELS_LUX")
            else ["Selected subdaily reader port is pending."]
        )
    if region == "CAMELS_US":
        return []
    try:
        profile = table_profile(sage_root, region)
        if "Daily" not in profile["meta"]["resolutions"]:
            return [f"{region} does not supply daily data."]
        from regional_special import SPECIAL

        if region in SPECIAL:
            return []
        if (
            profile["format"] != "csv"
            or profile["layout"] != "one_file_per_basin"
            or not profile["pattern"]
        ):
            return [
                f'{region}: specialized {profile["format"]} /'
                f' {profile["layout"]} reader is pending.'
            ]
        missing = set(("P", "T", "Ep", "Q")) - profile["variables"].keys()
        return (
            [
                f"{region}: unresolved reader variables:"
                f' {", ".join(sorted(missing))}.'
            ]
            if missing
            else []
        )
    except (ValueError, KeyError) as exc:
        return [str(exc)]


def _table(
    file,
    skiprows=0,
    columns=None,
    dtype=None,
    static=False,
    transient=False,
):
    """
    Read each file once; reuse unchanged shared attributes, never
    forcing arrays.
    """
    from local_data import require_local

    file = require_local(file)
    if not transient and (
        static
        or columns is not None
        or file.name == "gauge_information.txt"
    ):
        info = file.stat()
        return _static_table(
            str(file.resolve()),
            info.st_mtime_ns,
            info.st_size,
            skiprows,
            tuple(columns) if columns is not None else None,
            dtype,
        ).copy(deep=True)
    return _parse_table(
        file.read_bytes(),
        skiprows,
        columns,
        dtype,
        round_trip=file.name.endswith("_hourly_nldas.csv"),
    )


@lru_cache(maxsize=32)
def _static_table(path, mtime_ns, size, skiprows, columns, dtype):
    """Read a statically defined regional metadata table."""
    return _parse_table(Path(path).read_bytes(), skiprows, columns, dtype)


def _parse_table(
    raw, skiprows=0, columns=None, dtype=None, round_trip=False
):
    """Parse a regional metadata table."""
    from io import BytesIO, TextIOWrapper

    with TextIOWrapper(
        BytesIO(raw), encoding="utf-8-sig", errors="strict"
    ) as handle:
        for _ in range(skiprows):
            handle.readline()
        first = handle.readline()
    delimiter = (
        ";"
        if first.count(";") > first.count(",")
        else "," if "," in first else "\t" if "\t" in first else r"\s+"
    )
    options = dict(
        sep=delimiter,
        engine="c",
        encoding="utf-8-sig",
        skiprows=skiprows,
    )
    if columns is not None:
        options["usecols"] = lambda name: str(name).strip() in columns
    if dtype is not None:
        options["dtype"] = dtype
    options["low_memory"] = False
    options["float_precision"] = "round_trip"
    if round_trip:
        options["float_precision"] = "round_trip"
    frame = pd.read_csv(BytesIO(raw), **options)
    frame.columns = [str(c).strip() for c in frame.columns]
    return frame


def _dates(frame, profile):
    """Resolve regional reference dates."""
    if all(profile["ymd"]):
        cols = profile["ymd"]
        return pd.to_datetime(
            dict(
                year=frame[cols[0]],
                month=frame[cols[1]],
                day=frame[cols[2]],
            ),
            errors="raise",
        ).to_numpy(dtype="datetime64[D]")
    key = profile["date"]
    if not key or key not in frame:
        raise ValueError(f"Missing date column {key!r}")
    formats = {
        "yyyy-MM-dd": "%Y-%m-%d",
        "dd/MM/yyyy": "%d/%m/%Y",
        "yyyyMMdd": "%Y%m%d",
    }
    fmt = formats.get(profile["date_format"])
    return pd.to_datetime(frame[key], format=fmt, errors="raise").to_numpy(
        dtype="datetime64[D]"
    )


def _value(frame, spec, key):
    """Extract a regional metadata value."""
    if "source" in spec:
        if spec["source"] not in frame:
            raise ValueError(f'Missing {key} column: {spec["source"]}')
        values = pd.to_numeric(
            frame[spec["source"]], errors="coerce"
        ).to_numpy(dtype=float, copy=True)
        values[np.isin(values, [-999, -9999])] = np.nan
        return values
    sources = [s for s in spec.get("sources", []) if s in frame]
    if not sources:
        raise ValueError(f"Missing {key} columns: {spec.get((
                'sources'
            ),[])}")
    values = np.column_stack(
        [
            pd.to_numeric(frame[s], errors="coerce").to_numpy(float)
            for s in sources
        ]
    )
    values[np.isin(values, [-999, -9999])] = np.nan
    if spec.get("operation") == "mean":
        return pd.DataFrame(values).mean(axis=1).to_numpy(float)
    if spec.get("operation") == "coalesce":
        out = values[:, 0].copy()
        for col in values.T[1:]:
            out = np.where(np.isfinite(out), out, col)
        return out
    return values[:, 0]


def _file_root(root):
    """Resolve the root containing the regional source files."""
    root = Path(root)
    for candidate in (root / "daily" / "timeseries", root / "daily", root):
        if candidate.is_dir():
            return candidate
    return root


def catalogue(root, region, sage_root, **selections):
    """Return available basin identifiers and their display metadata.

    Resolve the selected resolution and regional directory layout before
    scanning files. No meteorological time series are loaded here.
    """
    from lamah_dataset import REGIONS, basin_ids, resolve_root

    if region in REGIONS:
        root = resolve_root(root, region, selections.get("resolution", 1))
        allowed = basin_ids(root, region)
        region = "LAMAH_CE"
    resolution = selections.get("resolution", 1)
    if resolution != 1:
        from regional_subdaily import inventory

        return [
            r
            for r in inventory(root, region, resolution)
            if region not in REGIONS or r[0] in allowed
        ]
    import xia_camels_us

    if xia_camels_us.project_region(root) == region:
        return xia_camels_us.catalogue(root)
    if region == "CAMELS_US":
        # Keep the existing forcing/discharge intersection semantics.
        product = selections.get("product", "nldas")
        base = Path(root) / "daily" / "v1p2"
        from camels_us import forcing_files

        force = {p.name[:8] for p in forcing_files(root, product)}
        flow = {
            p.name[:8]
            for p in (base / "streamflow").glob("*/*_streamflow_qc.txt")
        }
        return [(g, "", "") for g in sorted(force & flow)]
    profile = table_profile(
        sage_root,
        region,
        **{
            k: v
            for k, v in selections.items()
            if k in ("pet", "temperature", "precipitation")
        },
    )
    pattern = profile["pattern"]
    from regional_special import SPECIAL, inventory

    if region in SPECIAL:
        return inventory(root, region, profile)
    if (
        not pattern
        or profile["layout"] != "one_file_per_basin"
        or profile["format"] != "csv"
    ):
        raise ValueError(
            f'{region}: reader layout {profile["layout"]} needs a'
            " specialized adapter"
        )
    regex = (
        re.escape(pattern)
        .replace(re.escape("{gauge}"), r"(?P<gauge>.+?)")
        .replace(re.escape("{gauge02}"), r"(?P<gauge>.+?)")
        .replace(r"\*", ".*")
    )
    result = []
    folder = (
        Path(root) / "daily" / profile["file_folder"]
        if profile.get("file_folder")
        else _file_root(root)
    )
    for file in folder.glob(
        pattern.replace("{gauge}", "*").replace("{gauge02}", "*")
    ):
        match = re.fullmatch(
            regex, str(file.relative_to(folder)).replace("\\", "/")
        )
        if match:
            result.append((match.group("gauge"), "", ""))
    return sorted(
        {r for r in result if region not in REGIONS or r[0] in allowed}
    )


def _metadata(profile, root, gauge, forcing_file=None):
    """Assemble the basin metadata for a reader profile."""
    text = re.sub(
        r"\b(?:B|D|H)\.aux\.tables\.",
        "S.aux.tables.",
        profile.get("text", ""),
    )
    out = {}
    # MATLAB visits auxiliary tables in declaration order and keeps the
    # first finite value. A set made PET nondeterministic
    # between Python processes when both tables supplied elevation/area.
    tables = profile.get("auxiliary_tables")
    if tables is None:
        tables = []
        for name in dict.fromkeys(
            re.findall(r"S\.aux\.tables\.(\w+)\.file", text)
        ):
            prefix = "S.aux.tables." + name + "."
            table = {
                k: _literal(text, prefix + k)
                for k in (
                    "file",
                    "key",
                    "lat",
                    "elev",
                    "area",
                    "strip_prefix",
                )
            }
            table["alternates"] = re.findall(
                re.escape(prefix) + r"file\s*=\s*'([^']*)'", text
            )
            scale = re.search(
                re.escape(prefix + "area_scale") + r"\s*=\s*([0-9.eE+-]+)",
                text,
            )
            table["area_scale"] = float(scale[1]) / 1e6 if scale else 1.0
            tables.append(table)
    for table in tables:
        file = table.get("file")
        key = table.get("key")
        if not file or not key:
            continue
        candidates = [
            (_file_root(root) / file).resolve(),
            Path(root) / Path(file).name,
        ]
        for alternate in table.get("alternates", []):
            candidates.extend(
                [
                    Path(root) / Path(alternate).name,
                    (_file_root(root) / alternate).resolve(),
                ]
            )
        if forcing_file:
            candidates.insert(
                0, (Path(forcing_file).parent / file).resolve()
            )
        path = next((p for p in candidates if p.is_file()), candidates[0])
        if not path.is_file():
            continue
        wanted = [
            key,
            *[table.get(kind) for kind in ("lat", "elev", "area")],
        ]
        frame = _table(path, columns=wanted, dtype=str)
        if key not in frame:
            continue
        ids = (
            frame[key]
            .str.strip()
            .str.replace(r"\.0$", "", regex=True)
            .str.replace(r"^CAN_", "", regex=True)
            if profile["meta"]["acronym"] == "CA"
            else frame[key].str.strip().str.replace("\\.0$", "", regex=True)
        )
        if profile["meta"]["acronym"] == "CZ":
            ids = ids.str.replace(r"(?i)^camelscz_", "", regex=True)
        if profile["meta"]["acronym"] == "DE":
            ids = ids.str.replace(r"^DE", "", regex=True)
        strip = table.get("strip_prefix")
        if strip:
            ids = ids.str.replace("^" + re.escape(strip), "", regex=True)
        row = frame.loc[
            (ids == str(gauge))
            | (ids.str.lstrip("0") == str(gauge).lstrip("0"))
        ]
        if row.empty:
            continue
        for kind in ("lat", "elev", "area"):
            col = table.get(kind)
            if col and col in row:
                value = pd.to_numeric(row.iloc[0][col], errors="coerce")
                if np.isfinite(value) and kind not in out:
                    # Python discharge conversion uses km2; MATLAB aux uses
                    # m2.
                    out[kind] = float(value) * (
                        table.get("area_scale", 1.0)
                        if kind == "area"
                        else 1.0
                    )
    # Some schemas derive static latitude/elevation from time-series fields.
    fallback = Path(root) / "gauge_information.txt"
    if fallback.is_file() and (
        "lat" not in out or "area" not in out or "elev" not in out
    ):
        frame = _table(fallback)
        if "gauge_id" in frame:
            ids = (
                frame["gauge_id"]
                .astype(str)
                .str.replace(r"\.0$", "", regex=True)
            )
            row = frame.loc[
                (ids == str(gauge))
                | (ids.str.lstrip("0") == str(gauge).lstrip("0"))
            ]
            if not row.empty:
                for kind, col in [
                    ("lat", "gauge_lat"),
                    ("elev", "gauge_elev"),
                    ("area", "area_km2" if "area_km2" in row else "area"),
                ]:
                    if col in row and kind not in out:
                        value = pd.to_numeric(
                            row.iloc[0][col], errors="coerce"
                        )
                        if np.isfinite(value):
                            out[kind] = float(value)
    if profile["meta"]["acronym"] == "SE" and "lat" not in out:
        path = Path(root) / "catchments_physical_properties.csv"
        if path.exists():
            frame = _table(path, static=True)
            row = frame.loc[frame["ID"].astype(str) == str(gauge)]
            if not row.empty:
                for kind, col in [
                    ("lat", "Latitude_WGS84"),
                    ("area", "Area_km2"),
                    ("elev", "Elevation_mabsl"),
                ]:
                    out[kind] = float(row.iloc[0][col])
    if profile["meta"]["acronym"] == "PL" and "lat" not in out:
        path = Path(root) / "CAMELS_PL_topographic_attributes.csv"
        if path.exists():
            frame = _table(path, static=True)
            row = frame.loc[frame["gauge_id"].astype(str) == str(gauge)]
            if not row.empty:
                for kind, col in [
                    ("lat", "gauge_lon"),
                    ("area", "area_metadata"),
                    ("elev", "gauge_elev"),
                ]:
                    out[kind] = float(row.iloc[0][col])
    if profile["meta"]["acronym"] == "EU" and "lat" not in out:
        path = Path(root) / "CHAMELLION_topographic_attributes.csv"
        if path.exists():
            frame = _table(path, static=True)
            row = frame.loc[frame["gauge_id"].astype(str) == str(gauge)]
            if not row.empty:
                for kind, col in [
                    ("lat", "gauge_lat"),
                    ("area", "area"),
                    ("elev", "gauge_elev"),
                ]:
                    if col in row:
                        out[kind] = float(row.iloc[0][col])
    if profile["meta"]["acronym"] in ("AT", "IS") and "lat" not in out:
        path = Path(root) / "Gauge_attributes.csv"
        if path.exists():
            frame = _table(path, static=True)
            key = "ID" if "ID" in frame else "id"
            row = frame.loc[frame[key].astype(str) == str(gauge)]
            if not row.empty:
                from pyproj import Transformer

                record = row.iloc[0]
                longitude, latitude = Transformer.from_crs(
                    3035 if profile["meta"]["acronym"] == "AT" else 3057,
                    4326,
                    always_xy=True,
                ).transform(float(record["lon"]), float(record["lat"]))
                out["lat"] = latitude
                if "elev" not in out:
                    out["elev"] = float(
                        record["elev"]
                        if "elev" in record
                        else record["elevation"]
                    )
    if profile["meta"]["acronym"] == "KRH" and "area" not in out:
        path = Path(root) / "attributes_general.csv"
        if path.exists():
            frame = _table(path, static=True)
            row = frame.loc[frame["STAID"].astype(str) == str(gauge)]
            if not row.empty:
                out.update(
                    area=float(row.iloc[0]["Area"]),
                    lat=float(row.iloc[0]["Lat_gage"]),
                )
    if profile["meta"]["acronym"] == "BR":
        for name, columns in [
            (
                "camels_br_topography.txt",
                {"area": "area", "elev": "elev_mean"},
            ),
            ("camels_br_location.txt", {"lat": "gauge_lat"}),
        ]:
            # Shared metadata is independent of the basin time series.
            # Avoid reopening both full tables for every basin.
            if all(kind in out for kind in columns):
                continue
            path = Path(root) / name
            if not path.exists():
                continue
            frame = _table(path, columns=["gauge_id", *columns.values()])
            row = frame.loc[frame["gauge_id"].astype(str) == str(gauge)]
            if not row.empty:
                for kind, col in columns.items():
                    if kind not in out and col in row:
                        out[kind] = float(row.iloc[0][col])
    if profile["meta"]["acronym"] == "NZ" and (
        "lat" not in out or "area" not in out
    ):
        path = Path(root) / "1.CAMELS_NZ_Catchment_information.csv"
        if path.exists():
            frame = _table(path, static=True)
            row = frame.loc[frame["RID"].astype(str) == str(gauge)]
            if not row.empty:
                out.update(
                    lat=float(row.iloc[0]["Latitude (WGS 84)"]),
                    area=float(row.iloc[0]["uparea"]),
                    elev=float(row.iloc[0]["elevation"]),
                )
    return out


def extraterrestrial(dates, latitude):
    """Compute daily extraterrestrial radiation from dates and latitude."""
    doy = pd.DatetimeIndex(dates).dayofyear.to_numpy()
    phi = np.deg2rad(latitude)
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365)
    delta = 0.409 * np.sin(2 * np.pi * doy / 365 - 1.39)
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(delta), -1, 1))
    return (
        24
        * 60
        / np.pi
        * 0.0820
        * dr
        * (
            ws * np.sin(phi) * np.sin(delta)
            + np.cos(phi) * np.cos(delta) * np.sin(ws)
        )
    )


def read_region(
    root,
    gauge,
    *,
    region="CAMELS_US",
    sage_root=None,
    data_actions=None,
    **options,
):
    """Read canonical forcing and discharge for one basin.

    Return dates, forcing columns P/Ep/T, observed discharge and score
    masks. Track every source dependency for cache invalidation and apply
    registered data corrections once. Basin selection limits wide tables to
    the columns needed by the current run.
    """
    import xia_camels_us

    if xia_camels_us.project_region(root) == region:
        from timeseries_storage import store_timeseries

        data = xia_camels_us.read(root, gauge, **options)
        if not (options.get("allow_invalid") or data["train"].any()):
            raise ValueError("No valid Project 6 training discharge.")
        data["matlab_observed"] = data["observed"].copy()
        return store_timeseries(data)

    from lamah_dataset import REGIONS, basin_ids, resolve_root

    if region in REGIONS:
        root = resolve_root(root, region, options.get("resolution", 1))
        if str(gauge) not in basin_ids(root, region):
            raise ValueError(f"Gauge {gauge} does not belong to {region}")
        region = "LAMAH_CE"
    # Accept publisher-prefixed IDs as reader inputs as well as GUI IDs.
    if region == "CAMELS_CZ":
        gauge = re.sub(r"(?i)^camelscz_", "", str(gauge))
    if region == "CAMELS_DE":
        gauge = str(gauge).removeprefix("DE")
    from hydrologic_data_actions import load, apply, exclusion_reason

    actions = (
        load(region, options.get("resolution", 1), sage_root)
        if data_actions is None
        else data_actions
    )
    reason = exclusion_reason(actions, gauge)
    if reason and not options.get("allow_invalid", False):
        raise ValueError(f"Excluded gauge {gauge}: {reason}")
    from local_data import capture_sources
    from hydro_reader import read as read_hydro_timeseries

    with capture_sources(options.pop("basin_selection", ())) as sources:
        data = read_hydro_timeseries(
            root, gauge, region=region, sage_root=sage_root, **options
        )
    data["_source_files"] = sorted(sources)
    data = apply(data, actions, gauge)
    # Match MATLAB's final materialization, after curated data actions.
    for key in ("forcing", "observed", "_source_SWE", "_source_SM"):
        if key in data:
            data[key] = np.asarray(data[key], dtype=np.float32).astype(
                float
            )
    data["matlab_observed"] = data["observed"].copy()
    if (
        not options.get("allow_invalid", False)
        and not options.get("meteo_only", False)
        and not np.asarray(data["train"]).any()
    ):
        raise ValueError(
            "No valid training discharge after regional data corrections."
        )
    from timeseries_storage import store_timeseries

    return store_timeseries(data)


def _read_region(*args, **kwargs):
    """Compatibility entry point for the generic schema engine."""
    from hydro_reader import read

    return read(*args, **kwargs)


_legacy_region_registry = region_registry
_legacy_table_profile = table_profile


@lru_cache(maxsize=4)
def region_registry(sage_root):
    """Return installed region metadata keyed by canonical region code."""
    from runtime_metadata import regions

    result = regions(sage_root)
    return (
        result if result is not None else _legacy_region_registry(sage_root)
    )


def table_profile(sage_root, region, pet=1, temperature=1, precipitation=1):
    """Resolve the selected PET, temperature and precipitation schema.

    Prefer bundled runtime metadata; retain the legacy schema fallback.
    Reject unsupported forcing selections instead of substituting data.
    """
    from runtime_metadata import load
    from copy import deepcopy

    from lamah_dataset import reader_region

    region = reader_region(region)
    data = load(str(sage_root))
    if data is None:
        return _legacy_table_profile(
            sage_root, region, pet, temperature, precipitation
        )
    key = f"{region}:{int(pet)}:{int(temperature)}:{int(precipitation)}"
    if key not in data["table_profiles"]:
        raise ValueError(
            f"Invalid {region} forcing selection: {pet}, {temperature},"
            f" {precipitation}"
        )
    profile = deepcopy(data["table_profiles"][key])
    if "error" in profile:
        raise ValueError(profile["error"])
    profile["meta"] = region_registry(str(sage_root))[region]
    return profile
