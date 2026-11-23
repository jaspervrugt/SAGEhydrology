# ==========================================================================
# SAGEhydrology | Python source
# File: python/hydro_reader.py
#
# One schema-driven meteo/discharge reader, ported from MATLAB's engine.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""One schema-driven meteo/discharge reader, ported from MATLAB's engine."""

from pathlib import Path
import numpy as np
from hydro_schema import configuration, identifier
from hydro_io import read_raw
from hydro_aux import auxiliary
from hydro_derivations import derive, variable


def read(
    root,
    gauge,
    *,
    region="CAMELS_US",
    sage_root=None,
    product="nldas",
    pet=1,
    temperature=1,
    precipitation=1,
    resolution=1,
    allow_invalid=False,
    meteo_only=False,
    **period
):
    """Read basin forcing and discharge using MATLAB-compatible source
    rules.

    Parameters
    ----------
    root : path-like
        Regional dataset directory.
    gauge : str
        Basin identifier, normalized using the selected reader profile.
    region, sage_root, product : str or path-like
        Regional profile, SAGE metadata root and forcing product selection.
    pet, temperature, precipitation : int
        Source-column selections defined by the regional reader schema.
    resolution : int
        Samples per day: 1, 24, or 96.
    allow_invalid : bool
        Permit incomplete inputs for data-quality screening.
    meteo_only : bool
        Request forcing data without reading discharge.
    **period
        Training/evaluation start and end dates, plus spin-up days.

    Returns
    -------
    dict
        Aligned dates, forcing columns, discharge and source metadata.

    Notes
    -----
    Schema-specific clocks and missing-value rules follow the MATLAB reader.
    """
    meteo, discharge, paths, options = configuration(
        sage_root,
        region,
        resolution,
        product,
        pet,
        temperature,
        precipitation,
    )
    gauge = identifier(gauge, meteo["id"])
    step = int(meteo["timeline"]["step"])
    start = min(
        np.datetime64(period["train_start"], "D"),
        np.datetime64(period["eval_start"], "D"),
    )
    end = max(
        np.datetime64(period["train_end"], "D"),
        np.datetime64(period["eval_end"], "D"),
    )
    reference = np.datetime64(meteo["timeline"]["reference"], "s")
    # prepare_split offsets are relative to the schema's reference date;
    # preserve its clock hour (e.g. GB's 09:00 hydrological day).
    timeofday = reference - reference.astype("datetime64[D]")
    dates = np.arange(
        start - np.timedelta64(int(period["spinup"]), "D") + timeofday,
        # MATLAB builds complete day counts from the schema reference
        # clock. GB hydrological days run 09:00 through 08:00 next day.
        end + np.timedelta64(1, "D") + timeofday,
        np.timedelta64(step, "s"),
    ).astype("datetime64[s]")
    meteo_root = Path(root) / paths["meteo"].replace("\\", "/")
    aux = auxiliary(meteo_root, root, gauge, meteo)
    raw, used = read_raw(meteo_root, gauge, meteo, dates, options)
    raw = derive(raw, meteo, dates, aux)
    forcing = np.column_stack(
        [
            variable(raw[k], meteo["variables"][k], aux, step)
            for k in ("P", "Ep", "T")
        ]
    )
    # MATLAB's materialization boundary is single for every region.
    matlab_forcing = forcing.astype(np.float32).astype(float)
    if meteo_only:
        q = np.full(len(dates), np.nan)
    elif "Q" in raw and discharge.get("reuse_existing_q"):
        q = variable(raw["Q"], meteo["variables"]["Q"], aux, step)
    else:
        flow_root = Path(root) / paths["discharge"].replace("\\", "/")
        qraw, qused = read_raw(
            flow_root,
            identifier(gauge, discharge["id"]),
            discharge,
            dates,
            options,
        )
        q = variable(qraw["Q"], discharge["variables"]["Q"], aux, step)
        used.extend(qused)
    valid = np.isfinite(q) & (q >= 0)
    q[~valid] = np.nan
    day = dates.astype("datetime64[D]")
    train = (
        valid
        & (day >= np.datetime64(period["train_start"]))
        & (day <= np.datetime64(period["train_end"]))
    )
    evaluation = (
        valid
        & (day >= np.datetime64(period["eval_start"]))
        & (day <= np.datetime64(period["eval_end"]))
    )
    if meteo_only:
        train = (day >= np.datetime64(period["train_start"])) & (
            day <= np.datetime64(period["train_end"])
        )
        evaluation = (day >= np.datetime64(period["eval_start"])) & (
            day <= np.datetime64(period["eval_end"])
        )
    if not allow_invalid:
        if not np.all(np.isfinite(forcing)):
            raise ValueError("Missing or invalid meteorological forcing")
        if not train.any() and not meteo_only:
            raise ValueError("No valid training discharge")
    result = dict(
        gauge=gauge,
        dates=dates,
        forcing=forcing,
        observed=q,
        matlab_forcing=matlab_forcing,
        matlab_observed=q.astype(np.float32).astype(float),
        train=train,
        evaluation=evaluation,
        area=aux[2],
        forcing_file=str(used[0]),
        streamflow_file=str(used[-1]),
    )
    for field, key in [
        ("SWE", "_source_SWE"),
        ("SM", "_source_SM"),
        ("Radiation", "matlab_SW"),
        ("DayLength", "matlab_DayLength"),
    ]:
        if field in raw:
            value = variable(
                raw[field], meteo["variables"][field], aux, step
            )
            result[key] = value.astype(np.float32).astype(float)
    return result
