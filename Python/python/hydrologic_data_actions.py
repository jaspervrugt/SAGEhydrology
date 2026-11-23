# SITEhydrology | Python source
# ==========================================================================
# SAGEhydrology | Python source
# File: python/hydrologic_data_actions.py
#
# Curated SAGE data actions; source inputs are never modified.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Curated SAGE data actions; source inputs are never modified."""

import csv
from datetime import date
import hashlib
import json
from pathlib import Path
import re
import numpy as np

FIELDS = {
    "region",
    "resolution",
    "gauge_id",
    "variable",
    "action",
    "date",
    "reason",
    "evidence",
    "factor",
    "enabled",
}


def normalize_gauge(value):
    text = str(value).strip().upper()
    return (
        str(int(text.split(".")[0]))
        if re.fullmatch(r"\d+(?:\.0+)?", text)
        else text
    )


def load(region, resolution=1, sage_root=None):
    """
    Use the installed registry when present, otherwise the bundled
    copy.
    """
    file = (
        Path(sage_root)
        / "utils"
        / "metadata"
        / "hydrologic_data_actions.csv"
        if sage_root
        else Path(__file__).with_suffix(".csv")
    )
    if not file.is_file():
        file = Path(__file__).with_suffix(".csv")
    resolution = {1: "Daily", 24: "Hourly", 96: "15 minutes"}.get(
        resolution, str(resolution)
    )
    with file.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = FIELDS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "Data-action registry is missing columns: "
                + ", ".join(sorted(missing))
            )
        rows = list(reader)
    actions = []
    for row in rows:
        row = {k: v.strip() for k, v in row.items()}
        if (
            row["enabled"].lower() not in ("true", "1")
            or row["region"].casefold() != str(region).strip().casefold()
            or row["resolution"].casefold()
            not in ("*", resolution.casefold())
        ):
            continue
        row["gauge_id"] = normalize_gauge(row["gauge_id"])
        row["action"] = row["action"].lower()
        if row["action"] not in ("scale", "mask", "exclude"):
            raise ValueError(
                "Unsupported regional data action: " + row["action"]
            )
        row["factor"] = (
            float(row["factor"]) if row["action"] == "scale" else None
        )
        if row["action"] == "scale" and (
            not np.isfinite(row["factor"]) or row["factor"] <= 0
        ):
            raise ValueError(
                "Data-action scale factor must be finite and positive."
            )
        if row["action"] == "mask":
            date.fromisoformat(row["date"])
        row["registry_file"] = str(file.resolve())
        actions.append(row)
    return actions


def exclusion_reason(actions, gauge):
    return "; ".join(
        a["reason"]
        for a in actions
        if a["action"] == "exclude"
        and normalize_gauge(a["gauge_id"]) == normalize_gauge(gauge)
    )


def note(actions):
    lines = []
    for a in actions:
        if a["action"] == "exclude":
            text = f"Excluded gauge {a['gauge_id']}: {a['reason']}."
        elif a["action"] == "mask":
            text = (
                f"Masked {a['variable']} for gauge {a['gauge_id']} on"
                f" {a['date']}: {a['reason']}."
            )
        else:
            text = (
                f"Scaled {a['variable']} for gauge {a['gauge_id']} by"
                f" {a['factor']:.12g}: {a['reason']}."
            )
        lines.append(text)
    return "\n".join(lines)


def sync_notes(current, actions):
    current = re.sub(
        r"(?s)\s*Automatic data"
        r" actions\n.*?----------------------------------------\s*",
        "\n",
        str(current),
    ).strip()
    generated = note(actions)
    return (
        current
        + ("\n\n" if current and generated else "")
        + (
            "Automatic data actions\n"
            + generated
            + "\n----------------------------------------"
            if generated
            else ""
        )
    )


def configure(config):
    c = dict(config)
    actions = load(
        c.get("region", "CAMELS_US"),
        c.get("resolution", 1),
        c.get("sage_root"),
    )
    c.update(
        data_actions=actions,
        data_action_note=note(actions),
        notes=sync_notes(c.get("notes", ""), actions),
    )
    return c


def effective_actions(actions, gauges):
    selected = {normalize_gauge(g) for g in gauges}
    fields = (
        "region",
        "resolution",
        "gauge_id",
        "variable",
        "action",
        "date",
        "factor",
    )
    return [
        {
            key: (
                normalize_gauge(a[key]) if key == "gauge_id" else a.get(key)
            )
            for key in fields
        }
        for a in actions
        if normalize_gauge(a["gauge_id"]) in selected
    ]


def apply(data, actions, gauge):
    """
    Correct canonical arrays once; preserve inputs and repair score
    masks.
    """
    selected = [
        a
        for a in actions
        if a["action"] != "exclude"
        and normalize_gauge(a["gauge_id"]) == normalize_gauge(gauge)
    ]
    if not selected:
        return data
    signature = hashlib.sha256(
        json.dumps(
            effective_actions(selected, [gauge]), sort_keys=True
        ).encode()
    ).hexdigest()
    previous = data.get("_data_action_signature")
    if previous is not None:
        if str(previous) != signature:
            raise ValueError(
                "Data were already corrected with a different regional"
                " registry; reread the source."
            )
        return data
    result = dict(data)
    for a in selected:
        variable = a["variable"].upper()
        key = {
            "Q": "observed",
            "SWE": "observed_SWE",
            "SM": "observed_SM",
        }.get(variable)
        if key is None:
            raise ValueError(
                "Unsupported data-action variable: " + variable
            )
        if key not in result:
            continue
        values = np.array(result[key], float, copy=True)
        if a["action"] == "scale":
            values *= a["factor"]
        else:
            dates = np.asarray(result["dates"])
            target = np.datetime64(a["date"])
            matches = dates == target
            # The common readers expose the selected slice, whereas MATLAB
            # applies masks on its full timeline before slicing. Dates
            # outside
            # this slice cannot affect the returned series.
            if (
                matches.sum() != 1
                and len(dates)
                and dates.min() <= target <= dates.max()
            ):
                raise ValueError(
                    f"Mask date {a['date']} is unavailable or duplicated"
                    f" for gauge {gauge}."
                )
            values[matches] = np.nan
        result[key] = values
    result["_data_action_signature"] = np.array(signature)
    if "observed" in result:
        valid = np.isfinite(result["observed"]) & (result["observed"] >= 0)
        for key in ("train", "evaluation"):
            if key in result:
                result[key] = np.asarray(result[key], bool) & valid
    return result
