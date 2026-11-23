# ==========================================================================
# SAGEhydrology | Python source
# File: python/gauge_names.py
#
# Display-only port of MATLAB standardize_camels_gauge_names.m.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Display-only port of MATLAB standardize_camels_gauge_names.m.

National spelling conventions are deliberately different. Source
metadata and
gauge identifiers are never changed. US replacement tables are
bundled verbatim.
"""

import json, re
from functools import lru_cache
from pathlib import Path


@lru_cache(maxsize=1)
def rules():
    """Load the regional gauge-name formatting rules."""
    return json.loads(
        Path(__file__)
        .with_name("gauge_name_rules.json")
        .read_text(encoding="utf-8")
    )


def text(value):
    """Convert the supplied gauge-name value to text."""
    if value is None:
        return ""
    try:
        if value != value:
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def title_case(value, small=(), canadian=False, upper_only=False):
    """Apply the gauge-name title-case rules."""
    if upper_only and not (
        value == value.upper() and value != value.lower()
    ):
        return value

    def cap(w):
        """Capitalize a gauge-name token."""
        if canadian:
            if w.upper() in (
                "I",
                "II",
                "III",
                "IV",
                "V",
                "VI",
                "VII",
                "VIII",
                "IX",
                "X",
            ):
                return w.upper()
            if len(w) >= 3 and w.startswith("mc"):
                return "Mc" + w[2].upper() + w[3:]
            if len(w) >= 3 and w.startswith("o'"):
                return "O'" + w[2].upper() + w[3:]
        return w[:1].upper() + w[1:]

    return " ".join(
        (
            w
            if i and w in small
            else (
                "-".join(cap(p) for p in w.split("-"))
                if canadian
                else cap(w)
            )
        )
        for i, w in enumerate(value.lower().split())
    )


def us_name(value):
    """Format a United States gauge name."""
    value = re.sub(r"\s*@\s*", " at ", value)
    for pattern, replacement in rules()["us_replacements"]:
        pattern = (
            pattern.replace(r"\.?\>", r"\b\.?")
            .replace(r"\<", r"\b")
            .replace(r"\>", r"\b")
        )
        value = re.sub(pattern, replacement, value)
    for word in ("near", "at", "above", "below", "upstream", "downstream"):
        value = re.sub(r"(?i)\b" + word + r"\b", word, value)
    for short, long in [
        ("So", "South"),
        ("N", "North"),
        ("S", "South"),
        ("E", "East"),
        ("W", "West"),
        ("M", "Middle"),
    ]:
        value = re.sub(r"(?i)^" + short + r"\.?\s+", long + " ", value)
    value = re.sub(r"(?i)\bL\b\.?(?=\s+[A-Za-z])", "Little", value)
    water = (
        "\\b(River|Creek|Brook|Branch|Fork|Canal|Wash|Reservoir|Lake|Ba"
        "you|Arroyo)"
    )
    value = re.sub("(?i)" + water + r"\s+A\s+", r"\1 at ", value)
    value = re.sub("(?i)" + water + r"\s+Bel\s+", r"\1 below ", value)
    if re.search(r"(?i),?\s*(CA|FL)\.?\s*$", value):
        value = re.sub(
            r"(?i)(?<!Highway )(?<!-)\bA\b(?=\s+\S)", "at", value
        )
    code = ""
    for alias, postal in rules()["us_state_aliases"]:
        pattern = r"(?i)(?:\s*,\s*|\s+)" + re.escape(alias) + r"\.?\s*$"
        if re.search(pattern, value):
            value = re.sub(pattern, "", value)
            code = postal
            break
    if not code:
        match = re.search(r"(?i)(?:\s*,\s*|\s+)([A-Z]{2})\.?\s*$", value)
        postal = (
            "AL AK AZ AR CA CO CT DE FL GA HI ID IL IN IA KS KY LA ME MD MA"
            " MI MN MS MO MT NE NV NH NJ NM NY NC ND OH OK OR PA RI SC SD"
            " TN TX UT VT VA WA WV WI WY".split()
        )
        if match and match[1].upper() in postal:
            code = match[1].upper()
            value = value[: match.start()]
    value = re.sub(r"\s+", " ", value).strip()
    return value + ", " + code if code else value


def standardize(name, region="", gauge=""):
    """Standardize gauge names for the selected region."""
    value = text(name)
    gauge = text(gauge)
    reg = (
        str(region)
        .upper()
        .replace("CAMELSH_", "")
        .replace("CAMELS_", "")
        .replace("READ_ATTR_", "")
        .strip()
    )
    if reg in ("SPAT", "NAM"):
        return standardize(
            value, "CA" if gauge.upper().startswith("CAN_") else "US", gauge
        )
    if reg == "MACH_US":
        reg = "US"
    if reg in ("CHAMELLION", "EU") and gauge.startswith("DK_"):
        value = re.sub(r"^\s*\d+(?:\.\d+)?\s*,\s*", "", value)
    if reg == "COL" and gauge:
        value = re.sub(
            r"[_\s]*" + re.escape(gauge) + r"\s*0?\s*$", "", value
        )
    if reg == "ET":
        value = re.sub(r"Nr\.?", " near ", value)
        value = re.sub(r"(?<=[A-Za-z])At(?=[A-Z\s])", " at ", value)
        value = re.sub(r"([a-z])([A-Z])", r"\1 \2", value)
        value = re.sub(r"(?i)Finoteselam", "Finote Selam", value)
        value = re.sub(r"^L\.\s*", "Little ", value)
    value = value.replace("_", " ")
    if reg in ("CH", "CL", "COL", "DE", "DK", "ES", "GB", "IND", "SE"):
        value = value.replace("-", " ")
    value = re.sub(r"\s+", " ", value).strip()
    if reg in ("CHAMELLION", "EU"):
        value = re.sub(r"\s*,\s*", ", ", value).strip()
    if reg == "KR":
        value = re.sub(r"\s*\(", " (", value)
    if reg == "AU":
        value = re.sub(r"\s*\([A-Za-z0-9]+\)\s*$", "", value).strip()
    if not value:
        value = (
            (
                {
                    "AT": "AT_",
                    "COL": "COL_",
                    "DK": "DK_",
                    "ES": "ES_",
                    "ET": "ET_",
                    "PL": "PL_",
                }.get(reg, "Gauge ")
                + gauge
            )
            if gauge
            else "Gauge"
        )
    small = {
        "BR": "da das de do dos e",
        "CA": "at of the and in on by near above below under over to from",
        "CH": "au aux de des du d la le les l et en sur sous",
        "ET": "at near and",
        "GB": "at of the and in on by near",
        "US": "at of the and in on by near above below to from",
    }
    if reg in ("CL", "COL", "ES"):
        value = title_case(value, "a al de del el en la las los y".split())
    elif reg in (
        "BR",
        "CA",
        "CH",
        "DE",
        "ET",
        "GB",
        "IND",
        "SE",
        "US",
        "CHAMELLION",
        "EU",
    ):
        value = title_case(
            value,
            small.get(reg, "").split(),
            canadian=reg == "CA",
            upper_only=reg in ("US", "CHAMELLION", "EU"),
        )
    if reg == "US":
        value = us_name(value)
    if reg == "FI":
        value = re.sub(r"\s*-\s*", " - ", value)
        value = re.sub(r"\s*,\s*", ", ", value)
        value = re.sub(r"\s+", " ", value).strip()
        for word in (
            "luusua",
            "virtuaali",
            "kokonaisvirtaama",
            "pato",
            "ala",
        ):
            for before in (r"(?<= - )", r"(?<=, )"):
                value = re.sub(
                    before + word + r"(?=$|[ ,;/])", word, value, flags=re.I
                )
        for word in (
            "Kongas",
            "Roukkajankoski",
            "Jolma",
            "Lohikoski",
            "Hulttilanjoki",
            "Vaajakoski",
            "Kurittukoski",
            "Hanhikoski",
        ):
            for before in ("^", r"(?<= - )"):
                value = re.sub(
                    before + word.lower() + r"(?=$|[ ,;/])",
                    word,
                    value,
                    flags=re.I,
                )
    return value
