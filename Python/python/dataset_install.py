# SITEhydrology | Python source
# File: python/dataset_install.py
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
"""dataset install support for SITE Python."""

from http.client import IncompleteRead
from urllib.error import URLError, HTTPError
from time import sleep

# SITEhydrology | Python source
# ==========================================================================
# SAGEhydrology | Python source
# File: python/dataset_install.py
#
# Cancellable regional archive downloads and non-overwriting installations.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Cancellable regional archive downloads and non-overwriting
installations.
"""

import hashlib
import socket
from http.client import IncompleteRead
from urllib.error import URLError, HTTPError
from time import sleep
import json
from pathlib import Path
import re
import shutil
import sys
import tarfile
from urllib.parse import urlparse, unquote
from urllib.request import Request, urlopen
import uuid
import zipfile
from functools import lru_cache
from regional_data import (
    region_registry,
    table_profile,
    clean_source,
    strings,
)


@lru_cache(maxsize=64)
def _country_basins(sage_root, region):
    from regional_data import GRDC

    meta = region_registry(str(Path(sage_root).resolve()))[region]
    if meta["acronym"] not in GRDC:
        return None
    files = list(
        Path(meta["source"]).parent.glob(meta["acronym"] + "_*_basins.txt")
    )
    if not files:
        raise ValueError(
            "Country basin list is required to subset the GRDC archive"
        )
    ids = set()
    for file in files:
        ids.update(re.findall(r"(?m)^\s*(\d+)\s*$", file.read_text()))
    return ids


class InstallCancelled(Exception):
    pass


def reference_urls(site_root, region):
    root = Path(site_root)
    file = next(
        (
            p
            for p in (
                root / "gui" / "SITE_ui.m",
                root / "GUI" / "SITE_ui.m",
                root.parent / "Matlab" / "gui" / "SITE_ui.m",
                root / "gui" / "SAGE_ui.m",
                root.parent / "Matlab" / "gui" / "SAGE_ui.m",
            )
            if p.is_file()
        ),
        None,
    )
    if file is None:
        return []
    text = file.read_text(encoding="utf-8-sig")
    begin = text.find("function [txt,urls] = localDatasetReferenceNotes")
    if begin < 0:
        return []
    section = text[begin:]
    section = section[
        : section.find("function urls = localDatasetCitationURLs")
    ]
    for block in re.split(r"\bcase\s+", section)[1:]:
        labels = block.split("\n", 1)[0]
        if region not in strings(labels):
            continue
        dois = re.findall(r"doi:([^\s\'\]\}]+)", block)
        return list(
            dict.fromkeys(
                "https://doi.org/" + doi.rstrip(".,") for doi in dois
            )
        )
    return []


def download_options(sage_root, region, resolution=1):
    """
    Read the existing installer URLs and checksums without running
    MATLAB.
    """
    from lamah_dataset import reader_region

    region = reader_region(region)
    if region == "LAMAH_CE":
        if resolution not in (1, 24):
            return []
        name, checksum = (
            ("2_LamaH-CE_daily.tar.gz", "69fd2733e969513403f923ecc5eaa3dc")
            if resolution == 1
            else (
                "1_LamaH-CE_daily_hourly.tar.gz",
                "121b2292b288ef1c2ef4d96250be77ac",
            )
        )
        return [
            dict(
                name=name,
                md5=checksum,
                url="https://zenodo.org/records/5153305/files/"
                + name
                + "?download=1",
            )
        ]
    text = clean_source(Path(sage_root) / "utils" / "data_helpers.m")
    tag = region_registry(str(Path(sage_root).resolve()))[region][
        "acronym"
    ].lower()
    tag = {
        "hyd": "hyd_responses",
        "nam": "ca",
        "ru": "ru",
        "za": "grdc_africa",
        "na": "grdc_africa",
        "ar": "grdc_africa",
        "be": "grdc_africa",
        "ee": "grdc_africa",
        "ie": "grdc_africa",
        "jm": "grdc_africa",
        "no": "grdc_africa",
        "pr": "grdc_africa",
        "de": "de_daily",
        "deh": "de_hourly",
        "ush": "ush",
        "kr": "kr",
    }.get(tag, tag)
    if region == "HYD_RESPONSES":
        tag = "hyd_responses"
    if region == "CAMELS_DEH":
        tag = "de_hourly"
    if region == "CAMELSH_KR":
        tag = "krh"
    match = re.search(r"function out = local_download_" + tag + r"\(", text)
    block = (
        re.split(r"\nfunction ", text[match.end() :])[0] if match else ""
    )
    result = []
    helper = {
        "CAMELS_AU": "au",
        "CAMELS_BR": "br",
        "CAMELS_CL": "cl",
        "CAMELS_FR": "fr",
        "CAMELS_COL": "col",
        "CAMELS_NZ": "nz",
        "CAMELSH_KR": "kr",
    }.get(region)
    if helper:
        match = re.search(
            r"function files = local_camels_" + helper + r"_files\(", text
        )
        if match:
            section = re.split(r"\nfunction ", text[match.end() :])[0]
            base_match = re.search(
                r"(?:base|root)\s*=\s*(\[.*?\]|'[^']*');", section, re.S
            )
            base = (
                "".join(strings(base_match.group(1))) if base_match else ""
            )
            for call in re.finditer(
                r"files\(end\+1\)\s*=\s*local_\w+_file\((.*?)\);",
                section,
                re.S,
            ):
                args = strings(call.group(1))
                names = [
                    a
                    for a in args
                    if a.endswith((".zip", ".csv", ".txt", ".tar.gz"))
                ]
                if not names or not base:
                    continue
                name = names[0]
                checksum = next(
                    (a for a in args if re.fullmatch("[a-fA-F0-9]{32}", a)),
                    None,
                )
                from urllib.parse import quote

                url = (
                    base
                    + quote(name)
                    + (
                        "/content"
                        if "/api/records/" in base
                        else "?download=1" if "zenodo.org" in base else ""
                    )
                )
                result.append(dict(name=name, url=url, md5=checksum))
            if helper == "fr":
                for index in (1, 2):
                    key = "files(" + str(index) + ")."
                    name = re.search(
                        re.escape(key) + r"fileName\s*=\s*'([^']*)'",
                        section,
                    )
                    url = re.search(
                        re.escape(key) + r"url\s*=\s*(\[.*?\]|\'.*?\');",
                        section,
                        re.S,
                    )
                    md5 = re.search(
                        re.escape(key) + r"md5\s*=\s*'([^']*)'", section
                    )
                    if name and url:
                        result.append(
                            dict(
                                name=name.group(1),
                                url="".join(strings(url.group(1))),
                                md5=md5.group(1) if md5 else None,
                            )
                        )
            if helper in ("nz", "col"):
                for call in re.finditer(
                    r"files\(end\+1\)\s*=\s*struct\((.*?)\);", section, re.S
                ):
                    chunk = call.group(1)

                    def prop(key):
                        found = re.search(
                            "'" + key + r"'\s*,\s*(\[.*?\]|'[^']*')",
                            chunk,
                            re.S,
                        )
                        return (
                            "".join(strings(found.group(1)))
                            if found
                            else ""
                        )

                    name = prop("fileName")
                    url = prop("url")
                    if not name or not url:
                        continue
                    if (
                        "hourly" in name.lower()
                        and resolution == 1
                        or "daily" in name.lower()
                        and resolution != 1
                    ):
                        continue
                    result.append(
                        dict(name=name, url=url, md5=prop("md5") or None)
                    )
    for match in re.finditer(
        r"(?<!\w)url\s*=\s*(\[.*?\]|\'.*?\');", block, re.S
    ):
        expression = match.group(1)
        # Only complete literal URLs; dynamic URLs require their own recipe.
        if re.sub(r"'[^']*'|[\[\]\s]", "", expression):
            continue
        url = "".join(strings(expression))
        name = unquote(urlparse(url).path.split("/")[-1])
        if not name or not name.lower().endswith(
            (".zip", ".tar.gz", ".tgz", ".csv", ".txt")
        ):
            continue
        checksum = re.search(
            r"expectedMD5\s*=\s*'([a-fA-F0-9]{32})'", block
        )
        item = dict(
            url=url, name=name, md5=checksum.group(1) if checksum else None
        )
        if item not in result:
            result.append(item)
    if region == "LAMAH_CE":
        result = [
            item
            for item in result
            if ("daily_hourly" in item["name"]) == (resolution != 1)
        ]
    explicit = {
        "CAMELS_DK": [
            dict(
                name="dataverse_files.zip",
                url=(
                    "https://dataverse.geus.dk/api/access/dataset/"
                    ":persistentId/?persistentId=doi%3A10.22008%2FFK2%2"
                    "FAZXSYP"
                ),
                md5=None,
            )
        ],
        "BULL_ES": [
            dict(
                name=name,
                url="https://zenodo.org/records/10605646/files/"
                + name
                + "?download=1",
                md5=checksum,
            )
            for name, checksum in [
                ("attributes.zip", "2e0321e828217b52560fd8886ef9e5f7"),
                ("timeseries.zip", "3194058494c6f8170433642a4dddf1cf"),
            ]
        ],
        "CAMELS_IS": [
            dict(
                name="lamah_ice"
                + ("_hourly" if resolution == 24 else "")
                + ".zip",
                url=(
                    "https://www.hydroshare.org/resource/"
                    "705d69c0f77c48538d83cf383f8c63d6/data/contents/"
                    "lamah_ice"
                )
                + ("_hourly" if resolution == 24 else "")
                + ".zip",
                md5=None,
            )
        ],
        "CAMELSH_US": [
            dict(
                name=name,
                url="https://zenodo.org/records/15066778/files/"
                + name
                + "?download=1",
                md5=None,
            )
            for name in ("attributes.7z", "timeseries.7z", "info.csv")
        ],
        "CAMELS_CZ": [
            dict(
                name="CAMELS-CZ.zip",
                url=(
                    "https://zenodo.org/records/17769325/files/"
                    "CAMELS-CZ.zip?download=1"
                ),
                md5="f614413ba5c58630ce5cae4ed1fc7e3c",
            )
        ],
        "MACH_US": [
            dict(
                name=name,
                url="https://zenodo.org/records/18686475/files/"
                + name
                + "?download=1",
                md5=checksum,
            )
            for name, checksum in [
                ("MACH_ts.zip", "ad9f88795b383278af00da84ec4453a1"),
                ("attributes.zip", "f5a9c8d12c153690a25b6f210a47a0e2"),
            ]
        ],
        "CHAMELLION": [
            dict(
                name=name,
                url=(
                    "https://datapub.gfz.de/download/"
                    "10.5880.FIDGEO.2026.073-Wasfdsg/"
                    "2026-073_Maharjan-et-al_CHAMELLION-data/"
                )
                + name,
                md5=None,
            )
            for name in [
                "timeseries.zip",
                "CHAMELLION_topographic_attributes.csv",
                "CHAMELLION_climatic_attributes.csv",
                "CHAMELLION_hydrologic_attributes.csv",
                "CHAMELLION_landcover_attributes.csv",
                "CHAMELLION_soil_attributes.csv",
            ]
        ],
        "CAMELS_KR": [
            dict(
                name="CAMELS-KR.zip",
                url=(
                    "https://zenodo.org/api/records/21930882/files/"
                    "CAMELS-KR.zip/content"
                ),
                md5="5ffce545fc0b3b79d57382a33dc597ca",
            )
        ],
        "HYDRO_CIS": [
            dict(
                name="Russia_HydroMeteo_Database_v04.zip",
                url=(
                    "https://zenodo.org/api/records/8432070/files/"
                    "Russia_HydroMeteo_Database_v04.zip/content"
                ),
                md5=None,
            )
        ],
        "CAMELS_MX": [
            dict(
                name="Caravan-csv.tar.gz",
                url=(
                    "https://zenodo.org/api/records/15530022/files/"
                    "Caravan-csv.tar.gz/content"
                ),
                md5=None,
            )
        ],
        "CAMELS_SE": [
            dict(
                name=name + ".zip",
                url=(
                    "https://api.researchdata.se/dataset/2023-173/1/"
                    "file/data?filePath="
                )
                + name.replace(" ", "+")
                + ".zip",
                md5=None,
            )
            for name in ("catchment properties", "catchment time series")
        ],
    }
    if region == "CAMELS_US":
        items = [
            dict(
                name=name,
                url="https://zenodo.org/records/15529996/files/"
                + name
                + "?download=1",
                md5=None,
            )
            for name in (
                "camels_clim.txt",
                "camels_geol.txt",
                "camels_hydro.txt",
                "camels_name.txt",
                "camels_soil.txt",
                "camels_topo.txt",
                "camels_vege.txt",
            )
        ]
        if resolution == 1:
            items.append(
                dict(
                    name="basin_timeseries_v1p2_metForcing_obsFlow.zip",
                    url=(
                        "https://zenodo.org/records/15529996/files/"
                        "basin_timeseries_v1p2_metForcing_obsFlow.zip?d"
                        "ownload=1"
                    ),
                    md5="8e9a466710e8270b58f01d332a87184f",
                )
            )
        else:
            items.extend(
                dict(
                    name=name,
                    url="https://zenodo.org/records/4072701/files/"
                    + name
                    + "?download=1",
                    md5=None,
                )
                for name in (
                    "nldas_hourly_csv.tar.gz",
                    "usgs_streamflow_csv.tar.gz",
                )
            )
        return items
    if region in (
        "CAMELS_CA",
        "CAMELS_SPAT",
        "CAMELS_GB",
        "CAMELS_ET",
        "EStreams_PT",
        "CAMELS_IS",
    ):
        return [
            dict(
                name="Acquire " + region + " provider data",
                url="",
                provider=region,
            )
        ]
    if region == "CAMELS_CZ" and resolution == 24:
        return explicit[region] + [
            dict(
                name="CAMELS-CZ-Hourly.zip",
                url=(
                    "https://zenodo.org/records/17769325/files/"
                    "CAMELS-CZ-Hourly.zip?download=1"
                ),
                md5="182308706b70a616045ca12ecfd6466f",
            )
        ]
    return explicit.get(region, result)


def _check(stop):
    if stop is not None and stop.is_set():
        raise InstallCancelled("Installation canceled")


def installed_attribute(item, directory):
    (
        'Reuse complete local CAMELS attribute tables, not '
        'daily/hourly archives.'
    )
    name = item["name"]
    if not re.fullmatch(
        r"camels_(clim|geol|hydro|name|soil|topo|vege)\.txt", name
    ):
        return None
    candidates = [Path(directory) / name]
    if name == "camels_name.txt":
        candidates.append(Path(directory) / "camels_name_clean.txt")
    for path in candidates:
        if not path.is_file():
            continue
        try:
            import csv

            with path.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle, delimiter=";")
                if "gauge_id" not in (reader.fieldnames or []):
                    continue
                rows = list(reader)
            ids = [row.get("gauge_id", "") for row in rows]
            if len(ids) >= 500 and len(set(ids)) == len(ids) and all(ids):
                return path
        except (OSError, UnicodeError, csv.Error):
            continue
    return None


class IncompleteDownloadError(ValueError):
    """A response ended before its declared length; safe to retry."""


def download(
    item,
    directory,
    progress=lambda *args: None,
    stop=None,
    *,
    max_bytes=None,
):
    directory = Path(directory)
    name = Path(item["name"]).name
    if name != item["name"] or name in ("", ".", ".."):
        raise ValueError("Invalid download filename")
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / (name + ".lock")).open("a+b") as lock:
        if lock.seek(0, 2) == 0:
            lock.write(b"0")
            lock.flush()
        lock.seek(0)
        try:
            if __import__("os").name == "nt":
                import msvcrt

                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl

                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError(
                "Another installer is downloading " + name
            ) from exc
        try:
            return _download(
                item, directory, progress, stop, max_bytes=max_bytes
            )
        finally:
            lock.seek(0)
            if __import__("os").name == "nt":
                msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(lock, fcntl.LOCK_UN)


def _archive_identity(path):
    stat = Path(path).stat()
    return dict(
        size=stat.st_size,
        mtime_ns=stat.st_mtime_ns,
        ctime_ns=stat.st_ctime_ns,
    )


def _verified_archive(path, item=None):
    (
        'Reuse a verification receipt only while the local file is '
        'unchanged.'
    )
    path = Path(path)
    try:
        record = json.loads(
            path.with_name(path.name + ".download.json").read_text()
        )
        if record.get("identity") != _archive_identity(path):
            return None
        if not re.fullmatch(r"[0-9a-f]{64}", record.get("sha256", "")):
            return None
        if item is not None and (
            record.get("url") != item["url"]
            or (
                item.get("sha256")
                and record["sha256"] != item["sha256"].lower()
            )
            or (
                item.get("md5") and record.get("md5") != item["md5"].lower()
            )
        ):
            return None
        return record["sha256"]
    except (OSError, ValueError, TypeError):
        return None


def _replace_with_retry(source, destination):
    (
        'Tolerate brief Windows sharing locks while committing cache '
        'metadata.'
    )
    for attempt in range(6):
        try:
            return Path(source).replace(destination)
        except PermissionError:
            if attempt == 5:
                raise
            sleep(min(1.6, 0.1 * 2**attempt))


def _save_verified_archive(path, item, checksum):
    path = Path(path)
    record = dict(
        url=item["url"],
        sha256=checksum,
        md5=(item.get("md5") or "").lower(),
        size=path.stat().st_size,
        identity=_archive_identity(path),
    )
    receipt = path.with_name(path.name + ".download.json")
    staging = receipt.with_name(
        receipt.name + "." + uuid.uuid4().hex + ".tmp"
    )
    try:
        staging.write_text(json.dumps(record))
        _replace_with_retry(staging, receipt)
    finally:
        if staging.exists():
            staging.unlink()


def _download(
    item,
    directory,
    progress=lambda *args: None,
    stop=None,
    *,
    max_bytes=None,
):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    name = Path(item["name"]).name
    if name != item["name"] or name in ("", ".", ".."):
        raise ValueError("Invalid download filename")
    target = directory / name
    receipt = directory / (name + ".download.json")
    _check(stop)
    if target.is_file():
        progress("Checking completed download: " + name, 0, 0)
        checksum = _verified_archive(target, item)
        valid = checksum is not None
        if not valid:
            original_identity = _archive_identity(target)
            md5, sha = hashlib.md5(), hashlib.sha256()
            done = 0
            total = original_identity["size"]
            with target.open("rb") as handle:
                while True:
                    _check(stop)
                    chunk = handle.read(1024 * 1024)
                    if not chunk:
                        break
                    md5.update(chunk)
                    sha.update(chunk)
                    done += len(chunk)
                    progress(
                        "Checking completed download: " + name,
                        min(done, max(0, total - 1)),
                        total,
                    )
            checksum = sha.hexdigest()
            if item.get("sha256"):
                valid = checksum == item["sha256"].lower()
            elif item.get("md5"):
                valid = md5.hexdigest() == item["md5"].lower()
            elif receipt.is_file():
                try:
                    record = json.loads(receipt.read_text())
                    valid = (
                        record["url"] == item["url"]
                        and record["size"] == total
                        and record["sha256"] == checksum
                    )
                except (OSError, ValueError, KeyError):
                    pass
            valid = valid and original_identity == _archive_identity(target)
            if valid:
                _save_verified_archive(target, item, checksum)
        if valid and (
            max_bytes is None or target.stat().st_size <= max_bytes
        ):
            progress("Reusing completed download: " + name, 1, 1)
            return target
    temporary = directory / (name + ".part")
    resume = directory / (name + ".resume.json")

    def save_state(value):
        staging = directory / (
            resume.name + "." + uuid.uuid4().hex + ".tmp"
        )
        try:
            staging.write_text(json.dumps(value))
            _replace_with_retry(staging, resume)
        finally:
            if staging.exists():
                staging.unlink()

    def reset_partial():
        for path in (temporary, resume):
            if path.exists():
                path.unlink()

    max_attempts = 8
    for attempt in range(1, max_attempts + 1):
        _check(stop)
        try:
            state = (
                json.loads(resume.read_text()) if resume.exists() else {}
            )
        except (OSError, ValueError):
            state = {}
        offset = temporary.stat().st_size if temporary.is_file() else 0
        if offset and (
            state.get("url") != item["url"]
            or not (
                state.get("validator")
                or (item.get("md5") and state.get("md5") == item["md5"])
            )
        ):
            progress(
                "Saved transfer cannot be verified; restarting " + name,
                0,
                0,
            )
            reset_partial()
            state, offset = {}, 0
        try:
            digest, sha256 = hashlib.md5(), hashlib.sha256()
            if offset:
                progress("Checking saved download: " + name, 0, 0)
                with temporary.open("rb") as previous:
                    for chunk in iter(
                        lambda: previous.read(1024 * 1024), b""
                    ):
                        _check(stop)
                        digest.update(chunk)
                        sha256.update(chunk)
            if offset and state.get("total") == offset:
                done = total = offset
            else:
                headers = {
                    "User-Agent": "SITE-Python/1.0",
                    "Accept-Encoding": "identity",
                }
                if offset:
                    headers["Range"] = f"bytes={offset}-"
                    if state.get("validator"):
                        headers["If-Range"] = state["validator"]
                with urlopen(
                    Request(item["url"], headers=headers), timeout=120
                ) as response:
                    status = getattr(response, "status", None) or 200
                    if (
                        response.headers.get(
                            "x-amzn-waf-action", ""
                        ).lower()
                        == "challenge"
                    ):
                        raise ValueError(
                            "Provider requires browser verification for "
                            + name
                            + "; download the archive through the provider"
                            " website and use local archive installation"
                        )
                    length = int(response.headers.get("Content-Length", 0))
                    validator = response.headers.get("ETag", "")
                    if not validator or validator.startswith("W/"):
                        validator = response.headers.get(
                            "Last-Modified", ""
                        )
                    if status == 206:
                        match = re.fullmatch(
                            r"bytes (\d+)-(\d+)/(\d+)",
                            response.headers.get("Content-Range", ""),
                        )
                        if (
                            not match
                            or int(match[1]) != offset
                            or int(match[2]) < offset
                            or int(match[3]) <= int(match[2])
                            or (
                                length
                                and length != int(match[2]) - offset + 1
                            )
                            or (
                                state.get("total")
                                and state["total"] != int(match[3])
                            )
                            or (
                                offset
                                and validator
                                and validator != state.get("validator")
                            )
                        ):
                            reset_partial()
                            raise IncompleteDownloadError(
                                f"Invalid resume response for {name};"
                                " restarting safely"
                            )
                        total = int(match[3])
                    elif status == 200:
                        if offset:
                            progress(
                                "Server cannot resume this transfer;"
                                " restarting "
                                + name,
                                0,
                                0,
                            )
                        offset = 0
                        digest, sha256 = hashlib.md5(), hashlib.sha256()
                        total = length
                    else:
                        raise ValueError(
                            f"Unexpected HTTP status {status} for {name}"
                        )
                    if max_bytes is not None and (
                        total > max_bytes or offset > max_bytes
                    ):
                        raise ValueError(
                            "Download exceeds configured byte limit"
                        )
                    save_state(
                        dict(
                            url=item["url"],
                            validator=validator
                            or state.get("validator", ""),
                            total=total,
                            md5=item.get("md5"),
                        )
                    )
                    done = offset
                    if offset:
                        progress(
                            "Resuming " + name,
                            min(done, total - 1) if total else done,
                            total,
                        )
                    with temporary.open("ab" if offset else "wb") as handle:
                        while True:
                            _check(stop)
                            size = 1024 * 1024
                            if max_bytes is not None:
                                size = min(size, max_bytes - done + 1)
                            chunk = response.read(size)
                            if not chunk:
                                break
                            if (
                                max_bytes is not None
                                and done + len(chunk) > max_bytes
                            ):
                                raise ValueError(
                                    "Download exceeds configured byte limit"
                                )
                            handle.write(chunk)
                            digest.update(chunk)
                            sha256.update(chunk)
                            done += len(chunk)
                            progress(
                                name,
                                min(done, total - 1) if total else done,
                                total,
                            )
            if total and done != total:
                raise IncompleteDownloadError(
                    f"Incomplete download for {name}: {done}/{total} bytes"
                )
            if (
                item.get("md5")
                and digest.hexdigest().lower() != item["md5"].lower()
            ) or (
                item.get("sha256")
                and sha256.hexdigest() != item["sha256"].lower()
            ):
                reset_partial()
                raise ValueError(f"Checksum mismatch for {name}")
            if target.exists():
                with target.open("rb") as handle:
                    if (
                        hashlib.file_digest(handle, "sha256").hexdigest()
                        != sha256.hexdigest()
                    ):
                        raise FileExistsError(target)
                temporary.unlink()
                receipt.write_text(
                    json.dumps(
                        dict(
                            url=item["url"],
                            size=done,
                            sha256=sha256.hexdigest(),
                        )
                    )
                )
            else:
                # Write the receipt first so a shutdown immediately after
                # publishing the file still leaves a verifiable cache entry.
                receipt.write_text(
                    json.dumps(
                        dict(
                            url=item["url"],
                            size=done,
                            sha256=sha256.hexdigest(),
                        )
                    )
                )
                temporary.rename(target)
            _save_verified_archive(target, item, sha256.hexdigest())
            if resume.exists():
                resume.unlink()
            progress(name, max(1, done), max(1, done))
            return target
        except (
            TimeoutError,
            ConnectionError,
            URLError,
            IncompleteRead,
            IncompleteDownloadError,
        ) as exc:
            if isinstance(exc, HTTPError) and exc.code == 416:
                reset_partial()
            transient = not isinstance(exc, HTTPError) or exc.code in (
                408,
                416,
                429,
                500,
                502,
                503,
                504,
            )
            if not transient or attempt == max_attempts:
                if isinstance(exc, IncompleteDownloadError):
                    raise IncompleteDownloadError(
                        f"{exc} after {attempt} attempt(s); saved partial"
                        " retained"
                    ) from exc
                raise RuntimeError(
                    f"Download failed for {name} after {attempt}"
                    f" attempt(s): {exc}; saved partial retained"
                ) from exc
            delay = min(120, 10 * 2 ** (attempt - 1))
            if isinstance(exc, HTTPError) and exc.headers:
                retry_after = exc.headers.get("Retry-After", "")
                try:
                    from email.utils import parsedate_to_datetime
                    from datetime import datetime, timezone

                    requested = (
                        float(retry_after)
                        if retry_after.strip().isdigit()
                        else (
                            parsedate_to_datetime(retry_after)
                            - datetime.now(timezone.utc)
                        ).total_seconds()
                    )
                    delay = min(300, max(delay, requested))
                except (TypeError, ValueError, OverflowError):
                    pass
            progress(
                f"Retrying {name}: attempt {attempt + 1}/{max_attempts} in"
                f" {delay:g} s after {exc}",
                0,
                0,
            )
            if stop is not None:
                if stop.wait(delay):
                    raise InstallCancelled(
                        "Installation canceled; saved download retained"
                    )
            else:
                sleep(delay)


def extract(archive, destination, progress=lambda *args: None, stop=None):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)

    def target(name):
        if "\\" in name or ":" in name:
            raise ValueError(f"Unsafe archive member: {name}")
        file = (destination / name).resolve()
        if not file.is_relative_to(destination):
            raise ValueError(f"Archive member escapes destination: {name}")
        return file

    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as handle:
            members = handle.infolist()
            for index, member in enumerate(members):
                _check(stop)
                file = target(member.filename)
                if (member.external_attr >> 16) & 0o170000 == 0o120000:
                    raise ValueError("Archive symlinks are unsupported")
                if member.is_dir():
                    file.mkdir(parents=True, exist_ok=True)
                    continue
                file.parent.mkdir(parents=True, exist_ok=True)
                with handle.open(member) as source, file.open("xb") as dest:
                    shutil.copyfileobj(source, dest)
                progress(
                    "Extracting " + member.filename, index + 1, len(members)
                )
    elif str(archive).lower().endswith(".7z"):
        import py7zr

        with py7zr.SevenZipFile(archive) as handle:
            for member in handle.files:
                target(member.filename)
                if member.is_symlink:
                    raise ValueError("Archive symlinks are unsupported")
            _check(stop)
            progress("Extracting " + Path(archive).name, 0, 0)
            handle.extractall(path=destination)
            _check(stop)
    elif tarfile.is_tarfile(archive):
        total = Path(archive).stat().st_size
        archive_name = Path(archive).name

        class ArchiveReader:
            (
                'Report compressed bytes consumed during streaming '
                'extraction.'
            )

            def __init__(self, raw, message):
                self.raw = raw
                self.message = message
                self.high_water = 0
                progress(message, 0, total)

            def read(self, size=-1):
                _check(stop)
                data = self.raw.read(size)
                self.high_water = max(self.high_water, self.raw.tell())
                # Completion is reported only after the stage succeeds.
                progress(
                    self.message,
                    min(self.high_water, max(0, total - 1)),
                    total,
                )
                return data

            def __getattr__(self, name):
                return getattr(self.raw, name)

        extracting = "Extracting archive: " + archive_name
        with Path(archive).open("rb") as raw:
            # Validate each member immediately before writing
            # it; no index pass.
            with tarfile.open(
                fileobj=ArchiveReader(raw, extracting), mode="r|*"
            ) as handle:
                for member in handle:
                    _check(stop)
                    file = target(member.name)
                    if member.isdir():
                        file.mkdir(parents=True, exist_ok=True)
                        continue
                    if not member.isfile():
                        raise ValueError(
                            "Archive links and special files are"
                            " unsupported"
                        )
                    file.parent.mkdir(parents=True, exist_ok=True)
                    with handle.extractfile(member) as source, file.open(
                        "xb"
                    ) as dest:
                        while True:
                            _check(stop)
                            chunk = source.read(1024 * 1024)
                            if not chunk:
                                break
                            dest.write(chunk)
        progress(extracting, total, total)
    else:
        raise ValueError("Choose a ZIP or TAR archive")
    return destination


def packaged_elevation_band_name(region):
    if region == "CAMELS_US":
        return "elev_bands671.npz"
    if region in ("HYDRO_CIS", "HYD_RESPONSES", "CAMELS_BR",
                  "CAMELS_COL", "CAMELS_NZ", "CAMELS_SE"):
        return "cartplus_terrain_bands.csv"
    return None


def install_elevation_bands(directory, region, sage_root):
    """Copy model elevation data from the packaged regional resources."""
    name = packaged_elevation_band_name(region)
    if name is None:
        raise ValueError(f"No packaged elevation bands for {region}.")
    target = Path(directory) / name
    if target.is_file():
        return target
    relative = Path("elevation_bands") / region / name
    app_root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[1]))
    packaged = app_root / "resources" / "sage"
    for root in (Path(sage_root), packaged):
        source = root / relative
        if source.is_file():
            shutil.copy2(source, target)
            return target
    raise FileNotFoundError(f"Missing packaged elevation bands for {region}.")


def missing_camels_us_nldas(directory, sage_root):
    """Return CAMELS-US basins without a local, nonempty daily NLDAS file."""
    from camels_us import forcing_files

    universe = Path(sage_root) / "regions/US/US/US_671_basins.txt"
    expected = set(universe.read_text(encoding="utf-8-sig").split())
    if len(expected) != 671:
        raise ValueError("Expected 671 IDs in the CAMELS-US basin list.")
    present = set()
    for path in forcing_files(directory, "nldas"):
        info = path.stat()
        if info.st_size > 0 and not (
            getattr(info, "st_file_attributes", 0) & 0x400000
        ):
            present.add(path.name[:8])
    return sorted(expected - present)


def install_archives(
    archives,
    parent,
    region,
    sage_root,
    progress=lambda *args: None,
    stop=None,
    resolution=1,
    existing_directory=None,
    keep_source_archives=False,
):
    """
    Keep raw payloads and create a new normalized regional
    installation.
    """
    parent = Path(parent).resolve()
    parent.mkdir(parents=True, exist_ok=True)
    lamah = region in ("LAMAH_CE", "LAMAH_AT")
    additive = existing_directory is not None
    if additive and region != "CAMELS_US":
        raise ValueError(
            "Additive installation is only supported for CAMELS_US"
        )
    directory = (
        Path(existing_directory).resolve()
        if additive
        else parent
        / (region if lamah else region + "_" + uuid.uuid4().hex[:10])
    )
    if additive and directory.parent != parent:
        raise ValueError(
            "Existing region must be directly inside the selected data root"
        )
    if lamah:
        from lamah_dataset import available

        if available(directory, region, resolution):
            if (
                not (directory / "gauge_information.txt").is_file()
                or not (directory / "lamah_installation.json").is_file()
            ):
                from lamah_dataset import finalize

                finalize(directory, region, resolution)
            return directory
    directory.mkdir(exist_ok=lamah or additive)
    payload = directory / "source_archives"
    if lamah or additive:
        payload = payload / uuid.uuid4().hex[:10]
    payload.mkdir(parents=True)
    manifest = dict(
        region=region,
        resolution=resolution,
        status="installing",
        archives=[],
        files=0,
    )
    marker = directory / (
        f"python_installation_{resolution}_{uuid.uuid4().hex[:10]}.json"
        if additive
        else "python_installation.json"
    )
    marker.write_text(json.dumps(manifest, indent=2))
    try:
        for i, archive in enumerate(archives):
            _check(stop)
            archive = Path(archive)
            if archive.is_dir():
                extracted = archive
            elif archive.suffix.lower() in (".csv", ".txt", ".xlsx"):
                extracted = payload / str(i)
                extracted.mkdir()
                shutil.copy2(archive, extracted / archive.name)
            else:
                extracted = extract(
                    archive, payload / str(i), progress, stop
                )
            if archive.is_file():
                checksum = _verified_archive(archive)
                if checksum is None:
                    message = "Checking archive: " + archive.name
                    total = archive.stat().st_size
                    progress(message, 0, total)
                    digest = hashlib.sha256()
                    done = 0
                    with archive.open("rb") as handle:
                        while True:
                            _check(stop)
                            chunk = handle.read(1024 * 1024)
                            if not chunk:
                                break
                            digest.update(chunk)
                            done += len(chunk)
                            progress(
                                message, min(done, max(0, total - 1)), total
                            )
                    checksum = digest.hexdigest()
                    progress(message, total, total)
                manifest["archives"].append(
                    dict(name=archive.name, sha256=checksum)
                )
            else:
                manifest["archives"].append(
                    dict(
                        name=archive.name, source_kind="provider-directory"
                    )
                )
            for depth in range(3):
                nested = [
                    p
                    for p in extracted.rglob("*")
                    if p.is_file()
                    and p.suffix.lower() in (".zip", ".7z")
                    and not p.with_name(p.name + ".contents").exists()
                ]
                if not nested:
                    break
                for inner in nested:
                    extract(
                        inner,
                        inner.with_name(inner.name + ".contents"),
                        progress,
                        stop,
                    )
            from lamah_dataset import REGIONS

            allowed = None
            if region in REGIONS:
                import csv

                gauge_file = next(extracted.rglob("Gauge_attributes.csv"))
                catch_file = next(
                    p
                    for p in extracted.rglob("Catchment_attributes.csv")
                    if not any(
                        part.startswith(("B_basins", "C_basins"))
                        for part in p.parts
                    )
                )
                with gauge_file.open(encoding="utf-8-sig") as handle:
                    gauges = list(csv.DictReader(handle, delimiter=";"))
                with catch_file.open(encoding="utf-8-sig") as handle:
                    catchments = {
                        row["ID"]
                        for row in csv.DictReader(handle, delimiter=";")
                    }
                allowed = {
                    row["ID"]
                    for row in gauges
                    if row["ID"] in catchments
                    and (region != "LAMAH_AT" or row["country"] == "AUT")
                }
            for source in extracted.rglob("*"):
                _check(stop)
                if not source.is_file():
                    continue
                if allowed is not None:
                    match = re.fullmatch(r"ID_(\d+)\.csv", source.name)
                    if match and match.group(1) not in allowed:
                        continue
                relative = _destination(
                    source, region, sage_root, resolution
                )
                if relative is None:
                    continue
                target = directory / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.exists():
                    if lamah and target.name in (
                        "Gauge_attributes.csv",
                        "Catchment_attributes.csv",
                    ):
                        from lamah_dataset import rows

                        published = {row["ID"]: row for row in rows(source)}
                        expected = {
                            key: row
                            for key, row in published.items()
                            if key in allowed
                        }
                        actual = {row["ID"]: row for row in rows(target)}
                        if actual != expected and actual != published:
                            raise ValueError(
                                f"Conflicting LamaH metadata: {relative}"
                            )
                        continue
                    with target.open("rb") as a, source.open("rb") as b:
                        if (
                            hashlib.file_digest(a, "sha256").digest()
                            != hashlib.file_digest(b, "sha256").digest()
                        ):
                            raise ValueError(
                                f"Conflicting archive files for {relative}"
                            )
                    progress("Already installed: " + str(relative), 1, 1)
                    continue
                shutil.copy2(source, target)
                manifest["files"] += 1
                progress("Installed: " + str(relative), 1, 1)
        if region in ("LAMAH_CE", "LAMAH_AT"):
            from lamah_dataset import finalize

            finalize(directory, region, resolution)
        if region == "CAMELS_NZ":
            prepare_nz_gauges(directory)
        if region == "CAMELS_IS":
            prepare_is_gauges(directory)
        if region == "CAMELSH_KR":
            prepare_krh_gauges(directory)
        if region == "CAMELS_CL":
            prepare_cl_gauges(directory)
        if region == "HYD_RESPONSES":
            prepare_hyd_gauges(directory)
        if region == "CAMELS_COL":
            prepare_col_gauges(directory)
        if region == "CAMELS_DK":
            prepare_dk_gauges(directory)
        if region == "CAMELS_BR":
            prepare_br_gauges(directory)
        if region == "CAMELS_PL":
            prepare_pl_gauges(directory)
        if region == "CAMELS_ES":
            prepare_es_gauges(directory)
        if region == "CAMELS_CH":
            prepare_ch_gauges(directory)
        if region == "CAMELS_JP":
            series = list(
                (directory / "daily/timeseries").glob("varssim*.csv")
            )
            if len(series) != 87:
                raise ValueError(
                    "Expected 87 MERV-Jp version 2.1 files; found"
                    f" {len(series)}"
                )
            name = "MERV_Jp_135_HydroATLAS_attributes.xlsx"
            target = directory / name
            if not target.exists():
                source = Path(sage_root) / "regions/JP" / name
                if not source.is_file():
                    raise ValueError(
                        "Missing packaged MERV-Jp HydroATLAS attributes"
                    )
                shutil.copy2(source, target)
            manifest["dataset_version"] = "2.1"
        if region == "CAMELS_US":
            prepare_camels_names(directory)
            if resolution == 1:
                missing = missing_camels_us_nldas(directory, sage_root)
                if missing:
                    raise ValueError(
                        f"Daily NLDAS installation incomplete: {len(missing)}"
                        " of 671 CAMELS-US basins lack a local forcing file."
                    )
                manifest["daily_nldas_basins"] = 671
        from regional_data import catalogue

        gauges = catalogue(
            directory, region, sage_root, resolution=resolution
        )
        if not gauges:
            raise ValueError(
                "Archives do not contain recognized regional basin inputs;"
                " source payloads retained for diagnosis."
            )
        if packaged_elevation_band_name(region):
            bands = install_elevation_bands(directory, region, sage_root)
            manifest["cartplus_elevation_bands"] = bands.name
        manifest.update(status="installed", basins=len(gauges))
    except Exception as exc:
        manifest.update(
            status=(
                "canceled"
                if isinstance(exc, InstallCancelled)
                else "failed"
            ),
            error=str(exc),
        )
        raise
    finally:
        marker.write_text(json.dumps(manifest, indent=2))
    relative_payload = payload.relative_to(directory)
    directory = finalize_install_directory(directory, parent, region)
    if not keep_source_archives:
        from installer_cleanup import cleanup_extraction

        try:
            cleanup_extraction(directory, relative_payload)
            progress("Temporary extraction files removed", 1, 1)
        except OSError:
            progress(
                "Temporary extraction files retained: folder is in use",
                1,
                1,
            )
    return directory


def finalize_install_directory(directory, parent, region):
    """Promote a successful staging folder; never replace existing data."""
    directory, parent = Path(directory).resolve(), Path(parent).resolve()
    if (
        not re.fullmatch(r"[A-Za-z0-9_]+", region)
        or directory.parent != parent
    ):
        raise ValueError("Unexpected regional installation path")
    if directory.name == region:
        return directory
    if not re.fullmatch(
        re.escape(region) + r"_[0-9a-f]{10}", directory.name
    ):
        return directory
    marker = directory / "python_installation.json"
    if (
        not marker.is_file()
        or json.loads(marker.read_text()).get("status") != "installed"
    ):
        return directory
    target = parent / region
    if target.exists():
        return directory
    try:
        directory.rename(target)
    except FileExistsError:
        # Another installation may have completed while this one ran.
        return directory
    return target


def prepare_nz_gauges(directory):
    (
        'Create the MATLAB-compatible lookup using station IDs, not '
        'river IDs.'
    )
    import pandas as pd

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    frame = pd.read_csv(
        directory / "1.CAMELS_NZ_Catchment_information.csv",
        dtype={"Station_ID": str},
    )
    frame = frame.rename(
        columns={
            "Station_ID": "gauge_id",
            "Station Name": "gauge_name",
            "Latitude (WGS 84)": "gauge_lat",
            "Longitude(WGS 84)": "gauge_lon",
            "elevation": "gauge_elev",
            "uparea": "area_km2",
        }
    )
    with target.open("x", encoding="utf-8", newline="") as handle:
        frame[[
            "gauge_id",
            "gauge_name",
            "gauge_lat",
            "gauge_lon",
            "gauge_elev",
            "area_km2",
        ]].to_csv(handle, sep="\t", index=False)


def prepare_is_gauges(directory):
    (
        'Generate MATLAB Iceland basin metadata from published source '
        'tables.'
    )
    import pandas as pd
    from pyproj import Transformer

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    gauge = pd.read_csv(
        directory / "Gauge_attributes.csv", sep=";", dtype={"id": str}
    )
    catchment = pd.read_csv(
        directory / "Catchment_attributes.csv",
        sep=";",
        usecols=["id", "area_calc"],
        dtype={"id": str},
    )
    frame = gauge.merge(
        catchment, on="id", how="inner", validate="one_to_one"
    )
    lon, lat = Transformer.from_crs(3057, 4326, always_xy=True).transform(
        frame.lon.to_numpy(), frame.lat.to_numpy()
    )
    result = pd.DataFrame({
        "gauge_id": frame.id,
        "gauge_name": frame["name"],
        "gauge_lat": lat,
        "gauge_lon": lon,
        "gauge_elev": frame.elevation,
        "area_km2": frame.area_calc,
    })
    with target.open("x", encoding="utf-8", newline="") as handle:
        result.to_csv(handle, sep="\t", index=False)


def prepare_krh_gauges(directory):
    (
        'Provide the MATLAB Korea-hourly basin metadata to the master '
        'reader.'
    )
    import pandas as pd

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    frame = pd.read_csv(
        directory / "attributes_general.csv", dtype={"STAID": str}
    )
    elevation = pd.read_csv(
        directory / "attributes_HydroATLAS.csv",
        usecols=["STAID", "ele_mt_sav"],
        dtype={"STAID": str},
    )
    frame = frame.drop_duplicates("STAID", keep="first")
    elevation = elevation.drop_duplicates("STAID", keep="first")
    frame = frame.merge(
        elevation, on="STAID", how="left", validate="one_to_one"
    )
    frame = frame.rename(
        columns={
            "STAID": "gauge_id",
            "Area": "area_km2",
            "Lat_gage": "gauge_lat",
            "Lon_gage": "gauge_lon",
            "ele_mt_sav": "gauge_elev",
        }
    )
    with target.open("x", encoding="utf-8", newline="") as handle:
        frame[
            ["gauge_id", "area_km2", "gauge_lat", "gauge_lon", "gauge_elev"]
        ].to_csv(handle, sep="\t", index=False)


def prepare_cl_gauges(directory):
    """Transpose published Chile descriptors into the reader gauge table."""
    import pandas as pd

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    frame = pd.read_csv(
        directory / "1_CAMELScl_attributes.txt", sep="\t", dtype=str
    )
    frame = (
        frame.set_index("gauge_id").T.rename_axis("gauge_id").reset_index()
    )
    frame["gauge_id"] = frame.gauge_id.str.strip()
    frame = frame.rename(
        columns={"area": "area_km2", "elev_gauge": "gauge_elev"}
    )
    columns = [
        "gauge_id",
        "gauge_name",
        "gauge_lat",
        "gauge_lon",
        "area_km2",
        "gauge_elev",
    ]
    with target.open("x", encoding="utf-8", newline="") as handle:
        frame[columns].to_csv(handle, sep="\t", index=False)


def prepare_hyd_gauges(directory):
    """Materialize MATLAB HYD basin metadata for the master reader."""
    import pandas as pd
    from pyproj import Transformer

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    source = (
        directory
        / "HYDRESPONSES_descriptive_general_station_information.csv"
    )
    frame = pd.read_csv(
        source, sep=";", dtype={"catchmentnr": str}, na_values=["NA"]
    )
    lon, lat = Transformer.from_crs(21781, 4326, always_xy=True).transform(
        frame.x_coord.to_numpy(), frame.y_coord.to_numpy()
    )
    result = pd.DataFrame({
        "gauge_id": frame.catchmentnr.str.zfill(4),
        "gauge_name": (
            frame.water_name.fillna("") + " / " + frame.place.fillna("")
        ),
        "gauge_lat": lat,
        "gauge_lon": lon,
        "gauge_elev": frame.FOEN_stat_H,
        "area_km2": frame.FOEN_area,
    })
    if result.gauge_id.duplicated().any():
        raise ValueError("Duplicate HYD-RESPONSES gauge identifiers")
    with target.open("x", encoding="utf-8", newline="") as handle:
        result.to_csv(handle, sep="\t", index=False)


def prepare_col_gauges(directory):
    (
        'Normalize Colombian Web Mercator station coordinates and '
        'basin areas.'
    )
    import pandas as pd
    from pyproj import Transformer

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    location = pd.read_csv(
        directory / "02_CAMELS_COL_Catchment_information.csv",
        sep=";",
        encoding="latin-1",
        dtype={"gauge_id": str},
    )
    topo = pd.read_csv(
        directory / "10_CAMELS_COL_Physiograpic_characteristics.csv",
        sep=";",
        encoding="latin-1",
        dtype={"gauge_id": str},
    )
    frame = location.merge(
        topo[["gauge_id", "area"]], on="gauge_id", validate="one_to_one"
    )
    projection = Transformer.from_crs(3857, 4326, always_xy=True)
    lon, lat = projection.transform(
        frame.gauge_lon.to_numpy(), frame.gauge_lat.to_numpy()
    )
    frame["gauge_lat"], frame["gauge_lon"] = lat, lon
    frame = frame.rename(
        columns={"area": "area_km2", "gauge_department": "gauge_name"}
    )
    if len(frame) != len(location) or frame.gauge_lat.isna().any():
        raise ValueError("Incomplete Colombian gauge metadata")
    with target.open("x", encoding="utf-8", newline="") as handle:
        frame[[
            "gauge_id",
            "gauge_name",
            "gauge_lat",
            "gauge_lon",
            "gauge_elev",
            "area_km2",
        ]].to_csv(handle, sep="\t", index=False)


def prepare_dk_gauges(directory):
    """Match MATLAB's Danish UTM projection and catchment-area units."""
    import pandas as pd
    from pyproj import Transformer

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    source = pd.read_csv(
        directory / "CAMELS_DK_topography.csv", dtype={"catch_id": str}
    )
    available_ids = {
        p.stem.removeprefix("CAMELS_DK_obs_based_")
        for p in (directory / "daily/timeseries").glob(
            "CAMELS_DK_obs_based_*.csv"
        )
    }
    source = source[source.catch_id.isin(available_ids)]
    if len(source) != len(available_ids):
        raise ValueError(
            "Danish observed basins missing topographic metadata"
        )
    projection = Transformer.from_crs(25832, 4326, always_xy=True)
    lon, lat = projection.transform(
        source.catch_outlet_lon.to_numpy(),
        source.catch_outlet_lat.to_numpy(),
    )
    frame = pd.DataFrame({
        "gauge_id": source.catch_id,
        "gauge_name": "CAMELS-DK " + source.catch_id,
        "gauge_lat": lat,
        "gauge_lon": lon,
        "gauge_elev": source.elev_mean,
        "area_km2": source.catch_area / 1e6,
    })
    if frame.gauge_id.duplicated().any() or frame.gauge_lat.isna().any():
        raise ValueError("Invalid Danish gauge metadata")
    with target.open("x", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, sep="\t", index=False)


def prepare_br_gauges(directory):
    (
        'Join Brazilian station locations with the published '
        'modeled-basin areas.'
    )
    import pandas as pd

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    location = pd.read_csv(
        directory / "camels_br_location.txt",
        sep=r"\s+",
        dtype={"gauge_id": str},
    )
    topo = pd.read_csv(
        directory / "camels_br_topography.txt",
        sep=r"\s+",
        dtype={"gauge_id": str},
    )
    frame = topo[["gauge_id", "elev_gauge", "area"]].merge(
        location[["gauge_id", "gauge_name", "gauge_lat", "gauge_lon"]],
        on="gauge_id",
        validate="one_to_one",
    )
    frame = frame.rename(
        columns={"area": "area_km2", "elev_gauge": "gauge_elev"}
    )
    if len(frame) != len(topo) or frame.gauge_lat.isna().any():
        raise ValueError("Incomplete Brazilian gauge metadata join")
    with target.open("x", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, sep="\t", index=False)


def prepare_pl_gauges(directory):
    (
        "Match MATLAB's Polish coordinate aliases and auxiliary reader "
        'fields.'
    )
    import pandas as pd

    directory = Path(directory)
    target = directory / "gauge_information_PL.csv"
    if target.exists():
        return
    frame = pd.read_csv(
        directory / "CAMELS_PL_topographic_attributes.csv",
        dtype={"gauge_id": str},
    )
    frame = frame.rename(
        columns={
            "gauge_lon": "gauge_lat",
            "gauge_lat": "gauge_lon",
            "area_metadata": "area_km2",
        }
    )
    if frame.gauge_id.duplicated().any() or frame.gauge_lat.isna().any():
        raise ValueError("Invalid Polish gauge metadata")
    with target.open("x", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, index=False)


def prepare_es_gauges(directory):
    (
        "Create canonical Spanish gauge metadata from the publisher's "
        'attributes.'
    )
    import pandas as pd

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    frame = pd.read_csv(
        directory / "attributes_other_camelses.csv", dtype={"gauge_id": str}
    )
    frame["gauge_id"] = frame.gauge_id.str.replace(
        r"^camelses[_-]*", "", regex=True
    )
    frame = frame.rename(columns={"area": "area_km2"})
    atlas = pd.read_csv(
        directory / "attributes_hydroatlas_camelses.csv",
        dtype={"gauge_id": str},
    )
    atlas["gauge_id"] = atlas.gauge_id.str.replace(
        r"^camelses[_-]*", "", regex=True
    )
    frame["gauge_elev"] = frame.gauge_id.map(
        atlas.set_index("gauge_id")["ele_mt_sav"]
    )
    if frame.gauge_id.duplicated().any() or frame.gauge_lat.isna().any():
        raise ValueError("Invalid Spanish gauge metadata")
    with target.open("x", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, sep="\t", index=False)


def prepare_ch_gauges(directory):
    (
        'Materialize the Swiss metadata expected by the shared hydro '
        'schema.'
    )
    import pandas as pd

    directory = Path(directory)
    target = directory / "gauge_information.txt"
    if target.exists():
        return
    source = directory / "CAMELS_CH_topographic_attributes.csv"
    frame = pd.read_csv(
        source, comment="#", encoding="latin-1", dtype={"gauge_id": str}
    )
    columns = {
        "gauge_id": "gauge_id",
        "gauge_name": "gauge_name",
        "gauge_lat": "gauge_lat",
        "gauge_lon": "gauge_lon",
        "gauge_elevation": "gauge_elev",
        "area": "area_km2",
    }
    frame = frame[list(columns)].rename(columns=columns)
    if frame.gauge_id.duplicated().any() or frame.gauge_lat.isna().any():
        raise ValueError("Invalid Swiss gauge metadata")
    with target.open("x", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, sep="\t", index=False)


def prepare_camels_names(directory):
    (
        'Create required display metadata for a fresh daily or hourly '
        'install.'
    )
    import csv
    from gauge_names import standardize

    directory = Path(directory)
    target = directory / "camels_name_clean.txt"
    source = directory / "camels_name.txt"
    if target.exists() or not source.is_file():
        return
    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        if not {"gauge_id", "gauge_name"}.issubset(reader.fieldnames or []):
            raise ValueError(
                "CAMELS-US basin-name table has invalid columns"
            )
        rows = list(reader)
    if not rows or len({r["gauge_id"] for r in rows}) != len(rows):
        raise ValueError(
            "CAMELS-US basin-name table is empty or contains duplicate IDs"
        )
    with target.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=";")
        writer.writerow(["gauge_id", "huc_02", "gauge_name"])
        for row in rows:
            gauge = row["gauge_id"].strip().zfill(8)
            writer.writerow([
                gauge,
                row.get("huc_02", gauge[:2]),
                standardize(row["gauge_name"], "CAMELS_US", gauge),
            ])


@lru_cache(maxsize=128)
def attribute_destinations(region, resolution=1):
    (
        'Use the same exported attribute-table contract as the SAGE '
        'GUI reader.'
    )
    import sys

    base = Path(
        getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[1])
    )
    path = base / "gui_python/sage_region_metadata.json"
    if not path.is_file():
        return {}
    schema = (
        json.loads(path.read_text(encoding="utf-8"))
        .get(region, {})
        .get("attributes", {})
    )
    if resolution == 24:
        schema.update(
            schema.get("profiles", {}).get("hourly", {}).get("schema", {})
        )
    tables = schema.get("tables", [])
    if isinstance(tables, dict):
        tables = [tables]
    return {
        Path(t["file"].replace("\\", "/")).name: Path(
            t["file"].replace("\\", "/")
        )
        for t in tables
        if t.get("file")
    }


def _destination(source, region, sage_root, resolution=1):
    """
    Normalize published filenames into the paths used by regional
    readers.
    """
    from lamah_dataset import reader_region

    region = reader_region(region)
    name = source.name
    if name.endswith((".lock", ".download.json", ".resume.json", ".part")):
        return None
    # Basin-area alternatives have identically named attribute tables;
    # choose
    # total-upstream areas before applying the attribute manifest.
    if region in ("LAMAH_CE", "CAMELS_IS") and any(
        part.startswith(("B_basins", "C_basins")) for part in source.parts
    ):
        return None
    declared_attribute = attribute_destinations(region, resolution).get(
        name
    )
    if declared_attribute is not None:
        return declared_attribute
    parts = source.parts
    if region == "CAMELS_AU":
        if name == "Indices_MasterTable.csv":
            return Path("CAMELS_AUS_Attributes&Indices_MasterTable.csv")
        if name in ("precipitation_SILO.csv", "precipitation_AGCD.csv"):
            return Path("daily/precipitation") / name
        if name in (
            "tmin_SILO.csv",
            "tmax_SILO.csv",
            "tmin_AGCD.csv",
            "tmax_AGCD.csv",
        ):
            return Path("daily/temperature") / name
        if name in (
            "evap_syn_SILO.csv",
            "evap_pan_SILO.csv",
            "evap_morton_lake_SILO.csv",
            "et_tall_crop_SILO.csv",
            "et_short_crop_SILO.csv",
            "et_morton_wet_SILO.csv",
            "et_morton_point_SILO.csv",
            "et_morton_actual_SILO.csv",
        ):
            return Path("daily/evaporative_demand") / name
    if region == "CAMELS_JP" and re.fullmatch(r"varssim\d+\.csv", name):
        return (
            Path("daily/timeseries") / name if "ver2_1" in parts else None
        )
    if region == "CAMELS_SE" and re.fullmatch(
        r"catchment_id_.*\.csv", name
    ):
        return Path("daily") / name
    profile = table_profile(sage_root, region)
    if region in ("LAMAH_CE", "CAMELS_IS") and any(
        part.startswith(("B_basins", "C_basins")) for part in parts
    ):
        return None
    allowed = _country_basins(sage_root, region)
    if (
        allowed is not None
        and name.startswith("GRDC_")
        and source.stem.removeprefix("GRDC_") not in allowed
    ):
        return None
    if region in ("LAMAH_CE", "CAMELS_IS") and re.fullmatch(
        r"ID_\d+\.csv", name
    ):
        if region == "LAMAH_CE" and any(
            part.startswith("F_hydrol_model") for part in parts
        ):
            return None
        if region == "CAMELS_IS":
            if "D_gauges" in parts:
                if source.parent.name == "daily":
                    preferred = source.parent.parent / "daily_filtered"
                    if preferred.is_dir() and any(
                        preferred.glob("ID_*.csv")
                    ):
                        return None
                if source.parent.name not in (
                    "daily",
                    "daily_filtered",
                    "hourly",
                ):
                    return None
            elif not (
                "A_basins_total_upstrm" in parts
                and "meteorological_data" in parts
            ):
                return None
        with source.open(encoding="utf-8-sig") as handle:
            header = handle.readline().lower()
        discharge = "qobs" in header
        hourly = bool(re.search(r"(?:^|[;,])hh(?:[;,]|$)", header))
        if region == "CAMELS_IS" and discharge and not hourly:
            hourly = (
                any("hourly" in part.lower() for part in parts)
                or resolution == 24
            )
        stream = "hourly" if hourly else "daily"
        base = (
            Path(stream) / "timeseries"
            if region == "LAMAH_CE"
            else Path(stream)
        )
        return (
            base
            / (
                "streamflow"
                if discharge and region == "LAMAH_CE"
                else "discharge" if discharge else "forcing"
            )
            / name
        )
    if region == "CAMELS_GB" and name.startswith(
        "camels_gb_v2_hydromet_daily_timeseries_"
    ):
        return Path("daily/timeseries") / name.replace(
            "camels_gb_v2_hydromet_daily_timeseries_",
            "CAMELS_GB_hydromet_timeseries_",
        )
    if region == "CAMELS_GB" and name.startswith(
        "camels_gb_v2_hydromet_hourly_timeseries_"
    ):
        return Path("hourly/timeseries") / name
    if (
        region == "CAMELS_CZ"
        and name.startswith("camelscz_")
        and source.suffix.lower() == ".csv"
    ):
        with source.open(encoding="utf-8-sig") as handle:
            header = handle.readline().strip().split(",")
        # Daily companion files can occur during an hourly installation.
        hourly = any(column.strip().strip('"').lower() == "time" for column in header)
        return (
            Path("hourly/timeseries" if hourly else "daily/timeseries")
            / name
        )
    if region == "CAMELS_US":
        match = re.fullmatch(
            r"(\d{8})_lump_(cida|daymet|maurer|nldas)_forcing_leap\.txt",
            name,
        )
        if match:
            return (
                Path("daily/v1p2/forcing")
                / ("daymet" if match.group(2) == "cida" else match.group(2))
                / match.group(1)[:2]
                / name
            )
        if name.endswith("_streamflow_qc.txt"):
            return Path("daily/v1p2/streamflow") / name[:2] / name
        if name.endswith("_hourly_nldas.csv"):
            return Path("hourly/forcing") / name
        if name.endswith("-usgs-hourly.csv"):
            return Path("hourly/streamflow") / name
        if name.startswith("camels_"):
            return Path(name)
    if (
        region in ("CAMELSH_US", "CAMELSH_KR")
        and source.suffix in (".nc", ".csv")
        and re.fullmatch(r"\d+", source.stem)
    ):
        return Path("hourly/timeseries") / name
    if region == "CAMELS_DEH" and name.startswith(
        "CAMELS_DE_1h_hydromet_timeseries_"
    ):
        return Path("hourly/timeseries") / name
    if region == "CAMELS_LUX" and (
        "_hourly_ID_" in name or "_15min_ID_" in name
    ):
        return (
            Path(
                "hourly/timeseries"
                if "_hourly_ID_" in name
                else "15min/timeseries"
            )
            / name
        )
    if region == "CAMELS_IND":
        if re.fullmatch(r"\d+\.csv", name):
            return Path("daily/catchment_mean_forcings") / name
        if name == "streamflow_observed.csv":
            return Path("daily/streamflow_timeseries") / name
    pattern = profile["pattern"]
    if pattern:
        regex = (
            re.escape(Path(pattern).name)
            .replace(re.escape("{gauge}"), ".+")
            .replace(re.escape("{gauge02}"), ".+")
            .replace(r"\*", ".*")
        )
        if re.fullmatch(regex, name):
            if region == "CAMELS_KR":
                return Path("daily/forcing") / name
            return Path("daily/timeseries") / Path(pattern).parent / name
    if name.startswith("CAMELS_KR_Hydrological_timeseries_"):
        return Path("daily/discharge") / name
    if region in ("CAMELS_CA", "CAMELS_SPAT") and name.endswith(".nc"):
        stream = (
            "hourly"
            if "_rdrs_" in name or "_hourly_flow_" in name
            else "daily"
        )
        if "flow_observations" in name:
            return Path(stream) / "discharge" / name
        return (
            Path(stream)
            / "forcing"
            / ("daymet" if stream == "daily" else "rdrs")
            / name
        )
    if region == "HYDRO_CIS" and name.endswith(".nc"):
        return Path("daily/timeseries") / name
    if region == "CAMELS_CH":
        if name.startswith("CAMELS_CH_obs_based_"):
            return Path("daily/timeseries/observation_based") / name
        if name.startswith("CAMELS_CH_sim_based_"):
            return Path("daily/timeseries/simulation_based") / name
    if region == "CAMELS_NZ":
        for folder, token in [
            ("precipitation", "precipitation"),
            ("temperature", "temperature"),
            ("relative_humidity", "RH"),
            ("streamflow", "flow"),
        ]:
            if name.startswith("daily_" + token + "_station_id_"):
                return Path("daily/timeseries") / folder / name
        for folder, token in [
            ("precipitation", "precipitation"),
            ("temperature", "temperature"),
            ("pet", "PET"),
            ("streamflow", "flow"),
        ]:
            if name.startswith(token + "_station_id_"):
                return Path("hourly/timeseries") / folder / name
    if region == "CAMELS_BR":
        for folder, token in [
            ("precipitation", "_precipitation"),
            ("pet", "_potential_evapotransp"),
            ("temperature", "_temperature"),
            ("streamflow", "_streamflow"),
        ]:
            if re.fullmatch(r"\d+" + token + r"\.txt", name):
                return Path("daily") / folder / name
    if region == "CAMELS_CL" and name == "13_CAMELScl_swe.txt":
        return Path("daily/swe") / name
    if region in ("CAMELS_AU", "CAMELS_CL"):
        for folder in (
            "precipitation",
            "temperature",
            "evaporative_demand",
            "pet",
            "streamflow",
            "swe",
        ):
            if folder in parts:
                return Path("daily") / folder / name
        if (
            name == "streamflow_mmd.csv"
            or name == "3_CAMELScl_streamflow_mm.txt"
        ):
            return Path("daily/streamflow") / name
        if name.startswith(
            ("4_CAMELScl", "5_CAMELScl", "6_CAMELScl", "7_CAMELScl")
        ):
            return Path("daily/precipitation") / name
        if name.startswith(("8_CAMELScl", "9_CAMELScl", "10_CAMELScl")):
            return Path("daily/temperature") / name
        if name.startswith(("11_CAMELScl", "12_CAMELScl")):
            return Path("daily/pet") / name
    if region == "BULL_ES":
        for prefix, folder in [
            ("AEMET_", "AEMET"),
            ("streamflow_", "streamflow"),
            ("ERA5_Land_", "ERA5_Land"),
            ("EMO1_", "EMO1_arc"),
        ]:
            if name.startswith(prefix):
                return Path("daily/timeseries") / folder / name
    if source.suffix.lower() in (".csv", ".txt", ".xlsx") and (
        "attribute" in name.lower()
        or "gauge" in name.lower()
        or "station" in name.lower()
        or "catchment" in name.lower()
        or "topograph" in name.lower()
        or name.startswith("camels_br_")
        or name in ("static_data.csv", "camels-spat-metadata.csv")
    ):
        return Path(name)
    return None


_legacy_download_options = download_options


def download_options(sage_root, region, resolution=1):
    from runtime_metadata import load
    from copy import deepcopy

    from lamah_dataset import reader_region

    region = reader_region(region)
    data = load(str(sage_root))
    items = (
        deepcopy(data["downloads"][region + ":" + str(resolution)])
        if data
        else _legacy_download_options(sage_root, region, resolution)
    )

    if region == "CAMELS_NZ":
        # The public Figshare API serves the same file IDs without the
        # Canterbury website's interactive browser challenge.
        checksums = {
            "56902355": "00f877983df4be56fe490e98ed44fa84",
            "56902367": "4871d33a5d175feb89a9f81f96b3c95c",
            "56902370": "ddfff72bed22a6ad12cf52fd2d4ca196",
            "56902373": "089757d4b019487fefd8f20d7099403d",
            "56902382": "c3221428156d8798515f0d83c53f531d",
        }
        for item in items:
            match = re.search(r"/files/(\d+)", item.get("url", ""))
            if match:
                file_id = match.group(1)
                item["url"] = (
                    "https://api.figshare.com/v2/file/download/" + file_id
                )
                if file_id in checksums:
                    item["md5"] = checksums[file_id]

    if region == "CAMELS_CZ":
        # MATLAB selects one self-contained archive for the chosen
        # resolution.
        wanted = (
            "CAMELS-CZ-Hourly.zip" if resolution == 24 else "CAMELS-CZ.zip"
        )
        items = [item for item in items if item.get("name") == wanted]

    if region == "CAMELS_AU":
        for item in items:
            if item.get("name") == "Indices_MasterTable.csv":
                item["url"] = (
                    'https://zenodo.org/records/14289037/'
                    'files/CAMELS_AUS_Attribu'
                    'tes%26Indices_MasterTable.csv?download=1'
                )
                item["md5"] = "aa47ba598d0486d5ea4ccca6e132a7be"

    if region == "CAMELS_SE":
        # Current file URLs linked by the provider catalogue
        # (legacy query URLs return 404).
        from urllib.parse import quote

        for item in items:
            item["url"] = (
                "https://api.researchdata.se/dataset/2023-173/1/file/data/"
                + quote(item["name"], safe="")
            )

    if region == "CAMELS_US" and resolution == 24:
        checksums = {
            "nldas_hourly_csv.tar.gz": "c7e9e5d92af2153ca83114b70fe315b8",
            "usgs_streamflow_csv.tar.gz": (
                "9b1f154c1504b06bc4d1504a0afd4a08"
            ),
        }
        for item in items:
            url = urlparse(item.get("url", ""))
            if url.hostname == "zenodo.org" and url.path.startswith(
                "/records/4072701/files/"
            ):
                if item.get("name") in checksums:
                    item["md5"] = checksums[item["name"]]
    return items
