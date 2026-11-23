# SITEhydrology | Python source
# ==========================================================================
# SAGEhydrology | Python source
# File: python/provider_install.py
#
# Provider acquisition workflows ported from SAGE data_helpers.m.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Provider acquisition workflows ported from SAGE data_helpers.m."""

from pathlib import Path
import json
import re
import subprocess
import sys
import time
from urllib.parse import urljoin, quote
from urllib.request import Request, urlopen
import pandas as pd
from dataset_install import download, extract, _check, InstallCancelled

PROVIDERS = {
    "CAMELS_GB",
    "CAMELS_CA",
    "CAMELS_SPAT",
    "CAMELS_ET",
    "EStreams_PT",
    "CAMELS_IS",
}


def fetch(url):
    with urlopen(
        Request(url, headers={"User-Agent": "SITE-Python/1.0"}), timeout=120
    ) as response:
        return response.read()


def directory_files(url, pattern):
    text = fetch(url).decode("utf-8")
    links = re.findall(r'href=["\']([^"\']+)', text, re.I)
    return list(
        dict.fromkeys(
            urljoin(url, link) for link in links if re.search(pattern, link)
        )
    )


def acquire(
    region, stage, sage_root, resolution, progress, stop, options=None
):
    stage = Path(stage)
    stage.mkdir(parents=True, exist_ok=True)
    options = options or {}
    if region == "CAMELS_GB":
        return great_britain(stage, resolution, progress, stop)
    if region in ("CAMELS_CA", "CAMELS_SPAT"):
        return globus(stage, region, resolution, progress, stop, options)
    if region == "CAMELS_ET":
        return ethiopia(stage, progress, stop, options)
    if region == "EStreams_PT":
        return portugal(stage, sage_root, progress, stop)
    if region == "CAMELS_IS":
        return iceland(stage, resolution, progress, stop)
    raise ValueError(region)


def iceland(stage, resolution, progress, stop):
    """
    Use selected content, API fallback, then the resource bag like
    SAGE.
    """
    import zipfile

    resource = "705d69c0f77c48538d83cf383f8c63d6"
    name = "lamah_ice" + ("_hourly" if resolution == 24 else "") + ".zip"
    urls = [
        (
            "https://www.hydroshare.org/resource/"
            f"{resource}"
            "/data/contents/"
            f"{name}"
        ),
        (
            "https://www.hydroshare.org/hsapi/resource/"
            f"{resource}"
            "/files/"
            f"{name}"
        ),
    ]
    for i, url in enumerate(urls):
        _check(stop)
        try:
            archive = download(
                dict(name=name, url=url), stage / str(i), progress, stop
            )
            if zipfile.is_zipfile(archive):
                return archive
        except InstallCancelled:
            raise
        except (OSError, ValueError):
            pass
    archive = download(
        dict(
            name="resource-bag.zip",
            url=(
                "https://www.hydroshare.org/django_irods/download/bags/"
                f"{resource}"
                ".zip"
            ),
        ),
        stage / "bag",
        progress,
        stop,
    )
    selected = stage / name
    with zipfile.ZipFile(archive) as handle:
        matches = [
            m for m in handle.infolist() if Path(m.filename).name == name
        ]
        if len(matches) != 1:
            raise ValueError(
                "HydroShare resource bag does not contain the selected"
                " LamaH-Ice archive"
            )
        with handle.open(matches[0]) as source, selected.open(
            "xb"
        ) as target:
            done = 0
            while True:
                _check(stop)
                chunk = source.read(1024 * 1024)
                if not chunk:
                    break
                target.write(chunk)
                done += len(chunk)
                progress(name, done, matches[0].file_size)
    return selected


def great_britain(stage, resolution, progress, stop):
    base = (
        "https://catalogue.ceh.ac.uk/datastore/eidchub/"
        "9a46d428-958f-4ac1-86eb-94eee70c0955/"
    )
    urls = directory_files(base + "Catchment_Attributes/", r"\.csv(?:\?|$)")
    # Hourly forcing needs the companion daily temperature/PET tables.
    for stream in (["daily", "hourly"] if resolution == 24 else ["daily"]):
        urls += directory_files(
            base
            + "Catchment_Timeseries/hydro-meteorological/"
            + stream
            + "/",
            r"\.csv(?:\?|$)",
        )
    if not urls:
        raise ValueError("The provider returned no GB data-file links")
    for i, url in enumerate(urls):
        _check(stop)
        name = url.split("/")[-1].split("?")[0]

        def file_progress(message, done, total):
            if "attribute" in name.lower() or message.startswith(
                "Retrying "
            ):
                progress(message, done, total)
            fraction = min(1.0, done / total) if total else 0.0
            progress(
                "Dataset files",
                int((i + fraction) * 1000000),
                len(urls) * 1000000,
            )

        file_progress(name, 0, 1)
        download(dict(name=name, url=url), stage, file_progress, stop)
        file_progress(name, 1, 1)
    return stage


def ethiopia(stage, progress, stop, options):
    if not options.get("registered"):
        raise ValueError(
            "Complete the CAMELS-Eth provider registration and tick the"
            " registration checkbox before downloading."
        )
    base = "https://www.walris.wlrc-eth.org/camels-eth/"
    folder = stage / "daily" / "timeseries"
    folder.mkdir(parents=True)
    for kind in ("topographic", "soil", "landcover"):
        name = "CAMELS_ET_" + kind + "_attributes.csv"
        download(
            dict(name=name, url=base + "downloads/" + name),
            stage,
            progress,
            stop,
        )
    table = pd.read_csv(
        stage / "CAMELS_ET_topographic_attributes.csv",
        dtype={"gauge_id": str},
    )
    ids = table["gauge_id"].str.upper().str.strip().tolist()
    timeline = pd.date_range("1981-01-01", "2018-12-31")
    for i, gauge in enumerate(ids):
        _check(stop)
        frame = pd.DataFrame(index=timeline)
        for endpoint, key in [
            ("chirps", "precipitation"),
            ("chirts-tmin", "tmin"),
            ("chirts-tmax", "tmax"),
            ("gleam-ep", "pet"),
            ("runoff", "runoff"),
        ]:
            _check(stop)
            payload = json.loads(
                fetch(
                    base
                    + "api/public/"
                    + endpoint
                    + "?code="
                    + quote(gauge)
                )
            )
            embedded = next(
                (
                    value
                    for key, value in payload.items()
                    if "embedded" in key.lower()
                ),
                None,
            )
            if not isinstance(embedded, dict):
                raise ValueError(
                    "CAMELS-Eth response is missing embedded observations"
                )
            rows = next(
                (
                    value
                    for key, value in embedded.items()
                    if "link" not in key.lower()
                ),
                [],
            )
            values = pd.DataFrame(rows)
            if values.empty:
                frame[key] = float("nan")
                continue
            dates = pd.to_datetime(values["date"]).dt.normalize()
            numbers = pd.to_numeric(values["value"], errors="coerce")
            series = (
                pd.Series(numbers.to_numpy(), index=dates)
                .groupby(level=0)
                .first()
            )
            frame[key] = series.reindex(timeline)
        frame.index.name = "date"
        frame.to_csv(folder / ("CAMELS_ET_" + gauge + ".csv"))
        progress("Ethiopia " + gauge, i + 1, len(ids))
    return stage


def portugal(stage, sage_root, progress, stop):
    base = "https://zenodo.org/records/17598150/files/"
    raw = stage / "raw"
    raw.mkdir()
    for name in (
        "attributes.zip",
        "streamflow_gauges.zip",
        "meteorology.zip",
    ):
        _check(stop)
        archive = download(
            dict(name=name, url=base + name + "?download=1"),
            stage,
            progress,
            stop,
        )
        extract(archive, raw / name.removesuffix(".zip"), progress, stop)
    mapping = pd.read_csv(
        Path(sage_root) / "regions" / "PT" / "snirh_station_map.csv",
        dtype=str,
    )
    output = stage / "prepared"
    folder = output / "daily" / "timeseries"
    folder.mkdir(parents=True)
    for i, row in mapping.iterrows():
        _check(stop)
        basin = row["basin_id"]
        files = list(raw.rglob("estreams_meteorology_" + basin + ".csv"))
        if len(files) != 1:
            raise ValueError("Portugal meteorology missing for " + basin)
        frame = pd.read_csv(files[0])
        frame["date"] = pd.to_datetime(frame["date"])
        url = (
            "https://snirh.apambiente.pt/snirh/_dadosbase/site/"
            "paraCSV/dados_csv.php?sites="
            + quote(row["gauge_code"])
            + "&pars=1850&tmin=01/01/1972&tmax=31/12/2022&formato=csv"
        )
        try:
            raw_body = fetch(url)
            try:
                body = raw_body.decode("utf-8-sig")
            except UnicodeError:
                body = raw_body.decode("cp1252")
        except Exception as exc:
            raise RuntimeError(
                "SNIRH download failed. The provider may require a"
                " Portuguese IP address; retry after establishing provider"
                " access. "
                + str(exc)
            ) from exc
        import io

        flow = pd.read_csv(
            io.StringIO(body),
            skiprows=3,
            header=None,
            dtype=str,
            sep=None,
            engine="python",
        )
        dates = pd.to_datetime(
            flow.iloc[:, 0], format="%d/%m/%Y %H:%M", errors="coerce"
        ).dt.normalize()
        values = pd.to_numeric(
            flow.iloc[:, 1].str.replace(",", ".", regex=False),
            errors="coerce",
        )
        valid = dates.notna() & values.notna() & (values >= 0)
        series = (
            pd.Series(values[valid].to_numpy(), index=dates[valid])
            .groupby(level=0)
            .first()
        )
        frame["q_cms_obs"] = series.reindex(frame["date"]).to_numpy()
        frame.to_csv(
            folder / ("EStreams_PT_" + basin + ".csv"), index=False
        )
        progress("Portugal " + basin, i + 1, len(mapping))
    ids = set(mapping["basin_id"])
    for source in raw.rglob("estreams_*.csv"):
        if "meteorology_" in source.name:
            continue
        table = pd.read_csv(source, dtype={"basin_id": str})
        if "basin_id" not in table:
            continue
        table = table[table["basin_id"].isin(ids)]
        table.to_csv(output / (source.stem + "_PT.csv"), index=False)
    return output


def cli(arguments, stop=None):
    from globus_tls import environment

    command = [sys.executable, "-m", "globus_cli", *arguments]
    process = subprocess.Popen(
        command,
        env=environment(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        creationflags=(
            subprocess.CREATE_NO_WINDOW if sys.platform == "win32" else 0
        ),
    )
    while True:
        try:
            output, error = process.communicate(timeout=0.2)
            break
        except subprocess.TimeoutExpired:
            if stop is not None and stop.is_set():
                process.terminate()
                process.communicate()
                raise InstallCancelled()
    if process.returncode:
        raise RuntimeError("Globus command failed: " + error.strip())
    return (
        json.loads(output)
        if output.strip().startswith(("{", "["))
        else output
    )


def globus(stage, region, resolution, progress, stop, options):
    destination = options.get("endpoint", "").strip()
    import uuid

    try:
        destination = str(uuid.UUID(destination))
    except ValueError:
        raise ValueError(
            "Enter the Globus Connect Personal collection UUID for this"
            " computer."
        ) from None
    source = "f163c1b3-9c88-42f6-a7bb-5839ed6c4063"
    origin = "/1/published/publication_1301/submitted_data/"
    local = str(stage.resolve()).replace("\\", "/")
    local = (
        "/" + local[0].upper() + local[2:] + "/"
        if re.match(r"^[A-Za-z]:/", local)
        else local + "/"
    )
    cli(["ls", destination + ":" + local, "--format", "json"], stop)

    def transfer(lines, label):
        batch = stage / (label + ".batch")
        batch.write_text("\n".join(a + " " + b for a, b in lines) + "\n")
        result = cli(
            [
                "transfer",
                source + ":" + origin,
                destination + ":" + local,
                "--batch",
                str(batch),
                "--sync-level",
                "checksum",
                "--label",
                "SITE Python " + label,
                "--format",
                "json",
            ],
            stop,
        )
        task = result["task_id"]
        try:
            while True:
                _check(stop)
                state = cli(
                    ["task", "show", task, "--format", "json"], stop
                )
                status = state["status"]
                progress("Globus " + label + " · " + status, 0, 0)
                if status == "SUCCEEDED":
                    break
                if status in ("FAILED", "INACTIVE"):
                    raise RuntimeError(
                        "Globus transfer " + status + ": " + task
                    )
                for _ in range(25):
                    _check(stop)
                    time.sleep(0.2)
        except InstallCancelled:
            cli(["task", "cancel", task])
            raise

    transfer(
        [
            ("camels-spat-metadata.csv", "camels-spat-metadata.csv"),
            ("attributes/attributes-lumped.csv", "attributes-lumped.csv"),
        ],
        "metadata",
    )
    metadata = stage / "camels-spat-metadata.csv"
    if not metadata.exists():
        raise RuntimeError(
            "Globus destination is not accessible on this computer; check"
            " collection and path permissions."
        )
    table = pd.read_csv(metadata, dtype=str)
    countries = (
        ["CAN", "CANADA"]
        if region == "CAMELS_CA"
        else ["CAN", "CANADA", "USA", "US", "UNITED STATES"]
    )
    table = table[
        table["Country"].str.upper().str.strip().isin(countries)
    ].copy()
    country = (
        table["Country"]
        .str.upper()
        .str.strip()
        .replace({"CANADA": "CAN", "US": "USA", "UNITED STATES": "USA"})
    )
    table["Station_id"] = (
        country
        + "_"
        + table["Station_id"]
        .str.upper()
        .str.strip()
        .str.replace(r"^(CAN|USA)[-_]?", "", regex=True)
    )
    if table.empty:
        raise ValueError("No basins for this region in Globus metadata")
    if table["Station_id"].duplicated().any():
        raise ValueError("Duplicate basin IDs in Globus metadata")
    import shutil

    shutil.copyfile(metadata, stage / "camels-spat-metadata.csv.raw")
    table.to_csv(metadata, index=False)
    stream = "daily" if resolution == 1 else "hourly"
    lines = []
    for _, row in table.iterrows():
        basin = row["Station_id"]
        category = row["subset_category"].lower()
        if not re.fullmatch(r"[A-Z0-9_-]+", basin):
            raise ValueError("Invalid basin identifier in Globus metadata")
        if category not in ("headwater", "meso-scale", "macro-scale"):
            raise ValueError("Unknown CAMELS-SPAT subset category")
        product = "daymet" if resolution == 1 else "rdrs"
        lines.extend([
            (
                (
                    "forcing/"
                    f"{category}"
                    "/"
                    f"{product}"
                    "/"
                    f"{product}"
                    "-lumped/"
                    f"{basin}"
                    "_"
                    f"{product}"
                    "_lumped.nc"
                ),
                f"{stream}/forcing/{product}/{basin}_{product}_lumped.nc",
            ),
            (
                (
                    "observations/"
                    f"{category}"
                    "/obs-"
                    f"{stream}"
                    "/"
                    f"{basin}"
                    "_"
                    f"{stream}"
                    "_flow_observations.nc"
                ),
                f"{stream}/discharge/{basin}_{stream}_flow_observations.nc",
            ),
        ])
    transfer(lines, "basin-data")
    return stage
