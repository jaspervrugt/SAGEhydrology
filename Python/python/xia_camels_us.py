"""Install and read the published CaRT+ daily CAMELS-US dataset."""

from functools import lru_cache
from pathlib import Path

from shutil import copyfileobj
from tempfile import NamedTemporaryFile
from zipfile import ZipFile

import numpy as np


FORCING = "camels_us_daymet671_v4full.npz"
BANDS = "elev_bands671.npz"
FILES = {
    FORCING: "84482d304b59aad8717e97eb7a90784bb08a25256ec45b1d6dac4de223139a54",
    BANDS: "f980712c6a64107cb63e50eced3b3d520fe0805890fdf66cadb12a46c532e0ae",
}
GB_FORCING = "camels_gb_cart.npz"
GB_FILES = {GB_FORCING: "454187814591e3e57f9d10d79b7fabe997f9b7aa30fab4bb958502b2c5206238"}
GB_SERIES = ("P", "PET", "Tm", "q")
SERIES = ("P", "PET", "Tmax", "Tmin", "q")
RELEASE = "https://github.com/xiaxilin/cartplus/releases/download/v1.0/"


def directory(software_root, variant="us"):
    """Return the project folder, retaining legacy US installations."""
    root = Path(software_root) / "Data" / "Projects" / "6"
    if variant == "us" and ready(root):
        return root
    return root / variant.upper()


def project_region(root):
    """Identify the dataset from its reserved Project 6 path."""
    parts = tuple(p.lower() for p in Path(root).parts)
    if parts[-3:] == ("data", "projects", "6"):
        return "CAMELS_US"
    if parts[-4:-1] == ("data", "projects", "6"):
        return {"us": "CAMELS_US", "gb": "CAMELS_GB"}.get(parts[-1])
    return None


def matches(root):
    """Identify a Project 6 US data path."""
    return project_region(root) == "CAMELS_US"


def standard_root(root, region):
    """Locate the regular regional dataset for optional display metadata."""
    root = Path(root)
    depth = 2 if root.name == "6" else 3
    return root.parents[depth - 1] / region


def ready(root):
    """Check that the installer has completed both downloads."""
    root = Path(root)
    gb = project_region(root) == "CAMELS_GB"
    return (root / ".xia-v1.0").is_file() and all(
        (root / name).is_file() for name in (GB_FILES if gb else FILES)
    ) and all((root / f"{key}.npy").is_file() for key in (GB_SERIES if gb else SERIES))


def install(root, progress=lambda *args: None, stop=None,
            stage=lambda message: None):
    """Download and verify the two release assets."""
    from dataset_install import _check, download

    root = Path(root)
    if ready(root):
        stage("Project 6 files already installed.")
        return
    gb = project_region(root) == "CAMELS_GB"
    files, series = (GB_FILES, GB_SERIES) if gb else (FILES, SERIES)
    for name, checksum in files.items():
        stage("Downloading and verifying " + name)
        download(
            dict(name=name, url=RELEASE + name, sha256=checksum),
            root, progress, stop,
        )
    # Memory-map the arrays so basin workers share the OS file cache.
    stage("Preparing daily forcing and discharge arrays")
    with ZipFile(root / (GB_FORCING if gb else FORCING)) as archive:
        for index, key in enumerate(series, 1):
            _check(stop)
            name = f"{key}.npy"
            progress("Preparing " + name, index - 1, len(series))
            with archive.open(name) as source, NamedTemporaryFile(
                dir=root, delete=False
            ) as target:
                temp = Path(target.name)
                try:
                    copyfileobj(source, target)
                    target.close()
                    temp.replace(root / name)
                finally:
                    temp.unlink(missing_ok=True)
    _check(stop)
    (root / ".xia-v1.0").touch()


def _source(root):
    if not ready(root):
        raise ValueError("Use Data Install on the Project 6 page first.")
    gb = project_region(root) == "CAMELS_GB"
    path = Path(root) / (GB_FORCING if gb else FORCING)
    stat = path.stat()
    return path, stat.st_size, stat.st_mtime_ns


@lru_cache(maxsize=1)
def _metadata(path, size, mtime_ns):
    with np.load(path) as data:
        gb = path.name == GB_FORCING
        ids = [str(g) if gb else str(g).zfill(8) for g in data["ids"]]
        return (
            ids, data["days"].astype("datetime64[D]"),
            list(data["static_names"]), data["static" if gb else "static_raw"]
        )


@lru_cache(maxsize=1)
def _series(path, size, mtime_ns):
    return tuple(
        np.load(path.parent / f"{key}.npy", mmap_mode="r")
        for key in (GB_SERIES if path.name == GB_FORCING else SERIES)
    )


def catalogue(root):
    """Return the gauge IDs in the release."""
    ids, _, _, _ = _metadata(*_source(root))
    return [(g, "", "") for g in ids]


def attribute_names(root):
    """Return the attribute order supplied in the archive."""
    return list(map(str, _metadata(*_source(root))[2]))


def attributes(root):
    """Return the publisher's unnormalized static attributes."""
    ids, _, names, values = _metadata(*_source(root))
    return {
        gauge: dict(zip(names, row))
        for gauge, row in zip(ids, values.tolist())
    }


def read(root, gauge, *, train_start, train_end, eval_start, eval_end,
         spinup, resolution=1, allow_invalid=False, meteo_only=False,
         **options):
    """Read published P, PET, temperature and discharge for one basin."""
    if resolution != 1:
        raise ValueError("The Project 6 dataset has daily resolution.")
    source = _source(root)
    ids, days, names, static = _metadata(*source)
    gb = source[0].name == GB_FORCING
    gauge = str(gauge) if gb else str(gauge).zfill(8)
    col = ids.index(gauge)
    ts, te, es, ee = (
        np.datetime64(d, "D") for d in
        (train_start, train_end, eval_start, eval_end)
    )
    start = min(ts, es) - np.timedelta64(int(spinup), "D")
    end = max(te, ee)
    if (start < days[0] and not gb) or end > days[-1]:
        raise ValueError("Requested dates exceed the Project 6 record.")
    select = slice(days.searchsorted(start), days.searchsorted(end, "right"))
    day = days[select]
    values = [arr[select, col].astype(float) for arr in _series(*source)]
    if gb:
        p, pet, temperature, q = values
        forcing = np.column_stack((p, pet, temperature - 273.15))
        if start < days[0]:
            warm_days, warm_forcing = _gb_warmup(*source)
            take = (warm_days >= start) & (warm_days < days[0])
            prefix = warm_days[take]
            if not np.array_equal(prefix, np.arange(start, days[0])):
                raise ValueError("Requested warm-up exceeds the Project 6 GB record.")
            warm = warm_forcing[take, col][:, [0, 2, 1]].astype(float)
            warm[:, 2] -= 273.15
            forcing = np.concatenate((warm, forcing))
            q = np.concatenate((np.full(len(prefix), np.nan), q))
            day = np.concatenate((prefix, day))
    else:
        p, pet, tmax, tmin, q = values
        forcing = np.column_stack((p, pet, 0.5 * (tmax + tmin)))
    if meteo_only:
        q[:] = np.nan
    valid = np.isfinite(q) & (q >= 0)
    q[~valid] = np.nan
    train = (day >= ts) & (day <= te) & (valid | meteo_only)
    evaluation = (day >= es) & (day <= ee) & (valid | meteo_only)
    if not allow_invalid and not np.isfinite(forcing).all():
        raise ValueError("Missing or invalid Project 6 forcing.")
    path = str(source[0])
    return dict(
        gauge=gauge, dates=day.astype("datetime64[s]"), forcing=forcing,
        observed=q, matlab_forcing=forcing.astype(np.float32).astype(float),
        train=train, evaluation=evaluation,
        area=None if gb else float(static[col, names.index("area_gages2")]),
        forcing_file=path, streamflow_file=path, _source_files=[path],
    )

@lru_cache(maxsize=1)
def _gb_elevations(path, size, mtime_ns):
    with np.load(path) as data:
        return dict(zip(map(str, data["ids"]), data["elev_quantiles"]))


def gb_bands(root, gauge):
    """Return Xia's five equal-area elevation bands for a British basin."""
    values = _gb_elevations(*_source(root))[str(gauge)]
    return np.interp(
        [10, 30, 50, 70, 90], [0, 10, 50, 90, 100],
        np.maximum.accumulate(values),
    )


@lru_cache(maxsize=1)
def _gb_warmup(path, size, mtime_ns):
    with np.load(path) as data:
        return data["spin_days"].astype("datetime64[D]"), data["spin_forcing"]
