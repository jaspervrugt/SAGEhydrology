# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/cartplus_inputs.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
CAMELS-US inputs for the standalone CaRT+ CPU kernel.

Band grouping follows Xilin Xia's MIT-licensed CAMELS-US loader at revision
83682732728646669a190835d95ab3e60a8b5cb5.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np


SOURCES = {
    "CAMELS_GB": ("daily/CAMELS_GB_topographic_attributes.csv",
                  "CAMELS_GB_topographic_attributes.csv"),
    "CAMELS_AU": ('CAMELS_AUS_Attributes&Indices_' 'MasterTable.csv',),
    "CAMELS_ET": ("CAMELS_ET_topographic_attributes.csv",),
    "CAMELS_FI": ("CAMELS_FI_topographic_attributes.csv",),
    "CAMELS_DE": ("daily/CAMELS_DE_topographic_attributes.csv",
                  "CAMELS_DE_topographic_attributes.csv"),
    "CAMELS_PL": ("CAMELS_PL_topographic_attributes.csv",),
    "CHAMELLION": ("CHAMELLION_topographic_attributes.csv",),
    "CAMELS_FR": ('CAMELS_FR_topography_quantiles' '_attributes.csv',),
}
SUBDAILY_REGIONS = frozenset(("CAMELS_DEH", "CAMELSH_KR", "CAMELSH_US"))
REGIONS = ("CAMELS_US", *SOURCES)
SOURCES.update({
    "CAMELS_DEH": ('hourly/CAMELS_DE_1h_topographi' 'c_attributes.csv',
        "CAMELS_DE_1h_topographic_attributes.csv"),
    "CAMELSH_US": ("attributes/attributes_gageii_Topo.csv",
        "attributes_gageii_Topo.csv"),
})


@lru_cache(maxsize=16)
def _read_source(path, modified):
    """
    Read the installed elevation attribute table.

    SYNOPSIS:
        _read_source(...)

    INPUT ARGUMENTS:
    path : Input table, library, or output-file path.
    modified : Source modification time used as a cache key.

    OUTPUT ARGUMENTS:
    List of attribute dictionaries parsed from the source table.
    """
    import csv

    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter=";" if "quantiles" in
            path.name else ","))


def _regional_band_source(data_root, region):
    """
    Locate the installed regional elevation-band source.

    SYNOPSIS:
        _regional_band_source(...)

    INPUT ARGUMENTS:
    data_root : Root directory of the installed hydrological data.
    region : Dataset/region identifier used by the regional readers.

    OUTPUT ARGUMENTS:
    Path of the regional hypsometry or elevation-attribute source.
    """
    root = Path(data_root)
    path = next((root / name for name in SOURCES[region] if (root / name).
        is_file()), None)
    if path is None:
        raise FileNotFoundError(f'{region}'
            ' elevation attributes are missing un' 'der ' f'{root}' '.')
    return path


def regional_bands(data_root, region):
    """
    Read or generate all regional five-band elevation records.

    SYNOPSIS:
        regional_bands(...)

    INPUT ARGUMENTS:
    data_root : Root directory of the installed hydrological data.
    region : Dataset/region identifier used by the regional readers.

    OUTPUT ARGUMENTS:
    Mapping from basin identifiers to five elevation-band heights (m).
    """
    path = _regional_band_source(data_root, region)
    return _generate_bands(path, path.stat().st_mtime_ns, region)


def regional_band_count(data_root, region):
    """
    Count basins with usable regional elevation information.

    SYNOPSIS:
        regional_band_count(...)

    INPUT ARGUMENTS:
    data_root : Root directory of the installed hydrological data.
    region : Dataset/region identifier used by the regional readers.

    OUTPUT ARGUMENTS:
    Number of basin records with usable elevation information.
    """
    if region == "CAMELSH_US":
        path = _regional_band_source(data_root, region)
        return sum(
            elevation_summary(region, row) is not None
            for row in _read_source(path, path.stat().st_mtime_ns)
        )
    return len(regional_bands(data_root, region))


@lru_cache(maxsize=8)
def _generate_bands(path, modified, region):
    """
    Generate five elevation bands from regional attributes.

    SYNOPSIS:
        _generate_bands(...)

    INPUT ARGUMENTS:
    path : Input table, library, or output-file path.
    modified : Source modification time used as a cache key.
    region : Dataset/region identifier used by the regional readers.

    OUTPUT ARGUMENTS:
    Mapping from basin identifiers to five elevation-band heights (m).
    """
    rows = _read_source(path, modified)
    centres = (np.arange(5) + 0.5) / 5
    bands = {}
    if region == "CAMELS_FR":
        profiles = {}
        for row in rows:
            try:
                percentile = int(row["top_qnt_quant"])
                elevation = float(row["top_qnt_altitude"])
            except (ValueError, TypeError):
                continue
            if np.isfinite(elevation):
                profiles.setdefault(row["sta_code_h3"].strip(), {})[
                    percentile] = elevation
        for gauge, profile in profiles.items():
            if all(q in profile for q in (10, 30, 50, 70, 90)):
                bands[gauge] = np.array([profile[q] for q in (10, 30, 50,
                    70, 90)])
    elif region == "CAMELSH_US":
        for row in rows:
            heights = bands_from_attributes(region, row)
            if heights is not None:
                bands[row["STAID"].strip().zfill(8)] = heights
    elif region != "CAMELS_AU":
        lower = 5 if region in ("CAMELS_DE", "CAMELS_DEH", "CAMELS_PL",
            "CHAMELLION") else 10
        percentiles = [0, lower, 50, 100 - lower, 100]
        columns = ("elev_min", f"elev_{lower}", "elev_50",
                   f"elev_{100 - lower}", "elev_max")
        for row in rows:
            try:
                q = np.array([float(row[key]) for key in columns])
            except (ValueError, TypeError):
                continue
            if np.isfinite(q).all():
                bands[row["gauge_id"].strip()] = np.interp(
                    100 * centres, percentiles, np.maximum.accumulate(q)
                )
    else:
        from scipy.optimize import brentq
        from scipy.stats import truncnorm

        for row in rows:
            lo, mean, hi = (float(row[key]) for key in ("elev_min",
                "elev_mean", "elev_max"))
            if not np.isfinite([lo, mean, hi]).all() or hi <= mean:
                continue
            if mean <= lo:
                lo = mean - .5 * (hi - mean)
            span = hi - lo
            scale = .17 * span
            target = np.clip(mean, lo + .02 * span, hi - .02 * span)

            def distribution(loc):
                """
                Construct the truncated-normal elevation approximation.

                SYNOPSIS:
                    distribution(...)

                INPUT ARGUMENTS:
                loc : Location parameter of the elevation approximation (m).

                OUTPUT ARGUMENTS:
                Truncated-normal elevation distribution for the given
                location.
                """
                return truncnorm((lo - loc) / scale, (hi - loc) / scale,
                    loc=loc, scale=scale)

            loc = brentq(lambda loc: distribution(loc).mean() - target,
                         lo - 20 * span, hi + 20 * span)
            bands[row["station_id"].strip()] = distribution(loc).ppf(
                centres)
    return bands


# Regions with basin elevation summaries in the regional attribute tables.
ATTRIBUTE_COLUMNS = {
    **dict.fromkeys(
        ("CAMELS_AR", "CAMELS_BE", "CAMELS_CZ", "CAMELS_EE", "BULL_ES",
         "CAMELS_ES", "CAMELS_IE", "CAMELS_IL", "CAMELS_JM", "CAMELS_JP",
         "CAMELS_MX", "CAMELS_NA", "CAMELS_NO", "CAMELS_PR", "CAMELS_ZA",
         "CAMELSH_KR"),
        ("ele_mt_smn", "ele_mt_sav", "ele_mt_smx"),
    ),
    **dict.fromkeys(
        ("CAMELS_CA", "CAMELS_CH", "CAMELS_CL", "CAMELS_DK", "CAMELS_IND",
         "CAMELS_KR", "CAMELS_SPAT", "CAMELS_PE", "MACH_US"),
        ("elev_min", "elev_mean", "elev_max"),
    ),
    "CAMELS_LUX": ("Z_MIN", "Z_MEAN", "Z_MAX"),
    "EStreams_PT": ("ele_mt_min", "ele_mt_mean", "ele_mt_max"),
    "CAMELSH_US": ("ELEV_MIN_M_BASIN", "ELEV_MEAN_M_BASIN",
        "ELEV_MAX_M_BASIN"),
}
SPREAD_REGIONS = frozenset(("CAMELS_IS", "LAMAH_AT", "LAMAH_CE"))
ATTRIBUTE_REGIONS = (frozenset(ATTRIBUTE_COLUMNS) - {"CAMELSH_US"}) | \
    SPREAD_REGIONS
BAND_CENTRES = np.array([10, 30, 50, 70, 90])


def elevation_summary(region, row):
    """
    Extract available elevation summary statistics.

    SYNOPSIS:
        elevation_summary(...)

    INPUT ARGUMENTS:
    region : Dataset/region identifier used by the regional readers.
    row : Basin attribute record from the regional source.

    OUTPUT ARGUMENTS:
    Available elevation summary tuple (m), or None when missing.
    """
    if region in SPREAD_REGIONS:
        try:
            mean, spread, span = (float(row[k]) for k in ("elev_mean",
                "elev_std", "elev_ran"))
        except (KeyError, ValueError, TypeError):
            return None
        if np.isfinite([mean, spread, span]).all() and 0 <= spread <= span \
            / 2:
            return mean, spread
        return None
    try:
        lo, mean, hi = (float(row[key]) for key in ATTRIBUTE_COLUMNS[
            region])
    except (KeyError, ValueError, TypeError):
        return None
    if np.isfinite([lo, mean, hi]).all() and lo <= mean <= hi:
        return lo, mean, hi
    return None


def bands_from_attributes(region, row, *, with_method=False):
    """
    Reconstruct five elevations from published attributes.

    SYNOPSIS:
        bands_from_attributes(...)

    INPUT ARGUMENTS:
    region : Dataset/region identifier used by the regional readers.
    row : Basin attribute record from the regional source.
    with_method : Whether to return a description of the band-generation
        method.

    OUTPUT ARGUMENTS:
    Five elevations (m), or None when unavailable; with_method also returns
    a method description.
    """
    def result(heights, method):
        """
        Return reconstructed elevations and optional method metadata.

        SYNOPSIS:
            result(...)

        INPUT ARGUMENTS:
        heights : Five reconstructed elevation-band heights (m).
        method : Selected derivative or solver method.

        OUTPUT ARGUMENTS:
        (heights, method) if with_method else heights : computed result in
        the operand units.
        """
        return (heights, method) if with_method else heights

    if region in ("CAMELS_CH", "CAMELS_KR"):
        if region == "CAMELS_CH":
            percentiles = [0, 10, 25, 50, 75, 90, 100]
            keys = ["elev_min", *(f"elev_percentile{q}" for q in
                percentiles[1:-1]), "elev_max"]
        else:
            percentiles = [0, 5, 95, 100]
            keys = ["elev_min", "elev_5", "elev_95", "elev_max"]
        try:
            elevations = np.array([float(row[key]) for key in keys])
        except (KeyError, ValueError, TypeError):
            elevations = np.array([])
        if elevations.size and np.isfinite(elevations).all() and (np.diff(
            elevations) >= 0).all():
            return result(np.interp(BAND_CENTRES, percentiles, elevations),
                          'Interpolation of published bas'
                              'in elevation percentiles')
    summary = elevation_summary(region, row)
    if summary is None:
        return result(None, "")
    if region in SPREAD_REGIONS:
        from scipy.special import ndtri

        mean, spread = summary
        return result(mean + spread * ndtri(BAND_CENTRES / 100),
                      'Normal approximation from basi'
                          'n mean and elevation standard ' 'deviation')
    lo, mean, hi = summary
    if lo == hi:
        return result(np.full(5, lo),
            "Uniform elevation from basin summaries")
    from scipy.optimize import brentq
    from scipy.stats import truncnorm

    span = hi - lo
    scale = .17 * span
    target = np.clip(mean, lo + .02 * span, hi - .02 * span)

    def distribution(loc):
        """
        Construct the truncated-normal elevation approximation.

        SYNOPSIS:
            distribution(...)

        INPUT ARGUMENTS:
        loc : Location parameter of the elevation approximation (m).

        OUTPUT ARGUMENTS:
        Truncated-normal elevation distribution for the given location.
        """
        return truncnorm((lo - loc) / scale, (hi - loc) / scale, loc=loc,
            scale=scale)

    loc = brentq(lambda loc: distribution(loc).mean() - target,
                 lo - 20 * span, hi + 20 * span)
    return result(distribution(loc).ppf(BAND_CENTRES / 100),
                  'Truncated-normal approximation'
                      ' from basin minimum, mean and ' 'maximum elevations')


TERRAIN_REGIONS = frozenset((
    "HYDRO_CIS", "HYD_RESPONSES", "CAMELS_BR", "CAMELS_COL", "CAMELS_NZ",
        "CAMELS_SE",
))


def terrain_bands(data_root):
    """
    Load terrain-based elevations from the installed dataset.

    SYNOPSIS:
        terrain_bands(...)

    INPUT ARGUMENTS:
    data_root : Root directory of the installed hydrological data.

    OUTPUT ARGUMENTS:
    Mapping of installed basin identifiers to terrain elevation data.
    """
    path = Path(data_root) / "cartplus_terrain_bands.csv"
    return _terrain_bands(path, path.stat().st_mtime_ns)


@lru_cache(maxsize=8)
def _terrain_bands(path, modified):
    """
    Read cached terrain-attribute elevation records.

    SYNOPSIS:
        _terrain_bands(...)

    INPUT ARGUMENTS:
    path : Input table, library, or output-file path.
    modified : Source modification time used as a cache key.

    OUTPUT ARGUMENTS:
    Mapping of installed basin identifiers to terrain elevation data.
    """
    rows = _read_source(path, modified)
    groups = {}
    for row in rows:
        groups.setdefault(row["basin_id"], []).append(row)
    result = {}
    for gauge, items in groups.items():
        items.sort(key=lambda row: int(row["band"]))
        heights = np.array([float(row["elevation_m"]) for row in items])
        fractions = np.array([float(row["area_fraction"]) for row in items])
        if ([int(row["band"]) for row in items] != [1, 2, 3, 4, 5]
                or not np.isfinite(heights).all() or not np.isfinite(
                    fractions).all()
                or (fractions < 0).any() or not np.isclose(fractions.sum()
                    , 1)
                or (np.diff(heights) < 0).any()
                or not items[0]["source_method"]):
            raise ValueError(f"Invalid terrain bands for basin {gauge}.")
        result[gauge] = (heights, fractions, items[0]["source_method"])
    return result


def band_status(data_root, region, resolution):
    """
    Describe the elevation-band source available for a dataset.

    SYNOPSIS:
        band_status(...)

    INPUT ARGUMENTS:
    data_root : Root directory of the installed hydrological data.
    region : Dataset/region identifier used by the regional readers.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.

    OUTPUT ARGUMENTS:
    Human-readable availability or data-installation status.
    """
    if resolution != 1:
        return "CaRT+ uses daily forcing."
    if region not in REGIONS:
        return 'No CaRT+ band source configure' 'd for this dataset.'
    root = Path(data_root)
    if region == "CAMELS_US":
        return ('Available from installed CAMEL' 'S-US hypsometry'
                if (root / "elev_bands671.npz").is_file()
                else 'Missing CAMELS-US hypsometry (' 'elev_bands671.npz)')
    if any((root / name).is_file() for name in SOURCES[region]):
        method = "min, mean and max elevation" if region == "CAMELS_AU" \
            else "elevation quantiles"
        return 'Generated from installed ' f'{method}' ' when CaRT+ runs'
    return 'Elevation attributes required.' ' Install the regional dataset.'


def camels_us_band_data(data_root):
    """
    Load basin identifiers and CAMELS-US elevation fractions.

    SYNOPSIS:
        camels_us_band_data(...)

    INPUT ARGUMENTS:
    data_root : Root directory of the installed hydrological data.

    OUTPUT ARGUMENTS:
    Basin identifiers, five elevations (m), and band area fractions.
    """
    path = Path(data_root) / "elev_bands671.npz"
    stat = path.stat()
    return _camels_us_band_data(path, stat.st_mtime_ns, stat.st_size)


@lru_cache(maxsize=8)
def _camels_us_band_data(path, modified, size):
    """
    Read published CAMELS-US hypsometry from its source.

    SYNOPSIS:
        _camels_us_band_data(...)

    INPUT ARGUMENTS:
    path : Input table, library, or output-file path.
    modified : Source modification time used as a cache key.
    size : Requested state/array size.

    OUTPUT ARGUMENTS:
    Basin identifiers, five elevations (m), and band area fractions.
    """
    with np.load(path, allow_pickle=False) as data:
        ids = np.char.zfill(data["ids"].astype(str), 8)
        heights = data["band_elev"]
        fractions = data["band_frac"]
    if (ids.ndim != 1 or not len(ids) or heights.ndim != 2
            or heights.shape != fractions.shape or len(heights) != len(ids
                )):
        raise ValueError(f"Invalid CAMELS-US elevation archive: {path}.")
    return ids, heights, fractions


def elevation_bands(data_root, gauge, region="CAMELS_US"):
    """
    Read or reconstruct five basin elevation bands.

    SYNOPSIS:
        elevation_bands(...)

    INPUT ARGUMENTS:
    data_root : Root directory of the installed hydrological data.
    gauge : Basin identifier in the selected regional dataset.
    region : Dataset/region identifier used by the regional readers.

    OUTPUT ARGUMENTS:
    Five elevations (m) and normalized nonnegative area fractions.
    """
    if region in TERRAIN_REGIONS:
        try:
            heights, fractions, _ = terrain_bands(data_root)[str(gauge).
                strip()]
        except KeyError:
            raise ValueError('No terrain elevation bands for ' f'{region}'
                ' basin ' f'{gauge}' '.') from None
        return heights.copy(), fractions.copy()
    if region != "CAMELS_US":
        bands = regional_bands(data_root, region)
        try:
            return bands[str(gauge).strip()].copy(), np.full(5, .2)
        except KeyError:
            raise ValueError('No usable CaRT+ elevation attributes'
                ' for ' f'{region}' ' basin ' f'{gauge}' '.') from None
    ids, all_heights, all_fractions = camels_us_band_data(data_root)
    matches = np.flatnonzero(ids == str(gauge).zfill(8))
    if len(matches) != 1:
        raise ValueError('CaRT and CaRT+ need one elevation re'
            'cord for ' f'{gauge}' '.')
    heights = all_heights[matches[0]]
    fractions = all_fractions[matches[0]]
    positive = fractions > 0
    heights, fractions = heights[positive], fractions[positive]
    if (not len(heights) or not np.isfinite(heights).all()
            or not np.isfinite(fractions).all()):
        raise ValueError(f"Invalid CaRT+ elevation bands for {gauge}.")
    order = np.argsort(heights)
    heights, fractions = heights[order], fractions[order]
    groups = np.minimum((np.cumsum(fractions) - 1e-9) // 0.2, 4).astype(int)
    area = np.bincount(groups, weights=fractions, minlength=5)
    weighted = np.bincount(groups, weights=fractions * heights, minlength=5)
    elevation = np.divide(weighted, area, out=np.zeros(5), where=area > 0)
    return elevation, area / area.sum()


def attach_elevation_bands(dat, data_root, gauge, region="CAMELS_US"):
    """
    Attach validated elevation bands to basin data in place.

    SYNOPSIS:
        attach_elevation_bands(...)

    INPUT ARGUMENTS:
    dat : Basin data structure; modified to attach elevation bands.
    data_root : Root directory of the installed hydrological data.
    gauge : Basin identifier in the selected regional dataset.
    region : Dataset/region identifier used by the regional readers.

    OUTPUT ARGUMENTS:
    None. Attaches validated elevation/area arrays to dat in place.
    """
    elevation, area = elevation_bands(data_root, gauge, region)
    dat["cartplus_elevation_bands"] = {
        "elevation_m": elevation,
        "area_fraction": area,
    }


# Source group, upstream name, and multiplier for each native column.
# A None multiplier indicates an identity transformation.
PARAMETER_MAP = (
    (0, "TTs", None),
    (0, "lCF", 1.0),
    (0, "SFCF", None),
    (0, "lks", 1.0),
    (0, "lrfa", 1.0),
    (0, "lGam", 1.0),
    (0, "lCu", 1.0),
    (0, "ltu", 1.0),
    (0, "lC", 1.0),
    (0, "rKp", (4.0 / 9.0) ** 4),
    (0, "lbP", 1.0),
    (0, "lPs", 1.0),
    (0, "gie", None),
    (0, "wz", None),
    (0, "lCz", 1.0),
    (0, "ltz", 1.0),
    (0, "lCd", 1.0),
    (0, "lbd", 1.0),
    (0, "lkd", 1.0),
    (0, "ltau", 1.0),
    (0, "lamE", None),
    (0, "kv", None),
    (2, "Tv", None),
    (2, "lsv", 1.0),
    (2, "ltv", 1.0),
    (1, "lCt", 1.0),
    (1, "ltt", 1.0),
    (1, "Amin", None),
    (2, "thc", None),
    (2, "lsth", 1.0),
    (1, "lams", None),
    (1, "lamW", None),
    (0, "gR", None),
    (0, "X4", None),
    (0, "lX3", 1.0),
    (0, "X2", None),
    (0, "lJin", 1.0),
    (2, "TF", None),
    (2, "lsF", 1.0),
    (1, "aF", None),
    (1, "k", None),
    (0, "lLc", 1.0),
)


def physical_parameters(static, law, shared):
    """
    Map CaRT+ calibration groups to the 42 physical parameters.

    SYNOPSIS:
        physical_parameters(...)

    INPUT ARGUMENTS:
    static : Mapping from static parameter names to calibration values.
    law : Mapping from law-parameter names to calibration values.
    shared : CaRT+ shared/fixed physical parameter values.

    OUTPUT ARGUMENTS:
    Length-42 physical parameter vector in PARAMETER_MAP order.
    """
    groups = static, law, shared
    return np.array([
        groups[group][key] if multiplier is None else multiplier * np.exp(
            groups[group][key])
        for group, key, multiplier in PARAMETER_MAP
    ], dtype=np.float64)
