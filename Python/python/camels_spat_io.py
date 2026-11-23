# ==========================================================================
# SAGEhydrology | Python source
# File: python/camels_spat_io.py
#
# Use published paths directly, retaining nonstandard-layout discovery.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Use published paths directly, retaining nonstandard-layout discovery."""

from pathlib import Path


def inputs(root, region, gauge, resolution=1):
    """Read the selected forcing inputs and reference time window."""
    root = Path(root)
    full = "CAN_" + str(gauge) if region == "CAMELS_CA" else str(gauge)
    stream, product, suffix = (
        ("daily", "daymet", "daymet_lumped.nc")
        if resolution == 1
        else ("hourly", "rdrs", "rdrs_lumped.nc")
    )
    name = full + "_" + suffix
    forcing = root / stream / "forcing" / product / name
    if not forcing.is_file():
        candidates = list((root / stream / "forcing").glob("**/" + name))
        if len(candidates) != 1:
            raise FileNotFoundError(
                f"Expected one {name}; found {len(candidates)}"
            )
        forcing = candidates[0]
    discharge = (
        root
        / stream
        / "discharge"
        / (full + "_" + stream + "_flow_observations.nc")
    )
    if not discharge.is_file():
        raise FileNotFoundError(discharge)
    return forcing, discharge
