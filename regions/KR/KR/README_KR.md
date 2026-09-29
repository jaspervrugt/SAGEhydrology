# CAMELS-KR daily integration

Official release: Lee, S., Jung, G., and Ahn, K.-H. (2026), CAMELS-KR,
version 1.1, https://doi.org/10.5281/zenodo.21930882 (CC BY 4.0).
Paper: https://doi.org/10.5194/essd-2026-544 (ESSD preprint).
Archive CAMELS-KR.zip: 359530860 bytes; MD5 5ffce545fc0b3b79d57382a33dc597ca.

Select **South Korea (CAMELS, daily)** in the GUI and use **Install daily data**.
The installer verifies the archive checksum and the 282-gauge inventory.
Installed root: Data/CAMELS_KR; native CSV filenames are retained under
daily/forcing and daily/discharge. Seven native attribute tables sit at root.
The existing CAMELSH_KR / regions/KR hourly dataset and KR alias are preserved.
The internal module and basin-file prefix for the daily dataset is KR.

Gauge IDs (seven digits) identify files and joins; basin_id is not a gauge ID.
Precipitation uses prec (mm/day); temperature uses temp_avg (degrees C).
PET defaults to pet (FAO Penman-Monteith, mm/day), with pet_gleam selectable.
aet_gleam is actual evaporation and is not used as PET.
Observed discharge_spec is supplied in mm/day and is not area-normalized again.
discharge_vol (m3/s) is available in the source for independent checking.
Blank cells and -999/-9999 become NaN; negative fluxes become invalid, zero
remains valid, and negative temperatures are retained. No simulated discharge
is installed or substituted into observed gaps.

Forcing spans 1981-01-01 through 2025-12-31. Observed coverage varies by gauge.
Default scoring periods are October 2005–September 2015 for training and
October 2015–September 2025 for evaluation, with 365 days of spin-up.
All 282 source gauges enter standard check_basins/filter_basins screening;
the usable count depends on forcing completeness and observed Q coverage.
Location area is km2; basin_area is renamed area for SAGE metadata discovery.
Flow-derived static attributes are available but excluded from default predictors.
Source attributes otherwise retain their values, including source inconsistencies.

Validation: in MATLAB, add this region's tests folder to the path and run
`test_camels_kr(sourceRoot)`, where sourceRoot is the extracted CAMELS-KR
folder from the official archive. This installs into a temporary directory,
checks registry isolation, 282-gauge inventory, schemas, attributes, real
reader output, PET selection, screening, and synthetic missing-value cases.

Release audit: all 282 forcing and observed files have 16,436 unique daily
rows. Selected KMA precipitation, mean temperature and FAO PET are complete.
There are 2,344,192 blank observed discharge entries and 2,350 negative
entries (including converted -999 volumetric-discharge sentinels).
GLEAM PET has 60 negative entries, which remain invalid for screening.
The supplied specific discharge is used as published; recalculating from
rounded basin areas does not reproduce it exactly for every gauge.
