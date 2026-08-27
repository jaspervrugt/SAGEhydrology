# CAMELSH-KR integration

Dataset: CAMELSH-KR, Zenodo record 15073264.

- Resolution: hourly
- Basin count: 178
- Default modeling universe: 165 basins (140 training and 25 evaluation).
  All 178 source files remain installed. Eight basins have insufficient
  discharge coverage in at least one default period: 1018623, 1018680,
  1021650, 1022640, 1023670, 2001658, 2004625, and 2012625. Five basins
  contain severe source-data defects during a default period: 1002625,
  1003620, 1022670, 2002655, and 2012690. These exclusions are recorded
  explicitly rather than corrected silently by the reader. Full findings
  are in `Data/CAMELSH_KR/CAMELSH_KR_discharge_audit.md`; the excluded
  IDs are also listed in `KR_13_excluded_basins.txt`.
- Data period: 2000-2019
- Installed time-series directory: `Data/CAMELSH_KR/hourly/timeseries`
- Each basin CSV contains meteorological forcing, observed streamflow, and
  water level. The meteorological and discharge directories are therefore
  the same.
- Basin identifiers are preserved as seven-character strings.

`read_attr_KR.m`, `read_meteo_KR.m`, and `read_Q_KR.m` are implemented.
The hourly reader offers ERA5-Land or observed precipitation and air
temperature, uses the supplied signed ERA5-Land potential evaporation,
retains snow cover/depth as auxiliary forcing, and converts volumetric
streamflow to catchment-average runoff depth.

The distributed `snow_depth` values are retained without a unit-changing
conversion because version 1 supplies no variable dictionary. They are not
renamed to SWE.
