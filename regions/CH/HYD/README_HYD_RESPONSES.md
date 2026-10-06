# HYD-RESPONSES integration

SAGE displays this dataset as **HYD-RESPONSES** under the country
**Switzerland** and installs it under `Data/HYD_RESPONSES`. The underscore is
used for the internal directory and region code; the published hyphenated name
is retained in all user-facing text.

The official v3 archive contains daily catchment-average hydro-meteorological
time series and measured FOEN streamflow for 184 Swiss catchments. SAGE uses
MeteoSwiss RhiresD precipitation, MeteoSwiss TabsD temperature, ERA5-Land
potential evaporation, and measured FOEN discharge by default. Alternative
ERA5-Land precipitation and temperature products, and the MeteoSwiss daily
minimum/maximum temperature average, remain selectable.

The published base-variable metadata labels discharge `Q` globally as m³/s.
Seven small experimental catchments (`2206`, `2251`, `2252`, `2282`, `2283`,
`2414`, and `2437`) contain discharge values in L/s. SAGE applies the required
factor of `1e-3` to those seven series before converting discharge to mm/day;
the downloaded source files are not modified.

SAGE normalizes the published numbered archive folders to the same regional
layout used by the other datasets:

- `daily/timeseries` contains the 184 model-input files;
- `daily/metadata` and the other named `daily` subfolders retain additional
  published daily products;
- `climatology` retains the published climatological products;
- the six descriptor CSV files and `gauge_information.txt` sit directly in
  `Data/HYD_RESPONSES`;
- `shapefiles` contains the published catchment outlets and outlines.

All numeric descriptors in the six published descriptor tables are
represented in the SAGE attribute catalog. Physical catchment descriptors are
eligible default predictors. Streamflow-record properties, Q347, baseflow,
delayed-flow, runoff-ratio, flashiness, Pardé, and other discharge-derived
signatures remain available for explicit selection but are intentionally
excluded from the default predictor set. The published categorical dominant
delayed-flow class is shown but disabled because SAGE model inputs must be
numeric; its four numeric delayed-flow contributions remain selectable.

Data and documentation:

- von Matt, C., Stocker, B., and Martius, O. (2026), *HYD-RESPONSES: daily
  hydro-meteorological catchment-level time series to analyse hydrological
  drought dynamics in response to cumulative water deficits in Swiss
  catchments*, Earth System Science Data, 18, 4113–4145.
  <https://doi.org/10.5194/essd-18-4113-2026>
- HYD-RESPONSES v3 archive: <https://doi.org/10.5281/zenodo.14713274>
- License: Creative Commons Attribution 4.0 International (CC BY 4.0).
