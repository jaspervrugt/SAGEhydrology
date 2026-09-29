# HydroCIS Russia daily integration

Official release: Abramov, D. and Kurochkina, L. (2023),
*Hydro-meteorological database for watersheds across the CIS*, version 1.3,
https://doi.org/10.5281/zenodo.8432070 (CC BY 4.0).
Archive `Russia_HydroMeteo_Database_v04.zip`: 3,504,912,123 bytes;
MD5 `8d536a7b76d8fdd6004e321db03f4220`.

Select **Russia** in the GUI and use **Download daily**. The Summary tab
identifies the selection as **Russia (HydroCIS)**; its separate Resolution
line reports that the data are daily.
The installer verifies the archive checksum, joins the published NetCDF,
attribute, gauge, and catchment inventories by `gauge_id`, and installs the
1,886 Russian gauges retained by the source data-quality screen. The 278 GRDC
gauges elsewhere in the CIS that were added to version 1.3 are not mixed into
the Russia region. The archive is retained in Downloads after installation.

Installed root: `Data/HYDRO_CIS`. Combined forcing and observed-discharge
NetCDF records are stored under `daily/timeseries`; `static_data.csv`, derived
UTF-8 `gauge_information.txt`, and both published GeoPackages are stored at
the data root. The 1,886-basin inventory belongs to the regional code at
`regions/RU/RU_1886_basins.txt`. Daily records span 2008-01-01 through
2020-12-31.

Gauge names are displayed as cleaned English transliterations followed by
the official Russian name in parentheses. The English labels expand or
remove the abbreviated Russian feature and settlement prefixes used in the
source files; the original Cyrillic text remains unchanged in
`gauge_information.txt`.

Because the published static table does not contain climate summaries,
installation derives `climate_attributes.csv` once from the common
2008--2020 daily record. ERA5-Land is preferred, with ERA5 and MSWEP used
where available as fallbacks; Oudin potential evaporation provides the PET
summary. The resulting mean temperature, annual precipitation and PET,
aridity, and precipitation-weighted snow fraction support SAGE's
hydroclimatic basin-zone classification. Individual basins lacking the
required source series retain the classifier's documented fallback.

Selectable precipitation products are ERA5-Land, ERA5, MSWEP, GPCP, and
IMERG. Selectable temperature products are the daily Tmin/Tmax means from
ERA5-Land and ERA5. Potential evaporation defaults to the Oudin formulation
derived from the selected temperature and gauge latitude; published GLEAM
potential evaporation is the second choice. GLEAM actual evaporation and its
component fluxes are not offered as potential-evaporation inputs. Observed
specific discharge is read directly from `q_mm_day`.

The Attributes tab exposes 24 catchment predictors. The first 19 are complete
HydroATLAS fields supplied by the official HydroCIS release and remain the
default set; they cover topography, land cover, lakes and reservoirs, soil,
and karst. Five optional climate summaries derived during installation are
also selectable: long-term mean air temperature, mean annual precipitation,
mean annual Oudin potential evaporation, aridity, and precipitation-weighted
snow fraction. Mean annual precipitation is complete for all 1,886 Russian
gauges. The other four climate summaries are unavailable for gauges 1105,
2111, 2155, 2192, and 2241 because neither ERA5-Land nor ERA5 temperature is
present in those records. The separate gauge-elevation field is retained as
station metadata in `gauge_information.txt`; it is not offered as a catchment
predictor because nine Russian gauges lack a value. Native source field names
remain unchanged for reproducible joins; only the GUI labels are expanded.
When any of the four affected climate predictors is selected, SAGE excludes
the five incomplete gauges from the eligible sampling pool and reports their
identifiers; it does not impute climate values. Selecting the original 19
attributes, or adding mean annual precipitation alone, retains all 1,886
Russian gauges.
