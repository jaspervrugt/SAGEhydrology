function [dat,aux] = read_meteo(region,dirM,bas,split,meteo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_METEO Read regional CAMELS meteorological data selected watersheds.
%
%  Select the meteorological schema for REGION and read basin forcing data.
%  The schema maps regional files, variables, and units to the common SAGE
%  data structure; it can also supply observed discharge or snow water
%  equivalent when those variables are available in the source dataset.
%
% SYNOPSIS:
%   [dat,aux] = read_meteo(region,dirM,bas,split,meteo)
%
% INPUT ARGUMENTS:
%   region          CAMELS region identifier:
%                    'CAMELS_AR'  = Argentina
%                    'CAMELS_AT'  = Austria
%                    'CAMELS_AU'  = Australia
%                    'CAMELS_BE'  = Belgium
%                    'CAMELS_BR'  = Brazil
%                    'CAMELS_CA'  = Canada
%                    'CAMELS_CH'  = Switzerland
%                    'CAMELS_CL'  = Chile
%                    'CAMELS_COL' = Colombia
%                    'CAMELS_CZ'  = Czechia
%                    'CAMELS_DE'  = Germany, CAMELS-DE daily
%                    'CAMELS_DEH' = Germany, CAMELS-DE-1h hourly
%                    'CAMELS_DK'  = Denmark
%                    'CAMELS_EE'  = Estonia
%                    'CAMELS_ES'  = Spain
%                    'BULL_ES'    = Spain (BULL)
%                    'CAMELS_FI'  = Finland
%                    'CAMELS_FR'  = France
%                    'CAMELS_GB'  = Great Britain
%                    'CAMELS_IE'  = Ireland
%                    'CAMELS_IND' = India
%                    'CAMELS_IS'  = Iceland
%                    'CAMELS_JM'  = Jamaica
%                    'CAMELSH_KR' = South Korea (hourly)
%                    'CAMELS_KR'  = South Korea (daily)
%                    'CAMELS_LUX' = Luxembourg
%                    'CAMELS_MX'  = Mexico
%                    'CAMELS_NA'  = Namibia
%                    'CAMELS_NO'  = Norway
%                    'CAMELS_NZ'  = New Zealand
%                    'CAMELS_PE'  = Peru
%                    'CAMELS_PL'  = Poland
%                    'CAMELS_PR'  = Puerto Rico
%                    'CAMELS_SE'  = Sweden
%                    'CAMELS_US'  = United States
%                    'CAMELS_USH' = United States
%                    'CAMELS_ZA'  = South Africa
%   dirM            directory containing meteorological data
%   bas             structure with selected-basin information
%    .id_gauge       K-by-1 gauge/catchment identifiers
%   split           prepared training/evaluation time split
%   meteo           meteorological reading options
%    .data           selected forcing product
%    .pet            selected potential evaporation product
%    .progressFcn    optional progress callback
%
% OUTPUT ARGUMENTS:
%   dat{k}          basin data with meteorology and available observations
%    .meteo          meteorological forcing structure
%     .P              precipitation
%     .Ep             potential evaporation
%     .T              air temperature
%     .LAI            leaf-area index, when supplied by the data source
%     .bad            indices with invalid meteorological data
%    .obs            observed-variable structure
%     .SWE            snow water equivalent, when available
%      .value          observed SWE (mm); empty when unavailable
%      .bad            invalid-SWE mask
%      .units          SWE units (mm)
%      .source         source dataset description
%     .Q              discharge, when included in the meteo schema
%    .gauge          gauge/catchment identifier
%    .fname          source file names
%   aux             K-by-3 basin metadata matrix:
%                   column 1 = latitude (degrees)
%                   column 2 = elevation (m)
%                   column 3 = basin area (m^2)
%
% NOTES:
%   Region-specific file names, formats, and units are defined by the
%   meteorological schema selected by REGION. Fields marked "when available"
%   depend on the selected regional data source.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    schema = meteo_schema(region);
    [dat,aux] = read_meteo_data( ...
        dirM,bas,split,meteo,schema);

end
