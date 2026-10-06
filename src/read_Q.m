function dat = read_Q(region,dirQ,mdl,dat,bas,split,aux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_Q Read regional CAMELS discharge data for selected watersheds.
%
%  Select the discharge schema for REGION and read observed discharge for
%  the selected watersheds. The schema maps regional files, variables, and
%  units to the common SAGE data structure.
%
% SYNOPSIS:
%   dat = read_Q(region,dirQ,mdl,dat,bas,split,aux)
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
%                    'HYD_RESPONSES' = Switzerland HYD-RESPONSES
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
%   dirQ            directory containing discharge data
%   mdl             model settings; passed to the discharge reader
%   dat             K-by-1 cell array of existing basin data structures
%   bas             structure with selected-basin information
%    .K              number of selected watersheds
%    .id_gauge       K-by-1 gauge/catchment identifiers
%    .data_actions   optional actions for handling source data
%    .progressFcn    optional progress callback
%   split           prepared training/evaluation time split
%   aux             K-by-3 basin metadata matrix:
%                   column 1 = latitude (degrees)
%                   column 2 = elevation (m)
%                   column 3 = basin area (m^2)
%
% OUTPUT ARGUMENTS:
%   dat{k}          input basin data augmented with observed discharge
%    .obs            observed-variable structure
%     .Q              discharge structure
%      .value           observed discharge in canonical model units
%      .bad             invalid-discharge mask
%      .units           canonical discharge units
%      .source          source dataset description
%    .fname          source file names
%
% NOTES:
%   Region-specific file names, formats, and units are defined by the
%   discharge schema selected by REGION. Existing canonical discharge may
%   be reused when the schema permits it.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    schema = discharge_schema(region);
    dat = read_discharge_data( ...
        dirQ,mdl,dat,bas,split,aux,schema);

end
