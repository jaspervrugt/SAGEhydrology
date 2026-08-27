function [dat,aux] = read_meteo(region,dirM,bas,split,meteo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_METEO Read regional CAMELS meteorological data for each watershed
%
% SYNOPSIS: [dat,aux] = read_meteo(region,dirM,bas,split,meteo)
%   region      string with CAMELS region
%                'CAMELS_AR'  = Argentina
%                'CAMELS_AT'  = Austria
%                'CAMELS_AU'  = Australia
%                'CAMELS_BE'  = Belgium
%                'CAMELS_BR'  = Brazil
%                'CAMELS_CA'  = Canada
%                'CAMELS_CH'  = Switzerland
%                'CAMELS_CL'  = Chile
%                'CAMELS_COL' = Colombia
%                'CAMELS_CZ'  = Czechia
%                'CAMELS_DE'  = Germany
%                'CAMELS_DK'  = Denmark
%                'CAMELS_EE'  = Estonia
%                'CAMELS_ES'  = Spain
%                'BULL_ES'    = Spain (BULL)
%                'CAMELS_FI'  = Finland
%                'CAMELS_FR'  = France
%                'CAMELS_GB'  = Great Britain
%                'CAMELS_IE'  = Ireland
%                'CAMELS_IND' = India
%                'CAMELS_IS'  = Iceland
%                'CAMELS_JM'  = Jamaica 
%                'CAMELS_KR'  = South Korea 
%                'CAMELS_LUX' = Luxembourg
%                'CAMELS_MX'  = Mexico
%                'CAMELS_NA'  = Namibia 
%                'CAMELS_NO'  = Norway 
%                'CAMELS_NZ'  = New Zealand
%                'CAMELS_PE'  = Peru
%                'CAMELS_PL'  = Poland
%                'CAMELS_PR'  = Puerto Rico
%                'CAMELS_SE'  = Sweden
%                'CAMELS_US'  = United States
%                'CAMELS_USH' = United States
%                'CAMELS_ZA'  = South Africa
%   dirM        string with main directory of meteorological data
%   bas         structure basin information
%    .K          number of selected watersheds
%    .id_gauge   Kx1 vector gauge/catchment identifiers
%   split       structure with training/evaluation time split
%    .dt         temporal resolution [1 = daily, 24 = hourly]
%    .idx        indices of scored period after spin-up
%    .local      optional flag for basin-local split
%   meteo       structure basin-average meteorological data
%    .data       choice of meteorological forcing product
%    .pet        choice of potential evaporation product
%    .progressFcn optional GUI progress callback
%   dat         OUTPUT: Kx1 cell structure with meteorological information
%    {k}.meteo.P   precipitation
%    {k}.meteo.Ep  potential evaporation
%    {k}.meteo.T   air temperature
%    {k}.meteo.bad indices with invalid meteorological data
%    {k}.gauge     gauge/catchment identifier
%    {k}.fname     source file names
%   aux         OUTPUT: Kx3 matrix with basin metadata
%                column 1 = latitude  (degrees)
%                column 2 = elevation (m)
%                column 3 = basin area (m^2)
%
% DESCRIPTION:
%   This function dispatches meteorological-data reading to the appropriate
%   region-specific CAMELS reader. Each regional reader converts its native
%   file names, column names, units, and data layout to the common SAGE
%   structure used by the hydrologic models.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    useGeneric = true;
    
    if useGeneric
        schema = meteo_schema(region);
        [dat,aux] = read_meteo_data( ...
            dirM,bas,split,meteo,schema);
        return
    else
        suffix = region_helpers('short',region); %#ok
        reader = str2func(['read_meteo_' suffix]);
        [dat,aux] = reader(dirM,bas,split,meteo);
    end

end
