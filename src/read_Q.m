function dat = read_Q(region,dirQ,mdl,dat,bas,split,aux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_Q Read regional CAMELS discharge data for each watershed
%
% SYNOPSIS: dat = read_Q(region,dirQ,mdl,dat,bas,split,aux)
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
%   dirQ        string with main directory of discharge data
%   mdl         structure with model state/parameter information
%    .mode       assessment design
%    .id_train   training-period indices or basin-local masks
%    .id_eval    evaluation-period indices or basin-local masks
%   dat         Kx1 cell structure with meteorological/model information
%   bas         structure basin information
%    .K          number of selected watersheds
%    .K_t        number of training watersheds
%    .K_e        number of evaluation watersheds
%    .id_gauge   Kx1 vector gauge/catchment identifiers
%    .progressFcn optional GUI progress callback
%   split       structure with training/evaluation time split
%    .dt         temporal resolution [1 = daily, 24 = hourly]
%    .idx        indices of scored period after spin-up
%   aux         Kx3 matrix with basin metadata
%                column 1 = latitude  (degrees)
%                column 2 = elevation (m)
%                column 3 = basin area (m^2)
%   dat         OUTPUT: input dat structure with discharge info added
%    {k}.y_n      observed discharge transformed to model units
%    {k}.bad      indices with invalid discharge data
%    {k}.use      'training' or 'evaluation'
%    {k}.fname    source file names
%
% DESCRIPTION:
%   This function dispatches discharge-data reading to the appropriate
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
        schema = discharge_schema(region);
        dat = read_discharge_data( ...
            dirQ,mdl,dat,bas,split,aux,schema);
    else
        suffix = region_helpers('short',region);    %#ok
        fname = ['read_Q_' suffix];    
        if exist(fname,'file') ~= 2
            error('read_Q:unknownRegion', ...
                'Missing discharge reader: %s.m',fname);
        end    
        reader = str2func(fname);
        dat = reader(dirQ,mdl,dat,bas,split,aux);
    end

end
