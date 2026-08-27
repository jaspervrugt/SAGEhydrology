function out = data_helpers(action,varargin)
%DATA_HELPERS Region-aware data utilities for SAGE GUI
%
% Usage from SAGE_ui:
%   out = data_helpers('supports',region,'daily');
%   out = data_helpers('supports',region,'hourly');
%   out = data_helpers('supports',region,'15min');
%   out = data_helpers('status',cfg);
%   out = data_helpers('download',cfg,'daily',ui);
%   out = data_helpers('download',cfg,'hourly',ui);
%   out = data_helpers('download',cfg,'15min',ui);
%   out = data_helpers('mapstatus',cfg);
%   out = data_helpers('ensuremap',cfg,ui);
%
% The optional ui structure for 'download' may contain:
%   .fig          parent uifigure for progress dialog
%   .logFcn       function handle for logging text
%   .refreshFcn   function handle called after download
%   .confirmFcn   function handle: @(title,msg) true/false
%   .progressTitle custom title string
%
% This helper keeps region-specific data logic out of SAGE_ui.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 1 ...
            || isempty(action)
        error('data_helpers:missingAction', ...
            'Action must be provided.');
    end
    
    action = lower(strtrim(char(action)));
    
    switch action
        case 'supports'
            out = local_supports(varargin{:});
    
        case 'status'
            out = local_status(varargin{:});
    
        case 'hasdaily'
            out = local_has_stream( ...
                varargin{1},'daily');
    
        case 'hashourly'
            out = local_has_stream( ...
                varargin{1},'hourly');

        case {'has15min','has15minute'}
            out = local_has_stream( ...
                varargin{1},'15min');
    
        case {'hasmetadata','hasgauge', ...
                'hasattributes'}
            out = local_has_metadata( ...
                varargin{:});
    
        case 'download'
            out = local_download( ...
                varargin{:});
    
        case 'regionname'
            out = local_region_name( ...
                varargin{:});
    
        case 'datadir'
            out = local_data_dir( ...
                varargin{:});
    
        case 'streamdir'
            out = local_stream_dir( ...
                varargin{:});
    

        case 'ensuremap'
            out = local_ensure_natural_earth_map( ...
                varargin{:});

        case {'mapfile','mapstatus'}
            out = local_natural_earth_map_status( ...
                varargin{:});

        otherwise
            error('data_helpers:badAction', ...
                ['Unknown data_helpers ' ...
                'action: %s'],action);
    end
end

% =====================
% Public action helpers
% =====================
function tf = local_supports(region,stream)
    region = local_region_code(region);
    stream = lower(strtrim(char(stream)));
    
    switch region
        case {'CAMELS_AT', ...
                'CAMELS_US', ...
                'CAMELS_DE', ...
                'CAMELS_GB', ...
                'CAMELS_IS'}
            tf = any(strcmp(stream,{'daily','hourly'}));
        case 'CAMELSH_US'
            tf = strcmp(stream,'hourly');
        case 'CAMELS_CZ'
            tf = any(strcmp(stream,{'daily','hourly'}));
        case 'CAMELS_CA'
            tf = any(strcmp(stream,{'daily','hourly'}));
        case {'CAMELS_AU', ...
                'CAMELS_BR', ...
                'CAMELS_CL', ...
                'CAMELS_CH', ...
                'CAMELS_COL', ...
                'CAMELS_DK', ...
                'CAMELS_ES', ...
                'BULL_ES', ...
                'CAMELS_FR', ...
                'CAMELS_FI', ...
                'CAMELS_IND', ...
                'CAMELS_IL', ...
                'CAMELS_JP', ...
                'MACH_US', ...
                'CAMELS_MX', ...
                'CAMELS_PE', ...
                'CAMELS_PL', ...
                'CAMELS_SE'}
            tf = strcmp(stream,'daily');
        case 'CAMELS_NZ'
            tf = any(strcmp(stream, ...
                {'daily','hourly'}));
        case 'CAMELSH_KR'
            tf = strcmp(stream,'hourly');
        case {'CAMELS_ZA','CAMELS_NA','CAMELS_AR','CAMELS_BE', ...
                'CAMELS_EE','CAMELS_IE','CAMELS_JM','CAMELS_NO','CAMELS_PR'}
            tf = strcmp(stream,'daily');
        case 'CAMELS_LUX'
            tf = any(strcmp(stream, ...
                {'daily','hourly','15min'}));
        otherwise
            tf = false;
    end
end

function S = local_status(cfg)
    region = local_cfg_region(cfg);
    
    S = struct();
    S.region = region;
    S.regionName = local_region_name( ...
        region);
    S.supportsDaily = local_supports( ...
        region,'daily');
    S.supportsHourly = local_supports( ...
        region,'hourly');
    S.supports15min = local_supports( ...
        region,'15min');
    S.hasDaily = local_has_stream( ...
        cfg,'daily');
    S.hasHourly = local_has_stream( ...
        cfg,'hourly');
    S.has15min = local_has_stream( ...
        cfg,'15min');
    S.hasMetadata = local_has_metadata(cfg);
    
    S.dailyText = local_present_text( ...
        S.hasDaily);
    S.hourlyText = local_present_text( ...
        S.hasHourly);
    S.min15Text = local_present_text( ...
        S.has15min);
    S.metadataText = local_present_text( ...
        S.hasMetadata);
end

function ok = local_has_stream(cfg,stream)
    region = local_cfg_region(cfg);
    stream = lower(strtrim(char(stream)));

    if ~local_supports(region,stream)
        ok = false;
        return
    end

    dirD = local_cfg_dirD(cfg,region);

    switch region

        case 'CAMELS_AU' %1
            dDaily = fullfile(dirD,'daily');
            ok = isfile(fullfile(dDaily, ...
                'streamflow','streamflow_mmd.csv')) ...
                && isfile(fullfile(dDaily, ...
                'precipitation','precipitation_SILO.csv')) ...
                && isfile(fullfile(dDaily, ...
                'precipitation','precipitation_AWAP.csv')) ...
                && local_count_files(fullfile(dDaily, ...
                'evaporative_demand'),{'*.csv'}) >= 8 ...
                && isfile(fullfile(dDaily, ...
                'temperature','tmin_SILO.csv')) ...
                && isfile(fullfile(dDaily, ...
                'temperature','tmax_SILO.csv')) ...
                && isfile(fullfile(dDaily, ...
                'temperature','tmin_AWAP.csv')) ...
                && isfile(fullfile(dDaily, ...
                'temperature','tmax_AWAP.csv'));

        case 'CAMELS_AT' %2
            if strcmp(stream,'daily')
                dF = fullfile(dirD,'daily', ...
                    'timeseries','forcing');
                dQ = fullfile(dirD,'daily', ...
                    'timeseries','streamflow');
            else
                dF = fullfile(dirD,'hourly', ...
                    'timeseries','forcing');
                dQ = fullfile(dirD,'hourly', ...
                    'timeseries','streamflow');
            end        
            ok = isfolder(dF) ...
                && isfolder(dQ) ...
                && local_count_files(dF,{'ID_*.csv'}) >= 859 ...
                && local_count_files(dQ,{'ID_*.csv'}) >= 859;

        case 'CAMELS_BR' %3
            ok = local_br_component_ok(dirD, ...
                'precipitation',897) ...
                && local_br_component_ok(dirD, ...
                'pet',897) ...
                && local_br_component_ok(dirD, ...
                'temperature',897) ...
                && local_br_component_ok(dirD, ...
                'streamflow',897);

        case 'CAMELS_CA' %4
            nExpected = local_ca_expected_basin_count(dirD);
            if strcmp(stream,'daily')
                dF = fullfile(dirD,'daily','forcing','daymet');
                dQ = fullfile(dirD,'daily','discharge');
                nF = local_count_files(dF, ...
                    {'CAN_*_daymet_lumped.nc'});
                nQ = local_count_files(dQ, ...
                    {'CAN_*_daily_flow_observations.nc'});
            elseif strcmp(stream,'hourly')
                dF = fullfile(dirD,'hourly','forcing','rdrs');
                dQ = fullfile(dirD,'hourly','discharge');
                nF = local_count_files(dF, ...
                    {'CAN_*_rdrs_lumped.nc'});
                nQ = local_count_files(dQ, ...
                    {'CAN_*_hourly_flow_observations.nc'});
            else
                nF = 0;
                nQ = 0;
            end
            if isfinite(nExpected) ...
                    && nExpected > 0
                ok = nF >= nExpected ...
                    && nQ >= nExpected;
            else
                ok = nF > 0 ...
                    && nQ > 0;
            end

        case 'CAMELS_CH' %5
            dObs = fullfile(dirD,'daily','timeseries', ...
                'observation_based');
            dSim = fullfile(dirD,'daily','timeseries', ...
                'simulation_based');
            ok = strcmp(stream,'daily') ...
                && isfolder(dObs) ...
                && isfolder(dSim) ...
                && local_count_files(dObs, ...
                {'CAMELS_CH_obs_based_*.csv'}) > 0 ...
                && local_count_files(dSim, ...
                {'CAMELS_CH_sim_based_*.csv'}) > 0;

        case 'CAMELS_CL' %6
            dDaily = fullfile(dirD,'daily');
            ok = isfile(fullfile(dDaily, ...
                'streamflow','3_CAMELScl_streamflow_mm.txt')) ...
                && isfile(fullfile(dDaily, ...
                'precipitation','4_CAMELScl_precip_cr2met.txt')) ...
                && isfile(fullfile(dDaily, ...
                'precipitation','5_CAMELScl_precip_chirps.txt')) ...
                && isfile(fullfile(dDaily, ...
                'precipitation','6_CAMELScl_precip_mswep.txt')) ...
                && isfile(fullfile(dDaily, ...
                'precipitation','7_CAMELScl_precip_tmpa.txt')) ...
                && isfile(fullfile(dDaily, ...
                'temperature','8_CAMELScl_tmin_cr2met.txt')) ...
                && isfile(fullfile(dDaily, ...
                'temperature','9_CAMELScl_tmax_cr2met.txt')) ...
                && isfile(fullfile(dDaily, ...
                'temperature','10_CAMELScl_tmean_cr2met.txt')) ...
                && isfile(fullfile(dDaily, ...
                'pet','11_CAMELScl_pet_8d_modis.txt')) ...
                && isfile(fullfile(dDaily, ...
                'pet','12_CAMELScl_pet_hargreaves.txt')) ...
                && isfile(fullfile(dDaily, ...
                'swe','13_CAMELScl_swe.txt'));

        case 'CAMELS_COL' %7
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d, ...
                {'Hydromet_data_*.txt'}) > 0;

        case 'CAMELS_CZ' %8
            d = local_cz_timeseries_install_dir(dirD,stream);
            ok = any(strcmp(stream,{'daily','hourly'})) ...
                && ~isempty(d) ...
                && local_count_files(d,{'camelscz_*.csv'}) >= 249;

        case 'CAMELS_DE' %9
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily','timeseries');
                ok = isfolder(d) ...
                    && local_count_files(d, ...
                    {'CAMELS_DE_hydromet_timeseries_DE*.csv'}) >= 1484;
            elseif strcmp(stream,'hourly')
                d = fullfile(dirD,'hourly','timeseries');
                ok = isfolder(d) ...
                    && local_count_files(d, ...
                    {'CAMELS_DE_1h_hydromet_timeseries_*.csv'}) >= 1256;
            else
                ok = false;
            end

        case 'CAMELS_DK' %10
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d, ...
                {'CAMELS_DK_obs_based_*.csv'}) >= 304;

        case 'CAMELS_ES' %11
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d, ...
                {'camelses_*.csv'}) > 0;

        case 'BULL_ES'
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && local_count_files(fullfile(d,'streamflow'), ...
                {'streamflow_*.csv'}) >= 484 ...
                && local_count_files(fullfile(d,'AEMET'), ...
                {'AEMET_*.csv'}) >= 484 ...
                && local_count_files(fullfile(d,'ERA5_Land'), ...
                {'ERA5_Land_*.csv'}) >= 484 ...
                && local_count_files(fullfile(d,'EMO1_arc'), ...
                {'EMO1_*.csv'}) >= 484;

        case 'CAMELS_FI' %12
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d, ...
                {'CAMELS_FI_hydromet_timeseries_*.csv'}) > 0;

        case 'CAMELS_FR' %13
            d = fullfile(dirD,'daily','timeseries');
            ok = isfolder(d) ...
                && local_count_files(d, ...
                {'CAMELS_FR_tsd_*.csv'}) >= 50;

        case 'CAMELS_GB' %14
            dDaily = fullfile(dirD,'daily','timeseries');
            dHourly = fullfile(dirD,'hourly','timeseries');
            nDaily = local_count_files(dDaily,{'*.csv'});
            nHourly = local_count_files(dHourly,{'*.csv'});
            if strcmp(stream,'daily')
                ok = isfolder(dDaily) ...
                    && nDaily >= 669;
            elseif strcmp(stream,'hourly')
                % Hourly forcing uses daily PET and temperature, so both
                % hourly and daily time-series collections are required.
                ok = isfolder(dHourly) ...
                    && nHourly >= 664 ...
                    && isfolder(dDaily) ...
                    && nDaily >= 669;
            else
                ok = false;
            end

        case 'CAMELS_LUX' %15
            d = fullfile(dirD,stream,'timeseries');
            if strcmp(stream,'daily')
                pats = {['CAMELS_LUX_hydromet_' ...
                    'timeseries__daily_ID_*.csv']};
                nExpected = 50;
            elseif strcmp(stream,'hourly')
                pats = {['CAMELS_LUX_hydromet_' ...
                    'timeseries_hourly_ID_*.csv']};
                nExpected = 51;
            elseif strcmp(stream,'15min')
                % Keep this tolerant to small naming differences between
                % archive versions while requiring Luxembourg 15-minute CSVs.
                pats = {['CAMELS_LUX_hydromet_' ...
                    'timeseries_15min_ID_*.csv'], ...
                    ['CAMELS_LUX_hydromet_' ...
                    'timeseries__15min_ID_*.csv'], ...
                    'CAMELS_LUX*15min*ID_*.csv'};
                nExpected = 51;
            else
                pats = {};
                nExpected = inf;
            end
            ok = ~isempty(pats) ...
                && isfinite(nExpected) ...
                && isfolder(d) ...
                && local_count_files(d,pats) >= nExpected;

        case 'CAMELS_IL' %16
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d,{'il_*.csv'}) >= 95;

        case 'CAMELS_IND' %17
            d = fullfile(dirD,'daily', ...
                'catchment_mean_forcings');
            ok = isfile(fullfile(dirD,'daily', ...
                'streamflow_timeseries', ...
                'streamflow_observed.csv')) ...
                && isfolder(d) ...
                && local_count_files(d, ...
                {'*.csv'}) >= 242;

        case 'CAMELS_IS' %18
            dF = fullfile(dirD,stream,'forcing');
            dQ = fullfile(dirD,stream,'discharge');

            if strcmp(stream,'daily')
                nExpected = 111;
            elseif strcmp(stream,'hourly')
                nExpected = 76;
            else
                nExpected = inf;
            end

            ok = isfinite(nExpected) ...
                && isfolder(dF) ...
                && isfolder(dQ) ...
                && local_count_files( ...
                    dF,{'ID_*.csv'}) >= nExpected ...
                && local_count_files( ...
                    dQ,{'ID_*.csv'}) >= nExpected ...
                && local_is_csv_kind( ...
                    fullfile(dF,'ID_1.csv'), ...
                    'forcing') ...
                && local_is_csv_kind( ...
                    fullfile(dQ,'ID_1.csv'), ...
                    'discharge');

        case 'CAMELS_JP'
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d,{'varssim*.csv'}) >= 87;

        case 'CAMELS_MX'
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d,{'hysets_*.csv'}) >= 46;

        case 'MACH_US'
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d,{'basin_*_MACH.csv'}) >= 1014;

        case 'CAMELS_NZ' %19
            ok = local_has_nz_timeseries( ...
                dirD,stream);

        case 'CAMELSH_KR'
            d = fullfile(dirD,'hourly','timeseries');
            ok = strcmp(stream,'hourly') ...
                && isfolder(d) ...
                && local_count_files(d,{'*.csv'}) >= 178;

        case {'CAMELS_ZA','CAMELS_NA','CAMELS_AR','CAMELS_BE', ...
                'CAMELS_EE','CAMELS_IE','CAMELS_JM','CAMELS_NO','CAMELS_PR'}
            d=fullfile(dirD,'daily','timeseries');
            spec=local_grdc_africa_spec(extractAfter(region,'CAMELS_'));
            ok=strcmp(stream,'daily')&&isfolder(d) ...
                && local_count_files(d,{'*.csv'})>=spec.count;

        case 'CAMELS_PE' %20
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d,{'PE_*.csv'}) >= 136;

        case 'CAMELS_PL' %21
            d = fullfile(dirD,'daily','timeseries');
            ok = strcmp(stream,'daily') ...
                && isfolder(d) ...
                && local_count_files(d, ...
                {'CAMELS_PL_hydromet_timeseries_*.csv'}) >= 354;            

        case 'CAMELS_SE' %22
            d = fullfile(dirD,'daily');
            ok = isfolder(d) ...
                && local_count_files(d, ...
                {'catchment_id_*.csv'}) >= 50;

        case 'CAMELS_US' %23
            if strcmp(stream,'daily')
                d0 = fullfile(dirD,'daily');
                dF = fullfile(dirD,'daily', ...
                    'v1p2','forcing');
                dQ = fullfile(dirD,'daily', ...
                    'v1p2','streamflow');
                ok = isfolder(d0) && ...
                    (local_count_files(d0, ...
                    {'*.txt','*.csv','*.mat'}) > 0 ...
                    || local_count_files(dF, ...
                     {'*.txt','*.csv','*.mat'}) > 0 ...
                    || local_count_files(dQ, ...
                     {'*.txt','*.csv','*.mat'}) > 0);
            else
                d0 = fullfile(dirD,'hourly');
                dF = fullfile(dirD,'hourly', ...
                    'forcing');
                dQ = fullfile(dirD,'hourly', ...
                    'streamflow');
                ok = isfolder(d0) ...
                    && (local_count_files(d0, ...
                    {'*.txt','*.csv','*.mat'}) > 0 ...
                    || local_count_files(dF, ...
                    {'*.txt','*.csv','*.mat'}) > 0 ...
                    || local_count_files(dQ, ...
                     {'*.txt','*.csv','*.mat'}) > 0);
            end

        case 'CAMELSH_US'
            d=fullfile(dirD,'hourly','timeseries');
            ok=strcmp(stream,'hourly')&&isfolder(d) ...
                && local_count_files(d,{'*.nc'})>=3166;

        otherwise
            ok = false;
    end
end

function ok = local_has_metadata(cfg)
    region = local_cfg_region(cfg);
    dirD = local_cfg_dirD(cfg,region);

    switch region
        case 'CAMELS_AT' %1
            ok = isfile(fullfile(dirD, ...
                'Catchment_attributes.csv')) ...
                && isfile(fullfile(dirD, ...
                'Gauge_attributes.csv'));

        case 'CAMELS_AU' %2
            ok = isfile(fullfile(dirD, ...
                ['CAMELS_AUS_Attributes' ...
                '&Indices_MasterTable.csv'])) ...
                && (isfile(fullfile(dirD, ...
                'id_name_metadata.csv')) ...
                || isfile(fullfile(dirD, ...
                '01_id_name_metadata.csv')));

        case 'CAMELS_BR' %3
            ok = local_all_files_exist(dirD, ...
                local_br_metadata_files());

        case 'CAMELS_CA' %4
            ok = isfile(fullfile(dirD, ...
                'camels-spat-metadata.csv')) ...
                && isfile(fullfile(dirD, ...
                'attributes-lumped.csv'));

        case 'CAMELS_CH' %5
            ok = local_all_files_exist(dirD, ...
                local_ch_metadata_files());

        case 'CAMELS_CL' %6
            ok = isfile(fullfile(dirD, ...
                '1_CAMELScl_attributes.txt'));

        case 'CAMELS_COL' %7
            ok = local_all_files_exist(dirD, ...
                local_col_metadata_files());

        case 'CAMELS_CZ' %8
            ok = local_all_files_exist(dirD,{ ...
                'attributes_caravan_camelscz.csv', ...
                'attributes_hydroatlas_camelscz.csv', ...
                'attributes_other_camelscz.csv'});

        case 'CAMELS_DE' %9
            stream = local_cfg_stream(cfg);
            dirDE = fullfile(dirD,stream);
            ok = local_all_files_exist(dirDE, ...
                local_de_metadata_files(stream));

        case 'CAMELS_DK' %10
            ok = local_count_files(dirD,{'*.csv'}) >= 7;

        case 'CAMELS_ES' %11
            ok = local_all_files_exist(dirD, ...
                local_es_metadata_files());

        case 'BULL_ES'
            ok = local_all_files_exist(dirD,{ ...
                'attributes_other_.csv', ...
                'attributes_caravan_.csv', ...
                'attributes_hydroatlas_.csv'});

        case 'CAMELS_FI' %12
            ok = local_all_files_exist(dirD, ...
                local_fi_metadata_files());

        case 'CAMELS_FR' %13
            ok = local_all_files_exist(dirD, ...
                local_fr_metadata_files());

        case 'CAMELS_GB' %14
            stream = local_cfg_stream(cfg);
            dAttr = fullfile(dirD,stream);
            ok = local_count_files(dAttr,{'*.csv'}) >= 8 ...
                || isfile(fullfile(dAttr,'gauge_information.txt'));

        case 'CAMELS_IL' %15
            d = fullfile(dirD,'daily');
            ok = local_all_files_exist(d,{ ...
                'attributes_caravan_il.csv', ...
                'attributes_hydroatlas_il.csv', ...
                'attributes_other_il.csv'});

        case 'CAMELS_IND' %16
            ok = local_all_files_exist(dirD, ...
                local_ind_metadata_files());

        case 'CAMELS_IS' %17
            ok = isfile(fullfile(dirD, ...
                'Catchment_attributes.csv')) ...
                && isfile(fullfile(dirD, ...
                'Gauge_attributes.csv'));

        case 'CAMELS_LUX' %18
            ok = local_all_files_exist(dirD, ...
                local_lux_metadata_files());

        case 'CAMELS_NZ' %19
            ok = local_all_files_exist(dirD, ...
                local_nz_metadata_files());

        case 'CAMELSH_KR'
            ok = local_all_files_exist(dirD, ...
                local_kr_metadata_files());

        case 'CAMELS_JP'
            ok = isfile(fullfile(dirD, ...
                'MERV_Jp_135_HydroATLAS_attributes.xlsx'));

        case 'CAMELS_MX'
            ok = local_all_files_exist(dirD,{ ...
                'attributes_caravan_hysets.csv', ...
                'attributes_hydroatlas_hysets.csv', ...
                'attributes_other_hysets.csv'});

        case 'MACH_US'
            ok = local_all_files_exist(dirD,{ ...
                fullfile('attributes','site_info.csv'), ...
                fullfile('attributes','overall_climate.csv'), ...
                fullfile('attributes','soil.csv'), ...
                fullfile('attributes','geology.csv'), ...
                fullfile('attributes','hydrology.csv'), ...
                fullfile('attributes','anthropogenic.csv')});

        case {'CAMELS_ZA','CAMELS_NA','CAMELS_AR','CAMELS_BE', ...
                'CAMELS_EE','CAMELS_IE','CAMELS_JM','CAMELS_NO','CAMELS_PR'}
            ok=local_all_files_exist(dirD,{ ...
                'attributes_caravan_grdc.csv', ...
                'attributes_hydroatlas_grdc.csv', ...
                'attributes_other_grdc.csv', ...
                'attributes_additional_grdc.csv'});

        case 'CAMELS_PE' %20
            ok = local_all_files_exist(dirD,{ ...
                'stations.csv','data_dictionary.csv', ...
                'climatic_indices.csv','geologic_attributes.csv', ...
                'human_intervention_attributes.csv', ...
                'hydrological_signatures.csv','landcover_attributes.csv', ...
                'soil_attributes.csv','topographic_attributes.csv'});

        case 'CAMELS_PL' %21
            ok = local_all_files_exist(dirD, ...
                local_pl_metadata_files());

        case 'CAMELS_SE' %22
            ok = isfile(fullfile(dirD, ...
                'catchments_physical_properties.csv')) ...
                && isfile(fullfile(dirD, ...
                'catchments_landcover.csv')) ...
                && isfile(fullfile(dirD, ...
                'catchments_soil_classes.csv'));

        case 'CAMELS_US' %23
            % US normally has gauge_information.txt 
            % plus CAMELS attribute txt files.
            ok = isfile(fullfile(dirD, ...
                'gauge_information.txt')) ... 
                || local_count_files(dirD, ...
                {'*.txt'}) >= 5;
        case 'CAMELSH_US'
            ok=local_all_files_exist(dirD,{ ...
                'info.csv', ...
                fullfile('attributes','attributes_gageii_BasinID.csv'), ...
                fullfile('attributes','attributes_nldas2_climate.csv'), ...
                fullfile('attributes','attributes_gageii_Topo.csv')});

        otherwise
            ok = false;
    end
end

function out = local_download(cfg,stream,ui)
%LOCAL_DOWNLOAD Region-aware download/install dispatcher.
%
% US uses install_SAGEhydrology(root,'daily'/'hourly',opts).
% GB/BR hooks are included as dataset-specific placeholders so the GUI can
% route correctly without growing nested callbacks. Fill URLs/installers
% below once the final data sources are fixed.

    if nargin < 3 ...
            || isempty(ui)
        ui = struct();
    end
    
    region = local_cfg_region(cfg);
    stream = lower(strtrim(char(stream)));
    
    if ~local_supports(region,stream)
        error('data_helpers:unsupportedStream', ...
            '%s data are not available for %s.', ...
            upper(stream),local_region_name(region));
    end
    
    switch region

        case 'CAMELS_AT' %1
            out = local_download_at(cfg,stream,ui);

        case 'CAMELS_AU' %2
            out = local_download_au(cfg,stream,ui);
            
        case 'CAMELS_BR' %3
            out = local_download_br(cfg,stream,ui);
    
        case 'CAMELS_CA' %4
            out = local_download_ca(cfg,stream,ui);

        case 'CAMELS_CH' %5
            out = local_download_ch(cfg,stream,ui);

        case 'CAMELS_CL' %6
            out = local_download_cl(cfg,stream,ui);

        case 'CAMELS_COL' %7
            out = local_download_col(cfg,stream,ui);

        case 'CAMELS_CZ' %8
            out = local_download_cz(cfg,stream,ui);

        case 'CAMELS_DE' %9
            out = local_download_de(cfg,stream,ui);

        case 'CAMELS_DK' %10
            out = local_download_dk(cfg,stream,ui);

        case 'CAMELS_ES' %11
            out = local_download_es(cfg,stream,ui);

        case 'BULL_ES'
            out = local_download_bull(cfg,stream,ui);

        case 'CAMELS_FI' %12
            out = local_download_fi(cfg,stream,ui);

        case 'CAMELS_FR' %13
            out = local_download_fr(cfg,stream,ui);

        case 'CAMELS_GB' %14
            out = local_download_gb(cfg,stream,ui);

        case 'CAMELS_IL' %15
            out = local_download_il(cfg,stream,ui);
            
        case 'CAMELS_IND' %16
            out = local_download_ind(cfg,stream,ui);

        case 'CAMELS_IS' %17
            out = local_download_is(cfg,stream,ui);

        case 'CAMELS_LUX' %18
            out = local_download_lux(cfg,stream,ui);

        case 'CAMELS_JP'
            out = local_download_jp(cfg,stream,ui);

        case 'CAMELS_MX'
            out = local_download_mx(cfg,stream,ui);

        case 'MACH_US'
            out = local_download_mach(cfg,stream,ui);

        case 'CAMELS_NZ' %19
            out = local_download_nz(cfg,stream,ui);

        case 'CAMELSH_KR'
            out = local_download_kr(cfg,stream,ui);

        case 'CAMELS_ZA'
            out = local_download_za(cfg,stream,ui);

        case 'CAMELS_NA'
            out = local_download_na(cfg,stream,ui);
        case 'CAMELS_AR'
            out = local_download_grdc_africa(cfg,stream,ui,'AR');
        case 'CAMELS_BE'
            out = local_download_grdc_africa(cfg,stream,ui,'BE');
        case 'CAMELS_EE'
            out = local_download_grdc_africa(cfg,stream,ui,'EE');
        case 'CAMELS_IE'
            out = local_download_grdc_africa(cfg,stream,ui,'IE');
        case 'CAMELS_JM'
            out = local_download_grdc_africa(cfg,stream,ui,'JM');
        case 'CAMELS_NO'
            out = local_download_grdc_africa(cfg,stream,ui,'NO');
        case 'CAMELS_PR'
            out = local_download_grdc_africa(cfg,stream,ui,'PR');

        case 'CAMELS_PE' %20
            out = local_download_pe(cfg,stream,ui);

        case 'CAMELS_PL' %21
            out = local_download_pl(cfg,stream,ui);

        case 'CAMELS_SE' %22
            out = local_download_se(cfg,stream,ui);
            
        case 'CAMELS_US' %23
            out = local_download_us(cfg,stream,ui);   
        case 'CAMELSH_US'
            out = local_download_ush(cfg,stream,ui);
            
        otherwise
            error('data_helpers:unknownRegion', ...
                'Unknown data region: %s',region);
    end
end

function stream = local_cfg_stream(cfg)

    stream = 'daily';
    if isstruct(cfg) ...
            && isfield(cfg,'prd') ...
            && isfield(cfg.prd,'dt') ...
            && ~isempty(cfg.prd.dt)

        dt = double(cfg.prd.dt);

        if abs(dt - 96) < 1e-8 ...
                || abs(dt - 0.25) < 1e-8 ...
                || abs(dt - 15/1440) < 1e-8 ...
                || abs(dt - 900) < 1e-8
            stream = '15min';
        elseif abs(dt - 24) < 1e-8 ...
                || abs(dt - 1/24) < 1e-8 ...
                || abs(dt - 3600) < 1e-8
            stream = 'hourly';
        else
            stream = 'daily';
        end
    end
end

% ========================
% Download implementations
% ========================

% =============================================
function out = local_download_at(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_AT Download and install CAMELS-AT / LamaH-CE data.

    logFcn = local_ui_log(ui);

    if ~any(strcmp(stream,{'daily','hourly'}))
        error('data_helpers:ATBadStream', ...
            ['CAMELS-AT supports ' ...
            'stream = daily or hourly.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_AT');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,stream, ...
        'timeseries','forcing'));
    local_mkdir(fullfile(dirD,stream, ...
        'timeseries','streamflow'));

    msg = sprintf(['CAMELS-AT / LamaH-CE %s ' ...
        'is downloaded from Zenodo ' ...
        'as one tar.gz archive. The ' ...
        'installer copies catchment ' ...
        'attributes, gauge attributes, ' ...
        'basin forcing files, and ' ...
        ['gauge streamflow files into ' ...
        'Data/CAMELS_AT.'] newline ...
        'Continue?'],stream);

    if ~local_ui_confirm(ui, ...
        'Download Austria data',msg)
        out = struct('ok',false,'canceled',true);
        return
    end

    if strcmp(stream,'daily')
        archiveName = '2_LamaH_daily.tar';
        url = ['https://zenodo.org/records/5153305/files/' ...
            '2_LamaH-CE_daily.tar.gz?download=1'];
        expectedMD5 = '69fd2733e969513403f923ecc5eaa3dc';
    else
        archiveName = '1_LamaH-CE_daily_hourly.tar.gz';
        url = ['https://zenodo.org/records/5153305/files/' ...
            '1_LamaH-CE_daily_hourly.tar.gz?download=1'];
        expectedMD5 = '121b2292b288ef1c2ef4d96250be77ac';
    end
    
    downloadDir = local_default_download_dir();
    archiveFile = fullfile(downloadDir,archiveName);

    % downloadDir = local_default_download_dir();
    % gzName = '2_LamaH-CE_daily.tar.gz';
    % gzFile = fullfile(downloadDir,gzName);
    % url = ['https://zenodo.org/records/5153305/files/' ...
    %     '2_LamaH-CE_daily.tar.gz?download=1'];
    % expectedMD5 = '69fd2733e969513403f923ecc5eaa3dc';

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Austria data', ...
        'Starting CAMELS-AT / LamaH-CE download ...');

    try
        if local_at_install_complete(dirD,stream)
            logFcn(['CAMELS-AT appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true, ...
                'region','CAMELS_AT', ...
                'dirD',dirD, ...
                'dirM',fullfile(dirD,stream, ...
                'timeseries','forcing'), ...
                'dirQ',fullfile(dirD,stream, ...
                'timeseries','streamflow'), ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(archiveFile)
            try
                local_verify_md5(archiveFile, ...
                    expectedMD5,logFcn);
                logFcn(['Using verified ' ...
                    'existing file: ' archiveFile]);
            catch
                logFcn(['Existing CAMELS-AT archive ' ...
                    'failed MD5; deleting.']);
                delete(archiveFile);
            end
        end

        if ~isfile(archiveFile)
            logFcn(['Downloading CAMELS-AT ' ...
                'archive to: ' archiveFile]);
            local_download_file_retry(url,archiveFile, ...
                ['CAMELS-AT / LamaH-CE ' stream],d, ...
                @() userCanceled,@ui_progress,logFcn,3);
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        local_verify_md5(archiveFile, ...
            expectedMD5,logFcn);

        unzipDir = fullfile(downloadDir, ...
            'CAMELS_AT_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        try
            d.Indeterminate = true;
            d.Message = ['Untarring CAMELS-AT archive. ' ...
                'This can take a while ...'];
            drawnow limitrate nocallbacks
        catch
        end

        untar(archiveFile,unzipDir);

        local_install_at_files(unzipDir, ...
            dirD,logFcn,stream);

        ok = local_at_install_complete(dirD,stream);
        out = struct('ok',ok, ...
            'region','CAMELS_AT', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,stream, ...
                'timeseries','forcing'), ...
            'dirQ',fullfile(dirD,stream, ...
                'timeseries','streamflow'), ...
            'archive',archiveFile, ...
            'unzipDir',unzipDir);

        if ok
            nF = local_count_files(fullfile(dirD, ...
                stream,'timeseries','forcing'), ...
                {'ID_*.csv'});
            nQ = local_count_files(fullfile(dirD, ...
                stream,'timeseries','streamflow'), ...
                {'ID_*.csv'});
            logFcn(sprintf(['CAMELS-AT %s install ' ...
                'finished: %d forcing ' ...
                'files and %d streamflow files.'],stream, ...
                nF,nQ));
            local_cleanup_download_file(archiveFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-AT install ' ...
                'finished, but some required ' ...
                'files appear to be missing.']);
        end

    catch ME
        try
            if exist('archiveFile','var')
                local_cleanup_download_file(archiveFile, ...
                    logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder(unzipDir, ...
                    logFcn);
            end
        catch
        end
        logFcn('CAMELS-AT install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% ==================================================
function tf = local_at_install_complete(dirD,stream)
% ==================================================

    dF = fullfile(dirD,stream, ...
        'timeseries','forcing');
    dQ = fullfile(dirD,stream, ...
        'timeseries','streamflow');
    
    tf = isfile(fullfile(dirD, ...
        'Catchment_attributes.csv')) ...
        && isfile(fullfile(dirD, ...
        'Gauge_attributes.csv')) ...
        && isfolder(dF) ...
        && isfolder(dQ) ...
        && local_count_files(dF,{'ID_*.csv'}) >= 859 ...
        && local_count_files(dQ,{'ID_*.csv'}) >= 859;

end

% =========================================================
function local_install_at_files(srcRoot,dirD,logFcn,stream)
% =========================================================

    dataDir = local_find_child_dir(srcRoot, ...
        'CAMELS_AT');
    if isempty(dataDir)
        dataDir = local_find_child_dir(srcRoot, ...
            'CAMELS-AT');
    end
    if isempty(dataDir)
        dataDir = srcRoot;
    end

    dirA = fullfile(dataDir,'A_basins_total_upstrm');
    dirG = fullfile(dataDir,'D_gauges');
    if ~isfolder(dirA)
        error('data_helpers:ATMissingABasins', ...
            ['Could not find A_basins_total_upstrm ' ...
            'in CAMELS-AT archive.']);
    end
    if ~isfolder(dirG)
        error('data_helpers:ATMissingDGauges', ...
            ['Could not find D_gauges in ' ...
            'CAMELS-AT archive.']);
    end

    local_mkdir(dirD);

    dstF = fullfile(dirD,stream, ...
        'timeseries','forcing');
    dstQ = fullfile(dirD,stream, ...
        'timeseries','streamflow');
    local_mkdir(dstF);
    local_mkdir(dstQ);

    fCatch = fullfile(dirA,'1_attributes', ...
        'Catchment_attributes.csv');
    fGauge = fullfile(dirG,'1_attributes', ...
        'Gauge_attributes.csv');
    if ~isfile(fCatch)
        error('data_helpers:ATMissingCatchmentAttributes', ...
            ['Missing A_basins_total_upstrm/' ...
            '1_attributes/Catchment_attributes.csv.']);
    end
    if ~isfile(fGauge)
        error('data_helpers:ATMissingGaugeAttributes', ...
            ['Missing D_gauges/1_attributes/' ...
            'Gauge_attributes.csv.']);
    end

    copyfile(fCatch,fullfile(dirD, ...
        'Catchment_attributes.csv'),'f');
    copyfile(fGauge,fullfile(dirD, ...
        'Gauge_attributes.csv'),'f');

    srcF = fullfile(dirA,'2_timeseries',stream);
    srcQ = fullfile(dirG,'2_timeseries',stream);

    if ~isfolder(srcF)
        error('data_helpers:ATMissingForcing', ...
            ['Missing A_basins_total_upstrm/' ...
            '2_timeseries/%s.'],stream);
    end
    if ~isfolder(srcQ)
        error('data_helpers:ATMissingStreamflow', ...
            'Missing D_gauges/2_timeseries/%s.',stream);
    end

    DF = dir(fullfile(srcF,'ID_*.csv'));
    DQ = dir(fullfile(srcQ,'ID_*.csv'));
    for i = 1:numel(DF)
        copyfile(fullfile(DF(i).folder,DF(i).name), ...
            fullfile(dstF,DF(i).name),'f');
    end
    for i = 1:numel(DQ)
        copyfile(fullfile(DQ(i).folder,DQ(i).name), ...
            fullfile(dstQ,DQ(i).name),'f');
    end

    local_write_at_gauge_information(dirD,logFcn);

    nF = local_count_files(fullfile(dirD, ...
        stream,'timeseries','forcing'),{'ID_*.csv'});
    nQ = local_count_files(fullfile(dirD, ...
        stream,'timeseries','streamflow'),{'ID_*.csv'});

    logFcn(sprintf(['CAMELS-AT %s install finished: %d forcing ' ...
        'files and %d streamflow files.'],stream,nF,nQ));
end

% ====================================================
function local_write_at_gauge_information(dirD,logFcn)
% ====================================================
    fGauge = fullfile(dirD, ...
        'Gauge_attributes.csv');
    fCatch = fullfile(dirD, ...
        'Catchment_attributes.csv');
    if ~isfile(fGauge) ...
            || ~isfile(fCatch)
        return
    end
    try
        dg = local_detect_delimiter(fGauge);
        dc = local_detect_delimiter(fCatch);
        G = readtable(fGauge,'FileType','text', ...
            'Delimiter',dg, ...
            'VariableNamingRule','preserve');
        C = readtable(fCatch,'FileType','text', ...
            'Delimiter',dc, ...
            'VariableNamingRule','preserve');
        G.Properties.VariableNames = ...
            matlab.lang.makeValidName( ...
            string(G.Properties.VariableNames));
        C.Properties.VariableNames = ...
            matlab.lang.makeValidName( ...
            string(C.Properties.VariableNames));
        if ~ismember('ID',G.Properties.VariableNames) ...
                || ~ismember('ID', ...
                C.Properties.VariableNames)
            return
        end
        idG = string(G.ID);
        idG = regexprep(idG,'\.0+$','');
        idC = string(C.ID);
        idC = regexprep(idC,'\.0+$','');
        [tf,loc] = ismember(idG,idC);
        T = table();
        T.gauge_id = idG(:);
        T.gauge_name = "AT_" + idG(:);
        T.gauge_lat = local_table_var(G, ...
            {'gauge_lat','lat', ...
            'latitude','lat_wgs84'});
        T.gauge_lon = local_table_var(G, ...
            {'gauge_lon','lon', ...
            'longitude','lon_wgs84'});
        T.gauge_elev = local_table_var(G, ...
            {'gauge_elev','elev','elevation'});
        area = nan(height(G),1);
        if any(tf)
            a = local_table_var(C, ...
                {'area_calc', ...
                'area','area_gov'});
            area(tf) = a(loc(tf));
        end
        T.area_km2 = area;
        writetable(T,fullfile(dirD, ...
            'gauge_information.txt'), ...
            'Delimiter','\t', ...
            'FileType','text');
    catch ME
        logFcn(['Could not write CAMELS-AT ' ...
            'gauge_information.txt: ' ...
            ME.message]);
    end
end

% ============================================
function delim = local_detect_delimiter(fname)
% ============================================
%LOCAL_DETECT_DELIMITER Detect comma, semicolon, or tab delimited text.

    fid = fopen(char(fname),'r');
    if fid < 0
        error('data_helpers:openFailed', ...
            'Could not open %s',char(string(fname)));
    end
    cleanup = onCleanup(@() fclose(fid));

    line = fgetl(fid);
    if ~ischar(line)
        error('data_helpers:emptyFile', ...
            'Empty file: %s',char(string(fname)));
    end

    nTab = count(string(line),"\t");
    nSemicolon = count(string(line),";");
    nComma = count(string(line),",");

    if nTab >= nSemicolon ...
            && nTab >= nComma
        delim = '\t';
    elseif nSemicolon >= nComma
        delim = ';';
    else
        delim = ',';
    end
end

% ===================================
function x = local_table_var(T,cands)
% ===================================
    x = nan(height(T),1);
    vn = string(T.Properties.VariableNames);
    vnl = lower(regexprep(vn,'[^a-zA-Z0-9]',''));
    for i = 1:numel(cands)
        c = lower(regexprep(string(cands{i}), ...
            '[^a-zA-Z0-9]',''));
        j = find(vnl == c,1,'first');
        if isempty(j)
            j = find(contains(vnl,c),1,'first');
        end
        if ~isempty(j)
            y = T{:,j};
            if isnumeric(y)
                x = double(y(:));
            else
                x = str2double(strrep( ...
                    string(y(:)),',','.'));
            end
            return
        end
    end
end

% =============================================
function out = local_download_au(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_AU Download and install CAMELS-AU daily data.
    logFcn = local_ui_log(ui);
    if ~strcmp(stream,'daily')
        error('data_helpers:AUHourlyUnsupported', ...
            ['Hourly data are not ' ...
            'available for CAMELS-AU.']);
    end
    
    dirD = local_cfg_dirD(cfg, ...
        'CAMELS_AU');
    dirDaily = fullfile(dirD,'daily');
    local_mkdir(dirD);
    local_mkdir(dirDaily);
    local_mkdir(fullfile(dirDaily, ...
        'streamflow'));
    local_mkdir(fullfile(dirDaily, ...
        'precipitation'));
    local_mkdir(fullfile(dirDaily, ...
        'evaporative_demand'));
    local_mkdir(fullfile(dirDaily, ...
        'temperature'));
    
    if ~local_ui_confirm(ui, ...
            'Download Australia data', ...
            ['CAMELS-AU requires ' ...
            'several PANGAEA files. ' ...
            'The download may be large.' newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end
    
    files = local_camels_au_files();
    downloadDir = local_default_download_dir();
    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Australia data', ...
        'Starting CAMELS-AU download ...');
    
    try
        for i = 1:numel(files)
            if userCanceled
                out = struct('ok',false,'canceled',true);
                local_close_progress(d);
                return
            end
    
            f = files(i);
            targetFile = fullfile(downloadDir,f.fileName);
            if local_au_target_complete(dirD,f)
                logFcn(sprintf(['%s already complete; ' ...
                    'skipping download/install.'],f.label));
                files(i).install = false;
                files(i).localFile = targetFile;
                continue
            end
    
            files(i).install = true;
            if isfile(targetFile)
                logFcn(['Using existing file: ' targetFile]);
                files(i).localFile = targetFile;
            else
                logFcn(sprintf('Downloading %s to: %s', ...
                    f.label,targetFile));
                local_download_file_retry(f.url,targetFile, ...
                    f.label,d,@() userCanceled, ...
                    @(info) ui_progress(info,i,numel(files)), ...
                    logFcn,3);
                files(i).localFile = targetFile;
            end
        end
    
        installRoot = fullfile(downloadDir,'CAMELS_AU_download');
        if isfolder(installRoot)
            try
                rmdir(installRoot,'s');
            catch
            end
        end
        local_mkdir(installRoot);
    
        for i = 1:numel(files)
            f = files(i);
            if ~f.isZip ...
                    || ~isfield(files(i),'install') ...
                    || ~files(i).install
                continue
            end
            dst = fullfile(installRoot, ...
                erase(f.fileName,'.zip'));
            local_mkdir(dst);
            logFcn(sprintf('Unzipping %s ...', ...
                f.fileName));
            unzip(files(i).localFile,dst);
            files(i).unzipDir = dst;
        end
    
        % 1) id/name metadata -> Data/CAMELS_AU
        local_install_au_metadata( ...
            files,dirD,logFcn);
    
        % 2) streamflow -> Data/CAMELS_AU/streamflow
        local_install_au_streamflow( ...
            files,dirD,logFcn);
    
        % 3) hydrometeorology -> precipitation, 
        %                        evaporative_demand, temperature
        local_install_au_hydromet( ...
            files,dirD,logFcn);
    
        % 4) master attributes table -> Data/CAMELS_AU
        local_install_au_master( ...
            files,dirD,logFcn);

        % Clean downloaded archives/files and temporary unzip tree from
        % the user's Downloads folder. Installed data remain only under
        % Data/CAMELS_AU.
        local_cleanup_download_files(files,logFcn);
        local_cleanup_download_folder(installRoot,logFcn);
    
        out = struct('ok',true, ...
            'region','CAMELS_AU', ...
            'dirD',dirD, ...
            'downloadDir','', ...
            'installRoot','');
        logFcn('CAMELS-AU install finished.');
    
    catch ME
        try
            if exist('files','var')
                local_cleanup_download_files(files,logFcn);
            end
        catch
        end
        try
            if exist('installRoot','var')
                local_cleanup_download_folder(installRoot,logFcn);
            end
        catch
        end
        logFcn('CAMELS-AU install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end
    
    local_close_progress(d);

    function ui_progress(info,iFile,nFile)
        if nargin >= 2
            if isfield(info,'frac') ...
                    && isfinite(info.frac)
                info.frac = ((iFile - 1) ...
                    + info.frac) / nFile;
            else
                info.frac = (iFile - 1) / nFile;
            end
            if isfield(info,'label')
                info.label = sprintf('%s (%d/%d)', ...
                    info.label,iFile,nFile);
            else
                info.label = sprintf('File %d/%d', ...
                    iFile,nFile);
            end
        end
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% =============================================
function out = local_download_br(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_BR Download and install CAMELS-BR daily data.
    logFcn = local_ui_log(ui);
    if ~strcmp(stream,'daily')
        error('data_helpers:BRHourlyUnsupported', ...
            'Hourly data are not available for CAMELS-BR.');
    end
    
    root = local_cfg_root(cfg); %#ok<NASGU>
    dirD = local_cfg_dirD(cfg,'CAMELS_BR');
    local_mkdir(dirD);
    
    if ~local_ui_confirm(ui, ...
            'Download Brazil data', ...
            ['CAMELS-BR requires several Zenodo files. ' ...
            'The download may be large.' newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end
    
    files = local_camels_br_files();
    downloadDir = local_default_download_dir();
    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Brazil data', ...
        'Starting CAMELS-BR download ...');
    
    try
        for i = 1:numel(files)
            if userCanceled
                out = struct('ok', ...
                    false,'canceled',true);
                local_close_progress(d);
                return
            end
    
            f = files(i);
            targetFile = fullfile(downloadDir, ...
                f.fileName);
            files(i).install = local_br_needs_download(dirD,f);
            if ~files(i).install
                logFcn(sprintf(['%s already complete; ' ...
                    'skipping download/install.'],f.label));
                files(i).localFile = targetFile;
                continue
            end
            if isfile(targetFile) ...
                    && isfield(f,'md5') ...
                    && ~isempty(f.md5)
                try
                    local_verify_md5(targetFile, ...
                        f.md5,logFcn);
                    logFcn(sprintf(['Using verified ' ...
                        'existing file: %s'], ...
                        targetFile));
                    files(i).localFile = ...
                        targetFile;
                    continue
                catch
                    logFcn(sprintf(['Existing file ' ...
                        'failed MD5; deleting: %s'], ...
                        targetFile));
                    delete(targetFile);
                end
            end
            logFcn(sprintf('Downloading %s to: %s', ...
                f.label,targetFile));
            local_download_file_retry(f.url, ...
                targetFile,f.label,d, ...
                @() userCanceled, ...
                @(info) ui_progress(info, ...
                i,numel(files)), ...
                logFcn,3);
    
            if isfield(f,'md5') ...
                    && ~isempty(f.md5)
                local_verify_md5(targetFile, ...
                    f.md5,logFcn);
            end
            files(i).localFile = targetFile;
        end
    
        installRoot = fullfile(downloadDir, ...
            'CAMELS_BR_download');
        if isfolder(installRoot)
            try
                rmdir(installRoot,'s');
            catch
            end
        end
        local_mkdir(installRoot);
    
        for i = 1:numel(files)
            f = files(i);
            if ~f.isZip ...
                    || ~isfield(f,'install') ...
                    || ~f.install
                continue
            end
            dst = fullfile(installRoot, ...
                erase(f.fileName,'.zip'));
            local_mkdir(dst);
            logFcn(sprintf('Unzipping %s ...', ...
                f.fileName));
            unzip(f.localFile,dst);
            files(i).unzipDir = dst;
        end
    
        % 1) attributes -> Data/CAMELS_BR
        if local_br_should_install(files,'attributes')
            local_install_br_zip(files, ...
                'attributes',dirD, ...
                {'*.txt','*.csv', ...
                '*.xlsx','*.xls'}, ...
                logFcn);
        end
        % 2) streamflow -> Data/CAMELS_BR/streamflow
        % The streamflow and precipitation archives should not be mixed in one
        % folder; otherwise files can be overwritten or misread later.
        if local_br_should_install(files, ...
                'streamflow')
            local_install_br_zip(files, ...
                'streamflow', ...
                fullfile(dirD,'daily', ...
                'streamflow'), ...
                {'*.txt'}, ...
                logFcn);
        end
        % 3) precipitation -> Data/CAMELS_BR/precipitation
        if local_br_should_install(files, ...
                'precipitation')
            local_install_br_zip(files, ...
                'precipitation', ...
                fullfile(dirD,'daily', ...
                'precipitation'), ...
                {'*.txt'}, ...
                logFcn);
        end
        % 4) potential evapotranspiration -> Data/CAMELS_BR/pet
        if local_br_should_install(files, ...
                'pet')
            local_install_br_zip(files, ...
                'pet', ...
                fullfile(dirD,'daily','pet'), ...
                {'*.txt'}, ...
                logFcn);
        end
        % 5) temperature -> Data/CAMELS_BR/temperature
        if local_br_should_install(files, ...
                'temperature')
            local_install_br_zip(files, ...
                'temperature', ...
                fullfile(dirD,'daily', ...
                'temperature'), ...
                {'*.txt'}, ...
                logFcn);
        end
        % 6) soil moisture -> Data/CAMELS_BR/soil_moisture
        if local_br_should_install(files, ...
                'soil_moisture')
            local_install_br_zip(files, ...
                'soil_moisture', ...
                fullfile(dirD, ...
                'soil_moisture'), ...
                {'*.txt'}, ...
                logFcn);
        end
        % 7) readme -> Data/CAMELS_BR, optional
        idx = find(strcmp({files.key},'readme'),1);
        if ~isempty(idx) ...
                && isfield(files,'install') ...
                && files(idx).install ...
                && isfield(files,'localFile') ...
                && isfile(files(idx).localFile)
            copyfile(files(idx).localFile, ...
                fullfile(dirD,files(idx).fileName),'f');
            logFcn(['Copied CAMELS-BR ' ...
                'readme to: ' dirD]);
        end
    
        % Clean downloaded archives/files and temporary unzip tree from
        % the user's Downloads folder. Installed data remain only under
        % Data/CAMELS_BR.
        local_cleanup_download_files(files,logFcn);
        local_cleanup_download_folder(installRoot,logFcn);

        out = struct('ok',true, ...
            'region','CAMELS_BR', ...
            'dirD',dirD, ...
            'downloadDir','', ...
            'installRoot','');
        logFcn('CAMELS-BR install finished.');
    
    catch ME
        try
            if exist('files','var')
                local_cleanup_download_files( ...
                    files,logFcn);
            end
        catch
        end
        try
            if exist('installRoot','var')
                local_cleanup_download_folder( ...
                    installRoot,logFcn);
            end
        catch
        end
        logFcn('CAMELS-BR install failed:');
        logFcn(ME.message);
        out = struct('ok',false, ...
            'error',ME.message);
    end
    
    local_close_progress(d);

    function ui_progress(info,iFile,nFile)

        if nargin >= 2
            if isfield(info,'frac') ...
                    && isfinite(info.frac)
                info.frac = ((iFile - 1) ...
                    + info.frac) / nFile;
            else
                info.frac = (iFile - 1) / nFile;
            end
            if isfield(info,'label')
                info.label = sprintf('%s (%d/%d)', ...
                    info.label,iFile,nFile);
            else
                info.label = sprintf('File %d/%d', ...
                    iFile,nFile);
            end
        end
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% =============================================
function out = local_download_cl(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_CL Download and install CAMELS-CL daily data.
    logFcn = local_ui_log(ui);
    if ~strcmp(stream,'daily')
        error('data_helpers:CLHourlyUnsupported', ...
            'Hourly data are not available for CAMELS-CL.');
    end
    
    dirD = local_cfg_dirD(cfg,'CAMELS_CL');
    dirDaily = fullfile(dirD,'daily');
    local_mkdir(dirD);
    local_mkdir(dirDaily);
    local_mkdir(fullfile(dirDaily,'streamflow'));
    local_mkdir(fullfile(dirDaily,'precipitation'));
    local_mkdir(fullfile(dirDaily,'temperature'));
    local_mkdir(fullfile(dirDaily,'pet'));
    local_mkdir(fullfile(dirDaily,'swe'));
    
    if ~local_ui_confirm(ui, ...
            'Download Chile data', ...
            ['CAMELS-CL requires multiple PANGAEA files. ' ...
            'The download may be large.' newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end
    
    files = local_camels_cl_files();
    downloadDir = local_default_download_dir();
    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Chile data', ...
        'Starting CAMELS-CL download ...');
    
    try
        for i = 1:numel(files)
            if userCanceled
                out = struct('ok',false,'canceled',true);
                local_close_progress(d);
                return
            end
    
            f = files(i);
            targetFile = fullfile(downloadDir,f.fileName);
    
            if local_cl_target_complete(dirD,f)
                logFcn(sprintf(['%s already complete; ' ...
                    'skipping download/install.'],f.label));
                files(i).install = false;
                files(i).localFile = targetFile;
                continue
            end
    
            files(i).install = true;
    
            if isfile(targetFile)
                logFcn(['Using existing file: ' targetFile]);
                files(i).localFile = targetFile;
            else
                logFcn(sprintf('Downloading %s to: %s', ...
                    f.label,targetFile));
                local_download_file_retry(f.url,targetFile, ...
                    f.label,d,@() userCanceled, ...
                    @(info) ui_progress(info,i,numel(files)), ...
                    logFcn,3);
                files(i).localFile = targetFile;
            end
        end
    
        installRoot = fullfile(downloadDir, ...
            'CAMELS_CL_download');
        if isfolder(installRoot)
            try
                rmdir(installRoot,'s');
            catch
            end
        end
        local_mkdir(installRoot);
    
        for i = 1:numel(files)
            f = files(i);
            if ~f.isZip ...
                    || ~isfield(files(i),'install') ...
                    || ~files(i).install
                continue
            end
            dst = fullfile(installRoot, ...
                erase(f.fileName,'.zip'));
            local_mkdir(dst);
            logFcn(sprintf('Unzipping %s ...', ...
                f.fileName));
            unzip(files(i).localFile,dst);
            files(i).unzipDir = dst;
        end
    
        local_install_cl_files(files,dirD,logFcn);

        % Clean downloaded archives/files and temporary unzip tree from
        % the user's Downloads folder. Installed data remain only under
        % Data/CAMELS_CL.
        local_cleanup_download_files(files,logFcn);
        local_cleanup_download_folder(installRoot,logFcn);
    
        out = struct('ok',true, ...
            'region','CAMELS_CL', ...
            'dirD',dirD, ...
            'downloadDir','', ...
            'installRoot','');
        logFcn('CAMELS-CL install finished.');
    
    catch ME
        try
            if exist('files','var')
                local_cleanup_download_files(files,logFcn);
            end
        catch
        end
        try
            if exist('installRoot','var')
                local_cleanup_download_folder(installRoot,logFcn);
            end
        catch
        end
        logFcn('CAMELS-CL install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end
    
    local_close_progress(d);

    function ui_progress(info,iFile,nFile)
        if nargin >= 2
            if isfield(info,'frac') ...
                    && isfinite(info.frac)
                info.frac = ((iFile - 1) ...
                    + info.frac) / nFile;
            else
                info.frac = (iFile - 1) / nFile;
            end
            if isfield(info,'label')
                info.label = sprintf('%s (%d/%d)', ...
                    info.label,iFile,nFile);
            else
                info.label = sprintf('File %d/%d', ...
                    iFile,nFile);
            end
        end
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end


% =============================================
function out = local_download_ch(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_CH Download and install CAMELS-CH daily data.

    logFcn = local_ui_log(ui);

    if ~strcmp(stream,'daily')
        error('data_helpers:CHOnlyDaily', ...
            ['CAMELS-CH support is daily only. ' ...
            'Use stream = ''daily''.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_CH');
    dirTS = fullfile(dirD,'daily','timeseries');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(dirTS);
    local_mkdir(fullfile(dirTS, ...
        'observation_based'));
    local_mkdir(fullfile(dirTS, ...
        'simulation_based'));

    if ~local_ui_confirm(ui, ...
            'Download Switzerland data', ...
            ['CAMELS-CH is downloaded from Zenodo as one ZIP file. ' ...
            'The archive contains static attributes plus daily ' ...
            'observation-based and simulation-based time series.' newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipName = 'camels_ch.zip';
    zipFile = fullfile(downloadDir,zipName);
    url = ['https://zenodo.org/records/15025258/files/' ...
        'camels_ch.zip?download=1'];
    expectedMD5 = '04f909d9904375647d030c4ab8ddfdbe';

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Switzerland data', ...
        'Starting CAMELS-CH download ...');

    try
        if local_ch_install_complete(dirD)
            logFcn(['CAMELS-CH appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true, ...
                'region','CAMELS_CH', ...
                'dirD',dirD, ...
                'dirM',fullfile(dirD,'daily','timeseries'), ...
                'dirQ',fullfile(dirD,'daily','timeseries', ...
                'observation_based'), ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(zipFile)
            try
                local_verify_md5(zipFile, ...
                    expectedMD5,logFcn);
                logFcn(['Using verified existing file: ' ...
                    zipFile]);
            catch
                logFcn(['Existing CAMELS-CH ZIP failed MD5; ' ...
                    'deleting and downloading again.']);
                delete(zipFile);
            end
        end

        if ~isfile(zipFile)
            logFcn(['Downloading CAMELS-CH archive to: ' ...
                zipFile]);
            local_download_file_retry(url,zipFile, ...
                'CAMELS-CH',d, ...
                @() userCanceled, ...
                @(info) ui_progress(info), ...
                logFcn,3);
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        local_verify_md5(zipFile,expectedMD5,logFcn);

        unzipDir = fullfile(downloadDir, ...
            'CAMELS_CH_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        logFcn(['Unzipping CAMELS-CH archive to: ' ...
            unzipDir]);
        if ~isempty(d) ...
                && isvalid(d)
            try
                d.Indeterminate = true;
                d.Message = ['Unzipping CAMELS-CH archive. ' ...
                    'This can take a while ...'];
                drawnow limitrate nocallbacks
            catch
            end
        end
        unzip(zipFile,unzipDir);

        dataDir = local_find_ch_payload_dir(unzipDir);
        if isempty(dataDir) ...
                || ~isfolder(dataDir)
            error('data_helpers:CHMissingDataDir', ...
                ['Could not find the camels_ch directory ' ...
                'inside the archive.']);
        end

        local_install_ch_files(dataDir,dirD,logFcn);

        okAttr = local_all_files_exist(dirD, ...
            local_ch_metadata_files());
        okDaily = local_has_ch_timeseries(dirD);

        out = struct('ok',okAttr && okDaily, ...
            'region','CAMELS_CH', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,'daily','timeseries'), ...
            'dirQ',fullfile(dirD,'daily','timeseries', ...
            'observation_based'), ...
            'archive',zipFile, ...
            'unzipDir',unzipDir);

        if out.ok
            nObs = local_count_files(fullfile(dirD, ...
                'daily','timeseries','observation_based'), ...
                {'CAMELS_CH_obs_based_*.csv'});
            nSim = local_count_files(fullfile(dirD, ...
                'daily','timeseries','simulation_based'), ...
                {'CAMELS_CH_sim_based_*.csv'});
            logFcn(sprintf(['CAMELS-CH install finished: ' ...
                '%d observation-based and %d simulation-based ' ...
                'daily time-series CSV files.'],nObs,nSim));
            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-CH install finished, but ' ...
                'some required files appear to be missing.']);
        end

    catch ME
        try
            if exist('zipFile','var')
                local_cleanup_download_file(zipFile,logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder(unzipDir,logFcn);
            end
        catch
        end
        logFcn('CAMELS-CH install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end


% ==============================================
function out = local_download_col(cfg,stream,ui)
% ==============================================
%LOCAL_DOWNLOAD_COL Download and install CAMELS-COL daily data.

    logFcn = local_ui_log(ui);

    if ~strcmp(stream,'daily')
        error('data_helpers:COLOnlyDaily', ...
            ['CAMELS-COL support is daily only. ' ...
            'Use stream = ''daily''.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_COL');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily','timeseries'));

    if ~local_ui_confirm(ui, ...
            'Download Colombia data', ...
            ['CAMELS-COL is downloaded from Zenodo. ' ...
            'The installer downloads the attributes workbook, ' ...
            'the individual attribute CSV files, and the daily ' ...
            'hydrometeorological time-series archive.' newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    files = local_camels_col_files();
    downloadDir = local_default_download_dir();
    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Colombia data', ...
        'Starting CAMELS-COL download ...');

    try
        if local_col_install_complete(dirD)
            logFcn(['CAMELS-COL appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true, ...
                'region','CAMELS_COL', ...
                'dirD',dirD, ...
                'dirM',fullfile(dirD,'daily','timeseries'), ...
                'dirQ',fullfile(dirD,'daily','timeseries'), ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        for i = 1:numel(files)
            if userCanceled
                out = struct('ok',false,'canceled',true);
                local_close_progress(d);
                return
            end

            f = files(i);
            targetFile = fullfile(downloadDir,f.fileName);
            files(i).localFile = targetFile;
            files(i).install = true;

            if isfile(targetFile)
                try
                    local_verify_md5(targetFile,f.md5,logFcn);
                    logFcn(['Using verified existing file: ' ...
                        targetFile]);
                catch
                    logFcn(['Existing CAMELS-COL file failed MD5; ' ...
                        'deleting and downloading again: ' ...
                        targetFile]);
                    try
                        delete(targetFile);
                    catch
                    end
                end
            end

            if ~isfile(targetFile)
                logFcn(sprintf('Downloading %s to: %s', ...
                    f.label,targetFile));
                local_download_file_retry(f.url,targetFile, ...
                    f.label,d,@() userCanceled, ...
                    @(info) ui_progress(info,i,numel(files)), ...
                    logFcn,3);
            end

            local_verify_md5(targetFile,f.md5,logFcn);
            files(i).localFile = targetFile;
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        unzipDir = fullfile(downloadDir,'CAMELS_COL_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        for i = 1:numel(files)
            f = files(i);
            if ~f.isZip
                continue
            end
            dst = fullfile(unzipDir,erase(f.fileName,'.zip'));
            local_mkdir(dst);
            logFcn(sprintf('Unzipping %s ...',f.fileName));
            if ~isempty(d) && isvalid(d)
                try
                    d.Indeterminate = true;
                    d.Message = ['Unzipping ' f.fileName ' ...'];
                    drawnow limitrate nocallbacks
                catch
                end
            end
            unzip(files(i).localFile,dst);
            files(i).unzipDir = dst;
        end

        local_install_col_files(files,dirD,logFcn);

        okAttr = local_all_files_exist(dirD, ...
            local_col_metadata_files());
        okDaily = local_has_col_timeseries(dirD);

        out = struct('ok',okAttr && okDaily, ...
            'region','CAMELS_COL', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,'daily','timeseries'), ...
            'dirQ',fullfile(dirD,'daily','timeseries'), ...
            'archive','', ...
            'unzipDir','');

        if out.ok
            nTs = local_count_files(fullfile(dirD, ...
                'daily','timeseries'),{'Hydromet_data_*.txt'});
            logFcn(sprintf(['CAMELS-COL install finished: ' ...
                '%d daily hydrometeorological time-series files.'], ...
                nTs));
            local_cleanup_download_files(files,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
        else
            logFcn(['CAMELS-COL install finished, but ' ...
                'some required files appear to be missing.']);
        end

    catch ME
        try
            if exist('files','var')
                local_cleanup_download_files(files,logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder(unzipDir,logFcn);
            end
        catch
        end
        logFcn('CAMELS-COL install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info,iFile,nFile)
        if nargin >= 2
            if isfield(info,'frac') ...
                    && isfinite(info.frac)
                info.frac = ((iFile - 1) ...
                    + info.frac) / nFile;
            else
                info.frac = (iFile - 1) / nFile;
            end
            if isfield(info,'label')
                info.label = sprintf('%s (%d/%d)', ...
                    info.label,iFile,nFile);
            else
                info.label = sprintf('File %d/%d', ...
                    iFile,nFile);
            end
        end
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% =============================================
function out = local_download_dk(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_DK Download and install CAMELS-DK daily data.

    logFcn = local_ui_log(ui);

    if ~strcmp(stream,'daily')
        error('data_helpers:DKOnlyDaily', ...
            ['CAMELS-DK support is daily only. ' ...
            'Use stream = ''daily''.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_DK');
    dirT = fullfile(dirD,'daily','timeseries');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(dirT);

    if ~local_ui_confirm(ui, ...
            'Download Denmark data', ...
            ['CAMELS-DK is downloaded from GEUS Dataverse as ' ...
             'dataverse_files.zip. The installer copies the seven ' ...
             'attribute CSV files to Data\\CAMELS_DK and extracts the ' ...
             '304 gauged-catchment CSV files to ' ...
             'Data\\CAMELS_DK\\timeseries.' newline newline ...
             'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipName = 'dataverse_files.zip';
    zipFile = fullfile(downloadDir,zipName);
    url = ['https://dataverse.geus.dk/api/access/dataset/' ...
        ':persistentId/?persistentId=' ...
        'doi%3A10.22008%2FFK2%2FAZXSYP'];

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Denmark data', ...
        'Starting CAMELS-DK download ...');

    try
        if local_dk_install_complete(dirD)
            logFcn(['CAMELS-DK appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true, ...
                'region','CAMELS_DK', ...
                'dirD',dirD, ...
                'dirM',dirT, ...
                'dirQ',dirT, ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        if ~isfile(zipFile)
            logFcn(['Downloading CAMELS-DK archive to: ' zipFile]);
            local_download_file_retry(url,zipFile, ...
                'CAMELS-DK',d,@() userCanceled, ...
                @(info) ui_progress(info),logFcn,3);
        else
            logFcn(['Using existing CAMELS-DK archive: ' zipFile]);
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        unzipDir = fullfile(downloadDir,'CAMELS_DK_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        try
            d.Indeterminate = true;
            d.Message = ['Unzipping CAMELS-DK archive. ' ...
                'This can take a while ...'];
            drawnow limitrate nocallbacks
        catch
        end
        unzip(zipFile,unzipDir);

        local_install_dk_files(unzipDir,dirD,logFcn);

        ok = local_dk_install_complete(dirD);
        out = struct('ok',ok, ...
            'region','CAMELS_DK', ...
            'dirD',dirD, ...
            'dirM',dirT, ...
            'dirQ',dirT, ...
            'archive',zipFile, ...
            'unzipDir',unzipDir);

        if ok
            nA = local_count_files(dirD,{'*.csv'});
            nT = local_count_files(dirT, ...
                {'CAMELS_DK_obs_based_*.csv'});
            logFcn(sprintf(['CAMELS-DK install finished: ' ...
                '%d attribute CSV files and %d gauged-catchment ' ...
                'daily time-series CSV files.'],nA,nT));
            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-DK install finished, but some required ' ...
                'files appear to be missing.']);
        end

    catch ME
        % Preserve the archive and temporary extraction tree on failure.
        % They are valuable for diagnosis and allow the next run to retry
        % installation without downloading the full Dataverse archive again.
        logFcn('CAMELS-DK install failed:');
        logFcn(ME.message);
        if exist('zipFile','var') && isfile(zipFile)
            logFcn(['Preserved downloaded archive for retry: ' zipFile]);
        end
        if exist('unzipDir','var') && isfolder(unzipDir)
            logFcn(['Preserved temporary extraction folder: ' unzipDir]);
        end
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% =====================================================
function tf = local_dk_install_complete(dirD)
% =====================================================
    tf = local_count_files(dirD,{'*.csv'}) >= 7 ...
        && local_count_files(fullfile(dirD,'daily','timeseries'), ...
        {'CAMELS_DK_obs_based_*.csv'}) >= 304;
end

% =====================================================
function local_install_dk_files(srcRoot,dirD,logFcn)
%LOCAL_INSTALL_DK_FILES Install attributes and gauged daily time series.

    dirT = fullfile(dirD,'daily','timeseries');
    local_mkdir(dirD);
    local_mkdir(dirT);

    % ---------------------------------------------------------------
    % Attributes
    % ---------------------------------------------------------------
    % Dataverse ZIPs can preserve the Attributes directory in several
    % layouts, including an additional dataset folder or another directory
    % below Attributes. Search recursively and accept CSVs whose path
    % contains an Attributes directory.
    allCsv = dir(fullfile(srcRoot,'**','*.csv'));
    isAttr = false(numel(allCsv),1);

    for i = 1:numel(allCsv)
        p = strrep(fullfile(allCsv(i).folder,allCsv(i).name),'\','/');
        isAttr(i) = ~isempty(regexpi(p,'(^|/)Attributes(/|$)','once')) ...
            && isempty(regexpi(allCsv(i).name, ...
                '^CAMELS_DK_obs_based_.*\.csv$','once'));
    end

    A = allCsv(isAttr);

    % Fallback for archives that flatten the directory names. Attribute
    % files are the CAMELS-DK CSVs that are not observed-basin time series.
    if numel(A) < 7
        isCandidate = false(numel(allCsv),1);
        for i = 1:numel(allCsv)
            nam = allCsv(i).name;
            isCandidate(i) = startsWith(nam,'CAMELS_DK_', ...
                    'IgnoreCase',true) ...
                && isempty(regexpi(nam, ...
                    '^CAMELS_DK_obs_based_.*\.csv$','once'));
        end
        A = allCsv(isCandidate);
    end

    % Deduplicate by filename in case the archive contains repeated copies.
    if ~isempty(A)
        [~,ia] = unique(lower(string({A.name})),'stable');
        A = A(ia);
    end

    if numel(A) < 7
        % Log the archive layout to make any future archive change obvious.
        top = dir(srcRoot);
        top = top(~ismember({top.name},{'.','..'}));
        if isempty(top)
            topText = '(empty extraction folder)';
        else
            topText = strjoin(string({top.name}),', ');
        end
        error('data_helpers:DKMissingAttributes', ...
            ['Could not locate the seven CAMELS-DK attribute CSV files. ' ...
             'Found %d candidate attribute CSV files. Top-level archive ' ...
             'entries: %s'],numel(A),char(topText));
    end

    for i = 1:numel(A)
        copyfile(fullfile(A(i).folder,A(i).name), ...
            fullfile(dirD,A(i).name),'f');
    end

    % ---------------------------------------------------------------
    % Gauged catchments
    % ---------------------------------------------------------------
    % Locate the nested ZIP without assuming a fixed parent directory.
    Z = dir(fullfile(srcRoot,'**','*.zip'));
    keep = strcmpi(string({Z.name}),'Gauged_catchments.zip');
    Z = Z(keep);
    if isempty(Z)
        error('data_helpers:DKMissingGaugedArchive', ...
            ['Could not find Gauged_catchments.zip inside ' ...
             'dataverse_files.zip.']);
    end

    gaugedZip = fullfile(Z(1).folder,Z(1).name);
    gaugedDir = fullfile(srcRoot,'_gauged_catchments');
    if isfolder(gaugedDir)
        rmdir(gaugedDir,'s');
    end
    local_mkdir(gaugedDir);
    unzip(gaugedZip,gaugedDir);

    T = dir(fullfile(gaugedDir,'**','CAMELS_DK_obs_based_*.csv'));

    % Deduplicate by filename before validating and copying.
    if ~isempty(T)
        [~,ia] = unique(lower(string({T.name})),'stable');
        T = T(ia);
    end

    if numel(T) ~= 304
        error('data_helpers:DKTimeseriesCount', ...
            ['Expected exactly 304 unique gauged-catchment CSV files, ' ...
             'but found %d.'],numel(T));
    end

    for i = 1:numel(T)
        copyfile(fullfile(T(i).folder,T(i).name), ...
            fullfile(dirT,T(i).name),'f');
    end

    logFcn(sprintf(['Copied %d CAMELS-DK attribute CSV files and ' ...
        '%d gauged-catchment time-series CSV files.'], ...
        numel(A),numel(T)));
end

% =============================================

% =============================================
function out = local_download_cz(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_CZ Download and install CAMELS-CZ data.

    logFcn = local_ui_log(ui);
    if ~any(strcmp(stream,{'daily','hourly'}))
        error('data_helpers:CZUnsupportedResolution', ...
            'Unsupported CAMELS-CZ resolution: %s.',stream);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_CZ');
    dirTS = fullfile(dirD,stream,'timeseries');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,stream));
    local_mkdir(dirTS);

    % Migrate an older installation created under
    % Data/CAMELS_CZ/timeseries into the canonical daily/timeseries
    % location before deciding whether a new download is needed.
    legacyTS = fullfile(dirD,'timeseries');
    if strcmp(stream,'daily') ...
            && local_count_files(dirTS,{'camelscz_*.csv'}) < 249 ...
            && local_count_files(legacyTS,{'camelscz_*.csv'}) >= 249
        Dlegacy = dir(fullfile(legacyTS,'camelscz_*.csv'));
        for iLegacy = 1:numel(Dlegacy)
            copyfile(fullfile(Dlegacy(iLegacy).folder,Dlegacy(iLegacy).name), ...
                fullfile(dirTS,Dlegacy(iLegacy).name),'f');
        end
        logFcn(sprintf(['Migrated %d CAMELS-CZ time-series files from ' ...
            '%s to %s.'],numel(Dlegacy),legacyTS,dirTS));
    end

    archiveLabel = 'CAMELS-CZ.zip';
    archiveStem = 'CAMELS_CZ_download';
    expectedMD5 = 'f614413ba5c58630ce5cae4ed1fc7e3c';
    if strcmp(stream,'hourly')
        archiveLabel = 'CAMELS-CZ-Hourly.zip';
        archiveStem = 'CAMELS_CZ_hourly_download';
        expectedMD5 = '182308706b70a616045ca12ecfd6466f';
    end
    if ~local_ui_confirm(ui,'Download Czechia data', ...
            sprintf(['Download %s from Zenodo? The archive contains ' ...
            '249 %s catchment time series and three attribute files.'], ...
            archiveLabel,stream))
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipFile = fullfile(downloadDir,archiveLabel);
    unzipDir = fullfile(downloadDir,archiveStem);
    url = ['https://zenodo.org/records/17769325/files/' ...
        archiveLabel '?download=1'];

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui,'Downloading Czechia data', ...
        'Starting CAMELS-CZ download ...');

    try
        if local_cz_install_complete(dirD,stream)
            logFcn(sprintf(['CAMELS-CZ %s data appear complete; ' ...
                'skipping download.'],stream));
            out = struct('ok',true,'region','CAMELS_CZ', ...
                'dirD',dirD,'dirM',dirTS,'dirQ',dirTS, ...
                'archive','','unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(zipFile)
            try
                local_verify_md5(zipFile,expectedMD5,logFcn);
                logFcn(['Using verified existing file: ' zipFile]);
            catch
                delete(zipFile);
            end
        end

        if ~isfile(zipFile)
            local_download_file_retry(url,zipFile,'CAMELS-CZ',d, ...
                @() userCanceled,@(info) ui_progress(info),logFcn,3);
        end
        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end
        local_verify_md5(zipFile,expectedMD5,logFcn);

        if isfolder(unzipDir)
            rmdir(unzipDir,'s');
        end
        local_mkdir(unzipDir);
        logFcn('CAMELS-CZ archive verified. Extracting files ...');
        if ~isempty(d) && isvalid(d)
            d.Indeterminate = true;
            d.Message = 'Unzipping CAMELS-CZ archive ...';
            drawnow limitrate nocallbacks
        end
        unzip(zipFile,unzipDir);
        logFcn('CAMELS-CZ archive extracted. Installing files ...');

        payload = local_find_cz_payload(unzipDir);
        if isempty(payload)
            error('data_helpers:CZMissingPayload', ...
                'Could not locate the CAMELS-CZ payload inside the archive.');
        end

        srcAttr = fullfile(payload,'attributes','camelscz');
        srcTS = fullfile(payload,'timeseries','csv','camelscz');
        if ~isfolder(srcTS)
            Dsrc = dir(fullfile(unzipDir,'**','camelscz_*.csv'));
            if isempty(Dsrc)
                error('data_helpers:CZMissingTimeseries', ...
                    'Could not locate camelscz_*.csv files after unzipping.');
            end
            srcTS = Dsrc(1).folder;
        end
        attrs = {'attributes_caravan_camelscz.csv', ...
            'attributes_hydroatlas_camelscz.csv', ...
            'attributes_other_camelscz.csv'};
        for k = 1:numel(attrs)
            f = fullfile(srcAttr,attrs{k});
            if ~isfile(f)
                error('data_helpers:CZMissingAttribute', ...
                    'Missing archive file: %s',f);
            end
            copyfile(f,fullfile(dirD,attrs{k}),'f');
        end

        D = dir(fullfile(srcTS,'camelscz_*.csv'));
        if numel(D) < 249
            error('data_helpers:CZMissingTimeseries', ...
                'Expected 249 CAMELS-CZ CSV files; found %d.',numel(D));
        end
        if ~isempty(d) && isvalid(d)
            d.Indeterminate = false;
            d.Value = 0;
            d.Message = sprintf('Installing CAMELS-CZ files 0/%d ...', ...
                numel(D));
            drawnow limitrate nocallbacks
        end
        for k = 1:numel(D)
            copyfile(fullfile(D(k).folder,D(k).name), ...
                fullfile(dirTS,D(k).name),'f');
            if mod(k,5)==0 || k==numel(D)
                if ~isempty(d) && isvalid(d)
                    d.Value = k/numel(D);
                    d.Message = sprintf('Installing CAMELS-CZ files %d/%d ...', ...
                        k,numel(D));
                    drawnow limitrate nocallbacks
                end
            end
        end

        ok = local_cz_install_complete(dirD,stream);
        out = struct('ok',ok,'region','CAMELS_CZ', ...
            'dirD',dirD,'dirM',dirTS,'dirQ',dirTS, ...
            'archive',zipFile,'unzipDir',unzipDir);

        if ok
            logFcn(sprintf(['CAMELS-CZ install finished: %d %s ' ...
                'time-series files and three attribute files.'], ...
                numel(D),stream));
            if ~isempty(d) && isvalid(d)
                d.Indeterminate = false;
                d.Value = 1;
                d.Message = sprintf(['CAMELS-CZ installation complete: ' ...
                    '%d/%d files installed.'],numel(D),numel(D));
                drawnow
            end
            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        end
    catch ME
        logFcn('CAMELS-CZ install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end
    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

function tf = local_cz_install_complete(dirD,stream)
    dirTS = local_cz_timeseries_install_dir(dirD,stream);
    tf = local_all_files_exist(dirD,{ ...
        'attributes_caravan_camelscz.csv', ...
        'attributes_hydroatlas_camelscz.csv', ...
        'attributes_other_camelscz.csv'}) ...
        && ~isempty(dirTS) ...
        && local_count_files(dirTS,{'camelscz_*.csv'}) >= 249 ...
        && local_cz_stream_signature(dirTS,stream);
end

function tf = local_cz_stream_signature(dirTS,stream)
%LOCAL_CZ_STREAM_SIGNATURE Reject files installed at the wrong resolution.

    tf = false;
    D = dir(fullfile(dirTS,'camelscz_*.csv'));
    if isempty(D)
        return
    end

    fid = fopen(fullfile(D(1).folder,D(1).name),'r');
    if fid < 0
        return
    end
    cleaner = onCleanup(@()fclose(fid));
    header = fgetl(fid);
    if ~ischar(header)
        return
    end

    switch lower(char(stream))
        case 'daily'
            tf = contains(header, ...
                'potential_evaporation_sum_FAO_PENMAN_MONTEITH_CHMI') ...
                || contains(header, ...
                'potential_evaporation_sum_FAO_PENMAN_MONTEITH_ERA5_LAND');
        case 'hourly'
            tf = contains(header,'potential_evaporation_ERA5_LAND') ...
                && ~contains(header, ...
                'potential_evaporation_sum_FAO_PENMAN_MONTEITH_CHMI');
    end
end

function dirTS = local_cz_timeseries_install_dir(dirD,stream)
%LOCAL_CZ_TIMESERIES_INSTALL_DIR Return canonical CZ time-series path.

    candidate = fullfile(dirD,stream,'timeseries');
    if isfolder(candidate) ...
            && local_count_files(candidate,{'camelscz_*.csv'}) > 0
        dirTS = candidate;
    else
        dirTS = '';
    end
end

function payload = local_find_cz_payload(root)
    payload = '';
    candidates = {root, ...
        fullfile(root,'CAMELS-CZ'), ...
        fullfile(root,'CAMELS_CZ'), ...
        fullfile(root,'CAMELS-CZ-Hourly'), ...
        fullfile(root,'CAMELS_CZ_Hourly')};
    for k = 1:numel(candidates)
        if isfolder(fullfile(candidates{k},'attributes','camelscz')) ...
                && isfolder(fullfile(candidates{k}, ...
                'timeseries','csv','camelscz'))
            payload = candidates{k};
            return
        end
    end
    D = dir(fullfile(root,'**','attributes','camelscz'));
    if ~isempty(D)
        p = fileparts(fileparts(fullfile(D(1).folder,D(1).name)));
        if isfolder(fullfile(p,'timeseries','csv','camelscz'))
            payload = p;
        end
    end
end

function out = local_download_de(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_DE Download and install CAMELS-DE daily/hourly data.

    logFcn = local_ui_log(ui);

    stream = lower(strtrim(char(stream)));

    if ~any(strcmp(stream,{'daily','hourly'}))
        error('data_helpers:DEBadStream', ...
            'CAMELS-DE supports stream = daily or hourly.');
    end

    dirRoot = local_cfg_dirD(cfg,'CAMELS_DE');
    dirD = fullfile(dirRoot,stream);

    local_mkdir(dirRoot);
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily','timeseries'));

    switch stream
        case 'daily'
            out = local_download_de_daily(ui,dirD,logFcn);
        case 'hourly'
            out = local_download_de_hourly(ui,dirD,logFcn);
    end
end

% =====================================================
function out = local_download_de_hourly(ui,dirD,logFcn)
% =====================================================
%LOCAL_DOWNLOAD_DE Download and install CAMELS-DE hourly data.

    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily','timeseries'));

    if ~local_ui_confirm(ui, ...
            'Download Germany data', ...
            ['CAMELS-DE-1h is very large:' newline ...
             '  ZIP file: about 14 GB' newline ...
             '  Unzipped files: about 50 GB' newline ...
             'During installation, more than 100 ' ...
             'GB of free disk space ' ...
             'may be needed temporarily because ' ...
             'the archive is downloaded, ' ...
             ['unzipped, installed into Data\CAMELS_DE, ' ...
             'and then cleaned up.'] newline newline ...
             'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipName = 'CAMELS-DE-1h.zip';
    zipFile = fullfile(downloadDir,zipName);
    url = ['https://datapub.gfz.de/download/' ...
        '10.5880.FIDGEO.2026.045-Trfghbv/' ...
        '2026-045_Dolich-et-al_data/' ...
        'CAMELS-DE-1h.zip'];

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Germany data', ...
        'Starting CAMELS-DE-1h download ...');

    try
        logFcn(['Downloading CAMELS-DE-1h archive to: ' ...
            zipFile]);

        local_download_file_retry(url,zipFile, ...
            'CAMELS-DE-1h',d, ...
            @() userCanceled, ...
            @(info) ui_progress(info), ...
            logFcn,3);

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end
        unzipDir = fullfile(dirD,'_download_tmp');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        logFcn(['Unzipping CAMELS-DE-1h archive to: ' ...
            unzipDir]);
        if ~isempty(d) ...
                && isvalid(d)
            try
                d.Indeterminate = true;
                d.Message = ['Unzipping CAMELS-DE-1h ' ...
                    'archive. This can take a while ...'];
                drawnow limitrate nocallbacks
            catch
            end
        end
        unzip(zipFile,unzipDir);

        dataDir = local_find_child_dir(unzipDir, ...
            'CAMELS-DE-1h');
        if isempty(dataDir)
            dataDir = local_find_first_payload_dir( ...
                unzipDir,'CAMELS-DE-1h');
        end
        if isempty(dataDir) ...
                || ~isfolder(dataDir)
            error('data_helpers:DEMissingDataDir', ...
                ['Could not find the CAMELS-DE-1h ' ...
                'directory inside the archive.']);
        end

        local_install_de_files(dataDir,dirD,logFcn,'hourly');

        nAttr = local_count_files(dirD, ...
            {'CAMELS_DE_1h_*_attributes.csv', ...
             'CAMELS_DE_1h_simulation_benchmark.csv'});
        nTs = local_count_files(fullfile(dirD,'timeseries'), ...
            {'CAMELS_DE_1h_hydromet_timeseries_*.csv'});

        okAttr = local_all_files_exist(dirD, ...
            local_de_metadata_files('hourly'));
        okTs = nTs > 0;

        out = struct('ok',okAttr && okTs, ...
            'region','CAMELS_DE', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,'timeseries'), ...
            'archive',zipFile, ...
            'unzipDir',unzipDir, ...
            'n_attribute_files',nAttr, ...
            'n_timeseries_files',nTs);

        if out.ok
            logFcn(sprintf(['CAMELS-DE install finished: ' ...
                '%d attribute/benchmark CSV files and ' ...
                '%d hourly time-series CSV files.'], ...
                nAttr,nTs));
        
            % Remove downloaded archive and temporary unzip directory
            try
                if isfile(zipFile)
                    delete(zipFile);
                    logFcn(['Deleted downloaded ' ...
                        'archive: ' zipFile]);
                end
            catch MEdel
                logFcn(['Could not delete ' ...
                    'downloaded archive: ' ...
                    MEdel.message]);
            end
        
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-DE install finished, ' ...
                'but some required files ' ...
                'appear to be missing.']);
        end

    catch ME
        try
            if exist('zipFile','var')
                local_cleanup_download_file( ...
                    zipFile,logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder( ...
                    unzipDir,logFcn);
            end
        catch
        end
        logFcn('CAMELS-DE install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% ====================================================
function out = local_download_de_daily(ui,dirD,logFcn)
% ====================================================
%LOCAL_DOWNLOAD_DE_DAILY Download and install daily CAMELS-DE.

    if ~local_ui_confirm(ui, ...
            'Download Germany daily data', ...
            ['CAMELS-DE daily is downloaded from Zenodo as ' ...
             'camels_de.zip. The archive unpacks to camels_de, ' ...
             'with attribute files directly under camels_de and ' ...
             'hydrometeorological files under camels_de\timeseries.' ...
             newline newline ...
             'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipName = 'camels_de.zip';
    zipFile = fullfile(downloadDir,zipName);

    url = ['https://zenodo.org/records/13837553/files/' ...
           'camels_de.zip?download=1'];

    expectedMD5 = '438f25541d94a0db0337f6bd1cce3cd0';

    userCanceled = false;
    lastUI = tic;

    d = local_progress_dialog(ui, ...
        'Downloading Germany daily data', ...
        'Starting CAMELS-DE daily download ...');

    try
        if isfile(zipFile)
            try
                local_verify_md5(zipFile,expectedMD5,logFcn);
                logFcn(['Using verified existing file: ' zipFile]);
            catch
                logFcn(['Existing CAMELS-DE daily archive failed MD5; ' ...
                    'deleting and downloading again.']);
                delete(zipFile);
            end
        end

        if ~isfile(zipFile)
            logFcn(['Downloading CAMELS-DE daily archive to: ' zipFile]);

            local_download_file_retry(url,zipFile, ...
                'CAMELS-DE daily',d, ...
                @() userCanceled, ...
                @(info) ui_progress(info), ...
                logFcn,3);
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        local_verify_md5(zipFile,expectedMD5,logFcn);

        unzipDir = fullfile(dirD,'_download_tmp');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        logFcn(['Unzipping CAMELS-DE daily archive to: ' unzipDir]);

        if ~isempty(d) && isvalid(d)
            try
                d.Indeterminate = true;
                d.Message = ['Unzipping CAMELS-DE daily archive. ' ...
                    'This can take a while ...'];
                drawnow limitrate nocallbacks
            catch
            end
        end

        unzip(zipFile,unzipDir);

        dataDir = local_find_de_daily_payload_dir(unzipDir);

        if isempty(dataDir) ...
                || ~isfolder(dataDir)
            error('data_helpers:DEDailyMissingDataDir', ...
                ['Could not find CAMELS-DE daily payload inside ' ...
                 'the archive. Expected a folder containing ' ...
                 'CAMELS_DE_topographic_attributes.csv and a ' ...
                 'timeseries folder.']);
        end

        logFcn(['CAMELS-DE daily archive ' ...
            'payload folder: ' dataDir]);
        
        DattrProbe = dir(fullfile(dataDir, ...
            'CAMELS_DE_*attributes.csv'));
        DtsProbe = dir(fullfile(dataDir,'timeseries', ...
            'CAMELS_DE_hydromet_timeseries_DE*.csv'));
        
        logFcn(sprintf(['CAMELS-DE daily ' ...
            'payload contains %d attribute ' ...
            'CSV files and %d time-series ' ...
            'CSV files before install.'], ...
            numel(DattrProbe),numel(DtsProbe)));
        
        local_install_de_files(dataDir,dirD,logFcn,'daily');

        nAttr = local_count_files(dirD, ...
            {'CAMELS_DE_*_attributes.csv', ...
             'CAMELS_DE_simulation_benchmark.csv'});

        nTs = local_count_files(fullfile(dirD,'timeseries'), ...
            {'CAMELS_DE_hydromet_timeseries_DE*.csv'});

        okAttr = local_all_files_exist(dirD, ...
            local_de_metadata_files('daily'));

        okTs = nTs > 0;

        out = struct('ok',okAttr && okTs, ...
            'region','CAMELS_DE', ...
            'stream','daily', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,'daily','timeseries'), ...
            'dirQ',fullfile(dirD,'daily','timeseries'), ...
            'archive',zipFile, ...
            'unzipDir',unzipDir, ...
            'n_attribute_files',nAttr, ...
            'n_timeseries_files',nTs);

        if out.ok
            logFcn(sprintf(['CAMELS-DE daily install finished: ' ...
                '%d attribute/benchmark CSV files and ' ...
                '%d daily time-series CSV files.'],nAttr,nTs));

            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);

            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-DE daily install finished, but some ' ...
                'required files appear to be missing.']);
        end

    catch ME
        logFcn('CAMELS-DE daily install failed:');
        logFcn(ME.message);

        if exist('zipFile','var')
            logFcn(['Leaving downloaded archive for inspection: ' ...
                zipFile]);
        end
        if exist('unzipDir','var')
            logFcn(['Leaving temporary unzip folder for inspection: ' ...
                unzipDir]);
        end

        out = struct('ok',false, ...
            'error',ME.message, ...
            'archive',zipFile, ...
            'unzipDir',unzipDir);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% ==========================================================
function dataDir = local_find_de_daily_payload_dir(unzipDir)
% ==========================================================

    dataDir = '';

    if isempty(unzipDir) || ~isfolder(unzipDir)
        return
    end

    % Look for the topographic attributes file because it should be
    % present once per valid CAMELS-DE daily payload.
    q = dir(fullfile(unzipDir,'**', ...
        'CAMELS_DE_topographic_attributes.csv'));

    q = q(~[q.isdir]);

    if isempty(q)
        return
    end

    % The payload directory is the folder containing the attributes.
    cand = q(1).folder;

    % Confirm that the matching time-series folder exists below it.
    tsDir = fullfile(cand,'timeseries');
    if isfolder(tsDir) ...
            && ~isempty(dir(fullfile(tsDir, ...
            'CAMELS_DE_hydromet_timeseries_DE*.csv')))
        dataDir = cand;
        return
    end

    % Fallback: search below the candidate for the timeseries folder.
    qTs = dir(fullfile(cand,'**', ...
        'CAMELS_DE_hydromet_timeseries_DE*.csv'));

    qTs = qTs(~[qTs.isdir]);

    if ~isempty(qTs)
        dataDir = cand;
    end
end

% ===============================================
function out = local_download_bull(cfg,stream,ui)
% ===============================================
%LOCAL_DOWNLOAD_BULL Download and install the BULL database.

    logFcn = local_ui_log(ui);

    if ~strcmp(stream,'daily')
        error('data_helpers:BULLOnlyDaily', ...
            'The BULL database is daily only.');
    end

    dirD = local_cfg_dirD(cfg,'BULL_ES');
    dirTS = fullfile(dirD,'daily','timeseries');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(dirTS);

    if ~local_ui_confirm(ui,'Download Spain BULL data', ...
            ['BULL contains 484 Spanish basins. The official Zenodo ' ...
            'installation downloads attributes and approximately ' ...
            '2.2 GB of daily time series. Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    stageDir = fullfile(downloadDir,'BULL_ES_download');
    archives = { ...
        fullfile(downloadDir,'BULL_attributes.zip'), ...
        fullfile(downloadDir,'BULL_timeseries.zip')};
    urls = { ...
        ['https://zenodo.org/records/10605646/files/' ...
        'attributes.zip?download=1'], ...
        ['https://zenodo.org/records/10605646/files/' ...
        'timeseries.zip?download=1']};
    md5 = { ...
        '2e0321e828217b52560fd8886ef9e5f7', ...
        '3194058494c6f8170433642a4dddf1cf'};
    labels = {'BULL attributes','BULL daily time series'};
    d = local_progress_dialog(ui,'Downloading Spain BULL data', ...
        'Checking the BULL installation ...');
    userCanceled = false;
    lastUI = tic;

    try
        if local_bull_install_complete(dirD)
            logFcn(['BULL appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true,'region','BULL_ES', ...
                'dirD',dirD,'dirM',dirTS,'dirQ',dirTS, ...
                'archive','','unzipDir','');
            local_close_progress(d);
            return
        end

        for i = 1:numel(archives)
            keepPartial = false;
            if isfile(archives{i})
                try
                    local_verify_md5(archives{i},md5{i},logFcn);
                    logFcn(['Using verified existing file: ' ...
                        archives{i}]);
                catch
                    info = dir(archives{i});
                    keepPartial = i == 2 ...
                        && info.bytes < 2225187603 ...
                        && local_has_zip_signature(archives{i});
                    if keepPartial
                        logFcn(sprintf(['Resuming partial BULL time-' ...
                            'series archive at %.3f GB: %s'], ...
                            info.bytes/1e9,archives{i}));
                    else
                        logFcn(['Existing archive failed MD5; ' ...
                            'deleting and downloading again: ' ...
                            archives{i}]);
                        removed = local_remove_file_retry( ...
                            archives{i},logFcn);
                        if ~removed
                            error('data_helpers:BULLArchiveLocked', ...
                                ['The invalid BULL archive could not be ' ...
                                'removed. Close any program using the ' ...
                                'file and retry: %s'],archives{i});
                        end
                    end
                end
            end
            if i == 2 && (keepPartial || ~isfile(archives{i}))
                local_download_bull_timeseries( ...
                    urls{i},archives{i},d,logFcn);
            elseif ~isfile(archives{i})
                local_download_file_retry(urls{i},archives{i}, ...
                    labels{i},d,@() userCanceled,@ui_progress,logFcn,3);
            end
            if userCanceled
                out = struct('ok',false,'canceled',true);
                local_close_progress(d);
                return
            end
            local_verify_md5(archives{i},md5{i},logFcn);
        end

        if isfolder(stageDir)
            rmdir(stageDir,'s');
        end
        attrStage = fullfile(stageDir,'attributes');
        tsStage = fullfile(stageDir,'timeseries');
        local_mkdir(attrStage);
        local_mkdir(tsStage);
        try
            d.Indeterminate = true;
            d.Message = ['Extracting BULL archives. ' ...
                'This can take several minutes ...'];
            drawnow limitrate nocallbacks
        catch
        end
        logFcn('Extracting BULL attribute archive ...');
        unzip(archives{1},attrStage);
        logFcn('Extracting BULL time-series archive ...');
        unzip(archives{2},tsStage);

        attrRoot = local_find_bull_payload(attrStage, ...
            'attributes_other_.csv');
        tsRoot = local_find_bull_timeseries_payload(tsStage);
        if isempty(attrRoot) ...
                || isempty(tsRoot)
            error('data_helpers:BULLPayloadMissing', ...
                'Could not locate the BULL payload inside the archives.');
        end

        attrNames = {'attributes_other_.csv', ...
            'attributes_caravan_.csv', ...
            'attributes_hydroatlas_.csv'};
        for i = 1:numel(attrNames)
            copyfile(fullfile(attrRoot,attrNames{i}), ...
                fullfile(dirD,attrNames{i}),'f');
        end

        folders = {'streamflow','AEMET','ERA5_Land','EMO1_arc'};
        patterns = {'streamflow_*.csv','AEMET_*.csv', ...
            'ERA5_Land_*.csv','EMO1_*.csv'};
        for i = 1:numel(folders)
            src = fullfile(tsRoot,folders{i});
            dst = fullfile(dirTS,folders{i});
            if ~isfolder(src)
                error('data_helpers:BULLSourceMissing', ...
                    'Missing BULL source folder: %s',src);
            end
            local_copy_pattern(src,dst,patterns(i),logFcn);
        end

        ok = local_bull_install_complete(dirD);
        out = struct('ok',ok,'region','BULL_ES', ...
            'dirD',dirD,'dirM',dirTS,'dirQ',dirTS, ...
            'archive',archives{2},'unzipDir',stageDir);
        if ok
            logFcn(['BULL install finished: 484 basins with AEMET, ' ...
                'ERA5-Land, EMO-1, and observed streamflow.']);
            for i = 1:numel(archives)
                local_cleanup_download_file(archives{i},logFcn);
            end
            local_cleanup_download_folder(stageDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        else
            error('data_helpers:BULLInstallIncomplete', ...
                ['BULL extraction finished, but one or more required ' ...
                'collections contain fewer than 484 files.']);
        end
    catch ME
        logFcn('BULL install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message, ...
            'archive',archives{2},'unzipDir',stageDir);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

function local_download_bull_timeseries(url,targetFile,d,logFcn)
%LOCAL_DOWNLOAD_BULL_TIMESERIES Resume the 2.2 GB Zenodo transfer.

    expectedBytes = 2225187603;
    local_mkdir(fileparts(targetFile));

    try
        d.Indeterminate = true;
        d.Message = ['Downloading BULL daily time series (2.2 GB). ' ...
            'Interrupted transfers resume on the next attempt ...'];
        drawnow limitrate nocallbacks
    catch
    end

    if isfile(targetFile)
        info = dir(targetFile);
        if info.bytes > expectedBytes ...
                || ~local_has_zip_signature(targetFile)
            removed = local_remove_file_retry(targetFile,logFcn);
            if ~removed
                error('data_helpers:BULLArchiveLocked', ...
                    'Could not replace invalid archive: %s',targetFile);
            end
        elseif info.bytes == expectedBytes
            return
        end
    end

    if isfile(targetFile)
        info = dir(targetFile);
        logFcn(sprintf('Resuming BULL transfer at %.3f of %.3f GB ...', ...
            info.bytes/1e9,expectedBytes/1e9));
    else
        logFcn('Downloading BULL daily time series (2.225 GB) ...');
    end

    cmd = sprintf(['curl -L --fail --ssl-no-revoke --retry 10 ' ...
        '--retry-all-errors --retry-delay 5 --continue-at - ' ...
        '-o "%s" "%s"'],targetFile,url);
    [status,message] = system(cmd);
    if status ~= 0
        error('data_helpers:BULLDownloadInterrupted', ...
            ['The BULL transfer was interrupted. The partial archive ' ...
            'was retained and will resume on the next attempt. curl: %s'], ...
            strtrim(message));
    end

    info = dir(targetFile);
    if isempty(info) || info.bytes ~= expectedBytes
        error('data_helpers:BULLDownloadIncomplete', ...
            ['The BULL transfer ended at %.3f GB; expected %.3f GB. ' ...
            'The partial archive was retained for resumption.'], ...
            info.bytes/1e9,expectedBytes/1e9);
    end
end

function tf = local_bull_install_complete(dirD)

    dirTS = fullfile(dirD,'daily','timeseries');
    tf = local_all_files_exist(dirD,{ ...
        'attributes_other_.csv', ...
        'attributes_caravan_.csv', ...
        'attributes_hydroatlas_.csv'}) ...
        && local_count_files(fullfile(dirTS,'streamflow'), ...
        {'streamflow_*.csv'}) >= 484 ...
        && local_count_files(fullfile(dirTS,'AEMET'), ...
        {'AEMET_*.csv'}) >= 484 ...
        && local_count_files(fullfile(dirTS,'ERA5_Land'), ...
        {'ERA5_Land_*.csv'}) >= 484 ...
        && local_count_files(fullfile(dirTS,'EMO1_arc'), ...
        {'EMO1_*.csv'}) >= 484;
end

function root = local_find_bull_payload(stageDir,requiredName)

    root = '';
    if isfile(fullfile(stageDir,requiredName)) ...
            || isfolder(fullfile(stageDir,requiredName))
        root = stageDir;
        return
    end
    q = dir(fullfile(stageDir,'**',requiredName));
    if ~isempty(q)
        root = q(1).folder;
    end
end

function root = local_find_bull_timeseries_payload(stageDir)
%LOCAL_FIND_BULL_TIMESERIES_PAYLOAD Locate the CSV collection root.

    root = '';
    q = dir(fullfile(stageDir,'**','csv'));
    for i = 1:numel(q)
        if ~q(i).isdir
            continue
        end
        required = {'streamflow','AEMET','ERA5_Land','EMO1_arc'};
        candidates = unique({q(i).folder, ...
            fullfile(q(i).folder,q(i).name)},'stable');
        for j = 1:numel(candidates)
            present = cellfun(@(x) ...
                isfolder(fullfile(candidates{j},x)),required);
            if all(present)
                root = candidates{j};
                return
            end
        end
    end
end

% =============================================
function out = local_download_es(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_ES Download and install CAMELS-ES daily data.

    logFcn = local_ui_log(ui);

    if ~strcmp(stream,'daily')
        error('data_helpers:ESOnlyDaily', ...
            ['CAMELS-ES support is daily only. ' ...
            'Use stream = ''daily''.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_ES');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily','timeseries'));

    if ~local_ui_confirm(ui, ...
            'Download Spain data', ...
            ['CAMELS-ES is downloaded from Zenodo as one ZIP file. ' ...
            'The archive contains attribute CSV files and daily ' ...
            'hydrometeorological time series.' newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipName = 'CAMELS-ES.zip';
    zipFile = fullfile(downloadDir,zipName);
    url = ['https://zenodo.org/records/8428374/files/' ...
        'CAMELS-ES.zip?download=1'];
    expectedMD5 = '1a224276d5ce667ffbb6c3650e672fbc';

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Spain data', ...
        'Starting CAMELS-ES download ...');

    try
        if local_es_install_complete(dirD)
            logFcn(['CAMELS-ES appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true, ...
                'region','CAMELS_ES', ...
                'dirD',dirD, ...
                'dirM',fullfile(dirD,'daily','timeseries'), ...
                'dirQ',fullfile(dirD,'daily','timeseries'), ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(zipFile)
            try
                local_verify_md5(zipFile, ...
                    expectedMD5,logFcn);
                logFcn(['Using verified existing file: ' ...
                    zipFile]);
            catch
                logFcn(['Existing CAMELS-ES ZIP failed MD5; ' ...
                    'deleting and downloading again.']);
                try
                    delete(zipFile);
                catch
                end
            end
        end

        if ~isfile(zipFile)
            logFcn(['Downloading CAMELS-ES archive to: ' ...
                zipFile]);
            local_download_file_retry(url,zipFile, ...
                'CAMELS-ES',d, ...
                @() userCanceled, ...
                @(info) ui_progress(info), ...
                logFcn,3);
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        local_verify_md5(zipFile,expectedMD5,logFcn);

        unzipDir = fullfile(downloadDir, ...
            'CAMELS_ES_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        logFcn(['Unzipping CAMELS-ES archive to: ' ...
            unzipDir]);
        if ~isempty(d) ...
                && isvalid(d)
            try
                d.Indeterminate = true;
                d.Message = ['Unzipping CAMELS-ES archive. ' ...
                    'This can take a while ...'];
                drawnow limitrate nocallbacks
            catch
            end
        end

        unzip(zipFile,unzipDir);

        local_install_es_files(unzipDir,dirD,logFcn);

        okAttr = local_all_files_exist(dirD, ...
            local_es_metadata_files());
        okDaily = isfolder(fullfile(dirD, ...
            'daily','timeseries')) ...
            && local_count_files(fullfile(dirD, ...
            'daily','timeseries'), ...
            {'camelses_*.csv'}) > 0;

        out = struct('ok',okAttr && okDaily, ...
            'region','CAMELS_ES', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,'daily','timeseries'), ...
            'dirQ',fullfile(dirD,'daily','timeseries'), ...
            'archive',zipFile, ...
            'unzipDir',unzipDir);

        if out.ok
            nTs = local_count_files(fullfile(dirD, ...
                'daily','timeseries'), ...
                {'camelses_*.csv'});
            logFcn(sprintf(['CAMELS-ES install finished: ' ...
                '%d daily time-series CSV files.'],nTs));

            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);

            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-ES install finished, but ' ...
                'some required files appear to be missing.']);
        end

    catch ME
        try
            if exist('zipFile','var')
                local_cleanup_download_file(zipFile,logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder(unzipDir,logFcn);
            end
        catch
        end

        logFcn('CAMELS-ES install failed:');
        logFcn(ME.message);

        out = struct('ok',false, ...
            'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% =============================================
function out = local_download_gb(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_GB Download/install one CAMELS-GB v2 resolution.
%
% Daily installation:
%   attributes -> Data/CAMELS_GB/daily
%   hydromet   -> Data/CAMELS_GB/daily/timeseries
%
% Hourly installation:
%   attributes -> Data/CAMELS_GB/hourly
%   hydromet   -> Data/CAMELS_GB/hourly/timeseries
%
% The hourly installer does not download or write anything under daily.

    logFcn = local_ui_log(ui);
    stream = lower(strtrim(char(stream)));

    if ~any(strcmp(stream,{'daily','hourly'}))
        error('data_helpers:GBBadStream', ...
            'CAMELS-GB supports stream = daily or hourly.');
    end

    % Always recover the regional root first. cfg.dirD may currently point
    % to Data/CAMELS_GB/daily or Data/CAMELS_GB/hourly.
    dirRoot = local_cfg_dirD(cfg,'CAMELS_GB');
    dirStream = fullfile(dirRoot,stream);
    dirTS = fullfile(dirStream,'timeseries');

    local_mkdir(dirRoot);
    local_mkdir(dirStream);
    local_mkdir(dirTS);

    if strcmp(stream,'hourly')
        titleText = 'Download Great Britain hourly data';
        msg = [ ...
            'CAMELS-GB v2 hourly files will be installed as follows:' ...
            newline newline ...
            'Attributes:' newline ...
            '  Data\CAMELS_GB\hourly' newline ...
            'Hydromet CSV files:' newline ...
            '  Data\CAMELS_GB\hourly\timeseries' newline newline ...
            'No daily files will be downloaded or modified.' newline ...
            'Continue?'];
    else
        titleText = 'Download Great Britain daily data';
        msg = [ ...
            'CAMELS-GB v2 daily files will be installed as follows:' ...
            newline newline ...
            'Attributes:' newline ...
            '  Data\CAMELS_GB\daily' newline ...
            'Hydromet CSV files:' newline ...
            '  Data\CAMELS_GB\daily\timeseries' newline ...
            'Continue?'];
    end

    if ~local_ui_confirm(ui,titleText,msg)
        out = struct('ok',false,'canceled',true);
        return
    end

    base = [ ...
        'https://catalogue.ceh.ac.uk/datastore/eidchub/' ...
        '9a46d428-958f-4ac1-86eb-94eee70c0955/'];

    attrUrl = [base 'Catchment_Attributes/'];

    if strcmp(stream,'hourly')
        tsUrl = [base ...
            'Catchment_Timeseries/hydro-meteorological/hourly/'];
        prefix = 'camels_gb_v2_hydromet_hourly_timeseries_';
        expectedCount = 664;
    else
        tsUrl = [base ...
            'Catchment_Timeseries/hydro-meteorological/daily/'];
        prefix = 'camels_gb_v2_hydromet_daily_timeseries_';
        expectedCount = 669;
    end

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        ['Downloading Great Britain ' stream ' data'], ...
        ['Preparing CAMELS-GB v2 ' stream ' installation ...']);

    try
        attrFiles = local_gb_v2_attribute_files();
        tsFiles = local_gb_remote_csv_list(tsUrl,prefix);

        if numel(tsFiles) < expectedCount
            error('data_helpers:GBIndexIncomplete', ...
                ['The CAMELS-GB %s directory index returned %d CSV ' ...
                'files; expected at least %d.'], ...
                stream,numel(tsFiles),expectedCount);
        end

        nTotal = numel(attrFiles) + numel(tsFiles);
        iDone = 0;

        % Attributes always go to the selected stream root.
        for i = 1:numel(attrFiles)
            if userCanceled
                break
            end
            fname = attrFiles{i};
            dstFile = fullfile(dirStream,fname);

            iDone = local_gb_download_one( ...
                [attrUrl fname],dstFile, ...
                sprintf('GB %s attribute %s',stream,fname), ...
                iDone,nTotal,d,@() userCanceled,@ui_progress,logFcn);
        end

        % Hydromet files always go to the selected stream/timeseries folder.
        for i = 1:numel(tsFiles)
            if userCanceled
                break
            end
            fname = tsFiles{i};
            dstFile = fullfile(dirTS,fname);

            iDone = local_gb_download_one( ...
                [tsUrl fname],dstFile, ...
                sprintf('GB %s file %d/%d', ...
                    stream,i,numel(tsFiles)), ...
                iDone,nTotal,d,@() userCanceled,@ui_progress,logFcn);
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        okAttr = local_all_files_exist(dirStream,attrFiles);
        okTS = local_count_files(dirTS,{[prefix '*.csv']}) ...
            >= expectedCount;

        out = struct( ...
            'ok',okAttr && okTS, ...
            'region','CAMELS_GB', ...
            'dirD',dirStream, ...
            'dirM',dirTS, ...
            'dirQ',dirTS, ...
            'archive','', ...
            'unzipDir','');

        if out.ok
            logFcn(sprintf([ ...
                'CAMELS-GB v2 %s installation completed.' newline ...
                'Attributes: %s' newline ...
                'Hydromet: %s'], ...
                stream,dirStream,dirTS));
        else
            logFcn(sprintf([ ...
                'CAMELS-GB %s installation finished, but required ' ...
                'files are still missing under %s.'],stream,dirStream));
        end

    catch ME
        logFcn(['CAMELS-GB v2 ' stream ' installation failed:']);
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

function out = local_download_fr(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_FR Download and install CAMELS-FR daily data.
    logFcn = local_ui_log(ui);
    if ~strcmp(stream,'daily')
        error('data_helpers:FRHourlyUnsupported', ...
            ['Hourly data are not ' ...
            'available for CAMELS-FR.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_FR');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily','timeseries'));

    if ~local_ui_confirm(ui, ...
            'Download France data', ...
            ['CAMELS-FR requires two ' ...
            'data.gouv.fr ZIP files: ' ...
            'attributes and daily time series. ' ...
            'The download may be large.' newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    files = local_camels_fr_files();
    downloadDir = local_default_download_dir();
    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading France data', ...
        'Starting CAMELS-FR download ...');

    try
        if local_fr_install_complete(dirD)
            logFcn(['CAMELS-FR appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true, ...
                'region','CAMELS_FR', ...
                'dirD',dirD, ...
                'dirM',fullfile(dirD,'daily','timeseries'), ...
                'dirQ',fullfile(dirD,'daily','timeseries'), ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        for i = 1:numel(files)
            f = files(i);
            zipFile = fullfile(downloadDir,f.fileName);
            files(i).localFile = zipFile;

            if isfile(zipFile)
                try
                    local_verify_md5(zipFile,f.md5,logFcn);
                    logFcn(['Using verified ' ...
                        'existing file: ' zipFile]);
                catch
                    logFcn(['Existing CAMELS-FR ' ...
                        'ZIP failed MD5; ' ...
                        'deleting and ' ...
                        'downloading again.']);
                    delete(zipFile);
                end
            end

            if ~isfile(zipFile)
                logFcn(sprintf('Downloading %s to: %s', ...
                    f.label,zipFile));
                local_download_file_retry(f.url,zipFile, ...
                    f.label,d, ...
                    @() userCanceled, ...
                    @(info) ui_progress(info,i,numel(files)), ...
                    logFcn,3);
            end

            if userCanceled
                out = struct('ok',false,'canceled',true);
                local_close_progress(d);
                return
            end
            local_verify_md5(zipFile,f.md5,logFcn);
        end

        unzipDir = fullfile(downloadDir, ...
            'CAMELS_FR_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        for i = 1:numel(files)
            f = files(i);
            dst = fullfile(unzipDir, ...
                erase(f.fileName,'.zip'));
            local_mkdir(dst);
            logFcn(['Unzipping ' f.fileName ' to: ' dst]);
            if ~isempty(d) && isvalid(d)
                try
                    d.Indeterminate = true;
                    d.Message = ...
                        ['Unzipping ' f.fileName ' ...'];
                    drawnow limitrate nocallbacks
                catch
                end
            end
            unzip(files(i).localFile,dst);
            files(i).unzipDir = dst;
        end

        local_install_fr_files(unzipDir,dirD,logFcn);

        okAttr = local_all_files_exist(dirD, ...
            local_fr_metadata_files());
        okDaily = isfolder(fullfile(dirD, ...
            'daily','timeseries')) ...
            && local_count_files( ...
            fullfile(dirD,'daily','timeseries'), ...
            {'CAMELS_FR_tsd_*.csv'}) >= 50;

        out = struct('ok',okAttr && okDaily, ...
            'region','CAMELS_FR', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,'daily','timeseries'), ...
            'dirQ',fullfile(dirD,'daily','timeseries'), ...
            'archive','', ...
            'unzipDir','');

        if out.ok
            logFcn('CAMELS-FR install finished.');
            for i = 1:numel(files)
                local_cleanup_download_file( ...
                    files(i).localFile,logFcn);
            end
            local_cleanup_download_folder( ...
                unzipDir,logFcn);
        else
            logFcn(['CAMELS-FR install ' ...
                'finished, but some ' ...
                'required files appear ' ...
                'to be missing.']);
        end

    catch ME
        try
            if exist('files','var')
                for i = 1:numel(files)
                    if isfield(files(i), ...
                            'localFile')
                        local_cleanup_download_file( ...
                            files(i).localFile,logFcn);
                    end
                end
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder( ...
                    unzipDir,logFcn);
            end
        catch
        end
        logFcn('CAMELS-FR install failed:');
        logFcn(ME.message);
        out = struct('ok', ...
            false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info,iFile,nFile)
        if nargin >= 2
            if isfield(info,'frac') ...
                    && isfinite(info.frac)
                info.frac = ((iFile - 1) ...
                    + info.frac) / nFile;
            else
                info.frac = (iFile - 1) / nFile;
            end
            if isfield(info,'label')
                info.label = sprintf('%s (%d/%d)', ...
                    info.label,iFile,nFile);
            else
                info.label = sprintf('File %d/%d', ...
                    iFile,nFile);
            end
        end
        [userCanceled,lastUI] = ...
            local_update_progress_dialog(d, ...
            info,lastUI,userCanceled);
    end
end



% =============================================
function out = local_download_fi(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_FI Download and install CAMELS-FI daily data.
    logFcn = local_ui_log(ui);

    if ~strcmp(stream,'daily')
        error('data_helpers:FIOnlyDaily', ...
            ['CAMELS-FI support is daily only. ' ...
            'Use stream = ''daily''.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_FI');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily','timeseries'));

    if ~local_ui_confirm(ui, ...
            'Download Finland data', ...
            ['CAMELS-FI is downloaded from Zenodo as one ZIP file. ' ...
            'The archive contains static attributes, catchment ' ...
            'boundaries, and daily hydrometeorological time series.' newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipName = 'CAMELS-FI.zip';
    zipFile = fullfile(downloadDir,zipName);
    url = ['https://zenodo.org/records/20225368/files/' ...
        'CAMELS-FI.zip?download=1'];
    expectedMD5 = 'f50bf2d972f42b6fc4db690ce201482f';

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Finland data', ...
        'Starting CAMELS-FI download ...');

    try
        if local_fi_install_complete(dirD)
            logFcn(['CAMELS-FI appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true, ...
                'region','CAMELS_FI', ...
                'dirD',dirD, ...
                'dirM',fullfile(dirD,'daily','timeseries'), ...
                'dirQ',fullfile(dirD,'daily','timeseries'), ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(zipFile)
            try
                local_verify_md5(zipFile, ...
                    expectedMD5,logFcn);
                logFcn(['Using verified existing file: ' ...
                    zipFile]);
            catch
                logFcn(['Existing CAMELS-FI ZIP failed MD5; ' ...
                    'deleting and downloading again.']);
                delete(zipFile);
            end
        end

        if ~isfile(zipFile)
            logFcn(['Downloading CAMELS-FI archive to: ' ...
                zipFile]);
            local_download_file_retry(url,zipFile, ...
                'CAMELS-FI',d, ...
                @() userCanceled, ...
                @(info) ui_progress(info), ...
                logFcn,3);
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        local_verify_md5(zipFile,expectedMD5,logFcn);

        unzipDir = fullfile(downloadDir, ...
            'CAMELS_FI_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        logFcn(['Unzipping CAMELS-FI archive to: ' ...
            unzipDir]);
        if ~isempty(d) ...
                && isvalid(d)
            try
                d.Indeterminate = true;
                d.Message = ['Unzipping CAMELS-FI archive. ' ...
                    'This can take a while ...'];
                drawnow limitrate nocallbacks
            catch
            end
        end
        unzip(zipFile,unzipDir);

        local_install_fi_files(unzipDir,dirD,logFcn);

        okAttr = local_all_files_exist(dirD, ...
            local_fi_metadata_files());
        okDaily = isfolder(fullfile(dirD, ...
            'daily','timeseries')) ...
            && local_count_files(fullfile(dirD, ...
            'daily','timeseries'), ...
            {'CAMELS_FI_hydromet_timeseries_*.csv'}) > 0;

        out = struct('ok',okAttr && okDaily, ...
            'region','CAMELS_FI', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,'daily','timeseries'), ...
            'dirQ',fullfile(dirD,'daily','timeseries'), ...
            'archive',zipFile, ...
            'unzipDir',unzipDir);

        if out.ok
            nTs = local_count_files(fullfile(dirD, ...
                'daily','timeseries'), ...
                {'CAMELS_FI_hydromet_timeseries_*.csv'});
            logFcn(sprintf(['CAMELS-FI install finished: ' ...
                '%d daily time-series CSV files.'],nTs));
            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-FI install finished, but ' ...
                'some required files appear to be missing.']);
        end

    catch ME
        try
            if exist('zipFile','var')
                local_cleanup_download_file(zipFile,logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder(unzipDir,logFcn);
            end
        catch
        end
        logFcn('CAMELS-FI install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end


% ==============================================
function out = local_download_lux(cfg,stream,ui)
% ==============================================
%LOCAL_DOWNLOAD_LUX Download/install CAMELS-LUX daily/hourly/15-min data.
    logFcn = local_ui_log(ui);

    stream = lower(strtrim(char(stream)));

    if ~any(strcmp(stream,{'daily','hourly','15min'}))
        error('data_helpers:LUXBadStream', ...
            ['CAMELS-LUX supports stream = daily, ' ...
            'hourly, or 15min.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_LUX');
    dirT = fullfile(dirD,stream,'timeseries');
    local_mkdir(dirD);
    local_mkdir(dirT);

    allStreams = {'daily','hourly','15min'};
    installed = false(size(allStreams));
    for i = 1:numel(allStreams)
        installed(i) = local_lux_install_complete( ...
            dirD,allStreams{i});
    end
    if all(installed)
        logFcn(['All CAMELS-LUX temporal resolutions are already ' ...
            'installed.']);
        out = struct('ok',true,'region','CAMELS_LUX', ...
            'dirD',dirD,'dirM',dirT,'dirQ',dirT, ...
            'archive','','unzipDir','', ...
            'installedStreams',{allStreams});
        return
    end

    [installStreams,canceled] = local_select_lux_resolutions( ...
        ui,allStreams,installed,stream);
    if canceled || isempty(installStreams)
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipName = 'CAMELS-LUX.zip';
    zipFile = fullfile(downloadDir,zipName);
    url = ['https://zenodo.org/records/18776538/files/' ...
        'CAMELS-LUX.zip?download=1'];
    expectedMD5 = '6c4a14a0feed08382a6b565a798d8fdc';

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Luxembourg data', ...
        'Starting the shared CAMELS-LUX download ...');

    try
        if isfile(zipFile)
            try
                local_verify_md5(zipFile, ...
                    expectedMD5,logFcn);
                logFcn(['Using verified existing file: ' ...
                    zipFile]);
            catch
                logFcn(['Existing CAMELS-LUX ZIP failed MD5; ' ...
                    'deleting and downloading again.']);
                delete(zipFile);
            end
        end

        if ~isfile(zipFile)
            logFcn(['Downloading CAMELS-LUX archive to: ' ...
                zipFile]);
            local_download_file_retry(url,zipFile, ...
                ['CAMELS-LUX ' stream],d, ...
                @() userCanceled, ...
                @(info) ui_progress(info), ...
                logFcn,3);
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        local_verify_md5(zipFile,expectedMD5,logFcn);

        unzipDir = fullfile(downloadDir, ...
            'CAMELS_LUX_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        logFcn(['Unzipping CAMELS-LUX archive to: ' ...
            unzipDir]);
        if ~isempty(d) ...
                && isvalid(d)
            try
                d.Indeterminate = true;
                d.Message = ['Unzipping CAMELS-LUX archive. ' ...
                    'This can take a while ...'];
                drawnow limitrate nocallbacks
            catch
            end
        end
        unzip(zipFile,unzipDir);

        for i = 1:numel(installStreams)
            installStream = installStreams{i};
            try
                d.Message = ['Installing CAMELS-LUX ' ...
                    installStream ' data ...'];
                drawnow limitrate nocallbacks
            catch
            end
            local_install_lux_files( ...
                unzipDir,dirD,logFcn,installStream);
            if ~local_lux_install_complete(dirD,installStream)
                error('data_helpers:LUXIncompleteInstall', ...
                    'CAMELS-LUX %s installation is incomplete.', ...
                    installStream);
            end
        end

        ok = all(cellfun(@(s)local_lux_install_complete(dirD,s), ...
            installStreams));

        out = struct('ok',ok, ...
            'region','CAMELS_LUX', ...
            'dirD',dirD, ...
            'dirM',dirT, ...
            'dirQ',dirT, ...
            'archive',zipFile, ...
            'unzipDir',unzipDir, ...
            'installedStreams',{installStreams});

        if out.ok
            for i = 1:numel(installStreams)
                installStream = installStreams{i};
                installDir = fullfile( ...
                    dirD,installStream,'timeseries');
                pats = local_lux_timeseries_patterns(installStream);
                nTs = local_count_files(installDir,pats);
                logFcn(sprintf(['CAMELS-LUX %s install finished: ' ...
                    '%d time-series CSV files.'],installStream,nTs));
            end
            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-LUX install finished, but ' ...
                'some required files appear to be missing.']);
        end

    catch ME
        try
            if exist('zipFile','var')
                local_cleanup_download_file(zipFile,logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder(unzipDir,logFcn);
            end
        catch
        end
        logFcn('CAMELS-LUX install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

function [streams,canceled] = local_select_lux_resolutions( ...
        ui,allStreams,installed,requested)
%LOCAL_SELECT_LUX_RESOLUTIONS Select missing products in one shared ZIP.

    streams = {};
    canceled = true;
    hasUI = isstruct(ui) ...
        && isfield(ui,'fig') ...
        && ~isempty(ui.fig) ...
        && isvalid(ui.fig);
    if ~hasUI
        j = find(strcmp(allStreams,requested),1);
        if ~isempty(j) && ~installed(j)
            streams = {requested};
        end
        canceled = false;
        return
    end

    displayNames = {'Daily','Hourly','15-minute'};
    expected = [50 51 51];
    dlg = uifigure('Name','Install CAMELS-LUX data', ...
        'WindowStyle','modal','Position',[100 100 500 255], ...
        'CloseRequestFcn',@cancelSelection);
    g = uigridlayout(dlg,[6 2], ...
        'RowHeight',{44,28,28,28,22,34}, ...
        'ColumnWidth',{'1x',150}, ...
        'RowSpacing',5,'ColumnSpacing',10, ...
        'Padding',[16 12 16 12]);
    intro = uilabel(g,'Text',[ ...
        'The same Zenodo archive contains all three resolutions. ' ...
        'Select the missing products to install in one pass.'], ...
        'WordWrap','on');
    intro.Layout.Row = 1;
    intro.Layout.Column = [1 2];

    cb = gobjects(numel(allStreams),1);
    for i = 1:numel(allStreams)
        if installed(i)
            suffix = sprintf(' — installed (%d basins)',expected(i));
        else
            suffix = sprintf(' — not installed (%d basins)',expected(i));
        end
        cb(i) = uicheckbox(g,'Text',[displayNames{i} suffix], ...
            'Value',~installed(i));
        cb(i).Layout.Row = i + 1;
        cb(i).Layout.Column = [1 2];
        if installed(i)
            cb(i).Enable = 'off';
        end
    end

    note = uilabel(g,'Text', ...
        'Installed resolutions are retained and are not copied again.', ...
        'FontAngle','italic');
    note.Layout.Row = 5;
    note.Layout.Column = [1 2];
    bCancel = uibutton(g,'Text','Cancel', ...
        'ButtonPushedFcn',@cancelSelection);
    bCancel.Layout.Row = 6;
    bCancel.Layout.Column = 1;
    bInstall = uibutton(g,'Text','Install selected', ...
        'ButtonPushedFcn',@acceptSelection);
    bInstall.Layout.Row = 6;
    bInstall.Layout.Column = 2;
    movegui(dlg,'center');
    uiwait(dlg);

    function acceptSelection(~,~)
        streams = {};
        for q = 1:numel(allStreams)
            if ~installed(q) && cb(q).Value
                streams{end+1} = allStreams{q}; %#ok<AGROW>
            end
        end
        canceled = false;
        if isvalid(dlg)
            delete(dlg);
        end
    end

    function cancelSelection(~,~)
        streams = {};
        canceled = true;
        if isvalid(dlg)
            delete(dlg);
        end
    end
end

% ==============================================
function out = local_download_il(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_IL Download and install daily Israel Caravan extension.

    logFcn = local_ui_log(ui);
    if ~strcmp(stream,'daily')
        error('data_helpers:ILOnlyDaily', ...
            'CAMELS-IL is a daily-only data package.');
    end

    dirRoot = local_cfg_dirD(cfg,'CAMELS_IL');
    dirD = fullfile(dirRoot,'daily');
    dirTS = fullfile(dirD,'timeseries');
    local_mkdir(dirRoot);
    local_mkdir(dirD);
    local_mkdir(dirTS);

    if ~local_ui_confirm(ui,'Download Israel data', ...
            ['Download Caravan_extension_Israel_Ver4.zip from Zenodo? ' ...
             'The archive contains 95 daily catchment time series and ' ...
             'three attribute files.'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipFile = fullfile(downloadDir,'Caravan_extension_Israel_Ver4.zip');
    unzipDir = fullfile(downloadDir,'CAMELS_IL_download');
    url = ['https://zenodo.org/records/15181680/files/' ...
        'Caravan_extension_Israel_Ver4.zip?download=1'];
    expectedMD5 = '4b591229905677413ec5d305a84de127';

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui,'Downloading Israel data', ...
        'Starting CAMELS-IL download ...');

    try
        if local_il_install_complete(dirRoot)
            logFcn('CAMELS-IL daily data appear complete; skipping download.');
            out = struct('ok',true,'region','CAMELS_IL', ...
                'dirD',dirD,'dirM',dirTS,'dirQ',dirTS, ...
                'archive','','unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(zipFile)
            try
                local_verify_md5(zipFile,expectedMD5,logFcn);
                logFcn(['Using verified existing file: ' zipFile]);
            catch
                delete(zipFile);
            end
        end
        if ~isfile(zipFile)
            local_download_file_retry(url,zipFile,'CAMELS-IL',d, ...
                @() userCanceled,@(info) ui_progress(info),logFcn,3);
        end
        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end
        local_verify_md5(zipFile,expectedMD5,logFcn);

        if isfolder(unzipDir), rmdir(unzipDir,'s'); end
        local_mkdir(unzipDir);
        if ~isempty(d) && isvalid(d)
            d.Indeterminate = true;
            d.Message = 'Unzipping CAMELS-IL archive ...';
            drawnow limitrate nocallbacks
        end
        unzip(zipFile,unzipDir);

        payload = local_find_il_payload(unzipDir);
        if isempty(payload)
            error('data_helpers:ILMissingPayload', ...
                'Could not locate Caravan_extension_Israel_Ver4 payload.');
        end
        srcAttr = fullfile(payload,'attributes','il');
        srcTS = fullfile(payload,'timeseries','csv','il');

        attrs = {'attributes_caravan_il.csv', ...
            'attributes_hydroatlas_il.csv', ...
            'attributes_other_il.csv'};
        for k = 1:numel(attrs)
            f = fullfile(srcAttr,attrs{k});
            if ~isfile(f)
                error('data_helpers:ILMissingAttribute', ...
                    'Missing archive file: %s',f);
            end
            copyfile(f,fullfile(dirD,attrs{k}),'f');
        end

        D = dir(fullfile(srcTS,'il_*.csv'));
        if numel(D) < 95
            error('data_helpers:ILMissingTimeseries', ...
                'Expected 95 CAMELS-IL CSV files; found %d.',numel(D));
        end
        for k = 1:numel(D)
            copyfile(fullfile(D(k).folder,D(k).name), ...
                fullfile(dirTS,D(k).name),'f');
            if mod(k,5)==0 || k==numel(D)
                if ~isempty(d) && isvalid(d)
                    d.Message = sprintf( ...
                        'Installing CAMELS-IL files %d/%d ...',k,numel(D));
                    drawnow limitrate nocallbacks
                end
            end
        end

        ok = local_il_install_complete(dirRoot);
        out = struct('ok',ok,'region','CAMELS_IL', ...
            'dirD',dirD,'dirM',dirTS,'dirQ',dirTS, ...
            'archive',zipFile,'unzipDir',unzipDir);
        if ok
            logFcn(sprintf(['CAMELS-IL install finished: %d daily ' ...
                'time-series files and three attribute files.'],numel(D)));
            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        end
    catch ME
        logFcn('CAMELS-IL install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end
    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

function tf = local_il_install_complete(dirRoot)
    d = fullfile(dirRoot,'daily');
    tf = local_all_files_exist(d,{ ...
        'attributes_caravan_il.csv', ...
        'attributes_hydroatlas_il.csv', ...
        'attributes_other_il.csv'}) ...
        && local_count_files(fullfile(d,'timeseries'),{'il_*.csv'}) >= 95;
end

function payload = local_find_il_payload(root)
    payload = '';
    candidates = {root,fullfile(root,'Caravan_extension_Israel_Ver4')};
    for k = 1:numel(candidates)
        if isfolder(fullfile(candidates{k},'attributes','il')) ...
                && isfolder(fullfile(candidates{k}, ...
                'timeseries','csv','il'))
            payload = candidates{k};
            return
        end
    end
    D = dir(fullfile(root,'**','attributes','il'));
    if ~isempty(D)
        p = fileparts(fileparts(fullfile(D(1).folder,D(1).name)));
        if isfolder(fullfile(p,'timeseries','csv','il'))
            payload = p;
        end
    end
end

function out = local_download_ind(cfg,stream,ui)
% ==============================================
%LOCAL_DOWNLOAD_IND Download and install CAMELS-IND daily data.
    logFcn = local_ui_log(ui);
    if ~strcmp(stream,'daily')
        error('data_helpers:INDHourlyUnsupported', ...
            'Hourly data are not available for CAMELS-IND.');
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_IND');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily', ...
        'catchment_mean_forcings'));

    if ~local_ui_confirm(ui, ...
            'Download India data', ...
            ['CAMELS-IND is downloaded from Zenodo as one ZIP file. ' ...
            'The download may be large.' newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipName = 'CAMELS_IND_All_Catchments.zip';
    zipFile = fullfile(downloadDir,zipName);
    url = ['https://zenodo.org/records/14999580/files/' ...
        'CAMELS_IND_All_Catchments.zip?download=1'];
    expectedMD5 = 'b0b48fbcbc0a0e86e01bebcf152984a2';

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading India data', ...
        'Starting CAMELS-IND download ...');

    try
        if local_ind_install_complete(dirD)
            logFcn(['CAMELS-IND appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true, ...
                'region','CAMELS_IND', ...
                'dirD',dirD, ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(zipFile)
            try
                local_verify_md5(zipFile, ...
                    expectedMD5,logFcn);
                logFcn(['Using verified existing file: ' ...
                    zipFile]);
            catch
                logFcn(['Existing CAMELS-IND ZIP failed MD5; ' ...
                    'deleting and downloading again.']);
                delete(zipFile);
            end
        end

        if ~isfile(zipFile)
            logFcn(['Downloading CAMELS-IND archive to: ' ...
                zipFile]);
            local_download_file_retry(url,zipFile, ...
                'CAMELS-IND',d, ...
                @() userCanceled, ...
                @(info) ui_progress(info), ...
                logFcn,3);
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        local_verify_md5(zipFile,expectedMD5,logFcn);

        unzipDir = fullfile(downloadDir, ...
            'CAMELS_IND_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        logFcn(['Unzipping CAMELS-IND archive to: ' ...
            unzipDir]);
        if ~isempty(d) ...
                && isvalid(d)
            try
                d.Indeterminate = true;
                d.Message = ['Unzipping CAMELS-IND archive. ' ...
                    'This can take a while ...'];
                drawnow limitrate nocallbacks
            catch
            end
        end
        unzip(zipFile,unzipDir);

        dataDir = local_find_child_dir(unzipDir, ...
            'CAMELS_IND_All_Catchments');
        if isempty(dataDir)
            dataDir = local_find_first_payload_dir( ...
                unzipDir,'CAMELS_IND_All_Catchments');
        end
        if isempty(dataDir) ...
                || ~isfolder(dataDir)
            error('data_helpers:INDMissingDataDir', ...
                ['Could not find CAMELS_IND_All_Catchments ' ...
                'inside the archive.']);
        end

        local_install_ind_files(unzipDir,dirD,logFcn);

        okAttr = local_all_files_exist(dirD, ...
            local_ind_metadata_files());
        okDaily = isfile(fullfile(dirD,'daily', ...
            'streamflow_timeseries', ...
            'streamflow_observed.csv')) ...
            && local_count_files(fullfile(dirD,'daily', ...
            'catchment_mean_forcings'),{'*.csv'}) >= 200;

        out = struct('ok',okAttr && okDaily, ...
            'region','CAMELS_IND', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,'daily', ...
            'catchment_mean_forcings'), ...
            'dirQ',fullfile(dirD,'daily', ...
            'streamflow_timeseries'), ...
            'archive',zipFile, ...
            'unzipDir',unzipDir);

        if out.ok
            logFcn('CAMELS-IND install finished.');
            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-IND install finished, but ' ...
                'some required files appear to be missing.']);
        end

    catch ME
        try
            if exist('zipFile','var')
                local_cleanup_download_file(zipFile,logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder(unzipDir,logFcn);
            end
        catch
        end
        logFcn('CAMELS-IND install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end


% =============================================
function out = local_download_kr(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_KR Download and install hourly CAMELSH-KR from Zenodo.

    logFcn = local_ui_log(ui);
    if ~strcmp(stream,'hourly')
        error('data_helpers:KRBadStream', ...
            'CAMELSH-KR currently supports hourly data only.');
    end

    dirD = local_cfg_dirD(cfg,'CAMELSH_KR');
    dirT = fullfile(dirD,'hourly','timeseries');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'hourly'));

    if ~local_ui_confirm(ui,'Download South Korea data', ...
            sprintf(['CAMELSH-KR contains 178 hourly basin files. ' ...
            'The time-series ZIP is 3.06 GB; forcing and streamflow ' ...
            'are stored together in each CSV.' newline ...
            'Five smaller attribute files will also be downloaded. ' ...
            'Continue?']))
        out = struct('ok',false,'canceled',true);
        return
    end

    files = local_camels_kr_files();
    downloadDir = local_default_download_dir();
    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui,'Downloading South Korea data', ...
        'Starting CAMELSH-KR download ...');
    zipFile = '';
    try
        if local_kr_install_complete(dirD)
            out = struct('ok',true,'region','CAMELSH_KR', ...
                'dirD',dirD,'dirM',dirT,'dirQ',dirT, ...
                'archive','','unzipDir','');
            local_close_progress(d);
            return
        end

        for i = 1:numel(files)
            f = files(i);
            if f.isZip
                targetFile = fullfile(downloadDir,f.fileName);
                zipFile = targetFile;
            else
                targetFile = fullfile(dirD,f.fileName);
            end

            if isfile(targetFile)
                try
                    local_verify_md5(targetFile,f.md5,logFcn);
                catch
                    delete(targetFile);
                end
            end
            if ~isfile(targetFile)
                logFcn(sprintf('Downloading %s ...',f.label));
                local_download_file_retry(f.url,targetFile,f.label,d, ...
                    @() userCanceled,@ui_progress,logFcn,3);
            end
            local_verify_md5(targetFile,f.md5,logFcn);
        end

        if ~local_has_kr_timeseries(dirD)
            if isempty(zipFile) || ~isfile(zipFile)
                error('data_helpers:KRMissingArchive', ...
                    'The verified CAMELSH-KR time-series archive is missing.');
            end
            if ~isempty(d) && isvalid(d)
                d.Indeterminate = true;
                d.Message = ['Extracting 178 hourly CAMELSH-KR ' ...
                    'basin files ...'];
                drawnow limitrate nocallbacks
            end
            unzip(zipFile,fullfile(dirD,'hourly'));
        end

        if ~local_kr_install_complete(dirD)
            error('data_helpers:KRIncompleteInstall', ...
                ['CAMELSH-KR installation is incomplete. Expected five ' ...
                'attribute files and 178 hourly basin CSV files.']);
        end

        logFcn(['Installed CAMELSH-KR hourly forcing and streamflow ' ...
            'for 178 basins in: ' dirT]);
        out = struct('ok',true,'region','CAMELSH_KR', ...
            'dirD',dirD,'dirM',dirT,'dirQ',dirT, ...
            'archive',zipFile,'unzipDir',fullfile(dirD,'hourly'));
    catch ME
        logFcn(['CAMELSH-KR install failed: ' ME.message]);
        out = struct('ok',false,'error',ME.message);
    end
    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% =============================================
function out = local_download_za(cfg,stream,ui)
%LOCAL_DOWNLOAD_ZA Shared GRDC-Caravan installer, requested country ZA.
    out=local_download_grdc_africa(cfg,stream,ui,'ZA');
end

function out = local_download_na(cfg,stream,ui)
%LOCAL_DOWNLOAD_NA Shared GRDC-Caravan installer, requested country NA.
    out=local_download_grdc_africa(cfg,stream,ui,'NA');
end

function out = local_download_grdc_africa(cfg,stream,ui,requested)
%LOCAL_DOWNLOAD_GRDC_AFRICA Install selected GRDC country subsets.
    if ~strcmp(stream,'daily')
        error('data_helpers:GRDCAfricaBadStream', ...
            'The supported GRDC-Caravan country subsets are daily.');
    end
    logFcn=local_ui_log(ui);
    req=local_grdc_africa_spec(requested);
    allCodes={'AR','BE','EE','IE','JM','NA','NO','PR','ZA'};
    specs=cell(size(allCodes)); installed=false(size(allCodes));
    for ispec=1:numel(allCodes)
        specs{ispec}=local_grdc_africa_spec(allCodes{ispec});
        installed(ispec)=local_grdc_africa_installed(cfg,specs{ispec});
    end
    [installCodes,canceledSelection]=local_select_grdc_regions( ...
        ui,specs,installed,requested);
    if canceledSelection
        out=struct('ok',false,'canceled',true); return
    end
    if isempty(installCodes)
        dirD=local_cfg_dirD(cfg,req.region); dirT=fullfile(dirD,'daily','timeseries');
        out=struct('ok',true,'region',req.region,'dirD',dirD,'dirM',dirT,'dirQ',dirT);
        return
    end

    names=cell(size(installCodes));
    for iname=1:numel(installCodes)
        sname=local_grdc_africa_spec(installCodes{iname});
        names{iname}=sname.name;
    end
    logFcn(['Selected GRDC-Caravan regions: ' strjoin(names,', ')]);

    downloadDir=local_default_download_dir();
    canonical=fullfile(downloadDir,'GRDC_Caravan_extension_csv.zip');
    userCopy=fullfile(downloadDir,'GRDC_Caravan_extension_csv - Copy.zip');
    if isfile(canonical)
        zipFile=canonical; deleteArchive=true;
    elseif isfile(userCopy)
        zipFile=userCopy; deleteArchive=false;
        logFcn(['Reusing retained GRDC-Caravan archive: ' zipFile]);
    else
        zipFile=canonical; deleteArchive=true;
    end
    url=['https://zenodo.org/records/15349031/files/' ...
        'GRDC_Caravan_extension_csv.zip?download=1'];
    tmp=fullfile(local_cfg_root(cfg),'Data','_grdc_caravan_extract');
    % Remove obsolete country-local extraction folders left by older
    % versions of the South Africa/Namibia installer.
    for iclean=1:numel(specs)
        legacy=fullfile(local_cfg_dirD(cfg,specs{iclean}.region), ...
            '_grdc_caravan_extract');
        if isfolder(legacy), local_za_remove_tree(legacy,logFcn,true); end
    end
    d=local_progress_dialog(ui,'Installing GRDC-Caravan regions', ...
        'Preparing the GRDC-Caravan archive ...');
    canceled=false; lastUI=tic;
    try
        if ~isfile(zipFile)
            local_download_file_retry(url,zipFile,'GRDC-Caravan CSV',d, ...
                @()canceled,@progress,logFcn,3);
        end
        local_verify_md5(zipFile,'2689c2bff8807f53c3de127827a3cd16',logFcn);
        if isfolder(tmp)
            local_za_remove_tree(tmp,logFcn,false);
            if isfolder(tmp)
                error('data_helpers:GRDCStaleExtraction', ...
                    'Cannot clear previous extraction directory: %s',tmp);
            end
        end
        local_mkdir(tmp);
        if ~isempty(d)&&isvalid(d)
            d.Indeterminate=true;
            d.Message='Extracting GRDC-Caravan once for selected countries ...';
            drawnow;
        end
        unzip(zipFile,tmp);
        attrNames={'attributes_other_grdc.csv','attributes_caravan_grdc.csv', ...
            'attributes_hydroatlas_grdc.csv','attributes_additional_grdc.csv'};
        src=cell(size(attrNames));
        for j=1:numel(attrNames), src{j}=local_za_find_one(tmp,attrNames{j}); end
        O=readtable(src{1},'VariableNamingRule','preserve');
        key=local_za_var(O,{'gauge_id','gaugeid','station_id'});
        country=local_za_var(O, ...
            {'country','country_name','country_code','country_iso3','country_id'});
        countryValues=strip(string(O.(country)));
        allCsv=dir(fullfile(tmp,'**','*.csv'));

        for ic=1:numel(installCodes)
            spec=local_grdc_africa_spec(installCodes{ic});
            take=false(height(O),1);
            for ia=1:numel(spec.aliases)
                take=take|strcmpi(countryValues,string(spec.aliases{ia}));
            end
            sourceIds=strip(string(O.(key)(take)));
            sourceIds=sort(unique(regexprep(sourceIds,'\.0+$','')));
            if numel(sourceIds)~=spec.count
                error('data_helpers:GRDCCountryCount', ...
                    'Expected %d %s basins, found %d.', ...
                    spec.count,spec.name,numel(sourceIds));
            end
            ids=regexprep(upper(sourceIds),'^GRDC[_-]','');
            dirD=local_cfg_dirD(cfg,spec.region);
            dirT=fullfile(dirD,'daily','timeseries');
            local_mkdir(dirT);
            for j=1:numel(attrNames)
                T=readtable(src{j},'VariableNamingRule','preserve');
                kname=local_za_var(T,{'gauge_id','gaugeid','station_id'});
                kid=regexprep(strip(string(T.(kname))),'\.0+$','');
                Tout=T(ismember(kid,sourceIds),:);
                writetable(Tout,fullfile(dirD,attrNames{j}));
            end
            copied=0;
            for j=1:numel(allCsv)
                [~,base]=fileparts(allCsv(j).name);
                if ismember(string(base),sourceIds)
                    copyfile(fullfile(allCsv(j).folder,allCsv(j).name), ...
                        fullfile(dirT,[base '.csv']),'f');
                    copied=copied+1;
                end
            end
            if copied~=spec.count
                error('data_helpers:GRDCTimeSeriesCount', ...
                    'Expected %d %s time-series files, copied %d.', ...
                    spec.count,spec.name,copied);
            end
            basinName=sprintf('%s_%d_basins.txt',spec.short,spec.count);
            writelines(ids,fullfile(dirD,basinName));
            sage=''; if isfield(cfg,'SAGEhydro'), sage=char(string(cfg.SAGEhydro)); end
            if isempty(sage)&&isfield(cfg,'SAGEdir'), sage=char(string(cfg.SAGEdir)); end
            if ~isempty(sage)
                rd=fullfile(sage,'regions',spec.short); local_mkdir(rd);
                writelines(ids,fullfile(rd,basinName));
            end
            logFcn(sprintf('Installed %d %s GRDC-Caravan basins in: %s', ...
                copied,spec.name,dirT));
        end

        clear T Tout O allCsv src
        local_za_remove_tree(tmp,logFcn,true);
        if deleteArchive, local_za_delete_file(zipFile,logFcn); end
        dirD=local_cfg_dirD(cfg,req.region); dirT=fullfile(dirD,'daily','timeseries');
        out=struct('ok',true,'region',req.region,'dirD',dirD,'dirM',dirT,'dirQ',dirT);
    catch ME
        logFcn(['GRDC-Caravan region install failed: ' ME.message]);
        out=struct('ok',false,'error',ME.message);
    end
    local_close_progress(d);
    function progress(info)
        [canceled,lastUI]=local_update_progress_dialog(d,info,lastUI,canceled);
    end
end

function [codes,canceled]=local_select_grdc_regions(ui,specs,installed,requested)
%LOCAL_SELECT_GRDC_REGIONS Checkbox selector for one shared archive.
    codes={}; canceled=true;
    hasUI=isstruct(ui)&&isfield(ui,'fig')&&~isempty(ui.fig)&&isvalid(ui.fig);
    if ~hasUI
        j=find(strcmp(cellfun(@(s)s.short,specs,'UniformOutput',false),requested),1);
        if ~isempty(j)&&~installed(j), codes={requested}; end
        canceled=false; return
    end
    dlg=uifigure('Name','Install GRDC-Caravan regions', ...
        'WindowStyle','modal','Position',[100 100 500 390], ...
        'CloseRequestFcn',@cancelSelection);
    g=uigridlayout(dlg,[numel(specs)+3 2]);
    g.RowHeight=[38 repmat(25,1,numel(specs)) 18 34];
    g.ColumnWidth={'1x',130};
    g.RowSpacing=4; g.ColumnSpacing=10;
    g.Padding=[16 12 16 12];
    lab=uilabel(g,'Text',['Select the country subsets to install. ' ...
        'The 8.8-GB archive is extracted only once.'],'WordWrap','on');
    lab.Layout.Row=1; lab.Layout.Column=[1 2];
    cb=gobjects(numel(specs),1);
    for k=1:numel(specs)
        s=specs{k}; suffix='';
        if installed(k), suffix='  [installed]'; end
        cb(k)=uicheckbox(g,'Text',sprintf('%s (%d basins)%s', ...
            s.name,s.count,suffix),'Value',strcmp(s.short,requested)||installed(k));
        cb(k).Layout.Row=k+1; cb(k).Layout.Column=[1 2];
        if installed(k), cb(k).Enable='off'; end
    end
    note=uilabel(g,'Text','Installed regions are retained and are not copied again.', ...
        'FontAngle','italic'); note.Layout.Row=numel(specs)+2; note.Layout.Column=[1 2];
    bCancel=uibutton(g,'Text','Cancel','ButtonPushedFcn',@cancelSelection);
    bCancel.Layout.Row=numel(specs)+3; bCancel.Layout.Column=1;
    bInstall=uibutton(g,'Text','Install selected','ButtonPushedFcn',@acceptSelection);
    bInstall.Layout.Row=numel(specs)+3; bInstall.Layout.Column=2;
    movegui(dlg,'center'); uiwait(dlg);
    function acceptSelection(~,~)
        codes = {};
        for q = 1:numel(specs)
            if ~installed(q) ...
                    && cb(q).Value
                codes{end+1} = specs{q}.short; %#ok
            end
        end
        canceled = false; uiresume(dlg); delete(dlg);
    end
    function cancelSelection(~,~)
        canceled = true; codes = {}; uiresume(dlg); delete(dlg);
    end
end

function spec = local_grdc_africa_spec(code)
    switch upper(char(string(code)))
        case 'ZA'
            spec=struct('short','ZA','region','CAMELS_ZA', ...
                'name','South Africa','count',434, ...
                'aliases',{{'ZA','ZAF','SOUTH AFRICA'}});
        case 'NA'
            spec=struct('short','NA','region','CAMELS_NA', ...
                'name','Namibia','count',51, ...
                'aliases',{{'NA','NAM','NAMIBIA'}});
        case 'AR'
            spec=struct('short','AR','region','CAMELS_AR', ...
                'name','Argentina','count',59,'aliases', ...
                {{'AR','ARG','ARGENTINA'}});
        case 'BE'
            spec=struct('short','BE','region','CAMELS_BE', ...
                'name','Belgium','count',55,'aliases', ...
                {{'BE','BEL','BELGIUM'}});
        case 'EE'
            spec=struct('short','EE','region','CAMELS_EE', ...
                'name','Estonia','count',51,'aliases', ...
                {{'EE','EST','ESTONIA'}});
        case 'IE'
            spec=struct('short','IE','region','CAMELS_IE', ...
                'name','Ireland','count',43,'aliases', ...
                {{'IE','IRL','IRELAND'}});
        case 'JM'
            spec=struct('short','JM','region','CAMELS_JM', ...
                'name','Jamaica','count',12,'aliases', ...
                {{'JM','JAM','JAMAICA'}});
        case 'NO'
            spec=struct('short','NO','region','CAMELS_NO', ...
                'name','Norway','count',206,'aliases', ...
                {{'NO','NOR','NORWAY'}});
        case 'PR'
            spec=struct('short','PR','region','CAMELS_PR', ...
                'name','Puerto Rico','count',24, ...
                'aliases',{{'PR','PRI','PUERTO RICO'}});
        otherwise
            error('data_helpers:BadGRDCAfricaCode', ...
                'Unknown code: %s',code);
    end
end

function tf = local_grdc_africa_installed(cfg,spec)
    c = cfg; c.region = spec.region; c.dirD = '';
    tf = local_has_stream(c,'daily') ...
        && local_has_metadata(c);
end

function f = local_za_find_one(root,name)
    q = dir(fullfile(root,'**',name));
    if isempty(q)
        error('data_helpers:ZAMissingArchiveFile', ...
            'Archive lacks %s.',name); 
    end
    f = fullfile(q(1).folder,q(1).name);
end

function n = local_za_var(T,candidates)
    vn = T.Properties.VariableNames;
    for j = 1:numel(candidates)
        q = find(strcmpi(vn,candidates{j}),1); 
        if ~isempty(q)
            n = vn{q}; 
            return; 
        end
    end
    error('data_helpers:ZAMissingColumn', ...
        'Missing required column (%s).', ...
        strjoin(candidates,', '));
end

function local_za_remove_tree(folder,logFcn,nonfatal)
    if nargin < 3
        nonfatal = false; 
    end
    if ~isfolder(folder)
        return; 
    end
    last = '';
    for attempt = 1:8
        try
            [ok,msg] = rmdir(folder,'s');
            if ok ...
                    || ~isfolder(folder)
                return; 
            end
            last = msg;
        catch ME
            last = ME.message;
        end
        pause(0.5*attempt); drawnow;
    end
    if nonfatal
        logFcn(['South Africa data are installed, but the temporary ' ...
            'extraction directory is still locked and can be removed ' ...
            'later: ' folder '. ' last]);
    else
        error('data_helpers:ZACleanupFailed','%s',last);
    end
end

function local_za_delete_file(file,logFcn)
    if ~isfile(file)
        return; 
    end
    last = '';
    for attempt = 1:8
        try
            delete(file);
            if ~isfile(file)
                return; 
            end
        catch ME
            last=ME.message;
        end
        pause(0.5*attempt); drawnow;
    end
    logFcn(['South Africa data are installed, ' ...
        'but the downloaded archive ' ...
        ['is still locked and should be ' ...
        'deleted manually: '] file '. ' last]);
end

% =============================================
function out = local_download_nz(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_NZ Download and install CAMELS-NZ daily or hourly data.

    logFcn = local_ui_log(ui);

    if ~any(strcmp(stream, ...
            {'daily','hourly'}))
        error('data_helpers:NZBadStream', ...
            ['CAMELS-NZ supports daily and hourly data. ' ...
            'Use stream = ''daily'' or ''hourly''.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_NZ');
    dirT = fullfile(dirD,stream,'timeseries');

    local_mkdir(dirD);
    local_mkdir(dirT);
    local_mkdir(fullfile(dirT,'precipitation'));
    local_mkdir(fullfile(dirT,'temperature'));
    local_mkdir(fullfile(dirT,'streamflow'));

    if strcmp(stream,'hourly')
        local_mkdir(fullfile(dirT,'pet'));
        dataText = ['hourly precipitation, temperature, ' ...
            'PET, and streamflow'];
    else
        local_mkdir(fullfile(dirT,'relative_humidity'));
        dataText = ['daily precipitation, temperature, ' ...
            'relative humidity, and streamflow'];
    end

    if ~local_ui_confirm(ui, ...
            'Download New Zealand data', ...
            sprintf(['CAMELS-NZ %s requires five Figshare ZIP files: ' ...
            'the shared catchment attributes plus %s.' newline ...
            'Continue?'],stream,dataText))
        out = struct('ok',false,'canceled',true);
        return
    end

    files = local_camels_nz_files(stream);
    downloadDir = local_default_download_dir();
    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading New Zealand data', ...
        sprintf('Starting CAMELS-NZ %s download ...',stream));

    try
        if local_nz_install_complete(dirD,stream)
            logFcn(sprintf(['CAMELS-NZ %s appears complete; ' ...
                'skipping download/install.'],stream));
            out = struct('ok',true, ...
                'region','CAMELS_NZ', ...
                'dirD',dirD, ...
                'dirM',dirT, ...
                'dirQ',fullfile(dirT,'streamflow'), ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        for i = 1:numel(files)
            if userCanceled
                out = struct('ok',false,'canceled',true);
                local_close_progress(d);
                return
            end

            f = files(i);
            targetFile = fullfile(downloadDir,f.fileName);
            files(i).localFile = targetFile;
            files(i).install = true;

            if isfile(targetFile) ...
                    && ~isempty(f.md5)
                try
                    local_verify_md5(targetFile, ...
                        f.md5,logFcn);
                    logFcn(['Using verified existing file: ' ...
                        targetFile]);
                catch
                    logFcn(['Existing CAMELS-NZ ZIP failed MD5; ' ...
                        'deleting and downloading again: ' ...
                        targetFile]);
                    delete(targetFile);
                end
            elseif isfile(targetFile)
                logFcn(['Using existing file: ' ...
                    targetFile]);
            end

            if ~isfile(targetFile)
                logFcn(sprintf('Downloading %s to: %s', ...
                    f.label,targetFile));
                local_download_file_retry_nz(f.url,targetFile, ...
                    f.label,d,@() userCanceled, ...
                    @(info) ui_progress(info,i,numel(files)), ...
                    logFcn,3);
            end

            % Chrome may save Figshare files with .zip appended.
            altFile = [targetFile '.zip'];
            if ~isfile(targetFile) ...
                    && isfile(altFile)
                movefile(altFile,targetFile,'f');
            end

            files(i).localFile = targetFile;
            if ~isempty(f.md5)
                local_verify_md5(targetFile, ...
                    f.md5,logFcn);
            end
        end

        if userCanceled
            out = struct('ok',false,'canceled',true);
            local_close_progress(d);
            return
        end

        unzipDir = fullfile(downloadDir, ...
            ['CAMELS_NZ_' stream '_download']);
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        for i = 1:numel(files)
            f = files(i);
            dst = fullfile(unzipDir, ...
                erase(f.fileName,'.zip'));
            local_mkdir(dst);
            logFcn(sprintf('Unzipping %s ...', ...
                f.fileName));
            if ~isempty(d) ...
                    && isvalid(d)
                try
                    d.Indeterminate = true;
                    d.Message = ['Unzipping ' ...
                        f.fileName ' ...'];
                    drawnow limitrate nocallbacks
                catch
                end
            end
            unzip(files(i).localFile,dst);
            files(i).unzipDir = dst;
        end

        local_install_nz_files( ...
            files,dirD,stream,logFcn);

        okAttr = local_all_files_exist(dirD, ...
            local_nz_metadata_files());
        okStream = local_has_nz_timeseries( ...
            dirD,stream);

        out = struct('ok',okAttr && okStream, ...
            'region','CAMELS_NZ', ...
            'dirD',dirD, ...
            'dirM',dirT, ...
            'dirQ',fullfile(dirT,'streamflow'), ...
            'archive','', ...
            'unzipDir','');

        if out.ok
            if strcmp(stream,'hourly')
                pP = 'precipitation_station_id_*.csv';
                pT = 'temperature_station_id_*.csv';
                pX = 'PET_station_id_*.csv';
                pQ = 'flow_station_id_*.csv';
                xFolder = 'pet';
                xText = 'PET';
            else
                pP = 'daily_precipitation_station_id_*.csv';
                pT = 'daily_temperature_station_id_*.csv';
                pX = 'daily_RH_station_id_*.csv';
                pQ = 'daily_flow_station_id_*.csv';
                xFolder = 'relative_humidity';
                xText = 'relative humidity';
            end

            nP = local_count_files(fullfile(dirT, ...
                'precipitation'),{pP});
            nT = local_count_files(fullfile(dirT, ...
                'temperature'),{pT});
            nX = local_count_files(fullfile(dirT, ...
                xFolder),{pX});
            nQ = local_count_files(fullfile(dirT, ...
                'streamflow'),{pQ});

            logFcn(sprintf(['CAMELS-NZ %s install finished: ' ...
                '%d precipitation, %d temperature, ' ...
                '%d %s, and %d streamflow CSV files.'], ...
                stream,nP,nT,nX,xText,nQ));

            local_cleanup_download_files(files,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
        else
            logFcn(sprintf(['CAMELS-NZ %s install finished, but ' ...
                'some required files appear to be missing.'],stream));
        end

    catch ME
        try
            if exist('files','var')
                local_cleanup_download_files(files,logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder( ...
                    unzipDir,logFcn);
            end
        catch
        end
        logFcn(sprintf('CAMELS-NZ %s install failed:',stream));
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info,iFile,nFile)
        if nargin >= 2
            if isfield(info,'frac') ...
                    && isfinite(info.frac)
                info.frac = ((iFile - 1) ...
                    + info.frac) / nFile;
            else
                info.frac = (iFile - 1) / nFile;
            end
            if isfield(info,'label')
                info.label = sprintf('%s (%d/%d)', ...
                    info.label,iFile,nFile);
            else
                info.label = sprintf('File %d/%d', ...
                    iFile,nFile);
            end
        end
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end


% =============================================
function out = local_download_se(cfg,stream,ui)
% =============================================
%LOCAL_DOWNLOAD_SE Download and install CAMELS-SE files.

    if nargin < 2 ...
            || isempty(stream)
        stream = 'daily';
    end
    if nargin < 3
        ui = [];
    end
    logFcn = local_ui_log(ui);

    d = local_progress_dialog(ui, ...
        'Downloading Sweden data', ...
        'Starting CAMELS-SE download ...');
    try
        if ~isempty(d) && isvalid(d)
            d.Indeterminate = true;
        end
    catch
    end

    if isfield(cfg,'dirD') ...
            && ~isempty(cfg.dirD)
        dirSE = cfg.dirD;
    elseif isfield(cfg,'root') ...
            && ~isempty(cfg.root)
        dirSE = fullfile(cfg.root, ...
            'Data','CAMELS_SE');
    else
        dirSE = fullfile(pwd, ...
            'Data','CAMELS_SE');
    end
    
    if ~isfolder(dirSE)
        mkdir(dirSE);
    end
    
    dirDaily = fullfile(dirSE,'daily');
    if ~isfolder(dirDaily)
        mkdir(dirDaily);
    end
    
    tmpDir = fullfile(dirSE, ...
        '_download_tmp');
    if ~isfolder(tmpDir)
        mkdir(tmpDir);
    end
    
    urlProp = ['https://api.researchdata.se/' ...
        'dataset/2023-173/1/file/data?' ...
        'filePath=catchment+properties.zip'];
    urlTS = ['https://api.researchdata.se/' ...
        'dataset/2023-173/1/file/data?' ...
        'filePath=catchment+time+series.zip'];
    
    try
        if ~isempty(d) ...
                && isvalid(d)
            d.Message = ['Downloading CAMELS-SE ' ...
                'catchment properties ...'];
            drawnow limitrate
        end
    catch
    end
    logFcn(['Downloading CAMELS-SE ' ...
        'catchment properties ...']);
    zipProp = fullfile(tmpDir, ...
        'catchment_properties.zip');
    local_websave(zipProp,urlProp);

    try
        if ~isempty(d) ...
                && isvalid(d)
            d.Message = ['Unzipping CAMELS-SE ' ...
                'catchment properties ...'];
            drawnow limitrate
        end
    catch
    end
    logFcn(['Unzipping CAMELS-SE ' ...
        'catchment properties ...']);
    unzipDirProp = fullfile(tmpDir, ...
        'catchment_properties_unzip');
    local_reset_dir(unzipDirProp);
    unzip(zipProp,unzipDirProp);
    
    % The archive creates a folder named "catchment properties".
    D = dir(fullfile(unzipDirProp, ...
        '**','*.csv'));
    for k = 1:numel(D)
        src = fullfile(D(k).folder,D(k).name);
        dst = fullfile(dirSE,D(k).name);
        copyfile(src,dst,'f');
    end

    try
        if ~isempty(d) ...
                && isvalid(d)
            d.Message = ['Downloading CAMELS-SE ' ...
                'daily catchment time series ...'];
            drawnow limitrate
        end
    catch
    end
    logFcn(['Downloading CAMELS-SE ' ...
        'daily catchment time series ...']);
    zipTS = fullfile(tmpDir, ...
        'catchment_time_series.zip');
    local_websave(zipTS,urlTS);
    
    try
        if ~isempty(d) ...
                && isvalid(d)
            d.Message = ['Unzipping CAMELS-SE ' ...
                'daily catchment time series ...'];
            drawnow limitrate
        end
    catch
    end
    logFcn(['Unzipping CAMELS-SE ' ...
        'daily catchment time series ...']);
    unzipDirTS = fullfile(tmpDir, ...
        'catchment_time_series_unzip');
    local_reset_dir(unzipDirTS);
    unzip(zipTS,unzipDirTS);
    
    % The archive creates a folder named "catchment time series".
    D = dir(fullfile(unzipDirTS, ...
        '**','catchment_id_*.csv'));
    for k = 1:numel(D)
        src = fullfile(D(k).folder,D(k).name);
        dst = fullfile(dirDaily,D(k).name);
        copyfile(src,dst,'f');
    end
    
    out = struct();
    out.region = 'CAMELS_SE';
    out.dirD = dirSE;
    out.dirM = dirDaily;
    out.dirQ = dirDaily;
    out.n_daily_files = numel(dir(fullfile(dirDaily, ...
        'catchment_id_*.csv')));
    out.ok_attributes = isfile(fullfile(dirSE, ...
        'catchments_physical_properties.csv')) ...
        && isfile(fullfile(dirSE, ...
        'catchments_landcover.csv')) ...
        && isfile(fullfile(dirSE, ...
        'catchments_soil_classes.csv'));
    out.ok_daily = out.n_daily_files > 0;
    
    out.ok = out.ok_attributes ...
        && out.ok_daily;

    if out.ok_attributes ...
            && out.ok_daily
        logFcn(sprintf(['CAMELS-SE download ' ...
            'complete: %d daily ' ...
            'catchment files installed.'], ...
            out.n_daily_files));
    else
        logFcn(['CAMELS-SE download ' ...
            'finished, but some ' ...
            'files are missing.']);
    end
    
    if strcmpi(stream,'hourly')
        warning('data_helpers:CAMELS_SE_dailyOnly', ...
            ['CAMELS-SE support ' ...
            'is daily only; ' ...
            'hourly stream was ignored.']);
    end

    % Final clean-up block
    try
        local_close_progress(d);
        d = [];
    
        drawnow;
        fclose('all');
        pause(0.5);
    
        if isfolder(tmpDir)
            local_remove_dir_retry(tmpDir,logFcn);
        end
    catch ME
        logFcn(['Could not remove CAMELS-SE ' ...
            'temporary download folder: ' ...
            ME.message]);
    end
    local_close_progress(d);
end

% =============================================
function out = local_download_is(cfg,stream,ui)
%LOCAL_DOWNLOAD_IS Download/install LamaH-Ice daily or hourly data.
%
% Download order:
%   1. Reuse an existing requested inner archive in _download_tmp_IS.
%   2. Reuse an existing complete HydroShare bag.
%   3. Try the individual inner-archive HydroShare endpoints first.
%   4. Fall back to the complete HydroShare resource bag.
%
% The complete HydroShare bag contains:
%   <resource-id>/<resource-id>/data/contents/lamah_ice.zip
%   <resource-id>/<resource-id>/data/contents/lamah_ice_hourly.zip
%
% Archive selection:
%   daily  -> lamah_ice.zip
%   hourly -> lamah_ice_hourly.zip
%
% The hourly inner archive also contains daily forcing and discharge data.

    logFcn = local_ui_log(ui);
    dirD = local_cfg_dirD(cfg,'CAMELS_IS');
    tmpRoot = fullfile(local_default_download_dir(), ...
        'CAMELS_IS_download');

    if ~any(strcmp(stream,{'daily','hourly'}))
        error('data_helpers:ISBadStream', ...
            ['LamaH-Ice supports stream = ' ...
             '''daily'' or ''hourly''.']);
    end

    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,stream,'forcing'));
    local_mkdir(fullfile(dirD,stream,'discharge'));

    if local_has_stream(cfg,stream) ...
            && local_has_metadata(cfg)
        out = struct('ok',true, ...
            'region','CAMELS_IS', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,stream,'forcing'), ...
            'dirQ',fullfile(dirD,stream,'discharge'));
        logFcn(sprintf(['LamaH-Ice %s data are already ' ...
            'installed; skipping download.'],stream));
        if isfolder(tmpRoot)
            logFcn(['Removing residual LamaH-Ice download files from ' ...
                'an earlier installation.']);
            local_cleanup_download_folder(tmpRoot,logFcn);
        end
        return
    end

    resourceID = '705d69c0f77c48538d83cf383f8c63d6';

    local_mkdir(tmpRoot);

    dailyZip = fullfile(tmpRoot,'lamah_ice.zip');
    hourlyZip = fullfile(tmpRoot,'lamah_ice_hourly.zip');
    bagName = [resourceID '.zip'];
    bagFile = fullfile(tmpRoot,bagName);

    if strcmp(stream,'hourly')
        archiveName = 'lamah_ice_hourly.zip';
        zipFile = hourlyZip;
        %archiveContainsHourly = true;
    else
        % A daily request always uses the dedicated daily archive. This
        % prevents the much larger hourly package from being selected just
        % because it happens to be present from an earlier attempt.
        archiveName = 'lamah_ice.zip';
        zipFile = dailyZip;
        %archiveContainsHourly = false;
    end

    if strcmp(stream,'hourly')
        msg = [ ...
            'LamaH-Ice hourly installation downloads the 9.1 GB ' ...
            'lamah_ice_hourly.zip archive from HydroShare. It requires ' ...
            'approximately 30 GB while extracted and contains both hourly ' ...
            'and daily forcing and discharge data.' ...
            newline 'Continue?'];
    else
        msg = [ ...
            'LamaH-Ice daily installation uses lamah_ice.zip from the ' ...
            'complete HydroShare resource bag, or downloads that archive ' ...
            'directly. Only daily data are installed.' ...
            newline 'Continue?'];
    end

    if ~local_ui_confirm(ui,'Download Iceland data',msg)
        out = struct('ok',false,'canceled',true);
        return
    end

    d = local_progress_dialog(ui, ...
        'Downloading Iceland data', ...
        ['Locating ' archiveName ' ...']);

    okArchive = local_is_valid_archive(zipFile);
    lastErr = '';

    if okArchive
        logFcn(['Using existing LamaH-Ice archive: ' zipFile]);
    else
        % -------------------------------------------------------------
        % First request only the selected inner archive. In particular,
        % this prevents an hourly install from appearing to stall while
        % the much larger complete HydroShare resource bag is downloaded.
        % -------------------------------------------------------------
        innerURLs = { ...
            ['https://www.hydroshare.org/resource/' resourceID ...
             '/data/contents/' archiveName], ...
            ['https://www.hydroshare.org/hsapi/resource/' resourceID ...
             '/files/' archiveName]};

        for i = 1:numel(innerURLs)
            try
                logFcn(['Trying individual LamaH-Ice archive URL: ' ...
                    innerURLs{i}]);

                if ~isempty(d) && isvalid(d)
                    d.Message = ['Downloading 9.1 GB ' archiveName ...
                        '. This may take a long time ...'];
                    drawnow limitrate nocallbacks
                end

                if strcmp(stream,'hourly')
                    local_download_is_large_archive( ...
                        innerURLs{i},zipFile,archiveName,logFcn);
                else
                    local_download_file_retry( ...
                        innerURLs{i},zipFile,archiveName,d, ...
                        @() false,[],logFcn,3);
                end

                okArchive = local_is_valid_archive(zipFile);
                if okArchive
                    break
                end
            catch ME
                lastErr = ME.message;
            end
        end

        % -------------------------------------------------------------
        % If direct content access is unavailable, download the complete
        % HydroShare resource bag and extract the requested nested archive.
        % -------------------------------------------------------------
        bagURL = ['https://www.hydroshare.org/django_irods/' ...
            'download/bags/' resourceID '.zip'];

        okBag = false;
        if ~okArchive
            okBag = local_is_valid_archive(bagFile);
        end

        if ~okArchive && okBag
            logFcn(['Using existing HydroShare resource bag: ' bagFile]);
        elseif ~okArchive
            if isfile(bagFile)
                try
                    delete(bagFile);
                catch
                end
            end

            try
                logFcn(['Trying direct HydroShare resource-bag download: ' ...
                    bagURL]);

                local_download_file_retry( ...
                    bagURL,bagFile, ...
                    'Complete LamaH-Ice HydroShare resource', ...
                    d,@() false,[],logFcn,3);

                okBag = local_is_valid_archive(bagFile);
            catch ME
                lastErr = ME.message;
                okBag = false;
                logFcn(['Direct HydroShare bag download failed: ' ...
                    ME.message]);
            end
        end

        if ~okArchive && okBag
            try
                zipFile = local_extract_is_inner_archive( ...
                    bagFile,tmpRoot,archiveName,d,logFcn);
                okArchive = local_is_valid_archive(zipFile);
            catch ME
                lastErr = ME.message;
                okArchive = false;
                logFcn(['Could not extract ' archiveName ...
                    ' from the HydroShare bag: ' ME.message]);
            end
        end

    end

    if ~okArchive
        local_close_progress(d);
        error('data_helpers:ISDownloadFailed', ...
            ['Could not obtain %s automatically. Download either the ' ...
             'complete HydroShare bag (%s) or %s manually and place it ' ...
             'in %s. Last error: %s'], ...
             archiveName,bagName,archiveName,tmpRoot,lastErr);
    end

    unzipRoot = fullfile(tmpRoot, ...
        ['unzipped_' erase(archiveName,'.zip')]);
    local_reset_dir(unzipRoot);

    try
        if ~isempty(d) && isvalid(d)
            d.Message = ['Extracting ' archiveName ' ...'];
            drawnow limitrate nocallbacks
        end
    catch
    end

    try
        unzip(zipFile,unzipRoot);
    catch ME
        local_close_progress(d);
        error('data_helpers:ISUnzipFailed', ...
            ['Could not extract %s. The downloaded file may be an HTML ' ...
             'page rather than a ZIP archive. Cause: %s'], ...
             zipFile,ME.message);
    end

    if contains(lower(archiveName),'hourly')
        archiveContainsHourly = true;
        base = local_find_dir_recursive(unzipRoot,'lamah_ice_hourly');
    else
        archiveContainsHourly = false;
        base = local_find_dir_recursive(unzipRoot,'lamah_ice');
    end

    if isempty(base)
        base = unzipRoot;
    end

    % -----------------
    % Static attributes
    % -----------------
    fc = local_find_file_recursive_manual( ...
        base,'Catchment_attributes.csv');
    fg = local_find_file_recursive_manual( ...
        base,'Gauge_attributes.csv');

    if isempty(fc) || isempty(fg)
        local_close_progress(d);
        error('data_helpers:ISMissingAttributes', ...
            ['Could not find Catchment_attributes.csv and ' ...
             'Gauge_attributes.csv in %s.'],archiveName);
    end

    copyfile(fc,fullfile(dirD, ...
        'Catchment_attributes.csv'),'f');
    copyfile(fg,fullfile(dirD, ...
        'Gauge_attributes.csv'),'f');

    extras = { ...
        'hydro_indices_1981_2018.csv', ...
        'hydro_indices_1981_2018_unfiltered.csv', ...
        'meteorological_data_means_1981_to_2018.csv', ...
        'meteorological_data_means_1991_to_2018.csv', ...
        'water_balance.csv', ...
        'water_balance_unfiltered.csv'};

    for i = 1:numel(extras)
        fExtra = local_find_file_recursive_manual(base,extras{i});
        if ~isempty(fExtra)
            copyfile(fExtra, ...
                fullfile(dirD,extras{i}),'f');
        end
    end

    % ----------------
    % Daily time series
    % ----------------
    if strcmp(stream,'daily') || archiveContainsHourly
        srcDF = local_find_is_timeseries_dir( ...
            base,'daily','forcing');
        srcDQ = local_find_is_timeseries_dir( ...
            base,'daily','discharge');

        dstDF = fullfile( ...
            dirD,'daily','forcing');
        dstDQ = fullfile( ...
            dirD,'daily','discharge');

        local_assert_is_timeseries_count(srcDF,111,'daily forcing');
        local_assert_is_timeseries_count(srcDQ,111,'daily discharge');

        local_copy_pattern(srcDF, ...
            dstDF,{'ID_*.csv'},logFcn);
        local_copy_pattern(srcDQ, ...
            dstDQ,{'ID_*.csv'},logFcn);
    end

    % -----------------
    % Hourly time series
    % -----------------
    if strcmp(stream,'hourly')
        srcHF = local_find_is_timeseries_dir( ...
            base,'hourly','forcing');
        srcHQ = local_find_is_timeseries_dir( ...
            base,'hourly','discharge');

        dstHF = fullfile( ...
            dirD,'hourly','forcing');
        dstHQ = fullfile( ...
            dirD,'hourly','discharge');

        % The hourly meteorological branch contains forcing for all 111
        % LamaH-Ice basins, whereas hourly discharge is available for 76
        % gauges. Install only the forcing files whose gauge IDs occur in
        % the hourly discharge branch.
        local_assert_is_timeseries_min_count(srcHF,76,'hourly forcing');
        local_assert_is_timeseries_count(srcHQ,76,'hourly discharge');

        local_copy_is_hourly_pair(srcHF,srcHQ, ...
            dstHF,dstHQ,logFcn);
    end

    okRequested = local_has_stream(cfg,stream);
    okMetadata = local_has_metadata(cfg);

    out = struct('ok',okRequested && okMetadata, ...
        'region','CAMELS_IS', ...
        'dirD',dirD, ...
        'dirM',fullfile(dirD,stream,'forcing'), ...
        'dirQ',fullfile(dirD,stream,'discharge'), ...
        'archive',zipFile, ...
        'bag',bagFile, ...
        'unzipDir',unzipRoot);

    if out.ok
        nF = local_count_files(fullfile(dirD, ...
            stream,'forcing'),{'ID_*.csv'});
        nQ = local_count_files(fullfile(dirD, ...
            stream,'discharge'),{'ID_*.csv'});

        logFcn(sprintf(['LamaH-Ice %s installation complete: ' ...
            '%d forcing files and %d discharge files.'], ...
            stream,nF,nQ));

        if strcmp(stream,'hourly') ...
                && local_has_stream(cfg,'daily')
            logFcn(['Daily LamaH-Ice data were also installed from ' ...
                'lamah_ice_hourly.zip.']);
        end

        % Close the progress dialog before removing the extraction tree
        % and downloaded archives. On Windows, open UI/file handles can
        % otherwise prevent removal of the temporary directory.
        local_close_progress(d);
        d = [];
        drawnow limitrate nocallbacks

        try
            local_cleanup_download_folder(tmpRoot,logFcn);
            out.archive = '';
            out.bag = '';
            out.unzipDir = '';
        catch ME
            logFcn(['LamaH-Ice installation succeeded, but the temporary ' ...
                'download directory could not be removed: ' ME.message]);
        end
    else
        logFcn(['LamaH-Ice extraction finished, but required ' ...
            stream ' files or metadata are missing.']);
    end

    local_close_progress(d);
end

function local_assert_is_timeseries_count(folder,expected,label)
%LOCAL_ASSERT_IS_TIMESERIES_COUNT Validate LamaH-Ice before installation.

    actual = local_count_files(folder,{'ID_*.csv'});
    if actual ~= expected
        error('data_helpers:ISIncompleteArchive', ...
            ['Expected %d LamaH-Ice %s files in the extracted archive, ' ...
            'but found %d. Existing installed files were not changed.'], ...
            expected,label,actual);
    end
end

function local_assert_is_timeseries_min_count(folder,expected,label)
%LOCAL_ASSERT_IS_TIMESERIES_MIN_COUNT Require at least expected files.

    actual = local_count_files(folder,{'ID_*.csv'});
    if actual < expected
        error('data_helpers:ISIncompleteArchive', ...
            ['Expected at least %d LamaH-Ice %s files in the extracted ' ...
             'archive, but found %d. Existing installed files were not ' ...
             'changed.'],expected,label,actual);
    end
end

function local_copy_is_hourly_pair(srcF,srcQ,dstF,dstQ,logFcn)
%LOCAL_COPY_IS_HOURLY_PAIR Install the common 76 hourly gauge IDs.

    Q = dir(fullfile(srcQ,'ID_*.csv'));
    names = string({Q(~[Q.isdir]).name});
    names = sort(names(:));

    missing = strings(0,1);
    for i = 1:numel(names)
        if ~isfile(fullfile(srcF,names(i)))
            missing(end+1,1) = names(i); %#ok<AGROW>
        end
    end

    if ~isempty(missing)
        error('data_helpers:ISHourlyGaugeMismatch', ...
            ['Hourly meteorological data are missing %d gauge file(s) ' ...
             'that occur in hourly discharge, for example %s. Existing ' ...
             'installed files were not changed.'], ...
            numel(missing),missing(1));
    end

    local_mkdir(dstF);
    local_mkdir(dstQ);

    for i = 1:numel(names)
        copyfile(fullfile(srcF,names(i)), ...
            fullfile(dstF,names(i)),'f');
        copyfile(fullfile(srcQ,names(i)), ...
            fullfile(dstQ,names(i)),'f');
    end

    % Remove stale ID files left by an interrupted or older installer only
    % after every requested source pair has been validated and copied.
    local_remove_is_unmatched_hourly(dstF,names);
    local_remove_is_unmatched_hourly(dstQ,names);

    logFcn(sprintf(['Copied %d matching hourly forcing and discharge ' ...
        'files to %s and %s'],numel(names),dstF,dstQ));
end

function local_remove_is_unmatched_hourly(folder,names)
%LOCAL_REMOVE_IS_UNMATCHED_HOURLY Remove stale hourly ID files only.

    D = dir(fullfile(folder,'ID_*.csv'));
    for i = 1:numel(D)
        if D(i).isdir || any(strcmp(string(D(i).name),names))
            continue
        end
        delete(fullfile(D(i).folder,D(i).name));
    end
end

% =====================================================
function innerZip = local_extract_is_inner_archive( ...
        bagFile,tmpRoot,archiveName,d,logFcn)
% =====================================================
% Extract the HydroShare resource bag and copy the requested nested archive
% from <resource-id>/<resource-id>/data/contents into _download_tmp_IS.

    bagRoot = fullfile(tmpRoot,'unzipped_hydroshare_bag');
    local_reset_dir(bagRoot);

    try
        if ~isempty(d) && isvalid(d)
            d.Message = ['Extracting HydroShare resource bag to locate ' ...
                archiveName ' ...'];
            drawnow limitrate nocallbacks
        end
    catch
    end

    unzip(bagFile,bagRoot);

    found = local_find_file_recursive_manual(bagRoot,archiveName);
    if isempty(found)
        error('data_helpers:ISMissingInnerArchive', ...
            'Could not find %s inside HydroShare resource bag %s.', ...
            archiveName,bagFile);
    end

    innerZip = fullfile(tmpRoot,archiveName);
    copyfile(found,innerZip,'f');

    logFcn(['Extracted nested LamaH-Ice archive from HydroShare bag: ' ...
        innerZip]);

    try
        local_cleanup_download_folder(bagRoot,logFcn);
    catch
    end
end

% =========================================
function tf = local_is_valid_archive(fname)
% =========================================
    tf = false;
    if ~isfile(fname)
        return
    end

    info = dir(fname);
    largeEnough = ~isempty(info) ...
        && isfinite(info.bytes) ...
        && info.bytes > 1e6;
    if ~largeEnough
        return
    end

    % A partial multi-gigabyte transfer is nonempty but lacks a valid ZIP
    % central directory. Opening ZipFile validates that directory without
    % extracting the archive or reading its full contents into memory.
    try
        z = java.util.zip.ZipFile(java.io.File(fname));
        z.close();
        tf = true;
    catch
        tf = false;
    end
end

function local_download_is_large_archive(url,targetFile,label,logFcn)
%LOCAL_DOWNLOAD_IS_LARGE_ARCHIVE Resume the 9.1 GB LamaH-Ice transfer.

    local_mkdir(fileparts(targetFile));
    if local_is_valid_archive(targetFile)
        logFcn(['Using complete existing archive: ' targetFile]);
        return
    end

    if isfile(targetFile) ...
            && ~local_has_zip_signature(targetFile)
        logFcn(['Discarding non-ZIP response: ' targetFile]);
        delete(targetFile);
    end

    existingBytes = 0;
    if isfile(targetFile)
        info = dir(targetFile);
        existingBytes = info.bytes;
    end
    if existingBytes > 0
        logFcn(sprintf('Resuming partial %s download at %.2f GB: %s', ...
            label,existingBytes/1e9,targetFile));
    else
        logFcn(strcat('Downloading ',label,' to ',targetFile));
    end

    cmd = sprintf(['curl -L --fail --ssl-no-revoke --retry 5 ' ...
        '--retry-delay 5 --continue-at - -o "%s" "%s"'], ...
        targetFile,url);
    [status,message] = system(cmd);
    if status ~= 0
        error('data_helpers:ISLargeDownloadFailed', ...
            ['The resumable LamaH-Ice download stopped. Its partial file ' ...
             'was retained and will resume on the next attempt. curl: %s'], ...
            strtrim(message));
    end
    if ~local_is_valid_archive(targetFile)
        error('data_helpers:ISIncompleteArchive', ...
            ['The LamaH-Ice transfer ended without a valid ZIP central ' ...
             'directory. The partial file was retained for resumption: %s'], ...
            targetFile);
    end
end

function tf = local_has_zip_signature(fname)
%LOCAL_HAS_ZIP_SIGNATURE Distinguish a partial ZIP from an HTML response.

    tf = false;
    fid = fopen(fname,'r');
    if fid < 0
        return
    end
    cleanup = onCleanup(@() fclose(fid));
    signature = fread(fid,2,'*uint8');
    tf = numel(signature) == 2 ...
        && signature(1) == uint8('P') ...
        && signature(2) == uint8('K');
    clear cleanup
end

% =============================================================
function d = local_find_is_timeseries_dir( ...
        root,resolution,kind)
%LOCAL_FIND_IS_TIMESERIES_DIR Locate exact LamaH-Ice source folder.
%
% Forcing:
%   A_basins_total_upstrm/2_timeseries/<resolution>
%
% Discharge:
%   D_gauges/2_timeseries/<resolution>
%
% Daily discharge may be stored in daily_filtered. Hourly discharge must
% come specifically from D_gauges/2_timeseries/hourly.

    resolution = lower( ...
        char(string(resolution)));

    kind = lower( ...
        char(string(kind)));

    if strcmp(kind,'forcing')
        branchName = ...
            'A_basins_total_upstrm';
    elseif strcmp(kind,'discharge')
        branchName = 'D_gauges';
    else
        error('data_helpers:ISBadKind', ...
            'Unknown Iceland time-series kind: %s.', ...
            kind);
    end

    branch = local_find_dir_recursive( ...
        root,branchName);

    if isempty(branch)
        error('data_helpers:ISMissingBranch', ...
            ['Could not locate %s in the ' ...
             'LamaH-Ice archive.'], ...
            branchName);
    end

    tsRoot = fullfile( ...
        branch,'2_timeseries');

    if ~isfolder(tsRoot)
        error('data_helpers:ISMissingTimeseriesRoot', ...
            ['Could not locate %s/2_timeseries ' ...
             'in the LamaH-Ice archive.'], ...
            branchName);
    end

    %candidates = {};

    if strcmp(resolution,'daily') ...
            && strcmp(kind,'discharge')
        candidates = { ...
            fullfile(tsRoot, ...
                'daily_filtered'), ...
            fullfile(tsRoot,'daily')};
    else
        candidates = { ...
            fullfile(tsRoot,resolution)};
    end

    for i = 1:numel(candidates)
        candidate = candidates{i};

        if local_is_timeseries_folder( ...
                candidate,kind)
            d = candidate;
            return
        end
    end

    % Tolerate one extra release-specific folder below 2_timeseries.
    D = dir(fullfile(tsRoot,'**','ID_1.csv'));

    for i = 1:numel(D)
        folder = D(i).folder;
        folderKey = lower( ...
            strrep(folder,'\','/'));

        if ~contains(folderKey, ...
                ['/' resolution])
            continue
        end

        if local_is_timeseries_folder( ...
                folder,kind)
            d = folder;
            return
        end
    end

    error('data_helpers:ISMissingFolder', ...
        ['Could not locate LamaH-Ice %s %s ' ...
         'directory below %s.'], ...
        resolution,kind,tsRoot);
end

% =====================================================
function tf = local_is_timeseries_folder(folder,kind)
% =====================================================

    tf = false;

    if ~isfolder(folder)
        return
    end

    sample = fullfile(folder,'ID_1.csv');

    if ~isfile(sample)
        D = dir(fullfile(folder,'ID_*.csv'));

        if isempty(D)
            return
        end

        sample = fullfile( ...
            D(1).folder,D(1).name);
    end

    tf = local_is_csv_kind(sample,kind);
end

% ============================================
function tf = local_is_csv_kind(fname,kind)
% ============================================
% Verify file type from its header so forcing files cannot satisfy the
% discharge installation check and vice versa.

    tf = false;

    if ~isfile(fname)
        return
    end

    fid = fopen(fname,'r');

    if fid < 0
        return
    end

    cleanup = onCleanup(@() fclose(fid));
    header = fgetl(fid);

    if ~ischar(header)
        return
    end

    key = lower(regexprep( ...
        string(header),'[^a-z0-9]',''));

    if strcmpi(kind,'forcing')
        hasPrecip = contains(key,'prec') ...
            || contains(key,'precipitation');

        hasTemp = contains(key,'temp');

        tf = hasPrecip && hasTemp;

    elseif strcmpi(kind,'discharge')
        tf = contains(key,'qobs') ...
            || contains(key,'discharge');
    end
end

% =============================================
function out = local_download_ca(cfg,stream,ui)
%LOCAL_DOWNLOAD_CA Download/install CAMELS-SPAT Canada through Globus.
%
% Daily:
%   Daymet lumped forcing + daily observed discharge.
% Hourly:
%   RDRS lumped forcing + hourly observed discharge.
%
% Hourly files are transferred sequentially (discharge first, then RDRS)
% so that the C:\Users\<user>\Downloads staging folder never needs to hold
% both complete products at the same time.

    logFcn = local_ui_log(ui);
    stream = lower(strtrim(char(stream)));
    if ~any(strcmp(stream,{'daily','hourly'}))
        error('data_helpers:CABadStream', ...
            'CAMELS-CA supports stream = daily or hourly.');
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_CA');
    if strcmp(stream,'daily')
        dirF = fullfile(dirD,'daily','forcing','daymet');
        dirQ = fullfile(dirD,'daily','discharge');
        forcingLabel = 'Daymet';
    else
        dirF = fullfile(dirD,'hourly','forcing','rdrs');
        dirQ = fullfile(dirD,'hourly','discharge');
        forcingLabel = 'RDRS';
    end
    local_mkdir(dirD);
    local_mkdir(dirF);
    local_mkdir(dirQ);

    if local_ca_install_complete(dirD,stream)
        logFcn(sprintf(['CAMELS-CA %s data appear complete; ' ...
            'skipping download.'],stream));
        out = local_ca_result(true,dirD,'',stream);
        return
    end

    if strcmp(stream,'hourly')
        msg = [ ...
            'SAGE will use Globus to download Canadian CAMELS-SPAT ' ...
            'hourly discharge and lumped RDRS forcing.' newline newline ...
            'To limit temporary C: drive use, hourly discharge is ' ...
            'downloaded and installed first; its staging files are then ' ...
            'removed before RDRS forcing is downloaded.' newline newline ...
            'Globus CLI and Globus Connect Personal must be running. ' ...
            'Continue?'];
    else
        msg = [ ...
            'SAGE will use Globus to download Canadian CAMELS-SPAT ' ...
            'metadata, lumped attributes, lumped Daymet forcing, and ' ...
            'daily discharge.' newline newline ...
            'Globus CLI and Globus Connect Personal must be running. ' ...
            'Continue?'];
    end

    if ~local_ui_confirm(ui,'Download Canada data',msg)
        out = struct('ok',false,'canceled',true);
        return
    end

    sourceCollection = ...
        'f163c1b3-9c88-42f6-a7bb-5839ed6c4063';
    sourceRoot = ...
        '/1/published/publication_1301/submitted_data/';
    repoURL = ['https://globus.frdr.ca/file-manager?' ...
        'origin_id=' sourceCollection ...
        '&origin_path=%2F1%2Fpublished%2Fpublication_1301%2F' ...
        'submitted_data%2F'];

    stageRoot = fullfile(local_default_download_dir(), ...
        'CAMELS_CA_Globus_stage');
    local_mkdir(stageRoot);

    d = local_progress_dialog(ui,'Downloading Canada data', ...
        'Checking Globus installation ...');
    taskMetadata = '';
    taskDischarge = '';
    %taskForcing = '';

    try
        globusExe = local_ca_find_globus_cli();
        if isempty(globusExe)
            error('data_helpers:CAGlobusMissing', ...
                ['Globus CLI was not found. Install it and run ' ...
                'globus login once, then retry.']);
        end

        local_ca_ensure_globus_login(globusExe,logFcn);
        destinationCollection = local_ca_local_endpoint_id( ...
            globusExe,logFcn);
        destinationRoot = local_ca_globus_local_path(stageRoot);
        destinationRoot = local_ca_verify_destination( ...
            globusExe,destinationCollection,destinationRoot,stageRoot);
        logFcn(['Globus destination path: ' destinationRoot]);

        instructionFile = local_ca_write_transfer_instructions( ...
            dirD,sourceCollection,sourceRoot,destinationCollection, ...
            destinationRoot);
        logFcn(['Wrote CAMELS-CA Globus transfer details: ' ...
            instructionFile]);

        % Metadata are shared by daily and hourly. Reuse installed copies
        % when possible; otherwise obtain them first.
        fMetaInstalled = fullfile(dirD,'camels-spat-metadata.csv');
        fAttrInstalled = fullfile(dirD,'attributes-lumped.csv');
        if ~(isfile(fMetaInstalled) && isfile(fAttrInstalled))
            local_ca_set_progress(d,true,0, ...
                'Downloading metadata and attribute tables ...');
            batchMeta = fullfile(stageRoot, ...
                'CAMELS_CA_batch_metadata.txt');
            local_ca_write_metadata_batch(batchMeta);
            taskMetadata = local_ca_submit_globus_batch(globusExe, ...
                sourceCollection,sourceRoot,destinationCollection, ...
                destinationRoot,batchMeta,'CAMELS-CA metadata',logFcn);
            local_ca_wait_for_task(globusExe,taskMetadata,d,logFcn, ...
                'Downloading CAMELS-CA metadata');
            local_ca_install_metadata_from_stage(stageRoot,dirD,logFcn);
        end

        basinInfo = local_ca_read_basin_manifest(fMetaInstalled);

        if strcmp(stream,'daily')
            batchBoth = fullfile(stageRoot, ...
                'CAMELS_CA_batch_daily_timeseries.txt');
            local_ca_write_product_batch(batchBoth,basinInfo, ...
                'daily','both');
            local_ca_set_progress(d,true,0, ...
                sprintf(['Submitting %d Daymet and %d daily discharge ' ...
                'files ...'],height(basinInfo),height(basinInfo)));
            taskForcing = local_ca_submit_globus_batch(globusExe, ...
                sourceCollection,sourceRoot,destinationCollection, ...
                destinationRoot,batchBoth, ...
                'CAMELS-CA Daymet and daily discharge',logFcn);
            local_ca_wait_for_task(globusExe,taskForcing,d,logFcn, ...
                'Downloading CAMELS-CA daily files');
            local_ca_install_product_from_stage( ...
                stageRoot,dirD,'daily','both',d,logFcn);
        else
            % 1) Hourly discharge. Install and clear it from C: before RDRS.
            batchQ = fullfile(stageRoot, ...
                'CAMELS_CA_batch_hourly_discharge.txt');
            local_ca_write_product_batch(batchQ,basinInfo, ...
                'hourly','discharge');
            local_ca_set_progress(d,true,0, ...
                sprintf('Submitting %d hourly discharge files ...', ...
                height(basinInfo)));
            taskDischarge = local_ca_submit_globus_batch(globusExe, ...
                sourceCollection,sourceRoot,destinationCollection, ...
                destinationRoot,batchQ, ...
                'CAMELS-CA hourly discharge',logFcn);
            local_ca_wait_for_task(globusExe,taskDischarge,d,logFcn, ...
                'Downloading CAMELS-CA hourly discharge');
            local_ca_install_product_from_stage( ...
                stageRoot,dirD,'hourly','discharge',d,logFcn);
            local_ca_remove_stage_product(stageRoot,'discharge',logFcn);

            % 2) RDRS forcing.
            batchF = fullfile(stageRoot, ...
                'CAMELS_CA_batch_hourly_rdrs.txt');
            local_ca_write_product_batch(batchF,basinInfo, ...
                'hourly','forcing');
            local_ca_set_progress(d,true,0, ...
                sprintf('Submitting %d RDRS forcing files ...', ...
                height(basinInfo)));
            taskForcing = local_ca_submit_globus_batch(globusExe, ...
                sourceCollection,sourceRoot,destinationCollection, ...
                destinationRoot,batchF, ...
                'CAMELS-CA hourly RDRS forcing',logFcn);
            local_ca_wait_for_task(globusExe,taskForcing,d,logFcn, ...
                'Downloading CAMELS-CA RDRS forcing');
            local_ca_install_product_from_stage( ...
                stageRoot,dirD,'hourly','forcing',d,logFcn);
            local_ca_remove_stage_product(stageRoot,'forcing',logFcn);
        end

        ok = local_ca_install_complete(dirD,stream);
        out = local_ca_result(ok,dirD,stageRoot,stream);
        out.repository = repoURL;
        out.sourceCollection = sourceCollection;
        out.sourceRoot = sourceRoot;
        out.destinationCollection = destinationCollection;
        out.taskMetadata = taskMetadata;
        out.taskDischarge = taskDischarge;
        out.taskForcing = taskForcing;

        nExpected = local_ca_expected_basin_count(dirD);
        if strcmp(stream,'daily')
            nF = local_count_files(dirF,{'CAN_*_daymet_lumped.nc'});
            nQ = local_count_files(dirQ, ...
                {'CAN_*_daily_flow_observations.nc'});
        else
            nF = local_count_files(dirF,{'CAN_*_rdrs_lumped.nc'});
            nQ = local_count_files(dirQ, ...
                {'CAN_*_hourly_flow_observations.nc'});
        end

        if ok
            logFcn(sprintf(['CAMELS-CA %s download finished: %d/%d %s ' ...
                'forcing files and %d/%d discharge files.'], ...
                stream,nF,nExpected,forcingLabel,nQ,nExpected));
            local_ca_cleanup_stage(stageRoot,logFcn);
            out.stageRoot = '';
        else
            logFcn(sprintf(['CAMELS-CA %s download is incomplete: ' ...
                '%d forcing and %d discharge files found; expected %s.'], ...
                stream,nF,nQ,local_ca_count_text(nExpected)));
            logFcn(['The staging folder was preserved for retry: ' ...
                stageRoot]);
        end

    catch ME
        logFcn('CAMELS-CA automatic download failed:');
        logFcn(ME.message);
        logFcn(['Globus source: ' repoURL]);
        logFcn(['Staging folder preserved at: ' stageRoot]);
        out = struct('ok',false,'error',ME.message, ...
            'repository',repoURL,'stageRoot',stageRoot);
    end

    local_close_progress(d);
end

function exe = local_ca_find_globus_cli()
    exe = '';
    candidates = {'globus','globus.exe'};
    for i = 1:numel(candidates)
        [status,~] = system(sprintf('%s version',candidates{i}));
        if status == 0
            exe = candidates{i};
            return
        end
    end
end

function local_ca_ensure_globus_login(exe,logFcn)
    [status,txt] = system(sprintf('%s whoami',exe));
    if status == 0
        logFcn(['Globus login active for: ' strtrim(txt)]);
        return
    end
    logFcn('No active Globus CLI login; starting authentication.');
    [status,txt] = system(sprintf('%s login',exe));
    if status ~= 0
        error('data_helpers:CAGlobusLogin', ...
            'Globus login failed: %s',strtrim(txt));
    end
end

function id = local_ca_local_endpoint_id(exe,logFcn)
    [status,txt] = system(sprintf('%s endpoint local-id --quiet',exe));
    id = strtrim(txt);
    if status ~= 0 || isempty(regexp(id, ...
            '^[0-9a-fA-F-]{36}$','once'))
        error('data_helpers:CANoLocalGlobusEndpoint', ...
            ['Could not identify a local Globus Connect Personal ' ...
            'collection. Install and start Globus Connect Personal.']);
    end
    logFcn(['Local Globus collection: ' id]);
end

function p = local_ca_globus_local_path(localPath)
    localPath = char(localPath);
    a = strrep(localPath,'\','/');
    homeDirs = {getenv('USERPROFILE'),getenv('HOME')};
    for i = 1:numel(homeDirs)
        h = strrep(char(homeDirs{i}),'\','/');
        h = regexprep(h,'/+$','');
        if isempty(h)
            continue
        end
        if strcmpi(a,h)
            p = '/~/';
            return
        end
        prefix = [h '/'];
        if startsWith(a,prefix,'IgnoreCase',true)
            relativePath = a(numel(prefix)+1:end);
            p = ['/~/' relativePath];
            if p(end) ~= '/'
                p = sprintf('%s/',p);
            end
            return
        end
    end
    tok = regexp(a,'^([A-Za-z]):/(.*)$','tokens','once');
    if ~isempty(tok)
        p = ['/' upper(tok{1}) '/' tok{2}];
        if p(end) ~= '/'
            p = [p '/'];
        end
        return
    end
    error('data_helpers:CABadLocalPath', ...
        ['Could not convert the local staging folder to a Globus path: ' ...
         '%s'],localPath);
end

function selectedPath = local_ca_verify_destination( ...
    exe,collection,preferredPath,localPath)

    candidates = {preferredPath,local_ca_globus_drive_path(localPath)};
    candidates = unique(candidates,'stable');
    details = strings(numel(candidates),1);
    for i = 1:numel(candidates)
        pathName = candidates{i};
        cmd = sprintf('%s ls "%s:%s"',exe,collection,pathName);
        [status,txt] = system(cmd);
        if status == 0
            selectedPath = pathName;
            return
        end
        details(i) = string(strtrim(txt));
    end
    error('data_helpers:CAGlobusDestination', ...
        ['Globus Connect Personal cannot access the staging folder. ' ...
        'Tried %s. Confirm that Globus Connect Personal is running ' ...
        'under the current Windows account. In Options > Access, add ' ...
        'the folder %s with read and write access. Globus reported: %s'], ...
        strjoin(string(candidates),', '),localPath, ...
        strjoin(details,' | '));
end

function p = local_ca_globus_drive_path(localPath)
    a = strrep(char(localPath),'\','/');
    tok = regexp(a,'^([A-Za-z]):/(.*)$','tokens','once');
    if isempty(tok)
        p = a;
        return
    end
    p = ['/' upper(tok{1}) '/' tok{2}];
    if p(end) ~= '/'
        p = [p '/'];
    end
end

function local_ca_write_metadata_batch(fileName)
    fid = fopen(fileName,'w');
    if fid < 0
        error('data_helpers:CABatchWrite','Could not write %s.',fileName);
    end
    c = onCleanup(@() fclose(fid));
    fprintf(fid,'camels-spat-metadata.csv camels-spat-metadata.csv\n');
    fprintf(fid,['attributes/attributes-lumped.csv ' ...
        'attributes/attributes-lumped.csv\n']);
end

function T = local_ca_read_basin_manifest(metadataFile)
    M = readtable(metadataFile,'VariableNamingRule','preserve');
    vn = string(M.Properties.VariableNames);
    jCountry = find(strcmpi(vn,'Country'),1);
    jID = find(strcmpi(vn,'Station_id'),1);
    jCat = find(strcmpi(vn,'subset_category'),1);
    if isempty(jCountry) || isempty(jID) || isempty(jCat)
        error('data_helpers:CABadMetadata', ...
            ['Metadata must contain Country, Station_id, and ' ...
            'subset_category columns.']);
    end
    country = upper(strtrim(string(M{:,jCountry})));
    keep = country == "CAN" | country == "CANADA";
    station = upper(strtrim(string(M{keep,jID})));
    category = lower(strtrim(string(M{keep,jCat})));
    basin = regexprep(station,'^CAN[-_]?','');
    basin = "CAN_" + basin;
    validCat = ismember(category, ...
        ["headwater","meso-scale","macro-scale"]);
    if any(~validCat)
        bad = unique(category(~validCat));
        error('data_helpers:CABadCategory', ...
            'Unknown subset_category: %s',char(strjoin(bad,', ')));
    end
    T = unique(table(basin(:),category(:), ...
        'VariableNames',{'basin','category'}),'rows','stable');
    if isempty(T)
        error('data_helpers:CANoBasins', ...
            'No Canadian basins were found in metadata.');
    end
end

function local_ca_write_product_batch(fileName,T,stream,product)
    fid = fopen(fileName,'w');
    if fid < 0
        error('data_helpers:CABatchWrite','Could not write %s.',fileName);
    end
    c = onCleanup(@() fclose(fid));
    for i = 1:height(T)
        b = char(T.basin(i));
        cat = char(T.category(i));
        if strcmp(stream,'daily')
            if any(strcmp(product,{'forcing','both'}))
                src = sprintf(['forcing/%s/daymet/daymet-lumped/' ...
                    '%s_daymet_lumped.nc'],cat,b);
                dst = sprintf('forcing/%s_daymet_lumped.nc',b);
                fprintf(fid,'%s %s\n',src,dst);
            end
            if any(strcmp(product,{'discharge','both'}))
                src = sprintf(['observations/%s/obs-daily/' ...
                    '%s_daily_flow_observations.nc'],cat,b);
                dst = sprintf('discharge/%s_daily_flow_observations.nc',b);
                fprintf(fid,'%s %s\n',src,dst);
            end
        else
            if any(strcmp(product,{'forcing','both'}))
                src = sprintf(['forcing/%s/rdrs/rdrs-lumped/' ...
                    '%s_rdrs_lumped.nc'],cat,b);
                dst = sprintf('forcing/%s_rdrs_lumped.nc',b);
                fprintf(fid,'%s %s\n',src,dst);
            end
            if any(strcmp(product,{'discharge','both'}))
                src = sprintf(['observations/%s/obs-hourly/' ...
                    '%s_hourly_flow_observations.nc'],cat,b);
                dst = sprintf('discharge/%s_hourly_flow_observations.nc',b);
                fprintf(fid,'%s %s\n',src,dst);
            end
        end
    end
end

function taskID = local_ca_submit_globus_batch(exe,srcID,srcRoot, ...
    dstID,dstRoot,batchFile,label,logFcn)
    cmd = sprintf([ ...
        '%s transfer "%s:%s" "%s:%s" --batch "%s" ' ...
        '--sync-level size --label "%s" --notify off ' ...
        '--jmespath task_id --format UNIX'], ...
        exe,srcID,srcRoot,dstID,dstRoot,batchFile,label);
    [status,txt] = system(cmd);
    taskID = strtrim(txt);
    if status ~= 0 || isempty(regexp(taskID, ...
            '^[0-9a-fA-F-]{36}$','once'))
        error('data_helpers:CAGlobusSubmit', ...
            'Globus transfer submission failed: %s',strtrim(txt));
    end
    logFcn(sprintf('%s task submitted: %s',label,taskID));
end

function local_ca_wait_for_task(exe,taskID,d,logFcn,label)
    while true
        cmd = sprintf(['%s task show %s --jmespath ' ...
            'status --format UNIX'],exe,taskID);
        [status,txt] = system(cmd);
        if status ~= 0
            error('data_helpers:CAGlobusStatus', ...
                'Could not query Globus task %s: %s', ...
                taskID,strtrim(txt));
        end
        taskStatus = upper(strtrim(txt));
        local_ca_set_progress(d,true,0, ...
            sprintf('%s: %s',label,taskStatus));
        switch taskStatus
            case 'SUCCEEDED'
                logFcn([label ' completed successfully.']);
                return
            case 'FAILED'
                [~,detail] = system(sprintf('%s task show %s', ...
                    exe,taskID));
                error('data_helpers:CAGlobusFailed', ...
                    'Globus task failed: %s',strtrim(detail));
            otherwise
                pause(5);
        end
    end
end

function local_ca_set_progress(d,indeterminate,value,msg)
    try
        if ~isempty(d) && isvalid(d)
            d.Indeterminate = logical(indeterminate);
            if ~indeterminate
                d.Value = min(1,max(0,value));
            end
            d.Message = msg;
            drawnow limitrate
        end
    catch
    end
end

function local_ca_install_metadata_from_stage(srcRoot,dirD,logFcn)
    fMeta = fullfile(srcRoot,'camels-spat-metadata.csv');
    fAttr = fullfile(srcRoot,'attributes','attributes-lumped.csv');
    if ~isfile(fMeta)
        error('data_helpers:CAMissingMetadata', ...
            'Could not find camels-spat-metadata.csv in staging.');
    end
    if ~isfile(fAttr)
        error('data_helpers:CAMissingAttributes', ...
            'Could not find attributes-lumped.csv in staging: %s',fAttr);
    end
    copyfile(fMeta,fullfile(dirD,'camels-spat-metadata.csv'),'f');
    copyfile(fAttr,fullfile(dirD,'attributes-lumped.csv'),'f');
    logFcn('Installed CAMELS-SPAT metadata and lumped attributes.');
end

function local_ca_install_product_from_stage( ...
        srcRoot,dirD,stream,product,d,logFcn)
    if strcmp(stream,'daily')
        dirF = fullfile(dirD,'daily','forcing','daymet');
        dirQ = fullfile(dirD,'daily','discharge');
        forcePattern = 'CAN_*_daymet_lumped.nc';
        qPattern = 'CAN_*_daily_flow_observations.nc';
        forceMsg = 'Installing Daymet forcing';
        qMsg = 'Installing daily discharge';
    else
        dirF = fullfile(dirD,'hourly','forcing','rdrs');
        dirQ = fullfile(dirD,'hourly','discharge');
        forcePattern = 'CAN_*_rdrs_lumped.nc';
        qPattern = 'CAN_*_hourly_flow_observations.nc';
        forceMsg = 'Installing RDRS forcing';
        qMsg = 'Installing hourly discharge';
    end
    local_mkdir(dirF);
    local_mkdir(dirQ);

    forcing = {};
    discharge = {};
    if any(strcmp(product,{'forcing','both'}))
        forcing = local_ca_find_all(fullfile(srcRoot,'forcing'), ...
            forcePattern);
    end
    if any(strcmp(product,{'discharge','both'}))
        discharge = local_ca_find_all(fullfile(srcRoot,'discharge'), ...
            qPattern);
    end
    nAll = numel(forcing) + numel(discharge);
    nDone = 0;
    for i = 1:numel(forcing)
        movefile(forcing{i},fullfile(dirF, ...
            local_ca_file_name(forcing{i})),'f');
        nDone = nDone + 1;
        local_ca_progress(d,nDone,nAll,forceMsg);
    end
    for i = 1:numel(discharge)
        movefile(discharge{i},fullfile(dirQ, ...
            local_ca_file_name(discharge{i})),'f');
        nDone = nDone + 1;
        local_ca_progress(d,nDone,nAll,qMsg);
    end
    logFcn(sprintf(['Installed %d %s forcing files and %d %s ' ...
        'discharge files.'],numel(forcing),stream, ...
        numel(discharge),stream));
end

function local_ca_remove_stage_product(stageRoot,product,logFcn)
    d = fullfile(stageRoot,product);
    try
        if isfolder(d)
            rmdir(d,'s');
            logFcn(['Removed CAMELS-CA staging product: ' d]);
        end
    catch ME
        logFcn(['Could not remove staging product ' d ': ' ME.message]);
    end
end

function tf = local_ca_install_complete(dirD,stream)
    tf = isfile(fullfile(dirD,'camels-spat-metadata.csv')) ...
        && isfile(fullfile(dirD,'attributes-lumped.csv'));
    if ~tf, return, end
    nExpected = local_ca_expected_basin_count(dirD);
    if strcmp(stream,'daily')
        nF = local_count_files(fullfile(dirD,'daily','forcing','daymet'), ...
            {'CAN_*_daymet_lumped.nc'});
        nQ = local_count_files(fullfile(dirD,'daily','discharge'), ...
            {'CAN_*_daily_flow_observations.nc'});
    else
        nF = local_count_files(fullfile(dirD,'hourly','forcing','rdrs'), ...
            {'CAN_*_rdrs_lumped.nc'});
        nQ = local_count_files(fullfile(dirD,'hourly','discharge'), ...
            {'CAN_*_hourly_flow_observations.nc'});
    end
    tf = isfinite(nExpected) && nExpected > 0 ...
        && nF >= nExpected && nQ >= nExpected;
end

function n = local_ca_expected_basin_count(dirD)
    n = NaN;
    f = fullfile(dirD,'camels-spat-metadata.csv');
    if ~isfile(f), return, end
    try
        T = readtable(f,'VariableNamingRule','preserve');
        vn = string(T.Properties.VariableNames);
        j = find(strcmpi(vn,'Country'),1);
        if ~isempty(j)
            country = upper(strtrim(string(T{:,j})));
            n = sum(country == "CAN" | country == "CANADA");
        end
    catch
        n = NaN;
    end
end

function f = local_ca_write_transfer_instructions(dirD,srcID,srcRoot, ...
    dstID,dstRoot)
    f = fullfile(dirD,'CAMELS_CA_Globus_instructions.txt');
    fid = fopen(f,'w');
    if fid < 0
        error('data_helpers:CAInstructions','Could not write %s.',f);
    end
    c = onCleanup(@() fclose(fid));
    fprintf(fid,'CAMELS-CA automatic selective Globus transfer\n\n');
    fprintf(fid,'Source collection: %s\n',srcID);
    fprintf(fid,'Source root: %s\n',srcRoot);
    fprintf(fid,'Destination collection: %s\n',dstID);
    fprintf(fid,'Destination root: %s\n',dstRoot);
end

function files = local_ca_find_all(root,pattern)
    files = cell(0,1);
    if ~isfolder(root), return, end
    D = dir(fullfile(root,'**',pattern));
    for i = 1:numel(D)
        if ~D(i).isdir
            files{end+1,1} = fullfile(D(i).folder,D(i).name); %#ok<AGROW>
        end
    end
end

function local_ca_progress(d,k,n,msg)
    try
        if ~isempty(d) && isvalid(d)
            if n > 0
                d.Indeterminate = false;
                d.Value = min(1,max(0,k/n));
            else
                d.Indeterminate = true;
            end
            d.Message = sprintf('%s (%d/%d)',msg,k,n);
            drawnow limitrate
        end
    catch
    end
end

function name = local_ca_file_name(pathName)
    [~,n,e] = fileparts(pathName);
    name = [n e];
end

function s = local_ca_count_text(n)
    if isfinite(n), s = sprintf('%d',n); else, s = 'unknown'; end
end

function local_ca_cleanup_stage(stageRoot,logFcn)
    try
        if isfolder(stageRoot)
            rmdir(stageRoot,'s');
            logFcn(['Removed CAMELS-CA staging folder: ' stageRoot]);
        end
    catch ME
        logFcn(['Could not remove CAMELS-CA staging folder: ' ...
            ME.message]);
    end
end

function out = local_ca_result(ok,dirD,stageRoot,stream)
    if strcmp(stream,'daily')
        dirM = fullfile(dirD,'daily','forcing','daymet');
        dirQ = fullfile(dirD,'daily','discharge');
    else
        dirM = fullfile(dirD,'hourly','forcing','rdrs');
        dirQ = fullfile(dirD,'hourly','discharge');
    end
    out = struct('ok',logical(ok), ...
        'region','CAMELS_CA','dirD',dirD, ...
        'dirM',dirM,'dirQ',dirQ, ...
        'stream',stream,'stageRoot',stageRoot);
end

% =============================================
function out = local_download_us(cfg,stream,ui)
% =============================================
    root = local_cfg_root(cfg);
    logFcn = local_ui_log(ui);
    
    if strcmp(stream,'hourly')
        if ~local_ui_confirm(ui, ...
                'Download hourly data', ...
                ['Hourly CAMELS downloads ' ...
                'are very large (many GB).' newline ...
                 'Continue?'])
            out = struct('ok',false, ...
                'canceled',true);
            return
        end
    end
    
    try
        addpath(fileparts(mfilename('fullpath')));
    catch
    end
    
    if ~exist('install_SAGEhydrology','file')
        logFcn(['install_SAGEhydrology.m ' ...
            'not on path.']);
        logFcn(['Put install_SAGEhydrology.m ' ...
            'next to SAGE_ui ' ...
            'or add it to the MATLAB path.']);
        out = struct('ok',false,'reason', ...
            'missing installer');
        return
    end
    
    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        sprintf('Downloading %s data',stream), ...
        sprintf('Starting %s install ...',stream));
    
    opts = struct();
    opts.region = 'CAMELS_US';
    opts.force = false;
    % Do not leave downloaded archives behind after a successful install.
    opts.keepArchives = false;
    opts.addToPath = false;
    opts.logFcn = logFcn;
    opts.cancelFcn = @() userCanceled;
    opts.progressFcn = @(info) ui_progress(info);
    
    try
        logFcn(sprintf(['Downloading %s ' ...
            'CAMELS-US data ...'],stream));
        out = install_SAGEhydrology(root, ...
            stream,opts);
        ok = isstruct(out) ...
            && isfield(out,'ok') ...
            && out.ok;
        if ok
            try
                local_migrate_camels_us_data(root);
            catch ME
                logFcn(['US data migration ' ...
                    'warning: ' ME.message]);
            end
            logFcn(sprintf(['%s data install ' ...
                'finished.'],local_cap(stream)));
        else
            logFcn(sprintf(['%s data install ' ...
                'did not complete successfully.'], ...
                local_cap(stream)));
        end
    catch ME
        logFcn(sprintf('%s install failed:', ...
            local_cap(stream)));
        logFcn(ME.message);
        out = struct('ok',false,'error', ...
            ME.message);
    end
    
    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

% ===================
% other local helpers
% ===================
function files = local_camels_fr_files()

    files = struct([]);

    files(1).key = 'attributes';
    files(1).label = 'CAMELS-FR attributes';
    files(1).fileName = 'CAMELS_FR_attributes.zip';
    files(1).url = ['https://entrepot.recherche.data.gouv.fr/' ...
        'api/access/datafile/:persistentId?persistentId=' ...
        'doi:10.57745/GLLWWH&format=original'];
    files(1).md5 = '807c82e51d5380dc1f635b00711a3a79';

    files(2).key = 'timeseries';
    files(2).label = 'CAMELS-FR daily time series';
    files(2).fileName = 'CAMELS_FR_time_series.zip';
    files(2).url = ['https://entrepot.recherche.data.gouv.fr/' ...
        'api/access/datafile/:persistentId?persistentId=' ...
        'doi:10.57745/ZLLPKI&format=original'];
    files(2).md5 = 'dd48efe7cca89e86d8435a9888ebcdca';

end

function local_install_fr_files(unzipRoot,dirD,logFcn)
%LOCAL_INSTALL_FR_FILES Install CAMELS-FR files from the unzipped archives.
    if nargin < 3 ...
            || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n', ...
            char(string(s)));
    end

    if isempty(unzipRoot) ...
            || ~isfolder(unzipRoot)
        error(['data_helpers:' ...
            'FRBadUnzipRoot'], ...
            ['CAMELS-FR unzip root ' ...
            'does not exist: %s'],unzipRoot);
    end

    attrDir = local_find_dir_recursive( ...
        unzipRoot,'static_attributes');
    if isempty(attrDir)
        fProbe = local_find_file_recursive_manual( ...
            unzipRoot, ...
            ['CAMELS_FR_station_' ...
            'general_attributes.csv']);
        if ~isempty(fProbe)
            attrDir = fileparts(fProbe);
        end
    end
    if isempty(attrDir) ...
        || ~isfolder(attrDir)
        error(['data_helpers:' ...
            'FRMissingAttributes'], ...
            ['Could not find static_' ...
            'attributes or the CAMELS-FR ' ...
            'attribute CSV files inside ' ...
            'the attributes archive.']);
    end

    attrFiles = local_fr_metadata_files();
    for i = 1:numel(attrFiles)
        src = fullfile(attrDir, ...
            attrFiles{i});
        if ~isfile(src)
            src = local_find_file_recursive_manual( ...
                unzipRoot,attrFiles{i});
        end
        if isempty(src) ...
                || ~isfile(src)
            error(['data_helpers:' ...
                'FRMissingAttributeFile'], ...
                ['Could not find ' ...
                'required CAMELS-FR file: %s'], ...
                attrFiles{i});
        end
        copyfile(src,fullfile(dirD, ...
            attrFiles{i}),'f');
        logFcn(['Copied ' attrFiles{i} ' to: ' dirD]);
    end

    tsDir = local_find_dir_recursive(unzipRoot,'daily');
    if isempty(tsDir)
        fProbe = local_find_file_recursive_manual(unzipRoot, ...
            'CAMELS_FR_tsd_A105003001.csv');
        if ~isempty(fProbe)
            tsDir = fileparts(fProbe);
        end
    end
    if isempty(tsDir) ...
            || ~isfolder(tsDir)
        error(['data_helpers:' ...
            'FRMissingTimeseries'], ...
            ['Could not find the ' ...
            'daily CAMELS-FR time-series ' ...
            'folder inside the ' ...
            'time-series archive.']);
    end

    tsDst = fullfile(dirD,'daily','timeseries');
    local_mkdir(fileparts(tsDst));
    local_copy_folder_contents(tsDir, ...
        tsDst,logFcn);
end

function tf = local_fr_install_complete(dirD)
    tf = local_all_files_exist(dirD, ...
        local_fr_metadata_files()) ...
        && isfolder(fullfile(dirD, ...
        'daily','timeseries')) ...
        && local_count_files(fullfile(dirD, ...
        'daily','timeseries'), ...
        {'CAMELS_FR_tsd_*.csv'}) >= 50;
end



function files = local_fi_metadata_files()
    files = { ...
        'CAMELS_FI_climatic_attributes.csv', ...
        'CAMELS_FI_geology_attributes.csv', ...
        'CAMELS_FI_humaninfluence_attributes.csv', ...
        'CAMELS_FI_hydrologic_attributes.csv', ...
        'CAMELS_FI_landcover_attributes.csv', ...
        'CAMELS_FI_meta_attributes.csv', ...
        'CAMELS_FI_soil_attributes.csv', ...
        'CAMELS_FI_topographic_attributes.csv'};
end

function tf = local_fi_install_complete(dirD)
    tf = local_all_files_exist(dirD, ...
        local_fi_metadata_files()) ...
        && isfolder(fullfile(dirD, ...
        'daily','timeseries')) ...
        && local_count_files(fullfile(dirD, ...
        'daily','timeseries'), ...
        {'CAMELS_FI_hydromet_timeseries_*.csv'}) > 0;
end

function local_install_fi_files(unzipRoot,dirD,logFcn)
%LOCAL_INSTALL_FI_FILES Install CAMELS-FI files from the unzipped archive.
    if nargin < 3 ...
            || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n', ...
            char(string(s)));
    end

    if isempty(unzipRoot) ...
            || ~isfolder(unzipRoot)
        error(['data_helpers:' ...
            'FIBadUnzipRoot'], ...
            ['CAMELS-FI unzip root ' ...
            'does not exist: %s'], ...
            unzipRoot);
    end

    dataDir = fullfile(unzipRoot, ...
        'CAMELS-FI','data');
    if ~isfolder(dataDir)
        dataDir = local_find_dir_recursive( ...
            unzipRoot,'data');
    end
    if isempty(dataDir) ...
            || ~isfolder(dataDir)
        fProbe = local_find_file_recursive_manual( ...
            unzipRoot, ...
            'CAMELS_FI_meta_attributes.csv');
        if ~isempty(fProbe)
            dataDir = fileparts(fProbe);
        end
    end
    if isempty(dataDir) ...
            || ~isfolder(dataDir)
        error(['data_helpers:' ...
            'FIMissingDataDir'], ...
            ['Could not find CAMELS-FI/' ...
            'data or the required ' ...
            'attribute CSV files ' ...
            'inside the archive.']);
    end

    attrFiles = local_fi_metadata_files();
    for i = 1:numel(attrFiles)
        src = fullfile(dataDir,attrFiles{i});
        if ~isfile(src)
            src = local_find_file_recursive_manual( ...
                unzipRoot,attrFiles{i});
        end
        if isempty(src) ...
                || ~isfile(src)
            error(['data_helpers:' ...
                'FIMissingAttributeFile'], ...
                ['Could not find ' ...
                'required CAMELS-FI ' ...
                'attribute file: %s'],attrFiles{i});
        end
        copyfile(src,fullfile(dirD, ...
            attrFiles{i}),'f');
        logFcn(['Copied ' attrFiles{i} ' to: ' dirD]);
    end

    % Optional but useful for mapping/GIS workflows.
    gpkg = ['CAMELS_FI_catchment_' ...
        'boundaries.gpkg'];
    src = fullfile(dataDir,gpkg);
    if ~isfile(src)
        src = local_find_file_recursive_manual( ...
            unzipRoot,gpkg);
    end
    if ~isempty(src) ...
        && isfile(src)
        copyfile(src,fullfile(dirD,gpkg),'f');
        logFcn(['Copied ' gpkg ' to: ' dirD]);
    end

    tsSrc = fullfile(dataDir, ...
        'timeseries');
    if ~isfolder(tsSrc)
        tsSrc = local_find_dir_recursive( ...
            dataDir,'timeseries');
    end
    if isempty(tsSrc) ...
            || ~isfolder(tsSrc)
        fProbe = local_find_file_recursive_manual( ...
            unzipRoot, ...
            ['CAMELS_FI_hydromet_' ...
            'timeseries_1000_' ...
            '19610101-20231231.csv']);
        if ~isempty(fProbe)
            tsSrc = fileparts(fProbe);
        end
    end
    if isempty(tsSrc) || ~isfolder(tsSrc)
        error(['data_helpers:' ...
            'FIMissingTimeseries'], ...
            ['Could not find data/' ...
            'timeseries or the CAMELS-FI ' ...
            'daily time-series CSV ' ...
            'files inside the archive.']);
    end

    tsDst = fullfile(dirD,'daily','timeseries');
    local_mkdir(fileparts(tsDst));
    local_copy_folder_contents(tsSrc,tsDst,logFcn);
end


function files = local_lux_metadata_files()
    files = { ...
        'CAMELS_LUX_climatic_attributes.csv', ...
        'CAMELS_LUX_geologic_attributes.csv', ...
        'CAMELS_LUX_landuse_attributes.csv', ...
        'CAMELS_LUX_meta_attributes.csv', ...
        'CAMELS_LUX_topographic_attributes.csv'};
end

function tf = local_lux_install_complete(dirD,stream)
    stream = lower(strtrim(char(stream)));

    if ~any(strcmp(stream,{'daily','hourly','15min'}))
        tf = false;
        return
    end

    tsDir = fullfile(dirD,stream,'timeseries');
    pats = local_lux_timeseries_patterns(stream);
    expected = local_lux_expected_count(stream);

    tf = local_all_files_exist(dirD, ...
        local_lux_metadata_files()) ...
        && isfolder(tsDir) ...
        && local_count_files(tsDir,pats) >= expected;
end

function expected = local_lux_expected_count(stream)
%LOCAL_LUX_EXPECTED_COUNT Required basin files at each resolution.

    switch lower(strtrim(char(stream)))
        case 'daily'
            expected = 50;
        case {'hourly','15min'}
            expected = 51;
        otherwise
            expected = Inf;
    end
end

function pats = local_lux_timeseries_patterns(stream)
%LOCAL_LUX_TIMESERIES_PATTERNS Accepted CAMELS-LUX filenames by stream.
    stream = lower(strtrim(char(stream)));
    switch stream
        case 'daily'
            pats = {['CAMELS_LUX_hydromet_' ...
                'timeseries__daily_ID_*.csv']};
        case 'hourly'
            pats = {['CAMELS_LUX_hydromet_' ...
                'timeseries_hourly_ID_*.csv']};
        case '15min'
            pats = {['CAMELS_LUX_hydromet_' ...
                'timeseries_15min_ID_*.csv'], ...
                ['CAMELS_LUX_hydromet_' ...
                'timeseries__15min_ID_*.csv'], ...
                'CAMELS_LUX*15min*ID_*.csv'};
        otherwise
            pats = {};
    end
end

function local_install_lux_files(unzipRoot,dirD,logFcn,stream)
%LOCAL_INSTALL_LUX_FILES Install CAMELS-LUX files from the unzipped archive.
    if nargin < 3 ...
            || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n', ...
            char(string(s)));
    end
    if nargin < 4 ...
            || isempty(stream)
        stream = 'hourly';
    end
    stream = lower(strtrim(char(stream)));

    if ~any(strcmp(stream,{'daily','hourly','15min'}))
        error('data_helpers:LUXBadStream', ...
            ['CAMELS-LUX supports stream = daily, ' ...
            'hourly, or 15min.']);
    end

    if isempty(unzipRoot) ...
            || ~isfolder(unzipRoot)
        error(['data_helpers:' ...
            'LUXBadUnzipRoot'], ...
            ['CAMELS-LUX unzip root ' ...
            'does not exist: %s'],unzipRoot);
    end

    dataDir = local_find_child_dir( ...
        unzipRoot,'CAMELS-LUX');
    if isempty(dataDir)
        dataDir = local_find_first_payload_dir( ...
            unzipRoot,'CAMELS-LUX');
    end
    if isempty(dataDir) ...
            || ~isfolder(dataDir)
        fProbe = local_find_file_recursive_manual( ...
            unzipRoot, ...
            'CAMELS_LUX_meta_attributes.csv');
        if ~isempty(fProbe)
            dataDir = fileparts(fProbe);
        end
    end
    if isempty(dataDir) ...
            || ~isfolder(dataDir)
        error(['data_helpers:' ...
            'LUXMissingDataDir'], ...
            ['Could not find ' ...
            'CAMELS-LUX or the required ' ...
            'attribute CSV files ' ...
            'inside the archive.']);
    end

    attrFiles = local_lux_metadata_files();
    for i = 1:numel(attrFiles)
        src = fullfile(dataDir,attrFiles{i});
        if ~isfile(src)
            src = local_find_file_recursive_manual( ...
                unzipRoot,attrFiles{i});
        end
        if isempty(src) ...
                || ~isfile(src)
            error(['data_helpers:' ...
                'LUXMissingAttributeFile'], ...
                ['Could not find ' ...
                'required CAMELS-LUX ' ...
                'attribute file: %s'], ...
                attrFiles{i});
        end
        copyfile(src,fullfile(dirD,attrFiles{i}),'f');
        logFcn(['Copied ' attrFiles{i} ' to: ' dirD]);
    end

    tsSrc = fullfile(dataDir,'timeseries',stream);
    if ~isfolder(tsSrc)
        tsSrc = local_find_dir_recursive(unzipRoot,stream);
    end

    probeNames = local_lux_probe_names(stream);

    if isempty(tsSrc) ...
            || ~isfolder(tsSrc)
        for iProbe = 1:numel(probeNames)
            fProbe = local_find_file_recursive_manual( ...
                unzipRoot,probeNames{iProbe});
            if ~isempty(fProbe)
                tsSrc = fileparts(fProbe);
                break
            end
        end
    end
    if isempty(tsSrc) || ~isfolder(tsSrc)
        error(['data_helpers:' ...
            'LUXMissingTimeseries'], ...
            ['Could not find timeseries/' stream ...
            ' or the CAMELS-LUX ' stream ...
            ' time-series CSV files inside the archive.']);
    end

    tsDst = fullfile(dirD,stream,'timeseries');
    local_mkdir(fileparts(tsDst));
    local_mkdir(tsDst);
    local_mkdir(fileparts(tsDst));
    local_copy_folder_contents(tsSrc,tsDst,logFcn);
end


function names = local_lux_probe_names(stream)
%LOCAL_LUX_PROBE_NAMES Representative filenames used to locate a stream.
    stream = lower(strtrim(char(stream)));
    switch stream
        case 'daily'
            names = {['CAMELS_LUX_hydromet_' ...
                'timeseries__daily_ID_01.csv']};
        case 'hourly'
            names = {['CAMELS_LUX_hydromet_' ...
                'timeseries_hourly_ID_01.csv']};
        case '15min'
            names = {['CAMELS_LUX_hydromet_' ...
                'timeseries_15min_ID_01.csv'], ...
                ['CAMELS_LUX_hydromet_' ...
                'timeseries__15min_ID_01.csv']};
        otherwise
            names = {};
    end
end

function local_install_ind_files(unzipRoot,dirD,logFcn)
%LOCAL_INSTALL_IND_FILES Install CAMELS-IND files from the unzipped archive.
%
% Do not assume any fixed wrapper-folder depth. The Zenodo ZIP may unpack as
%
%   CAMELS_IND_download/CAMELS_IND_All_Catchments/...
%
% or with an additional nested folder on some MATLAB/Windows combinations.
% Therefore this routine searches from the unzip root for either the expected
% folders or the expected payload files.

    if nargin < 3 ...
            || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n', ...
            char(string(s)));
    end

    if isempty(unzipRoot) ...
            || ~isfolder(unzipRoot)
        error(['data_helpers:' ...
            'INDBadUnzipRoot'], ...
            ['CAMELS-IND unzip root ' ...
            'does not exist: %s'], ...
            unzipRoot);
    end

    % -----------------------
    % 1) Catchment attributes
    % -----------------------
    attrSrc = local_find_dir_recursive( ...
        unzipRoot,'attributes_csv');

    % Extra robust fallback: find one required attribute file and use its
    % parent folder. This catches cases where the folder name is different,
    % hidden by an extra wrapper, or not found by folder traversal.
    if isempty(attrSrc)
        fProbe = local_find_file_recursive_manual( ...
            unzipRoot, ...
            'camels_ind_topo.csv');
        if ~isempty(fProbe)
            attrSrc = fileparts(fProbe);
        end
    end

    if isempty(attrSrc) ...
            || ~isfolder(attrSrc)
        local_log_ind_unzip_tree( ...
            unzipRoot,logFcn);
        error(['data_helpers:' ...
            'INDMissingAttributes'], ...
            ['Could not find attributes_csv ' ...
            'or camels_ind_topo.csv ' ...
            'inside the CAMELS-IND archive.']);
    end

    attrFiles = local_ind_metadata_files();
    for i = 1:numel(attrFiles)
        src = fullfile(attrSrc,attrFiles{i});
        if ~isfile(src)
            src = local_find_file_recursive_manual( ...
                attrSrc,attrFiles{i});
        end
        if isempty(src) ...
                || ~isfile(src)
            error(['data_helpers:' ...
                'INDMissingAttributeFile'], ...
                ['Could not find ' ...
                'required CAMELS-IND ' ...
                'attribute file: %s'],attrFiles{i});
        end
        copyfile(src,fullfile( ...
            dirD,attrFiles{i}),'f');
        logFcn(['Copied ' attrFiles{i} ' to: ' dirD]);
    end

    % -----------------------------
    % 2) Daily meteorological files
    % -----------------------------
    forcSrc = local_find_dir_recursive( ...
        unzipRoot, ...
        'catchment_mean_forcings');

    if isempty(forcSrc)
        fProbe = local_find_file_recursive_manual( ...
            unzipRoot,'17001.csv');
        if ~isempty(fProbe)
            forcSrc = fileparts(fProbe);
        end
    end

    if isempty(forcSrc) ...
            || ~isfolder(forcSrc)
        local_log_ind_unzip_tree(unzipRoot,logFcn);
        error('data_helpers:INDMissingForcings', ...
            ['Could not find catchment_mean_forcings ' ...
            'inside the CAMELS-IND archive.']);
    end

    forcDst = fullfile(dirD,'daily', ...
        'catchment_mean_forcings');
    local_copy_folder_contents(forcSrc, ...
        forcDst,logFcn);

    % ---------------------------------
    % 3) Observed streamflow wide CSV
    % ---------------------------------
    qSrc = local_find_file_recursive_manual( ...
        unzipRoot, ...
        'streamflow_observed.csv');

    if isempty(qSrc) ...
            || ~isfile(qSrc)
        local_log_ind_unzip_tree( ...
            unzipRoot,logFcn);
        error(['data_helpers:' ...
            'INDMissingStreamflow'], ...
            ['Could not find ' ...
            'streamflow_observed.csv ' ...
            'inside the CAMELS-IND archive.']);
    end

    qDst = fullfile(dirD,'daily', ...
        'streamflow_timeseries');
    local_mkdir(qDst);
    copyfile(qSrc,fullfile(qDst, ...
        'streamflow_observed.csv'),'f');
    logFcn(['Copied streamflow_observed.csv to: ' qDst]);

end

function local_log_ind_unzip_tree(root,logFcn)
%LOCAL_LOG_IND_UNZIP_TREE Log a shallow directory/file listing for debugging.
    try
        logFcn(['CAMELS-IND unzip root was: ' root]);
        q = dir(root);
        q = q(~ismember({q.name},{'.','..'}));
        n = min(numel(q),40);
        for ii = 1:n
            if q(ii).isdir
                logFcn(['  [dir]  ' ...
                    fullfile(q(ii).folder,q(ii).name)]);
            else
                logFcn(['  [file] ' ...
                    fullfile(q(ii).folder,q(ii).name)]);
            end
        end
    catch
    end
end

function d = local_find_dir_recursive(root,targetName)
%LOCAL_FIND_DIR_RECURSIVE Manually find a subdirectory by exact name.
% This avoids depending on dir(fullfile(root,'**',name)), which can be
% inconsistent across older MATLAB releases and platforms.

    d = '';
    if isempty(root) ...
            || ~isfolder(root)
        return
    end

    targetName = char(string(targetName));

    q = dir(root);
    for i = 1:numel(q)
        if ~q(i).isdir
            continue
        end
        nm = q(i).name;
        if strcmp(nm,'.') ...
                || strcmp(nm,'..')
            continue
        end
        thisDir = fullfile(q(i).folder,nm);
        if strcmpi(nm,targetName)
            d = thisDir;
            return
        end
    end

    for i = 1:numel(q)
        if ~q(i).isdir
            continue
        end
        nm = q(i).name;
        if strcmp(nm,'.') || strcmp(nm,'..')
            continue
        end
        thisDir = fullfile(q(i).folder,nm);
        d = local_find_dir_recursive(thisDir,targetName);
        if ~isempty(d)
            return
        end
    end
end

function f = local_find_file_recursive_manual(root,targetName)
%LOCAL_FIND_FILE_RECURSIVE_MANUAL Manually find a file by exact name.

    f = '';
    if isempty(root) || ~isfolder(root)
        return
    end

    targetName = char(string(targetName));

    q = dir(root);
    for i = 1:numel(q)
        if q(i).isdir
            continue
        end
        if strcmpi(q(i).name,targetName)
            f = fullfile(q(i).folder,q(i).name);
            return
        end
    end

    for i = 1:numel(q)
        if ~q(i).isdir
            continue
        end
        nm = q(i).name;
        if strcmp(nm,'.') ...
                || strcmp(nm,'..')
            continue
        end
        f = local_find_file_recursive_manual( ...
            fullfile(q(i).folder,nm),targetName);
        if ~isempty(f)
            return
        end
    end
end

function tf = local_ind_install_complete(dirD)

    tf = local_all_files_exist(dirD, ...
        local_ind_metadata_files()) ...
        && isfile(fullfile(dirD,'daily', ...
        'streamflow_timeseries', ...
        'streamflow_observed.csv')) ...
        && local_count_files(fullfile(dirD,'daily', ...
        'catchment_mean_forcings'),{'*.csv'}) >= 200;

end

function local_install_de_files(dataDir,dirD,logFcn,stream)

    if nargin < 4 ...
            || isempty(stream)
            stream = 'hourly';
    end
    stream = lower(strtrim(char(stream)));

    attrFiles = local_de_metadata_files(stream);
    for i = 1:numel(attrFiles)
        src = fullfile(dataDir,attrFiles{i});
        if ~isfile(src)
            src = local_find_file_recursive( ...
                dataDir,attrFiles(i));
        end
        if isempty(src) ...
                || ~isfile(src)
            error('data_helpers:DEMissingAttributeFile', ...
                ['Could not find ' ...
                'required CAMELS-DE ' ...
                'file inside archive: %s'], ...
                attrFiles{i});
        end
        copyfile(src,fullfile(dirD, ...
            attrFiles{i}),'f');
        logFcn(['Copied ' attrFiles{i} ' to: ' dirD]);
    end

    switch stream
        case 'daily'
            pat = 'CAMELS_DE_hydromet_timeseries_DE*.csv';
            label = 'daily';
        case 'hourly'
            pat = 'CAMELS_DE_1h_hydromet_timeseries_*.csv';
            label = 'hourly';
        otherwise
            error('data_helpers:DEBadStream', ...
                'CAMELS-DE stream must be daily or hourly.');
    end

    tsSrc = fullfile(dataDir,'timeseries');

    if ~isfolder(tsSrc) ...
            || isempty(dir(fullfile(tsSrc,pat)))
        qTs = dir(fullfile(dataDir,'**',pat));
        qTs = qTs(~[qTs.isdir]);

        if ~isempty(qTs)
            tsSrc = qTs(1).folder;
        else
            tsSrc = '';
        end
    end
    if isempty(tsSrc) ...
            || ~isfolder(tsSrc)
        error('data_helpers:DEMissingTimeseries', ...
            ['Could not find the CAMELS-DE ' ...
            'timeseries directory ' ...
            'inside the archive.']);
    end

    D = dir(fullfile(tsSrc,pat));
    if isempty(D)
        error('data_helpers:DENoTimeseriesFiles', ...
            ['No CAMELS-DE %s time-series ' ...
            'CSV files found in: %s'],label,tsSrc);
    end

    tsDst = fullfile(dirD,'daily','timeseries');
    % Do not delete the destination until we know the source is valid.
    if ~isfolder(tsDst)
        mkdir(tsDst);
    end
    
    % Clear only old matching files, not the entire folder.
    Dold = dir(fullfile(tsDst,pat));
    for j = 1:numel(Dold)
        try
            delete(fullfile(Dold(j).folder,Dold(j).name));
        catch MEdel
            logFcn(['Could not delete old CAMELS-DE file: ' ...
                Dold(j).name ' -- ' MEdel.message]);
        end
    end
    
    nCopied = 0;
    for j = 1:numel(D)
        srcFile = fullfile(D(j).folder,D(j).name);
        dstFile = fullfile(tsDst,D(j).name);
        copyfile(srcFile,dstFile,'f');
        nCopied = nCopied + 1;
    end
    
    logFcn(sprintf(['Copied CAMELS-DE %s timeseries ' ...
        'to: %s (%d files).'], ...
        label,tsDst,nCopied));

    Dcheck = dir(fullfile(tsDst,pat));
    if isempty(Dcheck)
        error('data_helpers:DEInstallCopyFailed', ...
            ['CAMELS-DE %s time-series files were found in %s, ' ...
             'but none were copied to %s.'], ...
             label,tsSrc,tsDst);
    end
    % tsDst = fullfile(dirD,'timeseries');
    % 
    % if isfolder(tsDst)
    %     local_cleanup_download_folder(tsDst,logFcn);
    % end
    % 
    % movefile(tsSrc,tsDst,'f');
    % logFcn(sprintf(['Moved CAMELS-DE %s timeseries ' ...
    %     'folder to: %s (%d files).'], ...
    %     label,tsDst,numel(D)));

end

% ======================================
function files = local_camels_cl_files()
% ======================================
    base = ['https://store.pangaea.de/Publications/' ...
        'Alvarez-Garreton-etal_2018/'];
    
    files = struct('key',{},'label',{}, ...
        'fileName',{},'payload',{},'url',{},'isZip',{});
    
    files(end+1) = local_cl_file('attributes', ...
        'CAMELS-CL catchment attributes and metadata', ...
        '1_CAMELScl_attributes.zip', ...
        '1_CAMELScl_attributes.txt');
    
    files(end+1) = local_cl_file('streamflow', ...
        'CAMELS-CL streamflow in millimeters', ...
        '3_CAMELScl_streamflow_mm.zip', ...
        '3_CAMELScl_streamflow_mm.txt');
    
    files(end+1) = local_cl_file('precipitation', ...
        'CAMELS-CL precipitation from CR2MET', ...
        '4_CAMELScl_precip_cr2met.zip', ...
        '4_CAMELScl_precip_cr2met.txt');
    
    files(end+1) = local_cl_file('precipitation', ...
        'CAMELS-CL precipitation from CHIRPS', ...
        '5_CAMELScl_precip_chirps.zip', ...
        '5_CAMELScl_precip_chirps.txt');
    
    files(end+1) = local_cl_file('precipitation', ...
        'CAMELS-CL precipitation from MSWEP', ...
        '6_CAMELScl_precip_mswep.zip', ...
        '6_CAMELScl_precip_mswep.txt');
    
    files(end+1) = local_cl_file('precipitation', ...
        'CAMELS-CL precipitation from TMPA', ...
        '7_CAMELScl_precip_tmpa.zip', ...
        '7_CAMELScl_precip_tmpa.txt');
    
    files(end+1) = local_cl_file('temperature', ...
        'CAMELS-CL minimum air temperature', ...
        '8_CAMELScl_tmin_cr2met.zip', ...
        '8_CAMELScl_tmin_cr2met.txt');
    
    files(end+1) = local_cl_file('temperature', ...
        'CAMELS-CL maximum air temperature', ...
        '9_CAMELScl_tmax_cr2met.zip', ...
        '9_CAMELScl_tmax_cr2met.txt');
    
    files(end+1) = local_cl_file('temperature', ...
        'CAMELS-CL mean air temperature', ...
        '10_CAMELScl_tmean_cr2met.zip', ...
        '10_CAMELScl_tmean_cr2met.txt');
    
    files(end+1) = local_cl_file('pet', ...
        'CAMELS-CL potential evapotranspiration from MODIS', ...
        '11_CAMELScl_pet_8d_modis.zip', ...
        '11_CAMELScl_pet_8d_modis.txt');
    
    files(end+1) = local_cl_file('pet', ...
        'CAMELS-CL potential evapotranspiration from Hargreaves', ...
        '12_CAMELScl_pet_hargreaves.zip', ...
        '12_CAMELScl_pet_hargreaves.txt');
    
    files(end+1) = local_cl_file('swe', ...
        'CAMELS-CL snow water equivalent', ...
        '13_CAMELScl_swe.zip', ...
        '13_CAMELScl_swe.txt');

    function f = local_cl_file(key,label,fileName,payload)
        f = struct();
        f.key = key;
        f.label = label;
        f.fileName = fileName;
        f.payload = payload;
        f.url = [base fileName];
        f.isZip = true;
    end
end

function tf = local_cl_target_complete(dirD,f)

    switch f.key
        case 'attributes'
            dst = dirD;
        case 'streamflow'
            dst = fullfile(dirD,'daily','streamflow');
        case 'precipitation'
            dst = fullfile(dirD,'daily','precipitation');
        case 'temperature'
            dst = fullfile(dirD,'daily','temperature');
        case 'pet'
            dst = fullfile(dirD,'daily','pet');
        case 'swe'
            dst = fullfile(dirD,'daily','swe');
        otherwise
            dst = dirD;
    end
    
    tf = isfile(fullfile(dst,f.payload));

end

function local_install_cl_files(files,dirD,logFcn)
    
    for i = 1:numel(files)
        f = files(i);
    
        if ~isfield(f,'install') ...
                || ~f.install
            continue
        end
    
        if ~isfield(f,'unzipDir') ...
                || isempty(f.unzipDir)
            error('data_helpers:CLNotUnzipped', ...
                ['CAMELS-CL archive was not ' ...
                'unzipped: %s'],f.fileName);
        end
    
        src = local_find_file_recursive( ...
            f.unzipDir,{f.payload});
    
        if isempty(src)
            error('data_helpers:CLMissingPayload', ...
                ['Could not find %s inside ' ...
                '%s.'],f.payload,f.fileName);
        end
    
        switch f.key
            case 'attributes'
                dst = dirD;
            case 'streamflow'
                dst = fullfile(dirD,'daily','streamflow');
            case 'precipitation'
                dst = fullfile(dirD,'daily','precipitation');
            case 'temperature'
                dst = fullfile(dirD,'daily','temperature');
            case 'pet'
                dst = fullfile(dirD,'daily','pet');
            case 'swe'
                dst = fullfile(dirD,'daily','swe');
            otherwise
                dst = dirD;
        end
    
        local_mkdir(dst);
        copyfile(src,fullfile(dst,f.payload),'f');
    
        logFcn(['Copied ' f.payload ' to: ' dst]);
    end

end

% ======================================
function files = local_camels_au_files()
% ======================================
    base = ['https://download.pangaea.de/' ...
        'dataset/921850/files/'];
    files = struct('key',{},'label',{}, ...
        'fileName',{},'url',{},'isZip',{});
    files(end+1) = local_au_file('metadata', ...
        'CAMELS-AU id/name metadata', ...
        '01_id_name_metadata.zip',true);
    files(end+1) = local_au_file('streamflow', ...
        'CAMELS-AU streamflow', ...
        '03_streamflow.zip',true);
    files(end+1) = local_au_file('hydromet', ...
        'CAMELS-AU hydrometeorology', ...
        '05_hydrometeorology.zip',true);
    files(end+1) = local_au_file('master', ...
        'CAMELS-AU master attributes', ...
        ['CAMELS_AUS_Attributes&' ...
        'Indices_MasterTable.csv'],false);

    function f = local_au_file(key,label,fileName,isZip)
        f = struct();
        f.key = key;
        f.label = label;
        f.fileName = fileName;
        f.url = [base fileName];
        f.isZip = isZip;
    end
end

function tf = local_au_target_complete(dirD,f)
    switch f.key
        case 'metadata'
            tf = isfile(fullfile(dirD, ...
                'id_name_metadata.csv')) ...
                || isfile(fullfile(dirD, ...
                '01_id_name_metadata.csv'));
        case 'streamflow'
            tf = isfile(fullfile(dirD,'daily', ...
                'streamflow','streamflow_mmd.csv'));
        case 'hydromet'
            tf = isfile(fullfile(dirD,'daily', ...
                'precipitation', ...
                'precipitation_SILO.csv')) ...
                && isfile(fullfile(dirD,'daily', ...
                'precipitation', ...
                'precipitation_AWAP.csv')) ...
                && local_count_files( ...
                fullfile(dirD,'daily', ...
                'evaporative_demand'), ...
                {'*.csv'}) >= 8 ...
                && isfile(fullfile(dirD,'daily', ...
                'temperature', ...
                'tmin_SILO.csv')) ...
                && isfile(fullfile(dirD,'daily', ...
                'temperature', ...
                'tmax_SILO.csv')) ...
                && isfile(fullfile(dirD,'daily', ...
                'temperature', ...
                'tmin_AWAP.csv')) ...
                && isfile(fullfile(dirD,'daily', ...
                'temperature', ...
                'tmax_AWAP.csv'));
        case 'master'
            tf = isfile(fullfile(dirD, ...
                ['CAMELS_AUS_Attributes&' ...
                'Indices_MasterTable.csv']));
        otherwise
            tf = false;
    end
end

function local_install_au_metadata(files,dirD,logFcn)
    idx = find(strcmp({files.key},'metadata'),1);
    if isempty(idx) ...
            || ~isfield(files(idx),'install') ...
            || ~files(idx).install
        return
    end
    src = local_find_file_recursive(files(idx).unzipDir, ...
        {'01_id_name_metadata.csv', ...
        'id_name_metadata.csv'});
    if isempty(src)
        error('data_helpers:AUMissingMetadata', ...
            ['Could not find 01_id_name_metadata.csv ' ...
            'in CAMELS-AU metadata archive.']);
    end
    copyfile(src,fullfile(dirD, ...
        '01_id_name_metadata.csv'),'f');
    copyfile(src,fullfile(dirD, ...
        'id_name_metadata.csv'),'f');
    logFcn(['Copied CAMELS-AU id/name ' ...
        'metadata to: ' dirD]);
end

function local_install_au_streamflow(files,dirD,logFcn)
    idx = find(strcmp({files.key},'streamflow'),1);
    if isempty(idx) ...
            || ~isfield(files(idx),'install') ...
            || ~files(idx).install
        return
    end
    src = local_find_file_recursive( ...
        files(idx).unzipDir,{'streamflow_mmd.csv'});
    if isempty(src)
        error('data_helpers:AUMissingStreamflow', ...
            ['Could not find streamflow_mmd.csv ' ...
            'in CAMELS-AU streamflow archive.']);
    end
    dst = fullfile(dirD,'daily','streamflow');
    local_mkdir(dst);
    copyfile(src,fullfile(dst,'streamflow_mmd.csv'),'f');
    logFcn(['Copied CAMELS-AU streamflow to: ' dst]);
end

function local_install_au_hydromet(files, ...
    dirD,logFcn)
    idx = find(strcmp({files.key},'hydromet'),1);
    if isempty(idx) ...
            || ~isfield(files(idx),'install') ...
            || ~files(idx).install
        return
    end
    root = files(idx).unzipDir;
    
    precDst = fullfile(dirD,'daily', ...
        'precipitation');
    petDst = fullfile(dirD,'daily', ...
        'evaporative_demand');
    tmpDst = fullfile(dirD,'daily', ...
        'temperature');
    local_mkdir(precDst);
    local_mkdir(petDst);
    local_mkdir(tmpDst);
    
    local_copy_au_named_files(root,precDst, ...
        {'precipitation_SILO.csv', ...
        'precipitation_AWAP.csv'},logFcn);
    
    local_copy_au_named_files(root,petDst, ...
        {'evap_syn_SILO.csv', ...
         'evap_pan_SILO.csv', ...
         'evap_morton_lake_SILO.csv', ...
         'et_tall_crop_SILO.csv', ...
         'et_short_crop_SILO.csv', ...
         'et_morton_wet_SILO.csv', ...
         'et_morton_point_SILO.csv', ...
         'et_morton_actual_SILO.csv'},logFcn);
    
    local_copy_au_named_files(root,tmpDst, ...
        {'tmin_SILO.csv','tmax_SILO.csv', ...
         'tmin_AWAP.csv','tmax_AWAP.csv'},logFcn);
end

function local_install_au_master(files,dirD,logFcn)
    idx = find(strcmp({files.key},'master'),1);
    if isempty(idx) ...
            || ~isfield(files(idx),'install') ...
            || ~files(idx).install
        return
    end
    src = files(idx).localFile;
    if ~isfile(src)
        error('data_helpers:AUMissingMaster', ...
            ['CAMELS-AU master attributes ' ...
            'file was not downloaded.']);
    end
    copyfile(src,fullfile(dirD, ...
        'CAMELS_AUS_Attributes&Indices_MasterTable.csv'),'f');
    logFcn(['Copied CAMELS-AU ' ...
        'master attributes to: ' dirD]);
end

function local_copy_au_named_files(root,dst,names,logFcn)
    local_mkdir(dst);
    for i = 1:numel(names)
        src = local_find_file_recursive(root,names(i));
        if isempty(src)
            error('data_helpers:AUMissingFile', ...
                ['Could not find %s ' ...
                'in CAMELS-AU archive.'],names{i});
        end
        copyfile(src,fullfile(dst,names{i}),'f');
        logFcn(['Copied ' names{i} ' to: ' dst]);
    end
end

function f = local_find_file_recursive(root,names)
    f = '';
    if ~isfolder(root)
        return
    end
    if ischar(names) ...
            || isstring(names)
        names = cellstr(string(names));
    end
    for i = 1:numel(names)
        q = dir(fullfile(root,'**',names{i}));
        q = q(~[q.isdir]);
        if ~isempty(q)
            f = fullfile(q(1).folder,q(1).name);
            return
        end
    end
end

function files = local_camels_br_files()
    base = 'https://zenodo.org/records/15025488/files/';
    files = struct('key',{},'label',{},'fileName',{}, ...
        'url',{},'md5',{},'isZip',{});
    files(end+1) = local_br_file('attributes', ...
        'CAMELS-BR attributes', ...
        '01_CAMELS_BR_attributes.zip', ...
        '7de033200bdb25e4da9cf40addcd77aa',true);
    files(end+1) = local_br_file('streamflow', ...
        'CAMELS-BR streamflow', ...
        '02_CAMELS_BR_streamflow_all_catchments.zip', ...
        '01ac209ad3ca9a36f7f2510a6601a0a2',true);
    files(end+1) = local_br_file('precipitation', ...
        'CAMELS-BR precipitation', ...
        '05_CAMELS_BR_precipitation.zip', ...
        '8b3611c1fc0597836196fd3a41833855',true);
    files(end+1) = local_br_file('pet', ...
        'CAMELS-BR potential evapotranspiration', ...
        '07_CAMELS_BR_potential_evapotransp.zip', ...
        'a27b41843291af2a9d7dcdf996d575aa',true);
    files(end+1) = local_br_file('temperature', ...
        'CAMELS-BR temperature', ...
        '09_CAMELS_BR_temperature.zip', ...
        '9704d5ce24f414db65bbfc3196beaf6d',true);
    files(end+1) = local_br_file('soil_moisture', ...
        'CAMELS-BR soil moisture', ...
        '10_CAMELS_BR_soil_moisture.zip', ...
        '928ef8eef939e4574ce34aa8f8c89d1e',true);
    files(end+1) = local_br_file('readme', ...
        'CAMELS-BR readme', ...
        'CAMELS_BR_readme.txt', ...
        'fbb422eddabd323d1c7eb64452adb974',false);

    function f = local_br_file(key, ...
            label,fileName,md5,isZip)
        f = struct();
        f.key = key;
        f.label = label;
        f.fileName = fileName;
        f.url = [base fileName '?download=1'];
        f.md5 = md5;
        f.isZip = isZip;
    end
end

function local_install_br_zip(files, ...
    key,dst,patterns,logFcn)
    idx = find(strcmp({files.key},key),1);
    if isempty(idx)
        error('data_helpers:BRMissingSpec', ...
            ['Missing CAMELS-BR file ' ...
            'specification for %s.'],key);
    end
    f = files(idx);
    if ~isfield(f,'unzipDir') ...
            || isempty(f.unzipDir)
        error('data_helpers:BRNotUnzipped', ...
            ['CAMELS-BR archive was ' ...
            'not unzipped: %s'],f.fileName);
    end
    src = local_find_first_payload_dir(f.unzipDir, ...
        erase(f.fileName,'.zip'));
    if isempty(src)
        error('data_helpers:BRMissingPayloadDir', ...
            ['Could not find payload ' ...
            'directory for %s.'],f.fileName);
    end
    local_copy_pattern(src,dst,patterns,logFcn);
end

function local_download_file(url,targetFile, ...
    label,~,cancelFcn,progressFcn,logFcn)
%LOCAL_DOWNLOAD_FILE Robust downloader with Figshare/browser fallback.
% Some Figshare ndownloader links work in Chrome but return zero bytes
% through MATLAB websave/curl on some Windows systems.  In that case we
% open the URL in the default browser and wait until the browser-created
% file appears in Downloads.

    if nargin < 7 ...
            || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n', ...
            char(string(s)));
    end

    if nargin < 5 || isempty(cancelFcn)
        cancelFcn = @() false;
    end

    local_mkdir(fileparts(targetFile));
    targetFile = char(string(targetFile));
    url = char(string(url));

    % Reuse an existing nonempty file.  The caller will verify MD5 when
    % available, so this is safe.
    if isfile(targetFile)
        info = dir(targetFile);
        if info.bytes > 0
            logFcn(['Using existing file: ' targetFile]);
            if ~isempty(progressFcn)
                progressFcn(struct('label', ...
                    [label ' already downloaded'], ...
                    'frac',1,'speedMBs',NaN));
            end
            return
        else
            try
                delete(targetFile);
            catch
            end
        end
    end

    ok = false;
    msg = '';

    % 1) Try MATLAB websave first.
    try
        opts = weboptions('Timeout',Inf, ...
            'UserAgent',['Mozilla/5.0 ' ...
            '(Windows NT 10.0; Win64; x64)']);
        savedFile = websave(targetFile,url,opts);
        if isfile(savedFile) && ~strcmp(savedFile,targetFile)
            try
                movefile(savedFile,targetFile,'f');
            catch
            end
        end
        ok = isfile(targetFile) ...
            && dir(targetFile).bytes > 0;
    catch ME
        msg = ME.message;
    end

    % 2) Try Windows curl.  --ssl-no-revoke avoids common schannel
    % revocation failures on some Windows installations.
    if ~ok
        try
            if cancelFcn()
                error('data_helpers:DownloadCanceled', ...
                    'Download canceled by user.');
            end
            cmd = sprintf(['curl -L --fail --ssl-no-revoke ' ...
                '--retry 3 -A "Mozilla/5.0" ' ...
                '-o "%s" "%s"'],targetFile,url);
            [status,out] = system(cmd);
            ok = status == 0 ...
                && isfile(targetFile) ...
                && dir(targetFile).bytes > 0;
            if ~ok
                msg = [msg ' | curl: ' out];
            end
        catch ME2
            msg = [msg ' | curl: ' ME2.message];
        end
    end

    % 3) Try PowerShell.  Some machines block curl but allow .NET.
    if ~ok
        try
            if cancelFcn()
                error('data_helpers:DownloadCanceled', ...
                    'Download canceled by user.');
            end
            ps = sprintf(['powershell -NoProfile -' ...
                'ExecutionPolicy Bypass ' ...
                '-Command "[Net.ServicePointManager]:' ...
                ':SecurityProtocol = ' ...
                '[Net.SecurityProtocolType]::Tls12; ' ...
                '$wc = New-Object Net.WebClient; ' ...
                '$wc.Headers.Add(User-Agent,Mozilla/5.0); ' ...
                '$wc.DownloadFile(%s,%s)"'], ...
                url,targetFile);
            [status,out] = system(ps);
            ok = status == 0 ...
                && isfile(targetFile) ...
                && dir(targetFile).bytes > 0;
            if ~ok
                msg = [msg ' | powershell: ' out];
            end
        catch ME3
            msg = [msg ' | powershell: ' ME3.message];
        end
    end

    % 4) Last resort for Figshare: open the link in the browser and wait
    % for the browser download to finish.  This matches the behavior that
    % works when the same link is pasted into Chrome.
    if ~ok && contains(url,'figshare.canterbury.ac.nz')
        try
            logFcn(['Command-line download ' ...
                'failed; opening Figshare link ' ...
                'in the default browser.']);
            logFcn(['Please allow the ' ...
                'browser download to finish. ' ...
                'Expected file: ' targetFile]);
            web(url,'-browser');

            ok = local_wait_for_browser_download( ...
                targetFile,240,cancelFcn,progressFcn,label);
            if ~ok
                msg = [msg ' | browser ' ...
                    'fallback timed out.'];
            end
        catch ME4
            msg = [msg ' | browser fallback: ' ME4.message];
        end
    end

    if ~ok
        error('data_helpers:DownloadFailed', ...
            'Could not download %s to %s. %s', ...
            label,targetFile,msg);
    end

    if ~isempty(progressFcn)
        progressFcn(struct('label', ...
            [label ' downloaded'], ...
            'frac',1,'speedMBs',NaN));
    end
end

function ok = local_wait_for_browser_download( ...
    targetFile,maxSeconds,cancelFcn,progressFcn,label)
% Wait until Chrome/Edge/Firefox has completed a browser download.
%
% Figshare/Chrome can save files with slightly different names than the
% name requested by MATLAB, for example:
%   CAMELS_NZ_Catchment_Atrributes
%   CAMELS_NZ_Catchment_Atrributes (1)
%   CAMELS_NZ_Catchment_Atrributes.zip
%   CAMELS_NZ_Catchment_Atrributes (1).zip
% This helper accepts those browser-created variants, moves the valid
% nonzero file to targetFile, and ignores/deletes zero-byte placeholders.

    ok = false;
    t0 = tic;
    lastBytes = -1;
    stableCount = 0;
    targetFile = char(string(targetFile));
    [folder,base,ext] = fileparts(targetFile);

    if isempty(folder)
        folder = pwd;
    end

    while toc(t0) < maxSeconds
        if ~isempty(cancelFcn) && cancelFcn()
            return
        end

        % Candidate names produced by MATLAB/websave and by Chrome.
        cand = {targetFile};
        if isempty(ext)
            cand{end+1} = [targetFile '.zip'];                  %#ok
        else
            cand{end+1} = fullfile(folder,base);                %#ok
        end
        cand{end+1} = fullfile(folder,[base ' (1)' ext]);       %#ok
        if isempty(ext)
            cand{end+1} = fullfile(folder,[base ' (1).zip']);   %#ok
        else
            cand{end+1} = fullfile(folder,[base ' (1)']);       %#ok
        end

        % Also search the folder for additional browser duplicates.
        dd = dir(fullfile(folder,[base '*']));
        for jj = 1:numel(dd)
            if ~dd(jj).isdir
                cand{end+1} = fullfile(folder,dd(jj).name);     %#ok<AGROW>
            end
        end

        % Remove duplicate candidates while preserving order.
        [~,ia] = unique(cand,'stable');
        cand = cand(sort(ia));

        bestFile = '';
        bestBytes = -1;

        for jj = 1:numel(cand)
            f = cand{jj};
            if ~isfile(f)
                continue
            end

            % Ignore active/incomplete browser temp files.
            if endsWith(f,'.crdownload') || endsWith(f,'.tmp')
                continue
            end

            info = dir(f);
            if info.bytes == 0
                % Delete only zero-byte placeholders with matching basename.
                try
                    delete(f);
                catch
                end
                continue
            end

            if info.bytes > bestBytes
                bestFile = f;
                bestBytes = info.bytes;
            end
        end

        if ~isempty(bestFile)
            % Move/rename the browser-created file to the expected target.
            if ~strcmp(bestFile,targetFile)
                try
                    if isfile(targetFile)
                        delete(targetFile);
                    end
                    movefile(bestFile, ...
                        targetFile,'f');
                    logmsg = ['Detected ' ...
                        'browser-downloaded file: ' ...
                        bestFile ' -> ' targetFile];
                    % logFcn is not passed here, so use progress label below.
                catch
                    % If move fails, keep waiting.
                    pause(2)
                    continue
                end
            else
                logmsg = ['Detected browser-' ...
                    'downloaded file: ' targetFile];
            end

            % Require the final file to be stable for a few checks.
            info = dir(targetFile);
            bytes = info.bytes;
            tmp1 = [targetFile '.crdownload'];
            tmp2 = [targetFile '.tmp'];
            if bytes > 0 ...
                    && ~isfile(tmp1) ...
                    && ~isfile(tmp2)
                if bytes == lastBytes
                    stableCount = stableCount + 1;
                else
                    stableCount = 0;
                    lastBytes = bytes;
                end
                if stableCount >= 2
                    ok = true;
                    return
                end
            end
        else
            logmsg = ['Waiting for ' ...
                'browser download: ' label];
        end

        if ~isempty(progressFcn)
            progressFcn(struct('label',logmsg, ...
                'frac',NaN,'speedMBs',NaN));
        end
        pause(2)
    end
end

function local_cleanup_download_files(files,logFcn)
%LOCAL_CLEANUP_DOWNLOAD_FILES Delete downloaded files in file specs.
    if nargin < 2 || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n',char(string(s)));
    end
    if isempty(files)
        return
    end
    for i = 1:numel(files)
        if isfield(files(i),'localFile') ...
                && ~isempty(files(i).localFile)
            local_cleanup_download_file(files(i).localFile,logFcn);
        end
    end
end

function local_cleanup_download_file(fname,logFcn)
%LOCAL_CLEANUP_DOWNLOAD_FILE Delete one downloaded archive/file if present.
    if nargin < 2 ...
            || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n',char(string(s)));
    end
    if isempty(fname)
        return
    end
    fname = char(string(fname));
    if ~local_is_safe_download_cleanup_target(fname)
        logFcn(['Refused to delete a file outside the Downloads ' ...
            'directory: ' fname]);
        return
    end
    try
        if isfile(fname)
            delete(fname);
            logFcn(['Deleted downloaded file: ' fname]);
        end
    catch ME
        logFcn(['Could not delete downloaded file: ' ...
            fname ' | ' ME.message]);
    end
end

function local_cleanup_download_folder(folder,logFcn)
%LOCAL_CLEANUP_DOWNLOAD_FOLDER Delete one temporary extraction folder.
    if nargin < 2 ...
            || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n',char(string(s)));
    end
    if isempty(folder)
        return
    end
    folder = char(string(folder));
    if ~local_is_safe_download_cleanup_target(folder)
        logFcn(['Refused to delete a folder outside the Downloads ' ...
            'directory: ' folder]);
        return
    end
    try
        if isfolder(folder)
            removed = local_remove_dir_retry(folder,logFcn);
            if ~removed
                return
            end
            logFcn(['Deleted temporary download folder: ' folder]);
        end
    catch ME
        logFcn(['Could not delete temporary download folder: ' ...
            folder ' | ' ME.message]);
    end
end

function tf = local_is_safe_download_cleanup_target(target)
%LOCAL_IS_SAFE_DOWNLOAD_CLEANUP_TARGET Restrict installer cleanup scope.
% Cleanup helpers may remove downloaded archives and extraction folders,
% but must never remove a configured SAGE data directory.

    try
        downloadRoot = char(java.io.File( ...
            local_default_download_dir()).getCanonicalPath());
        targetPath = char(java.io.File(target).getCanonicalPath());
    catch
        downloadRoot = char(local_default_download_dir());
        targetPath = char(target);
    end

    downloadRoot = strrep(downloadRoot,'/',filesep);
    targetPath = strrep(targetPath,'/',filesep);
    downloadRoot = regexprep(downloadRoot, ...
        [regexptranslate('escape',filesep) '+$'],'');
    targetPath = regexprep(targetPath, ...
        [regexptranslate('escape',filesep) '+$'],'');

    prefix = [downloadRoot filesep];
    tf = startsWith(targetPath,prefix,'IgnoreCase',true) ...
        && ~strcmpi(targetPath,downloadRoot);
end

function local_verify_md5(fileName, ...
    expected,logFcn)
    actual = local_file_md5(fileName);
    if ~strcmpi(actual,expected)
        error('data_helpers:MD5Mismatch', ...
            ['MD5 mismatch for %s. Expected %s, ' ...
            'but found %s. Delete ' ...
            'the file and retry.'], ...
            fileName,expected,actual);
    end
    logFcn(sprintf('Verified MD5: %s', ...
        fileName));
end

function h = local_file_md5(fileName)
    md = java.security.MessageDigest.getInstance('MD5');

    fid = -1;
    for k = 1:20
        fid = fopen(fileName,'rb');
        if fid >= 0
            break
        end
        pause(0.25)
    end

    if fid < 0
        error('data_helpers:MD5OpenFailed', ...
            ['Could not open file for ' ...
            'MD5 check after waiting: %s'], ...
            fileName);
    end
    c = onCleanup(@() fclose(fid));
    while true
        data = fread(fid,1024*1024,'*uint8');
        if isempty(data)
            break
        end
        md.update(typecast(data(:),'int8'));
    end
    d = typecast(md.digest(),'uint8');
    h = lower(reshape(dec2hex(d,2).',1,[]));
end

function d = local_default_download_dir()
%LOCAL_DEFAULT_DOWNLOAD_DIR Return a writable download directory.

    if ispc
        home = getenv('USERPROFILE');
    else
        home = getenv('HOME');
    end

    if isempty(home)
        try
            home = char(java.lang.System.getProperty('user.home'));
        catch
            home = '';
        end
    end

    candidates = {};
    if ~isempty(home)
        candidates{end+1} = fullfile(home,'Downloads');
        candidates{end+1} = home;
    end
    candidates{end+1} = tempdir;

    %d = '';
    for i = 1:numel(candidates)
        candidate = candidates{i};
        try
            if ~isfolder(candidate)
                mkdir(candidate);
            end
            testFile = fullfile(candidate, ...
                sprintf('.sage_write_test_%d.tmp',matlabProcessID));
            fid = fopen(testFile,'w');
            if fid >= 0
                fclose(fid);
                delete(testFile);
                d = candidate;
                return
            end
        catch
        end
    end

    error('data_helpers:noWritableDownloadDirectory', ...
        'Could not locate a writable download directory.');
end

function local_mkdir(d)
    if ~isempty(d) ...
            && ~isfolder(d)
        mkdir(d);
    end
end

function d = local_find_child_dir(root,name)
    d = '';
    if ~isfolder(root)
        return
    end
    q = dir(fullfile(root,'**',name));
    for i = 1:numel(q)
        if q(i).isdir
            d = fullfile(q(i).folder, ...
                q(i).name);
            return
        end
    end
end

function d = local_find_first_payload_dir( ...
    root,preferredName)
    d = '';
    if ~isfolder(root)
        return
    end
    p = fullfile(root,preferredName);
    if isfolder(p)
        d = p;
        return
    end
    q = dir(root);
    q = q([q.isdir]);
    for i = 1:numel(q)
        nm = q(i).name;
        if strcmp(nm,'.') ...
                || strcmp(nm,'..')
            continue
        end
        d = fullfile(root,nm);
        return
    end
    d = root;
end

function local_copy_pattern(src,dst,patterns,logFcn)
    local_mkdir(dst);
    n = 0;
    for i = 1:numel(patterns)
        q = dir(fullfile(src,patterns{i}));
        for j = 1:numel(q)
            if q(j).isdir
                continue
            end
            copyfile(fullfile(q(j).folder,q(j).name), ...
                fullfile(dst,q(j).name),'f');
            n = n + 1;
        end
    end
    logFcn(sprintf(['Copied %d files ' ...
        'from %s to %s'], ...
        n,src,dst));
end

function local_copy_folder_contents(src,dst,logFcn)
    if isfolder(dst)
        try
            rmdir(dst,'s');
        catch
        end
    end
    local_mkdir(dst);
    q = dir(src);
    n = 0;
    for i = 1:numel(q)
        nm = q(i).name;
        if strcmp(nm,'.') ...
                || strcmp(nm,'..')
            continue
        end
        copyfile(fullfile(q(i).folder,nm), ...
            fullfile(dst,nm),'f');
        n = n + 1;
    end
    logFcn(sprintf(['Copied %d entries ' ...
        'from %s to %s'], ...
        n,src,dst));
end

% ================================
% Directories and region utilities
% ================================
function region = local_cfg_region(cfg)
    region = 'CAMELS_US';
    if isstruct(cfg) ...
            && isfield(cfg,'region') ...
            && ~isempty(cfg.region)
        region = cfg.region;
    elseif isstruct(cfg) ...
            && isfield(cfg,'dirD') ...
            && ~isempty(cfg.dirD)
        region = cfg.dirD;
    end
    region = local_region_code(region);
end

function root = local_cfg_root(cfg)
    if isstruct(cfg) ...
            && isfield(cfg,'root') ...
            && ~isempty(cfg.root)
        root = char(string(cfg.root));
    elseif isstruct(cfg) ...
            && isfield(cfg,'SAGEroot') ...
            && ~isempty(cfg.SAGEroot)
        root = char(string(cfg.SAGEroot));
    else
        error('data_helpers:missingRoot', ...
            ['No writable SAGE root was provided. ' ...
             'Select a root directory on the Paths tab.']);
    end
end

function dirD = local_cfg_dirD(cfg,region)
    if nargin < 2 ...
            || isempty(region)
        region = local_cfg_region(cfg);
    end

    if isstruct(cfg) ...
            && isfield(cfg,'dirD') ...
            && ~isempty(cfg.dirD)

        dirD0 = char(string(cfg.dirD));
        [parentDir,lastFolder] = fileparts(dirD0);

        % Normal regional root, for example Data/CAMELS_GB.
        if strcmpi(lastFolder,region)
            dirD = dirD0;
            return
        end

        % Prepared GUI configurations may point to
        % Data/CAMELS_GB/daily or Data/CAMELS_GB/hourly.
        % Recover Data/CAMELS_GB before constructing stream paths.
        if any(strcmpi(lastFolder,{'daily','hourly','15min'}))
            [~,parentName] = fileparts(parentDir);
            if strcmpi(parentName,region)
                dirD = parentDir;
                return
            end
        end
    end

    root = local_cfg_root(cfg);
    dirD = fullfile(root,'Data',region);
end

function dirD = local_data_dir(root,region)
    region = local_region_code(region);
    dirD = fullfile(char( ...
        string(root)),'Data',region);
end

function d = local_stream_dir(cfg,stream)
    region = local_cfg_region(cfg);
    dirD = local_cfg_dirD(cfg,region);
    stream = lower(strtrim(char(stream)));

    switch region
        case 'CAMELS_AU'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily', ...
                    'streamflow');
            else
                d = '';
            end
        case 'CAMELS_AT'
            if strcmp(stream,'daily')
                d = fullfile(dirD, ...
                    'daily','timeseries','forcing');
            elseif strcmp(stream,'hourly')
                d = fullfile(dirD, ...
                    'hourly','timeseries','forcing');
            else
                d = '';
            end            
        case 'CAMELS_BR'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily', ...
                    'streamflow');
            else
                d = '';
            end
        case 'CAMELS_CA'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily','forcing','daymet');
            elseif strcmp(stream,'hourly')
                d = fullfile(dirD,'hourly','forcing','rdrs');
            else
                d = '';
            end
        case 'CAMELS_CL'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily', ...
                    'streamflow');
            else
                d = '';
            end
        case 'CAMELS_CH'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily', ...
                    'timeseries');
            else
                d = '';
            end
        case 'CAMELS_COL'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily', ...
                    'timeseries');
            else
                d = '';
            end
        case 'CAMELS_DK'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily','timeseries');
            else
                d = '';
            end
        case 'CAMELS_CZ'
            if any(strcmp(stream,{'daily','hourly'}))
                d = fullfile(dirD,stream,'timeseries');
            else
                d = '';
            end
        case 'CAMELS_DE'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily','timeseries');
            elseif strcmp(stream,'hourly')
                d = fullfile(dirD,'hourly','timeseries');
            else
                d = '';
            end            
        case 'CAMELS_ES'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily', ...
                    'timeseries');
            else
                d = '';
            end
        case 'BULL_ES'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily','timeseries');
            else
                d = '';
            end
        case 'CAMELS_GB'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily','timeseries');
            elseif strcmp(stream,'hourly')
                d = fullfile(dirD,'hourly','timeseries');
            else
                d = '';
            end
        case 'CAMELS_FR'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily', ...
                    'timeseries');
            else
                d = '';
            end
        case 'CAMELS_FI'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily', ...
                    'timeseries');
            else
                d = '';
            end
        case 'CAMELS_LUX'
            if any(strcmp(stream, ...
                    {'daily','hourly','15min'}))
                d = fullfile(dirD,stream,'timeseries');
            else
                d = '';
            end
        case 'CAMELS_IS'
            if any(strcmp(stream,{'daily','hourly'}))
                d = fullfile(dirD,stream,'forcing');
            else
                d = '';
            end
        case 'CAMELS_IL'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily','timeseries');
            else
                d = '';
            end
        case 'CAMELS_IND'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily', ...
                    'catchment_mean_forcings');
            else
                d = '';
            end
        case {'CAMELS_MX','CAMELS_PE','MACH_US'}
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily','timeseries');
            else
                d = '';
            end
        case 'CAMELS_PL'
            if strcmp(stream,'daily')
                d = fullfile(dirD,'daily','timeseries');
            else
                d = '';
            end
        case 'CAMELS_SE'
            if strcmp(stream,'daily')
                d = fullfile(dirD, ...
                    'daily');
            else
                d = '';
            end
        case 'CAMELS_NZ'
            if any(strcmp(stream, ...
                    {'daily','hourly'}))
                d = fullfile(dirD,stream, ...
                    'timeseries');
            else
                d = '';
            end
        case 'CAMELSH_KR'
            if strcmp(stream,'hourly')
                d = fullfile(dirD,'hourly','timeseries');
            else
                d = '';
            end
        case {'CAMELS_ZA','CAMELS_NA','CAMELS_AR','CAMELS_BE', ...
                'CAMELS_EE','CAMELS_IE','CAMELS_JM','CAMELS_NO','CAMELS_PR'}
            if strcmp(stream,'daily'), d=fullfile(dirD,'daily','timeseries');
            else, d=''; end
        case 'CAMELS_US'
            d = fullfile(dirD,stream);
        case 'CAMELSH_US'
            if strcmp(stream,'hourly'), d=fullfile(dirD,'hourly','timeseries'); else, d=''; end
        otherwise
            d = '';
    end
end

function region = local_region_code(v)
    s = strtrim(char(string(v)));
    if any(strcmpi(s,{'CA','CANADA','CAMELS_CA','CAMELS-CA'}))
        region = 'CAMELS_CA';
        return
    end
    if any(strcmpi(s,{'IS','ICELAND','CAMELS_IS','CAMELS-IS'}))
        region = 'CAMELS_IS';
        return
    end
    if any(strcmpi(s,{'DK','DENMARK', ...
            'CAMELS_DK','CAMELS-DK'}))
        region = 'CAMELS_DK';
        return
    end
    region = region_helpers('code',v);
end

function name = local_region_name(region)
    r = local_region_code(region);
    if strcmp(r,'CAMELS_US')
        name = 'United States / CAMELS-US';
        return
    end
    if strcmp(r,'CAMELSH_US')
        name = 'United States / CAMELSH-US';
        return
    end
    if strcmp(r,'CAMELS_CA')
        name = 'Canada / CAMELS-SPAT';
        return
    end
    if strcmp(r,'CAMELS_IS')
        name = 'Iceland / LamaH-Ice';
        return
    end
    if strcmp(r,'CAMELS_DK')
        name = 'Denmark / CAMELS-DK';
        return
    end
    if strcmp(r,'CAMELS_MX')
        name = 'Mexico / Caravan-HYSETS';
        return
    end
    if strcmp(r,'MACH_US')
        name = 'United States / MACH';
        return
    end
    if any(strcmp(r,{'CAMELS_ZA','CAMELS_NA','CAMELS_AR','CAMELS_BE', ...
            'CAMELS_EE','CAMELS_IE','CAMELS_JM','CAMELS_NO','CAMELS_PR'}))
        name=[region_helpers('name',r) ' / GRDC-Caravan']; return
    end
    baseName = region_helpers('name',r);
    short = region_helpers('short',r);
    name = [baseName ' / CAMELS-' short];
end

function d = local_progress_dialog(ui,title,msg)
    d = [];
    if isstruct(ui) ...
            && isfield(ui,'fig') ...
            && ~isempty(ui.fig) ...
            && isvalid(ui.fig)
        try
            d = uiprogressdlg(ui.fig, ...
                'Title',title, ...
                'Message',msg, ...
                'Indeterminate',false, ...
                'Value',0, ...
                'Cancelable',true);
        catch
            d = [];
        end
    end
end

function local_close_progress(d)
    try
        if ~isempty(d) ...
                && isvalid(d)
            close(d);
        end
    catch
    end
end

function [userCanceled,lastUI] = ...
    local_update_progress_dialog( ...
    d,info,lastUI,userCanceled)

    if isempty(d) ...
            || ~isvalid(d)
        return
    end
    
    if isfield(info,'indeterminate') ...
            && logical(info.indeterminate)
        d.Indeterminate = true;
        d.Message = sprintf('%s: downloading ...', ...
            local_info_label(info));
        drawnow limitrate nocallbacks
        return
    end

    if toc(lastUI) < 0.25
        return
    end
    lastUI = tic;
    
    if d.CancelRequested
        userCanceled = true;
        d.Message = ['Cancel requested ... ' ...
            'stopping download'];
        drawnow limitrate nocallbacks
        return
    end
    
    if isfield(info,'frac') ...
            && isfinite(info.frac)
        d.Indeterminate = false;
        d.Value = info.frac;
        speed = local_info_speed(info);
        if isfield(info,'etaSec') ...
                && isfinite(info.etaSec)
            if isfinite(speed)
                d.Message = sprintf(['Downloaded: %5.1f%%  ' ...
                    '(%.2f MB/s, Time remaining %s)'], ...
                    100*info.frac, ...
                    speed,local_fmt_eta(info.etaSec));
            else
                d.Message = sprintf('Downloaded: %5.1f%%', ...
                    100*info.frac);
            end
        else
            if isfinite(speed)
                d.Message = sprintf(['Downloaded: %5.1f%%  ' ...
                    '(%.2f MB/s)'], ...
                    100*info.frac,speed);
            else
                d.Message = sprintf('Downloaded: %5.1f%%', ...
                    100*info.frac);
            end
        end
    else
        speed = local_info_speed(info);
        if isfinite(speed)
            d.Message = sprintf(['%s: downloading ... ' ...
                '(%.2f MB/s)'], ...
                local_info_label(info),speed);
        else
            d.Message = sprintf('%s: downloading ...', ...
                local_info_label(info));
        end
    end
    
    drawnow limitrate nocallbacks
end

function logFcn = local_ui_log(ui)
    if isstruct(ui) ...
            && isfield(ui,'logFcn') ...
            && isa(ui.logFcn,'function_handle')
        logFcn = ui.logFcn;
    else
        logFcn = @(s) fprintf('%s\n', ...
            char(string(s)));
    end
end

function tf = local_ui_confirm(ui,title,msg)
    if isstruct(ui) ...
            && isfield(ui,'confirmFcn') ...
            && isa(ui.confirmFcn,'function_handle')
        tf = ui.confirmFcn(title,msg);
        return
    end
    
    if isstruct(ui) ...
            && isfield(ui,'fig') ...
            && ~isempty(ui.fig) ...
            && isvalid(ui.fig)
        try
            choice = uiconfirm(ui.fig,msg,title, ...
                'Options',{'Continue','Cancel'}, ...
                'DefaultOption',2, ...
                'CancelOption',2);
            tf = strcmp(choice,'Continue');
            return
        catch
        end
    end
    
    tf = true;
end


% ===========================================
function files = local_gb_v2_attribute_files()
% ===========================================
    files = { ...
        'camels_gb_v2_climatic_attributes.csv', ...
        'camels_gb_v2_groundwaterwell_attributes.csv', ...
        'camels_gb_v2_humaninfluence_attributes.csv', ...
        'camels_gb_v2_hydrogeology_attributes.csv', ...
        'camels_gb_v2_hydrologic_attributes.csv', ...
        'camels_gb_v2_hydrometry_attributes.csv', ...
        'camels_gb_v2_landcover_attributes.csv', ...
        'camels_gb_v2_soil_attributes.csv', ...
        'camels_gb_v2_topographic_attributes.csv'};
end

% =====================================================
function files = local_gb_remote_csv_list(indexUrl,prefix)
% =====================================================
% Read an Apache-style directory index and return matching CSV hrefs.
    try
        html = webread(indexUrl);
    catch ME
        error('data_helpers:GBDirectoryIndexFailed', ...
            'Could not read CAMELS-GB directory index %s: %s', ...
            indexUrl,ME.message);
    end
    tok = regexp(char(html), ...
        'href\s*=\s*["'']([^"'']+\.csv)["'']', ...
        'tokens','ignorecase');
    files = cellfun(@(c) c{1},tok,'UniformOutput',false);
    files = cellfun(@local_gb_url_decode,files,'UniformOutput',false);
    files = files(startsWith(lower(string(files)),lower(string(prefix))));
    files = unique(cellstr(string(files)),'stable');
end

% ============================================
function s = local_gb_url_decode(s)
% ============================================
    s = char(string(s));
    s = strrep(s,'%20',' ');
end

% ==========================================================
function iDone = local_gb_download_one(url,dst,label, ...
    iDone,nTotal,d,cancelFcn,progressFcn,logFcn)
% ==========================================================
% Download one direct CSV, resuming an installation by skipping nonempty
% files. A partial/zero-byte file is deleted and downloaded again.
    if isfile(dst)
        info = dir(dst);
        name = lower(string(dst));
        minBytes = 100;
        if contains(name,'hydromet_daily_timeseries_')
            minBytes = 1e5;
        elseif contains(name,'hydromet_hourly_timeseries_')
            minBytes = 1e6;
        end
        if ~isempty(info) && info.bytes >= minBytes
            iDone = iDone + 1;
            progressFcn(struct('frac',iDone/nTotal, ...
                'label',[label ' (already present)']));
            return
        end
        try 
            delete(dst); 
        catch
        end
    end

    local_mkdir(fileparts(dst));
    local_download_file_retry(url,dst,label,d, ...
        cancelFcn, ...
        @file_progress,logFcn,3);
    iDone = iDone + 1;
    progressFcn(struct('frac',iDone/nTotal,'label',label));

    function file_progress(info)
        if isfield(info,'frac') && isfinite(info.frac)
            info.frac = (iDone + info.frac)/nTotal;
        else
            info.frac = iDone/nTotal;
        end
        info.label = label;
        progressFcn(info);
    end
end

% =====================
% File/status utilities
% =====================
function txt = local_present_text(tf)
    if tf
        txt = 'Present';
    else
        txt = 'Missing';
    end
end

function ok = local_all_files_exist(d,files)
    ok = isfolder(d);
    if ~ok
        return
    end
    for i = 1:numel(files)
        if ~isfile(fullfile(d,files{i}))
            ok = false;
            return
        end
    end
end

function n = local_count_files(d,patterns)
    n = 0;
    if ~isfolder(d)
        return
    end
    for i = 1:numel(patterns)
        q = dir(fullfile(d,patterns{i}));
        n = n + sum(~[q.isdir]);
    end
end

function s = local_fmt_eta(sec)
    sec = max(sec,0);
    mm = floor(sec/60);
    ss = floor(sec - 60*mm);
    if mm < 60
        s = sprintf('%02d:%02d',mm,ss);
    else
        hh = floor(mm/60);
        mm = mm - 60*hh;
        s = sprintf('%02d:%02d:%02d',hh,mm,ss);
    end
end

function s = local_info_label(info)
    if isfield(info,'label') ...
            && ~isempty(info.label)
        s = char(string(info.label));
    else
        s = 'Data';
    end
end

function speed = local_info_speed(info)
    if isfield(info,'speedMBs') ...
            && isfinite(info.speedMBs)
        speed = info.speedMBs;
    else
        speed = NaN;
    end
end

function s = local_cap(s)
    s = char(string(s));
    if isempty(s)
        return
    end
    s(1) = upper(s(1));
end


function files = local_camels_col_files()
    files = struct( ...
        'key',{},'label',{},'url',{}, ...
        'fileName',{},'isZip',{},'md5',{});

    files(end+1) = struct( ...
        'key','attributes_zip', ...
        'label','CAMELS-COL attributes workbook', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '01_CAMELS_COL_Attributes.zip?download=1'], ...
        'fileName','01_CAMELS_COL_Attributes.zip', ...
        'isZip',true, ...
        'md5','21945e75d1b0e4b04a854da79c8be218');

    files(end+1) = struct( ...
        'key','catchment_information', ...
        'label','CAMELS-COL catchment information', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '02_CAMELS_COL_Catchment_information.csv?download=1'], ...
        'fileName','02_CAMELS_COL_Catchment_information.csv', ...
        'isZip',false, ...
        'md5','0d50edc9dc1fb59ae5e7c045fbddeca9');

    files(end+1) = struct( ...
        'key','hydromet_zip', ...
        'label','CAMELS-COL hydrometeorological data', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '04_CAMELS_COL_Hydrometeorological_data.zip?download=1'], ...
        'fileName','04_CAMELS_COL_Hydrometeorological_data.zip', ...
        'isZip',true, ...
        'md5','2101d7cac0f87b2a8d871f7a622bf440');

    files(end+1) = struct( ...
        'key','geologic', ...
        'label','CAMELS-COL geologic characteristics', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '05_CAMELS_COL_Geologic_characteristics.csv?download=1'], ...
        'fileName','05_CAMELS_COL_Geologic_characteristics.csv', ...
        'isZip',false, ...
        'md5','ee1d3b3974dcf4ff82b5289b2fbe6b37');

    files(end+1) = struct( ...
        'key','land_cover', ...
        'label','CAMELS-COL land-cover characteristics', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '06_CAMELS_COL_Land_cover_characteristics.csv?download=1'], ...
        'fileName','06_CAMELS_COL_Land_cover_characteristics.csv', ...
        'isZip',false, ...
        'md5','e81e5583ff0bb67754feee5a0a157abd');

    files(end+1) = struct( ...
        'key','soil', ...
        'label','CAMELS-COL soil characteristics', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '07_CAMELS_COL_Soil_characteristics.csv?download=1'], ...
        'fileName','07_CAMELS_COL_Soil_characteristics.csv', ...
        'isZip',false, ...
        'md5','2c0a5a68f722b1dd8ecd87c632f0cbdd');

    files(end+1) = struct( ...
        'key','climatic_indices', ...
        'label','CAMELS-COL climatic indices', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '08_CAMELS_COL_Climatic_indices.csv?download=1'], ...
        'fileName','08_CAMELS_COL_Climatic_indices.csv', ...
        'isZip',false, ...
        'md5','83a9037a08612e905db96ebb16596196');

    files(end+1) = struct( ...
        'key','hydrological_signatures', ...
        'label','CAMELS-COL hydrological signatures', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '09_CAMELS_COL_Hydrological_signatures.csv?download=1'], ...
        'fileName','09_CAMELS_COL_Hydrological_signatures.csv', ...
        'isZip',false, ...
        'md5','ef22677856c2b103af6bffaca981c407');

    files(end+1) = struct( ...
        'key','physiographic', ...
        'label','CAMELS-COL physiographic characteristics', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '10_CAMELS_COL_Physiograpic_characteristics.csv?download=1'], ...
        'fileName','10_CAMELS_COL_Physiograpic_characteristics.csv', ...
        'isZip',false, ...
        'md5','5193a248d7ea2675675ed249d2b692c2');

    files(end+1) = struct( ...
        'key','land_use_capability', ...
        'label','CAMELS-COL land-use capability', ...
        'url',['https://zenodo.org/records/18794895/files/' ...
        '11_CAMELS_COL_Land_use_capability.csv?download=1'], ...
        'fileName','11_CAMELS_COL_Land_use_capability.csv', ...
        'isZip',false, ...
        'md5','bcc218c38b605e9068ed47a31920d7fb');
end

function files = local_col_metadata_files()
    files = { ...
        '01_CAMELS_COL_Attributes.xlsx', ...
        '02_CAMELS_COL_Catchment_information.csv', ...
        '05_CAMELS_COL_Geologic_characteristics.csv', ...
        '06_CAMELS_COL_Land_cover_characteristics.csv', ...
        '07_CAMELS_COL_Soil_characteristics.csv', ...
        '08_CAMELS_COL_Climatic_indices.csv', ...
        '09_CAMELS_COL_Hydrological_signatures.csv', ...
        '10_CAMELS_COL_Physiograpic_characteristics.csv', ...
        '11_CAMELS_COL_Land_use_capability.csv'};
end

function files = local_es_metadata_files()

    files = { ...
        'atributes_efas_hydrometeorology_camelses.csv'
        'atributes_efas_static_maps_camelses.csv'
        'attributes_caravan_camelses.csv'
        'attributes_efas_model_parameters_camelses.csv'
        'attributes_hydroatlas_camelses.csv'
        'attributes_other_camelses.csv'};
end

function local_install_es_files(unzipRoot,dirD,logFcn)
%LOCAL_INSTALL_ES_FILES Install CAMELS-ES files from unzipped archive.

    if nargin < 3 ...
            || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n', ...
            char(string(s)));
    end
    
    if isempty(unzipRoot) ...
            || ~isfolder(unzipRoot)
        error('data_helpers:ESBadUnzipRoot', ...
            ['CAMELS-ES unzip root ' ...
            'does not exist: %s'], ...
            unzipRoot);
    end
    
    attrFiles = local_es_metadata_files();
    
    for i = 1:numel(attrFiles)
        src = local_find_file_recursive_manual( ...
            unzipRoot,attrFiles{i});
    
        if isempty(src) ...
                || ~isfile(src)
            error('data_helpers:ESMissingAttributeFile', ...
                ['Could not find required ' ...
                'CAMELS-ES file: %s'], ...
                attrFiles{i});
        end
    
        copyfile(src,fullfile(dirD, ...
            attrFiles{i}),'f');
    
        logFcn(['Copied ' attrFiles{i} ...
            ' to: ' dirD]);
    end
    
    tsProbe = local_find_file_recursive_manual( ...
        unzipRoot,'camelses_1080.csv');
    
    if isempty(tsProbe)
        allTs = dir(fullfile(unzipRoot, ...
            '**','camelses_*.csv'));
        if ~isempty(allTs)
            tsProbe = fullfile(allTs(1).folder, ...
                allTs(1).name);
        end
    end
    
    if isempty(tsProbe) ...
            || ~isfile(tsProbe)
        error('data_helpers:ESMissingTimeseries', ...
            ['Could not find ' ...
            'CAMELS-ES daily ' ...
            'time-series CSV files ' ...
            'inside the archive.']);
    end
    
    tsSrc = fileparts(tsProbe);
    tsDst = fullfile(dirD,'daily','timeseries');
    
    local_copy_folder_contents(tsSrc, ...
        tsDst,logFcn);
end

function tf = local_es_install_complete(dirD)

    tf = local_all_files_exist(dirD, ...
        local_es_metadata_files()) ...
        && isfolder(fullfile(dirD, ...
        'daily','timeseries')) ...
        && local_count_files(fullfile(dirD, ...
        'daily','timeseries'), ...
        {'camelses_*.csv'}) > 0;
end

function tf = local_has_col_timeseries(dirD)
    d = fullfile(dirD,'daily','timeseries');
    tf = isfolder(d) ...
        && local_count_files(d,{'Hydromet_data_*.txt'}) > 0;
end

function tf = local_col_install_complete(dirD)
    tf = local_all_files_exist(dirD, ...
        local_col_metadata_files()) ...
        && local_has_col_timeseries(dirD);
end

function files = local_pl_metadata_files()
    files = { ...
        'CAMELS_PL_BDOT10K_land_cover_catchments.csv', ...
        'CAMELS_PL_Q_354_data.csv', ...
        'CAMELS_PL_climatic_attributes.csv', ...
        'CAMELS_PL_hydrologic_attributes.csv', ...
        'CAMELS_PL_landcover_attributes.csv', ...
        'CAMELS_PL_simulation_benchmark.csv', ...
        'CAMELS_PL_soil_attributes.csv', ...
        'CAMELS_PL_topographic_attributes.csv'};
end

function tf = local_pl_install_complete(dirD)
    tf = local_all_files_exist(dirD, ...
        local_pl_metadata_files()) ...
        && isfolder(fullfile(dirD,'daily','timeseries')) ...
        && local_count_files(fullfile(dirD,'daily','timeseries'), ...
        {'CAMELS_PL_hydromet_timeseries_*.csv'}) >= 354;
end

function local_install_pl_files(srcRoot,dirD,logFcn)

    dataDir = local_find_child_dir(srcRoot,'CAMELS-PL');
    if isempty(dataDir)
        dataDir = local_find_child_dir(srcRoot,'CAMELS_PL');
    end
    if isempty(dataDir)
        dataDir = srcRoot;
    end
    
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily','timeseries'));
    
    fAttr = local_pl_metadata_files();
    for i = 1:numel(fAttr)
        f0 = fullfile(dataDir,fAttr{i});
        if isfile(f0)
            copyfile(f0,fullfile(dirD,fAttr{i}),'f');
        else
            logFcn(['Missing CAMELS-PL file: ' fAttr{i}]);
        end
    end
    
    srcTS = fullfile(dataDir,'timeseries');
    if ~isfolder(srcTS)
        error('data_helpers:PLMissingTimeseries', ...
            'Could not find CAMELS-PL timeseries folder.');
    end
    
    D = dir(fullfile(srcTS, ...
        'CAMELS_PL_hydromet_timeseries_*.csv'));
    
    for i = 1:numel(D)
        copyfile(fullfile(D(i).folder,D(i).name), ...
            fullfile(dirD,'daily','timeseries',D(i).name),'f');
    end
    
    logFcn(sprintf(['Copied %d CAMELS-PL ' ...
        'hydrometeorological ' ...
        'time-series files.'],numel(D)));

end

function out = local_download_pe(cfg,stream,ui)
%LOCAL_DOWNLOAD_PE Download and install CAMELS-PE v1.0.1.
    if ~strcmpi(stream,'daily')
        error('data_helpers:PEdailyOnly','CAMELS-PE supports daily data only.');
    end
    dirD = local_cfg_dirD(cfg,'CAMELS_PE');
    dirTS = fullfile(dirD,'daily','timeseries');
    local_mkdir(dirD); local_mkdir(fullfile(dirD,'daily')); local_mkdir(dirTS);
    logFcn = local_ui_log(ui);
    if ~local_ui_confirm(ui,'Download Peru data', ...
            ['Download CAMELS-PE v1.0.1 from Zenodo? The archive contains ' ...
             '136 daily catchment time series, seven attribute files, and two metadata files.'])
        out = struct('ok',false,'canceled',true); return
    end
    downloadDir = local_default_download_dir();
    zipFile = fullfile(downloadDir,'CAMELS-PE_v1.0.1.zip');
    unzipDir = fullfile(downloadDir,'CAMELS_PE_download');
    url = ['https://zenodo.org/records/21195425/files/' ...
           'CAMELS-PE_v1.0.1.zip?download=1'];
    expectedMD5 = '13f127d381338eee0e35359c08dba199';
    userCanceled = false; lastUI = tic;
    d = local_progress_dialog(ui,'Downloading Peru data', ...
        'Starting CAMELS-PE download ...');
    try
        if local_pe_install_complete(dirD)
            logFcn('CAMELS-PE daily data appear complete; skipping download.');
            out = struct('ok',true,'region','CAMELS_PE','dirD',dirD, ...
                'dirM',dirTS,'dirQ',dirTS,'archive','','unzipDir','');
            local_close_progress(d); return
        end
        if isfile(zipFile)
            try
                local_verify_md5(zipFile,expectedMD5,logFcn);
            catch
                delete(zipFile);
            end
        end
        if ~isfile(zipFile)
            local_download_file_retry(url,zipFile,'CAMELS-PE',d, ...
                @() userCanceled,@(info) ui_progress(info),logFcn,3);
        end
        if userCanceled
            out = struct('ok',false,'canceled',true); local_close_progress(d); return
        end
        local_verify_md5(zipFile,expectedMD5,logFcn);
        if isfolder(unzipDir), rmdir(unzipDir,'s'); end
        local_mkdir(unzipDir); unzip(zipFile,unzipDir);
        payload = local_find_pe_payload(unzipDir);
        if isempty(payload)
            error('data_helpers:PEMissingPayload', ...
                'Could not locate CAMELS-PE/01_metadata through 03_timeseries.');
        end
        metaFiles = {'stations.csv','data_dictionary.csv'};
        attrFiles = {'climatic_indices.csv','geologic_attributes.csv', ...
            'human_intervention_attributes.csv','hydrological_signatures.csv', ...
            'landcover_attributes.csv','soil_attributes.csv', ...
            'topographic_attributes.csv'};
        for k=1:numel(metaFiles)
            copyfile(fullfile(payload,'01_metadata',metaFiles{k}), ...
                fullfile(dirD,metaFiles{k}),'f');
        end
        for k=1:numel(attrFiles)
            copyfile(fullfile(payload,'02_attributes',attrFiles{k}), ...
                fullfile(dirD,attrFiles{k}),'f');
        end
        srcTS = fullfile(payload,'03_timeseries','by_catchment');
        D = dir(fullfile(srcTS,'PE_*.csv'));
        if numel(D) < 136
            error('data_helpers:PEMissingTimeseries', ...
                'Expected at least 136 CAMELS-PE files; found %d.',numel(D));
        end
        for k=1:numel(D)
            copyfile(fullfile(D(k).folder,D(k).name), ...
                fullfile(dirTS,D(k).name),'f');
            if ~isempty(d) && isvalid(d) && (mod(k,5)==0 || k==numel(D))
                d.Message=sprintf('Installing CAMELS-PE files %d/%d ...',k,numel(D));
                drawnow limitrate nocallbacks
            end
        end
        ok = local_pe_install_complete(dirD);
        out = struct('ok',ok,'region','CAMELS_PE','dirD',dirD, ...
            'dirM',dirTS,'dirQ',dirTS,'archive',zipFile,'unzipDir',unzipDir);
        if ok
            logFcn(sprintf('CAMELS-PE install finished: %d daily files.',numel(D)));
            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive=''; out.unzipDir='';
        end
    catch ME
        logFcn('CAMELS-PE install failed:'); logFcn(ME.message);
        out=struct('ok',false,'error',ME.message);
    end
    local_close_progress(d);
    function ui_progress(info)
        [userCanceled,lastUI] = local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

function tf = local_pe_install_complete(dirD)
    tf = local_all_files_exist(dirD,{ ...
        'stations.csv','data_dictionary.csv','climatic_indices.csv', ...
        'geologic_attributes.csv','human_intervention_attributes.csv', ...
        'hydrological_signatures.csv','landcover_attributes.csv', ...
        'soil_attributes.csv','topographic_attributes.csv'}) ...
        && local_count_files(fullfile(dirD,'daily','timeseries'), ...
            {'PE_*.csv'}) >= 136;
end

function payload = local_find_pe_payload(root)
    payload = '';
    candidates = {fullfile(root,'CAMELS-PE_v1.0.1','CAMELS-PE'), ...
        fullfile(root,'CAMELS-PE'),root};
    for k=1:numel(candidates)
        c=candidates{k};
        if isfolder(fullfile(c,'01_metadata')) ...
                && isfolder(fullfile(c,'02_attributes')) ...
                && isfolder(fullfile(c,'03_timeseries','by_catchment'))
            payload=c; return
        end
    end
    D=dir(fullfile(root,'**','03_timeseries','by_catchment'));
    if ~isempty(D)
        c=fileparts(fileparts(fullfile(D(1).folder,D(1).name)));
        if isfolder(fullfile(c,'01_metadata')) && isfolder(fullfile(c,'02_attributes'))
            payload=c;
        end
    end
end

function out = local_download_mach(cfg,stream,ui)

    logFcn = local_ui_log(ui);
    if ~strcmp(stream,'daily')
        error('data_helpers:MACHOnlyDaily', ...
            'The MACH dataset is daily only.');
    end

    dirD = local_cfg_dirD(cfg,'MACH_US');
    dirTS = fullfile(dirD,stream,'timeseries');
    dirAttr = fullfile(dirD,'attributes');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,stream));
    local_mkdir(dirTS);
    local_mkdir(dirAttr);

    if ~local_ui_confirm(ui,'Download United States (MACH)', ...
            ['MACH version 4.0 provides 1,014 combined daily forcing ' ...
            'and streamflow files for 1980--2023, catchment ' ...
            'attributes, and basin boundaries. Approximately 520 MB ' ...
            'will be downloaded from Zenodo. Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    stageDir = fullfile(downloadDir,'MACH_US_download');
    files = { ...
        'MACH_ts.zip','ad9f88795b383278af00da84ec4453a1'; ...
        'attributes.zip','f5a9c8d12c153690a25b6f210a47a0e2'; ...
        'MACH_basins_all.gpkg','4733b6fa3efab7cab77ae110a49901c1'; ...
        'README.csv','c74872f3e235a591526020d36ef3bac2'};
    d = local_progress_dialog(ui,'Downloading MACH-US data', ...
        'Starting MACH version 4.0 download ...');
    userCanceled = false;
    lastUI = tic;

    try
        if local_mach_install_complete(dirD)
            logFcn(['MACH version 4.0 appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true,'region','MACH_US','dirD',dirD, ...
                'dirM',dirTS,'dirQ',dirTS,'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        for i = 1:size(files,1)
            name = files{i,1};
            target = fullfile(downloadDir,name);
            url = ['https://zenodo.org/records/18686475/files/' ...
                name '?download=1'];
            if isfile(target)
                try
                    local_verify_md5(target,files{i,2},logFcn);
                catch
                    delete(target);
                end
            end
            if ~isfile(target)
                local_download_file_retry(url,target, ...
                    ['MACH ' name],d,@() userCanceled, ...
                    @ui_progress,logFcn,3);
            end
            local_verify_md5(target,files{i,2},logFcn);
        end

        if isfolder(stageDir)
            rmdir(stageDir,'s');
        end
        local_mkdir(stageDir);
        try
            d.Indeterminate = true;
            d.Message = 'Extracting MACH time series and attributes ...';
            drawnow limitrate nocallbacks
        catch
        end
        unzip(fullfile(downloadDir,'MACH_ts.zip'),stageDir);
        unzip(fullfile(downloadDir,'attributes.zip'), ...
            fullfile(stageDir,'attributes'));

        ts = dir(fullfile(stageDir,'**','basin_*_MACH.csv'));
        if numel(ts) ~= 1014
            error('data_helpers:MACHBadArchive', ...
                'Expected 1,014 MACH time-series files, but found %d.', ...
                numel(ts));
        end
        for i = 1:numel(ts)
            copyfile(fullfile(ts(i).folder,ts(i).name), ...
                fullfile(dirTS,ts(i).name),'f');
        end

        attr = dir(fullfile(stageDir,'attributes','*.csv'));
        for i = 1:numel(attr)
            copyfile(fullfile(attr(i).folder,attr(i).name), ...
                fullfile(dirAttr,attr(i).name),'f');
        end
        copyfile(fullfile(downloadDir,'MACH_basins_all.gpkg'), ...
            fullfile(dirD,'MACH_basins_all.gpkg'),'f');
        copyfile(fullfile(downloadDir,'README.csv'), ...
            fullfile(dirD,'README.csv'),'f');
        local_write_mach_basin_universe(dirD);

        out = struct('ok',local_mach_install_complete(dirD), ...
            'region','MACH_US','dirD',dirD,'dirM',dirTS, ...
            'dirQ',dirTS,'archive','','unzipDir',stageDir);
        if out.ok
            logFcn(['MACH version 4.0 install finished: 1,014 ' ...
                'combined daily time-series files and attributes.']);
            for i = 1:size(files,1)
                local_cleanup_download_file( ...
                    fullfile(downloadDir,files{i,1}),logFcn);
            end
            local_cleanup_download_folder(stageDir,logFcn);
            out.unzipDir = '';
        end
    catch ME
        logFcn('MACH version 4.0 install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end
    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

function tf = local_mach_install_complete(dirD)
    tf = local_count_files(fullfile(dirD,'daily','timeseries'), ...
        {'basin_*_MACH.csv'}) >= 1014 ...
        && local_all_files_exist(dirD,{ ...
        fullfile('attributes','site_info.csv'), ...
        fullfile('attributes','overall_climate.csv'), ...
        fullfile('attributes','soil.csv'), ...
        fullfile('attributes','geology.csv'), ...
        fullfile('attributes','hydrology.csv'), ...
        fullfile('attributes','anthropogenic.csv'), ...
        'MACH_basins_all.gpkg','README.csv','MACH_1014_basins.txt'});
end

function local_write_mach_basin_universe(dirD)
    source = fullfile(dirD,'attributes','site_info.csv');
    opts = detectImportOptions(source,'VariableNamingRule','preserve');
    opts = setvartype(opts,'SITENO','string');
    T = readtable(source,opts);
    ids = strip(string(T.SITENO));
    ids = compose('%08d',str2double(ids));
    ids = sort(ids);
    if numel(unique(ids)) ~= 1014
        error('data_helpers:MACHBasinUniverse', ...
            'Expected 1,014 unique MACH gauge identifiers.');
    end
    writelines(ids,fullfile(dirD,'MACH_1014_basins.txt'));
end

function out = local_download_mx(cfg,stream,ui)

    logFcn = local_ui_log(ui);
    if ~strcmp(stream,'daily')
        error('data_helpers:MXOnlyDaily', ...
            'The Caravan-HYSETS Mexico subset is daily only.');
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_MX');
    dirTS = fullfile(dirD,stream,'timeseries');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,stream));
    local_mkdir(dirTS);

    if ~local_ui_confirm(ui,'Download Mexico data', ...
            ['Mexico is extracted from the official Caravan-HYSETS ' ...
            'CSV archive. Zenodo does not provide country-level files, ' ...
            'so the installer must download the complete 29 GB archive ' ...
            'but retains only 46 Mexican basins. The archive is deleted ' ...
            'after successful installation. Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    archiveFile = fullfile(downloadDir,'Caravan-csv.tar.gz');
    url = ['https://zenodo.org/api/records/15530022/files/' ...
        'Caravan-csv.tar.gz/content'];
    expectedMD5 = 'ede5b507bd7e7ab9558d8ab1cedf65b1';
    d = local_progress_dialog(ui,'Downloading Mexico data', ...
        'Starting the Caravan CSV download ...');
    userCanceled = false;
    lastUI = tic;

    try
        if local_count_files(dirTS,{'hysets_*.csv'}) >= 46 ...
                && local_all_files_exist(dirD,{ ...
                'attributes_caravan_hysets.csv', ...
                'attributes_hydroatlas_hysets.csv', ...
                'attributes_other_hysets.csv'})
            logFcn(['Caravan-HYSETS Mexico appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true,'region','CAMELS_MX', ...
                'dirD',dirD,'dirM',dirTS,'dirQ',dirTS, ...
                'archive','','unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(archiveFile)
            try
                local_verify_md5(archiveFile,expectedMD5,logFcn);
            catch
                delete(archiveFile);
            end
        end
        if ~isfile(archiveFile)
            local_download_file_retry(url,archiveFile, ...
                'Caravan CSV archive',d,@() userCanceled, ...
                @ui_progress,logFcn,3);
        end
        local_verify_md5(archiveFile,expectedMD5,logFcn);

        stageDir = fullfile(downloadDir,'CAMELS_MX_download');
        if isfolder(stageDir)
            rmdir(stageDir,'s');
        end
        local_mkdir(stageDir);
        try
            d.Indeterminate = true;
            d.Message = ['Extracting 46 Mexican HYSETS basins from ' ...
                'the Caravan archive ...'];
            drawnow limitrate nocallbacks
        catch
        end

        ids = local_mx_source_ids();
        members = { ...
            'Caravan-csv/attributes/hysets/attributes_other_hysets.csv', ...
            'Caravan-csv/attributes/hysets/attributes_caravan_hysets.csv', ...
            'Caravan-csv/attributes/hysets/attributes_hydroatlas_hysets.csv'};
        for i = 1:numel(ids)
            members{end+1} = ['Caravan-csv/timeseries/csv/hysets/' ...
                'hysets_' ids{i} '.csv']; %#ok<AGROW>
        end
        local_extract_tar_members(archiveFile,stageDir,members);

        sourceTS = fullfile(stageDir,'Caravan-csv','timeseries', ...
            'csv','hysets');
        for i = 1:numel(ids)
            name = ['hysets_' ids{i} '.csv'];
            copyfile(fullfile(sourceTS,name),fullfile(dirTS,name),'f');
        end

        sourceAttr = fullfile(stageDir,'Caravan-csv','attributes','hysets');
        attrNames = {'attributes_other_hysets.csv', ...
            'attributes_caravan_hysets.csv', ...
            'attributes_hydroatlas_hysets.csv'};
        sourceGauge = "hysets_" + string(ids(:));
        for i = 1:numel(attrNames)
            T = readtable(fullfile(sourceAttr,attrNames{i}), ...
                'VariableNamingRule','preserve');
            key = string(T{:,1});
            T = T(ismember(key,sourceGauge),:);
            if height(T) ~= 46
                error('data_helpers:MXAttributeSubset', ...
                    'Expected 46 Mexican rows in %s, but found %d.', ...
                    attrNames{i},height(T));
            end
            writetable(T,fullfile(dirD,attrNames{i}));
        end

        out = struct('ok', ...
            local_count_files(dirTS,{'hysets_*.csv'}) >= 46, ...
            'region','CAMELS_MX','dirD',dirD, ...
            'dirM',dirTS,'dirQ',dirTS, ...
            'archive',archiveFile,'unzipDir',stageDir);
        if out.ok
            logFcn(['Caravan-HYSETS Mexico install finished: ' ...
                '46 combined daily time-series files.']);
            local_cleanup_download_file(archiveFile,logFcn);
            local_cleanup_download_folder(stageDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        end
    catch ME
        logFcn('Caravan-HYSETS Mexico install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

function ids = local_mx_source_ids()

    ids = {'10034','10037','10063','11027','11040','12238', ...
        '12370','12371','12405','12436','12469','12484','12485', ...
        '12487','12488','12504','12520','12535','18277','18344', ...
        '18350','18538','18555','20017','20021','20022','20026', ...
        '20037','22008','22015','22016','26243','26249','26280', ...
        '26285','27001','27038','28001','28013','28015','28108', ...
        '28125','28136','30158','36071','9011'};
end

function local_extract_tar_members(archiveFile,targetDir,members)

    quotedMembers = cellfun(@local_shell_quote,members, ...
        'UniformOutput',false);
    command = ['tar -xzf ' local_shell_quote(archiveFile) ...
        ' -C ' local_shell_quote(targetDir) ' ' ...
        strjoin(quotedMembers,' ')];
    [status,message] = system(command);
    if status ~= 0
        error('data_helpers:TarExtractionFailed', ...
            'Selective Caravan extraction failed: %s',message);
    end
end

function value = local_shell_quote(value)

    value = char(string(value));
    if ispc
        value = ['"' strrep(value,'"','""') '"'];
    else
        value = ['''' strrep(value,'''','''"''"''') ''''];
    end
end

function out = local_download_jp(cfg,stream,ui)

    logFcn = local_ui_log(ui);

    if ~strcmp(stream,'daily')
        error('data_helpers:JPOnlyDaily', ...
            'MERV-Jp 2.0 support is daily only.');
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_JP');
    dirTS = fullfile(dirD,'daily','timeseries');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(dirTS);

    if ~local_ui_confirm(ui, ...
            'Download Japan data', ...
            ['MERV-Jp version 2.0 is downloaded from Zenodo. ' ...
            'The 579 MB archive contains combined daily forcing and ' ...
            'observed runoff for 87 Japanese basins. Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipFile = fullfile(downloadDir,'varssim.zip');
    url = ['https://zenodo.org/records/8176305/files/' ...
        'varssim.zip?download=1'];
    expectedMD5 = '8e538760fca9223a79d9a7c8da09f6ec';
    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Japan data', ...
        'Starting MERV-Jp 2.0 download ...');

    try
        if local_count_files(dirTS,{'varssim*.csv'}) >= 87
            local_install_jp_attributes(cfg,dirD,logFcn);
            logFcn(['MERV-Jp 2.0 appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true,'region','CAMELS_JP', ...
                'dirD',dirD,'dirM',dirTS,'dirQ',dirTS, ...
                'archive','','unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(zipFile)
            try
                local_verify_md5(zipFile,expectedMD5,logFcn);
                logFcn(['Using verified existing file: ' zipFile]);
            catch
                logFcn(['Existing MERV-Jp archive failed MD5; ' ...
                    'deleting it before downloading again.']);
                delete(zipFile);
            end
        end

        if ~isfile(zipFile)
            local_download_file_retry(url,zipFile, ...
                'MERV-Jp 2.0 archive',d,@() userCanceled, ...
                @ui_progress,logFcn,3);
        end
        local_verify_md5(zipFile,expectedMD5,logFcn);

        unzipDir = fullfile(downloadDir,'MERV_JP_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);
        try
            d.Indeterminate = true;
            d.Message = 'Extracting MERV-Jp version 2.0 ...';
            drawnow limitrate nocallbacks
        catch
        end
        unzip(zipFile,unzipDir);

        files = dir(fullfile(unzipDir,'**','ver2_0','varssim*.csv'));
        if numel(files) ~= 87
            error('data_helpers:JPBadArchive', ...
                ['Expected 87 version-2 files in varssim/ver2_0, ' ...
                'but found %d.'],numel(files));
        end
        for i = 1:numel(files)
            copyfile(fullfile(files(i).folder,files(i).name), ...
                fullfile(dirTS,files(i).name),'f');
        end
        local_install_jp_attributes(cfg,dirD,logFcn);

        out = struct('ok', ...
            local_count_files(dirTS,{'varssim*.csv'}) >= 87, ...
            'region','CAMELS_JP','dirD',dirD, ...
            'dirM',dirTS,'dirQ',dirTS, ...
            'archive',zipFile,'unzipDir',unzipDir);
        if out.ok
            logFcn(['MERV-Jp 2.0 install finished: ' ...
                '87 combined daily time-series files.']);
            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);
            out.archive = '';
            out.unzipDir = '';
        end
    catch ME
        logFcn('MERV-Jp 2.0 install failed:');
        logFcn(ME.message);
        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end
end

function local_install_jp_attributes(cfg,dirD,logFcn)
    name = 'MERV_Jp_135_HydroATLAS_attributes.xlsx';
    target = fullfile(dirD,name);
    if isfile(target)
        return
    end

    candidates = {};
    if isstruct(cfg) ...
            && isfield(cfg,'SAGEhydro') ...
            && ~isempty(cfg.SAGEhydro)
        candidates{end+1} = fullfile( ...
            char(cfg.SAGEhydro),'regions','JP',name);
    end
    sourceDir = fileparts(mfilename('fullpath'));
    candidates{end+1} = fullfile( ...
        sourceDir,'..','regions','JP',name);
    if isdeployed
        candidates{end+1} = fullfile( ...
            ctfroot,'regions','JP',name);
        candidates{end+1} = fullfile( ...
            ctfroot,'SAGEhydrology','regions','JP',name);
    end

    for i = 1:numel(candidates)
        if isfile(candidates{i})
            copyfile(candidates{i},target,'f');
            logFcn(['Installed packaged MERV-Jp ' ...
                'HydroATLAS attributes.']);
            return
        end
    end
    error('data_helpers:JPMissingAttributes', ...
        ['Could not locate the packaged MERV-Jp HydroATLAS ' ...
        'attribute workbook.']);
end

function out = local_download_pl(cfg,stream,ui)

    logFcn = local_ui_log(ui);

    if ~strcmp(stream,'daily')
        error('data_helpers:PLOnlyDaily', ...
            ['CAMELS-PL support is daily only. ' ...
            'Use stream = ''daily''.']);
    end

    dirD = local_cfg_dirD(cfg,'CAMELS_PL');
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily','timeseries'));

    if ~local_ui_confirm(ui, ...
            'Download Poland data', ...
            ['CAMELS-PL is downloaded from ' ...
            'Zenodo as one ZIP file. ' ...
            'The archive contains attribute ' ...
            'CSV files and daily ' ...
            ['hydrometeorological time series ' ...
            'for 354 catchments.'] newline ...
            'Continue?'])
        out = struct('ok',false,'canceled',true);
        return
    end

    downloadDir = local_default_download_dir();
    zipName = 'CAMELS-PL.zip';
    zipFile = fullfile(downloadDir,zipName);
    url = ['https://zenodo.org/records/20133183/files/' ...
        'CAMELS-PL.zip?download=1'];
    expectedMD5 = '0ff7d9635c941f4aac49c191376502c9';

    userCanceled = false;
    lastUI = tic;
    d = local_progress_dialog(ui, ...
        'Downloading Poland data', ...
        'Starting CAMELS-PL download ...');

    try
        if local_pl_install_complete(dirD)
            logFcn(['CAMELS-PL appears complete; ' ...
                'skipping download/install.']);
            out = struct('ok',true, ...
                'region','CAMELS_PL', ...
                'dirD',dirD, ...
                'dirM',fullfile(dirD,'daily','timeseries'), ...
                'dirQ',fullfile(dirD,'daily','timeseries'), ...
                'archive','', ...
                'unzipDir','');
            local_close_progress(d);
            return
        end

        if isfile(zipFile)
            try
                local_verify_md5(zipFile,expectedMD5,logFcn);
                logFcn(['Using verified ' ...
                    'existing file: ' zipFile]);
            catch
                logFcn(['Existing CAMELS-PL ' ...
                    'archive failed MD5; deleting.']);
                delete(zipFile);
            end
        end

        if ~isfile(zipFile)
            local_download_file_retry(url,zipFile, ...
                'CAMELS-PL archive',d, ...
                @() userCanceled, ...
                @ui_progress,logFcn,3);
        end

        local_verify_md5(zipFile,expectedMD5,logFcn);

        unzipDir = fullfile(downloadDir, ...
            'CAMELS_PL_download');
        if isfolder(unzipDir)
            try
                rmdir(unzipDir,'s');
            catch
            end
        end
        local_mkdir(unzipDir);

        try
            d.Indeterminate = true;
            d.Message = ['Unzipping ' ...
                'CAMELS-PL archive. ' ...
                'This can take a while ...'];
            drawnow limitrate nocallbacks
        catch
        end

        unzip(zipFile,unzipDir);

        local_install_pl_files(unzipDir,dirD,logFcn);

        okAttr = local_all_files_exist(dirD, ...
            local_pl_metadata_files());
        okDaily = isfolder(fullfile(dirD, ...
            'daily','timeseries')) ...
            && local_count_files(fullfile(dirD, ...
            'daily','timeseries'), ...
            {'CAMELS_PL_hydromet_timeseries_*.csv'}) >= 354;

        out = struct('ok',okAttr ...
            && okDaily,'region','CAMELS_PL', ...
            'dirD',dirD, ...
            'dirM',fullfile(dirD,'daily','timeseries'), ...
            'dirQ',fullfile(dirD,'daily','timeseries'), ...
            'archive',zipFile, ...
            'unzipDir',unzipDir);

        if out.ok
            nTs = local_count_files( ...
                fullfile(dirD,'daily','timeseries'), ...
                {'CAMELS_PL_hydromet_timeseries_*.csv'});
            logFcn(sprintf(['CAMELS-PL install finished: ' ...
                '%d daily time-series CSV files.'],nTs));

            local_cleanup_download_file(zipFile,logFcn);
            local_cleanup_download_folder(unzipDir,logFcn);

            out.archive = '';
            out.unzipDir = '';
        else
            logFcn(['CAMELS-PL install ' ...
                'finished, but ' ...
                'some required files ' ...
                'appear to be missing.']);
        end

    catch ME
        try
            if exist('zipFile','var')
                local_cleanup_download_file(zipFile,logFcn);
            end
        catch
        end
        try
            if exist('unzipDir','var')
                local_cleanup_download_folder(unzipDir,logFcn);
            end
        catch
        end

        logFcn('CAMELS-PL install failed:');
        logFcn(ME.message);

        out = struct('ok',false,'error',ME.message);
    end

    local_close_progress(d);

    function ui_progress(info)
        [userCanceled,lastUI] = ...
            local_update_progress_dialog( ...
            d,info,lastUI,userCanceled);
    end

end

function local_install_col_files(files,dirD,logFcn)
    local_mkdir(dirD);
    local_mkdir(fullfile(dirD,'daily'));
    local_mkdir(fullfile(dirD,'daily','timeseries'));

    % 1) Attribute workbook from the attributes ZIP.
    idx = find(strcmp({files.key}, ...
        'attributes_zip'),1);
    if isempty(idx) ...
            || ~isfield(files(idx),'unzipDir') ...
            || ~isfolder(files(idx).unzipDir)
        error('data_helpers:COLMissingAttributeZip', ...
            ['CAMELS-COL attribute archive ' ...
            'was not extracted.']);
    end

    q = dir(fullfile(files(idx).unzipDir, ...
        '**','01_CAMELS_COL_Attributes.xlsx'));
    if isempty(q)
        error('data_helpers:COLMissingAttributeWorkbook', ...
            ['Could not find ' ...
            '01_CAMELS_COL_Attributes.xlsx ' ...
            'inside the attributes archive.']);
    end
    copyfile(fullfile(q(1).folder,q(1).name), ...
        fullfile(dirD,q(1).name),'f');
    logFcn(['Copied CAMELS-COL attribute workbook to: ' dirD]);

    % 2) Individual attribute CSV files downloaded directly from Zenodo.
    metaFiles = local_col_metadata_files();
    for i = 1:numel(files)
        if files(i).isZip
            continue
        end
        if isfield(files(i),'localFile') ...
                && isfile(files(i).localFile)
            copyfile(files(i).localFile, ...
                fullfile(dirD,files(i).fileName),'f');
        end
    end
    nMeta = 0;
    for i = 1:numel(metaFiles)
        if isfile(fullfile(dirD,metaFiles{i}))
            nMeta = nMeta + 1;
        end
    end
    logFcn(sprintf(['Installed %d CAMELS-COL ' ...
        'metadata/attribute files.'],nMeta));

    % 3) Hydrometeorological TXT files -> Data/CAMELS_COL/daily/timeseries.
    idx = find(strcmp({files.key}, ...
        'hydromet_zip'),1);
    if isempty(idx) ...
            || ~isfield(files(idx),'unzipDir') ...
            || ~isfolder(files(idx).unzipDir)
        error('data_helpers:COLMissingHydrometZip', ...
            ['CAMELS-COL hydrometeorological ' ...
            'archive was not extracted.']);
    end

    q = dir(fullfile(files(idx).unzipDir, ...
        '**','Hydromet_data_*.txt'));
    if isempty(q)
        error('data_helpers:COLMissingHydrometFiles', ...
            ['Could not find Hydromet_data_*.txt ' ...
            'inside the CAMELS-COL hydrometeorological archive.']);
    end

    dst = fullfile(dirD,'daily','timeseries');
    local_mkdir(dst);
    n = 0;
    for i = 1:numel(q)
        if q(i).isdir
            continue
        end
        copyfile(fullfile(q(i).folder,q(i).name), ...
            fullfile(dst,q(i).name),'f');
        n = n + 1;
    end
    logFcn(sprintf(['Copied %d CAMELS-COL daily ' ...
        'hydrometeorological files to: %s'],n,dst));

end

function local_download_file_retry_nz(url, ...
    targetFile,label,d,cancelFcn, ...
    progressFcn,logFcn,nTry)
%LOCAL_DOWNLOAD_FILE_RETRY_NZ CAMELS-NZ/Figshare-only downloader.
% This avoids changing the generic downloader used by the other regions.

    if nargin < 8 || isempty(nTry)
        nTry = 3;
    end

    lastME = [];
    for k = 1:nTry
        if ~isempty(cancelFcn) && cancelFcn()
            error('data_helpers:DownloadCanceled', ...
                'Download canceled by user.');
        end

        try
            if k > 1
                logFcn(sprintf(['Retrying %s download ' ...
                    '(%d/%d) ...'],label,k,nTry));
                pause(3);
            end

            local_download_file_nz(url,targetFile,label, ...
                d,cancelFcn,progressFcn,logFcn);
            return

        catch ME
            lastME = ME;
            local_delete_nz_download_variants(targetFile);
            logFcn(sprintf(['%s download attempt %d ' ...
                'failed: %s'],label,k,ME.message));
        end
    end

    rethrow(lastME);
end

function local_download_file_nz(url,targetFile,label, ...
    d,cancelFcn,progressFcn,logFcn)
%LOCAL_DOWNLOAD_FILE_NZ Download CAMELS-NZ Figshare files automatically.
% Uses a Windows PowerShell/.NET downloader with certificate-revocation
% checking disabled.  This is used only for CAMELS-NZ/Figshare.

    if nargin < 7 || isempty(logFcn)
        logFcn = @(s) fprintf('%s\n',char(string(s)));
    end

    local_mkdir(fileparts(targetFile));
    local_delete_nz_download_variants(targetFile);

    if ~isempty(progressFcn)
        progressFcn(struct('label', ...
            ['Downloading ' label], ...
            'frac',NaN,'speedMBs',NaN));
    end

    ok = false;
    msg = '';

    % First try PowerShell/.NET. This avoids Windows curl schannel
    % revocation failures and follows Figshare redirects.
    if ispc
        try
            psUrl = local_ps_quote(url);
            psOut = local_ps_quote(targetFile);
            cmd = ['powershell -NoProfile -ExecutionPolicy Bypass ' ...
                '-Command "' ...
                '$ProgressPreference=''SilentlyContinue''; ' ...
                '[Net.ServicePointManager]::SecurityProtocol=' ...
                '[Net.SecurityProtocolType]::Tls12; ' ...
                '[Net.ServicePointManager]::CheckCertificateRevocationList=$false; ' ...
                '[System.Net.ServicePointManager]::ServerCertificateValidationCallback={ $true }; ' ...
                '$wc=New-Object System.Net.WebClient; ' ...
                '$wc.Headers.Add(''User-Agent'',''Mozilla/5.0''); ' ...
                '$wc.DownloadFile(' psUrl ',' psOut ')"'];

            [status,out] = system(cmd);
            ok = status == 0 ...
                && isfile(targetFile) ...
                && dir(targetFile).bytes > 0;
            if ~ok
                msg = ['powershell/.NET: ' out];
            end
        catch ME
            msg = ['powershell/.NET: ' ME.message];
        end
    end

    % Second try generic downloader, but with browser fallback disabled by
    % catching failures here.  This is mostly useful on Mac/Linux.
    if ~ok
        try
            local_download_file(url,targetFile,label, ...
                d,cancelFcn,progressFcn,logFcn);
            ok = local_accept_nz_download_variant(targetFile,logFcn);
        catch ME
            msg = [msg ' | generic: ' ME.message];
        end
    end

    % Last resort: open browser and wait for ANY valid Chrome/Edge name.
    % This is still automatic: the code detects " (1)" files and moves them.
    if ~ok
        try
            web(url,'-browser');
            ok = local_wait_for_nz_browser_file( ...
                targetFile,300,cancelFcn,progressFcn,label,logFcn);
            if ~ok
                msg = [msg ' | browser wait timed out'];
            end
        catch ME
            msg = [msg ' | browser: ' ME.message];
        end
    end

    if ~ok
        error('data_helpers:NZDownloadFailed', ...
            'Could not download %s to %s. %s', ...
            label,targetFile,msg);
    end

    if ~isempty(progressFcn)
        progressFcn(struct('label', ...
            [label ' downloaded'], ...
            'frac',1,'speedMBs',NaN));
    end
end

function ok = local_wait_for_nz_browser_file( ...
    targetFile,maxSeconds,cancelFcn,progressFcn,label,logFcn)

    ok = false;
    t0 = tic;
    lastBytes = -1;
    stableCount = 0;

    while toc(t0) < maxSeconds
        if ~isempty(cancelFcn) && cancelFcn()
            return
        end

        [cand,bytes] = local_find_nz_download_variant(targetFile);
        if ~isempty(cand) && bytes > 0
            if bytes == lastBytes
                stableCount = stableCount + 1;
            else
                lastBytes = bytes;
                stableCount = 0;
            end

            if stableCount >= 3
                ok = local_move_nz_variant(cand,targetFile,logFcn);
                return
            end
        end

        if ~isempty(progressFcn)
            progressFcn(struct('label', ...
                ['Waiting for browser download: ' label], ...
                'frac',NaN,'speedMBs',NaN));
        end
        pause(2)
    end
end

function ok = local_accept_nz_download_variant(targetFile,logFcn)
    [cand,bytes] = local_find_nz_download_variant(targetFile);
    ok = ~isempty(cand) && bytes > 0 ...
        && local_move_nz_variant(cand,targetFile,logFcn);
end

function [cand,bytes] = local_find_nz_download_variant(targetFile)
    cand = '';
    bytes = 0;

    [folder,name,ext] = fileparts(targetFile);
    bases = unique({[name ext],name},'stable');
    names = {};
    for b = 1:numel(bases)
        base = bases{b};
        if isempty(base)
            continue
        end
        names{end+1} = base; %#ok<AGROW>
        names{end+1} = [base '.zip']; %#ok<AGROW>
        for k = 1:20
            names{end+1} = sprintf('%s (%d)',base,k); %#ok<AGROW>
            names{end+1} = sprintf('%s (%d).zip',base,k); %#ok<AGROW>
            if endsWith(base,'.zip','IgnoreCase',true)
                nozip = extractBefore(base,strlength(base)-3);
                names{end+1} = sprintf('%s (%d).zip',nozip,k); %#ok<AGROW>
            end
        end
    end
    names = unique(names,'stable');

    bestBytes = -1;
    bestFile = '';
    for i = 1:numel(names)
        f = fullfile(folder,names{i});
        tmp = [f '.crdownload'];
        if isfile(f) && ~isfile(tmp)
            info = dir(f);
            if info.bytes > bestBytes
                bestBytes = info.bytes;
                bestFile = f;
            end
        end
    end

    if bestBytes > 0
        cand = bestFile;
        bytes = bestBytes;
    end
end

function ok = local_move_nz_variant(cand,targetFile,logFcn)
    ok = false;
    try
        if strcmp(cand,targetFile)
            ok = isfile(targetFile) && dir(targetFile).bytes > 0;
            return
        end
        if isfile(targetFile)
            try
                delete(targetFile);
            catch
            end
        end
        movefile(cand,targetFile,'f');
        logFcn(['Using browser-downloaded file: ' cand]);
        ok = isfile(targetFile) && dir(targetFile).bytes > 0;
    catch ME
        logFcn(['Could not move browser-downloaded file: ' ...
            ME.message]);
    end
end

function local_delete_nz_download_variants(targetFile)
    [folder,name,ext] = fileparts(targetFile);
    bases = unique({[name ext],name},'stable');
    for b = 1:numel(bases)
        base = bases{b};
        if isempty(base)
            continue
        end
        pats = {base,[base '.zip'],[base '*.crdownload']};
        for p = 1:numel(pats)
            q = dir(fullfile(folder,pats{p}));
            for i = 1:numel(q)
                if ~q(i).isdir
                    try
                        delete(fullfile(q(i).folder,q(i).name));
                    catch
                    end
                end
            end
        end
    end
end

function s = local_ps_quote(x)
    s = char(string(x));
    s = strrep(s,'''','''''');
    s = ['''' s ''''];
end


function local_download_file_retry(url, ...
    targetFile,label,d,cancelFcn, ...
    progressFcn,logFcn,nTry)

    if nargin < 8 ...
            || isempty(nTry)
        nTry = 3;
    end
    
    lastME = [];
    
    for k = 1:nTry
    
        if ~isempty(cancelFcn) ...
                && cancelFcn()
            error(['data_helpers:' ...
                'DownloadCanceled'], ...
                'Download canceled by user.');
        end
    
        try
            if k > 1
                logFcn(sprintf(['Retrying ' ...
                    '%s download (%d/%d) ...'], ...
                    label,k,nTry));
                pause(3);
            end
    
            local_download_file(url, ...
                targetFile,label,d, ...
                cancelFcn,progressFcn, ...
                logFcn);
            return
    
        catch ME
            lastME = ME;
    
            % Delete partial/corrupt file before retry
            if isfile(targetFile)
                try
                    delete(targetFile);
                catch
                end
            end
    
            logFcn(sprintf(['%s download ' ...
                'attempt %d failed: %s'], ...
                label,k,ME.message));
        end
    end
    
    rethrow(lastME);
end

function ok = local_br_component_ok(dirD,subdir,nMin)

    d = fullfile(dirD,'daily',subdir);
    ok = isfolder(d) ...
        && local_count_files(d, ...
        {'*.txt'}) >= nMin;
    
end

function files = local_camels_nz_files(stream)
    files = struct( ...
        'key',{}, ...
        'label',{}, ...
        'url',{}, ...
        'fileName',{}, ...
        'isZip',{}, ...
        'md5',{});

    files(end+1) = struct( ...
        'key','attributes', ...
        'label','CAMELS-NZ catchment attributes', ...
        'url',['https://figshare.canterbury.ac.nz/' ...
        'ndownloader/files/56902355'], ...
        'fileName',['CAMELS_NZ_Catchment_' ...
        'Atrributes.zip'], ...
        'isZip',true, ...
        'md5','00f877983df4be56fe490e98ed44fa84');

    if strcmp(stream,'hourly')
        files(end+1) = struct( ...
            'key','precipitation', ...
            'label','CAMELS-NZ hourly precipitation', ...
            'url',['https://figshare.canterbury.ac.nz/' ...
            'ndownloader/files/61527853'], ...
            'fileName',['CAMELS_NZ_hourly_' ...
            'Precipitation.zip'], ...
            'isZip',true, ...
            'md5','430c0203258ebe4e506a58eb4c13baee');

        files(end+1) = struct( ...
            'key','temperature', ...
            'label','CAMELS-NZ hourly temperature', ...
            'url',['https://figshare.canterbury.ac.nz/' ...
            'ndownloader/files/61527862'], ...
            'fileName',['CAMELS_NZ_hourly_' ...
            'Temperature.zip'], ...
            'isZip',true, ...
            'md5','f6a5cfd6766ba6a27b7eaec20842c3bc');

        files(end+1) = struct( ...
            'key','pet', ...
            'label','CAMELS-NZ hourly PET', ...
            'url',['https://figshare.canterbury.ac.nz/' ...
            'ndownloader/files/61527847'], ...
            'fileName',['CAMELS_NZ_hourly_' ...
            'PET.zip'], ...
            'isZip',true, ...
            'md5','5c717c2f31f8d1447b3289b6397f354b');

        files(end+1) = struct( ...
            'key','streamflow', ...
            'label','CAMELS-NZ hourly streamflow', ...
            'url',['https://figshare.canterbury.ac.nz/' ...
            'ndownloader/files/61527859'], ...
            'fileName',['CAMELS_NZ_hourly_' ...
            'Streamflow.zip'], ...
            'isZip',true, ...
            'md5','a9c5bcc6440198c9e00af051cb86d867');
    else
        files(end+1) = struct( ...
            'key','precipitation', ...
            'label','CAMELS-NZ daily precipitation', ...
            'url',['https://figshare.canterbury.ac.nz/' ...
            'ndownloader/files/56902367'], ...
            'fileName',['CAMELS_NZ_daily_' ...
            'Precipitation.zip'], ...
            'isZip',true, ...
            'md5','');

        files(end+1) = struct( ...
            'key','temperature', ...
            'label','CAMELS-NZ daily temperature', ...
            'url',['https://figshare.canterbury.ac.nz/' ...
            'ndownloader/files/56902382'], ...
            'fileName',['CAMELS_NZ_daily_' ...
            'Temperature.zip'], ...
            'isZip',true, ...
            'md5','');

        files(end+1) = struct( ...
            'key','relative_humidity', ...
            'label','CAMELS-NZ daily relative humidity', ...
            'url',['https://figshare.canterbury.ac.nz/' ...
            'ndownloader/files/56902370'], ...
            'fileName',['CAMELS_NZ_daily_' ...
            'Relative_Humidity.zip'], ...
            'isZip',true, ...
            'md5','');

        files(end+1) = struct( ...
            'key','streamflow', ...
            'label','CAMELS-NZ daily streamflow', ...
            'url',['https://figshare.canterbury.ac.nz/' ...
            'ndownloader/files/56902373'], ...
            'fileName',['CAMELS_NZ_daily_' ...
            'Streamflow.zip'], ...
            'isZip',true, ...
            'md5','');
    end
end

function files = local_nz_metadata_files()
    files = { ...
        '1.CAMELS_NZ_Catchment_information.csv', ...
        '2.CAMELS_NZ_Climatic_attribute.csv', ...
        '3.CAMELS_NZ_Landcover_attribute.csv', ...
        '4.CAMELS_NZ_Geology.csv', ...
        ['5.CAMELS_NZ_Anthropogenic_' ...
        'attribute.csv']};
end

function files = local_kr_metadata_files()
    files = {'attributes_climate_ERA5Land.csv', ...
        'attributes_climate_obs.csv','attributes_dam.csv', ...
        'attributes_general.csv','attributes_HydroATLAS.csv'};
end

function tf = local_has_kr_timeseries(dirD)
    d = fullfile(dirD,'hourly','timeseries');
    tf = isfolder(d) ...
        && local_count_files(d,{'*.csv'}) >= 178;
end

function tf = local_kr_install_complete(dirD)
    tf = local_all_files_exist(dirD,local_kr_metadata_files()) ...
        && local_has_kr_timeseries(dirD);
end

function files = local_camels_kr_files()
    root = 'https://zenodo.org/api/records/15073264/files/';
    files = struct('key',{},'label',{},'url',{}, ...
        'fileName',{},'isZip',{},'md5',{});
    files(end+1) = local_kr_file(root,'climate_era5', ...
        'CAMELSH-KR ERA5-Land climate attributes', ...
        'attributes_climate_ERA5Land.csv',false, ...
        '59e2b367819990f2c5b92cdf15fd1fd6');
    files(end+1) = local_kr_file(root,'climate_obs', ...
        'CAMELSH-KR observed climate attributes', ...
        'attributes_climate_obs.csv',false, ...
        '4cf4c0e7b44f7ac22fac2540e4f6b0fa');
    files(end+1) = local_kr_file(root,'dam', ...
        'CAMELSH-KR dam attributes','attributes_dam.csv',false, ...
        '3814b0a05c44397f1fc59ecdcdd5a326');
    files(end+1) = local_kr_file(root,'general', ...
        'CAMELSH-KR general attributes','attributes_general.csv',false, ...
        '2749de06bd195bd37bfbc325cae0700c');
    files(end+1) = local_kr_file(root,'hydroatlas', ...
        'CAMELSH-KR HydroATLAS attributes', ...
        'attributes_HydroATLAS.csv',false, ...
        '34cdabed3fa8e0e63822ebc3a548f113');
    files(end+1) = local_kr_file(root,'timeseries', ...
        'CAMELSH-KR hourly time series','timeseries.zip',true, ...
        '754f5de4fa828919a87c44b18f07969d');
end

function f = local_kr_file(root,key,label,fileName,isZip,md5)
    f = struct('key',key,'label',label, ...
        'url',[root fileName '/content'], ...
        'fileName',fileName,'isZip',isZip,'md5',md5);
end

function tf = local_has_nz_timeseries(dirD,stream)

    dirT = fullfile(dirD,stream, ...
        'timeseries');
    dP = fullfile(dirT, ...
        'precipitation');
    dT = fullfile(dirT, ...
        'temperature');
    dQ = fullfile(dirT, ...
        'streamflow');

    if strcmp(stream,'hourly')
        dX = fullfile(dirT,'pet');
        pP = {'precipitation_station_id_*.csv'};
        pT = {'temperature_station_id_*.csv'};
        pX = {'PET_station_id_*.csv'};
        pQ = {'flow_station_id_*.csv'};
    elseif strcmp(stream,'daily')
        dX = fullfile(dirT, ...
            'relative_humidity');
        pP = {'daily_precipitation_station_id_*.csv'};
        pT = {'daily_temperature_station_id_*.csv'};
        pX = {'daily_RH_station_id_*.csv'};
        pQ = {'daily_flow_station_id_*.csv'};
    else
        tf = false;
        return
    end

    tf = isfolder(dP) ...
        && isfolder(dT) ...
        && isfolder(dX) ...
        && isfolder(dQ) ...
        && local_count_files(dP,pP) > 0 ...
        && local_count_files(dT,pT) > 0 ...
        && local_count_files(dX,pX) > 0 ...
        && local_count_files(dQ,pQ) > 0;
end

function tf = local_nz_install_complete( ...
    dirD,stream)

    tf = local_all_files_exist(dirD, ...
        local_nz_metadata_files()) ...
        && local_has_nz_timeseries( ...
        dirD,stream);
end

function local_install_nz_files( ...
    files,dirD,stream,logFcn)

    local_mkdir(dirD);
    dirT = fullfile(dirD,stream, ...
        'timeseries');
    local_mkdir(dirT);

    % 1) Static attributes -> Data/CAMELS_NZ
    idx = find(strcmp({files.key}, ...
        'attributes'),1);
    if isempty(idx) ...
            || ~isfield(files(idx),'unzipDir') ...
            || ~isfolder(files(idx).unzipDir)
        error('data_helpers:NZMissingAttributeZip', ...
            ['CAMELS-NZ attribute ' ...
            'archive was not extracted.']);
    end

    attrFiles = local_nz_metadata_files();
    for i = 1:numel(attrFiles)
        q = dir(fullfile(files(idx).unzipDir, ...
            '**',attrFiles{i}));
        if isempty(q)
            error(['data_helpers:' ...
                'NZMissingAttributeFile'], ...
                ['Missing CAMELS-NZ ' ...
                'attribute file ' ...
                'inside archive: %s'], ...
                attrFiles{i});
        end
        copyfile(fullfile(q(1).folder, ...
            q(1).name), ...
            fullfile(dirD,attrFiles{i}),'f');
    end
    logFcn(sprintf(['Copied %d CAMELS-NZ ' ...
        'static attribute CSV files.'], ...
        numel(attrFiles)));

    if strcmp(stream,'hourly')
        pP = {'precipitation_station_id_*.csv'};
        pT = {'temperature_station_id_*.csv'};
        pX = {'PET_station_id_*.csv'};
        pQ = {'flow_station_id_*.csv'};
        xKey = 'pet';
        xFolder = 'pet';
    else
        pP = {'daily_precipitation_station_id_*.csv'};
        pT = {'daily_temperature_station_id_*.csv'};
        pX = {'daily_RH_station_id_*.csv'};
        pQ = {'daily_flow_station_id_*.csv'};
        xKey = 'relative_humidity';
        xFolder = 'relative_humidity';
    end

    % 2) Resolution-specific time series.
    local_install_nz_timeseries_zip(files, ...
        'precipitation', ...
        fullfile(dirT,'precipitation'), ...
        pP,logFcn);

    local_install_nz_timeseries_zip(files, ...
        'temperature', ...
        fullfile(dirT,'temperature'), ...
        pT,logFcn);

    local_install_nz_timeseries_zip(files, ...
        xKey,fullfile(dirT,xFolder), ...
        pX,logFcn);

    local_install_nz_timeseries_zip(files, ...
        'streamflow', ...
        fullfile(dirT,'streamflow'), ...
        pQ,logFcn);
end

function local_install_nz_timeseries_zip(files, ...
    key,dst,patterns,logFcn)

    idx = find(strcmp({files.key},key),1);
    if isempty(idx) ...
            || ~isfield(files(idx), ...
            'unzipDir') ...
            || ~isfolder(files(idx).unzipDir)
        error(['data_helpers:' ...
            'NZMissingTimeseriesZip'], ...
            ['CAMELS-NZ %s archive ' ...
            'was not extracted.'],key);
    end

    if isfolder(dst)
        try
            rmdir(dst,'s');
        catch
        end
    end
    local_mkdir(dst);

    n = 0;
    for ip = 1:numel(patterns)
        q = dir(fullfile(files(idx).unzipDir, ...
            '**',patterns{ip}));
        for j = 1:numel(q)
            if q(j).isdir
                continue
            end
            copyfile(fullfile(q(j).folder, ...
                q(j).name), ...
                fullfile(dst,q(j).name),'f');
            n = n + 1;
        end
    end

    if n == 0
        error(['data_helpers:' ...
            'NZNoTimeseriesFiles'], ...
            ['No CAMELS-NZ %s ' ...
            'files matching expected ' ...
            'patterns were found.'],key);
    end

    logFcn(sprintf(['Copied %d CAMELS-NZ ' ...
        '%s CSV files.'],n,key));
end

function files = local_ch_metadata_files()
    files = { ...
        'CAMELS_CH_climate_attributes_obs.csv', ...
        'CAMELS_CH_geology_attributes.csv', ...
        'CAMELS_CH_glacier_attributes.csv', ...
        'CAMELS_CH_humaninfluence_attributes.csv', ...
        'CAMELS_CH_hydrogeology_attributes.csv', ...
        'CAMELS_CH_hydrology_attributes_obs.csv', ...
        'CAMELS_CH_landcover_attributes.csv', ...
        'CAMELS_CH_soil_attributes.csv', ...
        'CAMELS_CH_topographic_attributes.csv'};
end

function tf = local_has_ch_timeseries(dirD)
    dObs = fullfile(dirD,'daily','timeseries', ...
        'observation_based');
    dSim = fullfile(dirD,'daily','timeseries', ...
        'simulation_based');
    tf = isfolder(dObs) ...
        && isfolder(dSim) ...
        && local_count_files(dObs, ...
        {'CAMELS_CH_obs_based_*.csv'}) > 0 ...
        && local_count_files(dSim, ...
        {'CAMELS_CH_sim_based_*.csv'}) > 0;
end

function tf = local_ch_install_complete(dirD)
    tf = local_all_files_exist(dirD, ...
        local_ch_metadata_files()) ...
        && local_has_ch_timeseries(dirD);
end

function local_install_ch_files(dataDir, ...
    dirD,logFcn)
    local_mkdir(dirD);

    attrSrc = fullfile(dataDir, ...
        'static_attributes');
    if ~isfolder(attrSrc)
        error(['data_helpers:' ...
            'CHMissingAttributes'], ...
            ['Could not find ' ...
            'static_attributes ' ...
            'inside CAMELS-CH archive.']);
    end

    attrFiles = local_ch_metadata_files();
    for i = 1:numel(attrFiles)
        src = fullfile(attrSrc,attrFiles{i});
        if ~isfile(src)
            error(['data_helpers:' ...
                'CHMissingAttributeFile'], ...
                ['Missing CAMELS-CH ' ...
                'attribute file ' ...
                'inside archive: %s'], ...
                attrFiles{i});
        end
        copyfile(src,fullfile(dirD, ...
            attrFiles{i}),'f');
    end
    logFcn(sprintf(['Copied %d CAMELS-CH ' ...
        'static attribute CSV files.'], ...
        numel(attrFiles)));

    tsSrc = fullfile(dataDir, ...
        'timeseries');
    obsSrc = fullfile(tsSrc, ...
        'observation_based');
    simSrc = fullfile(tsSrc, ...
        'simulation_based');
    if ~isfolder(obsSrc)
        error(['data_helpers:' ...
            'CHMissingObsTimeseries'], ...
            ['Could not find timeseries/' ...
            'observation_based ' ...
            'inside CAMELS-CH archive.']);
    end
    if ~isfolder(simSrc)
        error(['data_helpers:' ...
            'CHMissingSimTimeseries'], ...
            ['Could not find ' ...
            'timeseries/simulation_based ' ...
            'inside CAMELS-CH archive.']);
    end

    obsDst = fullfile(dirD,'daily','timeseries', ...
        'observation_based');
    simDst = fullfile(dirD,'daily','timeseries', ...
        'simulation_based');

    local_mkdir(fileparts(obsDst));
    local_copy_folder_contents(obsSrc, ...
        obsDst,logFcn);
    local_mkdir(fileparts(simDst));
    local_copy_folder_contents(simSrc, ...
        simDst,logFcn);
end

function dataDir = local_find_ch_payload_dir(unzipDir)
    dataDir = '';
    if ~isfolder(unzipDir)
        return
    end

    q = dir(fullfile(unzipDir,'**', ...
        'static_attributes'));
    for i = 1:numel(q)
        if q(i).isdir
            cand = q(i).folder;
            if isfolder(fullfile(cand,'timeseries'))
                dataDir = cand;
                return
            end
        end
    end

    q = dir(fullfile(unzipDir,'**','camels_ch'));
    for i = 1:numel(q)
        if q(i).isdir
            cand = fullfile(q(i).folder,q(i).name);
            if isfolder(fullfile(cand, ...
                    'static_attributes')) ...
                    && isfolder(fullfile(cand, ...
                    'timeseries'))
                dataDir = cand;
                return
            end
        end
    end

    q = dir(fullfile(unzipDir,'**','CAMELS_CH'));
    for i = 1:numel(q)
        if q(i).isdir
            cand = fullfile(q(i).folder,q(i).name);
            if isfolder(fullfile(cand, ...
                    'static_attributes')) ...
                    && isfolder(fullfile(cand, ...
                    'timeseries'))
                dataDir = cand;
                return
            end
        end
    end
end

function files = local_de_metadata_files(stream)

    if nargin < 1 ...
            || isempty(stream)
        stream = 'hourly';
    end
    stream = lower(strtrim(char(stream)));
    switch stream
        case 'daily'
            files = { ...
                'CAMELS_DE_climatic_attributes.csv', ...
                'CAMELS_DE_humaninfluence_attributes.csv', ...
                'CAMELS_DE_hydrogeology_attributes.csv', ...
                'CAMELS_DE_hydrologic_attributes.csv', ...
                'CAMELS_DE_landcover_attributes.csv', ...
                'CAMELS_DE_soil_attributes.csv', ...
                'CAMELS_DE_topographic_attributes.csv'};
        case 'hourly'
            files = { ...
                'CAMELS_DE_1h_climatic_attributes.csv', ...
                'CAMELS_DE_1h_humaninfluence_attributes.csv', ...
                'CAMELS_DE_1h_hydrogeology_attributes.csv', ...
                'CAMELS_DE_1h_hydrologic_attributes.csv', ...
                'CAMELS_DE_1h_landcover_attributes.csv', ...
                'CAMELS_DE_1h_soil_attributes.csv', ...
                'CAMELS_DE_1h_topographic_attributes.csv'};
        otherwise
            error('data_helpers:DEBadStream', ...
                'CAMELS-DE stream must be daily or hourly.');
    end
end

function files = local_fr_metadata_files()
    files = { ...
        'CAMELS_FR_geology_attributes.csv', ...
        'CAMELS_FR_human_influences_dams.csv', ...
        'CAMELS_FR_hydrogeology_attributes.csv', ...
        'CAMELS_FR_land_cover_attributes.csv', ...
        'CAMELS_FR_site_general_attributes.csv', ...
        'CAMELS_FR_soil_general_attributes.csv', ...
        'CAMELS_FR_soil_quantiles_attributes.csv', ...
        'CAMELS_FR_station_general_attributes.csv', ...
        'CAMELS_FR_topography_general_attributes.csv', ...
        'CAMELS_FR_topography_quantiles_attributes.csv', ...
        '00_description_land_cover_classes.txt', ...
        '00_description_geology_classes.txt'};
end


function files = local_ind_metadata_files()

    files = { ...
        'camels_ind_anth.csv', ...
        'camels_ind_clim.csv', ...
        'camels_ind_geol.csv', ...
        'camels_ind_hydro.csv', ...
        'camels_ind_land.csv', ...
        'camels_ind_name.csv', ...
        'camels_ind_soil.csv', ...
        'camels_ind_topo.csv'};

end

function files = local_br_metadata_files()

    files = { ...
        'camels_br_climate.txt', ...
        'camels_br_geology.txt', ...
        'camels_br_human_intervention.txt', ...
        'camels_br_hydrology.txt', ...
        'camels_br_land_cover.txt', ...
        'camels_br_location.txt', ...
        'camels_br_quality_check.txt', ...
        'camels_br_soil.txt', ...
        'camels_br_topography.txt'};

end

function tf = local_br_needs_download(dirD,f)

    tf = true;
    
    switch f.key
        case 'attributes'
            tf = ~local_all_files_exist(dirD, ...
                local_br_metadata_files());
    
        case 'streamflow'
            tf = ~local_br_component_ok(dirD, ...
                'streamflow',897);
    
        case 'precipitation'
            tf = ~local_br_component_ok(dirD, ...
                'precipitation',897);
    
        case 'pet'
            tf = ~local_br_component_ok(dirD, ...
                'pet',897);
    
        case 'temperature'
            tf = ~local_br_component_ok(dirD, ...
                'temperature',897);
    
        case {'soil_moisture','readme'}
            % Not required for current SAGE readiness.
            tf = false;
    end

end

function tf = local_br_should_install(files,key)

    idx = find(strcmp({files.key},key),1);
    tf = ~isempty(idx) ...
        && isfield(files,'install') ...
        && files(idx).install;

end

function local_websave(dst,url)

    try
        websave(dst,url);
    catch ME
        error('data_helpers:CAMELS_SE_downloadFailed', ...
            ['Could not download ' ...
            'CAMELS-SE file:\n%s\n\n' ...
             'MATLAB error:\n%s'], ...
             url,ME.message);
    end
    
end

function local_reset_dir(d)

    if isfolder(d)
        try
            rmdir(d,'s');
        catch
        end
    end
    mkdir(d);

end

function removed = local_remove_dir_retry(folder,logFcn)

    removed = true;

    if ~isfolder(folder)
        return
    end
    
    lastMsg = '';
    
    for i = 1:6
        try
            [ok,msg] = rmdir(folder,'s');
            if ok
                removed = true;
                return
            end
            lastMsg = msg;
        catch ME
            lastMsg = ME.message;
        end
        pause(0.5*i);
        drawnow;
    end
    logFcn(['Warning: temporary ' ...
        'folder remains: ' folder]);
    logFcn(['Reason: ' lastMsg]);
    removed = false;

end

function removed = local_remove_file_retry(fileName,logFcn)
%LOCAL_REMOVE_FILE_RETRY Delete a file and verify that it is gone.

    removed = true;
    if ~isfile(fileName)
        return
    end

    lastMsg = '';
    for i = 1:6
        try
            delete(fileName);
            if ~isfile(fileName)
                return
            end
            lastMsg = 'The file still exists after delete returned.';
        catch ME
            lastMsg = ME.message;
        end
        pause(0.5*i);
        drawnow;
    end

    logFcn(['Warning: downloaded file remains: ' fileName]);
    logFcn(['Reason: ' lastMsg]);
    removed = false;
end
% ==========================================================
function out = local_ensure_natural_earth_map(cfg,varargin)
%LOCAL_ENSURE_NATURAL_EARTH_MAP Explicitly install optional 10 m polygons.
%
% This helper is not called during GUI initialization. When invoked
% explicitly, downloaded 10 m files are written to a user-writable folder.
% Packaged 50 m map files are read from the compiled archive (ctfroot) or
% the development SAGEhydrology/maps folder. No writes are attempted below
% ctfroot.

    logFcn = @(x) fprintf('%s\n',char(string(x)));

    if nargin >= 2 && ~isempty(varargin{1})
        ui = varargin{1};
        if isstruct(ui) && isfield(ui,'logFcn') && ~isempty(ui.logFcn)
            logFcn = ui.logFcn;
        end
    end

    writableRoot = local_writable_map_root(cfg);
    packagedRoots = local_packaged_map_roots(cfg);

    dir10 = fullfile(writableRoot,'10');
    name10 = 'ne_10m_admin_0_countries';
    name50 = 'ne_50m_admin_0_countries';
    file10 = fullfile(dir10,[name10 '.shp']);

    if local_shape_complete(dir10,name10)
        out = local_map_result(true,file10,'10m',false,'');
        return
    end

    [file50,dir50] = local_find_packaged_shape( ...
        packagedRoots,name50);

    local_mkdir(writableRoot);
    local_mkdir(dir10);

    url = ['https://naciscdn.org/naturalearth/' ...
        '10m/cultural/ne_10m_admin_0_countries.zip'];
    zipFile = fullfile(tempdir,'ne_10m_admin_0_countries.zip');
    unzipDir = fullfile(tempdir,'SAGE_ne_10m_admin_0_countries');

    try
        if isfile(zipFile), delete(zipFile); end
        if isfolder(unzipDir), rmdir(unzipDir,'s'); end
        mkdir(unzipDir);

        logFcn(['Natural Earth 10m map not found; ' ...
            'downloading it once ...']);
        opts = weboptions('Timeout',60,'ContentType','binary');
        websave(zipFile,url,opts);
        unzip(zipFile,unzipDir);

        D = dir(fullfile(unzipDir,'**',[name10 '.*']));
        if isempty(D)
            error('data_helpers:MapArchiveMissingFiles', ...
                'The Natural Earth archive did not contain %s.*.',name10);
        end

        for i = 1:numel(D)
            copyfile(fullfile(D(i).folder,D(i).name), ...
                fullfile(dir10,D(i).name),'f');
        end

        if ~local_shape_complete(dir10,name10)
            error('data_helpers:MapInstallIncomplete', ...
                ['Natural Earth 10m installation is missing ' ...
                 'one or more required shapefile components.']);
        end

        logFcn(['Natural Earth 10m map installed in: ' dir10]);
        out = local_map_result(true,file10,'10m',true,'');

    catch ME
        try
            local_remove_partial_shape(dir10,name10);
        catch
        end

        if ~isempty(file50) && local_shape_complete(dir50,name50)
            logFcn(['Natural Earth 10m download failed; ' ...
                'using packaged 50m map.']);
            out = local_map_result(true,file50,'50m',false,ME.message);
        else
            out = local_map_result(false,'','none',false,ME.message);
        end
    end

    try
        if isfile(zipFile), delete(zipFile); end
    catch
    end
    try
        if isfolder(unzipDir), rmdir(unzipDir,'s'); end
    catch
    end
end

function out = local_natural_earth_map_status(cfg)
%LOCAL_NATURAL_EARTH_MAP_STATUS Resolve installed map without downloading.

    writableRoot = local_writable_map_root(cfg);
    dir10 = fullfile(writableRoot,'10');
    name10 = 'ne_10m_admin_0_countries';
    if local_shape_complete(dir10,name10)
        out = local_map_result(true,fullfile(dir10, ...
            [name10 '.shp']),'10m',false,'');
        return
    end

    name50 = 'ne_50m_admin_0_countries';
    [file50,dir50] = local_find_packaged_shape( ...
        local_packaged_map_roots(cfg),name50);
    if ~isempty(file50) ...
            && local_shape_complete(dir50,name50)
        out = local_map_result(true,file50,'50m',false,'');
    else
        out = local_map_result(false,'','none',false, ...
            'No Natural Earth map is installed.');
    end
end

function mapRoot = local_writable_map_root(cfg)
%LOCAL_WRITABLE_MAP_ROOT Return folder used for downloaded map files.

    try
        mapRoot = fullfile(prefdir,'SAGE','maps');
    catch
        root = local_cfg_root(cfg);
        mapRoot = fullfile(root,'SAGE_maps');
    end
end

function roots = local_packaged_map_roots(cfg)
%LOCAL_PACKAGED_MAP_ROOTS Return read-only/development map candidates.

    roots = strings(0,1);
    if isdeployed
        try
            roots(end+1) = string(fullfile(ctfroot,'maps'));
            roots(end+1) = string(fullfile( ...
                ctfroot,'SAGEhydrology','maps'));
        catch
        end
    end

    if isstruct(cfg) ...
            && isfield(cfg,'SAGEhydro') ...
            && ~isempty(cfg.SAGEhydro)
        roots(end+1) = string(fullfile( ...
            char(string(cfg.SAGEhydro)),'maps'));
    end

    try
        root = local_cfg_root(cfg);
        roots(end+1) = string(fullfile( ...
            root,'SAGEhydrology','maps'));
    catch
    end

    try
        here = fileparts(mfilename('fullpath'));
        roots(end+1) = string(fullfile(here,'..','..','maps'));
    catch
    end

    roots = unique(roots(strlength(roots)>0),'stable');
end

function [file,folder] = local_find_packaged_shape(roots,name)
%LOCAL_FIND_PACKAGED_SHAPE Find a complete shapefile in candidate roots.

    file = '';
    folder = '';
    for i = 1:numel(roots)
        root = char(roots(i));
        candidates = {fullfile(root,'50'),root};
        for j = 1:numel(candidates)
            if local_shape_complete(candidates{j},name)
                folder = candidates{j};
                file = fullfile(folder,[name '.shp']);
                return
            end
        end
    end
end

% ========================================
function tf = local_shape_complete(d,name)
% ========================================

    required = {'.shp','.shx','.dbf'};

    tf = isfolder(d);

    for i = 1:numel(required)
        tf = tf && isfile( ...
            fullfile(d,[name required{i}]));
    end
end

% ==============================================
function local_remove_partial_shape(d,name)
% ==============================================

    if ~isfolder(d)
        return
    end

    D = dir(fullfile(d,[name '.*']));

    for i = 1:numel(D)
        try
            delete(fullfile( ...
                D(i).folder,D(i).name));
        catch
        end
    end
end

% =======================================
function out = local_map_result( ...
        ok,file,scale,downloaded,message)
% =======================================

    out = struct();
    out.ok = logical(ok);
    out.file = char(string(file));
    out.scale = char(string(scale));
    out.downloaded = logical(downloaded);
    out.message = char(string(message));
end

% ==============================================
function out = local_download_ush(cfg,stream,ui)
%LOCAL_DOWNLOAD_USH Install observed CAMELSH-US data (3,166 basins).
% The non-observed 5,842-basin archive is intentionally excluded because
% SAGE calibration requires discharge observations.

    if ~strcmp(stream,'hourly')
        out = struct('ok',false, ...
            'reason','CAMELSH-US is hourly only'); 
        return
    end
    logFcn=local_ui_log(ui); dirD = local_cfg_dirD(cfg,'CAMELSH_US');
    if ~local_ui_confirm(ui,'Install CAMELSH-US', ...
            ['CAMELSH-US observed data require about 21 GB compressed ' ...
             'and considerably more space after extraction.' newline ...
             'Install attributes and the 3,166 observed-basin files?'])
        out = struct('ok',false, ...
            'canceled',true); 
        return
    end
    local_mkdir(dirD); stage = fullfile(dirD,'_camelsh_install');
    local_mkdir(stage); d = [];
    try
        d=local_progress_dialog(ui,'Installing CAMELSH-US', ...
            'Looking for existing downloads ...');
        dl=fullfile(getenv('USERPROFILE'),'Downloads');
        srcAttr=fullfile(dl,'attributes','attributes');
        srcTS=fullfile(dl,'timeseries','Data','CAMELSH','timeseries');
        srcInfo=fullfile(dl,'info.csv');
        if ~(isfolder(srcAttr) ...
                && local_count_files(srcAttr,{'*.csv'})>=28)
            a = local_ush_archive(dl,stage,'attributes.7z', ...
                'https://zenodo.org/records/15066778/files/attributes.7z?download=1', ...
                'CAMELSH-US attributes',d,logFcn);
            ex = fullfile(stage,'attributes_extract'); 
            local_7z_extract(a,ex);
            srcAttr = local_find_ush_dir(ex, ...
                'attributes_nldas2_climate.csv');
        end
        if ~(isfolder(srcTS) ...
                && local_count_files(srcTS,{'*.nc'})>=3166)
            a = local_ush_archive(dl,stage,'timeseries.7z', ...
                'https://zenodo.org/records/15066778/files/timeseries.7z?download=1', ...
                'CAMELSH-US observed time series',d,logFcn);
            ex = fullfile(stage,'timeseries_extract'); 
            local_7z_extract(a,ex);
            srcTS = local_find_ush_dir(ex,'*.nc');
        end
        if ~isfile(srcInfo)
            srcInfo=fullfile(stage,'info.csv');
            local_download_file_retry( ...
                'https://zenodo.org/records/15066778/files/info.csv?download=1', ...
                srcInfo,'CAMELSH-US availability table', ...
                d,@()false,[],logFcn);
        end
        dstAttr=fullfile(dirD,'attributes'); 
        dstTS = fullfile(dirD,'hourly','timeseries');
        local_mkdir(dstAttr); local_mkdir(dstTS);
        logFcn('Copying CAMELSH-US attributes ...'); 
        copyfile(fullfile(srcAttr,'*'),dstAttr,'f');
        logFcn('Copying 3,166 observed CAMELSH-US NetCDF files ...'); 
        copyfile(fullfile(srcTS,'*.nc'),dstTS,'f');
        copyfile(srcInfo,fullfile(dirD,'info.csv'),'f');
        ok=local_count_files(dstTS,{'*.nc'})>=3166 ...
            && isfile(fullfile(dstAttr, ...
            'attributes_nldas2_climate.csv'));
        if ~ok
            error('data_helpers:USHInstallIncomplete', ...
                'CAMELSH-US installation is incomplete.'); 
        end
        try
            rmdir(stage,'s'); 
        catch
        end
        for nm = {'attributes.7z','timeseries.7z'}
            f = fullfile(dl,nm{1}); 
            if isfile(f)
                try
                    delete(f); 
                catch
                end
            end
        end
        if ~isempty(d) ...
                && isvalid(d)
            close(d); 
        end
        out = struct('ok',true,'region', ...
            'CAMELSH_US','dirD',dirD, ...
            'dirM',dstTS,'dirQ',dstTS, ...
            'stream','hourly');
    catch ME
        if ~isempty(d) ...
                && isvalid(d)
            close(d);
        end
        logFcn(['CAMELSH-US install failed: ' ME.message]);
        out = struct('ok',false,'reason',ME.message);
    end
end

function a = local_ush_archive(downloads,stage,name,url,label,d,logFcn)
    a=fullfile(downloads,name);
    if ~isfile(a), a=fullfile(stage,name); local_download_file_retry(url,a,label,d,@()false,[],logFcn); end
end
function local_7z_extract(archive,dst)
    local_mkdir(dst); c={fullfile(getenv('ProgramFiles'),'7-Zip','7z.exe'), ...
        fullfile(getenv('ProgramFiles'),'AMD','CIM','Bin64','7z.exe'), ...
        fullfile(getenv('ProgramFiles'),'AMD','AMDInstallManager','7z.exe')}; exe='';
    for i=1:numel(c), if isfile(c{i}), exe=c{i}; break; end, end
    if isempty(exe), error('data_helpers:Missing7zip','7-Zip is required to extract CAMELSH-US archives.'); end
    cmd=sprintf('"%s" x -y -o"%s" "%s"',exe,dst,archive); [s,msg]=system(cmd);
    if s~=0, error('data_helpers:USHExtract','7-Zip extraction failed: %s',msg); end
end
function d=local_find_ush_dir(root,pattern)
    q=dir(fullfile(root,'**',pattern)); q=q(~[q.isdir]);
    if isempty(q), error('data_helpers:USHPayload','Could not find %s below %s.',pattern,root); end
    d=q(1).folder;
end
