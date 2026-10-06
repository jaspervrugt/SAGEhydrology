function R = region_config_HYD()
%REGION_CONFIG_HYD Regional defaults for HYD-RESPONSES Switzerland.
%
% The published archive is normalized to the standard SAGE layout under
% Data/HYD_RESPONSES. Descriptor tables sit at the regional root and the
% daily basin files sit in daily/timeseries.

    R = struct();
    R.code = 'HYD_RESPONSES';
    R.acronym = 'HYD';
    R.name = 'Switzerland';
    R.dataset = 'HYD-RESPONSES';
    R.data_root = 'HYD_RESPONSES';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'HYD_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 184;
    X.basins.training = 150;
    X.basins.evaluation = 34;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/01/2007';
    X.period.manual.train_end = '31/12/2014';
    X.period.manual.eval_start = '01/01/2015';
    X.period.manual.eval_end = '31/12/2022';
    X.period.common_start = '01/01/1992';
    X.period.common_end = '31/12/2022';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = X.paths.meteo;
    X.meteo.product.items = {'HYD-RESPONSES'};
    X.meteo.product.default = 'HYD-RESPONSES';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 MeteoSwiss RhiresD', ...
        '2 ERA5-Land total precipitation'};
    X.meteo.precipitation.default = '1 MeteoSwiss RhiresD';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 MeteoSwiss TabsD', ...
        '2 ERA5-Land 2-m temperature', ...
        '3 MeteoSwiss (TminD+TmaxD)/2'};
    X.meteo.temperature.default = '1 MeteoSwiss TabsD';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = {'1 ERA5-Land potential evaporation'};
    X.meteo.pet.default = '1 ERA5-Land potential evaporation';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Daily base-variable time series';
    C(1).type = 'folder';
    C(1).path = X.paths.meteo;
    C(1).pattern = 'HYDRESPONSES_basevars_catchment_*.csv';
    C(1).minimum_count = 184;
    C(2).name = 'Catchment descriptors';
    C(2).type = 'folder';
    C(2).path = '';
    C(2).pattern = 'HYDRESPONSES_descriptive_*.csv';
    C(2).minimum_count = 6;
    X.checks = C;
    R.by_resolution.Daily = X;

    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.variables.Q.area_normalize = true;
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();
end

function S = local_attribute_schema()
    S.name = 'HYD-RESPONSES catchment descriptors';
    S.root_candidates = {'.'};
    template = struct('file','','keys',{{'catchmentnr'}}, ...
        'make_valid_names',true,'make_unique_names',true, ...
        'join_type','left','delimiter',';', ...
        'encoding','UTF-8','comment_style','', ...
        'keep_columns',{{}},'row_filter',struct(), ...
        'prefix_columns','', ...
        'fill_missing_zero_prefixes',{{}});
    S.tables = repmat(template,6,1);
    S.tables(1).file = ...
        'HYDRESPONSES_descriptive_general_station_information.csv';
    S.tables(1).boolean_text = 'auto';
    S.tables(2).file = ...
        'HYDRESPONSES_descriptive_climatological_info.csv';
    S.tables(2).prefix_columns = 'clim_';
    S.tables(3).file = ...
        'HYDRESPONSES_descriptive_hydro_geological_terrestrial_info.csv';
    S.tables(3).prefix_columns = 'phys_';
    S.tables(3).fill_missing_zero_prefixes = { ...
        'REGION_','SOILDEPTH_','NUTRIENTCONT_', ...
        'SKELETALCONT_','WATERLOGGING_','PERMEABILITY_', ...
        'WATERSTORCAP_','ROCKGROUP_','AQUIFER_', ...
        'GENESE_','HYDROGEOL_'};
    S.tables(4).file = 'HYDRESPONSES_descriptive_Q347.csv';
    S.tables(4).prefix_columns = 'lowflow_';
    S.tables(5).file = 'HYDRESPONSES_descriptive_baseflowindex.csv';
    S.tables(5).row_filter = struct('column','time_scale', ...
        'value','yearly');
    S.tables(5).keep_columns = {'catchmentnr','Q','Q_baseflow','Q_BFI'};
    S.tables(5).prefix_columns = 'baseflow_';
    S.tables(6).file = 'HYDRESPONSES_descriptive_delayedflowindex.csv';
    S.tables(6).row_filter = struct('column','time_scale', ...
        'value','yearly');
    S.tables(6).keep_columns = {'catchmentnr','bp_1','bp_2','bp_3', ...
        'contr_1','contr_2','contr_3','contr_4','bias', ...
        'bfi_bp1_max','bfi_bp2_max','bfi_bp3_max','bfi_bp4_max','bfi_min'};
    S.tables(6).prefix_columns = 'delay_';
    S.id.column = 'catchmentnr';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.id.pad_width = 4;
    S.id.sort = 'numeric';
    S.metadata.name_components = {'water_name','place'};
    S.metadata.name_separator = ' at ';
    S.metadata.name_transform = '';
    S.metadata.standardize = false;
    S.aliases = repmat(struct('target','','sources',{{}}, ...
        'required',false,'default',NaN),4,1);
    S.aliases(1).target = 'area';
    S.aliases(1).sources = {'FOEN_area'};
    S.aliases(2).target = 'gauge_elev';
    S.aliases(2).sources = {'FOEN_stat_H'};
    S.aliases(3).target = 'p_mean';
    S.aliases(3).sources = {'clim_P_yearly_sum'};
    S.aliases(4).target = 'pet_mean';
    S.aliases(4).sources = {'clim_pev_yearly_mean'};
    S.projection.epsg = 21781;
    S.projection.x = 'x_coord';
    S.projection.y = 'y_coord';
    S.projection.latitude_target = 'gauge_lat';
    S.projection.longitude_target = 'gauge_lon';
    S.selection.default_file = '';
    S.selection.available_pattern = fullfile('daily','timeseries', ...
        'HYDRESPONSES_basevars_catchment_*.csv');
    S.selection.available_regex = 'catchment_(\d+)\.csv$';
    S.region = 'HYD_RESPONSES';
    S.zone.region = 'CH';
    S.progress.label = '... Reading HYD-RESPONSES attributes';
end

function S = local_meteo_schema()
    S.name = 'daily HYD-RESPONSES';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.id.pad_width = 4;
    S.timeline.reference = datetime(1950,1,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd';
    S.file.pattern = 'HYDRESPONSES_basevars_catchment_{gauge}.csv';
    S.file.delimiter = ';';
    S.missing_values = [-999 -99 -99.9 -99.99];
    S.variables.Q = local_variable('Q','m3/s','mm/day');
    S.variables.Q.area_normalize = true;
    S.variables.Q.source_description = [ ...
        'FOEN daily discharge with curated gauge-specific unit ' ...
        'corrections from the SAGE hydrologic data-action registry'];
    S.variables.P = local_variable('RhiresD_mean','mm/day','mm/day');
    S.variables.P.clip_min = 0;
    S.variables.Ep = local_variable('pev_mean_mean','mm/day','mm/day');
    S.variables.Ep.transform = @(x,context) abs(x);
    S.variables.Ep.clip_min = 0;
    S.variables.Ep.fill_missing = 0;
    S.variables.T = local_variable('TabsD_mean','degC','degC');
    S.variables.SWE = local_variable('SWECLQMD_mean','mm','mm');
    S.variables.SWE.source_description = 'SPASS snow-water equivalent';
    S.variables.SWE.clip_min = 0;
    S.variables.SM1 = local_variable( ...
        'swvl1_mean_mean','m3/m3','m3/m3');
    S.variables.SM2 = local_variable( ...
        'swvl2_mean_mean','m3/m3','m3/m3');
    S.variables.SM3 = local_variable( ...
        'swvl3_mean_mean','m3/m3','m3/m3');
    S.variables.SM = struct();
    S.variables.SM.derive = 'weighted_sum';
    S.variables.SM.inputs = {'SM1','SM2','SM3'};
    S.variables.SM.weights = [70 210 720];
    S.variables.SM.units = 'mm';
    S.variables.SM.target_units = 'mm';
    S.variables.SM.clip_min = 0;
    S.variables.SM.source_description = [ ...
        'ERA5-Land upper-1-m root-zone soil-water storage, derived ' ...
        'from volumetric layers 0-7, 7-28, and 28-100 cm'];
    Base = S;
    for precip = 1:2
        for temp = 1:3
            D = Base;
            if precip == 2
                D.variables.P = local_variable('tp_mean_mean', ...
                    'mm/day','mm/day');
                D.variables.P.clip_min = 0;
            end
            if temp == 2
                D.variables.T = local_variable('t2m_mean_mean', ...
                    'degC','degC');
            elseif temp == 3
                D.variables = rmfield(D.variables,'T');
                D.variables.Tmin = local_variable('TminD_mean', ...
                    'degC','degC');
                D.variables.Tmax = local_variable('TmaxD_mean', ...
                    'degC','degC');
                D.variables.T.source = '';
                D.variables.T.derive = 'mean_tmin_tmax';
            end
            name = sprintf('p%d_t%d',precip,temp);
            S.profiles.(name).match = struct( ...
                'dt',1,'precip',precip,'temp',temp,'pet',[0 1]);
            S.profiles.(name).schema = D;
        end
    end
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
