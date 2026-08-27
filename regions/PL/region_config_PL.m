function R = region_config_PL()
%REGION_CONFIG_PL Regional defaults for CAMELS-PL.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_PL';
    R.acronym = 'PL';
    R.name = 'Poland';
    R.dataset = 'CAMELS-PL';
    R.data_root = 'CAMELS_PL';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'PL_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 354;
    X.basins.training = 300;
    X.basins.evaluation = 54;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2005';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/1990';
    X.period.manual.eval_end = '30/09/2005';
    X.period.common_start = '01/10/1990';
    X.period.common_end = '30/09/2005';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-PL [daily]'};
    X.meteo.product.default = 'CAMELS-PL [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 Temperature'};
    X.meteo.temperature.default = '1 Temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 Hargreaves'};
    X.meteo.pet.default = '1 Hargreaves';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = '*';
    C(1).minimum_count = 1;
    X.checks = C;
    R.by_resolution.Daily = X;


    % Declarative generic-reader schemas owned by this region module.
    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_attribute_schema()
    files = { ...
        'CAMELS_PL_topographic_attributes.csv', ...
        'CAMELS_PL_climatic_attributes.csv', ...
        'CAMELS_PL_hydrologic_attributes.csv', ...
        'CAMELS_PL_landcover_attributes.csv', ...
        'CAMELS_PL_soil_attributes.csv', ...
        'CAMELS_PL_simulation_benchmark.csv', ...
        'CAMELS_PL_BDOT10K_land_cover_catchments.csv'};
    S.name = 'CAMELS-PL attributes';
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'encoding','UTF-8','make_valid_names',true, ...
        'make_unique_names',true,'key_type','char', ...
        'numeric_text','auto','numeric_text_exceptions',{{}}, ...
        'duplicate_policy','suffix','duplicate_suffix','_pl', ...
        'required',true,'prefix_columns',''),numel(files),1);
    exceptions = {'gauge_id','gauge_name','water_body_name', ...
        'high_prec_timing','low_prec_timing', ...
        'flow_period_start','flow_period_end'};
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id'};
        S.tables(i).numeric_text_exceptions = exceptions;
    end
    S.tables(7).required = false;
    S.tables(7).prefix_columns = 'bdot_';
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''; '[^0-9]',''};
    S.id.sort = 'numeric';
    S.metadata.name_sources = {'gauge_name','water_body_name'};
    S.metadata.name_transform = '';
    S.aliases = repmat(struct('target','','sources',{{}}, ...
        'required',true,'default',NaN),3,1);
    S.aliases(1).target = 'gauge_lat_dd';
    S.aliases(1).sources = {'gauge_lon'};
    S.aliases(2).target = 'gauge_lon_dd';
    S.aliases(2).sources = {'gauge_lat'};
    S.aliases(3).target = 'area_km2';
    S.aliases(3).sources = {'area_metadata'};
    S.region = 'CAMELS_PL';
    S.zone.region = 'PL';
    S.progress.label = '... Reading CAMELS-PL generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-PL Hargreaves';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'CAMELS_PL_hydromet_timeseries_{gauge}.csv';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd';
    S.variables.P = local_variable('precipitation_mean','mm/day','mm/day');
    S.variables.T = local_variable('temperature_mean','degC','degC');
    S.variables.Tmin = local_variable('minimum_temperature_mean','degC','degC');
    S.variables.Tmax = local_variable('maximum_temperature_mean','degC','degC');
    S.variables.Ep.derive = 'hargreaves';
    S.variables.Ep.latitude_fallback = 52;
    S.variables.Ep.fill_missing = 0;
    S.variables.Q = local_variable('discharge_spec_obs','mm/day','mm/day');
    S.aux.tables.gauges.file = '../../gauge_information_PL.csv';
    S.aux.tables.gauges.key = 'gauge_id';
    S.aux.tables.gauges.lat = 'gauge_lat';
    S.aux.tables.gauges.elev = 'gauge_elev';
    S.aux.tables.gauges.area = 'area_km2';
    S.aux.tables.gauges.area_scale = 1e6;
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
