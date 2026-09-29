function R = region_config_PT()
%REGION_CONFIG_PT EStreams-SNIRH Portugal daily configuration.

    R = struct();
    R.code = 'EStreams_PT';
    R.acronym = 'PT';
    R.name = 'Portugal';
    R.dataset = 'EStreams-SNIRH Portugal';
    R.data_root = 'EStreams_PT';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'PT_%d_basins.txt';

    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 280;
    X.basins.training = 224;
    X.basins.evaluation = 56;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2005';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/1990';
    X.period.manual.eval_end = '30/09/2005';
    X.period.common_start = '01/10/1990';
    X.period.common_end = '30/09/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    X.meteo.product.items = {'EStreams-SNIRH Portugal [daily]'};
    X.meteo.product.default = X.meteo.product.items{1};
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 E-OBS 0.1-degree precipitation'};
    X.meteo.precipitation.default = X.meteo.precipitation.items{1};
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = { ...
        '1 E-OBS 0.1-degree mean air temperature'};
    X.meteo.temperature.default = X.meteo.temperature.items{1};
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = { ...
        '1 Supplied Hargreaves potential evapotranspiration', ...
        '2 Supplied Makkink potential evapotranspiration', ...
        '3 Supplied Penman-Monteith potential evapotranspiration'};
    X.meteo.pet.default = X.meteo.pet.items{1};
    X.meteo.pet.enabled = true;
    C = struct('name','Daily basin time series','type','folder', ...
        'path',fullfile('daily','timeseries'),'pattern','EStreams_PT_*.csv', ...
        'minimum_count',280,'required',true);
    X.checks = C;
    R.by_resolution.Daily = X;

    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();
end

function S = local_meteo_schema()
    S.name = 'daily EStreams-SNIRH Portugal';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'EStreams_PT_{gauge}.csv';
    S.timeline.reference = datetime(1950,1,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd';
    S.variables.P = local_variable('p_mean','mm/day','mm/day');
    S.variables.T = local_variable('t_mean','degC','degC');
    S.variables.Tmin = local_variable('t_min','degC','degC');
    S.variables.Tmax = local_variable('t_max','degC','degC');
    S.variables.Q = local_variable('q_cms_obs','m3/s','mm/day');
    S.variables.Q.area_normalize = true;
    S.variables.Q.invalid_lt = 0;
    S.aux.tables.gauges.file = '../../estreams_gauging_stations_PT.csv';
    S.aux.tables.gauges.key = 'basin_id';
    S.aux.tables.gauges.lat = 'lat';
    S.aux.tables.gauges.elev = 'elevation';
    S.aux.tables.gauges.area = 'area_estreams';
    S.aux.tables.gauges.area_scale = 1e6;
    pet = {'pet_mean','petmk_mean','petpm_mean'};
    Base = S;
    for k=1:3
        D = Base;
        D.variables.Ep = local_variable(pet{k},'mm/day','mm/day');
        D.variables.Ep.clip_min = 0;
        name = sprintf('pet%d',k);
        S.profiles.(name).match = struct('dt',1,'pet',k);
        S.profiles.(name).schema = D;
    end
end

function S = local_attribute_schema()
    files = {'estreams_gauging_stations_PT.csv', ...
        'estreams_topography_attributes_PT.csv', ...
        'estreams_soil_attributes_PT.csv', ...
        'estreams_geology_attributes_PT.csv', ...
        'estreams_geologycontinental_attributes_PT.csv', ...
        'estreams_hydrology_attributes_PT.csv', ...
        'estreams_vegetation_attributes_PT.csv', ...
        'estreams_snowcover_attributes_PT.csv', ...
        'estreams_meteorology_density_PT.csv'};
    S.name = 'EStreams-SNIRH Portugal attributes';
    S.tables = repmat(struct('file','','keys',{{'basin_id'}}, ...
        'make_valid_names',true,'key_type','char','required',true), ...
        numel(files),1);
    for k=1:numel(files), S.tables(k).file=files{k}; end
    S.id.column = 'basin_id';
    S.id.uppercase = true;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {};
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_transform = '';
    S.metadata.name_override_file = 'station_names_PT.csv';
    S.aliases.target = 'area_km2';
    S.aliases.sources = {'area_estreams'};
    S.aliases.required = true;
    S.aliases.default = NaN;
    S.region = 'EStreams_PT';
    S.zone.region = 'PT';
    S.progress.label = '... Reading Portugal generic attributes';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
