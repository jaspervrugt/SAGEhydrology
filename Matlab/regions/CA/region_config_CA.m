function R = region_config_CA()
%REGION_CONFIG_CA Regional defaults for CAMELS-SPAT.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_CA';
    R.acronym = 'CA';
    R.name = 'Canada';
    R.dataset = 'CAMELS-SPAT';
    R.data_root = 'CAMELS_CA';
    R.resolutions = {'Daily','Hourly'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'CA_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 764;
    X.basins.training = 600;
    X.basins.evaluation = 164;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2009';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/1998';
    X.period.manual.eval_end = '30/09/2009';
    X.period.common_start = '01/10/1998';
    X.period.common_end = '30/09/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','forcing','daymet');
    X.paths.discharge = fullfile('daily','discharge');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-SPAT Daymet [daily]'};
    X.meteo.product.default = 'CAMELS-SPAT Daymet [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Daymet'};
    X.meteo.precipitation.default = '1 Daymet';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 Daymet (Tmin+Tmax)/2'};
    X.meteo.temperature.default = '1 Daymet (Tmin+Tmax)/2';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 Daymet Priestley-Taylor'};
    X.meteo.pet.default = '1 Daymet Priestley-Taylor';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'CAMELS-SPAT metadata';
    C(1).type = 'pattern';
    C(1).path = '';
    C(1).pattern = 'camels-spat-metadata.csv';
    C(1).minimum_count = 1;
    C(2).name = 'CAMELS-SPAT attributes';
    C(2).type = 'pattern';
    C(2).path = '';
    C(2).pattern = 'attributes-lumped.csv';
    C(2).minimum_count = 1;
    C(3).name = 'Daily Daymet forcing';
    C(3).type = 'pattern';
    C(3).path = fullfile('daily','forcing','daymet');
    C(3).pattern = 'CAN_*_daymet_lumped.nc';
    C(3).minimum_count = 764;
    C(4).name = 'Daily discharge';
    C(4).type = 'pattern';
    C(4).path = fullfile('daily','discharge');
    C(4).pattern = 'CAN_*_daily_flow_observations.nc';
    C(4).minimum_count = 764;
    X.checks = C;
    R.by_resolution.Daily = X;

    X = struct();
    X.label = 'Hourly';
    X.dt = 24;
    X.basins.universe = 764;
    X.basins.training = 600;
    X.basins.evaluation = 164;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2008';
    X.period.manual.train_end = '30/09/2018';
    X.period.manual.eval_start = '01/10/1998';
    X.period.manual.eval_end = '30/09/2008';
    X.period.common_start = '01/10/1998';
    X.period.common_end = '30/09/2018';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('hourly','forcing','rdrs');
    X.paths.discharge = fullfile('hourly','discharge');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-SPAT RDRS [hourly]'};
    X.meteo.product.default = 'CAMELS-SPAT RDRS [hourly]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 RDRS'};
    X.meteo.precipitation.default = '1 RDRS';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 RDRS air temperature'};
    X.meteo.temperature.default = '1 RDRS air temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 RDRS Penman-Monteith'};
    X.meteo.pet.default = '1 RDRS Penman-Monteith';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'CAMELS-SPAT metadata';
    C(1).type = 'pattern';
    C(1).path = '';
    C(1).pattern = 'camels-spat-metadata.csv';
    C(1).minimum_count = 1;
    C(2).name = 'CAMELS-SPAT attributes';
    C(2).type = 'pattern';
    C(2).path = '';
    C(2).pattern = 'attributes-lumped.csv';
    C(2).minimum_count = 1;
    C(3).name = 'Hourly RDRS forcing';
    C(3).type = 'pattern';
    C(3).path = fullfile('hourly','forcing','rdrs');
    C(3).pattern = 'CAN_*_rdrs_lumped.nc';
    C(3).minimum_count = 764;
    C(4).name = 'Hourly discharge';
    C(4).type = 'pattern';
    C(4).path = fullfile('hourly','discharge');
    C(4).pattern = 'CAN_*_hourly_flow_observations.nc';
    C(4).minimum_count = 764;
    X.checks = C;
    R.by_resolution.Hourly = X;

    % Declarative generic-reader schemas owned by this region module.
    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = local_discharge_schema();
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_attribute_schema()
    S.name = 'CAMELS-SPAT Canada attributes';
    S.tables = repmat(struct('file','','keys',{{'basin_id'}}, ...
        'layout','','make_valid_names',true, ...
        'make_unique_names',true,'column_renames',{{}}, ...
        'keep_columns',{{}},'row_filter',struct(), ...
        'attribute_column','','metadata_columns',0, ...
        'catalog_region','','exclude_catalog',{{}}, ...
        'derived',struct()),2,1);
    S.tables(1).file = 'camels-spat-metadata.csv';
    S.tables(1).keys = {'basin_id'};
    S.tables(1).column_renames = {'Station_id','basin_id'};
    S.tables(1).keep_columns = {'basin_id', ...
        'Station_name','Station_lat','Station_lon','Basin_area_km2'};
    S.tables(1).row_filter.column = 'Country';
    S.tables(1).row_filter.value = 'CAN';
    S.tables(2).file = 'attributes-lumped.csv';
    S.tables(2).layout = 'catalog_attribute_matrix';
    S.tables(2).make_valid_names = false;
    S.tables(2).make_unique_names = false;
    S.tables(2).attribute_column = 'Attribute';
    S.tables(2).metadata_columns = 4;
    S.tables(2).catalog_region = 'CAMELS_CA';
    S.tables(2).exclude_catalog = ...
        {'Station_lat','Station_lon','Basin_area_km2'};
    S.tables(2).derived.target = 'forest_fraction';
    S.tables(2).derived.operation = 'row_sum';
    S.tables(2).derived.rows = { ...
        'lc2_evergreen_needleleaf_forest_fraction', ...
        'lc2_evergreen_broadleaf_forest_fraction', ...
        'lc2_deciduous_needleleaf_forest_fraction', ...
        'lc2_deciduous_broadleaf_forest_fraction', ...
        'lc2_mixed_forest_fraction'};
    S.id.column = 'basin_id';
    S.id.uppercase = true;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'^CAN[-_]?',''};
    S.id.sort = 'text';
    S.metadata.name_sources = {'Station_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_CA';
    S.zone.region = 'CA';
    S.progress.label = '... Reading CAMELS-CA generic attributes';
end

function S = local_meteo_schema()
    S.name = 'CAMELS-SPAT Canada';
    S.format = 'netcdf';
    S.layout = 'one_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'variable';
    S.time.column = 'time';
    S.profiles.daily.match.dt = 1;
    S.profiles.daily.schema.file.pattern = 'CAN_{gauge}_daymet_lumped.nc';
    S.profiles.daily.schema.timeline.step = days(1);
    S.profiles.daily.schema.time.units = 'days';
    S.profiles.daily.schema.time.origin = datetime(1980,1,1,12,0,0);
    S.profiles.daily.schema.time.snap = 'day';
    S.profiles.daily.schema.variables.P = local_variable('prcp','mm/day','mm/day');
    S.profiles.daily.schema.variables.Ep = local_variable('pet','mm/day','mm/day');
    S.profiles.daily.schema.variables.Tmin = local_variable('tmin','degC','degC');
    S.profiles.daily.schema.variables.Tmax = local_variable('tmax','degC','degC');
    S.profiles.daily.schema.variables.T.derive = 'mean_Tmin_Tmax';
    S.profiles.hourly.match.dt = 24;
    S.profiles.hourly.schema.file.pattern = 'CAN_{gauge}_rdrs_lumped.nc';
    S.profiles.hourly.schema.timeline.step = hours(1);
    S.profiles.hourly.schema.time.units = 'hours';
    S.profiles.hourly.schema.time.origin = datetime(1980,1,1,9,0,0);
    S.profiles.hourly.schema.time.snap = 'hour';
    S.profiles.hourly.schema.variables.P = local_variable( ...
        'RDRS_v2.1_A_PR0_SFC','kg/m2/s','mm/hour');
    S.profiles.hourly.schema.variables.P.scale = 3600;
    S.profiles.hourly.schema.variables.Ep = local_variable('pet','kg/m2/s','mm/hour');
    S.profiles.hourly.schema.variables.Ep.scale = 3600;
    S.profiles.hourly.schema.variables.Ep.clip_min = 0;
    S.profiles.hourly.schema.variables.T = local_variable( ...
        'RDRS_v2.1_P_TT_1.5m','K','degC');
    S.profiles.hourly.schema.variables.T.offset = -273.15;
    S.aux.tables.stations.file = '../../../camels-spat-metadata.csv';
    S.aux.tables.stations.key = 'Station_id';
    S.aux.tables.stations.lat = 'Station_lat';
    S.aux.tables.stations.area = 'Basin_area_km2';
    S.aux.tables.stations.area_scale = 1e6;
    S.aux.tables.elevation.file = '../../../attributes-lumped.csv';
    S.aux.tables.elevation.layout = 'wide_attributes';
    S.aux.tables.elevation.attribute_column = 'Attribute';
    S.aux.tables.elevation.column_prefix = 'CAN_';
    S.aux.tables.elevation.strip_prefix = 'CAN_';
    S.aux.tables.elevation.elev_row = 'elev_mean';
end

function S = local_discharge_schema()
    S.name = 'CAMELS-CA discharge'; 
    S.format = 'netcdf';
    S.layout = 'one_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'variable'; 
    S.time.column = 'time';
    S.time.units = 'minutes'; 
    S.time.origin = datetime(1950,1,1);
    S.variables.Q = local_variable('q_obs','m3/s','mm/step');
    S.variables.Q.valid_min = 0; 
    S.variables.Q.area_normalize = true;
    S.variables.Q.quality.valid_min = 0;
    S.variables.Q.quality.valid_max = 0;
    S.profiles.daily.match.dt = 1;
    S.profiles.daily.schema.file.pattern =  ...
        'CAN_{gauge}_daily_flow_observations.nc';
    S.profiles.daily.schema.timeline.step = days(1);
    S.profiles.daily.schema.time.snap = 'day';
    S.profiles.daily.schema.variables.Q.quality.sources =  ...
        {'q_obs_is_ice_affected','q_obs_is_partial_day'};
    S.profiles.hourly.match.dt = 24;
    S.profiles.hourly.schema.file.pattern =  ...
        'CAN_{gauge}_hourly_flow_observations.nc';
    S.profiles.hourly.schema.timeline.step = hours(1);
    S.profiles.hourly.schema.time.snap = 'hour';
    S.profiles.hourly.schema.variables.Q.quality.sources =  ...
        {'q_obs_is_ice_affected','q_obs_is_below_sensor_level'};
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
