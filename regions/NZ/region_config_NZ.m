function R = region_config_NZ()
%REGION_CONFIG_NZ Regional defaults for CAMELS-NZ.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_NZ';
    R.acronym = 'NZ';
    R.name = 'New Zealand';
    R.dataset = 'CAMELS-NZ';
    R.data_root = 'CAMELS_NZ';
    R.resolutions = {'Daily','Hourly'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'NZ_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 341;
    X.basins.training = 300;
    X.basins.evaluation = 41;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2003';
    X.period.manual.train_end = '30/09/2013';
    X.period.manual.eval_start = '01/10/1993';
    X.period.manual.eval_end = '30/09/2003';
    X.period.common_start = '01/10/1993';
    X.period.common_end = '30/09/2013';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-NZ [daily]'};
    X.meteo.product.default = 'CAMELS-NZ [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 precipitation'};
    X.meteo.precipitation.default = '1 precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 air temperature'};
    X.meteo.temperature.default = '1 air temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 PET'};
    X.meteo.pet.default = '1 PET';
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

    X = struct();
    X.label = 'Hourly';
    X.dt = 24;
    X.basins.universe = 341;
    X.basins.training = 300;
    X.basins.evaluation = 41;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2003';
    X.period.manual.train_end = '30/09/2013';
    X.period.manual.eval_start = '01/10/1993';
    X.period.manual.eval_end = '30/09/2003';
    X.period.common_start = '01/10/1993';
    X.period.common_end = '30/09/2013';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('hourly','timeseries');
    X.paths.discharge = fullfile('hourly','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-NZ [hourly]'};
    X.meteo.product.default = 'CAMELS-NZ [hourly]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 precipitation'};
    X.meteo.precipitation.default = '1 precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 air temperature'};
    X.meteo.temperature.default = '1 air temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 PET'};
    X.meteo.pet.default = '1 PET';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('hourly','timeseries');
    C(1).pattern = '*';
    C(1).minimum_count = 1;
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
    S.name = 'CAMELS-NZ attributes';
    files = { ...
        '1.CAMELS_NZ_Catchment_information.csv', ...
        '2.CAMELS_NZ_Climatic_attribute.csv', ...
        '3.CAMELS_NZ_Landcover_attribute.csv', ...
        '4.CAMELS_NZ_Geology.csv', ...
        '5.CAMELS_NZ_Anthropogenic_attribute.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'delimiter',',','column_renames',{{}}, ...
        'boolean_text','auto','duplicate_policy','suffix', ...
        'duplicate_suffix','_nz'),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'station_id'};
    end
    S.tables(1).column_renames = { ...
        'Station_ID','station_id'; 'RID','rid'; ...
        'Station Name','station_name'; 'Latitude (WGS 84)','gauge_lat'; ...
        'Longitude(WGS 84)','gauge_lon'; 'uparea','area'; ...
        'Region','region_id'; 'UpStreamLakes','upstream_lakes'; ...
        'usLake','us_lake'; 'Stream_Order','stream_order'; ...
        'elevation','elev_mean'; 'usSteep','us_steep'; ...
        'usLowGrad','us_lowgrad'; 'usAveSlope','us_aveslope'; ...
        'DIST_SEA','dist_sea'; 'SRC_OF_FLW','src_of_flow'; ...
        'Records','records'};
    S.tables(2).column_renames = { ...
        'Station_ID','station_id'; 'RID','rid'; ...
        'StationName','station_name_clim'; 'latitude','latitude_clim'; ...
        'longitude','longitude_clim'; 'usDaysRainGT25','days_rain_gt25'; ...
        'usRainDays10','rain_days_10'; 'Mean Annual Rainfall','p_mean'; ...
        'usAnRainVar','p_var'; 'usPET','pet_mean'; ...
        'Climate Zone','climate_zone'};
    S.tables(3).column_renames = { ...
        'Station_ID','station_id'; 'RID','rid'; ...
        'StationName','station_name_land'; 'latitude','latitude_land'; ...
        'longitude','longitude_land'; 'LANDCOVER','landcover'};
    S.tables(4).column_renames = { ...
        'Station_ID','station_id'; 'RID','rid'; ...
        'StationName','station_name_geo'; 'latitude','latitude_geo'; ...
        'longitude','longitude_geo'; 'usParticleSize','particle_size'; ...
        'usHard','hardrock_index'; 'usCalc','carbonate_index'; ...
        'Geology','geology'};
    S.tables(5).column_renames = { ...
        'Station_ID','station_id'; 'RID','rid'; ...
        'StationName','station_name_ant'; 'latitude','latitude_ant'; ...
        'longitude','longitude_ant'; 'IsDam','is_dam'; ...
        'IsWeir','is_weir'; 'IsAbstracted','is_abstracted'; ...
        'IsEphemeral','is_ephemeral'; 'IsGWinfluenced','is_gw_influenced'; ...
        'IsRatingOk','is_rating_ok'; 'IsSnowInfluenced','is_snow_influenced'; ...
        'DSDamAffected','ds_dam_affected'; 'FlowRegime','flow_regime'; ...
        'Influence','influence'};
    S.id.column = 'station_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'["'']',''; '\.0+$',''};
    S.id.numeric_canonical = true;
    S.id.sort = 'numeric';
    S.metadata.name_sources = {'station_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_NZ';
    S.zone.region = 'NZ';
    S.progress.label = '... Reading CAMELS-NZ generic attributes';
end

function S = local_meteo_schema()
    S.name = 'CAMELS-NZ'; 
    S.format = 'csv'; 
    S.layout = 'multi_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'column'; 
    S.time.column = 'time';
    S.aux.tables.gauge.file = '../../gauge_information.txt';
    S.aux.tables.gauge.key = 'GAGE_ID'; 
    S.aux.tables.gauge.lat = 'LAT';
    S.aux.tables.gauge.area = 'DRAINAGE AREA (KM^2)';
    S.aux.tables.gauge.area_scale = 1e6;
    Base = S;
    D = Base;
    D.timeline.step = days(1); 
    D.time.input_format = 'yyyy-MM-dd';
    D.files.precip.pattern = 'precipitation/daily_precipitation_station_id_{gauge}.csv';
    D.files.temperature.pattern = 'temperature/daily_temperature_station_id_{gauge}.csv';
    D.files.humidity.pattern = 'relative_humidity/daily_RH_station_id_{gauge}.csv';
    D.variables.P = local_variable('precipitation','mm/day','mm/day'); 
    D.variables.P.file = 'precip';
    D.variables.TemperatureRaw = local_variable('temperature','','');
    D.variables.TemperatureRaw.file = 'temperature';
    D.variables.RH = local_variable('Relative_humidity','',''); 
    D.variables.RH.file = 'humidity';
    D.variables.T.derive = 'temperature_auto_celsius';
    D.variables.Ep.derive = 'oudin'; 
    D.variables.Ep.coefficient = 1/(100*2.45);
    D.variables.Ep.input_precision = 'double';
    S.profiles.daily.match.dt = 1; 
    S.profiles.daily.schema = D;
    H = Base;
    H.timeline.step = hours(1); 
    H.time.input_format = 'yyyy-MM-dd HH:mm:ss';
    H.files.precip.pattern = 'precipitation/precipitation_station_id_{gauge}.csv';
    H.files.temperature.pattern = 'temperature/temperature_station_id_{gauge}.csv';
    H.files.pet.pattern = 'pet/PET_station_id_{gauge}.csv';
    H.variables.P = local_variable('precipitation','mm/hour','mm/hour'); 
    H.variables.P.file = 'precip';
    H.variables.TemperatureRaw = local_variable('temperature','','');
    H.variables.TemperatureRaw.file = 'temperature';
    H.variables.T.derive = 'temperature_auto_celsius';
    H.variables.Ep = local_variable('PET','mm/hour','mm/hour'); 
    H.variables.Ep.file = 'pet';
    S.profiles.hourly.match.dt = 24; 
    S.profiles.hourly.schema = H;
end

function S = local_discharge_schema()
    S.name = 'CAMELS-NZ discharge';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'column';
    S.time.column = 'time';
    S.variables.Q = local_variable('flow','m3/s','mm/step');
    S.variables.Q.valid_min = 0;
    S.variables.Q.area_normalize = true;
    S.profiles.daily.match.dt = 1;
    S.profiles.daily.schema.file.pattern =  ...
        'streamflow/daily_flow_station_id_{gauge}.csv';
    S.profiles.daily.schema.timeline.step = days(1);
    S.profiles.daily.schema.time.input_format = 'yyyy-MM-dd';
    S.profiles.hourly.match.dt = 24;
    S.profiles.hourly.schema.file.pattern =  ...
        'streamflow/flow_station_id_{gauge}.csv';
    S.profiles.hourly.schema.timeline.step = hours(1);
    S.profiles.hourly.schema.time.input_format = 'yyyy-MM-dd HH:mm:ss';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
