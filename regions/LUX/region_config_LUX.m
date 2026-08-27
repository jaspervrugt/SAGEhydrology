function R = region_config_LUX()
%REGION_CONFIG_LUX Regional defaults for LamaH-Lux.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_LUX';
    R.acronym = 'LUX';
    R.name = 'Luxembourg';
    R.dataset = 'LamaH-Lux';
    R.data_root = 'CAMELS_LUX';
    R.resolutions = {'Daily','Hourly','15 minutes'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'LUX_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 50;
    X.basins.file = 'LUX_50_basins.txt';
    X.basins.training = 40;
    X.basins.evaluation = 10;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2013';
    X.period.manual.train_end = '30/09/2019';
    X.period.manual.eval_start = '01/10/2008';
    X.period.manual.eval_end = '30/09/2013';
    X.period.common_start = '01/10/2008';
    X.period.common_end = '30/09/2019';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-LUX [daily]'};
    X.meteo.product.default = 'CAMELS-LUX [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 Radar', ...
        '2 Station', ...
        '3 ERA5 total', ...
        '4 Radar minimum', ...
        '5 Radar maximum' ...
        };
    X.meteo.precipitation.default = '1 Radar';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 ERA5 2 m', ...
        '2 Station' ...
        };
    X.meteo.temperature.default = '1 ERA5 2 m';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 Penman-Monteith', ...
        '2 Oudin et al.', ...
        '0 Zero PET' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = ...
        'CAMELS_LUX_hydromet_timeseries__daily_ID_*.csv';
    C(1).minimum_count = 50;
    X.checks = C;
    R.by_resolution.Daily = X;

    X = struct();
    X.label = 'Hourly';
    X.dt = 24;
    X.basins.universe = 51;
    X.basins.file = 'LUX_51_basins.txt';
    X.basins.training = 40;
    X.basins.evaluation = 11;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2014';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/2008';
    X.period.manual.eval_end = '30/09/2014';
    X.period.common_start = '01/10/2008';
    X.period.common_end = '30/09/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('hourly','timeseries');
    X.paths.discharge = fullfile('hourly','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-LUX [hourly]'};
    X.meteo.product.default = 'CAMELS-LUX [hourly]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 Radar', ...
        '2 Station', ...
        '3 ERA5 total', ...
        '4 Radar minimum', ...
        '5 Radar maximum' ...
        };
    X.meteo.precipitation.default = '1 Radar';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 ERA5 2 m', ...
        '2 Station' ...
        };
    X.meteo.temperature.default = '1 ERA5 2 m';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 Penman-Monteith', ...
        '2 Oudin et al.', ...
        '0 Zero PET' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('hourly','timeseries');
    C(1).pattern = ...
        'CAMELS_LUX_hydromet_timeseries_hourly_ID_*.csv';
    C(1).minimum_count = 51;
    X.checks = C;
    R.by_resolution.Hourly = X;

    X = struct();
    X.label = '15 minutes';
    X.dt = 96;
    X.basins.universe = 51;
    X.basins.file = 'LUX_51_basins.txt';
    X.basins.training = 40;
    X.basins.evaluation = 11;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2014';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/2008';
    X.period.manual.eval_end = '30/09/2014';
    X.period.common_start = '01/10/2008';
    X.period.common_end = '30/09/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('15min','timeseries');
    X.paths.discharge = fullfile('15min','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-LUX [15 minutes]'};
    X.meteo.product.default = 'CAMELS-LUX [15 minutes]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 Radar', ...
        '2 Station', ...
        '3 ERA5 total', ...
        '4 Radar minimum', ...
        '5 Radar maximum' ...
        };
    X.meteo.precipitation.default = '1 Radar';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 ERA5 2 m', ...
        '2 Station' ...
        };
    X.meteo.temperature.default = '1 ERA5 2 m';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 Penman-Monteith', ...
        '2 Oudin et al.', ...
        '0 Zero PET' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('15min','timeseries');
    C(1).pattern = ...
        'CAMELS_LUX_hydromet_timeseries_15min_ID_*.csv';
    C(1).minimum_count = 51;
    X.checks = C;
    R.by_resolution.min15 = X;


    % Declarative generic-reader schemas owned by this region module.
    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_attribute_schema()
    S.name = 'CAMELS-LUX attributes';
    files = { ...
        'CAMELS_LUX_meta_attributes.csv', ...
        'CAMELS_LUX_climatic_attributes.csv', ...
        'CAMELS_LUX_topographic_attributes.csv', ...
        'CAMELS_LUX_geologic_attributes.csv', ...
        'CAMELS_LUX_landuse_attributes.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id'};
    end
    S.id.column = 'gauge_id';
    S.id.uppercase = true;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'["'']',''; '\.0+$',''; '^ID_?',''};
    S.id.numeric_canonical = true;
    S.metadata.name_sources = {'Station','stream'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_LUX';
    S.zone.region = 'LUX';
    S.progress.label = '... Reading CAMELS-LUX generic attributes';
end

function S = local_meteo_schema()
    S.name = 'CAMELS-LUX supplied PM PET';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'column';
    S.time.column = 'Date';
    S.variables.P = local_variable('RR_rad','mm/step','mm/step');
    S.variables.Ep = local_variable('PET_PM','mm/step','mm/step');
    S.variables.T = local_variable('t2m','degC','degC');
    S.variables.Q = local_variable('Qspec','mm/step','mm/step');
    S.profiles.daily.match.dt = 1;
    S.profiles.daily.schema.file.pattern = 'CAMELS_LUX_hydromet_timeseries__daily_{gauge}.csv';
    S.profiles.daily.schema.timeline.step = days(1);
    S.profiles.daily.schema.time.input_format = 'yyyy-MM-dd';
    S.profiles.hourly.match.dt = 24;
    S.profiles.hourly.schema.file.pattern = 'CAMELS_LUX_hydromet_timeseries_hourly_{gauge}.csv';
    S.profiles.hourly.schema.timeline.step = hours(1);
    S.profiles.hourly.schema.time.input_format = 'yyyy-MM-dd HH:mm:ss';
    S.profiles.quarter_hour.match.dt = 96;
    S.profiles.quarter_hour.schema.file.pattern = 'CAMELS_LUX_hydromet_timeseries_15min_{gauge}.csv';
    S.profiles.quarter_hour.schema.timeline.step = minutes(15);
    S.profiles.quarter_hour.schema.time.input_format = 'yyyy-MM-dd HH:mm:ss';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
