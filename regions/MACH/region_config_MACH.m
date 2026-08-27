function R = region_config_MACH()
%REGION_CONFIG_MACH Regional defaults for the daily MACH-US dataset.

    R = struct();
    R.code = 'MACH_US';
    R.acronym = 'MACH';
    R.name = 'United States (MACH)';
    R.dataset = 'MACH';
    R.data_root = 'MACH_US';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'MACH_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 1014;
    X.basins.training = 800;
    X.basins.evaluation = 214;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2003';
    X.period.manual.train_end = '30/09/2023';
    X.period.manual.eval_start = '01/10/1983';
    X.period.manual.eval_end = '30/09/2003';
    % MACH begins on 1 January 1980. The first complete water-year model
    % window, including the 365-day spin-up, therefore starts scoring on
    % 1 October 1981 (spin-up starts on 1 October 1980).
    X.period.common_start = '01/10/1981';
    X.period.common_end = '30/09/2023';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    X.meteo.product.items = {'MACH Daymet V4 + GLEAM 4.2a [daily]'};
    X.meteo.product.default = X.meteo.product.items{1};
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Daymet V4 precipitation'};
    X.meteo.precipitation.default = X.meteo.precipitation.items{1};
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 Daymet V4 mean temperature'};
    X.meteo.temperature.default = X.meteo.temperature.items{1};
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 GLEAM 4.2a potential evaporation'};
    X.meteo.pet.default = X.meteo.pet.items{1};
    X.meteo.pet.enabled = false;

    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Daily forcing and streamflow';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = 'basin_*_MACH.csv';
    C(1).minimum_count = 1014;
    C(2).name = 'MACH attributes';
    C(2).type = 'file';
    C(2).path = fullfile('attributes','site_info.csv');
    C(3).name = 'Basin boundaries';
    C(3).type = 'file';
    C(3).path = 'MACH_basins_all.gpkg';
    C(4).name = 'Basin universe';
    C(4).type = 'file';
    C(4).path = 'MACH_1014_basins.txt';
    X.checks = C;
    R.by_resolution.Daily = X;

    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();
end

function S = local_meteo_schema()
    S.name = 'daily MACH-US';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'basin_{gauge}_MACH.csv';
    S.file.delimiter = ',';
    S.id.pad_width = 8;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'DATE';
    S.time.input_format = '';
    S.variables.P = local_variable('PRCP','mm/day','mm/day');
    S.variables.P.clip_min = 0;
    S.variables.Ep = local_variable('PET','mm/day','mm/day');
    S.variables.Ep.clip_min = 0;
    S.variables.T = local_variable('TAIR','degC','degC');
    S.variables.Q = local_variable('OBSQ','mm/day','mm/day');
    S.reuse_existing_q = true;
    S.progress.label = '... Reading daily MACH-US generic data';
end

function S = local_attribute_schema()
    files = {fullfile('attributes','site_info.csv'), ...
        fullfile('attributes','overall_climate.csv'), ...
        fullfile('attributes','soil.csv'), ...
        fullfile('attributes','geology.csv'), ...
        fullfile('attributes','hydrology.csv'), ...
        fullfile('attributes','anthropogenic.csv')};
    S.name = 'MACH-US attributes';
    S.tables = repmat(struct('file','','keys',{{'SITENO'}}, ...
        'duplicate_policy','drop'),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
    end
    S.id.column = 'SITENO';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.id.numeric_canonical = true;
    S.id.pad_width = 8;
    S.id.sort = 'numeric';
    S.id.output_type = 'string';
    S.metadata.name_components = {'station_name','state'};
    S.metadata.name_separator = ', ';
    S.metadata.name_sources = {'station_name'};
    S.metadata.name_transform = '';
    S.region = 'MACH_US';
    S.zone.region = 'US';
    definitions = { ...
        'gauge_lat',{'dec_lat_va'},true; ...
        'gauge_lon',{'dec_long_va'},true; ...
        'area',{'area_sqkm'},true; ...
        'p_mean',{'PPT9120','PPT8110'},false; ...
        'pet_mean',{'PET9120','PET8110'},false; ...
        'aridity',{'mach_aridity'},false};
    S.aliases = repmat(struct('target','','sources',{{}}, ...
        'required',false,'default',NaN),size(definitions,1),1);
    for i = 1:size(definitions,1)
        S.aliases(i).target = definitions{i,1};
        S.aliases(i).sources = definitions{i,2};
        S.aliases(i).required = definitions{i,3};
    end
    S.progress.label = '... Reading MACH-US generic attributes';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
