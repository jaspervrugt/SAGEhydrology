function R = region_config_DK()
%REGION_CONFIG_DK Regional defaults for CAMELS-DK.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_DK';
    R.acronym = 'DK';
    R.name = 'Denmark';
    R.dataset = 'CAMELS-DK';
    R.data_root = 'CAMELS_DK';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'DK_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 304;
    X.basins.training = 250;
    X.basins.evaluation = 54;
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
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-DK [daily]'};
    X.meteo.product.default = 'CAMELS-DK [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 Mean temperature'};
    X.meteo.temperature.default = '1 Mean temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 Potential evapotranspiration'};
    X.meteo.pet.default = '1 Potential evapotranspiration';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = 'CAMELS_DK_obs_based_*.csv';
    C(1).minimum_count = 304;
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
    files = {'CAMELS_DK_topography.csv', ...
        'CAMELS_DK_climate.csv','CAMELS_DK_geology.csv', ...
        'CAMELS_DK_landuse.csv','CAMELS_DK_soil.csv', ...
        'CAMELS_DK_signature_obs_based.csv'};
    S.name = 'CAMELS-DK attributes';
    S.tables = repmat(struct('file','','keys',{{'catch_id'}}, ...
        'delimiter',',','make_valid_names',true, ...
        'make_unique_names',true,'numeric_text','threshold', ...
        'numeric_text_ratio',1, ...
        'numeric_text_exceptions',{{'catch_id'}}, ...
        'required',true),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
    end
    S.tables(6).required = false;
    S.root_candidates = {'','Attributes', ...
        'CAMELS_DK',fullfile('CAMELS_DK','Attributes')};
    S.id.column = 'catch_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.id.numeric_canonical = true;
    S.id.sort = 'numeric';
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_fallback_prefix = 'CAMELS-DK ';
    S.metadata.name_transform = '';
    S.metadata.standardize = false;
    S.aliases = repmat(struct('target','','sources',{{}}, ...
        'required',true,'default',NaN,'scale',[], ...
        'divisor',[]),2,1);
    S.aliases(1).target = 'area';
    S.aliases(1).sources = {'catch_area'};
    S.aliases(1).divisor = 1e6;
    S.aliases(2).target = 'gauge_elev';
    S.aliases(2).sources = {'elev_mean'};
    S.projection.epsg = 25832;
    S.projection.x = 'catch_outlet_lon';
    S.projection.y = 'catch_outlet_lat';
    S.projection.latitude_target = 'gauge_lat';
    S.projection.longitude_target = 'gauge_lon';
    S.selection.available_pattern = fullfile( ...
        'daily','timeseries','CAMELS_DK_obs_based_*.csv');
    S.selection.available_regex = ...
        '^CAMELS_DK_obs_based_(.+)\.csv$';
    S.region = 'CAMELS_DK';
    S.zone.region = 'DK';
    S.progress.label = '... Reading CAMELS-DK generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-DK (default products)';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'CAMELS_DK_obs_based_{gauge}.csv';
    S.file.delimiter = ',';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.column = 'time';
    S.variables.P = local_variable('precipitation','mm/day','mm/day');
    S.variables.Ep = local_variable('pet','mm/day','mm/day');
    S.variables.P.clip_min = 0;
    S.variables.Ep.clip_min = 0;
    S.variables.T = local_variable('temperature','degC','degC');
    S.variables.Q = local_variable('Qobs','m3/s','mm/day');
    S.variables.Q.area_normalize = true;
    S.progress.label = '... Reading daily CAMELS-DK generic data';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
