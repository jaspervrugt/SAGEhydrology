function R = region_config_COL()
%REGION_CONFIG_COL Regional defaults for CAMELS-COL.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_COL';
    R.acronym = 'COL';
    R.name = 'Colombia';
    R.dataset = 'CAMELS-COL';
    R.data_root = 'CAMELS_COL';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'COL_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    % Start from the complete local CAMELS-COL inventory. Data-quality
    % screening subsequently determines the active paired-period population.
    X.basins.universe = 346;
    X.basins.training = 280;
    X.basins.evaluation = 66;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/1990';
    X.period.manual.train_end = '30/09/1994';
    X.period.manual.eval_start = '01/10/1986';
    X.period.manual.eval_end = '30/09/1990';
    X.period.common_start = '01/10/1986';
    X.period.common_end = '30/09/1994';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-COL [daily]'};
    X.meteo.product.default = 'CAMELS-COL [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 (Tmin+Tmax)/2'};
    X.meteo.temperature.default = '1 (Tmin+Tmax)/2';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 ETP'};
    X.meteo.pet.default = '1 ETP';
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
    files = {'02_CAMELS_COL_Catchment_information.csv', ...
        '10_CAMELS_COL_Physiograpic_characteristics.csv', ...
        '08_CAMELS_COL_Climatic_indices.csv', ...
        '05_CAMELS_COL_Geologic_characteristics.csv', ...
        '06_CAMELS_COL_Land_cover_characteristics.csv', ...
        '07_CAMELS_COL_Soil_characteristics.csv', ...
        '11_CAMELS_COL_Land_use_capability.csv', ...
        '09_CAMELS_COL_Hydrological_signatures.csv', ...
        'CAMELS_COL_mean_temperature.csv'};
    S.name = 'CAMELS-COL attributes';
    S.tables = repmat(struct('file','','keys',{{'gauge_id'}}, ...
        'encoding','ISO-8859-1','make_valid_names',true, ...
        'make_unique_names',true,'drop_empty_rows',true, ...
        'numeric_text','threshold','numeric_text_ratio',0.9, ...
        'numeric_text_exceptions',{{}},'column_renames',{{}}, ...
        'duplicate_policy','suffix','duplicate_suffix','_col', ...
        'required',true,'keep_columns',{{}}),numel(files),1);
    exceptions = {'gauge_id','gauge_department','gauge_star','gauge_end'};
    commonRenames = {'gravelius_index_','gravelius_index'; ...
        'streng_chanel_','streng_chanel'; ...
        'Urban_area','urban_area_luc'; ...
        'No_classified','no_classified_luc'};
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).numeric_text_exceptions = exceptions;
        S.tables(i).column_renames = commonRenames;
    end
    S.tables(6).column_renames = [commonRenames; ...
        {'water_bodies_perc','soil_water_bodies_perc'}];
    S.tables(9).required = false;
    S.tables(9).encoding = '';
    S.tables(9).numeric_text = '';
    S.tables(9).keep_columns = {'gauge_id','temp_mean'};
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''; '[^0-9]',''};
    S.id.sort = 'numeric';
    S.metadata.name_sources = {'gauge_department'};
    S.metadata.name_fallback_prefix = 'COL_';
    S.metadata.name_transform = '';
    S.projection.epsg = 3857;
    S.projection.x = 'gauge_lon';
    S.projection.y = 'gauge_lat';
    S.projection.latitude_target = 'gauge_lat';
    S.projection.longitude_target = 'gauge_lon';
    S.region = 'CAMELS_COL';
    S.zone.region = 'COL';
    S.progress.label = '... Reading CAMELS-COL generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-COL';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'Hydromet_data_{gauge}.txt';
    S.file.delimiter = sprintf('\t');
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'Fecha';
    S.time.input_format = 'dd/MM/yyyy';
    S.variables.P.sources = {'Precipitacion','pr','precipitation'};
    S.variables.P.source_operation = 'aliases';
    S.variables.Ep.sources = {'ETP_','ETP','poten_evapo','pet'};
    S.variables.Ep.source_operation = 'aliases';
    S.variables.Tmin.sources = {'Temperatura_minima','t_min','tmin'};
    S.variables.Tmin.source_operation = 'aliases';
    S.variables.Tmax.sources = {'Temperatura_maxima','t_max','tmax'};
    S.variables.Tmax.source_operation = 'aliases';
    S.variables.T.derive = 'mean_tmin_tmax';
    S.variables.Q.sources = {'Caudal','streamflow','q'};
    S.variables.Q.source_operation = 'aliases';
    S.variables.Q.units = 'm3/s';
    S.variables.Q.target_units = 'mm/day';
    S.variables.Q.area_normalize = true;
    S.aux.tables.gauge.file = '../../gauge_information.txt';
    S.aux.tables.gauge.key = 'gauge_id';
    S.aux.tables.gauge.lat = 'gauge_lat';
    S.aux.tables.gauge.area = 'area_km2';
    S.aux.tables.gauge.area_scale = 1e6;
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
