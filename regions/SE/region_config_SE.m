function R = region_config_SE()
%REGION_CONFIG_SE Regional defaults for CAMELS-SE.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_SE';
    R.acronym = 'SE';
    R.name = 'Sweden';
    R.dataset = 'CAMELS-SE';
    R.data_root = 'CAMELS_SE';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'SE_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 50;
    X.basins.training = 40;
    X.basins.evaluation = 10;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2000';
    X.period.manual.train_end = '30/09/2010';
    X.period.manual.eval_start = '01/10/1990';
    X.period.manual.eval_end = '30/09/2000';
    X.period.common_start = '01/10/1990';
    X.period.common_end = '30/09/2010';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily');
    X.paths.discharge = fullfile('daily');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-SE [daily]'};
    X.meteo.product.default = 'CAMELS-SE [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Observed'};
    X.meteo.precipitation.default = '1 Observed';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 Observed air'};
    X.meteo.temperature.default = '1 Observed air';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 Oudin et al.'};
    X.meteo.pet.default = '1 Oudin et al.';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily');
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
    S.name = 'CAMELS-SE attributes';
    S.tables = repmat(struct('file','','pattern','', ...
        'keys',{{'ID'}},'make_valid_names',true, ...
        'valid_name_replacement','underscore', ...
        'drop_empty_columns',true,'unique_key',true, ...
        'required',true,'prefix_from_file',false, ...
        'prefix_columns','','column_renames',{{}}),4,1);
    S.tables(1).file = 'catchments_physical_properties.csv';
    S.tables(2).file = 'catchments_landcover.csv';
    S.tables(2).column_renames = { ...
        'Water_percentage','Water_percentage_T'};
    S.tables(3).file = 'catchments_soil_classes.csv';
    S.tables(3).column_renames = { ...
        'Water_percentage','Water_percentage_Tsoil'};
    S.tables(4).pattern = ...
        'catchments_hydrological_signatures*.csv';
    S.tables(4).required = false;
    S.tables(4).prefix_from_file = true;
    S.id.column = 'ID';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.id.numeric_canonical = true;
    S.id.sort = 'numeric';
    S.id.output_type = 'double';
    S.metadata.name_sources = {'Name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_SE';
    S.zone.region = 'SE';
    S.progress.label = '... Reading CAMELS-SE generic attributes';
end

function S = local_meteo_schema()

    S.name = 'daily CAMELS-SE';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    
    % Basin files have names such as:
    % catchment_id_123_....csv
    S.file.pattern = 'catchment_id_{gauge}_*.csv';
    S.file.delimiter = ',';
    
    S.id.pad_width = 0;
    
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    
    S.time.mode = 'ymd_columns';
    S.time.year_column = 'Year';
    S.time.month_column = 'Month';
    S.time.day_column = 'Day';
    
    S.variables.P = local_variable( ...
        'Pobs_mm','mm/day','mm/day');
    
    S.variables.T = local_variable( ...
        'Tobs_C','degC','degC');
    
    S.variables.Q = local_variable( ...
        'Qobs_mm','mm/day','mm/day');
    
    % PET is derived from T, date, and basin latitude.
    S.variables.Ep.derive = 'oudin';
    S.variables.Ep.units = 'mm/day';
    S.variables.Ep.target_units = 'mm/day';
    
    S.progress.label = ...
        '... Reading daily CAMELS-SE generic data';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
