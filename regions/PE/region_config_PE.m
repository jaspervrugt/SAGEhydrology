function R = region_config_PE()
%REGION_CONFIG_PE Regional defaults for CAMELS-PE.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_PE';
    R.acronym = 'PE';
    R.name = 'Peru';
    R.dataset = 'CAMELS-PE';
    R.data_root = 'CAMELS_PE';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'PE_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 57;
    X.basins.training = 45;
    X.basins.evaluation = 12;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2006';
    X.period.manual.train_end = '30/09/2016';
    X.period.manual.eval_start = '01/10/1997';
    X.period.manual.eval_end = '30/09/2006';
    X.period.common_start = '01/10/1997';
    X.period.common_end = '30/09/2016';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-PE [daily]'};
    X.meteo.product.default = 'CAMELS-PE [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = { ...
        '1 Mean temperature', ...
        '2 (Tmin+Tmax)/2' ...
        };
    X.meteo.temperature.default = '1 Mean temperature';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = {'1 Supplied PET'};
    X.meteo.pet.default = '1 Supplied PET';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = 'PE_*.csv';
    C(1).minimum_count = 136;
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
    S.name = 'CAMELS-PE attributes';
    files = { ...
        'stations.csv', ...
        'topographic_attributes.csv', ...
        'climatic_indices.csv', ...
        'geologic_attributes.csv', ...
        'soil_attributes.csv', ...
        'landcover_attributes.csv', ...
        'human_intervention_attributes.csv', ...
        'hydrological_signatures.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id'};
    end
    S.id.column = 'gauge_id';
    S.id.uppercase = true;
    S.id.strip = true;
    S.id.regex = {'["'']',''; '^PE_',''};
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_PE';
    S.zone.region = 'PE';
    S.progress.label = '... Reading CAMELS-PE generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-PE (default products)';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'PE_{gauge}.csv';
    S.file.delimiter = ',';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.column = 'date';
    S.variables.P = local_variable('prec','mm/day','mm/day');
    S.variables.Ep = local_variable('pet','mm/day','mm/day');
    S.variables.T = local_variable('tmean','degC','degC');
    S.variables.Q = local_variable('flow_obs','mm/day','mm/day');
    S.variables.P.clip_min = 0;
    S.variables.Ep.clip_min = 0;
    S.progress.label = '... Reading daily CAMELS-PE generic data';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
