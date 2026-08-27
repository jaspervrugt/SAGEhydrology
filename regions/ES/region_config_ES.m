function R = region_config_ES()
%REGION_CONFIG_ES Regional defaults for CAMELS-ES.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_ES';
    R.acronym = 'ES';
    R.name = 'Spain';
    R.dataset = 'CAMELS-ES';
    R.data_root = 'CAMELS_ES';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'ES_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 269;
    X.basins.training = 220;
    X.basins.evaluation = 49;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2006';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/1992';
    X.period.manual.eval_end = '30/09/2006';
    X.period.common_start = '01/10/1992';
    X.period.common_end = '30/09/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = { ...
        '1 EMO-1 [daily]', ...
        '2 ERA5-Land [daily]' ...
        };
    X.meteo.product.default = '1 EMO-1 [daily]';
    X.meteo.product.enabled = true;
    X.meteo.precipitation.items = {'1 EMO-1'};
    X.meteo.precipitation.default = '1 EMO-1';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 EMO-1'};
    X.meteo.temperature.default = '1 EMO-1';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 EMO-1'};
    X.meteo.pet.default = '1 EMO-1';
    X.meteo.pet.enabled = false;
    X.meteo.linked.source = 'product';
    X.meteo.linked.labels = { ...
        '1 EMO-1', ...
        '2 ERA5-Land' ...
        };
    X.meteo.linked.targets = { ...
        'precipitation', ...
        'temperature', ...
        'pet' ...
        };
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = '*';
    C(1).minimum_count = 1;
    C(2).name = 'Discharge data';
    C(2).type = 'folder';
    C(2).path = fullfile('daily','timeseries');
    C(2).pattern = '*';
    C(2).minimum_count = 1;
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
    S.name = 'CAMELS-ES attributes';
    files = { ...
        'attributes_other_camelses.csv', ...
        'attributes_caravan_camelses.csv', ...
        'atributes_efas_hydrometeorology_camelses.csv', ...
        'atributes_efas_static_maps_camelses.csv', ...
        'attributes_efas_model_parameters_camelses.csv', ...
        'attributes_hydroatlas_camelses.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true,'duplicate_policy','suffix', ...
        'duplicate_suffix','_es'),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id'};
    end
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = true;
    S.id.strip = true;
    S.id.regex = {'^camelses[_-]*',''; '\.0+$',''; '[^0-9]',''};
    S.id.sort = 'numeric';
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_ES';
    S.zone.region = 'ES';
    S.progress.label = '... Reading CAMELS-ES generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-ES';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'camelses_{gauge}.csv';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd';
    S.variables.Q = local_variable('streamflow','mm/day','mm/day');
    S.aux.tables.gauge.file = '../../gauge_information.txt';
    S.aux.tables.gauge.key = 'gauge_id';
    S.aux.tables.gauge.lat = 'lat';
    S.aux.tables.gauge.area = 'area';
    S.aux.tables.gauge.area_scale = 1e6;
    A.variables.P = local_variable('pr_emo1','mm/day','mm/day');
    A.variables.T = local_variable('ta_emo1','degC','degC');
    A.variables.Ep = local_variable('e0_emo1','mm/day','mm/day');
    S.profiles.emo1.match = struct('dt',1,'data',1);
    S.profiles.emo1.schema.variables = A.variables;
    B.variables.P = local_variable('total_precipitation_sum','mm/day','mm/day');
    B.variables.T = local_variable('temperature_2m_mean','degC','degC');
    B.variables.Ep = local_variable('potential_evaporation_sum','mm/day','mm/day');
    S.profiles.era5land.match = struct('dt',1,'data',2);
    S.profiles.era5land.schema.variables = B.variables;
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
