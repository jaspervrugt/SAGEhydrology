function R = region_config_FI()
%REGION_CONFIG_FI Regional defaults for CAMELS-FI.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_FI';
    R.acronym = 'FI';
    R.name = 'Finland';
    R.dataset = 'CAMELS-FI';
    R.data_root = 'CAMELS_FI';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'FI_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 320;
    X.basins.training = 270;
    X.basins.evaluation = 50;
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
    X.meteo.product.items = {'CAMELS-FI [daily]'};
    X.meteo.product.default = 'CAMELS-FI [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = { ...
        '1 Mean', ...
        '2 Minimum', ...
        '3 Maximum', ...
        '4 Grid minimum' ...
        };
    X.meteo.temperature.default = '1 Mean';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 native', ...
        '2 FMI', ...
        '3 Singer', ...
        '0 Zero PET' ...
        };
    X.meteo.pet.default = '1 native';
    X.meteo.pet.enabled = true;
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
    S.name = 'CAMELS-FI attributes';
    files = { ...
        'CAMELS_FI_meta_attributes.csv', ...
        'CAMELS_FI_climatic_attributes.csv', ...
        'CAMELS_FI_topographic_attributes.csv', ...
        'CAMELS_FI_soil_attributes.csv', ...
        'CAMELS_FI_geology_attributes.csv', ...
        'CAMELS_FI_landcover_attributes.csv', ...
        'CAMELS_FI_humaninfluence_attributes.csv', ...
        'CAMELS_FI_hydrologic_attributes.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true,'key_type','char', ...
        'drop_missing_key',true,'unique_key',true),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id'};
    end
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'["'']',''; '\.0+$',''};
    S.metadata.name_sources = {'gauge_name','basin_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_FI';
    S.zone.region = 'FI';
    S.progress.label = '... Reading CAMELS-FI generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-FI';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'CAMELS_FI_hydromet_timeseries_{gauge}_*.csv';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd';
    S.variables.P = local_variable('precipitation','mm/day','mm/day');
    S.variables.P.clip_min = 0;
    S.variables.P.fill_missing = 0;
    S.variables.Q = local_variable('discharge_spec','mm/day','mm/day');
    S.variables.Q.invalid_le = 0;
    S.aux.tables.meta.file = '../../CAMELS_FI_meta_attributes.csv';
    S.aux.tables.meta.key = 'gauge_id';
    S.aux.tables.meta.lat = 'gauge_lat';
    S.aux.tables.meta.area = 'area';
    S.aux.tables.meta.area_scale = 1e6;
    S.aux.tables.topo.file = '../../CAMELS_FI_topographic_attributes.csv';
    S.aux.tables.topo.key = 'gauge_id';
    S.aux.tables.topo.elev = 'elev_mean';
    tempPrefs = {{'temperature_mean','temperature_gmin'}, ...
        {'temperature_min','temperature_mean'}, ...
        {'temperature_max','temperature_mean'}, ...
        {'temperature_gmin','temperature_mean'}};
    petPrefs = {{'pet','pet_fmi','pet_singer'}, ...
        {'pet_fmi','pet','pet_singer'}, ...
        {'pet_singer','pet','pet_fmi'}};
    Base = S;
    for temp = 1:4
        for pet = 0:3
            D = Base;
            D.variables.T.sources = tempPrefs{temp};
            D.variables.T.source_operation = 'aliases';
            if pet == 0
                D.variables.Ep.derive = 'constant';
                D.variables.Ep.value = 0;
            else
                D.variables.Ep.sources = petPrefs{pet};
                D.variables.Ep.source_operation = 'aliases';
                D.variables.Ep.clip_min = 0;
                D.variables.Ep.fill_missing = 0;
            end
            name = sprintf('t%d_pet%d',temp,pet);
            S.profiles.(name).match = struct('dt',1,'temp',temp,'pet',pet);
            S.profiles.(name).schema = D;
        end
    end
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
