function R = region_config_BULL()
%REGION_CONFIG_BULL Regional defaults for the BULL database, Spain.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and complete water years.

    R = struct();
    R.code = 'BULL_ES';
    R.acronym = 'BULL';
    R.name = 'Spain (BULL)';
    R.dataset = 'BULL';
    R.data_root = 'BULL_ES';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'BULL_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 484;
    X.basins.training = 400;
    X.basins.evaluation = 84;
    X.period.spinup_days = 365;
    X.period.manual.eval_start = '01/10/1991';
    X.period.manual.eval_end = '30/09/2005';
    X.period.manual.train_start = '01/10/2005';
    X.period.manual.train_end = '30/09/2020';
    X.period.common_start = '01/10/1991';
    X.period.common_end = '30/09/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');

    X.meteo.product.items = {'BULL [daily]'};
    X.meteo.product.default = 'BULL [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 AEMET', ...
        '2 ERA5-Land', ...
        '3 EMO-1' ...
        };
    X.meteo.precipitation.default = '1 AEMET';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 AEMET (mean)', ...
        '2 ERA5-Land (mean)', ...
        '3 EMO-1 (mean)' ...
        };
    X.meteo.temperature.default = '1 AEMET (mean)';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 AEMET', ...
        '2 ERA5-Land', ...
        '3 EMO-1' ...
        };
    X.meteo.pet.default = '1 AEMET';
    X.meteo.pet.enabled = true;

    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    folders = {'streamflow','AEMET','ERA5_Land','EMO1_arc'};
    labels = {'Observed streamflow','AEMET forcing', ...
        'ERA5-Land forcing','EMO-1 forcing'};
    for i = 1:numel(folders)
        C(i).name = labels{i};
        C(i).type = 'folder';
        C(i).path = fullfile('daily','timeseries',folders{i});
        C(i).pattern = '*.csv';
        C(i).minimum_count = 484;
        C(i).required = true;
    end
    X.checks = C;
    R.by_resolution.Daily = X;

    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();
end

function S = local_attribute_schema()
    S.name = 'BULL attributes';
    files = { ...
        'attributes_other_.csv', ...
        'attributes_caravan_.csv', ...
        'attributes_hydroatlas_.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true,'duplicate_policy','suffix', ...
        'duplicate_suffix','_bull'),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id'};
    end
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = true;
    S.id.strip = true;
    S.id.regex = {'^bull[_-]*',''; '\.0+$',''};
    S.id.sort = 'numeric';
    S.metadata.name_sources = {};
    S.metadata.name_fallback_prefix = 'BULL_';
    S.metadata.name_transform = '';
    S.region = 'BULL_ES';
    S.zone.region = 'ES';
    S.progress.label = '... Reading BULL generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily BULL';
    S.format = 'csv';
    S.layout = 'multi_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd';
    S.progress.label = '... Reading daily BULL data';

    names = {'q','aemet','era5land','emo1'};
    patterns = { ...
        fullfile('streamflow','streamflow_{gauge}.csv'), ...
        fullfile('AEMET','AEMET_{gauge}.csv'), ...
        fullfile('ERA5_Land','ERA5_Land_{gauge}.csv'), ...
        fullfile('EMO1_arc','EMO1_{gauge}.csv')};
    for i = 1:numel(names)
        name = names{i};
        S.files.(name).pattern = patterns{i};
        S.files.(name).format = 'csv';
        S.files.(name).time.mode = 'column';
        S.files.(name).time.column = 'date';
        S.files.(name).time.input_format = 'yyyy-MM-dd';
    end

    % BULL follows the CARAVAN convention: streamflow is already
    % catchment-area-normalized runoff depth in mm/day.
    S.variables.Q = local_variable('streamflow','mm/day','mm/day');
    S.variables.Q.file = 'q';
    S.variables.Q.valid_min = 0;
    S.variables.Q.area_normalize = false;
    S.variables.P = local_variable( ...
        'total_precipitation','mm/day','mm/day');
    S.variables.P.file = 'aemet';
    S.variables.T = local_variable( ...
        'temperature_mean','degC','degC');
    S.variables.T.file = 'aemet';
    S.variables.Ep = local_variable( ...
        'potential_evapotranspiration','mm/day','mm/day');
    S.variables.Ep.file = 'aemet';
    S.variables.Ep.clip_min = 0;

    S.aux.tables.other.file = '../../attributes_other_.csv';
    S.aux.tables.other.key = 'gauge_id';
    S.aux.tables.other.lat = 'gauge_lat';
    S.aux.tables.other.area = 'area';
    S.aux.tables.other.area_scale = 1e6;
    S.aux.tables.other.strip_prefix = 'BULL_';
    S.aux.tables.hydro.file = '../../attributes_hydroatlas_.csv';
    S.aux.tables.hydro.key = 'gauge_id';
    S.aux.tables.hydro.elev = 'ele_mt_sav';
    S.aux.tables.hydro.strip_prefix = 'BULL_';

    sourceFiles = {'aemet','era5land','emo1'};
    for p = 1:3
        for t = 1:3
            for pet = 1:3
                D.variables.P = local_variable( ...
                    'total_precipitation','mm/day','mm/day');
                D.variables.P.file = sourceFiles{p};
                D.variables.T = local_variable( ...
                    'temperature_mean','degC','degC');
                D.variables.T.file = sourceFiles{t};
                D.variables.Ep = local_variable( ...
                    'potential_evapotranspiration','mm/day','mm/day');
                D.variables.Ep.file = sourceFiles{pet};
                D.variables.Ep.clip_min = 0;
                name = sprintf('p%d_t%d_pet%d',p,t,pet);
                S.profiles.(name).match = struct( ...
                    'dt',1,'precip',p,'temp',t,'pet',pet);
                S.profiles.(name).schema = D;
            end
        end
    end
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
