function R = region_config_GB()
%REGION_CONFIG_GB Regional defaults for CAMELS-GB.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_GB';
    R.acronym = 'GB';
    R.name = 'Great Britain';
    R.dataset = 'CAMELS-GB';
    R.data_root = 'CAMELS_GB';
    R.resolutions = {'Daily','Hourly'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'GB_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 669;
    X.basins.training = 500;
    X.basins.evaluation = 169;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2002';
    X.period.manual.train_end = '30/09/2015';
    X.period.manual.eval_start = '01/10/1990';
    X.period.manual.eval_end = '30/09/2002';
    X.period.common_start = '01/10/1990';
    X.period.common_end = '30/09/2015';
    X.paths.run_root = fullfile('daily');
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-GB [daily]'};
    X.meteo.product.default = 'CAMELS-GB [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 Temperature'};
    X.meteo.temperature.default = '1 Temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = { ...
        '1 PET', ...
        '2 PET with interception' ...
        };
    X.meteo.pet.default = '1 PET';
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

    X = struct();
    X.label = 'Hourly';
    X.dt = 24;
    X.basins.universe = 669;
    X.basins.training = 500;
    X.basins.evaluation = 169;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2002';
    X.period.manual.train_end = '30/09/2015';
    X.period.manual.eval_start = '01/10/1991';
    X.period.manual.eval_end = '30/09/2002';
    X.period.common_start = '01/10/1991';
    X.period.common_end = '30/09/2015';
    X.paths.run_root = fullfile('hourly');
    X.paths.meteo = fullfile('hourly','timeseries');
    X.paths.discharge = fullfile('hourly','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-GB v2 [hourly]'};
    X.meteo.product.default = 'CAMELS-GB v2 [hourly]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 CEH-GEAR', ...
        '2 Grid-to-Grid' ...
        };
    X.meteo.precipitation.default = '1 CEH-GEAR';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = {'1 Daily temperature'};
    X.meteo.temperature.default = '1 Daily temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = { ...
        '1 Daily PET', ...
        '2 Daily PET with interception' ...
        };
    X.meteo.pet.default = '1 Daily PET';
    X.meteo.pet.enabled = true;
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
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_attribute_schema()
    S = local_attribute_release('daily','CAMELS_GB_',false);
    H = local_attribute_release('hourly','camels_gb_v2_',true);
    S.profiles.hourly.match.stream = 'hourly';
    S.profiles.hourly.schema = H;
end

function S = local_attribute_release(folder,prefix,useFallback)
    kinds = {'topographic','climatic','hydrologic', ...
        'hydrogeology','soil','landcover', ...
        'humaninfluence','hydrometry'};
    S.name = ['CAMELS-GB ' folder ' attributes'];
    S.tables = repmat(struct('file','','keys',{{'gauge_id'}}, ...
        'make_valid_names',true,'duplicate_policy','suffix', ...
        'duplicate_suffix','_gb','fallback_file',''),numel(kinds),1);
    for i = 1:numel(kinds)
        S.tables(i).file = [prefix kinds{i} '_attributes.csv'];
        if useFallback
            S.tables(i).fallback_file = fullfile('..','daily', ...
                ['CAMELS_GB_' kinds{i} '_attributes.csv']);
        end
    end
    S.root_candidates = {folder,''};
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.id.numeric_canonical = true;
    S.id.sort = 'numeric';
    S.id.output_type = 'double';
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_transform = '';
    S.aliases = local_gb_aliases(folder);
    S.region = 'CAMELS_GB';
    S.zone.region = 'GB';
    S.progress.label = ['... Reading CAMELS-GB ' ...
        folder ' generic attributes'];
end

function A = local_gb_aliases(stream)
    bases = {'dwood','ewood','grass','shrub', ...
        'crop','urban','inwater','bares'};
    years = [2022 2021 2020 2019 2018 2017 2015 1990];
    if strcmpi(stream,'daily')
        years = [1990 years(years ~= 1990)];
    end
    A = repmat(struct('target','','sources',{{}}, ...
        'required',false,'default',NaN,'type','numeric'),11,1);
    for i = 1:numel(bases)
        A(i).target = [bases{i} '_perc'];
        A(i).sources = cellstr(string(bases{i}) + ...
            "_perc_" + string(years));
    end
    prefix = 'daily';
    if strcmpi(stream,'hourly')
        prefix = 'hourly';
    end
    targets = {'flow_period_start','flow_period_end', ...
        'flow_perc_complete'};
    for i = 1:numel(targets)
        j = numel(bases) + i;
        A(j).target = targets{i};
        A(j).sources = {[prefix '_' targets{i}]};
        if i < 3
            A(j).type = 'preserve';
        end
    end
end

function S = local_meteo_schema()
    S.name = 'CAMELS-GB';
    S.format = 'csv';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.id.strip_decimal = true;
    Base = S;
    for pet = 1:2
        D = Base;
        D.layout = 'multi_file_per_basin';
        D.timeline.step = days(1);
        D.files.daily.pattern = 'CAMELS_GB_hydromet_timeseries_{gauge}_*.csv';
        D.files.daily.format = 'csv';
        D.files.daily.time.mode = 'column';
        D.files.daily.time.column = 'date';
        D.files.daily.time.input_format = 'yyyy-MM-dd';
        D.variables.P = local_variable('precipitation','mm/day','mm/day');
        if pet == 1
            petSource = 'pet';
        else
            petSource = 'peti';
        end
        D.variables.Ep = local_variable(petSource,'mm/day','mm/day');
        D.variables.Ep.clip_min = 0;
        D.variables.T = local_variable('temperature','degC','degC');
        D.variables.Q = local_variable('discharge_spec','mm/day','mm/day');
        names = fieldnames(D.variables);
        for i = 1:numel(names)
            D.variables.(names{i}).file = 'daily';
        end
        D.aux.tables.topo.file = '../CAMELS_GB_topographic_attributes.csv';
        D.aux.tables.topo.key = 'gauge_id';
        D.aux.tables.topo.lat = 'gauge_lat';
        D.aux.tables.topo.elev = 'elev_mean';
        D.aux.tables.topo.area = 'area';
        D.aux.tables.topo.area_scale = 1e6;
        name = sprintf('daily_pet%d',pet);
        S.profiles.(name).match = struct('dt',1,'precip',1,'temp',1,'pet',pet);
        S.profiles.(name).schema = D;
    end
    for precip = 1:2
        for pet = 1:2
            H = Base;
            H.layout = 'multi_file_per_basin';
            H.timeline.reference = datetime(1950,10,1,9,0,0);
            H.timeline.step = hours(1);
            H.files.hourly.pattern = 'camels_gb_v2_hydromet_hourly_timeseries_{gauge}_*.csv';
            H.files.hourly.format = 'csv';
            H.files.hourly.time.mode = 'column';
            H.files.hourly.time.column = 'date';
            H.files.hourly.time.input_format = 'yyyy-MM-dd HH:mm:ss';
            H.files.daily.pattern = '../../daily/timeseries/CAMELS_GB_hydromet_timeseries_{gauge}_*.csv';
            H.files.daily.format = 'csv';
            H.files.daily.time.mode = 'column';
            H.files.daily.time.column = 'date';
            H.files.daily.time.input_format = 'yyyy-MM-dd';
            H.files.daily.time.join = 'calendar_day';
            if precip == 1
                pSource = 'precipitation_cehgear';
            else
                pSource = 'precipitation_gradgb';
            end
            H.variables.P = local_variable(pSource,'mm/hour','mm/hour');
            H.variables.P.file = 'hourly';
            if pet == 1
                petSource = 'pet';
            else
                petSource = 'peti';
            end
            H.variables.Ep = local_variable(petSource,'mm/day','mm/hour');
            H.variables.Ep.scale = 1/24;
            H.variables.Ep.clip_min = 0;
            H.variables.Ep.file = 'daily';
            H.variables.T = local_variable('temperature','degC','degC');
            H.variables.T.file = 'daily';
            H.variables.Q = local_variable('discharge_spec','mm/hour','mm/hour');
            H.variables.Q.file = 'hourly';
            H.aux.tables.topo.file = '../../daily/CAMELS_GB_topographic_attributes.csv';
            H.aux.tables.topo.key = 'gauge_id';
            H.aux.tables.topo.lat = 'gauge_lat';
            H.aux.tables.topo.elev = 'elev_mean';
            H.aux.tables.topo.area = 'area';
            H.aux.tables.topo.area_scale = 1e6;
            name = sprintf('hourly_p%d_pet%d',precip,pet);
            S.profiles.(name).match = struct('dt',24,'precip',precip,'temp',1,'pet',pet);
            S.profiles.(name).schema = H;
        end
    end
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
