function R = region_config_CZ()
%REGION_CONFIG_CZ Regional defaults for CAMELS-CZ.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_CZ';
    R.acronym = 'CZ';
    R.name = 'Czechia';
    R.dataset = 'CAMELS-CZ';
    R.data_root = 'CAMELS_CZ';
    R.resolutions = {'Daily','Hourly'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'CZ_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 249;
    X.basins.training = 200;
    X.basins.evaluation = 49;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2005';
    X.period.manual.train_end = '30/09/2015';
    X.period.manual.eval_start = '01/10/1995';
    X.period.manual.eval_end = '30/09/2005';
    X.period.common_start = '01/10/1995';
    X.period.common_end = '30/09/2015';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-CZ [daily]'};
    X.meteo.product.default = 'CAMELS-CZ [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 CHMI', ...
        '2 ERA5-Land', ...
        '3 Combined CHMI/ERA5-Land' ...
        };
    X.meteo.precipitation.default = '1 CHMI';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 CHMI', ...
        '2 ERA5-Land', ...
        '3 Combined CHMI/ERA5-Land' ...
        };
    X.meteo.temperature.default = '1 CHMI';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 CHMI (Penman-Monteith)', ...
        '2 ERA5-Land (Penman-Monteith)', ...
        '3 ERA5-Land native' ...
        };
    X.meteo.pet.default = '1 CHMI (Penman-Monteith)';
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

    H = X;
    H.label = 'Hourly';
    H.dt = 24;
    H.period.spinup_days = 365;
    H.period.manual.eval_start = '01/10/2017';
    H.period.manual.eval_end = '30/09/2021';
    H.period.manual.train_start = '01/10/2021';
    H.period.manual.train_end = '30/09/2024';
    H.period.common_start = '01/10/2017';
    H.period.common_end = '30/09/2024';
    H.paths.meteo = fullfile('hourly','timeseries');
    H.paths.discharge = fullfile('hourly','timeseries');
    H.meteo.product.items = {'CAMELS-CZ [hourly]'};
    H.meteo.product.default = 'CAMELS-CZ [hourly]';
    H.meteo.precipitation.items = { ...
        '1 CHMI', ...
        '2 ERA5-Land' ...
        };
    H.meteo.precipitation.default = '1 CHMI';
    H.meteo.temperature.items = { ...
        '1 CHMI', ...
        '2 ERA5-Land' ...
        };
    H.meteo.temperature.default = '1 CHMI';
    H.meteo.pet.items = { ...
        '1 ERA5-Land native (hourly)', ...
        '2 CHMI Penman-Monteith (daily/24)', ...
        '3 ERA5-Land Penman-Monteith (daily/24)', ...
        '4 ERA5-Land native (daily/24)' ...
        };
    H.meteo.pet.default = '1 ERA5-Land native (hourly)';
    C(1).name = 'Hourly hydrometeorological data';
    C(1).path = fullfile('hourly','timeseries');
    C(1).pattern = 'camelscz_*.csv';
    C(1).minimum_count = 249;
    C(2) = C(1);
    C(2).name = 'Daily PET data';
    C(2).path = fullfile('daily','timeseries');
    C(2).minimum_count = 249;
    C(2).required = false;
    H.checks = C;
    R.by_resolution.Hourly = H;


    % Declarative generic-reader schemas owned by this region module.
    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_attribute_schema()
    S.name = 'CAMELS-CZ attributes';
    files = { ...
        'attributes_other_camelscz.csv', ...
        'attributes_caravan_camelscz.csv', ...
        'attributes_hydroatlas_camelscz.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true,'exclude_columns_regex',''),3,1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id'};
    end
    S.tables(3).exclude_columns_regex = '_s(0[1-9]|1[0-2])$';
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = true;
    S.id.strip = true;
    S.id.regex = {'["'']',''};
    S.id.optional_prefix = 'camelscz_';
    S.id.output_regex = {'^camelscz_',''};
    S.id.output_uppercase = true;
    S.id.output_lowercase = false;
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_CZ';
    S.zone.region = 'CZ';
    S.progress.label = '... Reading CAMELS-CZ generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-CZ';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = '{gauge}.csv';
    S.id.prefix = 'camelscz_';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd';
    S.variables.Q = local_variable('streamflow','mm/day','mm/day');
    S.variables.Q.valid_min = -Inf;
    S.aux.tables.other.file = '../../attributes_other_camelscz.csv';
    S.aux.tables.other.key = 'gauge_id';
    S.aux.tables.other.lat = 'gauge_lat';
    S.aux.tables.other.area = 'area';
    S.aux.tables.other.area_scale = 1e6;
    pVars = {'total_precipitation_sum_CHMI', ...
        'total_precipitation_sum_ERA5_LAND', ...
        'total_precipitation_sum_CHMI_ERA5_LAND'};
    tVars = {'temperature_2m_mean_CHMI', ...
        'temperature_2m_mean_ERA5_LAND', ...
        'temperature_2m_mean_CHMI_ERA5_LAND'};
    eVars = {'potential_evaporation_sum_FAO_PENMAN_MONTEITH_CHMI', ...
        'potential_evaporation_sum_FAO_PENMAN_MONTEITH_ERA5_LAND', ...
        'potential_evaporation_sum_ERA5_LAND'};
    Base = S;
    for p = 1:3
        for t = 1:3
            for pet = 1:3
                D = Base;
                D.variables.P.sources = [pVars(p) pVars];
                D.variables.P.sources = unique(D.variables.P.sources,'stable');
                D.variables.P.source_operation = 'aliases';
                D.variables.T.sources = [tVars(t) tVars];
                D.variables.T.sources = unique(D.variables.T.sources,'stable');
                D.variables.T.source_operation = 'aliases';
                D.variables.Ep.sources = [eVars(pet) eVars];
                D.variables.Ep.sources = unique(D.variables.Ep.sources,'stable');
                D.variables.Ep.source_operation = 'aliases';
                D.variables.Ep.clip_min = 0;
                name = sprintf('p%d_t%d_pet%d',p,t,pet);
                S.profiles.(name).match = struct('dt',1,'precip',p,'temp',t,'pet',pet);
                S.profiles.(name).schema = D;
            end
        end
    end

    HourlyBase = Base;
    HourlyBase.name = 'hourly CAMELS-CZ';
    HourlyBase.layout = 'multi_file_per_basin';
    HourlyBase.progress.label = '... Reading hourly CAMELS-CZ data';
    HourlyBase.timeline.reference = datetime(1950,10,1,0,0,0);
    HourlyBase.timeline.step = hours(1);
    % schema.id.prefix already adds "camelscz_" to the basin identifier.
    HourlyBase.files.hourly.pattern = '{gauge}.csv';
    HourlyBase.files.hourly.format = 'csv';
    HourlyBase.files.hourly.time.mode = 'column';
    HourlyBase.files.hourly.time.column = 'time';
    HourlyBase.files.hourly.time.input_format = 'yyyy-MM-dd HH:mm:ss';
    HourlyBase.files.daily.pattern = ...
        '../../daily/timeseries/{gauge}.csv';
    HourlyBase.files.daily.format = 'csv';
    HourlyBase.files.daily.time.mode = 'column';
    HourlyBase.files.daily.time.column = 'date';
    HourlyBase.files.daily.time.input_format = 'yyyy-MM-dd';
    HourlyBase.files.daily.time.join = 'calendar_day';
    HourlyBase.aux.tables.other.file = ...
        '../../attributes_other_camelscz.csv';
    pHourly = {'total_precipitation_sum_CHMI', ...
        'total_precipitation_ERA5_LAND'};
    tHourly = {'temperature_2m_mean_CHMI', ...
        'temperature_2m_ERA5_LAND'};
    eDaily = {'potential_evaporation_sum_FAO_PENMAN_MONTEITH_CHMI', ...
        'potential_evaporation_sum_FAO_PENMAN_MONTEITH_ERA5_LAND', ...
        'potential_evaporation_sum_ERA5_LAND'};
    for p = 1:2
        for t = 1:2
            for pet = 1:4
                H = HourlyBase;
                H.variables.P = local_variable( ...
                    pHourly{p},'mm/hour','mm/hour');
                H.variables.P.file = 'hourly';
                H.variables.T = local_variable( ...
                    tHourly{t},'degC','degC');
                H.variables.T.file = 'hourly';
                H.variables.Q = local_variable( ...
                    'streamflow','mm/hour','mm/hour');
                H.variables.Q.file = 'hourly';
                if pet == 1
                    H.variables.Ep = local_variable( ...
                        'potential_evaporation_ERA5_LAND', ...
                        'mm/hour','mm/hour');
                    H.variables.Ep.file = 'hourly';
                else
                    H.variables.Ep = local_variable( ...
                        eDaily{pet - 1},'mm/day','mm/hour');
                    H.variables.Ep.scale = 1/24;
                    H.variables.Ep.file = 'daily';
                end
                H.variables.Ep.clip_min = 0;
                name = sprintf('hourly_p%d_t%d_pet%d',p,t,pet);
                S.profiles.(name).match = struct( ...
                    'dt',24,'precip',p,'temp',t,'pet',pet);
                S.profiles.(name).schema = H;
            end
        end
    end
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
