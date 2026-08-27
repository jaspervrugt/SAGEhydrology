function R = region_config_IL()
%REGION_CONFIG_IL Regional defaults for CAMELS-IL.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_IL';
    R.acronym = 'IL';
    R.name = 'Israel';
    R.dataset = 'CAMELS-IL';
    R.data_root = 'CAMELS_IL';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'IL_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 94;
    X.basins.training = 74;
    X.basins.evaluation = 20;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2004';
    X.period.manual.train_end = '30/09/2019';
    X.period.manual.eval_start = '01/10/1989';
    X.period.manual.eval_end = '30/09/2004';
    X.period.common_start = '01/10/1989';
    X.period.common_end = '30/09/2019';
    X.paths.run_root = fullfile('daily');
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-IL [daily]'};
    X.meteo.product.default = 'CAMELS-IL [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 ERA5-Land'};
    X.meteo.precipitation.default = '1 ERA5-Land';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 ERA5-Land'};
    X.meteo.temperature.default = '1 ERA5-Land';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = { ...
        '0 Zero PET', ...
        '1 Penman-Monteith', ...
        '2 Priestley-Taylor', ...
        '3 Makkink', ...
        '4 Supplied FAO Penman-Monteith', ...
        '5 Supplied ERA5-Land PET' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = 'il_*.csv';
    C(1).minimum_count = 95;
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
    S.name = 'CAMELS-IL attributes';
    files = { ...
        'attributes_other_il.csv', ...
        'attributes_caravan_il.csv', ...
        'attributes_hydroatlas_il.csv'};
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
    S.id.optional_prefix = 'il_';
    S.id.output_regex = {'^il_',''};
    S.id.output_uppercase = true;
    S.id.output_lowercase = false;
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_transform = 'title_case_words';
    S.region = 'CAMELS_IL';
    S.zone.region = 'IL';
    S.progress.label = '... Reading CAMELS-IL generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-IL';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'il_{gauge}.csv';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd';
    S.variables.P = local_variable('total_precipitation_sum','mm/day','mm/day');
    S.variables.T = local_variable('temperature_2m_mean','degC','degC');
    S.variables.Q = local_variable('streamflow','m3/s','mm/day');
    S.variables.Q.area_normalize = true;
    S.variables.Q.valid_min = 0;
    S.aux.tables.other.file = '../attributes_other_il.csv';
    S.aux.tables.other.key = 'gauge_id';
    S.aux.tables.other.strip_prefix = 'IL_';
    S.aux.tables.other.lat = 'gauge_lat';
    S.aux.tables.other.area = 'area';
    S.aux.tables.other.area_scale = 1e6;
    S.aux.series.elev.file_pattern = 'il_{gauge}.csv';
    S.aux.series.elev.source = 'surface_pressure_mean';
    S.aux.series.elev.reducer = 'median';
    S.aux.series.elev.derive = 'pressure_elevation';
    for method = 1:3
        D.variables.Tmin = local_variable('temperature_2m_min','degC','degC');
        D.variables.Tmax = local_variable('temperature_2m_max','degC','degC');
        D.variables.Tdew = local_variable('dewpoint_temperature_2m_mean','degC','degC');
        D.variables.WindU10 = local_variable('u_component_of_wind_10m_mean','m/s','m/s');
        D.variables.WindV10 = local_variable('v_component_of_wind_10m_mean','m/s','m/s');
        D.variables.Radiation = local_variable('surface_net_solar_radiation_mean','W/m2','W/m2');
        D.variables.VaporPressure.derive = 'vapor_pressure_dewpoint';
        D.variables.Wind2.derive = 'wind_uv_to_2m';
        D.variables.Ep.derive = 'fao56_daily';
        D.variables.Ep.method = method;
        name = sprintf('pet%d',method);
        S.profiles.(name).match.pet = method;
        S.profiles.(name).schema.variables = D.variables;
    end
    Z.variables.Ep.derive = 'constant';
    Z.variables.Ep.value = 0;
    S.profiles.zero.match.pet = 0;
    S.profiles.zero.schema.variables = Z.variables;
    F.variables.Ep = local_variable( ...
        'potential_evaporation_sum_FAO_PENMAN_MONTEITH','mm/day','mm/day');
    F.variables.Ep.clip_min = 0;
    S.profiles.supplied_fao.match.pet = 4;
    S.profiles.supplied_fao.schema.variables = F.variables;
    E.variables.Ep = local_variable( ...
        'potential_evaporation_sum_ERA5_LAND','mm/day','mm/day');
    E.variables.Ep.clip_min = 0;
    S.profiles.supplied_era5.match.pet = 5;
    S.profiles.supplied_era5.schema.variables = E.variables;
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
