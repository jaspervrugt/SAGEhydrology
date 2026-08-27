function R = region_config_AT()
%REGION_CONFIG_AT Regional defaults for CAMELS-AT.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_AT';
    R.acronym = 'AT';
    R.name = 'Austria';
    R.dataset = 'CAMELS-AT';
    R.data_root = 'CAMELS_AT';
    R.resolutions = {'Daily','Hourly'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'AT_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 859;
    X.basins.training = 700;
    X.basins.evaluation = 159;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2002';
    X.period.manual.train_end = '30/09/2017';
    X.period.manual.eval_start = '01/10/1987';
    X.period.manual.eval_end = '30/09/2002';
    X.period.common_start = '01/10/1987';
    X.period.common_end = '30/09/2017';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries','forcing');
    X.paths.discharge = fullfile('daily','timeseries','streamflow');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-AT [daily]'};
    X.meteo.product.default = 'CAMELS-AT [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = { ...
        '1 Mean', ...
        '2 (Tmin+Tmax)/2' ...
        };
    X.meteo.temperature.default = '1 Mean';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '0 Zero PET', ...
        '1 Penman-Monteith', ...
        '2 Priestley-Taylor', ...
        '3 Makkink', ...
        '4 LamaH ΣET' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries','forcing');
    C(1).pattern = '*';
    C(1).minimum_count = 1;
    C(2).name = 'Discharge data';
    C(2).type = 'folder';
    C(2).path = fullfile('daily','timeseries','streamflow');
    C(2).pattern = '*';
    C(2).minimum_count = 1;
    X.checks = C;
    R.by_resolution.Daily = X;

    X = struct();
    X.label = 'Hourly';
    X.dt = 24;
    X.basins.universe = 859;
    X.basins.training = 700;
    X.basins.evaluation = 159;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2005';
    X.period.manual.train_end = '30/09/2017';
    X.period.manual.eval_start = '01/10/1993';
    X.period.manual.eval_end = '30/09/2005';
    X.period.common_start = '01/10/1993';
    X.period.common_end = '30/09/2017';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('hourly','timeseries','forcing');
    X.paths.discharge = fullfile('hourly','timeseries','streamflow');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-AT [hourly]'};
    X.meteo.product.default = 'CAMELS-AT [hourly]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 Mean'};
    X.meteo.temperature.default = '1 Mean';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = { ...
        '0 Zero PET', ...
        '1 Penman-Monteith', ...
        '2 Priestley-Taylor', ...
        '3 Makkink', ...
        '4 LamaH ΣET' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('hourly','timeseries','forcing');
    C(1).pattern = '*';
    C(1).minimum_count = 1;
    C(2).name = 'Discharge data';
    C(2).type = 'folder';
    C(2).path = fullfile('hourly','timeseries','streamflow');
    C(2).pattern = '*';
    C(2).minimum_count = 1;
    X.checks = C;
    R.by_resolution.Hourly = X;


    % Declarative generic-reader schemas owned by this region module.
    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = local_discharge_schema();
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_attribute_schema()
    S.name = 'CAMELS-AT attributes';
    S.tables = repmat(struct('file','','keys',{{'ID'}}, ...
        'make_valid_names',true,'make_unique_names',true, ...
        'numeric_text','threshold','numeric_text_ratio',0.9, ...
        'numeric_text_exceptions',{{'ID'}},'join_type','left', ...
        'duplicate_policy','suffix','duplicate_suffix','_catch'),2,1);
    S.tables(1).file = 'Gauge_attributes.csv';
    S.tables(2).file = 'Catchment_attributes.csv';
    S.tables(2).join_type = 'inner';
    S.id.column = 'ID';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'^ID[_\-\s]*',''; '[^0-9]',''};
    S.id.sort = 'numeric';
    S.metadata.name_sources = {'name'};
    S.metadata.name_components = {'name','river'};
    S.metadata.name_separator = ': ';
    S.metadata.name_fallback_prefix = 'AT_';
    S.metadata.name_transform = '';
    S.metadata.standardize = false;
    S.aliases.target = 'area';
    S.aliases.sources = {'area_gov','area_calc'};
    S.aliases.required = false;
    S.aliases.default = NaN;
    S.projection.epsg = 3035;
    S.projection.x = 'lon';
    S.projection.y = 'lat';
    S.projection.latitude_target = 'gauge_lat';
    S.projection.longitude_target = 'gauge_lon';
    S.selection.available_pattern = ...
        fullfile('daily','timeseries','forcing','ID_*.csv');
    S.selection.available_regex = '^ID_(\d+)\.csv$';
    S.region = 'CAMELS_AT';
    S.zone.region = 'AT';
    S.progress.label = '... Reading CAMELS-AT generic attributes';
end

function S = local_meteo_schema()
    S.name = 'CAMELS-AT';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'ymd_columns';
    S.time.year_column = 'YYYY';
    S.time.month_column = 'MM';
    S.time.day_column = 'DD';
    S.variables.P = local_variable('prec','mm/step','mm/step');
    S.variables.T = local_variable('2m_temp_mean','degC','degC');
    S.aux.tables.catchments.file = '../../../Catchment_attributes.csv';
    S.aux.tables.catchments.key = 'ID';
    S.aux.tables.catchments.elev = 'elev_mean';
    S.aux.tables.catchments.area = 'area_calc';
    S.aux.tables.catchments.area_scale = 1e6;
    S.aux.tables.gauges.file = '../../../gauge_information.txt';
    S.aux.tables.gauges.key = 'gauge_id';
    S.aux.tables.gauges.lat = 'gauge_lat';
    S.aux.tables.gauges.elev = 'gauge_elev';
    S.aux.tables.gauges.area = 'area_km2';
    S.aux.tables.gauges.area_scale = 1e6;
    daily.file.pattern = 'ID_{gauge}.csv';
    daily.file.delimiter = ';';
    daily.timeline.step = days(1);
    hourly.file.pattern = 'ID_{gauge}.csv';
    hourly.file.delimiter = ';';
    hourly.timeline.step = hours(1);
    hourly.time.hour_column = 'hh';
    hourly.time.minute_column = 'mm';
    hourly.variables.T.source = '2m_temp';
    for method = 1:3
        D = daily;
        D.variables.Tmin = local_variable('2m_temp_min','degC','degC');
        D.variables.Tmax = local_variable('2m_temp_max','degC','degC');
        D.variables.Tdew = local_variable('2m_dp_temp_mean','degC','degC');
        D.variables.WindU10 = local_variable('10m_wind_u','m/s','m/s');
        D.variables.WindV10 = local_variable('10m_wind_v','m/s','m/s');
        D.variables.Radiation = local_variable( ...
            'surf_net_solar_rad_mean','W/m2','W/m2');
        D.variables.VaporPressure.derive = 'vapor_pressure_dewpoint';
        D.variables.Wind2.derive = 'wind_uv_to_2m';
        D.variables.Ep.derive = 'fao56_daily';
        D.variables.Ep.method = method;
        for temp = 1:2
            DT = local_daily_temperature(D,temp);
            dname = sprintf('daily_pet%d_t%d',method,temp);
            S.profiles.(dname).match = struct( ...
                'dt',1,'pet',method,'temp',temp);
            S.profiles.(dname).schema = DT;
        end

        H = hourly;
        H.variables.Tdew = local_variable('2m_dp_temp','degC','degC');
        H.variables.WindU10 = local_variable('10m_wind_u','m/s','m/s');
        H.variables.WindV10 = local_variable('10m_wind_v','m/s','m/s');
        H.variables.NetSolar = local_variable('surf_net_solar_rad','W/m2','W/m2');
        H.variables.NetThermal = local_variable('surf_net_therm_rad','W/m2','W/m2');
        H.variables.Pressure = local_variable('surf_press','Pa','Pa');
        H.variables.Ep.derive = 'lamah_hourly';
        H.variables.Ep.method = method;
        hname = sprintf('hourly_pet%d',method);
        S.profiles.(hname).match = struct( ...
            'dt',24,'pet',method,'temp',1);
        S.profiles.(hname).schema = H;
    end

    dailyZero = daily;
    dailyZero.variables.Ep.derive = 'constant';
    dailyZero.variables.Ep.value = 0;
    for temp = 1:2
        D = local_daily_temperature(dailyZero,temp);
        name = sprintf('daily_pet0_t%d',temp);
        S.profiles.(name).match = struct( ...
            'dt',1,'pet',0,'temp',temp);
        S.profiles.(name).schema = D;
    end

    hourlyZero = hourly;
    hourlyZero.variables.Ep.derive = 'constant';
    hourlyZero.variables.Ep.value = 0;
    S.profiles.hourly_pet0.match = struct( ...
        'dt',24,'pet',0,'temp',1);
    S.profiles.hourly_pet0.schema = hourlyZero;

    daily.variables.Ep = local_variable('total_et','mm/step','mm/step');
    daily.variables.Ep.transform = @(x,context) abs(x);
    for temp = 1:2
        D = local_daily_temperature(daily,temp);
        name = sprintf('daily_pet4_t%d',temp);
        S.profiles.(name).match = struct( ...
            'dt',1,'pet',4,'temp',temp);
        S.profiles.(name).schema = D;
    end
    hourly.variables.Ep = local_variable('total_et','mm/step','mm/step');
    hourly.variables.Ep.transform = @(x,context) abs(x);
    S.profiles.hourly_native.match = struct( ...
        'dt',24,'pet',4,'temp',1);
    S.profiles.hourly_native.schema = hourly;
end

function S = local_daily_temperature(S,temp)
    if temp == 2
        S.variables.T.source = [];
        S.variables.T.sources = {'2m_temp_min','2m_temp_max'};
        S.variables.T.source_operation = 'mean';
    end
end

function S = local_discharge_schema()
    S.name = 'CAMELS-AT discharge'; 
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'ymd_columns'; 
    S.time.year_column = 'YYYY';
    S.time.month_column = 'MM'; 
    S.time.day_column = 'DD';
    S.variables.Q = local_variable('qobs','m3/s','mm/step');
    S.variables.Q.invalid_le = -999; 
    S.variables.Q.area_normalize = true;
    S.profiles.daily.match.dt = 1;
    S.profiles.daily.schema.file.pattern = 'ID_{gauge}.csv';
    S.profiles.daily.schema.file.delimiter = ';';
    S.profiles.daily.schema.timeline.step = days(1);
    S.profiles.hourly.match.dt = 24;
    S.profiles.hourly.schema.file.pattern = 'ID_{gauge}.csv';
    S.profiles.hourly.schema.file.delimiter = ';';
    S.profiles.hourly.schema.timeline.step = hours(1);
    S.profiles.hourly.schema.time.hour_column = 'hh';
    S.profiles.hourly.schema.time.minute_column = 'mm';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
