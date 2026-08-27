function R = region_config_IS()
%REGION_CONFIG_IS Regional defaults for LamaH-Ice.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_IS';
    R.acronym = 'IS';
    R.name = 'Iceland';
    R.dataset = 'LamaH-Ice';
    R.data_root = 'CAMELS_IS';
    R.resolutions = {'Daily','Hourly'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'IS_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    % Start from the complete daily LamaH-Ice inventory. Basins that do
    % not satisfy the forcing and discharge requirements are excluded by
    % the common data-quality screening step before SAGE is initialized.
    X.basins.universe = 111;
    X.basins.file = 'IS_111_basins.txt';
    X.basins.training = 90;
    X.basins.evaluation = 21;
    X.basins.excluded_gauges = {'97'};
    X.basins.exclusion_reason = { ...
        'documented gauge-basin mapping inconsistency'};
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2012';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/2005';
    X.period.manual.eval_end = '30/09/2012';
    X.period.common_start = '01/10/2005';
    X.period.common_end = '30/09/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','forcing');
    X.paths.discharge = fullfile('daily','discharge');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'LamaH-Ice [daily]'};
    X.meteo.product.default = 'LamaH-Ice [daily]';
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
        '4 LamaH-Ice PET', ...
        '5 LamaH-Ice total ET' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','forcing');
    C(1).pattern = 'ID_*.csv';
    C(1).minimum_count = 111;
    C(2).name = 'Discharge data';
    C(2).type = 'folder';
    C(2).path = fullfile('daily','discharge');
    C(2).pattern = 'ID_*.csv';
    C(2).minimum_count = 111;
    X.checks = C;
    R.by_resolution.Daily = X;

    X = struct();
    X.label = 'Hourly';
    X.dt = 24;
    X.basins.universe = 76;
    X.basins.file = 'IS_76_basins.txt';
    X.basins.training = 65;
    X.basins.evaluation = 11;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2012';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/2005';
    X.period.manual.eval_end = '30/09/2012';
    X.period.common_start = '01/10/2005';
    X.period.common_end = '30/09/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('hourly','forcing');
    X.paths.discharge = fullfile('hourly','discharge');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'LamaH-Ice [hourly]'};
    X.meteo.product.default = 'LamaH-Ice [hourly]';
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
        '4 LamaH-Ice PET', ...
        '5 LamaH-Ice total ET' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('hourly','forcing');
    C(1).pattern = 'ID_*.csv';
    C(1).minimum_count = 76;
    C(2).name = 'Discharge data';
    C(2).type = 'folder';
    C(2).path = fullfile('hourly','discharge');
    C(2).pattern = 'ID_*.csv';
    C(2).minimum_count = 76;
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
    S.name = 'LamaH-Ice attributes';
    S.tables = repmat(struct('file','','keys',{{'id'}}, ...
        'delimiter',';','make_valid_names',true, ...
        'make_unique_names',true,'join_type','left', ...
        'duplicate_policy','suffix','duplicate_suffix','_catch'),2,1);
    S.tables(1).file = 'Gauge_attributes.csv';
    S.tables(2).file = 'Catchment_attributes.csv';
    S.tables(2).join_type = 'inner';
    S.id.column = 'id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'[^0-9]',''};
    S.id.sort = 'numeric';
    S.metadata.name_sources = {'name'};
    S.metadata.name_components = {'name','river'};
    S.metadata.name_separator = ': ';
    S.metadata.name_fallback_prefix = 'IS_';
    S.metadata.name_transform = '';
    S.aliases = repmat(struct('target','','sources',{{}}, ...
        'required',true,'default',NaN),2,1);
    S.aliases(1).target = 'gauge_elev';
    S.aliases(1).sources = {'elevation'};
    S.aliases(2).target = 'area';
    S.aliases(2).sources = {'area_calc'};
    S.projection.epsg = 3057;
    S.projection.x = 'lon';
    S.projection.y = 'lat';
    S.projection.latitude_target = 'gauge_lat';
    S.projection.longitude_target = 'gauge_lon';
    S.region = 'CAMELS_IS';
    S.zone.region = 'IS';
    S.progress.label = '... Reading CAMELS-IS generic attributes';
end

function S = local_meteo_schema()
    S.name = 'LamaH-Ice';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'ymd_columns';
    S.time.year_column = 'YYYY';
    S.time.month_column = 'MM';
    S.time.day_column = 'DD';
    S.variables.P = local_variable('prec','mm/step','mm/step');
    S.variables.T = local_variable('2m_temp_mean','degC','degC');
    S.aux.tables.catchments.file = fullfile( ...
        '..','..','Catchment_attributes.csv');
    S.aux.tables.catchments.key = 'id';
    S.aux.tables.catchments.elev = 'elev_mean';
    S.aux.tables.catchments.area = 'area_calc';
    S.aux.tables.catchments.area_scale = 1e6;
    S.aux.tables.gauges.file = fullfile( ...
        '..','..','gauge_information.txt');
    S.aux.tables.gauges.key = 'gauge_id';
    S.aux.tables.gauges.lat = 'gauge_lat';
    S.aux.tables.gauges.elev = 'gauge_elev';
    S.aux.tables.gauges.area = 'area';
    S.aux.tables.gauges.area_scale = 1e6;
    daily.file.pattern = 'ID_{gauge}.csv';
    daily.file.delimiter = ';';
    daily.timeline.step = days(1);
    hourly.file.pattern = 'ID_{gauge}.csv';
    hourly.file.delimiter = ';';
    hourly.timeline.step = hours(1);
    hourly.strict_time = false;
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
        dname = sprintf('daily_pet%d',method);
        S.profiles.(dname).match = struct('dt',1,'pet',method);
        S.profiles.(dname).schema = D;
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
        S.profiles.(hname).match = struct('dt',24,'pet',method);
        S.profiles.(hname).schema = H;
    end
    for dtValue = [1 24]
        if dtValue == 1
            base = daily;
            tag = 'daily';
        else
            base = hourly;
            tag = 'hourly';
        end
        Z = base;
        Z.variables.Ep.derive = 'constant';
        Z.variables.Ep.value = 0;
        name = [tag '_zero'];
        S.profiles.(name).match = struct('dt',dtValue,'pet',0);
        S.profiles.(name).schema = Z;
        P = base;
        P.variables.Ep = local_variable('pet','mm/step','mm/step');
        P.variables.Ep.valid_min = 0;
        name = [tag '_supplied'];
        S.profiles.(name).match = struct('dt',dtValue,'pet',4);
        S.profiles.(name).schema = P;
        E = base;
        E.variables.Ep = local_variable('total_et','mm/step','mm/step');
        E.variables.Ep.clip_min = 0;
        name = [tag '_total_et'];
        S.profiles.(name).match = struct('dt',dtValue,'pet',5);
        S.profiles.(name).schema = E;
    end
end

function S = local_discharge_schema()
    S.name = 'CAMELS-IS discharge'; 
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'ymd_columns'; 
    S.time.year_column = 'YYYY';
    S.time.month_column = 'MM'; 
    S.time.day_column = 'DD';
    S.variables.Q = local_variable('qobs','m3/s','mm/step');
    S.variables.Q.valid_min = 0; 
    S.variables.Q.area_normalize = true;
    S.variables.Q.quality.sources = {'qc_flag'};
    S.variables.Q.quality.valid_max = 100;
    S.profiles.daily.match.dt = 1;
    S.profiles.daily.schema.file.pattern = 'ID_{gauge}.csv';
    S.profiles.daily.schema.file.delimiter = ';';
    S.profiles.daily.schema.timeline.step = days(1);
    S.profiles.hourly.match.dt = 24;
    S.profiles.hourly.schema.file.pattern = 'ID_{gauge}.csv';
    S.profiles.hourly.schema.file.delimiter = ';';
    S.profiles.hourly.schema.timeline.step = hours(1);
    S.profiles.hourly.schema.time.mode = 'ymd_row_sequence';
    S.profiles.hourly.schema.time.steps_per_day = 24;
    S.profiles.hourly.schema.time.partial_first = 'final_slots';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
