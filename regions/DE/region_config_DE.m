function R = region_config_DE()
%REGION_CONFIG_DE Regional defaults for CAMELS-DE.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_DE';
    R.acronym = 'DE';
    R.name = 'Germany';
    R.dataset = 'CAMELS-DE';
    R.data_root = 'CAMELS_DE';
    R.resolutions = {'Daily','Hourly'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'DE_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 1490;
    X.basins.file = 'DE_1490_basins.txt';
    X.basins.training = 1200;
    X.basins.evaluation = 290;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2008';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/1996';
    X.period.manual.eval_end = '30/09/2008';
    X.period.common_start = '01/10/1996';
    X.period.common_end = '30/09/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'1 CAMELS-DE daily'};
    X.meteo.product.default = '1 CAMELS-DE daily';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 precipitation'};
    X.meteo.precipitation.default = '1 precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 air temperature'};
    X.meteo.temperature.default = '1 air temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = { ...
        '1 Penman-Monteith', ...
        '2 Priestley-Taylor', ...
        '3 Makkink' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
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
    X.basins.universe = 1259;
    X.basins.file = 'DE_1259_basins.txt';
    X.basins.training = 1000;
    X.basins.evaluation = 259;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2013';
    X.period.manual.train_end = '30/09/2021';
    X.period.manual.eval_start = '01/10/2005';
    X.period.manual.eval_end = '30/09/2013';
    X.period.common_start = '01/10/2005';
    X.period.common_end = '30/09/2021';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('hourly','timeseries');
    X.paths.discharge = fullfile('hourly','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'1 CAMELS-DE hourly'};
    X.meteo.product.default = '1 CAMELS-DE hourly';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 precipitation'};
    X.meteo.precipitation.default = '1 precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 air temperature'};
    X.meteo.temperature.default = '1 air temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = { ...
        '1 Penman-Monteith', ...
        '2 Priestley-Taylor', ...
        '3 Makkink' ...
        };
    X.meteo.pet.default = '1 Penman-Monteith';
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
    S = local_attribute_release('daily','CAMELS_DE_');
    H = local_attribute_release('hourly','CAMELS_DE_1h_');
    S.profiles.hourly.match.stream = 'hourly';
    S.profiles.hourly.schema = H;
end

function S = local_attribute_release(folder,prefix)
    suffixes = {'topographic_attributes.csv', ...
        'climatic_attributes.csv','hydrologic_attributes.csv', ...
        'hydrogeology_attributes.csv','soil_attributes.csv', ...
        'landcover_attributes.csv', ...
        'humaninfluence_attributes.csv'};
    S.name = ['CAMELS-DE ' folder ' attributes'];
    S.tables = repmat(struct('file','','keys',{{'gauge_id'}}, ...
        'make_valid_names',true,'duplicate_policy','suffix', ...
        'duplicate_suffix','_de'),numel(suffixes),1);
    for i = 1:numel(suffixes)
        S.tables(i).file = [prefix suffixes{i}];
    end
    S.root_candidates = {folder,''};
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'["'']',''; '\.0+$',''};
    S.id.sort = 'text';
    S.id.optional_prefix = 'DE';
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_transform = '';
    S.aliases = repmat(struct('target','','sources',{{}}, ...
        'required',false,'default',100),2,1);
    S.aliases(1).target = 'area_pct_in_germany';
    S.aliases(1).sources = {'area_pct_in_germany'};
    S.aliases(2).target = 'precipitation_perc_complete';
    S.aliases(2).sources = {'precipitation_perc_complete'};
    S.region = 'CAMELS_DE';
    S.zone.region = 'DE';
    S.progress.label = ['... Reading CAMELS-DE ' ...
        folder ' generic attributes'];
end

function S = local_meteo_schema()
    S.name = 'CAMELS-DE';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.timeline.reference = datetime(1950,10,1);
    S.time.mode = 'column';
    S.time.column = 'date';
    for hourly = [false true]
        if hourly
            tag = 'hourly';
            dt = 24;
            B.file.pattern = 'CAMELS_DE_1h_hydromet_timeseries_{gauge}.csv';
            B.timeline.step = hours(1);
            B.time.input_format = 'yyyy-MM-dd HH:mm:ss';
            B.variables.P = local_variable('precipitation_mean_gapfilled','mm/step','mm/step');
            B.variables.T = local_variable('air_temperature_mean','degC','degC');
            B.aux.tables.topo.file = '../CAMELS_DE_1h_topographic_attributes.csv';
        else
            tag = 'daily';
            dt = 1;
            B.file.pattern = 'CAMELS_DE_hydromet_timeseries_{gauge}.csv';
            B.timeline.step = days(1);
            B.time.input_format = 'yyyy-MM-dd';
            B.variables.P = local_variable('precipitation_mean','mm/step','mm/step');
            B.variables.T = local_variable('temperature_mean','degC','degC');
            B.aux.tables.topo.file = '../CAMELS_DE_topographic_attributes.csv';
        end
        B.id.prefix = 'DE';
        B.variables.Q = local_variable('discharge_spec_obs','mm/step','mm/step');
        B.aux.tables.topo.key = 'gauge_id';
        B.aux.tables.topo.lat = 'gauge_lat';
        B.aux.tables.topo.elev = 'elev_mean';
        B.aux.tables.topo.area = 'area';
        B.aux.tables.topo.area_scale = 1e6;
        Z = B;
        Z.variables.Ep.derive = 'constant';
        Z.variables.Ep.value = 0;
        name = [tag '_zero'];
        S.profiles.(name).match = struct('dt',dt,'pet',0);
        S.profiles.(name).schema = Z;
        for method = 1:3
            D = B;
            if hourly
                D.variables.Radiation = local_variable('global_radiation_mean','W/m2','W/m2');
                D.variables.Pressure = local_variable('air_pressure_surface_mean','hPa','Pa');
                D.variables.Pressure.scale = 100;
                D.variables.MixingRatio = local_variable('water_vapor_mixing_ratio_mean','g/kg','kg/kg');
                D.variables.MixingRatio.scale = 1/1000;
                D.variables.SpecificHumidity.derive = 'mixing_ratio_to_specific_humidity';
                D.variables.SpecificHumidity.input = 'MixingRatio';
                D.variables.SpecificHumidity.input_scale = 1/1000;
                D.variables.WindU10 = local_variable('wind_speed_eastward_mean','m/s','m/s');
                D.variables.WindV10 = local_variable('wind_speed_northward_mean','m/s','m/s');
                D.variables.Tdew = local_variable('dew_point_temperature_mean','degC','degC');
                D.variables.VaporPressure.derive = 'vapor_pressure_dewpoint';
                D.variables.CloudCover = local_variable('cloud_cover_mean','percent','percent');
                D.variables.Ep.derive = 'fao56_hourly_estimated_longwave';
                D.variables.Ep.pressure_scale = 100;
            else
                D.variables.Tmin = local_variable('temperature_min','degC','degC');
                D.variables.Tmax = local_variable('temperature_max','degC','degC');
                D.variables.Radiation = local_variable('radiation_global_mean','MJ/m2/day','W/m2');
                D.variables.Radiation.scale = 1/0.0864;
                D.variables.RelativeHumidity = local_variable('humidity_mean','percent','percent');
                D.variables.VaporPressure.derive = 'vapor_pressure_relative_humidity';
                D.variables.Wind2.derive = 'constant';
                D.variables.Wind2.value = 2;
                D.variables.Ep.derive = 'fao56_daily';
                D.variables.Ep.radiation_scale = 1/0.0864;
            end
            D.variables.Ep.method = method;
            D.variables.Ep.clip_min = 0;
            name = sprintf('%s_pet%d',tag,method);
            S.profiles.(name).match = struct('dt',dt,'pet',method);
            S.profiles.(name).schema = D;
        end
        clear B D Z
    end
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
