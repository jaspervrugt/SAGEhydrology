function R = region_config_AU()
%REGION_CONFIG_AU Regional defaults for CAMELS-AU.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_AU';
    R.acronym = 'AU';
    R.name = 'Australia';
    R.dataset = 'CAMELS-AU';
    R.data_root = 'CAMELS_AU';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'AU_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 222;
    X.basins.training = 170;
    X.basins.evaluation = 52;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/1999';
    X.period.manual.train_end = '30/09/2014';
    X.period.manual.eval_start = '01/10/1984';
    X.period.manual.eval_end = '30/09/1999';
    X.period.common_start = '01/10/1984';
    X.period.common_end = '30/09/2014';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily');
    X.paths.discharge = fullfile('daily','streamflow');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-AU [daily]'};
    X.meteo.product.default = 'CAMELS-AU [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 SILO precipitation', ...
        '2 AWAP precipitation' ...
        };
    X.meteo.precipitation.default = '1 SILO precipitation';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 Mean Tmin/Tmax (SILO)', ...
        '2 Mean Tmin/Tmax (AWAP)' ...
        };
    X.meteo.temperature.default = '1 Mean Tmin/Tmax (SILO)';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 Synthetic evaporation (SILO)', ...
        '2 Pan evaporation (SILO)', ...
        '3 Morton lake evaporation (SILO)', ...
        '4 Tall crop ET (SILO)', ...
        '5 Short crop ET (SILO)', ...
        '6 Morton wet ET (SILO)', ...
        '7 Morton point ET (SILO)', ...
        '8 Morton actual ET (SILO)' ...
        };
    X.meteo.pet.default = '1 Synthetic evaporation (SILO)';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Precipitation';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','precipitation');
    C(1).pattern = '*';
    C(1).minimum_count = 1;
    C(2).name = 'Temperature';
    C(2).type = 'folder';
    C(2).path = fullfile('daily','temperature');
    C(2).pattern = '*';
    C(2).minimum_count = 1;
    C(3).name = 'Evaporative demand';
    C(3).type = 'folder';
    C(3).path = fullfile('daily','evaporative_demand');
    C(3).pattern = '*';
    C(3).minimum_count = 1;
    C(4).name = 'Streamflow';
    C(4).type = 'folder';
    C(4).path = fullfile('daily','streamflow');
    C(4).pattern = '*';
    C(4).minimum_count = 1;
    X.checks = C;
    R.by_resolution.Daily = X;


    % Declarative generic-reader schemas owned by this region module.
    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = local_discharge_schema();
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_attribute_schema()
    S.name = 'CAMELS-AUS attributes';
    S.tables.file = 'CAMELS_AUS_Attributes&Indices_MasterTable.csv';
    S.tables.keys = {'station_id'};
    S.tables.make_valid_names = true;
    S.tables.key_type = 'char';
    S.id.column = 'station_id';
    S.id.uppercase = true;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.metadata.name_sources = {'station_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_AU';
    S.zone.region = 'AU';
    S.progress.label = '... Reading CAMELS-AUS generic attributes';
end

function S = local_meteo_schema()

    S.name = 'daily CAMELS-AU';
    S.format = 'csv';
    S.layout = 'wide_files';
    
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    
    S.time.mode = 'ymd_columns';
    S.time.year_column = 'year';
    S.time.month_column = 'month';
    S.time.day_column = 'day';
    
    % Precipitation
    S.variables.P.files = { ...
        'precipitation/precipitation_SILO.csv', ...
        'precipitation/precipitation_AWAP.csv'};
    S.variables.P.selector = 'precip';
    S.variables.P.default = 1;
    S.variables.P.units = 'mm/day';
    S.variables.P.target_units = 'mm/day';
    
    % PET
    S.variables.Ep.files = { ...
        'evaporative_demand/evap_syn_SILO.csv', ...
        'evaporative_demand/evap_pan_SILO.csv', ...
        'evaporative_demand/evap_morton_lake_SILO.csv', ...
        'evaporative_demand/et_tall_crop_SILO.csv', ...
        'evaporative_demand/et_short_crop_SILO.csv', ...
        'evaporative_demand/et_morton_wet_SILO.csv', ...
        'evaporative_demand/et_morton_point_SILO.csv', ...
        'evaporative_demand/et_morton_actual_SILO.csv'};
    S.variables.Ep.selector = 'pet';
    S.variables.Ep.default = 1;
    S.variables.Ep.units = 'mm/day';
    S.variables.Ep.target_units = 'mm/day';
    
    % Temperature
    S.variables.Tmin.files = { ...
        'temperature/tmin_SILO.csv', ...
        'temperature/tmin_AWAP.csv'};
    S.variables.Tmin.selector = 'temp';
    
    S.variables.Tmax.files = { ...
        'temperature/tmax_SILO.csv', ...
        'temperature/tmax_AWAP.csv'};
    S.variables.Tmax.selector = 'temp';
    
    % Native missing codes
    S.variables.P.invalid_le = -99;
    S.variables.Ep.invalid_le = -99;

    S.variables.Tmin.units = 'degC';
    S.variables.Tmin.target_units = 'degC';
    S.variables.Tmin.invalid_le = -99;

    S.variables.Tmax.units = 'degC';
    S.variables.Tmax.target_units = 'degC';
    S.variables.Tmax.invalid_le = -99;

    S.variables.T.units = 'degC';
    S.variables.T.target_units = 'degC';
    S.variables.T.derive = 'mean_Tmin_Tmax';

    S.progress.label = ...
        '... Reading daily CAMELS-AU generic data';
end

function S = local_discharge_schema()
    S.name = 'daily CAMELS-AU discharge';
    S.format = 'csv';
    S.layout = 'wide_files';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'ymd_columns';
    S.time.year_column = 'year';
    S.time.month_column = 'month';
    S.time.day_column = 'day';
    S.variables.Q.files = {'streamflow_mmd.csv'};
    S.variables.Q.units = 'mm/day';
    S.variables.Q.target_units = 'mm/day';
    S.variables.Q.invalid_le = -99;
    S.progress.label = '... Reading daily CAMELS-AU generic discharge';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
