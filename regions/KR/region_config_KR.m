function R = region_config_KR()
%REGION_CONFIG_KR Regional defaults for hourly CAMELSH-KR.

    R = struct();
    R.code = 'CAMELSH_KR';
    R.acronym = 'KR';
    R.name = 'South Korea';
    R.dataset = 'CAMELSH-KR';
    R.data_root = 'CAMELSH_KR';
    R.resolutions = {'Hourly'};
    R.default_resolution = 'Hourly';
    R.basin_file_pattern = 'KR_%d_basins.txt';

    X = struct();
    X.label = 'Hourly';
    X.dt = 24;
    % Thirteen source basins are excluded from the default periods:
    % eight have insufficient observed-Q coverage and five contain severe
    % source-data defects documented in CAMELSH_KR_discharge_audit.md.
    X.basins.universe = 165;
    X.basins.training = 140;
    X.basins.evaluation = 25;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2015';
    X.period.manual.train_end = '30/09/2019';
    X.period.manual.eval_start = '01/10/2010';
    X.period.manual.eval_end = '30/09/2015';
    X.period.common_start = '01/10/2010';
    X.period.common_end = '30/09/2019';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('hourly','timeseries');
    X.paths.discharge = fullfile('hourly','timeseries');

    X.meteo.product.items = {'CAMELSH-KR [hourly]'};
    X.meteo.product.default = 'CAMELSH-KR [hourly]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 ERA5-Land precipitation','2 observed precipitation'};
    X.meteo.precipitation.default = '1 ERA5-Land precipitation';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 ERA5-Land temperature','2 observed temperature'};
    X.meteo.temperature.default = '1 ERA5-Land temperature';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = {'1 ERA5-Land potential evaporation'};
    X.meteo.pet.default = '1 ERA5-Land potential evaporation';
    X.meteo.pet.enabled = false;

    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Hourly forcing and streamflow';
    C(1).type = 'folder';
    C(1).path = fullfile('hourly','timeseries');
    C(1).pattern = '*.csv';
    C(1).minimum_count = 178;
    C(2).name = 'General attributes';
    C(2).type = 'file';
    C(2).path = 'attributes_general.csv';
    C(3).name = 'Observed climate attributes';
    C(3).type = 'file';
    C(3).path = 'attributes_climate_obs.csv';
    C(4).name = 'ERA5-Land climate attributes';
    C(4).type = 'file';
    C(4).path = 'attributes_climate_ERA5Land.csv';
    C(5).name = 'HydroATLAS attributes';
    C(5).type = 'file';
    C(5).path = 'attributes_HydroATLAS.csv';
    C(6).name = 'Dam attributes';
    C(6).type = 'file';
    C(6).path = 'attributes_dam.csv';
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
    files = {'attributes_general.csv', ...
        'attributes_climate_obs.csv', ...
        'attributes_climate_ERA5Land.csv', ...
        'attributes_HydroATLAS.csv','attributes_dam.csv', ...
        'gauge_information.txt'};
    S.name = 'CAMELSH-KR attributes';
    S.tables = repmat(struct('file','','keys',{{'station_id'}}, ...
        'column_renames',{{}},'keep_columns',{{}}, ...
        'absolute_columns',{{}},'delimiter','', ...
        'duplicate_policy','drop'),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
    end
    S.tables(1).column_renames = { ...
        'STAID','station_id'; 'Lat_gage','gauge_lat'; ...
        'Lon_gage','gauge_lon'; 'Area','area'; ...
        'Area_HydroATLAS','area_hydroatlas'};
    S.tables(1).keep_columns = {'station_id','gauge_lat', ...
        'gauge_lon','area','area_hydroatlas'};
    S.tables(2).column_renames = { ...
        'STAID','station_id'; 'aridity_index','aridity'};
    S.tables(2).keep_columns = {'station_id','p_mean','aridity', ...
        'p_seasonality','frac_snow','high_prec_freq', ...
        'high_prec_dur','low_prec_freq','low_prec_dur'};
    S.tables(2).absolute_columns = {'aridity'};
    S.tables(3).column_renames = { ...
        'STAID','station_id'; 'p_mean','era5_p_mean'; ...
        'pet_mean','era5_pet_mean'; ...
        'aridity_index','era5_aridity'; ...
        'p_seasonality','era5_p_seasonality'; ...
        'frac_snow','era5_frac_snow'; ...
        'high_prec_freq','era5_high_prec_freq'; ...
        'high_prec_dur','era5_high_prec_dur'; ...
        'low_prec_freq','era5_low_prec_freq'; ...
        'low_prec_dur','era5_low_prec_dur'};
    S.tables(3).keep_columns = {'station_id','era5_p_mean', ...
        'era5_pet_mean','era5_aridity','era5_p_seasonality', ...
        'era5_frac_snow','era5_high_prec_freq', ...
        'era5_high_prec_dur','era5_low_prec_freq', ...
        'era5_low_prec_dur'};
    S.tables(3).absolute_columns = ...
        {'era5_pet_mean','era5_aridity'};
    S.tables(4).column_renames = {'STAID','station_id'};
    S.tables(5).column_renames = { ...
        'STAID','station_id'; 'vol_tot','dam_vol_tot'; ...
        'vol_eff','dam_vol_eff'; 'vol_low','dam_vol_low'; ...
        'vol_flood','dam_vol_flood'; 'WL_flood','dam_wl_flood'; ...
        'WL_normal','dam_wl_normal'; 'WL_low','dam_wl_low'};
    S.tables(5).keep_columns = {'station_id','dam_vol_tot', ...
        'dam_vol_eff','dam_vol_low','dam_vol_flood', ...
        'dam_wl_flood','dam_wl_normal','dam_wl_low'};
    S.tables(6).delimiter = sprintf('\t');
    S.tables(6).column_renames = {'gauge_id','station_id'};
    S.tables(6).keep_columns = ...
        {'station_id','gauge_name','gauge_name_kr'};
    S.id.column = 'station_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.id.sort = 'numeric';
    S.metadata.name_components = {'gauge_name'};
    S.metadata.parenthetical_components = {'gauge_name_kr'};
    S.metadata.name_transform = '';
    S.region = 'CAMELSH_KR';
    S.zone.region = 'KR';
    S.progress.label = '... Reading CAMELSH-KR generic attributes';
end

function S = local_meteo_schema()
    S.name = 'hourly CAMELSH-KR';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = '{gauge}.csv';
    S.file.delimiter = ',';
    S.id.pad_width = 0;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = hours(1);
    S.time.column = 'DateTime';
    S.time.input_format = 'dd-MMM-yyyy HH:mm:ss';
    S.variables.P = local_variable( ...
        'total_precipitation','mm/hour','mm/hour');
    S.variables.P.clip_min = 0;
    S.variables.Ep = local_variable( ...
        'potential_evaporation','mm/hour','mm/hour');
    S.variables.Ep.scale = -1;
    S.variables.Ep.clip_min = 0;
    S.variables.T = local_variable( ...
        'temperature_2m','degC','degC');
    S.variables.Q = local_variable( ...
        'streamflow','m3/s','mm/hour');
    S.variables.Q.area_normalize = true;
    S.progress.label = ...
        '... Reading hourly CAMELSH-KR generic data';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
