function R = region_config_MX()
%REGION_CONFIG_MX Defaults for the Mexican Caravan-HYSETS subset.

    R = struct();
    R.code = 'CAMELS_MX';
    R.acronym = 'MX';
    R.name = 'Mexico';
    R.dataset = 'Caravan-HYSETS';
    R.data_root = 'CAMELS_MX';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'MX_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 46;
    X.basins.training = 37;
    X.basins.evaluation = 9;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/1998';
    X.period.manual.train_end = '30/09/2018';
    X.period.manual.eval_start = '01/10/1978';
    X.period.manual.eval_end = '30/09/1998';
    X.period.common_start = X.period.manual.eval_start;
    X.period.common_end = X.period.manual.train_end;
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    X.meteo.product.items = {'Caravan-HYSETS ERA5-Land [daily]'};
    X.meteo.product.default = X.meteo.product.items{1};
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 ERA5-Land precipitation'};
    X.meteo.precipitation.default = X.meteo.precipitation.items{1};
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 ERA5-Land temperature'};
    X.meteo.temperature.default = X.meteo.temperature.items{1};
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 FAO Penman-Monteith potential evaporation'};
    X.meteo.pet.default = X.meteo.pet.items{1};
    X.meteo.pet.enabled = false;

    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Daily forcing and streamflow';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = 'hysets_*.csv';
    C(1).minimum_count = 46;
    C(2).name = 'Caravan attributes';
    C(2).type = 'file';
    C(2).path = 'attributes_caravan_hysets.csv';
    C(3).name = 'HydroATLAS attributes';
    C(3).type = 'file';
    C(3).path = 'attributes_hydroatlas_hysets.csv';
    C(4).name = 'Gauge metadata';
    C(4).type = 'file';
    C(4).path = 'attributes_other_hysets.csv';
    X.checks = C;
    R.by_resolution.Daily = X;

    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_meteo_schema()

    S.name = 'daily Caravan-HYSETS Mexico';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'hysets_{gauge}.csv';
    S.file.delimiter = ',';
    S.id.pad_width = 0;
    S.timeline.reference = datetime(1950,1,1);
    S.timeline.step = days(1);
    S.time.column = 'date';
    S.time.input_format = '';
    S.variables.P = local_variable( ...
        'total_precipitation_sum','mm/day','mm/day');
    S.variables.P.clip_min = 0;
    S.variables.Ep = local_variable( ...
        'potential_evaporation_sum_FAO_PENMAN_MONTEITH', ...
        'mm/day','mm/day');
    S.variables.Ep.clip_min = 0;
    S.variables.T = local_variable( ...
        'temperature_2m_mean','degC','degC');
    S.variables.Q = local_variable( ...
        'streamflow','mm/day','mm/day');
    S.progress.label = ...
        '... Reading daily Caravan-HYSETS Mexico generic data';

end

function S = local_attribute_schema()

    S.name = 'Caravan-HYSETS Mexico attributes';
    files = {'attributes_other_hysets.csv', ...
        'attributes_caravan_hysets.csv', ...
        'attributes_hydroatlas_hysets.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id','gaugeid','station_id'};
    end
    S.id.column = 'gauge_id';
    S.id.uppercase = true;
    S.id.strip = true;
    S.id.regex = {'^HYSETS[_-]',''; '\.0+$',''};
    S.metadata.name_sources = {'gauge_name','station_name','name'};
    S.metadata.name_transform = 'title_case';
    S.region = 'CAMELS_MX';
    S.zone.region = 'MX';
    definitions = { ...
        'gauge_lat',{'gauge_lat','gauge_latitude','lat','latitude'},true; ...
        'gauge_lon',{'gauge_lon','gauge_longitude','lon','longitude'},true; ...
        'area',{'area','area_km2','catchment_area','area_calc'},true; ...
        'p_mean',{'p_mean','precipitation_mean'},false; ...
        'pet_mean_fao_pm',{'pet_mean_fao_pm','pet_mean_FAO_PM', ...
            'potential_evaporation_mean_fao_pm'},false; ...
        'aridity_fao_pm',{'aridity_fao_pm','aridity_FAO_PM','aridity'},false; ...
        'p_seasonality',{'seasonality_FAO_PM','p_seasonality', ...
            'precipitation_seasonality'},false; ...
        'frac_snow',{'frac_snow','fraction_snow'},false; ...
        'high_prec_freq',{'high_prec_freq'},false; ...
        'high_prec_dur',{'high_prec_dur'},false; ...
        'low_prec_freq',{'low_prec_freq'},false; ...
        'low_prec_dur',{'low_prec_dur'},false};
    S.aliases = repmat(struct('target','','sources',{{}}, ...
        'required',false,'default',NaN),size(definitions,1),1);
    for i = 1:size(definitions,1)
        S.aliases(i).target = definitions{i,1};
        S.aliases(i).sources = definitions{i,2};
        S.aliases(i).required = definitions{i,3};
    end
    S.progress.label = ...
        '... Reading Caravan-HYSETS Mexico generic attributes';

end

function V = local_variable(source,units,targetUnits)

    V.source = source;
    V.units = units;
    V.target_units = targetUnits;

end
