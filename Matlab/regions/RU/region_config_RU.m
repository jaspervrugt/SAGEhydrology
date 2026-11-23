function R = region_config_RU()
%REGION_CONFIG_RU HydroCIS v1.3 daily data for the Russian Federation.
% Abramov and Kurochkina (2023), https://doi.org/10.5281/zenodo.8432070

    R = struct();
    R.code = 'HYDRO_CIS';
    R.acronym = 'RU';
    R.name = 'Russia';
    R.dataset = 'HydroCIS';
    R.data_root = 'HYDRO_CIS';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'RU_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    % The April 2023 area-quality screen retained 1,886 Russian basins.
    % The v1.3 archive also contains 278 GRDC gauges elsewhere in the CIS;
    % The consolidated data installer excludes those from this region.
    X.basins.universe = 1886;
    X.basins.training = 1500;
    X.basins.evaluation = 386;
    X.basins.file = 'RU_1886_basins.txt';
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2014';
    X.period.manual.train_end = '30/09/2020';
    X.period.manual.eval_start = '01/10/2009';
    X.period.manual.eval_end = '30/09/2014';
    X.period.common_start = '01/01/2008';
    X.period.common_end = '31/12/2020';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');

    X.meteo.product.items = {'HydroCIS v1.3 [daily]'};
    X.meteo.product.default = 'HydroCIS v1.3 [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 ERA5-Land precipitation', ...
        '2 ERA5 precipitation', ...
        '3 MSWEP precipitation', ...
        '4 GPCP precipitation', ...
        '5 IMERG precipitation'};
    X.meteo.precipitation.default = '1 ERA5-Land precipitation';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 ERA5-Land mean temperature', ...
        '2 ERA5 mean temperature'};
    X.meteo.temperature.default = '1 ERA5-Land mean temperature';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 Oudin potential evaporation', ...
        '2 GLEAM potential evaporation'};
    X.meteo.pet.default = '1 Oudin potential evaporation';
    X.meteo.pet.enabled = true;

    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),6,1);
    C(1).name = 'Daily forcing and discharge';
    C(1).type = 'folder';
    C(1).path = X.paths.meteo;
    C(1).pattern = '*.nc';
    C(1).minimum_count = 1886;
    C(2).name = 'HydroATLAS attributes';
    C(2).type = 'file'; C(2).path = 'static_data.csv';
    C(3).name = 'Derived climate attributes';
    C(3).type = 'file'; C(3).path = 'climate_attributes.csv';
    C(4).name = 'Gauge information';
    C(4).type = 'file'; C(4).path = 'gauge_information.txt';
    C(5).name = 'Gauge geometry';
    C(5).type = 'file'; C(5).path = fullfile('geometry','russia_gauges.gpkg');
    C(6).name = 'Catchment geometry';
    C(6).type = 'file'; C(6).path = fullfile('geometry','russia_ws.gpkg');
    X.checks = C;
    R.by_resolution.Daily = X;

    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = local_discharge_schema();
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();
end

function S = local_meteo_schema()
    S.name = 'daily HydroCIS v1.3';
    S.format = 'netcdf';
    S.layout = 'one_file_per_basin';
    S.file.pattern = '{gauge}.nc';
    S.id.pad_width = 0;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'variable';
    % The published NetCDF coordinate is named index and carries the
    % CF-style attribute "days since 2008-01-01".
    S.time.column = 'index';
    S.time.units = 'days';
    S.time.origin = datetime(2008,1,1);
    S.time.snap = 'day';
    S.time.implicit_index = true;
    S.variables.P = local_variable('prcp_e5l','mm/day');
    % HydroCIS contains nominal zero-precipitation values as small as
    % -3.4e-5 mm/day because of floating-point roundoff.  Accept only this
    % negligible negative range and clip it to physically meaningful zero;
    % larger negative values remain invalid.
    S.variables.P.valid_min = -1e-3;
    S.variables.P.clip_min = 0;
    S.variables.Tmin = local_variable('t_min_e5l','degC');
    S.variables.Tmax = local_variable('t_max_e5l','degC');
    S.variables.T = local_derived('mean_tmin_tmax','degC');
    S.variables.Ep = local_derived('oudin','mm/day');
    S.variables.Ep.valid_min = 0;
    S.resolve = @local_select_forcing;
    S.progress.label = '... Reading daily HydroCIS forcing';
end

function S = local_discharge_schema()
    S.name = 'daily observed HydroCIS discharge';
    S.format = 'netcdf';
    S.layout = 'one_file_per_basin';
    S.file.pattern = '{gauge}.nc';
    S.id.pad_width = 0;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'variable';
    S.time.column = 'index';
    S.time.units = 'days';
    S.time.origin = datetime(2008,1,1);
    S.time.snap = 'day';
    S.time.implicit_index = true;
    S.variables.Q = local_variable('q_mm_day','mm/day');
    S.variables.Q.valid_min = 0;
    S.progress.label = '... Reading daily HydroCIS discharge';
end

function S = local_attribute_schema()
    S.name = 'HydroCIS catchment attributes';
    S.tables = repmat(struct('file','','keys',{{'gauge_id'}}, ...
        'column_renames',{{}},'keep_columns',{{}}, ...
        'absolute_columns',{{}},'delimiter',',', ...
        'duplicate_policy','drop'),3,1);
    S.tables(1).file = 'gauge_information.txt';
    S.tables(1).delimiter = '\t';
    S.tables(1).keep_columns = {'gauge_id','gauge_name', ...
        'gauge_name_ru','gauge_lat','gauge_lon','gauge_elev','area_km2'};
    S.tables(2).file = 'static_data.csv';
    S.tables(3).file = 'climate_attributes.csv';
    S.id.column = 'gauge_id';
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.id.sort = 'numeric';
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_components = {'gauge_name'};
    S.metadata.parenthetical_components = {'gauge_name_ru'};
    S.metadata.name_transform = '';
    S.gauge_information.write = false;
    S.region = 'HYDRO_CIS';
    S.zone.region = 'RU';
    S.progress.label = '... Reading HydroCIS attributes';
end

function V = local_variable(source,units)
    V.source = source;
    V.units = units;
    V.target_units = units;
end

function V = local_derived(method,units)
    V.derive = method;
    V.units = units;
    V.target_units = units;
end

function S = local_select_forcing(S,options)
    precipitation = 1;
    temperature = 1;
    pet = 1;
    if isfield(options,'precipitation') && ~isempty(options.precipitation)
        precipitation = options.precipitation;
    end
    if isfield(options,'temperature') && ~isempty(options.temperature)
        temperature = options.temperature;
    end
    if isfield(options,'pet') && ~isempty(options.pet)
        pet = options.pet;
    end

    precipitationSources = { ...
        'prcp_e5l','prcp_e5','prcp_mswep','prcp_gpcp','prcp_imerg'};
    if ~isscalar(precipitation) || precipitation < 1 ...
            || precipitation > numel(precipitationSources)
        error('HYDRO_CIS:BadPrecipitation', ...
            'HydroCIS precipitation choice must be 1 through 5.');
    end
    S.variables.P.source = precipitationSources{precipitation};

    if isequal(temperature,1)
        suffix = 'e5l';
    elseif isequal(temperature,2)
        suffix = 'e5';
    else
        error('HYDRO_CIS:BadTemperature', ...
            'HydroCIS temperature choice must be 1 or 2.');
    end
    S.variables.Tmin.source = ['t_min_' suffix];
    S.variables.Tmax.source = ['t_max_' suffix];

    if isequal(pet,1)
        S.variables.Ep = local_derived('oudin','mm/day');
        S.variables.Ep.valid_min = 0;
    elseif isequal(pet,2)
        S.variables.Ep = local_variable('Ep','mm/day');
        S.variables.Ep.valid_min = 0;
    else
        error('HYDRO_CIS:BadPET', ...
            'HydroCIS potential-evaporation choice must be 1 or 2.');
    end
end
