function R = region_config_US()
%REGION_CONFIG_US Regional defaults for CAMELS-US.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_US';
    R.acronym = 'US';
    R.name = 'United States (CAMELS)';
    R.dataset = 'CAMELS-US';
    R.data_root = 'CAMELS_US';
    R.resolutions = {'Daily','Hourly'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'US_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = [531 671];
    X.basins.training = 400;
    X.basins.evaluation = 131;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/1999';
    X.period.manual.train_end = '30/09/2008';
    X.period.manual.eval_start = '01/10/1989';
    X.period.manual.eval_end = '30/09/1999';
    X.period.common_start = '01/10/1989';
    X.period.common_end = '30/09/2008';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','v1p2','forcing');
    X.paths.discharge = fullfile('daily','v1p2','streamflow');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = { ...
        '1 Daymet', ...
        '2 Maurer', ...
        '3 NLDAS' ...
        };
    X.meteo.product.default = '1 Daymet';
    X.meteo.product.enabled = true;
    X.meteo.precipitation.items = {'Controlled by meteo product'};
    X.meteo.precipitation.default = 'Controlled by meteo product';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'Controlled by meteo product'};
    X.meteo.temperature.default = 'Controlled by meteo product';
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
    C(1).path = fullfile('daily','v1p2','forcing');
    C(1).pattern = '*';
    C(1).minimum_count = 1;
    C(2).name = 'Discharge data';
    C(2).type = 'folder';
    C(2).path = fullfile('daily','v1p2','streamflow');
    C(2).pattern = '*';
    C(2).minimum_count = 1;
    X.checks = C;
    R.by_resolution.Daily = X;

    X = struct();
    X.label = 'Hourly';
    X.dt = 24;
    X.basins.universe = 499;
    X.basins.training = 400;
    X.basins.evaluation = 99;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2008';
    X.period.manual.train_end = '30/09/2015';
    X.period.manual.eval_start = '01/10/2000';
    X.period.manual.eval_end = '30/09/2008';
    X.period.common_start = '01/10/2000';
    X.period.common_end = '30/09/2015';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('hourly','forcing');
    X.paths.discharge = fullfile('hourly','streamflow');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'3 NLDAS'};
    X.meteo.product.default = '3 NLDAS';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'Controlled by meteo product'};
    X.meteo.precipitation.default = 'Controlled by meteo product';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'Controlled by meteo product'};
    X.meteo.temperature.default = 'Controlled by meteo product';
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
    C(1).path = fullfile('hourly','forcing');
    C(1).pattern = '*';
    C(1).minimum_count = 1;
    C(2).name = 'Discharge data';
    C(2).type = 'folder';
    C(2).path = fullfile('hourly','streamflow');
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
    files = {'camels_clim.txt','camels_geol.txt', ...
        'camels_hydro.txt','camels_soil.txt', ...
        'camels_topo.txt','camels_vege.txt', ...
        'camels_name_clean.txt'};
    S.name = 'CAMELS-US attributes';
    S.tables = repmat(struct('file','','keys',{{'gauge_id'}}, ...
        'delimiter',';','duplicate_policy','drop'),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
    end
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
    S.selection.use_requested = false;
    S.region = 'CAMELS_US';
    S.zone.region = 'US';
    S.progress.label = '... Reading CAMELS-US generic attributes';
end

function S = local_meteo_schema()
    S.name = 'hourly CAMELS-US';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = '{gauge}_hourly_nldas.csv';
    S.id.pad_width = 8;
    S.timeline.reference = datetime(1950,10,1); 
    S.timeline.step = hours(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd HH:mm:ss';
    S.variables.P = local_variable( ...
        'total_precipitation','mm/hour','mm/hour');
    S.variables.Ep = local_variable( ...
        'potential_evaporation','mm/hour','mm/hour');
    S.variables.Ep.clip_min = 0; 
    S.variables.Ep.fill_missing = 0;
    S.variables.T = local_variable('temperature','degC','degC');
    H = S;
    S.profiles.hourly.match = struct( ...
        'dt',24,'data',[1 2 3],'pet',[1 2 3]);
    S.profiles.hourly.schema = H;
    folders = {'daymet','maurer','nldas'};
    labels = {'cida','maurer','nldas'};
    for data = 1:3
        for pet = 1:3
            D = H;
            D.name = 'daily CAMELS-US';
            D.file.pattern = sprintf('%s/*/{gauge}_lump_%s_forcing_leap.txt', ...
                folders{data},labels{data});
            D.file.delimiter = {' ','\t'}; 
            D.file.header_lines = 3;
            D.timeline.step = days(1);
            D.time.mode = 'ymd_columns'; 
            D.time.year_column = 'Year';
            D.time.month_column = 'Mnth'; 
            D.time.day_column = 'Day';
            D.variables.P = local_variable('prcp(mm/day)','mm/day','mm/day');
            D.variables.Tmin = local_variable('tmin(C)','degC','degC');
            D.variables.Tmax = local_variable('tmax(C)','degC','degC');
            D.variables.Radiation = local_variable('srad(W/m2)','W/m2','W/m2');
            D.variables.VaporPressure = local_variable('vp(Pa)','Pa','Pa');
            if data > 1
                D.variables.P.source = 'PRCP(mm/day)';
                D.variables.Tmin.source = 'Tmin(C)';
                D.variables.Tmax.source = 'Tmax(C)';
                D.variables.Radiation.source = 'SRAD(W/m2)';
                D.variables.VaporPressure.source = 'Vp(Pa)';
            end
            D.variables.Wind2.derive = 'constant'; 
            D.variables.Wind2.value = 2;
            D.variables.T.source = ''; 
            D.variables.T.derive = 'mean_tmin_tmax';
            D.variables.Ep.source = ''; 
            D.variables.Ep.derive = 'fao56_daily';
            D.variables.Ep.method = pet;
            D.aux.file_header.lines = [1 2 3]; 
            D.aux.file_header.scale = [1 1 1];
            name = sprintf('daily_d%d_pet%d',data,pet);
            S.profiles.(name).match = struct('dt',1,'data',data,'pet',pet);
            S.profiles.(name).schema = D;
        end
    end
end

function S = local_discharge_schema()
    S.name = 'CAMELS-US discharge';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.id.pad_width = 8;
    S.timeline.reference = datetime(1950,10,1);
    S.file.pattern = '{gauge}-usgs-hourly.csv';
    S.timeline.step = hours(1);
    S.time.mode = 'column';
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd HH:mm:ss';
    S.variables.Q.sources = {'QObs(mm/h)','QObs_CAMELS(mm/h)'};
    S.variables.Q.source_operation = 'coalesce';
    S.variables.Q.units = 'mm/hour';
    S.variables.Q.target_units = 'mm/hour';
    H = S;
    S.profiles.hourly.match.dt = 24;
    S.profiles.hourly.schema = H;
    D = H;
    D.file.pattern = '*/{gauge}_streamflow_qc*.txt';
    D.file.delimiter = {' ','\t'};
    D.file.header_lines = 0;
    D.file.consecutive_delimiters = 'join';
    D.file.variable_names = ...
        {'gauge_id','year','month','day','Q_cfs','quality'};
    D.file.variable_types = ...
        {'double','double','double','double','double','string'};
    D.timeline.step = days(1);
    D.time.mode = 'ymd_columns';
    D.time.year_column = 'year';
    D.time.month_column = 'month';
    D.time.day_column = 'day';
    D.variables.Q.sources = {};
    D.variables.Q.source = 'Q_cfs';
    D.variables.Q.source_operation = 'single';
    D.variables.Q.units = 'ft3/s';
    D.variables.Q.target_units = 'mm/day';
    D.variables.Q.scale = 1/3.28084^3;
    D.variables.Q.area_normalize = true;
    S.profiles.daily.match.dt = 1;
    S.profiles.daily.schema = D;
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
