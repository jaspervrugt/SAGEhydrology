function R = region_config_IND()
%REGION_CONFIG_IND Regional defaults for CAMELS-IND.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_IND';
    R.acronym = 'IND';
    R.name = 'India';
    R.dataset = 'CAMELS-IND';
    R.data_root = 'CAMELS_IND';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'IND_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 242;
    X.basins.training = 210;
    X.basins.evaluation = 32;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2002';
    X.period.manual.train_end = '30/09/2012';
    X.period.manual.eval_start = '01/10/1992';
    X.period.manual.eval_end = '30/09/2002';
    X.period.common_start = '01/10/1992';
    X.period.common_end = '30/09/2012';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','catchment_mean_forcings');
    X.paths.discharge = fullfile('daily','streamflow_timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-IND [daily]'};
    X.meteo.product.default = 'CAMELS-IND [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = { ...
        '1 Mean', ...
        '2 (Tmin+Tmax)/2', ...
        '3 Tmax', ...
        '4 Tmin' ...
        };
    X.meteo.temperature.default = '1 Mean';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 native', ...
        '2 GLEAM' ...
        };
    X.meteo.pet.default = '1 native';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','catchment_mean_forcings');
    C(1).pattern = '*.csv';
    C(1).minimum_count = 242;
    C(2).name = 'Discharge data';
    C(2).type = 'folder';
    C(2).path = fullfile('daily','streamflow_timeseries');
    C(2).pattern = '*';
    C(2).minimum_count = 1;
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
    S.name = 'CAMELS-IND attributes';
    files = { ...
        'camels_ind_name.csv', ...
        'camels_ind_topo.csv', ...
        'camels_ind_clim.csv', ...
        'camels_ind_hydro.csv', ...
        'camels_ind_land.csv', ...
        'camels_ind_soil.csv', ...
        'camels_ind_geol.csv', ...
        'camels_ind_anth.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true,'key_type','char', ...
        'unique_key',true),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id'};
    end
    S.root_candidates = {'','attributes_csv'};
    S.id.column = 'gauge_id';
    S.id.uppercase = true;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'["'']',''; '^\D*([0-9]+).*$','$1'};
    S.id.numeric_canonical = true;
    S.id.pad_width = 5;
    S.metadata.name_sources = {'cwc_site_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_IND';
    S.zone.region = 'IND';
    S.progress.label = '... Reading CAMELS-IND generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-IND'; 
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = '{gauge}.csv';
    S.id.pad_width = 5;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'ymd_columns'; 
    S.time.year_column = 'year';
    S.time.month_column = 'month'; 
    S.time.day_column = 'day';
    S.variables.P = local_variable('prcp(mm/day)','mm/day','mm/day');
    S.aux.tables.topo.file = '../../camels_ind_topo.csv';
    S.aux.tables.topo.key = 'gauge_id'; 
    S.aux.tables.topo.lat = 'cwc_lat';
    S.aux.tables.topo.pad_width = 5;
    S.aux.tables.topo.elev = 'gauge_elevation';
    S.aux.tables.topo.area = 'cwc_area'; 
    S.aux.tables.topo.area_scale = 1e6;
    Base = S;
    for pet = 1:2
        for temp = 1:4
            D = Base;
            if pet == 1
                epSource = 'pet(mm/day)';
            else
                epSource = 'pet_gleam(mm/day)';
            end
            D.variables.Ep = local_variable(epSource,'mm/day','mm/day');
            switch temp
                case 1
                    D.variables.T = local_variable('tavg(C)','degC','degC');
                case 2
                    D.variables.Tmin = local_variable('tmin(C)','degC','degC');
                    D.variables.Tmax = local_variable('tmax(C)','degC','degC');
                    D.variables.T.derive = 'mean_tmin_tmax';
                case 3
                    D.variables.T = local_variable('tmax(C)','degC','degC');
                case 4
                    D.variables.T = local_variable('tmin(C)','degC','degC');
            end
            petMatch = pet;
            tempMatch = temp;
            if pet == 1
                petMatch = [0 1];
            end
            if temp == 1
                tempMatch = [0 1];
            end
            name = sprintf('pet%d_t%d',pet,temp);
            S.profiles.(name).match = struct('dt',1,'precip',[0 1], ...
                'pet',petMatch,'temp',tempMatch);
            S.profiles.(name).schema = D;
        end
    end
end

function S = local_discharge_schema()
    S.name = 'CAMELS-IND discharge';
    S.format = 'csv';
    S.layout = 'wide_files';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'ymd_columns';
    S.time.year_column = 'year';
    S.time.month_column = 'month';
    S.time.day_column = 'day';
    S.id.strip_leading_zeros = true;
    S.variables.Q.files = {'streamflow_observed.csv'};
    S.variables.Q.units = 'm3/s';
    S.variables.Q.target_units = 'mm/day';
    S.variables.Q.valid_min = 0;
    S.variables.Q.area_normalize = true;
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
