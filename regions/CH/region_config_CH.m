function R = region_config_CH()
%REGION_CONFIG_CH Regional defaults for CAMELS-CH.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_CH';
    R.acronym = 'CH';
    R.name = 'Switzerland';
    R.dataset = 'CAMELS-CH';
    R.data_root = 'CAMELS_CH';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'CH_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 221;
    X.basins.training = 180;
    X.basins.evaluation = 41;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2006';
    X.period.manual.train_end = '30/09/2014';
    X.period.manual.eval_start = '01/10/1998';
    X.period.manual.eval_end = '30/09/2006';
    X.period.common_start = '01/10/1998';
    X.period.common_end = '30/09/2014';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-CH [daily]'};
    X.meteo.product.default = 'CAMELS-CH [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 Observation-based', ...
        '2 Simulation-based' ...
        };
    X.meteo.precipitation.default = '1 Observation-based';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 Observed mean', ...
        '2 (Tmin+Tmax)/2' ...
        };
    X.meteo.temperature.default = '1 Observed mean';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = {'1 Simulation-based'};
    X.meteo.pet.default = '1 Simulation-based';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = '*';
    C(1).minimum_count = 1;
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
    S.name = 'CAMELS-CH attributes';
    definitions = { ...
        'CAMELS_CH_topographic_attributes.csv','ISO-8859-1'; ...
        'CAMELS_CH_climate_attributes_obs.csv','UTF-8'; ...
        'CAMELS_CH_hydrology_attributes_obs.csv','UTF-8'; ...
        'CAMELS_CH_hydrogeology_attributes.csv','UTF-8'; ...
        'CAMELS_CH_soil_attributes.csv','UTF-8'; ...
        'CAMELS_CH_landcover_attributes.csv','UTF-8'; ...
        'CAMELS_CH_geology_attributes.csv','UTF-8'; ...
        'CAMELS_CH_glacier_attributes.csv','UTF-8'; ...
        'CAMELS_CH_humaninfluence_attributes.csv','UTF-8'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true,'join_type','full', ...
        'delimiter',',','encoding','','comment_style','#'), ...
        size(definitions,1),1);
    for i = 1:size(definitions,1)
        S.tables(i).file = definitions{i,1};
        S.tables(i).keys = {'gauge_id'};
        S.tables(i).encoding = definitions{i,2};
    end
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'"',''; '\.0+$',''};
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_CH';
    S.zone.region = 'CH';
    S.progress.label = '... Reading CAMELS-CH generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-CH'; 
    S.format = 'csv';
    S.layout = 'multi_file_per_basin';
    S.id.pad_width = 4;
    S.timeline.reference = datetime(1950,10,1); 
    S.timeline.step = days(1);
    S.time.mode = 'column'; 
    S.time.column = 'date';
    S.time.input_format = 'yyyy-MM-dd';
    S.files.obs.pattern =  ...
        'observation_based/CAMELS_CH_obs_based_{gauge}.csv';
    S.files.sim.pattern =  ...
        'simulation_based/CAMELS_CH_sim_based_{gauge}.csv';
    S.missing_values = [-999 -99 -99.9 -99.99];
    S.variables.Q = local_variable('discharge_spec(mm/d)','mm/day','mm/day');
    S.variables.Q.file = 'obs';
    S.variables.Ep = local_variable('pet_sim(mm/d)','mm/day','mm/day');
    S.variables.Ep.file = 'sim'; 
    S.variables.Ep.clip_min = 0;
    S.variables.Ep.fill_missing = 0;
    S.variables.P = local_variable( ...
        'precipitation(mm/d)','mm/day','mm/day');
    S.variables.P.file = 'obs';
    S.variables.T = local_variable( ...
        'temperature_mean(degC)','degC','degC');
    S.variables.T.file = 'obs';
    S.aux.tables.gauge.file = '../../gauge_information.txt';
    S.aux.tables.gauge.key = 'gauge_id';
    S.aux.tables.gauge.lat = 'gauge_lat';
    Base = S;
    for data = 1:2
        for temp = 1:2
            D = Base;
            if data == 1
                D.variables.P = local_variable( ...
                    'precipitation(mm/d)','mm/day','mm/day');
                D.variables.P.file = 'obs';
            else
                D.variables.P = local_variable( ...
                    'precipitation_sim(mm/d)','mm/day','mm/day');
                D.variables.P.file = 'sim';
            end
            if temp == 1
                D.variables.T = local_variable( ...
                    'temperature_mean(degC)','degC','degC');
                D.variables.T.file = 'obs';
            else
                D.variables = rmfield(D.variables,'T');
                D.variables.Tmin = local_variable( ...
                    'temperature_min(degC)','degC','degC');
                D.variables.Tmin.file = 'obs';
                D.variables.Tmax = local_variable( ...
                    'temperature_max(degC)','degC','degC');
                D.variables.Tmax.file = 'obs';
                % Empty overrides remove the base profile's direct source
                % when the selected profile is recursively merged.
                D.variables.T.source = ''; 
                D.variables.T.file = '';
                D.variables.T.derive = 'mean_tmin_tmax';
            end
            name = sprintf('p%d_t%d',data,temp);
            S.profiles.(name).match = struct( ...
                'dt',1,'precip',data,'temp',temp,'pet',[0 1]);
            S.profiles.(name).schema = D;
        end
    end
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
