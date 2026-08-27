function R = region_config_FR()
%REGION_CONFIG_FR Regional defaults for CAMELS-FR.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_FR';
    R.acronym = 'FR';
    R.name = 'France';
    R.dataset = 'CAMELS-FR';
    R.data_root = 'CAMELS_FR';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'FR_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 654;
    X.basins.training = 500;
    X.basins.evaluation = 154;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/1999';
    X.period.manual.train_end = '30/09/2008';
    X.period.manual.eval_start = '01/10/1989';
    X.period.manual.eval_end = '30/09/1999';
    X.period.common_start = '01/10/1989';
    X.period.common_end = '30/09/2008';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-FR [daily]'};
    X.meteo.product.default = 'CAMELS-FR [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Precipitation'};
    X.meteo.precipitation.default = '1 Precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = { ...
        '1 Mean', ...
        '2 (Tmin+Tmax)/2', ...
        '3 Tmin', ...
        '4 Tmax' ...
        };
    X.meteo.temperature.default = '1 Mean';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 Penman-Monteith', ...
        '2 Penman', ...
        '3 Oudin et al.' ...
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


    % Declarative generic-reader schemas owned by this region module.
    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_attribute_schema()
    files = {'CAMELS_FR_station_general_attributes.csv', ...
        'CAMELS_FR_site_general_attributes.csv', ...
        'CAMELS_FR_topography_general_attributes.csv', ...
        'CAMELS_FR_topography_quantiles_attributes.csv', ...
        'CAMELS_FR_soil_general_attributes.csv', ...
        'CAMELS_FR_soil_quantiles_attributes.csv', ...
        'CAMELS_FR_land_cover_attributes.csv', ...
        'CAMELS_FR_hydrogeology_attributes.csv', ...
        'CAMELS_FR_geology_attributes.csv', ...
        'CAMELS_FR_human_influences_dams.csv', ...
        'CAMELS_FR_climate_zone_attributes.csv'};
    S.name = 'CAMELS-FR attributes';
    S.tables = repmat(struct('file','','keys',{{'sta_code_h3'}}, ...
        'delimiter',';','make_valid_names',true, ...
        'make_unique_names',true,'key_target','sta_code_h3', ...
        'join_key','','widen',struct(),'required',true, ...
        'keep_columns',{{}},'column_renames',{{}}, ...
        'duplicate_policy','drop','unique_key',true),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
    end
    S.tables(2).keys = {'sit_code_h3'};
    S.tables(2).key_target = 'sit_code_h3';
    S.tables(2).join_key = 'sit_code_h3';
    S.tables(4).widen.keys = {'top_qnt_quant'};
    S.tables(4).unique_key = false;
    S.tables(5).widen.keys = {'sol_stat','sol_agg_level'};
    S.tables(5).unique_key = false;
    S.tables(6).widen.keys = ...
        {'sol_qnt_agg_level','sol_qnt_quant'};
    S.tables(6).unique_key = false;
    S.tables(11).required = false;
    S.tables(11).delimiter = '';
    S.tables(11).keys = {'gauge_id','sta_code_h3'};
    S.tables(11).column_renames = {'gauge_id','sta_code_h3'};
    S.tables(11).keep_columns = ...
        {'sta_code_h3','aridity','frac_snow','temp_mean'};
    S.tables(11).duplicate_policy = 'replace_left';
    S.id.column = 'sta_code_h3';
    S.id.uppercase = true;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'["'']',''};
    S.id.sort = 'text';
    S.metadata.name_sources = {'sta_label','sit_label'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_FR';
    S.zone.region = 'FR';
    S.progress.label = '... Reading CAMELS-FR generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-FR';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'CAMELS_FR_tsd_{gauge}.csv';
    S.file.delimiter = ';';
    S.file.header_lines = 7;
    S.file.contiguous_time = true;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'tsd_date';
    S.time.input_format = 'yyyyMMdd';
    S.variables.P = local_variable('tsd_prec','mm/day','mm/day');
    S.variables.Q = local_variable('tsd_q_mm','mm/day','mm/day');
    S.aux.tables.site.file = '../../CAMELS_FR_site_general_attributes.csv';
    S.aux.tables.site.key = 'sit_code_h3';
    S.aux.tables.site.elev = 'sit_altitude';
    S.aux.tables.site.area = 'sit_area_topo';
    S.aux.tables.site.area_scale = 1e6;
    S.aux.tables.site.lat_x = 'sit_longitude';
    S.aux.tables.site.lat_y = 'sit_latitude';
    S.aux.tables.site.lat_derive = 'lambert93_to_wgs84';
    S.aux.tables.site.request_prefix_length = 8;
    petVars = {'tsd_pet_pm','tsd_pet_pe','tsd_pet_ou'};
    Base = S;
    for pet = 1:3
        for temp = 1:4
            D = Base;
            D.variables.Ep = local_variable(petVars{pet},'mm/day','mm/day');
            switch temp
                case 1
                    D.variables.T = local_variable('tsd_temp','degC','degC');
                case 2
                    D.variables.Tmin = local_variable('tsd_temp_min','degC','degC');
                    D.variables.Tmax = local_variable('tsd_temp_max','degC','degC');
                    D.variables.T.derive = 'mean_tmin_tmax';
                case 3
                    D.variables.T = local_variable('tsd_temp_min','degC','degC');
                case 4
                    D.variables.T = local_variable('tsd_temp_max','degC','degC');
            end
            name = sprintf('pet%d_t%d',pet,temp);
            petMatch = pet;
            tempMatch = temp;
            if pet == 1
                petMatch = [0 1];
            end
            if temp == 1
                tempMatch = [0 1];
            end
            S.profiles.(name).match = struct('dt',1,'pet',petMatch,'temp',tempMatch);
            S.profiles.(name).schema = D;
        end
    end
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
