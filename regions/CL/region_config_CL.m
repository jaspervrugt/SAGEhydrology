function R = region_config_CL()
%REGION_CONFIG_CL Regional defaults for CAMELS-CL.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_CL';
    R.acronym = 'CL';
    R.name = 'Chile';
    R.dataset = 'CAMELS-CL';
    R.data_root = 'CAMELS_CL';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'CL_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 516;
    X.basins.training = 434;
    X.basins.evaluation = 82;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2006';
    X.period.manual.train_end = '30/09/2014';
    X.period.manual.eval_start = '01/10/1998';
    X.period.manual.eval_end = '30/09/2006';
    X.period.common_start = '01/10/1998';
    X.period.common_end = '30/09/2014';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily');
    X.paths.discharge = fullfile('daily','streamflow');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-Chile [daily]'};
    X.meteo.product.default = 'CAMELS-Chile [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '1 CR2MET', ...
        '2 CHIRPS', ...
        '3 MSWEP', ...
        '4 TMPA' ...
        };
    X.meteo.precipitation.default = '1 CR2MET';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 CR2MET', ...
        '2 CR2MET (Tmin+Tmax)/2' ...
        };
    X.meteo.temperature.default = '1 CR2MET';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '1 Hargreaves', ...
        '2 MODIS 8-day' ...
        };
    X.meteo.pet.default = '1 Hargreaves';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Precipitation files';
    C(1).type = 'pattern';
    C(1).path = fullfile('daily','precipitation');
    C(1).pattern = '*.txt';
    C(1).minimum_count = 4;
    C(2).name = 'Temperature files';
    C(2).type = 'pattern';
    C(2).path = fullfile('daily','temperature');
    C(2).pattern = '*.txt';
    C(2).minimum_count = 3;
    C(3).name = 'PET files';
    C(3).type = 'pattern';
    C(3).path = fullfile('daily','pet');
    C(3).pattern = '*.txt';
    C(3).minimum_count = 2;
    C(4).name = 'Streamflow files';
    C(4).type = 'pattern';
    C(4).path = fullfile('daily','streamflow');
    C(4).pattern = '*.txt';
    C(4).minimum_count = 1;
    C(5).name = 'SWE files';
    C(5).type = 'pattern';
    C(5).path = fullfile('daily','swe');
    C(5).pattern = '*.txt';
    C(5).minimum_count = 1;
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
    S.name = 'CAMELS-CL attributes';
    S.tables.file = '1_CAMELScl_attributes.txt';
    S.tables.keys = {'gauge_id'};
    S.tables.layout = 'transposed_attributes';
    S.tables.delimiter = sprintf('\t');
    S.id.column = 'gauge_id';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'["'']',''; '\.0+$',''};
    S.id.sort = 'text';
    S.metadata.name_sources = {'gauge_name'};
    S.metadata.name_transform = '';
    S.region = 'CAMELS_CL';
    S.zone.region = 'CL';
    S.progress.label = '... Reading CAMELS-CL generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-CL'; 
    S.format = 'csv'; 
    S.layout = 'wide_files';
    S.file.delimiter = sprintf('\t');
    S.timeline.reference = datetime(1950,10,1); 
    S.timeline.step = days(1);
    S.time.mode = 'column'; 
    S.time.column = 'gauge_id';
    S.time.input_format = 'yyyy-MM-dd';
    pFiles = {'precipitation/4_CAMELScl_precip_cr2met.txt', ...
        'precipitation/5_CAMELScl_precip_chirps.txt', ...
        'precipitation/6_CAMELScl_precip_mswep.txt', ...
        'precipitation/7_CAMELScl_precip_tmpa.txt'};
    epFiles = {'pet/12_CAMELScl_pet_hargreaves.txt', ...
        'pet/11_CAMELScl_pet_8d_modis.txt'};
    S.variables.P.files = pFiles(1); 
    S.variables.P.units = 'mm/day';
    S.variables.P.target_units = 'mm/day'; 
    S.variables.P.invalid_le = -99;
    S.variables.Ep.files = epFiles(1); 
    S.variables.Ep.units = 'mm/day';
    S.variables.Ep.target_units = 'mm/day'; 
    S.variables.Ep.invalid_le = -99;
    S.variables.T.files = {'temperature/10_CAMELScl_tmean_cr2met.txt'};
    S.variables.T.units = 'degC'; 
    S.variables.T.target_units = 'degC';
    S.variables.T.invalid_le = -99;
    S.aux.tables.gauge.file = '../gauge_information.txt';
    S.aux.tables.gauge.key = 'gauge_id'; 
    S.aux.tables.gauge.lat = 'gauge_lat';
    Base = S;
    for data = 1:4
        for pet = 1:2
            for temp = 1:2
                D = Base;
                D.variables.P.files = pFiles(data);
                D.variables.Ep.files = epFiles(pet);
                if pet == 2
                    D.variables.Ep.time_mode = 'row_sequence_from_first';
                end
                if temp == 2
                    D.variables.T.files = {};
                    D.variables.T.derive = 'mean_tmin_tmax';
                    D.variables.Tmin.files = { ...
                        'temperature/8_CAMELScl_tmin_cr2met.txt'};
                    D.variables.Tmin.invalid_le = -99;
                    D.variables.Tmax.files = { ...
                        'temperature/9_CAMELScl_tmax_cr2met.txt'};
                    D.variables.Tmax.invalid_le = -99;
                end
                name = sprintf('p%d_pet%d_t%d',data,pet,temp);
                S.profiles.(name).match = struct( ...
                    'dt',1,'precip',data,'pet',pet,'temp',temp);
                S.profiles.(name).schema = D;
            end
        end
    end
end

function S = local_discharge_schema()
    S.name = 'CAMELS-CL discharge';
    S.format = 'csv';
    S.layout = 'wide_files';
    S.file.delimiter = sprintf('\t');
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'column';
    S.time.column = 'gauge_id';
    S.time.input_format = 'yyyy-MM-dd';
    S.variables.Q.files = {'3_CAMELScl_streamflow_mm.txt'};
    S.variables.Q.units = 'mm/day';
    S.variables.Q.target_units = 'mm/day';
    S.variables.Q.valid_min = 0;
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
