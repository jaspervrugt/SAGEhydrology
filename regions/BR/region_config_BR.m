function R = region_config_BR()
%REGION_CONFIG_BR Regional defaults for CAMELS-BR.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy and water-year boundaries.

    R = struct();
    R.code = 'CAMELS_BR';
    R.acronym = 'BR';
    R.name = 'Brazil';
    R.dataset = 'CAMELS-BR';
    R.data_root = 'CAMELS_BR';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'BR_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 685;
    X.basins.training = 500;
    X.basins.evaluation = 185;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2007';
    X.period.manual.train_end = '30/09/2017';
    X.period.manual.eval_start = '01/10/1996';
    X.period.manual.eval_end = '30/09/2007';
    X.period.common_start = '01/10/1996';
    X.period.common_end = '30/09/2017';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily');
    X.paths.discharge = fullfile('daily','streamflow');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CAMELS-BR [daily]'};
    X.meteo.product.default = 'CAMELS-BR [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = { ...
        '0 Mean of all products', ...
        '1 BR-GWGD', ...
        '2 CHIRPS', ...
        '3 CPC', ...
        '4 ERA5-Land', ...
        '5 MSWEP' ...
        };
    X.meteo.precipitation.default = '0 Mean of all products';
    X.meteo.precipitation.enabled = true;
    X.meteo.temperature.items = { ...
        '1 ERA5-Land', ...
        '2 CPC: (Tmin+Tmax)/2', ...
        '3 ERA5-Land: (Tmin+Tmax)/2', ...
        '4 BR-DWGD: (Tmin+Tmax)/2' ...
        };
    X.meteo.temperature.default = '1 ERA5-Land';
    X.meteo.temperature.enabled = true;
    X.meteo.pet.items = { ...
        '0 Mean of both products', ...
        '1 GLEAM', ...
        '2 ERA5-Land' ...
        };
    X.meteo.pet.default = '0 Mean of both products';
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
    C(3).name = 'PET';
    C(3).type = 'folder';
    C(3).path = fullfile('daily','pet');
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
    files = { ...
        'camels_br_location.txt', ...
        'camels_br_topography.txt', ...
        'camels_br_geology.txt', ...
        'camels_br_climate.txt', ...
        'camels_br_soil.txt', ...
        'camels_br_land_cover.txt', ...
        'camels_br_hydrology.txt', ...
        'camels_br_human_intervention.txt', ...
        'camels_br_quality_check.txt'};
    S.name = 'CAMELS-BR attributes';
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',true,'duplicate_policy','suffix', ...
        'duplicate_suffix','_br'),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = files{i};
        S.tables(i).keys = {'gauge_id'};
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
    S.selection.default_file = 'BR_685_basins.txt';
    S.region = 'CAMELS_BR';
    S.zone.region = 'BR';
    S.progress.label = '... Reading CAMELS-BR generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily CAMELS-BR';
    S.format = 'csv';
    S.layout = 'multi_file_per_basin';
    S.files.precip.pattern = ...
        'precipitation/{gauge}_precipitation.txt';
    S.files.precip.delimiter = ' ';
    S.files.precip.header_lines = 1;
    S.files.precip.contiguous_time = true;
    S.files.precip.variable_names = { ...
        'year','month','day','p_brdwgd','p_chirps', ...
        'p_cpc','p_era5land','p_mswep'};
    S.files.precip.variable_types = repmat({'double'},1,8);
    S.files.pet.pattern = ...
        'pet/{gauge}_potential_evapotransp.txt';
    S.files.pet.delimiter = ' ';
    S.files.pet.header_lines = 1;
    S.files.pet.contiguous_time = true;
    S.files.pet.variable_names = { ...
        'year','month','day','pet_gleam','pet_era5land'};
    S.files.pet.variable_types = repmat({'double'},1,5);
    S.files.temperature.pattern = ...
        'temperature/{gauge}_temperature.txt';
    S.files.temperature.delimiter = ' ';
    S.files.temperature.header_lines = 1;
    S.files.temperature.contiguous_time = true;
    S.files.temperature.variable_names = { ...
        'year','month','day','tmax_cpc','tmax_era5land', ...
        'tmax_brdwgd','tmean_era5land','tmin_cpc', ...
        'tmin_era5land','tmin_brdwgd'};
    S.files.temperature.variable_types = repmat({'double'},1,10);
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'ymd_columns';
    S.time.year_column = 'year';
    S.time.month_column = 'month';
    S.time.day_column = 'day';
    S.variables.P.sources = { ...
        'p_brdwgd','p_chirps','p_cpc','p_era5land','p_mswep'};
    S.variables.P.source_operation = 'mean';
    S.variables.P.units = 'mm/day';
    S.variables.P.target_units = 'mm/day';
    S.variables.P.file = 'precip';
    S.variables.Ep.sources = {'pet_gleam','pet_era5land'};
    S.variables.Ep.source_operation = 'mean';
    S.variables.Ep.units = 'mm/day';
    S.variables.Ep.target_units = 'mm/day';
    S.variables.Ep.file = 'pet';
    S.variables.T.sources = {'tmean_era5land'};
    S.variables.T.source_operation = 'coalesce';
    S.variables.T.units = 'degC';
    S.variables.T.target_units = 'degC';
    S.variables.T.file = 'temperature';
    S.aux.tables.gauges.file = '../gauge_information.txt';
    S.aux.tables.gauges.key = 'GAGE_ID';
    S.aux.tables.gauges.lat = 'LAT';
    S.aux.tables.gauges.area = 'DRAINAGE AREA (KM^2)';
    S.aux.tables.gauges.area_scale = 1e6;
    S.progress.label = ...
        '... Reading daily CAMELS-BR generic multi-file data';

    pSources = { ...
        S.variables.P.sources, ...
        {'p_brdwgd'}, ...
        {'p_chirps'}, ...
        {'p_cpc'}, ...
        {'p_era5land'}, ...
        {'p_mswep'}};
    petSources = { ...
        S.variables.Ep.sources, ...
        {'pet_gleam'}, ...
        {'pet_era5land'}};
    tempSources = { ...
        {'tmean_era5land'}, ...
        {'tmax_cpc','tmin_cpc'}, ...
        {'tmax_era5land','tmin_era5land'}, ...
        {'tmax_brdwgd','tmin_brdwgd'}};

    for precip = 0:5
        for pet = 0:2
            for temp = 1:4
                D.variables.P.sources = pSources{precip + 1};
                D.variables.P.source_operation = local_operation(precip);
                D.variables.Ep.sources = petSources{pet + 1};
                D.variables.Ep.source_operation = local_operation(pet);
                D.variables.T.sources = tempSources{temp};
                if temp == 1
                    D.variables.T.source_operation = 'coalesce';
                else
                    D.variables.T.source_operation = 'mean';
                end
                name = sprintf('p%d_pet%d_t%d',precip,pet,temp);
                S.profiles.(name).match = struct( ...
                    'dt',1,'precip',precip,'pet',pet,'temp',temp);
                S.profiles.(name).schema = D;
            end
        end
    end
end

function operation = local_operation(selection)
    if selection == 0
        operation = 'mean';
    else
        operation = 'coalesce';
    end
end

function S = local_discharge_schema()
    S.name = 'daily CAMELS-BR discharge';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = '{gauge}_streamflow.txt';
    S.file.delimiter = ' ';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'ymd_columns';
    S.time.year_column = 'year';
    S.time.month_column = 'month';
    S.time.day_column = 'day';
    S.variables.Q = local_variable('streamflow_m3s','m3/s','mm/day');
    S.variables.Q.valid_min = 0;
    S.variables.Q.area_normalize = true;
    S.variables.Q.quality.sources = ...
        {'qual_control_by_ana','qual_flag'};
    S.variables.Q.quality.valid_min = realmin;
    S.progress.label = '... Reading daily CAMELS-BR generic discharge';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
