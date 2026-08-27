function R = region_config_USH()
%REGION_CONFIG_USH Regional defaults for hourly CAMELSH-US.

    R = struct();
    R.code = 'CAMELSH_US';
    R.acronym = 'USH';
    R.name = 'United States (CAMELSH)';
    R.dataset = 'CAMELSH-US';
    R.data_root = 'CAMELSH_US';
    R.resolutions = {'Hourly'};
    R.default_resolution = 'Hourly';
    R.basin_file_pattern = 'USH_%d_basins.txt';

    X = struct();
    X.label = 'Hourly';
    X.dt = 24;
    % Basins have >=5% observed hourly discharge in both default periods.
    X.basins.universe = 2554;
    X.basins.training = 2000;
    X.basins.evaluation = 554;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/01/2010';
    X.period.manual.train_end = '31/12/2019';
    X.period.manual.eval_start = '01/01/2000';
    X.period.manual.eval_end = '31/12/2009';
    X.period.common_start = '01/01/1999';
    X.period.common_end = '31/12/2019';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('hourly','timeseries');
    X.paths.discharge = fullfile('hourly','timeseries');

    X.meteo.product.items = {'CAMELSH-US NLDAS-2 [hourly]'};
    X.meteo.product.default = 'CAMELSH-US NLDAS-2 [hourly]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'NLDAS-2 precipitation'};
    X.meteo.precipitation.default = 'NLDAS-2 precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'NLDAS-2 air temperature'};
    X.meteo.temperature.default = 'NLDAS-2 air temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'NLDAS-2 potential evaporation'};
    X.meteo.pet.default = 'NLDAS-2 potential evaporation';
    X.meteo.pet.enabled = false;

    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Observed hourly forcing and streamflow';
    C(1).type = 'folder';
    C(1).path = fullfile('hourly','timeseries');
    C(1).pattern = '*.nc';
    C(1).minimum_count = 3166;
    C(2).name = 'NLDAS-2 climate attributes';
    C(2).type = 'file';
    C(2).path = fullfile('attributes','attributes_nldas2_climate.csv');
    C(3).name = 'GAGES-II basin identifiers';
    C(3).type = 'file';
    C(3).path = fullfile('attributes','attributes_gageii_BasinID.csv');
    C(4).name = 'GAGES-II topography attributes';
    C(4).type = 'file';
    C(4).path = fullfile('attributes','attributes_gageii_Topo.csv');
    C(5).name = 'CAMELSH streamflow availability';
    C(5).type = 'file';
    C(5).path = 'info.csv';
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
    S.name = 'CAMELSH-US attributes';
    files = { ...
        'attributes_gageii_BasinID.csv', ...
        'attributes_nldas2_climate.csv', ...
        'attributes_gageii_Climate.csv', ...
        'attributes_gageii_Bas_Morph.csv', ...
        'attributes_gageii_Topo.csv', ...
        'attributes_gageii_Soils.csv', ...
        'attributes_gageii_Geology.csv', ...
        'attributes_gageii_Hydro.csv', ...
        'attributes_gageii_Landscape_Pat.csv', ...
        'attributes_gageii_LC06_Basin.csv', ...
        'attributes_gageii_HydroMod_Dams.csv', ...
        'attributes_gageii_HydroMod_Other.csv', ...
        'attributes_gageii_Nutrient_App.csv', ...
        'attributes_gageii_Pest_App.csv', ...
        'attributes_gageii_Pop_Infrastr.csv', ...
        'attributes_gageii_Prot_Areas.csv', ...
        'attributes_hydroATLAS.csv'};
    S.tables = repmat(struct('file','','keys',{{}}, ...
        'make_valid_names',false,'key_type','string'),numel(files),1);
    for i = 1:numel(files)
        S.tables(i).file = fullfile('attributes',files{i});
        S.tables(i).keys = {'STAID'};
    end
    S.id.column = 'STAID';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.id.pad_width = 8;
    S.metadata.name_sources = {'STANAME'};
    S.metadata.name_transform = '';
    S.region = 'CAMELSH_US';
    S.zone.region = 'USH';
    definitions = { ...
        'gauge_lat',{'LAT_GAGE'}; ...
        'gauge_lon',{'LNG_GAGE'}; ...
        'area',{'DRAIN_SQKM'}; ...
        'elev_mean',{'ELEV_MEAN_M_BASIN'}; ...
        'p_seasonality',{'p_seasonality'}; ...
        'frac_snow',{'frac_snow'}; ...
        'aridity',{'aridity_index'}; ...
        'pet_mean',{'pet_mean'}};
    S.aliases = repmat(struct('target','','sources',{{}}, ...
        'required',true,'default',NaN),size(definitions,1),1);
    for i = 1:size(definitions,1)
        S.aliases(i).target = definitions{i,1};
        S.aliases(i).sources = definitions{i,2};
    end
    S.progress.label = '... Reading CAMELSH-US generic attributes';
end

function S = local_meteo_schema()
    S.name = 'hourly CAMELSH-US NLDAS-2';
    S.format = 'netcdf';
    S.layout = 'one_file_per_basin';
    S.file.pattern = '{gauge}.nc';
    S.id.pad_width = 8;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = hours(1);
    S.time.mode = 'index';
    S.time.file_start = datetime(1980,1,1);
    S.variables.P = local_variable('Rainf','mm/hour','mm/hour');
    S.variables.P.valid_min = 0;
    S.variables.Ep = local_variable('PotEvap','mm/hour','mm/hour');
    S.variables.Ep.clip_min = 0;
    S.variables.T = local_variable('Tair','degC','degC');
    S.aux.tables.identity.file = ...
        '../../attributes/attributes_gageii_BasinID.csv';
    S.aux.tables.identity.key = 'STAID';
    S.aux.tables.identity.pad_width = 8;
    S.aux.tables.identity.lat = 'LAT_GAGE';
    S.aux.tables.identity.area = 'DRAIN_SQKM';
    S.aux.tables.identity.area_scale = 1e6;
    S.aux.tables.topography.file = ...
        '../../attributes/attributes_gageii_Topo.csv';
    S.aux.tables.topography.key = 'STAID';
    S.aux.tables.topography.pad_width = 8;
    S.aux.tables.topography.elev = 'ELEV_MEAN_M_BASIN';
    S.progress.label = ...
        '... Reading hourly CAMELSH-US generic data';
end

function S = local_discharge_schema()
    S.name = 'hourly CAMELSH-US discharge';
    S.format = 'netcdf';
    S.layout = 'one_file_per_basin';
    S.file.pattern = '{gauge}.nc';
    S.id.pad_width = 8;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = hours(1);
    S.time.mode = 'index';
    S.time.file_start = datetime(1980,1,1);
    S.variables.Q = local_variable('Streamflow','m3/s','mm/hour');
    S.variables.Q.invalid_le = -9990;
    S.variables.Q.valid_min = 0;
    S.variables.Q.area_normalize = true;
    S.progress.label = '... Reading hourly CAMELSH-US generic discharge';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
