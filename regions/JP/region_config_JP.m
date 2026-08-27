function R = region_config_JP()
%REGION_CONFIG_JP Regional defaults for MERV-Jp version 2.0.
%
% Time-series and HydroATLAS attribute support use the 87-basin version-2
% release. Dates reserve one complete year for model spin-up.

    R = struct();
    R.code = 'CAMELS_JP';
    R.acronym = 'JP';
    R.name = 'Japan';
    R.dataset = 'MERV-Jp 2.0';
    R.data_root = 'CAMELS_JP';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'JP_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 87;
    X.basins.file = 'JP_87_basins.txt';
    X.basins.training = 70;
    X.basins.evaluation = 17;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2002';
    X.period.manual.train_end = '30/09/2015';
    X.period.manual.eval_start = '01/10/1990';
    X.period.manual.eval_end = '30/09/2002';
    X.period.common_start = '01/10/1990';
    X.period.common_end = '30/09/2015';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');

    X.meteo.product.items = {'MERV-Jp 2.0 [daily]'};
    X.meteo.product.default = 'MERV-Jp 2.0 [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 Basin precipitation'};
    X.meteo.precipitation.default = '1 Basin precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 Mean air temperature'};
    X.meteo.temperature.default = '1 Mean air temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 Potential evapotranspiration'};
    X.meteo.pet.default = '1 Potential evapotranspiration';
    X.meteo.pet.enabled = false;

    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'MERV-Jp 2.0 forcing and streamflow';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = 'varssim*.csv';
    C(1).minimum_count = 87;
    C(2).name = 'MERV-Jp attributes';
    C(2).type = 'file';
    C(2).path = 'MERV_Jp_135_HydroATLAS_attributes.xlsx';
    C(2).required = true;
    X.checks = C;
    R.by_resolution.Daily = X;

    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = local_attribute_schema();

end

function S = local_attribute_schema()
    S.name = 'MERV-Jp HydroATLAS attributes';
    S.tables = struct( ...
        'file','MERV_Jp_135_HydroATLAS_attributes.xlsx', ...
        'sheet','HydroATLAS_all','keys',{{'basin_no'}}, ...
        'make_valid_names',true,'required',true);
    S.root_candidates = {'','attributes'};
    S.id.column = 'basin_no';
    S.id.uppercase = false;
    S.id.lowercase = false;
    S.id.strip = true;
    S.id.regex = {'\.0+$',''};
    S.id.numeric_canonical = true;
    S.id.sort = 'numeric';
    S.metadata.name_sources = {'station'};
    S.metadata.name_components = {'station','river'};
    S.metadata.parenthetical_components = {'station_jp','river_jp'};
    S.metadata.name_separator = ': ';
    S.metadata.name_fallback_prefix = 'MERV-Jp basin ';
    S.metadata.name_transform = 'title_case_words';
    S.metadata.standardize = false;
    S.aliases = repmat(struct('target','','sources',{{}}, ...
        'required',true,'default',NaN,'scale',[], ...
        'divisor',[]),7,1);
    S.aliases(1).target = 'gauge_lat';
    S.aliases(1).sources = {'outlet_lat'};
    S.aliases(2).target = 'gauge_lon';
    S.aliases(2).sources = {'outlet_lon'};
    S.aliases(3).target = 'gauge_elev';
    S.aliases(3).sources = {'ele_mt_sav'};
    S.aliases(4).target = 'area';
    S.aliases(4).sources = {'area_merv_km2'};
    S.aliases(5).target = 'p_mean';
    S.aliases(5).sources = {'pre_mm_syr'};
    S.aliases(6).target = 'pet_mean';
    S.aliases(6).sources = {'pet_mm_syr'};
    S.aliases(7).target = 'frac_snow';
    S.aliases(7).sources = {'snw_pc_syr'};
    S.aliases(7).divisor = 100;
    S.selection.available_pattern = fullfile( ...
        'daily','timeseries','varssim*.csv');
    S.selection.available_regex = '^varssim(\d+)\.csv$';
    S.region = 'CAMELS_JP';
    S.zone.region = 'JP';
    S.progress.label = '... Reading MERV-Jp generic attributes';
end

function S = local_meteo_schema()
    S.name = 'daily MERV-Jp 2.0';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'varssim{gauge}.csv';
    S.file.delimiter = ',';
    S.file.contiguous_time = true;
    S.id.pad_width = 3;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.mode = 'ymd_columns';
    S.time.year_column = 'Year';
    S.time.month_column = 'Month';
    S.time.day_column = 'Day';
    S.variables.P = local_variable('Precip','mm/day','mm/day');
    S.variables.P.clip_min = 0;
    S.variables.T = local_variable('Temp','degC','degC');
    S.variables.Ep = local_variable('PET','mm/day','mm/day');
    S.variables.Ep.clip_min = 0;
    S.variables.Q = local_variable('Obs flow','mm/day','mm/day');
    S.variables.Q.clip_min = 0;
    S.progress.label = '... Reading daily MERV-Jp generic data';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
