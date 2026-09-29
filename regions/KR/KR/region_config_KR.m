function R = region_config_KR()
%REGION_CONFIG_KR CAMELS-KR daily v1.1; independent of hourly CAMELSH-KR.
% Lee, Jung and Ahn (2026), https://doi.org/10.5281/zenodo.21930882
    R.code = 'CAMELS_KR';
    R.acronym = 'KR';
    R.name = 'South Korea (CAMELS) [daily]';
    R.dataset = 'CAMELS-KR';
    R.data_root = 'CAMELS_KR';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'KR_%d_basins.txt';
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 282;
    X.basins.training = 225;
    X.basins.evaluation = 57;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/10/2005';
    X.period.manual.train_end = '30/09/2015';
    X.period.manual.eval_start = '01/10/2015';
    X.period.manual.eval_end = '30/09/2025';
    X.period.common_start = '01/10/2005';
    X.period.common_end = '30/09/2025';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','forcing');
    X.paths.discharge = fullfile('daily','discharge');
    X.meteo.product.items = {'CAMELS-KR [daily]'};
    X.meteo.product.default = 'CAMELS-KR [daily]';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 KMA precipitation'};
    X.meteo.precipitation.default = '1 KMA precipitation';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 Mean temperature'};
    X.meteo.temperature.default = '1 Mean temperature';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 FAO Penman-Monteith PET','2 GLEAM PET'};
    X.meteo.pet.default = '1 FAO Penman-Monteith PET';
    X.meteo.pet.enabled = true;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),9,1);
    C(1).name = 'Daily meteorological forcing';
    C(1).type = 'folder'; C(1).path = X.paths.meteo;
    C(1).pattern = 'CAMELS_KR_Meteorological_timeseries_*.csv';
    C(1).minimum_count = 282;
    C(2).name = 'Daily observed discharge';
    C(2).type = 'folder'; C(2).path = X.paths.discharge;
    C(2).pattern = 'CAMELS_KR_Hydrological_timeseries_*.csv';
    C(2).minimum_count = 282;
    files = {'CAMELS_KR_location_attributes.csv', ...
        'CAMELS_KR_topography_attributes.csv', ...
        'CAMELS_KR_climate_attributes.csv', ...
        'CAMELS_KR_soil_attributes.csv', ...
        'CAMELS_KR_land cover_attributes.csv', ...
        'CAMELS_KR_human influence_attributes.csv', ...
        'CAMELS_KR_hydrology_attributes.csv'};
    for i=1:numel(files)
        C(i+2).name=files{i}; C(i+2).type='file'; C(i+2).path=files{i};
    end
    X.checks=C;
    R.by_resolution.Daily=X;
    S.name='daily CAMELS-KR';
    S.format='csv'; S.layout='one_file_per_basin';
    S.file.pattern='CAMELS_KR_Meteorological_timeseries_{gauge}.csv';
    S.file.delimiter=',';
    S.timeline.reference=datetime(1950,10,1); S.timeline.step=days(1);
    S.time.column='date'; S.time.input_format='yyyy-MM-dd';
    S.id.pad_width=0;
    % Empty source cells remain NaN. Negative fluxes are invalid, not zero.
    S.missing_values=[-999 -9999];
    S.variables.P=variable('prec','mm/day'); S.variables.P.valid_min=0;
    S.variables.Ep=variable('pet','mm/day'); S.variables.Ep.valid_min=0;
    S.variables.T=variable('temp_avg','degC');
    S.resolve=@select_pet;
    S.source_file=[mfilename('fullpath') '.m'];
    R.schema.meteo=S;
    Q=rmfield(S,{'variables','resolve'});
    Q.name='daily observed CAMELS-KR discharge';
    Q.file.pattern='CAMELS_KR_Hydrological_timeseries_{gauge}.csv';
    % Native discharge_spec is already mm/day; never use simulated Q.
    Q.variables.Q=variable('discharge_spec','mm/day');
    Q.variables.Q.valid_min=0;
    R.schema.discharge=Q;
    A.name='CAMELS-KR attributes';
    A.tables=repmat(struct('file','','keys',{{'gauge_id'}}, ...
        'column_renames',{{}},'delimiter',','),numel(files),1);
    for i=1:numel(files), A.tables(i).file=files{i}; end
    A.tables(1).column_renames={'basin_area','area'};
    A.id.column='gauge_id'; A.id.strip=true; A.id.sort='numeric';
    A.metadata.name_sources={'gauge_name'};
    A.metadata.name_transform='';
    A.region='CAMELS_KR'; A.zone.region='KR';
    A.progress.label='... Reading CAMELS-KR daily attributes';
    R.schema.attributes=A;
end
function V=variable(source,units)
    V.source=source; V.units=units; V.target_units=units;
end
function S=select_pet(S,options)
    choice=1;
    if isfield(options,'pet') && ~isempty(options.pet), choice=options.pet; end
    if isequal(choice,2), S.variables.Ep.source='pet_gleam';
    elseif ~isequal(choice,1)
        error('CAMELS_KR:BadPET','CAMELS-KR PET choice must be 1 or 2.');
    end
end
