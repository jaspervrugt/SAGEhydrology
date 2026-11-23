function R = region_config_ET()
%REGION_CONFIG_ET Regional defaults for daily CAMELS-Eth data.
    R=struct('code','CAMELS_ET','acronym','ET','name','Ethiopia', ...
        'dataset','CAMELS-Eth','data_root','CAMELS_ET', ...
        'resolutions',{{'Daily'}},'default_resolution','Daily', ...
        'basin_file_pattern','ET_%d_basins.txt');
    X=struct(); X.label='Daily'; X.dt=1;
    X.basins=struct('universe',83,'training',60,'evaluation',23);
    X.period.spinup_days=365;
    X.period.manual=struct('train_start','01/01/1996','train_end','31/12/1999', ...
        'eval_start','01/01/2000','eval_end','31/12/2003');
    X.period.common_start='01/01/1981'; X.period.common_end='31/12/2018';
    X.paths=struct('run_root','','meteo',fullfile('daily','timeseries'), ...
        'discharge',fullfile('daily','timeseries'));
    X.meteo.product=struct('items',{{'CAMELS-Eth [daily]'}},'default','CAMELS-Eth [daily]','enabled',false);
    X.meteo.precipitation=struct('items',{{'1 CHIRPS'}},'default','1 CHIRPS','enabled',false);
    X.meteo.temperature=struct('items',{{'1 CHIRTS: (Tmin+Tmax)/2'}},'default','1 CHIRTS: (Tmin+Tmax)/2','enabled',false);
    X.meteo.pet=struct('items',{{'1 GLEAM potential evaporation'}},'default','1 GLEAM potential evaporation','enabled',false);
    X.checks=struct('name','Daily time series','type','folder', ...
        'path',fullfile('daily','timeseries'),'pattern','CAMELS_ET_*.csv', ...
        'minimum_count',83,'required',true);
    R.by_resolution.Daily=X;
    R.schema.meteo=local_meteo_schema();
    R.schema.meteo.source_file=[mfilename('fullpath') '.m'];
    R.schema.discharge=derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file=R.schema.meteo.source_file;
    R.schema.attributes=local_attribute_schema();
end
function S=local_attribute_schema()
    files={'CAMELS_ET_topographic_attributes.csv','CAMELS_ET_soil_attributes.csv','CAMELS_ET_landcover_attributes.csv'};
    S.name='CAMELS-Eth physical attributes';
    S.tables=repmat(struct('file','','keys',{{'gauge_id'}},'make_valid_names',true, ...
        'make_unique_names',true,'duplicate_policy','suffix','duplicate_suffix','_et'),numel(files),1);
    for i=1:numel(files), S.tables(i).file=files{i}; end
    S.id=struct('column','gauge_id','uppercase',true,'lowercase',false,'strip',true, ...
        'regex',{{}},'sort','text','output_type','string');
    S.metadata=struct('name_sources',{{'gauge_name'}},'name_transform','');
    S.region='CAMELS_ET'; S.zone.region='ET';
    S.progress.label='... Reading CAMELS-Eth physical attributes';
end
function S=local_meteo_schema()
    S.name='daily CAMELS-Eth'; S.format='csv'; S.layout='one_file_per_basin';
    S.file.pattern='CAMELS_ET_{gauge}.csv'; S.file.delimiter=',';
    S.timeline.reference=datetime(1950,10,1); S.timeline.step=days(1);
    S.time=struct('mode','column','column','date','input_format','yyyy-MM-dd');
    S.variables.P=local_variable('precipitation','mm/day','mm/day');
    S.variables.Tmin=local_variable('tmin','degC','degC');
    S.variables.Tmax=local_variable('tmax','degC','degC'); S.variables.T.derive='mean_tmin_tmax';
    S.variables.Ep=local_variable('pet','mm/day','mm/day');
    S.variables.Q=local_variable('runoff','mm/day','mm/day'); S.variables.Q.valid_min=0;
    S.aux.tables.gauge=struct('file','../../gauge_information.txt','key','gauge_id', ...
        'lat','gauge_lat','area','area_km2','area_scale',1e6);
    S.progress.label='... Reading daily CAMELS-Eth data';
end
function V=local_variable(source,units,targetUnits)
    V=struct('source',source,'units',units,'target_units',targetUnits);
end