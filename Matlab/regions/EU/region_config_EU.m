function R = region_config_EU()
%REGION_CONFIG_EU Regional defaults for CHAMELLION historical.
%
% All X.paths fields are relative to Data/<R.data_root>.
% Dates use dd/MM/yyyy. Discharge coverage varies by catchment.

    R = struct();
    R.code = 'CHAMELLION';
    R.acronym = 'EU';
    R.name = 'Europe';
    R.dataset = 'CHAMELLION';
    R.data_root = 'CHAMELLION';
    R.resolutions = {'Daily'};
    R.default_resolution = 'Daily';
    R.basin_file_pattern = 'EU_%d_basins.txt';

    X = struct();
    X.label = 'Daily';
    X.dt = 1;
    X.basins.universe = 4560;
    X.basins.file = 'EU_4560_basins.txt';
    X.basins.training = 3648;
    X.basins.evaluation = 912;
    X.period.spinup_days = 365;
    X.period.manual.train_start = '01/01/2001';
    X.period.manual.train_end = '31/12/2019';
    X.period.manual.eval_start = '01/01/1981';
    X.period.manual.eval_end = '31/12/2000';
    X.period.common_start = '01/01/1981';
    X.period.common_end = '31/12/2019';
    X.paths.run_root = '';
    X.paths.meteo = fullfile('daily','timeseries');
    X.paths.discharge = fullfile('daily','timeseries');
    % Meteorological/PET controls shown by SAGE-GUI.
    X.meteo.product.items = {'CHAMELLION'};
    X.meteo.product.default = 'CHAMELLION';
    X.meteo.product.enabled = false;
    X.meteo.precipitation.items = {'1 E-OBS'};
    X.meteo.precipitation.default = '1 E-OBS';
    X.meteo.precipitation.enabled = false;
    X.meteo.temperature.items = {'1 E-OBS'};
    X.meteo.temperature.default = '1 E-OBS';
    X.meteo.temperature.enabled = false;
    X.meteo.pet.items = {'1 Hargreaves'};
    X.meteo.pet.default = '1 Hargreaves';
    X.meteo.pet.enabled = false;
    C = repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name = 'Meteorological data';
    C(1).type = 'folder';
    C(1).path = fullfile('daily','timeseries');
    C(1).pattern = '*';
    C(1).minimum_count = 4565;
    C(1).pattern = 'CHAMELLION_hydromet_timeseries_*.csv';
    C(2)=C(1);C(2).name='Topographic attributes';C(2).type='file';
    C(2).path='CHAMELLION_topographic_attributes.csv';C(2).pattern='';C(2).minimum_count=1;
    C(3)=C(2);C(3).name='Computed Hargreaves climate attributes';
    C(3).path='CHAMELLION_pet_attributes.csv';
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
    S.name='CHAMELLION historical attributes';
    files={'topographic','climatic','hydrologic','landcover','soil'};
    S.tables=repmat(struct('file','','keys',{{'gauge_id'}},'make_valid_names',true),numel(files),1);
    for j=1:numel(files),S.tables(j).file=['CHAMELLION_' files{j} '_attributes.csv'];end
    S.tables(6)=S.tables(1);S.tables(6).file='CHAMELLION_pet_attributes.csv';
    for j=1:6,S.tables(j).required=(j<=5);end
    S.id.column='gauge_id';S.id.sort='text';S.id.uppercase=false;S.id.lowercase=false;
    S.metadata.name_sources={'gauge_name'};S.metadata.name_transform='';
    S.region='CHAMELLION';S.zone.region='EU';
    S.progress.label='... Reading CHAMELLION historical attributes';
end
function S = local_meteo_schema()
    S.name='CHAMELLION historical daily E-OBS, Hargreaves PET';
    S.format='csv';S.layout='one_file_per_basin';
    S.file.pattern='CHAMELLION_hydromet_timeseries_{gauge}.csv';
    S.timeline.reference=datetime(1950,10,1);S.timeline.step=days(1);
    S.time.mode='column';S.time.column='date';S.time.input_format='yyyy-MM-dd';
    S.variables.P=variable('precipitation_mean','mm/day');S.variables.P.clip_min=0;
    S.variables.T=variable('air_temperature_mean','degC');
    S.variables.Tmin=variable('air_temperature_min','degC');
    S.variables.Tmax=variable('air_temperature_max','degC');
    S.variables.Ep.derive='hargreaves';
    % Do not replace missing meteorology/PET with zero.
    S.variables.Q=variable('discharge_spec_obs','mm/day');S.variables.Q.clip_min=0;
    S.aux.tables.gauges.file='../../CHAMELLION_topographic_attributes.csv';
    S.aux.tables.gauges.key='gauge_id';S.aux.tables.gauges.lat='gauge_lat';
    S.aux.tables.gauges.elev='gauge_elev';S.aux.tables.gauges.area='area';
    S.aux.tables.gauges.area_scale=1e6;
end
function V = variable(source,unit)
    V=struct('source',source,'units',unit,'target_units',unit);
end
