function R = region_config_NA()
%REGION_CONFIG_NA Regional defaults for daily Namibia GRDC-Caravan.
    R=struct(); R.code='CAMELS_NA'; R.acronym='NA';
    R.name='Namibia'; R.dataset='GRDC-Caravan'; R.data_root='CAMELS_NA';
    R.resolutions={'Daily'}; R.default_resolution='Daily';
    R.basin_file_pattern='NA_%d_basins.txt';
    X=struct(); X.label='Daily'; X.dt=1;
    % Of 51 source stations, 49 have at least 5% observed discharge in
    % both default twenty-water-year periods.
    X.basins.universe=49; X.basins.training=41; X.basins.evaluation=8;
    X.period.spinup_days=365;
    X.period.manual.train_start='01/10/1998'; X.period.manual.train_end='30/09/2018';
    X.period.manual.eval_start='01/10/1978'; X.period.manual.eval_end='30/09/1998';
    X.period.common_start='01/10/1978'; X.period.common_end='30/09/2018';
    X.paths.run_root=''; X.paths.meteo=fullfile('daily','timeseries');
    X.paths.discharge=fullfile('daily','timeseries');
    X.meteo.product.items={'GRDC-Caravan ERA5-Land [daily]'};
    X.meteo.product.default=X.meteo.product.items{1}; X.meteo.product.enabled=false;
    X.meteo.precipitation.items={'1 ERA5-Land precipitation'};
    X.meteo.precipitation.default=X.meteo.precipitation.items{1}; X.meteo.precipitation.enabled=false;
    X.meteo.temperature.items={'1 ERA5-Land temperature'};
    X.meteo.temperature.default=X.meteo.temperature.items{1}; X.meteo.temperature.enabled=false;
    X.meteo.pet.items={'1 FAO Penman-Monteith potential evaporation'};
    X.meteo.pet.default=X.meteo.pet.items{1}; X.meteo.pet.enabled=false;
    C=repmat(struct('name','','type','','path','','pattern','', ...
        'minimum_count',1,'required',true),0,1);
    C(1).name='Daily forcing and streamflow'; C(1).type='folder';
    C(1).path=fullfile('daily','timeseries'); C(1).pattern='*.csv'; C(1).minimum_count=51;
    C(2).name='Caravan attributes'; C(2).type='file'; C(2).path='attributes_caravan_grdc.csv';
    C(3).name='HydroATLAS attributes'; C(3).type='file'; C(3).path='attributes_hydroatlas_grdc.csv';
    C(4).name='Gauge metadata'; C(4).type='file'; C(4).path='attributes_other_grdc.csv';
    C(5).name='Additional attributes'; C(5).type='file'; C(5).path='attributes_additional_grdc.csv';
    X.checks=C; R.by_resolution.Daily=X;

    % Declarative generic-reader schemas owned by this region module.
    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;
    R.schema.attributes = grdc_attribute_schema('NA');

end

function S = local_meteo_schema()
    S.name = 'daily GRDC-Caravan';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'GRDC_{gauge}.csv';
    S.file.delimiter = ',';
    S.id.pad_width = 0;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.column = 'date';
    S.time.input_format = '';
    S.variables.P = local_variable( ...
        'total_precipitation_sum','mm/day','mm/day');
    S.variables.P.clip_min = 0;
    S.variables.Ep = local_variable( ...
        'potential_evaporation_sum_FAO_PENMAN_MONTEITH', ...
        'mm/day','mm/day');
    S.variables.Ep.clip_min = 0;
    S.variables.T = local_variable( ...
        'temperature_2m_mean','degC','degC');
    S.variables.Q = local_variable( ...
        'streamflow','mm/day','mm/day');
    S.progress.label = ...
        '... Reading daily GRDC-Caravan generic data';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
