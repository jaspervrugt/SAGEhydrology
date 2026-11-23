function R=region_config_GRDC(shortCode,regionName,sourceCount, ...
        modelCount,nTrain,nEval,evalStart,evalEnd,trainStart,trainEnd)
%REGION_CONFIG_GRDC Common configuration for GRDC-Caravan country subsets.
% Modeling universes use at least 5% observed discharge coverage in each
% default training/evaluation period; sourceCount is the installed total.
    shortCode=upper(char(string(shortCode)));
    R=struct(); R.code=['CAMELS_' shortCode]; R.acronym=shortCode;
    R.name=regionName; R.dataset='GRDC-Caravan'; R.data_root=R.code;
    R.resolutions={'Daily'}; R.default_resolution='Daily';
    R.basin_file_pattern=[shortCode '_%d_basins.txt'];
    X=struct(); X.label='Daily'; X.dt=1;
    X.basins.universe=modelCount; X.basins.training=nTrain;
    X.basins.evaluation=nEval;
    X.period.spinup_days=365;
    X.period.manual.train_start=trainStart;
    X.period.manual.train_end=trainEnd;
    X.period.manual.eval_start=evalStart;
    X.period.manual.eval_end=evalEnd;
    X.period.common_start=evalStart;
    X.period.common_end=trainEnd;
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
    C(1).path=fullfile('daily','timeseries'); C(1).pattern='*.csv'; C(1).minimum_count=sourceCount;
    C(2).name='Caravan attributes'; C(2).type='file'; C(2).path='attributes_caravan_grdc.csv';
    C(3).name='HydroATLAS attributes'; C(3).type='file'; C(3).path='attributes_hydroatlas_grdc.csv';
    C(4).name='Gauge metadata'; C(4).type='file'; C(4).path='attributes_other_grdc.csv';
    C(5).name='Additional attributes'; C(5).type='file'; C(5).path='attributes_additional_grdc.csv';
    X.checks=C; R.by_resolution.Daily=X;
    R.schema.attributes = grdc_attribute_schema(shortCode);
end
