function test_camels_kr(sourceRoot)
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(genpath(root));
cfg=region_helpers('config',root,'CAMELS_KR');
assert(strcmp(cfg.acronym,'KR') && cfg.by_resolution.Daily.dt==1);
assert(strcmp(region_helpers('code','KR'),'CAMELS_KR'));
assert(strcmp(region_helpers('code','KRH'),'CAMELSH_KR'));
assert(strcmp(region_helpers('code','CAMELS-KR'),'CAMELS_KR'));
hourly=region_helpers('config',root,'CAMELSH_KR');
assert(isequal(hourly.resolutions,{'Hourly'}) && hourly.by_resolution.Hourly.basins.universe==178);
[names,codes]=region_helpers('list');
assert(nnz(strcmp(codes,'CAMELS_KR'))==1 && numel(unique(names))==numel(names));
assert(data_helpers('supports','CAMELS_KR','daily'));
assert(~data_helpers('supports','CAMELS_KR','hourly'));
assert(data_helpers('supports','CAMELSH_KR','hourly'));
validate_hydro_schema(cfg.schema.meteo,'meteo');
validate_hydro_schema(cfg.schema.discharge,'Q');
validate_attribute_schema(cfg.schema.attributes);
temporaryRoot=tempname; mkdir(temporaryRoot); cleanup=onCleanup(@()rmdir(temporaryRoot,'s'));
target=fullfile(temporaryRoot,'Data','CAMELS_KR');
missing=struct('region','CAMELS_KR', ...
    'dirD',fullfile(target,'not-installed'),'root',temporaryRoot);
assert(~data_helpers('hasdaily',missing));
assert(data_helpers('install','CAMELS_KR',sourceRoot,target));
assert(data_helpers('hasdaily',struct('region','CAMELS_KR','dirD',target)));
c=struct('region','CAMELS_KR','dirD',target,'root',temporaryRoot);
status=data_helpers('status',c);
assert(status.hasDaily && status.hasMetadata && ~status.supportsHourly);
out=data_helpers('download',c,'daily',struct('confirmFcn',@(varargin)error('Unexpected prompt')));
assert(out.ok && strcmp(out.dirQ,fullfile(target,'daily','discharge')));
catalog=attr_catalog('CAMELS_KR');
assert(~any(ismember(catalog.names(catalog.default_ids),{'q_mean','flow_record'})));
[A,ids,~,zone]=read_attr('CAMELS_KR',target,struct('pr_attr',0));
assert(numel(ids)==282 && size(A,2)==282 && numel(unique(ids))==282);
bas=struct('K',2,'id_gauge',ids([1 282]));
bas.lat=[37.49;35]; bas.elev=[496.562;100]; bas.area=[160.9;100];
split=struct('dt',1,'mode',1,'dt0',datetime(2004,10,1), ...
 'dt_train0',datetime(2005,10,1),'dt_train1',datetime(2015,9,30), ...
 'dt_eval0',datetime(2015,10,1),'dt_eval1',datetime(2025,9,30),'dt_end',datetime(2025,9,30));
split=prepare_split(split);
[dat,aux]=read_meteo('CAMELS_KR',fullfile(target,'daily','forcing'),bas,split,struct('pet',1));
dat=read_Q('CAMELS_KR',fullfile(target,'daily','discharge'),struct(),dat,bas,split,aux);
for k=1:2
 id=string(bas.id_gauge(k));
 M=readtable(fullfile(target,'daily','forcing','CAMELS_KR_Meteorological_timeseries_'+id+'.csv'));
 H=readtable(fullfile(target,'daily','discharge','CAMELS_KR_Hydrological_timeseries_'+id+'.csv'));
 im=M.date>=split.dt0 & M.date<=split.dt_end;
 ih=H.date>=split.dt_train0 & H.date<=split.dt_end;
 assert(isequaln(dat{k}.meteo.P,single(M.prec(im))));
 assert(isequaln(dat{k}.meteo.T,single(M.temp_avg(im))));
 assert(isequaln(dat{k}.meteo.Ep,single(M.pet(im))));
 assert(isequaln(dat{k}.y_n,single(H.discharge_spec(ih))));
end
[gleam,~]=read_meteo('CAMELS_KR',fullfile(target,'daily','forcing'),bas,split,struct('pet',2));
M=readtable(fullfile(target,'daily','forcing','CAMELS_KR_Meteorological_timeseries_'+string(bas.id_gauge(2))+'.csv'));
assert(isequaln(gleam{2}.meteo.Ep,single(M.pet_gleam(M.date>=split.dt0 & M.date<=split.dt_end))));
mdl=struct('tout',numel(split.tout)-1,'id_train',split.id_train,'id_eval',split.id_eval);
eligibility=check_basins(dat,mdl,bas);
assert(numel(eligibility.valid)==2);
fprintf('PASS registry, daily/hourly isolation, schemas, 282-gauge installation/status, attributes, two real basin readers, PET selection and screening.\n');
test_camels_kr_edge();
end

function test_camels_kr_edge()

folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s'));
date=(datetime(2000,1,1):datetime(2000,1,8))';
prec=[0;1;NaN;-999;2;3;4;5]; pet=[1;1;1;1;1;1;1;1];
temp_avg=[-10;0;1;2;3;4;5;6]; pet_gleam=pet*2;
writetable(table(date,prec,pet,temp_avg,pet_gleam),fullfile(folder,'CAMELS_KR_Meteorological_timeseries_1001620.csv'));
discharge_spec=[0;1;NaN;-999;-0.003;2;3;4];
discharge_vol=ones(8,1)*1000; % Must never be substituted for missing Q.
writetable(table(date,discharge_spec,discharge_vol),fullfile(folder,'CAMELS_KR_Hydrological_timeseries_1001620.csv'));
split=struct('dt',1,'mode',1,'dt0',date(1),'dt_train0',date(1), ...
 'dt_train1',date(4),'dt_eval0',date(5),'dt_eval1',date(8),'dt_end',date(8));
bas=struct('K',1,'id_gauge',"1001620");
[dat,aux]=read_meteo('CAMELS_KR',folder,bas,split,struct('pet',1));
dat=read_Q('CAMELS_KR',folder,struct(),dat,bas,split,aux);
assert(dat{1}.meteo.T(1)==-10 && dat{1}.meteo.P(1)==0 && dat{1}.y_n(1)==0);
assert(all(isnan(dat{1}.meteo.P(3:4))) && all(isnan(dat{1}.y_n(3:5))));
assert(isequal(dat{1}.bad,[false false true true true false false false]));
fprintf('PASS missing/sentinel/negative runoff, zero runoff, negative temperature and no Q substitution.\n');
end
