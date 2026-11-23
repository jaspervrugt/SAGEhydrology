function result = test_canonical_cache(cacheBase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEST_CANONICAL_CACHE Check auxiliary inputs and NetCDF cache reuse.
%
% SYNOPSIS:
%   result = test_canonical_cache(cacheBase)
%
% INPUT:
%   cacheBase  OPTIONAL isolated directory for generated test fixtures
%
% OUTPUT:
%   result     Auxiliary-input and source-change cache validation checks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026 Jasper A. Vrugt                                       %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 1, cacheBase = tempname; end
    if ~isfolder(cacheBase), mkdir(cacheBase); end
    raw = fullfile(cacheBase,'raw');
    if ~isfolder(raw), mkdir(raw); end
    start = datetime(2020,1,1);
    Date = string(start+days((0:3).'),'yyyy-MM-dd');
    Q = (1:4).';
    writetable(table(Date,Q),fullfile(raw,'X.csv'));
    split = struct('dt',1,'mode',1,'dt0',start, ...
        'dt_train0',start,'dt_train1',start+days(1), ...
        'dt_eval0',start+days(2),'dt_eval1',start+days(3), ...
        'dt_end',start+days(3));
    split = prepare_split(split);
    bas = struct('id_gauge',"X",'K',1);
    S = struct('name','cache regression','format','csv', ...
        'layout','one_file_per_basin');
    S.file.pattern = '{gauge}.csv';
    S.time.column = 'Date';
    S.time.input_format = 'yyyy-MM-dd';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.variables.Q = struct('source','Q','units','m3/s', ...
        'target_units','mm/day','area_normalize',true);
    S.cache.root = fullfile(cacheBase,'csv_cache');
    aux = [45 200 2e6];
    first = read_discharge_data(raw,struct(),{struct()}, ...
        bas,split,aux,S);
    again = read_discharge_data(raw,struct(),{struct()}, ...
        bas,split,aux,S);
    assert(isequaln(first{1}.obs,again{1}.obs));
    doubled = aux; doubled(3) = 2*aux(3);
    changed = read_discharge_data(raw,struct(),{struct()}, ...
        bas,split,doubled,S);
    S.cache.enabled = false;
    expected = read_discharge_data(raw,struct(),{struct()}, ...
        bas,split,doubled,S);
    result.csv_aux_invalidated = isequaln( ...
        changed{1}.obs.Q.value,expected{1}.obs.Q.value);
    nc = fullfile(raw,'X.nc');
    if isfile(nc), delete(nc); end
    nccreate(nc,'Q','Dimensions',{'time',4});
    ncwrite(nc,'Q',Q);
    S.format = 'netcdf'; S.file.pattern = '{gauge}.nc';
    S.time = struct('mode','index','file_start',start);
    S.cache.enabled = true;
    S.cache.root = fullfile(cacheBase,'nc_cache');
    first = read_discharge_data(raw,struct(),{struct()}, ...
        bas,split,aux,S);
    again = read_discharge_data(raw,struct(),{struct()}, ...
        bas,split,aux,S);
    assert(isequaln(first{1}.obs,again{1}.obs));
    result.netcdf_cache_created = ~isempty(dir( ...
        fullfile(S.cache.root,'*.mat')));
    ncwrite(nc,'Q',2*Q);
    S.cache.enabled = true;
    changed = read_discharge_data(raw,struct(),{struct()}, ...
        bas,split,aux,S);
    S.cache.enabled = false;
    expected = read_discharge_data(raw,struct(),{struct()}, ...
        bas,split,aux,S);
    result.netcdf_source_invalidated = isequaln( ...
        changed{1}.obs.Q.value,expected{1}.obs.Q.value);
    A = (1:4).'; B = (11:14).'; C = (21:24).';
    writetable(table(Date,A,B,C),fullfile(raw,'wide.csv'));
    W = struct('name','wide cache regression','format','csv', ...
        'layout','wide_files');
    W.time.column = 'Date';
    W.time.input_format = 'yyyy-MM-dd';
    W.timeline.reference = datetime(1950,10,1);
    W.timeline.step = days(1);
    W.variables.Q = struct('files',{{'wide.csv'}}, ...
        'units','mm/day','target_units','mm/day');
    W.cache.enabled = true;
    W.cache.root = fullfile(cacheBase,'wide_cache');
    selected = struct('id_gauge',["A";"B"],'K',2);
    empty = repmat({struct()},2,1);
    first = read_discharge_data(raw,struct(),empty,selected,split,[],W);
    selected.id_gauge = ["B";"A"];
    reversed = read_discharge_data( ...
        raw,struct(),empty,selected,split,[],W);
    result.wide_order_preserved = isequaln( first{1}.obs.Q.value, ...
        reversed{2}.obs.Q.value) && ...
        isequaln(first{2}.obs.Q.value,reversed{1}.obs.Q.value);
    selected.id_gauge = ["B";"C"];
    subset = read_discharge_data( ...
        raw,struct(),empty,selected,split,[],W);
    result.wide_selection_preserved = isequal( subset{1}.obs.Q.value, ...
        single(B)) && ...
        isequal(subset{2}.obs.Q.value,single(C));
    result.passed = result.csv_aux_invalidated && ...
        result.netcdf_cache_created && ...
        result.netcdf_source_invalidated && ...
        result.wide_order_preserved && ...
        result.wide_selection_preserved;
    disp(result);
end
