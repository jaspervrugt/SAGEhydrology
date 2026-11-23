function result = test_numeric_csv_import(root)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEST_NUMERIC_CSV_IMPORT Check missing prefixes and invalid numeric text.
%
% SYNOPSIS:
%   result = test_numeric_csv_import(root)
%
% INPUT:
%   root       OPTIONAL isolated directory for generated test fixtures
%
% OUTPUT:
%   result     Numeric parsing, missing-value and invalid-text checks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026 Jasper A. Vrugt                                       %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 1, root = tempname; end
    if ~isfolder(root), mkdir(root); end
    file = fullfile(root,'X.csv');
    start = datetime(2020,1,1);
    split = struct('dt',1,'mode',1,'dt0',start, ...
        'dt_train0',start,'dt_train1',start+days(119), ...
        'dt_eval0',start+days(120),'dt_eval1',start+days(123), ...
        'dt_end',start+days(123));
    split = prepare_split(split);
    bas = struct('id_gauge',"X",'K',1);
    S = struct('name','numeric import regression','format','csv', ...
        'layout','one_file_per_basin');
    S.file.pattern = '{gauge}.csv';
    S.time.column = 'Date';
    S.time.input_format = 'yyyy-MM-dd';
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.variables.Q = struct('source','Q','units','mm/day', ...
        'target_units','mm/day');
    S.cache.enabled = false;
    tokens = [repmat("NA",120,1);"0.123456789"; ...
        "nUlL";"MiSsInG";"3.25"];
    clock = string(start+days((0:123).'),'yyyy-MM-dd');
    writetable(table(clock,tokens,'VariableNames',{'Date','Q'}),file);
    dat = read_discharge_data(root,struct(),{struct()}, ...
        bas,split,[],S);
    values = dat{1}.obs.Q.value;
    result.missing_prefix = all(~isfinite(values(1:120)));
    result.mixed_case_missing = all(~isfinite(values(122:123)));
    result.numeric_values = isequal(values([121 124]), ...
        single([0.123456789;3.25]));
    tokens(124) = "unexpected_text";
    writetable(table(clock,tokens,'VariableNames',{'Date','Q'}),file);
    result.invalid_text_rejected = false;
    try
        read_discharge_data(root,struct(),{struct()},bas,split,[],S);
    catch
        result.invalid_text_rejected = true;
    end
    result.passed = all(structfun(@(x)logical(x),result));
    assert(result.passed);
end
