function T = audit_pt_coverage(sourceDir,mapFile,outFile)
%AUDIT_PT_COVERAGE Summarize available daily SNIRH discharge by basin.
M = readtable(mapFile,'TextType','string');
n = height(M);
nValid = zeros(n,1); firstDate = NaT(n,1); lastDate = NaT(n,1);
nCommon = zeros(n,1);
commonA = datetime(1990,10,1); commonB = datetime(2020,9,30);
for i = 1:n
    stationFile = replace(M.gauge_id(i),'/','_') + ".csv";
    f = fullfile(sourceDir,stationFile);
    if ~isfile(f), continue, end
    try
        A = readtable(f,'Delimiter',',','NumHeaderLines',3, ...
            'ReadVariableNames',false,'TextType','string');
        if width(A) < 2, continue, end
        d = datetime(string(A{:,1}),'InputFormat','dd/MM/yyyy HH:mm');
        q = str2double(replace(string(A{:,2}),',','.'));
        good = ~isnat(d) & isfinite(q) & q >= 0;
        d = d(good);
        if isempty(d), continue, end
        [d,ia] = unique(d,'stable'); q = q(good); q = q(ia); %#ok<NASGU>
        nValid(i) = numel(d); firstDate(i)=min(d); lastDate(i)=max(d);
        nCommon(i)=sum(d>=commonA & d<=commonB);
    catch
    end
end
T = table(M.basin_id,M.gauge_code,nValid,firstDate,lastDate,nCommon, ...
    'VariableNames',{'basin_id','gauge_code','valid_days','first_date','last_date','common_days'});
writetable(T,outFile);
end
