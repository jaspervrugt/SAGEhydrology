function S = read_cartplus_project(C,bands_only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_CARTPLUS_PROJECT Read installed Project 6 inputs in MATLAB.
%   C contains the exported dates, attributes and prescribed basin roles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin > 1 && bands_only
        S = read_bands(C);
        return
    end

    gauges = string(C.gauges(:));
    dates = datetime(C.first):caldays(1):datetime(C.last);
    dates = dates(:);
    if strcmp(C.region,'CAMELS_AU')
        [P,T,E,Q,A,z,w] = read_aus(C,gauges,dates);
        floor = 1e-4;
    else
        [P,T,E,Q,A,z,w] = read_xia(C,gauges,dates);
        floor = 0;
    end

    train = dates >= datetime(C.train_start) ...
        & dates <= datetime(C.train_end);
    eval = dates >= datetime(C.eval_start) ...
        & dates <= datetime(C.eval_end) & C.evaluation_enabled;
    scored = true(size(dates));
    for i = 1:size(C.unscored,1)
        scored = scored & ~(dates >= datetime(C.unscored{i,1}) ...
            & dates <= datetime(C.unscored{i,2}));
    end
    assert(any(train & scored) && (~any(eval) || any(eval & scored)), ...
        'Unscored dates remove an entire scoring period.');
    Q(Q < 0) = NaN;
    eligible = all(isfinite(P) & isfinite(T) & isfinite(E),1) ...
        & quality(Q,train & scored) & ~ismember(gauges,string(C.excluded(:))).';
    if any(eval)
        eligible = eligible & quality(Q,eval & scored);
    end
    if any(~eligible)
        fprintf('Data screening excluded %d basins.\n',nnz(~eligible));
    end
    training = logical(C.training(eligible));
    assert(any(training),'Data screening left no training basins.');
    A = scale_attributes(A(eligible,:),training,C.preprocessing);
    n = nnz(eligible);
    forcing = cat(3,P(:,eligible),T(:,eligible),E(:,eligible), ...
        repmat(train & scored,1,n));
    S = C;
    S.A = A.';
    S.forcing = permute(forcing,[1 3 2]);
    S.observed = Q(:,eligible);
    S.elevation = z(:,eligible);
    S.area = w(:,eligible);
    S.floor = repmat(floor,n,1);
    S.training = training;
    S.train_mask = repmat(train & scored,1,n);
    S.eval_mask = repmat(eval & scored,1,n);
    S.gauges = cellstr(gauges(eligible));
    S.dates = dates;
end

function bands = read_bands(C)
%READ_BANDS Read model elevation bands without loading time series.

    gauges = string(C.gauges(:));
    folder = C.data_root;
    source = fullfile(folder,'cartplus_terrain_bands.csv');
    if isfile(source)
        t = readtable(source,'TextType','string');
        z = zeros(5,numel(gauges));
        w = z;
        for i = 1:numel(gauges)
            rows = sortrows(t(string(t.basin_id)==gauges(i),:),'band');
            assert(height(rows)==5,'Missing bands for %s.',gauges(i));
            z(:,i) = rows.elevation_m;
            w(:,i) = rows.area_fraction;
        end
    elseif isfile(fullfile(folder,'elev_bands671.npz'))
        temp = tempname;
        mkdir(temp);
        cleanup = onCleanup(@() rmdir(temp,'s')); %#ok<NASGU>
        unzip(fullfile(folder,'elev_bands671.npz'),temp);
        ids = read_npy(temp,'ids');
        cols = indices(pad(gauges,8,'left','0'),ids,'elevation gauges');
        [z,w] = us_bands(temp,cols);
    elseif strcmp(C.region,'CAMELS_AU')
        source = fullfile(folder, ...
            'CAMELS_AUS_Attributes&Indices_MasterTable.csv');
        opts = detectImportOptions(source,'VariableNamingRule','preserve');
        opts = setvartype(opts,'station_id','string');
        t = readtable(source,opts);
        cols = indices(gauges,t.station_id,'elevation gauges');
        z = zeros(5,numel(cols));
        for i = 1:numel(cols)
            j = cols(i);
            z(:,i) = aus_bands(t.elev_min(j),t.elev_mean(j),t.elev_max(j));
        end
        w = repmat(.2,size(z));
    else
        % Standard regional topographic percentiles, including GB and FI.
        files = dir(fullfile(folder,'*topographic_attributes.csv'));
        if isempty(files)
            files = dir(fullfile(folder,'daily','*topographic_attributes.csv'));
        end
        if isempty(files)
            files = dir(fullfile(folder,'hourly','*topographic_attributes.csv'));
        end
        assert(~isempty(files),'No CaRT elevation bands found under %s.', ...
            folder);
        source = fullfile(files(1).folder,files(1).name);
        opts = detectImportOptions(source,'VariableNamingRule','preserve');
        key = 'gauge_id';
        opts = setvartype(opts,key,'string');
        t = readtable(source,opts);
        cols = indices(gauges,t.(key),'elevation gauges');
        low = 10;
        if ismember('elev_5',t.Properties.VariableNames)
            low = 5;
        end
        keys = {'elev_min',sprintf('elev_%d',low),'elev_50', ...
            sprintf('elev_%d',100-low),'elev_max'};
        z = interp1([0 low 50 100-low 100],cummax(t{cols,keys},2).', ...
            [10 30 50 70 90].');
        w = repmat(.2,size(z));
    end
    assert(all(isfinite(z),'all') && all(isfinite(w),'all') && ...
        all(w>=0,'all') && all(abs(sum(w,1)-1)<1e-6), ...
        'Invalid CaRT elevation bands.');
    bands = struct('elevation',z,'area',w);
end

function ok = quality(Q,mask)
%QUALITY Require at least 5 percent coverage and variable discharge.

    q = Q(mask,:);
    q(~isfinite(q)) = NaN;
    ok = sum(isfinite(q),1) >= max(2,.05*nnz(mask)) ...
        & max(q,[],1,'omitnan') > min(q,[],1,'omitnan');
end

function A = scale_attributes(A,training,mode)
%SCALE_ATTRIBUTES Apply the selected standard SAGE preprocessing.

    if strcmp(mode,'matlab')
        assert(all(isfinite(A),'all'),'Selected attributes must be finite.');
        cohort = A;
        spread = std(cohort,0,1);
    else
        A(~isfinite(A)) = NaN;
        fill = median(A(training,:),1,'omitnan');
        assert(all(isfinite(fill)), ...
            'Every attribute requires finite training data.');
        replacement = repmat(fill,size(A,1),1);
        missing = isnan(A);
        A(missing) = replacement(missing);
        cohort = A(training,:);
        spread = std(cohort,1,1);
    end
    spread(spread == 0 | ~isfinite(spread)) = 1;
    A = (A-mean(cohort,1))./spread;
end

function [P,T,E,Q,A,z,w] = read_xia(C,gauges,dates)
%READ_XIA Read Xia's installed US or GB release and its elevation data.

    gb = strcmp(C.region,'CAMELS_GB');
    if gb
        archive = 'camels_gb_cart.npz';
        static = 'static';
        thermal = 'Tm';
    else
        archive = 'camels_us_daymet671_v4full.npz';
        static = 'static_raw';
        thermal = 'Tmax';
    end
    folder = tempname;
    mkdir(folder);
    cleanup = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
    file = fullfile(C.data_root,archive);
    assert(isfile(file),'Install Project 6 data first. Missing %s',file);
    unzip(file,folder);
    ids = string(read_npy(folder,'ids'));
    cols = indices(gauges,ids,'gauges');
    names = read_npy(folder,'static_names');
    values = read_npy(folder,static);
    A = values(cols,indices(string(C.attributes),names,'attributes'));
    days = datetime(read_npy(folder,'days'),'InputFormat','yyyy-MM-dd');
    [found,rows] = ismember(dates,days);
    if gb
        warm = datetime(read_npy(folder,'spin_days'), ...
            'InputFormat','yyyy-MM-dd');
        [is_warm,wrows] = ismember(dates(~found),warm);
        assert(all(is_warm),'Requested dates exceed the GB release.');
    else
        assert(all(found),'Requested dates exceed the US release.');
    end
    P = nan(numel(dates),numel(cols));
    T = P;
    E = P;
    Q = P;
    P(found,:) = subset(folder,'P',rows(found),cols);
    E(found,:) = subset(folder,'PET',rows(found),cols);
    Q(found,:) = subset(folder,'q',rows(found),cols);
    T(found,:) = subset(folder,thermal,rows(found),cols);
    if gb
        T(found,:) = T(found,:)-273.15;
        if any(~found)
            F = read_npy(folder,'spin_forcing');
            P(~found,:) = F(wrows,cols,1);
            T(~found,:) = F(wrows,cols,2)-273.15;
            E(~found,:) = F(wrows,cols,3);
        end
        elev = read_npy(folder,'elev_quantiles');
        z = interp1([0 10 50 90 100],cummax(elev(cols,:),2).', ...
            [10 30 50 70 90].');
        w = repmat(.2,size(z));
    else
        T(found,:) = .5*(T(found,:) ...
            + subset(folder,'Tmin',rows(found),cols));
        unzip(fullfile(C.data_root,'elev_bands671.npz'),folder);
        cols = indices(gauges,read_npy(folder,'ids'),'elevation gauges');
        [z,w] = us_bands(folder,cols);
    end
end

function data = subset(folder,name,rows,cols)
%SUBSET Select dates and gauges from a numeric archive member.

    data = read_npy(folder,name);
    data = data(rows,cols);
end

function [z,w] = us_bands(folder,cols)
%US_BANDS Aggregate elevation classes using Xia's cumulative-area rule.

    heights = read_npy(folder,'band_elev');
    weights = read_npy(folder,'band_frac');
    z = zeros(5,numel(cols));
    w = z;
    for i = 1:numel(cols)
        row = cols(i);
        valid = weights(row,:) > 0;
        [h,order] = sort(heights(row,valid));
        f = weights(row,valid);
        f = f(order);
        assert(~isempty(f) && all(isfinite(h)),'Invalid US bands.');
        g = min(5,1+floor((cumsum(f)-1e-9)/.2));
        area = accumarray(g(:),f(:),[5 1]);
        mass = accumarray(g(:),h(:).*f(:),[5 1]);
        z(area>0,i) = mass(area>0)./area(area>0);
        w(:,i) = area/sum(area);
    end
end

function [P,T,E,Q,A,z,w] = read_aus(C,gauges,dates)
%READ_AUS Read standard CAMELS-AUS daily tables and static attributes.

    root = fullfile(C.data_root,'daily');
    files = C.aus_files;
    P = wide(fullfile(root,files{1}),gauges,dates);
    E = wide(fullfile(root,files{2}),gauges,dates);
    T = .5*(wide(fullfile(root,files{3}),gauges,dates) ...
        + wide(fullfile(root,files{4}),gauges,dates));
    T = double(single(T));
    Q = wide(fullfile(root,'streamflow','streamflow_mmd.csv'),gauges,dates);
    file = fullfile(C.data_root,'CAMELS_AUS_Attributes&Indices_MasterTable.csv');
    opts = detectImportOptions(file,'VariableNamingRule','preserve');
    opts = setvartype(opts,'station_id','string');
    a = readtable(file,opts);
    cols = indices(gauges,a.station_id,'Australian attributes');
    if any(strcmp(C.attributes,'npp_diff'))
        monthly = a{:,compose('npp_%d',1:12)};
        a.npp_diff = max(monthly,[],2)-min(monthly,[],2);
        a.npp_diff(any(~isfinite(monthly),2)) = NaN;
    end
    A = a{cols,C.attributes};
    z = nan(5,numel(cols));
    for i = 1:numel(cols)
        row = cols(i);
        z(:,i) = aus_bands(a.elev_min(row),a.elev_mean(row),a.elev_max(row));
    end
    w = repmat(.2,size(z));
end

function value = wide(file,gauges,dates)
%WIDE Align an installed CAMELS-AUS wide table by day and gauge.

    fid = fopen(file,'rt');
    assert(fid >= 0,'Cannot open %s',file);
    cleanup = onCleanup(@() fclose(fid));
    header = strtrim(split(string(fgetl(fid)),','));
    clear cleanup
    cols = indices(gauges,erase(header,'"'),'Australian forcing gauges');
    x = readmatrix(file,'NumHeaderLines',1);
    days = datetime(x(:,1),x(:,2),x(:,3));
    [found,rows] = ismember(dates,days);
    value = nan(numel(dates),numel(cols));
    value(found,:) = double(single(x(rows(found),cols)));
    value(value <= -99) = NaN;
end

function z = aus_bands(lo,target,hi)
%AUS_BANDS Match Xia's truncated-normal elevation distribution.

    assert(all(isfinite([lo target hi])) && hi > target, ...
        'Invalid Australian elevation statistics.');
    if target <= lo
        lo = target-.5*(hi-target);
    end
    span = hi-lo;
    s = .17*span;
    target = min(max(target,lo+.02*span),hi-.02*span);
    mean_at = @(m) integral(@(x) x.*kernel(x,m,lo,hi,s),lo,hi) ...
        /integral(@(x) kernel(x,m,lo,hi,s),lo,hi);
    m = fzero(@(m) mean_at(m)-target,[lo-20*span hi+20*span]);
    mass = integral(@(x) kernel(x,m,lo,hi,s),lo,hi);
    z = zeros(5,1);
    for i = 1:5
        p = (i-.5)/5;
        z(i) = fzero(@(t) integral(@(x) kernel(x,m,lo,hi,s),lo,t) ...
            /mass-p,[lo hi]);
    end
end

function y = kernel(x,m,lo,hi,s)
%KERNEL Evaluate the normal density relative to its maximum on the range.

    edge = min(max(m,lo),hi);
    y = exp(-.5*((x-edge).*(x+edge-2*m))/s^2);
end

function at = indices(wanted,available,label)
%INDICES Resolve selected names without dropping or reordering basins.

    [found,at] = ismember(wanted,available);
    assert(all(found,'all'),'Missing %s: %s',label, ...
        strjoin(wanted(~found),', '));
end

function value = read_npy(folder,name)
%READ_NPY Read the numeric and Unicode arrays used by Xia's NPZ releases.

    file = fullfile(folder,[name,'.npy']);
    fid = fopen(file,'r','ieee-le');
    assert(fid >= 0,'Cannot read archive member %s',file);
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    magic = fread(fid,6,'uint8').';
    assert(isequal(magic,[147 double('NUMPY')]),'Invalid NPY file.');
    version = fread(fid,2,'uint8');
    if version(1) == 1
        count = fread(fid,1,'uint16');
    else
        count = fread(fid,1,'uint32');
    end
    header = char(fread(fid,count,'uint8').');
    desc = regexp(header,'''descr'': ''([^'']+)''','tokens','once');
    shape = regexp(header,'''shape'': \(([^)]*)\)','tokens','once');
    dims = str2double(regexp(shape{1},'\d+','match'));
    switch desc{1}
        case '<f4'
            value = fread(fid,prod(dims),'single=>double');
        case '<f8'
            value = fread(fid,prod(dims),'double');
        case '<i4'
            value = fread(fid,prod(dims),'int32=>double');
        otherwise
            assert(startsWith(desc{1},'<U'),'Unsupported NPY dtype.');
            width = str2double(desc{1}(3:end));
            chars = fread(fid,[width prod(dims)],'uint32=>double');
            value = strings(prod(dims),1);
            for i = 1:numel(value)
                column = chars(:,i);
                value(i) = string(char(column(column~=0).'));
            end
    end
    assert(numel(value)==prod(dims),'Incomplete NPY array: %s',name);
    if numel(dims) > 1
        if contains(header,'''fortran_order'': True')
            value = reshape(value,dims);
        else
            value = permute(reshape(value,fliplr(dims)),numel(dims):-1:1);
        end
    end
end
