function [dat,aux] = read_hydro_timeseries( ...
    dirData,bas,split,options,schema,dat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_HYDRO_TIMESERIES Read schema-described SAGE time-series data.
%
% SYNOPSIS:
%   [dat,aux] = read_hydro_timeseries(dirData,bas,split,options,schema)
%   [dat,aux] = read_hydro_timeseries( ...
%       dirData,bas,split,options,schema,dat)
%
% INPUT:
%   dirData    Directory containing basin time-series files
%   bas        Basin structure with field id_gauge and optional metadata
%   split      Prepared SAGE time-split structure
%   options    Runtime options, including optional progressFcn
%   schema     Validated declarative time-series schema
%   dat        OPTIONAL existing SAGE data cell array to augment
%
% OUTPUT:
%   dat        SAGE basin data with meteo and/or observed discharge
%   aux        K-by-3 matrix: latitude, elevation, and area in m2
%
% DESCRIPTION:
%   The reader supports one or several named CSV/NetCDF files per basin,
%   plus wide CSV files. Each named file is aligned independently to the
%   requested SAGE window. It applies schema-defined unit
%   transformations and quality control, optionally converts volumetric
%   discharge to runoff depth, and returns canonical SAGE fields.
%   Compact canonical caches store only processed hydrologic variables.
%   Full native-table and aggregate caches remain disabled by default.
%
% NOTES:
%   Dataset-specific PET derivations or unusual corrections can be supplied
%   through a variable specification's transform function handle. Existing
%   regional readers remain authoritative until numerical equivalence has
%   been demonstrated for their schemas.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6
        dat = [];
    end
    if nargin < 4 ...
            || isempty(options)
        options = struct();
    end
    if ~isstruct(bas) ...
            || ~isfield(bas,'id_gauge') ...
            || isempty(bas.id_gauge)
        error('read_hydro_timeseries:MissingBasins', ...
            'bas.id_gauge must be provided.');
    end

    split = prepare_split(split);
    schemaMode = 'both';
    if isfield(schema,'mode') ...
            && ~isempty(schema.mode)
        schemaMode = schema.mode;
    end
    schema = validate_hydro_schema(schema,schemaMode);
    ids = local_ids(bas.id_gauge,schema.id);
    K = numel(ids);
    if isempty(dat)
        dat = cell(K,1);
    elseif ~iscell(dat) ...
            || numel(dat) ~= K
        error('read_hydro_timeseries:BadDat', ...
            'Existing dat must be a K-element cell array.');
    end
    dat = dat(:);
    aux = local_aux(schema,bas,ids,dirData);

    request = build_requested_window(split,true);
    offsets = double(request.read_start_N:request.read_end_N).';
    targetTime = schema.timeline.reference ...
        + offsets.*schema.timeline.step;
    nTarget = numel(targetTime);
    scoreIdx = local_score_indices(split,nTarget);
    cacheRequestSignature = local_canonical_request_signature( ...
        targetTime,scoreIdx);

    [aggregateHit,dat,aux,aggregateFile,aggregateSignature] = ...
        local_aggregate_cache('load',dirData,ids,split,options,schema, ...
            targetTime,scoreIdx,dat,aux);
    if aggregateHit
        return
    end

    progress = reader_progress('start',schema.progress.label,K);

    % --------------------------------------------
    % Read source data according to storage layout
    % --------------------------------------------
    streamed = false;
    switch schema.layout
    
        case 'one_file_per_basin'
    
            streamed = true;
            csvPlan = [];
            for k = 1:K
                file = local_file(dirData, ...
                    schema.file.pattern,ids(k));
    
                if ~isfile(file)
                    error('read_hydro_timeseries:MissingFile', ...
                        'Missing time-series file: %s',file);
                end
                sourceFiles = {file};
                [cacheHit,dat{k}] = local_canonical_cache( ...
                    'load',dirData,ids(k),options,schema, ...
                    cacheRequestSignature,sourceFiles,dat{k});
                if ~cacheHit
                    switch schema.format
                        case 'csv'
                            [raw,csvPlan] = local_read_csv( ...
                                file,schema,targetTime,csvPlan);

                        case 'netcdf'
                            raw = local_read_netcdf( ...
                                file,schema,targetTime);

                        otherwise
                            error(['read_hydro_timeseries:' ...
                                'UnsupportedFormat'], ...
                                'Unsupported format: %s.',schema.format);
                    end

                    dat{k} = local_materialize_basin( ...
                        dat{k},raw,sourceFiles,schema,targetTime, ...
                        scoreIdx,k,aux,ids(k));
                    local_canonical_cache( ...
                        'save',dirData,ids(k),options,schema, ...
                        cacheRequestSignature,sourceFiles,dat{k});
                end
                progress = reader_progress('update',progress,k);
                local_progress_callback(options,k,K);
            end

        case 'multi_file_per_basin'

            streamed = true;
            namedPlans = struct();
            for k = 1:K
                sourceFiles = local_named_source_files( ...
                    dirData,ids(k),schema);
                [cacheHit,dat{k}] = local_canonical_cache( ...
                    'load',dirData,ids(k),options,schema, ...
                    cacheRequestSignature,sourceFiles,dat{k});
                if ~cacheHit
                    [raw,sourceFiles,namedPlans] = ...
                        local_read_named_files( ...
                        dirData,ids(k),schema,targetTime,namedPlans);
                    dat{k} = local_materialize_basin( ...
                        dat{k},raw,sourceFiles,schema,targetTime, ...
                        scoreIdx,k,aux,ids(k));
                    local_canonical_cache( ...
                        'save',dirData,ids(k),options,schema, ...
                        cacheRequestSignature,sourceFiles,dat{k});
                end
                progress = reader_progress('update',progress,k);
                local_progress_callback(options,k,K);
            end
    
        case 'wide_files'
    
            if ~strcmp(schema.format,'csv')
                error('read_hydro_timeseries:WideFormat', ...
                    'wide_files currently supports CSV input only.');
            end
    
            [rawAll,wideFiles] = local_read_wide_files( ...
                dirData,ids,targetTime,options,schema);
    
            filesAll = repmat({wideFiles},K,1);
    
        otherwise
            error('read_hydro_timeseries:UnsupportedLayout', ...
                'Unsupported layout: %s.',schema.layout);
    end
    
    % --------------------------------------------------
    % Convert raw source values to canonical SAGE fields
    % --------------------------------------------------
    if ~streamed
        for k = 1:K
            dat{k} = local_materialize_basin(dat{k},rawAll{k},filesAll{k}, ...
                schema,targetTime,scoreIdx,k,aux,ids(k));
            progress = reader_progress('update',progress,k);
            local_progress_callback(options,k,K);
        end
    end
    
    reader_progress('finish',progress,K);

    if isfield(split,'local') ...
            && split.local ...
            && isfield(split,'rainfall_block') ...
            && split.rainfall_block ...
            && any(strcmp(schema.mode,{'meteo','both'}))
        dat = rainfall_rank_split(dat,split);
    end

    local_aggregate_cache('save',dirData,ids,split,options,schema, ...
        targetTime,scoreIdx,dat,aux,aggregateFile,aggregateSignature);
end

function [hit,entry] = local_canonical_cache( ...
    action,dirData,id,options,schema,requestSignature,sourceFiles,entry)

    hit = false;
    eligible = isfield(schema,'cache') ...
        && schema.cache.enabled ...
        && schema.cache.canonical_enabled ...
        && ~any(endsWith(lower(string(sourceFiles)),'.nc'));
    if ~eligible
        return
    end

    profile = '';
    if isfield(schema,'selected_profile')
        profile = char(string(schema.selected_profile));
    end
    choices = struct();
    for name = {'data','pet','temp','precip','precipitation'}
        if isfield(options,name{1})
            choices.(name{1}) = options.(name{1});
        end
    end
    key = strjoin({char(string(schema.name)),profile, ...
        char(string(schema.mode)),jsonencode(choices)},char(30));
    cacheDir = fullfile(dirData,char(string(schema.cache.directory)));
    safeId = matlab.lang.makeValidName(char(string(id)));
    cacheFile = fullfile(cacheDir,sprintf('%s_%s.mat', ...
        safeId,local_signature_hash(key)));
    sourceSignature = local_canonical_source_signature( ...
        sourceFiles,schema);

    if strcmp(action,'load')
        if ~isfile(cacheFile)
            return
        end
        try
            C = load(cacheFile,'cacheVersion','requestSignature', ...
                'sourceSignature','payload');
            valid = isfield(C,'cacheVersion') ...
                && C.cacheVersion == 1 ...
                && isfield(C,'requestSignature') ...
                && strcmp(C.requestSignature,requestSignature) ...
                && isfield(C,'sourceSignature') ...
                && strcmp(C.sourceSignature,sourceSignature) ...
                && isfield(C,'payload');
            if valid
                entry = local_merge_canonical_payload(entry,C.payload);
                hit = true;
            end
        catch
            hit = false;
        end
        return
    end

    if ~strcmp(action,'save')
        error('read_hydro_timeseries:CanonicalCacheAction', ...
            'Unknown canonical cache action: %s.',action);
    end
    if ~isfolder(cacheDir)
        mkdir(cacheDir);
    end
    payload = local_canonical_payload(entry,schema);
    cacheVersion = 1;
    temporaryFile = [tempname(cacheDir) '.mat'];
    cleanup = onCleanup(@()local_delete_cache_file(temporaryFile));
    try
        save(temporaryFile,'cacheVersion','requestSignature', ...
            'sourceSignature','payload');
        [status,message] = movefile(temporaryFile,cacheFile,'f');
        if ~status
            warning('read_hydro_timeseries:CanonicalCacheWrite', ...
                'Could not write canonical cache %s: %s', ...
                cacheFile,message);
        end
    catch ME
        warning('read_hydro_timeseries:CanonicalCacheWrite', ...
            'Could not write canonical cache %s: %s', ...
            cacheFile,ME.message);
    end
    clear cleanup
end

function signature = local_canonical_request_signature(targetTime,scoreIdx)
    time = posixtime(targetTime(:));
    parts = {mat2str(time.'),mat2str(double(scoreIdx(:).'))};
    signature = local_signature_hash(strjoin(parts,char(30)));
end

function signature = local_canonical_source_signature(sourceFiles,schema)
    parts = strings(numel(sourceFiles)+2,1);
    for i = 1:numel(sourceFiles)
        info = dir(sourceFiles{i});
        if isempty(info)
            parts(i) = string(sourceFiles{i}) + "|missing";
        else
            parts(i) = string(sourceFiles{i}) + "|" ...
                + string(info.bytes) + "|" ...
                + compose('%.15g',info.datenum);
        end
    end
    codeInfo = dir(mfilename('fullpath'));
    if ~isempty(codeInfo)
        parts(end-1) = "reader|" + compose('%.15g',codeInfo.datenum);
    end
    if isfield(schema,'source_file') ...
            && ~isempty(schema.source_file)
        schemaInfo = dir(char(string(schema.source_file)));
        if ~isempty(schemaInfo)
            parts(end) = "schema|" + compose('%.15g',schemaInfo.datenum);
        end
    end
    signature = local_signature_hash(strjoin(parts,char(30)));
end

function payload = local_canonical_payload(entry,schema)
    payload = struct();
    if any(strcmp(schema.mode,{'meteo','both'})) ...
            && isfield(entry,'meteo')
        payload.meteo = entry.meteo;
    end
    if any(strcmp(schema.mode,{'Q','both'})) ...
            && isfield(schema.variables,'Q')
        if isfield(entry,'y_n')
            payload.y_n = entry.y_n;
        end
        if isfield(entry,'bad')
            payload.bad = entry.bad;
        end
    end
    if isfield(entry,'gauge')
        payload.gauge = entry.gauge;
    end
    if isfield(entry,'fname')
        payload.fname = entry.fname;
    end
end

function entry = local_merge_canonical_payload(entry,payload)
    if isempty(entry)
        entry = struct();
    end
    names = fieldnames(payload);
    for i = 1:numel(names)
        field = names{i};
        if strcmp(field,'fname') ...
                && isfield(entry,'fname') ...
                && ~isempty(entry.fname)
            files = cellstr(string(payload.fname(:)));
            for j = 1:numel(files)
                if ~any(strcmp(entry.fname,files{j}))
                    entry.fname{end+1} = files{j};
                end
            end
        else
            entry.(field) = payload.(field);
        end
    end
end

function local_delete_cache_file(file)
    if isfile(file)
        delete(file);
    end
end

function [hit,dat,aux,cacheFile,signature] = local_aggregate_cache( ...
        action,dirData,ids,split,options,schema,targetTime,scoreIdx,dat,aux, ...
        cacheFile,signature)
    hit = false;
    if nargin < 11
        cacheFile = '';
        signature = '';
    end
    eligible = isfield(schema,'cache') ...
        && schema.cache.enabled ...
        && schema.cache.aggregate_enabled ...
        && numel(ids) >= double(schema.cache.aggregate_min_basins) ...
        && any(strcmp(schema.layout,{'one_file_per_basin', ...
            'multi_file_per_basin'}));
    if ~eligible
        return
    end
    if strcmp(action,'load')
        signature = local_aggregate_signature(ids,split,options,schema, ...
            targetTime,scoreIdx);
        cacheDir = fullfile(dirData,char(string(schema.cache.directory)));
        cacheFile = fullfile(cacheDir,['_aggregate_' ...
            local_signature_hash(signature) '.mat']);
        if ~isfile(cacheFile)
            return
        end
        cacheInfo = dir(cacheFile);
        if isempty(cacheInfo) ...
            || cacheInfo.datenum < local_source_timestamp( ...
                dirData,schema)
            return
        end
        try
            C = load(cacheFile,'aggregateVersion','signature','dat','aux');
            if isfield(C,'aggregateVersion') ...
                    && C.aggregateVersion == 1 ...
                    && isfield(C,'signature') ...
                    && strcmp(C.signature,signature) ...
                    && isfield(C,'dat') ...
                    && isfield(C,'aux')
                dat = C.dat;
                aux = C.aux;
                hit = true;
            end
        catch
            hit = false;
        end
        return
    end
    if isempty(cacheFile)
        return
    end
    cacheDir = fileparts(cacheFile);
    if ~isfolder(cacheDir)
        mkdir(cacheDir);
    end
    aggregateVersion = 1;
    try
        save(cacheFile,'aggregateVersion','signature','dat','aux','-v7.3');
    catch ME
        warning('read_hydro_timeseries:AggregateCacheWrite', ...
            'Could not write aggregate cache %s: %s',cacheFile,ME.message);
    end
end

function signature = local_aggregate_signature( ...
        ids,split,options,schema,targetTime,scoreIdx)
    profile = '';
    if isfield(schema,'selected_profile')
        profile = schema.selected_profile;
    end
    choice = struct();
    for name = {'data','pet','temp','precipitation'}
        if isfield(options,name{1})
            choice.(name{1}) = options.(name{1});
        end
    end
    rain = false;
    if isfield(split,'local') ...
        && split.local ...
            && isfield(split,'rainfall_block')
        rain = logical(split.rainfall_block);
    end
    parts = {char(string(schema.name)),char(string(profile)), ...
        char(strjoin(string(ids(:)),'|')),char(string(targetTime(1))), ...
        char(string(targetTime(end))),mat2str(double(scoreIdx(:).')), ...
        jsonencode(choice),char(string(rain))};
    signature = strjoin(parts,char(30));
end

function hash = local_signature_hash(value)
    bytes = double(unicode2native(value,'UTF-8'));
    h = 2166136261;
    for i = 1:numel(bytes)
        h = mod((h+bytes(i))*16777619,2^32);
    end
    hash = lower(dec2hex(uint32(h),8));
end

function stamp = local_source_timestamp(dirData,schema)
    stamp = 0;
    patterns = {};
    if strcmp(schema.layout,'one_file_per_basin')
        patterns = {schema.file.pattern};
    elseif isfield(schema,'files')
        names = fieldnames(schema.files);
        for i = 1:numel(names)
            patterns{end+1} = schema.files.(names{i}).pattern; %#ok<AGROW>
        end
    end
    for i = 1:numel(patterns)
        pattern = strrep(char(string(patterns{i})),'{gauge02}','*');
        pattern = strrep(pattern,'{gauge}','*');
        info = dir(fullfile(dirData,pattern));
        if ~isempty(info)
            stamp = max(stamp,max([info.datenum]));
        end
    end
    codeInfo = dir(mfilename('fullpath'));
    if ~isempty(codeInfo)
        stamp = max(stamp,codeInfo.datenum);
    end
    if isfield(schema,'source_file') ...
            && ~isempty(schema.source_file)
        schemaInfo = dir(char(string(schema.source_file)));
    else
        schemaInfo = dir(which('meteo_schema'));
    end
    if ~isempty(schemaInfo)
        stamp = max(stamp,schemaInfo.datenum);
    end
end

function entry = local_materialize_basin(entry,raw,sourceFiles, ...
        schema,targetTime,scoreIdx,k,aux,id)
    if isempty(entry)
        entry = struct();
    end
    if isempty(sourceFiles)
        fileLabel = char(string(schema.name));
    else
        fileLabel = strjoin(sourceFiles,'; ');
    end
    if any(strcmp(schema.mode,{'meteo','both'}))
        P = local_variable(raw,schema.variables.P,'P', ...
            fileLabel,k,aux,schema.timeline.step);
        raw = local_derive_variables(raw,schema,targetTime,k,aux);
        Ep = local_variable(raw,schema.variables.Ep,'Ep', ...
            fileLabel,k,aux,schema.timeline.step);
        T = local_variable(raw,schema.variables.T,'T', ...
            fileLabel,k,aux,schema.timeline.step);
        entry.meteo.P = single(P);
        entry.meteo.Ep = single(Ep);
        entry.meteo.T = single(T);
        entry.meteo.bad = find(~isfinite(P) | ~isfinite(Ep) | ~isfinite(T)).';
    end
    if isfield(schema.variables,'Q') ...
        && ...
            any(strcmp(schema.mode,{'Q','both'}))
        Q = local_variable(raw,schema.variables.Q,'Q', ...
            fileLabel,k,aux,schema.timeline.step);
        entry.y_n = single(Q(scoreIdx));
        entry.bad = (~isfinite(Q(scoreIdx)) | Q(scoreIdx) < 0).';
    end
    entry.gauge = id;
    if ~isempty(sourceFiles)
        if ~isfield(entry,'fname') ...
            || isempty(entry.fname)
            entry.fname = sourceFiles;
        else
            for j = 1:numel(sourceFiles)
                if ~any(strcmp(entry.fname,sourceFiles{j}))
                    entry.fname{end+1} = sourceFiles{j};
                end
            end
        end
    end
end

function local_progress_callback(options,k,~)
    if isfield(options,'progressFcn') ...
            && ~isempty(options.progressFcn)
        try
            options.progressFcn(k);
        catch MEprogress
            warning('read_hydro_timeseries:ProgressFcn', ...
                'Progress callback failed: %s',MEprogress.message);
        end
    end
end

function [raw,plan] = local_read_csv(file,schema,targetTime,plan)
    if nargin < 4
        plan = [];
    end
    if ~isempty(plan)
        opts = plan.opts;
        sourceMap = plan.sourceMap;
    elseif isfield(schema.file,'variable_names') ...
        && ...
            ~isempty(schema.file.variable_names)
        names = cellstr(string(schema.file.variable_names));
        types = repmat({'double'},size(names));
        if isfield(schema.file,'variable_types') ...
            && ...
                ~isempty(schema.file.variable_types)
            types = cellstr(string(schema.file.variable_types));
        end
        opts = delimitedTextImportOptions( ...
            'Delimiter',schema.file.delimiter, ...
            'VariableNames',names,'VariableTypes',types, ...
            'VariableNamingRule','preserve');
        if isfield(schema.file,'consecutive_delimiters') ...
            && ...
                strcmpi(string(schema.file.consecutive_delimiters),'join')
            opts.ConsecutiveDelimitersRule = 'join';
        end
        opts.DataLines = [double(schema.file.header_lines)+1 Inf];
    else
        opts = detectImportOptions(file,'VariableNamingRule','preserve', ...
            'Delimiter',schema.file.delimiter, ...
            'NumHeaderLines',schema.file.header_lines);
    end
    if isempty(plan)
        available = opts.VariableNames;
        [sourceMap,needed] = local_resolve_sources(schema,available,file);
        switch schema.time.mode

        case 'column'
            needed = unique( ...
                [{char(schema.time.column)}; needed], ...
                'stable');

        case {'ymd_columns','ymd_row_sequence'}
            components = {char(schema.time.year_column); ...
                char(schema.time.month_column); ...
                char(schema.time.day_column)};
            for component = {'hour_column','minute_column','second_column'}
                if isfield(schema.time,component{1}) ...
                        && ~isempty(schema.time.(component{1}))
                    components{end+1,1} = ...
                        char(schema.time.(component{1})); %#ok<AGROW>
                end
            end
            needed = unique([components; ...
                needed], ...
                'stable');
        end

        missing = setdiff(needed,available);
        if ~isempty(missing)
            error('read_hydro_timeseries:MissingColumns', ...
                'Missing columns in %s: %s',file,strjoin(missing,', '));
        end
        opts.SelectedVariableNames = needed;
        firstDataLine = double(opts.DataLines(1));
        plan = struct('opts',opts,'sourceMap',sourceMap, ...
            'firstDataLine',firstDataLine);
    end
    [T,cacheUsed] = local_read_csv_native_cache( ...
        file,opts,schema,targetTime);
    readOpts = opts;
    if ~cacheUsed ...
        && isfield(schema.file,'contiguous_time') ...
            && schema.file.contiguous_time
        probeOpts = opts;
        probeOpts.DataLines = [plan.firstDataLine plan.firstDataLine];
        probe = readtable(file,probeOpts);
        sourceStart = local_csv_table_time(probe,schema.time,file);
        stepSeconds = seconds(schema.timeline.step);
        offset = round(seconds(targetTime(1)-sourceStart(1))/stepSeconds);
        if offset < 0
            error('read_hydro_timeseries:PeriodOutsideFile', ...
                'Requested period starts before %s.',file);
        end
        firstLine = plan.firstDataLine+offset;
        readOpts.DataLines = [firstLine firstLine+numel(targetTime)-1];
    end
    if ~cacheUsed
        T = readtable(file,readOpts);
    end
    alignedRows = isfield(schema.file,'contiguous_time') ...
        && schema.file.contiguous_time ...
        && height(T) == numel(targetTime);
    if alignedRows
        % Row slicing and native-cache indexing already used the declared
        % contiguous timeline to select exactly this window. Avoid parsing
        % and joining the same timestamps again for every basin.
        found = true(numel(targetTime),1);
        location = (1:numel(targetTime)).';
    else
        sourceTime = local_csv_table_time(T,schema.time,file);
        if schema.strict_time ...
                && any(diff(sourceTime) <= seconds(0))
            error('read_hydro_timeseries:NonIncreasingTime', ...
                'Timestamps are not strictly increasing in %s.',file);
        end
        if isfield(schema.time,'join') ...
                && strcmp(schema.time.join,'calendar_day')
            [found,location] = ismember( ...
                dateshift(targetTime,'start','day'), ...
                dateshift(sourceTime,'start','day'));
        else
            [found,location] = ismember(targetTime,sourceTime);
        end
    end
    raw = struct();
    raw.found = found;
    fields = fieldnames(schema.variables);
    for i = 1:numel(fields)
        field = fields{i};
        spec = schema.variables.(field);
        % Derived variables do not have a source column.
        hasSource = (isfield(spec,'source') ...
            && ~isempty(spec.source)) ...
            || (isfield(spec,'sources') ...
            && ~isempty(spec.sources));
        if ~hasSource
            continue
        end
        values = nan(numel(targetTime),1);
        columns = sourceMap.(field);
        X = nan(height(T),numel(columns));
        for j = 1:numel(columns)
            source = columns{j};
            X(:,j) = local_numeric_column(T.(source),source,file);
            X(:,j) = local_apply_native_missing(X(:,j),spec);
        end
        column = local_combine_sources(X,spec,field,file);
        if ~isempty(spec.quality)
            column(~local_csv_quality_mask(T,spec.quality,file)) = NaN;
        end
        values(found) = column(location(found));
        raw.(field) = values;
    end

end

function [T,used] = local_read_csv_native_cache( ...
        file,opts,schema,targetTime)
    T = table();
    used = false;
    if ~isfield(schema,'cache') ...
        || ~schema.cache.enabled ...
        || ~schema.cache.native_enabled
        return
    end
    info = dir(file);
    if isempty(info)
        return
    end
    cacheDir = fullfile(fileparts(file),char(string(schema.cache.directory)));
    [~,base,~] = fileparts(file);
    cacheFile = fullfile(cacheDir,[base '_native.mat']);
    requested = string(opts.SelectedVariableNames(:));
    returned = requested;
    if isfield(schema.file,'contiguous_time') ...
            && schema.file.contiguous_time ...
            && strcmp(schema.time.mode,'column')
        returned = requested(requested ~= string(schema.time.column));
    end
    cachedNames = strings(0,1);
    cacheVariables = strings(0,1);
    valid = false;
    if isfile(cacheFile)
        try
            % Cache column variable names are deterministic, so metadata
            % and the requested columns can be fetched with one MAT-file
            % open.  This matters for datasets containing hundreds of
            % basin files, especially on synchronized/network folders.
            wantedVariables = cellstr(local_cache_variable_names(returned));
            loadNames = [{'cacheVersion','cacheNames','cacheVariables', ...
                'cacheStart','cacheStep','cacheCount','cacheTime', ...
                'sourceBytes','sourceDatenum'},wantedVariables];
            missingVariableWarning = warning( ...
                'off','MATLAB:load:variableNotFound');
            restoreWarning = onCleanup(@() warning(missingVariableWarning));
            C = load(cacheFile,loadNames{:});
            valid = isfield(C,'cacheVersion') ...
                && any(C.cacheVersion == [2 3]) ...
                && isfield(C,'sourceBytes') ...
                && C.sourceBytes == info.bytes ...
                && isfield(C,'sourceDatenum') ...
                && abs(C.sourceDatenum-info.datenum) < 1e-8 ...
                && isfield(C,'cacheNames') ...
                && isfield(C,'cacheVariables') ...
                && isfield(C,'cacheStart') ...
                && isfield(C,'cacheStep') ...
                && isfield(C,'cacheCount') ...
                && isfield(C,'cacheTime');
            if valid
                cachedNames = string(C.cacheNames(:));
                cacheVariables = string(C.cacheVariables(:));
            end
        catch
            % A requested column may not exist yet when a runtime profile
            % augments an otherwise current cache.  Pay for a metadata-only
            % retry on that uncommon path; complete hits remain one load.
            try
                C = load(cacheFile,'cacheNames', ...
                    'cacheVariables','cacheStart','cacheStep','cacheCount', ...
                    'cacheTime','sourceBytes','sourceDatenum');
                versionOkay = ~isfield(C,'cacheVersion') ...
                    || any(C.cacheVersion == [2 3]);
                valid = versionOkay ...
                    && isfield(C,'sourceBytes') ...
                    && C.sourceBytes == info.bytes ...
                    && isfield(C,'sourceDatenum') ...
                    && abs(C.sourceDatenum-info.datenum) < 1e-8 ...
                    && all(isfield(C,{'cacheNames','cacheVariables', ...
                        'cacheStart','cacheStep','cacheCount','cacheTime'}));
                if valid
                    cachedNames = string(C.cacheNames(:));
                    cacheVariables = string(C.cacheVariables(:));
                    if all(ismember(returned,cachedNames))
                        [~,oldLoc] = ismember(returned,cachedNames);
                        oldVariables = cellstr(cacheVariables(oldLoc));
                        V = load(cacheFile,oldVariables{:});
                        for q = 1:numel(oldVariables)
                            C.(oldVariables{q}) = V.(oldVariables{q});
                        end
                    end
                end
            catch
                valid = false;
            end
        end
    end
    hitVariables = {};
    if valid ...
        && all(ismember(returned,cachedNames))
        [~,hitLoc] = ismember(returned,cachedNames);
        hitVariables = cellstr(cacheVariables(hitLoc));
    end
    completeHit = valid ...
        && ~isempty(hitVariables) ...
        && all(isfield(C,hitVariables));
    if completeHit
        [~,loc] = ismember(returned,cachedNames);
        rows = local_cached_time_rows('',C,targetTime,schema);
        variableNames = cellstr(cacheVariables(loc));
        values = cell(1,numel(loc));
        for j = 1:numel(loc)
            values{j} = C.(variableNames{j})(rows,:);
        end
        T = table(values{:},'VariableNames',cellstr(returned));
        used = true;
        return
    end
    importNames = requested;
    if valid
        importNames = unique([cachedNames; requested],'stable');
    end
    available = string(opts.VariableNames);
    importNames = importNames(ismember(importNames,available));
    fullOpts = opts;
    fullOpts.SelectedVariableNames = cellstr(importNames);
    % Native caches always cover the full source record. Requested-window
    % row slicing remains the uncached fast path and is used when caching
    % is explicitly disabled.
    fullOpts.DataLines = [double(opts.DataLines(1)) Inf];
    nativeTable = readtable(file,fullOpts);
    nativeTime = local_csv_table_time(nativeTable,schema.time,file);
    cacheTime = posixtime(nativeTime);
    cacheNames = string(nativeTable.Properties.VariableNames);
    cacheVariables = local_cache_variable_names(cacheNames);
    payload = struct();
    for j = 1:width(nativeTable)
        payload.(char(cacheVariables(j))) = ...
            nativeTable.(nativeTable.Properties.VariableNames{j});
    end
    cacheStart = cacheTime(1);
    cacheCount = numel(cacheTime);
    if cacheCount > 1
        cacheStep = median(diff(cacheTime));
    else
        cacheStep = NaN;
    end
    sourceBytes = info.bytes;
    sourceDatenum = info.datenum;
    cacheVersion = 3;
    if ~isfolder(cacheDir)
        mkdir(cacheDir);
    end
    payload.cacheNames = cacheNames;
    payload.cacheVariables = cacheVariables;
    payload.cacheTime = cacheTime;
    payload.cacheStart = cacheStart;
    payload.cacheStep = cacheStep;
    payload.cacheCount = cacheCount;
    payload.sourceBytes = sourceBytes;
    payload.sourceDatenum = sourceDatenum;
    payload.cacheVersion = cacheVersion;
    % Version 7.3 permits true selected-variable reads.  Version 7 MAT
    % files decompress the complete payload even when LOAD names only a
    % few columns, which is prohibitively slow for per-basin caches.
    save(cacheFile,'-struct','payload','-v7.3');
    C = payload;
    rows = local_cached_time_rows(cacheFile,C,targetTime,schema);
    T = nativeTable(rows,cellstr(requested));
    used = true;
end

function names = local_cache_variable_names(sourceNames)
    names = "column_"+string(matlab.lang.makeValidName( ...
        cellstr(string(sourceNames(:)))));
end

function rows = local_cached_time_rows(cacheFile,C,targetTime,schema)
    target = posixtime(targetTime(:));
    contiguous = isfield(schema.file,'contiguous_time') ...
        && schema.file.contiguous_time ...
        && isfinite(C.cacheStep);
    if contiguous
        rows = round((target-double(C.cacheStart))/double(C.cacheStep))+1;
        rows = unique(rows(rows >= 1 ...
            & rows <= double(C.cacheCount)),'stable');
        return
    end
    if isfield(C,'cacheTime')
        cacheTime = double(C.cacheTime(:));
    else
        V = load(cacheFile,'cacheTime');
        cacheTime = double(V.cacheTime(:));
    end
    if isfield(schema.time,'join') ...
        && strcmp(schema.time.join,'calendar_day')
        cacheKey = floor(cacheTime/86400);
        targetKey = floor(target/86400);
    else
        cacheKey = cacheTime;
        targetKey = target;
    end
    [found,loc] = ismember(targetKey,cacheKey);
    rows = unique(loc(found),'stable');
end

function sourceTime = local_csv_table_time(T,timeSpec,file)
    switch timeSpec.mode
        case 'column'
            sourceTime = local_datetime(T.(char(timeSpec.column)),timeSpec);
        case {'ymd_columns','ymd_row_sequence'}
            yr = local_numeric_column(T.(char(timeSpec.year_column)), ...
                char(timeSpec.year_column),file);
            mo = local_numeric_column(T.(char(timeSpec.month_column)), ...
                char(timeSpec.month_column),file);
            dy = local_numeric_column(T.(char(timeSpec.day_column)), ...
                char(timeSpec.day_column),file);
            if strcmp(timeSpec.mode,'ymd_row_sequence')
                sourceTime = local_ymd_row_sequence(yr,mo,dy,timeSpec,file);
            else
                hh = local_optional_time_component(T,timeSpec,'hour_column',file);
                mm = local_optional_time_component(T,timeSpec,'minute_column',file);
                ss = local_optional_time_component(T,timeSpec,'second_column',file);
                sourceTime = datetime(yr,mo,dy,hh,mm,ss);
            end
    end
end

function time = local_ymd_row_sequence(yr,mo,dy,spec,file)
    day = dateshift(datetime(yr,mo,dy),'start','day');
    n = numel(day);
    time = NaT(n,1);
    first = [true; diff(day) ~= days(0)];
    blockStart = find(first);
    blockEnd = [blockStart(2:end)-1; n];
    steps = double(spec.steps_per_day);
    step = days(1)/steps;
    for b = 1:numel(blockStart)
        idx = blockStart(b):blockEnd(b);
        count = numel(idx);
        if count > steps
            error('read_hydro_timeseries:TooManyRowsPerDay', ...
                '%s has %d rows on %s; schema permits %d.', ...
                file,count,string(day(idx(1))),steps);
        end
        startSlot = 0;
        if b == 1 ...
                && strcmp(spec.partial_first,'final_slots')
            startSlot = steps-count;
        end
        time(idx) = day(idx) + (startSlot:startSlot+count-1)'.*step;
    end
end

function value = local_optional_time_component(T,timeSpec,field,file)
    value = zeros(height(T),1);
    if isfield(timeSpec,field) ...
            && ~isempty(timeSpec.(field))
        name = char(timeSpec.(field));
        value = local_numeric_column(T.(name),name,file);
    end
end

function [sourceMap,sources] = local_resolve_sources( ...
    schema,available,file)
    sourceMap = struct();
    sources = cell(0,1);
    fields = fieldnames(schema.variables);
    for i = 1:numel(fields)
        field = fields{i};
        spec = schema.variables.(field);
        if isfield(spec,'source') ...
                && ~isempty(spec.source)
            requested = {char(string(spec.source))};
            operation = 'single';
        elseif isfield(spec,'sources') ...
                && ~isempty(spec.sources)
            requested = cellstr(string(spec.sources(:)));
            operation = char(string(spec.source_operation));
        else
            continue
        end
        present = requested(ismember(requested,available));
        if strcmp(operation,'single') ...
                && isempty(present)
            error('read_hydro_timeseries:MissingColumns', ...
                'Missing column in %s: %s',file,requested{1});
        elseif strcmp(operation,'aliases')
            if isempty(present)
                error(['read_hydro_timeseries:' ...
                    'MissingSourceAlternatives'], ...
                    ['None of the ordered aliases for ' ...
                    '%s exist in %s: %s'], ...
                    field,file,strjoin(requested,', '));
            end
            present = present(1);
        elseif isempty(present)
            error('read_hydro_timeseries:MissingSourceAlternatives', ...
                'None of the source columns for %s exist in %s: %s', ...
                field,file,strjoin(requested,', '));
        end
        sourceMap.(field) = present;
        sources = [sources; present(:)]; %#ok<AGROW>
        if isfield(spec,'quality') ...
                && ~isempty(spec.quality)
            qualitySources = cellstr(string(spec.quality.sources(:)));
            missingQuality = setdiff(qualitySources,available);
            if ~isempty(missingQuality)
                error('read_hydro_timeseries:MissingQualityColumns', ...
                    'Missing quality columns in %s: %s', ...
                    file,strjoin(missingQuality,', '));
            end
            sources = [sources; qualitySources(:)]; %#ok<AGROW>
        end
    end
    sources = unique(sources,'stable');
end

function valid = local_csv_quality_mask(T,quality,file)
    valid = true(height(T),1);
    sources = cellstr(string(quality.sources(:)));
    for j = 1:numel(sources)
        x = local_numeric_column(T.(sources{j}),sources{j},file);
        if quality.require_finite
            valid = valid ...
                & isfinite(x);
        end
        valid = valid & x >= double(quality.valid_min) ...
            & x <= double(quality.valid_max);
    end
end

function value = local_combine_sources(X,spec,field,file)
    if size(X,2) == 1 ...
            || ~isfield(spec,'sources') ...
            || isempty(spec.sources)
        value = X(:,1);
        return
    end
    operation = char(string(spec.source_operation));
    switch operation
        case 'aliases'
            value = X(:,1);
        case 'coalesce'
            value = nan(size(X,1),1);
            for j = 1:size(X,2)
                use = ~isfinite(value) ...
                    & isfinite(X(:,j));
                value(use) = X(use,j);
            end
        case 'mean'
            finite = isfinite(X);
            count = sum(finite,2);
            X(~finite) = 0;
            value = sum(X,2)./count;
            value(count == 0) = NaN;
        otherwise
            error('read_hydro_timeseries:BadSourceOperation', ...
                'Unsupported source operation %s for %s in %s.', ...
                operation,field,file);
    end
end

function [raw,filesUsed,plans] = local_read_named_files( ...
    dirData,id,schema,targetTime,plans)
%LOCAL_READ_NAMED_FILES Read and time-align named files for one basin.

    if nargin < 5 ...
            || isempty(plans)
        plans = struct();
    end
    raw = struct();
    filesUsed = cell(0,1);
    names = fieldnames(schema.files);
    variableNames = fieldnames(schema.variables);

    for i = 1:numel(names)
        name = names{i};
        selected = struct();
        for j = 1:numel(variableNames)
            field = variableNames{j};
            spec = schema.variables.(field);
            hasSource = (isfield(spec,'source') ...
                && ~isempty(spec.source)) ...
                || (isfield(spec,'sources') ...
                && ~isempty(spec.sources));
            if hasSource ...
                && isfield(spec,'file') ...
                    && strcmp(char(string(spec.file)),name)
                selected.(field) = spec;
            end
        end
        % A named file may be declared for provenance or future options
        % without being selected by this fixed schema.
        if isempty(fieldnames(selected))
            continue
        end

        fileSpec = schema.files.(name);
        file = local_file(dirData,fileSpec.pattern,id);
        if ~isfile(file)
            error('read_hydro_timeseries:MissingFile', ...
                'Missing named time-series file %s: %s', ...
            name,file);
        end

        subSchema = schema;
        subSchema.layout = 'one_file_per_basin';
        subSchema.variables = selected;
        subSchema.file = fileSpec;
        subSchema.format = char(string(fileSpec.format));
        if isfield(fileSpec,'time') ...
                && ~isempty(fieldnames(fileSpec.time))
            subSchema.time = local_overlay_struct(schema.time, ...
            fileSpec.time);
        end

        switch subSchema.format
            case 'csv'
                plan = [];
                if isfield(plans,name)
                    plan = plans.(name);
                end
                [part,plan] = local_read_csv( ...
                    file,subSchema,targetTime,plan);
                % A prepared plan is reusable only when the schema declares
                % the complete file columns. Otherwise headers may differ
                % legitimately between basins.
                if isfield(fileSpec,'variable_names') ...
                        && ~isempty(fileSpec.variable_names)
                    plans.(name) = plan;
                end
            case 'netcdf'
                part = local_read_netcdf(file,subSchema,targetTime);
            otherwise
                error('read_hydro_timeseries:UnsupportedFormat', ...
                    'Unsupported format for named file %s: %s.', ...
                    name,subSchema.format);
        end

        partFields = fieldnames(part);
        for j = 1:numel(partFields)
            field = partFields{j};
            if strcmp(field,'found')
                continue
            end
            if isfield(raw,field)
                error('read_hydro_timeseries:DuplicateVariable', ...
                    ['Variable %s was read from ' ...
                    'more than one named file.'],field);
            end
            raw.(field) = part.(field);
        end
        filesUsed{end+1,1} = file; %#ok<AGROW>
    end
end

function filesUsed = local_named_source_files(dirData,id,schema)
%LOCAL_NAMED_SOURCE_FILES Resolve selected named files without reading them.

    filesUsed = cell(0,1);
    names = fieldnames(schema.files);
    variableNames = fieldnames(schema.variables);
    for i = 1:numel(names)
        name = names{i};
        selected = false;
        for j = 1:numel(variableNames)
            spec = schema.variables.(variableNames{j});
            hasSource = (isfield(spec,'source') ...
                && ~isempty(spec.source)) ...
                || (isfield(spec,'sources') ...
                && ~isempty(spec.sources));
            if hasSource ...
                    && isfield(spec,'file') ...
                    && strcmp(char(string(spec.file)),name)
                selected = true;
                break
            end
        end
        if selected
            file = local_file( ...
                dirData,schema.files.(name).pattern,id);
            if ~isfile(file)
                error('read_hydro_timeseries:MissingFile', ...
                    'Missing named time-series file %s: %s', ...
                    name,file);
            end
            filesUsed{end+1,1} = file; %#ok<AGROW>
        end
    end
end

function out = local_overlay_struct(base,override)
    out = base;
    fields = fieldnames(override);
    for i = 1:numel(fields)
        out.(fields{i}) = override.(fields{i});
    end
end

function raw = local_read_netcdf(file,schema,targetTime)
    raw = struct();
    info = ncinfo(file);
    available = {info.Variables.Name};
    [sourceMap,~] = local_resolve_sources(schema,available,file);
    if strcmp(schema.time.mode,'index')
        stepSeconds = seconds(schema.timeline.step);
        start = round(seconds(targetTime(1)-schema.time.file_start) ...
            / stepSeconds) + 1;
        count = numel(targetTime);
        if start < 1
            error('read_hydro_timeseries:WindowBeforeFile', ...
                'Requested period starts before %s.',file);
        end
        fields = fieldnames(schema.variables);
        for i = 1:numel(fields)
            field = fields{i};
            spec = schema.variables.(field);
            if ~isfield(sourceMap,field)
                continue
            end
            columns = sourceMap.(field);
            X = nan(count,numel(columns));
            for j = 1:numel(columns)
                source = columns{j};
                try
                    x = double(ncread(file,source,start,count));
                catch ME
                    error('read_hydro_timeseries:NetCDFRead', ...
                        'Cannot read %s from %s: %s', ...
                    source,file,ME.message);
                end
                X(:,j) = local_apply_native_missing(x(:),spec);
            end
            raw.(field) = local_combine_sources(X,spec,field,file);
            if ~isempty(spec.quality)
                raw.(field)(~local_netcdf_quality_mask( ...
                    file,spec.quality,start,count)) = NaN;
            end
        end
    else
        sourceTime = local_netcdf_time(file,schema.time);
        [found,location] = ismember(targetTime,sourceTime);
        fields = fieldnames(schema.variables);
        
        for i = 1:numel(fields)
            field = fields{i};
            spec = schema.variables.(field);
            if ~isfield(sourceMap,field)
                continue
            end
            columns = sourceMap.(field);
            X = nan(numel(sourceTime),numel(columns));
            for j = 1:numel(columns)
                source = columns{j};
                try
                    x = double(ncread(file,source));
                catch ME
                    error('read_hydro_timeseries:NetCDFRead', ...
                        'Cannot read %s from %s: %s', ...
                    source,file,ME.message);
                end
                if numel(x) ~= numel(sourceTime)
                    error('read_hydro_timeseries:NetCDFLength', ...
                        ['%s and the time variable ' ...
                        'differ in length in %s.'], ...
                        source,file);
                end
                X(:,j) = local_apply_native_missing(x(:),spec);
            end
            values = nan(numel(targetTime),1);
            column = local_combine_sources(X,spec,field,file);
            if ~isempty(spec.quality)
                column(~local_netcdf_quality_mask( ...
                    file,spec.quality,1,numel(sourceTime))) = NaN;
            end
            values(found) = column(location(found));
            raw.(field) = values;
        end
    end
end

function valid = local_netcdf_quality_mask(file,quality,start,count)
    valid = true(count,1);
    sources = cellstr(string(quality.sources(:)));
    for j = 1:numel(sources)
        x = double(ncread(file,sources{j},start,count));
        x = x(:);
        if quality.require_finite
            valid = valid ...
                & isfinite(x);
        end
        valid = valid & x >= double(quality.valid_min) ...
            & x <= double(quality.valid_max);
    end
end

function raw = local_derive_variables( ...
    raw,schema,time,k,aux)

    fields = fieldnames(schema.variables);
    % Materialize constants first so later derivations may depend on them
    % regardless of structure field insertion order.
    for i = 1:numel(fields)
        field = fields{i};
        spec = schema.variables.(field);
        if isfield(spec,'derive') ...
            && strcmpi(string(spec.derive),'constant')
            if ~isfield(spec,'value') ...
                || ~isscalar(spec.value)
                error('read_hydro_timeseries:BadConstant', ...
                    'Constant derivation for %s requires scalar .value.',field);
            end
            raw.(field) = repmat(double(spec.value),numel(time),1);
        end
    end
    for i = 1:numel(fields)
        field = fields{i};
        spec = schema.variables.(field);
        if ~isfield(spec,'derive') ...
                || isempty(spec.derive)
            continue
        end
        % Layout-specific loaders may already materialize a declared
        % derived variable (for example a wide-file Tmin/Tmax mean).
        if isfield(raw,field)
            continue
        end
        method = lower(char(string(spec.derive)));
        switch method
            case 'constant'
                continue
            case 'mean_tmin_tmax'
                if ~isfield(raw,'Tmin') ...
                        || ~isfield(raw,'Tmax')
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        ['mean_Tmin_Tmax requires ' ...
                        'Tmin and Tmax.']);
                end
                raw.(field) = 0.5*(double(raw.Tmin)+double(raw.Tmax));
            case 'oudin'
                if ~isfield(raw,'T')
                    error( ...
                        'read_hydro_timeseries:MissingDerivedInput', ...
                        'Oudin PET requires temperature T.');
                end
                lat = aux(k,1);
                if ~isfinite(lat)
                    lat = 60;
                end
                if isfield(spec,'input_precision') ...
                        && strcmpi(string(spec.input_precision),'double')
                    Tair = double(raw.T);
                else
                    Tair = double(single(raw.T));
                end
                coefficient = 0.408/100;
                if isfield(spec,'coefficient')
                    coefficient = double(spec.coefficient);
                end
                raw.(field) = local_oudin_pet_generic( ...
                    time,Tair,lat,coefficient);
            case 'temperature_auto_celsius'
                input = 'TemperatureRaw';
                if isfield(spec,'input')
                    input = char(string(spec.input));
                end
                if ~isfield(raw,input)
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        'Automatic temperature conversion requires %s.',input);
                end
                value = double(raw.(input));
                finite = value(isfinite(value));
                if ~isempty(finite) ...
                    && median(finite) > 100
                    value = value-273.15;
                end
                raw.(field) = value;
            case 'hargreaves'
                if ~all(isfield(raw,{'T','Tmin','Tmax'}))
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        ['Hargreaves PET requires ' ...
                        'T, Tmin, and Tmax.']);
                end
                lat = aux(k,1);
                if ~isfinite(lat) ...
                        && isfield(spec,'latitude_fallback')
                    lat = double(spec.latitude_fallback);
                end
                raw.(field) = local_hargreaves_pet_generic( ...
                    time,double(raw.T),double(raw.Tmin), ...
                    double(raw.Tmax),lat);
                if isfield(spec,'fill_missing') ...
                        && ~isempty(spec.fill_missing)
                    raw.(field)(~isfinite(raw.(field))) = ...
                        double(spec.fill_missing);
                end
            case 'vapor_pressure_dewpoint'
                input = 'Tdew';
                if isfield(spec,'input')
                    input = char(string(spec.input));
                end
                if ~isfield(raw,input)
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        'Dew-point vapor pressure requires %s.',input);
                end
                Td = double(raw.(input));
                raw.(field) = 1000.*0.6108.*exp((17.27.*Td)./(Td+237.3));
            case 'vapor_pressure_relative_humidity'
                tInput = 'T';
                rhInput = 'RelativeHumidity';
                if isfield(spec,'temperature_input')
                    tInput = char(string(spec.temperature_input));
                end
                if isfield(spec,'humidity_input')
                    rhInput = char(string(spec.humidity_input));
                end
                if ~all(isfield(raw,{tInput,rhInput}))
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        ['Relative-humidity vapor pressure requires ' ...
                         '%s and %s.'],tInput,rhInput);
                end
                Ta = double(raw.(tInput));
                rh = min(max(double(raw.(rhInput)),0),100);
                raw.(field) = 1000.*(rh./100).*0.6108.* ...
                    exp((17.27.*Ta)./(Ta+237.3));
            case 'mixing_ratio_to_specific_humidity'
                input = 'MixingRatio';
                if isfield(spec,'input')
                    input = char(string(spec.input));
                end
                if ~isfield(raw,input)
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        'Mixing-ratio conversion requires %s.',input);
                end
                inputScale = 1;
                if isfield(spec,'input_scale')
                    inputScale = double(spec.input_scale);
                end
                ratio = double(raw.(input)).*inputScale;
                raw.(field) = ratio./(1+ratio);
            case 'wind_uv_to_2m'
                if ~all(isfield(raw,{'WindU10','WindV10'}))
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        'wind_uv_to_2m requires WindU10 and WindV10.');
                end
                raw.(field) = hypot(double(raw.WindU10), ...
                    double(raw.WindV10)).*(4.87/log(67.8*10-5.42));
            case 'fao56_daily'
                required = {'Tmin','Tmax','Radiation','VaporPressure','Wind2'};
                if ~all(isfield(raw,required))
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        ['fao56_daily requires Tmin, Tmax, Radiation, ' ...
                         'VaporPressure, and Wind2.']);
                end
                local_require_pet_method(spec,field);
                lat = aux(k,1);
                elev = aux(k,2);
                if ~isfinite(lat) ...
                    || ~isfinite(elev)
                    error('read_hydro_timeseries:MissingPetAux', ...
                        'fao56_daily requires finite latitude and elevation.');
                end
                radiationScale = 1;
                if isfield(spec,'radiation_scale')
                    radiationScale = double(spec.radiation_scale);
                end
                raw.(field) = et0_fao56_daily(double(raw.Tmax), ...
                    double(raw.Tmin),double(raw.Radiation).*radiationScale, ...
                    double(raw.VaporPressure),day(time,'dayofyear'), ...
                    lat,elev,double(raw.Wind2),double(spec.method));
            case 'lamah_hourly'
                required = {'T','Tdew','WindU10','WindV10', ...
                    'NetSolar','NetThermal','Pressure'};
                if ~all(isfield(raw,required))
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        ['lamah_hourly requires T, Tdew, WindU10, WindV10, ' ...
                         'NetSolar, NetThermal, and Pressure.']);
                end
                local_require_pet_method(spec,field);
                raw.(field) = et0_lamah_hourly(double(raw.T), ...
                    double(raw.Tdew),double(raw.WindU10), ...
                    double(raw.WindV10),double(raw.NetSolar), ...
                    double(raw.NetThermal),double(raw.Pressure), ...
                    double(spec.method));
            case 'fao56_hourly_estimated_longwave'
                required = {'T','Radiation','Pressure','SpecificHumidity', ...
                    'WindU10','WindV10','VaporPressure','CloudCover'};
                if ~all(isfield(raw,required))
                    error('read_hydro_timeseries:MissingDerivedInput', ...
                        ['fao56_hourly_estimated_longwave requires T, ' ...
                         'Radiation, Pressure, SpecificHumidity, WindU10, ' ...
                         'WindV10, VaporPressure, and CloudCover.']);
                end
                local_require_pet_method(spec,field);
                Ta = double(raw.T);
                Tk = Ta+273.15;
                ea = double(raw.VaporPressure)./1000;
                epsClear = 1.24.*max(ea./Tk,0).^(1/7);
                epsClear = min(max(epsClear,0.55),1.00);
                cloud = min(max(double(raw.CloudCover)./100,0),1);
                epsSky = min(max(epsClear.*(1+0.22.*cloud.^2),0.55),1.00);
                Rld = epsSky.*5.670374419e-8.*Tk.^4;
                pressureScale = 1;
                radiationScale = 1;
                if isfield(spec,'pressure_scale')
                    pressureScale = double(spec.pressure_scale);
                end
                if isfield(spec,'radiation_scale')
                    radiationScale = double(spec.radiation_scale);
                end
                raw.(field) = et0_fao56_hourly(Rld, ...
                    double(raw.Pressure).*pressureScale, ...
                    double(raw.Radiation).*radiationScale, ...
                    double(raw.SpecificHumidity), ...
                    Ta,double(raw.WindU10),double(raw.WindV10), ...
                    double(spec.method));
            otherwise
                error( ...
                    'read_hydro_timeseries:UnknownDerivation', ...
                    'Unknown derivation "%s".',method);
        end
    end
end

function local_require_pet_method(spec,field)
    if ~isfield(spec,'method') ...
        || ~isscalar(spec.method) ...
            || ~ismember(double(spec.method),[1 2 3])
        error('read_hydro_timeseries:BadPetMethod', ...
            '%s derivation requires method 1, 2, or 3.',field);
    end
end

function Ep = local_hargreaves_pet_generic(t,Tmean,Tmin,Tmax,lat)
    if ~isfinite(lat)
        error('read_hydro_timeseries:MissingLatitude', ...
            ['Hargreaves PET requires basin ' ...
            'latitude or latitude_fallback.']);
    end
    lat = max(min(double(lat),89.9),-89.9);
    phi = lat*pi/180;
    J = day(t,'dayofyear');
    Gsc = 0.0820;
    dr = 1 + 0.033*cos(2*pi*J/365);
    delta = 0.409*sin(2*pi*J/365-1.39);
    ws = acos(max(min(-tan(phi).*tan(delta),1),-1));
    Ra = (24*60/pi)*Gsc.*dr.*(ws.*sin(phi).*sin(delta) ...
        + cos(phi).*cos(delta).*sin(ws));
    dT = max(double(Tmax(:))-double(Tmin(:)),0);
    Ep = 0.0023.*(Ra(:)/2.45).*(double(Tmean(:))+17.8).*sqrt(dT);
    Ep(Ep < 0) = 0;
end

function Ep = local_oudin_pet_generic( ...
    t,Tair,latDeg,coefficient)

    doy = day(t,'dayofyear');
    phi = deg2rad(latDeg);
    Gsc = 0.0820;    
    dr = 1 + 0.033*cos(2*pi*doy/365);
    delta = ...
        0.409*sin(2*pi*doy/365 - 1.39);
    ws = acos(max(-1,min(1, ...
        -tan(phi).*tan(delta))));
    Ra = (24*60/pi) * Gsc .* dr .* ...
        (ws.*sin(phi).*sin(delta) + ...
        cos(phi).*cos(delta).*sin(ws));
    Ep = double(coefficient).*Ra(:).*max(Tair(:)+5,0);
    Ep(~isfinite(Ep)) = NaN;
    Ep = single(Ep);
    Ep(~isfinite(Ep) | Ep < 0) = 0;
end

function value = local_variable(raw,spec,field,file,k,aux,step)
    value = double(raw.(field)(:));
    missing = double(spec.missing_values(:));
    for i = 1:numel(missing)
        value(value == missing(i)) = NaN;
    end
    value(value < spec.valid_min | value > spec.valid_max) = NaN;
    value = value*double(spec.scale) + double(spec.offset);
    % Applying max/min with the default +/-Inf bounds is not a no-op in
    % MATLAB: max(NaN,-Inf) and min(NaN,Inf) replace the missing value.
    % Clip only finite observations and only when a real bound is declared.
    finiteValue = isfinite(value);
    if isfinite(spec.clip_min)
        value(finiteValue & value < spec.clip_min) = spec.clip_min;
    end
    if isfinite(spec.clip_max)
        value(finiteValue & value > spec.clip_max) = spec.clip_max;
    end
    if isfield(spec,'fill_missing') ...
        && ~isempty(spec.fill_missing)
        value(~isfinite(value)) = double(spec.fill_missing);
    end

    if spec.area_normalize
        area = aux(k,3);
        if ~isfinite(area) ...
            || area <= 0
            error('read_hydro_timeseries:MissingArea', ...
                'A positive basin area is required for %s.',file);
        end
        value = value*seconds(step)*1000/area;
    end
    if ~isempty(spec.transform)
        if ~isa(spec.transform,'function_handle')
            error('read_hydro_timeseries:BadTransform', ...
                'Transform for %s must be a function handle.',field);
        end
        context = struct('field',field,'file',file, ...
            'basin_index',k,'aux',aux(k,:),'step',step);
        value = spec.transform(value,context);
        value = double(value(:));
    end
end

function time = local_datetime(value,spec)
    if isdatetime(value)
        time = value(:);
    elseif isempty(spec.input_format)
        time = datetime(string(value(:)));
    else
        time = datetime(string(value(:)), ...
            'InputFormat',char(spec.input_format));
    end
    if ~isempty(spec.timezone)
        time.TimeZone = char(spec.timezone);
    end
end

function time = local_netcdf_time(file,spec)
    value = double(ncread(file,char(spec.column)));
    units = lower(string(spec.units));
    origin = spec.origin;
    % CF-compliant NetCDF files may use a different reference epoch for
    % every basin (for example, "days since 1997-01-01 12:00:00").  The
    % file metadata is authoritative; the schema values are retained as a
    % fallback for legacy files without a parseable units attribute.
    try
        fileUnits = string(ncreadatt(file,char(spec.column),'units'));
        token = regexp(char(fileUnits), ...
            '^\s*(seconds?|minutes?|hours?|days?)\s+since\s+(.+?)\s*$', ...
            'tokens','once','ignorecase');
        if ~isempty(token)
            units = lower(string(token{1}));
            origin = local_netcdf_time_origin(token{2});
        end
    catch
        % Use the validated schema definition when metadata is absent.
    end
    if units == "datetime"
        error('read_hydro_timeseries:NumericDatetime', ...
            'A numeric NetCDF time variable requires explicit units.');
    elseif startsWith(units,"hour")
        time = origin + hours(value(:));
    elseif startsWith(units,"day")
        time = origin + days(value(:));
    elseif startsWith(units,"minute")
        time = origin + minutes(value(:));
    elseif startsWith(units,"second")
        time = origin + seconds(value(:));
    else
        error('read_hydro_timeseries:BadTimeUnits', ...
            'Unsupported NetCDF time units: %s.',units);
    end
    if isfield(spec,'snap') ...
        && ~isempty(spec.snap)
        time = dateshift(time,'start',char(string(spec.snap)));
    end
end

function file = local_file(dirData,pattern,id)

    idText = char(string(id));
    idNumber = str2double(idText);
    if isfinite(idNumber) && idNumber == fix(idNumber)
        idText02 = sprintf('%02d',idNumber);
    else
        idText02 = idText;
    end
    relative = strrep(char(pattern),'{gauge02}',idText02);
    relative = strrep(relative,'{gauge}',idText);

    candidate = fullfile( ...
        char(string(dirData)),relative);

    if contains(candidate,'*') ...
            || contains(candidate,'?')

        D = dir(candidate);

        if isempty(D)
            error('read_hydro_timeseries:MissingFile', ...
                'No file matches pattern: %s',candidate);
        end

        if numel(D) > 1
            error('read_hydro_timeseries:AmbiguousFile', ...
                ['More than one file matches pattern %s. ' ...
                 'Schema must identify one source file.'], ...
                candidate);
        end

        file = fullfile(D(1).folder,D(1).name);

    else
        file = candidate;
    end
end

function ids = local_ids(value,spec)
    ids = strip(string(value(:)));
    if spec.strip_decimal
        ids = regexprep(ids,'\.0+$','');
    end
    if spec.pad_width > 0
        ids = arrayfun(@(x) pad(x,double(spec.pad_width),'left','0'),ids);
    end
    if strlength(string(spec.prefix)) > 0
        ids = string(spec.prefix) + ids;
    end
end

function idx = local_score_indices(split,n)
    if isfield(split,'idx') ...
        && ~isempty(split.idx)
        idx = round(double(split.idx(:)));
        idx = idx(isfinite(idx) & idx >= 1 & idx <= n);
    else
        idx = (1:n).';
    end
end

function aux = local_aux(schema,bas,ids,dirData)
    K = numel(ids);
    if ~isempty(schema.aux.provider)
        if ~isa(schema.aux.provider,'function_handle')
            error('read_hydro_timeseries:BadAuxProvider', ...
                'schema.aux.provider must be a function handle.');
        end
        aux = double(schema.aux.provider(ids,bas));
    elseif ~isempty(schema.aux.values)
        aux = double(schema.aux.values);
    elseif isfield(schema.aux,'tables') ...
        && ~isempty(schema.aux.tables)
        aux = local_aux_tables(schema.aux.tables,dirData,ids);
    else
        aux = nan(K,3);
        aux(:,1) = local_bas_vector(bas,{'lat','latitude'},K,1);
        aux(:,2) = local_bas_vector( ...
            bas,{'elev','elevation','elev_mean'},K,1);
        area = local_bas_vector(bas,{'area','area_m2'},K,1);
        area(isfinite(area) & area < 1e7) = ...
            area(isfinite(area) & area < 1e7)*1e6;
        aux(:,3) = area;
    end
    if isfield(schema.aux,'series') ...
        && ~isempty(schema.aux.series)
        aux = local_aux_series(aux,schema.aux.series,dirData,ids);
    end
    if isfield(schema.aux,'file_header') ...
        && ...
            ~isempty(schema.aux.file_header)
        aux = local_aux_file_header( ...
            aux,schema.aux.file_header,schema,dirData,ids);
    end
    if ~isequal(size(aux),[K 3])
        error('read_hydro_timeseries:BadAux', ...
            'Auxiliary metadata must be a K-by-3 numeric matrix.');
    end
end

function aux = local_aux_file_header(aux,spec,schema,dirData,ids)
    if ~strcmp(schema.layout,'one_file_per_basin')
        error('read_hydro_timeseries:BadAuxHeader', ...
            'aux.file_header requires one_file_per_basin layout.');
    end
    if ~isfield(spec,'lines') ...
            || numel(spec.lines) ~= 3
        error('read_hydro_timeseries:BadAuxHeader', ...
            'aux.file_header.lines must contain [lat elev area] line numbers.');
    end
    scale = ones(1,3);
    if isfield(spec,'scale')
        scale = double(spec.scale);
    end
    for k = 1:numel(ids)
        file = local_file(dirData,schema.file.pattern,ids(k));
        fid = fopen(file,'r');
        if fid < 0
            error('read_hydro_timeseries:OpenFile', ...
                'Cannot open auxiliary header file: %s',file);
        end
        cleanup = onCleanup(@() fclose(fid));
        lines = cell(max(double(spec.lines)),1);
        for j = 1:numel(lines)
            lines{j} = fgetl(fid);
        end
        for j = 1:3
            line = lines{double(spec.lines(j))};
            token = regexp(char(line), ...
                '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?','match','once');
            if ~isempty(token)
                aux(k,j) = str2double(token)*scale(j);
            end
        end
    end
end

function aux = local_aux_series(aux,series,dirData,ids)
    if ~isstruct(series) ...
        || ~isscalar(series)
        error('read_hydro_timeseries:BadAuxSeries', ...
            'schema.aux.series must be a scalar structure.');
    end
    fields = {'lat','elev','area'};
    for j = 1:3
        field = fields{j};
        if ~isfield(series,field) ...
                || isempty(series.(field))
            continue
        end
        spec = series.(field);
        for required = {'file_pattern','source','reducer'}
            if ~isfield(spec,required{1}) ...
                || isempty(spec.(required{1}))
                error('read_hydro_timeseries:BadAuxSeries', ...
                    'aux.series.%s requires .%s.',field,required{1});
            end
        end
        for k = 1:numel(ids)
            if isfinite(aux(k,j))
                continue
            end
            file = local_file(dirData,spec.file_pattern,ids(k));
            opts = detectImportOptions(file,'VariableNamingRule','preserve');
            opts.SelectedVariableNames = {char(string(spec.source))};
            T = readtable(file,opts);
            x = local_numeric_column(T.(char(string(spec.source))), ...
                char(string(spec.source)),file);
            switch lower(char(string(spec.reducer)))
                case 'median'
                    value = median(x(isfinite(x)),'omitnan');
                otherwise
                    error('read_hydro_timeseries:BadAuxReducer', ...
                        'Unsupported aux reducer %s.',spec.reducer);
            end
            if isfield(spec,'derive') ...
                && ~isempty(spec.derive)
                switch lower(char(string(spec.derive)))
                    case 'pressure_elevation'
                        if value > 20000
                            p = value/1000;
                        elseif value > 200
                            p = value/10;
                        else
                            p = value;
                        end
                        if ~isfinite(p) ...
                                || p <= 0 ...
                                || p > 110
                            value = NaN;
                        else
                            value = 44330*(1-(p/101.325)^0.1903);
                        end
                    otherwise
                        error('read_hydro_timeseries:BadAuxDerivation', ...
                            'Unsupported aux derivation %s.',spec.derive);
                end
            end
            scale = 1;
            if isfield(spec,'scale')
                scale = spec.scale;
            end
            aux(k,j) = double(value)*double(scale);
        end
    end
end

function aux = local_aux_tables(tables,dirData,ids)
    if ~isstruct(tables) ...
        || ~isscalar(tables)
        error('read_hydro_timeseries:BadAuxTables', ...
            'schema.aux.tables must be a scalar structure.');
    end
    aux = nan(numel(ids),3);
    names = fieldnames(tables);
    for i = 1:numel(names)
        spec = tables.(names{i});
        layout = 'rows';
        if isfield(spec,'layout') ...
            && ~isempty(spec.layout)
            layout = char(lower(string(spec.layout)));
        end
        if strcmp(layout,'wide_attributes')
            aux = local_aux_wide_attributes(aux,spec,dirData,ids,names{i});
            continue
        elseif ~strcmp(layout,'rows')
            error('read_hydro_timeseries:BadAuxTableLayout', ...
                'Unknown aux table layout %s for %s.',layout,names{i});
        end
        required = {'file','key'};
        for j = 1:numel(required)
            if ~isfield(spec,required{j}) ...
                || isempty(spec.(required{j}))
                error('read_hydro_timeseries:BadAuxTable', ...
                    'Aux table %s requires .%s.',names{i},required{j});
            end
        end
        file = fullfile(char(string(dirData)),char(string(spec.file)));
        opts = detectImportOptions(file, ...
            'VariableNamingRule','preserve');
        key = char(string(spec.key));
        opts = setvartype(opts,key,'string');
        T = readtable(file,opts);
        tableIds = local_metadata_ids(T.(char(string(spec.key))),spec);
        requestIds = upper(ids);
        if isfield(spec,'request_prefix_length') ...
                && ~isempty(spec.request_prefix_length)
            n = double(spec.request_prefix_length);
            requestIds = extractBefore(requestIds, ...
                min(strlength(requestIds)+1,n+1));
        end
        [found,loc] = ismember(requestIds,tableIds);
        fields = {'lat','elev','area'};
        for j = 1:3
            field = fields{j};
            if j == 1 ...
                    && isfield(spec,'lat_derive') ...
                    && strcmpi(string(spec.lat_derive),'lambert93_to_wgs84')
                required = {'lat_x','lat_y'};
                if ~all(isfield(spec,required))
                    error('read_hydro_timeseries:BadAuxDerivation', ...
                        'Lambert-93 latitude requires lat_x and lat_y.');
                end
                xCoord = local_numeric_column(T.(char(string(spec.lat_x))), ...
                    char(string(spec.lat_x)),file);
                yCoord = local_numeric_column(T.(char(string(spec.lat_y))), ...
                    char(string(spec.lat_y)),file);
                x = local_lambert93_latitude(xCoord,yCoord);
            else
                if ~isfield(spec,field) ...
                        || isempty(spec.(field))
                    continue
                end
                x = local_numeric_column(T.(char(string(spec.(field)))), ...
                    char(string(spec.(field))),file);
            end
            scaleField = [field '_scale'];
            scale = 1;
            if isfield(spec,scaleField)
                scale = spec.(scaleField);
            end
            use = found & ~isfinite(aux(:,j));
            aux(use,j) = x(loc(use))*double(scale);
        end
    end
end

function lat = local_lambert93_latitude(X,Y)
    X = double(X(:));
    Y = double(Y(:));
    lat = nan(size(X));
    good = isfinite(X) & isfinite(Y);
    degLike = good ...
        & abs(X) <= 20 ...
        & abs(Y) <= 60;
    lat(degLike) = Y(degLike);
    good = good & ~degLike;
    if ~any(good)
        return
    end
    e = sqrt(0.00669438002290);
    n = 0.7256077650532670;
    C = 11754255.426096;
    xs = 700000.0;
    ys = 12655612.049876;
    x = X(good);
    y = Y(good);
    R = sqrt((x-xs).^2+(y-ys).^2);
    L = -(1/n).*log(R./C);
    phi = 2*atan(exp(L))-pi/2;
    for it = 1:8
        phi = 2*atan(((1+e*sin(phi))./(1-e*sin(phi))).^(e/2).*exp(L))-pi/2;
    end
    lat(good) = phi*180/pi;
end

function aux = local_aux_wide_attributes(aux,spec,dirData,ids,name)
    required = {'file','attribute_column'};
    for j = 1:numel(required)
        if ~isfield(spec,required{j}) ...
            || isempty(spec.(required{j}))
            error('read_hydro_timeseries:BadAuxTable', ...
                'Wide aux table %s requires .%s.',name,required{j});
        end
    end
    file = fullfile(char(string(dirData)),char(string(spec.file)));
    T = readtable(file,'VariableNamingRule','preserve');
    attributeColumn = char(string(spec.attribute_column));
    columnNames = string(T.Properties.VariableNames);
    dataColumns = columnNames(columnNames ~= string(attributeColumn));
    if isfield(spec,'column_prefix') ...
        && ~isempty(spec.column_prefix)
        dataColumns = dataColumns(startsWith(dataColumns, ...
            string(spec.column_prefix),'IgnoreCase',true));
    end
    tableIds = local_metadata_ids(dataColumns,spec);
    [found,loc] = ismember(upper(ids),tableIds);
    fields = {'lat','elev','area'};
    attributes = string(T.(attributeColumn));
    for j = 1:3
        field = fields{j};
        rowField = [field '_row'];
        if ~isfield(spec,rowField) ...
                || isempty(spec.(rowField))
            continue
        end
        row = find(strcmpi(attributes,string(spec.(rowField))),1);
        if isempty(row)
            error('read_hydro_timeseries:MissingAuxAttribute', ...
                'Attribute %s is absent from %s.',spec.(rowField),file);
        end
        x = nan(numel(dataColumns),1);
        for c = 1:numel(dataColumns)
            x(c) = local_numeric_column(T{row,char(dataColumns(c))}, ...
                char(dataColumns(c)),file);
        end
        scaleField = [field '_scale'];
        scale = 1;
        if isfield(spec,scaleField)
            scale = spec.(scaleField);
        end
        use = found & ~isfinite(aux(:,j));
        aux(use,j) = x(loc(use))*double(scale);
    end
end

function ids = local_metadata_ids(value,spec)
    ids = upper(strip(string(value(:))));
    if ~isfield(spec,'strip_decimal') ...
        || spec.strip_decimal
        ids = regexprep(ids,'\.0+$','');
    end
    if isfield(spec,'strip_prefix') ...
        && ~isempty(spec.strip_prefix)
        ids = regexprep(ids,['^' regexptranslate('escape', ...
            char(string(spec.strip_prefix)))],'','ignorecase');
    end
    if isfield(spec,'pad_width') ...
        && ~isempty(spec.pad_width) ...
            && double(spec.pad_width) > 0
        ids = arrayfun(@(x) pad(x,double(spec.pad_width), ...
            'left','0'),ids);
    end
end

function value = local_bas_vector(bas,names,K,default) %#ok<INUSD>
    value = nan(K,1);
    for i = 1:numel(names)
        if isfield(bas,names{i}) ...
            && numel(bas.(names{i})) == K
            value = double(bas.(names{i})(:));
            return
        end
    end

    % SAMPLE_BASINS stores selected auxiliary metadata in bas.zone. Use
    % those vectors when the same fields are not duplicated at basin level.
    if isfield(bas,'zone') ...
            && isstruct(bas.zone) ...
            && isscalar(bas.zone)
        for i = 1:numel(names)
            if isfield(bas.zone,names{i}) ...
                    && numel(bas.zone.(names{i})) == K
                value = double(bas.zone.(names{i})(:));
                return
            end
        end
    end
end

function origin = local_netcdf_time_origin(value)
    text = strtrim(char(string(value)));
    text = regexprep(text,'T',' ');
    text = regexprep(text,'\s*(UTC|Z)\s*$','','ignorecase');
    formats = {'yyyy-MM-dd HH:mm:ss.SSSSSS', ...
        'yyyy-MM-dd HH:mm:ss.SSS','yyyy-MM-dd HH:mm:ss', ...
        'yyyy-MM-dd HH:mm','yyyy-MM-dd'};
    for i = 1:numel(formats)
        try
            origin = datetime(text,'InputFormat',formats{i});
            return
        catch
        end
    end
    error('read_hydro_timeseries:BadNetCDFTimeOrigin', ...
        'Cannot parse NetCDF time origin: %s.',value);
end

function x = local_numeric_column(x,source,file)

    if isnumeric(x) ...
            || islogical(x)
        x = double(x);
        return
    end
    
    if iscell(x) ...
            || isstring(x) ...
            || ischar(x) ...
            || iscategorical(x)
        s = string(x);
        s = strip(s);
    
        % Treat common textual missing-value representations as missing.
        missingMask = ismissing(s) ...
            | s == "" ...
            | strcmpi(s,"NA") ...
            | strcmpi(s,"NaN") ...
            | strcmpi(s,"null") ...
            | strcmpi(s,"missing");
    
        x = str2double(s);
        x(missingMask) = NaN;
    
        % Protect against silently converting genuinely nonnumeric data.
        bad = ~missingMask ...
            & isnan(x);
        if any(bad)
            examples = unique(s(bad));
            examples = examples(1:min(5,numel(examples)));
    
            error('read_hydro_timeseries:NonNumericColumn', ...
                ['Column "%s" in %s contains ' ...
                'nonnumeric values: %s'], ...
                source,file,strjoin(cellstr(examples),', '));
        end
    
        x = double(x);
        return
    end
    
    error('read_hydro_timeseries:UnsupportedColumnType', ...
        'Unsupported MATLAB type "%s" for column "%s" in %s.', ...
        class(x),source,file);
end

function [rawAll,filesUsed] = local_read_wide_files( ...
    dirData,ids,targetTime,options,schema)

    K = numel(ids);
    rawMatrix = struct();
    filesUsed = {};
    
    fields = fieldnames(schema.variables);
    
    % --------------------------------------------------
    % First read variables that come directly from files
    % --------------------------------------------------
    for i = 1:numel(fields)
    
        field = fields{i};
        spec = schema.variables.(field);
    
        if ~isfield(spec,'files') ...
            || isempty(spec.files)
            continue
        end
    
        file = local_selected_wide_file( ...
            dirData,spec,options,field);
    
        X = local_read_wide_csv( ...
            file,ids,targetTime,schema,spec);
    
        rawMatrix.(field) = X;
        filesUsed{end+1} = file; %#ok<AGROW>
    end
    
    filesUsed = unique(filesUsed,'stable');
    
    % --------------------------------
    % Then construct derived variables
    % --------------------------------
    for i = 1:numel(fields)
    
        field = fields{i};
        spec = schema.variables.(field);
    
        if ~isfield(spec,'derive') ...
                || isempty(spec.derive)
            continue
        end
    
        method = lower(char(string(spec.derive)));
    
        switch method
            case 'mean_tmin_tmax'
    
                if ~isfield(rawMatrix,'Tmin') ...
                        || ~isfield(rawMatrix,'Tmax')
                    error(['read_hydro_timeseries:' ...
                        'MissingDerivedInput'], ...
                        ['Derived temperature requires ' ...
                        'Tmin and Tmax.']);
                end
    
                rawMatrix.(field) = 0.5 * ...
                    (double(rawMatrix.Tmin) ...
                    + double(rawMatrix.Tmax));
    
            otherwise
                error('read_hydro_timeseries:UnknownDerivation', ...
                    'Unknown derivation "%s" for variable %s.', ...
                    method,field);
        end
    end
    
    % ---------------------------------------------------
    % Convert nTime-by-K matrices to one raw struct/basin
    % ---------------------------------------------------
    rawAll = cell(K,1);
    
    for k = 1:K
        R = struct();
    
        names = fieldnames(rawMatrix);
        for i = 1:numel(names)
            name = names{i};
            R.(name) = rawMatrix.(name)(:,k);
        end
    
        rawAll{k} = R;
    end
end

function file = local_selected_wide_file( ...
    dirData,spec,options,field)

    files = cellstr(string(spec.files));    
    choice = 1;
    
    if isfield(spec,'default') ...
            && ~isempty(spec.default)
        choice = round(double(spec.default));
    end
    
    if isfield(spec,'selector') ...
            && ~isempty(spec.selector)
    
        selector = char(string(spec.selector));   
        if isfield(options,selector) ...
                && ~isempty(options.(selector))
            choice = round(double(options.(selector)));
        end
    end
    
    if choice < 1 ...
            || choice > numel(files)
        error('read_hydro_timeseries:BadWideSelection', ...
            ['Selection %d for variable %s is outside ' ...
            'the valid range 1:%d.'], ...
            choice,field,numel(files));
    end
    
    file = fullfile( ...
        char(string(dirData)), ...
        files{choice});
    
    if ~isfile(file)
        error('read_hydro_timeseries:MissingFile', ...
            'Missing wide-file source: %s',file);
    end
end

function X = local_read_wide_csv( ...
        file,ids,targetTime,schema,spec)

    % --------------------------------------------------------
    % Read the native header directly. This is substantially
    % faster and more predictable than detectImportOptions for
    % very wide basin-by-column files.
    % --------------------------------------------------------
    fid = fopen(file,'r');
    if fid < 0
        error('read_hydro_timeseries:OpenFile', ...
            'Cannot open wide file: %s',file);
    end

    cleaner = onCleanup(@() fclose(fid));
    hdr = fgetl(fid);

    if ~ischar(hdr)
        error('read_hydro_timeseries:EmptyFile', ...
            'Wide file has no header: %s',file);
    end

    delimiter = schema.file.delimiter;
    rawNames = string(strsplit(hdr,delimiter));
    rawNames = strtrim(erase(rawNames,char(65279)));

    % Same normalization used by read_meteo_AU.
    names = matlab.lang.makeValidName( ...
        lower(regexprep(rawNames,'[^a-zA-Z0-9]','')), ...
        'ReplacementStyle','delete');

    if strcmp(schema.time.mode,'column')
        it = find(names == matlab.lang.makeValidName( ...
            lower(regexprep(string(schema.time.column), ...
            '[^a-zA-Z0-9]','')),'ReplacementStyle','delete'),1);
        if isempty(it)
            error('read_hydro_timeseries:MissingDateColumn', ...
                'Wide file %s lacks date column %s.', ...
                file,string(schema.time.column));
        end
        dateCols = it;
    else
        iy = find(ismember(names,{'year','yyyy','yr'}),1);
        im = find(ismember(names,{'month','mon','mm'}),1);
        id = find(ismember(names,{'day','dd'}),1);
        if isempty(iy) ...
                || isempty(im) ...
                || isempty(id)
            error('read_hydro_timeseries:MissingDateColumns', ...
                ['Wide file %s does not contain recognizable ' ...
                 'year/month/day columns.'],file);
        end
        dateCols = [iy im id];
    end

    % --------------------------------------------------------
    % Basin columns are every non-date column. Match using the
    % original header text because that contains station IDs.
    % --------------------------------------------------------
    dataCols = setdiff(1:numel(rawNames), ...
        dateCols,'stable');
    fileIDs = local_clean_wide_ids(rawNames(dataCols),schema.id);
    reqIDs = local_clean_wide_ids(ids,schema.id);
    [tf,jloc] = ismember(reqIDs,fileIDs);

    if ~all(tf)
        missing = ids(~tf);

        error('read_hydro_timeseries:MissingBasinColumns', ...
            ['Wide file %s is missing %d requested basin ' ...
             'column(s), e.g. %s.'], ...
            file,nnz(~tf),char(missing(1)));
    end

    [cacheHit,X] = local_wide_cache( ...
        'load',file,fileIDs,reqIDs,targetTime,schema,spec,[]);
    if cacheHit
        return
    end

    % Translate location within basin columns back to actual
    % CSV/table column indices.
    basinCols = dataCols(jloc);

    cacheEnabled = schema.cache.enabled ...
        && schema.cache.canonical_enabled;
    readBasinCols = basinCols;
    if cacheEnabled
        readBasinCols = dataCols;
    end

    selected = unique( ...
        [dateCols readBasinCols(:).'], ...
        'stable');

    % ---------------------------------------------------------
    % Define import behavior explicitly. CAMELS-AU wide files
    % are numeric throughout.
    % ---------------------------------------------------------
    opts = delimitedTextImportOptions( ...
        'Delimiter',schema.file.delimiter, ...
        'VariableNamingRule','preserve');

    opts.DataLines = [2 Inf];
    opts.VariableNames = cellstr(names);
    opts.VariableTypes = repmat( ...
        {'double'},1,numel(names));
    if strcmp(schema.time.mode,'column')
        opts.VariableTypes{it} = 'string';
    end
    opts.SelectedVariableNames = ...
        cellstr(names(selected));
    T = readtable(file,opts);

    % ----------------------
    % Construct native dates
    % ----------------------
    if strcmp(schema.time.mode,'column')
        parsedTime = datetime(T.(char(names(it))), ...
            'InputFormat',schema.time.input_format);
        if isfield(spec,'time_mode') ...
            && ...
                strcmpi(string(spec.time_mode), ...
                'row_sequence_from_first')
            sourceTime = parsedTime(1) + ...
                (0:height(T)-1)'.*schema.timeline.step;
        else
            sourceTime = parsedTime;
        end
    else
        yr = double(T.(char(names(iy))));
        mo = double(T.(char(names(im))));
        dy = double(T.(char(names(id))));
        sourceTime = datetime(yr,mo,dy);
    end

    if schema.strict_time ...
            && any(diff(sourceTime) <= seconds(0))

        error('read_hydro_timeseries:NonIncreasingTime', ...
            'Dates are not strictly increasing in %s.',file);
    end

    % Use single here to reproduce read_meteo_AU exactly.
    sourceValues = nan(height(T),numel(readBasinCols),'single');
    for k = 1:numel(readBasinCols)
        col = readBasinCols(k);
        name = char(names(col));
        column = double(T.(name));
        column = local_apply_native_missing( ...
            column,spec);
        sourceValues(:,k) = single(column);
    end

    if cacheEnabled
        local_wide_cache('save',file,fileIDs,reqIDs,targetTime, ...
            schema,spec,struct('time',sourceTime, ...
            'values',sourceValues));
        requestedValues = sourceValues(:,jloc);
    else
        requestedValues = sourceValues;
    end

    [found,row] = ismember(targetTime,sourceTime);
    X = nan(numel(targetTime),numel(ids),'single');
    X(found,:) = requestedValues(row(found),:);
end

function [hit,X] = local_wide_cache( ...
    action,file,fileIDs,requestIDs,targetTime,schema,spec,data)

    hit = false;
    X = [];
    eligible = schema.cache.enabled ...
        && schema.cache.canonical_enabled;
    if ~eligible
        return
    end
    cacheDir = fullfile(fileparts(file), ...
        char(string(schema.cache.directory)));
    [~,base] = fileparts(file);
    missing = [];
    invalid = [];
    if isfield(spec,'missing_values')
        missing = spec.missing_values;
    end
    if isfield(spec,'invalid_le')
        invalid = spec.invalid_le;
    end
    key = strjoin({file,char(string(schema.name)), ...
        mat2str(double(missing(:).')),mat2str(double(invalid(:).'))}, ...
        char(30));
    cacheFile = fullfile(cacheDir,sprintf('wide_%s_%s.mat', ...
        matlab.lang.makeValidName(base),local_signature_hash(key)));
    sourceSignature = local_canonical_source_signature({file},schema);

    if strcmp(action,'load')
        if ~isfile(cacheFile)
            return
        end
        try
            C = load(cacheFile,'cacheVersion','sourceSignature', ...
                'cacheIDs','cacheTime','cacheValues');
            valid = isfield(C,'cacheVersion') ...
                && C.cacheVersion == 1 ...
                && isfield(C,'sourceSignature') ...
                && strcmp(C.sourceSignature,sourceSignature) ...
                && all(isfield(C,{'cacheIDs','cacheTime','cacheValues'}));
            if ~valid
                return
            end
            [foundID,idLocation] = ismember(requestIDs,string(C.cacheIDs));
            [foundTime,timeLocation] = ismember( ...
                posixtime(targetTime(:)),double(C.cacheTime(:)));
            if ~all(foundID)
                return
            end
            X = nan(numel(targetTime),numel(requestIDs),'single');
            X(foundTime,:) = C.cacheValues( ...
                timeLocation(foundTime),idLocation);
            hit = true;
        catch
            hit = false;
            X = [];
        end
        return
    end

    if ~strcmp(action,'save')
        error('read_hydro_timeseries:WideCacheAction', ...
            'Unknown wide cache action: %s.',action);
    end
    if ~isfolder(cacheDir)
        mkdir(cacheDir);
    end
    cacheVersion = 1;
    cacheIDs = string(fileIDs(:));
    cacheTime = posixtime(data.time(:));
    cacheValues = single(data.values);
    temporaryFile = [tempname(cacheDir) '.mat'];
    cleanup = onCleanup(@()local_delete_cache_file(temporaryFile));
    try
        save(temporaryFile,'cacheVersion','sourceSignature', ...
            'cacheIDs','cacheTime','cacheValues');
        [status,message] = movefile(temporaryFile,cacheFile,'f');
        if ~status
            warning('read_hydro_timeseries:WideCacheWrite', ...
                'Could not write wide cache %s: %s',cacheFile,message);
        end
    catch ME
        warning('read_hydro_timeseries:WideCacheWrite', ...
            'Could not write wide cache %s: %s',cacheFile,ME.message);
    end
    clear cleanup
end

function ids = local_clean_wide_ids(ids,idSpec)

    ids = upper(strip(string(ids)));
    ids = erase(ids,char(65279));
    ids = erase(ids,'"');
    ids = regexprep(ids,'\.0+$','');
    if isfield(idSpec,'strip_leading_zeros') ...
            && idSpec.strip_leading_zeros
        numericID = ~cellfun('isempty', ...
            regexp(cellstr(ids),'^[0-9]+$','once'));
        ids(numericID) = regexprep(ids(numericID), ...
            '^0+(?=[0-9])','');
    end
end

function x = local_apply_native_missing(x,spec)
    
    x = double(x);
    
    if isfield(spec,'missing_values') ...
            && ~isempty(spec.missing_values)
    
        mv = double(spec.missing_values(:));
    
        for j = 1:numel(mv)
            x(x == mv(j)) = NaN;
        end
    end
    
    if isfield(spec,'invalid_le') ...
            && ~isempty(spec.invalid_le)
        x(x <= double(spec.invalid_le)) = NaN;
    end
end
