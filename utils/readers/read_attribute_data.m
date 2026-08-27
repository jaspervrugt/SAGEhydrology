function [A,id_gauge,gname,zone] = read_attribute_data( ...
    dirD,bas,schema)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_ATTRIBUTE_DATA Read attributes using a declarative regional schema.
%
% This reader is intentionally independent of src/read_attr.m while the
% regional migration gates are being completed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(bas)
        bas = struct();
    end
    schema = local_resolve_attribute_profile(schema,bas);
    schema = validate_attribute_schema(schema);
    dirA = local_attribute_root(dirD,schema);
    tableSchema = local_expand_tables(dirA,schema.tables);

    tables = cell(numel(tableSchema),1);
    for i = 1:numel(tableSchema)
        tables{i} = local_read_table(dirA,tableSchema(i),schema.id);
    end
    attributes = tables{1};
    for i = 2:numel(tables)
        if isempty(tables{i})
            continue
        end
        joinType = 'left';
        if isfield(tableSchema(i),'join_type') ...
                && ~isempty(tableSchema(i).join_type)
            joinType = tableSchema(i).join_type;
        end
        attributes = local_join( ...
            attributes,tables{i},local_join_key( ...
            tableSchema(i),schema.id.column),joinType, ...
            tableSchema(i));
    end

    for i = 1:numel(schema.aliases)
        attributes = local_add_alias(attributes,schema.aliases(i));
    end
    if isfield(schema,'projection') ...
            && ~isempty(schema.projection)
        attributes = local_project_coordinates( ...
            attributes,schema.projection);
    end

    idAll = local_normalize_id( ...
        attributes.(schema.id.column),schema.id);
    if strcmpi(schema.id.sort,'numeric')
        numericId = str2double(idAll);
        if all(isfinite(numericId))
            [~,order] = sort(numericId);
        else
            [~,order] = sort(idAll);
        end
    else
        [~,order] = sort(idAll);
    end
    idAll = idAll(order);
    attributes = attributes(order,:);

    if isfield(schema.metadata,'name_components') ...
            && ~isempty(schema.metadata.name_components)
        nameAll = local_composed_names( ...
            attributes,idAll,schema.metadata);
    else
        nameColumn = local_find_column( ...
            attributes,schema.metadata.name_sources,true);
        if isempty(nameColumn)
            prefix = '';
            if isfield(schema.metadata,'name_fallback_prefix')
                prefix = char(string( ...
                    schema.metadata.name_fallback_prefix));
            end
            nameAll = string(prefix) + idAll;
        else
            nameAll = string(attributes.(nameColumn));
        end
    end
    outputIdAll = local_output_id(idAll,schema.id);
    standardizeNames = true;
    if isfield(schema.metadata,'standardize') ...
            && ~schema.metadata.standardize
        standardizeNames = false;
    end
    if standardizeNames
        nameAll = standardize_camels_gauge_names( ...
            nameAll,schema.region,outputIdAll);
    end
    if strcmpi(schema.metadata.name_transform,'title_case')
        nameAll = local_title_case(nameAll,true);
    elseif strcmpi(schema.metadata.name_transform,'title_case_words')
        nameAll = local_title_case(nameAll,false);
    end

    local_write_gauge_information( ...
        dirD,attributes,outputIdAll,nameAll,schema);

    if (~isfield(bas,'id_gauge') ...
            || isempty(bas.id_gauge)) ...
            && isfield(schema.selection,'default_file') ...
            && ~isempty(schema.selection.default_file)
        selectionFile = fullfile( ...
            dirA,schema.selection.default_file);
        if isfile(selectionFile)
            bas.id_gauge = readmatrix(selectionFile);
        end
    end
    if (~isfield(bas,'id_gauge') ...
            || isempty(bas.id_gauge)) ...
            && isfield(schema.selection,'available_pattern') ...
            && ~isempty(schema.selection.available_pattern)
        files = dir(fullfile(dirA, ...
            schema.selection.available_pattern));
        selected = strings(0,1);
        expression = schema.selection.available_regex;
        for i = 1:numel(files)
            token = regexp(files(i).name,expression,'tokens','once');
            if ~isempty(token)
                selected(end+1,1) = string(token{1}); %#ok<AGROW>
            end
        end
        if ~isempty(selected)
            numericSelected = str2double(selected);
            if all(isfinite(numericSelected))
                [~,order] = sort(numericSelected);
                selected = selected(order);
            else
                selected = sort(selected);
            end
            bas.id_gauge = selected;
        end
    end
    useRequested = true;
    if isfield(schema.selection,'use_requested') ...
            && ~schema.selection.use_requested
        useRequested = false;
    end
    if useRequested ...
            && isfield(bas,'id_gauge') ...
            && ~isempty(bas.id_gauge)
        idRequested = local_normalize_id(bas.id_gauge(:),schema.id);
        idRequestedMatch = idRequested;
        idAllMatch = idAll;
        if isfield(schema.id,'optional_prefix') ...
                && ~isempty(schema.id.optional_prefix)
            prefix = string(schema.id.optional_prefix);
            if any(startsWith(idAllMatch,prefix)) ...
                    && ~all(startsWith(idRequestedMatch,prefix))
                idRequestedMatch(~startsWith( ...
                    idRequestedMatch,prefix)) = prefix + ...
                    idRequestedMatch(~startsWith( ...
                    idRequestedMatch,prefix));
            elseif any(startsWith(idRequestedMatch,prefix)) ...
                    && ~all(startsWith(idAllMatch,prefix))
                idAllMatch(~startsWith(idAllMatch,prefix)) = ...
                    prefix + idAllMatch(~startsWith(idAllMatch,prefix));
            end
        end
        [found,location] = ismember(idRequestedMatch,idAllMatch);
        if ~all(found)
            first = find(~found,1);
            error('read_attribute_data:MissingGauge', ...
                'Requested gauge is absent from %s: %s.', ...
                schema.name,idRequested(first));
        end
        attributes = attributes(location,:);
        id_gauge = local_output_id(idRequested,schema.id);
        gname = nameAll(location);
    else
        id_gauge = outputIdAll;
        gname = nameAll;
    end

    minPerZone = local_option(bas,'min_per_zone',5);
    maxZones = local_option(bas,'max_zones',8);
    zone = classify_camels_all_zones( ...
        schema.zone.region,attributes,minPerZone,maxZones);

    if ~isfield(bas,'id_attr') ...
            || isempty(bas.id_attr)
        error('read_attribute_data:MissingAttributeSelection', ...
            'bas.id_attr must contain selected attribute indices.');
    end
    printTable = local_option(bas,'pr_attr',1);
    [raw,~] = select_attr_from_catalog( ...
        attributes,bas.id_attr,schema.region,printTable);
    mu = mean(raw,1,'omitnan');
    sigma = std(raw,0,1,'omitnan');
    sigma(~isfinite(sigma) | sigma == 0) = 1;
    A = ((raw - mu)./sigma).';

    if isfield(schema.progress,'label') ...
            && ~isempty(schema.progress.label)
        fprintf('%s ... Done\n',schema.progress.label);
    end
end

function schema = local_resolve_attribute_profile(schema,bas)
    if ~isfield(schema,'profiles') ...
            || isempty(schema.profiles)
        return
    end
    stream = '';
    if isfield(bas,'stream') ...
            && ~isempty(bas.stream)
        stream = lower(char(string(bas.stream)));
    elseif isfield(bas,'dt') ...
            && ~isempty(bas.dt)
        dt = double(bas.dt);
        if abs(dt - 24) < 1e-8 ...
                || abs(dt - 1/24) < 1e-8 ...
                || abs(dt - 3600) < 1e-8
            stream = 'hourly';
        else
            stream = 'daily';
        end
    end
    if isempty(stream)
        return
    end
    names = fieldnames(schema.profiles);
    for i = 1:numel(names)
        profile = schema.profiles.(names{i});
        if isfield(profile,'match') ...
                && isfield(profile.match,'stream') ...
                && strcmpi(profile.match.stream,stream)
            schema = profile.schema;
            return
        end
    end
end

function expanded = local_expand_tables(dirD,tables)
    expanded = tables([]);
    for i = 1:numel(tables)
        S = tables(i);
        if ~isfield(S,'pattern') ...
                || isempty(S.pattern)
            expanded(end+1,1) = S; %#ok<AGROW>
            continue
        end
        files = dir(fullfile(dirD,S.pattern));
        [~,order] = sort({files.name});
        files = files(order);
        for j = 1:numel(files)
            item = S;
            item.file = files(j).name;
            if isfield(item,'prefix_from_file') ...
                    && item.prefix_from_file
                [~,base] = fileparts(files(j).name);
                base = regexprep(base, ...
                    '^catchments_hydrological_signatures_?','hs_');
                base = regexprep(base,'[^A-Za-z0-9]+','_');
                base = regexprep(base,'_+$','');
                item.prefix_columns = ...
                    [matlab.lang.makeValidName(base) '_'];
            end
            expanded(end+1,1) = item; %#ok<AGROW>
        end
    end
end

function dirA = local_attribute_root(dirD,schema)
    dirA = dirD;
    if ~isfield(schema,'root_candidates') ...
            || isempty(schema.root_candidates)
        return
    end
    candidates = cellstr(string(schema.root_candidates));
    for i = 1:numel(candidates)
        candidate = fullfile(dirD,candidates{i});
        found = true;
        for j = 1:numel(schema.tables)
            if ~isfile(fullfile(candidate,schema.tables(j).file))
                found = false;
                break
            end
        end
        if found
            dirA = candidate;
            return
        end
    end
end

function T = local_read_table(dirD,S,idSchema)
    file = fullfile(dirD,S.file);
    if ~isfile(file)
        if isfield(S,'required') ...
                && ~S.required
            T = [];
            return
        end
        error('read_attribute_data:MissingFile', ...
            'Missing attribute file: %s.',file);
    end
    if isfield(S,'layout') ...
            && strcmpi(S.layout,'transposed_attributes')
        T = local_read_transposed(file,S);
    elseif isfield(S,'layout') ...
            && strcmpi(S.layout,'catalog_attribute_matrix')
        T = local_read_catalog_matrix(file,S,idSchema);
    else
        T = local_read_native_table(file,S);
    end
    if isfield(S,'column_renames') ...
            && ~isempty(S.column_renames)
        renames = S.column_renames;
        for i = 1:size(renames,1)
            location = find(strcmp( ...
                T.Properties.VariableNames,renames{i,1}),1);
            if ~isempty(location)
                T.Properties.VariableNames{location} = renames{i,2};
            end
        end
    end
    if isfield(S,'make_valid_names') ...
            && S.make_valid_names
        replacement = 'delete';
        if isfield(S,'valid_name_replacement') ...
                && ~isempty(S.valid_name_replacement)
            replacement = S.valid_name_replacement;
        end
        T.Properties.VariableNames = matlab.lang.makeValidName( ...
            T.Properties.VariableNames,'ReplacementStyle',replacement);
    end
    if isfield(S,'make_unique_names') ...
            && S.make_unique_names
        T.Properties.VariableNames = matlab.lang.makeUniqueStrings( ...
            T.Properties.VariableNames);
    end
    if isfield(S,'row_filter') ...
            && isstruct(S.row_filter) ...
            && isfield(S.row_filter,'column') ...
            && ~isempty(S.row_filter.column)
        column = local_find_column( ...
            T,{S.row_filter.column},false);
        value = strtrim(string(T.(column)));
        keep = strcmpi(value,string(S.row_filter.value));
        T = T(keep,:);
    end
    if isfield(S,'keep_columns') ...
            && ~isempty(S.keep_columns)
        keep = cellstr(string(S.keep_columns));
        missing = setdiff(keep,T.Properties.VariableNames,'stable');
        if ~isempty(missing)
            error('read_attribute_data:MissingKeptColumn', ...
                'Missing requested attribute columns: %s.', ...
                strjoin(missing,', '));
        end
        T = T(:,keep);
    end
    if isfield(S,'numeric_text') ...
            && any(strcmpi(S.numeric_text,{'auto','threshold'}))
        exceptions = strings(0,1);
        if isfield(S,'numeric_text_exceptions')
            exceptions = string(S.numeric_text_exceptions);
        end
        for i = 1:width(T)
            name = T.Properties.VariableNames{i};
            if any(strcmpi(name,exceptions))
                continue
            end
            value = T.(name);
            if iscell(value) ...
                    || isstring(value) ...
                    || iscategorical(value)
                if strcmpi(S.numeric_text,'auto')
                    T.(name) = local_numeric(value);
                else
                    text = strtrim(string(value));
                    text = erase(text,'"');
                    text = strrep(text,',','.');
                    missingToken = ismissing(text) ...
                        | text == "" ...
                        | strcmpi(text,"NA") ...
                        | strcmpi(text,"NAN") ...
                        | strcmpi(text,"NULL") ...
                        | strcmpi(text,"NONE") ...
                        | text == "-9999";
                    text = regexprep(text, ...
                        '^(NA|NaN|null|NULL|nan|None|-9999)$', ...
                        'NaN','ignorecase');
                    converted = str2double(text);
                    present = ~missingToken;
                    finite = isfinite(converted);
                    ratio = 0.9;
                    if isfield(S,'numeric_text_ratio')
                        ratio = S.numeric_text_ratio;
                    end
                    if any(present) ...
                            && sum(finite(present)) ...
                            >= ratio*sum(present)
                        T.(name) = converted;
                    end
                end
            end
        end
    end
    if isfield(S,'boolean_text') ...
            && strcmpi(S.boolean_text,'auto')
        T = local_boolean_text(T);
    end
    if isfield(S,'drop_empty_columns') ...
            && S.drop_empty_columns
        keep = true(1,width(T));
        for i = 1:width(T)
            value = T.(i);
            if isnumeric(value)
                keep(i) = ~all(isnan(value));
            elseif iscellstr(value) ...
                    || isstring(value)
                keep(i) = any(strlength(strtrim(string(value))) > 0);
            end
        end
        T = T(:,keep);
    end
    if isfield(S,'drop_empty_rows') ...
            && S.drop_empty_rows ...
            && height(T) > 0
        empty = true(height(T),1);
        for i = 1:width(T)
            value = T.(i);
            if isnumeric(value) ...
                    || islogical(value)
                empty = empty & isnan(double(value));
            else
                value = string(value);
                empty = empty & (ismissing(value) ...
                    | strlength(strtrim(value)) == 0);
            end
        end
        T(empty,:) = [];
    end
    key = local_find_column(T,S.keys,false);
    targetKey = idSchema.column;
    if isfield(S,'key_target') ...
            && ~isempty(S.key_target)
        targetKey = S.key_target;
    end
    T.(targetKey) = local_normalize_id(T.(key),idSchema);
    if ~strcmp(key,targetKey)
        T(:,key) = [];
    end
    if isfield(S,'widen') ...
            && isstruct(S.widen) ...
            && isfield(S.widen,'keys') ...
            && ~isempty(S.widen.keys)
        T = local_widen_table(T,targetKey,S.widen.keys);
    end
    if isfield(S,'drop_missing_key') ...
            && S.drop_missing_key
        id = string(T.(targetKey));
        T(ismissing(id) | id == "",:) = [];
    end
    if isfield(S,'unique_key') ...
            && S.unique_key
        [~,keep] = unique(T.(targetKey),'stable');
        T = T(keep,:);
    end
    if isfield(S,'exclude_columns_regex') ...
            && ~isempty(S.exclude_columns_regex)
        variables = T.Properties.VariableNames;
        remove = false(size(variables));
        for i = 1:numel(variables)
            remove(i) = ~strcmp(variables{i},targetKey) ...
                && ~isempty(regexp(variables{i}, ...
                S.exclude_columns_regex,'once'));
        end
        T(:,remove) = [];
    end
    if isfield(S,'prefix_columns') ...
            && ~isempty(S.prefix_columns)
        names = T.Properties.VariableNames;
        for i = 1:numel(names)
            if strcmp(names{i},targetKey) ...
                    || startsWith(names{i},S.prefix_columns)
                continue
            end
            T.Properties.VariableNames{i} = ...
                [S.prefix_columns names{i}];
        end
    end
    if isfield(S,'absolute_columns') ...
            && ~isempty(S.absolute_columns)
        names = cellstr(string(S.absolute_columns));
        for i = 1:numel(names)
            if ismember(names{i},T.Properties.VariableNames)
                T.(names{i}) = abs(local_numeric(T.(names{i})));
            end
        end
    end
    if isfield(S,'fallback_file') ...
            && ~isempty(S.fallback_file) ...
            && isfile(fullfile(dirD,S.fallback_file))
        fallbackSchema = S;
        fallbackSchema.file = S.fallback_file;
        fallbackSchema = rmfield(fallbackSchema,'fallback_file');
        fallback = local_read_table(dirD,fallbackSchema,idSchema);
        T = local_coalesce_table(T,fallback,idSchema.column);
    end
end

function key = local_join_key(S,defaultKey)
    key = defaultKey;
    if isfield(S,'join_key') ...
            && ~isempty(S.join_key)
        key = S.join_key;
    end
end

function W = local_widen_table(T,idName,keyNames)
    keyNames = cellstr(string(keyNames));
    keyNames = keyNames(ismember(keyNames,T.Properties.VariableNames));
    ids = T.(idName);
    [uniqueIds,~,group] = unique(ids);
    if numel(uniqueIds) == height(T)
        W = T;
        return
    end
    variables = T.Properties.VariableNames;
    dataVariables = setdiff(variables, ...
        [{idName},keyNames],'stable');
    suffix = strings(height(T),1);
    for i = 1:height(T)
        parts = strings(1,numel(keyNames));
        for j = 1:numel(keyNames)
            value = T.(keyNames{j});
            parts(j) = string(value(i));
        end
        suffix(i) = local_clean_suffix(strjoin(parts,"_"));
    end
    raw = strings(height(T)*numel(dataVariables),1);
    count = 0;
    for i = 1:height(T)
        for j = 1:numel(dataVariables)
            count = count + 1;
            raw(count) = string(dataVariables{j}) + "_" + suffix(i);
        end
    end
    raw = raw(1:count);
    [uniqueRaw,~,map] = unique(raw);
    names = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName( ...
        cellstr(uniqueRaw),'ReplacementStyle','delete'));
    W = table(uniqueIds,'VariableNames',{idName});
    for i = 1:numel(names)
        W.(names{i}) = nan(numel(uniqueIds),1);
    end
    count = 0;
    for i = 1:height(T)
        for j = 1:numel(dataVariables)
            count = count + 1;
            column = names{map(count)};
            value = local_numeric(T.(dataVariables{j})(i));
            if ~isempty(value)
                W.(column)(group(i)) = value(1);
            end
        end
    end
end

function suffix = local_clean_suffix(value)
    suffix = lower(char(string(value)));
    suffix = regexprep(suffix,'[^a-zA-Z0-9]+','_');
    suffix = regexprep(suffix,'(^_+|_+$)','');
    if isempty(suffix)
        suffix = 'x';
    end
end

function T = local_coalesce_table(T,fallback,key)
    [found,location] = ismember(T.(key),fallback.(key));
    names = fallback.Properties.VariableNames;
    names(strcmp(names,key)) = [];
    for i = 1:numel(names)
        name = names{i};
        source = fallback.(name);
        if ~ismember(name,T.Properties.VariableNames)
            if isnumeric(source) ...
                    || islogical(source)
                value = nan(height(T),1);
                value(found) = double(source(location(found)));
            else
                value = strings(height(T),1);
                value(found) = string(source(location(found)));
            end
            T.(name) = value;
            continue
        end
        rows = find(found);
        if isnumeric(T.(name)) ...
                || islogical(T.(name))
            value = double(T.(name));
            sourceValue = double(source(location(found)));
            missing = ~isfinite(value(rows)) & isfinite(sourceValue);
            value(rows(missing)) = sourceValue(missing);
        else
            value = string(T.(name));
            sourceValue = string(source(location(found)));
            missing = ismissing(value(rows)) ...
                | strlength(strtrim(value(rows))) == 0;
            usable = ~ismissing(sourceValue) ...
                & strlength(strtrim(sourceValue)) > 0;
            use = missing & usable;
            value(rows(use)) = sourceValue(use);
        end
        T.(name) = value;
    end
end

function T = local_read_native_table(file,S)
    readArguments = {'VariableNamingRule','preserve'};
    if isfield(S,'sheet') ...
            && ~isempty(S.sheet)
        readArguments = [readArguments ...
            {'Sheet',S.sheet}];
    end
    if isfield(S,'delimiter') ...
            && ~isempty(S.delimiter)
        readArguments = [readArguments ...
            {'Delimiter',S.delimiter}];
    end
    if isfield(S,'encoding') ...
            && ~isempty(S.encoding)
        readArguments = [readArguments ...
            {'Encoding',S.encoding}];
    end
    try
        options = detectImportOptions(file,readArguments{:});
    catch
        options = detectImportOptions(file, ...
            'VariableNamingRule','preserve');
    end
    if isfield(S,'key_type') ...
            && ~isempty(S.key_type)
        optionKey = '';
        optionNames = options.VariableNames;
        candidates = cellstr(string(S.keys));
        for i = 1:numel(candidates)
            location = find(strcmpi(optionNames,candidates{i}),1);
            if ~isempty(location)
                optionKey = optionNames{location};
                break
            end
        end
        if ~isempty(optionKey)
            options = setvartype(options,optionKey,S.key_type);
        end
    end
    if isfield(S,'comment_style') ...
            && ~isempty(S.comment_style)
        try
            options.CommentStyle = S.comment_style;
        catch
        end
    end
    try
        T = readtable(file,options,readArguments{:});
    catch
        T = readtable(file,options);
    end
end

function T = local_read_transposed(file,S)
    delimiter = sprintf('\t');
    if isfield(S,'delimiter') ...
            && ~isempty(S.delimiter)
        delimiter = S.delimiter;
    end
    cells = readcell(file,'FileType','text', ...
        'Delimiter',delimiter,'TextType','string');
    names = string(cells(:,1));
    names = erase(strip(names),'"');
    names = matlab.lang.makeValidName(names, ...
        'ReplacementStyle','delete');
    names = matlab.lang.makeUniqueStrings(names);
    T = table();
    for i = 1:numel(names)
        value = string(cells(i,2:end).');
        value = erase(strip(value),'"');
        value(ismissing(value)) = "";
        numeric = str2double(value);
        present = strlength(value) > 0;
        if any(present) ...
                && all(isfinite(numeric(present)))
            T.(names{i}) = numeric;
        else
            T.(names{i}) = value;
        end
    end
end

function T = local_read_catalog_matrix(file,S,idSchema)
    source = readtable(file,'VariableNamingRule','preserve');
    attributeColumn = local_find_column( ...
        source,{S.attribute_column},false);
    rawNames = strtrim(string(source.(attributeColumn)));
    validNames = string(matlab.lang.makeValidName(rawNames));
    firstBasin = double(S.metadata_columns) + 1;
    basinColumns = source.Properties.VariableNames(firstBasin:end);
    id = local_normalize_id(string(basinColumns),idSchema);
    T = table(id,'VariableNames',{idSchema.column});
    catalog = attr_catalog(S.catalog_region);
    wanted = string(catalog.names(:));
    if isfield(S,'exclude_catalog') ...
            && ~isempty(S.exclude_catalog)
        wanted = wanted(~ismember(wanted,string(S.exclude_catalog)));
    end
    for i = 1:numel(wanted)
        name = wanted(i);
        derived = [];
        if isfield(S,'derived') ...
                && isstruct(S.derived) ...
                && ~isempty(S.derived)
            target = string({S.derived.target});
            location = find(target == name,1);
            if ~isempty(location)
                derived = S.derived(location);
            end
        end
        if ~isempty(derived) ...
                && strcmpi(derived.operation,'row_sum')
            components = string(derived.rows);
            values = nan(numel(components),numel(basinColumns));
            for j = 1:numel(components)
                row = find(validNames == components(j),1);
                if isempty(row)
                    error('read_attribute_data:MissingMatrixRow', ...
                        'Missing attribute-matrix row: %s.',components(j));
                end
                values(j,:) = local_matrix_row( ...
                    source{row,firstBasin:end});
            end
            value = sum(values,1,'omitnan');
            value(all(~isfinite(values),1)) = NaN;
        else
            row = find(validNames == name,1);
            if isempty(row)
                error('read_attribute_data:MissingMatrixRow', ...
                    'Missing attribute-matrix row: %s.',name);
            end
            value = local_matrix_row(source{row,firstBasin:end});
        end
        T.(char(name)) = value(:);
    end
end

function value = local_matrix_row(value)
    if iscell(value)
        output = nan(numel(value),1);
        for i = 1:numel(value)
            item = value{i};
            if isnumeric(item) ...
                    || islogical(item)
                output(i) = double(item(1));
            else
                output(i) = str2double(strrep( ...
                    string(item),',','.'));
            end
        end
        value = output;
    elseif isnumeric(value) ...
            || islogical(value)
        value = double(value(:));
    else
        value = str2double(strrep(string(value(:)),',','.'));
    end
    value = value(:).';
end

function T = local_boolean_text(T)
    for i = 1:width(T)
        value = T{:,i};
        if ~(iscell(value) ...
                || isstring(value) ...
                || iscategorical(value))
            continue
        end
        value = lower(strtrim(string(value)));
        present = ~ismissing(value) & value ~= "";
        if all(ismember(value(present),["true","false","0","1"]))
            T.(T.Properties.VariableNames{i}) = ...
                value == "true" | value == "1";
        end
    end
end

function left = local_join(left,right,key,joinType,S)
    duplicates = intersect( ...
        left.Properties.VariableNames,right.Properties.VariableNames);
    duplicates(strcmp(duplicates,key)) = [];
    policy = 'drop';
    if isfield(S,'duplicate_policy') ...
            && ~isempty(S.duplicate_policy)
        policy = lower(char(string(S.duplicate_policy)));
    end
    if strcmp(policy,'suffix')
        suffix = '_source';
        if isfield(S,'duplicate_suffix') ...
                && ~isempty(S.duplicate_suffix)
            suffix = char(string(S.duplicate_suffix));
        end
        for i = 1:numel(duplicates)
            old = duplicates{i};
            right.Properties.VariableNames{strcmp( ...
                right.Properties.VariableNames,old)} = [old suffix];
        end
    elseif strcmp(policy,'replace_left')
        left(:,duplicates) = [];
    else
        right(:,duplicates) = [];
    end
    if strcmpi(joinType,'inner')
        left = innerjoin(left,right,'Keys',key);
    else
        left = outerjoin(left,right,'Keys',key, ...
            'MergeKeys',true,'Type',joinType);
    end
end

function T = local_project_coordinates(T,S)
    x = local_numeric(T.(S.x));
    y = local_numeric(T.(S.y));
    geographic = abs(x) <= 180 & abs(y) <= 90;
    latitude = y;
    longitude = x;
    project = ~geographic;
    if any(project)
        crs = projcrs(S.epsg);
        [latitude(project),longitude(project)] = ...
            projinv(crs,x(project),y(project));
    end
    T.(S.latitude_target) = latitude;
    T.(S.longitude_target) = longitude;
end

function names = local_composed_names(T,id,metadata)
    components = cellstr(string(metadata.name_components));
    separator = ': ';
    if isfield(metadata,'name_separator')
        separator = char(string(metadata.name_separator));
    end
    prefix = '';
    if isfield(metadata,'name_fallback_prefix')
        prefix = char(string(metadata.name_fallback_prefix));
    end
    names = string(prefix) + string(id);
    first = local_find_column(T,components(1),true);
    if isempty(first)
        return
    end
    primary = strtrim(string(T.(first)));
    good = ~ismissing(primary) & primary ~= "";
    names(good) = primary(good);
    for i = 2:numel(components)
        column = local_find_column(T,components(i),true);
        if isempty(column)
            continue
        end
        value = strtrim(string(T.(column)));
        append = good & ~ismissing(value) & value ~= "";
        names(append) = names(append) + separator + value(append);
    end
    if isfield(metadata,'parenthetical_components') ...
            && ~isempty(metadata.parenthetical_components)
        secondary = local_component_names(T, ...
            metadata.parenthetical_components,separator);
        append = strlength(secondary) > 0;
        names(append) = names(append) + " (" ...
            + secondary(append) + ")";
    end
end

function names = local_component_names(T,components,separator)
    components = cellstr(string(components));
    names = strings(height(T),1);
    for i = 1:numel(components)
        column = local_find_column(T,components(i),true);
        if isempty(column)
            continue
        end
        value = strtrim(string(T.(column)));
        good = ~ismissing(value) & value ~= "";
        join = good & names ~= "";
        names(join) = names(join) + separator + value(join);
        first = good & names == "";
        names(first) = value(first);
    end
end

function T = local_add_alias(T,S)
    if ismember(S.target,T.Properties.VariableNames)
        return
    end
    source = local_find_column(T,S.sources,~S.required);
    if isempty(source)
        T.(S.target) = repmat(double(S.default),height(T),1);
    else
        if isfield(S,'type') ...
                && strcmpi(S.type,'preserve')
            T.(S.target) = T.(source);
        else
            T.(S.target) = local_numeric(T.(source));
        end
    end
    if isfield(S,'scale') ...
            && ~isempty(S.scale)
        T.(S.target) = T.(S.target)*S.scale;
    end
    if isfield(S,'divisor') ...
            && ~isempty(S.divisor)
        T.(S.target) = T.(S.target)/S.divisor;
    end
end

function name = local_find_column(T,candidates,optional)
    candidates = cellstr(string(candidates));
    variables = T.Properties.VariableNames;
    name = '';
    for i = 1:numel(candidates)
        location = find(strcmpi(variables,candidates{i}),1);
        if ~isempty(location)
            name = variables{location};
            return
        end
    end
    if ~optional
        error('read_attribute_data:MissingColumn', ...
            'Missing required column (%s).',strjoin(candidates,', '));
    end
end

function id = local_normalize_id(value,S)
    id = string(value(:));
    if S.strip
        id = strip(id);
    end
    if S.uppercase
        id = upper(id);
    end
    if S.lowercase
        id = lower(id);
    end
    for i = 1:size(S.regex,1)
        id = regexprep(id,S.regex{i,1},S.regex{i,2});
    end
    if S.numeric_canonical
        numericId = str2double(id);
        valid = isfinite(numericId) & numericId > 0;
        id(valid) = compose('%.0f',numericId(valid));
    end
    if S.pad_width > 0
        id = arrayfun(@(value) pad(value,S.pad_width, ...
            'left','0'),id);
    end
end

function id = local_output_id(id,S)
    id = string(id(:));
    for i = 1:size(S.output_regex,1)
        id = regexprep(id,S.output_regex{i,1},S.output_regex{i,2});
    end
    if S.output_uppercase
        id = upper(id);
    end
    if S.output_lowercase
        id = lower(id);
    end
    if strcmpi(S.output_type,'double')
        id = str2double(id);
    end
end

function value = local_numeric(value)
    if isnumeric(value) ...
            || islogical(value)
        value = double(value);
        return
    end
    text = string(value);
    text(ismissing(text) | text == "") = "NaN";
    value = str2double(text);
end

function local_write_gauge_information( ...
    dirD,attributes,id_gauge,gname,schema)

    outputFile = fullfile(dirD,'gauge_information.txt');
    enabled = true;
    if isfield(schema,'gauge_information') ...
            && isfield(schema.gauge_information,'write')
        enabled = logical(schema.gauge_information.write);
    end
    if ~enabled
        return
    end

    % Keep a complete existing lookup, but rebuild files created from an
    % earlier screened subset. This lets a region move to a universal basin
    % inventory without asking users to delete gauge_information.txt first.
    if isfile(outputFile)
        try
            existing = readtable(outputFile, ...
                'TextType','string','VariableNamingRule','preserve');
            idColumn = local_find_column(existing, ...
                {'gauge_id','gauge','id'},true);
            if ~isempty(idColumn)
                existingId = strip(string(existing.(idColumn)));
                existingId = regexprep(existingId,'\.0+$','');
                requestedId = strip(string(id_gauge(:)));
                requestedId = regexprep(requestedId,'\.0+$','');
                if all(ismember(requestedId,existingId))
                    return
                end
            end
        catch
            % A malformed or unreadable lookup is replaced below.
        end
    end

    latitude = local_optional_numeric_column(attributes,{ ...
        'gauge_lat_dd','gauge_lat','gauge_latitude','latitude','lat', ...
        'station_lat','station_latitude','LAT_GAGE'});
    longitude = local_optional_numeric_column(attributes,{ ...
        'gauge_lon_dd','gauge_lon','gauge_longitude','longitude','lon','long', ...
        'station_lon','station_longitude','LONG_GAGE'});
    elevation = local_optional_numeric_column(attributes,{ ...
        'gauge_elev','gauge_elevation','elevation','elev', ...
        'station_elev','station_elevation','ELEV_GAGE_M'});
    area = local_optional_numeric_column(attributes,{ ...
        'area_km2','area_calc','area','catchment_area_km2', ...
        'basin_area_km2','drainage_area_km2','DRAIN_SQKM'});
    if all(isnan(area))
        area = local_optional_numeric_column(attributes,{ ...
            'area_m2','catchment_area_m2','basin_area_m2'})/1e6;
    end

    gauge_id = string(id_gauge(:));
    gauge_name = string(gname(:));
    gauge_lat = latitude;
    gauge_lon = longitude;
    gauge_elev = elevation;
    area_km2 = area;
    overview = table(gauge_id,gauge_name,gauge_lat,gauge_lon, ...
        gauge_elev,area_km2);

    temporaryFile = [tempname(dirD) '.txt'];
    cleanup = onCleanup(@()local_delete_file(temporaryFile));
    try
        writetable(overview,temporaryFile, ...
            'Delimiter','\t','FileType','text');
        [status,message] = movefile(temporaryFile,outputFile,'f');
        if ~status
            warning('read_attribute_data:GaugeInformationWrite', ...
                'Could not write %s: %s',outputFile,message);
            return
        end
        fprintf('      Wrote gauge information file: %s\n', ...
            outputFile);
    catch ME
        warning('read_attribute_data:GaugeInformationWrite', ...
            'Could not write %s: %s',outputFile,ME.message);
    end
    clear cleanup
end

function value = local_optional_numeric_column(T,candidates)
    column = local_find_column(T,candidates,true);
    if isempty(column)
        value = nan(height(T),1);
    else
        value = local_numeric(T.(column));
        value = value(:);
    end
end

function local_delete_file(fileName)
    if isfile(fileName)
        delete(fileName);
    end
end

function value = local_option(S,field,default)
    value = default;
    if isfield(S,field) ...
            && ~isempty(S.(field))
        value = S.(field);
    end
end

function names = local_title_case(names,capitalizeAfterApostrophe)
    names = string(names(:));
    for i = 1:numel(names)
        value = char(lower(strip(names(i))));
        capitalize = true;
        for j = 1:numel(value)
            if isstrprop(value(j),'alpha')
                if capitalize
                    value(j) = upper(value(j));
                end
                capitalize = false;
            elseif any(value(j) == [' ' '-' '/' '('])
                capitalize = true;
            elseif value(j) == '''' ...
                    && capitalizeAfterApostrophe
                capitalize = true;
            end
        end
        names(i) = string(value);
    end
end
