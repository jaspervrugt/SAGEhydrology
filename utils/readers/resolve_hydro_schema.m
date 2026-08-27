function schema = resolve_hydro_schema(schema,split,options)
%RESOLVE_HYDRO_SCHEMA Select and merge a declarative runtime profile.
    if nargin < 3 || isempty(options), options = struct(); end
    if ~isfield(schema,'profiles') || isempty(schema.profiles), return; end
    if ~isstruct(schema.profiles) || ~isscalar(schema.profiles)
        error('resolve_hydro_schema:BadProfiles', ...
            'schema.profiles must be a scalar structure of named profiles.');
    end
    names = fieldnames(schema.profiles);
    matched = strings(0,1);
    for i = 1:numel(names)
        profile = schema.profiles.(names{i});
        if ~isstruct(profile) || ~isscalar(profile) ...
                || ~isfield(profile,'match') || ~isfield(profile,'schema')
            error('resolve_hydro_schema:BadProfile', ...
                'Profile %s requires scalar .match and .schema structures.',names{i});
        end
        if local_matches(profile.match,split,options)
            matched(end+1,1) = string(names{i}); %#ok<AGROW>
        end
    end
    if numel(matched) ~= 1
        error('resolve_hydro_schema:ProfileMatch', ...
            'Expected one matching profile; found %d (%s).', ...
            numel(matched),strjoin(cellstr(matched),', '));
    end
    selected = schema.profiles.(char(matched));
    schema = rmfield(schema,'profiles');
    schema = local_merge(schema,selected.schema);
    schema.selected_profile = char(matched);
end

function tf = local_matches(match,split,options)
    tf = true;
    fields = fieldnames(match);
    for i = 1:numel(fields)
        field = fields{i};
        if isfield(split,field)
            actual = split.(field);
        elseif isfield(options,field)
            actual = options.(field);
        else
            tf = false; return
        end
        expected = match.(field);
        tf = any(arrayfun(@(j) isequaln(actual,expected(j)),1:numel(expected)));
        if ~tf, return; end
    end
end

function out = local_merge(base,override)
    out = base;
    fields = fieldnames(override);
    for i = 1:numel(fields)
        field = fields{i};
        if isfield(out,field) && isstruct(out.(field)) ...
                && isscalar(out.(field)) && isstruct(override.(field)) ...
                && isscalar(override.(field))
            out.(field) = local_merge(out.(field),override.(field));
        else
            out.(field) = override.(field);
        end
    end
end
