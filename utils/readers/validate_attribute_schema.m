function schema = validate_attribute_schema(schema)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VALIDATE_ATTRIBUTE_SCHEMA Validate a declarative attribute schema.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isstruct(schema) ...
            || ~isscalar(schema)
        error('validate_attribute_schema:InvalidSchema', ...
            'Attribute schema must be a scalar structure.');
    end

    required = {'name','tables','id','metadata','region'};
    for i = 1:numel(required)
        field = required{i};
        if ~isfield(schema,field) ...
                || isempty(schema.(field))
            error('validate_attribute_schema:MissingField', ...
                'Attribute schema is missing required field: %s.',field);
        end
    end

    if ~isstruct(schema.tables) ...
            || isempty(schema.tables)
        error('validate_attribute_schema:BadTables', ...
            'schema.tables must be a nonempty structure array.');
    end

    for i = 1:numel(schema.tables)
        tableSchema = schema.tables(i);
        hasFile = isfield(tableSchema,'file') ...
            && ~isempty(tableSchema.file);
        hasPattern = isfield(tableSchema,'pattern') ...
            && ~isempty(tableSchema.pattern);
        if ~hasFile && ~hasPattern
            error('validate_attribute_schema:MissingFile', ...
                'Attribute table %d does not define .file or .pattern.',i);
        end
        if ~isfield(tableSchema,'keys') ...
                || isempty(tableSchema.keys)
            error('validate_attribute_schema:MissingKey', ...
                'Attribute table %d does not define .keys.',i);
        end
    end

    if ~isfield(schema.id,'column') ...
            || isempty(schema.id.column)
        schema.id.column = 'gauge_id';
    end
    if ~isfield(schema.id,'uppercase')
        schema.id.uppercase = false;
    end
    if ~isfield(schema.id,'lowercase')
        schema.id.lowercase = false;
    end
    if ~isfield(schema.id,'strip')
        schema.id.strip = true;
    end
    if ~isfield(schema.id,'regex')
        schema.id.regex = cell(0,2);
    end
    if ~iscell(schema.id.regex) ...
            || (~isempty(schema.id.regex) ...
            && size(schema.id.regex,2) ~= 2)
        error('validate_attribute_schema:BadIdRegex', ...
            'schema.id.regex must be an N-by-2 cell array.');
    end
    if ~isfield(schema.id,'output_uppercase')
        schema.id.output_uppercase = schema.id.uppercase;
    end
    if ~isfield(schema.id,'output_lowercase')
        schema.id.output_lowercase = schema.id.lowercase;
    end
    if ~isfield(schema.id,'output_regex')
        schema.id.output_regex = cell(0,2);
    end
    if ~isfield(schema.id,'pad_width')
        schema.id.pad_width = 0;
    end
    if ~isfield(schema.id,'numeric_canonical')
        schema.id.numeric_canonical = false;
    end
    if ~isfield(schema.id,'sort')
        schema.id.sort = 'text';
    end
    if ~isfield(schema.id,'output_type')
        schema.id.output_type = 'string';
    end
    if ~iscell(schema.id.output_regex) ...
            || (~isempty(schema.id.output_regex) ...
            && size(schema.id.output_regex,2) ~= 2)
        error('validate_attribute_schema:BadOutputIdRegex', ...
            'schema.id.output_regex must be an N-by-2 cell array.');
    end

    if ~isfield(schema,'aliases') ...
            || isempty(schema.aliases)
        schema.aliases = repmat(struct( ...
            'target','','sources',{{}},'required',false, ...
            'default',NaN),0,1);
    end
    if ~isstruct(schema.aliases)
        error('validate_attribute_schema:BadAliases', ...
            'schema.aliases must be a structure array.');
    end

    if ~isfield(schema.metadata,'name_sources') ...
            || isempty(schema.metadata.name_sources)
        schema.metadata.name_sources = {schema.id.column};
    end
    if ~isfield(schema.metadata,'name_transform')
        schema.metadata.name_transform = '';
    end
    if ~isfield(schema,'zone') ...
            || isempty(schema.zone)
        schema.zone = struct('region',schema.region);
    elseif ~isfield(schema.zone,'region') ...
            || isempty(schema.zone.region)
        schema.zone.region = schema.region;
    end
    if ~isfield(schema,'progress')
        schema.progress = struct();
    end
    if ~isfield(schema,'selection')
        schema.selection = struct();
    end
end
