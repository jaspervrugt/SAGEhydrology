function schema = validate_hydro_schema(schema,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VALIDATE_HYDRO_SCHEMA Validate a generic SAGE time-series schema.
%
% SYNOPSIS:
%   schema = validate_hydro_schema(schema)
%   schema = validate_hydro_schema(schema,mode)
%
% INPUT:
%   schema      Structure describing files, time, identifiers, variables,
%               transformations, auxiliary metadata, and progress text
%   mode        OPTIONAL reader mode: 'meteo', 'Q', or 'both'
%
% OUTPUT:
%   schema      Validated schema with defaults inserted
%
% DESCRIPTION:
%   This function validates the declarative configuration consumed by
%   read_hydro_timeseries. Dataset-specific behavior belongs in schema;
%   file reading, alignment, unit conversion, and quality control remain
%   centralized in the generic reader.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(mode)
        mode = 'both';
    end
    if ~isstruct(schema) ...
            || ~isscalar(schema)
        error('validate_hydro_schema:InvalidSchema', ...
            'schema must be a scalar structure.');
    end

    mode = validatestring(mode,{'meteo','Q','both'});
    schema.mode = mode;
    schema = local_default(schema,'name','generic dataset');
    schema = local_default(schema,'format','csv');
    schema = local_default(schema,'layout','one_file_per_basin');
    schema = local_default(schema,'missing_values',[-9999 -999]);
    schema = local_default(schema,'strict_time',true);
    schema = local_default(schema,'cache',struct());
    schema.cache = local_default(schema.cache, ...
    'enabled',true);
    schema.cache = local_default(schema.cache, ...
    'directory','_cache_generic_canonical');
    schema.cache = local_default(schema.cache, ...
    'canonical_enabled',true);
    schema.cache = local_default(schema.cache, ...
    'native_enabled',false);
    schema.cache = local_default(schema.cache, ...
    'aggregate_enabled',false);
    schema.cache = local_default(schema.cache, ...
    'aggregate_min_basins',50);

    schema.format = validatestring(lower(string(schema.format)), ...
        {'csv','netcdf'});
    schema.layout = validatestring(lower(string(schema.layout)), ...
        {'one_file_per_basin','multi_file_per_basin','wide_files'});

    if strcmp(schema.layout,'one_file_per_basin')
    
        local_require(schema,'file');
        local_require(schema.file,'pattern');
    
        schema.file = local_default(schema.file, ...
        'delimiter',',');
        schema.file = local_default(schema.file, ...
        'header_lines',0);
    
        if ~contains(string(schema.file.pattern),'{gauge}')
            error('validate_hydro_schema:MissingGaugeToken', ...
                'schema.file.pattern must contain the token {gauge}.');
        end
    
    elseif strcmp(schema.layout,'multi_file_per_basin')
        local_require(schema,'files');
        if ~isstruct(schema.files) ...
                || ~isscalar(schema.files) ...
                || isempty(fieldnames(schema.files))
            error('validate_hydro_schema:BadNamedFiles', ...
                'schema.files must be a nonempty scalar structure.');
        end
        names = fieldnames(schema.files);
        for iFile = 1:numel(names)
            name = names{iFile};
            fileSpec = schema.files.(name);
            if ~isstruct(fileSpec) ...
                    || ~isscalar(fileSpec)
                error('validate_hydro_schema:BadNamedFile', ...
                    'schema.files.%s must be a scalar structure.',name);
            end
            local_require(fileSpec,'pattern');
            if ~contains(string(fileSpec.pattern),'{gauge}')
                error('validate_hydro_schema:MissingGaugeToken', ...
                    'schema.files.%s.pattern must contain {gauge}.',name);
            end
            fileSpec = local_default(fileSpec,'format',schema.format);
            fileSpec.format = validatestring(lower(string(fileSpec.format)), ...
                {'csv','netcdf'});
            fileSpec = local_default(fileSpec,'delimiter',',');
            fileSpec = local_default(fileSpec,'header_lines',0);
            fileSpec = local_default(fileSpec,'time',struct());
            schema.files.(name) = fileSpec;
        end
        schema = local_default(schema,'file',struct());

    else
        % Wide-file datasets: individual variables declare their own files.
        schema = local_default(schema,'file',struct());
        schema.file = local_default(schema.file,'delimiter',',');
        schema.file = local_default(schema.file,'header_lines',0);
    end

    schema = local_default(schema,'id',struct());
    schema.id = local_default(schema.id,'pad_width',0);
    schema.id = local_default(schema.id,'prefix','');
    schema.id = local_default(schema.id,'strip_decimal',true);
    schema.id = local_default(schema.id,'strip_leading_zeros',false);

    local_require(schema,'timeline');
    local_require(schema.timeline,'reference');
    schema.timeline = local_default(schema.timeline,'step',days(1));
    if ~isdatetime(schema.timeline.reference) ...
            || ~isscalar(schema.timeline.reference)
        error('validate_hydro_schema:BadReference', ...
            'schema.timeline.reference must be a scalar datetime.');
    end
    if ~isduration(schema.timeline.step) ...
            || ~isscalar(schema.timeline.step) ...
            || seconds(schema.timeline.step) <= 0
        error('validate_hydro_schema:BadStep', ...
            'schema.timeline.step must be a positive scalar duration.');
    end

    schema = local_default(schema,'time',struct());
    schema.time = local_default(schema.time,'snap','');
    schema.time = local_default(schema.time,'join','exact');
    schema.time.join = validatestring(lower(string(schema.time.join)), ...
        {'exact','calendar_day'});
    if strcmp(schema.format,'csv')
    
        schema.time = local_default(schema.time,'mode','column');
    
        schema.time.mode = validatestring( ...
            lower(string(schema.time.mode)), ...
            {'column','ymd_columns','ymd_row_sequence'});
    
        switch schema.time.mode
    
            case 'column'
                local_require(schema.time,'column');
                schema.time = local_default( ...
                    schema.time,'input_format','');
                schema.time = local_default( ...
                    schema.time,'timezone','');
    
            case {'ymd_columns','ymd_row_sequence'}
                local_require(schema.time,'year_column');
                local_require(schema.time,'month_column');
                local_require(schema.time,'day_column');
                schema.time = local_default(schema.time,'hour_column','');
                schema.time = local_default(schema.time,'minute_column','');
                schema.time = local_default(schema.time,'second_column','');
                if strcmp(schema.time.mode,'ymd_row_sequence')
                    schema.time = local_default(schema.time, ...
                        'steps_per_day',round(days(1)/schema.timeline.step));
                    schema.time = local_default(schema.time, ...
                        'partial_first','final_slots');
                    if ~isscalar(schema.time.steps_per_day) ...
                            || schema.time.steps_per_day < 1 ...
                            || mod(schema.time.steps_per_day,1) ~= 0
                        error('validate_hydro_schema:BadStepsPerDay', ...
                            'time.steps_per_day must be a positive integer.');
                    end
                    schema.time.partial_first = validatestring( ...
                        lower(string(schema.time.partial_first)), ...
                        {'final_slots','first_slots'});
                end
        end
    
    else
        schema.time = local_default(schema.time,'mode','index');
    
        schema.time.mode = validatestring( ...
            lower(string(schema.time.mode)), ...
            {'index','variable'});
    
        if strcmp(schema.time.mode,'index')
            local_require(schema.time,'file_start');
    
            if ~isdatetime(schema.time.file_start) ...
                    || ~isscalar(schema.time.file_start)
                error('validate_hydro_schema:BadFileStart', ...
                    ['schema.time.file_start must be ' ...
                     'a scalar datetime.']);
            end
        else
            local_require(schema.time,'column');
            schema.time = local_default( ...
                schema.time,'units','datetime');
            schema.time = local_default( ...
                schema.time,'origin',NaT);
        end
    end

    if strcmp(schema.layout,'multi_file_per_basin')
        names = fieldnames(schema.files);
        for iFile = 1:numel(names)
            name = names{iFile};
            fileSpec = schema.files.(name);
            timeSpec = local_overlay(schema.time,fileSpec.time);
            local_validate_time_spec(timeSpec,fileSpec.format,name);
            fileSpec.time = timeSpec;
            schema.files.(name) = fileSpec;
        end
    end

    local_require(schema,'variables');
    required = strings(0,1);
    if any(strcmp(mode,{'meteo','both'}))
        required = [required;"P";"Ep";"T"];
    end
    if strcmp(mode,'Q')
        required = [required;"Q"];
    end
    for i = 1:numel(required)
        field = char(required(i));
        if ~isfield(schema.variables,field)
            error('validate_hydro_schema:MissingVariable', ...
                'schema.variables.%s is required in %s mode.',field,mode);
        end
    end
    fields = fieldnames(schema.variables);
    
    for i = 1:numel(fields)
    
        field = fields{i};
        spec = schema.variables.(field);
    
        if ~isstruct(spec) ...
                || ~isscalar(spec)
            error('validate_hydro_schema:BadVariable', ...
                'schema.variables.%s must be a scalar structure.',field);
        end
    
        hasScalarSource = isfield(spec,'source') ...
            && ~isempty(spec.source);
        hasSourceSet = isfield(spec,'sources') ...
            && ~isempty(spec.sources);
        if hasScalarSource && hasSourceSet
            error('validate_hydro_schema:ConflictingSources', ...
                ['schema.variables.%s must define .source or .sources, ' ...
                 'not both.'],field);
        end
        hasSource = hasScalarSource || hasSourceSet;

        hasDerive = isfield(spec,'derive') ...
            && ~isempty(spec.derive);

        if any(strcmp(schema.layout, ...
                {'one_file_per_basin','multi_file_per_basin'}))
    
            if ~hasSource ...
                    && ~hasDerive
                error('validate_hydro_schema:MissingVariableSource', ...
                    ['schema.variables.%s must define either ' ...
                     '.source or .derive.'],field);
            end
            if strcmp(schema.layout,'multi_file_per_basin') && hasSource
                if ~isfield(spec,'file') || isempty(spec.file)
                    error('validate_hydro_schema:MissingVariableFile', ...
                        ['schema.variables.%s must name a schema.files ' ...
                         'entry in .file.'],field);
                end
                fileName = char(string(spec.file));
                if ~isfield(schema.files,fileName)
                    error('validate_hydro_schema:UnknownVariableFile', ...
                        'schema.variables.%s references unknown file %s.', ...
                        field,fileName);
                end
            end
    
        else
            % Wide-file variables can either be read from one of several
            % files or derived from other variables.
            hasFiles = isfield(spec,'files') ...
                && ~isempty(spec.files);
            hasDerive = isfield(spec,'derive') ...
                && ~isempty(spec.derive);
    
            if ~hasFiles ...
                    && ~hasDerive
                error('validate_hydro_schema:MissingWideSource', ...
                    ['schema.variables.%s must define either ' ...
                     '.files or .derive for wide_files layout.'],field);
            end
    
            if hasFiles
                if ~iscell(spec.files) ...
                        && ~isstring(spec.files)
                    error('validate_hydro_schema:BadFiles', ...
                        'schema.variables.%s.files must be a cell/string array.', ...
                        field);
                end
    
                spec = local_default(spec,'selector','');
                spec = local_default(spec,'default',1);
            end
        end
    
        spec = local_default(spec,'units','');
        spec = local_default(spec,'target_units','');
        spec = local_default(spec,'scale',1);
        spec = local_default(spec,'offset',0);
        spec = local_default(spec,'missing_values',schema.missing_values);
        spec = local_default(spec,'valid_min',-Inf);
        spec = local_default(spec,'valid_max',Inf);
        spec = local_default(spec,'clip_min',-Inf);
        spec = local_default(spec,'clip_max',Inf);
        spec = local_default(spec,'area_normalize',false);
        spec = local_default(spec,'transform',[]);
        spec = local_default(spec,'invalid_le',[]);
        spec = local_default(spec,'fill_missing',[]);
        spec = local_default(spec,'quality',[]);
        if ~isempty(spec.quality)
            if ~isstruct(spec.quality) ...
                    || ~isscalar(spec.quality)
                error('validate_hydro_schema:BadQuality', ...
                    'schema.variables.%s.quality must be a scalar structure.',field);
            end
            local_require(spec.quality,'sources');
            qualitySources = string(spec.quality.sources(:));
            if isempty(qualitySources) ...
                    || any(strlength(qualitySources) == 0)
                error('validate_hydro_schema:BadQualitySources', ...
                    'schema.variables.%s.quality.sources must be nonempty.',field);
            end
            spec.quality.sources = cellstr(qualitySources);
            spec.quality = local_default(spec.quality,'valid_min',-Inf);
            spec.quality = local_default(spec.quality,'valid_max',Inf);
            spec.quality = local_default(spec.quality,'require_finite',true);
            if ~isscalar(spec.quality.valid_min) ...
                    || ~isscalar(spec.quality.valid_max) ...
                    || spec.quality.valid_min > spec.quality.valid_max
                error('validate_hydro_schema:BadQualityBounds', ...
                    'Quality bounds for %s must be ordered scalars.',field);
            end
        end
        if hasSourceSet
            if ~(iscellstr(spec.sources) ...
                    || isstring(spec.sources)) ...
                    || isempty(spec.sources)
                error('validate_hydro_schema:BadSources', ...
                    'schema.variables.%s.sources must be a string/cellstr list.', ...
                    field);
            end
            sourceNames = string(spec.sources(:));
            if any(strlength(sourceNames) == 0) ...
                    || numel(unique(sourceNames,'stable')) ~= numel(sourceNames)
                error('validate_hydro_schema:BadSources', ...
                    ['schema.variables.%s.sources must contain unique, ' ...
                     'nonempty names.'],field);
            end
            spec.sources = cellstr(sourceNames);
            spec = local_default(spec,'source_operation','aliases');
            spec.source_operation = validatestring( ...
                lower(string(spec.source_operation)), ...
                {'aliases','coalesce','mean'});
        end
    
        schema.variables.(field) = spec;
    end

    schema = local_default(schema,'file',struct());
    schema.file = local_default(schema.file,'delimiter',',');
    schema.file = local_default(schema.file,'header_lines',0);

    schema = local_default(schema,'aux',struct());
    schema.aux = local_default(schema.aux,'provider',[]);
    schema.aux = local_default(schema.aux,'values',[]);
    schema.aux = local_default(schema.aux,'tables',[]);
    schema.aux = local_default(schema.aux,'series',[]);
    hasProvider=~isempty(schema.aux.provider);
    hasValues=~isempty(schema.aux.values);
    hasDeclarative=~isempty(schema.aux.tables) || ~isempty(schema.aux.series);
    if nnz([hasProvider,hasValues,hasDeclarative]) > 1
        error('validate_hydro_schema:ConflictingAux', ...
            ['Configure provider, values, or declarative tables/series; ' ...
             'tables and series may be combined.']);
    end
    schema = local_default(schema,'progress',struct());
    schema.progress = local_default(schema.progress,'label', ...
        sprintf('... Reading %s data',char(string(schema.name))));
end

function S = local_default(S,field,value)
    if ~isfield(S,field) || isempty(S.(field))
        S.(field) = value;
    end
end

function local_require(S,field)
    if ~isfield(S,field) || isempty(S.(field))
        error('validate_hydro_schema:MissingField', ...
            'Required schema field is missing: %s.',field);
    end
end

function out = local_overlay(base,override)
    out = base;
    fields = fieldnames(override);
    for i = 1:numel(fields)
        out.(fields{i}) = override.(fields{i});
    end
end

function local_validate_time_spec(timeSpec,format,fileName)
    if strcmp(format,'csv')
        if ~isfield(timeSpec,'mode') || isempty(timeSpec.mode)
            error('validate_hydro_schema:MissingFileTime', ...
                'Named file %s requires a CSV time mode.',fileName);
        end
        mode = validatestring(lower(string(timeSpec.mode)), ...
            {'column','ymd_columns','ymd_row_sequence'});
        if strcmp(mode,'column')
            local_require(timeSpec,'column');
        else
            local_require(timeSpec,'year_column');
            local_require(timeSpec,'month_column');
            local_require(timeSpec,'day_column');
        end
    else
        if ~isfield(timeSpec,'mode') || isempty(timeSpec.mode)
            error('validate_hydro_schema:MissingFileTime', ...
                'Named file %s requires a NetCDF time mode.',fileName);
        end
        mode = validatestring(lower(string(timeSpec.mode)), ...
            {'index','variable'});
        if strcmp(mode,'index')
            local_require(timeSpec,'file_start');
        else
            local_require(timeSpec,'column');
        end
    end
end
