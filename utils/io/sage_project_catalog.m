function out = sage_project_catalog(action,varargin)
%SAGE_PROJECT_CATALOG Discover and load folder-backed SAGE projects.

% Each immediate child of <SAGEhydrology>/projects may contain:
%   project.json  required descriptive metadata
%   setup.json    optional GUI configuration overrides

% The JSON format deliberately contains data only.  Project discovery never
% executes MATLAB code, which makes locally added projects predictable and
% safe to inspect from the GUI.

    action = lower(char(string(action)));
    switch action
        case 'discover'
            out = local_discover(varargin{1});
        case 'configuration'
            out = local_configuration(varargin{1},varargin{2});
        otherwise
            error('SAGE:Projects:UnknownAction', ...
                'Unknown project-catalog action: %s',action);
    end
end

function projects = local_discover(projectRoot)
    projectRoot = char(string(projectRoot));
    projects = repmat(local_empty_project(),0,1);
    if ~isfolder(projectRoot)
        return
    end

    folders = dir(projectRoot);
    folders = folders([folders.isdir]);
    names = string({folders.name});
    folders = folders(~ismember(names,[".",".."])) ;

    for ii = 1:numel(folders)
        folder = fullfile(projectRoot,folders(ii).name);
        descriptor = fullfile(folder,'project.json');
        if ~isfile(descriptor)
            continue
        end
        try
            meta = jsondecode(fileread(descriptor));
            p = local_empty_project();
            p.folder = folder;
            p.descriptor = descriptor;
            p.id = local_text(meta,'id',folders(ii).name);
            p.name = local_text(meta,'name',p.id);
            p.kind = lower(local_text(meta,'kind','sage'));
            p.status = lower(local_text(meta,'status','draft'));
            p.summary = local_text(meta,'summary','');
            p.paper_title = local_text(meta,'paper_title','');
            p.authors = local_text(meta,'authors','');
            p.abstract = local_text(meta,'abstract','');
            p.reference = local_text(meta,'reference','');
            p.paper_url = local_text(meta,'paper_url','');
            p.notes = local_lines(meta,'notes');
            p.setup_file = local_resolve(folder, ...
                local_text(meta,'setup_file','setup.json'));
            p.source_export = local_resolve(folder, ...
                local_text(meta,'source_export',''));
            p.entry_point = local_resolve(folder, ...
                local_text(meta,'entry_point',''));
            p.ready = strcmp(p.status,'ready') && ...
                ((strcmp(p.kind,'sage') && isfile(p.setup_file)) || ...
                 (strcmp(p.kind,'standalone') && isfile(p.entry_point)));
            projects(end+1,1) = p; %#ok<AGROW>
        catch ME
            warning('SAGE:Projects:InvalidDescriptor', ...
                'Could not read %s: %s',descriptor,ME.message);
        end
    end

    if ~isempty(projects)
        [~,order] = sort(lower(string({projects.id})));
        projects = projects(order);
    end
end

function cfg = local_configuration(project,baseCfg)
    if ~strcmp(project.kind,'sage')
        error('SAGE:Projects:NotSageProject', ...
            '%s is a standalone project.',project.name);
    end
    if ~isfile(project.setup_file)
        error('SAGE:Projects:MissingSetup', ...
            'Project setup is not available: %s',project.setup_file);
    end
    overrides = jsondecode(fileread(project.setup_file));
    cfg = local_merge(baseCfg,overrides);

    % Preserve installation-specific locations.  Project files describe a
    % reproducible experiment, not the original author's computer.
    if isfield(baseCfg,'root'), cfg.root = baseCfg.root; end
    if isfield(baseCfg,'dirDroot'), cfg.dirDroot = baseCfg.dirDroot; end
    if isfield(baseCfg,'dirres'), cfg.dirres = baseCfg.dirres; end
    if isfield(baseCfg,'SAGEhydro'), cfg.SAGEhydro = baseCfg.SAGEhydro; end
    if isfield(baseCfg,'SAGEdir'), cfg.SAGEdir = baseCfg.SAGEdir; end
    if isfield(cfg,'dirDroot') && isfield(cfg,'region')
        cfg.dirD = fullfile(char(cfg.dirDroot),char(cfg.region));
    end
end

function target = local_merge(target,source)
    fields = fieldnames(source);
    for ii = 1:numel(fields)
        name = fields{ii};
        value = source.(name);
        if isstruct(value) && isscalar(value) && ...
                isfield(target,name) && isstruct(target.(name)) && ...
                isscalar(target.(name))
            target.(name) = local_merge(target.(name),value);
        else
            target.(name) = value;
        end
    end
end

function value = local_text(S,name,fallback)
    value = fallback;
    if isfield(S,name) && ~isempty(S.(name))
        value = char(string(S.(name)));
    end
end

function lines = local_lines(S,name)
    lines = strings(0,1);
    if isfield(S,name) && ~isempty(S.(name))
        lines = string(S.(name));
        lines = lines(:);
    end
end

function path = local_resolve(folder,value)
    value = char(string(value));
    if isempty(value)
        path = '';
    elseif isfile(value) || isfolder(value)
        path = value;
    else
        path = char(java.io.File(fullfile(folder,value)).getCanonicalPath());
    end
end

function p = local_empty_project()
    p = struct('id','','name','','kind','sage','status','draft', ...
        'summary','','paper_title','','authors','','abstract','', ...
        'reference','','paper_url','','notes',strings(0,1), ...
        'folder','','descriptor','','setup_file','','source_export','', ...
        'entry_point','','ready',false);
end
