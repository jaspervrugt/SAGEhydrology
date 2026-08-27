function sageRoot = ensureSAGEpath(rootFolder)
%ENSURESAGEPATH Add SAGEhydrology source folders to the MATLAB path.
%
%   sageRoot = ensureSAGEpath(rootFolder)
%
% Source MATLAB only. In deployed mode the MATLAB Compiler controls the
% runtime search path and this function returns immediately.
%
% The project is located first from this file itself:
%
%   SAGEhydrology/utils/ensureSAGEpath.m
%
% so a fresh MATLAB session does not depend on pwd, startup folders, or a
% correctly initialized root variable. rootFolder is retained as a fallback
% for exported scripts and alternative installations.

    sageRoot = '';

    if isdeployed
        return
    end

    if nargin < 1 ...
            || isempty(rootFolder)
        rootFolder = '';
    elseif ~(ischar(rootFolder) ...
            || (isstring(rootFolder) ...
            && isscalar(rootFolder)))
        error('ensureSAGEpath:InvalidInput', ...
            ['rootFolder must be empty, a character vector, ' ...
             'or a string scalar.']);
    end

    requiredFunctions = { ...
        'region_helpers', ...
        'attr_catalog', ...
        'read_attr', ...
        'read_meteo', ...
        'read_Q', ...
        'read_attribute_data', ...
        'read_meteo_data', ...
        'read_discharge_data', ...
        'attribute_schema', ...
        'meteo_schema', ...
        'discharge_schema', ...
        'read_model', ...
        'hymod', ...
        'hmodel', ...
        'sacsma', ...
        'xinanjiang', ...
        'gr4jA', ...
        'hbv', ...
        'cfe_nwm', ...
        'print_SAGE', ...
        'apply_parameter_override', ...
        'split_basin_ids'};

    rootFolder = char(string(rootFolder));

    % ---------------------------------------------------------------
    % 1. Preferred source: locate SAGEhydrology from this function.
    %
    % ensureSAGEpath.m lives in:
    %     <SAGEhydrology>/utils/ensureSAGEpath.m
    %
    % This makes startup independent of MATLAB's current directory.
    % ---------------------------------------------------------------
    thisFile = mfilename('fullpath');
    selfRoot = '';

    if ~isempty(thisFile)
        utilsFolder = fileparts(thisFile);
        candidate = fileparts(utilsFolder);

        if local_is_sage_root(candidate)
            selfRoot = candidate;
        end
    end

    % ---------------------------------------------------------------
    % 2. Fallback candidates from the caller-supplied root.
    % ---------------------------------------------------------------
    candidates = {};

    if ~isempty(rootFolder)
        candidates = { ...
            fullfile(rootFolder,'SAGEhydrology'), ...
            fullfile(rootFolder,'Software','SAGEhydrology'), ...
            rootFolder};
    end

    if ~isempty(selfRoot)
        % Prefer the known location of the function actually being run.
        candidates = [{selfRoot}, candidates];
    end

    % Remove duplicate candidate paths while preserving order.
    if ~isempty(candidates)
        keys = cellfun(@local_path_key,candidates, ...
            'UniformOutput',false);
        [~,ia] = unique(keys,'stable');
        candidates = candidates(sort(ia));
    end

    for ii = 1:numel(candidates)
        candidate = candidates{ii};

        if local_is_sage_root(candidate)
            sageRoot = candidate;
            break
        end
    end

    if isempty(sageRoot)
        if isempty(rootFolder)
            rootText = '<empty>';
        else
            rootText = rootFolder;
        end

        error('ensureSAGEpath:MissingRoot', ...
            ['Could not locate the SAGEhydrology project. ' ...
             'Caller root folder: %s\n' ...
             'ensureSAGEpath location: %s'], ...
            rootText,thisFile);
    end

    % ---------------------------------------------------------------
    % Add only code trees. genpath includes nested folders such as
    % utils/info, gui/diagnostics, model subfolders, and region readers.
    % ---------------------------------------------------------------
    codeRoots = { ...
        fullfile(sageRoot,'gui'), ...
        fullfile(sageRoot,'src'), ...
        fullfile(sageRoot,'utils'), ...
        fullfile(sageRoot,'models'), ...
        fullfile(sageRoot,'regions')};

    for ii = 1:numel(codeRoots)

        codeRoot = codeRoots{ii};

        if ~isfolder(codeRoot)
            warning('ensureSAGEpath:MissingFolder', ...
                'Expected code folder not found: %s', ...
                codeRoot);
            continue
        end

        addpath(genpath(codeRoot));
    end

    % ---------------------------------------------------------------
    % Validate functions required during startup / basin preparation.
    % ---------------------------------------------------------------
    missing = requiredFunctions( ...
        cellfun(@(name) isempty(which(name)),requiredFunctions));

    if ~isempty(missing)
        error('ensureSAGEpath:MissingFunctions', ...
            ['The following required functions ' ...
            'remain unavailable: %s'], ...
            strjoin(missing,', '));
    end
end


function tf = local_is_sage_root(candidate)
%LOCAL_IS_SAGE_ROOT True for a recognizable SAGEhydrology project root.

    tf = false;

    if isempty(candidate) ...
            || ~isfolder(candidate)
        return
    end

    tf = isfolder(fullfile(candidate,'gui')) ...
        && isfolder(fullfile(candidate,'src')) ...
        && isfolder(fullfile(candidate,'utils')) ...
        && isfolder(fullfile(candidate,'models')) ...
        && isfolder(fullfile(candidate,'regions'));
end


function key = local_path_key(p)
%LOCAL_PATH_KEY Normalize a path for duplicate detection.

    p = char(string(p));

    try
        p = char(java.io.File(p).getCanonicalPath());
    catch
    end

    if ispc
        p = lower(p);
    end

    key = p;
end
