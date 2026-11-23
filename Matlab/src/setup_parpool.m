function [p,mdl] = setup_parpool(mdl,bas,frestart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SETUP_PARPOOL Prepare basin-level parallel execution.
%
%  Reuse or restart a MATLAB pool when requested, choosing no more workers
%  than basins or available cores.
%
% SYNOPSIS:
%   [p,mdl] = setup_parpool(mdl,bas,frestart)
%   [p,mdl] = setup_parpool(mdl,bas)
%
% INPUT ARGUMENTS:
%   mdl             model execution settings
%    .calc           'seq', 'par', or 'parfeval'
%    .names          optional model names for unloading client MEX files
%    .numWorkers     optional positive integer worker count, capped by K
%                    An existing pool is resized only when this is set.
%   bas             selected-basin information
%    .K              number of selected basins
%   frestart        optional restart flag; default false
%
% OUTPUT ARGUMENTS:
%   p               pool object, or empty if sequential/unavailable
%   mdl             updated settings; .calc may fall back to 'seq'
%
% NOTES:
%   If Parallel Computing Toolbox is unavailable, execution falls back to
%   sequential mode.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Â© Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3 || ...
        isempty(frestart)
        frestart = 0;
    end
    p = [];

    % ----------------
    % Basic validation
    % ----------------
    if ~isstruct(mdl) || ...
        ~isfield(mdl,'calc') || ...
        isempty(mdl.calc)
        mdl.calc = 'seq';
        return
    end

    calc = lower(string(mdl.calc));
    if ~ismember(calc,["par","parfeval"])
        mdl.calc = 'seq';
        return
    end
    mdl.calc = char(calc);

    if ~isstruct(bas) || ...
        ~isfield(bas,'K') || ...
        isempty(bas.K)
        K = 1;
    else
        K = double(bas.K);
        if ~isfinite(K) || ...
            K < 1
            K = 1;
        end
    end

    % An explicit worker count is opt-in. Existing callers retain their
    % default pool reuse and automatic core-count behavior.
    requested = [];
    if isfield(mdl,'numWorkers') && ...
        ~isempty(mdl.numWorkers)
        requested = mdl.numWorkers;
        if ~isnumeric(requested) || ...
            ~isscalar(requested) || ...
            ~isreal(requested) || ...
            ~isfinite(requested) || ...
            requested < 1 || ...
            requested ~= fix(requested)
            error('setup_parpool:InvalidWorkers', ...
                'mdl.numWorkers must be a positive integer.');
        end
        requested = max(1,min(double(requested),K));
    end

    % ------------------------------------------
    % Check Parallel Computing Toolbox available
    % ------------------------------------------
    try
        hasPCT = license('test', ...
            'Distrib_Computing_Toolbox');
    catch
        hasPCT = false;
    end

    if ~hasPCT || ...
        exist('parpool','file') ~= 2 || ...
        exist('gcp','file') ~= 2
        mdl.calc = 'seq';
        warning('setup_parpool:noPCT', ...
            ['      Warning: setup_parpool: ' ...
            'Parallel Computing Toolbox not available ' ...
            'or incomplete -> forced to sequential mode']);
        return
    end

    % -------------------
    % Existing pool logic
    % -------------------
    try
        p = gcp('nocreate');
    catch
        p = [];
    end

    resize = ~isempty(p) && ...
        ~isempty(requested) && ...
        p.NumWorkers ~= requested;
    if ~isempty(p) && ...
        (frestart || ...
        resize)
        try
            delete(p);
        catch
        end
        p = [];
    end

    if ~isempty(p)
        return
    end

    % ------------------------------------
    % Unload model MEX functions on client
    % ------------------------------------
    if isfield(mdl,'names') && ...
        ~isempty(mdl.names)
        names = mdl.names;
        if isstring(names)
            names = cellstr(names);
        elseif ischar(names)
            names = {names};
        end

        patterns = {'crr_%s','crr_%s_mex'};

        for ii = 1:numel(names)
            base = strtrim(names{ii});
            if isempty(base)
                continue
            end

            for pp = 1:numel(patterns)
                fname = sprintf(patterns{pp},base);
                [~,fname] = fileparts(fname);

                try
                    clear(fname);
                catch
                end
            end
        end

        try
            rehash toolboxcache
        catch
        end

        disp(['      setup_parpool: ' ...
            'unloaded client crr_* MEX functions']);
    end

    % ----------------------
    % Determine worker count
    % ----------------------
    try
        nw = feature('numcores');
    catch
        nw = 1;
    end
    nw = max(1,min(double(nw),K));      %#ok
    if ~isempty(requested)
        % Used below inside evalc when opening the pool.
        nw = requested;                 %#ok<NASGU>
    end

    % ---------------------
    % Open new pool quietly
    % ---------------------
    fprintf('... Opening parallel pool for %s ... ',mdl.calc);

    try
        try
            evalc('p = parpool(''Processes'', nw);');
        catch
            evalc('p = parpool(nw);');
        end

        fprintf('Done (workers = %d)\n',p.NumWorkers);

    catch ME
        p = [];
        mdl.calc = 'seq';
        fprintf('Failed -> seq\n');
        warning('setup_parpool:parpoolFailed', ...
            '      Warning: setup_parpool: %s',ME.message);
    end

end
