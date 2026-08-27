function [p,mdl] = setup_parpool(mdl,bas,frestart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SETUP_PARPOOL  Initialize (and optionally restart) parallel pool quietly.
%
% SYNOPSIS: [p,mdl] = setup_parpool(mdl,bas,frestart)
%   mdl             structure with model settings
%    .calc           'seq', 'par', or 'parfeval'
%                    'par'      uses parfor
%                    'parfeval' uses asynchronous batched parfeval
%    .names          model names, e.g. ["hbv","xinanjiang",...]
%   bas             structure with basin info
%    .K              # basins (train + eval)
%   frestart        OPTIONAL: 0/1 close any existing pool and reopen
%   p               OUTPUT: parallel pool object (or [] if seq/no toolbox)
%   mdl             OUTPUT: may be updated (fallback to 'seq' if needed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 ...
        || isempty(frestart)
    frestart = 0;
end
p = [];

% ----------------
% Basic validation
% ----------------
if ~isstruct(mdl) ...
        || ~isfield(mdl,'calc') ...
        || isempty(mdl.calc)
    mdl.calc = 'seq';
    return
end

calc = lower(string(mdl.calc));
if ~ismember(calc,["par","parfeval"])
    mdl.calc = 'seq';
    return
end
mdl.calc = char(calc);

if ~isstruct(bas) ...
        || ~isfield(bas,'K') ...
        || isempty(bas.K)
    K = 1;
else
    K = double(bas.K);
    if ~isfinite(K) || K < 1
        K = 1;
    end
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

if ~hasPCT ...
        || exist('parpool','file') ~= 2 ...
        || exist('gcp','file') ~= 2
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

if ~isempty(p) && frestart
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
if isfield(mdl,'names') ...
        && ~isempty(mdl.names)
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
nw = max(1,min(double(nw),K)); %#ok

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
