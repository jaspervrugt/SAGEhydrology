function [mdl,d] = read_user_model_info(mdl,prd,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_MODEL Reads parameter information of model
% SYNOPSIS: [mdl,d] = read_user_model_info(mdl,prd,verbose)
%   mdl         structure with crr_model to use and settings
%    .model      choice of model
%                 1 hymod
%                 2 hmodel
%                 3 sacsma
%                 4 Xinanjiang
%                 5 gr4j
%                 6 hbv
%                 7 user-defined model
%    .mcode      scalar with numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation sacsma_ode C++
%    .calc       model execution
%      'seq'      sequential execution of watersheds
%      'par'      parallel execution of watersheds
%    .mode       assessment design
%                 1 = training basins | training period
%                 2 = training basins | evaluation period/mask
%                 3 = training + evaluation basins | training period
%                 4 = training + evaluation basins | evaluation period/mask
%    .names      list of model names
%    .eval_mode  evaluation design used during SAGE training
%      'te'       training basins evaluated on evaluation period
%      'et'       evaluation basins evaluated on training period
%      'ee'       evaluation basins evaluated on evaluation period
%   prd         structure with temporal resolution in field dt
%   verbose     OPTIONAL: print to screen
%   mdl         revised structure with crr_model to use and settings
%   d           number of model parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, May 2026                                  %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    verbose = 0;
end

if nargin < 2 ...
        || isempty(prd) ...
        || ~isfield(prd,'dt') ...
        || isempty(prd.dt)
    dt = 1;
else
    dt = double(prd.dt);
end

if ~ismember(dt,[1 24 96])
    error('read_user_model_info:InvalidResolution', ...
        'Unsupported temporal resolution dt = %g.',dt);
end

if isfield(mdl,'root')
    udir = user_model_dir_from_root(mdl.root);
else
    udir = user_model_dir_from_root(pwd);
end

% The installed manifest is authoritative for the authored model name.
% Refresh it whenever metadata is loaded so a long-lived GUI session cannot
% run a replacement AI model under the previous candidate's name.
manifestFile = fullfile(udir,'ai_model_manifest.json');
mdl.ai_requires_kosugi = false;
mdl.ai_requires_dual_kosugi = false;
if isfile(manifestFile)
    try
        manifest = jsondecode(fileread(manifestFile));
        authoredName = string(manifest.model_name);
        if isscalar(authoredName) && ~isempty(regexp(char(authoredName), ...
                '^[A-Za-z][A-Za-z0-9_]*$', 'once'))
            mdl.name = char(authoredName);
            if isfield(mdl,'names') && numel(mdl.names) >= 9
                if iscell(mdl.names)
                    mdl.names{9} = char(authoredName);
                else
                    mdl.names(9) = authoredName;
                end
            end
        end
    catch
        % Existing validation reports a malformed manifest elsewhere. Keep
        % metadata loading backward compatible when no authored name exists.
    end
end
specFile = fullfile(udir,'model_specification.json');
if isfile(specFile)
    try
        specification = jsondecode(fileread(specFile));
        forcingNames = string({specification.forcing_variables.name});
        mdl.ai_requires_kosugi = all(ismember( ...
            ["kosugi_a","kosugi_b","kosugi_c","kosugi_p0"], ...
            forcingNames));
    catch
        mdl.ai_requires_kosugi = false;
    end
end
prepareFile = fullfile(udir,'user_model_prepare.cpp');
if isfile(prepareFile)
    try
        prepareSource = fileread(prepareFile);
        mdl.ai_requires_dual_kosugi = contains(prepareSource, ...
            'data.hydro.fdc.kosugi_dual');
    catch
        mdl.ai_requires_dual_kosugi = false;
    end
end
infoFile = fullfile(udir, ...
    'user_model_info.mat');

if ~isfile(infoFile)
    error('read_user_model_info:MissingInfoFile', ...
        ['Cannot find user_model_info.mat ' ...
        'in user_model folder.']);
end

U = load(infoFile);
required = {'par_info','m','y0','pspace'};
for k = 1:numel(required)
    if ~isfield(U,required{k})
        error('read_user_model_info:MissingMetadata', ...
            'user_model_info.mat is missing %s.',required{k});
    end
end

par_info = U.par_info;

mdl.par_names = par_info(:,2)';
mdl.par_desc = par_info(:,3)';
par_units = par_info(:,4)';
th_min = cell2mat(par_info(:,5));
th_max = cell2mat(par_info(:,6));

mdl.th_min_daily = th_min;
mdl.th_max_daily = th_max;
mdl.par_units_daily = par_units;

[th_min,th_max,par_units] = local_convert_resolution( ...
    th_min,th_max,par_units,dt);

mdl.par_units = par_units;
mdl.th_min = th_min;
mdl.th_max = th_max;

mdl.m = U.m;
mdl.y0 = U.y0;
mdl.pspace = U.pspace;
if dt == 1
    mdl.L = 20;             % daily routing-memory length
else
    mdl.L = 100;            % cap subdaily routing state space
end
if isfield(U,'state_name')
    mdl.state_name = string(U.state_name(:));
    if numel(mdl.state_name) == 3
        mdl.state_name = [mdl.state_name; ...
            compose("unit_hydrograph_1_storage_%d",1:mdl.L)'; ...
            compose("unit_hydrograph_2_storage_%d",1:mdl.L)'];
    end
end
if isfield(U,'state_units')
    mdl.state_units = string(U.state_units(:));
    if numel(mdl.state_units) == 3
        mdl.state_units = [mdl.state_units; ...
            repmat("mm",2*mdl.L,1)];
    end
end
if isfield(U,'swe_ind')
    mdl.swe_ind = double(U.swe_ind(:)');
end
if isfield(U,'sm_ind')
    mdl.sm_ind = double(U.sm_ind(:)');
end

d = numel(mdl.par_names);

if verbose
    fprintf(['Loaded AI-assisted model metadata ' ...
        'with %d parameters.\n'],d);
end
end

function [th_min,th_max,par_units] = ...
    local_convert_resolution(th_min,th_max,par_units,dt)
%LOCAL_CONVERT_RESOLUTION Convert daily parameter definitions.

if dt == 1
    return
end

for j = 1:numel(par_units)
    unit = char(string(par_units{j}));

    if strcmp(unit,'1/d')
        th_min(j) = th_min(j)/dt;
        th_max(j) = th_max(j)/dt;
        par_units{j} = local_rate_unit(dt);

    elseif strcmp(unit,'mm/d')
        th_min(j) = th_min(j)/dt;
        th_max(j) = th_max(j)/dt;
        par_units{j} = local_flux_unit(dt);

    elseif strcmp(unit,'d')
        th_min(j) = th_min(j)*dt;
        th_max(j) = th_max(j)*dt;
        par_units{j} = local_time_unit(dt);

    elseif ~isempty(regexp(unit, ...
            '^mm/d/.+$|^mm/.+/d$','once'))
        th_min(j) = th_min(j)/dt;
        th_max(j) = th_max(j)/dt;
        suffix = regexprep(unit, ...
            '^mm/d/|^mm/|/d$','');
        par_units{j} = [local_flux_unit(dt) '/' suffix];
    end
end

end

function unit = local_rate_unit(dt)

if dt == 24
    unit = '1/h';
else
    unit = '1/15 min';
end

end

function unit = local_flux_unit(dt)

if dt == 24
    unit = 'mm/h';
else
    unit = 'mm/15 min';
end

end

function unit = local_time_unit(dt)

if dt == 24
    unit = 'h';
else
    unit = '15 min';
end

end

function udir = user_model_dir_from_root(rootDir)
%USER_MODEL_DIR_FROM_ROOT Locate external user_model directory.

if nargin < 1 ...
        || isempty(rootDir)
    rootDir = '';
end

candidates = {};
if isdeployed
    candidates{end+1,1}=sage_runtime_asset_folder('user_model',rootDir);
end

% 1. Prefer the installed AI-assisted model bound to the configured root.
% Path discovery may still point at the legacy GR4JB example and must not
% override an explicitly installed AI candidate.
if ~isdeployed && ~isempty(rootDir)
    candidates{end+1,1} = fullfile(rootDir, ...
        'AI_assisted_model','active');
    candidates{end+1,1} = fullfile(rootDir, ...
        'Software','AI_assisted_model','active');
end

% 2. If crr_user_model is already visible, retain it as a fallback.
try
    p = which(['crr_user_model.' mexext]);
    if ~isempty(p)
        candidates{end+1,1} = fileparts(p);
    end
catch
end

% 3. Deployed executable folder.
try
    if isdeployed
        exePath = matlab.internal.language. ...
            introspective.getExecutablePath;
        exeDir = fileparts(exePath);

        candidates{end+1,1} = fullfile(exeDir, ...
            'user_model');
        candidates{end+1,1} = exeDir;
    end
catch
end

% 4. Root-based source/development locations.
if ~isempty(rootDir)
    candidates{end+1,1} = fullfile(rootDir, ...
        'user_model');
    candidates{end+1,1} = fullfile(rootDir, ...
        'SAGEhydrology','user_model');
    candidates{end+1,1} = fullfile(rootDir, ...
        'Software','user_model');
end

% 5. MATLAB current-folder fallback.
try
    candidates{end+1,1} = fullfile( ...
        pwd,'user_model');
    candidates{end+1,1} = pwd;
catch
end

for k = 1:numel(candidates)
    c = candidates{k};
    if isfolder(c) && ...
            isfile(fullfile(c, ...
            'user_model_info.mat'))
        udir = c;
        return
    end
end

error('user_model_dir_from_root:NotFound', ...
    ['Cannot locate user_model ' ...
    'folder. Expected to find ' ...
     'user_model_info.mat in ' ...
     'the user_model folder.']);

end
