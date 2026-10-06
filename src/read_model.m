function [mdl,d] = read_model(mdl,prd,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_MODEL Read hydrologic model parameter information.
%
%  Sets model-specific states, parameter bounds, names, and time settings.
%
% SYNOPSIS:
%   [mdl,d] = read_model(mdl,prd)
%   [mdl,d] = read_model(mdl,prd,verbose)
%
% INPUT ARGUMENTS:
%   mdl             model, implementation, and assessment settings
%    .model          hydrologic model identifier
%    .mcode          numerical implementation
%    .calc           sequential or parallel model execution
%    .mode           assessment design code
%    .eval_mode      evaluation basin/period scenario
%   prd             time resolution and spin-up period
%    .dt             samples per day: 1, 24, or 96
%    .spinup         spin-up period in days
%   verbose         optional status-print switch; default true
%
% OUTPUT ARGUMENTS:
%   mdl             updated model structure
%    .y0             initial model states
%    .pspace         parameter-space code
%    .th_min         lower parameter bounds
%    .th_max         upper parameter bounds
%    .par_names      parameter names
%    .par_units      parameter units
%   d               number of hydrologic parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3
        verbose = true;
    end
    if nargin < 2
        dt = 1;
        prd = struct('dt',dt);
    else
        dt = prd.dt;
    end

    % Retain the active temporal resolution with the model metadata so
    % diagnostics and exported figures can identify it without requiring
    % a separate period structure.
    mdl.dt = dt;
    % Requirement flags belong to the freshly selected model, not a
    % previously loaded AI model. Its metadata reader sets them below.
    mdl.ai_requires_kosugi = false;
    mdl.ai_requires_dual_kosugi = false;
    
    if verbose 
        fprintf('... Reading model information');
    end
    
    switch dt
        case 1      % daily data
        case 24     % hourly data
        case 96     % 15-minute data
        otherwise
            fprintf('\n');
            error(['read_model: improper temporal ' ...
                'resolution dt = %g.'],dt);
    end
    
    % ---------------------------
    % 1) Daily master definitions
    % ---------------------------
    if mdl.model == 99
        [mdl,d] = read_user_model_info( ...
            mdl,prd,false);
        if verbose
            fprintf(' ... Done\n');
            fprintf(['      Loaded AI-assisted model metadata ' ...
                'with %d parameters\n'],d);
        end
        return
    end
    if mdl.model == 12
        plugin = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
            'private','mcp_salo');
        if isdeployed || ~isfile(fullfile(plugin,'read_mcp_salo_info.m'))
            error('read_model:MCPUnavailable','Private MN(DS)_SALO(5) is unavailable.');
        end
        addpath(plugin,'-begin');
        readMCPInfo = str2func('read_mcp_salo_info');
        [mdl,d] = readMCPInfo(mdl,prd,verbose);
        return
    end
    if mdl.model == 11
        if exist('read_gchm_info','file') ~= 2
            warning('read_model:GCHMUnavailable', ...
                ['gchm is a private development model and is not ' ...
                 'available in this SAGE distribution.']);
            error('read_model:GCHMUnavailable', ...
                ['Cannot run model 11 because the private gchm ' ...
                 'plug-in is unavailable.']);
        end
        readGCHMInfo = str2func('read_gchm_info');
        [mdl,d] = readGCHMInfo(mdl,prd,verbose);
        if d ~= 30 ...
                || numel(mdl.y0) ~= 4
            error('read_model:gchmVersion',['Model ID 11 gchm ' ...
                'requires 30 parameters and 4 numerical states.']);
        end
        return
    end
    
    [par_info,m,add_info] = ...
        local_get_model_par_info_daily(mdl.model,dt);
    
    par_names = par_info(:,2)';         % model parameter names
    par_desc = par_info(:,3)';          % short parameter descriptions
    par_units = par_info(:,4)';         % units of parameters
    th_min = cell2mat(par_info(:,5));   % lower bounds (daily master)
    th_max = cell2mat(par_info(:,6));   % upper bounds (daily master)
    
    % --------------------------------
    % Add snow module where applicable
    % --------------------------------
    if ismember(mdl.model,[1 2 3 4 5 7])
        % ==============================
        % Add snow to all models but HBV
        % ==============================
        par_snow = {'T_{\rm tr}','f_{\rm dd}'};
        par_snow_desc = {'snow/rain temperature threshold', ...
                         'degree-day melt factor'};
        
        par_names = [par_names, par_snow];
        par_desc  = [par_desc, par_snow_desc];
        par_units = [par_units, {'°C','mm/d/°C'}];
        th_min = [th_min ; -3 ; 0.10 ];
        th_max = [th_max ;  3 ; 10.0 ];
        m = m + 1;                      % # state variables with Swe
    end
    
    d = numel(par_names);               % # parameters
    
    mdl.th_min_daily = th_min;
    mdl.th_max_daily = th_max;
    mdl.par_units_daily = par_units;
    
    % ------------------------------------
    % 2) Daily -> selected time resolution
    % ------------------------------------
    if dt > 1
        [th_min,th_max,par_units] = ...
            convert_daily_to_subdaily( ...
            th_min,th_max,par_units,dt);
    end
    
    % ------------------------------------
    % 3) Resolution-specific overrides
    % ------------------------------------
    if dt == 24
        [th_min,th_max,par_units,add_info] = ...
            apply_hourly_overrides( ...
            mdl.model,par_names,th_min, ...
            th_max,par_units,add_info);
    
    elseif dt == 96
        [th_min,th_max,par_units,add_info] = ...
            apply_15min_overrides( ...
            mdl.model,par_names,th_min, ...
            th_max,par_units,add_info);
    end
    
    % ---------------------------------------------------------------
    % Define initial states, par ranges, parameter names, units, etc.
    % ---------------------------------------------------------------
    
    mdl.y0 = 1e-5*ones(m,1);    % Initial states
    mdl.pspace = 1;             % Treatment of parameter space
                                  % [0] hydrologic parameter values
                                  % [1] unit cube
                                  % [2] unconstrained (-∞,∞)
    mdl.par_names = par_names;  % Parameter names
    mdl.par_desc = par_desc;    % Parameter descriptions
    mdl.par_units = par_units;  % Parameter units
    mdl.th_min = th_min;        % Minimum parameter values
    mdl.th_max = th_max;        % Maximum parameter values
    
    if mdl.model == 5           % gr4j
        mdl.n1 = add_info.n1;   % Copy UH1 cascade order
        mdl.n2 = add_info.n2;   % Copy UH2 cascade order
    end
    if mdl.model == 7           % cfe_nwm
        mdl.giuh_ordnts = ...   % Ordinates unit hydrograph
            [0.1,0.35,0.2,0.14,0.1,0.06,0.05];
        % --> sum to 1
    end
    
    if verbose
        fprintf(' ... Done\n');
    end
end

% -------------
% Local helpers
% -------------
function [par_info,m,add_info] = ...
    local_get_model_par_info_daily(model,dt)
%LOCAL_GET_MODEL_PAR_INFO_DAILY Daily master parameter definitions
% par_info columns:
%   {idx, par_name, par_desc, par_unit, th_min, th_max}

    add_info = struct();
    
    switch model
    
        case 1  % HYMOD
            % Boyle, D. P. (2000). Multicriteria calibration of hydrologic
            %   models. PhD dissertation, University of Arizona.
            par_info = {
                1, 's_{\rm u,max}', 'maximum unsaturated storage', ...
                    'mm',    50,  5000;
                2, 'b',             'soil capacity variability', ...
                    '-',   0.10,  10.0;
                3, 'a',             'flow partition coefficient', ...
                    '-',   1e-5,   1.0;
                4, 'k_{\rm s}',     'slow reservoir recession', ...
                    '1/d', 1e-5,   1.0; % upper was 0.5 ** TODAY
                5, 'k_{\rm f}',     'fast reservoir recession', ...
                    '1/d',  0.1,  100}; % lower was 0.5, upper was 10 ** TODAY
            m = 6;
    
        case 2  % HMODEL
            % Schoups, G., Vrugt, J.A., Fenicia, F., & van de Giesen, N.C.
            %   (2010). Corruption of accuracy and efficiency of Markov 
            %   chain Monte Carlo simulation by inaccurate numerical 
            %   implementation of conceptual hydrologic models. Water 
            %   Resources Research, 46.        
            par_info = {
                1, 'I_{\rm max}',   'maximum interception storage', ...
                    'mm',     1,    20;
                2, 'S_{\rm u,max}', 'soil water capacity', ...
                    'mm',  10.0,  1000;
                3, 'Q_{\rm max}',   'maximum percolation rate', ...
                    'mm/d',1e-1,   100;
                4, 'a_{\rm E}',     'evaporation parameter', ...
                    '-',   1e-5,    10;
                5, 'a_{\rm F}',     'runoff parameter', ...
                    '-',    -10,    10;
                6, 'r_{\rm f}',     'fast reservoir time constant', ...
                    'd',    0.1,    10;
                7, 'r_{\rm s}',     'slow reservoir time constant', ...
                    'd',    1.0,   250};
            m = 5;
    
        case 3  % SACSMA
            % Burnash, R. J. C. (1973). A river basin simulation model for
            %   hydrology and flood operations. Technical Report, 
            %   U.S. National Weather Service.
            % Clark, M.P., Slater, A.G., Rupp, D.E., Woods, R.A.,
            %   Vrugt, J.A., Gupta, H.V., Wagener, T., & Hay, L.E. (2008).
            %   Framework for understanding structural errors (FUSE): A 
            %   modular framework to diagnose differences between 
            %   hydrological models. Water Resources Research, 44.        
            par_info = {
                1,  'u_{\rm f,max}',  'upper free-water storage', ...
                    'mm',     5,   500;
                2,  'u_{\rm t,max}',  'upper tension storage', ...
                    'mm',     1,   200;   % used to be 150 ** TODAY
                3,  'l_{\rm fp,max}', 'lower primary free storage', ...
                    'mm',     1,   1500;  % used to be 500 ** TODAY
                4,  'l_{\rm fs,max}', 'lower supplemental free storage', ...
                    'mm',     1,   500;   % used to be 750 ** TODAY
                5,  'l_{\rm t,max}',  'lower tension storage', ...
                    'mm',     5,  1000;
                6,  '\alpha',         'lower-layer percolation factor', ...
                    '-',   1e-2,   400;   % used to be 150 ** TODAY
                7,  '\psi',           'lower-layer percolation exponent', ...
                    '-',    0.1,    50;
                8,  'k_{\rm i}',      'interflow depletion rate', ...
                    '1/d', 1e-3,   100;   % used to be 50
                9,  '\kappa',         'percolation fraction to tension store',  ...
                    '-',   1e-5,  0.95;
                10, '\nu_{\rm p}',    'primary baseflow depletion', ...
                    '1/d', 1e-5,  0.50;   % lower was 1e-4 ** TODAY
                11, '\nu_{\rm s}',    'secondary baseflow depletion', ...
                    '1/d', 1e-5,  0.50;   % lower was 1e-4, upper 0.25 ** TODAY
                12, 'a_{\rm c,max}',  'maximum saturated area fraction', ...
                    '-',   1e-4,  0.20;   % lower was 1e-3, upper 0.15 ** TODAY
                13, 'k_{\rm f}',      'routing reservoir recession', ...
                    '1/d',   0.1,  100};  % used to be 10 ** TODAY
            m = 9;
    
        case 4  % XINANJIANG
            % Zhao, R. J. (1980). The Xinanjiang model applied in China.
            %   Journal of Hydrology, 135–155.
            par_info = {
                1,  'f_{\rm p}',    'PET-to-pan factor', ...
                    '-',      0.3,     1.7;
                2,  'A_{\rm im}',   'impervious area fraction', ...
                    '-',     1e-4,    0.20;
                3,  'a',            'tension-water inflection', ...
                    '-', -0.49999, 0.49999;
                4,  'b',            'tension-water shape', ...
                    '-',     1e-1,       5;    % used to be 2.0
                5,  'f_{\rm wm}',   'fraction of total as max tension', ...
                    '-',     1e-3,   0.999;    % used to be 0.95 ** TODAY
                6,  'f_{\rm lm}',   'first evaporation threshold fraction', ...
                    '-',     1e-3,    0.95;    % used to be 0.8 ** TODAY
                7,  'c',            'second evaporation threshold factor', ...
                    '-',     1e-3,       1;    % lower was 1e-2 ** TODAY
                8,  's_{\rm tot}',  'total soil moisture storage', ...
                    'mm',       25,    2e3;    % used to be 1,000
                9,  '\beta',        'free-water shape parameter', ...
                    '-',      0.1,       5;    % used to be 3
                10, 'k_{\rm i}',    'interflow parameter', ...
                    '1/d',    1e-5,    100;    % used to be 1 ** TODAY
                11, 'k_{\rm g}',    'groundwater parameter', ...
                    '1/d',    1e-5,    100;    % used to be 1 ** TODAY
                12, 'c_{\rm i}',    'interflow time coefficient', ...
                    '1/d',    1e-5,    100;    % used to be 1 ** TODAY
                13, 'c_{\rm g}',    'baseflow time coefficient', ...
                    '1/d',    1e-5,    100;    % used to be 1 ** TODAY
                14, 'k_{\rm f}',    'routing reservoir recession', ...
                    '1/d',     0.1,    100};   % used to be 10 ** TODAY
            m = 8;
    
        case 5  % GR4J
            % Perrin, C., Michel, C., & Andréassian, V. (2003). Improvement
            %   of a parsimonious model for streamflow simulation. 
            %   Journal of Hydrology, 279, 275–289.
            par_info = {
                1, 'x_{1}',        'production store capacity', ...    
                    'mm',     1,  4000;
                2, 'x_{2}',        'groundwater exchange coeff.', ...  
                    'mm/d', -20,    15;
                3, 'x_{3}',        'routing store capacity', ...
                    'mm',     1,  2000;
                4, 'x_{4}',        'routing time base', ...
                    'd',    0.1,    10;
                5, 'x_{5}',        'fast/slow flow partition', ...
                    '-',   1e-3,   1.0;
                6, 'f_{\rm p}',    'PET-to-pan factor', ...
                    '-',    0.3,   1.7};
            % JAV: added x5: flow partitioning factor
            % JAV: added f_p, to mitigate Ep bias/errors
    
            switch dt
                case 1
                    add_info.n1 = 2; add_info.n2 = 4;
                case 24
                    add_info.n1 = 2; add_info.n2 = 4;
            end
            m = add_info.n1 + add_info.n2 + 3;
    
        case 6  % HBV
            % Bergström, S. (1976). Development and application of a 
            %   conceptual runoff model for Scandinavian catchments. 
            %   SMHI Report RHO No. 7.
            % Bergström, S. (1995). The HBV model. In: Computer Models of
            %   Watershed Hydrology, Water Resources Publications.
            par_info = {
                1,  'f_{\rm c}',             'soil field capacity', ...
                    'mm',       25,  1500;  % used to be 50,800
                2,  '\beta',                 'soil nonlinearity exponent', ...
                    '-',       0.1,     6;
                3,  '\ell_{\rm lp}',         'ET limit parameter', ...
                    '-',       0.1,   1.0;
                4,  'k_{0}',                 'quickflow recession coefficient', ...
                    '1/d',    1e-3,   1.2;  % upper was 1.0 ** TODAY
                5,  'u_{\rm zl}',            'upper-zone quickflow threshold', ...
                    'mm',        0,   300;  % upper was 200 ** TODAY
                6,  'k_{1}',                 'upper-zone recession coefficient', ...
                    '1/d',    1e-5,   1.0;   % lower was 1e-3 ** TODAY
                7,  'k_{2}',                 'lower-zone recession coefficient', ...
                    '1/d',    1e-5,   1.0;   % lower 1e-4, upper 0.5 ** TODAY
                8,  '{\rm perc}',            'upper-to-lower percolation', ...
                    'mm/d',      0,    15;   % upper was 10 ** TODAY
                9,  'T_{\rm tr}',            'snow/rain temperature threshold', ...
                    '°C',       -3,     3;
                10, 'f_{\rm dd}',            'degree-day melt factor', ...
                    'mm/d/°C',  0.1,   10;
                11, '{\rm sf}_{\rm cf}',     'snowfall correction factor', ...
                    '-',       0.5,   1.5;
                12, 'f_{\rm r}',             'refreezing factor', ...
                    '-',       0.0,   0.2;
                13, 'b_{\rm rt}',            'routing width', ...
                    'd',       0.5,    25};
            m = 5;
            % renamed cf_{\rm max} to f_{\rm dd} --> rename in C++ code
    
        case 7  % CFE-NWM
            par_info = {
                1,  's_{\rm max}',   'maximum soil storage', ...
                'mm',    25,  2000;
                2,  's_{\rm fc}',    'field capacity fraction', ...
                '-',   0.50,  0.95;
                3,  's_{\rm wp}',    'wilting point fraction', ...
                '-',   1e-3,  0.50;
                4,  'k_{\rm sch}',   'soil scheme coefficient', ...
                '1/d', 1e-4,  25.0;
                5,  'a_{1}',         'surface runoff exponent', ...
                '-',   0.01,  10.0;
                6,  'k_{\rm perc}',  'percolation coefficient', ...
                '1/d', 1e-4,   100;
                7,  'lf_{\rm thr}',  'lateral-flow threshold', ...
                '-',   1e-5,   1.0;
                8,  'a_{2}',         'lateral-flow exponent', ...
                '1/d',  0.5,   6.0;
                9,  'k_{\rm lf}',    'lateral-flow coefficient', ...
                '1/d', 1e-4,   1e4;
                10, 'g_{\rm max}',   'maximum groundwater discharge', ...
                '1/d',    5,  1000;
                11, 'c_{\rm gw}',    'groundwater discharge coeff.', ... 
                '1/d', 1e-5,   1e3;
                12, '{\rm mm}',      'soil moisture exponent', ...       
                '-',   1e-2,    10;
                13, 'k_{\rm nsh}',   'Nash cascade coefficient', ...
                '1/d', 1e-4,  25.0};
            m = 13;
    
        otherwise
            fprintf('\n');
            fprintf(['      Error: read_model: ' ...
                'Unknown model id = %d\n'],model);
            error('read_model: demo_SAGE stops');
    end
end
