function ode = read_numsettings(ode_in)
%READ_NUMSETTINGS Defines solver settings for each watershed
% SYNOPSIS: ode = read_numsettings(ode_in)
%   ode_in      OPTIONAL: structure user-defined settings ODE solver
%   ode         Structure with settings for ODE solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('... Reading numerical ODE settings');
    
    % ---- Defaults ----
    ode_default = struct();
    ode_default.InitStep = 1e-2;    % initial time-step (T)    [1e-3]
    ode_default.MaxStep = 1.0;      % maximum time-step (T)
    ode_default.MinStep = 1e-4;     % minimum time-step (T)
    ode_default.RelTol = 1e-3;      % relative tolerance (-)   [1e-5]
    ode_default.AbsTol = 1e-3;      % absolute tolerances (mm) [1e-5]
    ode_default.Order = 2;          % 2nd order accurate (Heun) method
    ode_default.maxiter = 1e4;      % Maximum iterations numerical solver
    ode_default.mem = 0;            % store state variables? [0:no, 1:yes]
                                    % 0 should be faster for longer time series
    
    % ---- Merge user overrides (if provided) ----
    if nargin < 1 ...
            || isempty(ode_in)
        ode_in = struct();
    end
    if ~isstruct(ode_in)
        error(['      read_numsettings:BadInput', ...
            'ode_in must be a struct (or empty).']);
    end
    
    ode = ode_default;
    fn = fieldnames(ode_in);
    for k = 1:numel(fn)
        ode.(fn{k}) = ode_in.(fn{k});
    end
    
    % Optional: quick sanity checks (lightweight)
    if ode.MinStep <= 0 ...
            || ode.MaxStep <= 0 ...
            || ode.InitStep <= 0
        error(['      read_numsettings:BadSteps', ...
            'InitStep/MinStep/MaxStep must be > 0.']);
    end
    if ode.MinStep > ode.MaxStep
        error(['      read_numsettings:BadSteps', ...
            'MinStep cannot exceed MaxStep.']);
    end
    
    fprintf(' ... Done\n');

end