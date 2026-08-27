function [loss_value,out] = crr_model_cpp(x,mdl,dat,ode,loss,request)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRR_MODEL_CPP Run the C++ branch of crr_model and return the total loss
% and sensitivity vectors and/or matrices.
%
% SYNOPSIS: [loss_value,out] = crr_model_cpp(x,mdl,dat,ode,loss,request)
%   x           dx1 vector of parameter values [check pspace]
%   mdl         structure with model state/parameter info
%    .model      choice of model
%                 1 hymod
%                 2 hmodel
%                 3 sacsma
%                 4 xinanjiang
%                 5 gr4j
%                 6 hbv
%                 7 cfe_nwm
%                 8 user-specified model
%                21 gr4jB [analytic routing]
%                31 hymod data-assimilation implementation
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
%    .sp_method  string with split design
%                 'manual'
%                 'deterministic_block'
%                 'random_block'
%                 'random'
%                 'deterministic_kfold'
%                 'random_kfold'
%    .names      list of model names
%    .y0         mx1 vector of initial states
%    .pspace     0:hydrologic, 1:unit cube, 2:unconstrained parameters
%    .th_min     dx1 vector of lower parameter values [= in pspace]
%    .th_max     dx1 vector of upper parameter values [= in pspace]
%    .par_names  1xd cell parameter names
%    .tout       1x(n+m+1)+spinup*365 vector model output times
%    .idx        1x(n+m+1) vector of indices train&val periods
%    .id_train   1xn vector of indices training period
%    .id_eval    1xm vector of indices evaluation period
%   dat         structure with meteo, discharge and cache data
%   ode         structure with settings for ODE solver
%    .InitStep   Initial time step
%    .MaxStep    Maximum time step
%    .MinStep    Minimum time step
%    .RelTol     Relative tolerance
%    .AbsTol     Absolute tolerance
%    .Order      Order
%    .maxiter    Maximum number of iterations
%   loss        loss-function settings
%    .fnc        scalar choice of loss function
%                  1 sum of absolute residuals
%                  2 generalized least squares
%                  3 Nash-Sutcliffe efficiency
%                  4 Kling-Gupta efficiency
%                  5 Huber loss
%                  6 Flow duration curve loss
%                  7 Jawad Kling-Gupta efficiency
%    .n_win      for fnc = 7, moving-average window length in days
%    .method     scalar JKGE benchmark method
%                  1 = moving-average mean
%                  2 = section-wise mean
%                  3 = long-term mean
%                  4 = monthly climatology
%    .meta       metadata structure [method = 4]
%   request     logical output-selection structure
%    .q          return simulated discharge
%    .gradient   return the loss gradient
%    .jacobian   return the discharge Jacobian on the training mask
%    .metrics    return compact train/evaluation metrics
%    .attribution return total and net parameter attribution
%    .states     []/false, "all", or requested state-name strings
%   out         structure containing only the requested results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pspace = 2: vartheta = (-∞,∞)
% pspace = 1: ntheta = [0,1]
% pspace = 0: theta = tabulated ranges [= hydrologic parameter values]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6
        request = crr_request();
    elseif ~(isstruct(request) ...
            && isscalar(request) ...
            && isfield(request,'normalized') ...
            && isequal(request.normalized,true))
        request = crr_request(request);
    end

    nativeModels = 1:7;
    if ~ismember(mdl.model,nativeModels)
        error('crr_model_cpp:UnsupportedModel', ...
            ['The unified C++ backend does not support model %d. ' ...
             'Backend selection must be resolved by prepare_crr_backend ' ...
             'before crr_model_cpp is called.'],mdl.model);
    end
    if mdl.mcode ~= 4
        error('crr_model_cpp:InvalidMcode', ...
            ['The unified C++ backend requires mdl.mcode = 4; received %d. ' ...
             'Backend selection must be resolved by prepare_crr_backend.'], ...
            mdl.mcode);
    end
    if ~ismember(loss.fnc,1:7)
        error('crr_model_cpp:UnsupportedLoss', ...
            'The unified C++ backend does not support loss.fnc = %d.', ...
            loss.fnc);
    end
    if exist('crr_model_mex','file') ~= 3
        error('crr_model_cpp:MissingMex', ...
            ['The unified C++ backend was selected, but crr_model_mex is ' ...
             'not available. Run prepare_crr_backend before simulation.']);
    end

    req = request;
    if isempty(req.states)
        req.states = [];
    end
    [loss_value,out] = crr_model_mex(x,mdl,dat,ode,loss,req);
end
