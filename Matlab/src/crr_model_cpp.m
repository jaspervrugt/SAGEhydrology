function [loss_value,out] = crr_model_cpp(x,mdl,dat,ode,loss,request)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRR_MODEL_CPP Evaluate a CRR model with the C++ backend.
%
%  Routes simulated discharge, requested sensitivities, and loss evaluation
%  through the unified C++ MEX interface.
%
% SYNOPSIS:
%   [loss_value,out] = crr_model_cpp(x,mdl,dat,ode,loss,request)
%
% INPUT ARGUMENTS:
%   x               hydrologic parameters in mdl.pspace
%   mdl             model, parameter, and output-time settings
%    .model          supported built-in model identifier
%    .pspace         selected parameter space (0, 1, or 2)
%    .th_min         lower parameter bounds
%    .th_max         upper parameter bounds
%    .id_train       training observation indices
%    .id_eval        evaluation observation indices
%   dat             forcing and prepared observations
%   ode             numerical integration settings
%   loss            selected loss and benchmark settings
%   request         requested outputs and gradients
%
% OUTPUT ARGUMENTS:
%   loss_value      scalar objective value
%   out             requested simulations, sensitivities, and diagnostics
%
% NOTES:
%   User-defined model 99 is supported by the MATLAB backend only.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6
        request = crr_request();
    elseif ~(isstruct(request) && ...
        isscalar(request) && ...
        isfield(request,'normalized') && ...
        isequal(request.normalized,true))
        request = crr_request(request);
    end

    nativeModels = 1:7;
    if ~ismember(mdl.model,nativeModels)
        error('crr_model_cpp:UnsupportedModel', ...
            ['The unified C++ backend does ' ...
            'not support model %d. ' ...
             'Backend selection must be ' ...
             'resolved by crr_prepare_backend ' ...
             'before crr_model_cpp is called.'],mdl.model);
    end
    if mdl.mcode ~= 4
        error('crr_model_cpp:InvalidMcode', ...
            ['The unified C++ backend requires ' ...
            'mdl.mcode = 4; received %d. ' ...
             'Backend selection must be ' ...
             'resolved by crr_prepare_backend.'], ...
            mdl.mcode);
    end
    if ~ismember(loss.fnc,1:7)
        error('crr_model_cpp:UnsupportedLoss', ...
            ['The unified C++ backend does ' ...
            'not support loss.fnc = %d.'], ...
            loss.fnc);
    end
    if exist('crr_model_mex','file') ~= 3
        error('crr_model_cpp:MissingMex', ...
            ['The unified C++ backend was ' ...
            'selected, but crr_model_mex is ' ...
             'not available. Run crr_prepare_backend ' ...
             'before simulation.']);
    end

    mdl = apply_model_configuration(mdl,ode);
    [physical,slope] = parameter_scale('physical',mdl,x(:));
    normalized = parameter_scale('normalized',mdl,physical);
    mapping = mdl; mapping.pspace = 1;
    [~,normalizedSlope] = parameter_scale( ...
        'physical',mapping,normalized);
    mdl.attribution_scale = normalizedSlope./max(slope,realmin);

    [th,~] = parameter_scale('physical',mdl,x(:));
    modes = parameter_scale('modes',mdl);
    mdl.parameter_log = double(modes == "log");
    validate_initial_storage(mdl,th);
    req = request;
    if isempty(req.states)
        req.states = [];
    end
    standardizedQ = isempty(req.obs) && ...
        isempty(req.jac) && ...
        strcmp(observation_transform(loss,"Q"),'zscore');
    if standardizedQ
        req.q = true;
        req.jacobian = true;
    end
    [loss_value,out] = crr_model_mex(x,mdl,dat,ode,loss,req);
    if standardizedQ
        [loss_value,out] = zscore_q_result(out,dat,mdl,request);
    end
end
