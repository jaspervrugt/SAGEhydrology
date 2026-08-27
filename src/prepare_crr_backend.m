function [mdl,misc,status] = prepare_crr_backend(mdl,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREPARE_CRR_BACKEND Selects and prepares the CRR execution backend
%
% SYNOPSIS: [mdl,misc,status] = prepare_crr_backend(mdl,misc)
%
%   mdl          Structure with hydrologic-model and solver settings
%    .model       Choice of model
%                  1 hymod
%                  2 hmodel
%                  3 sacsma
%                  4 Xinanjiang
%                  5 gr4j
%                  6 hbv
%                  7 cfe_nwm
%                  8 user_model
%    .mcode       Numerical solution of the hydrologic model
%                  1 Runge-Kutta implementation in MATLAB
%                  2 ode45 implementation in MATLAB
%                  3 explicit Euler implementation in MATLAB
%                  4 Runge-Kutta implementation in C++
%    .names       List of model names
%   misc         Structure with runtime settings
%    .crr_backend Requested CRR routing backend
%                  'cpp'    unified crr_model_mex interface
%                  'matlab' MATLAB router with a standalone model MEX,
%                           or the pure MATLAB implementation
%
%   mdl          OUTPUT: model structure with the selected mcode and a
%                matching crr_backend field for downstream routing
%   misc         OUTPUT: runtime structure with the selected crr_backend
%   status       OUTPUT: text describing the selected backend and any
%                fallback that was required
%
% DESIGN
%   For built-in hydrologic models requested with mdl.mcode = 4, backend
%   selection follows this order:
%   1) Unified C++ MEX
%      Use crr_model_mex with misc.crr_backend = 'cpp'. If the binary is
%      absent in MATLAB mode, compile_crr_model_mex is called once.
%   2) Standalone model MEX
%      Use the MATLAB router with misc.crr_backend = 'matlab' and retain
%      mdl.mcode = 4. compile_model verifies or, in MATLAB mode, builds the
%      selected standalone model MEX.
%   3) Pure MATLAB Runge-Kutta
%      Use misc.crr_backend = 'matlab' and mdl.mcode = 1 when neither C++
%      option is available.
%
% NOTES
%   - Deployed applications never compile source code. They use packaged
%     MEX binaries and otherwise follow the available fallback path.
%   - user_model is external to the unified crr_model_mex library and
%     requires a compatible standalone crr_user_model MEX plug-in.
%   - A user-selected MATLAB solver (mdl.mcode ~= 4) is preserved and does
%     not trigger C++ backend preparation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(misc)
        misc = struct();
    end
    if ~isstruct(mdl) ...
            || ~isscalar(mdl)
        error('prepare_crr_backend:InvalidModel', ...
            'mdl must be a scalar structure.');
    end
    if ~isfield(mdl,'mcode') ...
            || isempty(mdl.mcode)
        mdl.mcode = 4;
    end

    modelName = local_model_name(mdl);
    isUserModel = strcmpi(modelName,'user_model') ...
        || (isfield(mdl,'model') ...
        && isequal(double(mdl.model),8));
    isCentralModel = isfield(mdl,'model') ...
        && ismember(double(mdl.model),1:7);

    if mdl.mcode ~= 4
        misc.crr_backend = 'matlab';
        mdl.crr_backend = misc.crr_backend;
        status = sprintf(['MATLAB hydrologic ' ...
            'solver selected (mcode = %d)'], ...
            mdl.mcode);
        return
    end

    if isCentralModel
        [centralOK,centralDetail] = local_prepare_central();
        if centralOK
            misc.crr_backend = 'cpp';
            mdl.crr_backend = misc.crr_backend;
            mdl.mcode = 4;
            status = ['unified C++ MEX (' centralDetail ')'];
            return
        end
    else
        centralDetail = sprintf( ...
            ['model %s uses the ' ...
            'standalone interface'],modelName);
    end

    misc.crr_backend = 'matlab';
    mdl.crr_backend = misc.crr_backend;
    mdl.mcode = 4;
    [standaloneFlag,standaloneDetail] = compile_model(mdl);
    if standaloneFlag >= 0
        status = sprintf(['MATLAB router with standalone %s ' ...
            'MEX (%s); unified backend unavailable: %s'], ...
            upper(modelName),standaloneDetail, ...
            centralDetail);
        return
    end

    if isUserModel
        error('prepare_crr_backend:UserModelUnavailable', ...
            ['The user_model requires a ' ...
            'compatible precompiled ' ...
             'crr_user_model MEX. %s'],standaloneDetail);
    end

    mdl.mcode = 1;
    misc.crr_backend = 'matlab';
    mdl.crr_backend = misc.crr_backend;
    status = sprintf(['MATLAB RK2; unified backend unavailable: ' ...
        '%s; standalone MEX unavailable: %s'], ...
        centralDetail,standaloneDetail);

    function [ok,detail] = local_prepare_central()
        ok = exist('crr_model_mex','file') == 3;
        if ok
            detail = 'existing binary found';
            return
        end
        if isdeployed
            detail = 'packaged crr_model_mex was not found';
            return
        end
        if exist('compile_crr_model_mex','file') ~= 2
            detail = ['compile_crr_model_mex ' ...
                'is not on the MATLAB path'];
            return
        end
        try
            clear crr_model_mex
            mexFile = compile_crr_model_mex();
            rehash toolboxcache
            ok = isfile(mexFile) ...
                && exist('crr_model_mex','file') == 3;
            if ok
                detail = 'compiled successfully';
            else
                detail = ['compilation did not ' ...
                    'produce a loadable MEX'];
            end
        catch ME
            ok = false;
            detail = ['compilation failed: ' ME.message];
        end
    end
end


function name = local_model_name(mdl)
%LOCAL_MODEL_NAME Resolve the selected model name for messages.

    name = 'selected model';
    if ~isfield(mdl,'model') ...
            || ~isfield(mdl,'names') ...
            || isempty(mdl.model) ...
            || isempty(mdl.names)
        return
    end
    try
        if iscell(mdl.names)
            name = char(mdl.names{mdl.model});
        else
            name = char(string(mdl.names(mdl.model)));
        end
        name = strtrim(name);
    catch
        name = 'selected model';
    end
end
