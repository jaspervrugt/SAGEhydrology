function [mdl,misc,status] = crr_prepare_backend(mdl,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRR_PREPARE_BACKEND Select and prepare the CRR backend.
%
%  Selects the unified C++ router or a MATLAB/model-specific fallback.
%
% SYNOPSIS:
%   [mdl,misc,status] = crr_prepare_backend(mdl,misc)
%
% INPUT ARGUMENTS:
%   mdl             hydrologic model and solver settings
%    .model          model identifier
%    .mcode          numerical implementation; 4 requests C++
%   misc            runtime settings
%    .crr_backend    requested backend: 'cpp' or 'matlab'
%
% OUTPUT ARGUMENTS:
%   mdl             model settings with selected backend and mcode
%   misc            runtime settings with selected backend
%   status          selection and fallback status message
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
        error('crr_prepare_backend:InvalidModel', ...
            'mdl must be a scalar structure.');
    end
    if ~isfield(mdl,'mcode') ...
            || isempty(mdl.mcode)
        mdl.mcode = 4;
    end

    modelName = local_model_name(mdl);
    requestedBackend = 'cpp';
    if isfield(misc,'crr_backend') && ~isempty(misc.crr_backend)
        requestedBackend = lower(char(string(misc.crr_backend)));
    end
    if ~ismember(requestedBackend,{'cpp','matlab'})
        error('crr_prepare_backend:InvalidBackend', ...
            'misc.crr_backend must be ''cpp'' or ''matlab''.');
    end
    isUserModel = strcmpi(modelName,'user_model') ...
        || (isfield(mdl,'model') ...
        && isequal(double(mdl.model),99));
    isCentralModel = isfield(mdl,'model') ...
        && ismember(double(mdl.model),1:7);
    isGCHM = strcmpi(modelName,'gchm') ...
        || (isfield(mdl,'model') ...
        && isequal(double(mdl.model),11));

    if mdl.mcode ~= 4
        misc.crr_backend = 'matlab';
        mdl.crr_backend = misc.crr_backend;
        status = sprintf(['MATLAB hydrologic ' ...
            'solver selected (mcode = %d)'], ...
            mdl.mcode);
        return
    end

    if isCentralModel && strcmp(requestedBackend,'cpp')
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
        if strcmp(requestedBackend,'matlab')
            status = sprintf('MATLAB router with standalone %s MEX (%s)', ...
                upper(modelName),standaloneDetail);
        else
            status = sprintf(['MATLAB router with standalone %s ' ...
                'MEX (%s); unified backend unavailable: %s'], ...
                upper(modelName),standaloneDetail,centralDetail);
        end
        return
    end

    if isGCHM
        warning('crr_prepare_backend:GCHMUnavailable', ...
            ['gchm is a private development model and is not ' ...
             'available in this SAGE distribution.']);
        error('crr_prepare_backend:GCHMUnavailable', ...
            ['gchm requires private/gchm source and a compatible standalone ' ...
             'crr_gchm MEX. %s'],standaloneDetail);
    end

    if isUserModel
        error('crr_prepare_backend:UserModelUnavailable', ...
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
    %LOCAL_PREPARE_CENTRAL Locate or build the unified native backend.

        ok = exist('crr_model_mex','file') == 3;
        if ok
            [isCurrent,~] = local_central_binary_is_current();
            if isCurrent
                detail = 'existing binary found';
                return
            end
            ok = false;
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

function [tf,detail] = local_central_binary_is_current()
%LOCAL_CENTRAL_BINARY_IS_CURRENT Detect source changes after MEX creation.

    tf = true;
    detail = '';
    mexPath = which('crr_model_mex');
    buildFcn = which('compile_crr_model_mex');
    if isempty(mexPath) || isempty(buildFcn)
        return
    end
    modelDir = fileparts(buildFcn);
    sources = {fullfile(modelDir,'crr_model_mex.cpp'), ...
        fullfile(modelDir,'delta_n.hpp'), ...
        fullfile(modelDir,'hymod','hymod.cpp'), ...
        fullfile(modelDir,'hymod','hymod.hpp'), ...
        fullfile(modelDir,'hmodel','hmodel.cpp'), ...
        fullfile(modelDir,'hmodel','hmodel.hpp'), ...
        fullfile(modelDir,'sacsma','sacsma.cpp'), ...
        fullfile(modelDir,'sacsma','sacsma.hpp'), ...
        fullfile(modelDir,'Xinanjiang','xinanjiang.cpp'), ...
        fullfile(modelDir,'Xinanjiang','xinanjiang.hpp'), ...
        fullfile(modelDir,'gr4jA','gr4jA.cpp'), ...
        fullfile(modelDir,'gr4jA','gr4jA.hpp'), ...
        fullfile(modelDir,'hbv','hbv.cpp'), ...
        fullfile(modelDir,'hbv','hbv.hpp'), ...
        fullfile(modelDir,'cfe_nwm','cfe_nwm.cpp'), ...
        fullfile(modelDir,'cfe_nwm','cfe_nwm.hpp'), ...
        buildFcn};
    mexInfo = dir(mexPath);
    if isempty(mexInfo)
        return
    end
    for i = 1:numel(sources)
        sourceInfo = dir(sources{i});
        if ~isempty(sourceInfo) ...
                && sourceInfo.datenum > mexInfo.datenum
            tf = false;
            detail = ['existing binary predates updated source code'];
            return
        end
    end
end


function name = local_model_name(mdl)
%LOCAL_MODEL_NAME Resolve the selected model name for messages.

    name = 'selected model';
    if ~isfield(mdl,'model') || isempty(mdl.model)
        return
    end
    try
        name = sage_model_name(mdl.model);
    catch
        name = 'selected model';
    end
end
