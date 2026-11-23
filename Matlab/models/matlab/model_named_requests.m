function [needSwe,needJswe,needSm,needJsm] = ...
    model_named_requests(structured,request)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL_NAMED_REQUESTS Decode named-state output requests.
%
%  Returns four logical flags for direct state outputs.
%
% SYNOPSIS:
%   [needSwe,needJswe,needSm,needJsm] = ...
%
% INPUT ARGUMENTS:
%   structured      whether the caller made a structured request
%   request         requested observations and Jacobians
%    .obs            observation names (SWE and/or SM)
%    .jac            Jacobian names (SWE and/or SM)
%
% OUTPUT ARGUMENTS:
%   needSwe         true if SWE observations are requested
%   needJswe        true if SWE Jacobians are requested
%   needSm          true if SM observations are requested
%   needJsm         true if SM Jacobians are requested
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    needSwe = false;
    needJswe = false;
    needSm = false;
    needJsm = false;
    if ~structured
        return
    end
    needSwe = any(request.obs == "SWE");
    needJswe = any(request.jac == "SWE");
    needSm = any(request.obs == "SM");
    needJsm = any(request.jac == "SM");
end
