function [needQ,needJq] = ...
    model_q_requests(structured,request,needJacobian)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL_Q_REQUESTS Decode discharge output requests.
%
%  Supports both legacy and structured discharge requests.
%
% SYNOPSIS:
%   [needQ,needJq] = model_q_requests(...)
%
% INPUT ARGUMENTS:
%   structured      whether the caller made a structured request
%   request         requested observations and Jacobians
%    .q              legacy discharge request flag
%    .obs            requested observation names
%    .jacobian       legacy Jacobian request flag
%    .jac            requested Jacobian names
%   needJacobian    legacy Jacobian requirement
%
% OUTPUT ARGUMENTS:
%   needQ           true if discharge output is requested
%   needJq          true if discharge Jacobian is requested
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~structured
        needQ = true;
        needJq = needJacobian;
        return
    end
    needQ = request.q || any(request.obs == "Q");
    needJq = request.jacobian || any(request.jac == "Q");
end
