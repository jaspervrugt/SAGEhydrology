function named = model_named_result(mem,Jth,needSwe,needJswe, ...
    needSm,needJsm,swe,Jswe,sm,Jsm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL_NAMED_RESULT Package direct named-state outputs.
%
%  Returns an empty structure when the full state history is retained.
%
% SYNOPSIS:
%   named = model_named_result(...)
%
% INPUT ARGUMENTS:
%   mem             state-history mode; direct outputs require zero
%   Jth             parameter-space scaling factors
%   needSwe         whether SWE output is requested
%   needJswe        whether SWE Jacobian is requested
%   needSm          whether SM output is requested
%   needJsm         whether SM Jacobian is requested
%   swe             simulated SWE
%   Jswe            SWE Jacobian before parameter scaling
%   sm              simulated SM
%   Jsm             SM Jacobian before parameter scaling
%
% OUTPUT ARGUMENTS:
%   named           requested named observations and Jacobians
%    .obs            SWE and/or SM observation arrays
%    .jac            scaled SWE and/or SM Jacobian arrays
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    named = struct();
    if mem ~= 0
        return
    end
    scale = reshape(Jth,1,[]);
    if needSwe
        named.obs.SWE = swe;
    end
    if needJswe
        named.jac.SWE = Jswe .* scale; 
    end
    if needSm
        named.obs.SM = sm; 
    end
    if needJsm
        named.jac.SM = Jsm .* scale; 
    end
end
