function [swe,Jswe,sm,Jsm] = ...
    model_named_arrays(n,d,needSwe,needJswe,needSm,needJsm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL_NAMED_ARRAYS Allocate requested named-state outputs.
%
%  Allocates only the state-output arrays requested by the caller.
%
% SYNOPSIS:
%   [swe,Jswe,sm,Jsm] = model_named_arrays(...)
%
% INPUT ARGUMENTS:
%   n               number of output times
%   d               number of model parameters
%   needSwe         whether SWE output is needed
%   needJswe        whether SWE Jacobian is needed
%   needSm          whether SM output is needed
%   needJsm         whether SM Jacobian is needed
%
% OUTPUT ARGUMENTS:
%   swe             n-by-1 SWE array, or empty
%   Jswe            n-by-d SWE Jacobian, or empty
%   sm              n-by-1 SM array, or empty
%   Jsm             n-by-d SM Jacobian, or empty
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    swe = [];
    Jswe = [];
    sm = [];
    Jsm = [];
    if needSwe, swe = nan(n,1); end
    if needJswe, Jswe = nan(n,d); end
    if needSm, sm = nan(n,1); end
    if needJsm, Jsm = nan(n,d); end
end
