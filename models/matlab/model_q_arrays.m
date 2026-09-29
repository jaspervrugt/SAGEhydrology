function [q,J] = model_q_arrays(n,d,needQ,needJq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL_Q_ARRAYS Allocate requested discharge outputs.
%
%  Allocates only the discharge arrays requested by the caller.
%
% SYNOPSIS:
%   [q,J] = model_q_arrays(n,d,needQ,needJq)
%
% INPUT ARGUMENTS:
%   n               number of output times
%   d               number of model parameters
%   needQ           whether discharge output is needed
%   needJq          whether discharge Jacobian is needed
%
% OUTPUT ARGUMENTS:
%   q               n-by-1 discharge array, or empty
%   J               n-by-d discharge Jacobian, or empty
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    q = [];
    J = [];
    if needQ, q = nan(n,1); end
    if needJq, J = nan(n,d); end
end
