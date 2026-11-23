function [q,J] = model_history_q( ...
    Z,rows,m,id,needQ,needJq,q,J)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL_HISTORY_Q Extract discharge from stored state history.
%
%  Takes first differences of selected state-history columns.
%
% SYNOPSIS:
%   [q,J] = model_history_q(Z,rows,m,id,needQ,needJq,q,J)
%
% INPUT ARGUMENTS:
%   Z               complete augmented state history
%   rows            output-time row indices
%   m               discharge accumulation state index
%   id              discharge sensitivity state indices
%   needQ           whether discharge was requested
%   needJq          whether discharge Jacobian was requested
%   q               existing discharge output
%   J               existing discharge Jacobian
%
% OUTPUT ARGUMENTS:
%   q               requested discharge increments
%   J               requested discharge sensitivity increments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if needQ, q = diff(Z(rows,m)); end
    if needJq, J = diff(Z(rows,id)); end
end
