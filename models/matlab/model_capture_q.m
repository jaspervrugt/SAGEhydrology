function [q,J] = model_capture_q( ...
    z,zprev,row,m,id,d,needQ,needJq,q,J)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL_CAPTURE_Q Store discharge at one output time.
%
%  Computes discharge from the state increment at a print time.
%
% SYNOPSIS:
%   [q,J] = model_capture_q(z,zprev,row,m,id,d,needQ,needJq,q,J)
%
% INPUT ARGUMENTS:
%   z               current augmented model state
%   zprev           previous augmented model state
%   row             output row index
%   m               number of physical model states
%   id              discharge sensitivity state indices
%   d               number of model parameters
%   needQ           whether discharge was requested
%   needJq          whether discharge Jacobian was requested
%   q               preallocated discharge output
%   J               preallocated discharge Jacobian
%
% OUTPUT ARGUMENTS:
%   q               discharge output with current row populated
%   J               discharge Jacobian with current row populated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if needQ
        q(row) = z(m) - zprev(m);
    end
    if needJq
        J(row,1:d) = z(id) - zprev(id);
    end
end
