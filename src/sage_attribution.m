function [A_t,A_n] = sage_attribution(J,delta,mdl,g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_ATTRIBUTION Returns time-resolved and net parameter attribution
% vectors from the Jacobian matrix and loss-sensitivity vector.
%
% SYNOPSIS:
%   [A_tw,A_net] = sage_attribution(J,delta,mdl)
%   [A_tw,A_net] = sage_attribution(J,delta,mdl,g)
%   J           nxd Jacobian matrix, J(t,j) = dq_t/dtheta_j
%   delta       nx1 loss-sensitivity vector, delta(t) = dloss/dq_t
%   mdl         structure with parameter bounds
%    .th_min     dx1 vector of lower parameter values
%    .th_max     dx1 vector of upper parameter values
%   g           OPTIONAL: dx1 gradient vector; if omitted, g = J'*delta
%   A_t         OUTPUT: dx1 time-resolved parameter attribution vector
%   A_n         OUTPUT: dx1 net gradient-based attribution vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4 ...
        || isempty(g)
    g = J' * delta;
end
d = size(J,2);
s = mdl.th_max(:) - mdl.th_min(:);
if numel(s) ~= d
    error(['      Error: sage_attribution: ' ...
        'parameter range length ' ...
        'does not match Jacobian dimension.']);
end
C = J .* delta;          % nxd matrix
C_scaled = C .* s.';     % nxd matrix

A_t = sum(abs(C_scaled),1).';
A_n = abs(s .* g);

end