function [A_t,A_n] = sage_attribution(J,delta,mdl,g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_ATTRIBUTION Compute total and net parameter attribution.
%
%  Scale time-resolved Jacobian contributions by each parameter's physical
%  range and aggregate them across time steps.
%
% SYNOPSIS:
%   [A_t,A_n] = sage_attribution(J,delta,mdl,g)
%   [A_t,A_n] = sage_attribution(J,delta,mdl)
%
% INPUT ARGUMENTS:
%   J               n-by-d discharge Jacobian
%   delta           n-by-1 loss sensitivity with respect to discharge
%   mdl             physical parameter bounds
%    .th_min         d-by-1 lower bounds
%    .th_max         d-by-1 upper bounds
%   g               optional d-by-1 loss gradient; default J' * delta
%
% OUTPUT ARGUMENTS:
%   A_t             d-by-1 sum of absolute scaled contributions
%   A_n             d-by-1 absolute scaled net gradient
%
% NOTES:
%   Parameter-range scaling uses mdl.th_max - mdl.th_min.
%
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
