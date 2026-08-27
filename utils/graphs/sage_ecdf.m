function [f,x] = sage_ecdf(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_ECDF Computes empirical cumulative distribution function values.
%
% SYNOPSIS:
%   [f,x] = sage_ecdf(z)
%
%   z           numeric vector
%   f           OUTPUT: ECDF values
%   x           OUTPUT: sorted support values of ECDF
%
% NOTES:
%   Input values are typically NSE values.
%   Non-finite filtering is usually handled before calling this routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Mar. 2026             %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z = z(:);
    z = z(isfinite(z) & isreal(z));
    if isempty(z)
        f = [];
        x = [];
        return
    end

    z = sort(z);
    n = numel(z);

    % unique support and counts (tie handling)
    [xu,~,ic] = unique(z,'sorted');
    cnt = accumarray(ic,1);
    fu = cumsum(cnt)/n;     % cdf at each unique value

    % prepend starting point at 0 (repeat xmin)
    x = [xu(1); xu];
    f = [0; fu];
end