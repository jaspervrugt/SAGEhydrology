function Fval = ecdf_value_from_stairs(xs,fs,x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ECDF_VALUE_FROM_STAIRS Returns ECDF value at a specified x-location.
%
% SYNOPSIS:
%   Fval = ecdf_value_from_stairs(xs,fs,x0)
%
%   xs          x-values of stair-step ECDF
%   fs          y-values of stair-step ECDF
%   x0          query x-location
%   Fval        OUTPUT: ECDF value corresponding to x0
%
% NOTES:
%   This helper is used to locate the y-position of the median marker on
%   the stair-step ECDF curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Mar. 2026             %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xs,fs are knots as produced by ecdf_to_stairs_fixed
xs = xs(:); fs = fs(:);
if numel(xs) ~= numel(fs)
    error('ecdf_value_from_stairs:SizeMismatch', ...
        'xs and fs must be same length.');
end
x0 = double(x0);

% Right-continuous: find last xs <= x0
idx = find(xs <= x0,1,'last');
if isempty(idx)
    Fval = 0;
else
    Fval = fs(idx);
end
end
