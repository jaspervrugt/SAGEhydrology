function [xs,fs] = ecdf_to_stairs_fixed(x,f,xL,xR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ECDF_TO_STAIRS_FIXED Converts ECDF support to stair-step plotting arrays.
%
% SYNOPSIS:
%   [xs,fs] = ecdf_to_stairs_fixed(x,f,xL,xR)
%
%   x           ECDF support values
%   f           ECDF values
%   xL          left x-limit of plot
%   xR          right x-limit of plot
%   xs          OUTPUT: x-values of stair-step representation
%   fs          OUTPUT: y-values of stair-step representation
%
% NOTES:
%   The returned arrays are suitable for line plotting and for building
%   filled polygons under the ECDF curve over a fixed x-range [xL,xR].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Mar. 2026             %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = x(:); f = f(:);

% Safety: ensure sorted
[x,ii] = sort(x,'ascend');
f = f(ii);

% Remove NaNs/Infs
ok = isfinite(x) & isfinite(f);
x = x(ok); f = f(ok);

if isempty(x)
    xs = [xL; xR];
    fs = [0; 0];
    return
end

% Clamp domain to [xL,xR] by adding boundary points
% Compute ECDF value at boundaries using right-continuous convention
fL = sage_ecdf_value_from_points(x,f,xL);
fR = sage_ecdf_value_from_points(x,f,xR);

% Keep only points within range
in = (x >= xL) & (x <= xR);
x_in = x(in);
f_in = f(in);

% Build stairs knots:
% Start at xL with f(xL)
xs = [xL; x_in; xR];
fs = [fL; f_in; fR];

% Ensure monotone
fs = max(fs, 0);
fs = min(fs, 1);

end

% -------------------------------------------------------------------------
function y = sage_ecdf_value_from_points(x,f,xq)
% Right-continuous ECDF value at xq given discrete ecdf points (x,f)
% ecdf returns f(i)=P(X<=x(i)). For xq between points, use last f with x<=xq.
idx = find(x <= xq, 1, 'last');
if isempty(idx)
    y = 0;
else
    y = f(idx);
end
end
