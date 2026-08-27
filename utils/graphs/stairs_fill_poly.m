function [xp,yp] = stairs_fill_poly(xs,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STAIRS_FILL_POLY Builds polygon coordinates for filled ECDF area.
%
% SYNOPSIS:
%   [xp,yp] = stairs_fill_poly(xs,fs)
%
%   xs          x-values of stair-step ECDF
%   fs          y-values of stair-step ECDF
%   xp          OUTPUT: x-coordinates of fill polygon
%   yp          OUTPUT: y-coordinates of fill polygon
%
% NOTES:
%   The returned polygon can be used directly with PATCH to fill the area
%   under the ECDF curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Mar. 2026             %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given stairs knots xs,fs (same length), produce (xp,yp) polygon coordinates
% that fill from y=0 up to the stairs curve.
xs = xs(:); fs = fs(:);
if numel(xs) ~= numel(fs)
    error('stairs_fill_poly:SizeMismatch','xs and fs must be same length.');
end
n = numel(xs);
if n < 2
    xp = xs;
    yp = zeros(size(xs));
    return
end

% Construct stairs vertices for patch:
% horizontal segments between xs(i)->xs(i+1) at height fs(i)
xv = zeros(2*(n-1),1);
yv = zeros(2*(n-1),1);
k = 1;
for i=1:n-1
    xv(k)   = xs(i);   yv(k)   = fs(i); k = k+1;
    xv(k)   = xs(i+1); yv(k)   = fs(i); k = k+1;
end

% Close polygon down to baseline y=0
xp = [xv; xs(end); xs(1)];
yp = [yv; 0; 0];
end
