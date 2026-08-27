function [xMark,yMark,x1m,x2m,xt,ha,mode] = ...
    median_marker_geom(medX,medF,xL,xR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MEDIAN_MARKER_GEOM Returns geometry of ECDF median annotation
%
% SYNOPSIS:
%   [xMark,yMark,x1m,x2m,xt,ha,mode] = ...
%       median_marker_geom(medX,medF,xL,xR)
%   medX        x-coordinate ECDF median
%   medF        ECDF value @ medX
%   xL          left x-limit ECDF panel
%   xR          right x-limit ECDF panel
%   xMark       OUTPUT: x-location square median marker
%   yMark       OUTPUT: y-location square median marker
%   x1m         OUTPUT: left x-location short horizontal marker
%   x2m         OUTPUT: right x-location short horizontal marker
%   xt          OUTPUT: x-location median text label
%   ha          OUTPUT: horizontal alignment median text label
%   mode        OUTPUT: structure with placement flags
%    .specialHalfAtLeft   logical flag for left-edge handling
%    .side                side of text label --> 'left' | 'right'
%
% NOTES:
%   1. Text placement rule is simple:
%        medX < 0  --> text right of marker
%        medX >= 0 --> text left of marker
%   2. Special handling remains for medians left of plotting window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode = struct('specialHalfAtLeft', ...
    false,'side','left');

% --------------
% Basic settings
% --------------
span = max(eps,xR - xL);
seg = 0.05 * span;      % short horizontal length
pad = 0.01 * span;      % text padding from segment end
inset = 0.0025 * span;  % keep off y-axis visually

% -----------------------------------
% Special case: median left of window
% -----------------------------------
if (medX <= xL) ...
        || ((medF >= 0.5) ...
        && (medX <= xL + inset))
    mode.specialHalfAtLeft = true;
    mode.side = 'right';

    xMark = xL;
    yMark = 0.5;

    x1m = xL + inset;
    x2m = min(xR,x1m + seg);

    xt = min(xR,x2m + pad);
    ha = 'left';
    return
end

% --------------------------
% Normal case: median inside
% --------------------------
xMark = min(max(medX,xL),xR);
yMark = min(max(medF,0),1);

% ------------------
% Text side rule:
% medX < 0 --> right
% medX >= 0 --> left
% ------------------
if medX < 0
    mode.side = 'right';

    x1m = xMark;
    x2m = min(xR - inset,xMark + seg);
    xt = min(xR - inset,x2m + pad);
    ha = 'left';
else
    mode.side = 'left';

    x2m = xMark;
    x1m = max(xL + inset,xMark - seg);
    xt = max(xL + inset,x1m - pad);
    ha = 'right';
end

end