function x = get_geo_column(T,cands,lo,hi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET_GEO_COLUMN Return the first valid geographic coordinate column.
%
% SYNOPSIS:
%   x = get_geo_column(T,cands,lo,hi)
%
% INPUT:
%   T       Table containing candidate geographic-coordinate columns
%   cands   Cell/string array with candidate variable names, test in order
%   lo      Lower admissible coordinate bound
%   hi      Upper admissible coordinate bound
%
% OUTPUT:
%   x       Numeric column vector with the first candidate column whose
%           finite values fall predominantly within [lo,hi]. If no valid
%           candidate is found, x is returned as NaN(height(T),1).
%
% DESCRIPTION:
%   This helper searches a table for candidate latitude, longitude, or
%   related geographic coordinate variables. Candidate columns are converted
%   to numeric values and screened against the supplied admissible range.
%   The first column passing the range check is returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, June 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = get_numeric_column(T,cands);
    bad = ~isfinite(x) | x < lo | x > hi;
    x(bad) = nan;

end
