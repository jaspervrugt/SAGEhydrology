function q = sage_quantile(x,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_QUANTILE Compute vector quantiles without a toolbox.
%
%  Use the midpoint interpolation of MATLAB QUANTILE for a finite numeric
%  vector. This also works before QUANTILE became part of base MATLAB.
%
% SYNOPSIS:
%   q = sage_quantile(x,p)
%
% INPUT ARGUMENTS:
%   x               finite numeric data vector
%   p               probabilities in [0,1]
%
% OUTPUT ARGUMENTS:
%   q               quantiles, with the same shape as p
%
% NOTES:
%   SAGE_PRCTILE implements the same midpoint rule on the 0-to-100 scale.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isnumeric(p) ...
            || any(~isfinite(p(:))) ...
            || any(p(:) < 0 | p(:) > 1)
        error('sage_quantile:InvalidProbabilities', ...
            'p must contain probabilities between 0 and 1.');
    end

    q = sage_prctile(x,100*p);

end
