function q = sage_prctile(x,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_PRCTILE Compute vector percentiles without a toolbox.
%
%  Match the midpoint interpolation used by MATLAB PRCTILE for a finite
%  numeric vector. This supports older MATLAB releases without Statistics
%  and Machine Learning Toolbox.
%
% SYNOPSIS:
%   q = sage_prctile(x,p)
%
% INPUT ARGUMENTS:
%   x               finite numeric data vector
%   p               percentiles in [0,100]
%
% OUTPUT ARGUMENTS:
%   q               percentiles, with the same shape as p
%
% NOTES:
%   Sorted value k is located at percentile 100*(k-0.5)/n. Results beyond
%   the first and last locations are clamped to the data minimum/maximum.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isnumeric(x) ...
            || (~isvector(x) && ~isempty(x)) ...
            || any(~isfinite(x(:)))
        error('sage_prctile:InvalidData', ...
            'x must be a finite numeric vector.');
    end
    if ~isnumeric(p) ...
            || any(~isfinite(p(:))) ...
            || any(p(:) < 0 | p(:) > 100)
        error('sage_prctile:InvalidPercentiles', ...
            'p must contain percentiles between 0 and 100.');
    end

    x = sort(x(:));
    if isempty(x)
        q = nan(size(p));
        return
    end

    n = numel(x);
    location = min(n,max(1,n * (p(:)/100) + 0.5));
    lower = floor(location);
    upper = ceil(location);
    weight = location - lower;
    values = (1-weight) .* x(lower) + weight .* x(upper);
    q = reshape(values,size(p));

end
