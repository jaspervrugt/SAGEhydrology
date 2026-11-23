function q = sage_percentile_linear(x,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_PERCENTILE_LINEAR NumPy linear/type-7 finite-vector percentiles.
%  Separate from SAGE_PRCTILE, which retains MATLAB midpoint interpolation.
%  Copyright (c) 2026 Jasper A. Vrugt, UCI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = sort(double(x(:)));
    if isempty(x) || ...
        any(~isfinite(x)) || ...
        any(~isfinite(p(:)) | ...
        p(:) < 0 | ...
        p(:) > 100)
        error('attributes:Percentile','Invalid percentile input.');
    end
    t = (numel(x)-1)*p/100+1;
    low = floor(t); high = ceil(t);
    q = reshape(x(low(:)).*(1-t(:)+low(:)) ...
        + x(high(:)).*(t(:)-low(:)),size(p));
end
