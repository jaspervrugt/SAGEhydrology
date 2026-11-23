function value=fdc_cdf_distance_sorted(qs,ys)
% FDC_CDF_DISTANCE_SORTED Stable integral of squared empirical CDF
% difference.
% Inputs have equal length; ties have zero integration width. This avoids
% subtracting large nearly equal pair-distance sums at a perfect/near fit.
qs=double(qs(:));ys=double(ys(:));n=numel(ys);
if n==0 || ...
    numel(qs)~=n || ...
    any(~isfinite(qs)) || ...
    any(~isfinite(ys))
    value=NaN;return
end
[x,order]=sort([qs;ys]);
balance=cumsum(2*double(order<=n)-1)/n;
value=sum(diff(x).*balance(1:end-1).^2);
end
