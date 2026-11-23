function [D_fdc,D_p,D_logp] = fdc_metrics_cached(q_n,fdc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FDC_METRICS_CACHED Three FDC divergences using one simulated-flow sort.
%
% SYNOPSIS:
%   [D_fdc,D_p,D_logp] = fdc_metrics_cached(q_n,fdc)
%
% INPUT:
%   q_n       nx1 simulated discharge vector
%   fdc       observed-record cache prepared by prep_stats
%
% OUTPUT:
%   D_fdc     physical-discharge CDF divergence
%   D_p       quantile-function divergence in probability space
%   D_logp    log-quantile divergence in probability space
%
% DESCRIPTION:
%   The simulated discharge is sorted only once. The three formulations
%   therefore share their dominant O(n log n) operation; D_p and D_logp
%   add only O(n) vector arithmetic to the existing D_fdc calculation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Sept. 2026                                %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = fdc.n;
    [D_fdc,D_p,D_logp] = deal(NaN);

    if n == 0 || ...
        isempty(q_n)
        return
    end

    q = double(q_n(:));
    if numel(q) ~= n
        error('fdc_metrics_cached:LengthMismatch', ...
            ['Simulated and observed FDC records must have equal ' ...
                'length.']);
    end
    if any(~isfinite(q))
        return
    end

    qs = sort(q);
    ys = double(fdc.ys(:));

    % Equation 20: squared CDF difference integrated over discharge.
    D_fdc = fdc_cdf_distance_sorted(qs,ys);
    D_fdc = max(D_fdc,0);               % roundoff safeguard at equality

    % Probability-space squared quantile difference.
    dq = qs - ys;
    D_p = mean(dq.^2);

    % Probability-space squared log-quantile difference.
    q0 = double(fdc.q0);
    if isfinite(q0) && ...
        q0 > 0 && ...
        all(qs + q0 > 0)
        if isfield(fdc,'log_ys') && ...
            numel(fdc.log_ys) == n
            log_ys = double(fdc.log_ys(:));
        else
            log_ys = log(ys + q0);
        end
        dlq = log(qs + q0) - log_ys;
        D_logp = mean(dlq.^2);
    end
end
