function [stats,fdc]=paired_observation_stats(y)
%PAIRED_OBSERVATION_STATS Rebuild caches after nonfinite simulations/J rows.
% Observation-only formulas are shared with the PREP_STATS implementation.
y=double(y(:));
stats=struct('mean',mean(y),'std',std(y),'TSS',sum((y-mean(y)).^2), ...
    'huber_scale',local_huber_scale(y));
fdc=local_fdc_cache(y);
end

function S_y = local_huber_scale(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL_HUBER_SCALE Compute robust observation scale for Huber loss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Remove missing and nonfinite observations
    y = double(y(:));
    y = y(isfinite(y));

    if isempty(y)
        S_y = NaN;
        return
    end

    % Normal-consistency factor:
    % 1 / norminv(0.75)
    kappa = 1.482602218505602;

    % Median absolute deviation about the median
    T_y = median(y);
    S_y = kappa * median(abs(y - T_y));

    % Fallback for ephemeral or zero-flow-dominated basins
    if ~isfinite(S_y) || ...
        S_y <= 0

        yp = y(y > 0);
        if ~isempty(yp)
            S_y = median(yp);
        end
    end

    % Numerical floor relative to observed discharge magnitude
    mean_abs_y = mean(abs(y));

    if ~isfinite(mean_abs_y)
        mean_abs_y = 1;
    end

    S_floor = max(1e-6, ...
        1e-6 * mean_abs_y);

    if ~isfinite(S_y) || ...
        S_y <= S_floor
        S_y = S_floor;
    end

end

function fdc = local_empty_fdc_cache()
%LOCAL_EMPTY_FDC_CACHE Initialize observed FDC cache.

    fdc = struct( ...
        'ys',[], ...
        'Py',[], ...
        'S_yy',NaN, ...
        'q0',NaN, ...
        'log_ys',[], ...
        'n',0, ...
        'D0',NaN, ...
        'D0_p',NaN, ...
        'D0_logp',NaN);

end

function fdc = local_fdc_cache(y)
%LOCAL_FDC_CACHE Precompute observation-only terms for FDC distance.
%
% The information-theoretic FDC distance contains
%
%   S_qy = sum_{i,j} |q_i - y_j|
%   S_qq = sum_{i,j} |q_i - q_j|
%   S_yy = sum_{i,j} |y_i - y_j|
%
% Because observed discharge y does not change during optimization,
% its sorted values, prefix sums and S_yy term are cached once.

    y = double(y(:));
    y = y(isfinite(y));

    n = numel(y);

    if n == 0
        fdc = local_empty_fdc_cache();
        return
    end

    % Sorted observed discharge
    ys = sort(y);

    % Prefix sums used by the O(n) cross-term calculation
    Py = [0; cumsum(ys)];

    % S_yy = sum_{i,j} |y_i - y_j|
    w = 2*(1:n)' - n - 1;
    S_yy = 2 * sum(w .* ys);

    % FDC divergence of the constant-median reference.
    % This defines the basin-specific zero-skill benchmark for S_fdc.
    ymed = median(ys);
    % Median-flow reference divergence with units of discharge.
    D0 = fdc_cdf_distance_sorted(repmat(ymed,n,1),ys);
    % Numerical safeguard only
    D0 = max(D0,0);

    % The probability-space formulations use the same constant-median
    % reference. The logarithmic offset equals one percent of the median
    % positive observed discharge, as documented in Appendix C.
    D0_p = mean((ymed - ys).^2);
    yp = ys(ys > 0);
    if isempty(yp)
        q0 = NaN;
        log_ys = nan(size(ys));
        D0_logp = NaN;
    else
        q0 = 0.01*median(yp);
        log_ys = log(ys + q0);
        D0_logp = mean((log(ymed + q0) - log_ys).^2);
    end

    fdc = struct( ...
        'ys',ys, ...
        'Py',Py, ...
        'S_yy',S_yy, ...
        'q0',q0, ...
        'log_ys',log_ys, ...
        'n',n, ...
        'D0',D0, ...
        'D0_p',D0_p, ...
        'D0_logp',D0_logp);
end

