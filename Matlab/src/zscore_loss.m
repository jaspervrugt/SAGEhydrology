function [value,delta] = zscore_loss(y,sim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZSCORE_LOSS Population z-score MSE and its complete trajectory derivative.
%  [value,delta] = zscore_loss(y,sim) returns mean((z(sim)-z(y)).^2).
%  The population scales (STD(...,1)) give value = 2*(1-corr(y,sim)).
%  Inputs contain paired valid samples from the scored period only.
%  delta includes derivatives of the simulated mean AND standard deviation.
%  Physical model outputs and their Jacobians remain unchanged.
%
% Copyright (c) 2026 Jasper A. Vrugt, University of California, Irvine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = double(y(:));
    sim = double(sim(:));
    n = numel(y);
    if n < 2 || ...
        numel(sim) ~= n || ...
        any(~isfinite(y)) || ...
        any(~isfinite(sim))
        error('SAGE:ZscoreSamples','Need two finite paired samples.');
    end
    yc=(y-y(1))-mean(y-y(1));
    sy = sqrt(mean(yc.^2));
    centered = (sim-sim(1))-mean(sim-sim(1));
    ss = sqrt(mean(centered.^2));
    if ~(sy > 0) || ...
        ~(ss > 0) || ...
        ~isfinite(sy) || ...
        ~isfinite(ss)
        error('SAGE:ZscoreConstant', ...
            ['Z-score loss requires nonzero observed/simulated ' ...
                'variability.']);
    end
    zs = centered/ss;
    zy = yc/sy;
    residual = zs-zy;
    value = mean(residual.^2);
    delta = 2/(n*ss)*(residual-mean(residual) ...
        - zs*mean(residual.*zs));
end
