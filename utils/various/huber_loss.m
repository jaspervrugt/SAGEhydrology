function [loss,delta] = huber_loss(res,S_y,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HUBER_LOSS Compute Huber loss using precomputed observation scale
%
% SYNOPSIS: [loss,delta] = huber_loss(res,S_y,c)
%   res     nx1 residual vector, res = y - q
%   S_y     scalar precomputed robust observation scale
%   c       scalar Huber threshold; default = 1.345
%   loss    OUTPUT: scalar Huber loss
%   delta   OUTPUT: nx1 loss-sensitivity vector
%
% DESCRIPTION:
%   This function evaluates the Huber loss, which combines quadratic and
%   linear penalties to provide robustness against outliers. Residuals
%   smaller than a threshold are treated with a squared penalty, while
%   larger residuals receive a linear penalty. This function also returns
%   the loss-sensitivity vector as this requires only 1 additional function
%   evaluation.
%
%   The residuals are normalized using the precomputed robust observation
%   scale S_y. This scale is prepared once from the observed discharge
%   record in prep_stats, thereby avoiding repeated median and MAD
%   calculations during model optimization.
%
%   Reference:
%     Huber, P. J. (1964), Robust estimation of a location parameter,
%     Annals of Mathematical Statistics, 35(1), 73–101.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delta = zeros(size(res));

    if isempty(res) ...
            || ~isscalar(S_y) ...
            || ~isfinite(S_y) ...
            || S_y <= 0
        loss = NaN;
        return
    end

    if nargin < 3 ...
            || isempty(c)
        c = 1.345;
    elseif ~(isscalar(c) ...
            && isnumeric(c) ...
            && isfinite(c) ...
            && c > 0)
        error(['      Error:huber_loss: ' ...
            'c must be a positive finite scalar.']);
    end

    good = isfinite(res);

    if ~any(good)
        loss = NaN;
        return
    end

    u = res(good) ./ S_y;
    au = abs(u);

    % min(|u|,c)
    v = min(au,c);

    % Equivalent to the piecewise Huber loss:
    % 0.5*u^2                     if |u| <= c
    % c*|u| - 0.5*c^2             if |u| > c
    loss = sum(v .* (au - 0.5*v));

    % res = y - q, therefore derivative with respect to q
    delta(good) = -(sign(u) .* v) ./ S_y;

end