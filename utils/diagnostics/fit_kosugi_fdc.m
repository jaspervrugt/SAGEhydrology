function fit = fit_kosugi_fdc(q,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIT_KOSUGI_FDC Robust fit of the 3-parameter Kosugi flow-duration curve
%
% SYNOPSIS:
%   fit = fit_kosugi_fdc(q)
%   fit = fit_kosugi_fdc(q,options)
%
% INPUT ARGUMENTS:
%   q                 observed discharge vector for ONE calibration period
%   options           optional settings structure
%    .zero_tol        values <= zero_tol are treated as zero [default: 0]
%    .min_positive    minimum number of positive observations [default: 20]
%    .n_c_grid        number of lower-shift starting values [default: 61]
%    .huber_delta     pseudo-Huber threshold in probability space [0.02]
%    .max_iter        fminsearch maximum iterations [default: 2000]
%    .max_fun_evals   fminsearch maximum evaluations [default: 5000]
%    .display         fminsearch display setting [default: 'off']
%
% OUTPUT ARGUMENT:
%   fit               structure with fitted parameters and diagnostics
%    .a               conditional positive-flow median parameter
%    .b               log-space spread parameter, b > 0
%    .c               lower shift/asymptote, c < min(q | q > zero_tol)
%    .p0              empirical zero-flow probability
%    .success         true if a finite Kosugi fit was obtained
%    .exitflag        fminsearch exit flag (0 if only initialization used)
%    .objective       final robust probability-space objective
%    .rmse_u          RMSE in exceedance-probability space
%    .mae_u           MAE in exceedance-probability space
%    .rmse_q          RMSE between empirical and fitted positive-flow FDC
%    .n               number of finite, nonnegative observations
%    .n_positive      number used to fit the continuous Kosugi component
%    .n_zero          number classified as zero flow
%    .n_negative      number of finite negative values excluded
%    .qmin_positive   minimum positive discharge used in the fit
%    .q50_positive    median positive discharge
%    .message         diagnostic text
%
% MODEL:
%   For positive flows, with exceedance probability u,
%
%       q(u) = c + (a-c) exp[b Phi^{-1}(1-u)]
%
%   or equivalently
%
%       U(q) = 0.5 erfc{ log[(q-c)/(a-c)] / (sqrt(2)b) }.
%
%   Exact zero flows are represented separately by
%
%       p0 = Pr(Q <= zero_tol).
%
% FITTING STRATEGY:
%   1. Remove nonfinite and negative observations.
%   2. Estimate p0 directly from the retained record.
%   3. Fit the continuous Kosugi distribution to positive flows only.
%   4. Search robustly over c. For every fixed c, estimate log(a-c) and b
%      from a robust IRLS regression of log(q-c) on Phi^{-1}(1-u).
%   5. Refine a,b,c simultaneously with fminsearch in unconstrained
%      coordinates. The transformation guarantees a>c, b>0 and c<qmin.
%   6. Minimize a smooth pseudo-Huber loss in exceedance-probability space.
%
% SAGE/SITE USE:
%   The optional preprocessing step is controlled in prep_stats with
%
%       loss.fdc.kosugi = true;
%
%   When enabled, prep_stats calls fit_kosugi_fdc using valid training-
%   period discharge only. Set loss.fdc.kosugi = false to deactivate fit.
%
% REFERENCE:
%   Sadegh, M., Vrugt, J.A., Gupta, H.V., and Xu, C. (2016), The soil
%   water characteristic as new class of closed-form parametric expressions
%   for the flow duration curve, Journal of Hydrology, 535, 438-456,
%   doi:10.1016/j.jhydrol.2016.01.027.
%
% NOTES:
%   * If p0 > 0, a is the median of the CONDITIONAL positive-flow Kosugi
%     component. p0 carries the discrete probability mass at zero.
%   * This routine uses only base MATLAB functions (no Optimization or
%     Statistics Toolbox functions are required).
%   * For SAGE/SITE, call this routine using TRAINING observations only.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Jasper A. Vrugt / OpenAI, Sep. 2026                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(options)
        options = struct();
    end
    if ~isstruct(options)
        error('fit_kosugi_fdc:Options', ...
            'options must be a structure.');
    end

    % ----------------
    % Default settings
    % ----------------
    zero_tol = local_get_option(options, ...
        'zero_tol',0);
    min_positive = local_get_option(options, ...
        'min_positive',20);
    n_c_grid = local_get_option(options, ...
        'n_c_grid',61);
    huber_delta = local_get_option(options, ...
        'huber_delta',0.02);
    max_iter = local_get_option(options, ...
        'max_iter',2000);
    max_fun_evals = local_get_option(options, ...
        'max_fun_evals',5000);
    display_opt = local_get_option(options, ...
        'display','off');

    if ~(isscalar(zero_tol) ...
            && isnumeric(zero_tol) ...
            && isfinite(zero_tol) ...
            && zero_tol >= 0)
        error('fit_kosugi_fdc:ZeroTolerance', ...
            ['options.zero_tol must be ' ...
            'a finite scalar >= 0.']);
    end
    if ~(isscalar(min_positive) ...
            && isnumeric(min_positive) ...
            && isfinite(min_positive) ...
            && min_positive >= 3)
        error('fit_kosugi_fdc:MinimumPositive', ...
            ['options.min_positive must be ' ...
            'a finite scalar >= 3.']);
    end
    min_positive = round(min_positive);
    if ~(isscalar(n_c_grid) ...
            && isnumeric(n_c_grid) ...
            && isfinite(n_c_grid) ...
            && n_c_grid >= 5)
        error('fit_kosugi_fdc:CGrid', ...
            ['options.n_c_grid must be ' ...
            'a finite scalar >= 5.']);
    end
    n_c_grid = round(n_c_grid);
    if ~(isscalar(huber_delta) ...
            && isnumeric(huber_delta) ...
            && isfinite(huber_delta) ...
            && huber_delta > 0)
        error('fit_kosugi_fdc:HuberDelta', ...
            ['options.huber_delta must be ' ...
            'a finite scalar > 0.']);
    end

    % -----------------------
    % Initialize output block
    % -----------------------
    fit = struct( ...
        'a',NaN, ...
        'b',NaN, ...
        'c',NaN, ...
        'p0',NaN, ...
        'success',false, ...
        'exitflag',0, ...
        'objective',NaN, ...
        'rmse_u',NaN, ...
        'mae_u',NaN, ...
        'rmse_q',NaN, ...
        'n',0, ...
        'n_positive',0, ...
        'n_zero',0, ...
        'n_negative',0, ...
        'qmin_positive',NaN, ...
        'q50_positive',NaN, ...
        'message',"");

    % ---------------------------------
    % Clean observations and estimate p0
    % ---------------------------------
    q = double(q(:));
    q = q(isfinite(q));
    if isempty(q)
        fit.message = "No finite " + ...
            "discharge observations.";
        return
    end

    fit.n_negative = sum(q < 0);
    q = q(q >= 0);
    fit.n = numel(q);
    if fit.n == 0
        fit.message = "No finite " + ...
            "nonnegative discharge observations.";
        return
    end

    is_zero = q <= zero_tol;
    fit.n_zero = sum(is_zero);
    fit.p0 = fit.n_zero / fit.n;

    qp = q(~is_zero);
    fit.n_positive = numel(qp);
    if fit.n_positive < min_positive
        fit.message = sprintf(['Only %d positive ' ...
            'discharge values; ' ...
            'at least %d are required.'], ...
            fit.n_positive,min_positive);
        return
    end

    % --------------------------------------------------------
    % Empirical positive-flow FDC and Gaussian latent coordinate
    % --------------------------------------------------------
    qp = sort(qp,'descend');
    n = numel(qp);
    u = ((1:n)' - 0.5) / n;                 % exceedance plotting position
    z = sqrt(2) * erfcinv(2*u);             % Phi^{-1}(1-u)

    qmin = min(qp);
    qmax = max(qp);
    q50 = median(qp);
    fit.qmin_positive = qmin;
    fit.q50_positive = q50;

    qrange = qmax - qmin;
    qscale = max([q50, qrange, qmax/10, realmin('double')]);
    if ~isfinite(qscale) ...
            || qscale <= 0
        fit.message = "Positive discharge record " + ...
            "has invalid scale.";
        return
    end

    % ------------------------------------------------------------------
    % Stage 1: robust grid search over c.  For fixed c,
    %          log(q-c) = beta0 + b*z, with a = c + exp(beta0).
    % ------------------------------------------------------------------
    dmin = max(1e-8*qscale,100*eps(qscale));
    dmax = max([10*qscale, qmin + 10*qscale, 100*dmin]);
    dgrid = logspace(log10(dmin),log10(dmax),n_c_grid)';

    % Include c = 0 explicitly whenever admissible.
    if qmin > 0
        d0 = qmin;
        if d0 >= dmin ...
                && d0 <= dmax
            dgrid = unique([dgrid; d0]);
        end
    end

    best_obj = Inf;
    best_abc = [NaN NaN NaN];

    for ii = 1:numel(dgrid)
        c_try = qmin - dgrid(ii);
        y = log(qp - c_try);

        [beta0,b_try,ok] = local_robust_line(z,y);
        if ~ok ...
                || ~(isfinite(b_try) ...
                && b_try > 0)
            continue
        end

        ac = exp(beta0);
        a_try = c_try + ac;
        if ~(isfinite(a_try) ...
                && a_try > c_try)
            continue
        end

        obj = local_objective_abc(a_try,b_try,c_try, ...
            qp,u,huber_delta);
        if obj < best_obj
            best_obj = obj;
            best_abc = [a_try,b_try,c_try];
        end
    end

    if ~all(isfinite(best_abc))
        fit.message = "Unable to construct " + ...
            "a valid Kosugi starting point.";
        return
    end

    % -----------------------------------------------------------
    % Stage 2: unconstrained nonlinear refinement.
    %
    %   c = qmin - qscale exp(eta3)
    %   a = c + qscale exp(eta1)
    %   b = exp(eta2)
    %
    % These transformations enforce all model-domain constraints.
    % -----------------------------------------------------------
    a0 = best_abc(1);
    b0 = best_abc(2);
    c0 = best_abc(3);
    eta0 = [ ...
        log(max((a0-c0)/qscale,realmin('double'))); ...
        log(max(b0,realmin('double'))); ...
        log(max((qmin-c0)/qscale,realmin('double')))];

    objective = @(eta) local_objective_eta(eta, ...
        qp,u,qmin,qscale,huber_delta);

    optim_options = optimset( ...
        'Display',char(string(display_opt)), ...
        'MaxIter',max_iter, ...
        'MaxFunEvals',max_fun_evals, ...
        'TolX',1e-9, ...
        'TolFun',1e-12);

    try
        [eta_star,obj_star,exitflag] = fminsearch( ...
            objective,eta0,optim_options);
    catch ME
        eta_star = eta0;
        obj_star = best_obj;
        exitflag = -99;
        fit.message = "fminsearch failed: " + string(ME.message);
    end

    [a_star,b_star,c_star,valid_star] = local_eta_to_abc( ...
        eta_star,qmin,qscale);

    % Keep the robust grid solution if refinement failed or became worse.
    if ~valid_star ...
            || ~isfinite(obj_star) ...
            || obj_star > best_obj
        a_star = a0;
        b_star = b0;
        c_star = c0;
        obj_star = best_obj;
        if exitflag > 0
            exitflag = 0;
        end
    end

    % ---------------
    % Fit diagnostics
    % ---------------
    log_ratio = log((qp-c_star)/(a_star-c_star));
    uhat = 0.5 * erfc(log_ratio/(sqrt(2)*b_star));
    qhat = c_star + (a_star-c_star) .* exp(b_star*z);

    fit.a = a_star;
    fit.b = b_star;
    fit.c = c_star;
    fit.exitflag = exitflag;
    fit.objective = obj_star;
    fit.rmse_u = sqrt(mean((u-uhat).^2));
    fit.mae_u = mean(abs(u-uhat));
    fit.rmse_q = sqrt(mean((qp-qhat).^2));
    fit.success = all(isfinite([fit.a fit.b fit.c fit.objective])) ...
        && fit.b > 0 ...
        && fit.a > fit.c ...
        && fit.c < qmin;

    if fit.success ...
            && strlength(fit.message) == 0
        if fit.p0 > 0
            fit.message = sprintf([ ...
                'Kosugi positive-flow ' ...
                'component fitted successfully; ' ...
                'p0 = %.6g.'],fit.p0);
        else
            fit.message = "Kosugi FDC " + ...
                "fitted successfully.";
        end
    elseif ~fit.success ...
            && strlength(fit.message) == 0
        fit.message = "Kosugi fit did not " + ...
            "produce a valid parameter set.";
    end

end

function value = local_get_option(S,name,defaultValue)
%LOCAL_GET_OPTION Return option value or default.
    if isfield(S,name) ...
            && ~isempty(S.(name))
        value = S.(name);
    else
        value = defaultValue;
    end
end

function [beta0,b,ok] = local_robust_line(x,y)
%LOCAL_ROBUST_LINE Robust IRLS fit y = beta0 + b*x.

    x = double(x(:));
    y = double(y(:));
    X = [ones(size(x)),x];
    ok = false;
    beta0 = NaN;
    b = NaN;

    if numel(x) < 3 ...
            || any(~isfinite(x)) ...
            || any(~isfinite(y))
        return
    end

    beta = X\y;
    if any(~isfinite(beta))
        return
    end

    for iter = 1:30
        r = y - X*beta;
        s = 1.482602218505602 * median(abs(r-median(r)));
        s_floor = max(1e-10,1e-10*max(1,median(abs(y))));
        if ~isfinite(s) ...
                || s < s_floor
            s = s_floor;
        end

        % Smooth pseudo-Huber IRLS weights; tuning ~= 1.345 robust sigma.
        t = r/(1.345*s);
        w = 1 ./ sqrt(1+t.^2);
        sw = sqrt(w);
        Xw = X .* sw;
        yw = y .* sw;
        beta_new = Xw\yw;

        if any(~isfinite(beta_new))
            return
        end
        if norm(beta_new-beta) <= 1e-10*(1+norm(beta))
            beta = beta_new;
            break
        end
        beta = beta_new;
    end

    beta0 = beta(1);
    b = beta(2);
    ok = isfinite(beta0) ...
        && isfinite(b);
end

function obj = local_objective_eta(eta,q,u,qmin,qscale,delta)
%LOCAL_OBJECTIVE_ETA Robust objective in unconstrained coordinates.

    [a,b,c,valid] = local_eta_to_abc(eta,qmin,qscale);
    if ~valid
        obj = realmax('double')/1e100;
        return
    end
    obj = local_objective_abc(a,b,c,q,u,delta);
end

function [a,b,c,valid] = local_eta_to_abc(eta,qmin,qscale)
%LOCAL_ETA_TO_ABC Map unconstrained coordinates to valid Kosugi parameters.

    eta = double(eta(:));
    valid = numel(eta)==3 ...
        && all(isfinite(eta)) ...
        && all(abs(eta) <= 40);
    if ~valid
        a = NaN; b = NaN; c = NaN;
        return
    end

    c = qmin - qscale*exp(eta(3));
    a = c + qscale*exp(eta(1));
    b = exp(eta(2));
    valid = all(isfinite([a b c])) ...
        && b > 0 ...
        && a > c ...
        && c < qmin;
end

function obj = local_objective_abc(a,b,c,q,u,delta)
%LOCAL_OBJECTIVE_ABC Smooth robust loss in exceedance-probability space.

    if ~(isfinite(a) ...
            && isfinite(b) ...
            && isfinite(c) ...
            && b > 0 ...
            && a > c ...
            && all(q > c))
        obj = Inf;
        return
    end

    log_ratio = log((q-c)/(a-c));
    uhat = 0.5 * erfc(log_ratio/(sqrt(2)*b));
    if any(~isfinite(uhat))
        obj = Inf;
        return
    end

    r = u - uhat;
    rho = delta^2 * (sqrt(1+(r/delta).^2)-1);
    obj = mean(rho);
end
