function fit = fit_dual_kosugi_fdc(q,options)
%FIT_DUAL_KOSUGI_FDC Fit a five-parameter dual-Kosugi flow-duration curve.
%
%   fit = fit_dual_kosugi_fdc(q)
%   fit = fit_dual_kosugi_fdc(q,options)
%
%   The conditional positive-flow exceedance distribution is
%
%     U(q) = w U1(q;a1,b1) + (1-w) U2(q;a2,b2),
%
%   where Ui(q) = 0.5*erfc(log(q/ai)/(sqrt(2)*bi)).  The parameters obey
%   0 < a1 < a2, b1 > 0, b2 > 0, and 0 < w < 1.  The empirical zero-flow
%   probability p0 is stored separately. Only finite, nonnegative flows are
%   used; the continuous mixture is fitted to positive flows.
%
%   Options (all optional):
%     .zero_tolerance  observations <= this value are zeros [0]
%     .min_positive    minimum number of positive observations [20]
%     .n_starts        deterministic multistart count [24]
%     .huber_delta     pseudo-Huber probability scale [0.02]
%     .max_iter        maximum fminsearch iterations [2000]
%     .max_fun_evals   maximum function evaluations [6000]
%     .inverse_points  points in cached inverse lookup [4097]

%   The returned structure contains a1, b1, a2, b2, w, p0, diagnostics,
%   sample counts, success, exitflag, and a human-readable message.

if nargin < 2 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error('fit_dual_kosugi_fdc:Options', ...
        'options must be a scalar structure.');
end

zeroTolerance = local_option(options,'zero_tolerance',0);
minPositive = local_option(options,'min_positive',20);
nStarts = local_option(options,'n_starts',24);
delta = local_option(options,'huber_delta',0.02);
maxIter = local_option(options,'max_iter',2000);
maxFunEvals = local_option(options,'max_fun_evals',6000);
inversePoints = local_option(options,'inverse_points',4097);

if ~(isscalar(zeroTolerance) && isfinite(zeroTolerance) ...
        && zeroTolerance >= 0)
    error('fit_dual_kosugi_fdc:ZeroTolerance', ...
        'zero_tolerance must be a finite nonnegative scalar.');
end
if ~(isscalar(minPositive) && isfinite(minPositive) ...
        && minPositive >= 5 && minPositive == fix(minPositive))
    error('fit_dual_kosugi_fdc:MinimumPositive', ...
        'min_positive must be an integer of at least five.');
end
if ~(isscalar(nStarts) && isfinite(nStarts) ...
        && nStarts >= 1 && nStarts == fix(nStarts))
    error('fit_dual_kosugi_fdc:Starts', ...
        'n_starts must be a positive integer.');
end
if ~(isscalar(delta) && isfinite(delta) && delta > 0)
    error('fit_dual_kosugi_fdc:HuberDelta', ...
        'huber_delta must be a finite positive scalar.');
end
if ~(isscalar(maxIter) && isfinite(maxIter) && maxIter >= 1 ...
        && maxIter == fix(maxIter))
    error('fit_dual_kosugi_fdc:MaxIterations', ...
        'max_iter must be a positive integer.');
end
if ~(isscalar(maxFunEvals) && isfinite(maxFunEvals) ...
        && maxFunEvals >= 1 && maxFunEvals == fix(maxFunEvals))
    error('fit_dual_kosugi_fdc:MaxFunctionEvaluations', ...
        'max_fun_evals must be a positive integer.');
end
if ~(isscalar(inversePoints) && isfinite(inversePoints) ...
        && inversePoints >= 257 && inversePoints == fix(inversePoints))
    error('fit_dual_kosugi_fdc:InversePoints', ...
        'inverse_points must be an integer of at least 257.');
end

fit = local_empty_fit();
q = double(q(:));
q = q(isfinite(q) & q >= 0);
fit.n = numel(q);
if isempty(q)
    fit.message = "No finite nonnegative discharge observations.";
    return
end

isZero = q <= zeroTolerance;
fit.n_zero = sum(isZero);
fit.p0 = fit.n_zero/fit.n;
qp = sort(q(~isZero),'descend');
fit.n_positive = numel(qp);
if fit.n_positive < minPositive
    fit.message = sprintf(['Only %d positive discharge values; ' ...
        'at least %d are required.'],fit.n_positive,minPositive);
    return
end

n = numel(qp);
u = ((1:n)'-0.5)/n;
fit.qmin_positive = min(qp);
fit.q50_positive = median(qp);
fit.qmax_positive = max(qp);

starts = local_starts(qp,nStarts);
if exist('dual_kosugi_objective_mex','file') == 3
    objective = @(eta) dual_kosugi_objective_mex(eta,qp,u,delta);
else
    objective = @(eta) local_objective(eta,qp,u,delta);
end
settings = optimset('Display','off','MaxIter',maxIter, ...
    'MaxFunEvals',maxFunEvals,'TolX',1e-9,'TolFun',1e-12);
bestValue = Inf;
bestEta = starts(1,:);
bestExitflag = 0;
for s = 1:size(starts,1)
    [eta,value,exitflag] = fminsearch(objective,starts(s,:),settings);
    if isfinite(value) && value < bestValue
        bestValue = value;
        bestEta = eta;
        bestExitflag = exitflag;
    end
end

[a1,b1,a2,b2,w,valid] = local_decode(bestEta);
if ~valid || ~isfinite(bestValue)
    fit.message = "Dual-Kosugi optimization did not return a valid fit.";
    return
end
uhat = local_exceedance(qp,a1,b1,a2,b2,w);
residual = u-uhat;
fit.a1 = a1;
fit.b1 = b1;
fit.a2 = a2;
fit.b2 = b2;
fit.w = w;
fit.objective = bestValue;
fit.rmse_u = sqrt(mean(residual.^2));
fit.mae_u = mean(abs(residual));
fit.exitflag = bestExitflag;
fit.n_starts = size(starts,1);
fit.inverse_u = linspace(1e-6,1-1e-6,inversePoints).';
fit.inverse_q = local_inverse_cache(fit.inverse_u, ...
    a1,b1,a2,b2,w);
fit.success = all(isfinite([a1 b1 a2 b2 w bestValue])) ...
    && a1 > 0 && a2 > a1 && b1 > 0 && b2 > 0 ...
    && w > 0 && w < 1;
if fit.success
    fit.message = "Dual-Kosugi FDC fitted successfully.";
else
    fit.message = "Dual-Kosugi fit did not satisfy validity checks.";
end
end

function value = local_option(options,name,defaultValue)
if isfield(options,name) && ~isempty(options.(name))
    value = double(options.(name));
else
    value = defaultValue;
end
end

function fit = local_empty_fit()
fit = struct('a1',NaN,'b1',NaN,'a2',NaN,'b2',NaN,'w',NaN, ...
    'p0',NaN,'success',false,'exitflag',0,'objective',NaN, ...
    'rmse_u',NaN,'mae_u',NaN,'n_starts',0,'n',0, ...
    'n_positive',0,'n_zero',0,'qmin_positive',NaN, ...
    'q50_positive',NaN,'qmax_positive',NaN, ...
    'inverse_u',zeros(0,1),'inverse_q',zeros(0,1), ...
    'message',"Dual-Kosugi fit was not attempted.");
end

function starts = local_starts(q,nStarts)
quantilePairs = [10 60; 20 80; 30 70; 10 90; 25 75; 40 85];
bPairs = [0.20 0.80; 0.35 1.00; 0.60 1.40; 0.15 1.80];
weights = [0.25 0.50 0.75];
candidates = zeros(size(quantilePairs,1)*size(bPairs,1)*numel(weights),5);
row = 0;
for i = 1:size(quantilePairs,1)
    a = prctile(q,quantilePairs(i,:));
    a = sort(max(a,realmin('double')));
    if a(2) <= a(1)
        a(2) = a(1)*(1+1e-6);
    end
    for j = 1:size(bPairs,1)
        for h = 1:numel(weights)
            row = row+1;
            candidates(row,:) = local_encode(a(1),bPairs(j,1), ...
                a(2),bPairs(j,2),weights(h));
        end
    end
end
index = unique(round(linspace(1,size(candidates,1), ...
    min(nStarts,size(candidates,1)))));
starts = candidates(index,:);
end

function eta = local_encode(a1,b1,a2,b2,w)
w = min(max(w,1e-8),1-1e-8);
eta = [log(a1),log(b1),log(max(a2-a1,realmin('double'))), ...
    log(b2),log(w/(1-w))];
end

function [a1,b1,a2,b2,w,valid] = local_decode(eta)
eta = double(eta(:).');
valid = numel(eta) == 5 && all(isfinite(eta)) && all(abs(eta) < 40);
if ~valid
    [a1,b1,a2,b2,w] = deal(NaN);
    return
end
a1 = exp(eta(1));
b1 = exp(eta(2));
a2 = a1+exp(eta(3));
b2 = exp(eta(4));
w = 1/(1+exp(-eta(5)));
valid = all(isfinite([a1 b1 a2 b2 w])) ...
    && a1 > 0 && a2 > a1 && b1 > 0 && b2 > 0 ...
    && w > 0 && w < 1;
end

function value = local_objective(eta,q,u,delta)
[a1,b1,a2,b2,w,valid] = local_decode(eta);
if ~valid
    value = realmax('double')/1e100;
    return
end
uhat = local_exceedance(q,a1,b1,a2,b2,w);
if any(~isfinite(uhat))
    value = realmax('double')/1e100;
    return
end
residual = u-uhat;
value = mean(delta^2*(sqrt(1+(residual/delta).^2)-1));
end

function u = local_exceedance(q,a1,b1,a2,b2,w)
q = max(double(q),realmin('double'));
u = w*0.5.*erfc(log(q/a1)/(sqrt(2)*b1)) ...
    +(1-w)*0.5.*erfc(log(q/a2)/(sqrt(2)*b2));
end

function q = local_inverse_cache(u,a1,b1,a2,b2,w)
% Invert the monotone mixture once; model evaluations only interpolate.
tail = 1e-9;
z = sqrt(2)*erfcinv(2*tail);
qlo = min([a1*exp(-b1*z),a2*exp(-b2*z)]);
qhi = max([a1*exp( b1*z),a2*exp( b2*z)]);
qlo = max(qlo,realmin('double'));
qhi = max(qhi,qlo*(1+1e-6));
qgrid = logspace(log10(qlo),log10(qhi),max(32769,8*numel(u))).';
ugrid = local_exceedance(qgrid,a1,b1,a2,b2,w);
[us,keep] = unique(flipud(ugrid),'stable');
qs = flipud(qgrid);
qs = qs(keep);
q = interp1(us,qs,u,'pchip');
if any(~isfinite(q)) || any(q <= 0) || any(diff(q) >= 0)
    error('fit_dual_kosugi_fdc:InverseCache', ...
        'Unable to construct a finite strictly decreasing inverse cache.');
end
end
