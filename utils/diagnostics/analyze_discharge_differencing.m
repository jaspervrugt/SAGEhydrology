function result = analyze_discharge_differencing(Q,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYZE_DISCHARGE_DIFFERENCING Estimate hourly discharge random errors.
%
% SYNOPSIS: result = analyze_discharge_differencing(Q,options)
%   Q         nx1 hourly discharge vector in consistent units
%   options   OPTIONAL analysis settings
%    .order          differencing order k (default: 3)
%    .half_window    observations on either side (default: 100)
%    .method         'moving' or 'all' (default: 'moving')
%    .min_discharge  smallest retained discharge (default: 0)
%    .min_sigma      smallest retained local sigma (default: 0)
%    .tolerance      paper-compatible alias for min_sigma
%    .samples_per_day samples in a complete day (default: [])
%   result    structure with raw, smoothed, and fitted diagnostics
%
% DESCRIPTION:
%   Implements the difference-based variance estimator of de Oliveira and
%   Vrugt (2022), Water Resources Research, doi:10.1029/2022WR032263:
%
%       sigma2(t) = diff(Q(t),k)^2 / nchoosek(2*k,k)
%
%   The valid pairs are sorted by discharge. For method 'moving', local
%   variances are averaged over a centered window containing up to
%   2*half_window+1 observations. The paper-consistent regression is
%
%       sigma = alpha*Q + beta*mean(Q).
%
%   An exploratory linear variance fit is returned separately. Missing or
%   invalid observations invalidate the complete k+1 observation window.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(options)
        options = struct();
    end
    options = local_options(options);

    if ~isnumeric(Q) ...
            || ~isvector(Q)
        error('analyze_discharge_differencing:InvalidDischarge', ...
            'Q must be a numeric vector.');
    end
    Q = double(Q(:));
    n = numel(Q);
    k = options.order;
    if n <= k
        error('analyze_discharge_differencing:RecordTooShort', ...
            'Q must contain more than %d observations.',k);
    end

    invalid = ~isfinite(Q) | Q < options.min_discharge;
    Qwork = Q;
    Qwork(invalid) = NaN;
    delta = diff(Qwork,k);
    normalization = nchoosek(2*k,k);
    sigma2 = delta.^2/normalization;
    dischargeWindow = nan(n-k,k+1);
    for j = 0:k
        dischargeWindow(:,j+1) = Qwork(k+1-j:n-j);
    end
    discharge = mean(dischargeWindow,2);
    sigma = sqrt(sigma2);

    valid = isfinite(discharge) ...
        & isfinite(sigma2) ...
        & sigma >= options.min_sigma;
    discharge = discharge(valid);
    sigma2 = sigma2(valid);
    sigma = sigma(valid);
    sourceIndex = find(valid) + k;

    if numel(discharge) < 3
        error('analyze_discharge_differencing:TooFewValidPairs', ...
            'At least three valid differenced discharge pairs are required.');
    end

    [dischargeSorted,order] = sort(discharge,'ascend');
    sigma2Sorted = sigma2(order);
    sigmaSorted = sigma(order);
    sourceIndexSorted = sourceIndex(order);
    window = 2*options.half_window + 1;

    switch options.method
        case 'moving'
            if numel(dischargeSorted) < window
                error('analyze_discharge_differencing:WindowTooWide', ...
                    ['The record requires at least %d valid pairs for ' ...
                    'a moving half-window of %d.'], ...
                    window,options.half_window);
            end
            first = options.half_window + 1;
            last = numel(dischargeSorted) - options.half_window;
            dischargeFit = dischargeSorted(first:last);
            sigma2Fit = movmean(sigma2Sorted,window,'Endpoints','discard');
            sigmaFit = sqrt(max(0,sigma2Fit));
        case 'all'
            dischargeFit = dischargeSorted;
            sigma2Fit = sigma2Sorted;
            sigmaFit = sigmaSorted;
        otherwise
            error('analyze_discharge_differencing:InvalidMethod', ...
                'Unknown method: %s.',options.method);
    end

    meanQ = mean(Qwork,'omitnan');
    if ~isfinite(meanQ) ...
            || meanQ <= 0
        error('analyze_discharge_differencing:InvalidMeanDischarge', ...
            'Mean valid discharge must be positive.');
    end

    Xsigma = [dischargeFit,meanQ*ones(size(dischargeFit))];
    coefficient = local_nnls2(Xsigma,sigmaFit);
    sigmaPredicted = Xsigma*coefficient;
    sigmaStats = local_fit_stats(sigmaFit,sigmaPredicted,2);

    Xvariance = [dischargeFit,ones(size(dischargeFit))];
    varianceCoefficient = local_nnls2(Xvariance,sigma2Fit);
    variancePredicted = Xvariance*varianceCoefficient;
    varianceStats = local_fit_stats(sigma2Fit,variancePredicted,2);

    result = struct();
    result.options = options;
    result.n_observations = n;
    result.n_valid_observations = sum(~invalid);
    result.n_pairs = numel(discharge);
    result.normalization = normalization;
    result.mean_discharge = meanQ;
    result.raw = struct( ...
        'discharge',discharge, ...
        'sigma',sigma, ...
        'sigma2',sigma2, ...
        'source_index',sourceIndex);
    result.sorted = struct( ...
        'discharge',dischargeSorted, ...
        'sigma',sigmaSorted, ...
        'sigma2',sigma2Sorted, ...
        'source_index',sourceIndexSorted);
    result.smooth = struct( ...
        'discharge',dischargeFit, ...
        'sigma',sigmaFit, ...
        'sigma2',sigma2Fit, ...
        'window',window);
    result.fit_sigma = struct( ...
        'slope',coefficient(1), ...
        'intercept',coefficient(2)*meanQ, ...
        'alpha',coefficient(1), ...
        'beta',coefficient(2), ...
        'predicted',sigmaPredicted, ...
        'r2',sigmaStats.r2, ...
        'rmse',sigmaStats.rmse);
    result.fit_variance = struct( ...
        'slope',varianceCoefficient(1), ...
        'intercept',varianceCoefficient(2), ...
        'predicted',variancePredicted, ...
        'paper_curve',(sigmaPredicted).^2, ...
        'r2',varianceStats.r2, ...
        'rmse',varianceStats.rmse);

    result.daily = struct();
    if ~isempty(options.samples_per_day)
        result.daily = local_daily_uncertainty( ...
            Qwork,coefficient,meanQ,options.samples_per_day);
    end

end

function options = local_options(options)

    defaults = struct( ...
        'order',3, ...
        'half_window',100, ...
        'method','moving', ...
        'min_discharge',0, ...
        'min_sigma',0, ...
        'samples_per_day',[]);
    if isfield(options,'tolerance') ...
            && ~isempty(options.tolerance) ...
            && (~isfield(options,'min_sigma') ...
            || isempty(options.min_sigma))
        options.min_sigma = options.tolerance;
    end
    names = fieldnames(defaults);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(options,name) ...
                || isempty(options.(name))
            options.(name) = defaults.(name);
        end
    end
    options.order = round(double(options.order));
    options.half_window = round(double(options.half_window));
    options.method = lower(strtrim(char(string(options.method))));
    options.min_discharge = double(options.min_discharge);
    options.min_sigma = double(options.min_sigma);
    if ~isempty(options.samples_per_day)
        options.samples_per_day = round(double(options.samples_per_day));
    end
    options.tolerance = options.min_sigma;
    if ~isscalar(options.order) ...
            || options.order < 1 ...
            || options.order > 10
        error('analyze_discharge_differencing:InvalidOrder', ...
            'options.order must be an integer from 1 through 10.');
    end
    if ~isscalar(options.half_window) ...
            || options.half_window < 0
        error('analyze_discharge_differencing:InvalidWindow', ...
            'options.half_window must be a nonnegative integer.');
    end
    if ~ismember(options.method,{'moving','all'})
        error('analyze_discharge_differencing:InvalidMethod', ...
            'options.method must be ''moving'' or ''all''.');
    end
    if ~isempty(options.samples_per_day) ...
            && (~isscalar(options.samples_per_day) ...
            || options.samples_per_day < 2)
        error('analyze_discharge_differencing:InvalidSamplesPerDay', ...
            'options.samples_per_day must be empty or at least two.');
    end

end

function daily = local_daily_uncertainty(Q,coefficient,meanQ,u)
%LOCAL_DAILY_UNCERTAINTY Apply de Oliveira and Vrugt (2022), Equation 5.

    nDay = floor(numel(Q)/u);
    Q = Q(1:nDay*u);
    hourly = reshape(Q,u,nDay);
    complete = all(isfinite(hourly),1);
    dailyDischarge = mean(hourly,1).';
    withinVariance = var(hourly,0,1).';
    hourlySigma = coefficient(1)*hourly + coefficient(2)*meanQ;
    randomVariance = sum(hourlySigma.^2,1).'/(u^2);
    frequencyVariance = withinVariance/u;
    dailyVariance = frequencyVariance + randomVariance;

    dailyDischarge(~complete) = NaN;
    frequencyVariance(~complete) = NaN;
    randomVariance(~complete) = NaN;
    dailyVariance(~complete) = NaN;
    daily = struct( ...
        'samples_per_day',u, ...
        'complete',complete(:), ...
        'discharge',dailyDischarge, ...
        'frequency_variance',frequencyVariance, ...
        'random_error_variance',randomVariance, ...
        'variance',dailyVariance, ...
        'sigma',sqrt(dailyVariance));
end

function coefficient = local_nnls2(X,y)

    coefficientLS = X\y;
    candidates = zeros(2,4);
    candidates(:,1) = max(0,coefficientLS);
    candidates(1,2) = max(0,X(:,1)'*y/(X(:,1)'*X(:,1)));
    candidates(2,3) = max(0,X(:,2)'*y/(X(:,2)'*X(:,2)));
    errorSum = inf(1,4);
    for j = 1:4
        if all(candidates(:,j) >= 0)
            residual = y - X*candidates(:,j);
            errorSum(j) = residual'*residual;
        end
    end
    [~,best] = min(errorSum);
    coefficient = candidates(:,best);

end

function stats = local_fit_stats(observed,predicted,nParameter)

    residual = observed - predicted;
    sse = sum(residual.^2);
    sst = sum((observed - mean(observed)).^2);
    if sst > 0
        r2 = 1 - sse/sst;
    else
        r2 = NaN;
    end
    denominator = max(1,numel(observed) - nParameter);
    stats = struct('r2',r2,'rmse',sqrt(sse/denominator));

end
