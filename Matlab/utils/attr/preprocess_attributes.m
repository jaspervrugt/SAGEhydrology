function [A,state] = preprocess_attributes(raw,training,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREPROCESS_ATTRIBUTES Fit one of nine base-MATLAB attribute transforms.
%  RAW is basins-by-attributes; TRAINING contains training row indices.
%  OPTIONS: mode (matlab/training/custom), method (1..9), cohort
%  (training/selected), missing (reject/median), methods (0 global, 1..9).
%  STATE contains all fitted statistics needed for independent inference.
%  MATLAB preset: selected cohort, sample SD, reject nonfinite values.
%  Training preset: training cohort, population SD, training median fill.
%  Copyright (c) 2026 Jasper A. Vrugt, University of California, Irvine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 3, options = struct(); end
    mode = local_option(options,'mode','matlab');
    method = local_option(options,'method',1);
    cohort = local_option(options,'cohort','training');
    missing = local_option(options,'missing','reject');
    population = false;
    if strcmp(mode,'matlab')
        method = 1; cohort = 'selected'; missing = 'reject';
    elseif strcmp(mode,'training')
        method = 1; cohort = 'training'; missing = 'median';
        population = true;
    elseif ~strcmp(mode,'custom')
        error('attributes:Mode','Unknown preprocessing preset.');
    end
    x = double(raw); training = double(training(:));
    if ~ismatrix(x) || ...
        isempty(training) || ...
        any(training < 1 | ...
        training > size(x,1) | ...
        training ~= fix(training))
        error('attributes:Rows','Invalid attribute matrix/training rows.');
    end
    if ~ismember(string(cohort),["training","selected"]) || ...
        ~ismember(string(missing),["reject","median"])
        error('attributes:Options','Invalid cohort or missing policy.');
    end
    if strcmp(missing,'reject') && ...
        any(~isfinite(x),'all')
        error('attributes:Missing','Nonfinite attributes are rejected.');
    end
    med = median(x(training,:),1,'omitmissing');
    for j = 1:size(x,2)
        v = x(training,j); v = v(isfinite(v));
        if isempty(v)
            error('attributes:Empty','No finite training attribute.');
        end
        med(j) = median(v);
        x(~isfinite(x(:,j)),j) = med(j);
    end
    fitted = x;
    if strcmp(cohort,'training'), fitted = x(training,:); end
    methods = repmat(double(method),1,size(x,2));
    if strcmp(mode,'custom') && ...
        isfield(options,'methods') && ...
        ~isempty(options.methods)
        overrides = double(options.methods(:).');
        if numel(overrides) ~= size(x,2)
            error('attributes:Methods','Attribute method count differs.');
        end
        if any(overrides < 0 | ...
            overrides > 9 | ...
            overrides ~= fix(overrides))
            error('attributes:Methods','Methods must be integers 0..9.');
        end
        methods(overrides ~= 0) = overrides(overrides ~= 0);
    end
    if any(~isfinite(methods) | ...
        methods < 1 | ...
        methods > 9 | ...
        methods ~= fix(methods))
        error('attributes:Method','Methods must be integers 1..9.');
    end
    columns = cell(1,size(x,2));
    for j = 1:size(x,2)
        v = fitted(:,j); m = methods(j);
        c = struct('method',m,'center',0,'scale',1);
        if m == 6
            c.lambda = 1;
            if max(v) > min(v)
                fun = @(p)local_power_loss(v,p);
                [c.lambda,~] = fminsearch(fun,1, ...
                    optimset('Display','off','TolX',1e-10, ...
                    'TolFun',1e-10,'MaxFunEvals',2000));
            end
        end
        if ismember(m,[7 8])
            [knots,~,group] = unique(v);
            counts = accumarray(group,1);
            first = cumsum([0;counts(1:end-1)]);
            probs = (first+(counts-1)/2)/max(1,numel(v)-1);
            if numel(knots) == 1, probs = 0.5; end
            c.knots = knots; c.probabilities = probs;
        end
        z = attribute_nonlinear(v,c);
        if ismember(m,[1 5 6])
            c.center = mean(z);
            if numel(z) > 1, c.scale = std(z,double(population)); end
        elseif ismember(m,[2 3])
            c.center = min(v); c.scale = max(v)-min(v);
        elseif m == 4
            c.center = median(v);
            q = sage_percentile_linear(v,[25 75]);
            c.scale = q(2)-q(1);
        end
        if ~isfinite(c.scale) || ...
            c.scale == 0
            c.scale = 1;
        end
        columns{j} = c;
    end
    state = struct('version',2,'mode',mode,'missing',missing, ...
        'cohort',cohort,'median',med,'columns',{columns});
    A = apply_attribute_preprocessing(raw,state);
end

function value = local_option(options,name,fallback)
%LOCAL_OPTION Resolve an optional preprocessing control.
    value = fallback;
    if isfield(options,name), value = options.(name); end
end

function loss = local_power_loss(v,p)
%LOCAL_POWER_LOSS Negative Yeo-Johnson Gaussian profile likelihood.
    z = attribute_nonlinear(v,struct('method',6,'lambda',p));
    variance = mean((z-mean(z)).^2);
    if ~isfinite(variance) || ...
        variance <= 0
        loss = realmax; return
    end
    loss = numel(v)/2*log(variance) ...
        -(p-1)*sum(sign(v).*log1p(abs(v)));
end
