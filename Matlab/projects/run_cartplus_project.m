function result = run_cartplus_project(S,root)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RUN_CARTPLUS_PROJECT Train and evaluate Project 6 in MATLAB.
%   S contains screened forcing, observations, and attributes.
%   MATLAB runs CaRT+, analytic sensitivities, and FFN optimization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2
        root = fileparts(fileparts(mfilename('fullpath')));
    end
    addpath(fullfile(root,'src'),fullfile(root,'models','cartplus'), ...
        fullfile(root,'utils','parameters'));
    S.lower = S.lower(:);
    S.upper = S.upper(:);
    S.fixed = S.fixed(:);
    S.multiplier = S.multiplier(:);
    net = S.net;
    net.r = size(S.A,1);
    net.d = numel(S.lower);
    [phi,~] = ffn_theta('init',net);
    opts = descent('init',S.alg,phi);
    mdl.mcode = S.mcode;
    useOde = isfield(S,'modelConfig') && isfield(S.modelConfig,'cart_formulation') && ...
        strcmp(S.modelConfig.cart_formulation,'ode');
    odeArgs = {};
    if useOde
        settings = rmfield(S.modelConfig,'cart_formulation');
        fields = fieldnames(settings);
        for k = 1:numel(fields)
            odeArgs(end+1:end+2) = {fields{k},settings.(fields{k})};
        end
    end
    [mapping,~] = read_cart_model(struct('model',18,'dt',1));
    mapping.th_min = S.lower;
    mapping.th_max = S.upper;
    if isfield(S,'parameter_scales')
        mapping.parameter_scales = S.parameter_scales;
    elseif isfield(S,'log_scale')
        modes = repmat("linear",numel(S.lower),1);
        modes(logical(S.log_scale)) = "log";
        mapping.parameter_scales = modes;
    end
    K = size(S.A,2);
    train = logical(S.training(:).');
    G = zeros(net.d,K);
    ntrain = nnz(train);
    assert(ntrain > 0,'Project 6 requires training basins.');
    data = cell(1,K);
    for k = 1:K
        data{k} = struct('forcing',S.forcing(:,:,k), ...
            'elevation',S.elevation(:,k),'area',S.area(:,k), ...
            'resolution',1,'floor',S.floor(k));
    end
    objective = nan(S.alg.i_max,1);
    Q = nan(size(S.observed));
    for i = 1:S.alg.i_max
        [theta,cache] = ffn_theta('eval',phi,S.A);
        [par,chain] = map_parameters(theta,S,mapping);
        objective(i) = 0;
        for k = find(train)
            if useOde
                [out,J] = cart_ode_discharge(data{k},par(:,k),18,odeArgs{:});
            else
                [out,J] = cartplus(par(:,k),mdl,data{k});
            end
            Q(:,k) = out(:,1);
            yobs = S.observed(:,k);
            use = S.train_mask(:,k) & isfinite(yobs) & yobs >= 0;
            assert(any(use),'A training basin has no scored discharge.');
            y = yobs(use);
            residual = out(use,1)-y;
            scale = 1;
            if S.loss == 3
                scale = sum((y-mean(y)).^2);
                assert(scale > 0,'NSE requires variable discharge.');
            end
            objective(i) = objective(i)+sum(residual.^2)/scale;
            dL = J(use,S.columns).'*residual;
            G(:,k) = (2/scale)*dL.*chain(:,k);
        end
        objective(i) = objective(i)/ntrain;
        fprintf('Iteration %d/%d   Objective %.8g\n', ...
            i,S.alg.i_max,objective(i));
        if i < S.alg.i_max
            gradient = ffn_theta('back',phi,cache,S.alg,G/ntrain);
            [phi,opts] = descent('dyn',S.alg,phi,gradient,i+1,opts);
        end
    end
    for k = find(~train)
        if useOde
            out = cart_ode_discharge(data{k},par(:,k),18,odeArgs{:});
        else
            out = cartplus(par(:,k),mdl,data{k});
        end
        Q(:,k) = out(:,1);
    end
    nse = nan(K,2);
    kge = nan(K,2);
    for k = 1:K
        [nse(k,1),kge(k,1)] = scores( ...
            S.observed(:,k),Q(:,k),S.train_mask(:,k));
        [nse(k,2),kge(k,2)] = scores( ...
            S.observed(:,k),Q(:,k),S.eval_mask(:,k));
    end
    result = struct('gauges',{S.gauges},'dates',S.dates, ...
        'training',train,'objective',objective,'theta',theta, ...
        'phi',phi,'Q',Q,'NSE',nse,'KGE',kge);
end

function [par,chain] = map_parameters(theta,S,mapping)
%MAP_PARAMETERS Convert normalized FFN outputs to native solver inputs.

    [x,chain] = parameter_scale('physical',mapping,theta);
    par = repmat(S.fixed,1,size(theta,2));
    par(S.columns,:) = x;
    j = logical(S.exp_scale);
    rows = S.columns(j);
    par(rows,:) = S.multiplier(j).*exp(x(j,:));
    chain(j,:) = chain(j,:).*par(rows,:);
end

function [nse,kge] = scores(y,q,mask)
%SCORES Compute NSE and KGE on one scored period.

    use = mask & isfinite(y) & y >= 0 & isfinite(q);
    y = y(use);
    q = q(use);
    nse = NaN;
    kge = NaN;
    if numel(y) < 2 || std(y) == 0
        return
    end
    nse = 1-sum((y-q).^2)/sum((y-mean(y)).^2);
    if mean(y) ~= 0 && std(q) > 0
        r = corrcoef(y,q);
        kge = 1-norm([r(1,2)-1,std(q)/std(y)-1,mean(q)/mean(y)-1]);
    end
end
