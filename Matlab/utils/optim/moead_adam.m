function out = moead_adam(objective,lower,upper,opts)
%MOEAD_ADAM Gradient-enhanced decomposition with neighborhood cooperation.
% objective(x) returns losses (2-by-1) and gradients (d-by-2).
% Experimental gradient variant of Zhang & Li (2007), not classical MOEA/D.
% G=tau*logsumexp(lambda.*(L-z)/tau);
% grad G = J_L * (softmax(lambda.*(L-z)/tau).*lambda).
% Bounds, shared ideal point, independent Adam states, and full archive.
% Reference: https://doi.org/10.1109/TEVC.2007.892759

    arguments
        objective (1,1) function_handle
        lower (:,1) double
        upper (:,1) double
        opts.population (1,1) double {mustBeInteger,mustBeGreaterThan( ...
            opts.population,1)} = 21
        opts.iterations (1,1) double {mustBeInteger,mustBePositive} = 150
        opts.neighbors (1,1) double {mustBeInteger,mustBePositive} = 5
        opts.learningRate (1,1) double {mustBePositive} = 0.01
        opts.tau (1,1) double {mustBePositive} = 0.02
        opts.seed (1,1) double = 1
        opts.parallel (1,1) logical = false
        opts.verbose (1,1) logical = true
    end

    assert(numel(lower)==numel(upper) && ...
        all(upper>lower),'Invalid bounds');
    oldRng = rng;
    cleanup = onCleanup(@()rng(oldRng));
    rng(opts.seed);
    d = numel(lower);
    K = opts.population;
    lambda = [linspace(0,1,K); linspace(1,0,K)];
    distance = sum((reshape(lambda,2,K,1) ...
        - reshape(lambda,2,1,K)).^2,1);
    [~,neighbor] = sort(reshape(distance,K,K),1);
    neighbor = neighbor(1:min(K,opts.neighbors),:);
    X=lower + (upper-lower).*(0.1+0.8*rand(d,K));
    [L,J] = evaluate(objective,X,opts.parallel);
    ideal = min(L,[],2);
    archiveX = X;
    archiveL = L;
    [archiveX,archiveL] = archive(archiveX,archiveL);
    m = zeros(d,K);
    v = m;
    age = zeros(1,K);
    history = zeros(opts.iterations,3);
    clock = tic;
    for iteration=1:opts.iterations
        % All candidates share this iteration's fixed ideal point.
        scaled = lambda .* (L-ideal)/opts.tau;
        scaled = scaled - max(scaled,[],1);
        probability = exp(scaled);
        probability = probability ./ sum(probability,1);
        grad = zeros(d,K);
        for k=1:K
            grad(:,k) = J(:,:,k) * (probability(:,k) .* lambda(:,k));
        end
        mNew = .9*m + .1*grad;
        vNew = .999*v + .001*grad.^2;
        ageNew = age + 1;
        proposal = X - opts.learningRate * (mNew./(1-.9.^ageNew))./ ...
            (sqrt(vNew./(1-.999.^ageNew)) + 1e-8);
        proposal = min(upper,max(lower,proposal));
        [Lp,Jp] = evaluate(objective,proposal,opts.parallel);
        ideal = min(ideal,min(Lp,[],2));
        [archiveX,archiveL] = archive([archiveX proposal], ...
            [archiveL Lp]);
        % Select cooperatively from the old incumbent and nearby proposals.
        % Compare every contender with the same updated ideal point.
        for k = 1:K
            candidates = neighbor(:,k)';
            values = [scalar(L(:,k),lambda(:,k),ideal,opts.tau), ...
                scalar(Lp(:,candidates),lambda(:,k),ideal,opts.tau)];
            [~,best] = min(values);
            if best==1
                m(:,k) = mNew(:,k);
                v(:,k) = vNew(:,k);
                age(k) = ageNew(k);
                continue
            end
            source = candidates(best-1);
            X(:,k) = proposal(:,source);
            L(:,k) = Lp(:,source);
            J(:,:,k) = Jp(:,:,source);
            if source==k
                m(:,k) = mNew(:,k);
                v(:,k) = vNew(:,k);
                age(k) = ageNew(k);
            else
                m(:,k) = 0;
                v(:,k) = 0;
                age(k) = 0;
            end
        end
        history(iteration,:) = [ideal' size(archiveL,2)];
        if opts.verbose && ...
            (iteration==1 || ...
            mod(iteration,10)==0 || ...
            iteration==opts.iterations)
            fprintf(['MOEA/D-Adam %d/%d ' ...
                '| ideal %.5g %.5g ' ...
                '| archive %d ' ...
                '| elapsed %.1f s\n'], ...
                iteration,opts.iterations, ...
                ideal,size(archiveL,2),toc(clock));
        end
    end
    out=struct('x',X, ...
        'losses',L, ...
        'weights',lambda, ...
        'ideal',ideal, ...
        'archiveX',archiveX, ...
        'archiveLosses',archiveL, ...
        'history',history, ...
        'elapsedSeconds',toc(clock), ...
        'evaluations',K*(opts.iterations+1), ...
        'options',opts);
end

function [L,J] = evaluate(fun,X,useParallel)
    K = size(X,2);
    d = size(X,1);
    L = zeros(2,K);
    J = zeros(d,2,K);
    if useParallel
        parfor k=1:K
            [lk,jk] = fun(X(:,k));
            L(:,k) = lk(:);
            J(:,:,k) = jk;
        end
    else
        for k=1:K
            [lk,jk] = fun(X(:,k));
            L(:,k) = lk(:);
            J(:,:,k) = jk;
        end
    end
    assert(all(isfinite(L),'all') && ...
        all(isfinite(J),'all'), ['Objective returned nonfinite ' ...
            'losses or gradients. Check data and model.']);
end

function value = scalar(L,lambda,ideal,tau)
    a = lambda.*(L-ideal)/tau;
    shift = max(a,[],1);
    value = tau*(shift + log(sum(exp(a-shift),1)));
end

function [X,L] = archive(X,L)
    [~,uniqueIndex] = unique(L','rows','stable');
    X = X(:,uniqueIndex);
    L = L(:,uniqueIndex);
    keep = true(1,size(L,2));
    for k = 1:size(L,2)
        keep(k) = ~any(all(L<=L(:,k),1) & ...
            any(L<L(:,k),1));
    end
    X = X(:,keep);
    L = L(:,keep);
end
