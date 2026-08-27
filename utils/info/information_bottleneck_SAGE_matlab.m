function info = information_bottleneck_SAGE_matlab(A,H,Y,opts)
%INFORMATION_BOTTLENECK_SAGE_MATLAB Numerical information-bottleneck diagnostics for SAGE.
%
% info = information_bottleneck_SAGE_matlab(A,H,Y,opts)
%
% Numerical computation only. This function creates no figures, axes,
% UI controls, dashboards, colormaps, or graphics. Visualization is
% handled exclusively by the SAGE GUI.
%
% A   : K x r basin-attribute matrix
% H   : 1 x L cell array; H{ell} is K x n_ell hidden activations
% Y   : K x p relevance/target matrix
%
% -------------------------------------------------------------------------
% CORE OPTIONS
% -------------------------------------------------------------------------
%   .nbins = 10
%   .target_nbins = 8
%   .quantization = 'quantile'      % 'equalwidth' or 'quantile'
%   .range = []                     % optional [lo hi] fixed activation range
%   .base = 2                       % information units; 2 = bits
%   .pairwise = true                % adjacent-neuron binned MI
%
% -------------------------------------------------------------------------
% KSG OPTIONS
% -------------------------------------------------------------------------
%   .ksg = true
%   .ksg_method = 'matlab'          % MATLAB reference implementation only
%   .ksg_k = 5                      % nearest-neighbor order
%   .ksg_max_dim = 4                % skip high-dimensional whole-layer KSG
%   .ksg_target_indices = 1         % targets for KSG neuron diagnostics
%   .ksg_standardize = true
%   .ksg_jitter = 1e-10             % deterministic tie-breaking jitter
%   .ksg_pairwise = false           % KSG for all adjacent neuron pairs is costly
%
%   This MATLAB implementation is fully self-contained. All KSG-1 mutual
%   information calculations, including neuron-target estimates, whole-layer
%   estimates, adjacent-neuron estimates, bootstrap calculations, and
%   permutation/Excess-MI calculations, use the local MATLAB KSG routine in
%   this file. No MEX function is called from this reference implementation.
%
% -------------------------------------------------------------------------
% UNCERTAINTY / FINITE-SAMPLE OPTIONS
% -------------------------------------------------------------------------
%   .bootstrap_n = 200              % 0 disables binned bootstrap
%   .permutation_n = 200            % 0 disables binned permutation baseline
%   .ci = [2.5 97.5]                % percentile interval
%   .random_seed = 1729
%
% KSG uncertainty is more expensive and is therefore opt-in:
%   .ksg_bootstrap_n = 0
%   .ksg_permutation_n = 0
%   .ksg_uncertainty_neurons = false
%   .ksg_uncertainty_target_indices = 1 % selected targets only; avoids huge cost
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% Binned / quantized diagnostics:
%   info.I_X_T
%   info.I_T_Y
%   info.I_T_Theta
%   info.I_T_theta
%   info.H_theta
%   info.NMI_T_theta
%   info.H_T
%   info.H_Y
%   info.H_neuron
%   info.I_neuron_theta
%   info.NMI_neuron_theta
%   info.I_adjacent
%   info.H_neuron_norm
%   info.effective_width
%   info.active_neuron_mask
%
% KSG continuous-variable diagnostics:
%   info.ksg.I_T_theta
%   info.ksg.I_neuron_theta
%   info.ksg.I_adjacent
%   info.ksg.valid_layer
%
% Finite-sample diagnostics:
%   info.uncertainty.binned.bootstrap_mean
%   info.uncertainty.binned.bootstrap_lo
%   info.uncertainty.binned.bootstrap_hi
%   info.uncertainty.binned.permutation_mean
%   info.uncertainty.binned.permutation_p
%   info.uncertainty.binned.excess
%
% Optional KSG uncertainty:
%   info.uncertainty.ksg.layer.*
%   info.uncertainty.ksg.neuron{ell}.*
%
% Sample-size metadata:
%   info.sample.K
%   info.sample.max_empirical_entropy
%   info.sample.target_bin_occupancy
%   info.sample.notes
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% 1. The binned information-plane estimator is retained for comparability
%    with early Information Bottleneck work.
%
% 2. Whole-layer joint-state MI can saturate when K is small relative to
%    dimensionality of the hidden representation. For that reason, this
%    function additionally provides low-dimensional KSG estimates.
%
% 3. KSG is most useful here for low-dimensional relations such as
%       neuron <-> NSE
%       neuron <-> KGE
%       neuron <-> one hydrologic parameter
%       adjacent neuron <-> adjacent neuron
%    Whole-layer KSG is skipped when dim(T_l)+dim(Y_j) exceeds
%    opts.ksg_max_dim.
%
% 4. Bootstrap intervals quantify resampling variability. Permutation tests
%    estimate a finite-sample/null MI baseline. The excess-information
%    diagnostic is
%
%       I_excess = I_observed - mean(I_permuted).
%
% 5. KSG bootstrap estimates can be sensitive to duplicated observations.
%    A tiny deterministic jitter is therefore applied after standardization.
%
% Example for NSE:
%   [nTheta,H] = ffn_theta('eval_info',phi,A);
%   Kt = bas.K_t;
%   Ht = H;
%   for ell = 1:numel(Ht)
%       Ht{ell} = Ht{ell}(1:Kt,:);
%   end
%   opts = struct();
%   opts.target_name = '\mathrm{NSE}';
%   opts.target_names = {'NSE'};
%   opts.quantization = 'quantile';
%   opts.ksg = true;
%   info = information_bottleneck_SAGE( ...
%       A(:,1:Kt).',Ht,NSE.tt(:),opts);

    if nargin < 4 ...
            || isempty(opts)
        opts = struct();
    end
    opts = local_defaults(opts);

    A = double(A);
    Y = double(Y);

    if ~iscell(H) ...
            || isempty(H)
        error('information_bottleneck_SAGE:H', ...
            'H must be a nonempty cell array.');
    end

    K0 = size(H{1},1);

    if size(A,1) ~= K0
        error('information_bottleneck_SAGE:A', ...
            'A must have K rows.');
    end
    if size(Y,1) ~= K0
        error('information_bottleneck_SAGE:Y', ...
            'Y must have K rows.');
    end
    for ell = 1:numel(H)
        if size(H{ell},1) ~= K0
            error('information_bottleneck_SAGE:HSize', ...
                'All H{ell} must have K rows.');
        end
    end

    % ----------------------
    % Finite-value screening
    % ----------------------
    keep = all(isfinite(A),2) & all(isfinite(Y),2);

    for ell = 1:numel(H)
        keep = keep & all(isfinite(H{ell}),2);
    end
    if nnz(keep) < 5
        error('information_bottleneck_SAGE:TooFewSamples', ...
            'Fewer than five finite samples remain.');
    end

    A = A(keep,:);
    Y = Y(keep,:);
    H = cellfun(@(x)double(x(keep,:)), ...
        H,'UniformOutput',false);

    K = size(A,1);
    p = size(Y,2);
    L = numel(H);

    % -------------------------
    % Quantize target variables
    % -------------------------
    Yq = local_quantize_matrix( ...
        Y,opts.target_nbins, ...
        opts.quantization,[]);

    yState = local_row_state(Yq);
    H_Y = local_entropy_discrete(yState,opts.base);

    H_theta = nan(1,p);
    for j = 1:p
        H_theta(j) = local_entropy_discrete( ...
            double(Yq(:,j)),opts.base);
    end

    % -----------------
    % Initialize output
    % -----------------
    info = struct();

    info.K = K;
    info.nLayers = L;
    info.nTargets = p;
    info.nAttributes = size(A,2);
    info.options = opts;
    info.keep = keep;

    info.I_X_T = nan(1,L);
    info.I_T_Y = nan(1,L);
    info.I_T_Theta = nan(1,L);
    info.I_T_theta = nan(L,p);
    info.NMI_T_theta = nan(L,p);

    info.H_T = nan(1,L);
    info.H_Y = H_Y;
    info.H_theta = H_theta;

    info.H_neuron = cell(1,L);
    info.H_neuron_norm = cell(1,L);
    info.active_neuron_mask = cell(1,L);
    info.effective_width = nan(1,L);
    info.I_neuron_theta = cell(1,L);
    info.NMI_neuron_theta = cell(1,L);
    info.I_adjacent = cell(1,max(L-1,0));

    info.quantized = struct( ...
        'H',{cell(1,L)}, ...
        'Y',Yq);

    % --------------------
    % Sample-size metadata
    % --------------------
    info.sample = struct();
    info.sample.K = K;
    info.sample.max_empirical_entropy = ...
        log(K)/log(opts.base);
    info.sample.target_bin_occupancy = ...
        K / opts.target_nbins;
    info.sample.hidden_bin_occupancy = ...
        K / opts.nbins;
    info.sample.notes = { ...
        sprintf('K = %d finite samples.',K), ...
        sprintf('Maximum empirical entropy = %.3f bits.', ...
            info.sample.max_empirical_entropy), ...
        sprintf('Nominal target occupancy = %.1f samples/bin.', ...
            info.sample.target_bin_occupancy), ...
        ['Whole-layer joint-state MI may ' ...
        'saturate when hidden dimension ' ...
         'is large relative to K; prefer ' ...
         'neuron-level KSG for inference.']};

    % ----------------------------------------
    % Binned layer-wise information quantities
    % ----------------------------------------
    for ell = 1:L

        Hq = local_quantize_matrix( ...
            H{ell},opts.nbins, ...
            opts.quantization,opts.range);

        info.quantized.H{ell} = Hq;

        tState = local_row_state(Hq);

        info.H_T(ell) = ...
            local_entropy_discrete(tState,opts.base);

        % Deterministic finite-sample information-plane coordinate.
        info.I_X_T(ell) = info.H_T(ell);

        info.I_T_Y(ell) = ...
            local_mutual_information_discrete( ...
            tState,yState,opts.base);

        info.I_T_Theta(ell) = info.I_T_Y(ell);

        for j = 1:p

            info.I_T_theta(ell,j) = ...
                local_mutual_information_discrete( ...
                tState,double(Yq(:,j)), ...
                opts.base);

            if isfinite(H_theta(j)) ...
                    && H_theta(j) > 0

                info.NMI_T_theta(ell,j) = ...
                    info.I_T_theta(ell,j) ...
                    / H_theta(j);
            end
        end

        nNeuron = size(Hq,2);

        hNeuron = nan(1,nNeuron);
        Mtheta = nan(nNeuron,p);

        for neuron = 1:nNeuron

            hNeuron(neuron) = ...
                local_entropy_discrete( ...
                double(Hq(:,neuron)), ...
                opts.base);

            for j = 1:p

                Mtheta(neuron,j) = ...
                    local_mutual_information_discrete( ...
                    double(Hq(:,neuron)), ...
                    double(Yq(:,j)), ...
                    opts.base);
            end
        end

        info.H_neuron{ell} = hNeuron;

        % Normalize scalar-neuron entropy by the largest entropy possible
        % under the selected hidden-activation quantizer. This yields a
        % compact [0,1]-style activity diagnostic. A neuron below
        % opts.active_entropy_threshold is treated as effectively inactive.
        HmaxNeuron = log(min(opts.nbins,K)) / log(opts.base);
        if isfinite(HmaxNeuron) ...
                && HmaxNeuron > 0
            hNeuronNorm = hNeuron / HmaxNeuron;
        else
            hNeuronNorm = nan(size(hNeuron));
        end
        hNeuronNorm = min(max(hNeuronNorm,0),1);
        activeMask = isfinite(hNeuronNorm) ...
            & hNeuronNorm >= opts.active_entropy_threshold;

        info.H_neuron_norm{ell} = hNeuronNorm;
        info.active_neuron_mask{ell} = activeMask;
        info.effective_width(ell) = sum(activeMask);

        info.I_neuron_theta{ell} = Mtheta;

        Mnmi = nan(size(Mtheta));

        for j = 1:p
            if isfinite(H_theta(j)) ...
                    && H_theta(j) > 0
                Mnmi(:,j) = ...
                    Mtheta(:,j) / H_theta(j);
            end
        end

        info.NMI_neuron_theta{ell} = Mnmi;
    end

    % -------------------------
    % Binned adjacent-neuron MI
    % -------------------------
    if opts.pairwise && L >= 2

        for ell = 1:L-1

            Qa = info.quantized.H{ell};
            Qb = info.quantized.H{ell+1};

            M = nan(size(Qa,2),size(Qb,2));

            for a = 1:size(Qa,2)
                for b = 1:size(Qb,2)

                    M(a,b) = ...
                        local_mutual_information_discrete( ...
                        double(Qa(:,a)), ...
                        double(Qb(:,b)), ...
                        opts.base);
                end
            end

            info.I_adjacent{ell} = M;
        end
    end

    % ==========================================
    % KSG continuous-variable mutual information
    % ==========================================
    info.ksg = struct();
    info.ksg.k = opts.ksg_k;
    info.ksg.method = 'matlab';
    info.ksg.max_dim = opts.ksg_max_dim;
    info.ksg.target_indices = ...
        local_valid_target_indices( ...
        opts.ksg_target_indices,p);

    info.ksg.uncertainty_target_indices = ...
        local_valid_target_indices( ...
        opts.ksg_uncertainty_target_indices,p);

    info.ksg.I_T_theta = nan(L,p);
    info.ksg.I_neuron_theta = cell(1,L);
    info.ksg.I_adjacent = cell(1,max(L-1,0));
    info.ksg.valid_layer = false(L,p);

    if opts.ksg

        ksgTargets = info.ksg.target_indices;

        for ell = 1:L

            nNeuron = size(H{ell},2);
            Mk = nan(nNeuron,p);

            % Whole-layer KSG only when dimensionality remains defensible.
            for jj = 1:numel(ksgTargets)

                j = ksgTargets(jj);

                totalDim = size(H{ell},2) + 1;

                if totalDim <= opts.ksg_max_dim

                    info.ksg.I_T_theta(ell,j) = ...
                        local_ksg_mi( ...
                        H{ell},Y(:,j),opts);

                    info.ksg.valid_layer(ell,j) = true;
                end
            end

            % Neuron-level KSG: scalar neuron vs scalar target.
            % This reference implementation is intentionally MATLAB-only.
            for neuron = 1:nNeuron
                for jj = 1:numel(ksgTargets)

                    j = ksgTargets(jj);

                    Mk(neuron,j) = ...
                        local_ksg_mi( ...
                        H{ell}(:,neuron), ...
                        Y(:,j),opts);
                end
            end

            info.ksg.I_neuron_theta{ell} = Mk;
        end

        % Optional adjacent-neuron KSG.
        if opts.ksg_pairwise ...
                && L >= 2

            for ell = 1:L-1
                Ha = H{ell};
                Hb = H{ell+1};
                M = nan(size(Ha,2),size(Hb,2));
                for a = 1:size(Ha,2)
                    for b = 1:size(Hb,2)
                        M(a,b) = ...
                            local_ksg_mi( ...
                            Ha(:,a),Hb(:,b),opts);
                    end
                end
                info.ksg.I_adjacent{ell} = M;
            end
        end
    end

    % =======================================
    % Binned bootstrap + permutation baseline
    % =======================================
    info.uncertainty = struct();
    info.uncertainty.binned = ...
        local_empty_uncertainty(L,p);

    if opts.bootstrap_n > 0 ...
            || opts.permutation_n > 0

        stream = RandStream( ...
            'mt19937ar','Seed',opts.random_seed);
        for ell = 1:L

            tState = local_row_state( ...
                info.quantized.H{ell});
            for j = 1:p
                yj = double(Yq(:,j));
                obs = info.I_T_theta(ell,j);
                U = local_mi_uncertainty_discrete( ...
                    tState,yj,obs,opts,stream);
                info.uncertainty.binned.bootstrap_mean(ell,j) = ...
                    U.bootstrap_mean;
                info.uncertainty.binned.bootstrap_lo(ell,j) = ...
                    U.bootstrap_lo;
                info.uncertainty.binned.bootstrap_hi(ell,j) = ...
                    U.bootstrap_hi;
                info.uncertainty.binned.permutation_mean(ell,j) = ...
                    U.permutation_mean;
                info.uncertainty.binned.permutation_p(ell,j) = ...
                    U.permutation_p;
                info.uncertainty.binned.excess(ell,j) = ...
                    U.excess;
            end
        end
    end

    % ========================
    % Optional KSG uncertainty
    % ========================
    info.uncertainty.ksg = struct();
    info.uncertainty.ksg.layer = ...
        local_empty_uncertainty(L,p);

    info.uncertainty.ksg.neuron = ...
        cell(1,L);

    if opts.ksg ...
            && (opts.ksg_bootstrap_n > 0 ...
            || opts.ksg_permutation_n > 0)

        streamK = RandStream( ...
            'mt19937ar', ...
            'Seed',opts.random_seed + 104729);
        ksgTargets = info.ksg.target_indices;
        % Whole-layer KSG uncertainty only for low-dimensional valid layers.
        for ell = 1:L

            for jj = 1:numel(ksgTargets)

                j = ksgTargets(jj);
                if ~info.ksg.valid_layer(ell,j)
                    continue
                end
                obs = info.ksg.I_T_theta(ell,j);
                U = local_mi_uncertainty_ksg( ...
                    H{ell},Y(:,j),obs,opts,streamK);
                info.uncertainty.ksg.layer = ...
                    local_store_uncertainty( ...
                    info.uncertainty.ksg.layer, ...
                    U,ell,j);
            end
        end

        % Neuron-level KSG uncertainty is opt-in because it can be costly.
        if opts.ksg_uncertainty_neurons

            uncertaintyTargets = ...
                info.ksg.uncertainty_target_indices;
            for ell = 1:L
                nNeuron = size(H{ell},2);
                UN = struct();
                UN.bootstrap_mean = nan(nNeuron,p);
                UN.bootstrap_lo = nan(nNeuron,p);
                UN.bootstrap_hi = nan(nNeuron,p);
                UN.permutation_mean = nan(nNeuron,p);
                UN.permutation_p = nan(nNeuron,p);
                UN.excess = nan(nNeuron,p);
                for neuron = 1:nNeuron
                    for jj = 1:numel(uncertaintyTargets)
                        j = uncertaintyTargets(jj);
                        obs = ...
                            info.ksg.I_neuron_theta{ell}(neuron,j);
                        U = local_mi_uncertainty_ksg( ...
                            H{ell}(:,neuron), ...
                            Y(:,j),obs,opts,streamK);
                        UN.bootstrap_mean(neuron,j) = ...
                            U.bootstrap_mean;
                        UN.bootstrap_lo(neuron,j) = ...
                            U.bootstrap_lo;
                        UN.bootstrap_hi(neuron,j) = ...
                            U.bootstrap_hi;
                        UN.permutation_mean(neuron,j) = ...
                            U.permutation_mean;
                        UN.permutation_p(neuron,j) = ...
                            U.permutation_p;
                        UN.excess(neuron,j) = ...
                            U.excess;
                    end
                end
                info.uncertainty.ksg.neuron{ell} = UN;
            end
        end
    end

end


% ==================================
function opts = local_defaults(opts)
% ==================================

    D = struct( ...
        'nbins',10, ...
        'target_nbins',8, ...
        'quantization','quantile', ...
        'range',[], ...
        'base',2, ...
        'pairwise',true, ...
        'ksg',true, ...
        'ksg_method','matlab', ...
        'ksg_k',5, ...
        'ksg_max_dim',4, ...
        'ksg_target_indices',1, ...
        'ksg_standardize',true, ...
        'ksg_jitter',1e-10, ...
        'ksg_pairwise',false, ...
        'bootstrap_n',200, ...
        'permutation_n',200, ...
        'ci',[2.5 97.5], ...
        'random_seed',1729, ...
        'ksg_bootstrap_n',0, ...
        'ksg_permutation_n',0, ...
        'ksg_uncertainty_neurons',false, ...
        'ksg_uncertainty_target_indices',1, ...
        'active_entropy_threshold',0.10);

    fn = fieldnames(D);

    for i = 1:numel(fn)

        if ~isfield(opts,fn{i}) ...
                || isempty(opts.(fn{i}))

            opts.(fn{i}) = D.(fn{i});
        end
    end

    opts.nbins = ...
        max(2,round(double(opts.nbins)));

    opts.target_nbins = ...
        max(2,round(double(opts.target_nbins)));

    opts.ksg_k = ...
        max(1,round(double(opts.ksg_k)));

    opts.ksg_max_dim = ...
        max(2,round(double(opts.ksg_max_dim)));

    opts.active_entropy_threshold = ...
        min(1,max(0,double(opts.active_entropy_threshold)));

    opts.bootstrap_n = ...
        max(0,round(double(opts.bootstrap_n)));

    opts.permutation_n = ...
        max(0,round(double(opts.permutation_n)));

    opts.ksg_bootstrap_n = ...
        max(0,round(double(opts.ksg_bootstrap_n)));

    opts.ksg_permutation_n = ...
        max(0,round(double(opts.ksg_permutation_n)));

    if numel(opts.ci) ~= 2 ...
            || any(~isfinite(opts.ci)) ...
            || opts.ci(1) < 0 ...
            || opts.ci(2) > 100 ...
            || opts.ci(1) >= opts.ci(2)

        error( ...
            'information_bottleneck_SAGE:ci', ...
            'opts.ci must be [lower upper] percentiles in [0,100].');
    end

    opts.ksg_method = lower(char(string(opts.ksg_method)));

    if ~strcmp(opts.ksg_method,'matlab')

        error( ...
            'information_bottleneck_SAGE_matlab:ksgMethod', ...
            ['information_bottleneck_SAGE_matlab is intentionally ' ...
             'self-contained and MATLAB-only. Use opts.ksg_method = ' ...
             '''matlab''. For the native implementation, call ' ...
             'information_bottleneck_SAGE_mex instead.']);
    end

    if ~ismember( ...
            lower(opts.quantization), ...
            {'equalwidth','quantile'})

        error( ...
            'information_bottleneck_SAGE:quantization', ...
            ['quantization must be ' ...
             'equalwidth or quantile.']);
    end
end


% =====================================
function Q = local_quantize_matrix( ...
    X,nbins,method,fixedRange)
% =====================================

    X = double(X);
    [n,p] = size(X);

    Q = zeros(n,p,'uint16');

    for j = 1:p

        x = X(:,j);
        ux = unique(x);

        if numel(ux) <= nbins ...
                && all(abs(ux-round(ux)) < 1e-12)

            [~,~,q] = unique(x,'sorted');
            Q(:,j) = uint16(q);
            continue
        end

        if strcmpi(method,'quantile')

            edges = quantile( ...
                x,linspace(0,1,nbins+1));

            edges(1) = -inf;
            edges(end) = inf;

            edges = unique(edges,'stable');

            if numel(edges) <= 2
                q = ones(n,1);
            else
                q = discretize(x,edges);
                q(~isfinite(q)) = 1;
            end

        else

            if isempty(fixedRange)
                lo = min(x);
                hi = max(x);
            else
                lo = fixedRange(1);
                hi = fixedRange(2);
            end

            if ~isfinite(lo) ...
                    || ~isfinite(hi) ...
                    || hi <= lo

                q = ones(n,1);

            else

                edges = linspace( ...
                    lo,hi,nbins+1);

                edges(1) = -inf;
                edges(end) = inf;

                q = discretize(x,edges);
                q(~isfinite(q)) = 1;
            end
        end

        Q(:,j) = uint16(q);
    end
end


% =================================
function state = local_row_state(Q)
% =================================

    [~,~,state] = ...
        unique(Q,'rows','stable');

    state = double(state);
end


% ==========================================
function Hx = local_entropy_discrete(x,base)
% ==========================================

    x = x(:);
    x = x(isfinite(x));

    if isempty(x)
        Hx = NaN;
        return
    end

    [~,~,g] = unique(x);

    counts = accumarray(g,1);

    p = counts / sum(counts);
    p = p(p > 0);

    Hx = -sum( ...
        p .* (log(p) / log(base)));
end


% =================================================
function I = local_mutual_information_discrete( ...
    x,y,base)
% =================================================

    x = x(:);
    y = y(:);

    keep = isfinite(x) ...
        & isfinite(y);

    x = x(keep);
    y = y(keep);

    if isempty(x)
        I = NaN;
        return
    end

    [~,~,ix] = unique(x);
    [~,~,iy] = unique(y);

    C = accumarray( ...
        [ix iy],1, ...
        [max(ix) max(iy)]);

    Pxy = C / sum(C(:));

    Px = sum(Pxy,2);
    Py = sum(Pxy,1);

    I = 0;

    for i = 1:size(Pxy,1)
        for j = 1:size(Pxy,2)

            pxy = Pxy(i,j);

            if pxy > 0

                I = I + ...
                    pxy * ...
                    (log( ...
                    pxy/(Px(i)*Py(j))) ...
                    / log(base));
            end
        end
    end

    if I < 0 ...
            && I > -1e-12

        I = 0;
    end
end


% =================================
function I = local_ksg_mi(X,Y,opts)
%LOCAL_KSG_MI Kraskov-Stoegbauer-Grassberger MI estimator (KSG-1).
%
% Uses the max/Chebyshev norm in joint space:
%
% I(X;Y) = psi(k) + psi(N)
%          - mean[psi(n_x+1) + psi(n_y+1)]
%
% and converts nats to opts.base units.
% =================================

    X = double(X);
    Y = double(Y);

    if isvector(X)
        X = X(:);
    end

    if isvector(Y)
        Y = Y(:);
    end

    good = all(isfinite(X),2) ...
        & all(isfinite(Y),2);

    X = X(good,:);
    Y = Y(good,:);

    N = size(X,1);
    k = opts.ksg_k;

    if N <= k + 1
        I = NaN;
        return
    end

    if size(X,2) + size(Y,2) ...
            > opts.ksg_max_dim

        I = NaN;
        return
    end

    if opts.ksg_standardize
        X = local_standardize_columns(X);
        Y = local_standardize_columns(Y);
    end

    if opts.ksg_jitter > 0
        X = local_deterministic_jitter( ...
            X,opts.ksg_jitter,1);
        Y = local_deterministic_jitter( ...
            Y,opts.ksg_jitter,101);
    end

    Z = [X Y];

    nx = zeros(N,1);
    ny = zeros(N,1);

    for i = 1:N

        dz = max(abs(Z - Z(i,:)),[],2);
        dz(i) = inf;

        ds = sort(dz,'ascend');
        epsilon = ds(k);

        if ~isfinite(epsilon)
            I = NaN;
            return
        end

        % Strict-radius count, matching KSG-1. Pull the radius inward by a
        % scale-aware amount so boundary ties are excluded.
        r = epsilon ...
            - max(1e-12, ...
            32*eps(max(1,epsilon)));

        r = max(r,0);

        dx = max(abs(X - X(i,:)),[],2);
        dy = max(abs(Y - Y(i,:)),[],2);

        dx(i) = inf;
        dy(i) = inf;

        nx(i) = sum(dx < r);
        ny(i) = sum(dy < r);
    end

    I_nats = psi(k) ...
        + psi(N) ...
        - mean( ...
        psi(nx + 1) ...
        + psi(ny + 1));

    I = I_nats / log(opts.base);

    % Finite-sample KSG estimates can be slightly negative.
    if I < 0 && I > -1e-6
        I = 0;
    end
end


% =======================================
function X = local_standardize_columns(X)
% =======================================

    mu = mean(X,1);
    sd = std(X,0,1);

    sd(~isfinite(sd) | sd <= 0) = 1;

    X = (X - mu) ./ sd;
end


% ===================================================
function X = local_deterministic_jitter(X,amp,offset)
% ===================================================

    if amp <= 0
        return
    end

    n = size(X,1);
    p = size(X,2);

    ii = (1:n).';

    for j = 1:p

        s = std(X(:,j),0,1);

        if ~isfinite(s) ...
                || s <= 0
            s = 1;
        end

        X(:,j) = X(:,j) ...
            + amp*s*sin( ...
            ii*(j + offset + sqrt(2)));
    end
end


% =============================================
function U = local_mi_uncertainty_discrete( ...
    x,y,obs,opts,stream)
% =============================================

    N = numel(x);

    U = local_empty_uncertainty_scalar();

    if opts.bootstrap_n > 0

        boot = nan(opts.bootstrap_n,1);

        for b = 1:opts.bootstrap_n

            idx = randi(stream,N,N,1);

            boot(b) = ...
                local_mutual_information_discrete( ...
                x(idx),y(idx),opts.base);
        end

        U.bootstrap_mean = ...
            mean(boot,'omitnan');

        U.bootstrap_lo = ...
            local_percentile( ...
            boot,opts.ci(1));

        U.bootstrap_hi = ...
            local_percentile( ...
            boot,opts.ci(2));
    end

    if opts.permutation_n > 0

        perm = nan(opts.permutation_n,1);

        for b = 1:opts.permutation_n

            idx = randperm(stream,N);

            perm(b) = ...
                local_mutual_information_discrete( ...
                x,y(idx),opts.base);
        end

        U.permutation_mean = ...
            mean(perm,'omitnan');

        U.permutation_p = ...
            (1 + sum(perm >= obs)) ...
            / (opts.permutation_n + 1);

        U.excess = ...
            obs - U.permutation_mean;
    end
end


% ========================================
function U = local_mi_uncertainty_ksg( ...
    X,Y,obs,opts,stream)
% ========================================

    N = size(X,1);

    U = local_empty_uncertainty_scalar();

    if opts.ksg_bootstrap_n > 0

        boot = nan( ...
            opts.ksg_bootstrap_n,1);

        for b = 1:opts.ksg_bootstrap_n

            idx = randi(stream,N,N,1);

            boot(b) = ...
                local_ksg_mi( ...
                X(idx,:),Y(idx,:),opts);
        end

        U.bootstrap_mean = ...
            mean(boot,'omitnan');

        U.bootstrap_lo = ...
            local_percentile( ...
            boot,opts.ci(1));

        U.bootstrap_hi = ...
            local_percentile( ...
            boot,opts.ci(2));
    end

    if opts.ksg_permutation_n > 0

        perm = nan( ...
            opts.ksg_permutation_n,1);

        for b = 1:opts.ksg_permutation_n

            idx = randperm(stream,N);

            perm(b) = ...
                local_ksg_mi( ...
                X,Y(idx,:),opts);
        end

        U.permutation_mean = ...
            mean(perm,'omitnan');

        U.permutation_p = ...
            (1 + sum(perm >= obs)) ...
            / (opts.ksg_permutation_n + 1);

        U.excess = ...
            obs - U.permutation_mean;
    end
end


% =======================================
function U = local_empty_uncertainty(L,p)
% =======================================

    U = struct();

    U.bootstrap_mean = nan(L,p);
    U.bootstrap_lo = nan(L,p);
    U.bootstrap_hi = nan(L,p);

    U.permutation_mean = nan(L,p);
    U.permutation_p = nan(L,p);

    U.excess = nan(L,p);
end


% ===========================================
function U = local_empty_uncertainty_scalar()
% ===========================================

    U = struct( ...
        'bootstrap_mean',NaN, ...
        'bootstrap_lo',NaN, ...
        'bootstrap_hi',NaN, ...
        'permutation_mean',NaN, ...
        'permutation_p',NaN, ...
        'excess',NaN);
end


% ===========================================
function S = local_store_uncertainty(S,U,i,j)
% ===========================================

    S.bootstrap_mean(i,j) = ...
        U.bootstrap_mean;

    S.bootstrap_lo(i,j) = ...
        U.bootstrap_lo;

    S.bootstrap_hi(i,j) = ...
        U.bootstrap_hi;

    S.permutation_mean(i,j) = ...
        U.permutation_mean;

    S.permutation_p(i,j) = ...
        U.permutation_p;

    S.excess(i,j) = ...
        U.excess;
end


% ===============================================
function q = local_percentile(x,pct)
%LOCAL_PERCENTILE Toolbox-free linear percentile.
% ===============================================

    x = sort(x(isfinite(x)));

    if isempty(x)
        q = NaN;
        return
    end

    if isscalar(x)
        q = x;
        return
    end

    pos = 1 ...
        + (numel(x)-1) ...
        * pct/100;

    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        q = x(lo);
    else
        w = pos - lo;
        q = (1-w)*x(lo) ...
            + w*x(hi);
    end
end


% ==============================================
function idx = local_valid_target_indices(idx,p)
% ==============================================

    if isempty(idx)
        idx = 1:p;
    end

    idx = unique( ...
        round(double(idx(:).')), ...
        'stable');

    idx = idx( ...
        isfinite(idx) ...
        & idx >= 1 ...
        & idx <= p);

    if isempty(idx)
        idx = 1:min(p,1);
    end
end




