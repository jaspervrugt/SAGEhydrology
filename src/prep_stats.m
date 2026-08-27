function [dat,loss] = prep_stats(dat,mdl,split,loss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREP_STATS Prepare observed discharge statistics for loss functions
%
% SYNOPSIS: [dat,loss] = prep_stats(dat,mdl,split,loss)
%   dat         cell structure with observed discharge records
%    {k}.y_n     nx1 measured discharge vector for basin k
%   mdl         structure with model and split information
%    .id_train  training-period indices
%    .id_eval   evaluation-period indices
%    .sp_method split design
%   split       structure with temporal information
%    .dt        data resolution: 1 = daily, 24 = hourly, 96 = 15-minute
%    .dt0       first datetime of record, required for JKGE method 4
%   loss        loss-function settings
%    .fnc       scalar loss function
%    .n_win     JKGE window length in days
%    .method    JKGE benchmark method
%                 1 = moving-average mean
%                 2 = section-wise mean
%                 3 = long-term mean
%                 4 = monthly climatology
%   dat         OUTPUT: updated cell structure
%    {k}.stats   observed discharge statistics for train/eval periods
%    {k}.jkge    JKGE benchmark vectors and cached bookkeeping
%    {k}.fdc     cached observed FDC quantities for train/eval periods
%    {k}.hydro   reserved for hydrologic metrics/signatures calculated
%                from observed discharge for the train/eval periods
%   loss        OUTPUT: updated loss-function settings
%    .meta      shared JKGE metadata, including month labels if needed
%    .fdc.D0t   Kx1 FDC references from training observations
%    .fdc.D0e   Kx1 FDC references from evaluation observations
%
% NOTES:
%   1. For JKGE, observed benchmark vectors are stored per basin.
%   2. Benchmark bookkeeping is cached for full, training, and evaluation
%      records to avoid repeated construction during optimization.
%   3. No dense or sparse n-by-n benchmark operators are formed.
%   4. Future observation-derived hydrologic metrics should be computed
%      here once during initialization and stored in dat{k}.hydro, with
%      separate training/evaluation fields. They should not be rebuilt
%      during optimization iterations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 4
        error(['      Error:prep_stats: ' ...
            'loss input argument is required.']);
    end
    
    if ~isstruct(loss)
        error(['      Error:prep_stats: ' ...
            'loss must be a structure.']);
    end
    
    if ~isfield(loss,'fnc') || isempty(loss.fnc)
        error(['      Error:prep_stats: ' ...
            'loss.fnc is required.']);
    end
    
    if ~(isscalar(loss.fnc) ...
            && isnumeric(loss.fnc) ...
            && isfinite(loss.fnc) ...
            && any(double(loss.fnc) == 1:7))
        error(['      Error:prep_stats: ' ...
            'loss.fnc must be one ' ...
            'of 1,2,3,4,5,6,7.']);
    end
        
    split = prepare_split(split);
    
    dt = double(split.dt);
    K = numel(dat);
    
    global_id_train = [];
    global_id_eval = [];
    
    if isfield(mdl,'id_train') ...
            && ~isempty(mdl.id_train)
        global_id_train = ...
            expand_index(mdl.id_train);
    end
    
    if isfield(mdl,'id_eval') ...
            && ~isempty(mdl.id_eval)
        global_id_eval = ...
            expand_index(mdl.id_eval);
    end
    
    doJKGE = isequal(double(loss.fnc),7);
    
    % ------------------------------------------------------
    % JKGE settings: only validate and construct if fnc == 7
    % ------------------------------------------------------
    method = [];
    n_win = [];
    
    if doJKGE
        % ---------------------------------------
        % JKGE benchmark scale in days
        % Default = 31 days if missing or invalid
        % ---------------------------------------
        nJKGE_days = 31;
    
        if ~isfield(loss,'n_win') || isempty(loss.n_win)
    
            fprintf(['      Warning:prep_stats: ' ...
                'loss.n_win not specified; ' ...
                'using default n_win = %d days.\n'], ...
                nJKGE_days);
    
        elseif ~(isscalar(loss.n_win) ...
                && isnumeric(loss.n_win) ...
                && isfinite(loss.n_win) ...
                && loss.n_win >= 1 ...
                && mod(loss.n_win,1) == 0)
    
            fprintf(['      Warning:prep_stats: ' ...
                'loss.n_win is invalid; ' ...
                'using default n_win = %d days.\n'], ...
                nJKGE_days);
    
        else
            nJKGE_days = double(loss.n_win);
        end
    
        % ---------------------------------
        % JKGE benchmark method
        % Default = 1 (moving-average mean)
        % ---------------------------------
        method = 1;
    
        if isfield(loss,'method') && ~isempty(loss.method)
    
            if isscalar(loss.method) ...
                    && isnumeric(loss.method) ...
                    && any(double(loss.method) ...
                    == [1 2 3 4])
    
                method = double(loss.method);
    
            else
                fprintf(['      Warning:prep_stats: ' ...
                    'loss.method invalid; ' ...
                    'using default method = 1 ' ...
                    '[moving-average mean].\n']);
            end
        end
    
        % --------------------------------
        % Enforce JKGE/split compatibility
        % --------------------------------
        if ~isfield(mdl,'sp_method') ...
                || isempty(mdl.sp_method)
            error(['      Error:prep_stats: ' ...
                'mdl.sp_method is required ' ...
                'to validate JKGE benchmark ' ...
                'compatibility.']);
        end
        
        sp_method = lower(string( ...
            mdl.sp_method));
        
        if any(strcmp(sp_method, ...
                ["random","random_kfold"])) ...
                && ~any(method == [3 4])
        
            error(['      Error:prep_stats: ' ...
                'invalid JKGE setup. For split ' ...
                'method "%s", JKGE only supports ' ...
                'method 3 [long-term mean] or ' ...
                'method 4 [monthly climatology]. ' ...
                'Methods 1 [moving-average mean] ' ...
                'and 2 [section-wise mean] ' ...
                'require contiguous ' ...
                'training/evaluation periods ' ...
                'and are not compatible with ' ...
                'random or random k-fold ' ...
                'splits.'],char(sp_method));
        end
        
        n_win = round(dt * nJKGE_days);
    
        % moving-average uses centered odd window
        if method == 1 ...
                && mod(n_win,2) == 0
            n_win = n_win + 1;
        end
    
        % ------------------------------------
        % initialize/update meta only for JKGE
        % only keep month-label information
        % ------------------------------------
        if ~isfield(loss,'meta') ...
                || isempty(loss.meta) ...
                || ~isstruct(loss.meta)
            loss.meta = struct();
        end
    
        if method == 4
            if ~isfield(split,'dt0') ...
                    || isempty(split.dt0)
                error(['      Error:prep_stats: ' ...
                    'split.dt0 is required ' ...
                    'for JKGE monthly ' ...
                    'climatology [method = 4].']);
            end
    
            n_all = numel(dat{1}.y_n);
            t0 = split.dt0;
    
            stepDuration = days(1/dt);
            t_all = (t0:stepDuration: ...
                t0 + stepDuration*(n_all-1))';

            mo_all = uint8(month(t_all));
    
            if numel(mo_all) ~= n_all
                error(['      Error:prep_stats: ' ...
                    'month vector length does ' ...
                    'not match full record.']);
            end
    
            loss.meta.mo_all = ...
                mo_all;
            loss.meta.mo_t = [];
            loss.meta.mo_e = [];
        else
            loss.meta.mo_all = [];
            loss.meta.mo_t = [];
            loss.meta.mo_e = [];
        end
    end
    
    % Basin- and period-specific FDC reference scales
    if ~isfield(loss,'fdc') ...
            || ~isstruct(loss.fdc)
        loss.fdc = struct();
    end
    loss.fdc.D0t = nan(K,1);
    loss.fdc.D0e = nan(K,1);

    fprintf(['... Preparing discharge ' ...
        'statistics %3d%%'],0);
    
    for k = 1:K
        y_n = double(dat{k}.y_n(:));
        
        if isfield(mdl,'local') ...
                && mdl.local == 1 ...
                && isfield(dat{k},'id_train') ...
                && ~isempty(dat{k}.id_train)
            id_train = double(dat{k}.id_train(:));
            if isfield(dat{k},'id_eval') ...
                    && ~isempty(dat{k}.id_eval)
                id_eval = double(dat{k}.id_eval(:));
            else
                id_eval = [];
            end
        else
            id_train = global_id_train(:);
            id_eval = global_id_eval(:);
        end
        hasEval = ~isempty(id_eval);
    
        if isempty(id_train)
            error(['      Error:prep_stats: ' ...
                'training indices are ' ...
                'empty for basin %d.'],k);
        end
    
        if max(id_train) > numel(y_n) ...
                || any(id_train < 1)
            error(['      Error:prep_stats: ' ...
                'training indices exceed ' ...
                'y_n length for basin %d.'],k);
        end
    
        if hasEval ...
                && (max(id_eval) > numel(y_n) ...
                || any(id_eval < 1))
            error(['      Error:prep_stats: ' ...
                'evaluation indices exceed ' ...
                'y_n length for basin %d.'],k);
        end
    
        % -----------------
        % initialize fields
        % -----------------
        dat{k}.stats = struct( ...
            'mut',[], ...
            'stdt',[], ...
            'TSSt',[], ...
            'Syt',[], ...
            'mue',[], ...
            'stde',[], ...
            'TSSe',[], ...
            'Sye',[]);
        dat{k}.jkge = struct( ...
            'm_y',[], ...
            'cache',[], ...
            'cache_t',[], ...
            'cache_e',[]);
    
        dat{k}.fdc = struct( ...
            't',local_empty_fdc_cache(), ...
            'e',local_empty_fdc_cache());

        % ---------------
        % training period
        % ---------------
        if isfield(dat{k},'bad') ...
                && ~isempty(dat{k}.bad)
            bad = dat{k}.bad(:);
            id_train = id_train(~bad(id_train));
        end
        yt = y_n(id_train);
    
        dat{k}.stats.mut = mean(yt, ...
            'omitnan');
        dat{k}.stats.stdt = std(yt, ...
            'omitnan');
        dat{k}.stats.TSSt = sum((yt - ...
            dat{k}.stats.mut).^2, ...
            'omitnan');
        dat{k}.stats.Syt = ...
            local_huber_scale(yt);

        % Cached observed quantities for fast FDC distance
        dat{k}.fdc.t = local_fdc_cache(yt);
        loss.fdc.D0t(k) = dat{k}.fdc.t.D0;

        % -----------------
        % evaluation period
        % -----------------
        if hasEval ...
                && ~isempty(id_eval)
            if isfield(dat{k},'bad') ...
                    && ~isempty(dat{k}.bad)
                bad = dat{k}.bad(:);
                id_eval = id_eval(~bad(id_eval));
            end
            ye = y_n(id_eval);
            
            dat{k}.stats.mue  = mean(ye, ...
                'omitnan');
            dat{k}.stats.stde = std(ye, ...
                'omitnan');
            dat{k}.stats.TSSe = sum((ye - ...
                dat{k}.stats.mue).^2, ...
                'omitnan');
            dat{k}.stats.Sye = ...
                local_huber_scale(ye);
            % Cached observed quantities for fast FDC distance
            dat{k}.fdc.e = local_fdc_cache(ye);
            loss.fdc.D0e(k) = dat{k}.fdc.e.D0;
        end
    
        % ----------------------
        % JKGE benchmark vectors
        % ----------------------
        if doJKGE
            m_y = jkge_benchmark(y_n,method, ...
                n_win,loss.meta.mo_all);
        
            dat{k}.jkge.m_y = single(m_y);
    
            dat{k}.jkge.cache = jkge_cache( ...
                y_n,method,n_win, ...
                loss.meta.mo_all);
        
        end
    
        if mod(k,20)==0 ...
                || k==K
            pct = floor(100*k/K);
            fprintf('\b\b\b\b%3d%%',pct);
        end
    end
    
    fprintf('\b\b\b\b... Done\n');

end

% ================
% helper functions
% ================

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
    if ~isfinite(S_y) ...
            || S_y <= 0

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

    if ~isfinite(S_y) ...
            || S_y <= S_floor
        S_y = S_floor;
    end

end

function fdc = local_empty_fdc_cache()
%LOCAL_EMPTY_FDC_CACHE Initialize observed FDC cache.

    fdc = struct( ...
        'ys',[], ...
        'Py',[], ...
        'S_yy',NaN, ...
        'n',0, ...
        'D0',NaN);

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
    A0 = mean(abs(ys - ymed));
    % Median-flow reference divergence with units of discharge.
    D0 = A0 - 0.5 * S_yy / n^2;
    % Numerical safeguard only
    D0 = max(D0,0);

    fdc = struct( ...
        'ys',ys, ...
        'Py',Py, ...
        'S_yy',S_yy, ...
        'n',n, ...
        'D0',D0);
end
