function [dat,loss] = prep_stats(dat,mdl,split,loss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREP_STATS Prepare observation-specific loss statistics.
%
%  Caches valid indices, reference statistics, JKGE benchmarks, and FDC
%  quantities for each selected observable.
%
% SYNOPSIS:
%   [dat,loss] = prep_stats(dat,mdl,split,loss)
%
% INPUT ARGUMENTS:
%   dat             basin records with named Q, SWE, or SM observations
%   mdl             training/evaluation indices and split method
%   split           resolution and time-origin information
%    .dt             samples per day: 1, 24, or 96
%    .dt0            record start time for monthly JKGE
%   loss            selected loss and observation settings
%    .fnc            loss-function identifier
%    .observed       selected Q, SWE, or SM names
%    .n_win          JKGE window length in days
%    .method         JKGE benchmark method (1 to 4)
%    .fdc            FDC formulation settings
%     .kosugi          optional Kosugi preprocessing switch [false]
%
% OUTPUT ARGUMENTS:
%   dat             basin records with prepared loss statistics
%    {k}.stats       named train/evaluation statistics and indices
%    {k}.jkge        named JKGE benchmark caches
%    {k}.fdc         named train/evaluation FDC caches
%    {k}.hydro.fdc.kosugi  optional training-derived Kosugi parameters
%   loss            loss settings with shared benchmark metadata
%    .meta           shared JKGE time metadata
%    .fdc            observation-specific FDC references
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
    
    if ~isfield(loss,'fnc') ...
            || isempty(loss.fnc)
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
    activeNames = "Q";
    if isfield(loss,'observed') ...
            && ~isempty(loss.observed)
        activeNames = upper(strtrim(string(loss.observed(:))));
        activeNames = unique(activeNames(strlength(activeNames)>0), ...
            'stable');
    end
    unknown = setdiff(activeNames,["Q";"SWE";"SM"]);
    if ~isempty(unknown)
        error('prep_stats:UnknownObservation', ...
            'Unknown observation(s): %s.', ...
            strjoin(cellstr(unknown),', '));
    end
    prepNames = activeNames(:);
    
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
    
        if ~isfield(loss,'n_win') ...
                || isempty(loss.n_win)
    
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
    
        if isfield(loss,'method') ...
                && ~isempty(loss.method)
    
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
    
            n_all = local_record_length(dat,prepNames);
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
    if ~isfield(loss.fdc,'formulation') ...
            || isempty(loss.fdc.formulation)
        loss.fdc.formulation = 1;
    end
    loss.fdc.formulation = local_fdc_formulation( ...
        loss.fdc.formulation);

    % --------------------------------------------------------------
    % Optional Kosugi FDC fit used by FDC-constrained dynamic models.
    % The fit is observation-derived and is performed ONCE using only
    % valid training-period Q. Set loss.fdc.kosugi = false to deactivate.
    % --------------------------------------------------------------
    if ~isfield(loss.fdc,'kosugi') ...
            || isempty(loss.fdc.kosugi)
        loss.fdc.kosugi = false;
    end
    if ~(isscalar(loss.fdc.kosugi) ...
            && (islogical(loss.fdc.kosugi) ...
            || isnumeric(loss.fdc.kosugi)) ...
            && isfinite(double(loss.fdc.kosugi)))
        error('prep_stats:Kosugi', ...
            'loss.fdc.kosugi must be a scalar logical/numeric flag.');
    end
    doKosugi = logical(loss.fdc.kosugi);

    emptyReferences = struct('D0t',nan(K,1),'D0e',nan(K,1), ...
        'D0pt',nan(K,1),'D0pe',nan(K,1), ...
        'D0logpt',nan(K,1),'D0logpe',nan(K,1));
    knownNames = {'Q','SWE','SM'};
    for j = 1:numel(knownNames)
        if isfield(loss.fdc,knownNames{j})
            loss.fdc = rmfield(loss.fdc,knownNames{j});
        end
    end
    for j = 1:numel(prepNames)
        loss.fdc.(char(prepNames(j))) = emptyReferences;
    end

    fprintf(['... Preparing observation ' ...
        'statistics %3d%%'],0);
    
    for k = 1:K
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
    
        % -----------------
        % initialize fields
        % -----------------
        dat{k}.jkge = struct();
    
        % -----------------------------------------------------------
        % Observation-specific preparation for future multi-data loss
        % -----------------------------------------------------------
        dat{k}.stats = struct();
        dat{k}.fdc = struct();
        if ~isfield(dat{k},'hydro') ...
                || ~isstruct(dat{k}.hydro)
            dat{k}.hydro = struct();
        end
        if ~isfield(dat{k}.hydro,'fdc') ...
                || ~isstruct(dat{k}.hydro.fdc)
            dat{k}.hydro.fdc = struct();
        end
        if ~doKosugi && isfield(dat{k}.hydro.fdc,'kosugi')
            dat{k}.hydro.fdc = rmfield(dat{k}.hydro.fdc,'kosugi');
        end

        for j = 1:numel(prepNames)
            name = char(prepNames(j));
            obs = [];
            if isfield(dat{k},'obs') ...
                    && isfield(dat{k}.obs,name)
                obs = dat{k}.obs.(name);
            end
            [stats,fdc] = local_observation_block( ...
                obs,id_train,id_eval);
            dat{k}.stats.(name) = stats;
            dat{k}.fdc.(name) = fdc;
            loss.fdc.(name) = local_store_fdc_reference( ...
                loss.fdc.(name),fdc,k);

            % ------------------------------------------------------
            % Robust three-parameter Kosugi fit for discharge only.
            % IMPORTANT: use TRAINING observations only. p0 is the
            % empirical zero-flow probability; a,b,c describe the
            % conditional positive-flow Kosugi component when p0>0.
            % ------------------------------------------------------
            if doKosugi && strcmp(name,'Q') ...
                    && stats.train.has_data
                q_train = double(obs.value(stats.train.indices));
                dat{k}.hydro.fdc.kosugi = fit_kosugi_fdc(q_train);
            end

            hasObs = ~isempty(obs) ...
                && isfield(obs,'value') ...
                && ~isempty(obs.value);
            if doJKGE && hasObs
                y_n = double(obs.value(:));
                if max(id_train) > numel(y_n) ...
                        || any(id_train < 1)
                    error('prep_stats:TrainingIndex', ...
                        ['Training indices exceed the %s record ' ...
                         'for basin %d.'],name,k);
                end
                if hasEval ...
                        && (max(id_eval) > numel(y_n) ...
                        || any(id_eval < 1))
                    error('prep_stats:EvaluationIndex', ...
                        ['Evaluation indices exceed the %s record ' ...
                         'for basin %d.'],name,k);
                end
                m_y = jkge_benchmark(y_n,method, ...
                    n_win,loss.meta.mo_all);
                dat{k}.jkge.(name) = struct( ...
                    'm_y',single(m_y), ...
                    'cache',jkge_cache(y_n,method,n_win, ...
                    loss.meta.mo_all));
                if strcmp(name,'Q')
                    dat{k}.jkge.m_y = ...
                        dat{k}.jkge.Q.m_y;
                    dat{k}.jkge.cache = ...
                        dat{k}.jkge.Q.cache;
                end
            elseif doJKGE
                dat{k}.jkge.(name) = struct( ...
                    'm_y',[],'cache',[]);
            end
        end
    
        if mod(k,20)==0 ...
                || k==K
            pct = floor(100*k/K);
            fprintf('\b\b\b\b%3d%%',pct);
        end
    end
    
    fprintf('\b\b\b\b... Done\n');

end

function n = local_record_length(dat,names)
%LOCAL_RECORD_LENGTH Find the common model-axis length.

    n = [];
    for k = 1:numel(dat)
        for j = 1:numel(names)
            name = char(names(j));
            if isfield(dat{k},'obs') ...
                    && isfield(dat{k}.obs,name) ...
                    && isfield(dat{k}.obs.(name),'value') ...
                    && ~isempty(dat{k}.obs.(name).value)
                candidate = numel(dat{k}.obs.(name).value);
                if isempty(n)
                    n = candidate;
                elseif candidate ~= n
                    error('prep_stats:ObservationLength', ...
                        ['Selected observation records must share ' ...
                         'one model time axis.']);
                end
            end
        end
    end
    if isempty(n)
        error('prep_stats:NoObservations', ...
            'No selected observation records are available.');
    end
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
        'q0',NaN, ...
        'log_ys',[], ...
        'n',0, ...
        'D0',NaN, ...
        'D0_p',NaN, ...
        'D0_logp',NaN);

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

    % The probability-space formulations use the same constant-median
    % reference. The logarithmic offset equals one percent of the median
    % positive observed discharge, as documented in Appendix C.
    D0_p = mean((ymed - ys).^2);
    yp = ys(ys > 0);
    if isempty(yp)
        q0 = NaN;
        log_ys = nan(size(ys));
        D0_logp = NaN;
    else
        q0 = 0.01*median(yp);
        log_ys = log(ys + q0);
        D0_logp = mean((log(ymed + q0) - log_ys).^2);
    end

    fdc = struct( ...
        'ys',ys, ...
        'Py',Py, ...
        'S_yy',S_yy, ...
        'q0',q0, ...
        'log_ys',log_ys, ...
        'n',n, ...
        'D0',D0, ...
        'D0_p',D0_p, ...
        'D0_logp',D0_logp);
end

function [stats,fdc] = local_observation_block(obs,id_train,id_eval)
%LOCAL_OBSERVATION_BLOCK Prepare one observation type independently.

    emptyPeriod = struct('mean',NaN,'std',NaN,'TSS',NaN, ...
        'huber_scale',NaN,'n',0,'indices',zeros(0,1), ...
        'has_data',false,'variable',false);
    stats = struct('available',false,'active',false, ...
        'units',"",'source',"",'n_total',0, ...
        'train',emptyPeriod,'eval',emptyPeriod);
    fdc = struct('t',local_empty_fdc_cache(), ...
        'e',local_empty_fdc_cache());

    if isempty(obs) ...
            || ~isstruct(obs) ...
            || ~isfield(obs,'value') ...
            || isempty(obs.value)
        return
    end
    y = double(obs.value(:));
    stats.n_total = numel(y);
    if isfield(obs,'units') ...
            && ~isempty(obs.units)
        stats.units = string(obs.units);
    end
    if isfield(obs,'source') ...
            && ~isempty(obs.source)
        stats.source = string(obs.source);
    end
    bad = false(size(y));
    if isfield(obs,'bad') ...
            && ~isempty(obs.bad)
        if numel(obs.bad) ~= numel(y)
            error('prep_stats:ObservationMaskLength', ...
                ['Observation value and ' ...
                'bad-mask lengths differ.']);
        end
        bad = logical(obs.bad(:));
    end
    bad = bad ...
        | ~isfinite(y);
    stats.available = any(~bad);
    stats.train = local_period_stats(y,bad,id_train);
    stats.eval = local_period_stats(y,bad,id_eval);
    stats.active = stats.train.has_data;
    fdc.t = local_fdc_cache(y(stats.train.indices));
    fdc.e = local_fdc_cache(y(stats.eval.indices));
end

function period = local_period_stats(y,bad,indices)
%LOCAL_PERIOD_STATS Statistics and valid indices for one period.

    period = struct('mean',NaN,'std',NaN,'TSS',NaN, ...
        'huber_scale',NaN,'n',0,'indices',zeros(0,1), ...
        'has_data',false,'variable',false);
    if isempty(indices)
        return
    end
    indices = double(indices(:));
    if any(indices < 1) ...
            || any(indices ~= fix(indices)) ...
            || any(indices > numel(y))
        error('prep_stats:ObservationIndex', ...
            ['Observation-period indices ' ...
            'are outside the data record.']);
    end
    indices = indices(~bad(indices));
    values = y(indices);
    period.indices = indices;
    period.n = numel(values);
    period.has_data = period.n > 0;
    if ~period.has_data
        return
    end
    period.mean = mean(values);
    period.std = std(values);
    period.TSS = sum((values-period.mean).^2);
    period.huber_scale = local_huber_scale(values);
    scale = max(1,max(abs(values)));
    period.variable = period.n >= 2 ...
        && isfinite(period.TSS) ...
        && period.TSS > eps(scale^2)*period.n;
end

function refs = local_store_fdc_reference(refs,fdc,k)
%LOCAL_STORE_FDC_REFERENCE Store per-observation FDC reference constants.

    refs.D0t(k) = fdc.t.D0;
    refs.D0e(k) = fdc.e.D0;
    refs.D0pt(k) = fdc.t.D0_p;
    refs.D0pe(k) = fdc.e.D0_p;
    refs.D0logpt(k) = fdc.t.D0_logp;
    refs.D0logpe(k) = fdc.e.D0_logp;
end

function formulation = local_fdc_formulation(value)
%LOCAL_FDC_FORMULATION Normalize the public 6a/6b/6c selection.

    if isnumeric(value) ...
            && isscalar(value) ...
            && isfinite(value)
        formulation = double(value);
    else
        key = lower(regexprep(char(string(value)),'[^a-z0-9]',''));
        switch key
            case {'1','a','fdc','dfdc','physical','physicalcdf'}
                formulation = 1;
            case {'2','b','p','dp','quantile'}
                formulation = 2;
            case {'3','c','logp','dlogp','logquantile'}
                formulation = 3;
            otherwise
                formulation = NaN;
        end
    end
    if ~ismember(formulation,1:3)
        error('prep_stats:BadFDCFormulation', ...
            ['loss.fdc.formulation must ' ...
            'be 1 (d_fdc), 2 (d_p), ' ...
            'or 3 (d_logp).']);
    end
end
