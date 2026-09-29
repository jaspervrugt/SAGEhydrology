function [split,mdl] = build_split(mdl,prd,bas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD_SPLIT Build temporal training and evaluation indices.
%
%  Constructs the scored window, split masks, and index ranges shared by
%  training and postprocessing.
%
% SYNOPSIS:
%   [split,mdl] = build_split(mdl,prd,bas)
%
% INPUT ARGUMENTS:
%   mdl             model and split settings
%   prd             period and split settings
%    .dt             samples per day: 1, 24, or 96
%    .spinup         optional spin-up in days
%    .dts/.dte       start/end of training period
%    .des/.dee       start/end of evaluation period
%    .ds/.de         start/end of complete scored period
%    .tr_frac        optional training fraction
%    .ev_frac        optional evaluation fraction
%    .block          optional block length
%    .seed           optional random seed
%    .n_folds        optional number of folds
%    .fold           optional evaluation fold
%    .method         optional split method
%   bas             basin counts and optional mode override
%
% OUTPUT ARGUMENTS:
%   split           temporal split and assessment design
%    .method         selected split method
%    .mode           assessment design code (1 to 4)
%    .idx            complete scored index range
%    .id_train       training indices
%    .id_eval        evaluation indices
%    .mask_train     logical training mask
%    .mask_eval      logical evaluation mask
%   mdl             model settings with downstream split indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3
        bas = struct();
    end
    
    if ~isfield(prd,'dt') ...
            || isempty(prd.dt)
        error(['      Error:build_split: ' ...
            'prd.dt is missing.']);
    end
    
    if ~isfield(prd,'spinup') ...
            || isempty(prd.spinup)
        spinup = 0;
    else
        spinup = double(prd.spinup);
    end
    
    if isfield(mdl,'mode') ...
            && ~isempty(mdl.mode)
        mode = local_parse_mode(mdl.mode);
    elseif isfield(bas,'mode') ...
            && ~isempty(bas.mode)
        mode = local_parse_mode(bas.mode);
    else
        mode = local_parse_mode( ...
            get_assess_mode(prd,bas));
    end
    
    if isfield(mdl,'sp_method') ...
            && ~isempty(mdl.sp_method)
        method = lower(string(mdl.sp_method));
    elseif isfield(prd,'method') ...
            && ~isempty(prd.method)
        method = lower(string(prd.method));
    else
        method = "manual";
    end
    method = char(method);
    
    dt = double(prd.dt);
    
    switch lower(method)
    
        case 'manual'
            if ~isfield(prd,'dts') ...
                    || isempty(prd.dts) ...
                    || ~isfield(prd,'dte') ...
                    || isempty(prd.dte)
                error(['      Error:build_split: ' ...
                    'Training period is incomplete ' ...
                    'for method = manual.']);
            end
    
            dts = local_to_ymd(prd.dts);
            dte = local_to_ymd(prd.dte);
    
            switch mode
                case {1,3}
                    des = dts;
                    dee = dte;
    
                case {2,4}
                    if ~isfield(prd,'des') ...
                            || isempty(prd.des) ...
                            || ~isfield(prd,'dee') ...
                            || isempty(prd.dee)
                        error(['      Error:build_split: ' ...
                            'Evaluation period is ' ...
                            'incomplete for mode = %d.'], ...
                            mode);
                    end
                    des = local_to_ymd(prd.des);
                    dee = local_to_ymd(prd.dee);
    
                otherwise
                    error(['      Error:build_split: ' ...
                        'Unknown mode = %d.'], mode);
            end
    
            split = struct();
            split.dt = dt;
            split.mode = mode;
            split.spinup = spinup;
            split.method = method;
            split.dts = dts;
            split.dte = dte;
            split.des = des;
            split.dee = dee;
    
            dtTrain0 = datetime(dts(1), ...
                dts(2),dts(3),0,0,0);
            dtTrain1 = datetime(dte(1), ...
                dte(2),dte(3),0,0,0);
            dtEval0  = datetime(des(1), ...
                des(2),des(3),0,0,0);
            dtEval1  = datetime(dee(1), ...
                dee(2),dee(3),0,0,0);
    
            if dtTrain1 < dtTrain0
                error(['      Error:build_split: ' ...
                    'Training end date must be ' ...
                    'after training start date.']);
            end
    
            if dtEval1 < dtEval0
                error(['      Error:build_split: ' ...
                    'Evaluation end date must be ' ...
                    'after evaluation start date.']);
            end
    
            dtFirst = min(dtTrain0,dtEval0);
            dtLast = max(dtTrain1,dtEval1);
    
            split.dt0 = dtFirst - days(spinup);
            split.dt_train0 = dtTrain0;
            split.dt_train1 = dtTrain1;
            split.dt_eval0 = dtEval0;
            split.dt_eval1 = dtEval1;
            split.dt_end = dtLast;
    
            [tout,idx,id_train,id_eval] = ...
                local_build_global_indices(split);
    
            split.tout = tout;
            split.idx = idx;
            split.id_train = id_train;
            split.id_eval = id_eval;
    
        case {'deterministic_block', ...
                'random_block', ...
                'random', ...
                'deterministic_kfold', ...
                'random_kfold', ...
                'rainfall_block'}
    
            if ~isfield(prd,'ds') ...
                    || isempty(prd.ds) ...
                    || ~isfield(prd,'de') ...
                    || isempty(prd.de)
                error(['      Error:build_split: ' ...
                    'ds/de are required ' ...
                    'for method = %s.'], ...
                    method);
            end
    
            d0 = local_to_ymd(prd.ds);
            d1 = local_to_ymd(prd.de);
    
            dtStart = datetime(d0(1), ...
                d0(2),d0(3),0,0,0);
            dtEnd   = datetime(d1(1), ...
                d1(2),d1(3),0,0,0);
    
            if dtEnd < dtStart
                error(['      Error:build_split: ' ...
                    'End day must be on/after ' ...
                    'start day.']);
            end
    
            split = struct();
            split.dt = dt;
            split.mode = mode;
            split.spinup = spinup;
            split.method = method;
            split.ds = d0;
            split.de = d1;
    
            split.dt0 = dtStart - days(spinup);
            split.dt_train0 = dtStart;
            split.dt_train1 = dtEnd;
            split.dt_eval0 = dtStart;
            split.dt_eval1 = dtEnd;
            split.dt_end = dtEnd;
    
            [id_all,tout,idx] = ...
                local_build_candidate_indices(split);
    
            switch lower(method)
                case 'deterministic_block'
                    [id_train,id_eval] = ...
                        local_build_deterministic_block_indices( ...
                        id_all,prd);
    
                case 'random_block'
                    [id_train,id_eval] = ...
                        local_build_random_block_indices( ...
                        id_all,prd);
    
                case 'random'
                    [id_train,id_eval] = ...
                        local_build_random_indices_from_pool( ...
                        id_all,prd);
    
                case 'deterministic_kfold'
                    [id_train,id_eval] = ...
                        local_build_deterministic_kfold_indices( ...
                        id_all,prd);
    
                case 'random_kfold'
                    [id_train,id_eval] = ...
                        local_build_random_kfold_indices( ...
                        id_all,prd);
                
                case 'rainfall_block'
                    id_train = [];
                    id_eval = [];
                    split.local = true;
                    split.split_source = 'meteo';
                    split.rainfall_block = true;
            end
    
            split.tout = tout;
            split.idx = idx;
            split.id_train = id_train;
            split.id_eval = id_eval;
    
        otherwise
            error(['      Error:build_split: ' ...
                'Unknown split method = %s.'], ...
                method);
    end
    
    n_scored = numel(split.idx);
    split.mask_train = false(1,n_scored);
    split.mask_eval = false(1,n_scored);
    
    if isfield(split,'local') ...
            && split.local
        split.id_train = [];
        split.id_eval = [];
        split.n_train = NaN;
        split.n_eval = NaN;
    else
        if ~isempty(split.id_train)
            split.mask_train( ...
                split.id_train) = true;
        end
        if ~isempty(split.id_eval)
            split.mask_eval( ...
                split.id_eval) = true;
        end
    
        split.n_train = ...
            numel(split.id_train);
        split.n_eval = ...
            numel(split.id_eval);
    end
    
    % ------------------------------------------
    % Reconcile mode with actual basin/time split
    % ------------------------------------------
    hasEvalBasins = isfield(bas,'K_e') ...
        && ~isempty(bas.K_e) && ...
        (double(bas.K_e) > 0);
    
    %hasEvalTime = ~isempty(id_eval);
    hasEvalTime = (~isempty(id_eval)) ...
        || (isfield(split,'local') ...
        && split.local);
    
    % mode codes:
    % 1 = training basins | training period
    % 2 = training basins | evaluation period/mask
    % 3 = training + evaluation basins | training period
    % 4 = training + evaluation basins | evaluation period/mask
    split.mode = 1 + 2*hasEvalBasins + hasEvalTime;
    
    % ------------------
    % Assign back to mdl
    % ------------------
    mdl.sp_method = lower(method);
    mdl.mode = split.mode;
    mdl.tout = tout(end);
    mdl.idx = [idx(1) idx(end)];
    
    if isfield(split,'local') ...
            && split.local
        mdl.id_train = [];
        mdl.id_eval = [];
        mdl.local = true;
    else
        mdl.local = false;
        if local_is_contiguous(id_train)
            mdl.id_train = [id_train(1) id_train(end)];
        else
            mdl.id_train = id_train(:).';
        end
        if local_is_contiguous(id_eval)
            mdl.id_eval = [id_eval(1) id_eval(end)];
        else
            mdl.id_eval = id_eval(:).';
        end
    end

end

function mode = local_parse_mode(x)
%LOCAL_PARSE_MODE Convert mode to canonical code 1..4

    if isnumeric(x)
        mode = double(x);
    elseif isstring(x) ...
            || ischar(x)
        s = lower(strtrim(char(x)));
    
        switch s
            case {'1','train_train','tt'}
                mode = 1;
            case {'2','train_eval','te'}
                mode = 2;
            case {'3','all_train','et'}
                mode = 3;
            case {'4','all_eval','ee'}
                mode = 4;
            otherwise
                error(['      Error:build_split: ' ...
                    'Unknown mode string: %s.'],s);
        end
    else
        error(['      Error:build_split: ' ...
            'mode must be numeric or text.']);
    end
    
    if ~isscalar(mode) ...
            || ~ismember(mode,[1 2 3 4])
        error(['      Error:build_split: ' ...
            'mode must be one of 1,2,3,4.']);
    end

end

function [id_all,tout,idx] = ...
    local_build_candidate_indices(split)
%LOCAL_BUILD_CANDIDATE_INDICES

    dt0 = split.dt0;
    dtStart = datetime(split.ds(1), ...
        split.ds(2),split.ds(3),0,0,0);
    dtEnd = datetime(split.de(1), ...
        split.de(2),split.de(3),0,0,0);
    dt = double(split.dt);
    
    switch dt
        case 1
            n0 = int64(days(dateshift(dtStart, ...
                'start','day') - dateshift(dt0, ...
                'start','day')));
            n1 = int64(days(dateshift(dtEnd, ...
                'start','day') - dateshift(dt0, ...
                'start','day')));
            nEnd = n1;
    
        case 24
            t0 = datetime(year(dt0), ...
                month(dt0),day(dt0),0,0,0);
            tStart = datetime(year(dtStart), ...
                month(dtStart),day(dtStart),0,0,0);
            tEnd = datetime(year(dtEnd), ...
                month(dtEnd),day(dtEnd),23,0,0);
    
            n0 = int64(round(hours(tStart - t0)));
            n1 = int64(round(hours(tEnd - t0)));
            nEnd = n1;
    
        case 96
            t0 = datetime(year(dt0), ...
                month(dt0),day(dt0),0,0,0);
            tStart = datetime(year(dtStart), ...
                month(dtStart),day(dtStart),0,0,0);
            tEnd = datetime(year(dtEnd), ...
                month(dtEnd),day(dtEnd),23,45,0);
        
            n0 = int64(round(minutes(tStart - t0) / 15));
            n1 = int64(round(minutes(tEnd - t0) / 15));
            nEnd = n1;
    
        otherwise
            error(['      Error:build_split: ' ...
                'Unsupported dt = %d.'],dt);
    end
    
    tout = 0:double(nEnd + 1);
    idx = (double(n0)+1):(double(n1)+2);
    id_all = 1:(double(n1-n0)+1);

end

function frac_train = ...
    local_get_train_frac(prd)
%LOCAL_GET_TRAIN_FRAC

    if isfield(prd,'tr_frac') ...
            && ~isempty(prd.tr_frac)
        frac_train = double(prd.tr_frac);
    elseif isfield(prd,'ev_frac') ...
            && ~isempty(prd.ev_frac)
        frac_train = 1 - double(prd.ev_frac);
    else
        frac_train = 0.8;
    end
    
    if ~isfinite(frac_train) ...
            || frac_train <= 0 ...
            || frac_train >= 1
        error(['      Error:build_split: ' ...
            'tr_frac must be strictly ' ...
            'between 0 and 1.']);
    end

end

function [id_train,id_eval] = ...
    local_build_deterministic_block_indices(id_all,prd)
%LOCAL_BUILD_DETERMINISTIC_BLOCK_INDICES

    frac_train = local_get_train_frac(prd);
    
    if isfield(prd,'block') ...
            && ~isempty(prd.block)
        block_size = round(double(prd.block));
    else
        block_size = 30;
    end
    
    if block_size < 1
        error(['      Error:build_split: ' ...
            'block must be >= 1.']);
    end
    
    nAll = numel(id_all);
    nBlocks = ceil(nAll / block_size);
    blocks = cell(nBlocks,1);
    
    for b = 1:nBlocks
        i0 = (b-1)*block_size + 1;
        i1 = min(b*block_size,nAll);
        blocks{b} = id_all(i0:i1);
    end
    
    nTrainBlocks = max(1,min(nBlocks-1, ...
        round(frac_train*nBlocks)));
    
    id_train = [];
    for b = 1:nTrainBlocks
        id_train = [id_train blocks{b}]; %#ok
    end
    
    id_eval = setdiff(id_all,id_train,'stable');
    
    if isempty(id_train) ...
            || isempty(id_eval)
        error(['      Error:build_split: ' ...
            'deterministic_block produced ' ...
            'empty training or evaluation indices.']);
    end
    
    id_train = sort(id_train);
    id_eval  = sort(id_eval);
end

function [id_train,id_eval] = ...
    local_build_random_indices_from_pool(id_all,prd)
%LOCAL_BUILD_RANDOM_INDICES_FROM_POOL

    nAll = numel(id_all);
    if nAll < 2
        error(['      Error:build_split: ' ...
            'Random split requires at ' ...
            'least two scored indices.']);
    end
    
    frac_train = local_get_train_frac(prd);
    
    if isfield(prd,'seed') ...
            && ~isempty(prd.seed)
        rng(double(prd.seed),'twister');
    end
    
    nTrain = max(1,min(nAll-1, ...
        round(frac_train*nAll)));
    rp = randperm(nAll);
    
    id_train = sort(id_all(rp(1:nTrain)));
    id_eval  = sort(id_all(rp(nTrain+1:end)));
end

function [id_train,id_eval] = ...
    local_build_random_block_indices(id_all,prd)
%LOCAL_BUILD_RANDOM_BLOCK_INDICES

    nAll = numel(id_all);
    if nAll < 2
        error(['      Error:build_split: ' ...
            'random_block requires at ' ...
            'least two scored indices.']);
    end
    
    frac_train = local_get_train_frac(prd);
    
    if isfield(prd,'block') ...
            && ~isempty(prd.block)
        block_size = round(double(prd.block));
    else
        block_size = 30;
    end
    
    if block_size < 1
        error(['      Error:build_split: ' ...
            'block must be >= 1.']);
    end
    
    if isfield(prd,'seed') ...
            && ~isempty(prd.seed)
        rng(double(prd.seed),'twister');
    end
    
    nBlocks = ceil(nAll / block_size);
    blocks = cell(nBlocks,1);
    
    for b = 1:nBlocks
        i0 = (b-1)*block_size + 1;
        i1 = min(b*block_size,nAll);
        blocks{b} = id_all(i0:i1);
    end
    
    rp = randperm(nBlocks);
    target_train = round(frac_train*nAll);
    
    id_train = [];
    id_eval  = [];
    nTrain = 0;
    
    for j = 1:nBlocks
        blk = blocks{rp(j)};
        if nTrain < target_train
            id_train = [id_train blk]; %#ok
            nTrain = nTrain + numel(blk);
        else
            id_eval = [id_eval blk]; %#ok
        end
    end
    
    if isempty(id_eval)
        blk = blocks{rp(end)};
        id_eval = blk;
        id_train = setdiff(id_train, ...
            id_eval,'stable');
    elseif isempty(id_train)
        blk = blocks{rp(1)};
        id_train = blk;
        id_eval = setdiff(id_eval, ...
            id_train,'stable');
    end
    
    id_train = sort(id_train);
    id_eval = sort(id_eval);
end

function [id_train,id_eval] = ...
    local_build_deterministic_kfold_indices( ...
    id_all,prd)
%LOCAL_BUILD_DETERMINISTIC_KFOLD_INDICES

    if ~isfield(prd,'n_folds') ...
            || isempty(prd.n_folds)
        n_folds = 5;
    else
        n_folds = round(double(prd.n_folds));
    end
    
    if ~isfield(prd,'fold') ...
            || isempty(prd.fold)
        fold = 1;
    else
        fold = round(double(prd.fold));
    end
    
    if n_folds < 2
        error(['      Error:build_split: ' ...
            'n_folds must be >= 2.']);
    end
    
    if fold < 1 ...
            || fold > n_folds
        error(['      Error:build_split: ' ...
            'fold must lie in [1,n_folds].']);
    end
    
    nAll = numel(id_all);
    edges = round(linspace(0, ...
        nAll,n_folds+1));
    
    i0 = edges(fold) + 1;
    i1 = edges(fold+1);
    
    id_eval = id_all(i0:i1);
    id_train = setdiff(id_all, ...
        id_eval,'stable');
    
    if isempty(id_train) ...
            || isempty(id_eval)
        error(['      Error:build_split: ' ...
            'deterministic_kfold produced ' ...
            'empty training or evaluation indices.']);
    end
end

function [id_train,id_eval] = ...
    local_build_random_kfold_indices( ...
    id_all,prd)
%LOCAL_BUILD_RANDOM_KFOLD_INDICES

    if ~isfield(prd,'n_folds') ...
            || isempty(prd.n_folds)
        n_folds = 5;
    else
        n_folds = round(double(prd.n_folds));
    end
    
    if ~isfield(prd,'fold') ...
            || isempty(prd.fold)
        fold = 1;
    else
        fold = round(double(prd.fold));
    end
    
    if n_folds < 2
        error(['      Error:build_split: ' ...
            'n_folds must be >= 2.']);
    end
    
    if fold < 1 ...
            || fold > n_folds
        error(['      Error:build_split: ' ...
            'fold must lie in [1,n_folds].']);
    end
    
    if isfield(prd,'seed') ...
            && ~isempty(prd.seed)
        rng(double(prd.seed),'twister');
    end
    
    nAll = numel(id_all);
    rp = randperm(nAll);
    id_perm = id_all(rp);
    
    edges = round(linspace(0, ...
        nAll,n_folds+1));
    i0 = edges(fold) + 1;
    i1 = edges(fold+1);
    
    id_eval = sort(id_perm(i0:i1));
    id_train = sort(setdiff(id_perm, ...
        id_eval,'stable'));
    
    if isempty(id_train) ...
            || isempty(id_eval)
        error(['      Error:build_split: ' ...
            'random_kfold produced empty ' ...
            'training or evaluation indices.']);
    end
end

function [tout,idx,id_train,id_eval] = ...
    local_build_global_indices(split)
%LOCAL_BUILD_GLOBAL_INDICES Build global time indices from split only

    dt0 = split.dt0;
    dtTrain0 = split.dt_train0;
    dtTrain1 = split.dt_train1;
    dtEval0 = split.dt_eval0;
    dtEval1 = split.dt_eval1;
    dtEnd = split.dt_end;
    dt = double(split.dt);
    
    switch dt
        case 1
            n0 = int64(days(dateshift( ...
                dtTrain0,'start','day') ...
                - dateshift(dt0,'start','day')));
            n1 = int64(days(dateshift( ...
                dtTrain1,'start','day') ...
                - dateshift(dt0,'start','day')));
            m0 = int64(days(dateshift( ...
                dtEval0,'start','day') ...
                - dateshift(dt0,'start','day')));
            m1 = int64(days(dateshift( ...
                dtEval1,'start','day') ...
                - dateshift(dt0,'start','day')));
            nEnd = int64(days(dateshift( ...
                dtEnd,'start','day') ...
                - dateshift(dt0,'start','day')));
    
        case 24
            t0 = datetime(year(dt0), ...
                month(dt0),day(dt0),0,0,0);
            tTrain0 = datetime(year(dtTrain0), ...
                month(dtTrain0),day(dtTrain0),0,0,0);
            tTrain1 = datetime(year(dtTrain1), ...
                month(dtTrain1),day(dtTrain1),23,0,0);
            tEval0 = datetime(year(dtEval0), ...
                month(dtEval0),day(dtEval0),0,0,0);
            tEval1 = datetime(year(dtEval1), ...
                month(dtEval1),day(dtEval1),23,0,0);
            tEnd = datetime(year(dtEnd), ...
                month(dtEnd),day(dtEnd),23,0,0);
    
            n0 = int64(round(hours(tTrain0 - t0)));
            n1 = int64(round(hours(tTrain1 - t0)));
            m0 = int64(round(hours(tEval0 - t0)));
            m1 = int64(round(hours(tEval1 - t0)));
            nEnd = int64(round(hours(tEnd - t0)));
    
        case 96
            t0 = datetime(year(dt0), ...
                month(dt0),day(dt0),0,0,0);
        
            tTrain0 = datetime(year(dtTrain0), ...
                month(dtTrain0),day(dtTrain0),0,0,0);
            tTrain1 = datetime(year(dtTrain1), ...
                month(dtTrain1),day(dtTrain1),23,45,0);
        
            tEval0 = datetime(year(dtEval0), ...
                month(dtEval0),day(dtEval0),0,0,0);
            tEval1 = datetime(year(dtEval1), ...
                month(dtEval1),day(dtEval1),23,45,0);
        
            tEnd = datetime(year(dtEnd), ...
                month(dtEnd),day(dtEnd),23,45,0);
        
            n0 = int64(round(minutes(tTrain0 - t0) / 15));
            n1 = int64(round(minutes(tTrain1 - t0) / 15));
            m0 = int64(round(minutes(tEval0 - t0) / 15));
            m1 = int64(round(minutes(tEval1 - t0) / 15));
            nEnd = int64(round(minutes(tEnd - t0) / 15));        
    
        otherwise
            error(['      Error:build_split: ' ...
                'Unsupported dt = %d. Use ' ...
                'dt = 1 [daily] or ' ...
                'dt = 24 [hourly].'], dt);
    end
    
    if n0 < 0 ...
            || m0 < 0 ...
            || n1 < n0 ...
            || m1 < m0 ...
            || nEnd < 0
        error(['      Error:build_split: ' ...
            'Invalid global time indexing ' ...
            'derived from split.']);
    end
    
    tout = 0:double(nEnd + 1);
    
    i0 = min(n0,m0);
    i1 = max(n1,m1);
    idx = (double(i0) + 1):(double(i1) + 2);
    
    id_train = (double(n0 - i0) + 1): ...
        (double(n1 - i0) + 1);
    id_eval  = (double(m0 - i0) + 1): ...
        (double(m1 - i0) + 1);
end

function flag = local_is_contiguous(id)
%LOCAL_IS_CONTIGUOUS True if id is a contiguous integer sequence

    id = double(id(:)).';
    
    if isempty(id)
        flag = false;
    elseif isscalar(id)
        flag = true;
    else
        flag = all(diff(id) == 1);
    end
end

function ymd = local_to_ymd(x)
%LOCAL_TO_YMD Convert [d m y] or [y m d] to [y m d]

    x = double(x(:)).';
    
    if numel(x) ~= 3
        error(['      Error:build_split: ' ...
            'Date must have three elements.']);
    end
    
    if x(1) > 1000
        ymd = x;
    elseif x(3) > 1000
        ymd = [x(3) x(2) x(1)];
    else
        error(['      Error:build_split: ' ...
            'Could not infer date format from [%g %g %g].'], ...
            x(1),x(2),x(3));
    end
end
