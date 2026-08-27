function file_mat = export_SAGE(mdl,dat,bas,prd,dirres,Q,prf,region,nTheta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXPORT_SAGE Export final SAGE training results to a single MAT file
%
% SYNOPSIS: file_mat = export_SAGE(mdl,dat,bas,prd,dirres, ...
%  Q,prf,region,nTheta)
%
%   mdl        structure model information
%    .local     train/eval handling
%                0 -> global train/eval split
%                1 -> basin-specific train/eval split
%    .id_train  training indices of global split
%    .id_eval   evaluation indices of global split
%   dat        1xK cell array basin data structures
%    {k}.y_n    measured discharge series basin k
%    {k}.id_train
%               training indices basin k (if mdl.local = 1)
%    {k}.id_eval
%               evaluation indices basin k (if mdl.local = 1)
%    {k}.use    basin usage flag
%                1 -> training basin
%                2 -> evaluation basin
%    {k}.id_USGS
%               USGS gauge identifier
%    {k}.gname  basin name
%   bas        structure basin information
%    .K_t       number of training basins
%    .K_e       number of evaluation basins
%   prd        structure training/evaluation period information
%    .dt        temporal data resolution
%                [1=daily, 24=hourly, 96=15-minute]
%    .ds        first day complete simulation period
%   dirres     directory for exported MAT file
%   Q          simulated discharge structure
%    .tt        training basins | training split
%    .te        training basins | evaluation split
%    .et        evaluation basins | training split
%    .ee        evaluation basins | evaluation split
%   prf        performance structure; prf.curr contains the final
%              basin-wise metrics and prf.iter the training histories
%    .curr      current basin-wise NSE, KGE, S_fdc and JKGE values;
%               each metric has .tt, .te, .et and .ee scenarios
%    .iter      scalar performance histories across SAGE iterations
%   region     Data region: US/GB/BR
%   nTheta        dxK matrix normalized model parameter values
%               d = number of model parameters
%               K = number of basins
%
%   file_mat   OUTPUT: full path exported MAT file
%
% The exported structure E contains:
%
%   E.basins
%    .row_id    basin row identifier
%    .basin_id  basin identifier
%    .basin_name
%                basin name
%    .basin_group
%                1 -> training basin
%                2 -> evaluation basin
%   E.performance
%    .current    final basin-wise metric structure copied from prf.curr
%    .iteration scalar training histories copied from prf.iter
%    .tt         performance metrics training basin | training split
%    .te         performance metrics training basin | evaluation split
%    .et         performance metrics evaluation basin | training split
%    .ee         performance metrics evaluation basin | evaluation split
%   E.theta_n    normalized parameter matrix
%   E.discharge
%    .split_code
%                KxN matrix train/eval flags
%                1 -> training output
%                2 -> evaluation output
%                NaN -> not used
%    .Q_measured
%                KxN matrix measured discharge
%    .Q_simulated
%                KxN matrix simulated discharge
%    .Q_scenario
%                original scenario-specific discharge structure Q
%    .time       Nx1 vector simulation time
%   E.metadata   metadata and configuration structures
%
% row_id links basins, performance metrics, parameters and discharge matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% © Written by Jasper A. Vrugt, Apr. 2026 / updated Aug. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    curr = prf.curr;
    NSE = curr.NSE;
    KGE = curr.KGE;
    S_fdc = curr.S_fdc;
    JKGE = curr.JKGE;

    if nargin < 9
        error('export_SAGE:NotEnoughInputs', ...
            ['Expected mdl, dat, bas, prd, ' ...
            'dirres, Q, prf, region, nTheta.']);
    end
    
    if ~isfolder(dirres)
        mkdir(dirres);
    end
    
    K = numel(dat);
    row_id = (1:K)';
    
    % ---------------------
    % 1. Basin master table
    % ---------------------
    basin_id = strings(K,1);
    basin_name = strings(K,1);
    basin_group = local_dat_use(dat);
    if all(isnan(basin_group))
        if isfield(bas,'K_t') ...
                && isfield(bas,'K_e')
            basin_group(1:bas.K_t) = 1;
            basin_group(bas.K_t+1: ...
                min(K,bas.K_t+bas.K_e)) = 2;
        end
    end
    
    for k = 1:K
        basin_id(k) = local_get_basin_id( ...
            dat{k},bas,k);
        basin_name(k) = local_get_basin_name( ...
            dat{k},bas,k);
    end
    
    E = struct();
    E.created = datetime('now');
    E.description = ['SAGE final export. ' ...
        'row_id links all ' ...
        'basin-level arrays and tables.'];
    E.codes.basin_group = ...
        ['1=training basin, ' ...
        '2=evaluation basin'];
    E.codes.split_code = ...
        ['1=training period/output, ' ...
        '2=evaluation period/output, ' ...
        'NaN=unused'];
    
    E.basins = table(row_id,basin_id, ...
        basin_name,basin_group, ...
        'VariableNames',{'row_id', ...
        'basin_id','basin_name', ...
        'basin_group'});
    
    % ----------------------
    % 2. Performance metrics
    % ----------------------
    E.performance = struct();
    E.performance.current = curr;
    E.performance.iteration = prf.iter;
    E.performance.tt = ...
        local_metric_table('tt', ...
        bas,NSE,KGE,S_fdc,JKGE);
    E.performance.te = ...
        local_metric_table('te', ...
        bas,NSE,KGE,S_fdc,JKGE);
    E.performance.et = ...
        local_metric_table('et', ...
        bas,NSE,KGE,S_fdc,JKGE);
    E.performance.ee = ...
        local_metric_table('ee', ...
        bas,NSE,KGE,S_fdc,JKGE);
    
    % ------------------------------
    % 3. Final normalized parameters
    % ------------------------------
    if isempty(nTheta)
        E.theta_n = table(row_id);
    else
        d = size(nTheta,1);
        theta = nan(K,d);
    
        % nTheta is usually d x K
        nk = min(K,size(nTheta,2));
        theta(1:nk,:) = nTheta(:,1:nk)';
    
        varNames = [{'row_id'}, ...
            arrayfun(@(j) ...
            sprintf('theta_%02d',j), ...
            1:d,'UniformOutput',false)];
    
        E.theta_n = array2table([row_id theta], ...
            'VariableNames',varNames);
    end
    
    % -------------------
    % 4. Discharge arrays
    % -------------------
    N = local_max_measured_length(dat);
    
    split_code  = nan(K,N);
    Q_measured  = nan(K,N);
    
    for k = 1:K
    
        y = local_get_measured(dat{k});
        if ~isempty(y)
            n = min(N,numel(y));
            Q_measured(k,1:n) = y(1:n);
        end
    
        idt = local_get_train_indices( ...
            mdl,dat{k},N);
        ide = local_get_eval_indices( ...
            mdl,dat{k},N);
    
        idt = idt(idt >= 1 ...
            & idt <= N);
        ide = ide(ide >= 1 ...
            & ide <= N);
    
        split_code(k,idt) = 1;
        split_code(k,ide) = 2;
    end
    
    Q_simulated = ...
        local_build_full_Qsim(Q,dat,mdl,N);
    
    E.discharge = struct();
    E.discharge.split_code = split_code;
    E.discharge.Q_measured = Q_measured;
    E.discharge.Q_simulated = Q_simulated;
    E.discharge.Q_scenario = Q;
    E.discharge.time = ...
        local_get_time(dat,prd,N);
    
    % -----------
    % 5. Metadata
    % -----------
    E.metadata = struct();
    E.metadata.mdl = mdl;
    E.metadata.bas = bas;
    E.metadata.prd = prd;
    E.metadata.size.K = K;
    E.metadata.size.N = N;
    if isempty(nTheta)
        E.metadata.size.d = 0;
    else
        E.metadata.size.d = size(nTheta,1);
    end
    
    file_mat = fullfile(dirres, ...
        local_export_filename(mdl, ...
        prd,region));
    save(file_mat,'E','-v7.3');

end

% =============
% Local helpers
% =============

% ---------------------------------------------------
function name = local_export_filename(mdl,prd,region)
% ---------------------------------------------------

    model = 'model';
    if isfield(mdl,'name') ...
            && ~isempty(mdl.name)
        model = char(string(mdl.name));
    elseif isfield(mdl,'model') ...
            && ~isempty(mdl.model)
        model = char(string(mdl.model));
    end
    
    dtName = 'dt';
    if isfield(prd,'dt')
        if prd.dt == 1
            dtName = 'daily';
        elseif prd.dt == 24
            dtName = 'hourly';
        elseif prd.dt == 96
            dtName = '15min';
        else
            dtName = sprintf('dt%g',prd.dt);
        end
    end
    
    regName = 'region';
    
    if nargin >= 3 ...
            && ~isempty(region)
        r = lower(strtrim(char(string(region))));
        switch r
            case {'united states', ...
                    'usa','us', ...
                    'conus','camels_us'}
                regName = 'CAMELS_US';
    
            case {'brazil','br','camels_br'}
                regName = 'CAMELS_BR';
    
            case {'great britain', ...
                    'uk','gb','camels_gb'}
                regName = 'CAMELS_GB';
    
            otherwise
                regName = char(string(region));
        end
    end
    
    stamp = char(datetime('now', ...
        'Format','yyyyMMdd_HHmmss'));
    
    name = sprintf('SAGE_export_%s_%s_%s_%s.mat', ...
        regName,model,dtName,stamp);
    
    name = regexprep(name,'[^A-Za-z0-9_.-]','_');
end

% ----------------------------------------
function id = local_get_basin_id(dk,bas,k)
% ----------------------------------------

    fields = {'id_USGS','USGS', ...
        'gauge_id','id_gauge','id'};
    for i = 1:numel(fields)
        f = fields{i};
        if isstruct(dk) ...
                && isfield(dk,f) ...
                && ~isempty(dk.(f))
            id = string(dk.(f));
            return
        end
    end
    for i = 1:numel(fields)
        f = fields{i};
        if isstruct(bas) ... 
                && isfield(bas,f) ...
                && numel(bas.(f)) >= k
            id = string(bas.(f)(k));
            return
        end
    end
    id = string(k);
end

% --------------------------------------------
function name = local_get_basin_name(dk,bas,k)
% --------------------------------------------

    name = "";
    fields = {'gname','name','basin_name'};
    for i = 1:numel(fields)
        f = fields{i};
        if isstruct(dk) ...
                && isfield(dk,f) ...
                && ~isempty(dk.(f))
            name = string(dk.(f));
            return
        end
    end
    for i = 1:numel(fields)
        f = fields{i};
        if isstruct(bas) ...
                && isfield(bas,f) ...
                && numel(bas.(f)) >= k
            name = string(bas.(f)(k));
            return
        end
    end
end

% ---------------------------------
function y = local_get_measured(dk)
% ---------------------------------

    y = [];
    fields = {'y_n','y','Qobs', ...
        'qobs','Q_meas','Q_measured'};
    for i = 1:numel(fields)
        f = fields{i};
        if isstruct(dk) ...
                && isfield(dk,f) ...
                && ~isempty(dk.(f))
            y = double(dk.(f)(:))';
            return
        end
    end
end

% ------------------------------------
function t = local_get_time(dat,prd,N)
% ------------------------------------

    if ~isempty(dat) ...
            && isstruct(dat{1})
        fields = {'t','time', ...
            'date','dates','Date'};
        for i = 1:numel(fields)
            f = fields{i};
            if isfield(dat{1},f) ...
                    && ~isempty(dat{1}.(f))
                t = dat{1}.(f)(:);
                return
            end
        end
    end
    
    if isfield(prd,'ds') ...
            && ~isempty(prd.ds)
        try
            dt0 = datetime(prd.ds, ...
                'InputFormat','dd/MM/yyyy');
            if isfield(prd,'dt')
                switch prd.dt
                    case 1
                        t = dt0 + days(0:N-1);
                    case 24
                        t = dt0 + hours(0:N-1);
                    case 96
                        t = dt0 + minutes(15*(0:N-1));
                    otherwise
                        % General interpretation:
                        % prd.dt = number of samples per day.
                        t = dt0 + days((0:N-1)/prd.dt);
                end
            else
                t = dt0 + days(0:N-1);
            end
        catch
            t = (1:N)';
        end
    else
        t = (1:N)';
    end
end

% ----------------------------------------------
function idt = local_get_train_indices(mdl,dk,N)
% ----------------------------------------------

    idt = [];
    
    if isfield(mdl,'local') ...
            && mdl.local == 1 ...
            && isstruct(dk) ...
            && isfield(dk,'id_train')
        idt = dk.id_train(:)';
    elseif isfield(mdl,'id_train') ...
            && ~isempty(mdl.id_train)
        idt = mdl.id_train(:)';
    end
    
    idt = local_logical_or_index(idt,N);
end

% ---------------------------------------------
function ide = local_get_eval_indices(mdl,dk,N)
% ---------------------------------------------

    ide = [];
    
    if isfield(mdl,'local') ...
            && mdl.local == 1 ...
            && isstruct(dk) ...
            && isfield(dk,'id_eval')
        ide = dk.id_eval(:)';
    elseif isfield(mdl,'id_eval') ...
            && ~isempty(mdl.id_eval)
        ide = mdl.id_eval(:)';
    end
    
    ide = local_logical_or_index(ide,N);
end

% ----------------------------------------
function id = local_logical_or_index(id,N)
% ----------------------------------------

    if isempty(id)
        id = [];
        return
    end
    
    if islogical(id)
        id = find(id);
    
    elseif isnumeric(id) ...
            && numel(id) == N ...
            && all(ismember(unique(id(~isnan(id))),[0 1]))
    
        id = find(id ~= 0);
    
    elseif isnumeric(id) ...
            && numel(id) == 2 ...
            && all(isfinite(id)) ...
            && id(1) >= 1 ...
            && id(2) <= N ...
            && id(2) >= id(1)
    
        id = double(id(1):id(2));
    
    else
        id = double(id(:)');
    end

end

% ----------------------------------------------------
function T = local_metric_table(scen,bas,NSE,KGE,S_fdc,JKGE)
% ----------------------------------------------------

    idx = local_scenario_rows(scen,bas);
    row_id = idx(:);
    
    n = numel(row_id);
    nse = local_metric_values(NSE,scen,n);
    kge = local_metric_values(KGE,scen,n);
    s_fdc = local_metric_values(S_fdc,scen,n);
    jkge = local_metric_values(JKGE,scen,n);
    kge_r = local_component_values(KGE,scen,'r',n);
    kge_alpha = local_component_values(KGE,scen,'alpha',n);
    kge_beta = local_component_values(KGE,scen,'beta',n);
    jkge_M = local_component_values(JKGE,scen,'M',n);
    jkge_V = local_component_values(JKGE,scen,'V',n);
    jkge_C = local_component_values(JKGE,scen,'C',n);
    
    T = table(row_id,nse,kge,kge_r,kge_alpha,kge_beta,s_fdc, ...
        jkge,jkge_M,jkge_V,jkge_C, ...
        'VariableNames',{'row_id', ...
        'NSE','KGE','KGE_r','KGE_alpha','KGE_beta','S_fdc', ...
        'JKGE','JKGE_M','JKGE_V','JKGE_C'});
end

% -----------------------------------------------------------
function v = local_component_values(M,scen,name,n)
% -----------------------------------------------------------

    v = nan(n,1);

    if ~isstruct(M) ...
            || ~isfield(M,'components') ...
            || ~isfield(M.components,scen) ...
            || ~isfield(M.components.(scen),name)
        return
    end

    x = M.components.(scen).(name);
    if isempty(x)
        return
    end

    x = double(x(:));
    m = min(n,numel(x));
    v(1:m) = x(1:m);
end

% ------------------------------------------
function idx = local_scenario_rows(scen,bas)
% ------------------------------------------

    Kt = local_getfield_default(bas,'K_t',0);
    Ke = local_getfield_default(bas,'K_e',0);
    
    switch scen
        case {'tt','te'}
            idx = 1:Kt;
        case {'et','ee'}
            idx = Kt + (1:Ke);
        otherwise
            idx = [];
    end
end

% ----------------------------------------
function v = local_metric_values(M,scen,n)
% ----------------------------------------

    v = nan(n,1);
    
    if isempty(M)
        return
    end
    
    if isstruct(M) ...
            && isfield(M,scen)
        x = M.(scen);
    elseif istable(M) ...
            && any(strcmp(M.Properties.VariableNames,scen))
        x = M.(scen);
    else
        x = M;
    end
    
    if isempty(x)
        return
    end
    
    x = double(x(:));
    
    if numel(x) >= n
        v = x(1:n);
    end
end

% --------------------------------------------
function val = local_getfield_default(S,f,def)
% --------------------------------------------

    if isstruct(S) ...
            && isfield(S,f) ...
            && ~isempty(S.(f))
        val = S.(f);
    else
        val = def;
    end
end

% ------------------------------------------------
function Qsim = local_build_full_Qsim(Q,dat,mdl,N)
% ------------------------------------------------

    K = numel(dat);
    Qsim = nan(K,N);
    
    if isempty(Q) ... 
            || ~isstruct(Q)
        return
    end
    
    % Scenario mapping:
    % tt = training basins     | training split
    % te = training basins     | evaluation split
    % et = evaluation basins   | training split
    % ee = evaluation basins   | evaluation split
    
    idxTrainBas = find(local_dat_use(dat) == 1);
    idxEvalBas = find(local_dat_use(dat) == 2);
    
    Qsim = local_insert_scenario(Qsim, ...
        Q,'tt',idxTrainBas,mdl,dat,N,1);
    Qsim = local_insert_scenario(Qsim, ...
        Q,'te',idxTrainBas,mdl,dat,N,2);
    Qsim = local_insert_scenario(Qsim, ...
        Q,'et',idxEvalBas,mdl,dat,N,1);
    Qsim = local_insert_scenario(Qsim, ...
        Q,'ee',idxEvalBas,mdl,dat,N,2);

end

% ----------------------------------------
function Qsim = local_insert_scenario( ...
    Qsim,Q,scen,basinRows,mdl, ...
    dat,N,periodCode)
% ----------------------------------------

    if ~isfield(Q,scen) ...
            || isempty(Q.(scen))
        return
    end
    
    X = Q.(scen);
    
    for ii = 1:numel(basinRows)
    
        k = basinRows(ii);
    
        q = local_extract_scenario_series( ...
            X,ii,k);
        if isempty(q)
            continue
        end
    
        if periodCode == 1
            id = local_get_train_indices( ...
                mdl,dat{k},N);
        else
            id = local_get_eval_indices( ...
                mdl,dat{k},N);
        end
    
        id = id(id >= 1 ...
            & id <= N);
        n = min(numel(id),numel(q));
    
        Qsim(k,id(1:n)) = q(1:n);
    end

end

% ------------------------------------------------
function q = local_extract_scenario_series(X,ii,k)
% ------------------------------------------------

    q = [];
    
    if iscell(X)
        if numel(X) >= ii ...
                && isnumeric(X{ii})
            q = double(X{ii}(:))';
        elseif numel(X) >= k ...
                && isnumeric(X{k})
            q = double(X{k}(:))';
        end
    
    elseif isnumeric(X)
    
        % Most pproc_SAGE output is time x basin:
        % Q.tt = nTrainTime x K_t
        % Q.te = nEvalTime  x K_t
        % Q.et = nTrainTime x K_e
        % Q.ee = nEvalTime  x K_e
        if size(X,2) >= ii ...
                && size(X,1) > size(X,2)
            q = double(X(:,ii))';
    
        % fallback for basin x time
        elseif size(X,1) >= ii
            q = double(X(ii,:));
    
        elseif size(X,2) >= ii
            q = double(X(:,ii))';
        end
    end

end

% -----------------------------
function u = local_dat_use(dat)
% -----------------------------

    K = numel(dat);
    u = nan(K,1);
    
    for k = 1:K
        if isfield(dat{k},'use') ...
                && ~isempty(dat{k}.use)
            val = dat{k}.use;
    
            if ischar(val) ...
                    || isstring(val)
                s = lower(string(val));
                if any(strcmp(s, ...
                        ["t","train","training"]))
                    u(k) = 1;
                elseif any(strcmp(s, ...
                        ["e","eval","evaluation"]))
                    u(k) = 2;
                end
            else
                u(k) = double(val);
            end
        end
    end

end

% -----------------------------------------
function N = local_max_measured_length(dat)
% -----------------------------------------

    N = 0;
    
    for k = 1:numel(dat)
        y = local_get_measured(dat{k});
        N = max(N,numel(y));
    end

end
