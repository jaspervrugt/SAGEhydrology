function varargout = camels(nTheta,mdl,dat,bas,ode,loss,misc,d,i,dirres)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CAMELS Runs the model for all K watersheds and returns the loss,
% gradients, compact performance metrics, and optional parameter
% attribution metrics.
%
% SYNOPSIS:
%  [ell,G,met] = camels(nTheta,mdl,dat,bas,ode,loss,misc,d,i,dirres)
%  [ell,G,met,At,An] = camels(nTheta,mdl,dat,bas,ode,loss,misc,d,i,dirres)
%  [ell,G,met,At,An,Qfdc] = camels(nTheta,mdl,dat,bas,ode,loss,misc,d, ...
%   i,dirres)
%   nTheta      d x K matrix normalized parameter values for K watersheds
%   mdl         structure with model state/parameter info
%    .model      choice of model
%                 1 hymod
%                 2 hmodel
%                 3 sacsma
%                 4 xinanjiang
%                 5 gr4j
%                 6 hbv
%                 7 gr4jB [analytic routing]
%    .mcode      scalar with numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation sacsma_ode C++
%    .calc       model execution
%      'seq'       sequential execution of watersheds
%      'par'       parallel execution of watersheds using parfor
%      'parfeval'  asynchronous parallel execution using batched parfeval
%    .mode       assessment design
%                 1 = training basins | training period
%                 2 = training basins | evaluation period/mask
%                 3 = training + evaluation basins | training period
%                 4 = training + evaluation basins | evaluation period/mask
%   dat         1 x K cell structure with forcing data for each watershed
%   bas         structure with basin information
%    .K          total number of CAMELS watersheds
%                [671 = all CAMELS / 531 = daily-restricted / 516 = hourly]
%    .K_t        number of training watersheds
%    .K_e        number of evaluation watersheds
%    .K          number of training + evaluation watersheds
%    .r          number of basin attributes
%    .id_t       K_t x 1 vector of training basin indices
%    .id_e       K_e x 1 vector of evaluation basin indices
%    .id_gauge   revised list of gauge basin codes
%   ode         structure with ODE solver settings
%    .InitStep   Initial time step
%    .MaxStep    Maximum time step
%    .MinStep    Minimum time step
%    .RelTol     Relative tolerance
%    .AbsTol     Absolute tolerance
%    .Order      Order
%    .maxiter    Maximum number of iterations
%    .mem        memory storage [=1] of states or not [=0]
%   loss        loss-function settings
%    .fnc        scalar choice of loss function 
%                  1 sum of absolute residuals
%                  2 generalized least squares
%                  3 Nash-Sutcliffe efficiency
%                  4 Kling-Gupta efficiency
%                  5 Huber loss
%                  6 Flow duration curve loss
%                  7 Jawad Kling-Gupta efficiency
%    .n_win      for fnc = 7, moving-average window length in days
%    .method     scalar JKGE benchmark method
%                  1 = moving-average mean
%                  2 = section-wise mean
%                  3 = long-term mean
%                  4 = monthly climatology
%   misc        remaining miscellaneous runtime settings
%    .io         output/logging
%     .plopt      printing individual (=1) or one (=2) progress figure
%     .file       0/1 write runtime+params file
%    .attr       0/1 compute parameter attribution
%    .crr_backend CRR execution backend
%                  'cpp'    central native C++ backend (default)
%                  'matlab' reference crr_model.m backend
%   d           number of parameters
%   i           descent iteration counter
%   dirres      results directory
%
% OUTPUT:
%   ell         1 x K vector with basin-wise training-period loss values
%               ell(1:K_t)   = losses for training basins on train period
%               ell(K_t+1:K) = losses for evaluation basins on train period
%   G           d x K_t matrix gradients of training watersheds
%   met         compact metrics for all K basins
%    .loss.t     training-period losses: SAR, GLS, Huber and RSS
%    .loss.e     evaluation-period losses: SAR, GLS, Huber and RSS
%    .performance.t training-period scores: NSE, KGE and JKGE
%    .performance.e evaluation-period scores: NSE, KGE and JKGE
%   At          OPTIONAL OUTPUT: d x K_t matrix time-weighted parameter
%                attribution values retained for training watersheds only
%   An          OPTIONAL OUTPUT: d x K_t matrix net gradient-based
%                attribution values retained for training watersheds only
%   Qfdc        OPTIONAL OUTPUT: retained discharge series for selected
%               postprocessor basins
%    .id         retained basin indices
%    .gauge       retained basin gauge identifiers
%    .req        requested scenario memberships (.tt,.te,.et,.ee)
%    .qy         1 x nKeep cell array with [q y] for each retained basin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Apr. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ---------------------
    % Unpack basin settings
    % ---------------------
    if ~isstruct(bas) ...
            || ~isfield(bas,'K') ...
            || ~isfield(bas,'K_t')
        error(['camels: ' ...
            'bas must contain ' ...
            'fields K and K_t.']);
    end
    
    K = bas.K;
    K_t = bas.K_t;
    
    % --------------------
    % Check model settings
    % --------------------
    if ~isstruct(mdl) ...
            || ~isfield(mdl,'calc') ...
            || isempty(mdl.calc)
        error(['camels: ' ...
            'mdl.calc must be ''seq'', ''par'', ' ...
            'or ''parfeval''.']);
    end
    
    if ~isfield(mdl,'mode') ...
            || isempty(mdl.mode) ...
            || ~ismember(mdl.mode,[1 2 3 4])
        error(['camels: ' ...
            'mdl.mode must be one of 1,2,3,4.']);
    end
    
    calc = char(lower(string(mdl.calc)));
    
    % -----------------
    % Check misc fields
    % -----------------
    if ~isstruct(loss) ...
            || ~isfield(loss,'fnc') ...
            || isempty(loss.fnc)
        error(['camels: ' ...
            'loss.fnc must be specified.']);
    end
    
    if ~isfield(misc,'io') ...
            || ~isstruct(misc.io) ...
            || ~isfield(misc.io,'file') ...
            || isempty(misc.io.file)
        prt_file = 0;
    else
        prt_file = double(misc.io.file ~= 0);
    end
    
    if ~isfield(misc,'attr') ...
            || isempty(misc.attr)
        attr = 0;
    else
        attr = double(misc.attr ~= 0);
    end
    
    
    % --------------------------------
    % Select CRR execution backend
    % --------------------------------
    % 'cpp'    : central native crr_model_mex via crr_model_cpp
    % 'matlab' : reference crr_model.m implementation
    %
    % The C++ router itself retains the MATLAB fallback for user_model and
    % any execution mode that is not supported by the native backend.
    if ~isfield(misc,'crr_backend') ...
            || isempty(misc.crr_backend)
        crr_backend = 'cpp';
    else
        crr_backend = char(lower(string(misc.crr_backend)));
    end

    if ~ismember(crr_backend,{'cpp','matlab'})
        error(['camels: unknown misc.crr_backend = ''%s''. ' ...
            'Use ''cpp'' or ''matlab''.'],crr_backend);
    end

    % ----------------
    % Check main input
    % ----------------
    if ~isnumeric(d) ...
            || ~isscalar(d) ...
            || d < 1 ...
            || mod(d,1) ~= 0
        error(['camels: ' ...
            'd must be a positive integer.']);
    end
    
    if ~isnumeric(nTheta) ...
            || size(nTheta,1) ~= d ...
            || size(nTheta,2) ~= K
        error(['camels: ' ...
            'nTheta must have size ' ...
            '%d x %d.'],d,K);
    end
    
    if ~iscell(dat) ...
            || numel(dat) ~= K
        error(['camels: ' ...
            'dat must be a cell array ' ...
            'with %d elements.'],K);
    end
    
    if prt_file
        if nargin < 10 ...
                || isempty(dirres)
            error(['camels: ' ...
                'dirres must be ' ...
                'provided when misc.io.file = 1.']);
        end
        if ~isfolder(dirres)
            mkdir(dirres);
        end
    end
    
    % ------------------------
    % Preallocate main outputs
    % ------------------------
    ell = nan(1,K);
    G = nan(d,K);
    
    fields = {'SARt','GLSt','NSEt', ...
              'KGEt','KGE_rt','KGE_alphat','KGE_betat', ...
              'Hubert','RSSt','Dfdct','JKGEt', ...
              'JKGE_Mt','JKGE_Vt','JKGE_Ct', ...
              'SARe','GLSe','NSEe', ...
              'KGEe','KGE_re','KGE_alphae','KGE_betae', ...
              'Hubere','RSSe','Dfdce','JKGEe', ...
              'JKGE_Me','JKGE_Ve','JKGE_Ce'};
    
    vals = repmat({nan(1,K)},1,numel(fields));
    met = cell2struct(vals,fields,2);
    
    if attr
        At = nan(d,K);
        An = nan(d,K);
    else
        At = [];
        An = [];
    end
    
    % -----------------------------------------
    % Determine which basins to retain for FDCs
    % -----------------------------------------
    [keep_idx,req_prt] = ...
        resolve_plot_basin_indices_local(misc,bas);
    
    storeQ = ~isempty(keep_idx);
    
    keep_mask = false(1,K);
    keep_mask(keep_idx) = true;
    
    if storeQ
        Qtmp = cell(1,K);
        Qfdc = struct();
        Qfdc.id = keep_idx(:)';
        Qfdc.gauge = bas.id_gauge(keep_idx);
        Qfdc.req = req_prt;
        Qfdc.qy = [];
    else
        Qtmp = {};
        Qfdc = [];
    end
    
    % -----------------------
    % Prepare runtime logging
    % -----------------------
    if prt_file
        runtime = nan(K,1);
        log_file = fullfile(dirres, ...
            'runtime_params.txt');
    
        if ~isnumeric(i) ...
                || ~isscalar(i) ...
                || ~isfinite(i)
            error(['camels: ' ...
                'iteration counter i must ' ...
                'be a finite scalar.']);
        end
    
        if i == 1
            fid = fopen(log_file,'w');
            if fid < 0
                error(['camels: ' ...
                    'could not open ' ...
                    'runtime log file for ' ...
                    'writing:\n      %s'], ...
                    log_file);
            end
            cObj = onCleanup(@() fclose(fid));
    
            fprintf(fid,'%-6s %-6s %-10s ', ...
                'it','k','runtime');
            for j = 1:d
                fprintf(fid,'%-10s', ...
                    sprintf('theta_%d',j));
            end
            fprintf(fid,'\n');
        end
    end
    
    % ---------------------------------------
    % Check requested execution configuration
    % ---------------------------------------
    req_crr = crr_request(struct('q',true, ...
        'gradient',true, ...
        'metrics',true, ...
        'attribution',logical(attr)));
    switch calc
        case 'seq'
            % ok
        case {'par','parfeval'}
            has_gcp = exist('gcp','file') == 2;
            if ~has_gcp
                fprintf(['      Warning: camels: ' ...
                    'gcp is not available. ' ...
                    'Falling back to sequential.\n']);
                calc = 'seq';
            else
                try
                    p = gcp('nocreate');
                catch
                    p = [];
                end
    
                if isempty(p)
                    fprintf(['      Warning: camels: ' ...
                        'no parallel pool found. ' ...
                        'Falling back to sequential.\n']);
                    calc = 'seq';
                end
            end
        otherwise
            error(['camels: ' ...
                'unknown mdl.calc = ''%s''. ' ...
                'Use ''seq'', ''par'', ' ...
                'or ''parfeval''.'],calc);
    end
    % -------------------
    % Sequential solution
    % -------------------
    switch calc
    
        case 'seq'

            % Optional GUI-yield callback. In sequential mode this runs on
            % the MATLAB client after each completed basin, allowing the
            % top-right GUI label to show live basin progress.
            yieldFcn = [];
            try
                if isfield(misc,'io') ...
                        && isfield(misc.io,'yieldFcn') ...
                        && isa(misc.io.yieldFcn,'function_handle')
                    yieldFcn = misc.io.yieldFcn;
                end
            catch
                yieldFcn = [];
            end

            tYield = tic;
            yieldEvery = 0.10;

            for k = 1:K
                if prt_file
                    t_start = tic;
                end
    
                if attr == 0
                    [ell(k),outk] = run_crr_local(crr_backend, ...
                        nTheta(:,k),mdl, ...
                        dat{k},ode,loss,req_crr);
                    qk = outk.q;
                    G(:,k) = outk.gradient;
                    metk = outk.metrics;
                else
                    [ell(k),outk] = run_crr_local(crr_backend, ...
                        nTheta(:,k),mdl, ...
                        dat{k},ode,loss,req_crr);
                    qk = outk.q;
                    G(:,k) = outk.gradient;
                    metk = outk.metrics;
                    At(:,k) = outk.attribution.total;
                    An(:,k) = outk.attribution.net;
                end
    
                if storeQ ...
                        && keep_mask(k)
                    Qtmp{k} = ...
                        single([qk(:) , ...
                        dat{k}.y_n(:)]);
                end
    
                for f = 1:numel(fields)
                    if isfield(metk,fields{f})
                        met.(fields{f})(k) = ...
                            metk.(fields{f});
                    end
                end
    
                if prt_file
                    runtime(k) = toc(t_start);
                end

                % Report progress to the GUI. The callback itself applies
                % additional renderer throttling, so this remains light.
                if ~isempty(yieldFcn) ...
                        && (toc(tYield) > yieldEvery ...
                        || k == 1 ...
                        || k == K)
                    try
                        yieldFcn(i,k,K);
                    catch
                        % Progress reporting must never interrupt SAGE.
                    end
                    drawnow limitrate
                    tYield = tic;
                end
            end
    
        case 'par'
            % parfor cannot conveniently write into met field-by-field,
            % so store metric arrays first and assemble met afterward.
            template = nan(1,K);
            SARt = template; GLSt = template; NSEt = template;
            KGEt = template; Hubert = template; RSSt = template;
            Dfdct = template;
            KGE_rt = template; KGE_alphat = template;
            KGE_betat = template; JKGEt = template;
            JKGE_Mt = template; JKGE_Vt = template;
            JKGE_Ct = template;
            SARe = template; GLSe = template; NSEe = template;
            KGEe = template; Hubere = template; RSSe = template;
            Dfdce = template;
            KGE_re = template; KGE_alphae = template;
            KGE_betae = template; JKGEe = template;
            JKGE_Me = template; JKGE_Ve = template;
            JKGE_Ce = template;
    
            if attr == 0
                parfor k = 1:K
                    t_start = tic;    
                    [ell(k),outk] = run_crr_local(crr_backend, ...
                        nTheta(:,k),mdl, ...
                        dat{k},ode,loss,req_crr);
                    qk = outk.q;
                    G(:,k) = outk.gradient;
                    metk = outk.metrics;
                    if storeQ ...
                            && keep_mask(k)
                        Qtmp{k} = ...
                            single([qk(:) , ...
                            dat{k}.y_n(:)]);
                    end
    
                    if isfield(metk,'SARt')   
                        SARt(k) = metk.SARt;   
                    end
                    if isfield(metk,'GLSt')   
                        GLSt(k) = metk.GLSt;   
                    end
                    if isfield(metk,'NSEt')   
                        NSEt(k) = metk.NSEt;   
                    end
                    if isfield(metk,'KGEt')   
                        KGEt(k) = metk.KGEt;   
                    end
                    KGE_rt(k) = metk.KGE_rt;
                    KGE_alphat(k) = metk.KGE_alphat;
                    KGE_betat(k) = metk.KGE_betat;
                    if isfield(metk,'Hubert') 
                        Hubert(k) = metk.Hubert; 
                    end
                    if isfield(metk,'RSSt')   
                        RSSt(k) = metk.RSSt;   
                    end
                    Dfdct(k) = metk.Dfdct;
                    if isfield(metk,'JKGEt')   
                        JKGEt(k) = metk.JKGEt;   
                    end
                    JKGE_Mt(k) = metk.JKGE_Mt;
                    JKGE_Vt(k) = metk.JKGE_Vt;
                    JKGE_Ct(k) = metk.JKGE_Ct;
    
                    if isfield(metk,'SARe')   
                        SARe(k) = metk.SARe;   
                    end
                    if isfield(metk,'GLSe')   
                        GLSe(k) = metk.GLSe;   
                    end
                    if isfield(metk,'NSEe')   
                        NSEe(k) = metk.NSEe;   
                    end
                    if isfield(metk,'KGEe')   
                        KGEe(k) = metk.KGEe;   
                    end
                    KGE_re(k) = metk.KGE_re;
                    KGE_alphae(k) = metk.KGE_alphae;
                    KGE_betae(k) = metk.KGE_betae;
                    if isfield(metk,'Hubere') 
                        Hubere(k) = metk.Hubere; 
                    end
                    if isfield(metk,'RSSe')   
                        RSSe(k) = metk.RSSe;   
                    end
                    Dfdce(k) = metk.Dfdce;
                    if isfield(metk,'JKGEe')   
                        JKGEe(k) = metk.JKGEe;   
                    end
                    JKGE_Me(k) = metk.JKGE_Me;
                    JKGE_Ve(k) = metk.JKGE_Ve;
                    JKGE_Ce(k) = metk.JKGE_Ce;
    
                    if prt_file
                        runtime(k) = toc(t_start);
                    end
                end
    
            else
                parfor k = 1:K
    
                    t_start = tic;    
                    [ell(k),outk] = run_crr_local(crr_backend, ...
                        nTheta(:,k),mdl, ...
                        dat{k},ode,loss,req_crr);
                    qk = outk.q;
                    G(:,k) = outk.gradient;
                    metk = outk.metrics;
                    At(:,k) = outk.attribution.total;
                    An(:,k) = outk.attribution.net;
                    if storeQ ...
                            && keep_mask(k)
                        Qtmp{k} = ...
                            single([qk(:) , ...
                            dat{k}.y_n(:)]);
                    end
    
                    if isfield(metk,'SARt')  
                        SARt(k) = metk.SARt;   
                    end
                    if isfield(metk,'GLSt')  
                        GLSt(k) = metk.GLSt;   
                    end
                    if isfield(metk,'NSEt')  
                        NSEt(k) = metk.NSEt;   
                    end
                    if isfield(metk,'KGEt')  
                        KGEt(k) = metk.KGEt;   
                    end
                    KGE_rt(k) = metk.KGE_rt;
                    KGE_alphat(k) = metk.KGE_alphat;
                    KGE_betat(k) = metk.KGE_betat;
                    if isfield(metk,'Hubert')
                        Hubert(k) = metk.Hubert; 
                    end
                    if isfield(metk,'RSSt')  
                        RSSt(k) = metk.RSSt;   
                    end
                    Dfdct(k) = metk.Dfdct;
                    if isfield(metk,'JKGEt')   
                        JKGEt(k) = metk.JKGEt;   
                    end
                    JKGE_Mt(k) = metk.JKGE_Mt;
                    JKGE_Vt(k) = metk.JKGE_Vt;
                    JKGE_Ct(k) = metk.JKGE_Ct;
    
                    if isfield(metk,'SARe')  
                        SARe(k) = metk.SARe;   
                    end
                    if isfield(metk,'GLSe')  
                        GLSe(k) = metk.GLSe;   
                    end
                    if isfield(metk,'NSEe')  
                        NSEe(k) = metk.NSEe;   
                    end
                    if isfield(metk,'KGEe')  
                        KGEe(k) = metk.KGEe;   
                    end
                    KGE_re(k) = metk.KGE_re;
                    KGE_alphae(k) = metk.KGE_alphae;
                    KGE_betae(k) = metk.KGE_betae;
                    if isfield(metk,'Hubere')
                        Hubere(k) = metk.Hubere; 
                    end
                    if isfield(metk,'RSSe')  
                        RSSe(k) = metk.RSSe;   
                    end
                    Dfdce(k) = metk.Dfdce;
                    if isfield(metk,'JKGEe')   
                        JKGEe(k) = metk.JKGEe;   
                    end
                    JKGE_Me(k) = metk.JKGE_Me;
                    JKGE_Ve(k) = metk.JKGE_Ve;
                    JKGE_Ce(k) = metk.JKGE_Ce;
    
                    if prt_file
                        runtime(k) = toc(t_start);
                    end
                end
            end
            values = {SARt,GLSt,NSEt,KGEt, ...
                      KGE_rt,KGE_alphat,KGE_betat, ...
                      Hubert,RSSt,Dfdct,JKGEt, ...
                      JKGE_Mt,JKGE_Vt,JKGE_Ct, ...
                      SARe,GLSe,NSEe,KGEe, ...
                      KGE_re,KGE_alphae,KGE_betae, ...
                      Hubere,RSSe,Dfdce,JKGEe, ...
                      JKGE_Me,JKGE_Ve,JKGE_Ce};
            for f = 1:numel(fields)
                met.(fields{f}) = values{f};
            end
    
        case 'parfeval'
    
            % -----------------------------------------------------------
            % parfeval execution:
            % Submit batches of basins asynchronously and collect results
            % as they finish. This keeps the GUI/client thread responsive
            % between completed batches.
            % -----------------------------------------------------------
            template = nan(1,K);
            SARt = template; GLSt = template; NSEt = template;
            KGEt = template; Hubert = template; RSSt = template;
            Dfdct = template;
            KGE_rt = template; KGE_alphat = template;
            KGE_betat = template; JKGEt = template;
            JKGE_Mt = template; JKGE_Vt = template;
            JKGE_Ct = template;
    
            SARe = template; GLSe = template; NSEe = template;
            KGEe = template; Hubere = template; RSSe = template;
            Dfdce = template;
            KGE_re = template; KGE_alphae = template;
            KGE_betae = template; JKGEe = template;
            JKGE_Me = template; JKGE_Ve = template;
            JKGE_Ce = template;
            if prt_file
                runtime = nan(K,1);
            end
            % Optional GUI-yield callback. This runs only on the client.
            yieldFcn = [];
            try
                if isfield(misc,'io') ...
                        && isfield(misc.io,'yieldFcn') ...
                        && isa(misc.io.yieldFcn,'function_handle')
                    yieldFcn = misc.io.yieldFcn;
                end
            catch
                yieldFcn = [];
            end
            pool = gcp('nocreate');
            if isempty(pool)
                error('camels:parfeval:noPool', ...
                    ['Parallel pool is required ' ...
                    'for parfeval execution.']);
            end
            % Use batches rather than one future per basin.
            % This reduces future-management overhead and closer to parfor
            nBatch = min(K,max(pool.NumWorkers, ...
                4*pool.NumWorkers));
            edges = round(linspace(1,K+1,nBatch+1));
            % Ensure strictly increasing batch edges.
            edges = unique(edges,'stable');
            if edges(1) ~= 1
                edges = [1 edges];
            end
            if edges(end) ~= K+1
                edges = [edges K+1];
            end
            nBatch = numel(edges) - 1;
            futures(nBatch,1) = parallel.FevalFuture;
            for b = 1:nBatch
                ids = edges(b):(edges(b+1)-1);
                storeMask_b = storeQ ...
                    & keep_mask(ids);
                % Outputs:
                % ids, ell_b, G_b, met_b, At_b, An_b, Q_b, runtime_b
                futures(b) = parfeval(pool, ...
                    @run_basin_batch_local,8, ...
                    ids,nTheta(:,ids),mdl,dat(ids),ode, ...
                    loss,attr,storeMask_b,crr_backend);
            end
    
            tYield = tic;
            yieldEvery = 0.10;
            nDone = 0;
    
            for b = 1:nBatch
                [~,ids,ell_b,G_b,met_b,At_b, ...
                    An_b,Q_b,runtime_b] = ...
                    fetchNext(futures);
                nIds = numel(ids);
                ell(ids) = ell_b;
                G(:,ids) = G_b;
                if attr
                    At(:,ids) = At_b;
                    An(:,ids) = An_b;
                end
                if prt_file
                    runtime(ids) = runtime_b;
                end
                % Store retained FDC discharge series
                if storeQ
                    for jj = 1:nIds
                        k = ids(jj);
                        if keep_mask(k)
                            Qtmp{k} = Q_b{jj};
                        end
                    end
                end
                % Copy metric arrays into global locations
                if isfield(met_b,'SARt')
                    SARt(ids) = met_b.SARt;
                end
                if isfield(met_b,'GLSt')
                    GLSt(ids) = met_b.GLSt;
                end
                if isfield(met_b,'NSEt')
                    NSEt(ids) = met_b.NSEt;
                end
                if isfield(met_b,'KGEt')
                    KGEt(ids) = met_b.KGEt;
                end
                KGE_rt(ids) = met_b.KGE_rt;
                KGE_alphat(ids) = met_b.KGE_alphat;
                KGE_betat(ids) = met_b.KGE_betat;
                if isfield(met_b,'Hubert')
                    Hubert(ids) = met_b.Hubert;
                end
                if isfield(met_b,'RSSt')
                    RSSt(ids) = met_b.RSSt;
                end
                if isfield(met_b,'Dfdct')
                    Dfdct(ids) = met_b.Dfdct;
                end
                if isfield(met_b,'JKGEt')
                    JKGEt(ids) = met_b.JKGEt;
                end
                JKGE_Mt(ids) = met_b.JKGE_Mt;
                JKGE_Vt(ids) = met_b.JKGE_Vt;
                JKGE_Ct(ids) = met_b.JKGE_Ct;
                
                if isfield(met_b,'SARe')
                    SARe(ids) = met_b.SARe;
                end
                if isfield(met_b,'GLSe')
                    GLSe(ids) = met_b.GLSe;
                end
                if isfield(met_b,'NSEe')
                    NSEe(ids) = met_b.NSEe;
                end
                if isfield(met_b,'KGEe')
                    KGEe(ids) = met_b.KGEe;
                end
                KGE_re(ids) = met_b.KGE_re;
                KGE_alphae(ids) = met_b.KGE_alphae;
                KGE_betae(ids) = met_b.KGE_betae;
                if isfield(met_b,'Hubere')
                    Hubere(ids) = met_b.Hubere;
                end
                if isfield(met_b,'RSSe')
                    RSSe(ids) = met_b.RSSe;
                end
                if isfield(met_b,'Dfdce')
                    Dfdce(ids) = met_b.Dfdce;
                end
                if isfield(met_b,'JKGEe')
                    JKGEe(ids) = met_b.JKGEe;
                end
                JKGE_Me(ids) = met_b.JKGE_Me;
                JKGE_Ve(ids) = met_b.JKGE_Ve;
                JKGE_Ce(ids) = met_b.JKGE_Ce;
                nDone = nDone + nIds;
                % Let the GUI process dropdown/tab/plot callbacks.
                if toc(tYield) > yieldEvery ...
                        || b == 1 ...
                        || b == nBatch
                    if ~isempty(yieldFcn)
                        try
                            yieldFcn(i,nDone,K);
                        catch
                        end
                    end
                    drawnow limitrate
                    tYield = tic;
                end
            end
            values = {SARt,GLSt,NSEt,KGEt, ...
                      KGE_rt,KGE_alphat,KGE_betat, ...
                      Hubert,RSSt,Dfdct,JKGEt, ...
                      JKGE_Mt,JKGE_Vt,JKGE_Ct, ...
                      SARe,GLSe,NSEe,KGEe, ...
                      KGE_re,KGE_alphae,KGE_betae, ...
                      Hubere,RSSe,Dfdce,JKGEe, ...
                      JKGE_Me,JKGE_Ve,JKGE_Ce};
            for f = 1:numel(fields)
                met.(fields{f}) = values{f};
            end
    end
    
    % -------------------
    % Append runtime file
    % -------------------
    if prt_file
        fid = fopen(log_file,'a');
        if fid < 0
            error(['camels: ' ...
                'could not open ' ...
                'runtime log file for ' ...
                'appending:\n      %s'], ...
                log_file);
        end
        cObj = onCleanup(@() fclose(fid));
    
        for k = 1:K
            fprintf(fid,'%-6d %-6d %-10.4f ', ...
                i,k,runtime(k));
            fprintf(fid,'%-10.4f',nTheta(:,k));
            fprintf(fid,'\n');
        end
    end
    
    % ----------------------------------------------
    % Retain gradients/attribution for training only
    % ----------------------------------------------
    G = G(:,1:K_t);
    
    if attr
        At = At(:,1:K_t);
        An = An(:,1:K_t);
    end
    
    % -----------------------
    % Return requested output
    % -----------------------
    if storeQ
        nq = numel(keep_idx);
        qy_keep = cell(1,nq);
        for j = 1:nq
            kk = keep_idx(j);
            qy_keep{j} = Qtmp{kk};
        end
        Qfdc.qy = qy_keep;
    end
    
    % Present metrics through the public semantic schema. Worker-local
    % arrays remain flat to keep parfor/parfeval aggregation inexpensive.
    met = local_package_metrics(met);
    met.loss.fnc = loss.fnc;
    met.loss.t.objective = ell;

    allout = {ell,G,met,At,An,Qfdc};
    varargout = allout(1:nargout);

end

% =============
% local helpers
% =============
function met = local_package_metrics(flat)
%LOCAL_PACKAGE_METRICS Organize public diagnostics by their meaning.

    met = struct();
    met.loss = struct();
    met.loss.t = struct( ...
        'SAR',flat.SARt, ...
        'GLS',flat.GLSt, ...
        'Huber',flat.Hubert, ...
        'RSS',flat.RSSt);
    met.loss.e = struct( ...
        'SAR',flat.SARe, ...
        'GLS',flat.GLSe, ...
        'Huber',flat.Hubere, ...
        'RSS',flat.RSSe);

    met.performance = struct();
    met.performance.t = struct( ...
        'NSE',flat.NSEt, ...
        'KGE',flat.KGEt, ...
        'KGE_components',struct( ...
            'r',flat.KGE_rt, ...
            'alpha',flat.KGE_alphat, ...
            'beta',flat.KGE_betat), ...
        'D_fdc',flat.Dfdct, ...
        'JKGE',flat.JKGEt, ...
        'JKGE_components',struct( ...
            'M',flat.JKGE_Mt, ...
            'V',flat.JKGE_Vt, ...
            'C',flat.JKGE_Ct));
    met.performance.e = struct( ...
        'NSE',flat.NSEe, ...
        'KGE',flat.KGEe, ...
        'KGE_components',struct( ...
            'r',flat.KGE_re, ...
            'alpha',flat.KGE_alphae, ...
            'beta',flat.KGE_betae), ...
        'D_fdc',flat.Dfdce, ...            
        'JKGE',flat.JKGEe, ...
        'JKGE_components',struct( ...
            'M',flat.JKGE_Me, ...
            'V',flat.JKGE_Ve, ...
            'C',flat.JKGE_Ce));
end

function [keep_idx,req] = ...
    resolve_plot_basin_indices_local(misc,bas)

    req = struct('tt',[],'te',[],'et',[],'ee',[]);
    flds = {'tt','te','et','ee'};
    
    if ~isfield(misc,'plot') ...
            || isempty(misc.plot)
        keep_idx = [];
        return
    end
    
    % ------------------------------------------------
    % Option 0: direct list of basin indices to retain
    % ------------------------------------------------
    if isfield(misc.plot,'k') ...
            && ~isempty(misc.plot.k)
    
        keep_idx = unique(double(misc.plot.k(:)))';
        keep_idx = keep_idx(isfinite(keep_idx) ...
            & keep_idx >= 1 ...
            & keep_idx <= bas.K);
    
        kt = keep_idx(keep_idx <= bas.K_t);
        ke = keep_idx(keep_idx >  bas.K_t);
    
        req.tt = kt;
        req.te = kt;
        req.et = ke;
        req.ee = ke;
    
        % no return
    end
    
    % ---------------------------------------
    % Option 1: already provided as scenarios
    % ---------------------------------------
    if isfield(misc.plot,'kscen') ...
            && isstruct(misc.plot.kscen)
    
        for ii = 1:numel(flds)
            f = flds{ii};
            if isfield(misc.plot.kscen,f) ...
                    && ~isempty(misc.plot.kscen.(f))
                req.(f) = unique(double(misc.plot.kscen.(f)(:)))';
            end
        end
    
    % ----------------------------------
    % Option 2: provided as gauge strings
    % ----------------------------------
    elseif isfield(misc.plot,'gaugescen') ...
            && isstruct(misc.plot.gaugescen)
    
        % Normalize all available basin IDs to 8-character strings
        gauge_all = string(bas.id_gauge);
        % remove Excel-style .0 if present
        gauge_all = regexprep(strtrim(gauge_all), '\.0$', '');   
        % left-pad with zeros to width 8
        gauge_all = compose("%08s", gauge_all);                 
        % compose pads with spaces first
        gauge_all = strrep(gauge_all, " ", "0");                
    
        for ii = 1:numel(flds)
            f = flds{ii};
            if isfield(misc.plot.gaugescen,f) ...
                    && ~isempty(misc.plot.gaugescen.(f))
    
                gauge_req = string(misc.plot.gaugescen.(f));
                gauge_req = regexprep(strtrim(gauge_req), '\.0$', '');
                gauge_req = compose("%08s", gauge_req);
                gauge_req = strrep(gauge_req, " ", "0");
    
                [tf,loc] = ismember(gauge_req,gauge_all);
                req.(f) = unique(loc(tf))';
            end
        end
    end
    
    % ----------------------------------------------------------
    % Convert requested global basin indices to retained storage
    % ----------------------------------------------------------
    req_global = req;
    
    keep_idx = unique([req_global.tt(:); ...
        req_global.te(:); ...
        req_global.et(:); ...
        req_global.ee(:)])';
    
    keep_idx = keep_idx(isfinite(keep_idx) ...
        & keep_idx >= 1 ...
        & keep_idx <= bas.K);
    
    % Now req.* becomes local indices into Qfdc.id / Qfdc.qy.
    % This means Qfdc.req.et may be 11:20 if the eval basins are
    % stored after the train basins in Qfdc.id.
    for ii = 1:numel(flds)
    
        f = flds{ii};
    
        r = req_global.(f)(:)';
        r = r(isfinite(r) ...
            & r >= 1 ...
            & r <= bas.K);
    
        [tf,loc] = ismember(r,keep_idx);
    
        req.(f) = loc(tf);
    end

end

function [ids,ell_b,G_b,met_b,At_b,An_b,Q_b,runtime_b] = ...
    run_basin_batch_local(ids,nTheta_b,mdl,dat_b, ...
    ode,loss,attr,storeMask_b,crr_backend)
%RUN_BASIN_BATCH_LOCAL Worker-side batch of basin model evaluations.
%
% This function runs on a parallel worker. It must not access GUI handles,
% nested GUI callbacks, or the S structure.

    ids = ids(:)';
    nb = numel(ids);

    d = size(nTheta_b,1);

    ell_b = nan(1,nb);
    G_b = nan(d,nb);

    fields = {'SARt','GLSt','NSEt', ...
              'KGEt','KGE_rt','KGE_alphat','KGE_betat', ...
              'Hubert','RSSt','Dfdct','JKGEt', ...
              'JKGE_Mt','JKGE_Vt','JKGE_Ct', ...
              'SARe','GLSe','NSEe', ...
              'KGEe','KGE_re','KGE_alphae','KGE_betae', ...
              'Hubere','RSSe','Dfdce','JKGEe', ...
              'JKGE_Me','JKGE_Ve','JKGE_Ce'};

    vals = repmat({nan(1,nb)},1,numel(fields));
    met_b = cell2struct(vals,fields,2);

    if attr
        At_b = nan(d,nb);
        An_b = nan(d,nb);
    else
        At_b = [];
        An_b = [];
    end

    req_crr = crr_request(struct('q',true,'gradient',true, ...
        'metrics',true,'attribution',logical(attr)));
    Q_b = cell(1,nb);
    runtime_b = nan(1,nb);

    for jj = 1:nb

        t0 = tic;

        try
            [ell_b(jj),outj] = run_crr_local( ...
                crr_backend,nTheta_b(:,jj),mdl, ...
                dat_b{jj},ode,loss,req_crr);
        catch ME
            gauge = local_worker_gauge(dat_b{jj});
            MEbasin = MException( ...
                'camels:BasinEvaluationFailed', ...
                ['Model evaluation failed for basin index %d ' ...
                 '(gauge %s; batch position %d of %d). ' ...
                 'Active parameter vector: %s'], ...
                ids(jj),gauge,jj,nb, ...
                mat2str(nTheta_b(:,jj)',8));
            MEbasin = addCause(MEbasin,ME);
            throw(MEbasin);
        end

        if attr == 0
            qj = outj.q;
            G_b(:,jj) = outj.gradient;
            metj = outj.metrics;
        else
            qj = outj.q;
            G_b(:,jj) = outj.gradient;
            metj = outj.metrics;
            At_b(:,jj) = outj.attribution.total;
            An_b(:,jj) = outj.attribution.net;
        end

        if storeMask_b(jj)
            Q_b{jj} = single([qj(:),dat_b{jj}.y_n(:)]);
        end

        for f = 1:numel(fields)
            if isfield(metj,fields{f})
                met_b.(fields{f})(jj) = metj.(fields{f});
            end
        end

        runtime_b(jj) = toc(t0);

    end

end

function [loss_value,out] = run_crr_local( ...
    crr_backend,x,mdl,dat,ode,loss,request)
%RUN_CRR_LOCAL Dispatch one basin to the requested CRR backend.
%
% Keeping the choice here avoids duplicating backend logic throughout the
% sequential, parfor, and parfeval execution paths.

    if ~local_complete_forcing(dat)
        [loss_value,out] = local_invalid_basin_result( ...
            numel(x),dat,request);
        return
    end

    switch crr_backend
        case 'cpp'
            [loss_value,out] = crr_model_cpp( ...
                x,mdl,dat,ode,loss,request);

        case 'matlab'
            [loss_value,out] = crr_model( ...
                x,mdl,dat,ode,loss,request);

        otherwise
            error('camels:UnknownCRRBackend', ...
                'Unknown CRR backend: %s',crr_backend);
    end
end


function tf = local_complete_forcing(dat)
%LOCAL_COMPLETE_FORCING Require finite P, Ep, and T over the full run.
%
% Conceptual hydrologic models are stateful. A missing forcing value at
% any integration step invalidates that step and every subsequent state;
% it is therefore unsafe to score any part of such a simulation.

    tf = false;
    if ~isstruct(dat) ...
            || ~isfield(dat,'meteo') ...
            || ~isstruct(dat.meteo)
        return
    end
    required = {'P','Ep','T'};
    n = [];
    for k = 1:numel(required)
        name = required{k};
        if ~isfield(dat.meteo,name) ...
                || isempty(dat.meteo.(name))
            return
        end
        value = dat.meteo.(name);
        if ~isnumeric(value) ...
                || ~isvector(value) ...
                || any(~isfinite(value(:)))
            return
        end
        if isempty(n)
            n = numel(value);
        elseif numel(value) ~= n
            return
        end
    end
    tf = true;
end


function [loss_value,out] = local_invalid_basin_result(d,dat,request)
%LOCAL_INVALID_BASIN_RESULT Return a shape-correct all-NaN evaluation.

    loss_value = NaN;
    out = struct();
    nQ = 0;
    if isstruct(dat) ...
            && isfield(dat,'y_n')
        nQ = numel(dat.y_n);
    end
    if request.q
        out.q = nan(nQ,1);
    end
    if request.gradient
        out.gradient = nan(d,1);
    end
    if request.jacobian
        out.jacobian = nan(0,d);
    end
    if ~isempty(request.states)
        out.states = [];
    end
    if request.metrics
        fields = {'SARt','GLSt','NSEt', ...
            'KGEt','KGE_rt','KGE_alphat','KGE_betat', ...
            'Hubert','RSSt','Dfdct','JKGEt', ...
            'JKGE_Mt','JKGE_Vt','JKGE_Ct', ...
            'SARe','GLSe','NSEe', ...
            'KGEe','KGE_re','KGE_alphae','KGE_betae', ...
            'Hubere','RSSe','Dfdce','JKGEe', ...
            'JKGE_Me','JKGE_Ve','JKGE_Ce'};
        out.metrics = cell2struct( ...
            repmat({NaN},size(fields)),fields,2);
    end
    if request.attribution
        out.attribution = struct( ...
            'total',nan(d,1),'net',nan(d,1));
    end
end


function gauge = local_worker_gauge(dat)
%LOCAL_WORKER_GAUGE Return a printable gauge identifier on a worker.

    gauge = 'unknown';

    if ~isstruct(dat) ...
            || ~isfield(dat,'gauge') ...
            || isempty(dat.gauge)
        return
    end

    value = dat.gauge;

    if isnumeric(value) ...
            && isscalar(value)
        gauge = char(string(value));
    elseif ischar(value)
        gauge = strtrim(value);
    elseif isstring(value) ...
            && isscalar(value)
        gauge = char(strtrim(value));
    end

end
