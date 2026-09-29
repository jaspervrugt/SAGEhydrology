function varargout = camels(nTheta,mdl,dat,bas,ode,loss,misc,d,i,dirres)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CAMELS Evaluate all selected basins and aggregate results.
%
%  Runs basin models and returns objective values, gradients, metrics, and
%  optional attribution diagnostics.
%
% SYNOPSIS:
%   [ell,G,met] = camels(nTheta,mdl,dat,bas,ode,loss,misc,d,i,dirres)
%   [ell,G,met,At,An] = camels(...)
%   [ell,G,met,At,An,Qfdc,loss,Aobs] = camels(...)
%
% INPUT ARGUMENTS:
%   nTheta          d-by-K normalized hydrologic parameters
%   mdl             model, execution, and assessment settings
%    .model          hydrologic model identifier
%    .mcode          numerical implementation
%    .calc           sequential or parallel execution
%    .mode           basin/period assessment design
%   dat             basin-specific forcing and observations
%   bas             training/evaluation basin counts and IDs
%   ode             numerical solver settings
%   loss            loss function and named-observation settings
%   misc            runtime and backend settings
%   d               number of hydrologic parameters
%   i               training iteration number
%   dirres          results directory
%
% OUTPUT ARGUMENTS:
%   ell             basin-wise objective values
%   G               hydrologic-parameter gradients
%   met             compact basin-wise performance metrics
%   At              optional attribution values
%   An              optional net attribution values
%   Qfdc            optional flow-duration-curve results
%   loss            updated loss settings
%   Aobs            optional observation attribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Aug. 2026             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % Backend preparation retains the MATLAB fallback for user_model,
    % optional models not compiled into crr_model_mex, and any execution
    % mode that is not supported by the native backend.
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
              'Hubert','RSSt','Dfdct','Dpt','Dlogpt','JKGEt', ...
              'JKGE_Mt','JKGE_Vt','JKGE_Ct', ...
              'SARe','GLSe','NSEe', ...
              'KGEe','KGE_re','KGE_alphae','KGE_betae', ...
              'Hubere','RSSe','Dfdce','Dpe','Dlogpe','JKGEe', ...
              'JKGE_Me','JKGE_Ve','JKGE_Ce'};
    jointFields = local_joint_metric_fields();
    fields = [fields jointFields];
    
    vals = repmat({nan(1,K)},1,numel(fields));
    met = cell2struct(vals,fields,2);
    
    if attr
        At = nan(d,K);
        An = nan(d,K);
        AobsTmp = cell(1,K);
    else
        At = [];
        An = [];
        AobsTmp = {};
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

    loss = local_prepare_loss_normalization( ...
        loss,nTheta,mdl,dat,bas,ode, ...
        crr_backend,calc,i);
    req_crr = local_crr_request(loss,attr);

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
                    qk = local_output_q(outk);
                    G(:,k) = outk.gradient;
                    metk = outk.metrics;
                else
                    [ell(k),outk] = run_crr_local(crr_backend, ...
                        nTheta(:,k),mdl, ...
                        dat{k},ode,loss,req_crr);
                    qk = local_output_q(outk);
                    G(:,k) = outk.gradient;
                    metk = outk.metrics;
                    At(:,k) = outk.attribution.total;
                    An(:,k) = outk.attribution.net;
                    AobsTmp{k} = local_observed_attribution( ...
                        outk.attribution, ...
                        local_observation_names(loss));
                end
    
                if storeQ && ~isempty(qk) ...
                        && keep_mask(k)
                    Qtmp{k} = ...
                        single([qk(:) , ...
                        dat{k}.obs.Q.value(:)]);
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
            Dfdct = template; Dpt = template; Dlogpt = template;
            KGE_rt = template; KGE_alphat = template;
            KGE_betat = template; JKGEt = template;
            JKGE_Mt = template; JKGE_Vt = template;
            JKGE_Ct = template;
            SARe = template; GLSe = template; NSEe = template;
            KGEe = template; Hubere = template; RSSe = template;
            Dfdce = template; Dpe = template; Dlogpe = template;
            KGE_re = template; KGE_alphae = template;
            KGE_betae = template; JKGEe = template;
            JKGE_Me = template; JKGE_Ve = template;
            JKGE_Ce = template;
            joint = nan(numel(jointFields),K);
    
            if attr == 0
                parfor k = 1:K
                    t_start = tic;    
                    [ell(k),outk] = run_crr_local(crr_backend, ...
                        nTheta(:,k),mdl, ...
                        dat{k},ode,loss,req_crr);
                    qk = local_output_q(outk);
                    G(:,k) = outk.gradient;
                    metk = outk.metrics;
                    joint(:,k) = local_joint_metric_vector( ...
                        metk,jointFields);
                    if storeQ && ~isempty(qk) ...
                            && keep_mask(k)
                        Qtmp{k} = ...
                            single([qk(:) , ...
                            dat{k}.obs.Q.value(:)]);
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
                    Dpt(k) = metk.Dpt;
                    Dlogpt(k) = metk.Dlogpt;
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
                    Dpe(k) = metk.Dpe;
                    Dlogpe(k) = metk.Dlogpe;
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
                    qk = local_output_q(outk);
                    G(:,k) = outk.gradient;
                    metk = outk.metrics;
                    joint(:,k) = local_joint_metric_vector( ...
                        metk,jointFields);
                    At(:,k) = outk.attribution.total;
                    An(:,k) = outk.attribution.net;
                    AobsTmp{k} = local_observed_attribution( ...
                        outk.attribution, ...
                        local_observation_names(loss));
                    if storeQ && ~isempty(qk) ...
                            && keep_mask(k)
                        Qtmp{k} = ...
                            single([qk(:) , ...
                            dat{k}.obs.Q.value(:)]);
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
                    Dpt(k) = metk.Dpt;
                    Dlogpt(k) = metk.Dlogpt;
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
                    Dpe(k) = metk.Dpe;
                    Dlogpe(k) = metk.Dlogpe;
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
                      Hubert,RSSt,Dfdct,Dpt,Dlogpt,JKGEt, ...
                      JKGE_Mt,JKGE_Vt,JKGE_Ct, ...
                      SARe,GLSe,NSEe,KGEe, ...
                      KGE_re,KGE_alphae,KGE_betae, ...
                      Hubere,RSSe,Dfdce,Dpe,Dlogpe,JKGEe, ...
                      JKGE_Me,JKGE_Ve,JKGE_Ce};
            for f = 1:numel(values)
                met.(fields{f}) = values{f};
            end
            for f = 1:numel(jointFields)
                met.(jointFields{f}) = joint(f,:);
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
            Dfdct = template; Dpt = template; Dlogpt = template;
            KGE_rt = template; KGE_alphat = template;
            KGE_betat = template; JKGEt = template;
            JKGE_Mt = template; JKGE_Vt = template;
            JKGE_Ct = template;
    
            SARe = template; GLSe = template; NSEe = template;
            KGEe = template; Hubere = template; RSSe = template;
            Dfdce = template; Dpe = template; Dlogpe = template;
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
                % ids, ell_b, G_b, met_b, At_b, An_b, Aobs_b, Q_b,
                % runtime_b
                futures(b) = parfeval(pool, ...
                    @run_basin_batch_local,9, ...
                    ids,nTheta(:,ids),mdl,dat(ids),ode, ...
                    loss,attr,storeMask_b,crr_backend);
            end
    
            tYield = tic;
            yieldEvery = 0.10;
            nDone = 0;
    
            for b = 1:nBatch
                [~,ids,ell_b,G_b,met_b,At_b, ...
                    An_b,Aobs_b,Q_b,runtime_b] = ...
                    fetchNext(futures);
                nIds = numel(ids);
                ell(ids) = ell_b;
                G(:,ids) = G_b;
                if attr
                    At(:,ids) = At_b;
                    An(:,ids) = An_b;
                    AobsTmp(ids) = Aobs_b;
                end
                if prt_file
                    runtime(ids) = runtime_b;
                end
                for jf = 1:numel(jointFields)
                    if isfield(met_b,jointFields{jf})
                        met.(jointFields{jf})(ids) = ...
                            met_b.(jointFields{jf});
                    end
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
                if isfield(met_b,'Dpt')
                    Dpt(ids) = met_b.Dpt;
                end
                if isfield(met_b,'Dlogpt')
                    Dlogpt(ids) = met_b.Dlogpt;
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
                if isfield(met_b,'Dpe')
                    Dpe(ids) = met_b.Dpe;
                end
                if isfield(met_b,'Dlogpe')
                    Dlogpe(ids) = met_b.Dlogpe;
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
                      Hubert,RSSt,Dfdct,Dpt,Dlogpt,JKGEt, ...
                      JKGE_Mt,JKGE_Vt,JKGE_Ct, ...
                      SARe,GLSe,NSEe,KGEe, ...
                      KGE_re,KGE_alphae,KGE_betae, ...
                      Hubere,RSSe,Dfdce,Dpe,Dlogpe,JKGEe, ...
                      JKGE_Me,JKGE_Ve,JKGE_Ce};
            for f = 1:numel(values)
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
        Aobs = local_package_observation_attribution( ...
            AobsTmp,local_observation_names(loss),d,K_t);
    else
        Aobs = [];
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
    met = local_package_metrics( ...
        met,local_observation_names(loss));
    met.info = struct('loss_fnc',loss.fnc, ...
        'observed',{cellstr(local_observation_names(loss))});
    if isfield(loss,'normalization')
        met.info.normalization = loss.normalization;
    end
    met.total.loss.t = ell;
    if attr
        met.attribution = struct( ...
            'total',At,'net',An,'observed',Aobs);
    end

    allout = {ell,G,met,At,An,Qfdc,loss,Aobs};
    varargout = allout(1:nargout);

end

% =============
% local helpers
% =============
function met = local_package_metrics(flat,names)
%LOCAL_PACKAGE_METRICS Organize public diagnostics by their meaning.

    met = struct();
    qDiagnostic.t = struct( ...
        'SAR',flat.SARt, ...
        'GLS',flat.GLSt, ...
        'Huber',flat.Hubert, ...
        'RSS',flat.RSSt);
    qDiagnostic.e = struct( ...
        'SAR',flat.SARe, ...
        'GLS',flat.GLSe, ...
        'Huber',flat.Hubere, ...
        'RSS',flat.RSSe);

    qPerformance.t = struct( ...
        'NSE',flat.NSEt, ...
        'KGE',flat.KGEt, ...
        'KGE_components',struct( ...
            'r',flat.KGE_rt, ...
            'alpha',flat.KGE_alphat, ...
            'beta',flat.KGE_betat), ...
        'D_fdc',flat.Dfdct, ...
        'D_p',flat.Dpt, ...
        'D_logp',flat.Dlogpt, ...
        'JKGE',flat.JKGEt, ...
        'JKGE_components',struct( ...
            'M',flat.JKGE_Mt, ...
            'V',flat.JKGE_Vt, ...
            'C',flat.JKGE_Ct));
    qPerformance.e = struct( ...
        'NSE',flat.NSEe, ...
        'KGE',flat.KGEe, ...
        'KGE_components',struct( ...
            'r',flat.KGE_re, ...
            'alpha',flat.KGE_alphae, ...
            'beta',flat.KGE_betae), ...
        'D_fdc',flat.Dfdce, ...            
        'D_p',flat.Dpe, ...
        'D_logp',flat.Dlogpe, ...
        'JKGE',flat.JKGEe, ...
        'JKGE_components',struct( ...
            'M',flat.JKGE_Me, ...
            'V',flat.JKGE_Ve, ...
            'C',flat.JKGE_Ce));

    met.variable = struct();
    met.total = struct('loss',struct( ...
        't',flat.JointLtot_t,'e',flat.JointLtot_e), ...
        'contribution',struct('t',struct(),'e',struct()));
    for j = 1:numel(names)
        name = char(names(j));
        met.total.contribution.t.(name) = ...
            flat.(['JointC' name 't']);
        met.total.contribution.e.(name) = ...
            flat.(['JointC' name 'e']);
        met.variable.(name) = struct();
        met.variable.(name).loss = struct( ...
            't',flat.(['JointL' name 't']), ...
            'e',flat.(['JointL' name 'e']));
        met.variable.(name).performance = struct( ...
            't',local_named_performance(flat,name,'t'), ...
            'e',local_named_performance(flat,name,'e'));
        if strcmp(name,'Q')
            met.variable.Q.diagnostic = qDiagnostic;
            if isscalar(names)
                % The Q-only backend returns the established compact
                % discharge metrics. Named Q_* worker fields are populated
                % only when the multi-observation evaluator is active.
                met.variable.Q.performance = qPerformance;
            end
        end
    end
end

function out = local_named_performance(flat,name,tag)
%LOCAL_NAMED_PERFORMANCE Package named performance metrics for one period.

    out = struct( ...
        'NSE',flat.([name '_NSE' tag]), ...
        'KGE',flat.([name '_KGE' tag]), ...
        'KGE_components',struct( ...
            'r',flat.([name '_KGE_r' tag]), ...
            'alpha',flat.([name '_KGE_alpha' tag]), ...
            'beta',flat.([name '_KGE_beta' tag])), ...
        'D_fdc',flat.([name '_Dfdc' tag]), ...
        'D_p',flat.([name '_Dp' tag]), ...
        'D_logp',flat.([name '_Dlogp' tag]), ...
        'JKGE',flat.([name '_JKGE' tag]), ...
        'JKGE_components',struct( ...
            'M',flat.([name '_JKGE_M' tag]), ...
            'V',flat.([name '_JKGE_V' tag]), ...
            'C',flat.([name '_JKGE_C' tag])));
end

function [keep_idx,req] = ...
    resolve_plot_basin_indices_local(misc,bas)
%RESOLVE_PLOT_BASIN_INDICES_LOCAL Resolve retained plotting-basin requests.

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
    
    % Gauge IDs are canonical. Numeric indices are only a legacy fallback
    % because data-quality screening can renumber the active population.
    if isfield(misc.plot,'gaugescen') ...
            && isstruct(misc.plot.gaugescen)
        gauge_all = local_plot_basin_key(bas.id_gauge);
        for ii = 1:numel(flds)
            f = flds{ii};
            if isfield(misc.plot.gaugescen,f) ...
                    && ~isempty(misc.plot.gaugescen.(f))
                gauge_req = local_plot_basin_key( ...
                    misc.plot.gaugescen.(f));
                [tf,loc] = ismember(gauge_req,gauge_all);
                req.(f) = unique(loc(tf),'stable')';
            end
        end
    elseif isfield(misc.plot,'kscen') ...
            && isstruct(misc.plot.kscen)
        for ii = 1:numel(flds)
            f = flds{ii};
            if isfield(misc.plot.kscen,f) ...
                    && ~isempty(misc.plot.kscen.(f))
                req.(f) = unique(double(misc.plot.kscen.(f)(:)))';
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

function key = local_plot_basin_key(values)
% Canonical, region-neutral key shared by GUI-selected and basin IDs.
    key = upper(strip(string(values(:))));
    key = regexprep(key,'\.0$','');
    key = regexprep(key,'[^A-Z0-9]','');
    numericOnly = ~cellfun('isempty', ...
        regexp(cellstr(key),'^[0-9]+$','once'));
    key(numericOnly) = regexprep(key(numericOnly),'^0+(?=[0-9])','');
end

function [ids,ell_b,G_b,met_b,At_b,An_b,Aobs_b,Q_b,runtime_b] = ...
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
              'Hubert','RSSt','Dfdct','Dpt','Dlogpt','JKGEt', ...
              'JKGE_Mt','JKGE_Vt','JKGE_Ct', ...
              'SARe','GLSe','NSEe', ...
              'KGEe','KGE_re','KGE_alphae','KGE_betae', ...
              'Hubere','RSSe','Dfdce','Dpe','Dlogpe','JKGEe', ...
              'JKGE_Me','JKGE_Ve','JKGE_Ce'};
    fields = [fields local_joint_metric_fields()];

    vals = repmat({nan(1,nb)},1,numel(fields));
    met_b = cell2struct(vals,fields,2);

    if attr
        At_b = nan(d,nb);
        An_b = nan(d,nb);
        Aobs_b = cell(1,nb);
    else
        At_b = [];
        An_b = [];
        Aobs_b = {};
    end

    req_crr = local_crr_request(loss,attr);
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
            qj = local_output_q(outj);
            G_b(:,jj) = outj.gradient;
            metj = outj.metrics;
        else
            qj = local_output_q(outj);
            G_b(:,jj) = outj.gradient;
            metj = outj.metrics;
            At_b(:,jj) = outj.attribution.total;
            An_b(:,jj) = outj.attribution.net;
            Aobs_b{jj} = local_observed_attribution( ...
                outj.attribution,local_observation_names(loss));
        end

        if storeMask_b(jj) && ~isempty(qj)
            Q_b{jj} = single([qj(:),dat_b{jj}.obs.Q.value(:)]);
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

    [loss_value,out] = local_combine_observations( ...
        loss_value,out,dat,mdl,loss,request);
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
            && isfield(dat,'obs') && isfield(dat.obs,'Q') ...
            && isfield(dat.obs.Q,'value')
        nQ = numel(dat.obs.Q.value);
    end
    if request.q
        out.q = nan(nQ,1);
    end
    if ~isempty(request.obs)
        out.obs = struct();
        for j = 1:numel(request.obs)
            name = char(request.obs(j));
            n = nQ;
            if isfield(dat,'obs') && isfield(dat.obs,name) ...
                    && isfield(dat.obs.(name),'value')
                n = numel(dat.obs.(name).value);
            end
            out.obs.(name) = nan(n,1);
        end
    end
    if ~isempty(request.jac)
        out.jac = struct();
        for j = 1:numel(request.jac)
            name = char(request.jac(j));
            n = nQ;
            if isfield(dat,'obs') && isfield(dat.obs,name) ...
                    && isfield(dat.obs.(name),'value')
                n = numel(dat.obs.(name).value);
            end
            out.jac.(name) = nan(n,d);
        end
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
            'Hubert','RSSt','Dfdct','Dpt','Dlogpt','JKGEt', ...
            'JKGE_Mt','JKGE_Vt','JKGE_Ct', ...
            'SARe','GLSe','NSEe', ...
            'KGEe','KGE_re','KGE_alphae','KGE_betae', ...
            'Hubere','RSSe','Dfdce','Dpe','Dlogpe','JKGEe', ...
            'JKGE_Me','JKGE_Ve','JKGE_Ce'};
        fields = [fields local_joint_metric_fields()];
        out.metrics = cell2struct( ...
            repmat({NaN},size(fields)),fields,2);
    end
    if request.attribution ...
            || local_wants_joint_attribution(request)
        out.attribution = struct( ...
            'total',nan(d,1),'net',nan(d,1));
        if local_wants_joint_attribution(request)
            out.attribution.observed = struct();
            for j = 1:numel(request.obs)
                name = char(request.obs(j));
                out.attribution.observed.(name) = ...
                    local_empty_attribution(d,NaN);
            end
        end
    end
end


function q = local_output_q(out)
%LOCAL_OUTPUT_Q Return Q only when the model explicitly supplied it.

    q = [];
    if isstruct(out) && isfield(out,'q')
        q = out.q;
    elseif isstruct(out) && isfield(out,'obs') ...
            && isfield(out.obs,'Q')
        q = out.obs.Q;
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

function request = local_crr_request(loss,attr)
%LOCAL_CRR_REQUEST Request legacy or named multi-observation results.

    names = local_observation_names(loss);
    named = numel(names) > 1 || names(1) ~= "Q";
    request = struct('q',any(names=="Q"),'gradient',true, ...
        'metrics',true, ...
        'attribution',logical(attr) && ~named);

    if named
        request.obs = names;
        request.jac = names;
    end

    request = crr_request(request);
    request.joint_attribution = logical(attr) && named;
end

function names = local_observation_names(loss)
%LOCAL_OBSERVATION_NAMES Normalize selected loss observation names.

    names = "Q";
    if isfield(loss,'observed') ...
            && ~isempty(loss.observed)
        names = upper(strtrim(string(loss.observed(:))));
        names = names(strlength(names) > 0);
    end
    names = unique(names,'stable');
    if isempty(names)
        names = "Q";
    end
    unknown = setdiff(names,["Q";"SWE";"SM"]);
    if ~isempty(unknown)
        error('camels:UnknownObservation', ...
            'Unknown loss observation(s): %s.', ...
            strjoin(cellstr(unknown),', '));
    end
end

function [lossValue,out] = local_combine_observations( ...
    lossValue,out,dat,mdl,loss,request)
%LOCAL_COMBINE_OBSERVATIONS Form one weighted basin loss and gradient.

    names = local_observation_names(loss);
    if isscalar(names) ...
            && names == "Q"
        return
    end
    if ~isfield(out,'obs') ...
            || ~isfield(out,'jac')
        error('camels:MissingNamedResults', ...
            ['The selected model backend did not return the named ' ...
             'observations and Jacobians required by the loss.']);
    end

    n = numel(names);
    weights = local_loss_vector(loss,'weight', ...
        names,ones(n,1)/n);
    scales = local_loss_scales(loss,names);
    baseCoefficients = weights./scales;
    hasCounts = isfield(loss,'normalization') ...
        && isfield(loss.normalization,'n_basin') ...
        && isfield(loss.normalization,'n_total');
    if hasCounts
        counts = local_loss_vector(loss.normalization, ...
            'n_basin',names,ones(n,1));
        baseCoefficients = baseCoefficients ...
            .* double(loss.normalization.n_total) ...
            ./ counts;
    end

    periods = {'train','eval'};
    tags = {'t','e'};
    result = struct();
    for p = 1:2
        values = nan(n,1);
        gradients = cell(n,1);
        details = cell(n,1);
        available = false(n,1);

        for j = 1:n
            name = char(names(j));
            gradients{j} = nan(size(out.gradient));
            if isfield(out.obs,name) ...
                    && isfield(out.jac,name)
                [values(j),gradients{j},details{j}] = ...
                    observation_loss(names(j),out.obs.(name), ...
                    out.jac.(name),dat,loss,periods{p});
                available(j) = details{j}.available;
            end
        end

        coefficients = baseCoefficients;
        if ~hasCounts && any(available)
            activeWeight = weights;
            activeWeight(~available) = 0;
            activeWeight = activeWeight/sum(activeWeight);
            coefficients = activeWeight./scales;
        end

        total = NaN;
        contribution = nan(n,1);
        if any(available)
            contribution(available) = coefficients(available) ...
                .* values(available);
            total = sum(contribution(available));
        end
        result.(tags{p}) = struct( ...
            'value',values,'gradient',{gradients}, ...
            'detail',{details},'available',available, ...
            'coefficient',coefficients, ...
            'contribution',contribution,'total',total);
    end

    lossValue = result.t.total;
    if ~any(result.t.available)
        out.gradient(:) = NaN;
    else
        totalGradient = zeros(size(out.gradient));
        for j = find(result.t.available(:))'
            totalGradient = totalGradient ...
                + result.t.coefficient(j) ...
                * result.t.gradient{j};
        end
        out.gradient = totalGradient;
    end

    out.loss = struct();
    out.gradient_observed = struct();
    for j = 1:n
        name = char(names(j));
        out.loss.(name) = result.t.value(j);
        out.gradient_observed.(name) = result.t.gradient{j};
    end
    out.loss_total = lossValue;
    out.loss_weight = local_named_struct(names,weights);
    out.loss_coefficient = local_named_struct( ...
        names,result.t.coefficient);
    out.loss_scale = local_named_struct(names,scales);
    out.loss_available = local_named_struct( ...
        names,result.t.available);

    if ~isfield(out,'metrics') ...
            || ~isstruct(out.metrics)
        out.metrics = struct();
    end
    out.metrics = local_store_joint_metrics( ...
        out.metrics,names,result);

    if local_wants_joint_attribution(request)
        out.attribution = local_joint_attribution( ...
            out,result,names,mdl);
    end
end

function tf = local_wants_joint_attribution(request)
%LOCAL_WANTS_JOINT_ATTRIBUTION Test the CAMELS-side attribution flag.

    tf = isfield(request,'joint_attribution') ...
        && logical(request.joint_attribution);
end

function attribution = local_joint_attribution( ...
    out,result,names,mdl)
%LOCAL_JOINT_ATTRIBUTION Compute aggregate and named attributions.

    d = numel(out.gradient);
    s = mdl.th_max(:)-mdl.th_min(:);
    if numel(s) ~= d
        error('camels:AttributionParameterSize', ...
            ['Parameter bounds do not match the ' ...
             'multiple-observation gradient.']);
    end

    observed = struct();
    contribution = [];
    for j = 1:numel(names)
        name = char(names(j));
        if ~result.t.available(j)
            observed.(name) = local_empty_attribution(d, ...
                result.t.coefficient(j));
            continue
        end
        J = double(out.jac.(name));
        delta = result.t.detail{j}.delta(:);
        if numel(delta) ~= size(J,1)
            error('camels:AttributionSensitivitySize', ...
                ['Loss sensitivity and Jacobian lengths ' ...
                 'differ for %s.'],name);
        end
        [At,An] = sage_attribution( ...
            J,delta,mdl,result.t.gradient{j});
        coefficient = result.t.coefficient(j);
        observed.(name) = struct( ...
            'total',At,'net',An, ...
            'weighted_total',abs(coefficient)*At, ...
            'weighted_net',abs(coefficient)*An, ...
            'coefficient',coefficient);
        C = coefficient*(J.*delta);
        if isempty(contribution)
            contribution = zeros(size(C));
        elseif ~isequal(size(contribution),size(C))
            error('camels:AttributionTrajectorySize', ...
                ['Named observation trajectories must have ' ...
                 'equal lengths for aggregate attribution.']);
        end
        contribution = contribution+C;
    end

    if isempty(contribution) || any(~isfinite(out.gradient))
        At = nan(d,1);
        An = nan(d,1);
    else
        At = sum(abs(contribution.*s.'),1).';
        An = abs(s.*out.gradient);
    end
    attribution = struct( ...
        'total',At,'net',An,'observed',observed);
end

function value = local_empty_attribution(d,coefficient)
%LOCAL_EMPTY_ATTRIBUTION Return one unavailable named attribution.

    value = struct( ...
        'total',nan(d,1),'net',nan(d,1), ...
        'weighted_total',nan(d,1), ...
        'weighted_net',nan(d,1), ...
        'coefficient',coefficient);
end

function observed = local_observed_attribution(attribution,names)
%LOCAL_OBSERVED_ATTRIBUTION Return a consistent named attribution.

    if isfield(attribution,'observed')
        observed = attribution.observed;
        return
    end
    observed = struct();
    if isscalar(names)
        name = char(names(1));
        observed.(name) = struct( ...
            'total',attribution.total, ...
            'net',attribution.net, ...
            'weighted_total',attribution.total, ...
            'weighted_net',attribution.net, ...
            'coefficient',1);
    end
end

function Aobs = local_package_observation_attribution( ...
    values,names,d,Kt)
%LOCAL_PACKAGE_OBSERVATION_ATTRIBUTION Assemble d x K named matrices.

    Aobs = struct();
    for j = 1:numel(names)
        name = char(names(j));
        Aobs.(name) = struct( ...
            'total',nan(d,Kt), ...
            'net',nan(d,Kt), ...
            'weighted_total',nan(d,Kt), ...
            'weighted_net',nan(d,Kt), ...
            'coefficient',nan(1,Kt));
    end
    for k = 1:Kt
        if isempty(values{k})
            continue
        end
        for j = 1:numel(names)
            name = char(names(j));
            if ~isfield(values{k},name)
                continue
            end
            item = values{k}.(name);
            Aobs.(name).total(:,k) = item.total;
            Aobs.(name).net(:,k) = item.net;
            Aobs.(name).weighted_total(:,k) = ...
                item.weighted_total;
            Aobs.(name).weighted_net(:,k) = ...
                item.weighted_net;
            Aobs.(name).coefficient(k) = item.coefficient;
        end
    end
end

function met = local_store_joint_metrics(met,names,result)
%LOCAL_STORE_JOINT_METRICS Flatten joint diagnostics for aggregation.

    fields = local_joint_metric_fields();
    for f = 1:numel(fields)
        met.(fields{f}) = NaN;
    end

    tags = {'t','e'};
    for p = 1:2
        tag = tags{p};
        R = result.(tag);
        met.(['JointLtot_' tag]) = R.total;
        for j = 1:numel(names)
            name = char(names(j));
            met.(['JointL' name tag]) = R.value(j);
            met.(['JointC' name tag]) = R.contribution(j);
            if ~isempty(R.detail{j})
                met.([name '_NSE' tag]) = R.detail{j}.NSE;
                met.([name '_KGE' tag]) = R.detail{j}.KGE;
                met.([name '_KGE_r' tag]) = ...
                    R.detail{j}.KGE_r;
                met.([name '_KGE_alpha' tag]) = ...
                    R.detail{j}.KGE_alpha;
                met.([name '_KGE_beta' tag]) = ...
                    R.detail{j}.KGE_beta;
                met.([name '_Dfdc' tag]) = R.detail{j}.D_fdc;
                met.([name '_Dp' tag]) = R.detail{j}.D_p;
                met.([name '_Dlogp' tag]) = R.detail{j}.D_logp;
                met.([name '_JKGE' tag]) = R.detail{j}.JKGE;
                met.([name '_JKGE_M' tag]) = R.detail{j}.JKGE_M;
                met.([name '_JKGE_V' tag]) = R.detail{j}.JKGE_V;
                met.([name '_JKGE_C' tag]) = R.detail{j}.JKGE_C;
            end
        end
    end
end

function fields = local_joint_metric_fields()
%LOCAL_JOINT_METRIC_FIELDS Flat worker fields for joint reporting.

    fields = {'JointLQt','JointLSWEt','JointLSMt','JointLtot_t', ...
        'JointCQt','JointCSWEt','JointCSMt', ...
        'JointLQe','JointLSWEe','JointLSMe','JointLtot_e', ...
        'JointCQe','JointCSWEe','JointCSMe', ...
        'Q_NSEt','Q_KGEt','Q_KGE_rt','Q_KGE_alphat', ...
        'Q_KGE_betat','Q_Dfdct','Q_Dpt','Q_Dlogpt', ...
        'Q_JKGEt','Q_JKGE_Mt','Q_JKGE_Vt','Q_JKGE_Ct', ...
        'Q_NSEe','Q_KGEe','Q_KGE_re','Q_KGE_alphae', ...
        'Q_KGE_betae','Q_Dfdce','Q_Dpe','Q_Dlogpe', ...
        'Q_JKGEe','Q_JKGE_Me','Q_JKGE_Ve','Q_JKGE_Ce', ...
        'SWE_NSEt','SWE_KGEt','SWE_KGE_rt','SWE_KGE_alphat', ...
        'SWE_KGE_betat','SWE_Dfdct', ...
        'SWE_Dpt','SWE_Dlogpt','SWE_JKGEt', ...
        'SWE_JKGE_Mt','SWE_JKGE_Vt','SWE_JKGE_Ct', ...
        'SWE_NSEe','SWE_KGEe','SWE_KGE_re','SWE_KGE_alphae', ...
        'SWE_KGE_betae','SWE_Dfdce','SWE_Dpe', ...
        'SWE_Dlogpe','SWE_JKGEe','SWE_JKGE_Me', ...
        'SWE_JKGE_Ve','SWE_JKGE_Ce', ...
        'SM_NSEt','SM_KGEt','SM_KGE_rt','SM_KGE_alphat', ...
        'SM_KGE_betat','SM_Dfdct','SM_Dpt','SM_Dlogpt', ...
        'SM_JKGEt','SM_JKGE_Mt','SM_JKGE_Vt','SM_JKGE_Ct', ...
        'SM_NSEe','SM_KGEe','SM_KGE_re','SM_KGE_alphae', ...
        'SM_KGE_betae','SM_Dfdce','SM_Dpe','SM_Dlogpe', ...
        'SM_JKGEe','SM_JKGE_Me','SM_JKGE_Ve','SM_JKGE_Ce'};
end

function values = local_joint_metric_vector(met,fields)
%LOCAL_JOINT_METRIC_VECTOR Collect flat joint fields for PARFOR.

    values = nan(numel(fields),1);
    for f = 1:numel(fields)
        if isfield(met,fields{f})
            values(f) = met.(fields{f});
        end
    end
end

function values = local_loss_scales(loss,names)
%LOCAL_LOSS_SCALES Return fixed positive scales for named losses.

    values = ones(numel(names),1);
    if ~isfield(loss,'normalization') ...
            || isempty(loss.normalization)
        return
    end
    normalization = loss.normalization;
    if isstruct(normalization) ...
            && isfield(normalization,'scale')
        values = local_loss_vector( ...
            normalization,'scale',names,values);
    end
    if any(~isfinite(values) | values <= 0)
        error('camels:InvalidLossScale', ...
            'Every multiple-observation loss scale must be positive.');
    end
end

function values = local_loss_vector(S,field,names,defaultValue)
%LOCAL_LOSS_VECTOR Read a numeric or named loss setting.

    values = double(defaultValue(:));
    if ~isstruct(S) ...
            || ~isfield(S,field) ...
            || isempty(S.(field))
        return
    end
    source = S.(field);
    if isnumeric(source)
        source = double(source(:));
        if isscalar(source)
            values(:) = source;
        elseif numel(source) == numel(names)
            values = source;
        else
            error('camels:LossSettingSize', ...
                'loss.%s must be scalar or match loss.observed.',field);
        end
    elseif isstruct(source)
        for j = 1:numel(names)
            name = char(names(j));
            if isfield(source,name)
                values(j) = double(source.(name));
            end
        end
    else
        error('camels:InvalidLossSetting', ...
            'loss.%s must be numeric or a named structure.',field);
    end
    if strcmp(field,'weight')
        if any(~isfinite(values) | values < 0) ...
                || sum(values) <= 0
            error('camels:InvalidLossWeight', ...
                ['Multiple-observation loss weights must be ' ...
                 'finite and nonnegative.']);
        end
        values = values/sum(values);
    elseif any(~isfinite(values) | values <= 0)
        error('camels:InvalidLossScale', ...
            ['Multiple-observation loss scales must be ' ...
             'finite and positive.']);
    end
end

function S = local_named_struct(names,values)
%LOCAL_NAMED_STRUCT Store one scalar under each observation name.

    S = struct();
    for j = 1:numel(names)
        S.(char(names(j))) = values(j);
    end
end

function loss = local_prepare_loss_normalization( ...
    loss,nTheta,mdl,dat,bas,ode,crrBackend,calc,iteration)
%LOCAL_PREPARE_LOSS_NORMALIZATION Establish fixed observable scales.

    names = local_observation_names(loss);
    if isscalar(names) ...
            && names == "Q"
        return
    end
    if isscalar(names)
        loss.normalization = struct( ...
            'method','none', ...
            'scale',local_named_struct(names,1), ...
            'floor',1e-12, ...
            'observed',cellstr(names));
        return
    end

    method = "auto";
    if isfield(loss,'normalization') ...
            && isstruct(loss.normalization) ...
            && isfield(loss.normalization,'method') ...
            && ~isempty(loss.normalization.method)
        method = lower(string(loss.normalization.method));
    elseif isfield(loss,'normalization') ...
            && isstruct(loss.normalization) ...
            && isfield(loss.normalization,'scale') ...
            && ~isempty(loss.normalization.scale)
        method = "manual";
    end
    if ~isscalar(method)
        error('camels:NormalizationMethod', ...
            ['loss.normalization.method ' ...
            'must be scalar text.']);
    end
    requestedMethod = method;
    if isfield(loss,'normalization') ...
            && isfield(loss.normalization,'requested_method') ...
            && ~isempty(loss.normalization.requested_method)
        requestedMethod = lower(string( ...
            loss.normalization.requested_method));
    end
    if method == "auto"
        method = local_auto_normalization_method(loss.fnc);
    end

    floorValue = 1e-12;
    if isfield(loss,'normalization') ...
            && isfield(loss.normalization,'floor') ...
            && ~isempty(loss.normalization.floor)
        floorValue = double(loss.normalization.floor);
    end
    if ~isscalar(floorValue) ...
            || ~isfinite(floorValue) ...
            || floorValue <= 0
        error('camels:NormalizationFloor', ...
            ['loss.normalization.floor ' ...
            'must be finite and positive.']);
    end

    counts = local_observation_basin_counts(dat,names,bas.K_t);
    if any(counts <= 0)
        missing = names(counts <= 0);
        error('camels:NoObservationBasins', ...
            ['No training basins contain ' ...
            'observations for: %s.'], ...
            strjoin(cellstr(missing),', '));
    end

    hasScale = isfield(loss,'normalization') ...
        && isstruct(loss.normalization) ...
        && isfield(loss.normalization,'scale') ...
        && ~isempty(loss.normalization.scale);

    switch method
        case {"manual","fixed"}
            if ~hasScale
                error('camels:MissingManualScale', ...
                    ['Manual normalization requires ' ...
                     'loss.normalization.scale.']);
            end
            scales = local_loss_vector( ...
                loss.normalization,'scale', ...
                names,ones(numel(names),1));

        case "none"
            scales = ones(numel(names),1);

        case "initial_loss"
            if hasScale
                scales = local_loss_vector( ...
                    loss.normalization,'scale',names, ...
                    ones(numel(names),1));
            else
                if iteration ~= 1
                    error('camels:MissingInitialScale', ...
                        ['Initial-loss scales are ' ...
                        'absent after iteration 1. ' ...
                         'Restore the saved ' ...
                         'loss.normalization structure.']);
                end
                scales = local_initial_loss_scales( ...
                    loss,nTheta,mdl,dat,bas.K_t,ode, ...
                    crrBackend,calc,names);
            end

        case "benchmark"
            if double(loss.fnc) == 4 ...
                    || double(loss.fnc) == 7
                error('camels:BenchmarkNormalization', ...
                    ['Constant-observation ' ...
                    'benchmark normalization is ' ...
                     'undefined for KGE and JKGE. ' ...
                     'Use initial_loss or ' ...
                     'manual normalization.']);
            end
            scales = local_benchmark_loss_scales( ...
                loss,dat,bas.K_t,names);

        otherwise
            error('camels:NormalizationMethod', ...
                ['Unknown normalization method ' ...
                '"%s". Use auto, initial_loss, ' ...
                 'benchmark, manual, or none.'],method);
    end

    scales = max(abs(scales),floorValue);
    loss.normalization.method = char(method);
    loss.normalization.requested_method = char(requestedMethod);
    loss.normalization.scale = local_named_struct(names,scales);
    loss.normalization.floor = floorValue;
    loss.normalization.n_basin = local_named_struct(names,counts);
    loss.normalization.n_total = bas.K_t;
    loss.normalization.observed = cellstr(names);
    if method == "initial_loss" ...
            && ~hasScale
        loss.normalization.iteration = iteration;
    end
end

function method = local_auto_normalization_method(lossFnc)
%LOCAL_AUTO_NORMALIZATION_METHOD Select normalization for a loss function.

    switch double(lossFnc)
        case {1,2,5,6}
            method = "benchmark";
        case {3,4,7}
            method = "none";
        otherwise
            error('camels:NormalizationLoss', ...
                'Cannot select normalization for loss function %g.', ...
                double(lossFnc));
    end
end

function counts = local_observation_basin_counts(dat,names,Kt)
%LOCAL_OBSERVATION_BASIN_COUNTS Count eligible training basins by type.

    counts = zeros(numel(names),1);
    for j = 1:numel(names)
        field = char(names(j));
        for k = 1:Kt
            if isfield(dat{k},'stats') ...
                    && isfield(dat{k}.stats,field) ...
                    && isfield(dat{k}.stats.(field),'train') ...
                    && dat{k}.stats.(field).train.has_data
                counts(j) = counts(j) + 1;
            end
        end
    end
end

function scales = local_initial_loss_scales( ...
    loss,nTheta,mdl,dat,Kt,ode,crrBackend,calc,names)
%LOCAL_INITIAL_LOSS_SCALES Evaluate the unweighted initial objectives.

    probeLoss = loss;
    probeLoss.normalization = struct( ...
        'method','none', ...
        'scale',local_named_struct(names,ones(numel(names),1)));
    request = crr_request(struct('q',false,'gradient',true, ...
        'metrics',false,'attribution',false, ...
        'obs',names,'jac',names));
    raw = nan(numel(names),Kt);

    if ismember(calc,{'par','parfeval'})
        parfor k = 1:Kt
            raw(:,k) = local_probe_observation_losses( ...
                crrBackend,nTheta(:,k),mdl,dat{k}, ...
                ode,probeLoss,request,names);
        end
    else
        for k = 1:Kt
            raw(:,k) = local_probe_observation_losses( ...
                crrBackend,nTheta(:,k),mdl,dat{k}, ...
                ode,probeLoss,request,names);
        end
    end

    scales = mean(raw,2,'omitnan');
    if any(~isfinite(scales))
        bad = names(~isfinite(scales));
        error('camels:InitialLossScale', ...
            'Cannot estimate initial loss scale for: %s.', ...
            strjoin(cellstr(bad),', '));
    end
end

function values = local_probe_observation_losses( ...
    crrBackend,x,mdl,dat,ode,loss,request,names)
%LOCAL_PROBE_OBSERVATION_LOSSES Return raw named losses for one basin.

    [~,out] = run_crr_local( ...
        crrBackend,x,mdl,dat,ode,loss,request);
    values = nan(numel(names),1);
    if ~isfield(out,'loss')
        return
    end
    for j = 1:numel(names)
        field = char(names(j));
        if isfield(out.loss,field)
            values(j) = out.loss.(field);
        end
    end
end

function scales = local_benchmark_loss_scales(loss,dat,Kt,names)
%LOCAL_BENCHMARK_LOSS_SCALES Use a constant-mean observation benchmark.

    raw = nan(numel(names),Kt);
    for j = 1:numel(names)
        field = char(names(j));
        for k = 1:Kt
            if ~isfield(dat{k},'stats') ...
                    || ~isfield(dat{k}.stats,field) ...
                    || ~dat{k}.stats.(field).train.has_data
                continue
            end
            y = dat{k}.obs.(field).value(:);
            mu = dat{k}.stats.(field).train.mean;
            sim = repmat(mu,numel(y),1);
            J = zeros(numel(y),1);
            raw(j,k) = observation_loss( ...
                names(j),sim,J,dat{k},loss);
        end
    end
    scales = mean(raw,2,'omitnan');
    if any(~isfinite(scales) | scales <= 0)
        bad = names(~isfinite(scales) | scales <= 0);
        error('camels:BenchmarkLossScale', ...
            'Cannot estimate benchmark loss scale for: %s.', ...
            strjoin(cellstr(bad),', '));
    end
end
