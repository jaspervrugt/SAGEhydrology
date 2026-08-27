function [Q,curr,nTheta_out] = pproc_SAGE(part,nTheta,mdl, ...
    dat,bas,prd,ode,loss,dirres,met,prf_check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PPROC_SAGE Postprocess SAGE or SITE simulations for explicit
% basin/split scenario groups and plot-ready discharge outputs.
%
% SYNOPSIS:
%   Q = pproc_SAGE( ...
%       part,nTheta,mdl,dat,bas,prd,ode,loss,dirres)
%   [Q,curr,nTheta_out] = pproc_SAGE( ...
%       part,nTheta,mdl,dat,bas,prd,ode,loss,dirres,met)
%   Q = pproc_SAGE( ...
%       part,nTheta,mdl,dat,bas,prd,ode,loss,dirres,met,prf)
%
%   part        which results are analyzed?
%     'site'     SITE single-basin calibration results
%     'sage'     SAGE regional/continental training results
%
%   nTheta      normalized hydrologic parameter values
%     'sage'     d x K x i_max or d x K matrix of submitted parameters
%     'site'     ignored; parameters are read from SITE result files
%
%   mdl         structure with model, solver, and split information
%    .model      model index
%                 1 hymod
%                 2 hmodel
%                 3 sacsma
%                 4 xinanjiang
%                 5 gr4jA
%                 6 hbv
%                 7 gr4jB [analytic routing]
%    .mcode      numerical model implementation
%                 1 Runge--Kutta MATLAB
%                 2 ode45 MATLAB
%                 3 explicit Euler MATLAB
%                 4 Runge--Kutta C++ MEX
%    .crr_backend backend selected by prepare_crr_backend
%                 'cpp' or 'matlab'
%    .calc       model execution
%      'seq'       sequential execution of watersheds
%      'par'       asynchronous parallel execution using rolling parfeval
%      'parfeval'  asynchronous parallel execution using rolling parfeval
%    .mode       assessment design
%                 1 = training basins | training split
%                 2 = training basins | training + evaluation splits
%                 3 = training + evaluation basins | training split
%                 4 = training + evaluation basins | training + eval splits
%    .sp_method  temporal split method
%                 'manual'
%                 'traditional_block'
%                 'block'
%                 'deterministic_block'
%                 'random_block'
%                 'random'
%                 'deterministic_kfold'
%                 'random_kfold'
%                 'rainfall_block'
%    .local      optional rainfall-split flag
%                 0 = global split (default)
%                     mdl.id_train and mdl.id_eval define common
%                     train/eval indices for all basins
%                 1 = rainfall-based local split
%                     each basin k uses its own rainfall-derived
%                     train/eval indices stored in
%                     dat{k}.id_train and dat{k}.id_eval
%
%                 For all non-rainfall split methods, mdl.local should
%                 remain 0.
%    .names      list of model names
%    .y0         m x 1 vector of initial states
%    .pspace     parameter-space flag
%                 0 physical parameter space
%                 1 unit-cube normalized parameter space
%                 2 unconstrained transformed parameter space
%    .th_min     d x 1 vector of lower parameter bounds
%    .th_max     d x 1 vector of upper parameter bounds
%    .par_names  1 x d cell array with parameter names
%    .tout       model output times
%    .idx        scored-window indices [compact pair or full vector]
%    .id_train   global training indices [compact pair or full vector]
%    .id_eval    global evaluation indices [compact pair or full vector]
%
%   dat         K x 1 cell array with basin data
%    {k}.gauge      gauge identifier
%    {k}.y_n        observed discharge series [normalized]
%    {k}.bad        logical vector with invalid observations
%    {k}.id_train   optional basin-specific training indices
%    {k}.id_eval    optional basin-specific evaluation indices
%    {k}.stats      split statistics used by NSE/KGE diagnostics
%      .TSSt        total sum of squares, training split
%      .TSSe        total sum of squares, evaluation split
%      .mut         mean discharge, training split
%      .stdt        standard deviation, training split
%      .mue         mean discharge, evaluation split
%      .stde        standard deviation, evaluation split
%    {k}.jkge.m_y   observed JKGE benchmark on the full model time axis
%
%   bas         basin metadata
%    .K          total number of basins
%    .K_t        number of training basins
%    .K_e        number of evaluation basins
%    .r          number of basin attributes
%    .id_gauge   gauge IDs in final basin order [train; eval]
%    .gname      gauge names in final basin order [train; eval]
%
%   prd         temporal settings
%    .dt         temporal resolution [1 daily, 24 hourly, 96 at 15-minute]
%    .spinup     spin-up period in days
%    .[others]   split-specific settings used upstream
%
%   ode         ODE solver settings
%
%   loss        loss-function settings
%    .fnc       scalar loss-function selector
%                   1 sum of absolute residuals
%                   2 generalized least squares
%                   3 Nash--Sutcliffe efficiency
%                   4 Kling--Gupta efficiency
%                   5 Huber loss
%                   6 flow-duration-curve loss
%                   7 Jawad--Kling--Gupta efficiency
%    .n_win     JKGE benchmark window length for methods 1 and 2
%    .method    JKGE benchmark method
%                   1 moving-average mean
%                   2 section-wise mean
%                   3 long-term mean
%                   4 monthly climatology
%    .M        JKGE benchmark-matching definition
%                   1 paper ratio-based M component
%                   2 revised norm-based M component
%                 If missing or empty, M = 2 is used.
%    .meta       metadata
%
%   dirres      results directory of SAGE/SITE
%
%   met         OPTIONAL metrics structure from camels.m
%               For SAGE, supplied values are used directly. Public fields
%               are
%               .loss         objective-related diagnostics
%               .performance  NSE/KGE/JKGE for train/eval periods
%               .hydro        reserved for final simulated hydrologic
%                             metrics/signatures for train/eval periods
%
%   prf_check   OPTIONAL performance structure. When supplied, metrics are
%               recomputed from Q and compared with both met and
%               prf_check.curr.
%
%   Q           simulated discharge structure
%    .tt         training basins | training split
%    .te         training basins | evaluation split
%    .et         evaluation basins | training split
%    .ee         evaluation basins | evaluation split
%
%   curr        optional grouped performance structure with fields
%               .NSE, .KGE, .S_fdc and .JKGE
%   nTheta_out  normalized hydrologic parameter values
%     'site'     d x K x n_m array, one layer per model
%     'sage'     d x K x 1 array for the submitted SAGE model
%
% NOTES:
%   1. Scenario notation is explicit throughout:
%        tt = training basins   | training split
%        te = training basins   | evaluation split
%        et = evaluation basins | training split
%        ee = evaluation basins | evaluation split
%
%   2. Global split mode:
%        mdl.id_train and mdl.id_eval define common train/eval indices
%        for all basins.
%
%   3. Local split mode:
%        if mdl.local = 1, each basin may define its own train/eval
%        indices through dat{k}.id_train and dat{k}.id_eval.
%
%   4. Block-type split methods store Q.tt/Q.te/Q.et/Q.ee on their
%        natural train/eval windows.
%
%   5. Mask-type and local split methods store Q fields on the full scored
%        window so that mask-aware and basin-specific plotting remain
%        possible.
%
%   6. JKGE is computed on the full model time axis with an explicit valid
%        index vector. Missing observations are excluded through the index,
%        but the simulated benchmark m_q is formed before reduction.
%
%   7. If met is supplied for SAGE, current metrics are taken directly
%        from met. Passing prf_check explicitly requests recomputation and
%        comparison against prf_check.curr.
%
%   8. Future simulation-derived hydrologic metrics belong in met.hydro.
%      Compute them here from the final Q simulation, once after training
%      has finished. Do not compute them in camels.m or crr_model.m during
%      every optimization iteration unless live GUI monitoring is enabled
%      deliberately and its additional runtime is accepted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Feb. 2026 / updated Aug. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 9
        dirres = '';
    end
    if nargin < 10
        met = [];
    end
    if nargin < 11
        prf_check = [];
    end
    verify_met = ~isempty(prf_check);
    if verify_met && (~isstruct(prf_check) ...
            || ~isfield(prf_check,'curr'))
        error('pproc_SAGE:InvalidPerformanceCheck', ...
            'The optional final argument must contain prf.curr.');
    end

    t_pproc_total = tic;
    
    if ~exist('loss','var') || isempty(loss) ...
            || ~isstruct(loss) || ~isfield(loss,'fnc') ...
            || isempty(loss.fnc) || ~isscalar(loss.fnc) ...
            || ~isnumeric(loss.fnc) ...
            || ~ismember(double(loss.fnc),1:7)
        error(['pproc_SAGE: loss.fnc must ' ...
            'be a numeric scalar ' ...
            'with value 1,2,3,4,5,6, or 7.']);
    end
    
    loss_fnc = double(loss.fnc);
    
    if ~exist('expand_index','file')
        error(['pproc_SAGE: ' ...
            'Required function expand_index ' ...
            'is not on the MATLAB path.']);
    end
    if ~exist('crr_model','file')
        error(['pproc_SAGE: ' ...
            'Required function crr_model ' ...
            'is not on the MATLAB path.']);
    end
    if ~isfield(mdl,'crr_backend') ...
            || isempty(mdl.crr_backend)
        error('pproc_SAGE:MissingCRRBackend', ...
            ['mdl.crr_backend is missing. Call prepare_crr_backend before ' ...
             'pproc_SAGE so postprocessing uses the prepared backend.']);
    end
    crr_backend = char(lower(string(mdl.crr_backend)));
    if ~ismember(crr_backend,{'cpp','matlab'})
        error('pproc_SAGE:UnknownCRRBackend', ...
            'Unknown mdl.crr_backend: %s.',crr_backend);
    end
    if ~exist('read_model','file')
        error(['pproc_SAGE: ' ...
            'Required function read_model ' ...
            'is not on the MATLAB path.']);
    end
    
    part = char(lower(string(part)));
    
    if ~ismember(part,{'site','sage'})
        error(['pproc_SAGE: ' ...
            'Unknown part = %s.'],part);
    end
    
    if nargin < 9 ...
            || isempty(dirres)
        error(['pproc_SAGE: ' ...
            'dirres must be ' ...
            'provided explicitly.']);
    end
    
    if ~isfolder(dirres)
        error(['pproc_SAGE: ' ...
            'Results directory does ' ...
            'not exist: %s'],dirres);
    end
    
    if ~isfield(mdl,'names') ...
            || isempty(mdl.names)
        error(['pproc_SAGE: ' ...
            'mdl.names is missing ' ...
            'or empty.']);
    end
    model_names = mdl.names;
    
    if ~isfield(mdl,'model') ...
            || isempty(mdl.model)
        if strcmp(part,'sage')
            error(['pproc_SAGE: ' ...
                'mdl.model is required for ' ...
                'SAGE postprocessing.']);
        else
            model = [];
        end
    else
        model = mdl.model;
    end
    
    if ~isfield(mdl,'calc') ...
            || isempty(mdl.calc)
        error(['pproc_SAGE: ' ...
            'mdl.calc must be ''seq'' ' ...
            'or ''par''.']);
    end
    calc = char(lower(string(mdl.calc)));
    if ~ismember(calc,{'seq','par','parfeval'})
        error(['pproc_SAGE: ' ...
            'Unknown mdl.calc = %s.'],calc);
    end
    
    if isfield(mdl,'sp_method') ...
            && ~isempty(mdl.sp_method)
        sp_method = char(lower( ...
            string(mdl.sp_method)));
    elseif isfield(mdl,'split_method') ...
            && ~isempty(mdl.split_method)
        sp_method = char(lower( ...
            string(mdl.split_method)));
    else
        sp_method = 'manual';
    end
    
    is_local_split = isfield(mdl,'local') ...
        && isequal(double(mdl.local),1);
    
    is_block_method = ~is_local_split ...
        && any(strcmpi(sp_method,{ ...
        'traditional_block', ...
        'manual', ...
        'block', ...
        'deterministic_block'}));
    
    if is_local_split
        id_q_tmp = expand_index(mdl.idx);
        nq_tmp = numel(id_q_tmp) - 1;
        id_tr = [];
        id_ev = [];
        id_all = 1:nq_tmp;
    else
        id_tr = expand_index(mdl.id_train);
        if isfield(mdl,'id_eval') ...
                && ~isempty(mdl.id_eval)
            id_ev = expand_index(mdl.id_eval);
        else
            id_ev = [];
        end
        id_all = unique([id_tr(:); id_ev(:)]).';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine temporal resolution --> daily | hourly for file naming    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt = [];
    if isfield(prd,'dt') ...
            && ~isempty(prd.dt)
        dt = prd.dt;
    end
    if isempty(dt)
        dt = 1;
    end
    
    if dt == 1
        dt_tag = 'daily';
    elseif dt == 24
        dt_tag = 'hourly';
    elseif dt == 96
        dt_tag = '15min';
    else
        fprintf('\n');
        error(['pproc_SAGE: ' ...
            'Unknown dt = %g. Use dt = 1 ' ...
            '(daily), 24 (hourly), or 96 (15-minute).'],dt);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Verify right order of data and basin IDs                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = bas.K;
    K_t = bas.K_t;
    K_e = bas.K_e;
    
    id_gauge = string(bas.id_gauge(:));
    
    id_list = strings(K,1);
    for k = 1:K
        if isfield(dat{k},'gauge') ...
                && ~isempty(dat{k}.gauge)
            id_list(k) = string(dat{k}.gauge);
        end
    end
    
    id_gauge = local_normalize_runtime_gauge_id(id_gauge);
    id_list = local_normalize_runtime_gauge_id(id_list);
    
    if numel(id_gauge) ~= numel(id_list)
        fprintf(['      pproc_SAGE: Mismatch of ' ...
            'gauge basin id vector lengths\n']);
    elseif all(id_gauge(:) == id_list(:))
        fprintf(['      pproc_SAGE: Perfect match of ' ...
            'gauge basin id codes\n']);
    else
        fprintf(['      pproc_SAGE: Mismatch of ' ...
            'gauge basin id codes\n']);
    end
    
    switch part
        case 'site'
            n_m = numel(model_names);
            id_model = 1:n_m;
        case 'sage'
            n_m = 1;
            id_model = model(1);
        otherwise
            error(['pproc_SAGE: ' ...
                'Unknown part = %s.'],part);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Allocate outputs                                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = numel(id_tr);
    m = numel(id_ev);
    nall = numel(id_all);
    
    % Discharge histories dominate RAM for hourly/15-minute runs.
    % Single precision halves the footprint and is sufficient for retained
    % postprocessing discharge and plotting.
    if is_block_method
        Qtt = nan(n,K_t,n_m,'single');
        Qte = nan(m,K_t,n_m,'single');
        Qet = nan(n,K_e,n_m,'single');
        Qee = nan(m,K_e,n_m,'single');
    else
        Qtt = nan(nall,K_t,n_m,'single');
        Qte = nan(nall,K_t,n_m,'single');
        Qet = nan(nall,K_e,n_m,'single');
        Qee = nan(nall,K_e,n_m,'single');
    end
    
    NSEtt = nan(K_t,n_m);
    NSEte = nan(K_t,n_m);
    NSEet = nan(K_e,n_m);
    NSEee = nan(K_e,n_m);
    
    KGEtt = nan(K_t,n_m);
    KGEte = nan(K_t,n_m);
    KGEet = nan(K_e,n_m);
    KGEee = nan(K_e,n_m);

    KGE_rtt = nan(K_t,n_m); KGE_rte = nan(K_t,n_m);
    KGE_ret = nan(K_e,n_m); KGE_ree = nan(K_e,n_m);
    KGE_alphatt = nan(K_t,n_m); KGE_alphate = nan(K_t,n_m);
    KGE_alphaet = nan(K_e,n_m); KGE_alphaee = nan(K_e,n_m);
    KGE_betatt = nan(K_t,n_m); KGE_betate = nan(K_t,n_m);
    KGE_betaet = nan(K_e,n_m); KGE_betaee = nan(K_e,n_m);
    
    JKGEtt = nan(K_t,n_m);
    JKGEte = nan(K_t,n_m);
    JKGEet = nan(K_e,n_m);
    JKGEee = nan(K_e,n_m);

    Dfdctt = nan(K_t,n_m);
    Dfdcte = nan(K_t,n_m);
    Dfdcet = nan(K_e,n_m);
    Dfdcee = nan(K_e,n_m);

    JKGE_Mtt = nan(K_t,n_m); JKGE_Mte = nan(K_t,n_m);
    JKGE_Met = nan(K_e,n_m); JKGE_Mee = nan(K_e,n_m);
    JKGE_Vtt = nan(K_t,n_m); JKGE_Vte = nan(K_t,n_m);
    JKGE_Vet = nan(K_e,n_m); JKGE_Vee = nan(K_e,n_m);
    JKGE_Ctt = nan(K_t,n_m); JKGE_Cte = nan(K_t,n_m);
    JKGE_Cet = nan(K_e,n_m); JKGE_Cee = nan(K_e,n_m);
    
    err_tt = nan(K_t,n_m);
    err_te = nan(K_t,n_m);
    err_et = nan(K_e,n_m);
    err_ee = nan(K_e,n_m);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read parameter values                                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch part
    
        case 'site'
            P = cell(1,n_m);      % physical parameters
            Pn = cell(1,n_m);     % normalized parameters
        
            for im = 1:n_m
                mname = local_model_name(model_names,im);
        
                [tf_file,file_used,file_kind] = ...
                    site_param_file_exists_local( ...
                    dirres,mname,dt_tag);
        
                if ~tf_file
                    fprintf(['      pproc_SAGE: ' ...
                        'no SITE parameter file for ' ...
                        'model %s --> outputs kept ' ...
                        'as NaN\n'],mname);
                    P{im} = [];
                    Pn{im} = [];
                    continue
                end
        
                [id_param,range_id,Pn_im] = ...
                    read_site_param_normalized_local( ...
                    dirres,mname,dt_tag,loss_fnc,part);
        
                [tf,loc] = ismember( ...
                    normalize_usgs(id_gauge), ...
                    normalize_usgs(id_param));
        
                if ~all(tf)
                    miss = string(id_gauge(~tf));
                    error(['pproc_SAGE: ' ...
                        'Parameter file exists ' ...
                        'for model %s (%s, %s), but ' ...
                        'is missing basin values. ' ...
                        'Example missing gauge id: %s.'], ...
                        mname,file_used,file_kind, ...
                        char(miss(1)));
                end
        
                range_id = range_id(loc);
                Pn_im = Pn_im(:,loc);
        
                [th_min_mat,th_max_mat] = ...
                    read_site_param_ranges_local( ...
                    dirres,mname,dt_tag, ...
                    range_id,size(Pn_im,1));
        
                % store normalized parameters
                Pn{im} = Pn_im;
        
                % convert to physical space for model simulation
                P{im} = th_min_mat + Pn_im .* ...
                    (th_max_mat - th_min_mat);
            end
    
        case 'sage'
            P = cell(1,1);    % physical/working parameter matrix used below
            Pn = cell(1,1);   % normalized parameter matrix to return
        
            if ndims(nTheta) == 3
                Pn{1} = nTheta(:,:,end);
            else
                Pn{1} = nTheta;
            end
            % In SAGE postprocessing, submitted nTheta is already normalized
            % parameter matrix expected by crr_model when mdl.pspace = 1.
            P{1} = Pn{1};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Package normalized parameter output                                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch part
        case 'site'
            d_out = NaN;
            for im = 1:n_m
                if ~isempty(Pn{im})
                    d_out = size(Pn{im},1);
                    break
                end
            end
    
            if isfinite(d_out)
                nTheta_out = nan(d_out,K,n_m);
                for im = 1:n_m
                    if ~isempty(Pn{im})
                        nTheta_out(:,:,im) = Pn{im};
                    end
                end
            else
                nTheta_out = [];
            end
    
        case 'sage'
            if ~isempty(Pn{1})
                d_out = size(Pn{1},1);
                nTheta_out = nan(d_out,K,1);
                nTheta_out(:,:,1) = Pn{1};
            else
                nTheta_out = [];
            end
    
        otherwise
            nTheta_out = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parallel fallback if no pool exists                                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch calc
        case 'seq'
            % ok
        case {'par','parfeval'}
            has_gcp = exist('gcp','file') == 2;
            if ~has_gcp
                fprintf(['      Warning: pproc_SAGE: ' ...
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
                    fprintf(['      Warning: pproc_SAGE: ' ...
                        'no parallel pool found. ' ...
                        'Falling back to sequential.\n']);
                    calc = 'seq';
                end
            end
        otherwise
            error(['pproc_SAGE: ' ...
                'unknown mdl.calc = ''%s''. ' ...
                'Use ''seq'', ''par'', ' ...
                'or ''parfeval''.'],calc);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main loop                                                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for im = 1:n_m
        model = id_model(im);
        mdl.model = model;
        mdl = read_model(mdl,prd,false);
        % Postprocessing needs discharge only. Stream q from the model
        % instead of retaining the complete augmented state trajectory.
        ode_pp = ode;
        if mdl.mcode ~= 2
            ode_pp.mem = 0;
        end

        if strcmp(part,'site') ...
                && isempty(P{im})
            fprintf(['      pproc_SAGE: model ' ...
                '%s skipped in postprocessing; ' ...
                'outputs left as NaN.\n'], ...
                local_model_name(model_names,im));
            continue
        end
    
        if strcmp(part,'site')
            mdl.pspace = 0;
        else
            mdl.pspace = 1;
        end
    
        id_q = expand_index(mdl.idx);
        nq = numel(id_q) - 1;
    
        if nq < 1
            error(['pproc_SAGE: ' ...
                'Expanded mdl.idx is empty.']);
        end

        Theta = P{im};
        if isempty(Theta)
            continue
        end

        % Decide whether parallel worker serialization is safe for this
        % model and requested time window. The estimate uses currently
        % available client memory and the expected temporary footprint.
        calc_model = calc;
        if ~strcmp(calc_model,'seq')
            [useParallel,memInfo] = local_parallel_memory_ok( ...
                nq,K,K_t,K_e,n,m,nall,n_m, ...
                is_block_method,dat,Theta);

            fprintf(['      pproc_SAGE memory check: available %.2f GB, ' ...
                'estimated additional peak %.2f GB, safety reserve ' ...
                '%.2f GB, workers %d.\n'], ...
                memInfo.availableGB,memInfo.estimatedPeakGB, ...
                memInfo.requiredHeadroomGB,memInfo.nWorkers);

            if ~useParallel
                fprintf(['      Warning: pproc_SAGE: insufficient memory ' ...
                    'for safe worker deserialization; using sequential ' ...
                    'postprocessing for this model.\n']);
                calc_model = 'seq';
            end
        end

        % Largest temporary array in postprocessing.
        Qfull = nan(nq,K,'single');

        if size(Theta,2) ~= K
            error(['pproc_SAGE: ' ...
                'Parameter matrix has %d columns, ' ...
                'but K = %d basins were expected.'], ...
                size(Theta,2),K);
        end

        req_crr = crr_request(struct('q',true));

        t_sim = tic;
        switch calc_model
            case 'seq'
                for k = 1:K
                    if ~local_complete_forcing(dat{k})
                        continue
                    end
                    [~,outk] = local_run_crr_backend( ...
                        crr_backend,Theta(:,k),mdl,dat{k}, ...
                        ode_pp,loss,req_crr);
                    qk = outk.q;
                    if numel(qk) ~= nq
                        error('pproc_SAGE:DischargeLengthMismatch', ...
                            ['Basin %s returned %d discharge values, ' ...
                             'but %d were expected from mdl.idx.'], ...
                            char(id_gauge(k)),numel(qk),nq);
                    end
                    Qfull(:,k) = single(qk(:));
                    clear qk
                end

            case {'par','parfeval'}
                % Genuine asynchronous execution. Keep one basin future
                % per worker in flight and launch the next basin as soon as
                % any worker finishes. This avoids the synchronization
                % barriers of fixed batches and limits queued input copies.
                Qfull = local_run_parfeval_basins( ...
                    Qfull,Theta,mdl,dat,ode_pp,loss, ...
                    nq,id_gauge,crr_backend);
        end
        if any(strcmp(calc_model,{'par','parfeval'}))
            execution_label = 'parfeval';
        else
            execution_label = calc_model;
        end
        fprintf(['      pproc_SAGE: basin simulations ' ...
            'finished in %.2f s (%s).\n'], ...
            toc(t_sim),execution_label);
    
        % ------------------------------------------------------------
        % Future hydrologic signatures for simulated discharge should be
        % computed once from the final Qfull in this postprocessing path
        % and stored under met.hydro (train/eval fields). Keep that work
        % outside the per-iteration camels/crr_model optimization path.
        %
        % Recompute metrics from Qfull only if needed for verification
        % or if met was not supplied.
        % ------------------------------------------------------------
        need_recompute = isempty(met) ...
            || verify_met;
    
        if need_recompute
            t_met = tic;
            [NSEtt(:,im),KGEtt(:,im),JKGEtt(:,im), ...
                KGE_rtt(:,im),KGE_alphatt(:,im), ...
                KGE_betatt(:,im),JKGE_Mtt(:,im), ...
                JKGE_Vtt(:,im),JKGE_Ctt(:,im),Dfdctt(:,im), ...
                err_tt(:,im)] = local_metric_block( ...
                Qfull,dat,loss,id_tr,1:K_t,'t', ...
                is_local_split,id_tr);
            
            [NSEte(:,im),KGEte(:,im),JKGEte(:,im), ...
                KGE_rte(:,im),KGE_alphate(:,im), ...
                KGE_betate(:,im),JKGE_Mte(:,im), ...
                JKGE_Vte(:,im),JKGE_Cte(:,im),Dfdcte(:,im), ...
                err_te(:,im)] = local_metric_block( ...
                Qfull,dat,loss,id_ev,1:K_t,'e', ...
                is_local_split,id_tr);
    
            if K_e > 0
                [NSEet(:,im),KGEet(:,im),JKGEet(:,im), ...
                    KGE_ret(:,im),KGE_alphaet(:,im), ...
                    KGE_betaet(:,im),JKGE_Met(:,im), ...
                    JKGE_Vet(:,im),JKGE_Cet(:,im),Dfdcet(:,im), ...
                    err_et(:,im)] = ...
                    local_metric_block(Qfull, ...
                    dat,loss,id_tr,K_t+1:K,'t', ...
                    is_local_split,id_tr);
            
                [NSEee(:,im),KGEee(:,im),JKGEee(:,im), ...
                    KGE_ree(:,im),KGE_alphaee(:,im), ...
                    KGE_betaee(:,im),JKGE_Mee(:,im), ...
                    JKGE_Vee(:,im),JKGE_Cee(:,im),Dfdcee(:,im), ...
                    err_ee(:,im)] = ...
                    local_metric_block(Qfull, ...
                    dat,loss,id_ev,K_t+1:K,'e', ...
                    is_local_split,id_tr);
            end

            if verify_met ...
                    && ~isempty(met) ...
                    && strcmpi(part,'sage')
                NSEtt_chk = NSEtt(:,im); KGEtt_chk = KGEtt(:,im);
                JKGEtt_chk = JKGEtt(:,im);
                NSEte_chk = NSEte(:,im); KGEte_chk = KGEte(:,im);
                JKGEte_chk = JKGEte(:,im);
                NSEet_chk = NSEet(:,im); KGEet_chk = KGEet(:,im);
                JKGEet_chk = JKGEet(:,im);
                NSEee_chk = NSEee(:,im); KGEee_chk = KGEee(:,im);
                JKGEee_chk = JKGEee(:,im);
            end
            fprintf(['      pproc_SAGE: metric ' ...
                'processing finished in %.2f s.\n'], ...
                toc(t_met));
        end
        
        % -----------------------------------------------
        % If met is supplied --> use it directly for SAGE
        % -----------------------------------------------
        if ~isempty(met) ...
                && strcmpi(part,'sage')
            NSEtt(:,im) = met.performance.t.NSE(1:K_t).';
            KGEtt(:,im) = met.performance.t.KGE(1:K_t).';
            JKGEtt(:,im) = met.performance.t.JKGE(1:K_t).';
            KGE_rtt(:,im) = met.performance.t.KGE_components.r(1:K_t).';
            KGE_alphatt(:,im) = ...
                met.performance.t.KGE_components.alpha(1:K_t).';
            KGE_betatt(:,im) = ...
                met.performance.t.KGE_components.beta(1:K_t).';
            JKGE_Mtt(:,im) = met.performance.t.JKGE_components.M(1:K_t).';
            JKGE_Vtt(:,im) = met.performance.t.JKGE_components.V(1:K_t).';
            JKGE_Ctt(:,im) = met.performance.t.JKGE_components.C(1:K_t).';
            Dfdctt(:,im) = met.performance.t.D_fdc(1:K_t).';
    
            if K_t > 0
                NSEte(:,im) = met.performance.e.NSE(1:K_t).';
                KGEte(:,im) = met.performance.e.KGE(1:K_t).';
                JKGEte(:,im) = met.performance.e.JKGE(1:K_t).';
                KGE_rte(:,im) = met.performance.e.KGE_components.r(1:K_t).';
                KGE_alphate(:,im) = ...
                    met.performance.e.KGE_components.alpha(1:K_t).';
                KGE_betate(:,im) = ...
                    met.performance.e.KGE_components.beta(1:K_t).';
                JKGE_Mte(:,im) = met.performance.e.JKGE_components.M(1:K_t).';
                JKGE_Vte(:,im) = met.performance.e.JKGE_components.V(1:K_t).';
                JKGE_Cte(:,im) = met.performance.e.JKGE_components.C(1:K_t).';
                Dfdcte(:,im) = met.performance.e.D_fdc(1:K_t).';
            end
    
            if K_e > 0
                NSEet(:,im) = met.performance.t.NSE(K_t+1:K).';
                KGEet(:,im) = met.performance.t.KGE(K_t+1:K).';
                JKGEet(:,im) = met.performance.t.JKGE(K_t+1:K).';
                KGE_ret(:,im) = met.performance.t.KGE_components.r(K_t+1:K).';
                KGE_alphaet(:,im) = ...
                    met.performance.t.KGE_components.alpha(K_t+1:K).';
                KGE_betaet(:,im) = ...
                    met.performance.t.KGE_components.beta(K_t+1:K).';
                JKGE_Met(:,im) = met.performance.t.JKGE_components.M(K_t+1:K).';
                JKGE_Vet(:,im) = met.performance.t.JKGE_components.V(K_t+1:K).';
                JKGE_Cet(:,im) = met.performance.t.JKGE_components.C(K_t+1:K).';
                Dfdcet(:,im) = met.performance.t.D_fdc(K_t+1:K).';
    
                NSEee(:,im) = met.performance.e.NSE(K_t+1:K).';
                KGEee(:,im) = met.performance.e.KGE(K_t+1:K).';
                JKGEee(:,im) = met.performance.e.JKGE(K_t+1:K).';
                KGE_ree(:,im) = met.performance.e.KGE_components.r(K_t+1:K).';
                KGE_alphaee(:,im) = ...
                    met.performance.e.KGE_components.alpha(K_t+1:K).';
                KGE_betaee(:,im) = ...
                    met.performance.e.KGE_components.beta(K_t+1:K).';
                JKGE_Mee(:,im) = met.performance.e.JKGE_components.M(K_t+1:K).';
                JKGE_Vee(:,im) = met.performance.e.JKGE_components.V(K_t+1:K).';
                JKGE_Cee(:,im) = met.performance.e.JKGE_components.C(K_t+1:K).';
                Dfdcee(:,im) = met.performance.e.D_fdc(K_t+1:K).';
            end
    
            if verify_met
                local_compare_metric_block('NSE tt', ...
                    NSEtt(:,im),NSEtt_chk);
                local_compare_metric_block('NSE te', ...
                    NSEte(:,im),NSEte_chk);
                local_compare_metric_block('KGE tt', ...
                    KGEtt(:,im),KGEtt_chk);
                local_compare_metric_block('KGE te', ...
                    KGEte(:,im),KGEte_chk);
                local_compare_metric_block('JKGE tt', ...
                    JKGEtt(:,im),JKGEtt_chk);
                local_compare_metric_block('JKGE te', ...
                    JKGEte(:,im),JKGEte_chk);
                if K_e > 0
                    local_compare_metric_block('NSE et', ...
                        NSEet(:,im),NSEet_chk);
                    local_compare_metric_block('NSE ee', ...
                        NSEee(:,im),NSEee_chk);
                    local_compare_metric_block('KGE et', ...
                        KGEet(:,im),KGEet_chk);
                    local_compare_metric_block('KGE ee', ...
                        KGEee(:,im),KGEee_chk);
                    local_compare_metric_block('JKGE et', ...
                        JKGEet(:,im),JKGEet_chk);
                    local_compare_metric_block('JKGE ee', ...
                        JKGEee(:,im),JKGEee_chk);
                end
            end
        end
    
        % -------------------------
        % Package discharge outputs
        % -------------------------
        if is_block_method
            if n > 0
                Qtt(:,:,im) = Qfull(id_tr,1:K_t);
            end
    
            if m > 0
                Qte(:,:,im) = Qfull(id_ev,1:K_t);
            end
    
            if K_e > 0
                if n > 0
                    Qet(:,:,im) = Qfull(id_tr,K_t+1:K);
                end
                if m > 0
                    Qee(:,:,im) = Qfull(id_ev,K_t+1:K);
                end
            end
        else
            Qtt(:,:,im) = Qfull(id_all,1:K_t);
            Qte(:,:,im) = Qfull(id_all,1:K_t);
    
            if K_e > 0
                Qet(:,:,im) = Qfull(id_all,K_t+1:K);
                Qee(:,:,im) = Qfull(id_all,K_t+1:K);
            end
        end
    
        % Release the full temporary simulation matrix once scenario
        % outputs have been retained.
        clear Qfull

        % ----------------------
        % Report TSS consistency
        % ----------------------
        if need_recompute
            report_err_block(local_tag_string( ...
                'training basins','training', ...
                sp_method),err_tt(:,im));
    
            report_err_block(local_tag_string( ...
                'training basins','evaluation', ...
                sp_method),err_te(:,im));
    
            report_err_block(local_tag_string( ...
                'evaluation basins','training', ...
                sp_method),err_et(:,im));
    
            report_err_block(local_tag_string( ...
                'evaluation basins','evaluation', ...
                sp_method),err_ee(:,im));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Package explicit outputs                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q = struct();
    Q.tt = Qtt;
    Q.te = Qte;
    Q.et = Qet;
    Q.ee = Qee;
    S_fdc = struct('tt',[],'te',[],'et',[],'ee',[]);
    if isfield(loss,'fdc') && isstruct(loss.fdc) ...
            && isfield(loss.fdc,'D0t') ...
            && isfield(loss.fdc,'D0e') ...
            && numel(loss.fdc.D0t) == K ...
            && numel(loss.fdc.D0e) == K
        D0t = double(loss.fdc.D0t(:));
        D0e = double(loss.fdc.D0e(:));
        S_fdc.tt = local_fdc_score_matrix(Dfdctt,D0t(1:K_t));
        S_fdc.te = local_fdc_score_matrix(Dfdcte,D0e(1:K_t));
        if K_e > 0
            S_fdc.et = local_fdc_score_matrix( ...
                Dfdcet,D0t(K_t+1:K));
            S_fdc.ee = local_fdc_score_matrix( ...
                Dfdcee,D0e(K_t+1:K));
        end
    end
    
    NSE = struct();
    NSE.tt = NSEtt;
    NSE.te = NSEte;
    NSE.et = NSEet;
    NSE.ee = NSEee;
    
    KGE = struct();
    KGE.tt = KGEtt;
    KGE.te = KGEte;
    KGE.et = KGEet;
    KGE.ee = KGEee;
    KGE.components.tt = struct('r',KGE_rtt, ...
        'alpha',KGE_alphatt,'beta',KGE_betatt);
    KGE.components.te = struct('r',KGE_rte, ...
        'alpha',KGE_alphate,'beta',KGE_betate);
    KGE.components.et = struct('r',KGE_ret, ...
        'alpha',KGE_alphaet,'beta',KGE_betaet);
    KGE.components.ee = struct('r',KGE_ree, ...
        'alpha',KGE_alphaee,'beta',KGE_betaee);
    
    JKGE = struct();
    JKGE.tt = JKGEtt;
    JKGE.te = JKGEte;
    JKGE.et = JKGEet;
    JKGE.ee = JKGEee;
    JKGE.components.tt = struct('M',JKGE_Mtt, ...
        'V',JKGE_Vtt,'C',JKGE_Ctt);
    JKGE.components.te = struct('M',JKGE_Mte, ...
        'V',JKGE_Vte,'C',JKGE_Cte);
    JKGE.components.et = struct('M',JKGE_Met, ...
        'V',JKGE_Vet,'C',JKGE_Cet);
    JKGE.components.ee = struct('M',JKGE_Mee, ...
        'V',JKGE_Vee,'C',JKGE_Cee);

    curr = struct('NSE',NSE,'KGE',KGE, ...
        'S_fdc',S_fdc,'JKGE',JKGE);

    if verify_met
        scenarios = {'tt','te','et','ee'};
        metricNames = {'NSE','KGE','S_fdc','JKGE'};
        for imetric = 1:numel(metricNames)
            metricName = metricNames{imetric};
            if ~isfield(prf_check.curr,metricName)
                error('pproc_SAGE:MissingPerformanceMetric', ...
                    ['prf.curr.%s is required ' ...
                    'for verification.'],metricName);
            end
            for iscenario = 1:numel(scenarios)
                scenario = scenarios{iscenario};
                local_compare_metric_block( ...
                    sprintf('%s %s',metricName,scenario), ...
                    curr.(metricName).(scenario), ...
                    prf_check.curr.(metricName).(scenario));
            end
        end
        fprintf(['      pproc_SAGE: ' ...
            'prf.curr verification passed.\n']);
    end

    fprintf(['      pproc_SAGE: total ' ...
        'postprocessing time %.2f s.\n'], ...
        toc(t_pproc_total));
end

function Qfull = local_run_parfeval_basins( ...
    Qfull,Theta,mdl,dat,ode_pp,loss,nq,id_gauge,crr_backend)
%LOCAL_RUN_PARFEVAL_BASINS Run basin simulations with batched parfeval.
%
% Match the stable execution architecture used by camels.m: submit a
% modest number of basin batches rather than creating one future per basin.
% The GUI yield callback is removed before worker submission because nested
% GUI callbacks can capture and serialize the complete GUI workspace.

    K = size(Theta,2);
    if K == 0
        return
    end

    p = gcp('nocreate');
    if isempty(p)
        error('pproc_SAGE:NoParallelPool', ...
            ['parfeval postprocessing requires an active ' ...
             'parallel pool.']);
    end

    % Use the same batching strategy as camels.m. With four workers this
    % creates sixteen futures, each processing about K/16 basins locally.
    nBatch = min(K,max(double(p.NumWorkers), ...
        4*double(p.NumWorkers)));
    edges = round(linspace(1,K+1,nBatch+1));
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
        futures(b) = parfeval(p,@local_pproc_basin_batch,2, ...
            ids,Theta(:,ids),mdl,dat(ids),ode_pp, ...
            loss,nq,id_gauge(ids),crr_backend);
    end

    nDone = 0;
    try
        for b = 1:nBatch
            [~,ids,Qb] = fetchNext(futures);

            if size(Qb,1) ~= nq ...
                    || size(Qb,2) ~= numel(ids)
                error('pproc_SAGE:DischargeBatchSizeMismatch', ...
                    ['Returned discharge batch has size %d x %d; ' ...
                     'expected %d x %d.'], ...
                    size(Qb,1),size(Qb,2),nq,numel(ids));
            end

            Qfull(:,ids) = Qb;
            nDone = nDone + numel(ids);
            clear Qb

            fprintf(['      pproc_SAGE: parfeval basins ' ...
                '%d/%d completed.\n'],nDone,K);
        end
    catch ME
        try
            cancel(futures);
        catch
        end
        rethrow(ME)
    end

    % Release future objects and any retained output metadata immediately.
    try
        delete(futures);
    catch
    end
end

function [ids,Qb] = local_pproc_basin_batch( ...
    ids,Theta_b,mdl,dat_b,ode_pp,loss,nq,gauges,crr_backend)
%LOCAL_PPROC_BASIN_BATCH Worker-side batch of discharge-only simulations.

    ids = ids(:).';
    nb = numel(ids);
    Qb = nan(nq,nb,'single');
    req_crr = crr_request(struct('q',true));

    for jj = 1:nb
        if ~local_complete_forcing(dat_b{jj})
            continue
        end
        [~,outj] = local_run_crr_backend( ...
            crr_backend,Theta_b(:,jj),mdl,dat_b{jj}, ...
            ode_pp,loss,req_crr);
        q = outj.q;

        if numel(q) ~= nq
            error('pproc_SAGE:DischargeLengthMismatch', ...
                ['Basin %s returned %d discharge values, ' ...
                 'but %d were expected from mdl.idx.'], ...
                char(string(gauges(jj))),numel(q),nq);
        end

        Qb(:,jj) = single(q(:));
        clear q
    end
end

function [tf,info] = local_parallel_memory_ok( ...
    nq,K,K_t,K_e,n,m,nall,n_m,is_block_method,dat,Theta)  %#ok
%LOCAL_PARALLEL_MEMORY_OK Estimate whether parallel postprocessing is safe.
%
% This is a conservative client-side estimate. It includes the complete
% discharge matrix, retained scenario arrays, worker return payloads,
% resident input data, and a multiplier for serialization copies and
% memory fragmentation.

    bytesSingle = 4;
    bytesDouble = 8;

    % Complete temporary discharge matrix on the MATLAB client.
    qfullBytes = double(nq) * double(K) * bytesSingle;

    % Exact retained-output allocation used near the top of this file.
    if is_block_method
        retainedElements = ...
            double(n) * double(K_t) + ...
            double(m) * double(K_t) + ...
            double(n) * double(K_e) + ...
            double(m) * double(K_e);
    else
        retainedElements = ...
            double(nall) * double(K_t) + ...
            double(nall) * double(K_t) + ...
            double(nall) * double(K_e) + ...
            double(nall) * double(K_e);
    end
    retainedBytes = retainedElements * double(n_m) * bytesSingle;

    % The helper must receive these variables explicitly. A local function
    % cannot inspect variables in its caller workspace with WHOS.
    wDat = whos('dat');
    wTheta = whos('Theta');
    residentInputBytes = double(wDat.bytes + wTheta.bytes);

    nWorkers = 1;
    try
        p = gcp('nocreate');
        if ~isempty(p)
            nWorkers = max(1,double(p.NumWorkers));
        end
    catch
    end

    % parfor may have several completed q vectors waiting to be
    % deserialized at once. Model output qk is normally double precision.
    returnedBytes = double(nq) * bytesDouble * nWorkers;

    % Conservative allowance for serialized copies of inputs, returned
    % values, Qfull assignment, model temporaries, and fragmentation.
    estimatedPeak = ...
        2.50 * qfullBytes + ...
        2.50 * returnedBytes + ...
        0.50 * residentInputBytes + ...
        0.35 * retainedBytes;

    availableBytes = local_available_memory_bytes();

    if isfinite(availableBytes) ...
            && availableBytes > 0
        % Preserve enough room for MATLAB, the GUI, Java, plotting, and
        % brief allocation spikes. Use at least 1.5 GB and otherwise 35%.
        requiredHeadroom = max(1.5 * 1024^3,0.35 * availableBytes);
        tf = estimatedPeak + requiredHeadroom < availableBytes;
    else
        requiredHeadroom = 1.5 * 1024^3;
        % Unknown-memory platforms use a deliberately cautious fallback.
        tf = estimatedPeak < 1.25 * 1024^3;
    end

    info = struct();
    info.availableGB = availableBytes / 1024^3;
    info.estimatedPeakGB = estimatedPeak / 1024^3;
    info.requiredHeadroomGB = requiredHeadroom / 1024^3;
    info.qfullGB = qfullBytes / 1024^3;
    info.retainedGB = retainedBytes / 1024^3;
    info.residentInputGB = residentInputBytes / 1024^3;
    info.workerReturnGB = returnedBytes / 1024^3;
    info.nWorkers = nWorkers;
end

function [loss_value,out] = local_run_crr_backend( ...
    crr_backend,x,mdl,dat,ode,loss,request)
%LOCAL_RUN_CRR_BACKEND Use the backend selected by prepare_crr_backend.

    switch crr_backend
        case 'cpp'
            [loss_value,out] = crr_model_cpp( ...
                x,mdl,dat,ode,loss,request);
        case 'matlab'
            [loss_value,out] = crr_model( ...
                x,mdl,dat,ode,loss,request);
        otherwise
            error('pproc_SAGE:UnknownCRRBackend', ...
                'Unknown CRR backend: %s.',crr_backend);
    end
end

function tf = local_complete_forcing(datk)
%LOCAL_COMPLETE_FORCING Require a complete finite meteorological record.

    tf = isstruct(datk) && isfield(datk,'meteo') ...
        && isstruct(datk.meteo);
    if ~tf
        return
    end
    names = {'P','Ep','T'};
    n = nan(1,numel(names));
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(datk.meteo,name) ...
                || ~isnumeric(datk.meteo.(name)) ...
                || ~isvector(datk.meteo.(name)) ...
                || isempty(datk.meteo.(name)) ...
                || any(~isfinite(datk.meteo.(name)(:)))
            tf = false;
            return
        end
        n(i) = numel(datk.meteo.(name));
    end
    tf = all(n == n(1));
end

function S = local_fdc_score_matrix(D,D0)
%LOCAL_FDC_SCORE_MATRIX Convert raw FDC divergence to skill.
%
% Use the same definition as pmetrics:
%   S_fdc = 1 - D_fdc/D0

    D = double(D);
    D0 = double(D0(:));
    S = nan(size(D));
    if isempty(D) ...
            || isempty(D0) ...
            || size(D,1) ~= numel(D0)
        return
    end
    for j = 1:size(D,2)
        ok = isfinite(D(:,j)) & isfinite(D0) ...
            & D(:,j) >= 0 & D0 > 0;
        S(ok,j) = 1 - D(ok,j) ./ D0(ok);
    end
end

function bytes = local_available_memory_bytes()
%LOCAL_AVAILABLE_MEMORY_BYTES Best-effort available-memory query.

    %bytes = NaN;

    if ispc
        try
            u = memory;
            if isfield(u,'MemAvailableAllArrays')
                bytes = double(u.MemAvailableAllArrays);
                return
            end
            if isfield(u,'MaxPossibleArrayBytes')
                bytes = double(u.MaxPossibleArrayBytes);
                return
            end
        catch
        end
    end

    if isunix && ~ismac
        try
            txt = fileread('/proc/meminfo');
            tok = regexp(txt,'MemAvailable:\s+(\d+)\s+kB', ...
                'tokens','once');
            if ~isempty(tok)
                bytes = str2double(tok{1}) * 1024;
                return
            end
        catch
        end
    end

    try
        rt = java.lang.Runtime.getRuntime();
        bytes = double(rt.maxMemory() - ...
            (rt.totalMemory() - rt.freeMemory()));
    catch
        bytes = NaN;
    end
end

% ----------------------
% local helper functions
% ----------------------
function report_err_block(tag,errv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REPORT_ERR_BLOCK Prints consistency diagnostics for total sum of squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    errv = errv(isfinite(errv));
    if isempty(errv)
        return
    end
    
    mx = 100*max(abs(errv));
    
    if mx < 0.1
        fprintf(['      pproc_SAGE: Stored and ' ...
            'recomputed TSS agree for %s\n'],tag);
    else
        fprintf(['      pproc_SAGE: Recomputed ' ...
            'TSS on pproc support differs from ' ...
            'stored TSS by at most %5.4f%% ' ...
            'for %s\n'],mx,tag);
    end
end


function id = local_normalize_runtime_gauge_id(id)
%LOCAL_NORMALIZE_RUNTIME_GAUGE_ID Normalize known dataset wrappers.
% This comparison is diagnostic only; it does not reorder basins.

    id = upper(strtrim(string(id(:))));
    id = regexprep(id,'\.0+$','');

    % CAMELS-CZ time-series files use camelscz_<code>, whereas basin lists
    % and some attribute sources return the bare code. Treat both forms 
    % identically.
    id = regexprep(id,'^CAMELSCZ[_-]?','');
end

function [NSE,KGE,JKGE,KGE_r,KGE_alpha,KGE_beta, ...
    JKGE_M,JKGE_V,JKGE_C,D_fdc,err] = local_metric_block( ...
    Qfull,dat,loss,id_split,id_bas,which_split, ...
    is_local_split,id_train_global)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL_METRIC_BLOCK Calculates NSE/KGE/JKGE block values for one scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 7
        is_local_split = false;
    end
    if nargin < 8
        id_train_global = [];
    end
    
    nb = numel(id_bas);
    NSE = nan(nb,1);
    KGE = nan(nb,1);
    JKGE = nan(nb,1);
    KGE_r = nan(nb,1);
    KGE_alpha = nan(nb,1);
    KGE_beta = nan(nb,1);
    JKGE_M = nan(nb,1);
    JKGE_V = nan(nb,1);
    JKGE_C = nan(nb,1);
    D_fdc = nan(nb,1);
    doJKGE = isfield(loss,'fnc') ...
        && double(loss.fnc) == 7;
    if isfield(loss,'M') ...
            && ~isempty(loss.M)
        Mdef = loss.M;
    else
        Mdef = 2;
    end
    
    if ~(isscalar(Mdef) ...
            && isnumeric(Mdef) ...
            && isfinite(Mdef) ...
            && any(double(Mdef) == [1 2]))
        Mdef = 2;
    end
    
    Mdef = double(Mdef);
    
    err = nan(nb,1);
    
    for j = 1:nb
        k = id_bas(j);
    
        if k < 1 ...
                || k > numel(dat)
            continue
        end
        if ~isfield(dat{k},'y_n') ...
                || isempty(dat{k}.y_n)
            continue
        end
    
        yk = dat{k}.y_n(:);
        nq = min(size(Qfull,1),numel(yk));
    
        if is_local_split
            switch lower(which_split)
                case 't'
                    if isfield(dat{k},'id_train') ...
                            && ~isempty(dat{k}.id_train)
                        id_use = dat{k}.id_train(:);
                    else
                        id_use = [];
                    end    
                case 'e'
                    if isfield(dat{k},'id_eval') ...
                            && ~isempty(dat{k}.id_eval)
                        id_use = dat{k}.id_eval(:);
                    else
                        id_use = [];
                    end
            end
        else
            id_use = id_split(:);
        end
        id_use = id_use(id_use >= 1 ...
            & id_use <= nq);
    
        if isempty(id_use)
            continue
        end
    
        if isfield(dat{k},'bad') ...
                && ~isempty(dat{k}.bad)
            badk = dat{k}.bad(:);
            if numel(badk) >= max(id_use)
                id_use = id_use(~badk(id_use));
            end
        end
    
        if isempty(id_use)
            continue
        end
    
        y_use = yk(id_use);
        q_use = Qfull(id_use,k);
    
        good = isfinite(y_use) ...
            & isfinite(q_use);
        y_use = y_use(good);
        q_use = q_use(good);
    
        if isempty(y_use)
            continue
        end

        try
            if isfield(dat{k},'fdc') ...
                    && isfield(dat{k}.fdc,lower(which_split))
                D_fdc(j) = fdc_loss_cached( ...
                    q_use,dat{k}.fdc.(lower(which_split)));
            end
        catch
            D_fdc(j) = NaN;
        end
    
        e_use = y_use - q_use;
        RSS_use = sum(e_use.^2,'omitnan');
        
        den = y_use - mean(y_use,'omitnan');
        TSS_recomp = sum(den.^2,'omitnan');
    
        switch lower(which_split)
    
            case 't'
                mu_use  = NaN;
                std_use = NaN;
    
                if isfield(dat{k},'stats') ...
                        && isfield(dat{k}.stats,'mut')
                    mu_use = dat{k}.stats.mut;
                end
                if isfield(dat{k},'stats') ...
                        && isfield(dat{k}.stats,'stdt')
                    std_use = dat{k}.stats.stdt;
                end
    
                [KGE(j),KGE_r(j),KGE_alpha(j), ...
                    KGE_beta(j)] = kge(y_use,q_use, ...
                    mu_use,std_use);
    
                if doJKGE ...
                        && isfield(dat{k},'jkge') ...
                        && isfield(dat{k}.jkge,'m_y') ...
                        && numel(dat{k}.jkge.m_y) >= nq
                    y_jkge = yk(1:nq);
                    q_jkge = Qfull(1:nq,k);
                    m_jkge = dat{k}.jkge.m_y(1:nq);
                    
                    id_jkge = id_use(:);
                    id_jkge = id_jkge(id_jkge >= 1 ...
                        & id_jkge <= nq);
                    
                    method_jkge = double(loss.method);
                    
                    switch method_jkge
                        case {1,2}
                            aux_jkge = loss.n_win;
                        case 3
                            aux_jkge = [];
                        case 4
                            aux_jkge = loss.meta.mo_all;
                    end
    
                    if numel(id_jkge) >= 2
                        try
                            if is_local_split
                                if isfield(dat{k},'id_train') ...
                                        && ~isempty(dat{k}.id_train)
                                    id_bench = dat{k}.id_train(:);
                                else
                                    id_bench = [];
                                end
                            else
                                id_bench = id_train_global(:);
                            end
                            id_bench = id_bench(id_bench >= 1 ...
                                & id_bench <= nq);
                            if isfield(dat{k},'bad') ...
                                    && ~isempty(dat{k}.bad)
                                badk = dat{k}.bad(:);
                                if numel(badk) >= max(id_bench)
                                    id_bench = id_bench(~badk(id_bench));
                                end
                            end
                            id_bench = id_bench(isfinite(y_jkge(id_bench)) ...
                                & isfinite(q_jkge(id_bench)) ...
                                & isfinite(m_jkge(id_bench)));
    
                            [JKGE_bench,~,m_q,M_bench, ...
                                V_bench,C_bench] = jkge_grad( ...
                                y_jkge,q_jkge,m_jkge, ...
                                id_bench,method_jkge,aux_jkge, ...
                                dat{k}.jkge.cache,Mdef);
                            if strcmpi(which_split,'t')
                                JKGE(j) = JKGE_bench;
                                JKGE_M(j) = M_bench;
                                JKGE_V(j) = V_bench;
                                JKGE_C(j) = C_bench;
                            else
                                [JKGE(j),JKGE_M(j), ...
                                    JKGE_V(j),JKGE_C(j)] = ...
                                    jkge_score_given_mq( ...
                                    y_jkge,q_jkge,m_jkge, ...
                                    m_q,id_jkge,loss);
                            end
                        catch ME
                            fprintf(['JKGE recompute failed ' ...
                                'basin %d split %s: %s\n'], ...
                                k,which_split,ME.message);
                            JKGE(j) = NaN;
                        end                        
                    end
                end
    
                TSSt_use = NaN;
                if isfield(dat{k},'stats') ...
                        && isfield(dat{k}.stats,'TSSt')
                    TSSt_use = dat{k}.stats.TSSt;
                end
    
                % ------------------------------------
                TSS_metric = TSSt_use;
                if ~(isfinite(TSS_metric) ...
                        && TSS_metric > 0)
                    TSS_metric = TSS_recomp;
                end
                
                if isfinite(TSS_metric) ...
                        && TSS_metric > 0
                    NSE(j) = 1 - RSS_use / TSS_metric;
                end
                % ------------------------------------
                if isfinite(TSSt_use) ...
                        && TSSt_use ~= 0
                    err(j) = (TSSt_use - TSS_recomp) / TSSt_use;
                end
    
            case 'e'
                mu_use = NaN;
                std_use = NaN;
    
                if isfield(dat{k},'stats') ...
                        && isfield(dat{k}.stats,'mue')
                    mu_use = dat{k}.stats.mue;
                end
                if isfield(dat{k},'stats') ...
                        && isfield(dat{k}.stats,'stde')
                    std_use = dat{k}.stats.stde;
                end
    
                [KGE(j),KGE_r(j),KGE_alpha(j), ...
                    KGE_beta(j)] = kge(y_use,q_use, ...
                    mu_use,std_use);
    
                if doJKGE ...
                        && isfield(dat{k},'jkge') ...
                        && isfield(dat{k}.jkge,'m_y') ...
                        && numel(dat{k}.jkge.m_y) >= nq
                    y_jkge = yk(1:nq);
                    q_jkge = Qfull(1:nq,k);
                    m_jkge = dat{k}.jkge.m_y(1:nq);
                    
                    id_jkge = id_use(:);
                    id_jkge = id_jkge(id_jkge >= 1 ...
                        & id_jkge <= nq);
                    
                    method_jkge = double(loss.method);
    
                    switch method_jkge
                        case {1,2}
                            aux_jkge = loss.n_win;
                        case 3
                            aux_jkge = [];
                        case 4
                            aux_jkge = loss.meta.mo_all;
                    end
    
                    if numel(id_jkge) >= 2
                        try
                            if is_local_split
                                if isfield(dat{k},'id_train') ...
                                        && ~isempty(dat{k}.id_train)
                                    id_bench = dat{k}.id_train(:);
                                else
                                    id_bench = [];
                                end
                            else
                                id_bench = id_train_global(:);
                            end
                            id_bench = id_bench(id_bench >= 1 ...
                                & id_bench <= nq);
                            if isfield(dat{k},'bad') ...
                                    && ~isempty(dat{k}.bad)
                                badk = dat{k}.bad(:);
                                if numel(badk) >= max(id_bench)
                                    id_bench = id_bench(~badk(id_bench));
                                end
                            end
                            id_bench = id_bench(isfinite(y_jkge(id_bench)) ...
                                & isfinite(q_jkge(id_bench)) ...
                                & isfinite(m_jkge(id_bench)));
                            
                            [JKGE_bench,~,m_q,M_bench, ...
                                V_bench,C_bench] = jkge_grad( ...
                                y_jkge,q_jkge,m_jkge, ...
                                id_bench,method_jkge,aux_jkge, ...
                                dat{k}.jkge.cache,Mdef);
                            
                            if strcmpi(which_split,'t')
                                JKGE(j) = JKGE_bench;
                                JKGE_M(j) = M_bench;
                                JKGE_V(j) = V_bench;
                                JKGE_C(j) = C_bench;
                            else
                                [JKGE(j),JKGE_M(j), ...
                                    JKGE_V(j),JKGE_C(j)] = ...
                                    jkge_score_given_mq( ...
                                    y_jkge,q_jkge,m_jkge, ...
                                    m_q,id_jkge,loss);
                            end
                        catch ME
                            fprintf(['JKGE recompute failed ' ...
                                'basin %d split %s: %s\n'], ...
                                k,which_split,ME.message);
                            JKGE(j) = NaN;
                        end                        
                    end
                end
                
                TSSe_use = NaN;
                if isfield(dat{k},'stats') ...
                        && isfield(dat{k}.stats,'TSSe')
                    TSSe_use = dat{k}.stats.TSSe;
                end
    
                % ------------------------------------
                TSS_metric = TSSe_use;
                if ~(isfinite(TSS_metric) ...
                        && TSS_metric > 0)
                    TSS_metric = TSS_recomp;
                end
                
                if isfinite(TSS_metric) ...
                        && TSS_metric > 0
                    NSE(j) = 1 - RSS_use / TSS_metric;
                end
                % ------------------------------------
    
                if isfinite(TSSe_use) ...
                        && TSSe_use ~= 0
                    err(j) = (TSSe_use - TSS_recomp) / TSSe_use;
                end
    
            otherwise
                error(['pproc_SAGE: which_split ' ...
                    'must be ''t'' or ''e''.']);
        end
    end
end

function local_compare_metric_block(tag, ...
    met_val,recomp_val)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL_COMPARE_METRIC_BLOCK Compares supplied and recomputed metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    good = isfinite(met_val) ...
        & isfinite(recomp_val);
    if ~any(good)
        return
    end
    
    mx = max(abs(met_val(good) ...
        - recomp_val(good)));
    
    if mx < 1e-8
        fprintf(['      pproc_SAGE: ' ...
            'Perfect match supplied vs ' ...
            'recomputed %s\n'],tag);
    else
        fprintf(['      pproc_SAGE: Maximum ' ...
            'difference of %10.4e in %s\n'], ...
            mx,tag);
    end
end

function tag = local_tag_string(basin_txt,split_txt,sp_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL_TAG_STRING Creates scenario label string for diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if any(strcmpi(sp_method,{ ...
            'traditional_block', ...
            'manual', ...
            'block', ...
            'deterministic_block'}))
        word = 'period';
    else
        word = 'mask';
    end
    
    tag = sprintf('%s / %s %s', ...
        basin_txt,split_txt,word);
end

function [id_param,range_id,Pn] = ...
    read_site_param_normalized_local(dirres,mname,dt_tag,loss_fnc,part)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_SITE_PARAM_NORMALIZED_LOCAL Read SITE normalized parameters
% from new xlsx or old csv format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 5 || isempty(part)
        part = 'site';
    end
    
    file_xlsx = fullfile(dirres, ...
        sprintf('param_%s_%s.xlsx',mname,dt_tag));
    file_csv = fullfile(dirres, ...
        sprintf('param_%s_%s.csv',mname,dt_tag));
    
    if isfile(file_xlsx)
    
        req_sheet = loss_sheet_name_local(loss_fnc);
    
        if strcmp(part,'site')
            [T,sheet_used] = ...
                sage_readtable_sheet_site_fallback( ...
                file_xlsx,req_sheet);
            if ~strcmpi(sheet_used,req_sheet)
                fprintf(['      pproc_SAGE: ' ...
                    'SITE workbook %s does not ' ...
                    'contain sheet %s --> using ' ...
                    'sheet %s instead\n'], ...
                    file_xlsx,req_sheet, ...
                    sheet_used);
            end
        else
            T = sage_readtable_sheet(file_xlsx,req_sheet);
            sheet_used = req_sheet;
        end
    
        vars = string(T.Properties.VariableNames);
    
        if ~ismember("gauge_ID",vars)
            error(['pproc_SAGE: ' ...
                'Missing gauge_ID in %s ' ...
                '(sheet %s).'],file_xlsx, ...
                sheet_used);
        end
        if ~ismember("range_id",vars)
            error(['pproc_SAGE: Missing ' ...
                'range_id in %s ' ...
                '(sheet %s).'],file_xlsx, ...
                sheet_used);
        end
    
        theta_vars = vars(startsWith(vars,"theta_"));
        theta_vars = sort_theta_vars_local(theta_vars);
    
        if isempty(theta_vars)
            error(['pproc_SAGE: No theta_j ' ...
                'columns found in %s ' ...
                '(sheet %s).'],file_xlsx, ...
                sheet_used);
        end
    
        id_param = T.gauge_ID;
        range_id = T.range_id;
        Pn = nan(numel(theta_vars),height(T));
        for j = 1:numel(theta_vars)
            vname = char(theta_vars(j));
            x = T.(vname);
            
            if isnumeric(x)
                x = double(x(:));
            elseif iscell(x)
                x = str2double(string(x(:)));
            elseif isstring(x) || ischar(x) || iscategorical(x)
                x = str2double(string(x(:)));
            else
                error(['read_site_param_normalized_local:' ...
                    'UnknownType'], ...
                    'Column %s has unsupported type %s.', ...
                    vname,class(x));
            end
            
            if any(~isfinite(x))
                error(['read_site_param_normalized_local:' ...
                    'NonNumericTheta'], ...
                    ['Column %s contains non-numeric or ' ...
                    'missing values.'],vname);
            end       
            Pn(j,:) = x.';
        end
    
    elseif isfile(file_csv)
        T = sage_readtable_any(file_csv);
    
        vars = string(T.Properties.VariableNames);
    
        if ~ismember("gauge_ID",vars)
            error(['pproc_SAGE: ' ...
                'Missing gauge_ID in %s.'], ...
                file_csv);
        end
        if ~ismember("range_id",vars)
            error(['pproc_SAGE: ' ...
                'Missing range_id in %s.'], ...
                file_csv);
        end
    
        theta_vars = vars(startsWith(vars,"theta_"));
        theta_vars = sort_theta_vars_local(theta_vars);
    
        if isempty(theta_vars)
            error(['pproc_SAGE: ' ...
                'No theta_j columns found in %s.'], ...
                file_csv);
        end
    
        id_param = T.gauge_ID;
        range_id = T.range_id;
        Pn = nan(numel(theta_vars),height(T));
        for j = 1:numel(theta_vars)
            Pn(j,:) = T.(theta_vars(j)).';
        end
    
    else
        error(['pproc_SAGE: ' ...
            'Missing SITE parameter file. ' ...
            'Neither %s nor %s exists.'], ...
            file_xlsx,file_csv);
    end
end

function T = sage_readtable_any(fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_READTABLE_ANY Read text/csv file as table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist('readtable','file') ~= 2
        error(['pproc_SAGE: ' ...
            'readtable is required to read %s.'], ...
            fname);
    end
    
    try
        if exist('detectImportOptions','file') == 2
            opts = detectImportOptions(fname);
            T = readtable(fname,opts);
        else
            T = readtable(fname);
        end
    catch ME
        error(['pproc_SAGE: ' ...
            'Could not read %s (%s).'],fname, ...
            ME.message);
    end
end

function T = sage_readtable_sheet(fname, ...
    sheet_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_READTABLE_SHEET Read one workbook sheet as table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist('readtable','file') ~= 2
        error(['pproc_SAGE: ' ...
            'readtable is required to read %s.'], ...
            fname);
    end
    
    try
        if exist('detectImportOptions','file') == 2
            opts = detectImportOptions(fname, ...
                'Sheet',sheet_name);
            T = readtable(fname,opts, ...
                'Sheet',sheet_name);
        else
            T = readtable(fname,'Sheet', ...
                sheet_name);
        end
    catch ME
        error(['pproc_SAGE: ' ...
            'Could not read sheet %s from %s (%s).'], ...
            sheet_name,fname,ME.message);
    end
end

function theta_vars = sort_theta_vars_local(theta_vars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SORT_THETA_VARS_LOCAL Sort theta_1, theta_2, ..., theta_d numerically
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(theta_vars)
        return
    end
    
    num = nan(size(theta_vars));
    for j = 1:numel(theta_vars)
        tok = regexp(theta_vars(j), ...
            '^theta_(\d+)$','tokens','once');
        if ~isempty(tok)
            num(j) = str2double(tok{1});
        end
    end
    
    [~,ord] = sort(num);
    theta_vars = theta_vars(ord);
end

function s = loss_sheet_name_local(loss_fnc)
%%%%%%%%%%%%%%%%%%%%%%
%LOSS_SHEET_NAME_LOCAL
%%%%%%%%%%%%%%%%%%%%%%

    switch loss_fnc
        case 1
            s = 'SAR';
        case 2
            s = 'GLS';
        case 3
            s = 'NSE';
        case 4
            s = 'KGE';
        case 5
            s = 'Huber';
        case 6
            s = 'FDC';
        case 7
            s = 'JKGE';
        otherwise
            error(['pproc_SAGE: ' ...
                'Unknown loss function: %g.'], ...
                loss_fnc);
    end
end

function mname = local_model_name(model_names,im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL_MODEL_NAME Return model name robustly for cell or string input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if iscell(model_names)
        mname = char(model_names{im});
    else
        mname = char(string(model_names(im)));
    end
end

function [th_min_mat,th_max_mat] = ...
    read_site_param_ranges_local(dirres, ...
    mname,dt_tag,range_id,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_SITE_PARAM_RANGES_LOCAL Read parameter ranges for each basin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    file_xlsx = fullfile(dirres, ...
        sprintf('param_ranges_%s.xlsx', ...
        dt_tag));
    
    file_csv_model = fullfile(dirres, ...
        sprintf('param_ranges_%s_%s.csv', ...
        mname,dt_tag));
    
    file_csv_global = fullfile(dirres, ...
        sprintf('param_ranges_%s.csv', ...
        dt_tag));
    
    if isfile(file_xlsx)
        T = sage_readtable_sheet(file_xlsx, ...
            mname);
        fmt = 'xlsx';
    
    elseif isfile(file_csv_model)
        T = sage_readtable_any(file_csv_model);
        fmt = 'csv_model';
    
    elseif isfile(file_csv_global)
        T = sage_readtable_any(file_csv_global);
        fmt = 'csv_global';
    
    else
        error(['pproc_SAGE: ' ...
            'Missing parameter range file for model %s. ' newline ...
            '      %s' newline ...
            '      %s' newline ...
            '      %s'], ...
            mname,file_xlsx,file_csv_model, ...
            file_csv_global);
    end
    
    vars = string(T.Properties.VariableNames);
    
    need_basic = ["range_id","th_min","th_max"];
    for i = 1:numel(need_basic)
        if ~ismember(need_basic(i),vars)
            error(['pproc_SAGE: ' ...
                'Parameter range file for %s ' ...
                'is missing column %s.'], ...
                mname,need_basic(i));
        end
    end
    
    nb = numel(range_id);
    th_min_mat = nan(d,nb);
    th_max_mat = nan(d,nb);
    
    for k = 1:nb
        rid = range_id(k);
    
        rows = T(T.range_id == rid,:);
        if isempty(rows)
            error(['pproc_SAGE: ' ...
                'range_id = %g not found in ' ...
                'parameter range file for %s.'], ...
                rid,mname);
        end
    
        switch fmt
            case {'csv_model','csv_global'}
                if ~ismember("n_par",vars)
                    error(['pproc_SAGE: ' ...
                        'Legacy parameter range file for %s ' ...
                        'is missing column n_par.'],mname);
                end
                [~,ord] = sort(rows.n_par);
                rows = rows(ord,:);
    
            case 'xlsx'
                % preserve workbook row order
    
            otherwise
                error(['pproc_SAGE: ' ...
                    'Unknown parameter range file format ' ...
                    'for model %s.'],mname);
        end
    
        if height(rows) ~= d
            error(['pproc_SAGE: ' ...
                'range_id = %g for %s has %d rows, ' ...
                'but model expects d = %d parameters.'], ...
                rid,mname,height(rows),d);
        end
    
        th_min_mat(:,k) = rows.th_min;
        th_max_mat(:,k) = rows.th_max;
    end
end

function [tf,file_used,file_kind] = ...
    site_param_file_exists_local(dirres,mname,dt_tag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SITE_PARAM_FILE_EXISTS_LOCAL Check whether SITE parameter file exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    file_xlsx = fullfile(dirres, ...
        sprintf('param_%s_%s.xlsx', ...
        mname,dt_tag));
    file_csv = fullfile(dirres, ...
        sprintf('param_%s_%s.csv', ...
        mname,dt_tag));
    
    if isfile(file_xlsx)
        tf = true;
        file_used = file_xlsx;
        file_kind = 'xlsx';
    elseif isfile(file_csv)
        tf = true;
        file_used = file_csv;
        file_kind = 'csv';
    else
        tf = false;
        file_used = '';
        file_kind = '';
    end
end

function [T,sheet_used] = ...
    sage_readtable_sheet_site_fallback(fname,req_sheet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_READTABLE_SHEET_SITE_FALLBACK Read requested SITE sheet if present
% otherwise fall back to a usable metric sheet.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cands = [string(req_sheet), ...
             "NSE","KGE","SAR", ...
             "GLS","Huber","FDC"];
    
    [~,ia] = unique(lower(cands),'stable');
    cands = cands(ia);
    
    for i = 1:numel(cands)
        try
            T_try = sage_readtable_sheet(fname, ...
                char(cands(i)));
            if local_is_valid_theta_table(T_try)
                T = T_try;
                sheet_used = char(cands(i));
                return
            end
        catch
        end
    end
    
    sheets = local_sheetnames_safe(fname);
    for i = 1:numel(sheets)
        try
            T_try = sage_readtable_sheet(fname, ...
                char(sheets{i}));
            if local_is_valid_theta_table(T_try)
                T = T_try;
                sheet_used = char(sheets{i});
                return
            end
        catch
        end
    end
    
    error(['pproc_SAGE: ' ...
        'Could not read a usable ' ...
        'parameter sheet from %s. ' ...
        'Tried requested sheet %s ' ...
        'and standard fallbacks.'], ...
        fname,req_sheet);
end

function tf = local_is_valid_theta_table(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL_IS_VALID_THETA_TABLE Check whether a table contains SITE parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    try
        vars = string(T.Properties.VariableNames);
        has_gauge  = ismember("gauge_ID",vars);
        has_range = ismember("range_id",vars);
        has_theta = any(startsWith(vars,"theta_"));
        tf = has_gauge && has_range && has_theta;
    catch
        tf = false;
    end
end

function sheets = local_sheetnames_safe(fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL_SHEETNAMES_SAFE Return workbook sheet names if supported
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sheets = {};
    
    try
        if exist('sheetnames','file') == 2
            s = sheetnames(fname);
            if isstring(s)
                sheets = cellstr(s(:));
            elseif iscell(s)
                sheets = s(:);
            end
        elseif exist('xlsfinfo','file') == 2
            [~,s] = xlsfinfo(fname);
            if iscell(s)
                sheets = s(:);
            end
        end
    catch
        sheets = {};
    end
end
