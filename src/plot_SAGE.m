function varargout = plot_SAGE(part,mdl,dat,bas,prd,Q,curr, ...
    Qfdc,region,tTheta,latlon,nTheta,At,An,gaugescen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT_SAGE Makes figures of SAGE results
%
% SYNOPSIS: plot_SAGE(part,mdl,dat,bas,prd,Q,curr,Qfdc,region, ...
%   tTheta,latlon,nTheta,At,An,gaugescen)
%  zoneTbl = plot_SAGE(part,mdl,dat,bas,prd,Q,curr,Qfdc, ...
%   region,tTheta,latlon,nTheta,At,An,gaugescen)
%   part        which data are analyzed?
%     'site'    single-site training results
%     'sage'    SAGE training results
%     'sage_multi' combined SAGE output from multiple model workspaces;
%                  used by the paper hydrograph composer
%   mdl         structure with model state/parameter information
%    .model      choice of model
%                 1 hymod
%                 2 hmodel
%                 3 sacsma
%                 4 xinanjiang
%                 5 gr4jA
%                 6 hbv
%                 7 gr4jB [analytic routing]
%    .mcode      scalar with numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler MATLAB
%                 4 Runge Kutta implementation C++
%    .calc       model execution
%      'seq'      sequential execution of watersheds
%      'par'      parallel execution of watersheds
%    .mode       assessment design
%                 1 = training basins only | training period only
%                 2 = training basins only | training and evaluation
%                     period/mask
%                 3 = training and evaluation basins | training period only
%                 4 = training and evaluation basins | training and
%                     evaluation period/mask
%    .sp_method  string with split design
%                 'manual'
%                 'deterministic_block'
%                 'traditional_block'
%                 'block'
%                 'random_block'
%                 'random'
%                 'deterministic_kfold'
%                 'random_kfold'
%    .names      list of model names
%    .y0         mx1 vector of initial states
%    .pspace     0 hydrologic, 1 unit cube, 2 unconstrained parameters
%    .th_min     dx1 vector of lower parameter values
%    .th_max     dx1 vector of upper parameter values
%    .par_names  1xd cell with parameter names
%    .tout       model output times
%    .idx        indices training + evaluation periods
%    .id_train   indices training period
%    .id_eval    indices evaluation period
%   dat         cell structure with basin data
%    {k}.gauge    gauge code of basin k
%    {k}.y_n     observed discharge series [normalized]
%    {k}.bad     logical vector with invalid observations
%   bas         structure with basin information
%    .K          total number of watersheds
%    .K_t        number of training watersheds
%    .K_e        number of evaluation watersheds
%    .r          number of basin attributes
%    .gname      gauge names in final basin order [train; eval]
%    .id_gauge    gauge codes in final basin order [train; eval]
%    .zone       structure with hydroclimatic basin classification
%     .id         (K_t+K_e)x1 string array with zone labels per basin
%                 e.g., "humid_rain", "subhumid_snow", "dry_rain"
%     .num        (K_t+K_e)x1 numeric vector with integer zone identifiers
%     .names      mx1 string array with unique zone names
%     .aridity    (K_t+K_e)x1 vector of aridity index values
%     .frac_snow  (K_t+K_e)x1 vector of snow-fraction values
%   prd         structure about training/evaluation/spin-up period
%    .dt         temporal data resolution [1 daily, 24 hourly, 96 at 15-minute]
%    .spinup     spin-up period in days
%    .[others]   depends on mdl.sp_method
%   Q           discharge structure
%    .tt         discharge, training basins | training period
%    .te         discharge, training basins | evaluation period/mask
%    .et         discharge, evaluation basins | training period
%    .ee         discharge, evaluation basins | evaluation period/mask
%   curr        current basin-wise performance structure
%    .NSE        NSE structure
%      .tt       training basins | training period/mask
%      .te       training basins | evaluation period/mask
%      .et       evaluation basins | training period/mask
%      .ee       evaluation basins | evaluation period/mask
%    .KGE        KGE structure
%      .tt/.te/.et/.ee as defined above
%    .S_fdc      flow-duration-curve skill structure
%      .tt/.te/.et/.ee as defined above
%    .JKGE       JKGE structure
%      .tt/.te/.et/.ee as defined above
%   Qfdc        Flow duration curve of gaugescen
%    .id
%    .gauge
%    .req       
%    .qy        cell structure simulated/measured q
%   region      CAMELS data region (US, GB, BZ, AU, etc.)
%   tTheta        OPTIONAL: parameter-trace summaries
%               preferred size: (i_max)xdx7: 5,15,25,50,75,85,95%
%               Set [] if no trace figure is desired
%   latlon      OPTIONAL: (K_t+K_e)x2 matrix lat/long (°) selected basins
%   nTheta         OPTIONAL: normalized hydrologic parameter values
%    'site'      dxKxn_m array, one layer per model
%    'sage'      dx(K_t+K_e)x1 array for the submitted SAGE model
%   At          OPTIONAL: d x K_t matrix time-weighted parameter
%                attribution values retained for training watersheds only
%   An          OPTIONAL: d x K_t matrix net gradient-based
%                attribution values retained for training watersheds only
%   gaugescen    OPTIONAL: basin selection structure
%    .tt         requested basins for training basins | training period
%    .te         requested basins for training basins | evaluation period
%    .et         requested basins for evaluation basins | training period
%    .ee         requested basins for evaluation basins | evaluation period
%    .train      applied to training-basin plots
%    .eval       applied to evaluation-basin plots
%
% OUTPUT:
%   zoneTbl     OPTIONAL table with zone/scenario/model summaries for the
%               GUI postprocessor table. Columns include Zone, ZoneName,
%               Scenario, ScenarioName, Model, n, T_NSE, T_KGE, T_JKGE, and
%               S_IB.
%
% KEY PHILOSOPHY: 
%  global split: 
%    use mdl.id_train / mdl.id_eval
%  local rainfall split:
%    use dat{k}.id_train / dat{k}.id_eval
%    store Q on full scored window 
%
% NOTES:
%   1. This version uses one notation only
%        tt = training basins | training period
%        te = training basins | evaluation period/mask
%        et = evaluation basins | training period
%        ee = evaluation basins | evaluation period/mask
%   2a. Figure 2:
%        1x4 ECDF panel of NSE for tt | te | et | ee
%   2b. Figure 3:
%        1x4 ECDF panel of KGE for tt | te | et | ee
%   2c. Figure 4:
%        1x4 ECDF panel of JKGE for tt | te | et | ee
%   2d. Figure 5:
%        1x4 ECDF panel of NSE for hydroclimatic zone
%   2e. Figure 6:
%        1x4 ECDF panel of KGE for hydroclimatic zone
%   2f. Figure 7:
%        1x4 ECDF panel of JKGE for hydroclimatic zone
%   3a. Figure 11:
%        parameter traces
%   3b. Figure 12:
%        regional parameter maps
%   3c. Figure 13:
%        regional variogram figures
%   3d. Figures 14-15
%        parameter attribution  
%   4.  Figure 20+: 
%        time-series figures generated separately for each scenario
%   5.  Figure 200+: 
%        flow duration curves generated separately for each scenario
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% © Written by Jasper A. Vrugt, Feb. 2026 / updated Aug. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NSE = curr.NSE;
    KGE = curr.KGE;
    S_fdc = curr.S_fdc;
    JKGE = curr.JKGE;

    % defaults for optional inputs
    if nargin < 8 ...
            || isempty(Qfdc)
        Qfdc = [];
    end
    if nargin < 9 ...
            || isempty(region)
        region = 'US';
    end
    try
        region = region_helpers('code',region);
    catch
    end
    if nargin < 10 ...
            || isempty(tTheta)
        tTheta = [];
    end
    if nargin < 11 ...
            || isempty(latlon)
        latlon = [];
    end
    if nargin < 12 ...
            || isempty(nTheta)
        nTheta = [];
    end
    if nargin < 13 ...
            || isempty(At)
        At = [];
    end
    if nargin < 14 ...
            || isempty(An)
        An = [];
    end
    if nargin < 15 ...
            || isempty(gaugescen)
        gaugescen = struct();
    end
    
    Bprt = local_default_Bprt_from_USGSscen_Qfdc( ...
        gaugescen,Qfdc);

    if strcmpi(part,'sage') ...
            && ~isempty(tTheta) ...
            && ndims(tTheta) ~= 3
        error(['      Error: plot_SAGE: ' ...
            'tTheta must be empty or a 3-D array ' ...
            'of size [i_max x d x 7].']);
    end
    if strcmpi(part,'site')
        mdl.model = [];
    end
    
    n_char = 50;        
    plot_map = 2;       % 1 = full map in gray
                        % 2 = country borders only parameter maps
                        
    % ---------- daily vs hourly labeling ----------
    dt = 1;
    if isfield(prd,'dt') ...
            && ~isempty(prd.dt)
        dt = prd.dt;
    end
    
    if dt == 24
        dt_str = 'hourly';
        q_unit = 'mm/h';
        t_unit = 'hours';
    elseif dt == 96
        dt_str = '15-minute';
        q_unit = 'mm/15\,min';
        t_unit = '15-minute intervals';
    else
        dt_str = 'daily';
        q_unit = 'mm/d';
        t_unit = 'days';
    end
    dt_fig_str = local_capfirst(dt_str);
    
    K = bas.K;
    K_t = bas.K_t;
    K_e = bas.K_e;
    
    if isfield(mdl,'sp_method') ...
            && ~isempty(mdl.sp_method)
        sp_method = lower( ...
            string(mdl.sp_method));
    else
        sp_method = "manual";
    end
    
    is_rainfall_block = strcmpi(sp_method, ...
        'rainfall_block');
    
    if is_rainfall_block
        mdl.local = 1;
    end
    
    if is_rainfall_block
        samp_word = 'split';
    elseif any(strcmpi(sp_method, ...
            {'manual', ...
            'traditional_block', ...
            'block', ...
            'deterministic_block'}))
        samp_word = 'period';
    else
        samp_word = 'mask';
    end
    
    % -------- resolution --------
    ppi = get(0,'ScreenPixelsPerInch');
    scr = get(0,'ScreenSize');
    figW = scr(3) / ppi;
    figH = scr(4) / ppi;            % 92.5% of screen height
    
    % ---------- colors ----------
    colors = [ ...
        214/255  39/255  40/255; ... % red    - hymod
         44/255 160/255  44/255; ... % green  - hmodel
        148/255 103/255 189/255; ... % purple - sacsma
        140/255  86/255  75/255; ... % brown  - Xinanjiang
         31/255 119/255 180/255; ... % blue   - gr4jA
        255/255 127/255  14/255; ... % orange - hbv
        0.0000 0.7000 0.7000;   ... % teal    — gr4jB
        0.9000 0.0000 0.9000;   ... % magenta — cfe_nwm
        0.6500 0.6500 0.6500];      % gray    — unknown
    
    fontsize_legend = 16;
    face_alpha = 0.50;
    marker_alpha = 0.70; 
    model_names = mdl.names;
    
    % % soft background shading
    % cShadeTrain = [0.90 0.94 0.98];   % blue-gray
    % cShadeEval = [0.99 0.94 0.88];   % warm orange-gray
    % 
    % % muted data-point colors
    % cPtTrain = [0.20 0.45 0.70];      % muted blue
    % cPtEval = [0.85 0.50 0.20];      % muted orange
    cShadeTrain = [0.94 0.94 0.94];
    cShadeEval = [0.90 0.92 0.96];
    
    cPtTrain = [0.25 0.25 0.25];
    cPtEval = [0.55 0.55 0.55];
    
    if strcmpi(part,'sage')
        n_m = 1;
    elseif strcmpi(part,'sage_multi')
        n_m = local_infer_nmodels(Q,NSE,KGE);
    else
        n_m = local_infer_nmodels(Q,NSE,KGE);
        assert(n_m <= size(colors,1), ...
            sprintf(['      plot_SAGE: ' ...
            'not enough rows in colors ' ...
            'for n_m = %d models.'],n_m));   
    end
    
    id_dat = cell(K,1);
    for k = 1:K
        id_dat{k} = local_basin_code(dat{k}.gauge,region);
    end
    
    switch lower(part)
        case 'site'
            id_model = 1:n_m;
        case 'sage'
            if isempty(mdl.model)
                id_model = 1;
            else
                id_model = mdl.model(1);
            end
        case 'sage_multi'
            id_model = 1:n_m;
        otherwise
            error(['      Error: plot_SAGE: ' ...
                'part must be ''site'', ''sage'', or ''sage_multi''.']);
    end
    
    % ---------- model name used in figure titles ----------
    if strcmpi(part,'sage')
        if isempty(mdl.model)
            model_name_fig = 'SAGE';
        else
            model_name_fig = local_model_display_name( ...
                model_names{mdl.model(1)});
        end
    elseif strcmpi(part,'sage_multi')
        model_name_fig = 'SAGE';
    else
        if n_m == 1 ...
                && ~isempty(id_model)
            model_name_fig = local_model_display_name( ...
                model_names{id_model(1)});
        else
            model_name_fig = 'SITE';
        end
    end
    
    if is_rainfall_block
        local_check_rainfall_indices(dat,K);
        id_tr = [];
        id_ev = [];
    else
        id_tr = expand_index(mdl.id_train);
        if isfield(mdl,'id_eval') ...
                && ~isempty(mdl.id_eval)
            id_ev = expand_index(mdl.id_eval);
        else
            id_ev = [];
        end
    end
    
    % ----------------------------------------------
    % Respect mdl.mode even if Q/NSE/KGE contain all
    % four scenario fields from a previous run
    % ----------------------------------------------
    [Qplot,NSEplot,KGEplot,JKGEplot] = ...
        local_apply_mode_filter( ...
        mdl.mode,Q,NSE,KGE,JKGE);

    % Paper-only hydrograph composition. This path still uses the native
    % plot_SAGE forcing, temperature, discharge, date-axis, and basin
    % resolution helpers, but skips unrelated ECDF/FDC/parameter figures.
    if isfield(gaugescen,'paper_hydrograph') ...
            && ~isempty(gaugescen.paper_hydrograph)
        hPaper = local_plot_paper_hydrograph( ...
            gaugescen.paper_hydrograph,mdl,dat,bas,prd,Qplot, ...
            q_unit,t_unit,n_m,id_model,colors,model_names,id_dat, ...
            sp_method,cShadeTrain,cShadeEval,face_alpha,region);
        if nargout >= 1
            varargout{1} = hPaper;
        end
        return
    end
    
    % =============================
    % Build canonical scenario list
    % =============================
    scenariosNSE = ...
        local_build_metric_scenarios( ...
        NSEplot,Qplot,dt_str,K_t,K_e,id_tr,id_ev, ...
        gaugescen,samp_word);
    
    scenariosKGE = ...
        local_build_metric_scenarios( ...
        KGEplot,Qplot,dt_str,K_t,K_e,id_tr,id_ev, ...
        gaugescen,samp_word);
    
    scenariosJKGE = ...
        local_build_metric_scenarios( ...
        JKGEplot,Qplot,dt_str,K_t,K_e,id_tr,id_ev, ...
        gaugescen,samp_word);

    Sfdcplot = struct('tt',[],'te',[],'et',[],'ee',[]);
    if isstruct(S_fdc)
        Sfdcplot = S_fdc;
        for field = {'tt','te','et','ee'}
            name = field{1};
            if ~isfield(Qplot,name) || isempty(Qplot.(name))
                Sfdcplot.(name) = [];
            end
        end
    end
    scenariosSfdc = ...
        local_build_metric_scenarios( ...
        Sfdcplot,Qplot,dt_str,K_t,K_e,id_tr,id_ev, ...
        gaugescen,samp_word);
    
    % --------------------------------------------------
    % Zone/scenario summary table for GUI postprocessor
    % --------------------------------------------------
    try
        zoneTbl = local_build_zone_summary_table( ...
            bas,scenariosNSE,scenariosKGE,scenariosJKGE, ...
            n_m,id_model,model_names);
    catch ME
        warning('plot_SAGE:zoneSummaryFailed', ...
            ['      Warning: zone summary table failed. ', ...
             'Continuing with figures.\n%s'],ME.message);
        zoneTbl = table();
    end
    
    % =============================
    % Figures: ECDF figure printing
    % =============================
    try 
        % ===========================================
        % Figure 2: 1x4 ECDF panel of NSE by scenario
        % ===========================================
        plot_metric_ecdf_figure(2, ...
            part,mdl, ...
            dt_str,samp_word, ...
            scenariosNSE,n_m,id_model, ...
            model_names, ...
            colors,figW,figH,model_name_fig,'NSE', ...
            'Nash-Sutcliffe Efficiency');
        
        % ===========================================
        % Figure 3: 1x4 ECDF panel of KGE by scenario
        % ===========================================
        plot_metric_ecdf_figure(3, ...
            part,mdl, ...
            dt_str,samp_word, ...
            scenariosKGE,n_m,id_model, ...
            model_names, ...
            colors,figW,figH,model_name_fig,'KGE', ...
            'Kling-Gupta Efficiency');
        
        % Figure 4: JKGE when selected; otherwise S_FDC by scenario
        hasJKGE = any(cellfun(@(s) ...
            local_has_data(s.metric), ...
            scenariosJKGE));
        if hasJKGE
            plot_metric_ecdf_figure(4, ...
                part,mdl, ...
                dt_str,samp_word, ...
                scenariosJKGE,n_m,id_model, ...
                model_names, ...
                colors,figW,figH,model_name_fig, ...
                'JKGE', ...
                'Jawad-Kling-Gupta Efficiency');
        else
            plot_metric_ecdf_figure(4, ...
                part,mdl, ...
                dt_str,samp_word, ...
                scenariosSfdc,n_m,id_model, ...
                model_names, ...
                colors,figW,figH,model_name_fig, ...
                'S_{\rm FDC}', ...
                'Flow-duration-curve skill');
        end
        
        % =======================================
        % Figure 5: 1x4 ECDF panel of NSE by zone
        % =======================================
        zoneIDsForLayout = [];
        if isfield(bas,'zone') && isfield(bas.zone,'num')
            zoneIDsForLayout = unique(double(bas.zone.num(:)));
            zoneIDsForLayout = zoneIDsForLayout( ...
                isfinite(zoneIDsForLayout) & zoneIDsForLayout > 0);
        end
        nZonePages = max(1,ceil(numel(zoneIDsForLayout)/4));
        plot_metric_ecdf_by_zone_figure(5, ...
            part,mdl, ...
            dt_str,samp_word, ...
            scenariosNSE,n_m,id_model, ...
            model_names,colors,figW,figH, ...
            model_name_fig, ...
            'NSE','Nash-Sutcliffe Efficiency',bas);
        
        % =======================================
        % Figure 6: 1x4 ECDF panel of KGE by zone
        % =======================================
        plot_metric_ecdf_by_zone_figure(5+nZonePages, ...
            part,mdl, ...
            dt_str,samp_word, ...
            scenariosKGE,n_m,id_model, ...
            model_names,colors,figW,figH, ...
            model_name_fig, ...
            'KGE','Kling-Gupta Efficiency',bas);
        
        % ========================================
        % Figure 7: zone ECDF of JKGE or raw D_FDC
        % =========================================
        if hasJKGE
            plot_metric_ecdf_by_zone_figure(5+2*nZonePages, ...
                part,mdl, ...
                dt_str,samp_word, ...
                scenariosJKGE,n_m,id_model, ...
                model_names,colors,figW,figH, ...
                model_name_fig, ...
                'JKGE', ...
                'Jawad-Kling-Gupta Efficiency',bas);
        else
            plot_metric_ecdf_by_zone_figure(5+2*nZonePages, ...
                part,mdl, ...
                dt_str,samp_word, ...
                scenariosSfdc,n_m,id_model, ...
                model_names,colors,figW,figH, ...
                model_name_fig, ...
                'S_{\rm FDC}', ...
                'Flow-duration-curve skill',bas);
        end
    
    catch ME
        warning('plot_SAGE:ecdfFailed', ...
            ['      Warning: plot_SAGE: ' ...
            'ECDF figures could not be ', ...
             'created and will be skipped.\n%s'], ...
             ME.message);
    end
    
    % ===========================
    % Figure 11: parameter traces
    % ===========================
    if strcmpi(part,'sage') ...
            && ~isempty(tTheta)
        pspace_trace = 1;   % default normalized trace space
        if isfield(mdl,'pspace') ...
                && ~isempty(mdl.pspace)
            pspace_trace = double(mdl.pspace);
            if ~ismember(pspace_trace,[0 1])
                pspace_trace = 1;
            end
        end   
        plot_parameter_traces_figure(11, ...
            mdl,tTheta,figW,figH,model_name_fig, ...
            pspace_trace);
    end
    
    % ==================================
    % Figure 12: regional parameter maps
    % ==================================
    if strcmpi(part,'sage') ...
            && ~isempty(nTheta) ...
            && ~isempty(latlon)
        try
            if plot_map == 1
                plot_parameter_maps_figure(12, ...
                    mdl,bas,latlon,nTheta,region, ...
                    model_name_fig,figW,figH, ...
                    gaugescen);
            elseif plot_map == 2
                plot_parameter_maps_country_figure(12, ...
                    mdl,bas,latlon,nTheta,region, ...
                    model_name_fig,figW,figH, ...
                    gaugescen);
            end
        catch ME
            warning('plot_SAGE:parameterMapFailed', ...
                ['      Warning: plot_SAGE: ' ...
                'parameter-map figure ' ...
                 'could not be created and ' ...
                 'will be skipped.\n%s'], ...
                ME.message);
        end
    end

    % ===============================
    % Figure 13: parameter variograms
    % ===============================
    if ~isempty(nTheta) ...
            && ~isempty(latlon)
        try
            if strcmpi(part,'sage')
                % SAGE: nTheta is d x K
                plot_parameter_variograms_figure(13, ...
                    mdl,latlon,nTheta,model_name_fig,region)
    
            elseif strcmpi(part,'site')
                % SITE: nTheta is d x K x n_m
                plot_parameter_variograms_SITE_figure(13, ...
                    mdl,latlon,nTheta,region);
            end
    
        catch ME
            warning('plot_SAGE:parameterVariogramFailed', ...
                ['      Warning: plot_SAGE: ' ...
                'parameter-variogram figure ' ...
                 'could not be created and ' ...
                 'will be skipped.\n%s'], ...
                ME.message);
        end
    end
    
    % ===================================
    % Figures 14-15: gradient attribution
    % ===================================
    if strcmpi(part,'sage') ...
            && ~isempty(At) ...
            && ~isempty(An)
        try
            % SAGE attribution:
            % At and An are expected as d x K x i_max arrays. They are
            % initialized with NaNs, so the latest iteration is the largest
            % third-dimension slice with at least one finite entry in At or An.
            %
            % The GUI reduces each d x K slice across basins and then plots
            % parameter x iteration. This reproduces that same convention here.
            plot_parameter_attribution_figure(14, ...
                mdl,At,An,figW,figH, ...
                model_name_fig,'abs');
            plot_parameter_attribution_figure(15, ...
                mdl,At,An,figW,figH, ...
                model_name_fig,'rel');
        catch ME
            warning('plot_SAGE:parameterAttributionFailed', ...
                ['      Warning: plot_SAGE: ' ...
                'gradient-attribution figures ' ...
                 'could not be created and ' ...
                 'will be skipped.\n%s'], ...
                ME.message);
        end
    end
    
    % =========================================
    % Figure 20+: time series by scenario group
    % =========================================
    
    try
        if local_use_four_scenario_timeseries(sp_method)
            local_plot_timeseries_fourgroups( ...
                20,part,mdl,dat,bas,prd,scenariosNSE, ...
                q_unit,t_unit,dt_fig_str,n_m,id_model, ...
                colors,model_names,id_dat,Qfdc,Bprt,n_char, ...
                figW,figH,model_name_fig,samp_word, ...
                fontsize_legend,region,sp_method, ...
                cShadeTrain,cShadeEval,face_alpha);
        else
            local_plot_timeseries_postprocessor( ...
                20,part,mdl,dat,bas,prd,Q,Qfdc,gaugescen, ...
                q_unit,t_unit,dt_fig_str,n_m,id_model, ...
                colors,model_names,id_dat,Bprt,n_char, ...
                figW,figH,model_name_fig,sp_method, ...
                fontsize_legend,cPtTrain,cPtEval, ...
                cShadeTrain,cShadeEval,face_alpha, ...
                marker_alpha,region);
        end
    catch ME
        warning('plot_SAGE:timeSeriesFailed', ...
            ['      Warning: plot_SAGE: ' ...
            'time-series figures ' ...
            'could not be created and ' ...
            'will be skipped.\n%s'], ...
            getReport(ME,'extended', ...
            'hyperlinks','off'));
    end
    
    % =================================
    % Figure 200+: flow duration curves
    % =================================
    try
        local_plot_fdc_postprocessor( ...
            200,part,dt_str,mdl,dat,bas,Qfdc, ...
            model_name_fig,q_unit,figH, ...
            fontsize_legend,samp_word,region);
    catch ME
        warning('plot_SAGE:FDCFailed', ...
            ['      Warning: plot_SAGE: ' ...
            'FDC figures could not be ' ...
            'created and will be skipped.\n%s'], ...
            ME.message);
    end
    % ========================================
    
    try
        setappdata(0,'SAGE_zone_summary_table',zoneTbl);
    catch
    end
    
    if nargout >= 1
        varargout{1} = zoneTbl;
    end

end

% ================
% helper functions
% ================
function plot_metric_ecdf_figure(figNo,part, ...
    mdl,dt_str,samp_word,scenarios,n_m,id_model, ...
    model_names,colors,figW,figH,model_name_fig, ...
    metric_tag,xlab_str)
%PLOT_METRIC_ECDF_FIGURE Plots scenario-specific ECDFs for one metric.

    figW = 0.98 * figW;
    figH = min(figH,7.2);
    
    figure(figNo); clf;
    fig = gcf;
    set(fig,'Name', ...
        sprintf(['%s: Empirical Cumulative ' ...
        'Distribution Functions: %s'], ...
        model_name_fig,metric_tag), ...
        'NumberTitle','off', ...
        'color','w', ...
        'Units','inches', ...
        'Position',[0.5 0.5 figW figH]); %[0.5 0.5 20 5.5]);
    
    % Explicit normalized locations are stable under exportgraphics and
    % PowerPoint rasterization. A tiled layout combined with `axis square`
    % shrank the axes when the wide figure was rendered off screen.
    marginLeft = 0.055;
    marginRight = 0.020;
    panelBottom = 0.180;
    panelHeight = 0.680;
    % Convert the normalized height to a normalized width using the actual
    % figure aspect ratio. This makes the axes physically square while the
    % explicit positions remain stable during PowerPoint export.
    panelWidth = panelHeight*(figH/figW);
    gap = (1-marginLeft-marginRight-4*panelWidth)/3;
    ax = gobjects(4,1);
    for ii = 1:4
        panelLeft = marginLeft + (ii-1)*(panelWidth+gap);
        ax(ii) = axes(fig,'Units','normalized', ...
            'Position',[panelLeft panelBottom panelWidth panelHeight]);
        box(ax(ii),'on');
        hold(ax(ii),'on');
    end
    
    [xL,xR] = local_metric_ecdf_limits(scenarios,metric_tag);
    fill_alpha = 0.20;
    line_width = 1.8;
    med_lw = 1.1;
    med_ms = 6;
    % One typography standard for every ECDF metric. These are the former
    % S_FDC sizes, retained as the preferred appearance for all figures.
    fnt_axis = 18;
    fnt_label = 18;
    fnt_title = 16;
    fnt_med = 16;
    
    panelTitles = ...
        {'train basins | train period', ...
         'train basins | eval period', ...
         'eval basins | train period', ...
         'eval basins | eval period'};
    
    for ii = 1:4
        set(ax(ii),'tickdir','out', ...
            'fontsize',fnt_axis, ...
            'linewidth',1, ...
            'box','off', ...                   %  remove full box
            'Layer','top', ...                 % keep ticks on top of fills
            'ticklength',[0.02 0.02]);
        
        xlim(ax(ii),[xL xR]);
        ylim(ax(ii),[0 1]);
    
        xlabel(ax(ii),xlab_str, ...
            'interpreter','latex', ...
            'fontsize',fnt_label);
    
        ttl = strrep(panelTitles{ii},'_','\_');
        title(ax(ii), ...
            ['\textbf{' ttl '}'], ...
            'interpreter','latex', ...
            'fontsize',fnt_title);
    
        if ii == 1
            ylabel(ax(ii),'ECDF', ...
                'interpreter','latex', ...
                'fontsize',fnt_label);
        else
            set(ax(ii),'YTickLabel',[]);
            ylabel(ax(ii),'');
        end
        % ---- manual top/right axes ----
        xl = xlim(ax(ii));
        yl = ylim(ax(ii));
        
        col = ax(ii).XColor;
        lw = ax(ii).LineWidth;
        
        % top horizontal line
        line(ax(ii),xl,[yl(2) yl(2)], ...
            'color',col,'linewidth',lw, ...
            'clipping','off', ...
            'handlevisibility','off');
        % right vertical line
        line(ax(ii),[xl(2) xl(2)],yl, ...
            'color',col,'linewidth',lw, ...
            'clipping','off', ...
            'handlevisibility','off');
    end
    
    for is = 1:4
        M = scenarios{is}.metric;
    
        if ~local_has_data(M)
            cla(ax(is));
            box(ax(is),'on');
            xlim(ax(is),[xL xR]);
            ylim(ax(is),[0 1]);
    
            if is == 1
                ylabel(ax(is),'ECDF', ...
                    'interpreter','latex', ...
                    'fontsize',fnt_label);
            else
                set(ax(is),'YTickLabel',[]);
            end
    
            xlabel(ax(is),xlab_str, ...
                'interpreter','latex', ...
                'fontsize',fnt_label);
            ht = title(ax(is),panelTitles{is});
            set(ht,'interpreter','none', ...
                'fontsize',fnt_title);
    
            text(ax(is),0,0.5,'Not available', ...
                'horizontalalignment','center', ...
                'verticalalignment','middle', ...
                'interpreter','latex', ...
                'fontsize',20);
            drawnow;
            continue
        end
    
        for mdl_i = 1:n_m
            c = colors(id_model( ...
                min(mdl_i,numel(id_model))),:);
            im = id_model(min(mdl_i,numel(id_model)));
            nam = local_model_display_name(model_names{im});
            nFinite = nnz(isfinite(M(:,mdl_i)));
            lbl = sprintf(['\\texttt{%s}\\; ' ...
                '($n = %d$)'], ...
                local_latex_escape(nam),nFinite);
    
            plot_ecdf_panel(ax(is), ...
                M(:,mdl_i),c,lbl, ...
                part,metric_tag, ...
                scenarios{is}.tag, ...
                xL,xR,fill_alpha, ...
                line_width,med_lw, ...
                med_ms,fnt_med);
        end
    
        legend(ax(is),'show', ...
            'interpreter','latex', ...
            'fontsize',16, ...
            'location','northwest', ...
            'box','off');
        drawnow;
    end
    
    part_name = local_part_name(part,mdl);
    part_name_tex = strrep(part_name,'_','\_');
    if any(strcmpi(metric_tag,{'NSE','KGE','JKGE'}))
        metric_title = sprintf('$\\mathrm{%s}$',metric_tag);
    elseif strcmpi(part,'sage_multi')
        model_name_fig = 'SAGE';
    else
        metric_title = sprintf('$%s$',metric_tag);
    end
    
    figureTitleSize = 21;
    % Keep the title clear of the upper export boundary; LaTeX ascenders
    % can otherwise be clipped in both the GUI and the PPTX rasterization.
    annotation(fig,'textbox',[0.02 0.935 0.96 0.055], ...
        'String',sprintf(['$\\texttt{%s}$: ' ...
        '%s %s performance summary for basin/period scenarios'], ...
        part_name_tex,dt_str,metric_title), ...
        'Interpreter','latex', ...
        'FontSize',figureTitleSize, ...
        'EdgeColor','none', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');

end

function plot_metric_ecdf_by_zone_figure( ...
    figNo,part,mdl,dt_str,samp_word, ...
    scenarios,n_m,id_model,model_names, ...
    colors,figW,figH,model_name_fig, ...
    metric_tag,~,bas)
%PLOT_METRIC_ECDF_BY_ZONE_FIGURE Plots metric ECDFs by hydroclimatic zone.

    if ~isfield(bas,'zone') ...
            || ~isfield(bas.zone,'num') ...
            || isempty(bas.zone.num)
        return
    end
    
    zoneNo = double( ...
        bas.zone.num(:));
    zoneIDs = unique(zoneNo( ...
        isfinite(zoneNo) ...
        & zoneNo > 0));
    
    if isempty(zoneIDs)
        return
    end
    
    zoneNames = ...
        local_zone_names_from_bas(bas);

    nScen = 4;
    zonesPerPage = 4;
    nPage = ceil(numel(zoneIDs)/zonesPerPage);
    [xL,xR] = local_metric_ecdf_limits(scenarios,metric_tag);

    for page = 1:nPage
        firstZone = (page-1)*zonesPerPage + 1;
        lastZone = min(page*zonesPerPage,numel(zoneIDs));
        pageZones = zoneIDs(firstZone:lastZone);
        nCol = numel(pageZones);

        f = figure(figNo+page-1); clf(f);
        pageWidth = min(0.95*figW,2.85*nCol+1.65);
        pageHeight = min(0.95*figH,10.8);
        set(f,'Name',sprintf( ...
            '%s: %s ECDF by hydroclimatic zone (%d/%d)', ...
            model_name_fig,metric_tag,page,nPage), ...
            'NumberTitle','off','Color','w','Units','inches', ...
            'Position',[0.2 0.2 pageWidth pageHeight]);

        tlo = tiledlayout(f,nScen,nCol, ...
            'TileSpacing','compact','Padding','compact');

        for izPage = 1:nCol
            zid = pageZones(izPage);
            for is = 1:nScen %#ok<ALIGN>
                axh = nexttile(tlo,(is-1)*nCol+izPage);
            hold(axh,'on');
            box(axh,'off');
    
            set(axh,'tickdir','out', ...
                'ticklength',[0.015 0.015], ...
                'fontsize',13, ...
                'linewidth',1);
    
            xlim(axh,[xL xR]);
            ylim(axh,[0 1]);
            axis(axh,'square');
            local_draw_top_right_frame_axis_color(axh);
    
            % Two-line titles prevent neighboring zone names colliding.
            if is == 1
                if zid <= numel(zoneNames)
                    title(axh,{sprintf('Zone %d',zid), ...
                        char(string(zoneNames(zid)))}, ...
                        'Interpreter','none','FontSize',12);
                else
                    title(axh,sprintf('Zone %d',zid), ...
                        'Interpreter','none','FontSize',12);
                end
            end
            
            % Row labels: scenarios, split across two lines
            if izPage == 1
                ytxt = local_scenario_ylabel(is,samp_word);
                ylabel(axh,ytxt, ...
                    'interpreter','latex', ...
                    'fontsize',13);
            else
                axh.YTickLabel = [];
            end
    
            if is == nScen
                xlabel(axh,['$' local_metric_symbol(metric_tag) '$'], ...
                    'interpreter','latex', ...
                    'fontsize',12);
            else
                axh.XTickLabel = [];
            end
    
            S = scenarios{is};
            M = S.metric;
    
            if ~local_has_data(M) ...
                    || isempty(S.id_global)
                text(axh,0,0.5,'Not available', ...
                    'horizontalalignment','center', ...
                    'verticalalignment','middle');
                continue
            end
    
            idAll = S.id_global(:);
            validRows = idAll >= 1 ...
                & idAll <= numel(zoneNo);
    
            rowKeep = false(numel(idAll),1);
            tmp = zoneNo(idAll(validRows)) == zid;
            rowKeep(validRows) = tmp(:);
    
            if ~any(rowKeep)
                text(axh,0,0.5,'No basins', ...
                    'horizontalalignment','center', ...
                    'verticalalignment','middle');
                continue
            end
    
            Mz = M(rowKeep,:);
            nz = sum(rowKeep);
    
            % Interbasin score for this zone/scenario/model:
            % S_IB = 1/nz * sum_i (1 - metric_i)
            SIBz = nan(1,n_m);
            for mdl_i = 1:n_m
                v = Mz(:,mdl_i);
                v = v(isfinite(v));
                if ~isempty(v)
                    SIBz(mdl_i) = mean(1 - v);
                end
            end
    
            for mdl_i = 1:n_m
                c = colors(id_model( ...
                    min(mdl_i,numel(id_model))),:);
                im = id_model( ...
                    min(mdl_i,numel(id_model)));
                nam = ...
                    local_model_display_name( ...
                    model_names{im});
                lbl = ['\texttt{' local_latex_escape(nam) '}'];
    
                plot_ecdf_panel(axh, ...
                    Mz(:,mdl_i),c,lbl, ...
                    part,metric_tag,S.tag, ...
                    xL,xR,0.12, ...
                    1.35,0.8,4,10,true);
            end
    
            % n in top-left
            text(axh,0.02,0.95, ...
                sprintf('$n = %d$',nz), ...
                'Units','normalized', ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','top', ...
                'interpreter','latex', ...
                'fontsize',13);
            
            SIBtxt = local_sib_text(SIBz);
            % S_IB in bottom-right
            if ~isempty(SIBtxt)
                text(axh,0.98,0.05, ...
                    SIBtxt, ...
                    'Units','normalized', ...
                    'HorizontalAlignment','right', ...
                    'VerticalAlignment','bottom', ...
                    'interpreter','latex', ...
                    'fontsize',13);
            end
    
            if n_m > 1 && izPage == nCol && is == 1
                legend(axh,'show', ...
                    'interpreter','latex', ...
                    'fontsize',10, ...
                    'location','northeast', ...
                    'box','off');
            end
        end
        end

        part_name = local_part_name(part,mdl);
        part_name_tex = strrep(part_name,'_','\_');
        if any(strcmpi(metric_tag,{'NSE','KGE','JKGE'}))
            metric_title = sprintf('$\\mathrm{%s}$',metric_tag);
        else
            metric_title = sprintf('$%s$',metric_tag);
        end
        titleLine1 = sprintf('$\\texttt{%s}$: %s %s ECDF', ...
            part_name_tex,dt_str,metric_title);
        titleLine2 = sprintf(['hydroclimatic zone and basin/%s ' ...
            'scenario --- %d/%d'],samp_word,page,nPage);
        sgtitle(tlo,{titleLine1,titleLine2}, ...
            'Interpreter','latex','FontSize',15);
    end

end

function tf = local_has_data(A)
% LOCAL_HAS_DATA

    tf = ~isempty(A);
    
    if tf
        try
            tf = any(isfinite(A(:)));
        catch
            tf = true;
        end
    end
end

function zoneTbl = local_build_zone_summary_table( ...
    bas,scenariosNSE,scenariosKGE,scenariosJKGE, ...
    n_m,id_model,model_names)
%LOCAL_BUILD_ZONE_SUMMARY_TABLE Summarize metrics by zone and scenario

    varNames = {'Zone','ZoneName','Scenario','ScenarioName', ...
        'Model','n','T_NSE','T_KGE','T_JKGE','S_IB'};
    zoneTbl = table('Size',[0 numel(varNames)], ...
        'VariableTypes',{'double','string','string','string', ...
        'string','double','double','double','double','double'}, ...
        'VariableNames',varNames);
    
    if ~isfield(bas,'zone') ...
            || ~isfield(bas.zone,'num') ...
            || isempty(bas.zone.num)
        return
    end
    
    zoneNo = double(bas.zone.num(:));
    zoneIDs = unique(zoneNo(isfinite(zoneNo) ...
        & zoneNo > 0));
    if isempty(zoneIDs)
        return
    end
    zoneIDs = (1:max(zoneIDs)).';
    zoneNames = local_zone_names_from_bas(bas);
    
    nScen = 4;
    rows = cell(0,numel(varNames));
    
    for iz = 1:numel(zoneIDs)
        zid = zoneIDs(iz);
        if iz <= numel(zoneNames)
            zname = string(zoneNames(iz));
        else
            zname = "Z" + string(zid);
        end
    
        for is = 1:nScen
            SN = scenariosNSE{is};
            SK = scenariosKGE{is};
            SJ = scenariosJKGE{is};
    
            idAll = SN.id_global(:);
            validRows = idAll >= 1 ...
                & idAll <= numel(zoneNo);
            rowKeep = false(numel(idAll),1);
            tmp = zoneNo(idAll(validRows)) == zid;
            rowKeep(validRows) = tmp(:);
    
            if ~any(rowKeep)
                continue
            end
    
            for mdl_i = 1:n_m
                im = id_model(min(mdl_i, ...
                    numel(id_model)));
                try
                    modelName = string( ...
                        local_model_display_name( ...
                        model_names{im}));
                catch
                    modelName = "model_" + ...
                        string(im);
                end
    
                vNSE = local_metric_col( ...
                    SN.metric,rowKeep,mdl_i);
                vKGE = local_metric_col( ...
                    SK.metric,rowKeep,mdl_i);
                vJKGE = local_metric_col( ...
                    SJ.metric,rowKeep,mdl_i);
    
                nVal = max([nnz(isfinite(vNSE)), ...
                    nnz(isfinite(vKGE)), ...
                    nnz(isfinite(vJKGE))]);
    
                % Interbasin score based on NSE:
                % S_IB = 1/n * sum_i (1 - NSE_i)
                v = vNSE(isfinite(vNSE));
                if isempty(v)
                    SIB = NaN;
                else
                    SIB = mean(1 - v);
                end
    
                rows(end+1,:) = {double(zid),zname, ...         
                    string(SN.tag), ...
                    string(SN.title),modelName, ...
                    double(nVal),median(vNSE,'omitnan'), ...
                    median(vKGE,'omitnan'), ...
                    median(vJKGE,'omitnan'),SIB};    %#ok
            end
        end
    end
    
    if isempty(rows)
        return
    end
    
    zoneTbl = cell2table(rows,'VariableNames',varNames);
end

function v = local_metric_col(M,rowKeep,mdl_i)
%LOCAL_METRIC_COL Return finite metric column for one model

    v = nan(0,1);
    if isempty(M)
        return
    end
    try
        if size(M,2) < mdl_i
            return
        end
        v = M(rowKeep,mdl_i);
        v = v(:);
    catch
        v = nan(0,1);
    end
end


function A = local_get_metric(S, ...
    fieldName,defaultVal)
%LOCAL_GET_METRIC Extracts one metric/scenario array with empty protection.

    A = defaultVal;
    
    if isstruct(S) ...
            && isfield(S,fieldName)
        A = S.(fieldName);
    end
end

function n_m = local_infer_nmodels(Q,NSE,KGE)
%LOCAL_INFER_NMODELS Infers the number of models represented in plot inputs.

    C = {Q,NSE,KGE};
    
    % Prefer 3D arrays: [time x basin x model]
    for i = 1:numel(C)
        X = C{i};
    
        if isstruct(X)
            fn = fieldnames(X);
            for j = 1:numel(fn)
                A = X.(fn{j});
                if ~isempty(A) ...
                        && ndims(A) == 3
                    n_m = size(A,3);
                    return
                end
            end
        end
    end
    
    % Then 2D metric matrices: [basin x model]
    for i = 1:numel(C)
        X = C{i};
    
        if isstruct(X)
            fn = fieldnames(X);
            for j = 1:numel(fn)
                A = X.(fn{j});
                if ~isempty(A) ...
                        && ismatrix(A) ...
                        && size(A,1) > 1 ...
                        && size(A,2) <= 20
                    n_m = size(A,2);
                    return
                end
            end
        end
    end
    
    n_m = 1;
end

function plot_ecdf_panel(axh,z,c, ...
    lbl,part,metric_tag, ...
    scenario_tag, ...
    xL,xR,fill_alpha,line_width, ...
    med_lw,med_ms,fnt_med,compactLabel)
%PLOT_ECDF_PANEL Draws and styles one empirical CDF panel.

    if nargin < 15
        compactLabel = false;
    end

    z = z(:);
    z = z(isfinite(z));
    % For plotting only: prevent extreme JKGE values from creating
    % pathological ECDF polygons or off-axis median geometry.
    z_raw = z;
    z = max(xL,min(xR,z));
    
    if isempty(z)
        return
    end
    
    [f,x] = sage_ecdf(z);
    [xs,fs] = ecdf_to_stairs_fixed(x,f,xL,xR);
    
    % -------------------------------------
    % SITE: line only
    % SAGE: line + fill + median annotation
    % -------------------------------------
    if strcmpi(part,'sage')
        [xp,yp] = stairs_fill_poly(xs,fs);
        p = patch(axh,xp,yp,c, ...
            'facealpha',fill_alpha, ...
            'edgealpha',0);
        set(p,'handlevisibility','off');
    end
    
    hLine = plot(axh,xs,fs, ...
        'color',c, ...
        'linewidth',line_width);
    set(hLine,'displayname',lbl);
    
    if strcmpi(part,'site')
        return
    end
    
    %medX = median(z);
    medX_raw = median(z_raw);
    medX = max(xL,min(xR,medX_raw));
    
    medF = ecdf_value_from_stairs(xs,fs,medX);
    
    [xMark,yMark,x1m,x2m,xt,ha,mode] = ...
        median_marker_geom(medX,medF,xL,xR);
    
    if ~mode.specialHalfAtLeft
        line(axh,[xMark xMark], ...
            [0 yMark], ...
            'color','k', ...
            'linewidth',med_lw, ...
            'handlevisibility','off');
    end
    
    line(axh,[x1m x2m], ...
        [yMark yMark], ...
        'color','k', ...
        'linewidth',med_lw, ...
        'handlevisibility','off');
    
    line(axh,xMark,yMark, ...
        'marker','s', ...
        'markerfacecolor','k', ...
        'markeredgecolor','k', ...
        'linestyle','none', ...
        'markersize',med_ms, ...
        'handlevisibility','off');
    
    % txt = sprintf(['$\\widehat{T}_' ...
    %     '{\\mathrm{%s}_{\\rm %s}} = %.3f$'], ...
    %     metric_tag,scenario_tag,medX);
    if medX_raw < xL
        medTxt = sprintf('< %.1f',xL);
    elseif medX_raw > xR
        medTxt = sprintf('> %.1f',xR);
    else
        medTxt = sprintf('%.3f',medX_raw);
    end
    
    if compactLabel
        txt = sprintf('$%s = %s$', ...
            local_metric_median_symbol(metric_tag),medTxt);
        text(axh,xt,yMark,txt, ...
            'Interpreter','latex', ...
            'FontWeight','bold', ...
            'HorizontalAlignment',ha, ...
            'VerticalAlignment','middle', ...
            'Color','k','FontSize',fnt_med, ...
            'HandleVisibility','off');
    else
        txt = sprintf('$%s_{\\rm %s} = %s$', ...
            local_metric_median_symbol(metric_tag),scenario_tag,medTxt);
        text(axh,xt,yMark,txt, ...
            'Interpreter','latex', ...
            'FontWeight','bold', ...
            'HorizontalAlignment',ha, ...
            'VerticalAlignment','middle', ...
            'Color','k','FontSize',fnt_med, ...
            'HandleVisibility','off');
    end
end

function scenarios = local_build_metric_scenarios( ...
    M,Q,dt_str,K_t,K_e,id_tr,id_ev,gaugescen,samp_word)
%LOCAL_BUILD_METRIC_SCENARIOS Builds metric scenarios.

    scenarios = {};
    
    % ---------------------------------
    % training basins | training period
    % ---------------------------------
    S = struct();
    S.tag = 'tt';
    S.title = sprintf(['%s training basin ' ...
        '| training %s'],dt_str,samp_word);
    S.metric = local_get_metric(M,'tt',[]);
    S.Q = local_get_metric(Q,'tt',[]);
    S.id_global = 1:K_t;
    S.idx = id_tr;
    S.codes = local_get_codes(gaugescen,'tt');
    S.basin_idx = local_get_idx(gaugescen,'tt');
    scenarios{end+1} = S;
    
    % -----------------------------------
    % training basins | evaluation period
    % -----------------------------------
    S = struct();
    S.tag = 'te';
    S.title = sprintf(['%s training basin ' ...
        '| evaluation %s'],dt_str,samp_word);
    S.metric = local_get_metric(M,'te',[]);
    S.Q = local_get_metric(Q,'te',[]);
    S.id_global = 1:K_t;
    S.idx = id_ev;
    S.codes = local_get_codes(gaugescen,'te');
    S.basin_idx = local_get_idx(gaugescen,'te');
    scenarios{end+1} = S;
    
    % -----------------------------------
    % evaluation basins | training period
    % -----------------------------------
    S = struct();
    S.tag = 'et';
    S.title = sprintf(['%s evaluation basin ' ...
        '| training %s'],dt_str,samp_word);
    S.metric = local_get_metric(M,'et',[]);
    S.Q = local_get_metric(Q,'et',[]);
    S.id_global = (K_t+1):(K_t+K_e);
    S.idx = id_tr;
    S.codes = local_get_codes(gaugescen,'et');
    S.basin_idx = local_get_idx(gaugescen,'et');
    scenarios{end+1} = S;
    
    % -------------------------------------
    % evaluation basins | evaluation period
    % -------------------------------------
    S = struct();
    S.tag = 'ee';
    S.title = sprintf(['%s evaluation basin ' ...
        '| evaluation %s'],dt_str,samp_word);
    S.metric = local_get_metric(M,'ee',[]);
    S.Q = local_get_metric(Q,'ee',[]);
    S.id_global = (K_t+1):(K_t+K_e);
    S.idx = id_ev;
    S.codes = local_get_codes(gaugescen,'ee');
    S.basin_idx = local_get_idx(gaugescen,'ee');
    scenarios{end+1} = S;

end


function pick = local_pick_qfdc_ids(Qfdc,id_global,tag)
%LOCAL_PICK_QFDC_IDS Return FDC-selected global dat{k} basin indices.
%
% Qfdc.id is interpreted as global dat{k} indices.  If Qfdc.req.(tag)
% exists, it is interpreted as indices into Qfdc.id for that scenario
% (tt/te/et/ee), matching the FDC postprocessor convention.  The result is
% filtered to the supplied id_global group and returned in Qfdc order.

    pick = [];
    
    if nargin < 1 ...
            || isempty(Qfdc) ...
            || nargin < 2 ...
            || isempty(id_global)
        return
    end
    
    if nargin < 3
        tag = '';
    end
    
    if ~isstruct(Qfdc) ...
            || ~isfield(Qfdc,'id') ...
            || isempty(Qfdc.id)
        return
    end
    
    idsAll = local_numeric_vector(Qfdc.id);
    
    if isempty(idsAll)
        return
    end
    
    ids = idsAll;
    
    % Scenario-specific FDC request. Qfdc.req.(tag) usually stores positions
    % into Qfdc.id, not the global dat{k} IDs themselves.
    try
        if ~isempty(tag) ...
                && isfield(Qfdc,'req') ...
                && isstruct(Qfdc.req) ...
                && isfield(Qfdc.req,tag) ...
                && ~isempty(Qfdc.req.(tag))
            req = local_numeric_vector(Qfdc.req.(tag));
            req = round(req(isfinite(req)));
            req = req(req >= 1 ...
                & req <= numel(idsAll));
            if ~isempty(req)
                ids = idsAll(req);
            end
        end
    catch
        ids = idsAll;
    end
    
    ids = ids(:);
    ids = ids(isfinite(ids));
    ids = round(ids);
    
    id_global = double(id_global(:));
    pick = intersect(ids,id_global,'stable');

end

function v = local_numeric_vector(x)
%LOCAL_NUMERIC_VECTOR Convert numeric/cell/string vectors to double column.

    if isempty(x)
        v = [];
        return
    end
    
    if isnumeric(x) ...
            || islogical(x)
        v = double(x(:));
        return
    end
    
    if iscell(x)
        try
            x = string(x(:));
        catch
            v = [];
            return
        end
    end
    
    if ischar(x)
        x = string({x});
    end
    
    if isstring(x)
        v = str2double(x(:));
    else
        v = [];
    end

end


function pick = pick_basins(codes, ...
    id_dat,id_global,Bpick,region)
%PICK_BASINS Resolves requested gauge codes to basin indices.

    K = numel(id_global);
    Bpick = min(Bpick,K);
    
    if isempty(codes)
        ii = randperm(K,Bpick);
        pick = id_global(ii);
        return
    end
    
    if isempty(id_dat)
        ii = randperm(K,Bpick);
        pick = id_global(ii);
        return
    end
    
    if isstring(id_dat)
        id_dat = cellstr(id_dat);
    elseif ischar(id_dat)
        id_dat = cellstr(string(id_dat));
    end
    
    if isnumeric(codes)
        codes = codes(:);
        codes = codes(isfinite(codes));
    
        if isempty(codes)
            ii = randperm(K,Bpick);
            pick = id_global(ii);
            return
        end
    
        codes = arrayfun(@(c) ...
            local_basin_code(c,region), ...
            codes,'UniformOutput',false);
    end
    
    if isstring(codes)
        codes = cellstr(codes);
    end
    if ischar(codes)
        codes = {codes};
    end
    if ~iscell(codes)
        codes = {codes};
    end
    
    for i = 1:numel(codes)
        c = codes{i};
    
        if isempty(c)
            codes{i} = '';
            continue
        end
    
        if isnumeric(c) ...
                && isscalar(c) ...
                && ~isnan(c)
            s = local_basin_code(c,region);
        elseif isstring(c)
            s = char(c);
        elseif ischar(c)
            s = c;
        else
            try
                s = char(string(c));
            catch
                s = '';
            end
        end
    
        s = strtrim(s);
        s = local_basin_code(s,region);
    
        codes{i} = s;
    end
    
    codes = codes(~cellfun(@isempty,codes));
    codes = unique(codes,'stable');
    
    if isempty(codes)
        ii = randperm(K,Bpick);
        pick = id_global(ii);
        return
    end
    
    id_dat = cellfun(@(x) ...
        local_basin_code(x,region), ...
        id_dat,'UniformOutput',false);
    
    id_sub = id_dat(id_global);
    
    id_sub = cellfun(@(x) ...
        local_basin_code(x,region), ...
        id_sub,'UniformOutput',false);

    % ---------------------------------------------------------
    % Robust matching: remove leading zeros for numeric IDs.
    % This fixes regions such as CAMELS-GB where requested codes
    % may be "00041004" but available basin IDs are "41004".
    % ---------------------------------------------------------
    codes_match = ...
        local_normalize_basin_code_for_match(codes);
    id_sub_match = ...
        local_normalize_basin_code_for_match(id_sub);
    
    pick = [];
    for i = 1:numel(codes_match)
        id = codes_match{i};
        j = find(strcmp(id_sub_match,id),1,'first');
        if ~isempty(j)
            pick(end+1,1) = id_global(j); %#ok
        end
    end
    
    pick = unique(pick,'stable');
    
    if isempty(pick)
        nshow = min(10,numel(id_sub));
        ex = strjoin(id_sub(1:nshow),', ');

        warning('plot_SAGE:pick_basins:fallback', ...
            ['plot_SAGE: pick_basins: ' ...
            'none of the requested ', ...
            'gauge codes were found. ' ...
            'Requested: %s. ' ...
            'Example subset: %s. Using ' ...
            'available basins instead.'], ...
            strjoin(codes,', '),ex);

        pick = id_global(:);
        pick = pick(1:min(numel(pick),Bpick));
    else
        pick = pick(1:min(numel(pick),Bpick));
    end
end

function out = local_normalize_basin_code_for_match(x)
%LOCAL_NORMALIZE_BASIN_CODE_FOR_MATCH Normalize basin IDs for robust matching.

    if ischar(x)
        x = {x};
    elseif isstring(x)
        x = cellstr(x);
    elseif isnumeric(x)
        x = cellstr(string(x(:)));
    end

    if ~iscell(x)
        x = {x};
    end

    out = cell(size(x));

    for i = 1:numel(x)
        s = string(x{i});
        s = strtrim(s);
        s = regexprep(s, ...
            '\.0+$','');
        s = regexprep(s, ...
            '[^0-9A-Za-z]','');

        % If purely numeric, remove leading zeros.
        if ~isempty(regexp(char(s), ...
                '^\d+$','once'))
            s = regexprep(s,'^0+','');
            if strlength(s) == 0
                s = "0";
            end
        end

        out{i} = char(s);
    end
end

function R = local_resolve_requested_basins(codes,id_dat,id_global,region)
%LOCAL_RESOLVE_REQUESTED_BASINS Convert requested gauge codes to plot indices.

    R = struct('req',{{}}, ...
        'kdat',[],'icol',[], ...
        'missing',{{}}, ...
        'outside',{{}});

    if isempty(codes) ...
            || isempty(id_dat) ...
            || isempty(id_global)
        return
    end

    if isstring(id_dat)
        id_dat = cellstr(id_dat);
    elseif ischar(id_dat)
        id_dat = cellstr(string(id_dat));
    end

    if isstring(codes)
        codes = cellstr(codes);
    elseif ischar(codes)
        codes = {codes};
    elseif isnumeric(codes)
        codes = num2cell(codes(:));
    elseif ~iscell(codes)
        codes = cellstr(string(codes));
    end

    id_global = id_global(:);

    id_all = cellfun(@(x) ...
        local_basin_code(x,region), ...
        id_dat(:),'UniformOutput',false);

    id_sub = id_all(id_global);

    id_all_match = ...
        local_normalize_basin_code_for_match(id_all);
    id_sub_match = ...
        local_normalize_basin_code_for_match(id_sub);

    for i = 1:numel(codes)

        req = local_basin_code(codes{i},region);
        if isempty(req)
            continue
        end

        req_match = ...
            local_normalize_basin_code_for_match({req});
        req_match = req_match{1};

        % Exact match first
        j = find(strcmp(id_sub,req),1,'first');
        
        % Fallback: normalized match, fixes leading-zero numeric IDs
        if isempty(j)
            j = find(strcmp(id_sub_match,req_match),1,'first');
        end

        if ~isempty(j)
            R.req{end+1,1} = req;
            R.kdat(end+1,1) = id_global(j);
            R.icol(end+1,1) = j;
            continue
        end

        % Exact match first
        kall = find(strcmp(id_all,req),1,'first');

        % Fallback: normalized match
        if isempty(kall)
            kall = find(strcmp(id_all_match,req_match), ...
                1,'first');
        end

        if isempty(kall)
            R.missing{end+1,1} = req;
        else
            R.outside{end+1,1} = ...
                sprintf('%s {dat=%d}',req,kall);
        end
    end
end

function us = local_get_usgs(id_dat,k)
%LOCAL_GET_USGS Returns the preferred gauge identifier for one basin.

    us = 'gauge: ?';
    if isempty(id_dat)
        return
    end
    
    try
        if isstring(id_dat)
            id_dat = cellstr(id_dat);
        end
    
        if k < 1 ...
                || k > numel(id_dat)
            return
        end
    
        s = id_dat{k};
        if isstring(s)
            s = char(s); 
        end
        if iscell(s)  
            s = s{1};   
        end
    
        if isempty(s)
            return
        end
    
        us = ['gauge: ' s];
    catch
    end
end

function s = local_gauge_name(basins,k)
%LOCAL_GAUGE_NAME Returns a display name for one basin.

    s = "";
    
    try
        if ~isfield(basins,'gname') ...
                || isempty(basins.gname) ...
                || k < 1
            return
        end
    
        g = basins.gname;
    
        if istable(g)
            c = min(3,width(g));
            if k <= height(g) ...
                    && c >= 1
                v = g{k,c};
                if iscell(v)
                    v = v{1};
                end
                s = string(v);
            end
    
        elseif iscell(g)
            if k <= numel(g)
                s = string(g{k});
            end
    
        elseif isstring(g)
            if k <= numel(g)
                s = g(k);
            end
    
        elseif ischar(g)
            if size(g,1) >= k
                s = string(strtrim(g(k,:)));
            else
                s = string(strtrim(g));
            end
    
        elseif iscategorical(g)
            if k <= numel(g)
                s = string(g(k));
            end
    
        else
            if k <= numel(g)
                s = string(g(k));
            else
                s = string(g);
            end
        end
    catch
        s = "";
    end
    
    s = strtrim(s);
    
    if s == "" ...
            || strcmpi(s,"<missing>")
        s = "";
    end
end

function name = local_part_name(part,mdl)
%LOCAL_PART_NAME Returns a display label for the selected analysis mode.

    if strcmpi(part,'sage')
        if isempty(mdl.model)
            name = 'SAGE';
        else
            name = char(string( ...
                mdl.names(mdl.model(1))));
        end
    else
        name = 'SITE';
    end
end

function codes = local_get_codes(S,fieldName)
%LOCAL_GET_CODES Returns normalized requested basin codes from a structure.

    codes = [];
    
    if isstruct(S) && isfield(S,fieldName)
        codes = S.(fieldName);
        return
    end
    
    if ~isstruct(S)
        return
    end
    
    switch lower(fieldName)
        case {'tt','te'}
            if isfield(S,'train')
                codes = S.train;
            end
    
        case {'et','ee'}
            if isfield(S,'eval')
                codes = S.eval;
            end
    end
end

function idx = local_get_idx(S,fieldName)
%LOCAL_GET_IDX Return selected global dat{k} basin indices.

    idx = [];
    
    if ~isstruct(S)
        return
    end
    
    if isfield(S,'idxscen') ...
            && isstruct(S.idxscen) ...
            && isfield(S.idxscen,fieldName)
        idx = S.idxscen.(fieldName);
    elseif isfield(S,'idx') ...
            && isstruct(S.idx) ...
            && isfield(S.idx,fieldName)
        idx = S.idx.(fieldName);
    end
    
    idx = local_numeric_vector(idx);
    idx = round(idx(isfinite(idx)));

end

function q = local_get_discharge(Q,ii,mdl_i)
%LOCAL_GETQ

    if ndims(Q) == 3
        q = Q(:,ii,mdl_i);
    else
        q = Q(:,ii);
    end
end

function plot_parameter_traces_figure(figNo, ...
    mdl,tTheta,figW,figH,model_name_fig,pspace)
%PLOT_PARAMETER_TRACES_FIGURE Plot parameter trace summaries in one figure
%
% SYNOPSIS:
%   plot_parameter_traces_figure(figNo,mdl,tTheta,model_name_fig)
%   plot_parameter_traces_figure(figNo,mdl,tTheta,model_name_fig,pspace)
%
% INPUT:
%   figNo          figure number
%   mdl            model structure
%   tTheta           [i_max x d x 7] parameter trace percentiles
%                  third dimension = [5 15 25 50 75 85 95]
%   figW           figure width in inches [= from screen]
%   figH           figure height in inches [= from screen]
%   model_name_fig figure title string
%   pspace         1 = normalized parameter space (default)
%                  0 = hydrologic parameter space
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --------------
    % Input handling
    % --------------
    if nargin < 7 ...
            || isempty(pspace)
        pspace = 1;
    end
    
    if isempty(tTheta)
        return
    end
    
    if ndims(tTheta) ~= 3 ...
            || size(tTheta,3) ~= 7
        error('plot_parameter_traces_figure:BadSize', ...
            'tTheta must have size [i_max x d x 7].');
    end
    
    [i_max,d,~] = size(tTheta); %#ok
    
    % ---------------------
    % Keep valid iterations
    % ---------------------
    keep = squeeze(any(isfinite(tTheta(:,:,4)),2));
    if ~any(keep)
        return
    end
    
    tTheta = tTheta(keep,:,:);
    x = find(keep);
    
    % ------
    % Layout
    % ------
    ncol = 2;
    nrow = ceil(d/ncol);
    
    leftMargin = 0.12;
    rightMargin = 0.04;
    topMargin = 0.04;
    bottomMargin = 0.06;
    
    colGap = 0.1;    
    rowGap = 0.025;
    
    axW = (1 - leftMargin - rightMargin ...
        - colGap) / ncol;
    axH = (1 - topMargin - bottomMargin ...
        - (nrow-1)*rowGap) / nrow;
    
    figW = min(0.95*figW,8.5);
    figH = min(0.95*figH,1.2*nrow + 0.8);
    
    figTrace = figure(figNo); 
    clf(figTrace);
    set(figTrace, ...
        'Units','inches', ...
        'Position',[0.1 0.1 figW figH], ...
        'NumberTitle','off', ...
        'color','w', ...
        'Name',sprintf(['%s: ' ...
        'Parameter traces'], ...
        model_name_fig));
    
    ax = gobjects(d,1);
    ylh = gobjects(d,1);
    
    % -----
    % Style
    % -----
    fsAxis = 14;
    fsLabel = 14;
    fsPanel = 13;
    
    lwMed = 1.5;
    lwBound = 1.0;
    
    % band1 = [0.90 0.90 0.90];
    % band2 = [0.75 0.75 0.75];
    % band3 = [0.55 0.55 0.55];
    band1 = [0.85 0.92 0.98];   % light blue (5–95)
    band2 = [0.65 0.80 0.93];   % medium (15–85)
    band3 = [0.40 0.65 0.85];   % darker (25–75)
    
    % -------------------------
    % Common x-axis range/ticks
    % -------------------------
    xmax_all = x(end);
    xt_all = local_trace_xticks(xmax_all);
    
    % -----------
    % Loop panels
    % -----------
    for r = 1:nrow
        for c = 1:ncol
    
            j = (r-1)*ncol + c;
            if j > d
                continue
            end
    
            left = leftMargin + ...
                (c-1)*(axW + colGap);
            bottom = 1 - topMargin ...
                - r*axH - (r-1)*rowGap;
    
            ax(j) = axes('Position', ...
                [left bottom axW axH]);
            hold(ax(j),'on');
            box(ax(j),'off');
    
            Q = reshape(tTheta(:,j,:), ...
                size(tTheta,1),size(tTheta,3));
    
            % -------------------------------------
            % Transform to physical space if needed
            % -------------------------------------
            if pspace == 0
                [thmin,thmax] = ...
                    local_trace_param_bounds(mdl,j);
                if isfinite(thmin) ...
                        && isfinite(thmax)
                    Q = thmin + Q.*(thmax - thmin);
                end
            end
    
            q05 = Q(:,1);
            q15 = Q(:,2);
            q25 = Q(:,3);
            q50 = Q(:,4);
            q75 = Q(:,5);
            q85 = Q(:,6);
            q95 = Q(:,7);
    
            good = isfinite(q50);
            if ~any(good)
                continue
            end
            xx = x(good);
    
            % ----------------
            % Percentile bands
            % ----------------
            fill(ax(j),[xx; flipud(xx)], ...
                [q05(good); flipud(q95(good))], ...
                band1,'EdgeColor','none');
    
            fill(ax(j),[xx; flipud(xx)], ...
                [q15(good); flipud(q85(good))], ...
                band2,'EdgeColor','none');
    
            fill(ax(j),[xx; flipud(xx)], ...
                [q25(good); flipud(q75(good))], ...
                band3,'EdgeColor','none');
    
            % median
            plot(ax(j),xx, ...
                q50(good),'k', ...
                'linewidth',lwMed);
    
            xmax = xmax_all;
    
            % --------------
            % Y-axis scaling
            % --------------
            if pspace == 1
                plot(ax(j),[1 xmax], ...
                    [0 0],'r--', ...
                    'linewidth',lwBound);
                plot(ax(j),[1 xmax], ...
                    [1 1],'r--', ...
                    'linewidth',lwBound);
    
                ylim(ax(j),[-0.1 1.1]);
                yt = 0:0.2:1;
    
            else
                [thmin,thmax] = ...
                    local_trace_param_bounds(mdl,j);
    
                if isfinite(thmin) ...
                        && isfinite(thmax)
                    plot(ax(j),[1 xmax], ...
                        [thmin thmin],'r--', ...
                        'linewidth',lwBound);
                    plot(ax(j),[1 xmax], ...
                        [thmax thmax],'r--', ...
                        'linewidth',lwBound);
                end
    
                yall = [q05(good); q15(good); ...
                        q25(good); q50(good); ...
                        q75(good); q85(good); ...
                        q95(good); thmin; thmax];
                yall = yall(isfinite(yall));
    
                if isempty(yall)
                    ylim(ax(j),[-0.1 1.1]);
                else
                    ylo = min(yall);
                    yhi = max(yall);
    
                    if yhi <= ylo
                        dy = max(1e-6, ...
                            0.01*max(abs([ylo yhi])));
                        ylo = ylo - dy;
                        yhi = yhi + dy;
                    end
    
                    pad = 0.08*(yhi - ylo);
                    ylim(ax(j),[ylo-pad yhi+pad]);
                end
    
                yt = get(ax(j),'YTick');
            end
    
            xlim(ax(j),[1 xmax]);
    
            % ------------
            % Axis styling
            % ------------
            set(ax(j), ...
                'fontname','Times', ...
                'fontsize',fsAxis, ...
                'tickdir','out', ...
                'linewidth',1, ...
                'XTick',xt_all, ...
                'YTick',yt, ...
                'XMinorTick','on', ...
                'YMinorTick','off');
    
            ax(j).XRuler.TickLabelGapOffset = -2;
            ax(j).YRuler.TickLabelGapOffset = -1;
    
            % -------------------------------------
            % X label (bottom panel in each column)
            % -------------------------------------
            if j + ncol > d
                xlabel(ax(j),'Iteration, $i$', ...
                    'interpreter','latex', ...
                    'fontsize',fsLabel);
            else
                set(ax(j),'XTickLabel',[]);
            end
    
            % -------------------------------------------------------
            % Y-axis: show ticks + labels in BOTH columns (left side)
            % -------------------------------------------------------
            set(ax(j), ...
                'YAxisLocation','left', ...
                'YTickLabelMode','auto');
            
            ax(j).YRuler.TickLabelGapOffset = -1;
    
            % -------
            % Y label
            % -------
            ylh(j) = ylabel(ax(j), ...
                trace_ylabel_general( ...
                mdl,j,pspace), ...
                'interpreter','latex', ...
                'fontsize',fsLabel);
    
            % panel label
            local_panel_label(ax(j),j,fsPanel);
    
            % top/right frame
            local_draw_top_right_box(ax(j));
        end
    end
    
    % ------------
    % Align labels
    % ------------
    drawnow;
    
    leftIdx = 1:2:d;
    rightIdx = 2:2:d;
    
    if ~isempty(leftIdx)
        local_align_ylabels(ax(leftIdx));
    end
    if ~isempty(rightIdx)
        local_align_ylabels(ax(rightIdx));
    end
    
    drawnow;
    
    for j = 1:d
        if ~isgraphics(ylh(j))
            continue
        end
        if mod(j,2) == 0
            ylh(j).Position(1) = ...
                ylh(j).Position(1) + 0.12;
        else
            ylh(j).Position(1) = ...
                ylh(j).Position(1) - 0.04;
        end
    end
    
    % -----
    % Title
    % -----
    ttl = strrep(model_name_fig,'_','\_');
    annotation(figTrace,'textbox',[0 0.96 1 0.03], ...
        'String',sprintf(['$\\texttt{%s}$: ' ...
        'Parameter trace summary'],ttl), ...
        'interpreter','latex', ...
        'fontsize',16, ...
        'horizontalalignment','center', ...
        'verticalalignment','bottom', ...
        'EdgeColor','none');

end

function xt = local_trace_xticks(n)
%LOCAL_TRACE_XTICKS Clean iteration ticks starting at 1

    n = max(1,round(n));
    
    if n <= 10
        xt = 1:n;
    
    elseif n <= 14
        xt = [1 2:2:n];
    
    elseif n < 20
        xt = [1 5 10 15];
        xt = xt(xt <= n);
    
    elseif n < 50
        step = 5;
        xt = [1 step:step:n];
    
    elseif n < 100
        step = 10;
        xt = [1 step:step:n];
    
    else
        step = local_nice_trace_step(n);
        xt = [1 step:step:n];
    end
    
    xt = unique(xt);
end

function step = local_nice_trace_step(n)
%LOCAL_NICE_TRACE_STEP Choose readable trace tick spacing

    raw = n/5;
    pow10 = 10^floor(log10(raw));
    m = raw/pow10;
    
    if m <= 1
        nice = 1;
    elseif m <= 2
        nice = 2;
    elseif m <= 5
        nice = 5;
    else
        nice = 10;
    end
    
    step = nice*pow10;
    step = max(1,round(step));
end

function [thmin,thmax] = local_trace_param_bounds(mdl,j)
%LOCAL_TRACE_PARAM_BOUNDS Get physical parameter bounds if available

    thmin = NaN;
    thmax = NaN;
    
    try
        if isfield(mdl,'th_min') ...
                && numel(mdl.th_min) >= j
            thmin = double(mdl.th_min(j));
        end
        if isfield(mdl,'th_max') ...
                && numel(mdl.th_max) >= j
            thmax = double(mdl.th_max(j));
        end
    catch
    end
end

function local_panel_label(ax,j,fs)
%LOCAL_PANEL_LABEL Panel label: (a), (b), ..., (z), (aa), ...

    xl = xlim(ax);
    yl = ylim(ax);
    
    off = 0.015;
    x = xl(1) + off*(xl(2)-xl(1));
    y = yl(2) - off*(yl(2)-yl(1));
    
    txt = local_panel_label_string(j);
    
    text(ax,x,y,txt, ...
        'horizontalalignment','left', ...
        'verticalalignment','top', ...
        'fontname','Times', ...
        'fontsize',fs, ...
        'interpreter','none');
end

function s = local_panel_label_string(j)
%LOCAL_PANEL_LABEL_STRING Return (a), (b), ..., (z), (aa), ...

    letters = 'abcdefghijklmnopqrstuvwxyz';
    
    if j <= 26
        s = ['(' letters(j) ')'];
        return
    end
    
    q = j;
    out = '';
    while q > 0
        r = mod(q-1,26);
        out = [letters(r+1) out]; %#ok<AGROW>
        q = floor((q-1)/26);
    end
    s = ['(' out ')'];
end

function local_draw_top_right_frame_axis_color(ax)
%LOCAL_DRAW_TOP_RIGHT_FRAME_AXIS_COLOR Draw top/right lines without box ticks

    if ~isgraphics(ax)
        return
    end
    
    xl = xlim(ax);
    yl = ylim(ax);
    lw = ax.LineWidth;
    
    % Delete old frame lines, if the axes is refreshed
    if isappdata(ax,'TopRightFrameAxisColorHandles')
        hOld = getappdata(ax,'TopRightFrameAxisColorHandles');
        try
            delete(hOld(isgraphics(hOld)));
        catch
        end
    end
    
    h1 = line(ax,[xl(1) xl(2)],[yl(2) yl(2)], ...
        'color',ax.XColor, ...
        'linestyle','-', ...
        'linewidth',lw, ...
        'clipping','off', ...
        'handlevisibility','off');
    h2 = line(ax,[xl(2) xl(2)],[yl(1) yl(2)], ...
        'color',ax.YColor, ...
        'linestyle','-', ...
        'linewidth',lw, ...
        'clipping','off', ...
        'handlevisibility','off');
    
    setappdata(ax,'TopRightFrameAxisColorHandles',[h1 h2]);
end

function local_draw_top_right_box(ax)
%LOCAL_DRAW_TOP_RIGHT_BOX Draw top/right frame lines manually

    xl = xlim(ax);
    yl = ylim(ax);
    lw = ax.LineWidth;
    
    line(ax,[xl(1) xl(2)],[yl(2) yl(2)], ...
        'color','k', ...
        'linewidth',lw, ...
        'clipping','off');
    
    line(ax,[xl(2) xl(2)],[yl(1) yl(2)], ...
        'color','k', ...
        'linewidth',lw, ...
        'clipping','off');
end

function local_align_ylabels(axs)
%LOCAL_ALIGN_YLABELS Align y-labels for a vector of axes handles

    axs = axs(isgraphics(axs));
    if isempty(axs)
        return
    end
    
    drawnow;
    
    xpos = nan(numel(axs),1);
    for i = 1:numel(axs)
        yl = get(axs(i),'YLabel');
        if isgraphics(yl)
            pos = yl.Position;
            xpos(i) = pos(1);
        end
    end
    
    xpos = xpos(isfinite(xpos));
    if isempty(xpos)
        return
    end
    
    xmin = min(xpos);
    
    for i = 1:numel(axs)
        yl = get(axs(i),'YLabel');
        if isgraphics(yl)
            pos = yl.Position;
            pos(1) = xmin;
            yl.Position = pos;
        end
    end
end

function plot_parameter_maps_figure(figNo, ...
    mdl,bas,latlon,nTheta,region, ...
    model_name_fig,figW,figH,gaugescen)
%PLOT_PARAMETER_MAPS_FIGURE Plot final normalized parameter values on 
% regional map
%   latlon : K x 2 matrix with [lat lon]
%   nTheta    : d x K matrix of final normalized parameter values in [0,1]
%
% Uses same logic as SAGE_ui:
%   - geoaxes
%   - geobasemap('grayland')
%   - training basins: blue circles
%   - evaluation basins: orange squares
%   - lighter color = smaller normalized value
%   - darker color = larger normalized value

    if nargin < 10
        gaugescen = struct();
    elseif isempty(gaugescen)
        gaugescen = struct();
    end
    
    if isempty(latlon) ...
            || isempty(nTheta)
        return
    end
    
    if size(latlon,2) < 2
        warning('plot_parameter_maps_figure:badLatLon', ...
            ['      Warning: plot_parameter_maps_figure: ' ...
             'latlon must have size [K x 2]. ' ...
             'Skipping figure.']);
        return
    end
    
    if size(nTheta,2) ~= size(latlon,1)
        warning('plot_parameter_maps_figure:badSize', ...
            ['      Warning: plot_parameter_maps_figure: ' ...
             'size(nTheta,2) must equal size(latlon,1). ' ...
             'Skipping figure.']);
        return
    end
    
    if ~isfield(bas,'K_t') ...
            || ~isfield(bas,'K_e')
        warning('plot_parameter_maps_figure:missingBasinFields', ...
            ['      Warning: plot_parameter_maps_figure: ' ...
             'basins.K_t and/or basins.K_e missing. ' ...
             'Skipping figure.']);
        return
    end
    
    if ~isfield(bas,'id_gauge') ...
            || isempty(bas.id_gauge)
        warning('plot_parameter_maps_figure:missinggauge', ...
            ['      Warning: plot_parameter_maps_figure: ' ...
             'basins.id_gauge missing. Skipping figure.']);
        return
    end
    
    if bas.K_t + bas.K_e ~= size(latlon,1)
        warning('plot_parameter_maps_figure:badBasinCount', ...
            ['      Warning: plot_parameter_maps_figure: ' ...
             'bas.K_t + bas.K_e does not match size(latlon,1). ' ...
             'Skipping figure.']);
        return
    end
    
    if exist('geoaxes','file') ~= 2 ...
            || exist('geoscatter','file') ~= 2
        warning('plot_parameter_maps_figure:noGeoAxes', ...
            ['      Warning: plot_parameter_maps_figure: ' ...
             'geoaxes/geoscatter not available in this MATLAB ' ...
             'installation. Skipping figure.']);
        return
    end
    
    d = size(nTheta,1);
    
    if d < 1
        return
    end
    
    lat = latlon(:,1);
    lon = latlon(:,2);
    
    K_t = bas.K_t;
    K_e = bas.K_e;
    
    % [latLim,lonLim,mapName] = ...
    %     region_map_limits(region,lat,lon);
    [latLim,lonLim,mapName] = ...
        local_region_helpers_plot('maplimits', ...
        region,lat,lon);
    
    [nrow,ncol] = local_parameter_map_layout(d);
    
    % wider than tall, but not too wide
    figW = min(0.98*figW,16.3);                 % fixed width for all models
    figH = min(0.925*figH,3.0 + 2.45*nrow);     % calibrated from Xinanjiang template
    
    figure(figNo); clf;
    set(gcf,'Name',sprintf(['%s: ' ...
        '%s parameter maps'], ...
        model_name_fig,mapName), ...
        'NumberTitle','off', ...
        'color','w', ...
        'Units','inches', ...
        'Position',[0.1 0.1 figW figH]);
    
    left = 0.025;
    right = 0.020;
    top = 0.055;
    bottom = 0.075;

    hgap = 0.006;
    vgap = 0.012;
    
    panelW = (1 - left - right - ...
        (ncol-1)*hgap)/ncol;
    panelH = (1 - bottom - top - ...
        (nrow-1)*vgap)/nrow;
    
    blueBase = [0 0.4470 0.7410];
    orangeBase = [0.8500 0.3250 0.0980];
    
    % Line 2783-2784: change marker sizes
    msT = 20;   % was 65
    msV = 20;   % was 70  (won't be used if you hide eval basins)
    lwU = 0.40;
    lwS = 0.75;
    
    selTrain = false(K_t,1);
    selEval = false(K_e,1);
    
    try
        [selTrain,selEval] = ...
            local_get_parameter_map_selection( ...
            bas,gaugescen,region);
    catch
    end
    
    gax = gobjects(d,1);
    plot_eval = 0; % do not plot evaluation basins

    for j = 1:d
    
        row = ceil(j/ncol);
        col = mod(j-1,ncol) + 1;
    
        x0 = left + (col-1)*(panelW + hgap);
        y0 = 1 - top - row*panelH - (row-1)*vgap;
    
        gax(j) = geoaxes('Units','normalized', ...
            'Position',[x0 y0 panelW panelH]);
    
        hold(gax(j),'on');
        geobasemap(gax(j),'grayland');
        geolimits(gax(j),latLim,lonLim);
    
        % remove box / frame / graticule look 
        % as much as geoaxes allows
        try
            gax(j).Toolbar.Visible = 'off';
        catch
        end
        try
            gax(j).Grid = 'off';
        catch
        end
        try
            gax(j).MeridianLabel = 'off';
            gax(j).ParallelLabel = 'off';
        catch
        end
        try
            gax(j).LatitudeAxis.Visible = 'off';
            gax(j).LongitudeAxis.Visible = 'off';
        catch
        end
        try
            gax(j).Box = 'off';
        catch
        end
        try
            gax(j).FontSize = 10;
        catch
        end
    
        z = double(nTheta(j,:)).';
    
        % ---------------
        % training basins
        % ---------------
        if K_t > 0
            idxT = (1:K_t).';
            cT = zeros(K_t,3);
    
            for ii = 1:K_t
                cT(ii,:) = local_tint_color( ...
                    z(idxT(ii)),0,1,blueBase);
            end
    
            iu = find(~selTrain);
            is = find(selTrain);
    
            if ~isempty(iu)
                h = geoscatter(gax(j), ...
                    lat(idxT(iu)), ...
                    lon(idxT(iu)), ...
                    msT,'o','filled', ...
                    'MarkerEdgeColor', ...
                    [1 1 1], ...
                    'linewidth',lwU);
                h.CData = cT(iu,:);
            end
    
            if ~isempty(is)
                h = geoscatter(gax(j), ...
                    lat(idxT(is)), ...
                    lon(idxT(is)), ...
                    msT,'o','filled', ...
                    'MarkerEdgeColor', ...
                    [0 0 0], ...
                    'linewidth',lwS);
                h.CData = cT(is,:);
            end
        end
    
        % -----------------
        % evaluation basins
        % -----------------
        if K_e > 0 && (plot_eval == 1)
            idxE = (K_t+1:K_t+K_e).';
            cE = zeros(K_e,3);
    
            for kk = 1:K_e
                cE(kk,:) = local_tint_color( ...
                    z(idxE(kk)),0,1,orangeBase);
            end
    
            iu = find(~selEval);
            is = find(selEval);
    
            if ~isempty(iu)
                h = geoscatter(gax(j), ...
                    lat(idxE(iu)), ...
                    lon(idxE(iu)), ...
                    msV,'s','filled', ...
                    'MarkerEdgeColor', ...
                    [1 1 1], ...
                    'linewidth',lwU);
                h.CData = cE(iu,:);
            end
    
            if ~isempty(is)
                h = geoscatter(gax(j), ...
                    lat(idxE(is)), ...
                    lon(idxE(is)), ...
                    msV,'s','filled', ...
                    'MarkerEdgeColor', ...
                    [0 0 0], ...
                    'linewidth',lwS);
                h.CData = cE(is,:);
            end
        end
    
        % panel label slightly lower than before
        annotation('textbox', ...
            [x0 + 0.78*panelW, y0 + 0.9*panelH, ...
            0.08*panelW, 0.05*panelH], ...
        'String',local_panel_label_string(j), ...
        'interpreter','latex', ...
        'fontsize',20, ...
        'EdgeColor','none', ...
        'horizontalalignment','right', ...
        'verticalalignment','middle');
    
        % larger parameter symbol
        annotation('textbox', ...
            [x0 + 0.84*panelW, y0 + 0.88*panelH, ...
            0.12*panelW, 0.06*panelH], ...
        'String',trace_ylabel_general(mdl,j,1), ...
        'interpreter','latex', ...
        'fontsize',21, ...
        'EdgeColor','none', ...
        'horizontalalignment','left', ...
        'verticalalignment','middle');
    
        try
            if j == 1              % or use j == nPar for last panel
                gax(j).Scalebar.Visible = "on";
            else
                gax(j).Scalebar.Visible = "off";
            end
        catch
        end
    end
    
    ttl = local_latex_escape(model_name_fig);
    mapNameTex = local_latex_escape(mapName);
    
    annotation('textbox',[0 0.965 1 0.03], ...
        'String',sprintf(['\\texttt{%s}: ' ...
        '%s parameter maps'],ttl,mapNameTex), ...
        'interpreter','latex', ...
        'EdgeColor','none', ...
        'horizontalalignment','center', ...
        'verticalalignment','middle', ...
        'fontsize',18);
    
    nEmpty = nrow*ncol - d;
    if nEmpty >= 2
        local_add_parameter_colorbar_in_grid( ...
            figNo,blueBase, ...
            left,hgap,vgap, ...
            panelW,panelH,nrow,ncol,d);
    else
        local_add_parameter_colorbar( ...
            figNo,blueBase);
    end

end

function local_add_parameter_colorbar(figNo,blueBase)
%LOCAL_ADD_PARAMETER_COLORBAR Add horizontal normalized-parameter color bar

    figure(figNo);
    
    cbH = 0.018;
    cbW = 0.42;
    cbX = 0.5 - cbW/2;
    cbY = 0.055;

    n = 256;
    v = linspace(0,1,n);
    
    CB = zeros(n,3);
    for i = 1:n
        CB(i,:) = local_tint_color(v(i),0,1,blueBase);
    end
    
    ax1 = axes('Position',[cbX cbY cbW cbH]);
    hold(ax1,'on');
    
    for i = 1:n-1
        patch(ax1,[v(i) v(i+1) v(i+1) v(i)], ...
            [0 0 1 1],CB(i,:), ...
            'EdgeColor','none');
    end
    
    set(ax1, ...
        'XLim',[0 1], ...
        'YLim',[0 1], ...
        'YTick',[], ...
        'XTick',0:0.2:1, ...
        'XTickLabel',{'0','','','','','1'}, ...
        'tickdir','out', ...
        'XAxisLocation','bottom', ...
        'fontname','Calibri', ...
        'fontsize',17, ...
        'box','off', ...
        'Layer','top');
    
    % Manually draw a box
    set(ax1,'Box','off','TickDir','out','XAxisLocation','bottom');
    line(ax1,[0 1],[1 1],'Color','k','LineWidth',0.25);
    % End manually draw a box
    
    ax1.XRuler.TickLabelGapOffset = -3;
    hx1 = xlabel(ax1,'normalized parameter value', ...
        'fontname','Calibri', ...
        'fontsize',21);    
    hx1.Units = 'normalized';
    hx1.Position(2) = -0.75;

end

function plot_parameter_maps_country_figure(figNo, ...
    mdl,bas,latlon,nTheta,region, ...
    model_name_fig,figW,figH,gaugescen)
%PLOT_PARAMETER_MAPS_COUNTRY_FIGURE Parameter maps with country boundary only

    plot_parameter_maps_figure(figNo, ...
        mdl,bas,latlon,nTheta,region, ...
        model_name_fig,figW,figH, ...
        gaugescen);
    
    fig = figure(figNo);
    ax = findall(fig,'Type','geoaxes');
    
    [latLim,lonLim] = local_parameter_map_limits(latlon,region);

    for i = 1:numel(ax)
        try
            geobasemap(ax(i),'none');
        catch
        end
    
        try
            local_add_country_boundary(ax(i),region);
        catch ME
            warning('plot_SAGE:countryBoundaryFailed', ...
                ['      Warning: country boundary ' ...
                'could not be plotted.\n%s'], ...
                ME.message);
        end
        % Force the map back to the basin extent after plotting boundary
        try
            geolimits(ax(i),latLim,lonLim);
        catch
        end
    end
end

function local_add_country_boundary(axh,region)
%LOCAL_ADD_COUNTRY_BOUNDARY Add Natural Earth country outline to geoaxes

    persistent G
    persistent Gfile

    shp = local_find_natural_earth_countries();

    if isempty(shp)
        warning('plot_SAGE:noCountryBoundaryFile', ...
            ['No Natural Earth 10m or 50m ' ...
             'country shapefile was found.']);
        return
    end

    % Reload when the selected map file changes. This lets a newly
    % installed 10m map be used without restarting MATLAB.
    if isempty(G) ...
            || isempty(Gfile) ...
            || ~strcmpi(Gfile,shp)

        G = readgeotable(shp);
        Gfile = shp;
    end

    country = region_helpers('country_name',region);
    
    if country == ""
        return
    end
    
    vn = string(G.Properties.VariableNames);
    
    if any(strcmpi(vn,'ADMIN'))
        names = string(G.ADMIN);
    elseif any(strcmpi(vn,'NAME'))
        names = string(G.NAME);
    elseif any(strcmpi(vn,'SOVEREIGNT'))
        names = string(G.SOVEREIGNT);
    else
        warning('plot_SAGE:countryNameField', ...
            ['Could not identify country-name ' ...
            'field in Natural Earth table.']);
        return
    end
    
    row = strcmpi(names,country);
    
    if ~any(row)
        warning('plot_SAGE:countryNotFound', ...
            ['Country boundary not found ' ...
            'for region %s.'],string(region));
        return
    end
    
    hold(axh,'on');

    geoplot(axh,G.Shape(row), ...
        'FaceColor','none', ...
        'EdgeColor',[0 0 0], ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');
end

function shp = local_find_natural_earth_countries()
%LOCAL_FIND_NATURAL_EARTH_COUNTRIES Resolve export-map boundaries.
%
% An explicitly installed 10 m map is preferred for high-resolution
% postprocessing output. Downloading is never performed here. The bundled
% 50 m map keeps plot_SAGE fully functional offline.

    shp = '';

    pathNames = { ...
        'ne_10m_admin_0_countries.shp', ...
        'ne_50m_admin_0_countries.shp'};

    for i = 1:numel(pathNames)
        p = which(pathNames{i});

        if ~isempty(p) ...
                && isfile(p)
            shp = p;
            return
        end
    end

    relativeFiles = { ...
        fullfile('10', ...
            'ne_10m_admin_0_countries.shp'), ...
        fullfile('50', ...
            'ne_50m_admin_0_countries.shp'), ...
        'ne_50m_admin_0_countries.shp'};

    thisDir = fileparts( ...
        mfilename('fullpath'));

    mapRoots = { ...
        fullfile(prefdir,'SAGE','maps'), ...
        fullfile(thisDir,'maps'), ...
        fullfile(thisDir,'..','maps'), ...
        fullfile(thisDir,'..','..','maps')};

    if isdeployed
        mapRoots = [mapRoots, { ...
            fullfile(ctfroot,'maps'), ...
            fullfile(ctfroot, ...
                'SAGEhydrology','maps')}];
    end

    for r = 1:numel(mapRoots)
        for i = 1:numel(relativeFiles)
            candidate = fullfile( ...
                mapRoots{r},relativeFiles{i});

            if isfile(candidate)
                shp = candidate;
                return
            end
        end
    end
end

function [selTrain,selEval] = ...
    local_get_parameter_map_selection( ...
    basins,gaugescen,region)
%LOCAL_GET_PARAMETER_MAP_SELECTION Selects training/evaluation basins for maps.
    
    selTrain = false(basins.K_t,1);
    selEval = false(basins.K_e,1);
    
    if nargin < 2 ...
            || isempty(gaugescen) ...
            || ~isstruct(gaugescen)
        return
    end
    
    id_all = strings(numel(basins.id_gauge),1);
    for k = 1:numel(id_all)
        id_all(k) = local_basin_code( ...
            basins.id_gauge(k),region);
    end
    
    reqT = {};
    reqE = {};
    
    if isfield(gaugescen,'train') ...
            && ~isempty(gaugescen.train)
        reqT = [reqT,cellstr( ...
            string(gaugescen.train))];
    end
    if isfield(gaugescen,'tt') ...
            && ~isempty(gaugescen.tt)
        reqT = [reqT,cellstr( ...
            string(gaugescen.tt))];
    end
    if isfield(gaugescen,'te') ...
            && ~isempty(gaugescen.te)
        reqT = [reqT,cellstr( ...
            string(gaugescen.te))];
    end
    
    if isfield(gaugescen,'eval') ...
            && ~isempty(gaugescen.eval)
        reqE = [reqE,cellstr( ...
            string(gaugescen.eval))];
    end
    if isfield(gaugescen,'et') ...
            && ~isempty(gaugescen.et)
        reqE = [reqE,cellstr( ...
            string(gaugescen.et))];
    end
    if isfield(gaugescen,'ee') ...
            && ~isempty(gaugescen.ee)
        reqE = [reqE,cellstr( ...
            string(gaugescen.ee))];
    end
    
    reqT = unique(string(reqT(:)), ...
        'stable');
    reqE = unique(string(reqE(:)), ...
        'stable');
    
    for i = 1:numel(reqT)
        s = string(local_basin_code(reqT(i),region));
        jj = find(id_all(1:basins.K_t) == s, ...
            1,'first');
        if ~isempty(jj)
            selTrain(jj) = true;
        end
    end
    
    for i = 1:numel(reqE)
        s = string(local_basin_code(reqE(i),region));
        jj = find(id_all(basins.K_t+1: ...
            basins.K_t+basins.K_e) == s, ...
            1,'first');
        if ~isempty(jj)
            selEval(jj) = true;
        end
    end
end

function c = local_tint_color( ...
    v,vmin,vmax,baseColor)
%LOCAL_TINT_COLOR Blends a base color toward white by a specified fraction.

    if ~isfinite(v)
        c = [0.75 0.75 0.75];
        return
    end
    
    if ~isfinite(vmin) ...
            || ~isfinite(vmax) ...
            || vmax <= vmin
        t = 0.75;
    else
        t = (v - vmin) / (vmax - vmin);
        t = max(0,min(1,t));
        t = t.^0.7;
    end
    
    a = 0.20 + 0.80*t;
    c = (1-a)*[1 1 1] + a*baseColor;
end

function runs = local_contiguous_runs(id)
%LOCAL_CONTIGUOUS_RUNS Return contiguous runs as an N x 2 matrix

    id = double(id(:)).';   % force row vector
    
    if isempty(id)
        runs = zeros(0,2);
        return
    end
    
    ibreak = [1, find(diff(id) > 1) + 1, numel(id) + 1];
    nrun = numel(ibreak) - 1;
    
    runs = zeros(nrun,2);
    for k = 1:nrun
        i1 = ibreak(k);
        i2 = ibreak(k+1) - 1;
        runs(k,1) = id(i1);
        runs(k,2) = id(i2);
    end

end


function plot_parameter_attribution_figure(figNo, ...
    mdl,At,An,figW,figH,model_name_fig,viewMode)
%PLOT_PARAMETER_ATTRIBUTION_FIGURE Plot SAGE attribution like the GUI.
%
% At and An are expected as d x K x i_max arrays, initialized with NaNs.
% The latest iteration is inferred as the last third-dimension slice with
% at least one finite entry in At or An.  Each slice is reduced over basins
% using the same default statistic used by the GUI:
%
%     median(abs(A(:,:,i)),2,'omitnan')
%
% The resulting matrices are d x i, so the x-axis is iteration number,
% not basin number.

    if nargin < 6 ...
            || isempty(figH)
        figH = 8;
    end
    if nargin < 7 ...
            || isempty(model_name_fig)
        model_name_fig = 'SAGE';
    end
    if nargin < 8 ...
            || isempty(viewMode)
        viewMode = 'abs';
    end
    
    if isempty(At) ...
            || isempty(An)
        return
    end
    
    [AtPlot,AnPlot,last_i] = ...
        local_build_attribution_history(At,An);
    
    if last_i < 1 ...
            || isempty(AtPlot) ...
            || isempty(AnPlot)
        return
    end
    
    isRel = strcmpi(viewMode,'rel') ...
        || strcmpi(viewMode,'relative');
    
    if isRel
        AtPlot = local_attr_relative_rows(AtPlot);
        AnPlot = local_attr_relative_rows(AnPlot);
        viewTxt = 'relative';
        ylabTxt = 'Relative attribution';
    else
        viewTxt = 'absolute';
        ylabTxt = 'Attribution';
    end
    
    d = size(AtPlot,1);
    
    namesLatex = cell(d,1);
    for j = 1:d
        try
            namesLatex{j} = ...
                trace_ylabel_general(mdl,j,0);
        catch
            namesLatex{j} = ...
                sprintf('$\\theta_{%d}$',j);
        end
    end
    
    figW = 0.7*figW; figH = 0.9*figH;
    fig = figure(figNo);
    clf(fig);
    set(fig,'color','w', ...
        'Name',sprintf(['%s: ' ...
        'SAGE attribution (%s)'], ...
        model_name_fig,viewTxt), ...
        'NumberTitle','off');
    try
        set(fig,'Units','inches', ...
                'Position',[0.1 0.1 figW figH]);
    catch
    end
    
    tl = tiledlayout(fig,2,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');
    
    xt = local_attr_xticks_plot(last_i);
    
    % --------------------------
    % Time-wise attribution, A_t
    % --------------------------
    ax1 = nexttile(tl,1);
    imagesc(ax1,AtPlot);
    axis(ax1,'tight');
    ax1.YDir = 'normal';
    ax1.Toolbar.Visible = 'off';
    ax1.Box = 'on';
    ax1.FontSize = 13;
    ax1.YTick = 1:d;
    ax1.YTickLabel = namesLatex;
    ax1.TickLabelInterpreter = 'latex';
    ax1.XTick = xt;
    xlabel(ax1,'Iteration', ...
        'interpreter','latex');
    title(ax1, ...
        tern_local(strcmpi(viewTxt,'relative'), ...
        ['Time-wise attribution, ' ...
         '$A_{\rm t}$ (relative)'], ...
        ['Time-wise attribution, ' ...
         '$A_{\rm t}$ (absolute)']), ...
        'interpreter','latex');
    cb = colorbar(ax1,'southoutside');
    cb.Label.String = tern_local(isRel, ...
        'Relative attribution', ...
        'Absolute attribution');
    cb.Label.Interpreter = 'latex';
    
    % --------------------
    % Net attribution, A_n
    % --------------------
    ax2 = nexttile(tl,2);
    imagesc(ax2,AnPlot);
    axis(ax2,'tight');
    ax2.YDir = 'normal';
    ax2.Toolbar.Visible = 'off';
    ax2.Box = 'on';
    ax2.FontSize = 13;
    ax2.YTick = 1:d;
    ax2.YTickLabel = namesLatex;
    ax2.TickLabelInterpreter = 'latex';
    ax2.XTick = xt;
    xlabel(ax2,'Iteration', ...
        'interpreter','latex');
    title(ax2, ...
        tern_local(strcmpi(viewTxt,'relative'), ...
        ['Net attribution, ' ...
         '$A_{\rm n}$ (relative)'], ...
        ['Net attribution, ' ...
         '$A_{\rm n}$ (absolute)']), ...
        'interpreter','latex');
    cb = colorbar(ax2,'southoutside');
    cb.Label.String = tern_local(isRel, ...
        'Relative attribution', ...
        'Absolute attribution');
    cb.Label.Interpreter = 'latex';
    
    % ---------------------
    % Latest iteration bars
    % ---------------------
    ax3 = nexttile(tl,3);
    bar(ax3,AtPlot(:,last_i));
    grid(ax3,'on');
    ax3.Toolbar.Visible = 'off';
    ax3.Box = 'on';
    ax3.FontSize = 13;
    ax3.XTick = 1:d;
    ax3.XTickLabel = namesLatex;
    ax3.TickLabelInterpreter = 'latex';
    ax3.XTickLabelRotation = 45;
    ylabel(ax3,ylabTxt,'interpreter','latex');
    title(ax3,sprintf(['Latest iteration: ' ...
        '$A_{\\rm t}$ ($i = %d$)'],last_i), ...
        'interpreter','latex');
    
    ax4 = nexttile(tl,4);
    bar(ax4,AnPlot(:,last_i));
    grid(ax4,'on');
    ax4.Toolbar.Visible = 'off';
    ax4.Box = 'on';
    ax4.FontSize = 13;
    ax4.XTick = 1:d;
    ax4.XTickLabel = namesLatex;
    ax4.TickLabelInterpreter = 'latex';
    ax4.XTickLabelRotation = 45;
    ylabel(ax4,ylabTxt,'interpreter','latex');
    title(ax4,sprintf(['Latest iteration: ' ...
        '$A_{\\rm n}$ ($i = %d$)'],last_i), ...
        'interpreter','latex');
    
    try
        sgtitle(fig,sprintf(['\\texttt{%s}: ' ...
            'gradient attribution (%s view)'], ...
            strrep(model_name_fig,'_','\_'),viewTxt), ...
            'interpreter','latex', ...
            'fontweight','bold', ...
            'fontsize',18);
    catch
    end

end

function [AtHist,AnHist,last_i] = ...
    local_build_attribution_history(At,An)
%LOCAL_BUILD_ATTRIBUTION_HISTORY Convert d x K x i arrays to d x i.
%
% The latest iteration is the final third-dimension slice for which At or
% An contains at least one finite value.

    last_i = 0;
    AtHist = [];
    AnHist = [];
    
    if isempty(At) ...
            || isempty(An)
        return
    end
    
    szAt = size(At);
    szAn = size(An);
    
    if numel(szAt) < 3
        szAt(3) = 1;
    end
    if numel(szAn) < 3
        szAn(3) = 1;
    end
    
    d = szAt(1);
    K = szAt(2);
    iMax = min(szAt(3),szAn(3));
    
    if d < 1 ...
            || K < 1 ...
            || iMax < 1
        return
    end
    
    if szAn(1) ~= d ...
            || szAn(2) ~= K
        warning('plot_SAGE:attributionSizeMismatch', ...
            ['      Warning: plot_SAGE: At and An must have ' ...
             'the same d x K x i size. Skipping attribution.']);
        return
    end
    
    valid = false(1,iMax);
    for i = 1:iMax
        Ai = At(:,:,i);
        Ni = An(:,:,i);
        valid(i) = any(isfinite(Ai(:))) ...
            || any(isfinite(Ni(:)));
    end
    
    idx = find(valid,1,'last');
    if isempty(idx)
        return
    end
    
    last_i = idx;
    
    AtHist = NaN(d,last_i);
    AnHist = NaN(d,last_i);
    
    for i = 1:last_i
        AtHist(:,i) = local_attr_reduce_default(At(:,:,i));
        AnHist(:,i) = local_attr_reduce_default(An(:,:,i));
    end

end

function q = local_attr_reduce_default(A)
%LOCAL_ATTR_REDUCE_DEFAULT GUI default reduction across basins.
%
% Equivalent to attrReduce(A,'median(abs(.))') in SAGE_ui.

    if isempty(A)
        q = NaN(0,1);
        return
    end
    
    q = median(abs(A),2,'omitnan');

end

function Arel = local_attr_relative_rows(A)
%LOCAL_ATTR_RELATIVE_ROWS Match GUI row-wise relative scaling.

    Arel = abs(A);
    
    if isempty(Arel)
        return
    end
    
    rowmax = max(Arel,[],2,'omitnan');
    rowmax(~isfinite(rowmax) ...
        | rowmax <= 0) = 1;
    Arel = Arel ./ rowmax;

end

function xt = local_attr_xticks_plot(n)
%LOCAL_ATTR_XTICKS_PLOT Selects iteration ticks for attribution histories.

    if n <= 10
        xt = 1:n;
    elseif n <= 20
        xt = unique([1 5:5:n]);
    elseif n <= 50
        xt = unique([1 10:10:n]);
    elseif n <= 100
        xt = unique([1 20:20:n]);
    elseif n <= 500
        xt = unique([1 50:50:n]);
    else
        xt = unique([1 100:100:n]);
    end
    
    if isempty(xt)
        xt = 1;
    elseif xt(end) ~= n
        xt = unique([xt n]);
    end

end

function out = tern_local(tf,a,b)
%TERN_LOCAL Returns one of two values according to a logical condition.

    if tf
        out = a;
    else
        out = b;
    end

end

function [npairs,VG] = plot_parameter_variograms_figure(figNo, ...
    mdl,latlon,nTheta,model_name_fig,region)
%PLOT_PARAMETER_VARIOGRAMS_FIGURE Plot empirical spatial variograms
% of final parameter values and fit spherical variogram models.
%
% SYNOPSIS:
%   [npairs,VG] = plot_parameter_variograms_figure( ...
%       figNo,mdl,latlon,nTheta,model_name_fig)
%
% INPUT:
%   figNo           figure number
%   mdl             model structure
%    .par_names     1xd cell/string array of parameter names
%   latlon          K x 2 matrix with [lat lon] in degrees
%   nTheta             d x K matrix of parameter values
%   model_name_fig  string used in figure title
%
% OUTPUT:
%   npairs          1 x nBins vector with number of pairs per lag bin
%   VG              structure with variogram information
%    .binEdges      1 x (nBins+1) vector of lag-bin edges [km]
%    .binCenter     1 x nBins vector of lag-bin centers [km]
%    .gamma         d x nBins matrix of empirical semivariances
%    .npairs        1 x nBins vector of number of pairs
%    .fitPar        d x 3 matrix [nugget, partial sill, range]
%    .fitOK         d x 1 logical fit flag
%    .sill          d x 1 vector = nugget + partial sill
%    .label_plain   d x 1 cell with plain-text labels
%    .label_latex   d x 1 cell with LaTeX labels
%    .hmax_fit      scalar fit cutoff distance [km]
%    .xmax          scalar x-axis maximum [km]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt / adapted for SAGE UI                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6 ...
            || isempty(region)
        region = 'CAMELS_FR';
    end    
    try
        region = region_helpers('code',region);
    catch
        region = 'CAMELS_FR';
    end
    if nargin < 5 ...
            || isempty(model_name_fig)
        model_name_fig = 'SAGE';
    end
    
    VG = struct([]);
    
    if isempty(latlon) ...
            || isempty(nTheta)
        npairs = [];
        return
    end
    
    if size(latlon,2) < 2
        warning('plot_parameter_variograms_figure:badLatLon', ...
            ['      Warning: plot_parameter_variograms_figure: ' ...
             'latlon must have size [K x 2]. ' ...
             'Skipping figure.']);
        npairs = [];
        return
    end
    
    if size(nTheta,2) ~= size(latlon,1)
        warning('plot_parameter_variograms_figure:badSize', ...
            ['      Warning: plot_parameter_variograms_figure: ' ...
             'size(nTheta,2) must equal size(latlon,1). ' ...
             'Skipping figure.']);
        npairs = [];
        return
    end
    
    d = size(nTheta,1);
    K = size(nTheta,2);
    
    if d < 1 ...
            || K < 2
        npairs = [];
        return
    end
    
    % --------------------------------
    % User settings for variogram plot
    % --------------------------------
    nBins = 30;         % number of empirical lag bins
    hmax_fit = region_helpers('variogramrange',region);

    % ------
    % Colors
    % ------
    colors = hsv(d);
    
    % ----------------------------
    % Compute empirical variograms
    % ----------------------------
    Dkm = local_great_circle_distance_km( ...
        latlon(:,1),latlon(:,2));
    
    iu = triu(true(K),1);
    distVec = Dkm(iu);
    
    maxDistKm = max(distVec);
    if ~isfinite(maxDistKm) ...
            || maxDistKm <= 0
        warning('plot_parameter_variograms_figure:badDistance', ...
            ['      Warning: plot_parameter_variograms_figure: ' ...
             'Invalid basin distances. Skipping figure.']);
        npairs = [];
        return
    end
    
    binEdges = linspace(0,maxDistKm,nBins+1);
    binCenter = 0.5*(binEdges(1:end-1) + binEdges(2:end));
    binID = discretize(distVec,binEdges);
    
    gamma = nan(d,nBins);
    npairs = zeros(1,nBins);
    
    for b = 1:nBins
        mask = (binID == b);
        npairs(b) = sum(mask);
    
        if npairs(b) == 0
            continue
        end
    
        for j = 1:d
            z = double(nTheta(j,:)).';
            DZ = z - z.';
            dzVec = DZ(iu);
            gamma(j,b) = 0.5 * mean(dzVec(mask).^2,'omitnan');
        end
    end
    
    % ------------------------------------
    % Fit spherical model to each variogram
    % ------------------------------------
    fitPar = nan(d,3);   % [nugget, partial sill, range]
    fitOK = false(d,1);
    
    for j = 1:d
        gj = gamma(j,:);
    
        good = isfinite(gj) ...
            & isfinite(binCenter) ...
            & (npairs > 0) ...
            & (binCenter <= hmax_fit);
    
        h = binCenter(good);
        g = gj(good);
        w = npairs(good);
    
        if numel(h) < 4
            continue
        end
    
        % Initial guesses
        gmin = min(g);
        gmax = max(g);
        if ~isfinite(gmin)
            %gmin = 0; 
        end
        if ~isfinite(gmax)
            gmax = 1; 
        end
    
        nug0 = max(0,0.05*gmax);
        psill0 = max(1e-10,gmax - nug0);
        range0 = min(hmax_fit,0.5*max(h));
        if ~isfinite(range0) ...
                || range0 <= 0
            range0 = 1000;
        end
    
        p0 = [nug0, psill0, range0];
    
        objfun = @(p) local_variogram_objective( ...
            p,h,g,w,hmax_fit);
    
        opts = optimset('Display','off', ...
                        'MaxFunEvals',5000, ...
                        'MaxIter',5000, ...
                        'TolX',1e-8, ...
                        'TolFun',1e-8);
    
        try
            phat = fminsearch(objfun,p0,opts);
    
            phat(1) = max(0,phat(1));               % nugget
            phat(2) = max(0,phat(2));               % partial sill
            phat(3) = max(1,min(hmax_fit,phat(3))); % range
    
            if all(isfinite(phat))
                fitPar(j,:) = phat;
                fitOK(j) = true;
            end
        catch
            % leave as NaN/false
        end
    end

    % ------------------------
    % Prepare output structure
    % ------------------------
    label_plain = cell(d,1);
    label_latex = cell(d,1);
    for j = 1:d
        label_plain{j} = ...
            local_plain_label(mdl,j);
            % Normalized parameter label:
        % underline the main parameter symbol, but do not print subscripts.
        label_latex{j} = ...
            local_variogram_label_normalized_latex(mdl,j);
    end
    
    VG = struct();
    VG.binEdges = binEdges;
    VG.binCenter = binCenter;
    VG.gamma = gamma;
    VG.npairs = npairs;
    VG.fitPar = fitPar;
    VG.fitOK = fitOK;
    VG.sill = fitPar(:,1) + fitPar(:,2);
    VG.label_plain = label_plain;
    VG.label_latex = label_latex;
    VG.hmax_fit = hmax_fit;
    
    % ------
    % Figure
    % ------   
    figure(figNo); clf;
    set(gcf, ...
        'Name',sprintf('%s: Parameter variograms', ...
        model_name_fig), ...
        'NumberTitle','off', ...
        'Color','w', ...
        'Units','inches', ...
        'Position',[0.4 0.4 10.5 9.5]);
    
    % Use the regional fit cutoff also as the plotted x-axis cutoff.
    % Example: Germany hmax_fit = 600 km -> x-axis [0 600].
    xmaxPlot = hmax_fit;
    
    % Defensive fallback.
    if ~isfinite(xmaxPlot) ...
            || xmaxPlot <= 0
        xmaxPlot = maxDistKm;
    end
    xmaxPlot = max(1,xmaxPlot);
    
    VG.xmax = xmaxPlot;
    
    % Start with a uniformly spaced, region-independent tick set. It is
    % thinned below after the rendered subplot width is known. Never append
    % xmaxPlot to an existing sequence: that creates a short final interval
    % (for example 0, 1000, 1200) and overlapping endpoint labels.
    maxXTicks = 6;
    [xt,xtlbl] = local_variogram_xticks(xmaxPlot,maxXTicks);
    % ---- manual layout ----
    if d <= 12
        ncol = 3;
    else
        ncol = 4;
    end
    nrow = ceil(d/ncol);

    figW = 17.0;
    figH = 9.0;

    figure(figNo); clf;
    set(gcf,'Name',sprintf('%s: Parameter variograms', ...
        model_name_fig), ...
        'NumberTitle','off', ...
        'Color','w', ...
        'Units','inches', ...
        'Position',[0.6 0.6 figW figH]);

    fig = gcf;
    ttl = strrep(char(string(model_name_fig)),'_','\_');
    annotation(fig,'textbox',[0.00 0.955 1.00 0.035], ...
        'String',sprintf(['$\\texttt{%s}$: Variograms of final ' ...
        'normalized parameter values'],ttl), ...
        'Interpreter','latex', ...
        'FontSize',20, ...
        'FontWeight','bold', ...
        'EdgeColor','none', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');

    left = 0.080;
    right = 0.050;
    bottom = 0.1;

    top = 0.125;

    % Keep sufficient separation for neighboring y labels and for the
    % titles/x labels of vertically adjacent variogram panels.
    hgap = 0.085;
    vgap = 0.065;

    panelW = (1-left-right-(ncol-1)*hgap)/ncol;
    panelH = (1-bottom-top-(nrow-1)*vgap)/nrow;

    % A region-specific tick interval is useful only while its labels fit.
    % Thin the common tick set using the actual on-screen subplot width so
    % four-column figures remain readable on narrower laptop displays.
    maxXTicksFit = local_variogram_max_xticks( ...
        fig,panelW,15,xmaxPlot);
    if numel(xt) > maxXTicksFit
        [xt,xtlbl] = local_variogram_xticks( ...
            xmaxPlot,maxXTicksFit);
    end

    for j = 1:d

        row = ceil(j/ncol);
        col = mod(j-1,ncol) + 1;
        
        x0 = left + (col-1)*(panelW+hgap);
        y0 = 1 - top - row*panelH - (row-1)*vgap;
        
        ax = axes('Units','normalized', ...
            'Position',[x0 y0 panelW panelH]);
        hold(ax,'on');
        box(ax,'on');
    
        gj = gamma(j,:);
    
        YS = local_variogram_axis_scale( ...
            gj,fitOK(j),fitPar(j,:),hmax_fit, ...
            @local_spherical_variogram);
        
        scaleFac = YS.scaleFac;

        good = isfinite(gj) ...
            & isfinite(binCenter) ...
            & (npairs > 0) ...
            & (binCenter <= hmax_fit);
    
        % empirical variogram
        if any(good)
            plot(ax,binCenter(good),scaleFac*gj(good), ...
                's', ...
                'Color',colors(j,:), ...
                'MarkerFaceColor',colors(j,:), ...
                'MarkerEdgeColor',colors(j,:), ...
                'MarkerSize',5, ...
                'LineWidth',1.0);
        end
    
        % fitted spherical curve
        if fitOK(j)
            hfit = linspace(0,hmax_fit,400);
            gfit = local_spherical_variogram(hfit,fitPar(j,:));

            plot(ax,hfit,scaleFac*gfit, ...
                '-', ...
                'Color',colors(j,:), ...
                'LineWidth',2.0);
        end

        xlim(ax,[0 xmaxPlot]);
        xticks(ax,xt);
        
        set(ax, ...
            'TickDir','out', ...
            'FontSize',15, ...
            'LineWidth',1, ...
            'TickLabelInterpreter','latex');
        ax.XRuler.TickLabelGapOffset = -2;
        ax.XTickLabelRotation = 0;
        
        isBottomPanel = (j > d - ncol);

        if isBottomPanel
            xticklabels(ax,xtlbl);
            xlabel(ax,'Lag distance, $h$ (km)', ...
                'Interpreter','latex', ...
                'FontSize',16);
        else
            xticklabels(ax,repmat({''},size(xt)));
        end
    
        % only left column gets y-label
        % y-label for every panel
        ylab = YS.ylabelLatex;
        
        axpos = ax.Position;
        xlabpos = axpos(1) - 0.045;
        
        annotation('textbox', ...
            [xlabpos axpos(2)+0.4*axpos(4) 0.035 0.001*axpos(4)], ...
            'String',ylab, ...
            'Interpreter','latex', ...
            'FontSize',15, ...
            'EdgeColor','none', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'Rotation',90);

        % scaled y-limits
        ylim(ax,YS.ylim);
        yticks(ax,YS.yticks);
        yticklabels(ax,YS.yticklabels);

        panelStr = sprintf('%s %s', ...
            local_panel_label_string(j),label_latex{j});

        text(ax,0.01,1.03,panelStr, ...
            'Units','normalized', ...
            'Interpreter','latex', ...
            'FontSize',16, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top');
        box(ax,'off');
        set(ax,...
            'TickDir','out',...
            'TickLength',[0.03 0.015], ...
            'XAxisLocation','bottom',...
            'YAxisLocation','left',...
            'LineWidth',1.0);
    end
end

% ======================
% local helper functions
% ======================
function Dkm = local_great_circle_distance_km(latDeg,lonDeg)
%LOCAL_GREAT_CIRCLE_DISTANCE_KM Pairwise great-circle distances in km.

    R = 6371.0;
    
    lat = deg2rad(double(latDeg(:)));
    lon = deg2rad(double(lonDeg(:)));
    
    n = numel(lat);
    Dkm = zeros(n,n);
    
    for i = 1:n
        dlat = lat - lat(i);
        dlon = lon - lon(i);
    
        a = sin(dlat/2).^2 + cos(lat(i)) ...
            .* cos(lat) .* sin(dlon/2).^2;
        c = 2 * atan2(sqrt(a),sqrt(max(0,1-a)));
        Dkm(i,:) = R * c;
    end
end

function val = local_variogram_objective(p,h,g,w,hmax_fit)
%LOCAL_VARIOGRAM_OBJECTIVE Weighted least-squares objective for spherical fit.

    if any(~isfinite(p)) ...
            || numel(p) ~= 3
        val = 1e30;
        return
    end
    
    c0 = p(1);
    c = p(2);
    a = p(3);
    
    if c0 < 0 ...
            || c < 0 ...
            || a <= 0 ...
            || a > hmax_fit
        val = 1e30 + 1e6*(abs(min(c0,0)) ...
            + abs(min(c,0)) + abs(a));
        return
    end    
    gfit = local_spherical_variogram(h,p);
    res = g(:) - gfit(:);
    ww = w(:) ./ max(1,sum(w));    
    val = sum(ww .* res.^2);    
    % slight regularization to discourage extremely tiny ranges
    if a < 50
        val = val + 1e-6*(50 - a)^2;
    end
end

function g = local_spherical_variogram(h,p)
%LOCAL_SPHERICAL_VARIOGRAM Spherical variogram model.
%
% p = [nugget, partial_sill, range]

    c0 = p(1);
    c = p(2);
    a = p(3);    
    h = double(h);
    g = nan(size(h));    
    if ~isfinite(c0) ...
            || ~isfinite(c) ...
            || ~isfinite(a) ...
            || a <= 0
        g(:) = nan;
        return
    end    
    r = h ./ a;    
    inside = (h <= a);
    outside = ~inside;    
    g(inside) = c0 + c .* (1.5*r(inside) - 0.5*r(inside).^3);
    g(outside) = c0 + c;
end

function s = local_variogram_label_latex(mdl,j)
%LOCAL_VARIOGRAM_LABEL_LATEX Build latex parameter symbol from mdl.par_names.

    pname = '';
    try
        if isfield(mdl,'par_names') ...
                && numel(mdl.par_names) >= j
            if iscell(mdl.par_names)
                pname = char(string(mdl.par_names{j}));
            else
                pname = char(string(mdl.par_names(j)));
            end
        end
    catch
    end    
    if isempty(pname)
        s = sprintf('$\\theta_{%d}$',j);
        return
    end    
    [mainSym,subSym] = local_split_symbol(pname);    
    if isempty(subSym)
        s = ['$' mainSym '$'];
    else
        s = ['$' mainSym '_{' subSym '}$'];
    end
end

function [mainSym,subSym] = local_split_symbol(p)
%LOCAL_SPLIT_SYMBOL Split parameter name into main and sub symbols.

    p = strtrim(char(p));
    k = find(p == '_',1,'first');
    
    if isempty(k)
        mainRaw = p;
        subRaw = '';
    else
        mainRaw = strtrim(p(1:k-1));
        subRaw = strtrim(p(k+1:end));
    end
    
    mainSym = local_latex_main(mainRaw);
    subSym = local_latex_sub(subRaw);
    end
    
    function s = local_latex_main(raw)
    %LOCAL_LATEX_MAIN Latex formatting for main parameter symbol.
    
    raw = strtrim(char(raw));
    rl = lower(raw);
    
    if isempty(raw)
        s = '';
        return
    end
    
    if startsWith(raw,'\')
        s = raw;
        return
    end
    
    switch rl
        case 'alpha'
            s = '\alpha'; return
        case 'beta'
            s = '\beta'; return
        case 'gamma'
            s = '\gamma'; return
        case 'delta'
            s = '\delta'; return
        case 'epsilon'
            s = '\epsilon'; return
        case 'lambda'
            s = '\lambda'; return
        case 'omega'
            s = '\omega'; return
        case 'theta'
            s = '\theta'; return
        case 'phi'
            s = '\phi'; return
        case 'psi'
            s = '\psi'; return
        case 'mu'
            s = '\mu'; return
        case 'sigma'
            s = '\sigma'; return
        case 'ell'
            s = '\ell'; return
    end
    
    tok = regexp(raw,'^([A-Za-z]+)(\d+)$','tokens','once');
    if ~isempty(tok)
        head = tok{1};
        num = tok{2};
        if isscalar(head)
            s = [head '_{' num '}'];
        else
            s = ['{\rm ' head '}_{' num '}'];
        end
        return
    end
    
    switch rl
        case 'stot'
            s = 's'; return
        case 'smax'
            s = 's'; return
        case 'smin'
            s = 's'; return
        case 'cmax'
            s = 'c'; return
        case 'cmin'
            s = 'c'; return
        case 'ks'
            s = 'k'; return
        case 'ki'
            s = 'k'; return
        case 'kb'
            s = 'k'; return
        case 'kq'
            s = 'k'; return
        case 'uzl'
            s = 'u'; return
        case 'bexp'
            s = 'b'; return
        case 'qb'
            s = 'q'; return
        case 'qs'
            s = 'q'; return
        case 'qt'
            s = 'q'; return
        case 'ep'
            s = 'e'; return
        case 'pet'
            s = 'ET'; return
        case 'fp'
            s = 'f'; return
    end    
    if isscalar(raw)
        s = raw;
    else
        s = ['{\rm ' raw '}'];
    end
end

function s = local_latex_sub(raw)
%LOCAL_LATEX_SUB Latex formatting for parameter subscripts.

    raw = strtrim(char(raw));
    if isempty(raw)
        s = '';
        return
    end
    if startsWith(raw,'\')
        s = raw;
        return
    end
    if ~isempty(regexp(raw,'^\d+$','once'))
        s = raw;
        return
    end
    
    switch lower(raw)
        case 'tot'
            s = '\rm tot';
        case 'max'
            s = '\rm max';
        case 'min'
            s = '\rm min';
        case 's'
            s = '\rm s';
        case 'i'
            s = '\rm i';
        case 'b'
            s = '\rm b';
        case 'q'
            s = '\rm q';
        case 'zl'
            s = '\rm zl';
        case 'exp'
            s = '\rm exp';
        case 't'
            s = '\rm t';
        case 'p'
            s = '\rm p';
        otherwise
            s = ['\rm ' raw];
    end
end

function s = local_format_tick_with_commas(x)
%LOCAL_FORMAT_TICK_WITH_COMMAS Format integer tick labels with commas.

    if ~isfinite(x)
        s = '';
        return
    end
    x = round(x);
    neg = x < 0;
    x = abs(x);
    s = sprintf('%d',x);
    n = length(s);
    if n <= 3
        if neg
            s = ['-' s];
        end
        return
    end
    k = mod(n,3);
    if k == 0
        k = 3;
    end
    parts = {s(1:k)};
    for i = k+1:3:n
        parts{end+1} = s(i:min(i+2,n)); %#ok<AGROW>
    end    
    s = strjoin(parts,',');
    if neg
        s = ['-' s];
    end
end

function s = local_plain_label(mdl,j)
%LOCAL_PLAIN_LABEL Plain-text parameter label for command-window output.

    s = sprintf('theta_%d',j);    
    try
        if isfield(mdl,'par_names') ...
                && numel(mdl.par_names) >= j
            if iscell(mdl.par_names)
                s = char(string(mdl.par_names{j}));
            else
                s = char(string(mdl.par_names(j)));
            end
        end
    catch
    end
end

function YS = local_variogram_axis_scale(gj,fitOK,fitPar,hmax_fit,modelFcn)
%LOCAL_VARIOGRAM_AXIS_SCALE Choose pleasant y-axis scaling for variograms.
%
% This helper chooses a power-of-ten scale so that small variogram values
% are shown with integer-like y-axis ticks.
%
% Example:
%   ymax = 0.1201
%   scaleExp = 2
%   plotted y = 10^2 * gamma(h)
%   y-axis label = gamma(h) (x 10^{-2})
%   ticks = 0,2,4,6,8,10,12,14
% Thus, underlying gamma(h) = displayed tick value x 10^{-2}.
%
% INPUT
%   gj        empirical variogram values
%   fitOK     true/false, whether fitted model is available
%   fitPar    fitted variogram parameters
%   hmax_fit  maximum lag distance for fitted curve
%   modelFcn  function handle, e.g. @local_spherical_variogram
%
% OUTPUT
%   YS.scaleExp
%   YS.scaleFac
%   YS.ymax
%   YS.ymaxScaled
%   YS.ylim
%   YS.yticks
%   YS.yticklabels
%   YS.ylabelLatex

    if nargin < 5 ...
            || isempty(modelFcn)
        modelFcn = [];
    end

    yall = gj(:);

    if fitOK ...
            && ~isempty(modelFcn)
        try
            hfit = linspace(0,hmax_fit,200);
            yfit = modelFcn(hfit,fitPar);
            yall = [yall; yfit(:)];
        catch
        end
    end

    yall = yall(isfinite(yall) ...
        & yall >= 0);

    if isempty(yall)
        YS = local_empty_variogram_axis_scale();
        return
    end

    ymax = max(yall);

    if ~isfinite(ymax) ...
            || ymax <= 0
        YS = local_empty_variogram_axis_scale();
        return
    end

    % ---------------------------------------------------------
    % Choose exponent.
    %
    % We want ymax*10^scaleExp to land preferably around 2-20.
    % The extra decade when mantissa < 2 gives:
    %   0.1201 -> 12.01, not 1.201
    %   1.201  -> 12.01, not 1.201
    % ---------------------------------------------------------
    p = floor(log10(ymax));
    m = ymax / 10^p;

    if m < 2
        scaleExp = -p + 1;
    else
        scaleExp = -p;
    end

    scaleExp = max(0,scaleExp);
    scaleFac = 10^scaleExp;
    ymaxScaled = ymax * scaleFac;

    % -----------------------------------------------------
    % Nice integer ticks with maximum number of tick labels
    % -----------------------------------------------------
    maxTicks = 6;
    
    tickStep = local_nice_tick_step_limited( ...
        ymaxScaled,maxTicks);
    ytop = tickStep * ceil( ...
        ymaxScaled/tickStep);
    
    % Make sure the top is not zero.
    if ytop <= 0 ...
            || ~isfinite(ytop)
        ytop = 1;
    end
    yt = 0:tickStep:ytop;
    % If roundoff still gives too many ticks, increase the step.
    while numel(yt) > maxTicks
        tickStep = local_next_nice_step(tickStep);
        ytop = tickStep * ceil(ymaxScaled/tickStep);
        yt = 0:tickStep:ytop;
    end
    YS = struct();
    YS.scaleExp = scaleExp;
    YS.scaleFac = scaleFac;
    YS.ymax = ymax;
    YS.ymaxScaled = ymaxScaled;
    YS.ylim = [0 ytop];
    YS.yticks = yt;
    YS.yticklabels = compose('%g',yt);
    if scaleExp == 0
        YS.ylabelLatex = '$\gamma(h)$';
    else
        multiplierExp = -scaleExp;
        YS.ylabelLatex = sprintf( ...
            '$\\gamma(h)\\;(\\times 10^{%d})$',multiplierExp);
    end
end

function YS = local_empty_variogram_axis_scale()
%LOCAL_EMPTY_VARIOGRAM_AXIS_SCALE Returns default scaling for an empty variogram panel.

    YS = struct();
    YS.scaleExp = 0;
    YS.scaleFac = 1;
    YS.ymax = 0;
    YS.ymaxScaled = 0;
    YS.ylim = [0 1];
    YS.yticks = 0:0.2:1;
    YS.yticklabels = compose('%g',YS.yticks);
    YS.ylabelLatex = '$\gamma(h)$';

end

function step = local_nice_tick_step_limited(ymaxScaled,maxTicks)
%LOCAL_NICE_TICK_STEP_LIMITED Selects a readable tick interval subject to a tick limit.

    if nargin < 2 ...
            || isempty(maxTicks)
        maxTicks = 6;
    end
    if ~isfinite(ymaxScaled) ...
            || ymaxScaled <= 0
        step = 1;
        return
    end
    % We need at most maxTicks labels including zero.
    % Thus the maximum number of intervals is maxTicks - 1.
    raw = ymaxScaled / max(1,maxTicks - 1);
    p = floor(log10(raw));
    m = raw / 10^p;
    if m <= 1
        nice = 1;
    elseif m <= 2
        nice = 2;
    elseif m <= 5
        nice = 5;
    else
        nice = 10;
    end
    step = nice * 10^p;
end

function step2 = local_next_nice_step(step)
%LOCAL_NEXT_NICE_STEP Advances to the next conventional numeric tick interval.

    if ~isfinite(step) ...
            || step <= 0
        step2 = 1;
        return
    end
    p = floor(log10(step));
    m = step / 10^p;
    if m < 2
        m2 = 2;
    elseif m < 5
        m2 = 5;
    elseif m < 10
        m2 = 10;
    else
        m2 = 20;
    end
    step2 = m2 * 10^p;
end

function [xt,xtlbl] = local_variogram_xticks(xmaxPlot,maxTicks)
%LOCAL_VARIOGRAM_XTICKS Nice sparse x ticks for variogram subplots.

    if nargin < 2 ...
            || isempty(maxTicks)
        maxTicks = 4;
    end

    if ~isfinite(xmaxPlot) ...
            || xmaxPlot <= 0
        xt = [0 1];
    else
        rawStep = xmaxPlot / max(1,maxTicks - 1);
        step = local_nice_distance_step(rawStep);
        xt = 0:step:xmaxPlot;
        % If still too many labels, increase step.
        while numel(xt) > maxTicks
            step = local_next_distance_step(step);
            xt = 0:step:xmaxPlot;
        end
    end
    xtlbl = cell(size(xt));
    for i = 1:numel(xt)
        xtlbl{i} = local_format_tick_with_commas(xt(i));
    end
end

function maxTicks = local_variogram_max_xticks( ...
        fig,panelWidthNormalized,fontSize,xmaxPlot)
%LOCAL_VARIOGRAM_MAX_XTICKS Tick capacity from rendered subplot width.

    try
        figPixels = getpixelposition(fig,true);
        panelPixels = max(120, ...
            double(figPixels(3))*double(panelWidthNormalized));

        endpointLabels = { ...
            local_format_tick_with_commas(0), ...
            local_format_tick_with_commas(xmaxPlot)};
        maxChars = max(cellfun(@numel,endpointLabels));

        % Approximate the widest label at the active axes font. The extra
        % 20 pixels provide a visible gap rather than merely avoiding an
        % exact bounding-box collision.
        labelPixels = 0.62*double(fontSize)*maxChars;
        maxTicks = floor(panelPixels/max(labelPixels + 20,1));
        maxTicks = max(2,min(6,maxTicks));
    catch
        maxTicks = 4;
    end
end

function step = local_nice_distance_step(rawStep)
%LOCAL_NICE_DISTANCE_STEP Selects a readable distance-axis tick interval.

    if ~isfinite(rawStep) ...
            || rawStep <= 0
        step = 1000;
        return
    end
    exponent = floor(log10(rawStep));
    fraction = rawStep / 10^exponent;
    niceFractions = [1 2 2.5 4 5 10];
    k = find(niceFractions >= fraction,1,'first');
    step = niceFractions(k) * 10^exponent;
end

function step2 = local_next_distance_step(step)
%LOCAL_NEXT_DISTANCE_STEP Advances to the next conventional distance interval.

    exponent = floor(log10(step));
    fraction = step / 10^exponent;
    niceFractions = [1 2 2.5 4 5 10];
    k = find(niceFractions > fraction*(1 + 10*eps),1,'first');
    if isempty(k)
        step2 = 2 * 10^(exponent + 1);
    else
        step2 = niceFractions(k) * 10^exponent;
    end
end

function s = local_variogram_label_normalized_latex(mdl,j)
%LOCAL_VARIROGRAM_LABEL_NORMALIZED_LATEX Normalized parameter label.
%
% Underlines the main parameter symbol and preserves the subscript.
%
% Examples:
%   s_tot   -> \underline{s}_{\rm tot}
%   k_s     -> \underline{k}_{\rm s}
%   sf_cf   -> \underline{\rm sf}_{\rm cf}
%   alpha   -> \underline{\alpha}

    pname = '';

    try
        if isfield(mdl,'par_names') ...
                && numel(mdl.par_names) >= j
            if iscell(mdl.par_names)
                pname = char(string(mdl.par_names{j}));
            else
                pname = char(string(mdl.par_names(j)));
            end
        end
    catch
    end

    if isempty(pname)
        s = sprintf('$\\underline{\\theta}_{%d}$',j);
        return
    end

    [mainSym,subSym] = local_split_symbol(pname);

    if isempty(subSym)
        s = ['$\underline{' mainSym '}$'];
    else
        s = ['$\underline{' mainSym '}_{' subSym '}$'];
    end

end

function VG = plot_parameter_variograms_SITE_figure( ...
    figNo,mdl,latlon,nTheta,region)
%PLOT_PARAMETER_VARIOGRAMS_SITE_FIGURE Plots SITE parameter variograms.
% Plot empirical spatial variograms for all SITE models in one figure.
%
% SYNOPSIS:
%   VG = plot_parameter_variograms_SITE_figure( ...
%       figNo,mdl,latlon,nTheta)
%
% INPUT:
%   figNo   figure number
%   mdl     model structure with .names and .par_names
%   latlon  K x 2 [lat lon]
%   nTheta     d x K x n_m normalized parameter values
%
% OUTPUT:
%   VG      1 x n_m struct array with fields:
%             .binEdges
%             .binCenter
%             .gamma
%             .npairs
%             .fitPar
%             .fitOK
%             .sill
%             .label_plain
%             .label_latex
%             .hmax_fit
%             .xmax
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 5 ...
            || isempty(region)
        region = 'US';
    end
    
    if isempty(latlon) ...
            || isempty(nTheta)
        VG = struct([]);
        return
    end
    
    if ndims(nTheta) ~= 3
        error(['plot_parameter_variograms_SITE_figure: ' ...
            'nTheta must be d x K x n_m for SITE.']);
    end
    
    [d,K,n_m] = size(nTheta);
    
    if size(latlon,1) ~= K ...
            || size(latlon,2) < 2
        error(['plot_parameter_variograms_SITE_figure: ' ...
            'latlon must be K x 2 and match size(nTheta,2).']);
    end
    
    % --------------------------------
    % User settings for variogram plot
    % --------------------------------
    nBins = 30;         % number of empirical lag bins
    hmax_fit = local_region_helpers_plot( ...
        'variogramrange',region);
    
    % ------
    % Colors
    % ------
    colors = hsv(d);
    
    % -------------------------------
    % Compute pairwise distances once
    % -------------------------------
    Dkm = local_great_circle_distance_km( ...
        latlon(:,1),latlon(:,2));
    iu = triu(true(K),1);
    distVec = Dkm(iu);
    
    maxDistKm = max(distVec);
    if ~isfinite(maxDistKm) ...
            || maxDistKm <= 0
        error(['plot_parameter_variograms_SITE_figure: ' ...
            'Invalid basin distances.']);
    end
    
    binEdges = linspace(0,maxDistKm,nBins+1);
    binCenter = 0.5*(binEdges(1:end-1) ...
        + binEdges(2:end));
    binID = discretize(distVec,binEdges);
    
    npairs = zeros(1,nBins);
    for b = 1:nBins
        npairs(b) = sum(binID == b);
    end
    
    xmax = min(1.15*hmax_fit,max(binCenter));
    
    % ---------------------------
    % Compute/fits for each model
    % ---------------------------
    VG = repmat(struct( ...
        'binEdges',[], ...
        'binCenter',[], ...
        'gamma',[], ...
        'npairs',[], ...
        'fitPar',[], ...
        'fitOK',[], ...
        'sill',[], ...
        'label_plain',[], ...
        'label_latex',[], ...
        'hmax_fit',[], ...
        'xmax',[]), 1,n_m);
    
    for im = 1:n_m
        gamma = nan(d,nBins);
    
        for b = 1:nBins
            mask = (binID == b);
            if ~any(mask)
                continue
            end
    
            for j = 1:d
                z = double(nTheta(j,:,im)).';
                DZ = z - z.';
                dzVec = DZ(iu);
                gamma(j,b) = 0.5 * ...
                    mean(dzVec(mask).^2,'omitnan');
            end
        end
    
        fitPar = nan(d,3);   % [nugget, partial sill, range]
        fitOK = false(d,1);
    
        for j = 1:d
            gj = gamma(j,:);
    
            good = isfinite(gj) ...
                & isfinite(binCenter) ...
                & (npairs > 0) ...
                & (binCenter <= hmax_fit);
    
            h = binCenter(good);
            g = gj(good);
            w = npairs(good);
    
            if numel(h) < 4
                continue
            end
    
            gmax = max(g);
            if ~isfinite(gmax)
                gmax = 1; 
            end
    
            nug0 = max(0,0.05*gmax);
            psill0 = max(1e-10,gmax - nug0);
            range0 = min(hmax_fit,0.5*max(h));
            if ~isfinite(range0) || range0 <= 0
                range0 = 1000;
            end
    
            p0 = [nug0, psill0, range0];
    
            objfun = @(p) local_variogram_objective( ...
                p,h,g,w,hmax_fit);
    
            opts = optimset('Display','off', ...
                            'MaxFunEvals',5000, ...
                            'MaxIter',5000, ...
                            'TolX',1e-8, ...
                            'TolFun',1e-8);
    
            try
                phat = fminsearch(objfun,p0,opts);
                phat(1) = max(0,phat(1));
                phat(2) = max(0,phat(2));
                phat(3) = ...
                    max(1,min(hmax_fit,phat(3)));
    
                if all(isfinite(phat))
                    fitPar(j,:) = phat;
                    fitOK(j) = true;
                end
            catch
            end
        end
    
        label_plain = cell(d,1);
        label_latex = cell(d,1);
        for j = 1:d
            label_plain{j} = ...
                local_plain_label(mdl,j);
            label_latex{j} = ...
                local_variogram_label_latex(mdl,j);
        end
    
        VG(im).binEdges = binEdges;
        VG(im).binCenter = binCenter;
        VG(im).gamma = gamma;
        VG(im).npairs = npairs;
        VG(im).fitPar = fitPar;
        VG(im).fitOK = fitOK;
        VG(im).sill = fitPar(:,1) + fitPar(:,2);
        VG(im).label_plain = label_plain;
        VG(im).label_latex = label_latex;
        VG(im).hmax_fit = hmax_fit;
        VG(im).xmax = xmax;
    end
    
    % -------------
    % Figure layout
    % -------------
    ncol_top = 2;
    nrow_top = ceil(n_m / ncol_top);
    
    figure(figNo); clf;
    set(gcf, ...
        'Name',['SITE: Spatial variograms ' ...
        'of final parameter values'], ...
        'NumberTitle','off', ...
        'color','w', ...
        'Units','inches', ...
        'Position',[0.3 0.3 15 3.4*nrow_top + 2.0]);
    
    tlo = tiledlayout(nrow_top+1,ncol_top, ...
        'TileSpacing','compact', ...
        'Padding','compact');
    
    % ----------------------------------------
    % Top tiles: one variogram panel per model
    % ----------------------------------------
    for im = 1:n_m
        ax = nexttile(tlo,im);
        hold(ax,'on');
        box(ax,'on');
    
        gamma = VG(im).gamma;
        fitPar = VG(im).fitPar;
        fitOK = VG(im).fitOK;
    
        for j = 1:d
            gj = gamma(j,:);
    
            % empirical points only up to hmax_fit
            good = isfinite(gj) ...
                & isfinite(binCenter) ...
                & (npairs > 0) ...
                & (binCenter <= hmax_fit);
    
            if any(good)
                plot(ax,binCenter(good),gj(good), ...
                    's', ...
                    'color',colors(j,:), ...
                    'MarkerFaceColor',colors(j,:), ...
                    'MarkerEdgeColor',colors(j,:), ...
                    'markersize',6, ...
                    'linewidth',1.2, ...
                    'handlevisibility','off');
            end
    
            % fitted spherical curve
            if fitOK(j)
                hfit = linspace(0,hmax_fit,300);
                gfit = ...
                    local_spherical_variogram( ...
                    hfit,fitPar(j,:));
    
                plot(ax,hfit,gfit, ...
                    '-', ...
                    'color',colors(j,:), ...
                    'linewidth',2.0, ...
                    'handlevisibility','off');
            end
        end
    
        set(ax, ...
            'tickdir','out', ...
            'fontsize',11, ...
            'linewidth',1);
    
        ylabel(ax,'$\gamma(h)$', ...
            'interpreter','latex', ...
            'fontsize',14);
    
        if im > (nrow_top-1)*ncol_top
            xlabel(ax,'$h$ (km)', ...
                'interpreter','latex', ...
                'fontsize',14);
        else
            set(ax,'XTickLabel',[]);
        end
    
        grid(ax,'off');
        xlim(ax,[0 xmax]);
    
        tickStep = local_region_helpers_plot( ...
            'variogramtick',region);
        xt = 0:tickStep:ceil(xmax/tickStep)*tickStep;
        xtlbl = cell(size(xt));
        for i = 1:numel(xt)
            xtlbl{i} = local_format_tick_with_commas(xt(i));
        end
        set(ax, ...
            'XTick',xt, ...
            'XTickLabel',xtlbl, ...
            'TickLabelInterpreter','latex');
    
        ax.XRuler.TickLabelGapOffset = -2;
        ax.YRuler.TickLabelGapOffset = -2;
    
        % robust y-limits per panel
        yall = gamma(:);
        for j = 1:d
            if fitOK(j)
                yall = [yall; ...
                    local_spherical_variogram( ...
                    linspace(0,hmax_fit,50), ...
                    fitPar(j,:)).']; %#ok
            end
        end
        yall = yall(isfinite(yall));
    
        if isempty(yall)
            ylim(ax,[0 1]);
        else
            ylo = min(yall);
            yhi = max(yall);
            if yhi <= ylo
                dy = max(1e-6,0.01*max(abs([ylo yhi])));
                ylo = ylo - dy;
                yhi = yhi + dy;
            end
            pad = 0.1*(yhi - ylo);
            ylim(ax,[max(0,ylo-pad) yhi+pad]);
        end
    
        % title per panel
        if isfield(mdl,'names') ...
                && numel(mdl.names) >= im
            ttl = upper(char(string(mdl.names(im))));
        else
            ttl = sprintf('Model %d',im);
        end
        title(ax,ttl, ...
            'interpreter','none', ...
            'fontsize',12);
    
        % manual legend in white block
        % xL1 = 3820;
        % xL2 = 3960;
        % xT = 4010;
        xL1 = 0.90*xmax;
        xL2 = 0.94*xmax;
        xT = 0.955*xmax;
        
        yl = ylim(ax);
        yminL = yl(1) + 0.10*(yl(2)-yl(1));
        ymaxL = yl(1) + 0.92*(yl(2)-yl(1));
    
        if d == 1
            yLegend = 0.5*(yminL+ymaxL);
        else
            yLegend = linspace(ymaxL,yminL,d);
        end
    
        for j = 1:d
            plot(ax,[xL1 xL2],[yLegend(j) yLegend(j)], ...
                '-', ...
                'color',colors(j,:), ...
                'linewidth',2.0, ...
                'clipping','off', ...
                'handlevisibility','off');
    
            text(ax,xT,yLegend(j),VG(im).label_latex{j}, ...
                'interpreter','latex', ...
                'color',colors(j,:), ...
                'fontsize',14, ...
                'horizontalalignment','left', ...
                'verticalalignment','middle', ...
                'clipping','off');
        end
    end
    
    % turn off unused tiles in top block
    for im = n_m+1:nrow_top*ncol_top
        ax = nexttile(tlo,im);
        axis(ax,'off');
    end
    
    % --------------------------------------------------
    % Bottom row: common histogram spanning both columns
    % --------------------------------------------------
    axH = nexttile(tlo, nrow_top*ncol_top + 1, [1 2]);
    hold(axH,'on'); box(axH,'on');
    
    bar(axH,binCenter,npairs,0.7, ...
        'FaceColor',[0.4 0.4 0.4], ...
        'EdgeColor','none');
    
    set(axH, ...
        'tickdir','out', ...
        'fontsize',16, ...
        'linewidth',1);
    
    ylabel(axH,'Number of pairs', ...
        'interpreter','latex', ...
        'fontsize',18);
    
    xlabel(axH,'Lag distance, $h$ (km)', ...
        'interpreter','latex', ...
        'fontsize',18);
    
    xlim(axH,[0 xmax]);
    
    xt = 0:tickStep:ceil(xmax/tickStep)*tickStep;
    xtlbl = cell(size(xt));
    for i = 1:numel(xt)
        xtlbl{i} = local_format_tick_with_commas(xt(i));
    end
    
    set(axH, ...
        'XTickMode','manual', ...
        'XTick',xt, ...
        'XTickLabelMode','manual', ...
        'XTickLabel',xtlbl, ...
        'TickLabelInterpreter','latex');
    
    axH.XRuler.TickLabelGapOffset = -2;
    axH.YRuler.TickLabelGapOffset = -2;
    
    ymax_pairs = max(npairs);
    if isempty(ymax_pairs) ...
            || ~isfinite(ymax_pairs) ...
            || ymax_pairs <= 0
        ymax_pairs = 1;
    end
    ylim(axH,[0 ymax_pairs + 0.10*ymax_pairs]);
    
    % -----
    % Title
    % -----
    annotation('textbox',[0 0.965 1 0.03], ...
        'String',['\texttt{SITE}: ' ...
        'Spatial variograms of final ' ...
        'parameter values'], ...
        'interpreter','latex', ...
        'EdgeColor','none', ...
        'horizontalalignment','center', ...
        'verticalalignment','middle', ...
        'fontsize',15);

end

function local_plot_timeseries_postprocessor( ...
    figBase,part,mdl,dat,bas,prd,Q,Qfdc,gaugescen, ...
    q_unit,t_unit,dt_str,n_m,id_model,colors, ...
    model_names,id_dat,Bprt,n_char,figW,figH, ...
    model_name_fig,sp_method,fontsize_legend, ...
    cPtTrain,cPtEval,cShadeTrain,cShadeEval, ...
    face_alpha,marker_alpha,region)

%LOCAL_PLOT_TIMESERIES_POSTPROCESSOR Unified postprocessor for all split
% methods. Training observations are always training data (id_train) and
% evaluation observations are always evaluation data (id_eval).
% Visual rule:
%   - use background shading if either train/eval contains a contiguous
%     block of at least 30 points
%   - otherwise use colored observation markers
%   - force point-style for pure random sampling

    K_t = bas.K_t;
    K_e = bas.K_e;
    showForcingPanels = local_show_forcing_panels();
    figW = 0.95*figW;
    figH = 0.95*figH;
    id_train = expand_index(mdl.id_train);
    
    if isfield(mdl,'id_eval') ...
            && ~isempty(mdl.id_eval)
        id_eval = expand_index(mdl.id_eval);
    else
        id_eval = [];
    end
    
    groups = {};
    
    S = struct();
    S.tag = 'train_mask';
    S.title = sprintf(['%s training basin ' ...
        '| full scored window'],dt_str);
    S.id_global = 1:K_t;
    S.codes = local_collect_group_codes( ...
        gaugescen,'train');
    S.Qtr = local_get_metric(Q,'tt',[]);
    S.Qev = local_get_metric(Q,'te',[]);
    groups{end+1} = S;
    
    S = struct();
    S.tag = 'eval_mask';
    S.title = sprintf(['%s evaluation basin ' ...
        '| full scored window'],dt_str);
    S.id_global = (K_t+1):(K_t+K_e);
    S.codes = local_collect_group_codes( ...
        gaugescen,'eval');
    S.Qtr = local_get_metric(Q,'et',[]);
    S.Qev = local_get_metric(Q,'ee',[]);
    groups{end+1} = S;
    
    for is = 1:numel(groups)
        G = groups{is};
    
        if isempty(G.id_global)
            continue
        end
    
        if isempty(G.codes)
            % No explicit gaugescen request was supplied for this group.
            % Before random Bprt selection, reuse the FDC-selected dat{k}
            % indices when Qfdc.id exists.  This keeps time-series and FDC
            % figures focused on the same basins.
            pick = local_pick_qfdc_ids(Qfdc,G.id_global,'');
            if isempty(pick)
                pick = pick_basins([],id_dat, ...
                    G.id_global, ...
                    min(Bprt,numel(G.id_global)),region);
            end
        else
            % Explicit gaugescen/user selection has highest priority.
            pick = pick_basins(G.codes,id_dat, ...
                G.id_global, ...
                min(Bprt,numel(G.id_global)),region);
        end
    
        if isempty(pick)
            continue
        end
    
        nPerFig = 2;   % full-width panels
        nFig = ceil(numel(pick)/nPerFig);
    
        for jf = 1:nFig
            f = figure(figBase);
            clf(f,'reset');
            figBase = figBase + 1;
    
            nThisFig = min(nPerFig, ...
                numel(pick) - (jf-1)*nPerFig);
            nRows = max(1,nThisFig);
            figH_page = figH * (0.56 + 0.44*(nRows == 2));
    
            set(f, ...
                'Name',sprintf('%s: Time series - %s (%d/%d)', ...
                model_name_fig,G.title,jf,nFig), ...
                'NumberTitle','off', ...
                'color','w', ...
                'Units','inches', ...
                'Position',[0.3 0.3 figW figH_page], ...
                'Visible','on');
    
            axh = gobjects(nRows,1);
            axForcing = gobjects(nRows,1);
            axTemperature = gobjects(nRows,1);
    
            for jp = 1:nRows
                kk = (jf-1)*nPerFig + jp;
                if showForcingPanels
                    [pF,pT,pQ] = local_timeseries_axes_positions( ...
                        nRows,1,jp);
                    axForcing(jp) = axes(f, ...
                        'Units','normalized','Position',pF, ...
                        'PositionConstraint','innerposition');
                    axTemperature(jp) = axes(f, ...
                        'Units','normalized','Position',pT, ...
                        'PositionConstraint','innerposition');
                    axh(jp) = axes(f, ...
                        'Units','normalized','Position',pQ, ...
                        'PositionConstraint','innerposition');
                else
                    [~,~,pQ] = local_timeseries_axes_positions( ...
                        nRows,1,jp);
                    axh(jp) = axes(f, ...
                        'Units','normalized','Position',pQ, ...
                        'PositionConstraint','innerposition');
                end
                box(axh(jp),'on');
                hold(axh(jp),'on');
                set(axh(jp), ...
                    'Layer','top', ...
                    'TickLength',[0.004 0.004]);
    
                if kk > numel(pick)
                    axis(axh(jp),'off');
                    continue
                end
    
                kdat = pick(kk);
    
                switch lower(G.tag)
                    case 'train_mask'
                        qCol = kdat;
                    case 'eval_mask'
                        qCol = kdat - K_t;
                    otherwise
                        qCol = kdat;
                end
    
                if strcmpi(sp_method,'rainfall_block')
                    id_train_k = local_safe_index( ...
                        dat{kdat},'id_train');
                    id_eval_k = local_safe_index( ...
                        dat{kdat},'id_eval');
                else
                    id_train_k = id_train;
                    id_eval_k = id_eval;
                end
                
                local_plot_basin_postprocessor(axh(jp), ...
                    dat{kdat},prd,mdl,qCol, ...
                    G.Qtr,G.Qev,id_train_k, ...
                    id_eval_k,sp_method, ...
                    n_m,id_model,colors, ...
                    q_unit,t_unit, ...
                    cPtTrain, ...
                    cPtEval, ...
                    cShadeTrain, ...
                    cShadeEval, ...
                    face_alpha, ...
                    marker_alpha);

                % The discharge panel defines the displayed time window.
                % Preserve its scored-window limits so that plotting the
                % full forcing vectors cannot reintroduce the spin-up.
                qXLim = xlim(axh(jp));
                qXTick = get(axh(jp),'XTick');

                if showForcingPanels
                    local_plot_forcing_postprocessor( ...
                        axForcing(jp),dat{kdat}, ...
                        id_train_k,id_eval_k,q_unit, ...
                        sp_method,cShadeTrain,cShadeEval, ...
                        face_alpha);
                    local_plot_temperature_strip( ...
                        axTemperature(jp),dat{kdat}, ...
                        id_train_k,id_eval_k,sp_method, ...
                        cShadeTrain,cShadeEval,face_alpha, ...
                        qXLim,axh(jp).LineWidth);
                    xlim(axForcing(jp),qXLim);
                    linkaxes([axForcing(jp), ...
                        axTemperature(jp),axh(jp)],'x');
                    xlim(axh(jp),qXLim);
                    set(axForcing(jp), ...
                        'XTick',qXTick, ...
                        'XTickLabel',[]);
                end

                gauge = local_gauge_name(bas,kdat);
                us = local_get_usgs(id_dat,kdat);
                zoneTxt = local_zone_label_short(bas,kdat);
                invalidForcing = local_incomplete_forcing(dat{kdat});
                labelColor = [0.15 0.15 0.15];
                if invalidForcing
                    labelColor = [0.80 0 0];
                end
                
                if showForcingPanels
                    add_basin_label(axForcing(jp),us, ...
                        zoneTxt,gauge,n_char,0.5,labelColor);
                else
                    add_basin_label(axh(jp),us, ...
                        zoneTxt,gauge,n_char,0,labelColor);
                end
                if invalidForcing
                    local_add_invalid_forcing_notice(axh(jp));
                end
            end

            drawnow limitrate nocallbacks

            goodAx = axh(isgraphics(axh));
            set(goodAx, ...
                'tickdir','out', ...
                'TickLength',[0.004 0.004], ...
                'fontsize',12);

            if showForcingPanels
                goodForcing = axForcing(isgraphics(axForcing));
                set(goodForcing, ...
                    'tickdir','out', ...
                    'TickLength',[0.004 0.004], ...
                    'fontsize',12);
                goodTemperature = axTemperature( ...
                    isgraphics(axTemperature));
                set(goodTemperature, ...
                    'XTick',[], ...
                    'YTick',[]);
            end
    
            % ----- proxy legend so marker sizes can be controlled -----
            axLeg = [];
            if ~isempty(goodAx)
                vis = get(goodAx,'Visible');
                if iscell(vis)
                    iv = find(strcmp(vis,'on'), ...
                        1,'first');
                else
                    if strcmp(vis,'on')
                        iv = 1;
                    else
                        iv = [];
                    end
                end
            
                if ~isempty(iv)
                    axLeg = goodAx(iv);
                else
                    axLeg = goodAx(1);
                end
            end
    
            if ~isempty(axLeg)
                [hLeg,lblLeg] = local_make_postproc_legend( ...
                    axLeg,n_m,id_model,colors,cPtTrain, ...
                    cPtEval,model_names);
            
                legend(axLeg,hLeg,lblLeg, ...
                    'Location','northeast', ...
                    'interpreter','latex', ...
                    'fontsize',fontsize_legend, ...
                    'box','off');
            end
    
            if nRows == 1
                local_position_single_basin_ylabels( ...
                    axForcing(1),axh(1));
            end

            model_name_tex = local_latex_escape(model_name_fig);
            local_timeseries_figure_title(f,sprintf( ...
                '\\texttt{%s}: %s --- figure %d/%d', ...
                model_name_tex,G.title,jf,nFig));
    
            drawnow limitrate nocallbacks
            figure(f);
            drawnow limitrate nocallbacks
        end
    end
end

function local_plot_basin_postprocessor(axh, ...
    datk,prd,mdl,qCol,Qtr,Qev,id_train,id_eval, ...
    sp_method,n_m,id_model,colors,q_unit, ...
    t_unit,cPtTrain,cPtEval,cShadeTrain, ...
    cShadeEval,face_alpha,marker_alpha)
%LOCAL_PLOT_BASIN_POSTPROCESSOR Plots discharge for one basin and scenario.

    fontsize_axis = 18;
    lw_sim = 1.25;
    
    cObsEdge = [0 0 0]; %#ok
    
    minRunForShade = 30;
    
    y = datk.y_n(:);
    bad = datk.bad(:);
    
    if isempty(y)
        axis(axh,'off');
        return
    end
    
    xT = id_train(:).';
    xE = id_eval(:).';
    
    xT = xT(xT >= 1 & xT <= numel(y));
    xE = xE(xE >= 1 & xE <= numel(y));
    
    idx_all = unique([xT(:); xE(:)]).';
    if isempty(idx_all)
        axis(axh,'off');
        return
    end
    
    if ~isempty(bad)
        if ~isempty(xT) ...
                && numel(bad) >= max(xT)
            xT_plot = xT(~bad(xT));
        else
            xT_plot = xT;
        end
    
        if ~isempty(xE) ...
                && numel(bad) >= max(xE)
            xE_plot = xE(~bad(xE));
        else
            xE_plot = xE;
        end
    else
        xT_plot = xT;
        xE_plot = xE;
    end
    
    % decide whether to use shading
    useShading = local_should_use_shading( ...
        sp_method,xT_plot,xE_plot, ...
        minRunForShade);
    
    % y-limits from obs + sims
    yy_all = [];
    if ~isempty(xT_plot)
        yy_all = [yy_all; y(xT_plot)]; 
    end
    if ~isempty(xE_plot)
        yy_all = [yy_all; y(xE_plot)]; 
    end
    
    for mdl_i = 1:n_m
        q = [];
        if local_has_data(Qtr)
            q = local_get_discharge(Qtr,qCol,mdl_i);
        elseif local_has_data(Qev)
            q = local_get_discharge(Qev,qCol,mdl_i);
        end
    
        q = q(:);
        if isempty(q)
            continue
        end
    
        if numel(q) == numel(idx_all)
            q_sim = q;
        elseif numel(q) >= max(idx_all)
            q_sim = q(idx_all);
        else
            n = min(numel(q),numel(idx_all));
            q_sim = q(1:n);
        end
    
        q_sim = q_sim(isfinite(q_sim) ...
            & isreal(q_sim));
        yy_all = [yy_all; q_sim(:)]; %#ok
    end
    
    yy_all = yy_all(isfinite(yy_all));
    
    if isempty(yy_all)
        yl = [0 1];
    else
        ymin = min(yy_all);
        ymax = max(yy_all);
        yr = ymax - ymin;
    
        if yr <= 0
            yr = max(1e-6,abs(ymax));
        end
    
        ylo = ymin - 0.03*yr;
        yhi = ymax + 0.03*yr;
    
        if ylo < 0
            ylo = 0.5*ylo;
        end
    
        yl = [ylo yhi];
    end
    
    hold(axh,'on');
    
    % ----------------------------
    % background shading if useful
    % ----------------------------
    if useShading
        runsT = local_contiguous_runs(xT_plot);
        runsE = local_contiguous_runs(xE_plot);
    
        for r = 1:size(runsT,1)
            x1 = runsT(r,1) - 0.5;
            x2 = runsT(r,2) + 0.5;
            patch(axh,[x1 x2 x2 x1], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                cShadeTrain, ...
                'FaceAlpha',face_alpha, ...
                'EdgeColor','none', ...
                'handlevisibility','off');
        end
    
        for r = 1:size(runsE,1)
            x1 = runsE(r,1) - 0.5;
            x2 = runsE(r,2) + 0.5;
            patch(axh,[x1 x2 x2 x1], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                cShadeEval, ...
                'FaceAlpha',face_alpha, ...
                'EdgeColor','none', ...
                'handlevisibility','off');
        end
    
        % observations as open filled circles
        if ~isempty(xT_plot)
            scatter(axh,xT_plot,y(xT_plot), ...
                6, ...
                cPtTrain, ...
                'filled', ...
                'MarkerFaceAlpha',marker_alpha, ...
                'MarkerEdgeAlpha',marker_alpha, ...
                'handlevisibility','off');
        end
    
        if ~isempty(xE_plot)
            scatter(axh,xE_plot,y(xE_plot), ...
                6, ...
                cPtEval, ...
                'filled', ...
                'MarkerFaceAlpha',marker_alpha, ...
                'MarkerEdgeAlpha',marker_alpha, ...
                'handlevisibility','off');
        end
    
    else
        % no shading -> use colored obs points
        if ~isempty(xT_plot)
            scatter(axh,xT_plot,y(xT_plot), ...
                6, ...
                cPtTrain, ...
                'filled', ...
                'MarkerFaceAlpha',marker_alpha, ...
                'MarkerEdgeAlpha',marker_alpha, ...
                'handlevisibility','off');
        end
    
        if ~isempty(xE_plot)
            scatter(axh,xE_plot,y(xE_plot), ...
                6, ...
                cPtEval, ...
                'filled', ...
                'MarkerFaceAlpha',marker_alpha, ...
                'MarkerEdgeAlpha',marker_alpha, ...
                'handlevisibility','off');
        end
    end
    
    % simulations on top
    for mdl_i = 1:n_m
        q = [];
    
        if local_has_data(Qtr)
            q = local_get_discharge(Qtr,qCol,mdl_i);
        elseif local_has_data(Qev)
            q = local_get_discharge(Qev,qCol,mdl_i);
        end
    
        q = q(:);
        if isempty(q)
            continue
        end
    
        if numel(q) == numel(idx_all)
            x_sim = idx_all;
            q_sim = q;
        elseif numel(q) >= max(idx_all)
            x_sim = idx_all;
            q_sim = q(idx_all);
        else
            n = min(numel(q),numel(idx_all));
            x_sim = idx_all(1:n);
            q_sim = q(1:n);
        end
    
        good_sim = isfinite(q_sim) ...
            & isreal(q_sim);
        x_sim = x_sim(good_sim);
        q_sim = q_sim(good_sim);
    
        if ~isempty(x_sim)
            im = id_model(min(mdl_i, ...
                numel(id_model)));
            plot(axh,x_sim,q_sim, ...
                'linewidth',lw_sim, ...
                'color',colors(im,:), ...
                'handlevisibility','off');
        end
    end
    
    x1 = idx_all(1);
    x2 = idx_all(end);
    
    xlim(axh,[x1-0.5 x2+0.5]);
    useDateAxis = local_use_dates_on_xaxis(prd,mdl);    
    
    if useDateAxis
        idx_all = unique([id_train(:); id_eval(:)]);
        idx_all = idx_all(isfinite(idx_all));
    
        if ~isempty(idx_all)
            xmin = min(idx_all);
            xmax = max(idx_all);
    
            xlim(axh,[xmin xmax]);
    
            local_set_date_ticks_fourgroups( ...
                axh,'tt',prd,idx_all,xmin,xmax);
            set(axh,'XTickLabelRotation',0, ...
                'fontsize',12);
            xlabel(axh,'Date', ...
                'interpreter','latex', ...
                'fontsize',18);
        end
    else
        local_set_integer_ticks(axh,x1,x2);
        xlabel(axh,sprintf('${\\rm time\\; (%s)}$',t_unit), ...
            'interpreter','latex', ...
            'fontsize',fontsize_axis);
    end
    
    ylim(axh,yl);
    grid(axh,'off');
    
    hY = ylabel(axh,sprintf('$Q$ (%s)',q_unit), ...
        'interpreter','latex', ...
        'fontsize',fontsize_axis);
    local_position_left_ylabel(hY);

end

function [xL,xR] = local_metric_ecdf_limits(scenarios,metricTag)
%LOCAL_METRIC_ECDF_LIMITS Fixed efficiency limits or rounded FDC limits.

    if ~startsWith(strtrim(metricTag),'D','IgnoreCase',true)
        xL = -1;
        xR = 1;
        return
    end

    values = [];
    for i = 1:numel(scenarios)
        z = double(scenarios{i}.metric(:));
        values = [values; z(isfinite(z) & z >= 0)]; %#ok<AGROW>
    end
    xL = 0;
    if isempty(values) || max(values) <= 0
        xR = 1;
        return
    end
    peak = max(values);
    scale = 10^floor(log10(peak));
    choices = scale*[1 2 5 10];
    xR = choices(find(choices >= peak,1,'first'));
end

function showForcingPanels = local_show_forcing_panels()
%LOCAL_SHOW_FORCING_PANELS Select the time-series slide layout.
% Set this value to false to restore the original discharge-only panels.

    showForcingPanels = true;
end

function local_plot_forcing_postprocessor( ...
    axh,datk,id_train,id_eval,forcingUnit, ...
    sp_method,cShadeTrain,cShadeEval,face_alpha)
%LOCAL_PLOT_FORCING_POSTPROCESSOR Plot precipitation and PET above Q.

    if ~isfield(datk,'meteo') ...
            || ~isstruct(datk.meteo) ...
            || ~isfield(datk.meteo,'P')
        axis(axh,'off');
        return
    end

    P = double(datk.meteo.P(:));
    if isfield(datk.meteo,'Ep')
        Ep = double(datk.meteo.Ep(:));
    else
        Ep = nan(size(P));
    end

    n = numel(P);
    if n == 0
        axis(axh,'off');
        return
    end
    if numel(Ep) < n
        Ep(end+1:n,1) = nan;
    else
        Ep = Ep(1:n);
    end
    x = (1:n).';

    idx = unique([id_train(:); id_eval(:)]);
    idx = idx(isfinite(idx) & idx >= 1 & idx <= n);
    if isempty(idx)
        idx = x;
    end
    x1 = idx(1);
    x2 = idx(end);

    xT = id_train(:).';
    xE = id_eval(:).';
    xT = xT(xT >= 1 & xT <= n);
    xE = xE(xE >= 1 & xE <= n);
    if isfield(datk,'bad') && ~isempty(datk.bad)
        bad = logical(datk.bad(:));
        if ~isempty(xT) && numel(bad) >= max(xT)
            xT = xT(~bad(xT));
        end
        if ~isempty(xE) && numel(bad) >= max(xE)
            xE = xE(~bad(xE));
        end
    end
    useShading = local_should_use_shading( ...
        sp_method,xT,xE,30);

    cP = [0 0 1];
    cEp = [0.46 0.55 0.89];

    yyaxis(axh,'left');
    bar(axh,x,P,1, ...
        'FaceColor',cP, ...
        'EdgeColor',cP, ...
        'LineWidth',0.15, ...
        'HandleVisibility','off');
    local_set_nonnegative_ylim(axh,P(idx));
    ylP = ylim(axh);
    if useShading
        hShade = [ ...
            local_plot_forcing_shading(axh,xT,ylP, ...
            cShadeTrain,face_alpha); ...
            local_plot_forcing_shading(axh,xE,ylP, ...
            cShadeEval,face_alpha)];
        if ~isempty(hShade)
            uistack(hShade(isgraphics(hShade)),'bottom');
        end
    end
    hY = ylabel(axh,sprintf('$P$ (%s)',forcingUnit), ...
        'Interpreter','latex', ...
        'FontSize',18, ...
        'Color',cP);
    local_position_left_ylabel(hY);
    axh.YAxis(1).Color = cP;

    yyaxis(axh,'right');
    plot(axh,x,Ep, ...
        'Color',cEp, ...
        'LineWidth',0.5, ...
        'HandleVisibility','off');
    ylabel(axh,sprintf('$E_{\\rm p}$ (%s)',forcingUnit), ...
        'Interpreter','latex', ...
        'FontSize',18, ...
        'Color',cEp);
    axh.YAxis(2).Color = cEp;
    local_set_nonnegative_ylim(axh,Ep(idx));
    axh.YAxis(2).Direction = 'reverse';

    xlim(axh,[x1-0.5 x2+0.5]);
    xlabel(axh,'');
    grid(axh,'off');
    box(axh,'on');
    set(axh,'Layer','top');
end

function local_plot_temperature_strip( ...
    axh,datk,id_train,id_eval,sp_method, ...
    cShadeTrain,cShadeEval,faceAlpha,qXLim,axisLineWidth)
%LOCAL_PLOT_TEMPERATURE_STRIP Mark time steps below freezing.

    hold(axh,'on');
    xlim(axh,qXLim);
    ylim(axh,[0 1]);

    xT = id_train(:).';
    xE = id_eval(:).';
    if isfield(datk,'bad') && ~isempty(datk.bad)
        bad = logical(datk.bad(:));
        if ~isempty(xT) && numel(bad) >= max(xT)
            xT = xT(~bad(xT));
        end
        if ~isempty(xE) && numel(bad) >= max(xE)
            xE = xE(~bad(xE));
        end
    end

    useShading = local_should_use_shading( ...
        sp_method,xT,xE,30);
    if useShading
        local_plot_forcing_shading( ...
            axh,xT,[0 1],cShadeTrain,faceAlpha);
        local_plot_forcing_shading( ...
            axh,xE,[0 1],cShadeEval,faceAlpha);
    end

    if isfield(datk,'meteo') ...
            && isstruct(datk.meteo) ...
            && isfield(datk.meteo,'T')
        T = double(datk.meteo.T(:));
        idCold = find(isfinite(T) & T < 0);
        idCold = idCold(idCold >= qXLim(1) ...
            & idCold <= qXLim(2));
        if ~isempty(idCold)
            coldColor = [0 0 0];
            bar(axh,idCold,ones(size(idCold)),1, ...
                'FaceColor',coldColor, ...
                'EdgeColor',coldColor, ...
                'LineWidth',0.10, ...
                'HandleVisibility','off');
        end
    end

    set(axh, ...
        'XLim',qXLim, ...
        'YLim',[0 1], ...
        'XTick',[], ...
        'YTick',[], ...
        'XColor','k', ...
        'YColor','k', ...
        'LineWidth',axisLineWidth, ...
        'Box','on', ...
        'Layer','top');

    hY = ylabel(axh,'$T<0$', ...
        'Interpreter','latex', ...
        'Rotation',0, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontSize',14, ...
        'Color','k');
    hY.Units = 'normalized';
    hY.Position = [-0.020 0.5 0];
end

function [pF,pT,pQ] = local_timeseries_axes_positions( ...
    nRows,nCols,index)
%LOCAL_TIMESERIES_AXES_POSITIONS Fixed normalized time-series geometry.

    if nCols == 1
        leftMargin = 0.090;
        rightMargin = 0.105;
        columnGap = 0;
    else
        leftMargin = 0.060;
        rightMargin = 0.060;
        columnGap = 0.110;
    end

    topMargin = 0.075;
    bottomMargin = 0.105;
    rowGap = 0.115*double(nRows > 1);
    panelWidth = (1-leftMargin-rightMargin ...
        -(nCols-1)*columnGap)/nCols;
    groupHeight = (1-topMargin-bottomMargin ...
        -(nRows-1)*rowGap)/nRows;

    row = floor((index-1)/nCols);
    column = mod(index-1,nCols);
    x = leftMargin + column*(panelWidth+columnGap);
    y = 1-topMargin-(row+1)*groupHeight-row*rowGap;

    % Reduce the equal clearances above and below the temperature strip
    % to 60% of their former size. Split the recovered height equally
    % between the forcing and discharge panels.
    forcingHeight = 0.322*groupHeight;
    dischargeHeight = 0.582*groupHeight;
    interPanelGap = groupHeight-forcingHeight-dischargeHeight;
    temperatureHeight = 0.030*groupHeight;
    temperatureBottom = y+dischargeHeight ...
        + 0.5*(interPanelGap-temperatureHeight);

    pQ = [x y panelWidth dischargeHeight];
    pT = [x temperatureBottom panelWidth temperatureHeight];
    pF = [x y+dischargeHeight+interPanelGap ...
        panelWidth forcingHeight];
end

function local_timeseries_figure_title(f,titleText)
%LOCAL_TIMESERIES_FIGURE_TITLE Add a stable figure-level title.

    annotation(f,'textbox',[0.02 0.955 0.96 0.035], ...
        'String',titleText, ...
        'Interpreter','latex', ...
        'FontSize',17, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'EdgeColor','none', ...
        'Margin',0);
end

function handles = local_plot_forcing_shading( ...
    axh,id,yl,colorValue,faceAlpha)
%LOCAL_PLOT_FORCING_SHADING Shade contiguous valid forcing intervals.

    runs = local_contiguous_runs(id);
    handles = gobjects(size(runs,1),1);
    for i = 1:size(runs,1)
        x1 = runs(i,1) - 0.5;
        x2 = runs(i,2) + 0.5;
        handles(i) = patch(axh,[x1 x2 x2 x1], ...
            [yl(1) yl(1) yl(2) yl(2)], ...
            colorValue, ...
            'FaceAlpha',faceAlpha, ...
            'EdgeColor','none', ...
            'HandleVisibility','off');
    end
end

function local_position_left_ylabel(hY,desiredX)
%LOCAL_POSITION_LEFT_YLABEL Align forcing and discharge label centers.

    if isempty(hY) || ~isgraphics(hY)
        return
    end
    if nargin < 2 || isempty(desiredX)
        desiredX = -0.0315;
    end
    hY.Units = 'normalized';
    pos = hY.Position;
    minFigureX = 0.012;

    ax = ancestor(hY,'axes');
    if ~isempty(ax) && isgraphics(ax)
        oldUnits = ax.Units;
        cleanup = onCleanup(@()set(ax,'Units',oldUnits));
        ax.Units = 'normalized';
        axPos = ax.Position;
        minimumX = (minFigureX-axPos(1)) ...
            / max(axPos(3),eps);
        pos(1) = max(desiredX,minimumX);
        clear cleanup
    else
        pos(1) = desiredX;
    end
    pos(2) = 0.5;
    hY.Position = pos;
end

function local_position_single_basin_ylabels(axForcing,axDischarge)
%LOCAL_POSITION_SINGLE_BASIN_YLABELS Move full-width labels slightly in.

    desiredX = -0.0175;
    if isgraphics(axForcing)
        try
            local_position_left_ylabel( ...
                axForcing.YAxis(1).Label,desiredX);
        catch
        end
    end
    if isgraphics(axDischarge)
        local_position_left_ylabel( ...
            axDischarge.YLabel,desiredX);
    end
end

function local_set_nonnegative_ylim(axh,y)
%LOCAL_SET_NONNEGATIVE_YLIM Apply a compact nonnegative forcing range.

    y = y(isfinite(y));
    if isempty(y)
        ylim(axh,[0 1]);
        return
    end

    ymax = max(y);
    if ymax <= 0
        ymax = 1;
    end
    ylim(axh,[0 1.05*ymax]);
end

function useShading = local_should_use_shading( ...
    sp_method,xTrain,xEval,minRunForShade)
%LOCAL_SHOULD_USE_SHADING Tests whether forcing should be rendered as shading.

    if nargin < 4 ...
            || isempty(minRunForShade)
        minRunForShade = 30;
    end
    
    sm = lower(string(sp_method));

    % Manual-date scenarios already separate training and evaluation
    % periods into distinct panels, so background shading is redundant.
    if strcmp(sm,'manual') ...
            || strcmp(sm,'manual_dates')
        useShading = false;
        return
    end
    
    % pure random sampling: always use points
    if strcmp(sm,'random') ...
            || strcmp(sm,'random_sampling')
        useShading = false;
        return
    end
    
    maxRunTrain = ...
        local_max_contiguous_run_length(xTrain);
    maxRunEval = ...
        local_max_contiguous_run_length(xEval);
    
    useShading = (maxRunTrain >= minRunForShade) ...
              || (maxRunEval  >= minRunForShade);
end

function L = local_max_contiguous_run_length(id)
%LOCAL_MAX_CONTIGUOUS_RUN_LENGTH Returns the longest consecutive index run.

    id = unique(double(id(:)).');
    id = id(isfinite(id));
    
    if isempty(id)
        L = 0;
        return
    end
    
    runs = local_contiguous_runs(id);
    if isempty(runs)
        L = 0;
    else
        L = max(runs(:,2) - runs(:,1) + 1);
    end
end

function [hLeg,lblLeg] = ...
    local_make_postproc_legend(axh, ...
    n_m,id_model,colors, ...
    cPtTrain,cPtEval,model_names)
% LOCAL_MAKE_POSTPROC_LEGEND

    hLeg = gobjects(0);
    lblLeg = {};
    
    hold(axh,'on');
    
    % model proxies
    for ii = 1:n_m
        c = colors(id_model(min(ii, ...
            numel(id_model))),:);
    
        h = plot(axh,nan,nan,'-', ...
            'color',c, ...
            'linewidth',2.5);
        hLeg(end+1) = h;            %#ok
    
        name_i = local_model_display_name( ...
            model_names{id_model(min(ii, ...
            numel(id_model)))});
    
        lblLeg{end+1} = sprintf('$\\ \\mathrm{%s}$', ...
            local_latex_escape(name_i)); %#ok
    end
    
    % training-data proxy
    h = plot(axh,nan,nan,'o', ...
        'MarkerFaceColor',cPtTrain, ...
        'MarkerEdgeColor',cPtTrain, ...
        'markersize',6, ...
        'linewidth',1.4, ...
        'linestyle','none');
    hLeg(end+1) = h;
    lblLeg{end+1} = '$\ \mathrm{training\; data}$';
    
    % evaluation-data proxy
    h = plot(axh,nan,nan,'o', ...
        'MarkerFaceColor',cPtEval, ...
        'MarkerEdgeColor',cPtEval, ...
        'markersize',6, ...
        'linewidth',1.4, ...
        'linestyle','none');
    hLeg(end+1) = h;
    lblLeg{end+1} = '$\ \mathrm{evaluation\; data}$';

end

function codes = local_collect_group_codes( ...
    gaugescen,groupName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL_COLLECT_GROUP_CODES Collect requested basin codes for group plots
%
% SYNOPSIS:
%   codes = local_collect_group_codes(gaugescen,groupName)
%   gaugescen     structure with optional basin selections
%    .tt         requested basins, training basins|training period
%    .te         requested basins, training basins|evaluation period/mask
%    .et         requested basins, evaluation basins|training period
%    .ee         requested basins, evaluation basins|evaluation period/mask
%    .train      requested basins, all training-basin plots
%    .eval       requested basins, all evaluation-basin plots
%   groupName   character string with group name
%     'train'    collect codes for training-basin plots
%     'eval'     collect codes for evaluation-basin plots
%   codes       OUTPUT: unique basin codes in stable order
%
% NOTES:
%   1. This function uses one notation only
%        tt = training basins | training period
%        te = training basins | evaluation period/mask
%        et = evaluation basins | training period
%        ee = evaluation basins | evaluation period/mask
%   2. Field .train applies to all training-basin plots
%   3. Field .eval  applies to all evaluation-basin plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    codes = [];
    
    if nargin < 2 ...
            || isempty(gaugescen) ...
            || ~isstruct(gaugescen)
        return
    end
    
    switch lower(groupName)
        case 'train'
            C = {};
    
            if isfield(gaugescen,'train') ...
                    && ~isempty(gaugescen.train)
                C = [C,cellstr(string(gaugescen.train))];
            end
            if isfield(gaugescen,'tt') ...
                    && ~isempty(gaugescen.tt)
                C = [C,cellstr(string(gaugescen.tt))];
            end
            if isfield(gaugescen,'te') ...
                    && ~isempty(gaugescen.te)
                C = [C,cellstr(string(gaugescen.te))];
            end
            codes = unique(string(C(:)),'stable');
    
        case 'eval'
            C = {};
    
            if isfield(gaugescen,'eval') ...
                    && ~isempty(gaugescen.eval)
                C = [C,cellstr(string(gaugescen.eval))];
            end
            if isfield(gaugescen,'et') ...
                    && ~isempty(gaugescen.et)
                C = [C,cellstr(string(gaugescen.et))];
            end
            if isfield(gaugescen,'ee') ...
                    && ~isempty(gaugescen.ee)
                C = [C,cellstr(string(gaugescen.ee))];
            end
            codes = unique(string(C(:)),'stable');
    
        otherwise
            codes = [];
    end

end

function tf = ...
    local_use_four_scenario_timeseries(sp_method)
%LOCAL_USE_FOUR_SCENARIO_TIMESERIES Tests whether four-scenario time-series layout is required.

    sm = string(sp_method);
    
    tf = any(strcmpi(sm,[ ...
        "manual", ...
        "manual_dates", ...
        "traditional_block", ...
        "deterministic_block", ...
        "block"]));

end

function hfig = local_plot_paper_hydrograph(cfg,mdl,dat,bas,prd,Q, ...
    q_unit,t_unit,n_m,id_model,colors,model_names,id_dat,sp_method, ...
    cShadeTrain,cShadeEval,faceAlpha,region)
%LOCAL_PLOT_PAPER_HYDROGRAPH Compose selected paper hydrograph panels.
% The graphics are delegated to the same plot_SAGE helpers used by the GUI:
% forcing bars/PET, temperature strip, discharge, observation treatment,
% and date conversion. This function only selects basins/water-year ranges
% and arranges two requested scenarios side by side for the manuscript.

    required = {'panels','figure_number'};
    for i = 1:numel(required)
        if ~isfield(cfg,required{i})
            error('plot_SAGE:PaperHydrographConfig', ...
                'paper_hydrograph.%s is required.',required{i});
        end
    end
    panels = cfg.panels;
    nPanels = numel(panels);
    if nPanels < 1 || nPanels > 2
        error('plot_SAGE:PaperHydrographPanels', ...
            'Paper hydrographs currently support one or two panels.');
    end

    hfig = figure(double(cfg.figure_number));
    clf(hfig,'reset');
    set(hfig,'Color','w','Visible','on','NumberTitle','off', ...
        'Units','pixels','Position',[20 90 1400 800]);
    if isfield(cfg,'name') && ~isempty(cfg.name)
        hfig.Name = char(string(cfg.name));
    end

    axF = gobjects(nPanels,1);
    axT = gobjects(nPanels,1);
    axQ = gobjects(nPanels,1);

    for ip = 1:nPanels
        P = panels(ip);
        tag = lower(char(string(P.tag)));
        if ~ismember(tag,{'tt','te','et','ee'})
            error('plot_SAGE:PaperHydrographTag', ...
                'Unsupported paper hydrograph scenario: %s',tag);
        end
        if ~isfield(Q,tag) || isempty(Q.(tag))
            error('plot_SAGE:PaperHydrographData', ...
                'Discharge scenario Q.%s is unavailable.',tag);
        end

        if ismember(tag,{'tt','te'})
            idGlobal = 1:bas.K_t;
        else
            idGlobal = (bas.K_t+1):(bas.K_t+bas.K_e);
        end
        R = local_resolve_requested_basins( ...
            string(P.gauge),id_dat,idGlobal,region);
        if isempty(R.kdat)
            error('plot_SAGE:PaperHydrographGauge', ...
                'Gauge %s is not available for scenario %s.', ...
                string(P.gauge),tag);
        end
        kdat = R.kdat(1);
        qCol = R.icol(1);

        if ismember(tag,{'tt','et'})
            idxFull = expand_index(mdl.id_train);
        else
            idxFull = expand_index(mdl.id_eval);
        end
        idxFull = idxFull(:);
        if isfield(P,'water_years') && ~isempty(P.water_years)
            waterYears = double(P.water_years(:).');
        elseif isfield(P,'water_year') && ~isempty(P.water_year)
            waterYears = double(P.water_year);
        else
            error('plot_SAGE:PaperHydrographWaterYears', ...
                'Each panel requires water_year or water_years.');
        end
        [idxPlot,rowMask] = local_select_paper_water_years( ...
            tag,prd,idxFull,waterYears);
        if isempty(idxPlot)
            error('plot_SAGE:PaperHydrographWaterYear', ...
                'Water Years %d--%d are unavailable for scenario %s.', ...
                min(waterYears),max(waterYears),tag);
        end

        Qscenario = Q.(tag);
        if size(Qscenario,1) == numel(idxFull)
            Qpanel = Qscenario(rowMask,:,:);
        elseif size(Qscenario,1) >= max(idxPlot)
            Qpanel = Qscenario(idxPlot,:,:);
        else
            error('plot_SAGE:PaperHydrographLength', ...
                'Q.%s cannot be aligned with its time indices.',tag);
        end

        [pF,pT,pQ] = local_timeseries_axes_positions(1,nPanels,ip);
        axF(ip) = axes(hfig,'Units','normalized','Position',pF, ...
            'PositionConstraint','innerposition');
        axT(ip) = axes(hfig,'Units','normalized','Position',pT, ...
            'PositionConstraint','innerposition');
        axQ(ip) = axes(hfig,'Units','normalized','Position',pQ, ...
            'PositionConstraint','innerposition');

        isTrainingData = local_is_training_data_scenario(tag);
        local_plot_basin_fourgroups(axQ(ip),dat{kdat},qCol,idxPlot, ...
            Qpanel,isTrainingData,n_m,id_model,colors,q_unit,t_unit, ...
            tag,prd,true,sp_method,cShadeTrain,cShadeEval,faceAlpha);

        if isTrainingData
            forcingTrain = idxPlot;
            forcingEval = [];
        else
            forcingTrain = [];
            forcingEval = idxPlot;
        end
        local_plot_forcing_postprocessor(axF(ip),dat{kdat}, ...
            forcingTrain,forcingEval,q_unit,sp_method, ...
            cShadeTrain,cShadeEval,faceAlpha);
        qXLim = xlim(axQ(ip));
        local_plot_temperature_strip(axT(ip),dat{kdat}, ...
            forcingTrain,forcingEval,sp_method,cShadeTrain,cShadeEval, ...
            faceAlpha,qXLim,axQ(ip).LineWidth);
        linkaxes([axF(ip),axT(ip),axQ(ip)],'x');
        xlim(axQ(ip),qXLim);
        local_set_multiyear_ticks_paper(axQ(ip),tag,prd,idxPlot);
        set(axF(ip),'XTick',get(axQ(ip),'XTick'),'XTickLabel',[]);
        grid(axF(ip),'off'); grid(axT(ip),'off'); grid(axQ(ip),'off');

        if isfield(P,'title') && ~isempty(P.title)
            panelTitle = char(string(P.title));
        else
            panelTitle = sprintf('USGS ID %s',string(P.gauge));
        end
        text(axF(ip),0.01,1.03,panelTitle,'Units','normalized', ...
            'FontName','Times New Roman','FontSize',18, ...
            'Interpreter','none','HorizontalAlignment','left', ...
            'VerticalAlignment','bottom','Clipping','off');
        if numel(waterYears) == 1
            periodLabel = sprintf('Water Year %d',waterYears(1));
        else
            periodLabel = sprintf('Water Years %d--%d', ...
                min(waterYears),max(waterYears));
        end
        xlabel(axQ(ip),periodLabel, ...
            'FontName','Times New Roman','FontSize',18, ...
            'Interpreter','none');
        % Keep paper labels clear of their tick values. Normalized offsets
        % remain stable when the discharge or forcing ranges change.
        axQ(ip).XLabel.Units = 'normalized';
        axQ(ip).XLabel.Position(2) = -0.155;
        % Place P and Q at 70% of their former horizontal offset so both
        % labels sit slightly closer to their respective left y-axes.
        leftLabelX = 0.70 * -0.053;
        rightLabelX = 1.035;
        axQ(ip).YLabel.Units = 'normalized';
        axQ(ip).YLabel.FontSize = 18;
        qLabelPos = axQ(ip).YLabel.Position;
        qLabelPos(1:2) = [leftLabelX 0.5];
        axQ(ip).YLabel.Position = qLabelPos;
        try
            axF(ip).YAxis(1).Label.Units = 'normalized';
            axF(ip).YAxis(1).Label.FontSize = 18;
            leftPos = axF(ip).YAxis(1).Label.Position;
            leftPos(1:2) = [leftLabelX 0.5];
            axF(ip).YAxis(1).Label.Position = leftPos;
            axF(ip).YAxis(2).Label.Units = 'normalized';
            axF(ip).YAxis(2).Label.FontSize = 18;
            rightPos = axF(ip).YAxis(2).Label.Position;
            rightPos(1:2) = [rightLabelX 0.5];
            axF(ip).YAxis(2).Label.Position = rightPos;
        catch
            % Older MATLAB releases can expose yyaxis labels differently.
        end
    end

    set([axF;axQ],'FontName','Times New Roman','FontSize',15, ...
        'TickDir','out','Layer','top');
    for ia = 1:numel(axQ)
        % Bring month names closer to the discharge axis while leaving the
        % water-year label farther below them.
        axQ(ia).XRuler.TickLabelGapOffset = -2;
    end
    % Use outward ticks that are 1.5 times MATLAB's current default length.
    paperAxes = [axF;axT;axQ];
    for ia = 1:numel(paperAxes)
        paperAxes(ia).TickDir = 'out';
        paperAxes(ia).TickLength = 0.75*paperAxes(ia).TickLength;
    end
    set(axT,'XTick',[],'YTick',[]);
    local_add_paper_timeseries_legend(axQ(1),n_m,id_model, ...
        colors,model_names);

    if isfield(cfg,'save_pdf') && logical(cfg.save_pdf) ...
            && isfield(cfg,'output_pdf') && ~isempty(cfg.output_pdf)
        fprintf('Saving %s\n',char(string(cfg.output_pdf)));
        exportgraphics(hfig,char(string(cfg.output_pdf)), ...
            'ContentType','vector');
    end
end

function [idxPlot,rowMask] = local_select_paper_water_years( ...
    tag,prd,idxFull,waterYears)
%LOCAL_SELECT_PAPER_WATER_YEARS Select a contiguous range of water years.

    idxPlot = [];
    rowMask = false(size(idxFull));
    if isempty(idxFull)
        return
    end
    dtDays = local_time_step_days(prd);
    t = local_index_to_datetime_fourgroups( ...
        tag,prd,idxFull(1),idxFull,dtDays);
    if all(isnat(t))
        return
    end
    waterYears = sort(waterYears(:));
    tStart = datetime(waterYears(1)-1,10,1);
    tEnd = datetime(waterYears(end),10,1);
    rowMask = t >= tStart & t < tEnd;
    idxPlot = idxFull(rowMask);
end

function local_set_multiyear_ticks_paper(ax,tag,prd,idx)
%LOCAL_SET_MULTIYEAR_TICKS_PAPER Label multi-year plots at six-month steps.

    if isempty(idx)
        return
    end
    dtDays = local_time_step_days(prd);
    t = local_index_to_datetime_fourgroups( ...
        tag,prd,idx(1),idx,dtDays);
    months = dateshift(t(1),'start','month'):calmonths(1): ...
        dateshift(t(end),'start','month');
    xt = nan(size(months));
    for i = 1:numel(months)
        [~,j] = min(abs(t-months(i)));
        xt(i) = idx(j);
    end
    labels = cellstr(string(months,'MMM'));
    set(ax,'XTick',xt,'XTickLabel',labels,'XTickLabelRotation',0);
end

function local_add_paper_timeseries_legend(ax,n_m,id_model, ...
    colors,model_names)
%LOCAL_ADD_PAPER_TIMESERIES_LEGEND Manual legend with colored text.

    xl = xlim(ax); yl = ylim(ax);
    % Place the manual legend beside the left y-axis. Keeping it away from
    % the right side avoids the principal hydrograph peaks and observations.
    x1 = xl(1) + 0.138*(xl(2)-xl(1));
    dx = 0.04*(xl(2)-xl(1));
    y1 = yl(1) + 0.92*(yl(2)-yl(1));
    dy = 0.085*(yl(2)-yl(1));
    cObs = [0.25 0.25 0.25];
    % Keep hydrograph peaks and observation points from crossing the
    % legend. The nearly opaque backing is restricted to the legend area.
    legendTop = y1 + 0.045*(yl(2)-yl(1));
    legendBottom = y1 - (n_m+0.55)*dy;
    patch(ax,[x1-0.018*(xl(2)-xl(1)) ...
              x1+0.245*(xl(2)-xl(1)) ...
              x1+0.245*(xl(2)-xl(1)) ...
              x1-0.018*(xl(2)-xl(1))], ...
             [legendBottom legendBottom legendTop legendTop], ...
             'w','EdgeColor','none','FaceAlpha',0.94, ...
             'HandleVisibility','off');
    plot(ax,x1+0.5*dx,y1,'o','MarkerSize',5, ...
        'MarkerFaceColor',cObs,'MarkerEdgeColor',cObs, ...
        'HandleVisibility','off');
    text(ax,x1+dx,y1,' Observed','Color',cObs, ...
        'FontName','Times New Roman','FontSize',14, ...
        'FontWeight','bold', ...
        'VerticalAlignment','middle');
    for i = 1:n_m
        y = y1-i*dy;
        im = id_model(min(i,numel(id_model)));
        c = colors(im,:);
        line(ax,[x1 x1+dx],[y y],'Color',c,'LineWidth',2.2, ...
            'HandleVisibility','off');
        label = local_model_display_name(model_names{im});
        text(ax,x1+dx,y,[' ' label],'Color',c, ...
            'FontName','Times New Roman','FontSize',14, ...
            'FontWeight','bold', ...
            'VerticalAlignment','middle');
    end
end

function local_plot_timeseries_fourgroups( ...
    figBase,part,mdl,dat,bas,prd,scenarios, ...
    q_unit,t_unit,dt_str,n_m,id_model, ...
    colors,model_names,id_dat,Qfdc,Bprt, ...
    n_char,figW,figH,model_name_fig,samp_word, ...
    fontsize_legend,region,sp_method, ...
    cShadeTrain,cShadeEval,face_alpha)
%LOCAL_PLOT_TIMESERIES_FOURGROUPS Plot four explicit basin/period groups
% SYNOPSIS:
%   local_plot_timeseries_fourgroups( ...
%       figBase,part,mdl,dat,bas,prd,scenarios, ...
%       q_unit,t_unit,dt_str,n_m,id_model,colors,model_names, ...
%       id_dat,Qfdc,Bprt,n_char,figW,figH,model_name_fig,samp_word, ...
%       fontsize_legend,region)
% This helper is used for train/eval methods with separable training and
% evaluation windows, such as:
%   - manual
%   - manual_dates
%   - traditional_block
%   - block
% The four explicit scenario groups are:
%   1) training basins   | training period
%   2) training basins   | evaluation period
%   3) evaluation basins | training period
%   4) evaluation basins | evaluation period
% Each figure uses the original 2x2 layout.

    figW = 0.95*figW;   % was 16
    figH = 0.95*figH;   % was 8
    
    allSelected = false;
    for is = 1:4
        allSelected = allSelected ...
            || ~isempty(scenarios{is}.codes);
    end
    
    useDateAxis = local_use_dates_on_xaxis(prd, ...
        mdl);
    
    for is = 1:4
        S = scenarios{is};
    
        % fprintf('\n--- %s ---\n',S.tag);
        % fprintf('S.codes = ');
        % disp(S.codes);

        % if isfield(S,'basin_idx') ...
        %         && ~isempty(S.basin_idx)
        %     fprintf('S.basin_idx = ');
        %     disp(S.basin_idx(:).');
        % elseif ~isempty(S.codes)
        %     Rtest = local_resolve_requested_basins( ...
        %         S.codes,id_dat,S.id_global,region);
        %     fprintf('resolved kdat = ');
        %     disp(Rtest.kdat.');
        %     fprintf('resolved qCols = ');
        %     disp(Rtest.icol.');
        % else
        %     fprintf(['S.codes is EMPTY -> code ' ...
        %         'will fall back to first basins.\n']);
        % end

        if isempty(S.id_global) ...
                || ~local_has_data(S.Q)
            continue
        end
            
        if isfield(S,'basin_idx') ...
                && ~isempty(S.basin_idx)

            pick = S.basin_idx(:).';
            pick = pick(ismember(pick,S.id_global));

            qCols = nan(size(pick));
            for jj = 1:numel(pick)
                qCols(jj) = find(S.id_global(:) == pick(jj),1,'first');
            end

        elseif ~isempty(S.codes)

            R = local_resolve_requested_basins( ...
                S.codes,id_dat,S.id_global,region);

            if ~isempty(R.outside)
                warning('plot_SAGE:requestedBasinsOutsideScenario', ...
                    ['Requested basin(s) exist in ' ...
                    'dat but are not in %s: %s. ', ...
                    'They will be skipped.'], ...
                    S.tag,strjoin(R.outside,', '));
            end

            if ~isempty(R.missing)
                warning('plot_SAGE:missingRequestedBasins', ...
                    ['Requested basin(s) do not ' ...
                    'exist in dat for %s: %s. ', ...
                    'They will be skipped.'], ...
                    S.tag,strjoin(R.missing,', '));
            end

            if isempty(R.kdat)
                continue
            end

            pick = R.kdat;
            qCols = R.icol;

        else
            % No explicit gaugescen request was supplied.  Before random Bprt
            % selection, reuse FDC-selected dat{k} indices when Qfdc.id exists.
            % If Qfdc.req has scenario-specific entries, use them first;
            % otherwise filter Qfdc.id by this scenario's global basin IDs.
            pick = local_pick_qfdc_ids(Qfdc,S.id_global,S.tag);
            if isempty(pick)
                pick = pick_basins([], ...
                    id_dat,S.id_global, ...
                    min(Bprt,numel(S.id_global)),region);
            end
            qCols = nan(size(pick));
            for jj = 1:numel(pick)
                qCols(jj) = find(S.id_global(:) ...
                    == pick(jj),1,'first');
            end
        end
    
        if isempty(pick)
            continue
        end
    
        nPerFig = 4;
        nFig = ceil(numel(pick)/nPerFig);
    
        for jf = 1:nFig
            f = figure(figBase);
            clf(f,'reset');
            figBase = figBase + 1;
    
            nThisFig = min(nPerFig, ...
                numel(pick) - (jf-1)*nPerFig);
            if nThisFig <= 2
                nTileRows = 1;
                nTileCols = max(1,nThisFig);
            else
                nTileRows = 2;
                nTileCols = 2;
            end
            figH_page = figH * (0.56 + 0.44*(nTileRows == 2));
    
            set(f, ...
                'Name',sprintf('%s: Time series - %s (%d/%d)', ...
                model_name_fig,S.title,jf,nFig), ...
                'NumberTitle','off', ...
                'color','w', ...
                'Units','inches', ...
                'Position',[0.3 0.3 figW figH_page], ...
                'Visible','on');
    
            axh = gobjects(nThisFig,1);
            axForcing = gobjects(nThisFig,1);
            axTemperature = gobjects(nThisFig,1);
    
            for jp = 1:nThisFig
                kk = (jf-1)*nPerFig + jp;
                [pF,pT,pQ] = local_timeseries_axes_positions( ...
                    nTileRows,nTileCols,jp);
                axForcing(jp) = axes(f, ...
                    'Units','normalized','Position',pF, ...
                    'PositionConstraint','innerposition');
                axTemperature(jp) = axes(f, ...
                    'Units','normalized','Position',pT, ...
                    'PositionConstraint','innerposition');
                axh(jp) = axes(f, ...
                    'Units','normalized','Position',pQ, ...
                    'PositionConstraint','innerposition');
                box(axh(jp),'on');
                hold(axh(jp),'on');
                set(axh(jp),'Layer','top');
    
                if kk > numel(pick)
                    axis(axh(jp),'off');
                    continue
                end
    
                kdat = pick(kk);
                qCol = qCols(kk);
    
                if ndims(S.Q) == 3
                    if qCol < 1 ...
                            || qCol > size(S.Q,2)
                        axis(axh(jp),'off');
                        text(axh(jp),0.5,0.5, ...
                            'qCol out of range', ...
                            'Units','normalized', ...
                            'horizontalalignment', ...
                            'center');
                        continue
                    end
                else
                    if qCol < 1 ...
                            || qCol > size(S.Q,2)
                        axis(axh(jp),'off');
                        text(axh(jp),0.5,0.5, ...
                            'qCol out of range', ...
                            'Units','normalized', ...
                            'horizontalalignment', ...
                            'center');
                        continue
                    end
                end
    
                isTrainData = ...
                    local_is_training_data_scenario(S.tag);
    
                idxPlot = S.idx;
                if isempty(idxPlot)
                    try
                        iq = find(double(Qfdc.id(:)) ...
                            == double(kdat),1,'first');
                        if ~isempty(iq)
                            idxPlot = ...
                                local_fdc_indices_for_basin( ...
                                Qfdc,dat,kdat,iq,S.tag,mdl, ...
                                numel(dat{kdat}.y_n));
                        end
                    catch
                        idxPlot = [];
                    end
                end
                
                if isempty(idxPlot)
                    axis(axh(jp),'off');
                    text(axh(jp),0.5,0.5,'No time indices', ...
                        'Units','normalized', ...
                        'HorizontalAlignment','center');
                    continue
                end
    
                local_plot_basin_fourgroups(axh(jp), ...
                    dat{kdat},qCol,idxPlot,S.Q, ...
                    isTrainData,n_m,id_model,colors, ...
                    q_unit,t_unit,S.tag, ...
                    prd,useDateAxis,sp_method,cShadeTrain, ...
                    cShadeEval,face_alpha);

                % The discharge scenario contains the authoritative
                % scored-window limits, including exclusion of spin-up.
                qXLim = xlim(axh(jp));
                qXTick = get(axh(jp),'XTick');

                if local_is_training_data_scenario(S.tag)
                    forcingTrain = idxPlot;
                    forcingEval = [];
                else
                    forcingTrain = [];
                    forcingEval = idxPlot;
                end
                local_plot_forcing_postprocessor( ...
                    axForcing(jp),dat{kdat}, ...
                    forcingTrain,forcingEval,q_unit, ...
                    sp_method,cShadeTrain,cShadeEval, ...
                    face_alpha);
                local_plot_temperature_strip( ...
                    axTemperature(jp),dat{kdat}, ...
                    forcingTrain,forcingEval,sp_method, ...
                    cShadeTrain,cShadeEval,face_alpha, ...
                    qXLim,axh(jp).LineWidth);
                xlim(axForcing(jp),qXLim);
                linkaxes([axForcing(jp), ...
                    axTemperature(jp),axh(jp)],'x');
                xlim(axh(jp),qXLim);
                set(axForcing(jp), ...
                    'XTick',qXTick, ...
                    'XTickLabel',[]);

                gauge = local_gauge_name(bas,kdat);
                us = local_get_usgs(id_dat,kdat);
                zoneTxt = local_zone_label_short( ...
                    bas,kdat);
                invalidForcing = local_incomplete_forcing(dat{kdat});
                labelColor = [0.15 0.15 0.15];
                if invalidForcing
                    labelColor = [0.80 0 0];
                end
                add_basin_label(axForcing(jp),us, ...
                    zoneTxt,gauge,n_char,0.5,labelColor);
                if invalidForcing
                    local_add_invalid_forcing_notice(axh(jp));
                end
            end


            drawnow limitrate nocallbacks
    
            goodAx = axh(isgraphics(axh));
            set(goodAx, ...
                'tickdir','out', ...
                'TickLength',[0.004 0.004], ...
                'fontsize',12);
            goodForcing = axForcing(isgraphics(axForcing));
            set(goodForcing, ...
                'tickdir','out', ...
                'TickLength',[0.004 0.004], ...
                'fontsize',12);
            goodTemperature = axTemperature( ...
                isgraphics(axTemperature));
            set(goodTemperature, ...
                'XTick',[], ...
                'YTick',[]);
            axLeg = [];
            if ~isempty(goodAx)
                vis = get(goodAx,'Visible');
                if iscell(vis)
                    iv = find(strcmp(vis,'on'), ...
                        1,'last');
                else
                    iv = strcmp(vis,'on');
                    iv = find(iv,1,'last');
                end
    
                if ~isempty(iv)
                    axLeg = goodAx(iv);
                else
                    axLeg = goodAx(1);
                end
            end
    
            if ~isempty(axLeg)
                [hLeg,lblLeg] = ...
                    local_make_fourgroup_legend( ...
                    axLeg,S.tag,n_m,id_model, ...
                    colors,model_names);
    
                legend(axLeg,hLeg,lblLeg, ...
                    'Location','northeast', ...
                    'interpreter','latex', ...
                    'fontsize',fontsize_legend, ...
                    'box','off');
            end
    
            if nThisFig == 1
                local_position_single_basin_ylabels( ...
                    axForcing(1),axh(1));
            end

            model_name_tex = local_latex_escape(model_name_fig);
            local_timeseries_figure_title(f,sprintf( ...
                '\\texttt{%s}: %s --- figure %d/%d', ...
                model_name_tex,S.title,jf,nFig));
            
            drawnow limitrate nocallbacks
            figure(f);
            drawnow limitrate nocallbacks
        end
    end

end

function local_plot_basin_fourgroups(axh,datk,qCol,idx,Qmat, ...
    isTrainData,n_m,id_model,colors,q_unit,t_unit,tag,prd, ...
    useDateAxis,sp_method,cShadeTrain,cShadeEval,faceAlpha)
%LOCAL_PLOT_BASIN_FOURGROUPS Plot one basin for one explicit scenario
%
% Observation color depends on whether the plotted observations belong to
% the training period or the evaluation period.

    fontsize_axis = 18;
    lw_sim = 1.10;
    
    cPtTrain = [0.25 0.25 0.25];
    cPtEval = [0.55 0.55 0.55];
    
    %cObsEdge = [0 0 0];
    cObsEdge = 'none'; %#ok
    msObs = 2.5;
    
    if isTrainData
        colObs = cPtTrain;
    else
        colObs = cPtEval;
    end
    
    idx = idx(:);
    y = datk.y_n(:);
    bad = datk.bad(:);
    
    if isempty(idx)
        axis(axh,'off');
        return
    end
    
    if numel(bad) < max(idx)
        good_obs = true(size(idx));
    else
        good_obs = ~bad(idx);
    end
    
    idx_obs = idx(good_obs);
    y_obs = y(idx_obs);
    
    if ~isempty(idx_obs)
        plot(axh,idx_obs,y_obs,'o', ...
            'color',colObs, ...
            'MarkerFaceColor',colObs, ...
            'MarkerEdgeColor',colObs, ...
            'markersize',msObs, ...
            'linestyle','none', ...
            'handlevisibility','off');
        hold(axh,'on');
    end
    
    idx_sim = idx(:);
    yminSim = inf;
    ymaxSim = -inf;
    
    for mdl_i = 1:n_m
        q = local_get_discharge(Qmat, ...
            qCol,mdl_i);
        q = q(:);
    
        if numel(q) == numel(idx_sim)
            x_sim = idx_sim;
            q_sim = q;
    
        elseif numel(q) >= max(idx_sim)
            x_sim = idx_sim;
            q_sim = q(idx_sim);
    
        else
            n = min(numel(q), ...
                numel(idx_sim));
            x_sim = idx_sim(1:n);
            q_sim = q(1:n);
        end
    
        good_sim = isfinite(q_sim) ...
            & isreal(q_sim);
    
        x_sim = x_sim(good_sim);
        q_sim = q_sim(good_sim);

        if ~isempty(q_sim)
            yminSim = min(yminSim,min(q_sim));
            ymaxSim = max(ymaxSim,max(q_sim));
        end
    
        if ~isempty(x_sim)
            im = id_model(min(mdl_i, ...
                numel(id_model)));
            plot(axh,x_sim,q_sim, ...
                'linewidth',lw_sim, ...
                'color',colors(im,:), ...
                'handlevisibility','off');
        end
    end
    
    ymin = inf;
    ymax = -inf;
       
    % ----------------------------------------------
    % x-limits: always use full simulated/scenario
    % period. Observations are added where available
    % ----------------------------------------------
    xmin = min(idx);
    xmax = max(idx);

    % ----------------------------------------
    % y-limits: use everything that is plotted
    % ----------------------------------------
    if ~isempty(y_obs)
        ytmp = y_obs(isfinite(y_obs) ...
            & isreal(y_obs));
        if ~isempty(ytmp)
            ymin = min(ymin,min(ytmp));
            ymax = max(ymax,max(ytmp));
        end
    end

    if isfinite(yminSim) && isfinite(ymaxSim)
        ymin = min(ymin,yminSim);
        ymax = max(ymax,ymaxSim);
    end
    
    h = findobj(axh,'Type','line');
    for k = 1:numel(h)
        yd = double(h(k).YData(:));
        yd = yd(isfinite(yd) ...
            & isreal(yd));
        if ~isempty(yd)
            ymin = min(ymin,min(yd));
            ymax = max(ymax,max(yd));
        end
    end
    
    if ~isfinite(xmin) ...
        || ~isfinite(xmax)
        grid(axh,'off');
        return
    end
    
    if ~isfinite(ymin) ...
            || ~isfinite(ymax)
        ymin = 0;
        ymax = 1;
    end
        
    if xmax <= xmin
        xmax = xmin + 1;
    end
    
    if ymax <= ymin
        dy = max(1e-6,0.01*max(abs([ymin ymax])));
        ymin = ymin - dy;
        ymax = ymax + dy;
    end
    
    pad = 0.05;
    dy = ymax - ymin;
    ylimv = [ymin-pad*dy, ...
        ymax+pad*dy];
    
    ylimv(1) = min(0,ylimv(1));
    
    xlim(axh,[xmin xmax]);
    ylim(axh,ylimv);

    useShading = local_should_use_shading( ...
        sp_method,idx_obs,[],30);
    if useShading
        if isTrainData
            shadeColor = cShadeTrain;
        else
            shadeColor = cShadeEval;
        end
        hShade = local_plot_forcing_shading( ...
            axh,idx_obs,ylimv,shadeColor,faceAlpha);
        if ~isempty(hShade)
            try
                uistack(hShade(isgraphics(hShade)),'bottom');
            catch
                % Some graphics backends do not permit child reordering.
            end
        end
    end
    
    set(axh, ...
        'XLimMode','manual', ...
        'YLimMode','manual');
    
    if useDateAxis
        local_set_date_ticks_fourgroups(axh, ...
            tag,prd,idx,xmin,xmax);
        xlabel(axh,'', ...
            'interpreter','none', ...
            'fontsize',fontsize_axis);
    else
        local_set_integer_ticks(axh,xmin,xmax);
        xlabel(axh,sprintf('${\\rm time\\; (%s)}$', ...
            t_unit), ...
            'interpreter','latex', ...
            'fontsize',fontsize_axis);
    end
    hY = ylabel(axh,sprintf('$Q$ (%s)',q_unit), ...
        'interpreter','latex', ...
        'fontsize',fontsize_axis);
    local_position_left_ylabel(hY);
    
    grid(axh,'off');
    box(axh,'on');
    set(axh,'Layer','top');
end

function tf = ...
    local_is_training_data_scenario(tag)
%LOCAL_IS_TRAINING_DATA_SCENARIO Tests whether a scenario uses training data.

    tf = any(strcmpi(tag,{ ...
        'tt', ...
        'et'}));

end

function ttl = local_group_title_fourgroups( ...
    tag,part,mdl,dt_str,samp_word,kdat)
%LOCAL_GROUP_TITLE_FOURGROUPS Creates a title for a four-scenario panel group.

    mname = local_part_name(part,mdl);
    
    if strcmpi(part,'sage')
        tprefix = upper(char(mname));
    else
        tprefix = 'SITE';
    end
    
    switch lower(tag)
        case 'tt'
            ttl = sprintf(['%s %s training basin ' ...
                '| training %s (dat = %d)'], ...
                tprefix,dt_str,samp_word,kdat);
    
        case 'te'
            ttl = sprintf(['%s %s training basin ' ...
                '| evaluation %s (dat = %d)'], ...
                tprefix,dt_str,samp_word,kdat);
    
        case 'et'
            ttl = sprintf(['%s %s evaluation basin ' ...
                '| training %s (dat = %d)'], ...
                tprefix,dt_str,samp_word,kdat);
    
        case 'ee'
            ttl = sprintf(['%s %s evaluation basin ' ...
                '| evaluation %s (dat = %d)'], ...
                tprefix,dt_str,samp_word,kdat);
    
        otherwise
            ttl = char(string(tag));
    end
end

function [hLeg,lblLeg] = ...
    local_make_fourgroup_legend( ...
    axh,tag,n_m,id_model,colors, ...
    model_names)
%LOCAL_MAKE_FOURGROUP_LEGEND Builds handles and labels for four-scenario legends.

    cPtTrain = [0.25 0.25 0.25];
    cPtEval = [0.55 0.55 0.55];
    
    hLeg = gobjects(0);
    lblLeg = {};
    
    hold(axh,'on');
    
    for ii = 1:n_m
        c = colors(id_model( ...
            min(ii,numel(id_model))),:);
    
        h = plot(axh,nan,nan,'-', ...
            'color',c, ...
            'linewidth',2.5);
        hLeg(end+1) = h; %#ok
    
        name_i = local_model_display_name( ...
            model_names{id_model(min(ii,numel(id_model)))});
    
        lblLeg{end+1} = sprintf('$\\ \\mathrm{%s}$', ...
            local_latex_escape(name_i)); %#ok
    
    end
    
    if local_is_training_data_scenario(tag)
        h = plot(axh,nan,nan,'o', ...
            'MarkerFaceColor',cPtTrain, ...
            'MarkerEdgeColor',cPtTrain, ...
            'markersize',6, ...
            'linewidth',1.4, ...
            'linestyle','none');
        hLeg(end+1) = h;
        lblLeg{end+1} = '$\ \mathrm{training\; data}$';
    else
        h = plot(axh,nan,nan,'o', ...
            'MarkerFaceColor',cPtEval, ...
            'MarkerEdgeColor',cPtEval, ...
            'markersize',6, ...
            'linewidth',1.4, ...
            'linestyle','none');
        hLeg(end+1) = h;
        lblLeg{end+1} = '$\ \mathrm{evaluation\; data}$';
    end

end

function t = ...
    local_index_to_datetime_fourgroups(tag, ...
    prd,idx0,idx,dt_days)
%LOCAL_INDEX_TO_DATETIME_FOURGROUPS Convert plotted absolute index to date
%
% For rainfall_block the plotted x-values are absolute indices on the full
% model time axis. Therefore index 1 corresponds to prd.ds, index 2 to
% prd.ds + one time step, etc. Do not reset the clock at idx(1); otherwise
% spinup/scored-window indices such as 366 are labeled as prd.ds.

    t = NaT;
    
    try
        if local_is_rainfall_block_prd(prd)
            d0 = local_make_datetime(prd.ds);
            idx_ref = 1;
        else
            switch lower(tag)
                case {'tt','et'}
                    if isfield(prd,'dts') ...
                            && ~isempty(prd.dts)
                        d0 = ...
                            local_make_datetime(prd.dts);
                    else
                        d0 = ...
                            local_make_datetime(prd.ds);
                    end
                    idx_ref = idx0;
                
                case {'te','ee'}
                    if isfield(prd,'des') ...
                            && ~isempty(prd.des)
                        d0 = ...
                            local_make_datetime(prd.des);
                    elseif isfield(prd,'dts') ...
                            && ~isempty(prd.dts)
                        d0 = ...
                            local_make_datetime(prd.dts);
                    else
                        d0 = ...
                            local_make_datetime(prd.ds);
                    end
                    idx_ref = idx0;
            end
        end
    
        if isnat(d0)
            return
        end
    
        t = d0 + days((idx - idx_ref) * dt_days);
    
    catch
        t = NaT;
    end
end

function d = local_make_datetime(v)
%LOCAL_MAKE_DATETIME Convert [d m y] to datetime

    d = NaT;
    
    try
        if isempty(v) ...
                || numel(v) < 3
            return
        end
    
        d = datetime(double(v(3)), ...
            double(v(2)),double(v(1)));
    catch
        d = NaT;
    end

end

function dt_days = local_time_step_days(prd)
%LOCAL_TIME_STEP_DAYS Return time step in days

    dt_days = 1;
    
    try
        if isfield(prd,'dt') ...
                && ~isempty(prd.dt)
            if prd.dt == 24
                dt_days = 1/24;
            elseif prd.dt == 96
                dt_days = 1/96;
            else
                dt_days = 1;
            end
        end
    catch
        dt_days = 1;
    end

end

function s = local_datetime_label(t)
%LOCAL_DATETIME_LABEL Format datetime label for axis

    s = "";
    
    if isnat(t)
        return
    end
    
    try
        s = string(datetime(t, ...
            'Format','dd/MM/yyyy'));
    catch
        try
            s = string(t);
        catch
            s = "";
        end
    end

end

function local_set_date_ticks_fourgroups(axh, ...
    tag,prd,idx,xmin,xmax)
%LOCAL_SET_DATE_TICKS_FOURGROUPS Put horizontal date labels on x-axis
% and adapt their number to available axis width

    if nargin < 6 ...
            || ~isfinite(xmin) ...
            || ~isfinite(xmax) ...
            || isempty(idx)
        return
    end
    
    if xmax < xmin
        tmp = xmin;
        xmin = xmax;
        xmax = tmp;
    end
    
    idx0 = idx(1);
    dt_days = local_time_step_days(prd);
    
    t1 = local_index_to_datetime_fourgroups( ...
        tag,prd,idx0,xmin,dt_days);
    t2 = local_index_to_datetime_fourgroups( ...
        tag,prd,idx0,xmax,dt_days);
    
    if isnat(t1) ...
            || isnat(t2)
        return
    end
    
    if t2 < t1
        tmp = t1;
        t1 = t2;
        t2 = tmp;
    end
       
    max_labels = min(5, ...
        local_max_date_labels(axh,prd));
    max_labels = max(4, ...
        max_labels);

    % Water-year-end ticks only: 30 September.
    yrStart = year(t1);
    yrEnd = year(t2);

    % If t2 is after 30 Sep, last visible WY tick is current year.
    % If t2 is before 30 Sep, last visible WY tick is previous year.
    if t2 < datetime(yrEnd,9,30)
        yrEnd = yrEnd - 1;
    end

    if t1 > datetime(yrStart,9,30)
        yrStart = yrStart + 1;
    end

    if yrEnd < yrStart
        return
    end

    spanYears = yrEnd - yrStart;
    stepYears = max(1,ceil(spanYears/(max_labels-1)));

    yrs = yrEnd:-stepYears:yrStart;
    yrs = sort(yrs(:));

    tick_dates = datetime(yrs,9,30);

    % Guarantee final water-year date is printed.
    tEnd = datetime(yrEnd,9,30);
    if ~ismember(tEnd,tick_dates)
        tick_dates = [tick_dates; tEnd];
    end

    tick_dates = unique(sort(tick_dates),'stable');

    xt = zeros(numel(tick_dates),1);
    for i = 1:numel(tick_dates)
        xt(i) = local_datetime_to_index_fourgroups( ...
            tag,prd,idx0,tick_dates(i),dt_days);
    end

    good = isfinite(xt) ...
        & xt >= xmin ...
        & xt <= xmax;

    xt = xt(good);
    tick_dates = tick_dates(good);

    if isempty(xt)
        return
    end

    xt = round(xt(:));

    labs = cell(numel(xt),1);
    for i = 1:numel(xt)
        labs{i} = ...
            char(local_datetime_label(tick_dates(i)));
    end
    
    set(axh, ...
        'XTick',xt, ...
        'XTickLabel',labs, ...
        'XTickLabelRotation',0);

end

function nlab = local_max_date_labels(axh,prd)
%LOCAL_MAX_DATE_LABELS Estimate how many horizontal date labels fit

    try
        oldUnits = axh.Units;
        axh.Units = 'pixels';
        pos = axh.Position;
        axh.Units = oldUnits;
    
        axw = pos(3);   % axis width in pixels
    
        if isfield(prd,'dt') ...
                && ~isempty(prd.dt) ...
                && any(prd.dt == [24 96])
            % Subdaily labels include time: dd-MMM-yyyy HH:mm
            px_per_label = 120;
        else
            % Daily labels: dd-MMM-yyyy
            px_per_label = 85;
        end
    
        nlab = floor(axw / px_per_label);
    
        % keep sensible bounds
        nlab = max(4,min(5,nlab));
    
    catch
        nlab = 5;
    end

end

function x = local_datetime_to_index_fourgroups( ...
    tag,prd,idx0,t,dt_days)
%LOCAL_DATETIME_TO_INDEX_FOURGROUPS Maps datetimes to indices on four-scenario axes.

    x = NaN;

    if isnat(t)
        return
    end

    try
        if local_is_rainfall_block_prd(prd)
            t0 = local_make_datetime(prd.ds);
            idx_ref = 1;
        else
            switch lower(tag)
                case {'tt','et'}
                    if isfield(prd,'dts') ...
                            && ~isempty(prd.dts)
                        t0 = local_make_datetime(prd.dts);
                    else
                        t0 = local_make_datetime(prd.ds);
                    end
                    idx_ref = idx0;

                case {'te','ee'}
                    if isfield(prd,'des') ...
                            && ~isempty(prd.des)
                        t0 = local_make_datetime(prd.des);
                    elseif isfield(prd,'dts') && ~isempty(prd.dts)
                        t0 = local_make_datetime(prd.dts);
                    else
                        t0 = local_make_datetime(prd.ds);
                    end
                    idx_ref = idx0;

                otherwise
                    return
            end
        end

        if isnat(t0)
            return
        end

        x = idx_ref + days(t - t0)/dt_days;

    catch
        x = NaN;
    end
end

function tf = local_is_rainfall_block_prd(prd)
%LOCAL_IS_RAINFALL_BLOCK_PRD True when period object is rainfall block

    try
        if isfield(prd,'method') ...
                && ~isempty(prd.method) ...
                && strcmpi(prd.method,'rainfall_block')
            tf = true;
            return
        end
    
        % Defensive fallback: rainfall_block uses full-period ds/de rather
        % than separate train/eval date fields dts/dte/des/dee.
        tf = isfield(prd,'ds') ...
            && ~isempty(prd.ds) ...
            && isfield(prd,'de') ...
            && ~isempty(prd.de) ...
            && ~isfield(prd,'dts');
    catch
        tf = false;
    end

end

function tf = local_use_dates_on_xaxis(prd,mdl)
%LOCAL_USE_DATES_ON_XAXIS Tests whether a plot should use calendar-date ticks.

    try
        if isfield(mdl,'sp_method') ...
                && ~isempty(mdl.sp_method)
            sm = lower(string(mdl.sp_method));
        else
            sm = "";
        end

        if strcmp(sm,"rainfall_block")
            tf = isfield(prd,'ds') ...
                && ~isempty(prd.ds) ...
                && isfield(prd,'de') ...
                && ~isempty(prd.de);
            return
        end

        % Any global train/eval split with explicit date fields should
        % use date labels, including deterministic_block/block.
        tf = isfield(prd,'dts') ...
            && ~isempty(prd.dts) ...
            && isfield(prd,'des') ...
            && ~isempty(prd.des);

        if ~tf
            tf = isfield(prd,'ds') ...
                && ~isempty(prd.ds);
        end

    catch
        tf = false;
    end

end

function local_set_integer_ticks(axh,xmin,xmax)
%LOCAL_SET_INTEGER_TICKS Set sensible integer ticks for non-2x2 plots

    if ~isfinite(xmin) ...
            || ~isfinite(xmax)
        return
    end
    
    if xmax < xmin
        tmp = xmin;
        xmin = xmax;
        xmax = tmp;
    end
    
    if xmax == xmin
        set(axh,'XTick',xmin,'XTickLabel', ...
            {sprintf('%d',round(xmin))});
        return
    end
    
    span = xmax - xmin;
    
    if span <= 10
        step = 1;
    elseif span <= 20
        step = 2;
    elseif span <= 50
        step = 5;
    elseif span <= 100
        step = 10;
    elseif span <= 200
        step = 20;
    elseif span <= 500
        step = 50;
    elseif span <= 1000
        step = 100;
    else
        step = 200;
    end
    
    xt = xmin:step:xmax;
    
    if isempty(xt) ...
            || xt(1) ~= xmin
        xt = [xmin xt];
    end
    if xt(end) ~= xmax
        xt = [xt xmax];
    end
    
    xt = unique(round(xt),'stable');
    
    labs = arrayfun(@(x) ...
        sprintf('%d',x),xt, ...
        'UniformOutput',false);
    
    set(axh, ...
        'XTick',xt, ...
        'XTickLabel',labs, ...
        'XTickLabelRotation',0);

end

function [Qout,NSEout,KGEout,JKGEout] = ...
    local_apply_mode_filter(mode,Qin,NSEin,KGEin,JKGEin)
%LOCAL_APPLY_MODE_FILTER Hide scenarios that are not active for mdl.mode
%
% mode 1 -> tt only
% mode 2 -> tt and te
% mode 3 -> tt and et
% mode 4 -> tt, te, et, ee

    Qout = Qin;
    NSEout = NSEin;
    KGEout = KGEin;
    JKGEout = JKGEin;
    
    % start by clearing all non-tt scenarios
    fields = {'te','et','ee'};
    for i = 1:numel(fields)
        f = fields{i};
    
        if isfield(Qout,f)
            Qout.(f) = [];
        end
        if isfield(NSEout,f)
            NSEout.(f) = [];
        end
        if isfield(KGEout,f)
            KGEout.(f) = [];
        end
        if isfield(JKGEout,f)
            JKGEout.(f) = [];
        end
    end
    
    switch mode
        case 1
            % keep tt only
    
        case 2
            if isfield(Qin,'te')
                Qout.te = Qin.te;   
            end
            if isfield(NSEin,'te')
                NSEout.te = NSEin.te; 
            end
            if isfield(KGEin,'te')
                KGEout.te = KGEin.te; 
            end
            if isfield(JKGEin,'te')
                JKGEout.te = JKGEin.te; 
            end
    
        case 3
            if isfield(Qin,'et')
                Qout.et = Qin.et;
            end
            if isfield(NSEin,'et')
                NSEout.et = NSEin.et;
            end
            if isfield(KGEin,'et')
                KGEout.et = KGEin.et;
            end
            if isfield(JKGEin,'et')
                JKGEout.et = JKGEin.et; 
            end
    
        case 4
            if isfield(Qin,'te')
                Qout.te = Qin.te;
            end
            if isfield(NSEin,'te')
                NSEout.te = NSEin.te;
            end
            if isfield(KGEin,'te')
                KGEout.te = KGEin.te;
            end
            if isfield(JKGEin,'te')
                JKGEout.te = JKGEin.te; 
            end
    
            if isfield(Qin,'et')
                Qout.et = Qin.et;
            end
            if isfield(NSEin,'et')
                NSEout.et = NSEin.et; 
            end
            if isfield(KGEin,'et')
                KGEout.et = KGEin.et; 
            end
            if isfield(JKGEin,'et')
                JKGEout.et = JKGEin.et; 
            end
    
            if isfield(Qin,'ee')
                Qout.ee = Qin.ee;   
            end
            if isfield(NSEin,'ee')
                NSEout.ee = NSEin.ee;
            end
            if isfield(KGEin,'ee')
                KGEout.ee = KGEin.ee;
            end
            if isfield(JKGEin,'ee')
                JKGEout.ee = JKGEin.ee; 
            end
    
        otherwise
            error(['      Error: plot_SAGE: ' ...
                'mdl.mode must be one of 1,2,3,4.']);
    end
end

function s = local_model_display_name(nameIn)
%LOCAL_MODEL_DISPLAY_NAME Return model name for display only

    s = char(string(nameIn));   % keep original capitalization
    s = strtrim(s);
    
    if strcmpi(s,'xinanjiang')
        s = 'Xinanjiang';
    end

end

function ztxt = local_zone_label_short(bas,k)
%LOCAL_ZONE_LABEL_SHORT Return compact hydroclimatic zone label

    ztxt = '';
    
    try
        if ~isfield(bas,'zone') ...
                || ~isfield(bas.zone,'num') ...
                || isempty(bas.zone.num) ...
                || k > numel(bas.zone.num)
            return
        end
    
        zn = double(bas.zone.num(k));
    
        zoneNames = ...
            local_zone_names_from_bas(bas);
    
        if isfinite(zn) ...
                && zn >= 1 ...
                && zn <= numel(zoneNames)
            ztxt = sprintf('Z%d: %s', ...
                zn,char(zoneNames(zn)));
        elseif isfinite(zn)
            ztxt = sprintf('Z%d',zn);
        end
    
    catch
        ztxt = '';
    end
end

function add_basin_label(ax,us,zoneTxt,gauge,n_char,topOffsetCm,textColor)
%ADD_BASIN_LABEL Add compact basin metadata label inside time-series panel

    if nargin < 5
        n_char = inf;
    end
    if nargin < 6
        topOffsetCm = 0;
    end
    if nargin < 7 || isempty(textColor)
        textColor = [0.15 0.15 0.15];
    end
    
    txt = strings(0,1);
    
    us = string(us);
    us = erase(us,"gauge: ");
    if strlength(us) > 0
        txt(end+1,1) = "gauge " + us;
    end
    
    zoneTxt = string(zoneTxt);
    if strlength(zoneTxt) > 0
        txt(end+1,1) = zoneTxt;
    end
    
    gauge = string(gauge);
    if strlength(gauge) > n_char
        gauge = extractBefore(gauge,n_char) + "...";
    end
    sameGauge = strcmpi(strtrim(gauge),strtrim(us));
    if strlength(gauge) > 0 && ~sameGauge
        txt(end+1,1) = gauge;
    end
    
    if isempty(txt)
        return
    end
    
    drawnow limitrate;
    
    offsetFraction = local_axes_vertical_fraction(ax,topOffsetCm);
    x0 = 0.02;
    y0 = 1-max(0.02,offsetFraction);
    
    text(ax,x0,y0,strjoin(txt,newline), ...
        'Units','normalized', ...
        'horizontalalignment','left', ...
        'verticalalignment','top', ...
        'interpreter','none', ...
        'fontsize',12, ...
        'color',textColor, ...
        'BackgroundColor','w', ...
        'Margin',3, ...
        'clipping','on');
    
    % The label is created after the forcing and discharge graphics and is
    % therefore already the uppermost ordinary child. Do not call uistack
    % here: axes configured with yyaxis contain internal ruler children for
    % which uistack can fail with "Children may only be set to a permutation
    % of itself" and abort creation of the remaining time-series panels.
end

function tf = local_incomplete_forcing(datk)
%LOCAL_INCOMPLETE_FORCING True when P, Ep, or T is absent or nonfinite.

    tf = true;
    if ~isstruct(datk) || ~isfield(datk,'meteo') ...
            || ~isstruct(datk.meteo)
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
            return
        end
        n(i) = numel(datk.meteo.(name));
    end
    tf = ~all(n == n(1));
end

function local_add_invalid_forcing_notice(axh)
%LOCAL_ADD_INVALID_FORCING_NOTICE Mark an excluded basin prominently.

    text(axh,0.5,0.5, ...
        '\textbf{Basin excluded: incomplete forcing record}', ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Interpreter','latex', ...
        'FontSize',20, ...
        'FontWeight','bold', ...
        'Color',[0.80 0 0], ...
        'BackgroundColor','white', ...
        'Margin',8, ...
        'Clipping','on');
end

function fraction = local_axes_vertical_fraction(ax,distanceCm)
%LOCAL_AXES_VERTICAL_FRACTION Convert a physical distance to axis units.

    fraction = 0;
    if distanceCm <= 0 || ~isgraphics(ax)
        return
    end

    oldUnits = ax.Units;
    cleanup = onCleanup(@()set(ax,'Units',oldUnits));
    ax.Units = 'centimeters';
    heightCm = ax.Position(4);
    fraction = min(0.25,distanceCm/max(heightCm,eps));
    clear cleanup
end

function zoneNames = local_zone_names_from_bas(bas)
%LOCAL_ZONE_NAMES_FROM_BAS Return compact display names for zones

    zoneNames = strings(0,1);

    try
        if isfield(bas,'zone') ...
                && isfield(bas.zone,'names_short') ...
                && ~isempty(bas.zone.names_short)

            zoneNames = string(bas.zone.names_short(:));

        elseif isfield(bas,'zone') ...
                && isfield(bas.zone,'names') ...
                && ~isempty(bas.zone.names)

            zoneNames = string(bas.zone.names(:));
            zoneNames = local_zone_display_name(zoneNames);
        end
    catch
        zoneNames = strings(0,1);
    end

end

function z = local_zone_display_name(z)
%LOCAL_ZONE_DISPLAY_NAME Shorten long region-specific zone labels

    z = string(z(:));
    z = strtrim(strrep(z,'_',' '));
    zl = lower(z);

    for i = 1:numel(z)
        s = zl(i);
        s = regexprep(s,'\s+',' ');
        s = strtrim(s);

        % First normalize rain/snow wording
        s = replace(s,"rain dominated","rain");
        s = replace(s,"snow dominated","snow");
        s = replace(s,"rainfall dominated","rain");
        s = replace(s,"snowfall dominated","snow");

        % Compact hydroclimatic labels
        s = replace(s,"very humid rain","VH rain");
        s = replace(s,"very humid snow","VH snow");
        s = replace(s,"sub humid rain","SH rain");
        s = replace(s,"sub humid snow","SH snow");
        s = replace(s,"subhumid rain","SH rain");
        s = replace(s,"subhumid snow","SH snow");
        s = replace(s,"humid rain","H rain");
        s = replace(s,"humid snow","H snow");
        s = replace(s,"dry rain","D rain");
        s = replace(s,"dry snow","D snow");
        s = replace(s,"semi arid rain","SA rain");
        s = replace(s,"semi arid snow","SA snow");
        s = replace(s,"arid rain","D rain");
        s = replace(s,"arid snow","D snow");

        % Brazil-style labels
        s = replace(s,"low runoff / forest","Low Q / forest");
        s = replace(s,"moderate runoff / forest","Mod Q / forest");
        s = replace(s,"high runoff / forest","High Q / forest");
        s = replace(s,"low runoff / mixed","Low Q / mixed");
        s = replace(s,"moderate runoff / mixed","Mod Q / mixed");
        s = replace(s,"high runoff / mixed","High Q / mixed");

        z(i) = strtrim(s);
    end
end

function txt = local_sib_text(SIBz)
%LOCAL_SIB_TEXT Formats an integrated basin loss score for display.

    txt = '';
    
    if isempty(SIBz) ...
            || ~isfinite(SIBz)
        return
    end
    
    txt = sprintf(['$\\widehat{S}_' ...
        '{\\rm ib} = %.3f$'],SIBz);

end

function symbol = local_metric_symbol(metricTag)
%LOCAL_METRIC_SYMBOL Format metric abbreviations consistently for axes.

    tag = char(string(metricTag));
    if contains(lower(tag),'fdc')
        symbol = 'S_{\mathrm{fdc}}';
        return
    end
    switch upper(strtrim(tag))
        case {'NSE','KGE','JKGE'}
            symbol = sprintf('\\mathrm{%s}',upper(strtrim(tag)));
        otherwise
            symbol = strrep(tag,'_','\\_');
    end
end

function symbol = local_metric_median_symbol(metricTag)
%LOCAL_METRIC_MEDIAN_SYMBOL Sample-median estimator for a metric.

    tag = upper(strtrim(char(string(metricTag))));
    if contains(lower(tag),'fdc')
        symbol = '\\widehat{T}_{S_{\\rm fdc}}';
        return
    end
    switch tag
        case 'NSE'
            symbol = '\\widehat{T}_{\\rm nse}';
        case 'KGE'
            symbol = '\\widehat{T}_{\\rm kge}';
        case 'JKGE'
            symbol = '\\widehat{T}_{\\rm jkge}';
        otherwise
            symbol = sprintf('\\widehat{T}_{%s}', ...
                local_metric_symbol(metricTag));
    end
end

function local_check_rainfall_indices(dat,K)
%LOCAL_CHECK_RAINFALL_INDICES Verify basin-specific rainfall split indices.

    for k = 1:K
        if ~isfield(dat{k},'id_train') ...
                || isempty(dat{k}.id_train) ...
                || ~isfield(dat{k},'id_eval') ...
                || isempty(dat{k}.id_eval)
            error(['plot_SAGE: rainfall_block requires ' ...
                'dat{%d}.id_train and dat{%d}.id_eval.'],k,k);
        end
    end
end

function s = local_latex_escape(s)
%LOCAL_LATEX_ESCAPE Escapes special characters for LaTeX-rendered text.

    s = char(s);
    
    s = strrep(s,'_','\_');
    s = strrep(s,'%','\%');
    s = strrep(s,'&','\&');
    s = strrep(s,'#','\#');
    
    end
    
    function local_plot_fdc_postprocessor( ...
        fig0,part,dt_str,mdl,dat,bas,Qfdc, ...
        model_name_fig,q_unit,figH, ...
        fontsize_legend,samp_word,region)
    
    if isempty(Qfdc) ...
            || ~isstruct(Qfdc)
        return
    end
    if ~isfield(Qfdc,'gauge') ...
            || isempty(Qfdc.gauge)
        return
    end
    if ~isfield(Qfdc,'qy') ...
            || isempty(Qfdc.qy)
        return
    end
    
    % Match retained Qfdc basins back to dat{k}
    usQ = strings(numel(Qfdc.gauge),1);
    for iq = 1:numel(Qfdc.gauge)
        usQ(iq) = string(local_basin_code( ...
            Qfdc.gauge(iq),region));
    end
    
    usD = strings(numel(dat),1);
    for k = 1:numel(dat)
        try
            usD(k) = string(local_basin_code( ...
                dat{k}.gauge,region));
        catch
            usD(k) = "";
        end
    end
    
    kRet = [];
    iqRet = [];
    
    usQmatch = string(local_normalize_basin_code_for_match( ...
        cellstr(usQ)));
    usDmatch = string(local_normalize_basin_code_for_match( ...
        cellstr(usD)));

    for iq = 1:numel(usQmatch)
        k = find(usDmatch == usQmatch(iq),1,'first');
        if ~isempty(k)
            kRet(end+1,1) = k;   %#ok<AGROW>
            iqRet(end+1,1) = iq;  %#ok<AGROW>
        end
    end
    
    if isempty(kRet)
        warning('plot_SAGE:FDCNoMatch', ...
            ['No Qfdc.gauge entries ' ...
            'matched dat{k}.gauge.']);
        return
    end
    
    % Basin groups: first K_t are 
    % training, next K_e are evaluation
    idTrainBas = (1:bas.K_t).';
    idEvalBas = (bas.K_t+1:bas.K_t+bas.K_e).';
    
    isTrain = ismember(kRet,idTrainBas);
    isEval = ismember(kRet,idEvalBas);
    
    pTrain = [];
    pEval = [];
    
    if isfield(Qfdc,'req') ...
            && isstruct(Qfdc.req)
        ipT = [];
        if isfield(Qfdc.req,'tt') ...
                && ~isempty(Qfdc.req.tt)
            ipT = [ipT; Qfdc.req.tt(:)];
        end
        if isfield(Qfdc.req,'te') ...
                && ~isempty(Qfdc.req.te)
            ipT = [ipT; Qfdc.req.te(:)];
        end
        ipT = unique(ipT,'stable');
        ipT = ipT(ipT >= 1 ...
            & ipT <= numel(Qfdc.id));
        if ~isempty(ipT)
            pTrain = double(Qfdc.id(ipT));
        end
    
        ipE = [];
        if isfield(Qfdc.req,'et') ...
                && ~isempty(Qfdc.req.et)
            ipE = [ipE; Qfdc.req.et(:)];
        end
        if isfield(Qfdc.req,'ee') ...
                && ~isempty(Qfdc.req.ee)
            ipE = [ipE; Qfdc.req.ee(:)];
        end
        ipE = unique(ipE,'stable');
        ipE = ipE(ipE >= 1 ...
            & ipE <= numel(Qfdc.id));
        if ~isempty(ipE)
            pEval = double(Qfdc.id(ipE));
        end
    end
    
    
    if isempty(pTrain)
        pTrain = kRet(isTrain);
    end
    if isempty(pEval)
        pEval = kRet(isEval);
    end
    
    figNo = fig0;
    
    for j0 = 1:2:numel(pTrain)
        ids = pTrain(j0:min(j0+1,numel(pTrain)));
    
        local_plot_one_fdc_group_page( ...
            figNo,part,dt_str,model_name_fig, ...
            'train',mdl,dat,bas,Qfdc, ...
            ids,q_unit,figH, ...
            fontsize_legend, ...
            samp_word,region);
    
        figNo = figNo + 1;
    end
    
    for j0 = 1:2:numel(pEval)
        ids = pEval(j0:min(j0+1,numel(pEval)));
    
        local_plot_one_fdc_group_page( ...
            figNo,part,dt_str,model_name_fig, ...
            'eval',mdl,dat,bas,Qfdc, ...
            ids,q_unit,figH, ...
            fontsize_legend, ...
            samp_word,region);
    
        figNo = figNo + 1;
    end
end

function local_plot_one_fdc_group_page( ...
    figNo,part,dt_str,model_name_fig, ...
    groupName,mdl,dat,bas,Qfdc,kList, ...
    q_unit,figH,fontsize_legend, ...
    samp_word,region)
%LOCAL_PLOT_ONE_FDC_GROUP_PAGE Plots one page of scenario-specific FDC panels.

    if isempty(kList)
        basinRange = '';
    else
        basinRange = sprintf(' basins %d-%d', ...
            min(kList),max(kList));
    end
    
    if strcmpi(groupName,'train')
        scTags = {'tt','te','tt','te'};
        scNames = { ...
            sprintf(['training basin ' ...
            '| training %s'],samp_word), ...
            sprintf(['training basin ' ...
            '| evaluation %s'],samp_word), ...
            sprintf(['training basin ' ...
            '| training %s'],samp_word), ...
            sprintf(['training basin ' ...
            '| evaluation %s'],samp_word)};
    else
        scTags = {'et','ee','et','ee'};
        scNames = { ...
            sprintf(['evaluation basin ' ...
            '| training %s'],samp_word), ...
            sprintf(['evaluation basin ' ...
            '| evaluation %s'],samp_word), ...
            sprintf(['evaluation basin ' ...
            '| training %s'],samp_word), ...
            sprintf(['evaluation basin ' ...
            '| evaluation %s'],samp_word)};
    end
    
    figName = sprintf('%s: FDC %s%s', ...
        model_name_fig,groupName,basinRange);
    
    figW = 0.925*figH;
    figH0 = 0.925*figH;
    
    % One FDC row contains the two periods for one basin.  If only one
    % basin is present on this page, use a one-row page instead of filling
    % the lower half with empty "Not available" panels.
    nActiveRows = min(2,numel(kList));
    if nActiveRows < 1
        nActiveRows = 1;
    end
    figH = figH0 * (0.56 + 0.44*(nActiveRows == 2));
    
    figure(figNo); clf;
    set(gcf,'color','w', ...
        'NumberTitle','off', ...
        'Units','inches', ...
        'Position',[0.2 0.2 figW figH], ...
        'Name',figName);
    
    xleft = 0.090;
    xgap = 0.455;
    dx = 0.375;
    if nActiveRows == 1
        yleft = 0.180;
        dy = 0.700;
        axPos = [ ...
            xleft        , yleft , dx , dy; ...
            xleft + xgap , yleft , dx , dy];
    else
        yleft = 0.080;
        ygap = 0.470;
        dy = 0.375;
        axPos = [ ...
            xleft        ,  yleft + ygap , dx , dy; ...
            xleft + xgap ,  yleft + ygap , dx , dy; ...
            xleft        ,  yleft        , dx , dy; ...
            xleft + xgap ,  yleft        , dx , dy];
    end
    
    nPanels = 2*nActiveRows;
    
    for is = 1:nPanels
        ax = axes('Position',axPos(is,:));
        hold(ax,'on');
        box(ax,'off');
        grid(ax,'off');
    
        ib = ceil(is/2);
        showXLabel = ib == nActiveRows;
        showYLabel = mod(is,2) == 1;
        title(ax, ['\textbf{' scNames{is} '}'], ...
            'interpreter','latex', ...
            'fontsize',13);
            
        if ib > numel(kList)
            text(ax,0.5,0.5,'Not available', ...
                'Units','normalized', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'fontsize',13);
            local_style_fdc_axes(ax,q_unit,showXLabel,showYLabel);
            ax.Title.FontSize = 15;
            continue
        end
    
        k = double(kList(ib));
        iq = local_qfdc_lookup(Qfdc,dat,bas,k,region);
        tag = scTags{is};
        isSelected = local_qfdc_has_requested_scenario( ...
            Qfdc,iq,tag);
    
        if ~isSelected
            text(ax,0.5,0.5,'Not selected', ...
                'Units','normalized', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'fontsize',13);
            ok = false;
        else
            [yObs,qSim,ok] = local_get_fdc_series( ...
                Qfdc,dat,k,iq,tag,mdl);
    
            if ok
                [cObs,cSim] = local_fdc_colors(tag);
                local_plot_fdc_pair(ax,yObs,qSim,cObs,cSim);
            else
                text(ax,0.5,0.5,'Not available', ...
                    'Units','normalized', ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'fontsize',13);
            end
        end
    
        local_style_fdc_axes(ax,q_unit,showXLabel,showYLabel);
        if ~isSelected
            set(ax,'XTick',[], ...
                'XTickLabel',{}, ...
                'XMinorTick','off');
            xlabel(ax,'');
        end
        ax.Title.FontSize = 15;
    
        us = local_fdc_usgs_string(dat,bas,k,region);
        ztxt = local_zone_label_short(bas,k);
    
        infoTxt = us;
        if ~isempty(ztxt)
            infoTxt = sprintf('%s\n%s',us,ztxt);
        end
        
        text(ax,0.02,0.02,infoTxt, ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom', ...
            'fontsize',12, ...
            'interpreter','none');
    
        if any(is == [1 2]) ...
                && ok
            lg = legend(ax,'show', ...
                'Box','off', ...
                'interpreter','latex', ...
                'fontsize',fontsize_legend);
            lg.Units = 'normalized';
            lg.Location = 'none';
            if nActiveRows == 2 ...
                    && is == 1
                % Top-left panel: shift right, but not too far.
                lg.Position(1) = ax.Position(1) ...
                    + 0.10*ax.Position(3);            
            elseif nActiveRows == 2 ...
                    && is == 2
                % Top-right panel: also shift right, away from y-axis.
                lg.Position(1) = ax.Position(1) ...
                    + 0.10*ax.Position(3);            
            else
                % One-row / lower panels: original placement.
                lg.Position(1) = ax.Position(1) ...
                    + 0.02*ax.Position(3);
            end
            lg.Position(2) = ax.Position(2) ...
                + 0.28*ax.Position(4);
        end
    end
    
    part_name = local_part_name(part,mdl);
    part_name_tex = strrep(part_name,'_','\_');
    
    annotation(gcf,'textbox',[0 0.96 1 0.03], ...
        'String',sprintf(['\\texttt{%s}: ' ...
        '%s Flow Duration Curve'], ...
        part_name_tex,dt_str), ...
        'interpreter','latex', ...
        'fontsize',18, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'EdgeColor','none');
end

function tf = local_qfdc_has_requested_scenario( ...
    Qfdc,iq,tag)
%LOCAL_QFDC_HAS_REQUESTED_SCENARIO Tests whether Qfdc contains a requested scenario.

    tf = true;
    
    if isempty(Qfdc) ...
            || ~isstruct(Qfdc) ...
            || ~isfield(Qfdc,'req') ...
            || ~isstruct(Qfdc.req)
        return
    end
    
    if ~isfield(Qfdc.req,tag)
        tf = false;
        return
    end
    
    r = Qfdc.req.(tag);
    r = double(r(:));
    
    tf = any(r == double(iq));

end

function [yObs,qSim,ok] = local_get_fdc_series( ...
    Qfdc,dat,k,iq,tag,mdl)
%LOCAL_GET_FDC_SERIES Retrieves observed and simulated FDC source series.

    ok = false;
    yObs = [];
    qSim = [];
    
    % Simulated discharge retained by camels.m/Qfdc.
    % Simulated/measured discharge retained by camels.m/Qfdc.
    try
        if isfield(Qfdc,'qy') ...
                && ~isempty(Qfdc.qy)
    
            if iscell(Qfdc.qy)
                qy = double(Qfdc.qy{iq});
    
                if isempty(qy) ...
                        || size(qy,2) < 2
                    return
                end
    
                qFull = qy(:,1);
                yFull = qy(:,2);
    
            else
                A = Qfdc.qy;
    
                if ndims(A) == 3
                    yFull = A(:,iq,1);
                    qFull = A(:,iq,2);
                elseif size(A,3) >= 2
                    yFull = A(:,iq,1);
                    qFull = A(:,iq,2);
                else
                    qFull = A(:,iq);
                    yFull = dat{k}.y_n(:);
                end
            end
    
        elseif isfield(Qfdc,'q') ...
                && ~isempty(Qfdc.q)
            A = Qfdc.q;
            qFull = A(:,iq);
            yFull = dat{k}.y_n(:);
        else
            return
        end
    
        yFull = yFull(:);
        qFull = qFull(:);
    
    catch
        return
    end
    
    n = min(numel(yFull),numel(qFull));
    yFull = yFull(1:n);
    qFull = qFull(1:n);
    
    idx = local_fdc_indices_for_basin(Qfdc,dat,k,iq,tag,mdl,n);
    
    if isempty(idx)
        return
    end
    
    idx = idx(idx >= 1 ...
        & idx <= n);
    
    if isfield(dat{k},'bad') ...
            && ~isempty(dat{k}.bad)
        badk = dat{k}.bad(:);
        if numel(badk) >= max(idx)
            idx = idx(~badk(idx));
        end
    end
    
    if isempty(idx)
        return
    end
    
    yObs = yFull(idx);
    qSim = qFull(idx);
    
    good = isfinite(yObs) ...
        & isfinite(qSim);
    yObs = yObs(good);
    qSim = qSim(good);
    
    ok = numel(yObs) >= 2 ...
        && numel(qSim) >= 2;
end


function idx = local_fdc_indices_for_basin( ...
    Qfdc,dat,k,iq,tag,mdl,n)
%LOCAL_FDC_INDICES_FOR_BASIN Resolves valid FDC indices for one basin and scenario.

    idx = [];
    
    % Best case: Qfdc.idx has cell arrays attached by scenario and basin.
    try
        if isfield(Qfdc,'idx') ...
                && isstruct(Qfdc.idx) ...
                && isfield(Qfdc.idx,tag)
            Z = Qfdc.idx.(tag);
            if iscell(Z) ...
                    && numel(Z) >= iq
                idx = Z{iq}(:);
            elseif isnumeric(Z)
                idx = Z(:);
            end
        end
    catch
        idx = [];
    end
    
    if ~isempty(idx)
        return
    end
    
    % Local rainfall split: use basin-specific train/eval indices.
    isLocal = false;
    try
        isLocal = isfield(mdl,'local') ...
            && mdl.local == 1;
    catch
    end
    try
        isLocal = isLocal ...
            || strcmpi(mdl.sp_method, ...
            'rainfall_block');
    catch
    end
    
    if isLocal
        switch lower(tag)
            case {'tt','et'}
                if isfield(dat{k},'id_train')
                    idx = dat{k}.id_train(:);
                end
            case {'te','ee'}
                if isfield(dat{k},'id_eval')
                    idx = dat{k}.id_eval(:);
                end
        end
        return
    end
    
    % Global split fallback.
    switch lower(tag)
        case {'tt','et'}
            if isfield(mdl,'id_train')
                idx = expand_index(mdl.id_train);
            end
        case {'te','ee'}
            if isfield(mdl,'id_eval')
                idx = expand_index(mdl.id_eval);
            end
    end
    
    idx = idx(:);
    idx = idx(idx >= 1 ...
        & idx <= n);
end

function [cObs,cSim] = local_fdc_colors(tag)
%LOCAL_FDC_COLORS Returns observed and simulated colors for an FDC scenario.

    switch lower(tag)
        case 'tt'
            cObs = [0.50 0.50 0.50];
            cSim = [0.00 0.00 0.00];
    
        case 'te'
            cObs = [0.93 0.74 0.56];
            cSim = [0.85 0.33 0.10];
    
        case 'et'
            cObs = [0.62 0.80 0.95];
            cSim = [0.00 0.45 0.74];
    
        case 'ee'
            cObs = [0.78 0.69 0.86];
            cSim = [0.49 0.18 0.56];
    
        otherwise
            cObs = [0.50 0.50 0.50];
            cSim = [0.00 0.00 0.00];
    end
end

function local_plot_fdc_pair(ax,yObs,qSim,cObs,cSim)
%LOCAL_PLOT_FDC_PAIR Plots paired observed and simulated flow-duration curves.

    [yS,pY] = local_fdc_curve(yObs);
    [qS,pQ] = local_fdc_curve(qSim);

    setappdata(ax,'SAGE_fdc_xobs',yS);
    setappdata(ax,'SAGE_fdc_xsim',qS);
    setappdata(ax,'SAGE_fdc_psim',pQ);

    if ~isempty(yS)
        plot(ax,yS,pY,'o', ...
            'color',cObs, ...
            'MarkerFaceColor','w', ...
            'MarkerEdgeColor',cObs, ...
            'MarkerSize',4, ...
            'linewidth',1, ...
            'LineStyle','none', ...
            'handlevisibility','off');
    end

    if ~isempty(qS)
        plot(ax,qS,pQ,'-', ...
            'color',cSim, ...
            'linewidth',1.8, ...
            'handlevisibility','off');
    end

    plot(ax,nan,nan,'o', ...
        'color',cObs, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor',cObs, ...
        'MarkerSize',7, ...
        'linewidth',2, ...
        'LineStyle','none', ...
        'DisplayName','$Q_{\rm obs}$');

    plot(ax,nan,nan,'-', ...
        'color',cSim, ...
        'linewidth',2.2, ...
        'DisplayName','$Q_{\rm sim}$');
end

function [qSort,p] = local_fdc_curve(q)
%LOCAL_FDC_CURVE Sorts discharge and returns exceedance probabilities.

    q = q(:);
    q = q(isfinite(q));
    q = q(q >= 0);
    
    if isempty(q)
        qSort = [];
        p = [];
        return
    end
    
    qSort = sort(q,'descend');
    n = numel(qSort);
    p = ((1:n)' - 0.5) ./ n;
    p = max(p,1e-6);
end

function local_style_fdc_axes(ax,q_unit,showXLabel,showYLabel)
%LOCAL_STYLE_FDC_AXES Applies common flow-duration-curve axis styling.

    set(ax,'XScale','log', ...
        'YScale','linear', ...
        'TickDir','out', ...
        'TickLength',[0.02 0.02], ...
        'fontname','Times', ...
        'fontsize',15, ...       % tick-label font size
        'linewidth',1, ...
        'XMinorTick','on', ...
        'YMinorTick','off');
    
    axis(ax,'square');
    
    if nargin < 3 ...
            || isempty(showXLabel)
        showXLabel = true;
    end
    if nargin < 4 ...
            || isempty(showYLabel)
        showYLabel = true;
    end
    
    if showXLabel
        hX = xlabel(ax, ...
            ['Discharge (' q_unit ')'], ...
            'interpreter','latex', ...
            'fontsize',14);
        hX.Units = 'normalized';
        posX = hX.Position;
        posX(1) = 0.5;
        posX(2) = -0.095;
        hX.Position = posX;
    else
        xlabel(ax,'');
    end
    
    if showYLabel
        hY = ylabel(ax, ...
            'Exceedance probability', ...
            'interpreter','latex', ...
                'fontsize',14);
        hY.Units = 'normalized';
        posY = hY.Position;
        posY(1) = posY(1)-0.025;
        hY.Position = posY;
    else
        ylabel(ax,'');
        set(ax,'YTickLabel',[]);
    end
    
    ylim(ax,[-0.02 1.02]);
    
    % Smart x-limits:
    % observed curve keeps its full low-flow range;
    % simulated curve may trim the flat tail near p = 1.
    xObs = [];
    xSim = [];
    pSim = [];

    try
        if isappdata(ax,'SAGE_fdc_xobs')
            xObs = getappdata(ax,'SAGE_fdc_xobs');
        end
        if isappdata(ax,'SAGE_fdc_xsim')
            xSim = getappdata(ax,'SAGE_fdc_xsim');
        end
        if isappdata(ax,'SAGE_fdc_psim')
            pSim = getappdata(ax,'SAGE_fdc_psim');
        end
    catch
        xObs = [];
        xSim = [];
        pSim = [];
    end

    % Fallback for older plots.
    if isempty(xObs) ...
            && isempty(xSim)
        xraw = [];
        h = findobj(ax,'Type','line');
        for i = 1:numel(h)
            x = h(i).XData(:);
            xraw = [xraw; x(:)]; %#ok<AGROW>
        end
        xraw = double(xraw(:));
        xraw = xraw(isfinite(xraw) ...
            & xraw > 0);
        xObs = xraw;
        xSim = [];
        pSim = [];
    end

    if ~isempty(xObs) ...
            || ~isempty(xSim)

        [xmin,xmax] = local_fdc_smart_xlimits( ...
            xObs,xSim,pSim,true);

        xlim(ax,[xmin xmax]);

        [xt,xtlbl] = ...
            local_fdc_major_xticks(xmin,xmax);

        set(ax, ...
            'XTickMode','manual', ...
            'XTick',xt, ...
            'XTickLabelMode','manual', ...
            'XTickLabel',xtlbl, ...
            'XTickLabelRotation',0, ...
            'TickLabelInterpreter','latex');

        try
            ax.XRuler.TickLabelGapOffset = -2;
            ax.YRuler.TickLabelGapOffset = -1;
        catch
        end
    end
    
    local_draw_top_right_box_no_legend(ax);

end

function [xmin,xmax] = local_fdc_smart_xlimits( ...
    xObs,xSim,pSim,useLog)
%LOCAL_FDC_SMART_XLIMITS Selects log-scale discharge limits from plotted curves.

    xObs = double(xObs(:));
    xSim = double(xSim(:));
    pSim = double(pSim(:));
    
    xObs = xObs(isfinite(xObs));
    
    goodSim = isfinite(xSim) & isfinite(pSim);
    xSim = xSim(goodSim);
    pSim = pSim(goodSim);
    
    if useLog
        xObs = xObs(xObs > 0);
    
        goodSim = xSim > 0;
        xSim = xSim(goodSim);
        pSim = pSim(goodSim);
    end
    
    if ~isempty(xObs)
        xminObs = min(xObs);
    else
        xminObs = inf;
    end
    
    if ~isempty(xSim)
        keepSim = pSim <= 0.98;
        if any(keepSim)
            xminSim = min(xSim(keepSim));
        else
            xminSim = min(xSim);
        end
    else
        xminSim = inf;
    end
    
    xmin = min([xminObs xminSim]);
    xmax = max([xObs(:); xSim(:)]);
    
    if ~isfinite(xmin) || ~isfinite(xmax)
        xmin = 1e-5;
        xmax = 1;
    end
    
    if useLog
        xmin = max(xmin,1e-5);
        if xmax <= xmin
            xmax = xmin * 10;
        end
    else
        if xmax <= xmin
            xmax = xmin + 1;
        end
    end
end

function [xt,xtlbl] = local_fdc_major_xticks(xmin,xmax)
%LOCAL_FDC_MAJOR_XTICKS Major log ticks for FDC panels.
%
% Limits are controlled by the data domain in local_style_fdc_axes.
% Here we select readable major-decade labels. If the visible range
% spans many decades, use an integer decade spacing and shift the first
% labeled decade upward when needed so the labels remain equally spaced
% and do not exceed about five labels.

    xmin = max(realmin,xmin);
    xmax = max(xmin*(1+eps),xmax);
    
    pLo = ceil(log10(xmin));
    pHi = floor(log10(xmax));
    
    % If fewer than three full decades are strictly inside the data limits,
    % use the nearest enclosing decades. This keeps the labels as powers of
    % ten instead of mixed numeric labels.
    if pHi < pLo || (pHi - pLo + 1) < 3
        pLo = floor(log10(xmin));
        pHi = ceil(log10(xmax));
    end
    
    % Guarantee at least three candidate decades where possible.
    while (pHi - pLo + 1) < 3
        pLo = pLo - 1;
        if (pHi - pLo + 1) >= 3
            break
        end
        pHi = pHi + 1;
    end
    
    pAll = pLo:pHi;
    maxLabels = 5;
    
    if numel(pAll) <= maxLabels
        p = pAll;
    else
        span = pHi - pLo;
    
        % Start with the smallest integer spacing that can keep the number
        % of labels near maxLabels. Then increase if needed.
        step = max(1,floor(span/(maxLabels-1)));
        while floor(span/step) + 1 > maxLabels
            step = step + 1;
        end
    
        % Anchor the highest label at pHi and move the first label upward
        % until all labels lie in the available decade range. Example:
        % pLo=-8, pHi=1 gives -7,-5,-3,-1,1 instead of uneven labels.
        pStart = pHi - step*(maxLabels-1);
        while pStart < pLo
            pStart = pStart + step;
        end
    
        p = pStart:step:pHi;
    end
    
    % If the above produced only 2 labels for a narrow domain, fall back to
    % three enclosing decades.
    if numel(p) < 3
        pMid = round((pLo + pHi)/2);
        p = unique([pLo pMid pHi],'stable');
    end
    
    xt = 10.^p;
    
    % Keep ticks finite/positive. We intentionally allow enclosing decade
    % ticks just outside raw data range for very narrow domains; this is
    % necessary to keep major-decade labels readable.
    xt = xt(:).';
    xt = xt(isfinite(xt) & xt > 0);
    
    xtlbl = cell(size(xt));
    for i = 1:numel(xt)
        xtlbl{i} = sprintf('$10^{%d}$',round(log10(xt(i))));
    end

end

function local_draw_top_right_box_no_legend(ax)
%LOCAL_DRAW_TOP_RIGHT_BOX_NO_LEGEND Completes an axes frame without adding legend objects.

    xl = xlim(ax);
    yl = ylim(ax);
    lw = ax.LineWidth;
    
    line(ax,[xl(1) xl(2)],[yl(2) yl(2)], ...
        'color','k', ...
        'linewidth',lw, ...
        'clipping','off', ...
        'handlevisibility','off');
    
    line(ax,[xl(2) xl(2)],[yl(1) yl(2)], ...
        'color','k', ...
        'linewidth',lw, ...
        'clipping','off', ...
        'handlevisibility','off');
end

function us = local_fdc_usgs_string(dat,bas,k,region)
%LOCAL_FDC_USGS_STRING Formats the regional gauge identifier for an FDC panel.

    us = sprintf('gauge %d',k);
    
    try
        us = ['gauge ' local_basin_code(dat{k}.gauge, ...
            region)];
        return
    catch
    end
    
    try
        us = ['gauge ' char(string(bas.id_gauge(k)))];
    catch
    end
end

function iq = local_qfdc_lookup(Qfdc,dat,bas,k,region)
%LOCAL_QFDC_LOOKUP Finds a basin in the stored flow-duration-curve selection.

    iq = [];
    
    try
        usTarget = string(local_basin_code( ...
            dat{k}.gauge,region));
    catch
        try
            usTarget = string(local_basin_code( ...
                bas.id_gauge(k),region));
        catch
            return
        end
    end
    
    try
        usQ = strings(numel(Qfdc.gauge),1);
        for ii = 1:numel(Qfdc.gauge)
            usQ(ii) = string(local_basin_code( ...
                Qfdc.gauge(ii),region));
        end
        iq = find(usQ == usTarget,1,'first');
    catch
        iq = [];
    end
end

function Bprt = local_default_Bprt_from_USGSscen_Qfdc( ...
    gaugescen,Qfdc)
    % Priority 1: explicit gaugescen entries
    u = strings(0,1);
    if isstruct(gaugescen) ...
            && ~isempty(fieldnames(gaugescen))
        fn = fieldnames(gaugescen);
        for j = 1:numel(fn)
            v = gaugescen.(fn{j});
            if isempty(v)
                continue
            end
            if iscell(v)
                v = string(v(:));
            elseif isnumeric(v)
                v = string(v(:));
            elseif ischar(v)
                v = string({v});
            elseif isstring(v)
                v = v(:);
            else
                continue
            end
            v = strtrim(v);
            v = v(v ~= "");
            u = [u; v]; %#ok<AGROW>
        end
    end
    if ~isempty(u)
        Bprt = numel(unique(u,'stable'));
        return
    end
    % Priority 2: Qfdc.id
    if isstruct(Qfdc) ...
            && isfield(Qfdc,'id') ...
            && ~isempty(Qfdc.id)
        Bprt = numel(Qfdc.id);
        return
    end
    % Priority 3: fallback
    Bprt = 10;
end

function idx = local_safe_index(S,field)
%LOCAL_SAFE_INDEX Returns a normalized index vector from a structure field.

    idx = [];
    
    try
        if ~isstruct(S) ...
                || ~isfield(S,field) ...
                || isempty(S.(field))
            return
        end
    
        idx = S.(field);
    
        if iscell(idx)
            idx = idx{1};
        end
    
        idx = double(idx(:));
        idx = idx(isfinite(idx) ...
            & idx >= 1);
        idx = unique(round(idx),'stable');
    
    catch
        idx = [];
    end

end

function s = local_capfirst(s)
%LOCAL_CAPFIRST Capitalize first character only.

    s = char(string(s));
    
    if isempty(s)
        return
    end
    
    s = [upper(s(1)) s(2:end)];
end

function s = local_latex_escape_title(s)
%LOCAL_LATEX_ESCAPE_TITLE Escapes a title for MATLAB's LaTeX interpreter.

    s = char(string(s));
    s = strrep(s,'\','\textbackslash{}');
    s = strrep(s,'_','\_');
    s = strrep(s,'%','\%');
    s = strrep(s,'&','\&');
    s = strrep(s,'#','\#');

end


% ----------------------------------------------------------------
function varargout = local_region_helpers_plot(op,region,varargin)
% ----------------------------------------------------------------
% Wrapper around region_helpers for plot_SAGE, with local fallbacks for
% newer CAMELS regions if an older region_helpers.m is still on the path.

    try
        [varargout{1:nargout}] = ...
            region_helpers(op,region,varargin{:});
        return
    catch
        % Fall through to local fallback below.
    end

end

% -------------------------------------
function s = local_basin_code(u,region)
% -------------------------------------
%LOCAL_BASIN_CODE Standardize basin codes for plotting/selection.

    if nargin < 2 ...
            || isempty(region)
        region = 'US';
    end

    reg = upper(strtrim(char(string(region))));

    % Convert input to clean string
    if isnumeric(u)
        if isempty(u) ...
                || isnan(u)
            s = '';
            return
        end
        s = sprintf('%.0f',double(u));
    else
        s = char(string(u));
    end

    s = strtrim(s);
    s = regexprep(s,'\.0+$','');

    if isempty(s)
        return
    end

    % CAMELS-CZ: use the bare uppercase gauge code internally.
    % Accept both O4263000 and camelscz_o4263000.
    if any(strcmp(reg,{'CAMELS_CZ','CZ','CZECHIA', ...
            'CZECH REPUBLIC'}))
        s = upper(strtrim(s));
        s = regexprep(s,'^CAMELSCZ[_-]?','');
        s = regexprep(s,'[^0-9A-Z]','');
        return
    end

    if any(strcmp(reg,{'CAMELS_CH', ...
            'CH','SWITZERLAND'}))

        s = char(string(s));
        s = strtrim(s);
        s = regexprep(s,'\D','');   % keep digits only
        % CAMELS-CH basin IDs are 4-digit codes.
        % Remove accidental US-style zero padding.
        if strlength(string(s)) > 4
            s = char(extractAfter(string(s), ...
                strlength(string(s))-4));
        end
        return
    end

    % CAMELS-DE: add/preserve DE prefix
    if any(strcmp(reg,{'CAMELS_DE', ...
        'DE','GERMANY'}))
    
        s = upper(strtrim(s));
        s = regexprep(s,'^DE','');
    
        % CAMELS-DE has mixed IDs:
        %   alphabetic-prefix IDs: G10410, C10160, A13080, ...
        %   numeric IDs:          00210550, 00710230, ...
        % Do not add DE. Match dat{k}.gauge exactly after removing DE.
        return
    end

    % CAMELS-IND: preserve 5-digit gauge IDs
    if any(strcmp(reg,{'CAMELS_IND', ...
            'IND','INDIA'}))
        n = str2double(s);
        if ~isnan(n)
            s = sprintf('%05d',round(n));
        end
        return
    end

    % CAMELS-LUX: use unpadded numeric IDs internally.  This accepts
    % 1, 01, ID01, and ID_01 as the same basin code.  The time-series
    % reader can still map this back to files named ID_01.
    if any(strcmp(reg,{'CAMELS_LUX', ...
            'LUX','LUXEMBOURG'}))
        s = upper(s);
        s = regexprep(s,'^ID_?','');
        n = str2double(s);
        if ~isnan(n)
            s = sprintf('%.0f',round(n));
        end
        return
    end

    % CAMELS-FI: preserve Finnish gauge IDs as strings.  Numeric gauges
    % stay unpadded, while virtual or compound gauge IDs with hyphens are
    % left unchanged.
    if any(strcmp(reg,{'CAMELS_FI', ...
            'FI','FIN','FINLAND'}))
        s = upper(s);
        n = str2double(s);
        if ~isnan(n)
            s = sprintf('%.0f',round(n));
        end
        return
    end

    % CAMELS-COL: preserve 8-digit IDEAM basin/gauge IDs.
    if any(strcmp(reg,{'CAMELS_COL', ...
            'COL','CO','COLOMBIA','COLUMBIA'}))
        s = upper(s);
        s = regexprep(s,'[^0-9A-Z]','');
        n = str2double(s);
        if ~isnan(n)
            s = sprintf('%.0f',round(n));
        end
        return
    end

    % CAMELS-AU: remove leading zeros only for purely numeric IDs.
    % Preserve alphanumeric IDs such as A2390531 and 422394A.
    if any(strcmp(reg,{'CAMELS_AU','AU','AUS','AUSTRALIA'}))
        s = upper(s);
        if ~isempty(regexp(s,'^\d+$','once'))
            s = regexprep(s,'^0+','');
            if isempty(s)
                s = '0';
            end
        end
        return
    end

    % CAMELS-US: preserve 8-digit USGS IDs
    if any(strcmp(reg,{'CAMELS_US','MACH_US','US', ...
            'USA','UNITED STATES'})) ...
            && all(isstrprop(s,'digit')) ...
            && numel(s) < 8
        s = [repmat('0',1,8-numel(s)) s];
        return
    end

    % CAMELS-SE: remove leading zeros
    if any(strcmp(reg,{'CAMELS_SE','SE','SWEDEN'}))
        s = regexprep(s,'^0+','');
        return
    end

end

function ytxt = local_scenario_ylabel(is,~)
%LOCAL_SCENARIO_YLABEL Returns the y-axis label for a scenario panel.

    switch is
        case 1
            ytxt = '\textbf{train basins | train period}';
        case 2
            ytxt = '\textbf{train basins | eval period}';
        case 3
            ytxt = '\textbf{eval basins | train period}';
        case 4
            ytxt = '\textbf{eval basins | eval period}';
    end

end

function [nrow,ncol] = local_parameter_map_layout(d)
%LOCAL_PARAMETER_MAP_LAYOUT Choose compact layout for parameter maps

    if d <= 4
        ncol = d;
    elseif d <= 8
        ncol = 4;
    elseif d <= 10
        ncol = 5;
    elseif d <= 15
        ncol = 5;     % HBV: 13 -> 3 x 5
    else
        ncol = 5;
    end
    
    nrow = ceil(d/ncol);
end

function local_add_parameter_colorbar_in_grid(figNo,blueBase, ...
    leftMargin,colGap,rowGap,axW,axH,nrow,ncol,d)
%LOCAL_ADD_PARAMETER_COLORBAR_IN_GRID Thin colorbar in empty map panels

    figure(figNo);

    topMargin = 0.055;

    firstEmpty = d + 1;
    row = ceil(firstEmpty/ncol);
    col = mod(firstEmpty-1,ncol) + 1;

    nEmpty = nrow*ncol - d;

    % Use at most two empty panels
    span = min([2,nEmpty,ncol-col+1]);

    x0 = leftMargin + (col-1)*(axW + colGap);
    y0 = 1 - topMargin - row*axH - (row-1)*rowGap;

    cbW = span*axW + (span-1)*colGap;
    cbH = 0.018;              % thin colorbar
    cbX = x0 + 0.08*cbW;
    cbY = y0 + 0.48*axH;
    cbW = 0.84*cbW;

    n = 256;
    v = linspace(0,1,n);

    CB = zeros(n,3);
    for i = 1:n
        CB(i,:) = local_tint_color(v(i),0,1,blueBase);
    end

    ax1 = axes('Position',[cbX cbY cbW cbH]);
    hold(ax1,'on');

    for i = 1:n-1
        patch(ax1,[v(i) v(i+1) v(i+1) v(i)], ...
            [0 0 1 1],CB(i,:), ...
            'EdgeColor','none');
    end

    set(ax1, ...
        'XLim',[0 1], ...
        'YLim',[0 1], ...
        'YTick',[], ...
        'XTick',0:0.2:1, ...
        'XTickLabel',{'0','','','','','1'}, ... %{'0','0.2','0.4','0.6','0.8','1.0'}, ...
        'tickdir','out', ...
        'XAxisLocation','bottom', ...
        'fontname','Calibri', ...
        'fontsize',17, ...
        'box','off', ...
        'Layer','top');

    % Manually draw a box
    set(ax1,'Box','off','TickDir','out','XAxisLocation','bottom');
    line(ax1,[0 1],[0 0],'Color','k','LineWidth',0.75);
    line(ax1,[0 1],[1 1],'Color','k','LineWidth',0.75);
    line(ax1,[0 0],[0 1],'Color','k','LineWidth',0.75);
    line(ax1,[1 1],[0 1],'Color','k','LineWidth',0.75);
    % End manually draw a box

    xlabel(ax1,'normalized parameter value', ...
        'fontname','Calibri', ...
            'fontsize',12);
end

function [latLim,lonLim] = local_parameter_map_limits(latlon,region)
%LOCAL_PARAMETER_MAP_LIMITS Return tight map limits around selected basins

    lat = latlon(:,1);
    lon = latlon(:,2);
    
    good = isfinite(lat) ...
        & isfinite(lon);
    lat = lat(good);
    lon = lon(good);
    
    if isempty(lat)
        latLim = [-90 90];
        lonLim = [-180 180];
        return
    end
    
    padLat = 0.04 * range(lat);
    padLon = 0.04 * range(lon);
    
    if padLat == 0
        padLat = 0.5;
    end
    if padLon == 0
        padLon = 0.5; 
    end
    
    latLim = [min(lat)-padLat max(lat)+padLat];
    lonLim = [min(lon)-padLon max(lon)+padLon];
    
    % Optional CONUS protection for US maps
    r = upper(string(region));
    if any(r == ["US","CAMELS_US","MACH_US"])
        latLim = [max(latLim(1),24) min(latLim(2),50)];
        lonLim = [max(lonLim(1),-125) min(lonLim(2),-66)];
    end
end
