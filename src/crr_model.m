function [loss_value,out] = crr_model(x,mdl,dat,ode,loss,request)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRR_MODEL Run the crr_model and return the total loss and sensitivity
% vectors and/or matrices.
%
% SYNOPSIS: [loss_value,out] = crr_model(x,mdl,dat,ode,loss,request)
%   x           dx1 vector of parameter values [check pspace]
%   mdl         structure with model state/parameter info
%    .model      choice of model 
%                 1 hymod
%                 2 hmodel
%                 3 sacsma
%                 4 xinanjiang
%                 5 gr4j
%                 6 hbv
%                 7 cfe_nwm
%                 8 user-specified model
%                21 gr4jB [analytic routing]
%                31 hymod data-assimilation implementation
%    .mcode      scalar with numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation sacsma_ode C++
%    .calc       model execution
%      'seq'      sequential execution of watersheds
%      'par'      parallel execution of watersheds
%    .mode       assessment design
%                 1 = training basins | training period
%                 2 = training basins | evaluation period/mask
%                 3 = training + evaluation basins | training period
%                 4 = training + evaluation basins | evaluation period/mask
%    .sp_method  string with split design
%                 'manual'
%                 'deterministic_block'
%                 'random_block'
%                 'random'
%                 'deterministic_kfold'
%                 'random_kfold'
%    .names      list of model names
%    .y0         mx1 vector of initial states
%    .pspace     0:hydrologic, 1:unit cube, 2:unconstrained parameters
%    .th_min     dx1 vector of lower parameter values [= in pspace]
%    .th_max     dx1 vector of upper parameter values [= in pspace]
%    .par_names  1xd cell parameter names
%    .tout       1x(n+m+1)+spinup*365 vector model output times
%    .idx        1x(n+m+1) vector of indices train&val periods
%    .id_train   1xn vector of indices training period
%    .id_eval    1xm vector of indices evaluation period
%   dat         structure with meteo, discharge and cache data
%   ode         structure with settings for ODE solver
%    .InitStep   Initial time step
%    .MaxStep    Maximum time step
%    .MinStep    Minimum time step
%    .RelTol     Relative tolerance
%    .AbsTol     Absolute tolerance
%    .Order      Order
%    .maxiter    Maximum number of iterations
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
%    .meta       metadata structure [method = 4]
%   request     logical output-selection structure
%    .q          return simulated discharge
%    .gradient   return the loss gradient
%    .jacobian   return the discharge Jacobian on the training mask
%    .metrics    return compact train/evaluation metrics
%    .attribution return total and net parameter attribution
%    .states     []/false, "all", or requested state-name strings
%   out         structure containing only the requested results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pspace = 2: vartheta = (-∞,∞)
% pspace = 1: ntheta = [0,1]
% pspace = 0: theta = tabulated ranges [= hydrologic parameter values]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6
        request = crr_request();
    elseif ~(isstruct(request) && isscalar(request) ...
            && isfield(request,'normalized') ...
            && isequal(request.normalized,true))
        request = crr_request(request);
    end
    out = struct();
    
    d = numel(x);                       % # parameters
    
    if isfield(mdl,'local') ...         % based on yearly rainfall ranking
            && mdl.local == 1
        id_tra = dat.id_train;          % all training-period entries
        id_eva = dat.id_eval;           % all evaluation-period entries
    else
        id_tra = expand_index( ...      % all training-period entries
            mdl.id_train);
        id_eva = expand_index( ...      % all evaluation-period entries
            mdl.id_eval);
    end
    
    good_tra = ~dat.bad(id_tra);        % good entries [=missing values remvd]
    good_eva = ~dat.bad(id_eva);        % good entries [=missing values remvd]
    id_tr = id_tra(good_tra);           % valid training entries
    id_ev = id_eva(good_eva);           % valid evaluation entries
    
    model = mdl.model;                  % choice of model [= integer]
    loss_fnc = loss.fnc;                % loss function [= integer]
    
    needGradient = request.gradient ...
        || request.attribution;
    needJ = request.jacobian ...
        || needGradient;
    needStates = ~isempty(request.states);
    Z = [];
    ode_model = ode;
    if needStates
        ode_model.mem = 1;
    end

    switch model
        case 1
            if needStates
                [q,J,~,Z] = hymod(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = hymod(x,mdl,dat.meteo,ode_model);
            else
                q = hymod(x,mdl,dat.meteo,ode_model);
            end
        case 2
            if needStates
                [q,J,~,Z] = hmodel(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = hmodel(x,mdl,dat.meteo,ode_model);
            else
                q = hmodel(x,mdl,dat.meteo,ode_model);
            end
        case 3
            if needStates
                [q,J,~,Z] = sacsma(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = sacsma(x,mdl,dat.meteo,ode_model);
            else
                q = sacsma(x,mdl,dat.meteo,ode_model);
            end
        case 4
            if needStates
                [q,J,~,Z] = xinanjiang(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = xinanjiang(x,mdl,dat.meteo,ode_model);
            else
                q = xinanjiang(x,mdl,dat.meteo,ode_model);
            end
        case 5
            if needStates
                [q,J,~,Z] = gr4jA(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = gr4jA(x,mdl,dat.meteo,ode_model);
            else
                q = gr4jA(x,mdl,dat.meteo,ode_model);
            end
        case 6
            if needStates
                [q,J,~,Z] = hbv(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = hbv(x,mdl,dat.meteo,ode_model);
            else
                q = hbv(x,mdl,dat.meteo,ode_model);
            end
        case 7
            if needStates
                [q,J,~,Z] = cfe_nwm(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = cfe_nwm(x,mdl,dat.meteo,ode_model);
            else
                q = cfe_nwm(x,mdl,dat.meteo,ode_model);
            end
        case 8
            if needStates
                [q,J,~,Z] = user_model(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = user_model(x,mdl,dat.meteo,ode_model);
            else
                q = user_model(x,mdl,dat.meteo,ode_model);
            end
        case 21
            if needStates
                [q,J,~,Z] = gr4jB(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = gr4jB(x,mdl,dat.meteo,ode_model);
            else
                q = gr4jB(x,mdl,dat.meteo,ode_model);
            end
        case 31
            if needStates
                [q,J,~,Z] = hymod2(x,mdl,dat.meteo,dat.y_n,ode_model);
            elseif needJ
                [q,J] = hymod2(x,mdl,dat.meteo,dat.y_n,ode_model);
            else
                q = hymod2(x,mdl,dat.meteo,dat.y_n,ode_model);
            end
        otherwise
            error('I do not know this model');
    end
    % ----------------------------
    % Training and evaluation data
    % ----------------------------
    if needJ
        J = J(id_tr,1:d);               % Jacobian on training mask only
    end
    y_t = dat.y_n(id_tr);               % nx1 observed discharge, training
    q_t = q(id_tr);                     % nx1 simulated discharge, training
    res_t = y_t - q_t;                  % nx1 residuals, training
    
    y_e = dat.y_n(id_ev);               % mx1 observed discharge, evaluation
    q_e = q(id_ev);                     % mx1 simulated discharge, evaluation
    res_e = y_e - q_e;                  % mx1 residuals, evaluation
    
    % --------------------------------
    % Optimization loss: training only
    % --------------------------------
    [JKGE_Mt,JKGE_Vt,JKGE_Ct] = deal(NaN);
    Dfdct = NaN;
    Dfdce = NaN; %#ok

    if isempty(res_t)
        loss_value = NaN;
        delta = zeros(size(res_t));
        JKGEt = NaN;
        m_q = NaN(size(dat.y_n));
    else
        switch loss_fnc
            case 1 % SAR
                loss_value = sum(abs(res_t));
            case 2 % GLS / RSS
                loss_value = sum(res_t.^2);
            case 3 % 1 - NSE
                if isfinite(dat.stats.TSSt) ...
                        && dat.stats.TSSt > 0
                    loss_value = sum(res_t.^2) / ...
                        dat.stats.TSSt;
                else
                    loss_value = NaN;
                end
            case 4 % 1 - KGE
                [KGEt,~,~,~] = kge(y_t, ...
                    q_t,dat.stats.mut, ...
                    dat.stats.stdt);
                loss_value = 1 - KGEt;
            case 5 % Huber loss
                [loss_value,delta] = ...
                    huber_loss(res_t, ...
                    dat.stats.Syt);
            case 6 % FDC
%                loss_value = fdc_loss(y_t,q_t);
                if numel(q_t) ~= dat.fdc.t.n
                    error(['      Error:crr_model: ' ...
                        'FDC training cache length mismatch.']);
                end
                Dfdct = fdc_loss_cached(q_t,dat.fdc.t);
                loss_value = Dfdct;
            case 7 % 1 - JKGE
                % Prepare input arguments
                if isfield(loss,'M') ...
                        && ~isempty(loss.M)
                    Mdef = loss.M;
                else
                    Mdef = 2;   % or 1 if you want the paper version as default
                end
                args_t = jkge_args( ...
                    loss,'all');
                if ~needGradient
                    [JKGEt,m_q,JKGE_Mt, ...
                        JKGE_Vt,JKGE_Ct] = jkge( ...
                        dat.y_n, ...
                        q, ...
                        dat.jkge.m_y, ...
                        id_tr, ...
                        args_t{:});
                else
                    switch loss.method
                        case {1,2}
                            aux = loss.n_win;
                        case 3
                            aux = [];
                        case 4
                            aux = loss.meta.mo_all;
                        otherwise
                            error('Unknown JKGE method.');
                    end
                    [JKGEt,delta_raw,m_q, ...
                        JKGE_Mt,JKGE_Vt, ...
                        JKGE_Ct] = ...
                        jkge_grad( ...
                        dat.y_n, ...
                        q, ...
                        dat.jkge.m_y, ...
                        id_tr, ...
                        loss.method, ...
                        aux, ...
                        dat.jkge.cache, ...
                        Mdef);
                    if numel(delta_raw) ...
                            == numel(dat.y_n)
                        delta = delta_raw(id_tr);
                    elseif numel(delta_raw) ...
                            == numel(id_tr)
                        delta = delta_raw;
                    else
                        error(['      ' ...
                            'Error: crr_model: ' ...
                            'JKGE gradient ' ...
                            'length mismatch.']);
                    end
                end
                loss_value = 1 - JKGEt; 

            otherwise
                error('Unknown loss function choice');
        end
    end

    % --------------------
    % Gradient computation
    % --------------------
    switch loss_fnc
        case 2 % GLS: Sigma_eps;
            args = {1};      
        otherwise 
            args = {};
    end

    if needGradient
        if isempty(res_t)
            delta = zeros(size(res_t));
            g = nan(d,1);
        else
            if loss_fnc == 5 ...
                    || loss_fnc == 7
                % delta already computed
            else
                delta = delta_n(loss_fnc, ...
                    y_t,q_t,args{:});
            end
            g = J' * delta;
        end
    end
    % ----------------
    % Return arguments
    % ----------------
    if request.q
        out.q = q; 
    end
    if request.gradient
        out.gradient = g;
    end
    if request.jacobian
        out.jacobian = J;
    end
    if ~isempty(request.states)
        out.states = model_states(mdl, ...
            Z,request.states,numel(q));
    end
    if request.metrics
        % -----------------------------------------------------
        % FDC distance: compute for current simulated discharge
        % -----------------------------------------------------
        if ~isfinite(Dfdct)
            Dfdct = fdc_loss_cached(q_t,dat.fdc.t);
        end

        if ~isempty(q_e) ...
                && isfield(dat,'fdc') ...
                && isfield(dat.fdc,'e') ...
                && dat.fdc.e.n > 0
            Dfdce = fdc_loss_cached(q_e,dat.fdc.e);
        else
            Dfdce = NaN;
        end

        met = local_metrics_struct(y_t,q_t, ...
            res_t,dat.stats.mut, ...
            dat.stats.stdt,dat.stats.TSSt, ...
            dat.stats.Syt,y_e,q_e,res_e, ...
            dat.stats.mue,dat.stats.stde, ...
            dat.stats.TSSe,dat.stats.Sye);

        met.Dfdct = Dfdct;
        met.Dfdce = Dfdce;
        
        if loss_fnc == 7
            met.JKGEt = JKGEt;
            met.JKGE_Mt = JKGE_Mt;
            met.JKGE_Vt = JKGE_Vt;
            met.JKGE_Ct = JKGE_Ct;
            if all(isfinite(m_q(id_ev)))
                [met.JKGEe,met.JKGE_Me, ...
                    met.JKGE_Ve,met.JKGE_Ce] = ...
                    jkge_score_given_mq( ...
                    dat.y_n,q,dat.jkge.m_y,m_q, ...
                    id_ev,loss);
            else
                met.JKGEe = NaN;
                [met.JKGE_Me,met.JKGE_Ve, ...
                    met.JKGE_Ce] = deal(NaN);
            end
        end
        out.metrics = met;
    end
    if request.attribution
        [At,An] = sage_attribution(J,delta,mdl,g);
        out.attribution = struct('total',At,'net',An);
    end
end

% =============
% local helpers
% =============
function met = local_metrics_struct( ...
    y_t,q_t,res_t,mu_t,std_t,TSSt,S_y_t, ...
    y_e,q_e,res_e,mu_e,std_e,TSSe,S_y_e)
%LOCAL_METRICS_STRUCT Compute compact train/evaluation metrics structure.
%
% Computes standard scalar performance metrics for the training and
% evaluation periods and stores them in a single structure. Metrics are
% computed separately for the training data and evaluation data.
%
% INPUT
%   y_t     nt x 1 observed discharge, training period
%   q_t     nt x 1 simulated discharge, training period
%   res_t   nt x 1 residuals, training period, y_t - q_t
%   mu_t    scalar mean of observed training discharge
%   std_t   scalar standard deviation of observed training discharge
%   TSSt    scalar total sum of squares, training period
%   S_y_t   scalar Huber scale training observations
%   y_e     ne x 1 observed discharge, evaluation period
%   q_e     ne x 1 simulated discharge, evaluation period
%   res_e   ne x 1 residuals, evaluation period, y_e - q_e
%   mu_e    scalar mean of observed evaluation discharge
%   std_e   scalar standard deviation of observed evaluation discharge
%   TSSe    scalar total sum of squares, evaluation period
%   S_y_e   scalar Huber scale evaluation observations
%
% OUTPUT
%   met     structure with fields:
%           SARt, GLSt, NSEt, KGEt, Hubert, RSSt, JKGEt
%           SARe, GLSe, NSEe, KGEe, Hubere, RSSe, JKGEe
%
% NOTES
%   - Suffix "t" denotes training-period metrics.
%   - Suffix "e" denotes evaluation-period metrics.
%   - JKGE is initialized by local_metric_block and may be overwritten
%     later when JKGE-specific benchmark vectors are available.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, April 2026                                %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [SARt,GLSt,NSEt,KGEt,KGE_rt,KGE_alphat, ...
        KGE_betat,Hubert,RSSt,JKGEt] = ...
        local_metric_block(y_t,q_t, ...
        res_t,mu_t,std_t,TSSt,S_y_t);
    
    [SARe,GLSe,NSEe,KGEe,KGE_re,KGE_alphae, ...
        KGE_betae,Hubere,RSSe,JKGEe] = ...
        local_metric_block(y_e,q_e, ...
        res_e,mu_e,std_e,TSSe,S_y_e);
    
    met = struct( ...
        'SARt',SARt, ...
        'GLSt',GLSt, ...
        'NSEt',NSEt, ...
        'KGEt',KGEt, ...
        'KGE_rt',KGE_rt, ...
        'KGE_alphat',KGE_alphat, ...
        'KGE_betat',KGE_betat, ...
        'Hubert',Hubert, ...
        'RSSt',RSSt, ...
        'JKGEt',JKGEt, ...
        'JKGE_Mt',NaN, ...
        'JKGE_Vt',NaN, ...
        'JKGE_Ct',NaN, ...
        'SARe',SARe, ...
        'GLSe',GLSe, ...
        'NSEe',NSEe, ...
        'KGEe',KGEe, ...
        'KGE_re',KGE_re, ...
        'KGE_alphae',KGE_alphae, ...
        'KGE_betae',KGE_betae, ...
        'Hubere',Hubere, ...
        'RSSe',RSSe, ...
        'JKGEe',JKGEe, ...
        'JKGE_Me',NaN, ...
        'JKGE_Ve',NaN, ...
        'JKGE_Ce',NaN);
end

function [SAR,GLS,NSE,KGE,KGE_r,KGE_alpha, ...
    KGE_beta,Huber,RSS,JKGE] = ...
    local_metric_block(y,q,res,mu_y,std_y, ...
    TSS,S_y)
%LOCAL_METRIC_BLOCK Compute standard discharge performance metrics.
%
% Computes a compact set of scalar performance metrics for one basin and
% one selected period using observed discharge y, simulated discharge q,
% and residuals res = y - q.
%
% INPUT
%   y       nx1 observed discharge vector
%   q       nx1 simulated discharge vector
%   res     nx1 residual vector, y - q
%   mu_y    scalar mean of observed discharge
%   std_y   scalar standard deviation of observed discharge
%   TSS     scalar total sum of squares of observed discharge
%   S_y     scalar robust Huber scale of observed discharge
%
% OUTPUT
%   SAR     sum of absolute residuals
%   GLS     generalized least-squares loss proxy; currently RSS
%   NSE     Nash-Sutcliffe efficiency
%   KGE     Kling-Gupta efficiency
%   Huber   Huber loss
%   RSS     residual sum of squares
%   JKGE    Jawad Kling-Gupta efficiency; NaN unless computed elsewhere
%
% NOTES
%   - Nonfinite paired values are removed before computing metrics.
%   - If TSS is missing, nonfinite or nonpositive, it is recomputed from y.
%   - JKGE is returned as NaN in this helper because JKGE requires the
%     benchmark vectors m_y and m_q.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, April 2026                                %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(res)
        [SAR,GLS,NSE,KGE,KGE_r,KGE_alpha, ...
            KGE_beta,Huber,RSS,JKGE] = ...
            deal(NaN);
        return
    end
    
    % Ensure column vectors
    y = y(:);
    q = q(:);
    res = res(:);
    
    % Keep only finite paired values for metric calculations
    good = isfinite(y) ...
        & isfinite(q) ...
        & isfinite(res);
    y = y(good);
    q = q(good);
    res = res(good);
    
    if isempty(y)
        [SAR,GLS,NSE,KGE,KGE_r,KGE_alpha, ...
            KGE_beta,Huber,RSS,JKGE] = ...
            deal(NaN);
        return
    end
    
    RSS = sum(res.^2);
    SAR = sum(abs(res));
    GLS = RSS;
    
    % Robust scalar TSS
    if ~(isscalar(TSS) ...
            && isfinite(TSS) ...
            && (TSS > 0))
        ybar = mean(y);
        TSS = sum((y - ybar).^2);
    end
    
    if isfinite(TSS) ...
            && (TSS > 0)
        NSE = 1 - RSS / TSS;
    else
        NSE = NaN;
    end

    [KGE,KGE_r,KGE_alpha,KGE_beta] = ...
        kge(y,q,mu_y,std_y);
    [Huber,~] = huber_loss(res,S_y);
    JKGE = NaN;
end
