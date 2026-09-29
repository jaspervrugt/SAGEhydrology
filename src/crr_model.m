function [loss_value,out] = crr_model(x,mdl,dat,ode,loss,request)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRR_MODEL Evaluate one rainfall-runoff model and requested outputs.
%
%  Run the selected conceptual rainfall-runoff model for one basin.
%  Compute its training-period loss and return only the requested
%  simulated trajectories, derivatives, metrics, or attributions.
%
% SYNOPSIS:
%   [loss_value,out] = crr_model(x,mdl,dat,ode,loss,request)
%   [loss_value,out] = crr_model(x,mdl,dat,ode,loss)
%
% INPUT ARGUMENTS:
%   x               d-by-1 parameter vector in the space selected by
%                   mdl.pspace
%   mdl             model configuration and assessment indices
%    .model          model code:
%                    1 HYMOD; 2 HMODEL; 3 SAC-SMA; 4 Xinanjiang
%                    5 GR4J; 6 HBV; 7 CFE-NWM; 11 GCHM; 99 user model
%    .pspace         parameter space: 0 physical, 1 normalized, 2
%                    unconstrained
%    .id_train       training-period indices or masks
%    .id_eval        evaluation-period indices or masks
%    .local          optional flag for basin-local indices in dat
%    .variant        optional GCHM implementation variant
%    .swe_ind        optional SWE state mapping for named requests
%    .sm_ind         optional soil-moisture state mapping
%   dat             basin data and precomputed observation summaries
%    .meteo          meteorological forcing passed to the selected model
%    .obs            observed variables
%     .Q              discharge observations
%      .value          observed discharge series
%      .bad            invalid-discharge mask
%    .stats          training/evaluation discharge summary statistics
%    .fdc            cached training/evaluation discharge distributions
%    .jkge           JKGE benchmark data when loss.fnc = 7
%   ode             numerical solver settings passed to the model
%   loss            loss-function configuration
%    .fnc            1 SAR; 2 squared residuals; 3 (1-NSE); 4 (1-KGE)
%                    5 Huber; 6 FDC divergence; 7 JKGE
%    .fdc            FDC settings when loss.fnc = 6
%     .formulation    1 physical-Q; 2 probability; 3 log-probability
%    .method         JKGE benchmark method when loss.fnc = 7
%    .n_win          JKGE moving-average window when required
%    .meta           JKGE metadata when required
%    .M              optional JKGE norm definition
%   request         optional output-selection structure; default is empty
%    .q              return the complete simulated-discharge series
%    .gradient       return the training-loss gradient
%    .jacobian       return the discharge Jacobian on valid training data
%    .metrics        return training/evaluation performance metrics
%    .attribution    return total and net gradient attribution
%    .states         return selected model-state trajectories
%    .obs            request named trajectories: Q, SWE, or SM
%    .jac            request named Jacobians: Q, SWE, or SM
%
% OUTPUT ARGUMENTS:
%   loss_value      scalar training loss; NaN for named-channel requests
%   out             structure of requested outputs
%    .q              complete simulated-discharge series, when requested
%    .gradient       d-by-1 gradient in the parameter space of x
%    .jacobian       valid-training discharge Jacobian, n-by-d
%    .metrics        training/evaluation performance metrics
%    .attribution    parameter-gradient attribution
%     .total          d-by-1 total attribution
%     .net            d-by-1 net attribution
%    .states         requested model-state trajectories
%    .obs            requested named simulated trajectories
%     .Q              discharge trajectory, when requested
%     .SWE            snow water equivalent, when requested
%     .SM             soil moisture, when requested
%    .jac            requested named Jacobians
%     .Q              discharge Jacobian, when requested
%     .SWE            SWE Jacobian, when requested
%     .SM             soil-moisture Jacobian, when requested
%
% NOTES:
%   The discharge Jacobian and gradient are scaled to the parameter
%   space supplied in x. Missing discharge observations are excluded
%   from the training loss and its Jacobian.
%   A nonempty request.obs or request.jac takes the named-channel path:
%   the training loss is not evaluated, and loss_value is NaN.
%   That path also returns zero gradient and empty metrics placeholders.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6
        request = crr_request();
    elseif ~(isstruct(request) ...
            && isscalar(request) ...
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

    % Named requests use the structured model interface directly. Q is
    % not requested implicitly when another observation is trained.
    if ~isempty(request.obs) || ~isempty(request.jac)
        [loss_value,out] = local_named_evaluation( ...
            x,mdl,dat,ode,request,d);
        return
    end
    
    good_tra = ~dat.obs.Q.bad(id_tra);  % good entries [=missing vls remvd]
    good_eva = ~dat.obs.Q.bad(id_eva);  % good entries [=missing vls remvd]
    id_tr = id_tra(good_tra);           % valid training entries
    id_ev = id_eva(good_eva);           % valid evaluation entries
    
    model = mdl.model;                  % choice of model [= integer]
    loss_fnc = loss.fnc;                % loss function [= integer]
    
    needGradient = request.gradient ...
        || request.attribution;
    namedObservables = upper(string(request.obs(:)));
    namedJacobians = upper(string(request.jac(:)));
    needNamedQJacobian = any(namedJacobians == "Q");
    needNamedSWE = any(namedObservables == "SWE") ...
        || any(namedJacobians == "SWE");
    needNamedSM = any(namedObservables == "SM") ...
        || any(namedJacobians == "SM");
    needJ = request.jacobian ...
        || needGradient ...
        || needNamedQJacobian ...
        || any(namedJacobians == "SWE") ...
        || any(namedJacobians == "SM");
    needStates = ~isempty(request.states) ...
        || needNamedSWE ...
        || needNamedSM;
    Z = [];
    Jth = [];
    ode_model = ode;
    if needStates
        ode_model.mem = 1;
    end

    switch model
        case 1
            if needStates
                [q,J,Jth,Z] = hymod(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = hymod(x,mdl,dat.meteo,ode_model);
            else
                q = hymod(x,mdl,dat.meteo,ode_model);
            end
        case 2
            if needStates
                [q,J,Jth,Z] = hmodel(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = hmodel(x,mdl,dat.meteo,ode_model);
            else
                q = hmodel(x,mdl,dat.meteo,ode_model);
            end
        case 3
            if needStates
                [q,J,Jth,Z] = sacsma(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = sacsma(x,mdl,dat.meteo,ode_model);
            else
                q = sacsma(x,mdl,dat.meteo,ode_model);
            end
        case 4
            if needStates
                [q,J,Jth,Z] = xinanjiang(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = xinanjiang(x,mdl,dat.meteo,ode_model);
            else
                q = xinanjiang(x,mdl,dat.meteo,ode_model);
            end
        case 5
            if needStates
                [q,J,Jth,Z] = gr4jA(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = gr4jA(x,mdl,dat.meteo,ode_model);
            else
                q = gr4jA(x,mdl,dat.meteo,ode_model);
            end
        case 6
            if needStates
                [q,J,Jth,Z] = hbv(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = hbv(x,mdl,dat.meteo,ode_model);
            else
                q = hbv(x,mdl,dat.meteo,ode_model);
            end
        case 7
            if needStates
                [q,J,Jth,Z] = cfe_nwm(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = cfe_nwm(x,mdl,dat.meteo,ode_model);
            else
                q = cfe_nwm(x,mdl,dat.meteo,ode_model);
            end
        case 11
            isOde = isfield(mdl,'variant') ...
                && strcmpi(string(mdl.variant),'gchm_ode');
            if isOde
                modelFcn = @gchm_ode;
            else
                modelFcn = @gchm;
            end
            if exist(func2str(modelFcn),'file') ~= 2
                error('crr_model:GCHMUnavailable', ...
                    ['gchm (model 11) is not installed. Add ' ...
                    'private/gchm to the MATLAB path.']);
            end
            
            if needStates
                [q,J,Jth,Z] = modelFcn(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = modelFcn(x,mdl,dat.meteo,ode_model);
            else
                q = modelFcn(x,mdl,dat.meteo,ode_model);
            end
        case 99
            if needStates
                [q,J,Jth,Z] = user_model(x,mdl,dat.meteo,ode_model);
            elseif needJ
                [q,J] = user_model(x,mdl,dat.meteo,ode_model);
            else
                q = user_model(x,mdl,dat.meteo,ode_model);
            end
        otherwise
            error('I do not know this model');
    end
    % ----------------------------
    % Training and evaluation data
    % ----------------------------
    if needJ
        Jfull = J(:,1:d);
        J = J(id_tr,1:d);               % Jacobian on training mask only
    end
    y_t = dat.obs.Q.value(id_tr);       % nx1 observed discharge, training
    q_t = q(id_tr);                     % nx1 simulated discharge, training
    res_t = y_t - q_t;                  % nx1 residuals, training
    
    y_e = dat.obs.Q.value(id_ev);       % mx1 observed discharge, evaluation
    q_e = q(id_ev);                     % mx1 simulated discharge, evaluation
    res_e = y_e - q_e;                  % mx1 residuals, evaluation
    
    % --------------------------------
    % Optimization loss: training only
    % --------------------------------
    [JKGE_Mt,JKGE_Vt,JKGE_Ct] = deal(NaN);
    [Dfdct,Dpt,Dlogpt] = deal(NaN);
    
    if isempty(res_t)
        loss_value = NaN;
        delta = zeros(size(res_t));
        JKGEt = NaN;
        m_q = NaN(size(dat.obs.Q.value));
    else
        switch loss_fnc
            case 1 % SAR
                loss_value = sum(abs(res_t));
            case 2 % GLS / RSS
                loss_value = sum(res_t.^2);
            case 3 % 1 - NSE
                if isfinite(dat.stats.Q.train.TSS) ...
                        && dat.stats.Q.train.TSS > 0
                    loss_value = sum(res_t.^2) / ...
                        dat.stats.Q.train.TSS;
                else
                    loss_value = NaN;
                end
            case 4 % 1 - KGE
                [KGEt,~,~,~] = kge(y_t, ...
                    q_t,dat.stats.Q.train.mean, ...
                    dat.stats.Q.train.std);
                loss_value = 1 - KGEt;
            case 5 % Huber loss
                [loss_value,delta] = ...
                    huber_loss(res_t, ...
                    dat.stats.Q.train.huber_scale);
            case 6 % FDC divergence: 6a d_fdc, 6b d_p, or 6c d_logp
                if numel(q_t) ~= dat.fdc.Q.t.n
                    error(['      Error:crr_model: ' ...
                        'FDC training cache length mismatch.']);
                end
                [Dfdct,Dpt,Dlogpt] = ...
                    fdc_metrics_cached(q_t,dat.fdc.Q.t);
                fdcForm = local_fdc_formulation(loss);
                values = [Dfdct,Dpt,Dlogpt];
                loss_value = values(fdcForm);
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
                        dat.obs.Q.value, ...
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
                        dat.obs.Q.value, ...
                        q, ...
                        dat.jkge.m_y, ...
                        id_tr, ...
                        loss.method, ...
                        aux, ...
                        dat.jkge.cache, ...
                        Mdef);
                    if numel(delta_raw) ...
                            == numel(dat.obs.Q.value)
                        delta = delta_raw;
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
        case 6 % selected FDC formulation
            args = {local_fdc_formulation(loss)};
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
            if loss_fnc == 7 ...
                    && numel(delta) == size(Jfull,1)
                g = Jfull' * delta;
            else
                g = J' * delta;
            end
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
    if ~isempty(namedObservables)
        out.obs = struct();
        if any(namedObservables == "Q")
            out.obs.Q = q;
        end
        if any(namedObservables == "SWE")
            out.obs.SWE = local_swe_from_history(mdl,Z,d,numel(q));
        end
        if any(namedObservables == "SM")
            out.obs.SM = local_sm_from_history(mdl,Z,d,numel(q));
        end
    end
    if ~isempty(namedJacobians)
        out.jac = struct();
        if needNamedQJacobian
            out.jac.Q = Jfull;
        end
        if any(namedJacobians == "SWE")
            [~,Jswe] = local_swe_from_history(mdl,Z,d,numel(q),Jth);
            out.jac.SWE = Jswe;
        end
        if any(namedJacobians == "SM")
            [~,Jsm] = local_sm_from_history(mdl,Z,d,numel(q),Jth);
            out.jac.SM = Jsm;
        end
    end
    if ~isempty(request.states)
        out.states = model_states(mdl, ...
            Z,request.states,numel(q));
    end
    if request.metrics
        % -----------------------------------------------------
        % FDC distance: compute for current simulated discharge
        % -----------------------------------------------------
        if ~all(isfinite([Dfdct,Dpt,Dlogpt]))
            [Dfdct,Dpt,Dlogpt] = ...
                fdc_metrics_cached(q_t,dat.fdc.Q.t);
        end

        if ~isempty(q_e) ...
                && isfield(dat,'fdc') ...
                && isfield(dat.fdc,'Q') ...
                && isfield(dat.fdc.Q,'e') ...
                && dat.fdc.Q.e.n > 0
            [Dfdce,Dpe,Dlogpe] = ...
                fdc_metrics_cached(q_e,dat.fdc.Q.e);
        else
            [Dfdce,Dpe,Dlogpe] = deal(NaN);
        end

        met = local_metrics_struct(y_t,q_t, ...
            res_t,dat.stats.Q.train.mean, ...
            dat.stats.Q.train.std,dat.stats.Q.train.TSS, ...
            dat.stats.Q.train.huber_scale,y_e,q_e,res_e, ...
            dat.stats.Q.eval.mean,dat.stats.Q.eval.std, ...
            dat.stats.Q.eval.TSS,dat.stats.Q.eval.huber_scale);

        met.Dfdct = Dfdct;
        met.Dfdce = Dfdce;
        met.Dpt = Dpt;
        met.Dpe = Dpe;
        met.Dlogpt = Dlogpt;
        met.Dlogpe = Dlogpe;
        
        if loss_fnc == 7
            met.JKGEt = JKGEt;
            met.JKGE_Mt = JKGE_Mt;
            met.JKGE_Vt = JKGE_Vt;
            met.JKGE_Ct = JKGE_Ct;
            if all(isfinite(m_q(id_ev)))
                [met.JKGEe,met.JKGE_Me, ...
                    met.JKGE_Ve,met.JKGE_Ce] = ...
                    jkge_score_given_mq( ...
                    dat.obs.Q.value,q,dat.jkge.m_y,m_q, ...
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
        if loss_fnc == 7 ...
                && numel(delta) == size(Jfull,1)
            [At,An] = sage_attribution(Jfull,delta,mdl,g);
        else
            [At,An] = sage_attribution(J,delta,mdl,g);
        end
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

function formulation = local_fdc_formulation(loss)
%LOCAL_FDC_FORMULATION Return 1/2/3 for loss 6a/6b/6c.

    formulation = 1;
    if isfield(loss,'fdc') ...
            && isstruct(loss.fdc) ...
            && isfield(loss.fdc,'formulation') ...
            && ~isempty(loss.fdc.formulation)
        formulation = double(loss.fdc.formulation);
    end
    if ~isscalar(formulation) ...
            || ~isfinite(formulation) ...
            || ~ismember(formulation,1:3)
        error('crr_model:BadFDCFormulation', ...
            ['loss.fdc.formulation ' ...
            'must be 1 (d_fdc), 2 (d_p), ' ...
            'or 3 (d_logp).']);
    end
end

function [loss_value,out] = local_named_evaluation( ...
    x,mdl,dat,ode,request,d)
%LOCAL_NAMED_EVALUATION Run only explicitly requested model channels.

    modelRequest = struct('obs',request.obs, ...
        'jac',request.jac,'states',request.states);
    switch mdl.model
        case 1
            result = hymod(x,mdl,dat.meteo,ode,0,modelRequest);
        case 2
            result = hmodel(x,mdl,dat.meteo,ode,0,modelRequest);
        case 3
            result = sacsma(x,mdl,dat.meteo,ode,0,modelRequest);
        case 4
            result = xinanjiang(x,mdl,dat.meteo,ode,0,modelRequest);
        case 5
            result = gr4jA(x,mdl,dat.meteo,ode,0,modelRequest);
        case 6
            result = hbv(x,mdl,dat.meteo,ode,0,modelRequest);
        case 7
            result = cfe_nwm(x,mdl,dat.meteo,ode,0,modelRequest);
        case 11
            isOde = isfield(mdl,'variant') ...
                && strcmpi(string(mdl.variant),'gchm_ode');
            if isOde
                modelFcn = @gchm_ode;
            else
                modelFcn = @gchm;
            end
            if exist(func2str(modelFcn),'file') ~= 2
                error('crr_model:GCHMUnavailable', ...
                    'Private gchm is not on the MATLAB path.');
            end
            result = modelFcn(x,mdl,dat.meteo,ode,0,modelRequest);
        case 99
            result = user_model(x,mdl,dat.meteo,ode,0,modelRequest);
        otherwise
            error('crr_model:UnknownModel', ...
                'Unknown model number %d.',mdl.model);
    end

    loss_value = NaN;
    out = result;
    out.gradient = zeros(d,1);
    out.metrics = struct();
    if isfield(result,'obs') && isfield(result.obs,'Q')
        out.q = result.obs.Q;
    end
end

function [swe,Jswe] = local_swe_from_history(mdl,Z,d,nq,Jth)
%LOCAL_SWE_FROM_HISTORY Extract SWE and its parameter sensitivities.

    if isempty(Z)
        error('crr_model:SWEUnavailable', ...
            ['The selected model did not ' ...
            'return augmented state history.']);
    end
    if mdl.model == 6
        m = numel(mdl.y0);
    else
        m = size(Z,2)/(d+1);
    end
    dz = size(Z,2)/m - 1;
    if m ~= fix(m) || dz ~= fix(dz)
        error('crr_model:InvalidAugmentedHistory', ...
            ['Cannot infer model-state ' ...
            'count from augmented history.']);
    end
    if isfield(mdl,'swe_ind') ...
            && ~isempty(mdl.swe_ind)
        stateIds = double(mdl.swe_ind(:)');
    elseif isfield(mdl,'state_name') ...
            && ~isempty(mdl.state_name)
        stateIds = find(strcmpi(string(mdl.state_name(:)), ...
            "snow_water_equivalent"))';
    elseif mdl.model ~= 99
        stateIds = 1;
    else
        error('crr_model:UserModelSWEUndefined', ...
            ['Define mdl.swe_ind or include a state named ' ...
             '"snow_water_equivalent" in mdl.state_name.']);
    end
    if isempty(stateIds) ...
            || any(stateIds ~= fix(stateIds)) ...
            || any(stateIds < 1) ...
            || any(stateIds >= m)
        error('crr_model:InvalidSWEStates', ...
            ['SWE state indices must ' ...
            'identify physical model states.']);
    end
    rows = mdl.idx(1)+1:mdl.idx(2);
    swe = sum(Z(rows,stateIds),2);
    if numel(swe) ~= nq
        error('crr_model:SWELengthMismatch', ...
            ['SWE and discharge histories ' ...
            'have different lengths.']);
    end
    if nargout < 2
        return
    end
    if nargin < 5 ...
            || numel(Jth) ~= d
        error('crr_model:SWEJacobianUnavailable', ...
            ['The model must return its ' ...
            'd-parameter transformation derivative.']);
    end
    Jswe = zeros(nq,d);
    for j = 1:min(d,dz)
        sensitivityCols = j*m + stateIds;
        Jswe(:,j) = sum(Z(rows,sensitivityCols),2) .* Jth(j);
    end
end

function [sm,Jsm] = local_sm_from_history(mdl,Z,d,nq,Jth)
%LOCAL_SM_FROM_HISTORY Extract soil-water storage and sensitivities.

    if isfield(mdl,'sm_ind') && ~isempty(mdl.sm_ind)
        stateIds = double(mdl.sm_ind(:)');
    else
        switch mdl.model
            case 2
                stateIds = 3;
            case 3
                stateIds = 2:6;
            case 4
                stateIds = 2:3;
            case {5,6,7}
                stateIds = 2;
            otherwise
                stateIds = [];
        end
    end
    if isempty(stateIds)
        error('crr_model:SMUnavailable', ...
            'SM is not defined for model %d.',mdl.model);
    end
    if isempty(Z)
        error('crr_model:SMUnavailable', ...
            'The model did not return augmented state history.');
    end
    if mdl.model == 6
        m = numel(mdl.y0);
    else
        m = size(Z,2)/(d+1);
    end
    dz = size(Z,2)/m - 1;
    if m ~= fix(m) || dz ~= fix(dz) ...
            || any(stateIds < 1) || any(stateIds >= m)
        error('crr_model:InvalidSMStates', ...
            'Invalid SM state mapping for model %d.',mdl.model);
    end
    rows = mdl.idx(1)+1:mdl.idx(2);
    sm = sum(Z(rows,stateIds),2);
    if numel(sm) ~= nq
        error('crr_model:SMLengthMismatch', ...
            'SM and discharge histories have different lengths.');
    end
    if nargout < 2, return, end
    Jsm = zeros(nq,d);
    for j = 1:min(d,dz)
        Jsm(:,j) = sum(Z(rows,j*m+stateIds),2).*Jth(j);
    end
end
