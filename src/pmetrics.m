function prf = pmetrics(bas,loss,L,met,i,prf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PMETRICS Update SAGE performance histories.
%
%  Converts basin metrics and losses into scenario values and iteration
%  summaries.
%
% SYNOPSIS:
%   prf = pmetrics(bas,loss,L,met,i,prf)
%
% INPUT ARGUMENTS:
%   bas             training and evaluation basin counts
%    .K_t            number of training basins
%    .K_e            number of evaluation basins
%   loss            selected loss and observed variables
%    .fnc            loss-function identifier
%    .observed       names in a joint objective
%   L               basin-wise objective values
%   met             compact metrics returned by CAMELS
%   i               current SAGE iteration
%   prf             prior performance-history structure
%
% OUTPUT ARGUMENTS:
%   prf             updated performance structure
%    .curr           current basin-wise scenario metrics
%    .iter           scalar metric and loss histories
%
% NOTES:
%   Scenario suffixes tt, te, et, and ee denote training or evaluation
%   basins crossed with training or evaluation periods.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Aug. 2026             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K = bas.K;
    K_t = bas.K_t;
    K_e = bas.K_e;
    loss_fnc = loss.fnc;
    names = local_observation_names(loss);
    base = local_base_metrics(met,names,K);
    
    idx_t = 1:K_t;
    if K_e > 0
        idx_e = K_t+1:K;
    else
        idx_e = [];
    end
    
    % -----------------------------------------
    % Scenario availability implied by mdl.mode
    % -----------------------------------------
    has = prf.has;
    
    % -------------------------------------------------------------
    % Build basin-wise performance for the current iteration.
    % These arrays are overwritten rather than stored by iteration.
    % -------------------------------------------------------------
    NSE = struct('tt',[],'te',[],'et',[],'ee',[]);
    KGE = struct('tt',[],'te',[],'et',[],'ee',[]);
    JKGE = struct('tt',[],'te',[],'et',[],'ee',[]);
    D_fdc = struct('tt',[],'te',[],'et',[],'ee',[]);
    D_p = struct('tt',[],'te',[],'et',[],'ee',[]);
    D_logp = struct('tt',[],'te',[],'et',[],'ee',[]);
    S_fdc = struct('tt',[],'te',[],'et',[],'ee',[]);
    S_p = struct('tt',[],'te',[],'et',[],'ee',[]);
    S_logp = struct('tt',[],'te',[],'et',[],'ee',[]);
    KGE.components = local_component_scenarios( ...
        {'r','alpha','beta'});
    JKGE.components = local_component_scenarios( ...
        {'M','V','C'});
    
    % -----------------------------
    % Always-available scenario: tt
    % -----------------------------
    NSE.tt = base.performance.t.NSE(idx_t).';
    KGE.tt = base.performance.t.KGE(idx_t).';
    JKGE.tt = base.performance.t.JKGE(idx_t).';
    D_fdc.tt = base.performance.t.D_fdc(idx_t).';
    D_p.tt = base.performance.t.D_p(idx_t).';
    D_logp.tt = base.performance.t.D_logp(idx_t).';
    KGE.components.tt = local_kge_components( ...
        base.performance.t,idx_t);
    JKGE.components.tt = local_jkge_components( ...
        base.performance.t,idx_t);
    
    % ------------------------------------
    % Optional scenarios based on mdl.mode
    % ------------------------------------
    if has.te
        NSE.te = base.performance.e.NSE(idx_t).';
        KGE.te = base.performance.e.KGE(idx_t).';
        JKGE.te = base.performance.e.JKGE(idx_t).';
        D_fdc.te = base.performance.e.D_fdc(idx_t).';
        D_p.te = base.performance.e.D_p(idx_t).';
        D_logp.te = base.performance.e.D_logp(idx_t).';
        KGE.components.te = local_kge_components( ...
            base.performance.e,idx_t);
        JKGE.components.te = local_jkge_components( ...
            base.performance.e,idx_t);
    end
    
    if has.et
        if isempty(idx_e)
            NSE.et = [];
            KGE.et = [];
            D_fdc.et = [];
            D_p.et = [];
            D_logp.et = [];
            JKGE.et = [];
        else
            NSE.et = base.performance.t.NSE(idx_e).';
            KGE.et = base.performance.t.KGE(idx_e).';
            JKGE.et = base.performance.t.JKGE(idx_e).';
            D_fdc.et = base.performance.t.D_fdc(idx_e).';
            D_p.et = base.performance.t.D_p(idx_e).';
            D_logp.et = base.performance.t.D_logp(idx_e).';
            KGE.components.et = local_kge_components( ...
                base.performance.t,idx_e);
            JKGE.components.et = local_jkge_components( ...
                base.performance.t,idx_e);
        end
    end
    
    if has.ee
        if isempty(idx_e)
            NSE.ee = [];
            KGE.ee = [];
            D_fdc.ee = [];
            D_p.ee = [];
            D_logp.ee = [];
            JKGE.ee = [];        
        else
            NSE.ee = base.performance.e.NSE(idx_e).';
            KGE.ee = base.performance.e.KGE(idx_e).';
            JKGE.ee = base.performance.e.JKGE(idx_e).';
            D_fdc.ee = base.performance.e.D_fdc(idx_e).';
            D_p.ee = base.performance.e.D_p(idx_e).';
            D_logp.ee = base.performance.e.D_logp(idx_e).';
            KGE.components.ee = local_kge_components( ...
                base.performance.e,idx_e);
            JKGE.components.ee = local_jkge_components( ...
                base.performance.e,idx_e);
        end
    end
    
    % --------------
    % Loss histories
    % --------------
    prf.iter.L.tt(i) = mean(L(idx_t),'omitnan');
    
    if has.te
        prf.iter.L.te(i) = local_block_loss(loss, ...
            base.diagnostic.e.SAR(idx_t), ...
            base.diagnostic.e.GLS(idx_t), ...
            base.performance.e.NSE(idx_t), ...
            base.performance.e.KGE(idx_t), ...
            base.diagnostic.e.Huber(idx_t), ...
            base.performance.e.D_fdc(idx_t), ...
            base.performance.e.D_p(idx_t), ...
            base.performance.e.D_logp(idx_t), ...
            base.performance.e.JKGE(idx_t));
    else
        prf.iter.L.te(i) = NaN;
    end
    
    if has.et ...
            && ~isempty(idx_e)
        prf.iter.L.et(i) = mean(L(idx_e), ...
            'omitnan');
    else
        prf.iter.L.et(i) = NaN;
    end
    
    if has.ee ...
            && ~isempty(idx_e)
        prf.iter.L.ee(i) = local_block_loss(loss, ...
            base.diagnostic.e.SAR(idx_e), ...
            base.diagnostic.e.GLS(idx_e), ...
            base.performance.e.NSE(idx_e), ...
            base.performance.e.KGE(idx_e), ...
            base.diagnostic.e.Huber(idx_e), ...
            base.performance.e.D_fdc(idx_e), ...
            base.performance.e.D_p(idx_e), ...
            base.performance.e.D_logp(idx_e), ...
            base.performance.e.JKGE(idx_e));
    else
        prf.iter.L.ee(i) = NaN;
    end
    
    % ----------------
    % Aggregate totals
    % ----------------
    prf.iter.SAR.tt(i) = mean(base.diagnostic.t.SAR(idx_t), ...
        'omitnan');
    prf.iter.GLS.tt(i) = mean(base.diagnostic.t.GLS(idx_t), ...
        'omitnan');
    prf.iter.RSS.tt(i) = mean(base.diagnostic.t.RSS(idx_t), ...
        'omitnan');
    
    prf.iter.NSE.tt(i) = mean(base.performance.t.NSE(idx_t), ...
        'omitnan');
    prf.iter.KGE.tt(i) = mean(base.performance.t.KGE(idx_t), ...
        'omitnan');
    
    prf.iter.Huber.tt(i) = mean(base.diagnostic.t.Huber(idx_t), ...
        'omitnan');
    prf.iter.JKGE.tt(i) = mean(base.performance.t.JKGE(idx_t), ...
        'omitnan');
    
    if has.te
        prf.iter.SAR.te(i) = mean(base.diagnostic.e.SAR(idx_t), ...
            'omitnan');
        prf.iter.GLS.te(i) = mean(base.diagnostic.e.GLS(idx_t), ...
            'omitnan');
        prf.iter.RSS.te(i) = mean(base.diagnostic.e.RSS(idx_t), ...
            'omitnan');
    
        prf.iter.NSE.te(i) = mean(base.performance.e.NSE(idx_t), ...
            'omitnan');
        prf.iter.KGE.te(i) = mean(base.performance.e.KGE(idx_t), ...
            'omitnan');
        
        prf.iter.Huber.te(i) = mean(base.diagnostic.e.Huber(idx_t), ...
            'omitnan');
        prf.iter.JKGE.te(i) = mean(base.performance.e.JKGE(idx_t), ...
            'omitnan');
    else
        prf.iter.SAR.te(i) = NaN;
        prf.iter.GLS.te(i) = NaN;
        prf.iter.RSS.te(i) = NaN;
        prf.iter.NSE.te(i) = NaN;
        prf.iter.KGE.te(i) = NaN;   
        prf.iter.Huber.te(i) = NaN;
        prf.iter.JKGE.te(i) = NaN;    
    end
    
    if has.et ...
            && ~isempty(idx_e)
        prf.iter.SAR.et(i) = mean(base.diagnostic.t.SAR(idx_e), ...
            'omitnan');
        prf.iter.GLS.et(i) = mean(base.diagnostic.t.GLS(idx_e), ...
            'omitnan');
        prf.iter.RSS.et(i) = mean(base.diagnostic.t.RSS(idx_e), ...
            'omitnan');
    
        prf.iter.NSE.et(i) = mean(base.performance.t.NSE(idx_e), ...
            'omitnan');
        prf.iter.KGE.et(i) = mean(base.performance.t.KGE(idx_e), ...
            'omitnan');
        
        prf.iter.Huber.et(i) = mean(base.diagnostic.t.Huber(idx_e), ...
            'omitnan');
        prf.iter.JKGE.et(i) = mean(base.performance.t.JKGE(idx_e), ...
            'omitnan');
    else
        prf.iter.SAR.et(i) = NaN;
        prf.iter.GLS.et(i) = NaN;
        prf.iter.RSS.et(i) = NaN;
    
        prf.iter.NSE.et(i) = NaN;
        prf.iter.KGE.et(i) = NaN;
        
        prf.iter.Huber.et(i) = NaN;
        prf.iter.JKGE.et(i) = NaN;
    end
    
    if has.ee ...
            && ~isempty(idx_e)
        prf.iter.SAR.ee(i) = mean(base.diagnostic.e.SAR(idx_e), ...
            'omitnan');
        prf.iter.GLS.ee(i) = mean(base.diagnostic.e.GLS(idx_e), ...
            'omitnan');
        prf.iter.RSS.ee(i) = mean(base.diagnostic.e.RSS(idx_e), ...
            'omitnan');
    
        prf.iter.NSE.ee(i) = mean(base.performance.e.NSE(idx_e), ...
            'omitnan');   
        prf.iter.KGE.ee(i) = mean(base.performance.e.KGE(idx_e), ...
            'omitnan');   
    
        prf.iter.Huber.ee(i) = mean(base.diagnostic.e.Huber(idx_e), ...
            'omitnan');
        prf.iter.JKGE.ee(i) = mean(base.performance.e.JKGE(idx_e), ...
            'omitnan');   
    else
        prf.iter.SAR.ee(i) = NaN;
        prf.iter.GLS.ee(i) = NaN;
        prf.iter.RSS.ee(i) = NaN;
    
        prf.iter.NSE.ee(i) = NaN;   
        prf.iter.KGE.ee(i) = NaN;   
    
        prf.iter.Huber.ee(i) = NaN;
        prf.iter.JKGE.ee(i) = NaN;   
    end
    
    % ----------------------------------------------------------------------
    % Basin- and period-specific dimensionless FDC skill scores
    % ----------------------------------------------------------------------
    if local_has_observation(loss,"Q")
    if ~isfield(loss,'fdc') ...
            || ~isstruct(loss.fdc) ...
            || ~isfield(loss.fdc,'Q') ...
            || ~isstruct(loss.fdc.Q) ...
            || ~isfield(loss.fdc.Q,'D0t') ...
            || ~isfield(loss.fdc.Q,'D0e') ...
            || ~isfield(loss.fdc.Q,'D0pt') ...
            || ~isfield(loss.fdc.Q,'D0pe') ...
            || ~isfield(loss.fdc.Q,'D0logpt') ...
            || ~isfield(loss.fdc.Q,'D0logpe') ...
            || numel(loss.fdc.Q.D0t) ~= K ...
            || numel(loss.fdc.Q.D0e) ~= K ...
            || numel(loss.fdc.Q.D0pt) ~= K ...
            || numel(loss.fdc.Q.D0pe) ~= K ...
            || numel(loss.fdc.Q.D0logpt) ~= K ...
            || numel(loss.fdc.Q.D0logpe) ~= K
    
        error(['      Error:pmetrics: ' ...
            'one or more FDC reference vectors are missing ' ...
            'or has incorrect size.']);
    end
    
    D0t = double(loss.fdc.Q.D0t(:));
    D0e = double(loss.fdc.Q.D0e(:));
    D0pt = double(loss.fdc.Q.D0pt(:));
    D0pe = double(loss.fdc.Q.D0pe(:));
    D0logpt = double(loss.fdc.Q.D0logpt(:));
    D0logpe = double(loss.fdc.Q.D0logpe(:));
    D0t_t = D0t(idx_t);
    D0e_t = D0e(idx_t);
    D0pt_t = D0pt(idx_t);
    D0pe_t = D0pe(idx_t);
    D0logpt_t = D0logpt(idx_t);
    D0logpe_t = D0logpe(idx_t);
    if ~isempty(idx_e)
        D0t_e = D0t(idx_e);
        D0e_e = D0e(idx_e);
        D0pt_e = D0pt(idx_e);
        D0pe_e = D0pe(idx_e);
        D0logpt_e = D0logpt(idx_e);
        D0logpe_e = D0logpe(idx_e);
    else
        D0t_e = [];
        D0e_e = [];
        D0pt_e = [];
        D0pe_e = [];
        D0logpt_e = [];
        D0logpe_e = [];
    end
    
    S_fdc.tt = local_fdc_score(D_fdc.tt,D0t_t);
    S_p.tt = local_fdc_score(D_p.tt,D0pt_t);
    S_logp.tt = local_fdc_score(D_logp.tt,D0logpt_t);
    if has.te
        S_fdc.te = local_fdc_score(D_fdc.te,D0e_t);
        S_p.te = local_fdc_score(D_p.te,D0pe_t);
        S_logp.te = local_fdc_score(D_logp.te,D0logpe_t);
    end
    if has.et ...
            && ~isempty(idx_e)
        S_fdc.et = local_fdc_score(D_fdc.et,D0t_e);
        S_p.et = local_fdc_score(D_p.et,D0pt_e);
        S_logp.et = local_fdc_score(D_logp.et,D0logpt_e);
    end
    if has.ee ...
            && ~isempty(idx_e)
        S_fdc.ee = local_fdc_score(D_fdc.ee,D0e_e);
        S_p.ee = local_fdc_score(D_p.ee,D0pe_e);
        S_logp.ee = local_fdc_score(D_logp.ee,D0logpe_e);
    end
    end
    
    % -----------------
    % Summary histories
    % -----------------
    % Support runs initialized before the S_ib,JKGE history was introduced.
    % Preallocate with NaN so MATLAB does not backfill earlier iterations
    % with misleading zeros when this field is first assigned mid-run.
    if ~isfield(prf.iter,'Sib_JKGE')
        prf.iter.Sib_JKGE = struct();
        sibScenarios = {'tt','te','et','ee'};
        for iScenario = 1:numel(sibScenarios)
            scn = sibScenarios{iScenario};
            prf.iter.Sib_JKGE.(scn) = ...
                nan(size(prf.iter.mJKGE.(scn)));
        end
    end

    prf.iter.mNSE.tt(i) = local_median(NSE.tt);
    prf.iter.mKGE.tt(i) = local_median(KGE.tt);
    prf.iter.mJKGE.tt(i) = local_median(JKGE.tt);
    prf.iter.Sib_NSE.tt(i) = local_mean_one_minus(NSE.tt);
    prf.iter.Sib_KGE.tt(i) = local_mean_one_minus(KGE.tt);
    if loss_fnc == 7
        prf.iter.Sib_JKGE.tt(i) = local_mean_one_minus(JKGE.tt);
    else
        prf.iter.Sib_JKGE.tt(i) = NaN;
    end
    prf.iter.Sib_S_fdc.tt(i) = local_mean_one_minus(S_fdc.tt);
    prf.iter.S_fdc.tt(i) = local_mean(S_fdc.tt);
    prf.iter.mS_fdc.tt(i) = local_median(S_fdc.tt);
    prf.iter.Sib_S_p.tt(i) = local_mean_one_minus(S_p.tt);
    prf.iter.S_p.tt(i) = local_mean(S_p.tt);
    prf.iter.mS_p.tt(i) = local_median(S_p.tt);
    prf.iter.Sib_S_logp.tt(i) = local_mean_one_minus(S_logp.tt);
    prf.iter.S_logp.tt(i) = local_mean(S_logp.tt);
    prf.iter.mS_logp.tt(i) = local_median(S_logp.tt);
    
    prf.iter.mNSE.te(i) = local_median(NSE.te);
    prf.iter.mKGE.te(i) = local_median(KGE.te);
    prf.iter.mJKGE.te(i) = local_median(JKGE.te);
    prf.iter.Sib_NSE.te(i) = local_mean_one_minus(NSE.te);
    prf.iter.Sib_KGE.te(i) = local_mean_one_minus(KGE.te);
    if loss_fnc == 7
        prf.iter.Sib_JKGE.te(i) = local_mean_one_minus(JKGE.te);
    else
        prf.iter.Sib_JKGE.te(i) = NaN;
    end
    prf.iter.Sib_S_fdc.te(i) = local_mean_one_minus(S_fdc.te);
    prf.iter.S_fdc.te(i) = local_mean(S_fdc.te);
    prf.iter.mS_fdc.te(i) = local_median(S_fdc.te);
    prf.iter.Sib_S_p.te(i) = local_mean_one_minus(S_p.te);
    prf.iter.S_p.te(i) = local_mean(S_p.te);
    prf.iter.mS_p.te(i) = local_median(S_p.te);
    prf.iter.Sib_S_logp.te(i) = local_mean_one_minus(S_logp.te);
    prf.iter.S_logp.te(i) = local_mean(S_logp.te);
    prf.iter.mS_logp.te(i) = local_median(S_logp.te);
    
    prf.iter.mNSE.et(i) = local_median(NSE.et);
    prf.iter.mKGE.et(i) = local_median(KGE.et);
    prf.iter.mJKGE.et(i) = local_median(JKGE.et);
    prf.iter.Sib_NSE.et(i) = local_mean_one_minus(NSE.et);
    prf.iter.Sib_KGE.et(i) = local_mean_one_minus(KGE.et);
    if loss_fnc == 7
        prf.iter.Sib_JKGE.et(i) = local_mean_one_minus(JKGE.et);
    else
        prf.iter.Sib_JKGE.et(i) = NaN;
    end
    prf.iter.Sib_S_fdc.et(i) = local_mean_one_minus(S_fdc.et);
    prf.iter.S_fdc.et(i) = local_mean(S_fdc.et);
    prf.iter.mS_fdc.et(i) = local_median(S_fdc.et);
    prf.iter.Sib_S_p.et(i) = local_mean_one_minus(S_p.et);
    prf.iter.S_p.et(i) = local_mean(S_p.et);
    prf.iter.mS_p.et(i) = local_median(S_p.et);
    prf.iter.Sib_S_logp.et(i) = local_mean_one_minus(S_logp.et);
    prf.iter.S_logp.et(i) = local_mean(S_logp.et);
    prf.iter.mS_logp.et(i) = local_median(S_logp.et);
    
    prf.iter.mNSE.ee(i) = local_median(NSE.ee);
    prf.iter.mKGE.ee(i) = local_median(KGE.ee);
    prf.iter.mJKGE.ee(i) = local_median(JKGE.ee);
    prf.iter.Sib_NSE.ee(i) = local_mean_one_minus(NSE.ee);
    prf.iter.Sib_KGE.ee(i) = local_mean_one_minus(KGE.ee);
    if loss_fnc == 7
        prf.iter.Sib_JKGE.ee(i) = local_mean_one_minus(JKGE.ee);
    else
        prf.iter.Sib_JKGE.ee(i) = NaN;
    end
    prf.iter.Sib_S_fdc.ee(i) = local_mean_one_minus(S_fdc.ee);
    prf.iter.S_fdc.ee(i) = local_mean(S_fdc.ee);
    prf.iter.mS_fdc.ee(i) = local_median(S_fdc.ee);
    prf.iter.Sib_S_p.ee(i) = local_mean_one_minus(S_p.ee);
    prf.iter.S_p.ee(i) = local_mean(S_p.ee);
    prf.iter.mS_p.ee(i) = local_median(S_p.ee);
    prf.iter.Sib_S_logp.ee(i) = local_mean_one_minus(S_logp.ee);
    prf.iter.S_logp.ee(i) = local_mean(S_logp.ee);
    prf.iter.mS_logp.ee(i) = local_median(S_logp.ee);

    if isfield(met,'total')
        prf = local_joint_histories( ...
            prf,met,loss,i,idx_t,idx_e,has);
    end
    
    % -------------------------------------------------------
    % Regional KGE and JKGE component histories per scenario
    % -------------------------------------------------------
    scenarios = {'tt','te','et','ee'};
    for j = 1:numel(scenarios)
        sc = scenarios{j};
        kc = KGE.components.(sc);
        jc = JKGE.components.(sc);
    
        prf.iter.KGE_r.(sc)(i) = local_mean(kc.r);
        prf.iter.mKGE_r.(sc)(i) = local_median(kc.r);
        prf.iter.KGE_alpha.(sc)(i) = local_mean(kc.alpha);
        prf.iter.mKGE_alpha.(sc)(i) = local_median(kc.alpha);
        prf.iter.KGE_beta.(sc)(i) = local_mean(kc.beta);
        prf.iter.mKGE_beta.(sc)(i) = local_median(kc.beta);
    
        prf.iter.JKGE_M.(sc)(i) = local_mean(jc.M);
        prf.iter.mJKGE_M.(sc)(i) = local_median(jc.M);
        prf.iter.JKGE_V.(sc)(i) = local_mean(jc.V);
        prf.iter.mJKGE_V.(sc)(i) = local_median(jc.V);
        prf.iter.JKGE_C.(sc)(i) = local_mean(jc.C);
        prf.iter.mJKGE_C.(sc)(i) = local_median(jc.C);
    end

    % Basin-wise performance for the most recently evaluated parameters.
    % Unlike prf.iter, these arrays are replaced at every iteration.
    prf.curr.NSE = NSE;
    prf.curr.KGE = KGE;
    prf.curr.S_fdc = S_fdc;
    prf.curr.S_p = S_p;
    prf.curr.S_logp = S_logp;
    prf.curr.JKGE = JKGE;
    prf.curr.loss_fnc = loss_fnc;
    prf.curr.fdc_formulation = local_fdc_formulation(loss);
    % Named-observation current values are populated by
    % local_joint_histories for exactly the selected variables.

end

function names = local_observation_names(loss)
%LOCAL_OBSERVATION_NAMES Return normalized requested observation names.

    names = "Q";
    if isfield(loss,'observed') && ~isempty(loss.observed)
        names = unique(upper(strtrim( ...
            string(loss.observed(:)))),'stable');
        names = names(strlength(names) > 0);
    end
end

function base = local_base_metrics(met,names,K)
%LOCAL_BASE_METRICS Select data for the original single-variable reports.

    if isfield(met.variable,'Q')
        base = met.variable.Q;
    else
        base = met.variable.(char(names(1)));
    end
    if ~isfield(base,'diagnostic')
        z = nan(1,K);
        block = struct('SAR',z,'GLS',z,'Huber',z,'RSS',z);
        base.diagnostic = struct('t',block,'e',block);
    end
    periods = {'t','e'};
    for j = 1:2
        per = periods{j};
        P = base.performance.(per);
        if ~isfield(P,'KGE_components')
            z = nan(1,K);
            P.KGE_components = struct( ...
                'r',z,'alpha',z,'beta',z);
        end
        base.performance.(per) = P;
    end
end

function tf = local_has_observation(loss,name)
%LOCAL_HAS_OBSERVATION Test whether a variable is a training objective.

    names = "Q";
    if isfield(loss,'observed') ...
            && ~isempty(loss.observed)
        names = unique(upper(strtrim( ...
            string(loss.observed(:)))),'stable');
        names = names(strlength(names) > 0);
    end
    tf = any(names == upper(string(name)));
end

function prf = local_joint_histories( ...
    prf,met,loss,i,idx_t,idx_e,has)
%LOCAL_JOINT_HISTORIES Retain joint losses and per-observation diagnostics.

    names = local_observation_names(loss);
    scenarios = {'tt','te','et','ee'};
    periods = {'t','e','t','e'};
    indices = {idx_t,idx_t,idx_e,idx_e};
    enabled = [true,has.te,has.et,has.ee];
    formulation = local_fdc_formulation(loss);

    dFields = {'D_fdc','D_p','D_logp'};
    d0Train = {'D0t','D0pt','D0logpt'};
    d0Eval = {'D0e','D0pe','D0logpe'};

    for j = 1:4
        sc = scenarios{j};
        per = periods{j};
        id = indices{j};
        if ~enabled(j) ...
                || isempty(id)
            prf = local_store_empty_joint(prf,sc,i,names);
            continue
        end

        prf.iter.joint.loss.total.(sc)(i) = ...
            local_sum(met.total.loss.(per)(id));
        for k = 1:numel(names)
            name = char(names(k));
            prf.iter.joint.loss.(name).(sc)(i) = ...
                local_sum(met.variable.(name).loss.(per)(id));
            prf.iter.joint.contribution.(name).(sc)(i) = ...
                local_sum(met.total.contribution. ...
                (per).(name)(id));

            D = met.variable.(name).diagnostic.(per);
            diagnosticNames = {'SAR','GLS','RSS','Huber'};
            for d = 1:numel(diagnosticNames)
                diagnosticName = diagnosticNames{d};
                values = D.(diagnosticName)(id);
                prf.curr.joint.(name).(diagnosticName).(sc) = ...
                    values.';
                prf.iter.joint.(name).mean. ...
                    (diagnosticName).(sc)(i) = local_mean(values);
                prf.iter.joint.(name).median. ...
                    (diagnosticName).(sc)(i) = local_median(values);
            end

            S = met.variable.(name).performance.(per);
            prf.curr.joint.(name).NSE.(sc) = S.NSE(id).';
            prf.curr.joint.(name).KGE.(sc) = S.KGE(id).';
            prf.iter.joint.(name).mean.NSE.(sc)(i) = ...
                local_mean(S.NSE(id));
            prf.iter.joint.(name).mean.KGE.(sc)(i) = ...
                local_mean(S.KGE(id));
            prf.iter.joint.(name).median.NSE.(sc)(i) = ...
                local_median(S.NSE(id));
            prf.iter.joint.(name).median.KGE.(sc)(i) = ...
                local_median(S.KGE(id));
            prf.iter.joint.(name).integrated.NSE.(sc)(i) = ...
                local_mean_one_minus(S.NSE(id));
            prf.iter.joint.(name).integrated.KGE.(sc)(i) = ...
                local_mean_one_minus(S.KGE(id));
            components = {'r','alpha','beta'};
            for c = 1:numel(components)
                component = components{c};
                values = S.KGE_components.(component)(id);
                prf.curr.joint.(name).KGE_components. ...
                    (component).(sc) = values.';
                prf.iter.joint.(name).KGE_components. ...
                    (component).(sc)(i) = local_median(values);
                prf.iter.joint.(name).mean.KGE_components. ...
                    (component).(sc)(i) = local_mean(values);
                prf.iter.joint.(name).median.KGE_components. ...
                    (component).(sc)(i) = local_median(values);
            end
            prf.curr.joint.(name).JKGE.(sc) = S.JKGE(id).';
            prf.iter.joint.(name).NSE.(sc)(i) = ...
                local_median(S.NSE(id));
            prf.iter.joint.(name).KGE.(sc)(i) = ...
                local_median(S.KGE(id));
            prf.iter.joint.(name).JKGE.(sc)(i) = ...
                local_median(S.JKGE(id));
            prf.iter.joint.(name).mean.JKGE.(sc)(i) = ...
                local_mean(S.JKGE(id));
            prf.iter.joint.(name).median.JKGE.(sc)(i) = ...
                local_median(S.JKGE(id));
            prf.iter.joint.(name).integrated.JKGE.(sc)(i) = ...
                local_mean_one_minus(S.JKGE(id));
            components = {'M','V','C'};
            for c = 1:numel(components)
                component = components{c};
                values = S.JKGE_components.(component)(id);
                prf.iter.joint.(name).JKGE_components. ...
                    (component).(sc)(i) = local_median(values);
                prf.iter.joint.(name).mean.JKGE_components. ...
                    (component).(sc)(i) = local_mean(values);
                prf.iter.joint.(name).median.JKGE_components. ...
                    (component).(sc)(i) = local_median(values);
            end

            durationNames = {'S_fdc','S_p','S_logp'};
            for d = 1:numel(durationNames)
                D = S.(dFields{d})(id).';
                score = nan(size(D));
                if isfield(loss,'fdc') ...
                        && isfield(loss.fdc,name)
                    refs = loss.fdc.(name);
                    if per == 't'
                        refField = d0Train{d};
                    else
                        refField = d0Eval{d};
                    end
                    if isfield(refs,refField) ...
                            && numel(refs.(refField)) >= max(id)
                        D0 = double(refs.(refField)(id));
                        score = local_fdc_score(D,D0);
                    end
                end
                durationName = durationNames{d};
                prf.curr.joint.(name).(durationName).(sc) = score;
                prf.iter.joint.(name).duration_metrics. ...
                    (durationName).mean.(sc)(i) = local_mean(score);
                prf.iter.joint.(name).duration_metrics. ...
                    (durationName).median.(sc)(i) = local_median(score);
                prf.iter.joint.(name).duration_metrics. ...
                    (durationName).integrated.(sc)(i) = ...
                    local_mean_one_minus(score);
                if d == formulation
                    prf.iter.joint.(name).duration.(sc)(i) = ...
                        local_median(score);
                end
            end
        end
    end
end

function prf = local_store_empty_joint(prf,sc,i,names)
%LOCAL_STORE_EMPTY_JOINT Store unavailable joint scenario values.

    prf.iter.joint.loss.total.(sc)(i) = NaN;
    for k = 1:numel(names)
        name = char(names(k));
        prf.iter.joint.loss.(name).(sc)(i) = NaN;
        prf.iter.joint.contribution.(name).(sc)(i) = NaN;
        prf.iter.joint.(name).NSE.(sc)(i) = NaN;
        prf.iter.joint.(name).KGE.(sc)(i) = NaN;
        diagnosticNames = {'SAR','GLS','RSS','Huber'};
        for d = 1:numel(diagnosticNames)
            diagnosticName = diagnosticNames{d};
            prf.iter.joint.(name).mean. ...
                (diagnosticName).(sc)(i) = NaN;
            prf.iter.joint.(name).median. ...
                (diagnosticName).(sc)(i) = NaN;
            prf.curr.joint.(name).(diagnosticName).(sc) = [];
        end
        summaries = {'mean','median','integrated'};
        metrics = {'NSE','KGE','JKGE'};
        for s = 1:numel(summaries)
            summary = summaries{s};
            for m = 1:numel(metrics)
                metric = metrics{m};
                prf.iter.joint.(name).(summary). ...
                    (metric).(sc)(i) = NaN;
            end
        end
        components = {'r','alpha','beta'};
        for c = 1:numel(components)
            component = components{c};
            prf.iter.joint.(name).KGE_components. ...
                (component).(sc)(i) = NaN;
            prf.iter.joint.(name).mean.KGE_components. ...
                (component).(sc)(i) = NaN;
            prf.iter.joint.(name).median.KGE_components. ...
                (component).(sc)(i) = NaN;
            prf.curr.joint.(name).KGE_components. ...
                (component).(sc) = [];
        end
        prf.iter.joint.(name).JKGE.(sc)(i) = NaN;
        components = {'M','V','C'};
        for c = 1:numel(components)
            component = components{c};
            prf.iter.joint.(name).JKGE_components. ...
                (component).(sc)(i) = NaN;
            prf.iter.joint.(name).mean.JKGE_components. ...
                (component).(sc)(i) = NaN;
            prf.iter.joint.(name).median.JKGE_components. ...
                (component).(sc)(i) = NaN;
        end
        prf.iter.joint.(name).duration.(sc)(i) = NaN;
        durationNames = {'S_fdc','S_p','S_logp'};
        durationSummaries = {'mean','median','integrated'};
        for d = 1:numel(durationNames)
            durationName = durationNames{d};
            for s = 1:numel(durationSummaries)
                summary = durationSummaries{s};
                prf.iter.joint.(name).duration_metrics. ...
                    (durationName).(summary).(sc)(i) = NaN;
            end
        end
        prf.curr.joint.(name).NSE.(sc) = [];
        prf.curr.joint.(name).KGE.(sc) = [];
        prf.curr.joint.(name).JKGE.(sc) = [];
    end
end

function S = local_component_scenarios(names)
%LOCAL_COMPONENT_SCENARIOS Initialize empty component structures.

    empty = cell2struct(repmat({[]},1,numel(names)), ...
        names,2);
    S = struct('tt',empty,'te',empty, ...
        'et',empty,'ee',empty);
end

function C = local_kge_components(P,idx)
%LOCAL_KGE_COMPONENTS Select basin-wise KGE components.

    C = struct( ...
        'r',P.KGE_components.r(idx).', ...
        'alpha',P.KGE_components.alpha(idx).', ...
        'beta',P.KGE_components.beta(idx).');
end

function C = local_jkge_components(P,idx)
%LOCAL_JKGE_COMPONENTS Select basin-wise JKGE components.

    C = struct( ...
        'M',P.JKGE_components.M(idx).', ...
        'V',P.JKGE_components.V(idx).', ...
        'C',P.JKGE_components.C(idx).');
end

function Lblk = local_block_loss(loss, ...
    SAR,GLS,NSE,KGE,Huber,D_fdc,D_p,D_logp,JKGE)
%LOCAL_BLOCK_LOSS Basin-aggregated loss over one scenario block.

    loss_fnc = loss.fnc;

    if isempty(SAR) ...
            && isempty(GLS) ...
            && isempty(NSE) ...
            && isempty(KGE) ...
            && isempty(Huber) ...
            && isempty(D_fdc) ...
            && isempty(D_p) ...
            && isempty(D_logp) ...
            && isempty(JKGE)
        Lblk = NaN;
        return
    end
    
    switch loss_fnc
        case 1
            Lblk = mean(SAR,'omitnan');
        case 2
            Lblk = mean(GLS,'omitnan');
        case 3
            Lblk = mean(1 - NSE,'omitnan');
        case 4
            Lblk = mean(1 - KGE,'omitnan');
        case 5
            Lblk = mean(Huber,'omitnan');
        case 6
            switch local_fdc_formulation(loss)
                case 1
                    Lblk = mean(D_fdc,'omitnan');
                case 2
                    Lblk = mean(D_p,'omitnan');
                case 3
                    Lblk = mean(D_logp,'omitnan');
            end
        case 7
            Lblk = mean(1 - JKGE,'omitnan');
        otherwise
            error(['      Error:pmetrics: ' ...
                'Unknown loss choice loss = %g.'], ...
                loss_fnc);
    end
end

function formulation = local_fdc_formulation(loss)
%LOCAL_FDC_FORMULATION Return the selected 6a/6b/6c formulation.
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
        error('pmetrics:BadFDCFormulation', ...
            ['loss.fdc.formulation must be 1 (d_fdc), 2 (d_p), ' ...
             'or 3 (d_logp).']);
    end
end

function m = local_median(x)
%LOCAL_MEDIAN Median with empty-vector protection.
    if isempty(x)
        m = NaN;
    else
        m = median(x,'omitnan');
    end
end

function s = local_mean_one_minus(x)
%LOCAL_MEAN_ONE_MINUS Mean(1-x) with empty-vector protection.
    if isempty(x)
        s = NaN;
    else
        s = mean(1 - x,'omitnan');
    end
end

function m = local_mean(x)
%LOCAL_MEAN Mean with empty-vector protection.
    if isempty(x)
        m = NaN;
    else
        m = mean(x,'omitnan');
    end
end

function s = local_sum(x)
%LOCAL_SUM Sum finite values with empty/all-missing protection.
    if isempty(x) ...
            || ~any(isfinite(x))
        s = NaN;
    else
        s = sum(x,'omitnan');
    end
end

function S = local_fdc_score(D,D0)
%LOCAL_FDC_SCORE Convert raw FDC divergence to dimensionless skill score.
%
%   S_fdc = 1 - D_fdc/D0
%
%   where D0 is the FDC divergence of the constant-median discharge
%   benchmark derived from observations for the period being scored.
%
%   D_fdc = 0  -> S_fdc = 1
%   D_fdc = D0 -> S_fdc = 0
%   D_fdc > D0 -> S_fdc < 0
%
%   Range: (-Inf,1]

    D  = double(D(:));
    D0 = double(D0(:));
    
    S = nan(size(D));
    
    if isempty(D) ...
            || isempty(D0)
        return
    end
    
    if numel(D) ~= numel(D0)
        error(['      Error:pmetrics: ' ...
            'D_fdc and D0 dimensions do not match.']);
    end
    
    ok = isfinite(D) ...
        & isfinite(D0) ...
        & D >= 0 ...
        & D0 > 0;
    
    S(ok) = 1 - D(ok) ./ D0(ok);

end
