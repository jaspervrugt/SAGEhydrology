function prf = pmetrics(bas,loss,L,met,i,prf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PMETRICS Calculates basin-aggregated performance metrics during SAGE
% descent and updates the performance-history structure using the compact
% basin-wise metrics returned by CAMELS/CRR_MODEL.
%
% This function translates the raw basin-wise metric arrays in 'met' into
% current basin-wise values in prf.curr and scalar iteration histories in
% prf.iter.
%
% The notation used throughout is:
%   tt = training basins | training period
%   te = training basins | evaluation period/mask
%   et = evaluation basins | training period
%   ee = evaluation basins | evaluation period/mask
%
% SYNOPSIS:
%  prf = pmetrics(bas,loss,L,met,i,prf)
%
%   bas         structure with basin information
%    .K          total number of watersheds (= K_t + K_e)
%    .K_t        number of training watersheds
%    .K_e        number of evaluation watersheds
%   loss        structure of loss function
%    .fnc        scalar loss function
%                 1 = sum of absolute residuals
%                 2 = generalized least squares / residual sum of squares
%                 3 = 1 - Nash-Sutcliffe efficiency
%                 4 = 1 - Kling-Gupta efficiency
%                 5 = Huber loss
%                 6 = flow-duration-curve loss
%                 7 = 1 - Jawad Kling-Gupta efficiency
%    .n_win      for fnc = 7, moving-average window length in days
%   L           1xK vector of basin-wise loss values returned by CAMELS
%   met         basin-wise output returned by CAMELS for all K basins:
%                .loss.t/.loss.e
%                   SAR, GLS, Huber and RSS values
%                .performance.t/.performance.e
%                   NSE, KGE, D_fdc and JKGE values, plus the nested KGE
%                   components r/alpha/beta and JKGE components M/V/C
%   i           descent iteration counter
%   prf         structure with performance histories
%
% OUTPUT:
%   prf         performance structure with two data lifetimes:
%                .curr.NSE/KGE/S_fdc/JKGE
%                   basin-wise values for the current iteration, organized
%                   as .tt, .te, .et and .ee scenarios
%                .iter
%                   scalar histories stored across SAGE iterations
%
% NOTES:
%   1. The four explicit scenarios tt/te/et/ee are stored in 'prf'.
%   2. No active-comparison aliases with suffix 'a' are used.
%   3. Unavailable scenarios are returned as [] and stored as NaN in 'prf'.
%   4. The integrated basin loss score is stored for NSE, KGE and S_fdc
%      in prf.iter.Sib_NSE, prf.iter.Sib_KGE and prf.iter.Sib_S_fdc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Aug. 2026             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K = bas.K;
    K_t = bas.K_t;
    K_e = bas.K_e;
    loss_fnc = loss.fnc;
    
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
    
    % ---------------------------------------------------------
    % Build basin-wise performance for the current iteration.
    % These arrays are overwritten rather than stored by iteration.
    % ---------------------------------------------------------
    NSE = struct('tt',[],'te',[],'et',[],'ee',[]);
    KGE = struct('tt',[],'te',[],'et',[],'ee',[]);
    JKGE = struct('tt',[],'te',[],'et',[],'ee',[]);
    D_fdc = struct('tt',[],'te',[],'et',[],'ee',[]);
    S_fdc = struct('tt',[],'te',[],'et',[],'ee',[]);
    KGE.components = local_component_scenarios( ...
        {'r','alpha','beta'});
    JKGE.components = local_component_scenarios( ...
        {'M','V','C'});
    
    % -----------------------------
    % Always-available scenario: tt
    % -----------------------------
    NSE.tt = met.performance.t.NSE(idx_t).';
    KGE.tt = met.performance.t.KGE(idx_t).';
    JKGE.tt = met.performance.t.JKGE(idx_t).';
    D_fdc.tt = met.performance.t.D_fdc(idx_t).';
    KGE.components.tt = local_kge_components( ...
        met.performance.t,idx_t);
    JKGE.components.tt = local_jkge_components( ...
        met.performance.t,idx_t);
    
    % ------------------------------------
    % Optional scenarios based on mdl.mode
    % ------------------------------------
    if has.te
        NSE.te = met.performance.e.NSE(idx_t).';
        KGE.te = met.performance.e.KGE(idx_t).';
        JKGE.te = met.performance.e.JKGE(idx_t).';
        D_fdc.te = met.performance.e.D_fdc(idx_t).';
        KGE.components.te = local_kge_components( ...
            met.performance.e,idx_t);
        JKGE.components.te = local_jkge_components( ...
            met.performance.e,idx_t);
    end
    
    if has.et
        if isempty(idx_e)
            NSE.et = [];
            KGE.et = [];
            D_fdc.et = [];
            JKGE.et = [];
        else
            NSE.et = met.performance.t.NSE(idx_e).';
            KGE.et = met.performance.t.KGE(idx_e).';
            JKGE.et = met.performance.t.JKGE(idx_e).';
            D_fdc.et = met.performance.t.D_fdc(idx_e).';
            KGE.components.et = local_kge_components( ...
                met.performance.t,idx_e);
            JKGE.components.et = local_jkge_components( ...
                met.performance.t,idx_e);
        end
    end
    
    if has.ee
        if isempty(idx_e)
            NSE.ee = [];
            KGE.ee = [];
            D_fdc.ee = [];
            JKGE.ee = [];        
        else
            NSE.ee = met.performance.e.NSE(idx_e).';
            KGE.ee = met.performance.e.KGE(idx_e).';
            JKGE.ee = met.performance.e.JKGE(idx_e).';
            D_fdc.ee = met.performance.e.D_fdc(idx_e).';
            KGE.components.ee = local_kge_components( ...
                met.performance.e,idx_e);
            JKGE.components.ee = local_jkge_components( ...
                met.performance.e,idx_e);
        end
    end
    
    % --------------
    % Loss histories
    % --------------
    prf.iter.L.tt(i) = mean(L(idx_t),'omitnan');
    
    if has.te
        prf.iter.L.te(i) = local_block_loss(loss_fnc, ...
            met.loss.e.SAR(idx_t),met.loss.e.GLS(idx_t), ...
            met.performance.e.NSE(idx_t),met.performance.e.KGE(idx_t), ...
            met.loss.e.Huber(idx_t),met.performance.e.D_fdc(idx_t), ...
            met.performance.e.JKGE(idx_t));
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
        prf.iter.L.ee(i) = local_block_loss(loss_fnc, ...
            met.loss.e.SAR(idx_e),met.loss.e.GLS(idx_e), ...
            met.performance.e.NSE(idx_e),met.performance.e.KGE(idx_e), ...
            met.loss.e.Huber(idx_e),met.performance.e.D_fdc(idx_e), ...
            met.performance.e.JKGE(idx_e));
    else
        prf.iter.L.ee(i) = NaN;
    end
    
    % ----------------
    % Aggregate totals
    % ----------------
    prf.iter.SAR.tt(i) = mean(met.loss.t.SAR(idx_t), ...
        'omitnan');
    prf.iter.GLS.tt(i) = mean(met.loss.t.GLS(idx_t), ...
        'omitnan');
    prf.iter.RSS.tt(i) = mean(met.loss.t.RSS(idx_t), ...
        'omitnan');
    
    prf.iter.NSE.tt(i) = mean(met.performance.t.NSE(idx_t), ...
        'omitnan');
    prf.iter.KGE.tt(i) = mean(met.performance.t.KGE(idx_t), ...
        'omitnan');
    
    prf.iter.Huber.tt(i) = mean(met.loss.t.Huber(idx_t), ...
        'omitnan');
    prf.iter.JKGE.tt(i) = mean(met.performance.t.JKGE(idx_t), ...
        'omitnan');
    
    if has.te
        prf.iter.SAR.te(i) = mean(met.loss.e.SAR(idx_t), ...
            'omitnan');
        prf.iter.GLS.te(i) = mean(met.loss.e.GLS(idx_t), ...
            'omitnan');
        prf.iter.RSS.te(i) = mean(met.loss.e.RSS(idx_t), ...
            'omitnan');
    
        prf.iter.NSE.te(i) = mean(met.performance.e.NSE(idx_t), ...
            'omitnan');
        prf.iter.KGE.te(i) = mean(met.performance.e.KGE(idx_t), ...
            'omitnan');
        
        prf.iter.Huber.te(i) = mean(met.loss.e.Huber(idx_t), ...
            'omitnan');
        prf.iter.JKGE.te(i) = mean(met.performance.e.JKGE(idx_t), ...
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
        prf.iter.SAR.et(i) = mean(met.loss.t.SAR(idx_e), ...
            'omitnan');
        prf.iter.GLS.et(i) = mean(met.loss.t.GLS(idx_e), ...
            'omitnan');
        prf.iter.RSS.et(i) = mean(met.loss.t.RSS(idx_e), ...
            'omitnan');
    
        prf.iter.NSE.et(i) = mean(met.performance.t.NSE(idx_e), ...
            'omitnan');
        prf.iter.KGE.et(i) = mean(met.performance.t.KGE(idx_e), ...
            'omitnan');
        
        prf.iter.Huber.et(i) = mean(met.loss.t.Huber(idx_e), ...
            'omitnan');
        prf.iter.JKGE.et(i) = mean(met.performance.t.JKGE(idx_e), ...
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
        prf.iter.SAR.ee(i) = mean(met.loss.e.SAR(idx_e), ...
            'omitnan');
        prf.iter.GLS.ee(i) = mean(met.loss.e.GLS(idx_e), ...
            'omitnan');
        prf.iter.RSS.ee(i) = mean(met.loss.e.RSS(idx_e), ...
            'omitnan');
    
        prf.iter.NSE.ee(i) = mean(met.performance.e.NSE(idx_e), ...
            'omitnan');   
        prf.iter.KGE.ee(i) = mean(met.performance.e.KGE(idx_e), ...
            'omitnan');   
    
        prf.iter.Huber.ee(i) = mean(met.loss.e.Huber(idx_e), ...
            'omitnan');
        prf.iter.JKGE.ee(i) = mean(met.performance.e.JKGE(idx_e), ...
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
    % Basin- and period-specific dimensionless FDC skill score
    % ----------------------------------------------------------------------
    if ~isfield(loss,'fdc') ...
            || ~isstruct(loss.fdc) ...
            || ~isfield(loss.fdc,'D0t') ...
            || ~isfield(loss.fdc,'D0e') ...
            || numel(loss.fdc.D0t) ~= K ...
            || numel(loss.fdc.D0e) ~= K
    
        error(['      Error:pmetrics: ' ...
            'loss.fdc.D0t or loss.fdc.D0e is missing ' ...
            'or has incorrect size.']);
    end
    
    D0t = double(loss.fdc.D0t(:));
    D0e = double(loss.fdc.D0e(:));
    D0t_t = D0t(idx_t);
    D0e_t = D0e(idx_t);
    if ~isempty(idx_e)
        D0t_e = D0t(idx_e);
        D0e_e = D0e(idx_e);
    else
        D0t_e = [];
        D0e_e = [];
    end
    
    S_fdc.tt = local_fdc_score(D_fdc.tt,D0t_t);
    if has.te
        S_fdc.te = local_fdc_score(D_fdc.te,D0e_t);
    end
    if has.et ...
            && ~isempty(idx_e)
        S_fdc.et = local_fdc_score(D_fdc.et,D0t_e);
    end
    if has.ee ...
            && ~isempty(idx_e)
        S_fdc.ee = local_fdc_score(D_fdc.ee,D0e_e);
    end
    
    % -----------------
    % Summary histories
    % -----------------
    prf.iter.mNSE.tt(i) = local_median(NSE.tt);
    prf.iter.mKGE.tt(i) = local_median(KGE.tt);
    prf.iter.mJKGE.tt(i) = local_median(JKGE.tt);
    prf.iter.Sib_NSE.tt(i) = local_mean_one_minus(NSE.tt);
    prf.iter.Sib_KGE.tt(i) = local_mean_one_minus(KGE.tt);
    prf.iter.Sib_S_fdc.tt(i) = local_mean_one_minus(S_fdc.tt);
    prf.iter.S_fdc.tt(i) = local_mean(S_fdc.tt);
    prf.iter.mS_fdc.tt(i) = local_median(S_fdc.tt);
    
    prf.iter.mNSE.te(i) = local_median(NSE.te);
    prf.iter.mKGE.te(i) = local_median(KGE.te);
    prf.iter.mJKGE.te(i) = local_median(JKGE.te);
    prf.iter.Sib_NSE.te(i) = local_mean_one_minus(NSE.te);
    prf.iter.Sib_KGE.te(i) = local_mean_one_minus(KGE.te);
    prf.iter.Sib_S_fdc.te(i) = local_mean_one_minus(S_fdc.te);
    prf.iter.S_fdc.te(i) = local_mean(S_fdc.te);
    prf.iter.mS_fdc.te(i) = local_median(S_fdc.te);
    
    prf.iter.mNSE.et(i) = local_median(NSE.et);
    prf.iter.mKGE.et(i) = local_median(KGE.et);
    prf.iter.mJKGE.et(i) = local_median(JKGE.et);
    prf.iter.Sib_NSE.et(i) = local_mean_one_minus(NSE.et);
    prf.iter.Sib_KGE.et(i) = local_mean_one_minus(KGE.et);
    prf.iter.Sib_S_fdc.et(i) = local_mean_one_minus(S_fdc.et);
    prf.iter.S_fdc.et(i) = local_mean(S_fdc.et);
    prf.iter.mS_fdc.et(i) = local_median(S_fdc.et);
    
    prf.iter.mNSE.ee(i) = local_median(NSE.ee);
    prf.iter.mKGE.ee(i) = local_median(KGE.ee);
    prf.iter.mJKGE.ee(i) = local_median(JKGE.ee);
    prf.iter.Sib_NSE.ee(i) = local_mean_one_minus(NSE.ee);
    prf.iter.Sib_KGE.ee(i) = local_mean_one_minus(KGE.ee);
    prf.iter.Sib_S_fdc.ee(i) = local_mean_one_minus(S_fdc.ee);
    prf.iter.S_fdc.ee(i) = local_mean(S_fdc.ee);
    prf.iter.mS_fdc.ee(i) = local_median(S_fdc.ee);
    
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
    prf.curr.JKGE = JKGE;

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

function Lblk = local_block_loss(loss_fnc, ...
    SAR,GLS,NSE,KGE,Huber,D_fdc,JKGE)
%LOCAL_BLOCK_LOSS Basin-aggregated loss over one scenario block.

    if isempty(SAR) ...
            && isempty(GLS) ...
            && isempty(NSE) ...
            && isempty(KGE) ...
            && isempty(Huber) ...
            && isempty(D_fdc) ...
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
            Lblk = mean(D_fdc,'omitnan');
        case 7
            Lblk = mean(1 - JKGE,'omitnan');
        otherwise
            error(['      Error:pmetrics: ' ...
                'Unknown loss choice loss = %g.'], ...
                loss_fnc);
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
