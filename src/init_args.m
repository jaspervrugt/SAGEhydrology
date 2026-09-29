function [prf,ax,tTheta,At,An,nTheta] = init_args(bas,mdl,alg,attr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INIT_ARGS Initialize SAGE diagnostics and result arrays.
%
%  Allocates performance histories, parameter traces, and optional
%  attribution arrays.
%
% SYNOPSIS:
%   [prf,ax,tTheta,At,An,nTheta] = init_args(bas,mdl,alg,attr)
%
% INPUT ARGUMENTS:
%   bas             selected basin counts
%    .K_t            number of training basins
%    .K_e            number of evaluation basins
%    .K              total number of basins
%   mdl             assessment settings
%    .mode           basin/period assessment design (1 to 4)
%   alg             optimization settings
%    .i_max          maximum number of training iterations
%   attr            scalar switch for attribution arrays
%
% OUTPUT ARGUMENTS:
%   prf             current metrics and iteration histories
%    .curr           basin-wise metrics for the current iteration
%    .iter           scalar histories across iterations
%   ax              graphics-handle structure
%   tTheta          normalized-parameter percentile traces
%   At              training-basin attribution array
%   An              training-basin net attribution array
%   nTheta          empty normalized-parameter array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K_t = bas.K_t;
    K_e = bas.K_e;
    K = bas.K;
    i_max = alg.i_max;
    
    if isfield(mdl,'th_max') ...
            && ~isempty(mdl.th_max)
        d = numel(mdl.th_max);
    elseif isfield(mdl,'th_max_daily') ...
            && ~isempty(mdl.th_max_daily)
        d = numel(mdl.th_max_daily);
    else
        error(['      Error:init_args: ' ...
            'Could not determine ' ...
            'parameter dimension d.']);
    end
    
    switch mdl.mode
        case 1
            has = struct('tt',true,  ...
                         'te',false, ...
                         'et',false, ...
                         'ee',false);
    
        case 2
            has = struct('tt',true,  ...
                         'te',true,  ...
                         'et',false, ...
                         'ee',false);
    
        case 3
            has = struct('tt',true,  ...
                         'te',false, ...
                         'et',true,  ...
                         'ee',false);
    
        case 4
            has = struct('tt',true, ...
                         'te',true, ...
                         'et',true, ...
                         'ee',true);
    
        otherwise
            error(['      Error:init_args: ' ...
                'mdl.mode must be one of 1,2,3,4.']);
    end
    
    tTheta = nan(i_max,d,7);
    
    ax = struct();
    ax.inited = false;
    
    prf = struct();
    prf.mode = mdl.mode;
    prf.K_t = K_t;
    prf.K_e = K_e;
    prf.K = K;
    prf.has = has;
    prf.curr = struct();
    prf.iter = struct();
    
    stats = {'L','SAR','GLS', ...
        'RSS','Huber','NSE','KGE','JKGE','S_fdc','S_p','S_logp', ...
        'mNSE','mKGE','mJKGE','mS_fdc','mS_p','mS_logp', ...
        'Sib_NSE','Sib_KGE','Sib_JKGE','Sib_S_fdc', ...
        'Sib_S_p','Sib_S_logp', ...
        'KGE_r','mKGE_r', ...
        'KGE_alpha','mKGE_alpha', ...
        'KGE_beta','mKGE_beta', ...
        'JKGE_M','mJKGE_M', ...
        'JKGE_V','mJKGE_V', ...
        'JKGE_C','mJKGE_C'};
    sets  = {'tt','te','et','ee'};
    
    for s = 1:numel(stats)
        for q = 1:numel(sets)
            prf.iter.(stats{s}).(sets{q}) = nan(1,i_max);
        end
    end
    
    prf.iter.cpuT = nan(1,i_max);
    % Learning rate used to create the phi evaluated at iteration i.
    % Entry 1 remains NaN because phi_1 is created by ffn_theta('init').
    % For i > 1, assign prf.iter.lr(i) = opts.lr_current after descent('dyn').
    prf.iter.lr = nan(1,i_max);
    prf.iter.gradNormPre = nan(1,i_max);
    prf.iter.gradNormPost = nan(1,i_max);
    prf.iter.updateNorm = nan(1,i_max);
    prf.iter.relativeUpdateNorm = nan(1,i_max);
    prf.iter.moment1Norm = nan(1,i_max);
    prf.iter.moment2Norm = nan(1,i_max);
    prf.iter.weightDecayNorm = nan(1,i_max);
    prf.iter.gradientClipped = false(1,i_max);
    
    if ~isempty(attr) && isequal(attr,1)
        At = nan(d,K_t,i_max);
        An = nan(d,K_t,i_max);
    else
        At = [];
        An = [];
    end    
    nTheta = [];

end
