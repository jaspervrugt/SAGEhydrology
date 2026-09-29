function met = hydro_metrics(qsim,qobs,idx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HYDRO_METRICS Compute discharge performance and flow metrics.
%
%   met = hydro_metrics(qsim,qobs)
%   met = hydro_metrics(qsim,qobs,idx)
%
% qsim is the n-by-m discharge returned by crr_model (one simulation per
% column); qobs is an n-vector of measured discharge (for example y_n).
% idx is an optional logical mask or vector of indices for a train/eval
% period. The 1-by-m output structure contains n, ME, MAE, RMSE, NRMSE,
% PBIAS, R, R2, NSE, KGE and its r/alpha/beta components, and observed,
% simulated, and percent-bias values for Q05, Q50, and Q95 exceedance flow.
%
% Nonfinite pairs are removed independently for each simulation. Undefined
% ratios and efficiencies are NaN. JKGE is deliberately excluded because
% it requires benchmark vectors and method metadata prepared elsewhere.
% This standalone function is not currently called by camels or crr_model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    narginchk(2,3);
    validateattributes(qsim,{'numeric'}, ...
        {'real','2d'}, ...
        mfilename,'qsim',1);
    validateattributes(qobs,{'numeric'}, ...
        {'real','vector'}, ...
        mfilename,'qobs',2);
    yall = double(qobs(:));
    qall = double(qsim);
    if isvector(qall)
        qall = qall(:);
    elseif size(qall,1) ~= numel(yall) ...
            && size(qall,2) == numel(yall)
        qall = qall.';
    end
    if size(qall,1) ~= numel(yall)
        error('hydro_metrics:LengthMismatch', ...
            ['qsim must have one ' ...
            'row per qobs value.']);
    end

    if nargin < 3 ...
            || isempty(idx)
        idx = (1:numel(yall)).';
    elseif islogical(idx)
        if numel(idx) ~= numel(yall)
            error('hydro_metrics:MaskLengthMismatch', ...
                ['Logical idx must have ' ...
                'the same length as qobs.']);
        end
        idx = find(idx(:));
    else
        validateattributes(idx,{'numeric'}, ...
            {'real','vector', ...
            'integer','positive'}, ...
            mfilename,'idx',3);
        idx = double(idx(:));
        if any(idx > numel(yall))
            error('hydro_metrics:IndexOutOfRange', ...
                'idx exceeds numel(qobs).');
        end
    end
    yall = yall(idx);
    qall = qall(idx,:);

    met = repmat(local_empty(),1,size(qall,2));
    for j = 1:size(qall,2)
        y = yall; q = qall(:,j);
        use = isfinite(y) ...
            & isfinite(q);
        y = y(use); 
        q = q(use);
        met(j).n = numel(y);
        if isempty(y)
            continue
        end

        e = q-y;
        my = mean(y); 
        mq = mean(q);
        sy = std(y,0); 
        sq = std(q,0);
        met(j).ME = mean(e);
        met(j).MAE = mean(abs(e));
        met(j).RMSE = sqrt(mean(e.^2));
        met(j).NRMSE = local_ratio(met(j).RMSE,my);
        met(j).PBIAS = 100*local_ratio(sum(e),sum(y));
        if numel(y) >= 2 ...
                && sy > 0 ...
                && sq > 0
            C = corrcoef(y,q);
            met(j).R = C(1,2);
            met(j).R2 = met(j).R^2;
        end
        tss = sum((y-my).^2);
        if tss > 0
            met(j).NSE = 1-sum(e.^2)/tss; 
        end

        met(j).KGE_r = met(j).R;
        met(j).KGE_alpha = local_ratio(sq,sy);
        met(j).KGE_beta = local_ratio(mq,my);
        c = [met(j).KGE_r met(j).KGE_alpha met(j).KGE_beta];
        if all(isfinite(c))
            met(j).KGE = 1-sqrt(sum((c-1).^2));
        end

        % Q05 is high flow and Q95 is low flow (exceedance convention).
        met(j).Q05obs = local_pct(y,95); 
        met(j).Q05sim = local_pct(q,95);
        met(j).Q50obs = local_pct(y,50); 
        met(j).Q50sim = local_pct(q,50);
        met(j).Q95obs = local_pct(y,5);
        met(j).Q95sim = local_pct(q,5);
        for nm = {'Q05','Q50','Q95'}
            s = nm{1};
            met(j).([s 'bias']) = 100*local_ratio( ...
                met(j).([s 'sim'])-met(j).([s 'obs']),met(j).([s 'obs']));
        end
    end
end

function s = local_empty()
    names = {'ME','MAE','RMSE','NRMSE', ...
        'PBIAS','R','R2','NSE','KGE', ...
        'KGE_r','KGE_alpha','KGE_beta', ...
        'Q05obs','Q05sim','Q05bias', ...
        'Q50obs','Q50sim','Q50bias', ...
        'Q95obs','Q95sim','Q95bias'};
    s = struct('n',0);
    for k = 1:numel(names)
        s.(names{k}) = NaN; 
    end
end

function z = local_ratio(a,b)
    if isfinite(a) ...
            && isfinite(b) ...
            && b ~= 0
        z = a/b; 
    else
        z = NaN; 
    end
end

function z = local_pct(x,p)
    x = sort(x(:)); 
    n = numel(x);
    if n == 1
        z = x; 
        return
    end
    pos = 1+(n-1)*p/100; 
    lo = floor(pos); 
    hi = ceil(pos);
    z = x(lo)+(pos-lo)*(x(hi)-x(lo));
end
