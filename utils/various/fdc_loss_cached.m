function loss = fdc_loss_cached(q_n,fdc)

    % Number of valid observations
    n = fdc.n;
    
    if n == 0 ...
            || isempty(q_n)
        loss = NaN;
        return
    end
    
    % Cached observed-discharge quantities
    % ys = fdc.ys;
    % Py = fdc.Py;
    % S_yy = fdc.S_yy;
    
    % Only simulated discharge changes between SAGE iterations
    qs = sort(q_n);
    
    % Cross term:
    % S_qy = sum_{i,j} |q_i - y_j|
    S_qy = sumAbsPairsCross_sorted(qs,fdc.ys,fdc.Py);
    
    % Simulated-discharge within-sample term:
    % S_qq = sum_{i,j} |q_i - q_j|
    S_qq = sumAbsPairsSame_sorted(qs);
    
    % Information-theoretic FDC distance
    loss = (S_qy - 0.5*(S_qq + fdc.S_yy)) / n^2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Secondary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = sumAbsPairsSame_sorted(xs)
% xs is sorted
    n = numel(xs);
    w = (2*(1:n)' - n - 1);
    S = 2 * sum(w .* xs);        % sum_{i,j} |x_i-x_j|
end

function S = sumAbsPairsCross_sorted(qs, ys, Py)
% qs, ys sorted. Py = [0; cumsum(ys)]
    n = numel(ys);
    j = 0;                       % number of ys <= current q
    
    S = 0;
    sumY = Py(end);
    
    for i = 1:n
        qi = qs(i);
        while (j < n) ...
                && (ys(j+1) <= qi)
            j = j + 1;
        end
        % left: y <= q
        left = j * qi - Py(j+1);
        % right: y > q
        right = (sumY - Py(j+1)) - (n - j) * qi;
        S = S + left + right;
    end
end