function loss = fdc_loss(y_n,q_n)

    n = numel(y_n);

    method = 3;

    % method = 2; % for automatic differentiation!

    switch method
        case 1 % Original implementation: O(n^2) scaling
            part1 = 0; part2 = 0;
            for i = 1:n
                part1 = part1 + sum(abs(q_n(i) - y_n));
            end
            for i = 1:n
                part2 = part2 + sum((abs(q_n(i) - q_n) + abs(y_n(i) - y_n)));
            end
            loss = 1/n^2 * part1 - 1/2 * 1/n^2 * part2;

        case 2 % AD implementation: O(n^2) but 3 matrices (= not large n)

            Dyq = abs(q_n - y_n.');     % n x n: Pairwise absolute differences |q_i - y_j|
            Dqq = abs(q_n - q_n.');     % n x n: Pairwise |q_i - q_j| and |y_i - y_j|: 
            Dyy = abs(y_n - y_n.');     % n x n
            % Use sums so it works both for double and dlarray
            part1 = sum(Dyq(:));        % sum_{i,j} |q_i - y_j|
            part2 = sum(Dqq(:)) + sum(Dyy(:));  % sum_{i,j}(|q_i - q_j| + |y_i - y_j|)
            loss = (1/n^2) * part1 - 0.5 * (1/n^2) * part2;

        case 3 % O(nlog(n))
            % Sort once
            ys = sort(y_n); qs = sort(q_n);
            % Prefix sums
            Py = [0; cumsum(ys)]; %Pq = [0; cumsum(qs)];       % (n+1)x1
            % Cross term S_qy = sum_{i,j} |q_i - y_j|
            S_qy = sumAbsPairsCross_sorted(qs, ys, Py);
            % Within-sample terms S_qq and S_yy = sum_{i,j}|x_i-x_j|
            S_qq = sumAbsPairsSame_sorted(qs);
            S_yy = sumAbsPairsSame_sorted(ys);
            loss = (1/n^2) * S_qy - 0.5 * (1/n^2) * (S_qq + S_yy);
    end

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
        while (j < n) && (ys(j+1) <= qi)
            j = j + 1;
        end
        % left: y <= q
        left  = j * qi - Py(j+1);
        % right: y > q
        right = (sumY - Py(j+1)) - (n - j) * qi;
        S = S + left + right;
    end
end