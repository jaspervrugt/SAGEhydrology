function loss = fdc_loss(y_n,q_n)

    n = numel(y_n);

    method = 3;

    % method = 2; % for automatic differentiation!

    switch method
        % Original implementation: O(n^2) scaling
        case 1
            part1 = 0; part2 = 0;
            for i = 1:n
                part1 = part1 + sum(abs(q_n(i) - y_n));
            end
            for i = 1:n
                part2 = part2 + sum((abs(q_n(i) - q_n) + abs(y_n(i) - ...
                    y_n)));
            end
            loss = 1/n^2 * part1 - 1/2 * 1/n^2 * part2;

        % AD implementation: O(n^2) but 3 matrices (= not large n)
        case 2

            % n x n: Pairwise absolute differences |q_i - y_j|
            Dyq = abs(q_n - y_n.');
            % n x n: Pairwise |q_i - q_j| and |y_i - y_j|:
            Dqq = abs(q_n - q_n.');
            Dyy = abs(y_n - y_n.');     % n x n
            % Use sums so it works both for double and dlarray
            part1 = sum(Dyq(:));        % sum_{i,j} |q_i - y_j|
            % sum_{i,j}(|q_i - q_j| + |y_i - y_j|)
            part2 = sum(Dqq(:)) + sum(Dyy(:));
            loss = (1/n^2) * part1 - 0.5 * (1/n^2) * part2;

        case 3                          % O(nlog(n))
            % Sort once
            ys = sort(y_n); qs = sort(q_n);
            loss = fdc_cdf_distance_sorted(qs,ys);
    end

end
