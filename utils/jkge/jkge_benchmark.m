function m_y = jkge_benchmark(y,method,n_win,mo,cache)
%JKGE_BENCHMARK Matrix-free JKGE benchmark computation
%
% SYNOPSIS:  b = jkge_benchmark(y,method,n_win,mo)
%            b = jkge_benchmark(y,method,n_win,mo,cache)
%   y        nx1 vector observed or simulated discharge
%   method   scalar JKGE benchmark definition method
%             1 = moving-average mean
%             2 = section-wise mean
%             3 = long-term mean
%             4 = monthly climatology
%   n_win   scalar centered moving window in days: methods 1&2
%   mo      nx1 month labels: method 4, empty, otherwise
%   cache   OPTIONMAL: structure with info to speed-up jkge computation
%   m_y     OUTPUT: nx1 mean benchmark vector
%
% NaN gap shorter than/near window  -> benchmark uses nearby finite values
% NaN gap much longer than window   -> interior benchmark = NaN
% JKGE scoring                      -> those NaN locations are skipped
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    cache = [];
end

y = double(y(:));
n = numel(y);

m_y = nan(n,1);

switch method
    case 1
        if isempty(cache)
            cache = jkge_cache(y, ...
                method,n_win,mo);
        end
        good = isfinite(y);
        y0 = y;
        y0(~good) = 0;
        cs = [0; cumsum(y0)];
        num = cs(double(cache.R) + 1) ...
            - cs(double(cache.L));
        den = double(cache.den);
        ok = den > 0;
        m_y(ok) = num(ok) ./ den(ok);

    case 2
        if isempty(cache)
            cache = jkge_cache(y, ...
                method,n_win,mo);
        end
        good = isfinite(y);
        gid = double(cache.gid);
        G = cache.G;
        sum_g = accumarray(gid(good), ...
            y(good),[G 1],@sum,0);
        mu_g = nan(G,1);
        ok_g = cache.cnt_g > 0;
        mu_g(ok_g) = sum_g(ok_g) ./ cache.cnt_g(ok_g);
        m_y = mu_g(gid);

    case 3
        good = isfinite(y);
        ng = nnz(good);
        if ng > 0
            m_y(:) = sum(y(good)) / ng;
        end

    case 4
        if isempty(cache)
            cache = jkge_cache(y, ...
                method,n_win,mo);
        end
        good = cache.good;
        gid = double(cache.gid);
        G = cache.G;
        sum_g = accumarray(gid(good), ...
            y(good),[G 1],@sum,0);
        mu_g = nan(G,1);
        ok_g = cache.cnt_g > 0;
        mu_g(ok_g) = sum_g(ok_g) ./ cache.cnt_g(ok_g);
        ok = good;
        m_y(ok) = mu_g(gid(ok));

    otherwise
        error(['      Error:jkge_benchmark: ' ...
            'unknown JKGE method = %g.'], ...
            method);
end

end