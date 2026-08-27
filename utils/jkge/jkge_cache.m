function cache = jkge_cache(y,method,n_win,mo)
%JKGE_CACHE Precompute JKGE benchmark bookkeeping

y = double(y(:));
n = numel(y);

cache = struct();
cache.method = method;
cache.n = n;
cache.good = logical(isfinite(y));

switch method
    case 1
        h = floor((n_win - 1)/2);
        idx = (1:n)';

        cache.L = uint32(max(1,idx - h));
        cache.R = uint32(min(n,idx + h));

        cnt = [0; cumsum(double(cache.good))];
        cache.den = single(cnt( ...
            double(cache.R) + 1) - ...
            cnt(double(cache.L)));

    case 2
        gid = ceil((1:n)' / n_win);
        cache.gid = uint32(gid);
        cache.G = double(gid(end));
        cache.cnt_g = single( ...
            accumarray(gid(cache.good), ...
            1,[cache.G 1],@sum,0));

    case 3
        cache.ng = nnz(cache.good);

    case 4
        mo = double(mo(:));

        valid_mo = isfinite(mo) ...
            & mo >= 1 ...
            & mo <= 12 ...
            & mod(mo,1) == 0;

        gid = double(mo);
        cache.gid = uint8(gid);
        cache.G = 12;
        cache.good = logical(cache.good ...
            & valid_mo);
        cache.cnt_g = single( ...
            accumarray(gid(cache.good), ...
            1,[cache.G 1],@sum,0));

    otherwise
        error(['      Error:jkge_cache: ' ...
            'unknown JKGE method = %g.'], ...
            method);
end
end