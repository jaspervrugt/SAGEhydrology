function z = jkge_transpose(v,q,method,n_win,mo,cache)
%JKGE_TRANSPOSE Apply transpose of JKGE benchmark operator without forming 
% the operator explicitly.
%  SYNOPSIS: z = jkge_transpose(v,q,method,n_win,mo)
%   v       nx1 vector to which B' is applied
%   q       nx1 simulated discharge on the same support
%   method  scalar JKGE benchmark method
%   n_win   scalar window/block size for methods 1 and 2
%   mo      nx1 month labels for method 4
%   cache   OPTIONAL: structure with info to speed-up jkge computation
%   z       OUTPUT: nx1 result of B' * v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, April 2026                                %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    cache = [];
end

q = double(q(:));

if isempty(cache)
    cache = jkge_cache(q,method,n_win,mo);
end

v = double(v(:));
n = numel(q);

if numel(v) ~= n
    error(['      Error:jkge_transpose: ' ...
        'v and q must have the same length.']);
end

z = zeros(n,1);

switch method
    case 1
        den = double(cache.den);
        L = double(cache.L);
        R = double(cache.R);
        good = cache.good;
        ok = isfinite(v) ...
            & den > 0;
        w = zeros(n,1);
        w(ok) = v(ok) ./ ...
            den(ok);
        acc = zeros(n + 1,1);
        ii = find(ok);
        for k = 1:numel(ii)
            j = ii(k);
            acc(L(j)) = ...
                acc(L(j)) + w(j);
            acc(R(j) + 1) = ...
                acc(R(j) + 1) - w(j);
        end
        z = cumsum(acc(1:n));
        z(~good) = 0;

    case 2
        gid = double(cache.gid);
        G = double(cache.G);
        good = cache.good;
        cnt_g = double(cache.cnt_g);
        use = isfinite(v) ...
            & good;
        sumv_g = accumarray(gid(use), ...
            v(use),[G 1],@sum,0);
        ok = good & cnt_g(gid) > 0;
        z(ok) = sumv_g(gid(ok)) ./ ...
            cnt_g(gid(ok));

    case 3
        good = cache.good;
        ng = nnz(good);
        if ng > 0
            s = sum(v(isfinite(v) ...
                & good));
            z(good) = s / ng;
        end

    case 4
        gid = double(cache.gid);
        G = double(cache.G);
        good = cache.good;
        cnt_g = double(cache.cnt_g);
        use = isfinite(v) ...
            & good;
        sumv_g = accumarray(gid(use), ...
            v(use),[G 1],@sum,0);
        ok = good & cnt_g(gid) > 0;
        z(ok) = sumv_g(gid(ok)) ./ ...
            cnt_g(gid(ok));

    otherwise
        error(['      Error:jkge_transpose: ' ...
            'unknown JKGE method = %g.'],method);
end

end