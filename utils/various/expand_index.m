function id = expand_index(id_in)
%EXPAND_INDEX Expand compact [first last] index pair or accept full vector
%
% Accepted forms:
%   [i0 i1]              -> i0:i1
%   [i1 i2 i3 ... in]    -> unchanged, provided strictly increasing
%
% Output:
%   1xn row vector of integer indices

if nargin < 1 || isempty(id_in)
    id = [];
    return
end

id_in = double(id_in(:)).';

if any(~isfinite(id_in)) ...
        || any(abs(id_in - round(id_in)) > 0)
    error(['      Error:expand_index: ' ...
        'indices must be finite ' ...
        'integers.']);
end

if numel(id_in) == 2
    if id_in(2) < id_in(1)
        error(['      Error:expand_index: ' ...
            'invalid index range ' ...
            '[%d %d].'],id_in(1),id_in(2));
    end
    id = id_in(1):id_in(2);
else
    if any(diff(id_in) <= 0)
        error(['      Error:expand_index: ' ...
            'full index vector must be ' ...
            'strictly increasing.']);
    end
    id = id_in;
end

end