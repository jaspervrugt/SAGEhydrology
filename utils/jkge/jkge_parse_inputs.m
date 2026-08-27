function [idx,method,n_win,mo,Mdef] = jkge_parse_inputs(n,varargin)
%JKGE_PARSE_INPUTS Parse optional valid-index vector and JKGE inputs
%
% SYNOPSIS: [idx,method,n_win,mo,Mdef] = jkge_parse_inputs(n,varargin)
%
% Optional trailing M selector:
%   M = 1 original ratio form from the paper
%   M = 2 norm-based formulation

method = NaN;
n_win = [];
mo = [];
Mdef = 2;

if isempty(varargin)
    error(['      Error:jkge: ' ...
        'missing JKGE benchmark inputs.']);
end

first = varargin{1};

% Optional idx input
isIdxInput = false;

if islogical(first) ...
        && numel(first) == n
    isIdxInput = true;

elseif isnumeric(first)

    if ~isscalar(first)
        isIdxInput = true;

    elseif numel(varargin) >= 2 ...
            && isnumeric(varargin{2}) ...
            && isscalar(varargin{2}) ...
            && any(double(varargin{2}) == [1 2 3 4])
        % scalar valid index followed by method
        isIdxInput = true;
    end
end

if isIdxInput

    if islogical(first)
        idx = find(first(:));
    else
        idx = double(first(:));
    end

    if isempty(idx)
        idx = zeros(0,1);
    end

    if any(~isfinite(idx)) ...
            || any(idx < 1) ...
            || any(idx > n) ...
            || any(mod(idx,1) ~= 0)
        error(['      Error:jkge: ' ...
            'idx must contain valid ' ...
            'integer indices.']);
    end

    idx = unique(idx(:),'stable');

    if isempty(idx)
        return
    end

    if numel(varargin) < 2
        error(['      Error:jkge: ' ...
            'method missing after idx.']);
    end

    method = varargin{2};
    aux = varargin(3:end);

else
    idx = (1:n)';
    method = varargin{1};
    aux = varargin(2:end);
end

if ~(isscalar(method) ...
        && isnumeric(method) ...
        && isfinite(method) ...
        && any(method == [1 2 3 4]))
    error(['      Error:jkge: ' ...
        'method must be 1, 2, 3, or 4.']);
end

method = double(method);

switch method
    case {1,2}
        if numel(aux) < 1
            error(['      Error:jkge: ' ...
                'method %d requires n_win.'], ...
                method);
        end

        n_win = aux{1};

        if ~(isscalar(n_win) ...
                && isnumeric(n_win) ...
                && isfinite(n_win) ...
                && n_win >= 1 ...
                && mod(n_win,1) == 0)
            error(['      Error:jkge: ' ...
                'n_win must be a ' ...
                'positive integer.']);
        end

        n_win = double(n_win);

        if method == 1 ...
                && mod(n_win,2) == 0
            n_win = n_win + 1;
        end

        if numel(aux) >= 2
            Mdef = local_parse_M(aux{2});
        end

    case 3
        n_win = [];

        if numel(aux) >= 1
            Mdef = local_parse_M(aux{1});
        end

    case 4
        if numel(aux) < 1
            error(['      Error:jkge: ' ...
                'method 4 requires ' ...
                'month labels.']);
        end

        mo = double(aux{1}(:));

        if numel(mo) ~= n
            error(['      Error:jkge: ' ...
                'month vector length ' ...
                'does not match y/q length.']);
        end

        if any(~isfinite(mo)) ...
                || any(mo < 1 ...
                | mo > 12 ...
                | mod(mo,1) ~= 0)
            error(['      Error:jkge: ' ...
                'month labels must be ' ...
                'integers in {1,...,12}.']);
        end

        if numel(aux) >= 2
            Mdef = local_parse_M(aux{2});
        end
end
end

function Mdef = local_parse_M(Mdef)
if ~(isscalar(Mdef) ...
        && isnumeric(Mdef) ...
        && isfinite(Mdef) ...
        && any(Mdef == [1 2]))
    error(['      Error:jkge: ' ...
        'M must be 1 or 2.']);
end

Mdef = double(Mdef);
end
