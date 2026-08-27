function s = normalize_usgs(x)
%NORMALIZE_USGS Normalize USGS basin IDs to 8-character strings
%
% Rules:
%   - numeric inputs are rounded and formatted as 8 digits
%   - string/cell inputs are trimmed
%   - Excel-style trailing ".0" is removed
%   - purely numeric strings shorter than 8 are left-padded with zeros
%   - empty/missing entries become ""
%   - suspicious nonnumeric text is left as-is after trimming, except that
%     internal spaces are removed only if the remaining text is all digits

if isnumeric(x) ...
        || islogical(x)
    x = double(x(:));

    s = strings(size(x));
    good = isfinite(x);

    if any(good)
        xr = round(x(good));

        if any(abs(x(good) - xr) > 1e-9)
            warning(['      Warning: normalize_usgs:' ...
                '    NonIntegerNumericID: ' ...
                'Some numeric USGS IDs are ' ...
                'not integers; rounding used.']);
        end

        s(good) = string(compose('%08d',xr));
    end

    s(~good) = "";
    return
end

s = string(x(:));
s = strtrim(s);
s(ismissing(s)) = "";

for jj = 1:numel(s)
    sj = s(jj);

    if strlength(sj) == 0
        continue
    end

    % remove Excel-style trailing ".0", ".00", etc.
    sj = regexprep(sj,'\.0+$','');

    % if spaces exist inside, remove only if remaining content is digits
    sj_nospace = regexprep(sj,'\s+','');
    if ~isempty(regexp(sj_nospace, ...
            '^\d+$','once'))
        sj = sj_nospace;
    end

    % purely numeric string -> left pad to 8 if needed
    if ~isempty(regexp(sj,'^\d+$', ...
            'once'))
        if strlength(sj) < 8
            sj = string([repmat('0',1, ...
                8-strlength(sj)) char(sj)]);
        end
        % if strlength(sj) < 8
        %     sj = pad(sj,8,'left','0');
        % end
    end

    s(jj) = sj;
end

s = s(:);
end