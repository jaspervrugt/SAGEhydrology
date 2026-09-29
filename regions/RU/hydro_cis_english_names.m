function names = hydro_cis_english_names(sourceNames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HYDRO_CIS_ENGLISH_NAMES Format HydroCIS gauge names for display.
%
% SYNOPSIS:
%   names = hydro_cis_english_names(sourceNames)
%
% INPUT:
%   sourceNames   Array of official transliterated HydroCIS gauge names
%
% OUTPUT:
%   names         English-readable names without Russian abbreviations
%
% DESCRIPTION:
%   HydroCIS transliterations use Russian feature and settlement prefixes,
%   for example r. (river), ruch. (stream), s. (village), and pos.
%   (settlement). This function preserves the transliterated proper names
%   while converting the source notation to concise English display names.
%   The unmodified Cyrillic source name is retained separately by the
%   installer in gauge_information.txt.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Sept. 2026                               %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sourceNames = string(sourceNames(:));
    names = sourceNames;
    locationPrefix = [ ...
        '^(?:zh\.-d\.st\.|gm\.st\.|gm\.p\.|met\.st\.|' ...
        'm\.st\.|k\.p\.|r\.p\.|pgt\.?|pos\.|slob\.|' ...
        'uroch\.|fakt\.|svh\.|rzd\.|st\.|s\.|d\.|g\.|' ...
        'h\.|p\.|sl\.|z\.|por\.|gm\.|met\.|c\.|pldp\.|' ...
        'm\.|k\.|uch\.)\s*'];

    for i = 1:numel(sourceNames)
        source = strtrim(char(sourceNames(i)));
        separator = regexp(source,'\s+-\s+','once');
        if isempty(separator)
            continue
        end

        feature = strtrim(source(1:separator-1));
        location = regexprep(source(separator:end),'^\s+-\s+','');
        if startsWith(feature,'ruch.','IgnoreCase',true)
            feature = [strtrim(feature(6:end)) ' Stream'];
        elseif startsWith(feature,'r.','IgnoreCase',true)
            feature = [strtrim(feature(3:end)) ' River'];
        elseif startsWith(feature,'balka ','IgnoreCase',true)
            feature = [strtrim(feature(7:end)) ' Gully'];
        end

        % Remove source-language administrative prefixes. The location
        % proper name is sufficient after the English preposition "at".
        location = regexprep(location,locationPrefix,'', ...
            'ignorecase');
        location = regexprep(location, ...
            '^(?:GMS|GMP|GP|GES|st-tsa)\s+','', ...
            'ignorecase');

        % Translate the recurring descriptive locations whose meaning is
        % unambiguous. Less common source descriptions remain transliterated
        % rather than receiving a speculative translation.
        location = regexprep(location, ...
            '^v\s+([0-9.]+)\s+km\s+ot\s+ust''ja$', ...
            '$1 km upstream from the mouth','ignorecase');
        location = regexprep(location, ...
            '^([0-9.]+)\s+km\s+ot\s+ust''ja$', ...
            '$1 km upstream from the mouth','ignorecase');
        location = regexprep(location, ...
            '^v\s+([0-9.]+)\s+m\s+ot\s+ust''ja$', ...
            '$1 m upstream from the mouth','ignorecase');
        location = regexprep(location,'^ust''e$', ...
            'the mouth','ignorecase');
        location = regexprep(location,'^istok$', ...
            'the source','ignorecase');
        location = regexprep(location, ...
            '^vyshe\s+ust''ja\s+ruch\.(.+)$', ...
            'upstream from the confluence with $1 Stream','ignorecase');
        location = regexprep(location, ...
            '^vyshe\s+ust''ja\s+r\.(.+)$', ...
            'upstream from the confluence with $1 River','ignorecase');
        location = regexprep(location, ...
            '^ust''e\s+r\.(.+)$', ...
            'the confluence with $1 River','ignorecase');

        % Expand or remove the same source abbreviations when they occur
        % inside a descriptive location or a parenthetical former name.
        location = regexprep(location,'met\.st\.([A-Za-z0-9''-]+)', ...
            '$1 meteorological station','ignorecase');
        location = regexprep(location,'g\.st\.','gauging station', ...
            'ignorecase');
        location = regexprep(location,'zh\.d\.','railway ', ...
            'ignorecase');
        location = regexprep(location,'\br\.([A-Za-z0-9''-]+)', ...
            '$1 River','ignorecase');
        location = regexprep(location, ...
            '\b(?:pos|pgt|slob|uroch|fakt|svh|rzd|pldp|uch|st|s|d|g|h)\.', ...
            '','ignorecase');
        location = regexprep(location,'\s+',' ');

        if any(startsWith(string(location), ...
                ["upstream ","downstream "],'IgnoreCase',true))
            names(i) = string([feature ' ' strtrim(location)]);
        else
            names(i) = string([feature ' at ' strtrim(location)]);
        end
    end

    names = reshape(names,size(sourceNames));
end
