function s = standardize_camels_gauge_names(s,region,id_gauge)
%STANDARDIZE_CAMELS_GAUGE_NAMES Standardize CAMELS gauge/station names.
%
% SYNOPSIS:
%   s = standardize_camels_gauge_names(s,region)
%   s = standardize_camels_gauge_names(s,region,id_gauge)
%
% INPUT:
%   s          Raw gauge/station names as string/cell/char/categorical
%   region     Region code, e.g. 'CAMELS_AU', 'AU', 'CAMELS_CH', 'CH'
%   id_gauge   Optional gauge IDs, used only for region-specific fallback
%              names and for removing appended IDs in some datasets.
%
% OUTPUT:
%   s          Cleaned string column vector.
%
% NOTES:
%   This helper replaces the small local name-cleaning functions in the
%   regional read_attr_XX readers. It intentionally does only display-name
%   cleaning. It does not choose which source column supplies the name.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, May 2026                                  %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(region)
        region = "";
    end

    if ischar(s) && (isrow(s) ...
            || isempty(s))
        s = string({s});
    else
        s = string(s);
    end
    s = s(:);

    if nargin < 3
        id_gauge = strings(numel(s),1);
    end

    if ischar(id_gauge) && (isrow(id_gauge) ...
            || isempty(id_gauge))
        id_gauge = string({id_gauge});
    else
        id_gauge = string(id_gauge);
    end
    id_gauge = id_gauge(:);

    if isscalar(id_gauge) && numel(s) > 1
        id_gauge = repmat(id_gauge,numel(s),1);
    elseif numel(id_gauge) ~= numel(s)
        id_gauge = strings(numel(s),1);
    end

    reg = upper(string(region));
    reg = erase(reg,"CAMELSH_");
    reg = erase(reg,"CAMELS_");
    reg = erase(reg,"READ_ATTR_");
    reg = strip(reg);
    if reg == "MACH_US"
        reg = "US";
    end

    % -------------------------
    % Region-specific defaults
    % -------------------------
    opt.title_case = false;
    opt.title_case_upper_only = false;
    opt.spanish_title = false;
    opt.canadian_title = false;
    opt.small_words = strings(0,1);
    opt.replace_hyphen = false;
    opt.remove_parenthetical_id = false;
    opt.remove_appended_id = false;
    opt.fallback_prefix = "Gauge ";

    switch char(reg)
        case 'AT'
            opt.fallback_prefix = "AT_";

        case 'AU'
            opt.remove_parenthetical_id = true;

        case 'BR'
            opt.title_case = true;
            opt.small_words = ["da","das","de","do","dos","e"];

        case 'CA'
            opt.title_case = true;
            opt.canadian_title = true;
            opt.small_words = ["at","of","the","and","in","on","by", ...
                "near","above","below","under","over","to","from"];

        case 'CH'
            opt.title_case = true;
            opt.replace_hyphen = true;
            opt.small_words = ["au","aux","de","des", ...
                "du","d","la","le","les","l", ...
                "et","en","sur","sous"];
        
        case 'CL'
            opt.title_case = true;
            opt.spanish_title = true;
            opt.replace_hyphen = true;

        case 'COL'
            opt.title_case = true;
            opt.spanish_title = true;
            opt.remove_appended_id = true;
            opt.replace_hyphen = true;
            opt.fallback_prefix = "COL_";

        case 'DE'
            opt.title_case = true;
            opt.replace_hyphen = true;

        case 'DK'
            % Preserve source capitalization where supplied, but normalize
            % separators and use a Denmark-specific fallback for missing
            % station names.
            opt.replace_hyphen = true;
            opt.fallback_prefix = "DK_";

        case 'ES'
            opt.title_case = true;
            opt.spanish_title = true;
            opt.replace_hyphen = true;
            opt.fallback_prefix = "ES_";

        case 'LUX'
            % Preserve source capitalization; 
            % only clean spaces/underscores.
        
        case {'FI','FR','NZ'}
            % Preserve source capitalization; 
            % only clean spaces/underscores.
        
        case 'PL'
            % Preserve source capitalization; 
            % use PL_ fallback if name missing.
            opt.fallback_prefix = "PL_";

        case 'GB'
            opt.title_case = true;
            opt.replace_hyphen = true;
            opt.small_words = ["at","of","the", ...
                "and","in","on","by","near"];
        
        case {'IND','SE'}
            opt.title_case = true;
            opt.replace_hyphen = true;

        case 'US'
            % The two US products combine station-name sources with mixed
            % capitalization and USGS abbreviations. Preserve names already
            % in mixed case, but title-case fully capitalized entries.
            opt.title_case = true;
            opt.title_case_upper_only = true;
            opt.small_words = ["at","of","the","and","in","on", ...
                "by","near","above","below","to","from"];

        otherwise
            % Conservative default: clean separators but preserve case.
    end

    % -----------------
    % Generic cleaning
    % -----------------
    s = strip(s);
    s(ismissing(s)) = "";

    % Remove appended gauge ID before general separator cleanup. This covers
    % cases such as "Caldas_26157020 0" in CAMELS-COL.
    if opt.remove_appended_id
        for k = 1:numel(s)
            id = strip(id_gauge(k));
            if id == "" || ismissing(id)
                continue
            end
            idpat = regexptranslate('escape',char(id));
            s(k) = regexprep(s(k), ...
                ['[_\s]*' idpat '\s*0?\s*$'],'');
        end
    end

    s = replace(s,"_"," ");

    if opt.replace_hyphen
        s = replace(s,"-"," ");
    end

    s = regexprep(s,'\s+',' ');
    s = strip(s);

    % Remove trailing station ID in parentheses, e.g. "River Name (136208A)".
    if opt.remove_parenthetical_id
        s = regexprep(s,'\s*\([A-Za-z0-9]+\)\s*$','');
        s = strip(s);
    end

    bad = ismissing(s) ...
        | strlength(s) == 0;
    if any(bad)
        id = id_gauge;
        id(ismissing(id) ...
            | strlength(strip(id)) == 0) = "";
        has_id = bad & strlength(strip(id)) > 0;
        s(has_id) = opt.fallback_prefix + strip(id(has_id));
        s(bad & ~has_id) = "Gauge";
    end

    if opt.title_case
        if opt.canadian_title
            s = local_title_case_canadian(s,opt.small_words);
        elseif opt.spanish_title
            s = local_title_case_spanish(s);
        else
            if opt.title_case_upper_only
                isUpper = (s == upper(s)) & (s ~= lower(s));
                s(isUpper) = local_title_case_simple( ...
                    s(isUpper),opt.small_words);
            else
                s = local_title_case_simple(s,opt.small_words);
            end
        end
    end

    if strcmp(char(reg),'US')
        s = local_standardize_us_display_names(s);
    end
    % CAMELS-FI display-name cleanup.
    % Preserve source capitalization generally, but standardize common
    % Finnish gauge descriptors after hyphens.
    if strcmp(char(reg),'FI')
        s = local_standardize_fi_display_names(s);
    end

    s = s(:);

end

% =====================================================
function s = local_standardize_us_display_names(s)
% =====================================================
% Expand established USGS station-name abbreviations and standardize the
% terminal postal state code. Ambiguous one-letter tokens are expanded only
% where the hydrologic context is clear.
    s = string(s(:));

    replacements = { ...
        '(?i)\<W\s+Fk\.?\>', 'West Fork'; ...
        '(?i)\<M\s+Fk\.?\>', 'Middle Fork'; ...
        '(?i)\<So\.?\s+Fk\.?\>', 'South Fork'; ...
        '(?i)\<No\.?\s+Fk\.?\>', 'North Fork'; ...
        '(?i)\<M\.?\s*F\.?\>', 'Middle Fork'; ...
        '(?i)\<W\.?\s*F\.?\>', 'West Fork'; ...
        '(?i)\<E\.?\s*F\.?\>', 'East Fork'; ...
        '(?i)\<N\.?\s*F\.?\>', 'North Fork'; ...
        '(?i)\<S\.?\s*F\.?\>', 'South Fork'; ...
        '(?i)\<W\s+Br\.?\>', 'West Branch'; ...
        '(?i)\<E\s+Br\.?\>', 'East Branch'; ...
        '(?i)\<N\s+Br\.?\>', 'North Branch'; ...
        '(?i)\<S\s+Br\.?\>', 'South Branch'; ...
        '(?i)\<W\s+B\.?\>',  'West Branch'; ...
        '(?i)\<E\s+B\.?\>',  'East Branch'; ...
        '(?i)\<N\s+B\.?\>',  'North Branch'; ...
        '(?i)\<S\s+B\.?\>',  'South Branch'; ...
        '(?i)\<So\.?\s+Br\.?\>', 'South Branch'; ...
        '(?i)\<E\s+Fk\.?\>', 'East Fork'; ...
        '(?i)\<N\s+Fk\.?\>', 'North Fork'; ...
        '(?i)\<S\s+Fk\.?\>', 'South Fork'; ...
        '(?i)\<SF\.?\>',     'South Fork'; ...
        '(?i)\<NF\.?\>',     'North Fork'; ...
        '(?i)\<EF\.?\>',     'East Fork'; ...
        '(?i)\<WF\.?\>',     'West Fork'; ...
        '(?i)\<MF\.?\>',     'Middle Fork'; ...
        '(?i)\<NB\.?\>',     'North Branch'; ...
        '(?i)\<SB\.?\>',     'South Branch'; ...
        '(?i)\<EB\.?\>',     'East Branch'; ...
        '(?i)\<WB\.?\>',     'West Branch'; ...
        '(?i)\<Fk\.?\>',     'Fork'; ...
        '(?i)\<Br\.?\>',     'Branch'; ...
        '(?i)\<Rv\.?\>',     'River'; ...
        '(?i)\<R\.?\>',      'River'; ...
        '(?i)\<Ck\.?\>',     'Creek'; ...
        '(?i)\<Crk\.?\>',    'Creek'; ...
        '(?i)\<Cr\.?\>',     'Creek'; ...
        '(?i)\<C\.?\>',      'Creek'; ...
        '(?i)\<Riv\.?\>',    'River'; ...
        '(?i)\<Bk\.?\>',     'Brook'; ...
        '(?i)\<Rn\.?\>',     'Run'; ...
        '(?i)\<Lk\.?\>',     'Lake'; ...
        '(?i)\<Res\.?\>',    'Reservoir'; ...
        '(?i)\<Nr\.?\>',     'near'; ...
        '(?i)\<Bl\.?\>',     'below'; ...
        '(?i)\<Blw\.?\>',    'below'; ...
        '(?i)\<Ab\.?\>',     'above'; ...
        '(?i)\<Abv\.?\>',    'above'; ...
        '(?i)\<Ds\.?\>',     'downstream'; ...
        '(?i)\<Us\.?\>',     'upstream'; ...
        '(?i)\<Hts\.?\>',    'Heights'; ...
        '(?i)\<Pk\.?\s+Dr\.?\>', 'Park Drive'; ...
        '(?i)\<Dr\.?\>',     'Drive'; ...
        '(?i)\<Rd\.?\>',     'Road'; ...
        '(?i)\<Ave\.?\>',    'Avenue'; ...
        '(?i)\<Blvd\.?\>',   'Boulevard'; ...
        '(?i)\<Hwy\.?\>',    'Highway'; ...
        '(?i)\<Hiway\.?\>',  'Highway'; ...
        '(?i)\<Highwy\.?\>', 'Highway'; ...
        '(?i)\<Cnty\.?\>',   'County'; ...
        '(?i)\<Co\.?\s+Road\>', 'County Road'; ...
        '(?i)\<Nat\.?\>',    'National'; ...
        '(?i)\<Hist\.?\>',   'Historic'; ...
        '(?i)\<Stn\.?\>',    'Station'; ...
        '(?i)\<Gen\.?\>',    'Generating'; ...
        '(?i)\<Cyn\.?\>',    'Canyon'; ...
        '(?i)\<RR\.?\>',     'Railroad'; ...
        '(?i)\<D\.?\>',      'Dam'; ...
        '(?i)\<Trib\.?\>',   'Tributary'; ...
        '(?i)\<Jct\.?\>',    'Junction'; ...
        '(?i)\<Mtn\.?\>',    'Mountain'; ...
        '(?i)\<Mt\.?(?=\s+\S)', 'Mount'; ...
        '(?i)\<Spgs\.?\>',   'Springs'; ...
        '(?i)\<Spr\.?\>',    'Spring'; ...
        '(?i)\<Reser\.?\>',  'Reservoir'; ...
        '(?i)\<Sta\.?\>',    'Station'; ...
        '(?i)\<Xing\.?\>',   'Crossing'; ...
        '(?i)\<Div\.?\>',    'Diversion'; ...
        '(?i)\<Confl\.?\>',  'confluence'; ...
        '(?i)\<Lwr\.?\>',    'Lower'; ...
        '(?i)\<Pkwy\.?\>',   'Parkway'; ...
        '(?i)\<Pky\.?\>',    'Parkway'; ...
        '(?i)\<Rte\.?\>',    'Route'; ...
        '(?i)\<Shwy\.?\>',   'State Highway'; ...
        '(?i)\<Bndy\.?\>',   'Boundary'; ...
        '(?i)\<Upst\.?\>',   'upstream'; ...
        '(?i)\<Pwrplnt\.?\>', 'Power Plant'; ...
        '(?i)\<Proj\.?\>',   'Project'; ...
        '(?i)\<Ctr\.?\>',    'Center'; ...
        '(?i)\<Rnch\.?\>',   'Ranch'; ...
        '(?i)\<Bch\.?\>',    'Beach'; ...
        '(?i)\<Ft\.?\>',     'Fort'};

    s = regexprep(s,'\s*@\s*',' at ');
    for j = 1:size(replacements,1)
        s = regexprep(s,replacements{j,1},replacements{j,2});
    end
    s = regexprep(s,'(?i)\<near\>','near');
    s = regexprep(s,'(?i)\<at\>','at');
    s = regexprep(s,'(?i)\<above\>','above');
    s = regexprep(s,'(?i)\<below\>','below');
    s = regexprep(s,'(?i)\<upstream\>','upstream');
    s = regexprep(s,'(?i)\<downstream\>','downstream');

    % Expand leading directional qualifiers used by USGS. At this stage,
    % compound forms such as "N F" and "S BR" have already become
    % "North Fork" and "South Branch" through the table above.
    s = regexprep(s,'(?i)^So\.?\s+','South ');
    s = regexprep(s,'(?i)^N\.?\s+','North ');
    s = regexprep(s,'(?i)^S\.?\s+','South ');
    s = regexprep(s,'(?i)^E\.?\s+','East ');
    s = regexprep(s,'(?i)^W\.?\s+','West ');
    s = regexprep(s,'(?i)^M\.?\s+','Middle ');

    % "L" is consistently "Little" in the selected CAMELSH-US names.
    % Requiring a following letter protects the lock-and-dam form "L & D".
    s = regexprep(s,'(?i)\<L\.?\>(?=\s+[A-Za-z])','Little');

    % In western USGS station names, the isolated token "A" is an
    % abbreviation for "at" when it follows a waterbody descriptor.
    % Keep route designations such as "Highway A" unchanged.
    s = regexprep(s,[ ...
        '(?i)\<(River|Creek|Brook|Branch|Fork|Canal|Wash|' ...
        'Reservoir|Lake|Bayou|Arroyo)\s+A\s+'], '$1 at ');
    s = regexprep(s,[ ...
        '(?i)\<(River|Creek|Brook|Branch|Fork|Canal|Wash|' ...
        'Reservoir|Lake|Bayou|Arroyo)\s+Bel\s+'], '$1 below ');

    % California (and a small number of Florida records) use "A" for
    % "at" even after a full place or facility name, for example
    % "Lake Arrowhead A Lake Arrowhead". Restrict this convention to
    % those terminal state codes and protect lettered/highway routes.
    atStyle = ~cellfun(@isempty,regexp(cellstr(s), ...
        '(?i),?\s*(CA|FL)\.?\s*$','once'));
    s(atStyle) = regexprep(s(atStyle), ...
        '(?i)(?<!Highway )(?<!-)\<A\>(?=\s+\S)','at');

    % Older USGS files use several non-postal state abbreviations.
    aliases = { ...
        'DISTRICT OF COLUMBIA','DC'; ...
        'ALABAMA','AL'; 'ALASKA','AK'; 'ARIZONA','AZ'; ...
        'ARKANSAS','AR'; 'CALIFORNIA','CA'; 'COLORADO','CO'; ...
        'CONNECTICUT','CT'; 'DELAWARE','DE'; 'FLORIDA','FL'; ...
        'GEORGIA','GA'; 'HAWAII','HI'; 'IDAHO','ID'; ...
        'ILLINOIS','IL'; 'INDIANA','IN'; 'IOWA','IA'; ...
        'KANSAS','KS'; 'KENTUCKY','KY'; 'LOUISIANA','LA'; ...
        'MAINE','ME'; 'MARYLAND','MD'; 'MASSACHUSETTS','MA'; ...
        'MICHIGAN','MI'; 'MINNESOTA','MN'; 'MISSISSIPPI','MS'; ...
        'MISSOURI','MO'; 'MONTANA','MT'; 'NEBRASKA','NE'; ...
        'NEVADA','NV'; 'NEW HAMPSHIRE','NH'; 'NEW JERSEY','NJ'; ...
        'NEW MEXICO','NM'; 'NEW YORK','NY'; ...
        'NORTH CAROLINA','NC'; 'NORTH DAKOTA','ND'; 'OHIO','OH'; ...
        'OKLAHOMA','OK'; 'OREGON','OR'; 'PENNSYLVANIA','PA'; ...
        'RHODE ISLAND','RI'; 'SOUTH CAROLINA','SC'; ...
        'SOUTH DAKOTA','SD'; 'TENNESSEE','TN'; 'TEXAS','TX'; ...
        'UTAH','UT'; 'VERMONT','VT'; 'VIRGINIA','VA'; ...
        'WASHINGTON','WA'; 'WEST VIRGINIA','WV'; ...
        'WISCONSIN','WI'; 'WYOMING','WY'; ...
        'ALA','AL'; 'ARIZ','AZ'; 'ARK','AR'; 'CALIF','CA'; ...
        'COLO','CO'; 'CONN','CT'; 'DEL','DE'; 'FLA','FL'; ...
        'ILL','IL'; 'IND','IN'; 'KANS','KS'; 'MASS','MA'; ...
        'MICH','MI'; 'MINN','MN'; 'MISS','MS'; 'MONT','MT'; ...
        'NEBR','NE'; 'NEV','NV'; 'N MEX','NM'; 'N Y','NY'; ...
        'N C','NC'; 'N DAK','ND'; 'OKLA','OK'; 'OREG','OR'; ...
        'PENN','PA'; 'S C','SC'; 'S DAK','SD'; 'TENN','TN'; ...
        'TEX','TX'; 'VT','VT'; 'W VA','WV'; 'WASH','WA'; ...
        'WIS','WI'; 'WYO','WY'};
    postal = ["AL","AK","AZ","AR","CA","CO","CT","DE","FL", ...
        "GA","HI","ID","IL","IN","IA","KS","KY","LA","ME", ...
        "MD","MA","MI","MN","MS","MO","MT","NE","NV","NH", ...
        "NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI", ...
        "SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY"];

    for k = 1:numel(s)
        name = strip(s(k));
        code = "";
        for j = 1:size(aliases,1)
            % Require a comma or whitespace before the state. Without this
            % boundary, full names such as Georgia and Maine were wrongly
            % interpreted from their final two letters as IA and NE.
            pat = ['(?i)(?:\s*,\s*|\s+)' ...
                regexptranslate('escape',aliases{j,1}) '\.?\s*$'];
            if ~isempty(regexp(char(name),pat,'once'))
                name = string(regexprep(char(name),pat,''));
                code = string(aliases{j,2});
                break
            end
        end
        if code == ""
            tok = regexp(char(name),'(?i)(?:\s*,\s*|\s+)([A-Z]{2})\.?\s*$', ...
                'tokens','once');
            if ~isempty(tok) && any(upper(string(tok{1})) == postal)
                code = upper(string(tok{1}));
                name = string(regexprep(char(name), ...
                    '(?i)(?:\s*,\s*|\s+)[A-Z]{2}\.?\s*$',''));
            end
        end
        name = strip(regexprep(name,'\s+',' '));
        if code ~= ""
            s(k) = name + ", " + code;
        else
            s(k) = name;
        end
    end
end
% ================================================
function s = local_title_case_simple(s,smallWords)
% ================================================
    s = string(s(:));

    if nargin < 2 || isempty(smallWords)
        smallWords = strings(0,1);
    end
    smallWords = lower(string(smallWords(:)));

    for k = 1:numel(s)
        if ismissing(s(k)) || strlength(s(k)) == 0
            continue
        end

        words = split(lower(s(k)));

        for j = 1:numel(words)
            w = string(words(j));

            if strlength(w) == 0
                continue
            end

            if j > 1 && any(w == smallWords)
                words(j) = w;
            else
                c = char(w);
                words(j) = string([upper(c(1)) c(2:end)]);
            end
        end

        s(k) = strjoin(words,' ');
    end
end

% ==================================================
function s = local_title_case_canadian(s,smallWords)
% ==================================================
% Title case for Canadian station names. Keeps common linking words in
% lower case (except at the start), preserves hyphenated names, and fixes
% common Mc-/O'-style names and Roman numerals.

    s = string(s(:));
    if nargin < 2 || isempty(smallWords)
        smallWords = strings(0,1);
    end
    smallWords = lower(string(smallWords(:)));

    for k = 1:numel(s)
        if ismissing(s(k)) || strlength(s(k)) == 0
            continue
        end

        words = split(lower(s(k)));

        for j = 1:numel(words)
            w = string(words(j));
            if strlength(w) == 0
                continue
            end

            if j > 1 && any(w == smallWords)
                words(j) = w;
                continue
            end

            parts = split(w,"-");
            for p = 1:numel(parts)
                parts(p) = local_ca_cap_word(parts(p));
            end
            words(j) = strjoin(parts,"-");
        end

        s(k) = strjoin(words,' ');
    end
end

% =========================================
function w = local_ca_cap_word(w)
% =========================================
    w = string(w);
    if strlength(w) == 0
        return
    end

    if any(upper(w) == ["I","II","III","IV","V","VI","VII","VIII","IX","X"])
        w = upper(w);
        return
    end

    c = char(lower(w));

    if numel(c) >= 3 && strcmp(c(1:2),'mc')
        c = ['Mc' upper(c(3)) c(4:end)];
        w = string(c);
        return
    end

    if numel(c) >= 3 && c(1) == 'o' && c(2) == ''''
        c = ['O''' upper(c(3)) c(4:end)];
        w = string(c);
        return
    end

    c(1) = upper(c(1));
    w = string(c);
end

% ======================================
function s = local_title_case_spanish(s)
% ======================================
    s = string(s(:));
    smallWords = ["a","al","de","del","el","en", ...
        "la","las","los","y"];

    for k = 1:numel(s)
        if ismissing(s(k)) ...
                || strlength(s(k)) == 0
            continue
        end

        words = split(lower(s(k)));
        for j = 1:numel(words)
            w = string(words(j));
            if strlength(w) == 0
                continue
            end

            if j > 1 && any(w == smallWords)
                words(j) = w;
            else
                c = char(w);
                words(j) = string([upper(c(1)) c(2:end)]);
            end
        end
        s(k) = strjoin(words,' ');
    end
end

% ================================================
function s = local_standardize_fi_display_names(s)
% ================================================
    s = string(s(:));
    s(ismissing(s)) = "";

    % Normalize spacing around separators:
    % "Pallasjarvi-luusua" -> "Pallasjarvi - luusua"
    s = regexprep(s,'\s*-\s*',' - ');

    % Normalize spacing after commas:
    % "Karhijarvenpato,ala" -> "Karhijarvenpato, ala"
    s = regexprep(s,'\s*,\s*',', ');

    % Normalize repeated spaces
    s = regexprep(s,'\s+',' ');
    s = strip(s);

    % Lowercase descriptive Finnish hydrologic/station qualifiers
    % when they occur after " - " or after a comma.
    lowerWords = ["luusua","virtuaali", ...
        "kokonaisvirtaama","pato","ala"];

    % Proper names that should be capitalized if source is inconsistent
    properWords = ["Kongas","Roukkajankoski","Jolma", ...
        "Lohikoski","Hulttilanjoki","Vaajakoski", ...
        "Kurittukoski","Hanhikoski"];

    for k = 1:numel(s)

        txt = char(s(k));

        % Descriptors: lowercase after dash or comma
        for j = 1:numel(lowerWords)
            w = char(lowerWords(j));

            txt = regexprep(txt, ...
                ['(?<= - )' w '(?=$|[ ,;/])'], ...
                w,'ignorecase');

            txt = regexprep(txt, ...
                ['(?<=, )' w '(?=$|[ ,;/])'], ...
                w,'ignorecase');
        end

        % Proper station/place names: capitalize after dash or at start
        for j = 1:numel(properWords)
            w = char(properWords(j));

            txt = regexprep(txt, ...
                ['^' lower(w) '(?=$|[ ,;/])'], ...
                w,'ignorecase');

            txt = regexprep(txt, ...
                ['(?<= - )' lower(w) '(?=$|[ ,;/])'], ...
                w,'ignorecase');
        end

        s(k) = string(txt);
    end
end
