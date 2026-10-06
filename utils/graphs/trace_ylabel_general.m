function s = trace_ylabel_general(mdl,j,pspace)
%LOCAL_TRACE_YLABEL_GENERAL
% pspace = 1: underline main symbol only
% pspace = 0: no underline

    if nargin < 3
        pspace = 0;
    end

    if ischar(j) ...
            || isstring(j)
        name = strtrim(char(j));
    else
        name = local_get_par_name(mdl,j);
    end
    
    [mainSym,subSym] = local_split_latex_symbol(name);
    
    if pspace == 1
        if isempty(subSym)
            s = sprintf('$\\underline{%s}$',mainSym);
        else
            s = sprintf('$\\underline{%s}_{%s}$',mainSym,subSym);
        end
    else
        if isempty(subSym)
            s = sprintf('$%s$',mainSym);
        else
            s = sprintf('$%s_{%s}$',mainSym,subSym);
        end
    end
end

function name = local_get_par_name(mdl,j)
%LOCAL_GET_PAR_NAME Return raw parameter name

    name = '';
    
    try
        if isfield(mdl,'par_names') ...
                && numel(mdl.par_names) >= j
            p = mdl.par_names{j};
            if isstring(p)
                p = char(p);
            end
            name = strtrim(char(p));
        end
    catch
    end

    if isempty(name)
        try
            if isfield(mdl,'par') ...
                    && isfield(mdl.par,'name') ...
                    && numel(mdl.par.name) >= j
                p = mdl.par.name{j};
                if isstring(p)
                    p = char(p);
                end
                name = strtrim(char(p));
            end
        catch
        end
    end
    
    if isempty(name)
        name = sprintf('theta_%d',j);
    end
end

function [mainSym,subSym] = local_split_latex_symbol(name)
%LOCAL_SPLIT_LATEX_SYMBOL Split parameter name into main symbol + subscript

    name = strtrim(char(name));
    subSym = '';

    % GCHM uses descriptive implementation names that are much too long to
    % typeset literally.  Keep this map here (rather than in an individual
    % GUI) so the on-screen traces and exported SAGE figures use exactly the
    % same compact, valid LaTeX notation.
    [isGchm,mainSym,subSym] = local_gchm_symbol(name);
    if isGchm
        return
    end
    
    % Case 1: something_{...}
    tok = regexp(name,'^(.+)_\{(.+)\}$','tokens','once');
    if ~isempty(tok)
        mainSym = local_format_symbol(tok{1});
        subSym = local_format_subsymbol(tok{2});
        return
    end
    
    % Case 2: something_x
    tok = regexp(name,'^([^_]+)_(.+)$','tokens','once');
    if ~isempty(tok)
        mainSym = local_format_symbol(tok{1});
        subSym = local_format_subsymbol(tok{2});
        return
    end
    
    % Case 3: no subscript
    mainSym = local_format_symbol(name);
end

function [matched,mainSym,subSym] = local_gchm_symbol(name)
%LOCAL_GCHM_SYMBOL Compact mathematical labels for all 30 GCHM parameters.

    matched = true;
    mainSym = '';
    subSym = '';

    switch lower(strtrim(char(name)))
        case 's_p_max_theta'
            mainSym = 'S'; subSym = '\mathrm{P,max},\theta';
        case 's_p_max_b_parameter_theta'
            mainSym = 'S'; subSym = '\mathrm{P,max},b\theta';
        case 's_b1_theta'
            mainSym = 'S'; subSym = '\mathrm{b1},\theta';
        case 's_b2_theta'
            mainSym = 'S'; subSym = '\mathrm{b2},\theta';
        case 's_b3_theta'
            mainSym = 'S'; subSym = '\mathrm{b3},\theta';
        case 's_a1_theta'
            mainSym = 'S'; subSym = '\mathrm{a1},\theta';
        case 's_max_theta'
            mainSym = 'S'; subSym = '\mathrm{max},\theta';
        case 's_max_b_parameter_theta'
            mainSym = 'S'; subSym = '\mathrm{max},b\theta';
        case 's_b_seepage_x'
            mainSym = 'S'; subSym = '\mathrm{b,seepage},X';
        case 's_a_seepage'
            mainSym = 'S'; subSym = '\mathrm{a,seepage}';
        case 's_c_seepage'
            mainSym = 'S'; subSym = '\mathrm{c,seepage}';
        case 's_b_out_x'
            mainSym = 'S'; subSym = '\mathrm{b,out},X';
        case 's_a_out'
            mainSym = 'S'; subSym = '\mathrm{a,out}';
        case 's_c_out'
            mainSym = 'S'; subSym = '\mathrm{c,out}';
        case 's_c_et'
            mainSym = 'S'; subSym = '\mathrm{c,ET}';
        case 's_c_remember'
            mainSym = 'S'; subSym = '\mathrm{c,remember}';
        case 'b_b1_theta'
            mainSym = 'B'; subSym = '\mathrm{b1},\theta';
        case 'b_b2_theta'
            mainSym = 'B'; subSym = '\mathrm{b2},\theta';
        case 'b_b3_theta'
            mainSym = 'B'; subSym = '\mathrm{b3},\theta';
        case 'b_a1_theta'
            mainSym = 'B'; subSym = '\mathrm{a1},\theta';
        case 'b_b_out_x'
            mainSym = 'B'; subSym = '\mathrm{b,out},X';
        case 'b_a_out'
            mainSym = 'B'; subSym = '\mathrm{a,out}';
        case 'b_c_out'
            mainSym = 'B'; subSym = '\mathrm{c,out}';
        case 'b_c_et'
            mainSym = 'B'; subSym = '\mathrm{c,ET}';
        case 'b_c_remember'
            mainSym = 'B'; subSym = '\mathrm{c,remember}';
        case 'b1_melting_soil_theta'
            mainSym = 'b'; subSym = '1,\mathrm{melt},\theta';
        case 'b2_melting_soil_theta'
            mainSym = 'b'; subSym = '2,\mathrm{melt},\theta';
        case 'c_melting_soil_theta'
            mainSym = 'c'; subSym = '\mathrm{melt},\theta';
        case 'b_rain_fraction_theta'
            mainSym = 'b'; subSym = '\mathrm{rain},\theta';
        case 'c_rain_fraction_theta'
            mainSym = 'c'; subSym = '\mathrm{rain},\theta';
        otherwise
            matched = false;
    end
end

function s = local_format_symbol(txt)
%LOCAL_FORMAT_SYMBOL Format main symbol for LaTeX

    txt = strtrim(char(txt));
    
    if isempty(txt)
        s = '';
        return
    end
    
    % already LaTeX
    if startsWith(txt,'\')
        s = txt;
        return
    end
    
    % Greek shortcuts
    switch lower(txt)
        case 'alpha'
            s = '\alpha'; return
        case 'beta'
            s = '\beta';  return
        case 'gamma'
            s = '\gamma'; return
        case 'delta'
            s = '\delta'; return
        case 'epsilon'
            s = '\epsilon'; return
        case 'eta'
            s = '\eta'; return
        case 'theta'
            s = '\theta'; return
        case 'lambda'
            s = '\lambda'; return
        case 'mu'
            s = '\mu'; return
        case 'phi'
            s = '\phi'; return
        case 'psi'
            s = '\psi'; return
        case 'sigma'
            s = '\sigma'; return
        case 'omega'
            s = '\omega'; return
        case 'tau'
            s = '\tau'; return
    end
    
    % single latin letter -> italic math
    if isscalar(txt) && isletter(txt)
        s = txt;
        return
    end
    
    % compact names like k1 -> k_{1}
    tok = regexp(txt,'^([A-Za-z])([0-9]+)$','tokens','once');
    if ~isempty(tok)
        s = sprintf('%s_{%s}',tok{1},tok{2});
        return
    end
    
    % fallback upright roman
    txt = strrep(txt,'_','\_');
    txt = regexprep(txt,'\s+','\\,');
    s = sprintf('\\mathrm{%s}',txt);
end

function s = local_format_subsymbol(txt)
%LOCAL_FORMAT_SUBSYMBOL Format subscript for LaTeX

    txt = strtrim(char(txt));
    
    if isempty(txt)
        s = '';
        return
    end
    
    % already LaTeX-like
    if contains(txt,'\') ...
            || contains(txt,'{') ...
            || contains(txt,'}')
        s = txt;
        return
    end
    
    % purely numeric
    if all(isstrprop(txt,'digit'))
        s = txt;
        return
    end
    
    % single letter -> upright
    if isscalar(txt) && isletter(txt)
        s = sprintf('\\mathrm{%s}',txt);
        return
    end
    
    % compact names like k1 -> k_{1}
    tok = regexp(txt,'^([A-Za-z])([0-9]+)$', ...
        'tokens','once');
    if ~isempty(tok)
        s = sprintf('%s_{%s}',tok{1},tok{2});
        return
    end
    
    % fallback upright roman
    txt = strrep(txt,'_','\_');
    txt = regexprep(txt,'\s+','\\,');
    s = sprintf('\\mathrm{%s}',txt);
end
