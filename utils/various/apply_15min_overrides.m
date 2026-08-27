function [th_min,th_max,par_units,add_info] = ...
    apply_15min_overrides( ...
    model,par_names,th_min,th_max,par_units,add_info)
%APPLY_15MIN_OVERRIDES
% Optional explicit 15-minute overrides
%
% Use this only where the literature/experience suggests that the
% 15-minute ranges should differ from simple daily->15-minute scaling.

switch model

    case 4  % XINANJIANG
        % setpar('k_{\rm i}',    1e-5,     2.0,   '1/15 min');
        % setpar('k_{\rm g}',    1e-5,     2.5,   '1/15 min');
        % setpar('c_{\rm i}',    1e-4,     2.5,   '1/15 min');
        % setpar('c_{\rm g}',    1e-6,     2.0,   '1/15 min');
        % setpar('k_{\rm f}',    0.01,    15.0,   '1/15 min');

    case 6  % HBV
        % setpar('k_{0}',        1e-4,     0.5,   '1/15 min');
        % setpar('k_{1}',        1e-4,     0.1,   '1/15 min');
        % setpar('k_{2}',        1e-6,     0.1,   '1/15 min');
        % setpar('{\rm perc}',      0,      10,   'mm/15 min');
        % setpar('cf_{\rm max}',    0,     0.5,   'mm/15 min/°C');
        % setpar('f_{\rm r}',       0,    0.25,   '-');
        % setpar('b_{\rm rt}',     10,     600,   '15 min');

    case 7  % GR4JB
        % Explicit 15-minute routing override from Mathias et al. logic
        setpar('x_{4}', 0.01/96, 2/96, '15 min');

    otherwise
        % no explicit 15-minute overrides
end

    function setpar(pname,lb,ub,unit)
        ii = find(strcmp(par_names,pname),1,'first');
        if ~isempty(ii)
            th_min(ii) = lb;
            th_max(ii) = ub;
            par_units{ii} = unit;
        end
    end

end