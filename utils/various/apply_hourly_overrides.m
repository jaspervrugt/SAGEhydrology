function [th_min,th_max,par_units,add_info] = ...
    apply_hourly_overrides( ...
    model,par_names,th_min,th_max,par_units,add_info)
%APPLY_HOURLY_OVERRIDES 
% Optional explicit hourly overrides
%
% Use this only where the literature/experience suggests that the hourly
% ranges should differ from simple daily->hourly scaling.

switch model

    case 4  % XINANJIANG
        % setpar('k_{\rm i}',    1e-5,      2.0,   '1/h');
        % setpar('k_{\rm g}',    1e-5,      2.5,   '1/h');
        % setpar('c_{\rm i}',    1e-4,      2.5,   '1/h');
        % setpar('c_{\rm g}',    1e-6,      2.0,   '1/h');
        % setpar('k_{\rm f}',    0.01,     15.0,   '1/h');

    case 6  % HBV
        % setpar('k_{0}',        1e-4,      0.5,   '1/h');
        % setpar('k_{1}',        1e-4,      0.1,   '1/h');
        % setpar('k_{2}',        1e-6,      0.1,   '1/h');
        % setpar('{\rm perc}',      0,       10,   'mm/h');
        % setpar('cf_{\rm max}',    0,      0.5,   'mm/h/°C');
        % setpar('f_{\rm r}',       0,     0.25,   '-');
        % setpar('b_{\rm rt}',     10,      600,   'h');

    case 7  % GR4JB
        % Explicit hourly routing override from Mathias et al. logic
        setpar('x_{4}', 0.01/24, 2/24, 'h');

    otherwise
        % no explicit hourly overrides
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