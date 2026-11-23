function [q,J] = cart_model_discharge(p,slope,mdl,dat,train,needJ,ode)
%CART_MODEL_DISCHARGE Adapt the shared basin data to CaRT model solvers.

    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'models','cart'), ...
        fullfile(root,'models','cartplus'));
    meteo = dat.meteo;
    mask = false(numel(meteo.P),1);
    mask(train) = true;
    data = dat.cart_bands;
    data.forcing = [double(meteo.P(:)),double(meteo.T(:)), ...
        double(meteo.Ep(:)),double(mask)];
    data.resolution = mdl.dt;
    J = [];
    useOde = isfield(mdl,'config') && ...
        isfield(mdl.config,'cart_formulation') && ...
        strcmp(mdl.config.cart_formulation,'ode');
    if useOde
        names = intersect(fieldnames(mdl.config), ...
            {'cart_solver','cart_refinement','cart_rtol','cart_atol', ...
            'cart_max_step','cart_min_step'},'stable');
        args = cell(1,2*numel(names));
        for k = 1:numel(names)
            args(2*k-1:2*k) = {names{k},mdl.config.(names{k})};
        end
    end
    if mdl.model == 19
        if useOde
            if needJ
                [q,J] = cart_ode_discharge(data,p,19,args{:});
            else
                q = cart_ode_discharge(data,p,19,args{:});
            end
        elseif needJ
            [q,J] = cart(p,mdl,data);
        else
            q = cart(p,mdl,data);
        end
    else
        par = mdl.cart_fixed;
        j = mdl.cart_columns;
        logmask = mdl.cart_multiplier ~= 0;
        scale = ones(size(p));
        scale(logmask) = mdl.cart_multiplier(logmask).*exp(p(logmask));
        par(j) = p;
        par(j(logmask)) = scale(logmask);
        if useOde
            if needJ
                [out,J] = cart_ode_discharge(data,par,18,args{:});
                J = J(:,j).*scale.';
            else
                out = cart_ode_discharge(data,par,18,args{:});
            end
        elseif needJ
            [out,J] = cartplus(par,mdl,data,ode);
            J = J(:,j).*scale.';
        else
            out = cartplus(par,mdl,data,ode);
        end
        q = out(:,1);
    end
    if needJ
        J = J.*slope.';
    end
end
