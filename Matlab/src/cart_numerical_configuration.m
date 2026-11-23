function cfg = cart_numerical_configuration(cfg)
%CART_NUMERICAL_CONFIGURATION Validate and identify numerical choices.
    if nargin == 0, cfg = struct(); end
    if ~isfield(cfg,'cart_formulation')
        cfg.cart_formulation = 'sequential';
    end
    if ~ismember(string(cfg.cart_formulation),["sequential","ode"])
        error('cart:formulation','Unknown CaRT formulation.');
    end
    if strcmp(cfg.cart_formulation,'ode')
        if ~isfield(cfg,'cart_solver'), cfg.cart_solver='dop853'; end
        methods = ["dop853","rk2","rk23","rk45", ...
            "heun_fixed","euler"];
        if ~ismember(string(cfg.cart_solver),methods)
            error('cart:solver','Unknown CaRT ODE method.');
        end
        [rtol,atol] = cart_ode_tolerances(cfg.cart_solver);
        defaults = struct('cart_refinement',16,'cart_rtol',rtol, ...
            'cart_atol',atol,'cart_min_step',0);
        for name = string(fieldnames(defaults)).'
            if ~isfield(cfg,name), cfg.(name)=defaults.(name); end
        end
        numbers = [cfg.cart_refinement,cfg.cart_rtol,cfg.cart_atol];
        if numel(numbers) ~= 3 || any(~isfinite(numbers)) || ...
                any(numbers <= 0) || ...
                cfg.cart_refinement ~= fix(cfg.cart_refinement)
            error('cart:settings','Invalid ODE stages or tolerances.');
        end
        fixed = ismember(string(cfg.cart_solver), ...
            ["heun_fixed","euler","rk45"]);
        if ~isfield(cfg,'cart_max_step') || isempty(cfg.cart_max_step)
            cfg.cart_max_step = 1e100;
            if fixed, cfg.cart_max_step = 1/96; end
        end
        limits = [cfg.cart_min_step,cfg.cart_max_step];
        if numel(limits) ~= 2 || any(~isfinite(limits)) || ...
                limits(1) < 0 || limits(2) <= 0 || limits(1) > limits(2)
            error('cart:steps','Invalid step limits in days.');
        end
        cfg.cart_implementation_revision = 'cart-numerics-v1';
    else
        if ~isfield(cfg,'cart_uniform_input')
            cfg.cart_uniform_input = false;
        end
        if ~isscalar(cfg.cart_uniform_input) || ...
                ~ismember(cfg.cart_uniform_input,[0 1])
            error('cart:scheme','Scheme flag must be logical.');
        end
        cfg.cart_upstream_revision = ...
            '83682732728646669a190835d95ab3e60a8b5cb5';
        if cfg.cart_uniform_input
            cfg.cart_upstream_revision = ...
                'f986af265633a7a980836c2e0f2124077d632420';
        end
    end
end
