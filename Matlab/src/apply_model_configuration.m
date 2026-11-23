function mdl = apply_model_configuration(mdl,ode)
% APPLY_MODEL_CONFIGURATION Resolve initial states and model-specific
% constants.
% ode.y0 overrides the legacy solver vector mdl.y0. Missing overrides retain
% the source defaults, including the cumulative-discharge bookkeeping state.
    if ismember(mdl.model,[18 19])
        % Both models compute initial stores from the forcing.
        if ~isfield(mdl,'config'), mdl.config=struct(); end
        mdl.config = cart_numerical_configuration(mdl.config);
        if strcmp(mdl.config.cart_formulation,'sequential') && ...
                mdl.config.cart_uniform_input
            mdl.th_min(2) = 1/24;
            mdl.th_min_daily(2) = 1/24;
        end
        return
    end
    if nargin<2, ode=struct(); end
    if mdl.model==7
        config=resolve_model_config(mdl);
        oldL=numel(mdl.giuh_ordnts); oldK=numel(mdl.y0)-4-oldL;
        newL=numel(config.giuh_ordnts); newK=config.K;
        if oldL~=newL || ...
            oldK~=newK
            old=mdl.y0(:); resized=1e-5*ones(4+newL+newK,1);
            resized(1:3)=old(1:3);
            n=min(oldL,newL); resized(3+(1:n))=old(3+(1:n));
            n=min(oldK,newK); resized(3+newL+(1:n))=old(3+oldL+(1:n));
            resized(end)=old(end); mdl.y0=resized;
            mdl.initial_state_defaults=1e-5*ones(size(resized));
            if isfield(ode,'y0') && ...
                numel(ode.y0)==numel(old)
                initial=ode.y0(:); resized(1:3)=initial(1:3);
                n=min(oldL,newL); resized(3+(1:n))=initial(3+(1:n));
                n=min(oldK,newK); resized(3+newL+(1:n))=initial(3+oldL+( ...
                    1:n));
                resized(end)=initial(end); ode.y0=resized;
            end
        end
        mdl.giuh_ordnts=config.giuh_ordnts;
        mdl.config=config;
    end
    if ~isfield(mdl,'initial_state_defaults')
        mdl.initial_state_defaults=mdl.y0(:);
    end
    if isfield(ode,'y0') && ...
        ~isempty(ode.y0)
        initial=ode.y0;
        if ~isnumeric(initial) || ...
            ~isreal(initial) || ...
            ~isvector(initial) || ...
            numel(initial)~=numel(mdl.y0) || ...
            any(~isfinite( initial)) || ...
            any(initial<0)
            error('apply_model_configuration:InitialStates', ...
                ['Initial states must contain %d finite ' ...
                    'nonnegative values.'],numel(mdl.y0));
        end
        physical=initial(1:end-1);
        if any(physical<=0 | ...
            physical>100)
            error('apply_model_configuration:InitialStorageRange', ...
                ['Editable initial storage values must be greater ' ...
                    'than zero and at most 100 mm.']);
        end
        mdl.y0=double(initial(:));
    end
    if ismember(mdl.model,1:7) || ...
        isfield(mdl,'configuration_definition')
        mdl.config=resolve_model_config(mdl);
    end
end
