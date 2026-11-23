function record = model_run_configuration(mdl)
% MODEL_RUN_CONFIGURATION Actual settings, source defaults and filename
% identity.
% Editable is not evidence of output independence. Only explicit
% AffectsOutput=false metadata excludes a value from the identity.
base=mdl;
if isfield(mdl,'initial_state_defaults')
    base.y0=mdl.initial_state_defaults;
elseif ismember(mdl.model,1:7)
    base.y0=1e-5*ones(size(mdl.y0));
end
meta=model_configuration(base);
if ismember(mdl.model,[18 19])
    config=struct();
else
    config=resolve_model_config(mdl);
end
rows={'Section','Symbol','Description','Unit','Default','Used','Changed'};
identity=''; influentialChange=false;
for j=1:height(meta.states)
    def=meta.states.Default(j); used=mdl.y0(j); changed=~isequaln(def,used);
    rows(end+1,:)={'Initial state',meta.states.Symbol{j}, ...
        meta.states.Description{j}, ...
        meta.states.Unit{j},def,used,changed};
    identity=[identity sprintf('state:%s=%.17g;',meta.states.Symbol{j}, ...
        used)];                         %#ok<AGROW>
    influentialChange=influentialChange|| ...
        changed;
end
[~,order]=sort(string(meta.values.Symbol));
for j=order(:)'
    key=meta.values.Symbol{j}; def=meta.values.Default(j); used=config.( ...
        key);
    if strcmp(key,'giuh_ordnts'), def=[0.1 0.35 0.2 0.14 0.1 0.06 0.05]; end
    changed=~isequaln(def,used);
    rows(end+1,:)={'Additional value',key,meta.values.Name{j}, ...
        meta.values.Unit{j},def,used,changed};
    affects=true;
    if ismember('AffectsOutput',meta.values.Properties.VariableNames)
        affects=meta.values.AffectsOutput(j);
    end
    if affects
        if isscalar(used)
            identity=[identity sprintf('value:%s=%.17g;',key,used)];
        else
            identity=[identity sprintf('value:%s=%s;',key,mat2str(used, ...
                17))];
        end                             %#ok<AGROW>
        influentialChange=influentialChange|| ...
            changed;
    end
end
if mdl.model==99 && ...
    isfield(mdl,'user_model_folder')
    identity=[identity 'definition:' sage_user_model_fingerprint( ...
        mdl.user_model_folder) ';'];
end
if ismember(mdl.model,[18 19])
    used = struct();
    if isfield(mdl,'config'), used=mdl.config; end
    used = cart_numerical_configuration(used);
    keys = {'cart_formulation','cart_uniform_input', ...
        'cart_upstream_revision'};
    if strcmp(used.cart_formulation,'ode')
        keys = {'cart_formulation','cart_solver','cart_refinement', ...
            'cart_rtol','cart_atol','cart_max_step','cart_min_step', ...
            'cart_implementation_revision'};
    end
    keys = sort(keys);
    for j = 1:numel(keys)
        key = keys{j};
        value = used.(key);
        token = char(string(value));
        if isnumeric(value) || islogical(value)
            token=sprintf('%.17g',value);
        end
        identity = [identity sprintf('numerics:%s=%s;',key,token)];
        rows(end+1,:) = {'Numerics',key,key,'-',[],value,false};
    end
    influentialChange = true;
end
hash=0;
for c=double(unicode2native(identity,'UTF-8'))
    hash=mod(hash*65599+c,2^32);
end
record.suffix='';
if influentialChange, record.suffix=['_cfg-' lower(dec2hex(hash,8))]; end
record.rows=rows;
record.identity=identity;
record.changed=influentialChange;
end
