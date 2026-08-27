function Q = derive_discharge_schema(M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DERIVE_DISCHARGE_SCHEMA Derive a shared-file discharge schema.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isfield(M,'profiles') ...
            || isempty(M.profiles)
        Q = local_projection(M);
        return
    end

    base = rmfield(M,'profiles');
    names = fieldnames(M.profiles);
    candidates = struct('dt',{},'schema',{});
    for i = 1:numel(names)
        profile = M.profiles.(names{i});
        D = local_merge(base,profile.schema);
        if ~isfield(D,'variables') ...
                || ~isfield(D.variables,'Q')
            continue
        end
        if isfield(profile.match,'dt') ...
                && ~isempty(profile.match.dt)
            values = double(profile.match.dt(:));
        else
            values = local_dt_from_step(D.timeline.step);
        end
        for j = 1:numel(values)
            if isempty(candidates) ...
                    || ~any([candidates.dt] == values(j))
                candidates(end+1).dt = values(j); %#ok<AGROW>
                candidates(end).schema = local_projection(D);
            end
        end
    end
    if isempty(candidates)
        error('derive_discharge_schema:MissingQ', ...
            'The meteorological schema does not declare Q.');
    end
    if isscalar(candidates)
        Q = candidates(1).schema;
        return
    end

    Q = base;
    Q.mode = 'Q';
    Q.variables = struct();
    if isfield(Q,'aux')
        Q = rmfield(Q,'aux');
    end
    for i = 1:numel(candidates)
        name = matlab.lang.makeValidName( ...
            sprintf('q_dt_%g',candidates(i).dt));
        Q.profiles.(name).match = struct('dt',candidates(i).dt);
        Q.profiles.(name).schema = candidates(i).schema;
    end
end

function Q = local_projection(S)
    if ~isfield(S,'variables') ...
            || ~isfield(S.variables,'Q')
        error('derive_discharge_schema:MissingQ', ...
            'The meteorological schema does not declare Q.');
    end
    Q = S;
    if isfield(Q,'profiles')
        Q = rmfield(Q,'profiles');
    end
    Q.variables = struct('Q',S.variables.Q);
    Q.mode = 'Q';
    Q.reuse_existing_q = true;
    if isfield(Q,'aux')
        Q = rmfield(Q,'aux');
    end
    Q.name = [char(string(S.name)) ' discharge'];
end

function dt = local_dt_from_step(step)
    stepHours = hours(step);
    if abs(stepHours-24) < 1e-9
        dt = 1;
    elseif abs(stepHours-1) < 1e-9
        dt = 24;
    elseif abs(stepHours-0.25) < 1e-9
        dt = 96;
    else
        error('derive_discharge_schema:UnknownResolution', ...
            'Cannot infer SAGE dt from a %g-hour timeline step.',stepHours);
    end
end

function out = local_merge(base,override)
    out = base;
    fields = fieldnames(override);
    for i = 1:numel(fields)
        field = fields{i};
        if isfield(out,field) ...
                && isstruct(out.(field)) ...
                && isscalar(out.(field)) ...
                && isstruct(override.(field)) ...
                && isscalar(override.(field))
            out.(field) = local_merge(out.(field),override.(field));
        else
            out.(field) = override.(field);
        end
    end
end
