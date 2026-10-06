function available = regional_observation_inventory(region,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REGIONAL_OBSERVATION_INVENTORY Report observations declared by a dataset.
%
% Availability is determined entirely from the installed SAGE regional
% schema. It does not inspect or read downloaded basin data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    available = struct('Q',true,'SWE',false,'SM',false);
    schema = meteo_schema(region);
    available = local_add_variables(available,schema);

    if ~isfield(schema,'profiles') ...
            || ~isstruct(schema.profiles)
        return
    end
    names = fieldnames(schema.profiles);
    for kk = 1:numel(names)
        profile = schema.profiles.(names{kk});
        if ~local_resolution_matches(profile,dt)
            continue
        end
        if isfield(profile,'schema') ...
                && isstruct(profile.schema)
            available = local_add_variables(available,profile.schema);
        end
    end
end

function available = local_add_variables(available,schema)
    if ~isfield(schema,'variables') ...
            || ~isstruct(schema.variables)
        return
    end
    available.Q = available.Q || isfield(schema.variables,'Q');
    available.SWE = available.SWE || isfield(schema.variables,'SWE');
    available.SM = available.SM || isfield(schema.variables,'SM');
end

function tf = local_resolution_matches(profile,dt)
    tf = true;
    if ~isfield(profile,'match') ...
            || ~isstruct(profile.match) ...
            || ~isfield(profile.match,'dt')
        return
    end
    expected = double(profile.match.dt(:));
    tf = any(expected == double(dt));
end
