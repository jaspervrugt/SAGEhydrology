function request = crr_request(request)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRR_REQUEST Validate and normalize a CRR_MODEL output request.
%
%  Fill omitted request fields with defaults and reject unknown fields or
%  observable names.
%
% SYNOPSIS:
%   request = crr_request(request)
%   request = crr_request()
%
% INPUT ARGUMENTS:
%   request         optional scalar structure of requested model outputs
%    .q              return simulated discharge
%    .gradient       return the selected-loss gradient
%    .jacobian       return the discharge Jacobian
%    .metrics        return compact performance metrics
%    .attribution    return parameter attribution
%    .states         false/empty, true/"all", or requested state names
%    .obs            named trajectories: Q, SWE, or SM
%    .jac            named Jacobians: Q, SWE, or SM
%
% OUTPUT ARGUMENTS:
%   request         normalized structure with every supported field
%    .normalized     true; marks the validated request contract
%
% NOTES:
%   Omitted logical fields are false; omitted name lists are empty.
%   Unknown fields and observable names cause an error.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaults = struct('q',false,'gradient',false,'jacobian',false, ...
        'metrics',false,'attribution',false,'states',strings(0,1), ...
        'obs',strings(0,1),'jac',strings(0,1), ...
        'normalized',true);
    if nargin == 0 || isempty(request), request = defaults; return, end
    if ~isstruct(request) || ~isscalar(request)
        error('crr_request:InvalidRequest','Request must be a scalar structure.');
    end
    known = fieldnames(defaults);
    unknown = setdiff(fieldnames(request),known);
    if ~isempty(unknown)
        error('crr_request:UnknownField','Unknown request field(s): %s.', ...
            strjoin(unknown,', '));
    end
    flags = {'q','gradient','jacobian','metrics','attribution'};
    for k = 1:numel(flags)
        f = flags{k};
        if ~isfield(request,f) || isempty(request.(f))
            request.(f) = false;
        elseif ~(islogical(request.(f)) && isscalar(request.(f)))
            error('crr_request:InvalidFlag', ...
                'request.%s must be a scalar logical.',f);
        end
    end
    if ~isfield(request,'states') || isempty(request.states) ...
            || (islogical(request.states) && isscalar(request.states) ...
            && ~request.states)
        request.states = strings(0,1);
    elseif islogical(request.states) && isscalar(request.states)
        request.states = "all";
    else
        if ischar(request.states)
            request.states = string({request.states});
        else
            request.states = string(request.states(:));
        end
        request.states = request.states(strlength(request.states) > 0);
    end
    listFields = {'obs','jac'};
    for k = 1:numel(listFields)
        f = listFields{k};
        if ~isfield(request,f) || isempty(request.(f)) ...
                || (islogical(request.(f)) && isscalar(request.(f)) ...
                && ~request.(f))
            request.(f) = strings(0,1);
        elseif islogical(request.(f))
            error('crr_request:InvalidList', ...
                'request.%s must contain observable names.',f);
        else
            if ischar(request.(f))
                names = string({request.(f)});
            else
                names = string(request.(f)(:));
            end
            names = upper(strtrim(names));
            names = names(strlength(names) > 0);
            allowed = ["Q";"SWE";"SM"];
            unknownNames = setdiff(names,allowed);
            if ~isempty(unknownNames)
                error('crr_request:UnknownObservable', ...
                    'Unknown request.%s name(s): %s.',f, ...
                    strjoin(cellstr(unknownNames),', '));
            end
            request.(f) = unique(names,'stable');
        end
    end
    request.normalized = true;
end
