function request = crr_request(request)
%CRR_REQUEST Normalize output requests for CRR_MODEL.
    defaults = struct('q',false,'gradient',false,'jacobian',false, ...
        'metrics',false,'attribution',false,'states',strings(0,1), ...
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
        request.states = string(request.states(:));
        request.states = request.states(strlength(request.states) > 0);
    end
    request.normalized = true;
end