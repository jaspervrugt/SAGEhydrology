function method = observation_transform(loss,name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBSERVATION_TRANSFORM Validate and return a response objective option.
%  loss.transform.SM = 'zscore' selects population z-score MSE = 2(1-r).
%  Unspecified responses use 'none' and retain loss.fnc unchanged.
%
% Copyright (c) 2026 Jasper A. Vrugt, University of California, Irvine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    method = 'none';
    if ~isfield(loss,'transform') || ...
        isempty(loss.transform)
        return
    end
    if ~isstruct(loss.transform) || ...
        ~isscalar(loss.transform)
        error('SAGE:Transform','loss.transform must be a named structure.');
    end
    fields = fieldnames(loss.transform);
    for k = 1:numel(fields)
        key = fields{k};
        value = string(loss.transform.(key));
        if ~ismember(key,{'Q','SWE','SM'}) || ...
            ~isscalar(value) || ...
            ~ismember(value,["none","zscore"])
            error('SAGE:Transform', ...
                'Use none or zscore for loss.transform.Q/SWE/SM.');
        end
    end
    key = char(name);
    if isfield(loss.transform,key)
        method = char(loss.transform.(key));
    end
end
