function dat = read_discharge_data(dirQ,mdl,dat,bas,split,aux,schema) %#ok<INUSD>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_DISCHARGE_DATA Read observed discharge using a declarative schema.
%
% SYNOPSIS:
%   dat = read_discharge_data(dirQ,mdl,dat,bas,split,aux,schema)
%
% INPUT:
%   dirQ       Directory containing discharge time-series files
%   dat        Existing meteorological SAGE data
%   bas        Basin structure
%   split      Prepared SAGE time-split structure
%   aux        Existing K-by-3 metadata; area is in column 3 [m2]
%   schema     Generic discharge schema containing variables.Q
%
% OUTPUT:
%   dat        Input data augmented with y_n, bad, and discharge filename
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 7
        error('read_discharge_data:MissingSchema', ...
            'A discharge schema must be provided.');
    end
    if nargin >= 6 ...
            && ~isempty(aux)
        schema.aux.values = aux;
    end
    schema = resolve_hydro_schema(schema,split,struct());
    schema = validate_hydro_schema(schema,'Q');
    if isfield(schema,'reuse_existing_q') ...
            && schema.reuse_existing_q ...
            && local_has_existing_q(dat)
        return
    end
    options = struct();
    if isfield(bas,'progressFcn') ...
            && ~isempty(bas.progressFcn)
        options.progressFcn = bas.progressFcn;
    end
    [dat,~] = read_hydro_timeseries( ...
        dirQ,bas,split,options,schema,dat);
end

function tf = local_has_existing_q(dat)
    tf = iscell(dat) ...
        && ~isempty(dat);
    for k = 1:numel(dat)
        entry = dat{k};
        if ~isstruct(entry) ...
                || ~isfield(entry,'y_n') ...
                || ~isfield(entry,'bad') ...
                || numel(entry.y_n) ~= numel(entry.bad)
            tf = false;
            return
        end
    end
end
