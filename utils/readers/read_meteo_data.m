function [dat,aux] = read_meteo_data(dirM,bas,split,meteo,schema)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_METEO_DATA Read meteorology using a declarative dataset schema.
%
% SYNOPSIS:
%   [dat,aux] = read_meteo_data(dirM,bas,split,meteo,schema)
%
% INPUT:
%   dirM       Directory containing meteorological time-series files
%   bas        Basin structure
%   split      Prepared SAGE time-split structure
%   meteo      Runtime meteorological choices and optional progressFcn
%   schema     Generic time-series schema
%
% OUTPUT:
%   dat        SAGE basin meteorology and optional combined discharge
%   aux        K-by-3 latitude, elevation, and area matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 4 ...
            || isempty(meteo)
        meteo = struct();
    end
    if isfield(schema,'resolve') ...
            && ~isempty(schema.resolve)
        schema = schema.resolve(schema,meteo);
    end
    schema = resolve_hydro_schema(schema,split,meteo);
    mode = 'meteo';
    if isfield(schema,'variables') ...
            && isfield(schema.variables,'Q')
        mode = 'both';
    end
    schema = validate_hydro_schema(schema,mode);
    [dat,aux] = read_hydro_timeseries( ...
        dirM,bas,split,meteo,schema);
end
