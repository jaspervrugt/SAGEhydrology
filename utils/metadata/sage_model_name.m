function name = sage_model_name(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_MODEL_NAME Return a canonical SAGE model name.
%
%  Map supported numeric model identifiers to names. Unknown codes are
%  returned as "model N".
%
% SYNOPSIS:
%   name = sage_model_name(model)
%
% INPUT ARGUMENTS:
%   model           scalar or array of numeric model identifiers
%
% OUTPUT ARGUMENTS:
%   name            char for scalar input; string array otherwise
%
% NOTES:
%   Codes: 1 HYMOD; 2 HMODEL; 3 SAC-SMA; 4 Xinanjiang; 5 GR4J; 6 HBV; 7
%   CFE-NWM; 11 GCHM; 99 user model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ids = [1 2 3 4 5 6 7 11 99];
    names = ["hymod","hmodel","sacsma","Xinanjiang", ...
        "gr4jA","hbv","cfe_nwm","gchm","user_model"];

    model = double(model);
    name = strings(size(model));
    for i = 1:numel(model)
        j = find(ids == model(i),1,'first');
        if isempty(j)
            name(i) = "model " + string(model(i));
        else
            name(i) = names(j);
        end
    end

    if isscalar(name)
        name = char(name);
    end
end
