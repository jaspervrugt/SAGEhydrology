function x = get_numeric_column(T,cands)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET_NUMERIC_COLUMN Return the first available numeric table column.
%
% SYNOPSIS:
%   x = get_numeric_column(T,cands)
%
% INPUT:
%   T       Table containing candidate attribute variables
%   cands   Cell/string array with candidate variable names, tested in order
%
% OUTPUT:
%   x       Numeric column vector with one value per row of T. If none of
%           the candidate variables are found, x is returned as
%           NaN(height(T),1).
%
% DESCRIPTION:
%   This helper searches a table for the first available candidate variable
%   name and converts column to numeric form using get_numeric_attribute.
%   It is used by regional CAMELS attribute readers to retrieve common
%   fields such as area, elevation, precipitation, PET, latitude, 
%   longitude, or other attributes that may have region-specific names.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, June 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = nan(height(T),1);
    vars = string(T.Properties.VariableNames);
    varClean = lower(regexprep(vars,'[^a-zA-Z0-9]',''));
    cands = string(cands);
    
    for i = 1:numel(cands)
        c = lower(regexprep(cands(i),'[^a-zA-Z0-9]',''));
        j = find(varClean == c,1,'first');
        if isempty(j)
            j = find(contains(varClean,c),1,'first');
        end
        if ~isempty(j)
            y = T{:,j};
            if isnumeric(y) ...
                    || islogical(y)
                x = double(y(:));
            else
                x = str2double(strrep(string(y(:)),',','.'));
            end
            return
        end
    end

end
