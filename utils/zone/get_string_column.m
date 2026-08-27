function s = get_string_column(T,names)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET_STRING_COLUMN Return the first available string table column.
%
% SYNOPSIS:
%   s = get_string_column(T,names)
%
% INPUT:
%   T       Table containing candidate text variables
%   names   Cell/string array with candidate variable names, test in order
%
% OUTPUT:
%   s       String column vector with one value per row of T. If none of 
%           candidate variables are found, s returns strings(height(T),1).
%
% DESCRIPTION:
%   This helper searches a table for the first available candidate variable
%   name and converts that column to a string array. It is used by regional
%   CAMELS attribute readers to retrieve gauge names, station IDs, site
%   labels, or other text metadata whose variable names differ by region.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, June 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    s = strings(height(T),1);
    names = string(names);
    for j = 1:numel(names)
        hit = find(strcmpi( ...
            string(T.Properties.VariableNames), ...
            names(j)),1);
        if ~isempty(hit)
            s = string(T.( ...
                T.Properties.VariableNames{hit}));
            s = s(:);
            return
        end
    end

end
