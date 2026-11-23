function [A,loc] = select_basin_attributes(A_reg,ID,id_selected)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SELECT_BASIN_ATTRIBUTES Select attributes for prepared basins.
%
%  Match prepared basin identifiers against the regional attribute table and
%  retain columns in the requested run order.
%
% SYNOPSIS:
%   A = select_basin_attributes(A_reg,ID,id_selected)
%   [A,loc] = select_basin_attributes(A_reg,ID,id_selected)
%
% INPUT ARGUMENTS:
%   A_reg           r-by-N matrix of regional basin attributes
%   ID              N identifiers corresponding to columns of A_reg
%   id_selected     K selected identifiers in run order
%
% OUTPUT ARGUMENTS:
%   A               r-by-K matrix in selected order; empty if A_reg is empty
%   loc             K positions of selected basins in ID
%
% NOTES:
%   Missing identifiers or mismatched attribute dimensions cause an error.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Apr. 2026             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    idSelected = normalize_usgs(id_selected);
    idAvailable = normalize_usgs(ID);
    idSelected = idSelected(:);
    idAvailable = idAvailable(:);
    [found,loc] = ismember(idSelected,idAvailable);
    if ~all(found)
        missing = strjoin(cellstr(idSelected(~found).'),', ');
        error('select_basin_attributes:MissingBasin', ...
            ['Prepared basin selection ' ...
            'contains basin(s) not found ' ...
             'in read_attr output: %s'],missing);
    end

    if isempty(A_reg)
        A = [];
        return
    end

    if size(A_reg,2) ~= numel(idAvailable)
        error('select_basin_attributes:AttributeCountMismatch', ...
            ['The number of attribute ' ...
            'columns (%d) does not match ' ...
             'the number of basin identifiers ' ...
             'returned by read_attr (%d).'], ...
            size(A_reg,2),numel(idAvailable));
    end

    A = A_reg(:,loc);
end
