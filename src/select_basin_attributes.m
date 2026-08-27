function [A,loc] = select_basin_attributes(A_reg,ID,id_selected)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SELECT_BASIN_ATTRIBUTES Select and order attributes for prepared basins.
%
% SYNOPSIS:
%   A = select_basin_attributes(A_reg,ID,id_selected)
%   [A,loc] = select_basin_attributes(A_reg,ID,id_selected)
%
% INPUT:
%   A_reg       r-by-N attribute matrix returned by read_attr
%   ID          N basin identifiers corresponding to columns of A_reg
%   id_selected identifiers of the prepared basins in their run order
%
% OUTPUT:
%   A           r-by-K attributes ordered as id_selected
%   loc         K indices of the selected basins in ID
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
