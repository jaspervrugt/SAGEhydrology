function [Araw,labels,attrVars] = select_attr_from_catalog(AA, ...
    id_attr,region,pr_table)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SELECT_ATTR_FROM_CATALOG Select catchment attributes using attr_catalog.
%
% SYNOPSIS:
%   [Araw,labels,attrVars] = select_attr_from_catalog(AA, ...
%    id_attr,region,pr_table)
%
% DESCRIPTION:
%   This helper selects attribute columns from the regional catchment
%   attribute table AA using the variable-name order defined in
%   attr_catalog(region). This ensures that GUI attribute IDs, printed
%   labels, and extracted columns all refer to the same attributes.
%
% INPUT:
%   AA          Table with regional catchment attributes
%   id_attr     Selected attribute IDs, matching attr_catalog(region).names
%   region      Region identifier, e.g. 'CAMELS_FR', 'CAMELS_AT'
%   pr_table    Print selected attribute table [1/0]
%
% OUTPUT:
%   Araw        K x r raw attribute matrix before standardization
%   labels      1 x r cell array with selected attribute labels
%   attrVars    r_all x 1 string array with catalog variable names
%
% NOTES:
%   This function does not standardize attributes. Standardization is done
%   in the corresponding read_attr_* function after Araw is returned.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, June 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 4 ...
            || isempty(pr_table)
        pr_table = 1;
    end
    
    region = region_helpers('code',region);
    
    ATTR = attr_catalog(region);
    
    attrVars = string(ATTR.names(:));
    labelsAll = string(ATTR.labels(:));
    
    if numel(labelsAll) ~= numel(attrVars)
        error('select_attr_from_catalog:catalogMismatch', ...
            ['ATTR.labels and ATTR.names ' ...
            'have different lengths for %s.'], ...
            region);
    end
    
    id_attr = id_attr(:).';
    
    if isempty(id_attr)
        error('select_attr_from_catalog:missingIdAttr', ...
            'id_attr is empty for %s.',region);
    end
    
    if any(id_attr < 1) || any(id_attr > numel(attrVars))
        error('select_attr_from_catalog:badIdAttr', ...
            ['id_attr contains indices ' ...
            'outside [1,%d] for %s.'], ...
            numel(attrVars),region);
    end

    tableVars = string(AA.Properties.VariableNames);
    missing = attrVars(~ismember(attrVars,tableVars));
    
    if ~isempty(missing)
        error('select_attr_from_catalog:missingCatalogVars', ...
            ['These attr_catalog variables ' ...
            'are missing from AA for %s: %s'], ...
            region,strjoin(missing,", "));
    end
    
    labels = cellstr(labelsAll(id_attr));
    
    if pr_table == 1
        fprintf('\n');
        fprintf('      %s selectable attributes \n',region);
        fprintf('      --------------------------------\n');
    
        for j = 1:numel(attrVars)
            mark = ' ';
            if any(id_attr == j)
                mark = '*';
            end
            fprintf('      %s %3d  %s\n',mark, ...
                j,char(labelsAll(j)));
        end
    end

    Araw = nan(height(AA),numel(id_attr));
    
    for j = 1:numel(id_attr)
        v = attrVars(id_attr(j));
        Araw(:,j) = local_numeric_attribute(AA.(char(v)),char(v));
    end

end

% =============
% local helpers
% =============
function x = local_numeric_attribute(v,vname)

    if isnumeric(v) ...
            || islogical(v)
        x = double(v(:));
        return
    end
    
    s = string(v);
    s = strip(s);
    
    s(ismissing(s) ...
        | strcmpi(s,"NA") ...
        | strcmpi(s,"NAN") ...
        | strcmpi(s,"NULL") ...
        | s == "") = missing;
    
    x = str2double(s);
    nonmiss = ~ismissing(s);
    numericOK = ~nonmiss | isfinite(x);
    
    if all(numericOK)
        x = x(:);
        return
    end
    
    error('select_attr_from_catalog:notNumeric', ...
        ['Attribute %s is not numeric. ' ...
        'Mark it unselectable in attr_catalog.'], ...
        vname);
end