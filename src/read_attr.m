function [A,id_gauge,gname,zone] = read_attr(region,dirD,bas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_ATTR Read regional catchment attributes.
%
%  Uses the regional schema to join, transform, and standardize basin
%  attributes.
%
% SYNOPSIS:
%   [A,id_gauge,gname,zone] = read_attr(region,dirD,bas)
%
% INPUT ARGUMENTS:
%   region          region supported by a region_config_* function
%   dirD            regional CAMELS data directory
%   bas             optional basin-selection and attribute settings
%    .id_attr        selected attribute indices
%    .id_gauge       requested basin IDs and output order
%    .dt             optional temporal resolution
%    .stream         optional release name
%
% OUTPUT ARGUMENTS:
%   A               standardized attribute matrix, r-by-K
%   id_gauge        K-by-1 gauge or catchment IDs
%   gname           K-by-1 basin names
%   zone            hydroclimatic basin classification
%    .id             K-by-1 zone labels
%    .num            K-by-1 zone numbers
%    .names          unique zone names
%    .aridity        K-by-1 aridity index
%    .frac_snow      K-by-1 snow fraction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3 ...
            || isempty(bas)
        bas = struct();
    end

    catalog = attr_catalog(region);
    if ~isfield(bas,'id_attr') ...
            || isempty(bas.id_attr)
        if isfield(catalog,'default_ids') ...
                && ~isempty(catalog.default_ids)
            bas.id_attr = catalog.default_ids;
        elseif isfield(catalog,'selectable') ...
                && ~isempty(catalog.selectable)
            bas.id_attr = find(catalog.selectable);
        else
            bas.id_attr = 1:numel(catalog.names);
        end
    end

    schema = attribute_schema(region);
    [A,id_gauge,gname,zone] = read_attribute_data( ...
        dirD,bas,schema);

    % Some regional catalogs intentionally expose useful optional fields
    % with incomplete source coverage. When that catalog opts in, exclude
    % affected gauges from the eligible sampling pool rather than
    % imputing predictors or failing later after a random split happens to
    % include one. Complete/default selections retain the full inventory.
    excludeIncomplete = isfield(catalog, ...
        'exclude_nonfinite_selected') ...
        && ~isempty(catalog.exclude_nonfinite_selected) ...
        && logical(catalog.exclude_nonfinite_selected(1));
    if excludeIncomplete && ~isempty(A)
        incomplete = any(~isfinite(A),1);
        if any(incomplete)
            removed = id_gauge(incomplete);
            keep = ~incomplete;
            A = A(:,keep);
            id_gauge = id_gauge(keep);
            gname = gname(keep);
            zone = local_subset_zone(zone,keep);
            warning('read_attr:IncompleteSelectedAttributes', ...
                ['Excluded %d %s gauge(s) from the eligible basin pool ' ...
                 'because at least one selected optional attribute is ' ...
                 'unavailable. Gauge IDs: %s.'], ...
                nnz(incomplete),char(region_helpers('name',region)), ...
                strjoin(cellstr(string(removed(:)).'),', '));
        end
    end

end

function subset = local_subset_zone(zone,keep)
%LOCAL_SUBSET_ZONE Keep zone vectors aligned with retained gauge rows.

    subset = zone;
    if ~isstruct(zone) || isempty(fieldnames(zone))
        return
    end
    n = numel(keep);
    fields = fieldnames(zone);
    for i = 1:numel(fields)
        name = fields{i};
        value = zone.(name);
        if isvector(value) && numel(value) == n
            subset.(name) = value(keep);
        end
    end
    if isfield(subset,'id') && ~isempty(subset.id)
        [names,~,number] = unique(string(subset.id(:)),'stable');
        subset.names = names;
        subset.num = number;
    end
end
