function [A,id_gauge,gname,zone] = read_attr(region,dirD,bas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ_ATTR Read regional catchment attributes using a declarative schema
%
% SYNOPSIS: [A,id_gauge,gname,zone] = read_attr(region,dirD,bas)
%   region      string with a region defined by region_config_XX.m
%   dirD        string with main directory of regional CAMELS data
%   bas         OPTIONAL: structure basin information
%    .K          number of selected watersheds
%    .K_t        number of training watersheds
%    .K_e        number of evaluation watersheds
%    .id_attr    rx1 vector of integers of basin attributes
%    .pr_attr    optional print attribute table [1/0]
%    .id_gauge   optional requested basin identifiers and output order
%    .dt         optional temporal resolution used by schema profiles
%    .stream     optional release name such as 'daily' or 'hourly'
%   A           OUTPUT: rxK matrix of standardized catchment attributes
%   id_gauge    OUTPUT: Kx1 vector gauge/catchment identifiers
%   gname       OUTPUT: Kx1 vector gauge/catchment names
%   zone        OUTPUT: structure with hydroclimatic basin classification
%    .id         Kx1 string array with zone labels
%    .num        Kx1 numeric vector with integer zone identifiers
%    .names      mx1 string array with unique zone names
%    .aridity    Kx1 vector of aridity index values
%    .frac_snow  Kx1 vector of snow-fraction values
%
% DESCRIPTION:
%   The regional configuration owns the attribute schema. The generic
%   reader applies its file layouts, joins, identifier normalization,
%   metadata rules, transformations, and auxiliary classifications to
%   produce the common SAGE attribute structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3 ...
            || isempty(bas)
        bas = struct();
    end

    if ~isfield(bas,'id_attr') ...
            || isempty(bas.id_attr)
        catalog = attr_catalog(region);
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

end
