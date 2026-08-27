function schema = attribute_schema(region)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ATTRIBUTE_SCHEMA Load a region-owned catchment-attribute schema.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    config = region_helpers('config',root,region);
    if ~isfield(config,'schema') ...
            || ~isfield(config.schema,'attributes')
        error('attribute_schema:MissingSchema', ...
            'Region %s does not define R.schema.attributes.', ...
            string(region));
    end
    schema = config.schema.attributes;
end
