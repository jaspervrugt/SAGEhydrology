function schema = discharge_schema(region)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DISCHARGE_SCHEMA Load a region-owned discharge schema.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    config = region_helpers('config',root,region);
    if ~isfield(config,'schema') ...
            || ~isfield(config.schema,'discharge')
        error('discharge_schema:MissingSchema', ...
            'Region %s does not define R.schema.discharge.',string(region));
    end
    schema = config.schema.discharge;
end
