function schema = meteo_schema(region)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%METEO_SCHEMA Load a region-owned meteorological schema.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    config = region_helpers('config',root,region);
    if ~isfield(config,'schema') ...
            || ~isfield(config.schema,'meteo')
        error('meteo_schema:MissingSchema', ...
            'Region %s does not define R.schema.meteo.',string(region));
    end
    schema = config.schema.meteo;
end
