function R=region_config_EE()
%REGION_CONFIG_EE Defaults for Estonia GRDC-Caravan.
    R=region_config_GRDC('EE','Estonia',51,51,42,9, ...
        '01/10/1972','30/09/1992','01/10/1992','30/09/2012');

    % Declarative generic-reader schemas owned by this region module.
    R.schema.meteo = local_meteo_schema();
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge = derive_discharge_schema(R.schema.meteo);
    R.schema.discharge.source_file = R.schema.meteo.source_file;

end

function S = local_meteo_schema()
    S.name = 'daily GRDC-Caravan';
    S.format = 'csv';
    S.layout = 'one_file_per_basin';
    S.file.pattern = 'GRDC_{gauge}.csv';
    S.file.delimiter = ',';
    S.id.pad_width = 0;
    S.timeline.reference = datetime(1950,10,1);
    S.timeline.step = days(1);
    S.time.column = 'date';
    S.time.input_format = '';
    S.variables.P = local_variable( ...
        'total_precipitation_sum','mm/day','mm/day');
    S.variables.P.clip_min = 0;
    S.variables.Ep = local_variable( ...
        'potential_evaporation_sum_FAO_PENMAN_MONTEITH', ...
        'mm/day','mm/day');
    S.variables.Ep.clip_min = 0;
    S.variables.T = local_variable( ...
        'temperature_2m_mean','degC','degC');
    S.variables.Q = local_variable( ...
        'streamflow','mm/day','mm/day');
    S.progress.label = ...
        '... Reading daily GRDC-Caravan generic data';
end

function V = local_variable(source,units,targetUnits)
    V.source = source;
    V.units = units;
    V.target_units = targetUnits;
end
