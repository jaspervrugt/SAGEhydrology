function R = region_config_NAM()
%REGION_CONFIG_NAM CAMELS-SPAT basin-averaged US and Canadian data.
% NAM avoids the existing NA acronym for Namibia. Dataset paths remain flat.
    sageRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    R = region_helpers('config',sageRoot,'CAMELS_CA');
    R.code = 'CAMELS_SPAT';
    R.acronym = 'NAM';
    R.name = 'North America';
    R.data_root = 'CAMELS_SPAT';
    R.default_resolution = 'Hourly';
    R.basin_file_pattern = 'NAM_%d_basins.txt';
    for field = {'Daily','Hourly'}
        name = field{1}; X = R.by_resolution.(name);
        X.basins.universe = 1426;
        X.basins.file = 'NAM_1426_basins.txt';
        X.basins.training = 1141; X.basins.evaluation = 285;
        if strcmp(name,'Hourly')
            % Full water years avoid the truncated final local day of RDRS.
            % Hourly flow coverage is basin-specific and screened after reading.
            X.period.manual.train_start = '01/10/2014';
            X.period.manual.train_end = '30/09/2018';
            X.period.manual.eval_start = '01/10/2010';
            X.period.manual.eval_end = '30/09/2014';
            X.period.common_start = '01/10/2010';
            X.period.common_end = '30/09/2018';
        end
        for j=1:numel(X.checks)
            X.checks(j).pattern = strrep(X.checks(j).pattern,'CAN_*','*');
            if X.checks(j).minimum_count==764,X.checks(j).minimum_count=1426;end
        end
        R.by_resolution.(name) = X;
    end
    A = R.schema.attributes;
    A.name = 'CAMELS-SPAT North America attributes';
    % The installed metadata retain a country prefix in Station_id.
    A.tables(1).row_filter = struct();
    A.id.regex = {};
    A.region = 'CAMELS_SPAT'; A.zone.region = 'NAM';
    A.progress.label = '... Reading CAMELS-SPAT North America attributes';
    R.schema.attributes = A;
    R.schema.meteo.name = 'CAMELS-SPAT North America';
    R.schema.meteo.profiles.daily.schema.file.pattern = '{gauge}_daymet_lumped.nc';
    R.schema.meteo.profiles.hourly.schema.file.pattern = '{gauge}_rdrs_lumped.nc';
    R.schema.meteo.aux.tables.elevation.column_prefix = {'CAN_','USA_'};
    R.schema.meteo.aux.tables.elevation.strip_prefix = '';
    R.schema.discharge.name = 'CAMELS-SPAT North America discharge';
    R.schema.discharge.profiles.daily.schema.file.pattern = '{gauge}_daily_flow_observations.nc';
    R.schema.discharge.profiles.hourly.schema.file.pattern = '{gauge}_hourly_flow_observations.nc';
    % Publisher quality flags differ between USGS and Canadian observations.
    for field={'daily','hourly'}
        name=field{1};Q=R.schema.discharge.profiles.(name).schema.variables.Q;
        Q.quality.file_prefix_sources.CAN = Q.quality.sources;
        Q.quality.file_prefix_sources.USA = {'q_obs_is_ice_affected', ...
            'q_obs_is_malfunction_affected','q_obs_is_backwater_affected'};
        R.schema.discharge.profiles.(name).schema.variables.Q=Q;
    end
    R.schema.meteo.source_file = [mfilename('fullpath') '.m'];
    R.schema.discharge.source_file = R.schema.meteo.source_file;
end
