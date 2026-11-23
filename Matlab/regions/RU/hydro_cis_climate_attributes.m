function T = hydro_cis_climate_attributes(seriesDir,ids,outputFile,latitude)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HYDRO_CIS_CLIMATE_ATTRIBUTES Derive Russia climate-zone attributes.
%
% SYNOPSIS:
%   T = hydro_cis_climate_attributes(seriesDir,ids)
%   T = hydro_cis_climate_attributes(seriesDir,ids,outputFile)
%
% INPUT:
%   seriesDir   Folder containing one consolidated HydroCIS NetCDF file
%               per gauge
%   ids         Gauge identifiers to process
%   outputFile  OPTIONAL CSV file receiving the derived table
%   latitude    OPTIONAL gauge-latitude vector; otherwise read from the
%               installed gauge_information.txt
%
% OUTPUT:
%   T           Table with mean temperature, annual precipitation and PET,
%               aridity, and precipitation-weighted snow fraction
%
% DESCRIPTION:
%   The published HydroCIS static table omits the climate summaries needed
%   by SAGE's global basin-zone classifier. They are derived once from the
%   common 2008--2020 daily record. ERA5-Land is preferred, followed by
%   ERA5; MSWEP supplies precipitation where both ERA products are absent.
%   Potential evaporation uses the Oudin formulation selected by default
%   in the GUI. Source gaps remain NaN. The classifier applies its
%   documented per-basin fallback, while read_attr excludes an affected
%   gauge from the eligible pool if the user selects an incomplete climate
%   summary as an FFN predictor.
%   Snow fraction is the fraction of ERA5-Land precipitation occurring on
%   days with mean air temperature at or below 0 degrees Celsius.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Sept. 2026                               %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3
        outputFile = '';
    end
    ids = strip(string(ids(:)));
    K = numel(ids);
    if nargin < 4 || isempty(latitude)
        dataRoot = fileparts(fileparts(seriesDir));
        G = readtable(fullfile(dataRoot,'gauge_information.txt'), ...
            'TextType','string','Delimiter','\t');
        sourceIds = strip(string(G.gauge_id));
        [found,location] = ismember(ids,sourceIds);
        if ~all(found)
            error('HYDRO_CIS:ClimateLatitude', ...
                'Gauge information is missing at least one basin ID.');
        end
        latitude = double(G.gauge_lat(location));
    else
        latitude = double(latitude(:));
    end
    if numel(latitude) ~= K || any(~isfinite(latitude))
        error('HYDRO_CIS:ClimateLatitude', ...
            'Finite latitude is required for every HydroCIS gauge.');
    end
    t_mean = nan(K,1);
    p_mean = nan(K,1);
    pet_mean = nan(K,1);
    aridity = nan(K,1);
    frac_snow = nan(K,1);
    daysPerYear = 365.2425;

    for k = 1:K
        file = fullfile(seriesDir,char(ids(k) + ".nc"));
        if ~isfile(file)
            error('HYDRO_CIS:MissingClimateSeries', ...
                'Missing HydroCIS time-series file: %s',file);
        end
        P = local_best_series(file, ...
            {'prcp_e5l','prcp_e5','prcp_mswep'});
        Tmin = local_best_series(file,{'t_min_e5l','t_min_e5'});
        Tmax = local_best_series(file,{'t_max_e5l','t_max_e5'});
        if ~isequal(numel(P),numel(Tmin),numel(Tmax))
            error('HYDRO_CIS:ClimateSeriesLength', ...
                'Climate variables differ in length for gauge %s.',ids(k));
        end

        % Tiny negative precipitation and evaporation values arise from
        % floating-point encoding of nominal zeros in the source files.
        P(isfinite(P)) = max(P(isfinite(P)),0);
        Tday = 0.5*(Tmin + Tmax);
        time = datetime(2008,1,1) + days((0:numel(Tday)-1).');
        doy = day(time,'dayofyear');
        phi = deg2rad(latitude(k));
        dr = 1 + 0.033*cos(2*pi*doy/365);
        delta = 0.409*sin(2*pi*doy/365 - 1.39);
        ws = acos(max(-1,min(1,-tan(phi).*tan(delta))));
        Ra = (24*60/pi)*0.0820.*dr.*(ws.*sin(phi).*sin(delta) ...
            + cos(phi).*cos(delta).*sin(ws));
        Ep = (0.408/100).*Ra.*max(Tday + 5,0);
        Ep(~isfinite(Tday)) = NaN;
        validP = isfinite(P);
        validT = isfinite(Tday);
        validEp = isfinite(Ep);
        if nnz(validT) >= 365
            t_mean(k) = mean(Tday(validT));
        end
        if nnz(validP) >= 365
            p_mean(k) = daysPerYear*mean(P(validP));
        end
        if nnz(validEp) >= 365
            pet_mean(k) = daysPerYear*mean(Ep(validEp));
        end
        if isfinite(p_mean(k)) && p_mean(k) > 0 && isfinite(pet_mean(k))
            aridity(k) = pet_mean(k)/p_mean(k);
        end
        validSnow = validP & validT;
        denominator = sum(P(validSnow));
        if denominator > 0
            frac_snow(k) = sum(P(validSnow & Tday <= 0))/denominator;
        end

        if mod(k,100) == 0 || k == K
            fprintf('... Derived HydroCIS climate attributes: %d/%d\n',k,K);
        end
    end

    gauge_id = ids;
    T = table(gauge_id,t_mean,p_mean,pet_mean,aridity,frac_snow);
    if strlength(string(outputFile)) > 0
        writetable(T,char(string(outputFile)), ...
            'FileType','text','Encoding','UTF-8');
    end
end

function value = local_best_series(file,candidates)
    value = [];
    for i = 1:numel(candidates)
        candidate = double(ncread(file,candidates{i}));
        candidate = candidate(:);
        if isempty(value) || nnz(isfinite(candidate)) > nnz(isfinite(value))
            value = candidate;
        end
        if nnz(isfinite(value)) >= 365
            return
        end
    end
end
