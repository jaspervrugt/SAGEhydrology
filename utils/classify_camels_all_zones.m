function zone = classify_camels_all_zones(region,AA,min_per_zone,max_zones)
%CLASSIFY_CAMELS_ALL_ZONES Dispatch hydroclimatic zones by region.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CLASSIFY_CAMELS_ZONES Classify CAMELS basins into hydroclimatic zones
%
% SYNOPSIS: zone = classify_camels_zones(AA)
%   AA          table with CAMELS catchment attributes for all basins
%               (rows = basins, columns = attributes). Must contain at least:
%                 .aridity      aridity index (PET / P or similar)
%                 .frac_snow    fraction of precipitation falling as snow
%   zone        OUTPUT: structure with hydroclimatic classification
%    .id         Kx1 string array with zone labels per basin
%                (e.g., "humid_rain", "dry_snow")
%    .num        Kx1 numeric vector with integer zone identifiers
%    .names      Mx1 string array with unique zone names
%    .aridity    Kx1 vector of aridity index values (copied from AA)
%    .frac_snow  Kx1 vector of snow fraction values (copied from AA)
%
% DESCRIPTION:
%   This function partitions CAMELS basins into hydroclimatic zones based
%   on a combination of aridity and snow dominance. The classification is
%   intentionally simple, interpretable, and robust, making it suitable for:
%
%     - stratified sampling of training/evaluation basins,
%     - regional performance diagnostics (ECDFs, maps),
%     - sensitivity and attribution analysis within SAGE.
%
%   The current implementation defines zones using threshold-based rules:
%
%     Aridity classes:
%       aridity < 1        → humid
%       1 ≤ aridity < 2    → subhumid
%       aridity ≥ 2        → dry
%
%     Snow regime:
%       frac_snow ≥ 0.30   → snow-dominated
%       frac_snow <  0.30  → rain-dominated
%
%   These are combined into composite labels such as:
%       "humid_rain", "subhumid_snow", "dry_rain", etc.
%
%   The thresholds can be adjusted or replaced by alternative schemes
%   (e.g., Köppen classification, clustering, or PCA-based grouping)
%   depending on the application.
%
% NOTES:
%   - Basins with missing or non-finite attributes are assigned to a
%     "missing" class.
%   - The returned zone numbering is consistent with the order of first
%     occurrence of each unique class (stable labeling).
%   - This function is designed to operate on the full CONUS dataset prior
%     to basin sampling to ensure consistent stratification.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    method = 1;     % 1: global using common hydroclimatic scheme.
                    % 2: local using region-specific legacy rules.
    switch method
        case 1
            zone = classify_camels_global_zones(region,AA, ...
                min_per_zone,max_zones);
        case 2
            zone = classify_camels_local_zones(region,AA);
        otherwise
            error('Unknown zone classification method: %d',method);
    end

end

% =========================================================
function zone = classify_camels_global_zones(region,AA, ...
    min_per_zone,max_zones)
% =========================================================
% Global hydroclimatic classifier (including CAMELS-PE)
%
% Classification based on four independent descriptors:
%   (1) Moisture regime (aridity index)
%   (2) Latitude (climate zone)
%   (3) Topography (mean elevation)
%   (4) Snow influence
%
% These descriptors are combined into a limited number
% of globally consistent hydroclimatic classes.
% =====================================================

    region = lower(string(region));
    K = height(AA);

    switch lower(region)

        case "fr"
            [lat,lon] = local_fr_coordinates(AA);

        otherwise
            lat = get_geo_column(AA, ...
                {'gauge_lat_dd','gauge_lat','lat','latitude', ...
                'lat_outlet','lat_wgs84','Lat','cwc_lat','ghi_lat', ...
                'centroid_lat','gauge_latitude','latitude_deg', ...
                'sta_y_w84','sta_y_w84_snap','sit_latitude', ...
                'Station_lat','station_lat'}, ...
                -90,90);

            lon = get_geo_column(AA, ...
                {'gauge_lon_dd','gauge_lon','lon','longitude', ...
                'long_outlet','lon_wgs84','Lon','cwc_lon','ghi_lon', ...
                'centroid_lon','gauge_longitude','longitude_deg', ...
                'sta_x_w84','sta_x_w84_snap','sit_longitude', ...
                'Station_lon','station_lon'}, ...
                -180,180);

    end
    
    area = get_numeric_column(AA, ...
        {'area','area_km2','catchment_area','area_calc', ...
         'area_gov','area_gages2','cwc_area','ghi_area', ...
         'sit_area_topo','sit_area_hydro','basin_area', ...
         'catch_area','catchment_area_km2','area_total', ...
         'Basin_area_km2','basin_area_km2'});
    
    elev = get_numeric_column(AA, ...
        {'elev_mean','mean_ele','gauge_elev','elev_med', ...
         'Elevation_mabsl','Z_MEAN','elv_mean','ele_mt_sav', ...
         'top_altitude_mean','sit_altitude','elevation_mean', ...
         'mean_elevation','elev','elevation','elev_m', ...
         'elev_mean'});
    
    p_mean = get_numeric_column(AA, ...
        {'p_mean','P_stn_sum','p_mean_mm','pre_mean', ...
         'precip_mean','precipitation_mean','mean_precip', ...
         'p_mean_annual','precip_mm_syr','annual_precipitation', ...
         'PR0_mean','prec_mean','Pmean_mm_year'});
    
    pet_mean = get_numeric_column(AA, ...
        {'pet_mean','PET_PM_sum','et0_mean','pet_mean_mm', ...
         'potential_evaporation_mean','pet_mean_annual', ...
         'pet_mm_syr','eto_mean','mean_pet','annual_pet', ...
         'pet1_mean','pet2_mean'});
    
    aridity = get_numeric_column(AA, ...
        {'aridity','aridity_index','aridity_catch', ...
         'arid_1','arid_2','aridity_p_pet','ai_mean', ...
         'aridity_ratio','pet_p_ratio','P_PET_PM','p_pet_ratio', ...
         'aridity1_mean','aridity2_mean'});
    
    ix = ~isfinite(aridity) & isfinite(p_mean) & ...
        isfinite(pet_mean) & p_mean > 0;
    aridity(ix) = pet_mean(ix) ./ p_mean(ix);

    % CAMELS-PL provides mean daily precipitation and snow fraction, but
    % no static PET or conventional PET/P aridity attribute. Preserve the
    % established CAMELS-PL moisture thresholds used by the legacy local
    % classifier and map them onto the common global aridity classes:
    %   p_mean < 1.6 mm d^-1       -> Dry
    %   1.6 <= p_mean < 2.3        -> Humid
    %   p_mean >= 2.3              -> Wet
    % The representative values below are class codes on the common
    % aridity scale; they are not measured PET/P ratios.
    if region == "pl"
        ixPL = ~isfinite(aridity) & isfinite(p_mean);
        aridity(ixPL & p_mean < 1.6) = 2.0;
        aridity(ixPL & p_mean >= 1.6 & p_mean < 2.3) = 1.0;
        aridity(ixPL & p_mean >= 2.3) = 0.5;
    end

    % CAMELS-FR does not provide a static climate/aridity attribute file,
    % nor a precipitation/PET pair in the supplied catchment attributes.
    % The legacy France zoning is geographic/topographic and treats the
    % moisture regime as humid. Encode that documented fallback explicitly
    % so France does not emit a misleading missing-attribute warning.
    % The value 1.0 is a class code on the common aridity scale; it is not
    % claimed to be a measured PET/P ratio.
    if region == "fr"
        ixFR = ~isfinite(aridity);
        aridity(ixFR) = 1.0;
    end

    % CAMELS-SE supplies annual precipitation and mean temperature, but no
    % static PET or aridity attribute. Its legacy regional classifier uses
    % temperature and elevation and does not distinguish moisture classes.
    % Retain that behavior within the global scheme by explicitly assigning
    % the humid class. The value is a class code, not a measured PET/P ratio.
    if region == "se"
        ixSE = ~isfinite(aridity);
        aridity(ixSE) = 1.0;
    end

    % CAMELS-DE provides p_mean and frac_snow, but no static PET or
    % conventional PET/P aridity attribute in either the daily or hourly
    % climatic table. Use precipitation terciles, as in the legacy
    % CAMELS-DE classifier, to define relative moisture classes. Because
    % terciles are calculated within the active attribute table, this works
    % consistently for daily p_mean [mm d^-1] and hourly p_mean [mm h^-1].
    % The assigned values are class codes on the common aridity scale; they
    % are not measured PET/P ratios.
    if region == "de"
        ixDE = ~isfinite(aridity) & isfinite(p_mean);
        pDE = p_mean(ixDE);
        if ~isempty(pDE)
            qDE = quantile(pDE,[1/3 2/3]);
            aridity(ixDE & p_mean <= qDE(1)) = 2.0;  % Dry
            aridity(ixDE & p_mean > qDE(1) ...
                & p_mean <= qDE(2)) = 1.3;           % Semi-arid
            aridity(ixDE & p_mean > qDE(2)) = 0.8;  % Humid
        end
    end
    % CAMELS-PL uses the exact field name frac_snow.
    
    frac_snow = get_numeric_column(AA, ...
        {'frac_snow','snow_frac','snow_fraction', ...
         'fraction_snow','snowfall_fraction','frac_precip_snow', ...
         'fracsnow1_mean','fracsnow2_mean'});
    
    Tmean = get_numeric_column(AA, ...
        {'Tmean_C','t_mean','temp_mean','temperature_mean', ...
         'mean_temperature','tmean','air_temperature_mean', ...
         'temperature_2m_mean','TT_mean','T_stn_mean'});
    
    ix = ~isfinite(frac_snow) ...
        & isfinite(Tmean);
    frac_snow(ix) = max(0,min(1,(2 - Tmean(ix))/8));

    if all(~isfinite(lat))
        warning('classify_camels_all_zones:latitude', ...
            'No usable latitude attribute found for region %s.',region);
    end

    if all(~isfinite(lon))
        warning('classify_camels_all_zones:longitude', ...
            'No usable longitude attribute found for region %s.',region);
    end

    if all(~isfinite(elev))
        warning('classify_camels_all_zones:elevation', ...
            ['No usable elevation attribute found for region %s. ' ...
             'Lowland will be used as the topographic fallback.'],region);
    end

    if all(~isfinite(aridity))
        warning('classify_camels_all_zones:aridity', ...
            ['No usable aridity attribute, or precipitation/PET pair, ' ...
             'was found for region %s. Humid will be used as the ' ...
             'moisture-regime fallback.'],region);
    end

    if all(~isfinite(frac_snow))
        warning('classify_camels_all_zones:snow', ...
            ['No usable snow-fraction or mean-temperature attribute ' ...
             'was found for region %s. Basins will default to ' ...
             'rain-dominated unless mountain rules apply.'],region);
    end

    classID = nan(K,1);
    for k = 1:K
    
        %% 1. Moisture regime
        if isfinite(aridity(k))
            if aridity(k) >= 2.0
                moisture = "Dry";
            elseif aridity(k) >= 1.2
                moisture = "Semi-arid";
            elseif aridity(k) >= 0.65
                moisture = "Humid";
            else
                moisture = "Wet";
            end
        else
            moisture = "Humid";
        end
    
        %% 2. Latitude band
        if isfinite(lat(k))
            alat = abs(lat(k));
            if alat < 23
                climate = "Tropical";
            elseif alat < 35
                climate = "Subtropical";
            elseif alat < 55
                climate = "Temperate";
            else
                climate = "Boreal";
            end
        else
            climate = "Temperate";
        end
    
        %% 3. Topography
        if isfinite(elev(k))
            if elev(k) >= 1500
                topo = "Mountain";
            elseif elev(k) >= 500
                topo = "Upland";
            else
                topo = "Lowland";
            end
        else
            topo = "Lowland";
        end
    
        %% 4. Snow influence
        snow = false;
        if isfinite(frac_snow(k))
            snow = frac_snow(k) > 0.25;
        end
        if isfinite(Tmean(k))
            snow = snow ...
                || Tmean(k) < 1.0;
        end

        %% 5. Integer classification
        if topo == "Mountain"    
            if moisture == "Dry"
                classID(k) = 4;                 % Dry Mountain
            elseif snow
                classID(k) = 2;                 % Cold Mountain
            else
                classID(k) = 3;                 % Mountain
            end
        elseif snow
            classID(k) = 1;
        elseif topo == "Upland"
    
            switch moisture
                case "Dry"
                    classID(k) = 8;             % Dry Upland
                case "Semi-arid"
                    classID(k) = 7;             % Semi-arid Upland
                case "Wet"
                    classID(k) = 9;             % Wet Upland
                otherwise
                    classID(k) = 6;             % Humid Upland
            end    
        else  % Lowland    
            switch moisture
                case "Dry"
                    if climate == "Tropical"
                        classID(k) = 10;        % Hot Dry
                    else
                        classID(k) = 11;        % Dry
                    end    
                case "Semi-arid"
                    classID(k) = 12;            % Semi-arid    
                case "Wet"
                    if climate == "Tropical"
                        classID(k) = 13;        % Tropical Wet
                    else
                        classID(k) = 14;        % Wet
                    end    
                otherwise  % Humid
                    switch climate
                        case "Tropical"
                            classID(k) = 15;    % Tropical Humid
                        case "Subtropical"
                            classID(k) = 16;    % Subtropical Humid
                        case "Temperate"
                            classID(k) = 17;    % Temperate Humid
                        otherwise
                            classID(k) = 18;    % Boreal Humid
                    end
            end
        end
    end

    %% 6. Labeling
    labelMap = strings(18,1);
    labelMap(1) = "Cold Snow";
    labelMap(2) = "Cold Mountain";
    labelMap(3) = "Mountain";
    labelMap(4) = "Dry Mountain";
    labelMap(5) = "Reserved";
    labelMap(6) = "Humid Upland";
    labelMap(7) = "Semi-arid Upland";
    labelMap(8) = "Dry Upland";
    labelMap(9) = "Wet Upland";
    labelMap(10) = "Hot Dry";
    labelMap(11) = "Dry";
    labelMap(12) = "Semi-arid";
    labelMap(13) = "Tropical Wet";
    labelMap(14) = "Wet";
    labelMap(15) = "Tropical Humid";
    labelMap(16) = "Subtropical Humid";
    labelMap(17) = "Temperate Humid";
    labelMap(18) = "Boreal Humid";

    zone = finalize_zone_from_classID(region,classID,labelMap, ...
        lat,lon,area,elev,p_mean,pet_mean,aridity,frac_snow);
    zone = merge_small_zones(zone,min_per_zone);
    zone = compress_zones(zone,max_zones);

end

% =====================================================================
function zone = finalize_zone_from_classID(region,classID,labelMap, ...
    lat,lon,area,elev,p_mean,pet_mean,aridity,frac_snow)
% =====================================================================

    K = numel(classID);

    bad = ~isfinite(classID) ...
        | classID < 1 ...
        | classID > numel(labelMap);

    id = strings(K,1);
    id(~bad) = labelMap(classID(~bad));
    id(bad) = "Unknown";

    % This is the actual merging step
    key = lower(strtrim(id));
    key = regexprep(key,'\s+',' ');

    usedKey = unique(key,'stable');

    num = nan(K,1);
    usedNames = strings(numel(usedKey),1);

    for j = 1:numel(usedKey)
        ix = key == usedKey(j);
        num(ix) = j;
        usedNames(j) = id(find(ix,1,'first'));
    end

    zone = struct();
    zone.id = id(:);
    zone.num = num(:);
    zone.codes = usedNames(:);
    zone.names = usedNames(:);
    zone.long = usedNames(:);

    zone.lat = lat(:);
    zone.lon = lon(:);
    zone.region = repmat(lower(string(region)),K,1);
    zone.area = area(:);
    zone.elev_mean = elev(:);
    zone.p_mean = p_mean(:);
    zone.pet_mean = pet_mean(:);
    zone.aridity = aridity(:);
    zone.frac_snow = frac_snow(:);

end

% ====================================================
function zone = classify_camels_local_zones(region,AA)
% ====================================================

    region = lower(string(region));
    f = str2func("classify_camels_" + region + "_zones");
    
    try
        zone = f(AA);
    catch ME
        error('classify_camels_all_zones:localRegion', ...
            ['No usable local zone classifier was found for region %s. ' ...
             'Expected function: %s.m\n%s'], ...
            region,string(func2str(f)),ME.message);
    end

end

% ==========================================
function zone = classify_camels_us_zones(AA)
% ==========================================

    aridity = AA.aridity;
    frac_snow = AA.frac_snow;

    K = height(AA);

    zone_id = strings(K,1);
    zone_name = strings(K,1);

    for k = 1:K
        if ~isfinite(aridity(k)) ...
                || ~isfinite(frac_snow(k))
            zone_id(k) = "missing";
            zone_name(k) = "Missing attributes";
            continue
        end

        if aridity(k) < 1
            aridClass = "humid";
        elseif aridity(k) < 2
            aridClass = "subhumid";
        else
            aridClass = "dry";
        end

        if frac_snow(k) >= 0.30
            snowClass = "snow";
        else
            snowClass = "rain";
        end

        zone_id(k) = aridClass + "_" + snowClass;
    end

    u = unique(zone_id,'stable');
    zone_num = nan(K,1);
    for j = 1:numel(u)
        zone_num(zone_id == u(j)) = j;
    end

    zone = struct();
    zone.id = zone_id;
    zone.num = zone_num;
    zone.names = u;
    zone.aridity = aridity;
    zone.frac_snow = frac_snow;
end

% ==========================================
function zone = classify_camels_at_zones(AA)
% ==========================================

    n = height(AA);
    z = strings(n,1);

    AI = get_numeric_column(AA, ...
        {'arid_1','aridity', ...
        'aridity_catch','arid_2'});
    Z = get_numeric_column(AA, ...
        {'elev_mean','elev_med', ...
        'gauge_elev'});

    for i = 1:n
        if isfinite(AI(i))
            if AI(i) < 0.65
                clim = "dry";
            elseif AI(i) < 1.5
                clim = "humid";
            else
                clim = "wet";
            end
        else
            clim = "climate_unknown";
        end

        if isfinite(Z(i))
            if Z(i) < 600
                topo = "lowland";
            elseif Z(i) < 1400
                topo = "montane";
            else
                topo = "alpine";
            end
        else
            topo = "topography_unknown";
        end

        z(i) = clim + "_" + topo;
    end

    zone.id = z(:);
    zone.name = z(:);
    zone.label = z(:);
    zone.group = z(:);

    zone.lat = get_numeric_column(AA, ...
        {'gauge_lat','lat', ...
        'latitude','lat_wgs84'});
    zone.lon = get_numeric_column(AA, ...
        {'gauge_lon','lon', ...
        'longitude','lon_wgs84'});
    zone.area = get_numeric_column(AA, ...
        {'area_calc','area_gov', ...
        'area','area_km2'});
    zone.elev = get_numeric_column(AA, ...
        {'gauge_elev','elev','elevation','elev_mean'});
end

% ==========================================
function zone = classify_camels_au_zones(AA)
% ==========================================

    K = height(AA);
    
    division = strings(K,1);
    river_region = strings(K,1);
    aridity = nan(K,1);
    frac_snow = nan(K,1);
    lat = nan(K,1);
    lon = nan(K,1);
    area = nan(K,1);
    
    if ismember('drainage_division', ...
            AA.Properties.VariableNames)
        division = string( ...
            AA.drainage_division);
    end
    if ismember('river_region', ...
            AA.Properties.VariableNames)
        river_region = string( ...
            AA.river_region);
    end
    if ismember('aridity', ...
            AA.Properties.VariableNames)
        aridity = double(AA.aridity);
    end
    if ismember('frac_snow', ...
            AA.Properties.VariableNames)
        frac_snow = double(AA.frac_snow);
    end
    if ismember('lat_outlet', ...
            AA.Properties.VariableNames)
        lat = double(AA.lat_outlet);
    end
    if ismember('long_outlet', ...
            AA.Properties.VariableNames)
        lon = double(AA.long_outlet);
    end
    if ismember('catchment_area', ...
            AA.Properties.VariableNames)
        area = double(AA.catchment_area);
    end
    
    id = strings(K,1);
    long = strings(K,1);
    
    for k = 1:K
        if isfinite(aridity(k))
            if aridity(k) < 0.80
                aLong = "humid";
            elseif aridity(k) < 1.50
                aLong = "subhumid";
            else
                aLong = "dry";
            end
        else
            aLong = "unknown aridity";
        end
    
        if isfinite(frac_snow(k)) ...
                && frac_snow(k) >= 0.05
            sLong = "snow";
        else
            sLong = "rain";
        end
    
        id(k) = zone_title_AU(aLong + "_" + sLong);
        long(k) = aLong + " / " + sLong;
    end
    
    codes = unique(id,'stable');
    names = strings(size(codes));
    num = nan(K,1);
    
    for j = 1:numel(codes)
        tf = id == codes(j);
        num(tf) = j;
        jj = find(tf,1,'first');
        names(j) = zone_title_AU(long(jj));
    end
    
    zone.id = id;
    zone.long = long;
    zone.num = num;
    zone.codes = codes;
    zone.names = names;
    zone.lat = lat;
    zone.lon = lon;
    zone.area = area;
    zone.division = division;
    zone.river_region = river_region;
    zone.aridity = aridity;
    zone.frac_snow = frac_snow;

end

% ==========================================
function zone = classify_camels_br_zones(AA)
% ==========================================

    K = height(AA);
    
    runoff_ratio = nan(K,1);
    forest_perc = nan(K,1);
    lat = nan(K,1);
    
    if ismember('runoff_ratio', ...
            AA.Properties.VariableNames)
        runoff_ratio = ...
            double(AA.runoff_ratio);
    end
    
    if ismember('forest_perc', ...
            AA.Properties.VariableNames)
        forest_perc = ...
            double(AA.forest_perc);
    end
    
    if ismember('gauge_lat', ...
            AA.Properties.VariableNames)
        lat = ...
            double(AA.gauge_lat);
    end
    
    id = strings(K,1);
    long = strings(K,1);
    
    for k = 1:K
    
        if ~isfinite(runoff_ratio(k))
            rLong = "unknown runoff";
        elseif runoff_ratio(k) < 0.30
            rLong = "low runoff";
        elseif runoff_ratio(k) < 0.60
            rLong = "moderate runoff";
        else
            rLong = "high runoff";
        end
    
        if ~isfinite(forest_perc(k))
            cLong = "unknown cover";
        elseif forest_perc(k) >= 50
            cLong = "forest";
        else
            cLong = "mixed";
        end
    
        id(k) = zone_title_BR(rLong + " " + cLong);
        long(k) = rLong + " / " + cLong;
    end
    
    codes = unique(id,'stable');
    names = strings(size(codes));
    num = nan(K,1);
    
    for j = 1:numel(codes)
        tf = id == codes(j);
        num(tf) = j;
    
        % Use the first descriptive label for this class
        jj = find(tf,1,'first');
        names(j) = zone_title_BR(long(jj));
    end
    
    zone.id = id;                  % compact class code, e.g. L-F
    zone.long = long;              % descriptive label, e.g. low runoff / forest
    zone.num = num;
    zone.codes = codes;            % compact legend code if needed
    zone.names = names;            % legend labels, e.g. Low runoff forest
    zone.lat = lat;
    zone.runoff_ratio = runoff_ratio;
    zone.forest_perc = forest_perc;

end

% ==========================================
function zone = classify_camels_ch_zones(AA)
% ==========================================

    K = height(AA);
    lat = get_numeric_column(AA, ...
        ["gauge_lat","lat", ...
        "lat_outlet"]);
    lon = get_numeric_column(AA, ...
        ["gauge_lon","lon","long", ...
        "lon_outlet"]);
    area = get_numeric_column(AA, ...
        ["area","catchment_area"]);
    aridity = get_numeric_column(AA, ...
        ["aridity","aridity_index"]);
    frac_snow = get_numeric_column(AA, ...
        ["frac_snow","snow_frac", ...
        "snow_fraction"]);

    id = strings(K,1);
    long = strings(K,1);

    for k = 1:K
        if isfinite(aridity(k))
            if aridity(k) < 0.80
                aLong = "humid";
            elseif aridity(k) < 1.50
                aLong = "subhumid";
            else
                aLong = "dry";
            end
        else
            aLong = "unknown";
        end

        if isfinite(frac_snow(k)) ...
                && frac_snow(k) >= 0.15
            sLong = "snow";
        else
            sLong = "rain";
        end

        id(k) = zone_title_CH( ...
            aLong + "_" + sLong);
        long(k) = aLong + " / " + sLong;
    end

    codes = unique(id,'stable');
    names = strings(size(codes));
    num = nan(K,1);

    for j = 1:numel(codes)
        tf = id == codes(j);
        num(tf) = j;
        jj = find(tf,1,'first');
        names(j) = zone_title_CH(long(jj));
    end

    zone.id = id;
    zone.long = long;
    zone.num = num;
    zone.codes = codes;
    zone.names = names;
    zone.lat = lat;
    zone.lon = lon;
    zone.area = area;
    zone.aridity = aridity;
    zone.frac_snow = frac_snow;

end

% ==========================================
function zone = classify_camels_cl_zones(AA)
% ==========================================

    K = height(AA);
    
    lat = get_numeric_column(AA, ...
        ["gauge_lat","lat","lat_outlet"]);
    lon = get_numeric_column(AA, ...
        ["gauge_lon","lon","long","lon_outlet"]);
    area = get_numeric_column(AA, ...
        ["area","catchment_area","area_gages2"]);
    aridity = get_numeric_column(AA, ...
        ["aridity","aridity_index"]);
    frac_snow = get_numeric_column(AA, ...
        ["frac_snow","snow_frac","snow_fraction"]);
    
    id = strings(K,1);
    long = strings(K,1);
    
    for k = 1:K
        if isfinite(aridity(k))
            if aridity(k) < 0.80
                aLong = "humid";
            elseif aridity(k) < 1.50
                aLong = "subhumid";
            else
                aLong = "dry";
            end
        else
            % Use latitude fallback when aridity is unavailable.
            if isfinite(lat(k)) ...
                    && lat(k) <= -40
                aLong = "southern";
            elseif isfinite(lat(k)) ...
                    && lat(k) <= -30
                aLong = "central";
            elseif isfinite(lat(k))
                aLong = "northern";
            else
                aLong = "unknown";
            end
        end
    
        if isfinite(frac_snow(k)) ...
                && frac_snow(k) >= 0.05
            sLong = "snow";
        else
            sLong = "rain";
        end
    
        id(k) = zone_title_CL(aLong + "_" + sLong);
        long(k) = aLong + " / " + sLong;
    end
    
    codes = unique(id,'stable');
    names = strings(size(codes));
    num = nan(K,1);
    
    for j = 1:numel(codes)
        tf = id == codes(j);
        num(tf) = j;
        jj = find(tf,1,'first');
        names(j) = zone_title_CL(long(jj));
    end
    
    zone.id = id;
    zone.long = long;
    zone.num = num;
    zone.codes = codes;
    zone.names = names;
    zone.lat = lat;
    zone.lon = lon;
    zone.area = area;
    zone.aridity = aridity;
    zone.frac_snow = frac_snow;

end

% ===========================================
function zone = classify_camels_col_zones(AA)
% ===========================================

    zone = struct();
    AI = double(AA.aridity);
    Z = double(AA.mean_ele);
    z = strings(height(AA),1);

    for i = 1:height(AA)
        if AI(i) < 0.65
            clim = "dry";
        elseif AI(i) < 1.5
            clim = "humid";
        else
            clim = "wet";
        end
        if Z(i) < 1000
            topo = "lowland";
        elseif Z(i) < 2500
            topo = "montane";
        else
            topo = "highland";
        end
        z(i) = clim + "_" + topo;
    end

    zone.id = z(:);
    zone.name = z(:);
    zone.label = z(:);
    zone.group = z(:);

    if ismember('gauge_lat', ...
            AA.Properties.VariableNames)
        zone.lat = double(AA.gauge_lat(:));
    else
        zone.lat = nan(height(AA),1);
    end
    if ismember('gauge_lon', ...
            AA.Properties.VariableNames)
        zone.lon = double(AA.gauge_lon(:));
    else
        zone.lon = nan(height(AA),1);
    end
    if ismember('area', ...
            AA.Properties.VariableNames)
        zone.area = double(AA.area(:));
    else
        zone.area = nan(height(AA),1);
    end
end

% ==========================================
function zone = classify_camels_de_zones(AA)
% ==========================================

    K = height(AA);

    lat = get_numeric_column(AA, ...
        "gauge_lat");
    lon = get_numeric_column(AA, ...
        "gauge_lon");
    area = get_numeric_column(AA, ...
        "area");
    elev_mean = get_numeric_column(AA, ...
        "elev_mean");
    p_mean = get_numeric_column(AA, ...
        "p_mean");
    frac_snow = get_numeric_column(AA, ...
        "frac_snow");
    forest_perc = get_numeric_column(AA, ...
        "forests_and_seminatural_areas_perc");

    id = strings(K,1);
    long = strings(K,1);

    % 3 precipitation classes x 2 snow classes = 6 zones
    qP = quantile(p_mean(isfinite(p_mean)), ...
        [1/3 2/3]);

    for k = 1:K

        if ~isfinite(p_mean(k))
            pClass = "unknown";
        elseif p_mean(k) <= qP(1)
            pClass = "dry";
        elseif p_mean(k) <= qP(2)
            pClass = "subhumid";
        else
            pClass = "humid";
        end

        if isfinite(frac_snow(k)) ...
                && frac_snow(k) >= 0.05
            sClass = "mixed";
        else
            sClass = "rain";
        end

        id(k) = pClass + "_" + sClass;
        long(k) = pClass + " " + sClass;
    end

    codes = unique(id,'stable');
    names = strings(size(codes));
    num = nan(K,1);

    for j = 1:numel(codes)
        tf = id == codes(j);
        num(tf) = j;
        jj = find(tf,1,'first');
        names(j) = zone_title_DE(long(jj));
    end

    zone.id = id;
    zone.long = long;
    zone.num = num;
    zone.codes = codes;
    zone.names = names;
    zone.lat = lat;
    zone.lon = lon;
    zone.area = area;
    zone.elev_mean = elev_mean;
    zone.p_mean = p_mean;
    zone.frac_snow = frac_snow;
    zone.forest_perc = forest_perc;

end

% ==========================================
function zone = classify_camels_es_zones(AA)
% ==========================================
    zone = struct();

    if ismember('aridity', ...
            AA.Properties.VariableNames)
        AI = double(AA.aridity);
    elseif ismember('pet_mean', ...
            AA.Properties.VariableNames) ...
            && ismember('p_mean', ...
            AA.Properties.VariableNames)
        AI = double(AA.pet_mean) ...
            ./ max(double(AA.p_mean),eps);
    else
        AI = nan(height(AA),1);
    end

    if ismember('elv_mean', ...
            AA.Properties.VariableNames)
        Z = double(AA.elv_mean);
    elseif ismember('ele_mt_sav', ...
            AA.Properties.VariableNames)
        Z = double(AA.ele_mt_sav);
    else
        Z = nan(height(AA),1);
    end

    z = strings(height(AA),1);
    for i = 1:height(AA)
        if isfinite(AI(i)) ...
                && AI(i) > 1.5
            clim = "dry";
        elseif isfinite(AI(i)) ...
                && AI(i) > 0.8
            clim = "transitional";
        else
            clim = "humid";
        end

        if isfinite(Z(i)) ...
                && Z(i) >= 1000
            topo = "mountain";
        elseif isfinite(Z(i)) ...
                && Z(i) >= 400
            topo = "upland";
        else
            topo = "lowland";
        end
        z(i) = clim + "_" + topo;
    end

    zone.id = z(:);
    zone.name = z(:);
    zone.label = z(:);
    zone.group = z(:);
    zone.lat = double(AA.gauge_lat(:));
    zone.lon = double(AA.gauge_lon(:));
    zone.area = double(AA.area(:));
end

% ==========================================
function zone = classify_camels_fi_zones(AA)
% ==========================================
    K = height(AA);

    aridity = get_numeric_column(AA, ...
        {"aridity"});
    frac_snow = get_numeric_column(AA, ...
        {"frac_snow"});
    lat = get_numeric_column(AA, ...
        {"gauge_lat"});
    lon = get_numeric_column(AA, ...
        {"gauge_lon"});
    area = get_numeric_column(AA, ...
        {"area"});
    elev = get_numeric_column(AA, ...
        {"elev_mean","elev_gauge"});

    z = strings(K,1);
    for k = 1:K
        if isfinite(frac_snow(k)) ...
                && frac_snow(k) >= 0.35
            snow = "snow_dominated";
        elseif isfinite(frac_snow(k)) ...
                && frac_snow(k) >= 0.20
            snow = "mixed_snow_rain";
        else
            snow = "rain_dominated";
        end

        if isfinite(lat(k)) ...
                && lat(k) >= 66
            region = "north";
        elseif isfinite(lat(k)) ...
                && lat(k) <= 62
            region = "south";
        else
            region = "central";
        end    
        z(k) = region + "_" + snow;
    end

    names = unique(z,'stable');
    num = nan(K,1);
    for j = 1:numel(names)
        num(z == names(j)) = j;
    end

    zone = struct();
    zone.id = z(:);
    zone.long = z(:);
    zone.num = num(:);
    zone.codes = names(:);
    zone.names = replace(names(:),"_"," ");
    zone.aridity = aridity(:);
    zone.frac_snow = frac_snow(:);
    zone.lat = lat(:);
    zone.lon = lon(:);
    zone.area = area(:);
    zone.elev_mean = elev(:);
    zone.p_mean = get_numeric_column( ...
        AA,{"p_mean"});
    zone.pet_mean = get_numeric_column( ...
        AA,{"pet_mean"});
end

% ==========================================
function zone = classify_camels_fr_zones(AA)
% ==========================================
    K = height(AA);
    [lat,lon] = local_fr_coordinates(AA);
    elev = get_numeric_column(AA,{"sit_altitude", ...
        "top_altitude_mean"});

    z = strings(K,1);
    for i = 1:K
        if isfinite(elev(i)) ...
                && elev(i) >= 800
            z(i) = "mountain";
        elseif isfinite(lat(i)) ...
                && lat(i) >= 48.5
            z(i) = "northern_france";
        elseif isfinite(lat(i)) ...
                && lat(i) <= 45.0
            z(i) = "southern_france";
        else
            z(i) = "central_france";
        end
    end

    [names,~,num] = unique(z);
    zone = struct();
    zone.id = z(:);
    zone.long = z(:);
    zone.num = num(:);
    zone.codes = names(:);
    zone.names = replace(names(:),"_"," ");
    zone.lat = lat(:);
    zone.lon = lon(:);
    zone.area = get_numeric_column( ...
        AA,{"sit_area_topo", ...
        "sit_area_hydro"});
    zone.elev_mean = elev(:);
    zone.p_mean = nan(K,1);
    zone.pet_mean = nan(K,1);
    zone.aridity = nan(K,1);
    zone.frac_snow = nan(K,1);
end

% ==========================================
function zone = classify_camels_gb_zones(AA)
% ==========================================
    K = height(AA);

    aridity = nan(K,1);
    frac_snow = nan(K,1);

    if ismember('aridity', ...
            AA.Properties.VariableNames)
        aridity = double(AA.aridity);
    end
    if ismember('frac_snow', ...
            AA.Properties.VariableNames)
        frac_snow = double(AA.frac_snow);
    end

    id = strings(K,1);
    for k = 1:K
        if ~isfinite(aridity(k))
            moist = "unknown";
        elseif aridity(k) < 0.75
            moist = "humid";
        elseif aridity(k) < 1.25
            moist = "subhumid";
        else
            moist = "dry";
        end

        if ~isfinite(frac_snow(k))
            snow = "unknown";
        elseif frac_snow(k) >= 0.10
            snow = "snow";
        else
            snow = "rain";
        end

        id(k) = moist + "_" + snow;
    end

    names = unique(id,'stable');
    num = nan(K,1);
    for j = 1:numel(names)
        num(id == names(j)) = j;
    end

    zone.id = id;
    zone.num = num;
    zone.names = names;
    zone.aridity = aridity;
    zone.frac_snow = frac_snow;
end

% ===========================================
function zone = classify_camels_ind_zones(AA)
% ===========================================

    K = height(AA);

    if ismember('river_basin', ...
            AA.Properties.VariableNames)
        zlong = string(AA.river_basin);
    else
        zlong = repmat("India",K,1);
    end

    % zlong = strip(zlong);
    % zlong(ismissing(zlong) ...
    %     | zlong == "") = "India";
    % 
    % [names,~,num] = unique(zlong,'stable'); 
    % --> leads to 15 zones - based on river basins
    zlong = strip(lower(zlong));
    zlong(ismissing(zlong) ...
        | zlong == "") = "india";

    zgroup = strings(size(zlong));

    for i = 1:numel(zlong)
        switch char(zlong(i))

            case {'wfrn','wfrs'}
                zgroup(i) = "west_coast";

            case {'efrn','efrs'}
                zgroup(i) = "east_coast";

            case {'cauvery','krishna','pennar'}
                zgroup(i) = "peninsular_south";

            case {'godavari','mahanadi','subernar'}
                zgroup(i) = "central_east";

            case {'narmada','tapi','mahi','sabarmat'}
                zgroup(i) = "central_west";

            case {'bb______'}
                zgroup(i) = "northeast";

            otherwise
                zgroup(i) = "other";
        end
    end

    [names,~,num] = unique(zgroup,'stable');

    codes = replace(names,"_"," ");

    zid = codes(num);

    zone = struct();
    zone.id = zid(:);
    zone.long = zgroup(:);
    zone.river_basin = zlong(:);
    zone.num = num(:);
    zone.codes = codes(:);
    zone.names = names(:);

    zone.lat = get_numeric_column(AA, ...
        ["cwc_lat","ghi_lat"]);
    zone.lon = get_numeric_column(AA, ...
        ["cwc_lon","ghi_lon"]);
    zone.area = get_numeric_column(AA, ...
        ["cwc_area","ghi_area"]);
    zone.elev_mean = get_numeric_column(AA, ...
        "elev_mean");
    zone.p_mean = get_numeric_column(AA, ...
        "p_mean");
    zone.pet_mean = get_numeric_column(AA, ...
        "pet_mean");
    zone.aridity = get_numeric_column(AA, ...
        ["aridity_p_pet","ai_mean"]);
    zone.frac_snow = nan(K,1);

    if ismember('flow_availability', ...
            AA.Properties.VariableNames)
        zone.flow_availability = ...
            get_numeric_attribute( ...
            AA.flow_availability, ...
            'flow_availability');
    end

    if ismember('cwc_river', ...
            AA.Properties.VariableNames)
        zone.river = string(AA.cwc_river);
    end

end

% ===========================================
function zone = classify_camels_lux_zones(AA)
% ===========================================
    K = height(AA);

    lat = get_numeric_column(AA, ...
        {"Lat"});
    lon = get_numeric_column(AA, ...
        {"Lon"});
    elev = get_numeric_column(AA, ...
        {"Z_MEAN"});
    perm = get_numeric_column(AA, ...
        {"permeable_formations"});
    area = get_numeric_column(AA, ...
        {"area_km2"});

    z = strings(K,1);
    for i = 1:K
        if isfinite(elev(i)) ...
                && elev(i) >= 400
            z(i) = "northern_uplands";
        elseif isfinite(perm(i)) ...
                && perm(i) >= 60
            z(i) = "permeable_sedimentary";
        elseif isfinite(lat(i)) ...
                && lat(i) < 49.65
            z(i) = "southern_lowlands";
        else
            z(i) = "central_luxembourg";
        end
    end

    names = unique(z,'stable');
    num = nan(K,1);
    for j = 1:numel(names)
        num(z == names(j)) = j;
    end

    zone = struct();
    zone.id = z(:);
    zone.long = z(:);
    zone.num = num(:);
    zone.codes = names(:);
    zone.names = replace(names(:),"_"," ");
    zone.lat = lat(:);
    zone.lon = lon(:);
    zone.area = area(:);
    zone.elev_mean = elev(:);
    zone.permeable = perm(:);
    zone.p_mean = get_numeric_column(AA,{"P_stn_sum"});
    zone.pet_mean = get_numeric_column(AA,{"PET_PM_sum"});
    zone.aridity = get_numeric_column(AA,{"P_PET_PM"});
    zone.frac_snow = nan(K,1);
end

% ==========================================
function zone = classify_camels_nz_zones(AA)
% ==========================================

    K = height(AA);
    lat = get_numeric_column(AA, ...
        "gauge_lat");
    lon = get_numeric_column(AA, ...
        "gauge_lon");
    area = get_numeric_column(AA, ...
        "area");
    p_mean = get_numeric_column(AA, ...
        "p_mean");
    pet_mean = get_numeric_column(AA, ...
        "pet_mean");
    elev = get_numeric_column(AA, ...
        "elev_mean");

    climate = get_string_column(AA, ...
        "climate_zone");
    landcover = get_string_column(AA, ...
        "landcover");
    geology = get_string_column(AA, ...
        "geology");
    influence = get_string_column(AA, ...
        "influence");

    snow = false(K,1);
    if ismember('is_snow_influenced', ...
            AA.Properties.VariableNames)
        snow = logical(AA.is_snow_influenced);
    end

    id = strings(K,1);
    long = strings(K,1);
    for k = 1:K
        c = climate(k);
        if ismissing(c) ...
                || strlength(strtrim(c)) == 0
            c = "Unknown";
        end
        if snow(k)
            s = "Snow";
        else
            s = "Rain";
        end
        id(k) = zone_title_NZ(c + "_" + s);
        long(k) = c + " / " + s;
    end

    codes = unique(id,'stable');
    names = strings(size(codes));
    num = nan(K,1);
    for j = 1:numel(codes)
        tf = id == codes(j);
        num(tf) = j;
        jj = find(tf,1,'first');
        names(j) = zone_title_NZ(long(jj));
    end

    zone.id = id;
    zone.long = long;
    zone.num = num;
    zone.codes = codes;
    zone.names = names;
    zone.lat = lat;
    zone.lon = lon;
    zone.area = area;
    zone.p_mean = p_mean;
    zone.pet_mean = pet_mean;
    zone.aridity = pet_mean ./ p_mean;
    zone.elev = elev;
    zone.climate_zone = climate;
    zone.landcover = landcover;
    zone.geology = geology;
    zone.influence = influence;
    zone.is_snow_influenced = snow;

end

% ===========================================
function zone = classify_camels_pl_zones(AA)
% ==========================================
    zone = struct();

    P = nan(height(AA),1);
    Z = nan(height(AA),1);

    if ismember('p_mean', ...
            AA.Properties.VariableNames)
        P = double(AA.p_mean(:));
    end
    if ismember('elev_mean', ...
            AA.Properties.VariableNames)
        Z = double(AA.elev_mean(:));
    elseif ismember('gauge_elev', ...
            AA.Properties.VariableNames)
        Z = double(AA.gauge_elev(:));
    end

    z = strings(height(AA),1);
    for i = 1:height(AA)
        if ~isfinite(P(i))
            clim = "unknown";
        elseif P(i) < 1.6
            clim = "dry";
        elseif P(i) < 2.3
            clim = "moderate";
        else
            clim = "wet";
        end

        if ~isfinite(Z(i))
            topo = "unknown";
        elseif Z(i) < 200
            topo = "lowland";
        elseif Z(i) < 500
            topo = "upland";
        else
            topo = "mountain";
        end
        z(i) = clim + "_" + topo;
    end

    zone.id = z(:);
    zone.name = z(:);
    zone.label = z(:);
    zone.group = z(:);
    zone.lat = double(AA.gauge_lat_dd(:));
    zone.lon = double(AA.gauge_lon_dd(:));
    zone.area = double(AA.area_km2(:));
end

% ==========================================
function zone = classify_camels_se_zones(AA)
% ==========================================

K = height(AA);
zone = struct();
zone.id = strings(K,1);
zone.num = nan(K,1);
zone.names = strings(0,1);
zone.aridity = nan(K,1);
zone.frac_snow = nan(K,1);

if ismember('Tmean_C', ...
        AA.Properties.VariableNames)
    Tm = double(AA.Tmean_C(:));
else
    Tm = nan(K,1);
end

if ismember('Elevation_mabsl', ...
        AA.Properties.VariableNames)
    elev = double(AA.Elevation_mabsl(:));
else
    elev = nan(K,1);
end

% Approximate snow tendency for plotting/grouping only. The CAMELS-SE
% attribute files do not provide a snow-fraction attribute directly.
fs = max(0,min(1,(2 - Tm)/8));
fs(~isfinite(fs)) = nan;
zone.frac_snow = fs;

for k = 1:K
    if isfinite(Tm(k)) ...
            && Tm(k) < 0
        zid = "cold snow dominated";
    elseif isfinite(Tm(k)) ...
            && Tm(k) < 3
        zid = "cool mixed snow rain";
    elseif isfinite(elev(k)) ...
            && elev(k) > 700
        zid = "upland cool catchment";
    else
        zid = "temperate rain dominated";
    end
    zone.id(k) = zid;
end

zone.names = unique(zone.id,'stable');
for k = 1:K
    zone.num(k) = ...
        find(zone.names == zone.id(k),1,'first');
end

end

% ===========================================
function [lat,lon] = local_fr_coordinates(AA)
% ===========================================

    lat = get_numeric_column(AA, ...
        {"sit_latitude","sta_y_w84", ...
        "sta_y_w84_snap"});

    lon = get_numeric_column(AA, ...
        {"sit_longitude","sta_x_w84", ...
        "sta_x_w84_snap"});

    bad = lat < 40 ...
        | lat > 52 ...
        | lon < -6 ...
        | lon > 10;

    if any(bad)
        lat2 = get_numeric_column(AA, ...
            {"sta_y_w84","sta_y_w84_snap", ...
            "sit_latitude"});

        lon2 = get_numeric_column(AA, ...
            {"sta_x_w84","sta_x_w84_snap", ...
            "sit_longitude"});

        ok2 = lat2 >= 40 ...
            & lat2 <= 52 ...
            & lon2 >= -6 ...
            & lon2 <= 10;

        lat(bad & ok2) = lat2(bad & ok2);
        lon(bad & ok2) = lon2(bad & ok2);
    end

end

% ===========================
function s = zone_title_AU(s)
% ===========================

    s = replace(string(s),' / ',' ');
    s = replace(s,'_',' ');
    s = lower(strtrim(s));
    
    if strlength(s) > 0
        w = char(s);
        w(1) = upper(w(1));
        s = string(w);
    end

end

% =============================
function out = zone_title_CH(s)
% =============================

    s = string(s);
    s = replace(s,"/","_");
    s = replace(s," ","_");
    s = lower(s);
    parts = split(s,"_");
    parts(parts == "") = [];
    if isempty(parts)
        out = "unknown";
        return
    end
    for i = 1:numel(parts)
        p = char(parts(i));
        parts(i) = ...
            upper(extractBefore(string(p),2)) ...
            + extractAfter(string(p),1);
    end
    out = strjoin(parts,"_");

end

% =============================
function out = zone_title_NZ(s)
% =============================

    s = string(s);
    s = replace(s,"/","_");
    s = replace(s," ","_");
    s = regexprep(s,'[^A-Za-z0-9_]+','');
    s = lower(s);
    parts = split(s,"_");
    parts(parts == "") = [];
    if isempty(parts)
        out = "unknown";
        return
    end
    for i = 1:numel(parts)
        p = char(parts(i));
        parts(i) = upper(extractBefore( ...
            string(p),2)) ...
            + extractAfter(string(p),1);
    end
    out = strjoin(parts,"_");

end

% =============================
function out = zone_title_DE(s)
% =============================

    s = string(s);
    s = replace(s,"_"," ");
    s = replace(s,"/"," ");
    s = regexprep(s,'\s+',' ');
    s = strip(lower(s));

    parts = split(s);
    for k = 1:numel(parts)
        if strlength(parts(k)) > 0
            parts(k) = upper( ...
                extractBefore(parts(k),2)) + ...
                extractAfter(parts(k),1);
        end
    end
    out = strjoin(parts," ");

end

% =============================
function out = zone_title_CL(s)
% =============================

    s = string(s);
    s = replace(s,"_"," ");
    s = replace(s,"/"," ");
    s = regexprep(s,'\s+',' ');
    s = strip(lower(s));

    parts = split(s);
    for k = 1:numel(parts)
        if strlength(parts(k)) > 0
            parts(k) = upper(extractBefore(parts(k),2)) + ...
                extractAfter(parts(k),1);
        end
    end
    out = strjoin(parts," ");

end

% ===========================
function s = zone_title_BR(s)
% ===========================

    s = replace(string(s),' / ',' ');
    s = lower(strtrim(s));

    if strlength(s) > 0
        w = char(s);
        w(1) = upper(w(1));
        s = string(w);
    end

end

% ================================
function out = local_zone_title(s)
% ================================

    s = string(s);
    s = replace(s,"/","_");
    s = replace(s," ","_");
    s = regexprep(s,'[^A-Za-z0-9_]+','');
    s = lower(strip(s));
    parts = split(s,"_");
    parts(parts == "") = [];
    if isempty(parts)
        out = "Unknown";
        return
    end
    for ii = 1:numel(parts)
        p = char(parts(ii));
        parts(ii) = upper(extractBefore(string(p),2)) + ...
            extractAfter(string(p),1);
    end
    out = strjoin(parts,"_");

end
