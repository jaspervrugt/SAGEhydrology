function ET0 = et0_lamah_hourly(T,Tdew,u10,v10,Rns,Rnt,Pa,method)
%ET0_LAMAH_HOURLY Calculates hourly reference ET0 from LamaH-CE variables
%
% SYNOPSIS: ET0 = et0_lamah_hourly(T_C,Tdew,u10,v10,Rns,Rnt,Pa,method)
%   T          2 m air temperature (°C)
%   Tdew       2 m dew-point temperature (°C)
%   u10        10 m eastward wind component (m/s)
%   v10        10 m northward wind component (m/s)
%   Rns        net shortwave radiation at surface (W/m2), positive downward
%   Rnt        net thermal radiation at surface (W/m2), positive upward
%   Pa         surface pressure (Pa)
%   method     PET method
%                1 = FAO-56 Penman-Monteith
%                2 = Priestley-Taylor
%                3 = Makkink
%   ET0        hourly reference evapotranspiration (mm/h)
%
% Notes:
%   LamaH-CE hourly timestamps are left-edge interval timestamps. A value
%   at 00:00 represents the mean or sum from 00:00 to 01:00 UTC.
%   Radiation is converted from W/m2 to MJ/m2/h using 0.0036.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Jun. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 8 ...
            || isempty(method)
        method = 1;
    end

    T = double(T(:));
    Tdew = double(Tdew(:));
    u10 = double(u10(:));
    v10 = double(v10(:));
    Rns = double(Rns(:));
    Rnt = double(Rnt(:));
    Pa = double(Pa(:));

    n = max([numel(T),numel(Tdew), ...
        numel(u10),numel(v10), ...
        numel(Rns),numel(Rnt),numel(Pa)]);

    T = local_expand(T,n);
    Tdew = local_expand(Tdew,n);
    u10 = local_expand(u10,n);
    v10 = local_expand(v10,n);
    Rns = local_expand(Rns,n);
    Rnt = local_expand(Rnt,n);
    Pa = local_expand(Pa,n);

    Tk = T + 273.15;
    p_kPa = Pa ./ 1000;
    gamma = 0.000665 .* p_kPa;

    % Wind speed at 10 m, adjusted to 2 m using FAO-56 correction.
    u10mag = hypot(u10,v10);
    u2 = u10mag .* (4.87 ./ ...
        log(67.8*10 - 5.42));

    % Actual vapor pressure from dew-point temperature.
    ea_kPa = 0.6108 .* exp((17.27 .* Tdew) ./ ...
        (Tdew + 237.3));

    % Saturation vapor pressure and slope at air temperature.
    es_kPa = 0.6108 .* exp((17.27 .* T) ./ ...
        (T + 237.3));
    Delta = 4098 .* es_kPa ./ (T + 237.3).^2;

    % LamaH-CE convention:
    %   surf_net_solar_rad: positive downward/to surface
    %   surf_net_therm_rad: positive upward/from surface
    % Available net radiation is therefore shortwave minus thermal.
    Rn = (Rns - Rnt) .* 0.0036; % MJ/m2/h

    % Hourly soil heat flux approximation.
    isDay = Rns > 0;
    G = zeros(size(Rn));
    G(isDay) = 0.10 .* Rn(isDay);
    G(~isDay) = 0.50 .* Rn(~isDay);

    switch double(method)
        case 1
            VPD = max(es_kPa - ea_kPa,0);
            ET0 = (0.408 .* Delta .* (Rn - G) + ...
                gamma .* (37 ./ Tk) .* u2 .* VPD) ./ ...
                (Delta + gamma .* (1 + 0.34 .* u2));

        case 2
            alphaPT = 1.26;
            ET0 = alphaPT .* (Delta ./ (Delta + gamma)) .* ...
                (0.408 .* (Rn - G));

        case 3
            c = 0.61;
            b = -0.12 / 24;
            Rs_MJph = max(Rns,0) .* 0.0036;
            ET0 = c .* (Delta ./ (Delta + gamma)) .* ...
                (0.408 .* Rs_MJph) + b;

        otherwise
            error(['      Error:et0_lamah_hourly: ' ...
                'Unknown method. ' ...
                'Use 1 = Penman-Monteith, ' ...
                '2 = Priestley-Taylor, ' ...
                'or 3 = Makkink.']);
    end

    ET0 = max(ET0,0);
    ET0(~isfinite(ET0)) = NaN;
end

% ============================
function x = local_expand(x,n)
% ============================
    if numel(x) == n
        return
    elseif isscalar(x)
        x = repmat(x,n,1);
    else
        error(['      Error:et0_lamah_hourly: ' ...
            'Input sizes are ' ...
            'not compatible.']);
    end
end
