function ET0 = et0_fao56_hourly(Rld,p,Rs,q,T,u10,v10,method)
%ET0_FAO56_HOURLY Calculates hourly ET0 using standard metereological
% variables.
% SYNOPSIS: ET0 = et0_fao56_hourly(Rld,p,Rs,q,T,u10,v10)
%   Rld     Longwave_radiation (W/m2)   [downwelling LW]
%   p       Air pressure (Pa)
%   Rs      Shortwave_radiation (W/m2)  [downwelling SW]
%   q       Specific_humidity (kg/kg)
%   T       Temperature (°C)     
%   u       Horizontal windspeed (m/s)  [assumed at 10 m]
%   v       Vertical windspeed (m/s)    [assumed at 10 m]
%   method  Calculation method for ET0
%    [1] Penman-Monteith
%    [2] Priestley Taylor
%    [3] Makkink
%   ET0     OUTPUT: hourly reference ET0 (mm/h) FAO-56 Penman–Monteith
%
% Notes:
%  - This uses net radiation estimated from downwelling SW/LW and surface
%    thermal emission eps*sigma*T^4 (no cloud/LAI details).
%  - Ground heat flux G is parameterized as a fraction of Rn (day/night).
%  - Wind speed is adjusted from 10 m to 2 m (FAO-56).
%  - Vapor pressures computed from q and p.
% 
% FAO-56 hourly reference ET0 (mm/h)
% Inputs can be vectors (n x 1). lat_deg,z_m,u2 can be scalar or n x 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If your temperature is actually in Kelvin, uncomment:
% T = T - 273.15;

Tk = T + 273.15;            % K
p_kPa = p / 1000;           % kPa
gamma = 0.000665 .* p_kPa;  % Psychrometric constant, kPa/°C
z = 10;                     % m
sigma = 5.670374419e-8;     % W/m2/K4
eps_sfc = 0.98;             % reference grass surface emissivity (approx.)
albedo  = 0.23;             % albedo of Earth surface

% Wind speed at 10 m and adjust to 2 m
u10mag = hypot(u10,v10);  % sqrt(u^2+v^2)
% FAO-56 wind adjustment: u2 = u_z * 4.87 / ln(67.8*z - 5.42)
u2 = u10mag .* (4.87 ./ log(67.8*z - 5.42));  % ~0.748 * u10

% Actual vapor pressure from specific humidity
% e_a in Pa (exact relation)
ea_Pa = (q .* p) ./ (0.622 + 0.378*q);
ea_kPa = ea_Pa / 1000;

% Saturation vapor pressure and slope of curve
% FAO-56 saturation vapor pressure es(T) in kPa
es_kPa = 0.6108 .* exp((17.27 .* T) ./ (T + 237.3));
% Slope Delta in kPa/°C:
Delta = 4098 .* es_kPa ./ (T + 237.3).^2;

% Net radiation Rn (MJ/m2/h)
% Convert W/m2 -> MJ/m2/h : multiply by 0.0036
Rs_MJph = Rs * 0.0036;
Rld_MJph = Rld * 0.0036;

% Net shortwave: (1 - albedo)*Rs
Rns = (1 - albedo) .* Rs_MJph;

% Estimate net longwave using downwelling minus upwelling thermal emission
% Rlu = eps*sigma*T^4  (W/m2) -> MJ/m2/h
Rlu_Wm2 = eps_sfc .* sigma .* (Tk.^4);
Rlu_MJph = Rlu_Wm2 * 0.0036;
% Net longwave (positive downward): Rnl = Rld - Rlu
Rnl = Rld_MJph - Rlu_MJph;
% Net radiation:
Rn = Rns + Rnl;  % MJ/m2/h

% Ground heat flux G (MJ/m2/h)
% Simple FAO-style hourly parameterization:
%  day: G ~ 0.1 Rn ; night: G ~ 0.5 Rn
isDay = Rs > 0;     % crude proxy (sun up if SW > 0)
G = zeros(size(Rn));
G(isDay) = 0.10 * Rn(isDay);
G(~isDay) = 0.50 * Rn(~isDay);

switch method
    case 1  % Penman–Monteith (FAO-56 hourly)
            % Hourly form uses constant 37 (instead of 900 for daily)
        VPD = max(es_kPa - ea_kPa,0);
        ET0 = (0.408 .* Delta .* (Rn - G) + gamma .* ...
            (37./Tk) .* u2 .* VPD) ./ ...
            (Delta + gamma .* (1 + 0.34 .* u2));
    case 2  % Priestley–Taylor (hourly adaptation)
        alphaPT = 1.26;  % often tuned; may differ at hourly scale
        ET0 = alphaPT .* (Delta./(Delta + gamma)) .* (0.408 .* (Rn - G));
    case 3  % Makkink (hourly adaptation; coeffs usually daily-calibrated)
        c = 0.61;
        b_day = -0.12;          % mm/day (classic)
        b = b_day / 24;         % convert to mm/h (rough)
        ET0 = c .* (Delta./(Delta + gamma)) .* (0.408 .* Rs_MJph) + b;
    otherwise
        error(['      Error:read_meteo:' ...
            'et0_fao56_hourly: Unknown method. ' ...
            'Use 1=PM, 2=Priestley-Taylor, 3=Makkink.']);
end
ET0 = max(ET0,0);               % Physical clipping (optional)

end