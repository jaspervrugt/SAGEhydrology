function ET0 = et0_fao56_daily(TmaxC,TminC,Rs_Wm2,Vp_Pa,DOY,lat_deg, ...
    z_m,u2,method)
%ET0_FAO56_DAILY Calculates daily ET0 using standard metereological
% variables.
% SYNOPSIS: ET0 = et0_fao56_daily(TmaxC,TminC,Rs_Wm2,Vp_Pa,DOY,lat_deg, ...
%    z_m,u2,method)
%   TmaxC       Maximum daily temperature in C
%   TminC       Minimum daily temperature in C
%   Rs_Wm2      Global radition in W/m2
%   Vp_Pa       Vapor pressure in Pascal
%   DOY         Day of Year
%   lat_deg     Latitude in degrees
%   z_m         Elevation in meters
%   u2          Windspeed at 2 m height in m/s
%   method      Calculation method for ET0
%    [1] Penman-Monteith
%    [2] Priestley Taylor
%    [3] Makkink
%   ET0         ET0 in mm/day
%
% FAO-56 daily reference ET0 (mm/day)
% Inputs can be vectors (n x 1). lat_deg,z_m,u2 can be scalar or n x 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tmean = 0.5*(TmaxC+TminC);          % Daily mean temperature in C
Rs = Rs_Wm2 * 0.0864;               % Radiation W/m2 -> MJ/m2/day
ea = Vp_Pa / 1000;                  % Actual vapor pressure (kPa)
es_Tmax = 0.6108 .* exp(17.27* ...  % Saturation vapor pressure (kPa)
    TmaxC./(TmaxC+237.3));
es_Tmin = 0.6108 .* exp(17.27* ...  % Saturation vapor pressure (kPa)
    TminC./(TminC+237.3));
es = 0.5*(es_Tmax + es_Tmin);       % Saturation vapor pressure (kPa)
Delta = 4098*(0.6108*exp( ...       % Slope of saturation vapor pressure curve (kPa/°C)
    17.27*Tmean./(Tmean + ...
    237.3))) ./ (Tmean+237.3).^2;
P = 101.3 * ((293 - 0.0065* ...     % Pressure from elevation (kPa)
    z_m)/293).^5.26;
gamma = 0.000665 * P;               % Psychrometric constant (kPa/°C)

%% Extraterrestrial radiation Ra (MJ/m2/day) for PM and PT
phi = lat_deg*pi/180;
dr  = 1 + 0.033*cos(2*pi*DOY/365);
delta = 0.409*sin(2*pi*DOY/365 - 1.39);
ws = acos(-tan(phi).*tan(delta));
Gsc = 0.0820;                       % MJ/m2/min
Ra = (24*60/pi) * Gsc .* dr .* ...  % Global radiation? W/m2
    (ws.*sin(phi).*sin(delta) + ...
    cos(phi).*cos(delta).*sin(ws));
Rso = (0.75 + 2e-5*z_m) .* Ra;      % Clear-sky radiation Rso (MJ/m2/day)
albedo = 0.23;                      % Albedo
Rns = (1-albedo).*Rs;               % Net shortwave
sigma = 4.903e-9;                   % MJ/K^4/m2/day
TmaxK = TmaxC + 273.16;             % Minimum Temperature in Kelvin
TminK = TminC + 273.16;             % Maximum Temperature in Kelvin
Rs_Rso = min(1.0, max(0.0, ...      % clamp [0,1]
    Rs./max(Rso,eps)));
Rnl = sigma .* ((TmaxK.^4 + ...     % Net longwave radiation in W/m2
    TminK.^4)/2) .* (0.34 - ...
    0.14*sqrt(max(ea,0))).* ...
    (1.35*Rs_Rso - 0.35);
Rn = Rns - Rnl;                     % Net radiation in W/m2
G = 0;                              % Soil heat flux (daily) in W/m2

switch method
    case 1  % Penman–Monteith (FAO-56)
        ET0 = (0.408*Delta.*(Rn-G) + gamma.*(900./(Tmean+273)).* ...
            u2.*(es-ea)) ./ (Delta + gamma.*(1 + 0.34*u2));
    case 2  % Priestley–Taylor
        alphaPT = 1.26;
        ET0 = alphaPT .* (Delta./(Delta + gamma)) .* (0.408*(Rn - G));
    case 3  % Makkink: ET0 = 0.61*(Delta/(Delta+gamma))*(Rs/lambda) - 0.12
        c = 0.61;
        b = -0.12; % mm/day
        ET0 = c*(Delta./(Delta + gamma)) .* (0.408*Rs) + b;
    otherwise
        error(['      Error:read_meteo:' ...
            'et0_fao56_daily:Unknown method. ' ...
            'Use 1 = PM, 2 = Priestley-Taylor, ' ...
            '3 = Makkink.']);
end

ET0 = max(0,ET0);       % Clean negatives (can occur in cold/dark days)

end
