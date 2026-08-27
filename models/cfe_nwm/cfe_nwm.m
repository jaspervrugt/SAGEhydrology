function varargout = cfe_nwm(par,mdl,data,ode,check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CFE_NWM: Runge Kutta implementation of cfe_nwm conceptual watershed model
% SYNOPSIS: varargout = cfe_nwm(par,mdl,data,ode,check)
%  par          dx1 vector of parameter values
%   s_max:par(1) maximum soil storage (mm)
%   s_fc:par(2)  field-capacity storage threshold (-)
%   s_wp:par(3)  wilting-point storage (-)
%   k_sch:par(4) Schaake adjusted magic constant (-)
%   a1:par(5)    soil primary outlet exponent (percolation)
%   k_perc:par(6) soil percolation coefficient
%   lf_thr:par(7) lateral-flow threshold storage (-)
%   a2:par(8)    soil secondary outlet exponent
%   k_lf:par(9)  soil lateral-flow coefficient
%   g_max:par(10) groundwater maximum storage (mm)
%   c_gw:par(11) groundwater discharge coefficient
%   mm:par(12)   groundwater exponent
%   k_nsh:par(13) fast-routing coefficient
%   T_tr:par(14) temperature threshold (°C)
%   f_dd:par(15) degree-day factor (mm/°C/T)
%  mdl          structure with model state/parameter info
%   .mcode       scalar numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation ode_cfe_nwm C++
%   .y0          mx1 vector of initial states
%   .pspace      0:hydrologic, 1:unit cube, 2:unconstrained parameters
%   .th_min      dx1 vector of lower parameter values [= in pspace]
%   .th_max      dx1 vector of upper parameter values [= in pspace]
%   .par_names   1xd cell parameter names
%   .id_train    1x2 vector of start and end index training period
%   .id_eval     1x2 vector of start and end index evaluation period
%   .eval_mode   evaluation design used during SAGE training
%     'per'       training basins evaluated on evaluation period
%     'bas'       evaluation basins evaluated on training period
%     'basper'    evaluation basins evaluated on evaluation period
%     'none'      no evaluation
%   .tout        final model print time (scalar)
%   .idx         1x2 vector of indices train&val periods
%  data         structure with meteorological data and other info
%   .P           (n+m)x1 record of precipitation (mm/T)
%   .Ep          (n+m)x1 record of potential evapotranspiration (mm/T)
%   .T           (n+m)x1 record of air temperature (°C)
%  ode          structure with numerical settings ODE solver
%   .InitStep    Initial time step
%   .MaxStep     Maximum time step
%   .MinStep     Minimum time step
%   .RelTol      Relative tolerance
%   .AbsTol      Absolute tolerance
%   .Order       Order
%   .maxiter     Maximum number of iterations
%   .mem         storage of state variables [0: no, 1: yes]
%  check        numerical check of J(x)_f and J(x)_th matrices (or not)
%   0            do not check
%   1            check Jacobian matrices of states and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    check = 0;          % no check of J(x)_f and J(x)_th
end
mcode = mdl.mcode;      % Formulation/language
                        % 1: Runge Kutta implementation MATLAB
                        % 2: ode45 implementation MATLAB
                        % 3: Explicit Euler int_steps MATLAB
                        % 4: Runge Kutta implementation ode_cfe_nwm C++
mem = ode.mem;          % state variable storage or not
if mcode == 2 && mem == 0
    warning(['cfe_nwm: ' ...
        'built-in ode45 ' ...
        'solver stores ' ...
        'states: mem = 1'])
    mem = 1;
end
T_sm = 1.0;                 % smoothing width (°C) partition/positive-part
eps_m = 1e-6;               % smoothing for min()
K = 3;                      % # Nash stores for routing
giuh = mdl.giuh_ordnts;     % GIUH ordinates
L = numel(giuh);            % GIUH length
dT = 1;                     % Forcing time step Schaake: 1=days, 1/24=hours
m = 3 + L + K + 1;          % # states: [Swe S G Q_1,...,Q_L N_1,...,N_K Q]
ns = mdl.tout + 1;          % # print times
d = numel(par);             % # parameters
n = mdl.idx(2)-mdl.idx(1);  % # elements discharge vector
nvar = m*(d+1);             % # number of variables
rho = 0.01;                 % Dimensionless smoothing coefficient
eps = 5;                    % Dimensionless smoothing coefficient
fail = false;               % Default: model completes run
id = m + (1:d)*m;           % Indices of sensitivity state variables
if mem == 0
    q_n = nan(n,1);
    J = nan(n,d); 
    ipr = mdl.idx(1);       % --> C++ code
else
    Z = nan(ns,nvar); 
    ipr = 0;
end
Z(1,1:m) = mdl.y0;          % Initialize state variables at time 0
Z(1,m+1:nvar) = 0;          % Initialize sensitivity at time 0

switch mdl.pspace
    case 0 % hydrologic parameter values
        th = par;
        if (any(th<mdl.th_min) || any(th>mdl.th_max))
            varargout = {nan(n,1),nan(n,d),nan(d,1),Z}; return
        end
        Jth = ones(d,1);                    % return dq_n/dth
    case 1 % normalized hydrologic parameter values
        nth = par;
        if (any(nth<0) || any(nth>1))
            varargout = {nan(n,1),nan(n,d),nan(d,1),Z}; return
        end
        dth_dnth = mdl.th_max - mdl.th_min; % dth/dnth
        th = mdl.th_min + nth.*dth_dnth;    % hydrologic parameter values
        Jth = dth_dnth;                     % return dq_n/dnth
    case 2 % unconstrained parameters (for training)
        varth = par;
        nth = 1./(1 + exp(-varth));         % normalized parameter values
        dth_dnth = mdl.th_max - mdl.th_min; % dth/dnth
        th = mdl.th_min + nth.*dth_dnth;    % hydrologic parameter values
        dnth_dvarth = nth.*(1-nth);         % dnth/dvarth
        Jth = dth_dnth .* dnth_dvarth;      % return dq_n/dvarth
end

s_max = th(1);                      % maximum soil storage (mm)
s_fc = th(2)*s_max;                 % field-capacity storage threshold (mm)
s_wp = th(3)*s_max;                 % wilting-point storage (mm)
lf_thr = s_wp + th(7)*(s_fc-s_wp);  % lateral-flow threshold storage (mm)

mdl.th = th;

% Initialization
switch mcode
    case {1,2,3}
        %% Unpack parameters
        k_sch = th(4);  % Schaake adjusted magic constant (-)
        a1 = th(5);     % soil primary outlet exponent (percolation)
        k_perc = th(6); % soil percolation coefficient
        a2 = th(8);     % soil secondary outlet exponent
        k_lf = th(9);   % soil lateral-flow coefficient
        g_max = th(10); % groundwater maximum storage (mm)
        c_gw = th(11);  % groundwater discharge coefficient
        mm = th(12);    % groundwater exponent
        k_nsh = th(13); % fast-routing coefficient
        T_tr = th(14);  % Temperature threshold (°C)
        f_dd = th(15);  % Degree-day factor (mm/°C/T)

    case 4
        data.s_max = s_max;     % maximum soil water storage (mm)
        data.s_fc = s_fc;       % field-capacity storage threshold (mm)
        data.s_wp = s_wp;       % wilting-point storage (mm)
        data.k_sch = th(4);     % Schaake adjusted magic constant (-)
        data.a1 = th(5);        % soil primary outlet exp. (percolation)
        data.k_perc = th(6);    % soil percolation coefficient
        data.lf_thr = lf_thr;   % lateral-flow threshold storage (mm)
        data.a2 = th(8);        % soil secondary outlet exponent
        data.k_lf = th(9);      % soil lateral-flow coefficient
        data.g_max = th(10);    % groundwater maximum storage (mm)
        data.c_gw = th(11);     % groundwater discharge coefficient
        data.mm = th(12);       % groundwater exponent
        data.k_nsh = th(13);    % fast-routing coefficient
        data.T_tr = th(14);     % Temperature threshold (°C)
        data.f_dd = th(15);     % Degree-day factor (mm/°C/T)

        data.K = K;             % # Nash-Cascade routing reservoirs
        data.L = L;             % # GIUH states
        data.dT = dT;           % forcing timestep Schaake [1=d, 1/24 = h]
        data.giuh = giuh;       % GIUH ordinates
        data.T_sm = T_sm;       % smoothing width (°C) partition/positive
        data.eps_m = eps_m;     % smoothing for min (°C)
        data.eps = eps;         % Dimensionless smoothing coef. [inactive]
        data.rho = rho;         % Dimensionless smoothing coef. [inactive]
        data.ipr = ipr;         % Time to print

end

% Execute model
switch mcode

    case 1 %% MATLAB: Runge Kutta (= model similar to cfe_nwm_odefcn)
        hin = ode.InitStep;     % Initial time step
        hmax_ = ode.MaxStep;    % Maximum time step
        hmin_ = ode.MinStep;    % minimum time step
        reltol = ode.RelTol;    % Relative tolerance
        abstol = ode.AbsTol;    % Absolute tolerance
        order = ode.Order;      % Order
        maxiter = ode.maxiter;  % Maximum iterations
        iterCount = 0; flag = 0;

        hCarry = max(hmin_,min(hin,hmax_)); % carry adaptive recommendation

        for s = 2:ns            % Integrate from tprint(1) to tprint(end)
            t1 = s-2; t2 = s-1;                 % Set start, end times
            h = min(hCarry,t2-t1);              % reuse prior recommendation
            if mem == 1
                Z(s,1:nvar) = Z(s-1,1:nvar);
                z = Z(s,1:nvar);                    % row vector
            else
                z = Z(1,1:nvar);
            end
            t = t1;                                 % Initial time
            % Integrate from t1 to t2
            while (t < t2)
                [ztmp,LTE] = rk2(t,z,h,s_max,s_fc, ...   % Evaluate rk2
                    s_wp,k_sch,a1,k_perc,lf_thr,a2, ...
                    k_lf,g_max,c_gw,mm,k_nsh,giuh,T_tr,f_dd,K,L, ...
                    dT,data,T_sm,eps_m,eps,rho,m,d);
                if any(~isfinite(ztmp)) || any(abs(ztmp) > 1e12)
                    fail = true; break;
                end
                w = 1 ./ (reltol*abs(ztmp(1:nvar)) + abstol);   % Weights
                wrms = sqrt(sum((w.*LTE).^2)/nvar);             % CORRECTED, Nov. 2022
                if (wrms <= 1) || (h <= hmin_)                  % Accept if error is small enough
                    z = ztmp; t = t + h;
                end
                hNext = h*max(0.2,min(5.0,0.9*wrms^(-1/order)));    % Compute new step
                hNext = max(hmin_,min(hNext,hmax_));
                hCarry = hNext;                         % retain before boundary clipping
                h = min(hNext,t2-t);   % Another turn
                if (iterCount >= maxiter)
                    if (flag == 0)
                        warning("WARNING: Max step limit reached at t = %.5f\n", t);
                        warning("WARNING: Parameter s_max  %8.5f\n", s_max);
                        warning("WARNING: Parameter s_wp   %8.5f\n", s_fc);
                        warning("WARNING: Parameter s_fc   %8.5f\n", s_wp);
                        warning("WARNING: Parameter k_sch  %8.5f\n", k_sch);
                        warning("WARNING: Parameter a1     %8.5f\n", a1);
                        warning("WARNING: Parameter k_perc %8.5f\n", k_perc);
                        warning("WARNING: Parameter lf_thr %8.5f\n", lf_thr);
                        warning("WARNING: Parameter a2     %8.5f\n", a2);
                        warning("WARNING: Parameter k_lf   %8.5f\n", k_lf);
                        warning("WARNING: Parameter g_max  %8.5f\n", g_max);
                        warning("WARNING: Parameter c_gw   %8.5f\n", c_gw);
                        warning("WARNING: Parameter mm     %8.5f\n", mm);
                        warning("WARNING: Parameter k_nsh  %8.5f\n", k_nsh);
                        warning("WARNING: Parameter T_tr    %8.5f\n", T_tr);
                        warning("WARNING: Parameter f_dd   %8.5f\n", f_dd);
                        flag = 1; fail = true; break;
                    end
                end
                if flag == 1
                    fprintf('Time: %7.5f  h: %7.5f  iterCount: %d', ...
                        t,h,iterCount);
                end                
            end
            if fail
                break;
            else
                iterCount = 0;
            end
            if mem == 1
                Z(s,1:nvar) = z;
            else
                if s >= ipr+1
                    q_n(s-ipr) = z(m) - Z(m);
                    J(s-ipr,1:d) = z(id) - Z(id);
                end
                Z(1,1:nvar) = z;
            end
        end

    case 2 %% MATLAB: ode45 implementation of cfe_nwm_odefcn
        hin    = ode.InitStep;   % Initial time step
        hmax_  = ode.MaxStep;    % Maximum time step
        reltol = ode.RelTol;     % Relative tolerance
        abstol = ode.AbsTol;     % Absolute tolerance
        mdl.tout = mdl.tout-1e-10;    % Loop ode goes to maxT and then still one try
        % Then, time index goes out of bound of forcing data
         ode_options = odeset('InitialStep',hin,...  % initial time-step (T)
            'MaxStep',hmax_, ...                    % maximum time-step (T)
            'RelTol',reltol, ...                    % relative tolerance
            'AbsTol',abstol);                       % absolute tol (mm)
        [~,Z] = ode45(@(t,z) cfe_nwm_aug_ode(t,z, ...
            s_max,s_fc,s_wp,k_sch,a1,k_perc,lf_thr, ...
            a2,k_lf,g_max,c_gw,mm,k_nsh,giuh,T_tr, ...
            f_dd,K,L,dT,data,T_sm,eps_m,eps,rho,m,d), ...
            0:mdl.tout,Z(1,1:nvar),ode_options);
        s = size(Z,1); if s < ns, fail = true; s = s+1; end

    case 3 %% MATLAB: Explicit Euler cfe_nwm_odefcn with int_steps steps
        int_steps = 500;                        % number of integration steps
        dt = 1/int_steps;                       % time of each int. step
        for s = 2:ns                            % Start time loop
            if mem == 1
                z = Z(s-1,1:nvar);              % Initialize state variables
            else
                z = Z(1,1:nvar);
            end
            for it = 1:int_steps                % Do integration in int_steps steps
                dzdt = cfe_nwm_aug_ode(s-2,z,s_max, ...
                    s_fc,s_wp,k_sch,a1,k_perc,lf_thr, ...
                    a2,k_lf,g_max,c_gw,mm,k_nsh, ...
                    giuh,T_tr,f_dd,K,L,dT,data,T_sm, ...
                    eps_m,eps,rho,m,d)';
                if s == 150 && check
                    [~,~,Jth_f,Jx_f] = rhs_only(s-2,z,th,K,L,dT, ...
                        data,T_sm,eps_m,eps,rho,m,d);
                    % --- Jth_f via 2-sided finite differences ---
                    Jth_fn = zeros(m,d);
                    for j = 1:d
                        h = 1e-6*max(1,abs(th(j)));
                        thp = th; thm = th;
                        thp(j) = thp(j) + h; thm(j) = thm(j) - h;
                        fp = rhs_only(s-2,z,thp,K,L,dT, ...
                            data,T_sm,eps_m,eps,rho,m,d);
                        fm = rhs_only(s-2,z,thm,K,L,dT, ...
                            data,T_sm,eps_m,eps,rho,m,d);
                        Jth_fn(:,j) = (fp - fm)/(2*h);
                    end
                    % --- Jx_f via 2-sided finite differences ---
                    Jx_fn = zeros(m,m);
                    for i = 1:m
                        h = 1e-6*max(1,abs(z(i)));
                        zp = z; zm = z;
                        zp(i) = zp(i) + h; zm(i) = zm(i) - h;
                        fp = rhs_only(s-2,zp,th,K,L,dT, ...
                            data,T_sm,eps_m,eps,rho,m,d);
                        fm = rhs_only(s-2,zm,th,K,L,dT, ...
                            data,T_sm,eps_m,eps,rho,m,d);
                        Jx_fn(:,i) = (fp - fm)/(2*h);
                    end
                    err_x = Jx_f(1:m,1:m) - Jx_fn(1:m,1:m);
                    err_th = Jth_f(1:m,1:d) - Jth_fn(1:m,1:d);
                    disp(err_th); disp(err_x); pause
                end
                z = z + dzdt * dt;              % Update states
            end
            if any(~isfinite(z)) || any(abs(z) > 1e12)
                fail = true; break;
            else
                % -----------------------------
                if mem == 1
                    Z(s,1:nvar) = z;            % State at t
                else
                    if s >= ipr+1
                        q_n(s-ipr) = z(m) - Z(m);
                        J(s-ipr,1:d) = z(id) - Z(id);
                    end
                    Z(1,1:nvar) = z;
                end
            end
        end                                     % End of time loop

    case 4 %% C++: Runge Kutta implementation (= similar to cfe_nwm_odefcn)
        
        if nargout == 1
            if mem == 1
                [Z,~] = crr_cfe_nwm(mdl.tout,Z(1,1:nvar)',data,ode);
                q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
            else
                [~,q_n] = crr_cfe_nwm(mdl.tout,Z(1,1:nvar)',data,ode);
            end
            J = [];
        else
            [Z,q_n,J] = crr_cfe_nwm(mdl.tout,Z(1,1:nvar)',data,ode);
        end

end

if fail == true && (mem == 1)
    Z(s:ns,:) = repmat(Z(s-1,:), ns-s+1, 1);
end

if mem == 1
    q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
    switch nargout
        case {2,3,4}
            % diff appropriate elements of sensitivity state variables
            J = diff(Z(mdl.idx(1):mdl.idx(2),id));
   end
end
if nargout == 1
    varargout = {q_n}; return
else
    J = J .* reshape(Jth,1,[]);
    if nargout == 2
        varargout = {q_n,J}; return
    elseif nargout == 3
        varargout = {q_n,J,Jth}; return
    else
        varargout = {q_n,J,Jth,Z}; return
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                   Secondary functions listed below                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% 1. Runge Kutta integration
function [z,LTE] = rk2(t,z,h,s_max,s_fc,s_wp,k_sch,a1,k_perc,lf_thr, ...
    a2,k_lf,g_max,c_gw,mm,k_nsh,giuh,T_tr,f_dd,K,L,dT,data,T_sm,eps_m, ...
    eps,rho,m,d)

dzdtE = cfe_nwm_aug_ode(t,z,s_max,s_fc,s_wp,k_sch,a1,k_perc,lf_thr, ...
    a2,k_lf,g_max,c_gw,mm,k_nsh,giuh,T_tr,f_dd,K,L,dT,data,T_sm,eps_m, ...
    eps,rho,m,d)';                                      % Euler
zE = z + h*dzdtE;
dzdtH = cfe_nwm_aug_ode(t,zE,s_max,s_fc,s_wp,k_sch,a1,k_perc,lf_thr, ...
    a2,k_lf,g_max,c_gw,mm,k_nsh,giuh,T_tr,f_dd,K,L,dT,data,T_sm,eps_m, ...
    eps,rho,m,d)';                                      % Heun
z = z + 0.5*h*(dzdtE + dzdtH);                          % New z
LTE = abs(zE - z);                                      % LTE

end

%% 2. CFE_NWM augmented ode with sensitivities as state variables
function dzdt = cfe_nwm_aug_ode(t,z,s_max,s_fc,s_wp,k_sch,a1,k_perc, ...
    lf_thr,a2,k_lf,g_max,c_gw,mm,k_nsh,giuh,T_tr,f_dd,K,L,dT,data,T_sm, ...
    eps_m,eps,rho,m,d)

x = z(1:m);                                         % x = [Su Ss Sf1 Sf2 
                                                    %      Sf3 Q]
S = reshape(z(m+1:m*(d+1)),m,d);                    % sensitivity matrix
[dxdt,dSdt] = cfe_nwm_odefcn(t,x,S,s_max,s_fc, ...  % compute dxdt & dSdt
    s_wp,k_sch,a1,k_perc,lf_thr,a2,k_lf,g_max, ...
    c_gw,mm,k_nsh,giuh,T_tr,f_dd,K,L,dT,data,T_sm, ...
    eps_m,eps,rho,m,d);
dzdt = [dxdt; dSdt(:)];                             % repack single vector

end

%% 3. CFE_NWM augmented ode with sensitivities as state variables
function [dxdt,dSdt,Jth_f,Jx_f] = cfe_nwm_odefcn(t,x,Ssens,s_max,s_fc, ...
    s_wp,k_sch,a1,k_perc,lf_thr,a2,k_lf,g_max,c_gw,mm,k_nsh,giuh,T_tr, ...
    f_dd,K,L,dT,data,T_sm,eps_m,~,~,m,d)
% Matches the provided C implementation:
% - SWE smooth positivity (smooth_pos)
% - Rainfall-first ET uses smooth_min2(Pliq,Ep)
% - Schaake uses smooth_pos on Pe and deficit (fully smoothed)
% - Percolation uses smooth_pos activation + smooth_min2 cap
% - Lateral flow uses smooth_pos activation + smooth min cap (with Q_perc)
% - Soil ET extraction remains piecewise (as in C)
% - eps = ~ and rho = ~

% State order:
%   x = [ Swe_u; S; G; Q_1..Q_L; F1..F_K; Qcum ]
% Sensitivities: Ssens is (m x d), columns dx/dtheta_j

dxdt = nan(m,1);

% -----------------
% Smooth primitives
% -----------------
smooth_pos     = @(a,ep) 0.5*(a + sqrt(a.^2 + ep.^2));
dsmooth_pos_da = @(a,ep) 0.5*(1 + a./sqrt(a.^2 + ep.^2));

smooth_min2    = @(a,b,ep) 0.5*(a + b - sqrt((a-b).^2 + ep.^2));
dmin2_da       = @(a,b,ep) 0.5*(1 - (a-b)./sqrt((a-b).^2 + ep.^2));  % d/da

% ----------------
% 0) Unpack states
% ----------------
Swe_u = x(1);             % RAW SWE
S     = x(2);             % Soil storage
G     = x(3);             % GW storage
Q     = x(4:3+L);         % GIUH queue (mm/T)
F     = x(4+L:3+L+K);     % Nash storages (mm)
% x(m) is Qcum

% Forcing index
it = floor(t) + 1;
P = data.P(it,1);
Ep = data.Ep(it,1);
T = data.T(it,1);

% -------------------------
% 0b) SWE smooth positivity
% -------------------------
Swe        = smooth_pos(Swe_u, eps_m);
dSwe_dSweu = dsmooth_pos_da(Swe_u, eps_m);

% ---------------------------------------------------
% 1) Snow module (smooth partition + smooth melt min)
% ---------------------------------------------------
T_smeps_m = max(T_sm, eps_m);
uT = (T - T_tr) / T_smeps_m;

snow_fr = 0.5*(1 - tanh(uT));
rain_fr = 1 - snow_fr;

sech2 = 1/(cosh(uT)^2);

P_snow = P * snow_fr;
P_rain = P * rain_fr;

aT = (T - T_tr);
denom_pos = sqrt(aT^2 + T_smeps_m^2);
posT = 0.5*(aT + denom_pos);
dpos_daT = 0.5*(1 + aT/denom_pos);
dpos_dT_tr = -dpos_daT;

M_pot = f_dd * posT;

% smooth min(Swe, M_pot)
dxy   = (Swe - M_pot);
sqrtm = sqrt(dxy^2 + eps_m^2);
M     = 0.5*(Swe + M_pot - sqrtm);

dM_dSwe  = 0.5*(1 - dxy/sqrtm);
dM_dMpot = 0.5*(1 + dxy/sqrtm);

% chain to RAW Swe_u
dM_dSweu = dM_dSwe * dSwe_dSweu;

% precip partition derivatives wrt T_tr
dsnow_dT_tr = 0.5 * sech2 * (1/T_smeps_m);
drain_dT_tr = -dsnow_dT_tr;

dP_snow_dT_tr = P * dsnow_dT_tr;
dP_rain_dT_tr = P * drain_dT_tr;

% melt potential derivatives
dMpot_dT_tr  = f_dd * dpos_dT_tr;
dMpot_df_dd = posT;

dM_dT_tr   = dM_dMpot * dMpot_dT_tr;
dM_df_dd  = dM_dMpot * dMpot_df_dd;

% liquid input
Pliq = P_rain + M;

dPliq_dSweu = dM_dSweu;
dPliq_dT_tr   = dP_rain_dT_tr + dM_dT_tr;
dPliq_df_dd  = dM_df_dd;

% ---------------------------------
% 2) Rainfall-first ET (SMOOTH MIN)
% ---------------------------------
E_r = smooth_min2(Pliq, Ep, eps_m);
dEr_dPliq = dmin2_da(Pliq, Ep, eps_m);

P_e     = Pliq - E_r;
Ep_star = Ep   - E_r;

dPe_dPliq     = 1 - dEr_dPliq;
dEpstar_dPliq = -dEr_dPliq;

dPe_dSweu  = dPe_dPliq * dPliq_dSweu;
dPe_dT_tr    = dPe_dPliq * dPliq_dT_tr;
dPe_df_dd   = dPe_dPliq * dPliq_df_dd;

dEpstar_dSweu = dEpstar_dPliq * dPliq_dSweu;
dEpstar_dT_tr   = dEpstar_dPliq * dPliq_dT_tr;
dEpstar_df_dd  = dEpstar_dPliq * dPliq_df_dd;

% ----------------------------------------
% 3) Schaake partitioning (FULLY SMOOTHED)
% ----------------------------------------
%I = 0; R_s = 0;

%dI_dS = 0; dI_dSmax = 0; dI_dKsch = 0; dI_dPe = 0;
%dRs_dS = 0; dRs_dSmax = 0; dRs_dKsch = 0; dRs_dPe = 0;

D_soil = s_max - S;
eTerm  = exp(-k_sch * dT);

eps_pe = 1e-6;
eps_D  = 1e-6;

Pe_pos = smooth_pos(P_e, eps_pe);
dPePos_dPe = dsmooth_pos_da(P_e, eps_pe);

D_pos  = smooth_pos(D_soil, eps_D);
dDpos_dD = dsmooth_pos_da(D_soil, eps_D);

Ic = D_pos * (1 - eTerm);   % mm
Px = Pe_pos * dT;           % mm

denPI = Px + Ic;
I_amt = Px * (Ic/denPI);
R_amt = Px - I_amt;

I   = I_amt / dT;
R_s = R_amt / dT;

dI_dPx_amt = (Ic^2) / (denPI^2);
dI_dIc_amt = (Px^2) / (denPI^2);

% Ic derivatives
dIc_dDpos = (1 - eTerm);
dIc_dS    = dIc_dDpos * dDpos_dD * (-1);
dIc_dSmax = dIc_dDpos * dDpos_dD * (+1);
dIc_dKsch = D_pos * (dT * eTerm);

dI_dS    = (dI_dIc_amt * dIc_dS)    / dT;
dI_dSmax = (dI_dIc_amt * dIc_dSmax) / dT;
dI_dKsch = (dI_dIc_amt * dIc_dKsch) / dT;

dRs_dS    = -dI_dS;
dRs_dSmax = -dI_dSmax;
dRs_dKsch = -dI_dKsch;

% chain through Pe_pos(Pe)
dI_dPe  = (dI_dPx_amt) * dPePos_dPe;
dRs_dPe = (1 - dI_dPx_amt) * dPePos_dPe;

% chain from snow variables into I and Rs through Pe
dI_dSweu  = dI_dPe  * dPe_dSweu;
dI_dT_tr    = dI_dPe  * dPe_dT_tr;
dI_df_dd   = dI_dPe  * dPe_df_dd;

dRs_dSweu = dRs_dPe * dPe_dSweu;
dRs_dT_tr   = dRs_dPe * dPe_dT_tr;
dRs_df_dd  = dRs_dPe * dPe_df_dd;

% -------------
% 4) GIUH queue
% -------------
Qshift = [Q(2:end), 0];   % horizontal concat
dQdt   = -Q + Qshift + giuh * R_s;
Q_giuh = Q(1);

% -----------------------------------------------------------------
% 5) Soil outlets: percolation (SMOOTHED) + lateral flow (SMOOTHED)
% -----------------------------------------------------------------

% ---- Percolation: smooth activation + smooth cap (matches C) ----
% Q_perc = 0;
% dQperc_dS = 0; dQperc_dSmax = 0; dQperc_dSfc = 0; 
% dQperc_da1 = 0; dQperc_dkperc = 0;

den1 = (s_max - s_fc);         % >0 by bounds
above_raw = S - s_fc;
eps_ab = 1e-6;

above_pos = smooth_pos(above_raw, eps_ab);
dAbove_dS   = dsmooth_pos_da(above_raw, eps_ab);
dAbove_dSfc = -dAbove_dS;

r1  = above_pos / den1;
eps_r1 = 1e-12;
r1e = r1 + eps_r1;

Qunc = k_perc * (r1e^a1);

dQunc_dkperc = (r1e^a1);
dQunc_da1    = Qunc * log(r1e);
dQunc_dr1    = k_perc * a1 * r1e^(a1-1);

dr1_dS    = dAbove_dS / den1;
dr1_dSmax = -above_pos / (den1^2);
dr1_dSfc  = (dAbove_dSfc*den1 + above_pos) / (den1^2);

dQunc_dS1    = dQunc_dr1 * dr1_dS;
dQunc_dSmax1 = dQunc_dr1 * dr1_dSmax;
dQunc_dSfc   = dQunc_dr1 * dr1_dSfc;

eps_cap_perc = 1e-6;
Q_perc = smooth_min2(Qunc, above_pos, eps_cap_perc);

w_unc_perc = dmin2_da(Qunc, above_pos, eps_cap_perc); % dQperc/dQunc
w_cap_perc = 1 - w_unc_perc;                          % dQperc/dAbove

dQperc_dkperc = w_unc_perc * dQunc_dkperc;
dQperc_da1    = w_unc_perc * dQunc_da1;

dQperc_dS    = w_unc_perc * dQunc_dS1    + w_cap_perc * dAbove_dS;
dQperc_dSmax = w_unc_perc * dQunc_dSmax1;
dQperc_dSfc  = w_unc_perc * dQunc_dSfc   + w_cap_perc * dAbove_dSfc;

% ---- Lateral flow: your smooth block but corrected to match C ----
eps_on  = 1e-6;
eps_cap = 1e-6;
eps_r   = 1e-12;

a2raw = S - lf_thr;
R_on  = sqrt(a2raw^2 + eps_on^2);
above2 = 0.5*(a2raw + R_on);

dabove2_dS     = 0.5*(1 + a2raw/R_on);
dabove2_dlfthr = -dabove2_dS;

den2 = (s_max - lf_thr);           % >0 by bounds
r2   = above2 / den2;

dr2_dS     = dabove2_dS / den2;
dr2_dSmax  = -above2 / (den2^2);
dr2_dlfthr = (dabove2_dlfthr*den2 + above2) / (den2^2);

r2e = r2 + eps_r;

Qunc2 = k_lf * (r2e^a2);

dQunc2_dklf = (r2e^a2);
dQunc2_da2  = Qunc2 * log(r2e);
dQunc2_dr2  = k_lf * a2 * r2e^(a2-1);

dQunc2_dS     = dQunc2_dr2 * dr2_dS;
dQunc2_dSmax  = dQunc2_dr2 * dr2_dSmax;
dQunc2_dlfthr = dQunc2_dr2 * dr2_dlfthr;

cap2 = above2 - Q_perc;

dcap2_dS     = dabove2_dS - dQperc_dS;
dcap2_dSmax  = -dQperc_dSmax;
dcap2_dlfthr = dabove2_dlfthr;

Dcap = Qunc2 - cap2;
Rcap = sqrt(Dcap^2 + eps_cap^2);
Q_lf = 0.5*(Qunc2 + cap2 - Rcap);

w_unc = 0.5*(1 - Dcap/Rcap);  % dQlf/dQunc2
w_cap = 1 - w_unc;            % dQlf/dcap2

dQlf_dklf   = w_unc * dQunc2_dklf;
dQlf_da2    = w_unc * dQunc2_da2;
dQlf_dS     = w_unc * dQunc2_dS     + w_cap * dcap2_dS;
dQlf_dSmax  = w_unc * dQunc2_dSmax  + w_cap * dcap2_dSmax;
dQlf_dlfthr = w_unc * dQunc2_dlfthr + w_cap * dcap2_dlfthr;

% ---------------------
% 6) Soil ET extraction
% ---------------------
E_s = 0;
dEs_dS = 0; dEs_dSfc = 0; dEs_dSwp = 0;
dEs_dEpstar = 0;

if Ep_star > 0
    if S <= s_wp
        E_s = 0;
    elseif S >= s_fc
        E_s = Ep_star;
        dEs_dEpstar = 1;
    else
        den = (s_fc - s_wp);
        phi = (S - s_wp)/den;
        E_s = Ep_star * phi;

        dEs_dS      = Ep_star / den;
        dEs_dSfc    = Ep_star * (-(S - s_wp)/den^2);
        dEs_dSwp    = Ep_star * ((S - s_fc)/den^2);
        dEs_dEpstar = phi;
    end
end

dEs_dSweu = dEs_dEpstar * dEpstar_dSweu;
dEs_dT_tr   = dEs_dEpstar * dEpstar_dT_tr;
dEs_df_dd  = dEs_dEpstar * dEpstar_df_dd;

% -------------------
% 7) Groundwater flux
% -------------------
expterm  = exp(mm * G / g_max);
flux_exp = expterm - 1;
Q_gw     = c_gw * flux_exp;

dQgw_dG    = c_gw * expterm * (mm/g_max);
dQgw_dC    = flux_exp;
dQgw_dmm   = c_gw * expterm * (G/g_max);
dQgw_dGmax = c_gw * expterm * (-mm*G/(g_max^2));

% ---------------
% 8) Nash routing
% ---------------
dF = zeros(K,1);
for k = 1:K
    if k == 1
        dF(k) = Q_lf - k_nsh * F(k);
    else
        dF(k) = k_nsh * F(k-1) - k_nsh * F(k);
    end
end
Q_fast = k_nsh * F(K);

% ----------------
% 9) Assemble dxdt
% ----------------
dxdt(:) = 0;
dxdt(1) = P_snow - M;
dxdt(2) = I - Q_perc - Q_lf - E_s;
dxdt(3) = Q_perc - Q_gw;

dxdt(4:3+L) = dQdt;
dxdt(4+L:3+L+K) = dF;         % these are indices 4+L ... 3+L+K
dxdt(m) = Q_fast + Q_giuh + Q_gw;

% -------------------------------------
% 1) State Jacobian Jx_f = df/dx  (m×m)
% -------------------------------------
Jx_f = zeros(m,m);

% Indices in state x
iSwe  = 1;
iS    = 2;
iG    = 3;
iq1   = 4;
iF1   = 4 + L;
iQcum = m;

% -------------------------------------------------------------------------
% SWE' = P_snow(T_tr) - M(Swe_u,T_tr,f_dd)
% IMPORTANT: M was computed with Swe = smooth_pos(Swe_u), so df/dSwe_u uses
% dM_dSweu (NOT dM_dSwe)
% -------------------------------------------------------------------------
Jx_f(iSwe,iSwe) = -dM_dSweu;

% -----------------------------------------------------------------------
% Soil' = I - Q_perc - Q_lf - E_s
% Soil depends on SWE through I and E_s via P_e/Ep_star (depends on Pliq)
% -----------------------------------------------------------------------
Jx_f(iS,iSwe) = dI_dSweu - dEs_dSweu;

% Soil' wrt Soil
Jx_f(iS,iS) = dI_dS - dQperc_dS - dQlf_dS - dEs_dS;

% -------------------
% GW' = Q_perc - Q_gw
% -------------------
Jx_f(iG,iS) = dQperc_dS;
Jx_f(iG,iG) = -dQgw_dG;

% --------------------------------------
% GIUH: dQ/dt = -Q + shift(Q) + giuh*R_s
% --------------------------------------
Jqq = -eye(L) + diag(ones(L-1,1), +1);
Jx_f(iq1:iq1+L-1, iq1:iq1+L-1) = Jqq;
Jx_f(iq1:iq1+L-1, iS)   = giuh(:) * dRs_dS;
Jx_f(iq1:iq1+L-1, iSwe) = giuh(:) * dRs_dSweu;

% -------------------
% Nash cascade states
% -------------------
Jx_f(iF1, iS)  = dQlf_dS;
Jx_f(iF1, iF1) = -k_nsh;
for kk = 2:K
    Jx_f(iF1+kk-1, iF1+kk-2) =  k_nsh;
    Jx_f(iF1+kk-1, iF1+kk-1) = -k_nsh;
end

% --------------------------------
% Qcum' = Q(1) + k_nsh*F(K) + Q_gw
% --------------------------------
Jx_f(iQcum, iq1) = 1;           % Q_giuh = Q(1)
Jx_f(iQcum, iF1+K-1) = k_nsh;   % Q_fast = k_nsh*F(K)
Jx_f(iQcum, iG) = dQgw_dG;      % Q_gw(G)

% ---------------------------------------------------------------------
% 2) Physical-parameter Jacobian Jp_f = df/dp  (m×15)
%    p = [s_max,s_fc,s_wp,k_sch,a1,k_perc,lf_thr,a2,k_lf,g_max,c_gw,mm,
%         k_nsh,T_tr,f_dd]
% ---------------------------------------------------------------------
Jp_f = zeros(m,15);

% indices in p
pSmax = 1; pSfc = 2; pSwp = 3; pKsch = 4; pa1 = 5; pkperc = 6; pLfthr = 7;
pa2 = 8; pkLf = 9; pGmax = 10; pCgw = 11; pmm = 12; pkNsh = 13; pT_tr = 14;
pf_dd = 15;

% ----------------------------------------------------
% SWE equation: SWE' = P_snow(T_tr) - M(T_tr,f_dd,...)
% ----------------------------------------------------
Jp_f(iSwe,pT_tr)  = dP_snow_dT_tr - dM_dT_tr;
Jp_f(iSwe,pf_dd) = -dM_df_dd;

% -------------------------------------------
% Soil equation: S' = I - Q_perc - Q_lf - E_s
% -------------------------------------------
Jp_f(iS,pSmax)  =  dI_dSmax - dQperc_dSmax - dQlf_dSmax;
Jp_f(iS,pSfc)   = -dQperc_dSfc - dEs_dSfc;
Jp_f(iS,pSwp)   = -dEs_dSwp;
Jp_f(iS,pKsch)  =  dI_dKsch;

Jp_f(iS,pa1)    = -dQperc_da1;
Jp_f(iS,pkperc) = -dQperc_dkperc;
Jp_f(iS,pLfthr) = -dQlf_dlfthr;
Jp_f(iS,pa2)    = -dQlf_da2;
Jp_f(iS,pkLf)   = -dQlf_dklf;

% snow params into soil through I and E_s
Jp_f(iS,pT_tr)  = Jp_f(iS,pT_tr)  + dI_dT_tr  - dEs_dT_tr;
Jp_f(iS,pf_dd) = Jp_f(iS,pf_dd) + dI_df_dd - dEs_df_dd;

% ----------------------------------------
% Groundwater equation: G' = Q_perc - Q_gw
% ----------------------------------------
Jp_f(iG,pSmax)  = +dQperc_dSmax;
Jp_f(iG,pSfc)   = +dQperc_dSfc;
Jp_f(iG,pa1)    = +dQperc_da1;
Jp_f(iG,pkperc) = +dQperc_dkperc;

Jp_f(iG,pGmax) = -dQgw_dGmax;
Jp_f(iG,pCgw)  = -dQgw_dC;
Jp_f(iG,pmm)   = -dQgw_dmm;

% -------------------------
% GIUH block depends on R_s
% -------------------------
Jp_f(iq1:iq1+L-1,pSmax) = giuh(:) * dRs_dSmax;
Jp_f(iq1:iq1+L-1,pKsch) = giuh(:) * dRs_dKsch;
Jp_f(iq1:iq1+L-1,pT_tr)   = giuh(:) * dRs_dT_tr;
Jp_f(iq1:iq1+L-1,pf_dd)  = giuh(:) * dRs_df_dd;

% ---------------------------------------------------------------------
% Nash cascade depends on Q_lf parameters through first reservoir input
% ---------------------------------------------------------------------
Jp_f(iF1,pSmax)  = Jp_f(iF1,pSmax)  + dQlf_dSmax;
Jp_f(iF1,pLfthr) = Jp_f(iF1,pLfthr) + dQlf_dlfthr;
Jp_f(iF1,pa2)    = Jp_f(iF1,pa2)    + dQlf_da2;
Jp_f(iF1,pkLf)   = Jp_f(iF1,pkLf)   + dQlf_dklf;

% Nash cascade parameter k_nsh
Jp_f(iF1, pkNsh) = -F(1);
for kk = 2:K
    Jp_f(iF1+kk-1,pkNsh) = (F(kk-1) - F(kk));
end

% --------------------------
% Qcum parameter derivatives
% --------------------------
Jp_f(iQcum,pCgw)  = dQgw_dC;
Jp_f(iQcum,pmm)   = dQgw_dmm;
Jp_f(iQcum,pGmax) = dQgw_dGmax;
Jp_f(iQcum,pkNsh) = F(K);

% ----------------------------------------------------
% 3) Chain rule: Jth_f = df/dth = Jp_f * dp/dth  (m×d)
% ----------------------------------------------------
Jth_f = zeros(m,d);

if s_max > 0
    th2 = s_fc/s_max;
    th3 = s_wp/s_max;
else
    th2 = 0; th3 = 0;
end
th7 = 0; denom = s_fc - s_wp;
if abs(denom) > 0
    th7 = (lf_thr - s_wp)/denom;
end

% mapping derivatives
dsmax_dth1 = 1;

dsfc_dth1 = th2;
dsfc_dth2 = s_max;

dswp_dth1 = th3;
dswp_dth3 = s_max;

dlf_dth7 = (s_fc - s_wp);
dlf_dth2 = th7 * s_max;
dlf_dth3 = (1-th7) * s_max;
dlf_dth1 = (1-th7)*th3 + th7*th2;

% th1
Jth_f(:,1) = Jp_f(:,pSmax)*dsmax_dth1 + Jp_f(:,pSfc)*dsfc_dth1 + ...
    Jp_f(:,pSwp)*dswp_dth1 + Jp_f(:,pLfthr)*dlf_dth1;

% th2
Jth_f(:,2) = Jp_f(:,pSfc)*dsfc_dth2 + Jp_f(:,pLfthr)*dlf_dth2;

% th3
Jth_f(:,3) = Jp_f(:,pSwp)*dswp_dth3 + Jp_f(:,pLfthr)*dlf_dth3;

% direct mappings
Jth_f(:,4)  = Jp_f(:,pKsch);
Jth_f(:,5)  = Jp_f(:,pa1);
Jth_f(:,6)  = Jp_f(:,pkperc);
Jth_f(:,7)  = Jp_f(:,pLfthr) * dlf_dth7;

Jth_f(:,8)  = Jp_f(:,pa2);
Jth_f(:,9)  = Jp_f(:,pkLf);
Jth_f(:,10) = Jp_f(:,pGmax);
Jth_f(:,11) = Jp_f(:,pCgw);
Jth_f(:,12) = Jp_f(:,pmm);
Jth_f(:,13) = Jp_f(:,pkNsh);
Jth_f(:,14) = Jp_f(:,pT_tr);
Jth_f(:,15) = Jp_f(:,pf_dd);

% ------------------
% 4) Sensitivity RHS
% ------------------
dSdt = Jx_f * Ssens + Jth_f;

end

%% 4. CFE_NWM: RHS-only for J(x)_f and J(th)_f check
function [dxdt,dSdt,Jth_f,Jx_f] = rhs_only(t,z,th,K,L,dT,data,T_sm, ...
    eps_m,eps,rho,m,d)

% Unpack once
s_max = th(1);                      % maximum soil storage (mm)
s_fc = th(2)*s_max;                 % field-capacity storage threshold (mm)
s_wp = th(3)*s_max;                 % wilting-point storage (mm)
lf_thr = s_wp + th(7)*(s_fc-s_wp);  % lateral-flow threshold storage (mm)
k_sch = th(4);                      % Schaake adjusted magic constant (-)
a1 = th(5);                         % soil primary outlet exp. (percolatin)
k_perc = th(6);                     % soil percolation coefficient
a2 = th(8);                         % soil secondary outlet exponent
k_lf = th(9);                       % soil lateral-flow coefficient
g_max = th(10);                     % groundwater maximum storage (mm)
c_gw = th(11);                      % groundwater discharge coefficient
mm = th(12);                        % groundwater exponent
k_nsh = th(13);                     % fast-routing coefficient
T_tr = th(14);                      % Temperature threshold (°C)
f_dd = th(15);                      % Degree-day factor (mm/°C/T)

Smat = ones(m,d);
% Call your scalar-argument RHS
[dxdt,dSdt,Jth_f,Jx_f] = cfe_nwm_odefcn(t,z,Smat,s_max,s_fc,s_wp,k_sch, ...
    a1,k_perc,lf_thr,a2,k_lf,g_max,c_gw,mm,k_nsh,giuh,T_tr,f_dd,K,L,dT, ...
    data,T_sm,eps_m,eps,rho,m,d);

end
