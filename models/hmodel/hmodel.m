function varargout = hmodel(par,mdl,data,ode,check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HMODEL: Runge Kutta implementation of Hmodel conceptual watershed model
% SYNOPSIS: varargout = hmodel(par,mdl,data,ode,check)
%  par          dx1 vector of parameter values
%   I_max:par(1) maximum storage interception reservoir
%   Su_max:par(2) maximum storage unsaturated zone
%   Q_max:par(3) maximum percolation flux
%   a_E:par(4)   evaporation parameter
%   a_F:par(5)   runoff parameter
%   r_f:par(6)   residence time fast reservoir
%   r_s:par(7)   residence time slow reservoir
%   T_tr:par(8)  temperature threshold (°C)
%   f_dd:par(9)  degree-day factor (mm/°C/T)
%  mdl          structure with model state/parameter info
%   .mcode       scalar numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation ode_hmodel C++
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
                        % 4: Runge Kutta implementation ode_hmodel C++
mem = ode.mem;          % state variable storage or not
if mcode == 2 && mem == 0
    warning(['hmodel: ' ...
        'built-in ode45 ' ...
        'solver stores ' ...
        'states: mem = 1'])
    mem = 1;
end
T_sm = 1.0;                 % smoothing width (°C) partition/positive-part
eps_m = 1e-6;               % smoothing for min()
m = 6;                      % # state variables [+ snow: 1 state = SWE]
ns = mdl.tout + 1;          % # print times
d = numel(par);             % # parameters
n = mdl.idx(2)-mdl.idx(1);  % # elements discharge vector
nvar = m*(d+1);             % # number of variables
rho = 0.01;                 % Dimensionless smoothing coefficient 
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
Z(1,1:m) = mdl.y0;          % initial state
Z(1,m+1:nvar) = 0;          % initial sensitivities

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

switch mcode

    case 1 %% MATLAB: Runge Kutta implementation
        hin = ode.InitStep;     % Initial time step
        hmax_ = ode.MaxStep;    % Maximum time step
        hmin_ = ode.MinStep;    % minimum time step
        reltol = ode.RelTol;    % Relative tolerance
        abstol = ode.AbsTol;    % Absolute tolerance
        order = ode.Order;      % Order
        maxiter = ode.maxiter;  % Maximum iterations
        iterCount = 0; flag = 0;
        
        hCarry = max(hmin_,min(hin,hmax_)); % carry adaptive recommendation

        for s = 2:ns            % Loop over elements
            t1 = s-2; t2 = s-1;                 % Set start, end times
            h = min(hCarry,t2-t1);              % reuse prior recommendation
            if mem == 1
                Z(s,1:nvar) = Z(s-1,1:nvar);
                z = Z(s,1:nvar);                    % row vector
            else
                z = Z(1,1:nvar);
            end
            t = t1;                             % Initial time
            % Integrate from t1 to t2
            while (t < t2)
                [ztmp,LTE] = rk2(t,z,h,th,...       % Evaluate rk2
                    data,T_sm,eps_m,m,d);
                if any(~isfinite(ztmp)) || any(abs(ztmp) > 1e12)
                    fail = true; break;
                end
                w = 1 ./ (reltol*abs(ztmp) + abstol);       % Weights                
                wrms = sqrt(sum((w.*LTE).^2)/nvar);         % CORRECTED, Nov. 2022
                if (wrms <= 1) || (h <= hmin_)              % Accept if error is small enough
                    z = ztmp; t = t + h;
                end
                hNext = h*max(0.2,min(5.0,0.9*wrms^(-1/order)));    % Compute new step
                hNext = max(hmin_,min(hNext,hmax_));
                hCarry = hNext;                         % retain before boundary clipping
                h = min(hNext,t2-t);   % Another turn
                if (iterCount >= maxiter)
                    if (flag == 0)
                        warning("WARNING: Max step limit " + ...
                            "reached at t = %.5f\n", t);
                        warning("WARNING: Parameter I_max  %8.5f\n", th(1));
                        warning("WARNING: Parameter Su_max %8.5f\n", th(2));
                        warning("WARNING: Parameter Q_max  %8.5f\n", th(3));
                        warning("WARNING: Parameter a_E    %8.5f\n", th(4));
                        warning("WARNING: Parameter a_F    %8.5f\n", th(5));
                        warning("WARNING: Parameter r_f    %8.5f\n", th(6));
                        warning("WARNING: Parameter r_s    %8.5f\n", th(7));                       
                        warning("WARNING: Parameter T_tr   %8.5f\n", th(8));
                        warning("WARNING: Parameter f_dd   %8.5f\n", th(9));
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

    case 2 %% MATLAB: ode45 implementation + augmented sensitivities
        hin = ode.InitStep;   % Initial time step
        hmax_ = ode.MaxStep;    % Maximum time step
        reltol = ode.RelTol;     % Relative tolerance
        abstol = ode.AbsTol;     % Absolute tolerance
        mdl.tout = mdl.tout-1e-10;    % Loop ode goes to maxT and then still one try
        % Then, time index goes out of bound of forcing data
        ode_options = odeset('InitialStep',hin,...  % initial time-step (T)
            'MaxStep',hmax_, ...                    % maximum time-step (T)
            'RelTol',reltol, ...                    % relative tolerance
            'AbsTol',abstol);                       % absolute tol (mm)
        [~,Z] = ode45(@(t,z) hmodel_aug_ode(t,z,th,data,T_sm,eps_m, ...
            m,d),0:mdl.tout,Z(1,1:nvar),ode_options);
        s = size(Z,1); if s < ns, fail = true; s = s+1; end

    case 3 %% MATLAB: Explicit Euler hymod_odefcn with int_steps steps
        int_steps = 500;                    % # integration steps
        dt = 1/int_steps;                   % integration timestep
        for s = 2:ns                        % Start time loop
            if mem == 1
                z = Z(s-1,1:nvar);          % Initialize state variables
            else
                z = Z(1,1:nvar);
            end
            for it = 1:int_steps                % integrate int_steps steps
                dzdt = hmodel_aug_ode(s-2, ...      % compute dzdt based on
                    z,th,data,T_sm,eps_m,m,d)';     % current state, x,P,Ep    
                if s == 3126 && check
                    [~,~,Jth_f,Jx_f] = rhs_only(s-2,z,th,data,T_sm, ...
                        eps_m,m,d);
                    % --- Jth_f via 2-sided finite differences ---
                    Jth_fn = zeros(m,d);
                    for j = 1:d
                        h = 1e-6*max(1,abs(th(j)));
                        thp = th; thm = th;
                        thp(j) = thp(j) + h; thm(j) = thm(j) - h;
                        fp = rhs_only(s-2,z,thp,data,T_sm, ...
                            eps_m,m,d);
                        fm = rhs_only(s-2,z,thm,data,T_sm, ...
                            eps_m,m,d);
                        Jth_fn(:,j) = (fp - fm)/(2*h);
                    end
                    % --- Jx_f via 2-sided finite differences ---
                    Jx_fn = zeros(m,m);
                    for i = 1:m
                        h = 1e-6*max(1,abs(z(i)));
                        zp = z; zm = z;
                        zp(i) = zp(i) + h; zm(i) = zm(i) - h;
                        fp = rhs_only(s-2,zp,th,data,T_sm, ...
                            eps_m,m,d);
                        fm = rhs_only(s-2,zm,th,data,T_sm, ...
                            eps_m,m,d);
                        Jx_fn(:,i) = (fp - fm)/(2*h);
                    end
                    err_x = Jx_f(1:m,1:m) - Jx_fn(1:m,1:m);
                    err_th = Jth_f(1:m,1:d) - Jth_fn(1:m,1:d);
                    disp(err_th); disp(err_x); pause
                end
                z = z + dzdt * dt;                  % update states
            end
            if any(~isfinite(z)) || any(abs(z) > 1e12)
                fail = true; break;
            else
                if mem == 1
                    Z(s,1:nvar) = z;                % State at t
                else
                    if s >= ipr+1
                        q_n(s-ipr) = z(m) - Z(m);
                        J(s-ipr,1:d) = z(id) - Z(id);
                    end
                    Z(1,1:nvar) = z;
                end
            end
        end                                         % End of time loop

    case 4 %% C++: Runge Kutta implementation ( = similar to hmodel_odefcn)
        data.I_max = th(1);     % Maximum storage of interception zone
        data.Su_max = th(2);    % Maximum storage of unsaturated zone
        data.Q_max = th(3);     % Maximum percolation rate
        data.a_E = th(4);       % Evaporation parameter
        data.a_F = th(5);       % Runoff parameter
        data.a_S = 1e-6;        % Percolation coefficient
        data.a_I = 50;          % Infiltration coefficient
        data.a_P = -50;         % Effective precipitation coefficient
        data.r_f = th(6);       % Time constant, fast reservoir
        data.r_s = th(7);       % Time constant, slow reservoir
        data.T_tr = th(8);      % Temperature threshold (°C)
        data.f_dd = th(9);      % Degree-day factor (mm/°C/T)
        
        data.T_sm = T_sm;       % smoothing width (°C) partition/positive
        data.eps_m = eps_m;     % smoothing for min (°C)    
        data.rho = rho;         % Smoothing parameter [inactive]
        data.ipr = ipr;         % Time to print
        
        if nargout == 1
            if mem == 1
                [Z,~] = crr_hmodel(mdl.tout,Z(1,1:nvar)',data,ode);
                q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
            else
                [~,q_n] = crr_hmodel(mdl.tout,Z(1,1:nvar)',data,ode);
            end
            J = [];
        else
            [Z,q_n,J] = crr_hmodel(mdl.tout,Z(1,1:nvar)',data,ode);
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
function [z,LTE] = rk2(t,z,h,th,data,T_sm,eps_m,m,d)

dzdtE = hmodel_aug_ode(t,z,th,data,T_sm,eps_m,m,d)';    % Euler
zE = z + h*dzdtE;
dzdtH = hmodel_aug_ode(t,zE,th,data,T_sm,eps_m,m,d)';   % Heun
z = z + 0.5*h*(dzdtE + dzdtH);                          % New z
LTE = abs(zE - z);                                      % LTE

end

%% 2. HMODEL augmented ode with sensitivities as state variables
function dzdt = hmodel_aug_ode(t,z,th,data,T_sm,eps_m,m,d)

x = z(1:m);                                     % x = [Swe Si Su Sf Ss Sq]
S = reshape(z(m+1:m*(d+1)),m,d);                % mxd sensitivity matrix
[dxdt,dSdt] = hmodel_odefcn(t,x,th,S,data, ...  % compute dxdt & dSdt
    T_sm,eps_m,m,d);
dzdt = [dxdt; dSdt(:)];                         % Repack single vector

end

%% 3. HMODEL: Secondary function, ODE solver
function [dxdt,dSdt,Jth_f,Jx_f] = hmodel_odefcn(t,x,th,S,data,T_sm, ...
    eps_m,m,d)

dxdt = nan(m,1);    % Initialize return argument
I_max = th(1);      % Maximum storage interception reservoir
Su_max = th(2);     % Maximum storage unsaturated zone
Q_max = th(3);      % Maximum percolation flux
a_E = th(4);        % Evaporation parameter
a_F = th(5);        % Runoff parameter
r_f = th(6);        % Residence time fast reservoir
r_s = th(7);        % Residence time slow reservoir
T_tr = th(8);      % Temperature threshold (°C)
f_dd = th(9);       % Degree-day factor (mm/°C/T)

a_S = 1e-6;         % Percolation coefficient
a_I = 50;           % Infiltration coefficient
a_P = -50;          % Effective precipitation coefficient

Swe = x(1);         % Snow water equivalent (mm)
Si = x(2);          % Interception storage (mm)
Su = x(3);          % Surface storage (mm)
Sf = x(4);          % Storage of fast/quick reservoirs (mm)
Ss = x(5);          % Storage of slow reservoir (mm)
it = floor(t) + 1;  % Truncate time to current time index
P = data.P(it,1);   % Get current rainfall (mm/T)
Ep = data.Ep(it,1); % Get current Ep (mm/T)
T = data.T(it,1);   % Get current temperature (°C)

% Keep storages nonnegative where needed
smooth_pos = @(a,eps) 0.5*(a + sqrt(a.^2 + eps.^2));
d_smooth_pos_da = @(a,eps) 0.5*(1 + a./sqrt(a.^2 + eps.^2));

Swe_u = Swe;
Swe = smooth_pos(Swe_u, eps_m);
dSwe_dx1 = d_smooth_pos_da(Swe_u, eps_m);

smooth_min = @(A,B,eps) 0.5*(A + B - sqrt((A-B).^2 + eps.^2));
smooth_max = @(A,B,eps) 0.5*(A + B + sqrt((A-B).^2 + eps.^2));
d_smooth_min_dA = @(A,B,eps) 0.5*(1 - (A-B)./sqrt((A-B).^2 + eps.^2));
d_smooth_max_dA = @(A,B,eps) 0.5*(1 + (A-B)./sqrt((A-B).^2 + eps.^2));

% clamp01(z) = min(max(z,0),1) smoothly
clamp01 = @(z,eps) smooth_min(smooth_max(z,0,eps),1,eps);
d_clamp01_dz = @(z,eps) ...
    d_smooth_min_dA(smooth_max(z,0,eps),1,eps) .* d_smooth_max_dA(z,0,eps);

% --------------------------------------------
% 3) Snow module (smooth HBV-style degree-day)
% --------------------------------------------
T_smeps_m = max(T_sm,eps_m);
uT = (T - T_tr)/T_smeps_m;
% Smooth snow fraction: ~1 when T<T_tr, ~0 when T>T_tr
snow_fr = 0.5*(1 - tanh(uT));
rain_fr = 1 - snow_fr;
P_snow = P * snow_fr;
P_rain = P * rain_fr;

% Smooth positive part: posT = max(T-T_tr,0) (smooth)
a = (T - T_tr);
posT = 0.5*(a + sqrt(a^2 + T_smeps_m^2));   % >=0, smooth
M_pot = f_dd * posT;                        % potential melt (mm/T)

% Smooth min(SWE, Mpot)
dxy = (Swe - M_pot);
sqrtm = sqrt(dxy^2 + eps_m^2);
M = 0.5*(Swe + M_pot - sqrtm);              % actual melt (mm/T)
Pliq = P_rain + M;                          % into unsat reservoir (mm/T)

% ---------------------------------------------------
%  Snow derivatives needed for coupling (x and theta)
% ---------------------------------------------------
% dsnow_fr/dT_tr
sech2 = 1/(cosh(uT)^2);
dsnow_dT_tr = 0.5 * sech2 / T_smeps_m;              % d(snow_fr)/dT_tr

% dposT/dT_tr, with posT = 0.5*(a + sqrt(a^2 + T_sm^2)), a=T-T_tr
dposT_da = 0.5*(1 + a/sqrt(a^2 + T_smeps_m^2));
dposT_dT_tr = -dposT_da;                          % negative

% M_pot derivatives
dMpot_dT_tr = f_dd * dposT_dT_tr;
dMpot_f_dd = posT;

% M derivatives (smooth min)
dM_dSwe = 0.5*(1 - dxy/sqrtm);
dM_dMpot = 0.5*(1 + dxy/sqrtm);

% Pliq derivatives
dPliq_dSwe = dM_dSwe;                              % only via melt
dPliq_dT_tr = (-P * dsnow_dT_tr) + dM_dMpot*dMpot_dT_tr; % Prain + M
dPliq_f_dd = dM_dMpot*dMpot_f_dd;

% Snow equation derivative pieces
dPsnow_dT_tr = P * dsnow_dT_tr;                    % Psnow = P*snow_fr
dM_dT_tr = dM_dMpot*dMpot_dT_tr;
dM_f_dd = dM_dMpot*dMpot_f_dd;

% ----------------------
% 4) Interception module
% ----------------------
dEvapI_dSi = 0;
dPe_dSi = 0;
dEp_e_dSi = 0;

dEvapI_dImax = 0;
dPe_dImax = 0;
dEp_e_dImax = 0;

if I_max > 0
    % --- Interception relative storage (MATCH expFlux clamp) ---
    Sir_raw = Si / I_max;
    Sir = clamp01(Sir_raw, eps_m);
    dSir_dSirraw = d_clamp01_dz(Sir_raw, eps_m);
    dSir_dSi = dSir_dSirraw * (1/I_max);
    dSir_dImax = dSir_dSirraw * (-Si/I_max^2);

    F_I = expFlux(Sir, a_I);
    G_I = expFlux(Sir, a_P);
    
    % Fluxes
    EvapI = Ep * F_I;
    P_e = Pliq * G_I;
    
    % Derivatives of F_I, G_I wrt Sir
    dF_I_dz = d_expFlux_dz(Sir, a_I);
    dG_I_dz = d_expFlux_dz(Sir, a_P);
    
    % Now derivatives of EvapI and Pe wrt Si/Imax
    dEvapI_dSi = Ep   * dF_I_dz * dSir_dSi;
    dEvapI_dImax = Ep   * dF_I_dz * dSir_dImax;
    
    dPe_dSi = Pliq * dG_I_dz * dSir_dSi;
    dPe_dImax = Pliq * dG_I_dz * dSir_dImax;
    
    % Now Ep_e and its derivatives (must come AFTER dEvapI_*)
    arg = Ep - EvapI;
    Ep_e = smooth_pos(arg, eps_m);
    dEp_e_darg = d_smooth_pos_da(arg, eps_m);
    
    dEp_e_dSi = dEp_e_darg * (-dEvapI_dSi);
    dEp_e_dImax = dEp_e_darg * (-dEvapI_dImax);
    
    % snow coupling
    dPe_dSwe = G_I * dPliq_dSwe;

else
    EvapI = 0;
    P_e = Pliq;
    Ep_e = Ep;

    dPe_dSwe = dPliq_dSwe;   % since P_e = Pliq
end

% --------------------------
% 5) Unsaturated zone module
% --------------------------
dEa_dSu = 0;
dprc_dSu = 0;
drnf_dSu = 0;

dEa_dSumax = 0;
dprc_dSumax = 0;
drnf_dSumax = 0;

if Su_max > 0

    % --- Interception relative storage (MATCH expFlux clamp) ---
    Sur_raw = Su / Su_max;
    Sur = clamp01(Sur_raw, eps_m);
    dSur_dSurraw = d_clamp01_dz(Sur_raw, eps_m);
    dSur_dSu = dSur_dSurraw * (1/Su_max);
    dSur_dSumax = dSur_dSurraw * (-Su/Su_max^2);

    dF_E_dz = d_expFlux_dz(Sur, a_E);
    dF_S_dz = d_expFlux_dz(Sur, a_S);
    dF_F_dz = d_expFlux_dz(Sur, a_F);

    dEa_dSu = Ep_e * dF_E_dz * dSur_dSu;
    dprc_dSu = Q_max * dF_S_dz * dSur_dSu;
    drnf_dSu = P_e  * dF_F_dz * dSur_dSu;

    dEa_dSumax = Ep_e * dF_E_dz * dSur_dSumax;
    dprc_dSumax = Q_max * dF_S_dz * dSur_dSumax;
    drnf_dSumax = P_e  * dF_F_dz * dSur_dSumax;
else
    Sur = 0;
end

F_E = expFlux(Sur, a_E);    % Ea factor
F_S = expFlux(Sur, a_S);    % percolation factor
F_F = expFlux(Sur, a_F);    % runoff factor

Ea = Ep_e * F_E;           % Actual evaporation (mm/T)
prc = Q_max * F_S;          % Percolation (mm/T)
rnf = P_e  * F_F;           % Runoff (mm/T)

qf = Sf/r_f;               % Outflow fast (mm/T)
qs = Ss/r_s;               % Outflow slow (mm/T)

% --------------------------
% 6) ODE system (m=6 states)
% --------------------------
dxdt(1) = P_snow - M;               % SWE
dxdt(2) = Pliq - EvapI - P_e;       % Si (liquid input to interception)
dxdt(3) = P_e - Ea - prc - rnf;     % Su
dxdt(4) = rnf - qf;                 % Sf
dxdt(5) = prc - qs;                 % Ss
dxdt(6) = qf + qs;                  % Q (streamflow from differencing)

% ---------------------------------
% 7) Jacobian Jx_f = df/dx  (m x m)
% ---------------------------------
Jx_f = zeros(m,m);

% Helpful: dEa/dSi through Ep_e(Si)
% dEa_dSi = dEp_e_dSi * F_E;

% rnf dependence on Si via P_e(Si) and on Swe via Pliq(Swe)->P_e
drnf_dSi = dPe_dSi  * F_F;
drnf_dSwe = dPe_dSwe * F_F;

% qf, qs derivatives
dqf_dSf = 1 / r_f;
dqs_dSs = 1 / r_s;

% ---- Row 1: f_Swe = Psnow - M
Jx_f(1,1) = - dM_dSwe;

% ---- Row 2: f_Si = Pliq - EvapI - P_e
% P_e depends on Swe through Pliq
Jx_f(2,1) = dPliq_dSwe - dPe_dSwe;
Jx_f(2,2) = -dEvapI_dSi - dPe_dSi;

% ---- Row 3: f_Su = P_e - Ea - prc - rnf
Jx_f(3,1) = dPe_dSwe - drnf_dSwe;               % via P_e(Swe)
%Jx_f(3,2) = dPe_dSi  - dEa_dSi - drnf_dSi;     % via Si
% since rnf = P_e * F_F and Ea = Ep_e * F_E
Jx_f(3,2) = dPe_dSi*(1 - F_F) - F_E*dEp_e_dSi;
Jx_f(3,3) = -dEa_dSu - dprc_dSu - drnf_dSu;     % via Su

% ---- Row 4: f_Sf = rnf - qf
Jx_f(4,1) = drnf_dSwe;
Jx_f(4,2) = drnf_dSi;
Jx_f(4,3) = drnf_dSu;
Jx_f(4,4) = -dqf_dSf;

% ---- Row 5: f_Ss = prc - qs
Jx_f(5,3) = dprc_dSu;
Jx_f(5,5) = -dqs_dSs;

% ---- Row 6: f_Q = qf + qs
Jx_f(6,4) = dqf_dSf;
Jx_f(6,5) = dqs_dSs;

% --- PATCH: map df/d(Swe_smooth) back to df/d(Swe_raw) ---
Jx_f(:,1) = Jx_f(:,1) * dSwe_dx1;
% ---------------------------------------------
% 8) Parameter Jacobian Jth_f = df/dth  (m x d)
% ---------------------------------------------
Jth_f = zeros(m,d);

% (i) Common pieces reused below
% Sensitivities Ea and rnf to I_max already via dEp_e_dImax, dPe_dImax
dEa_dImax = dEp_e_dImax * F_E;
drnf_dImax = dPe_dImax   * F_F;

% ---- theta(1) = I_max
j = 1;
Jth_f(2,j) = -dEvapI_dImax - dPe_dImax;                 % Si eq
Jth_f(3,j) =  dPe_dImax - dEa_dImax - drnf_dImax;       % Su eq
Jth_f(4,j) =  drnf_dImax;                               % Sf eq

% ---- theta(2) = Su_max
j = 2;
Jth_f(3,j) = -dEa_dSumax - dprc_dSumax - drnf_dSumax;   % Su eq
Jth_f(4,j) =  drnf_dSumax;                              % Sf eq
Jth_f(5,j) =  dprc_dSumax;                              % Ss eq

% ---- theta(3) = Q_max
dprc_dQmax = F_S;
j = 3;
Jth_f(3,j) = -dprc_dQmax;  % Su eq
Jth_f(5,j) =  dprc_dQmax;  % Ss eq

% ---- theta(4) = a_E
dF_E_da = d_expFlux_da(Sur, a_E);
dEa_daE = Ep_e * dF_E_da;

j = 4;
Jth_f(3,j) = -dEa_daE;     % Su eq

% ---- theta(5) = a_F
dF_F_da = d_expFlux_da(Sur, a_F);
drnf_daF = P_e * dF_F_da;

j = 5;
Jth_f(3,j) = -drnf_daF;    % Su eq
Jth_f(4,j) =  drnf_daF;    % Sf eq

% ---- theta(6) = r_f
j = 6;
Jth_f(4,j) =  Sf / r_f^2;  % f_Sf = rnf - Sf/r_f
Jth_f(6,j) = -Sf / r_f^2;  % f_Q =  Sf/r_f + Ss/r_s

% ---- theta(7) = r_s
j = 7;
Jth_f(5,j) =  Ss / r_s^2;  % f_Ss = prc - Ss/r_s
Jth_f(6,j) = -Ss / r_s^2;  % f_Q = Sf/r_f + Ss/r_s

% ---- theta(8) = T_tr  (snow threshold)
% Effects enter through Pliq and snow mass balance.
% dPliq/dT_tr already computed; 
% dPe/dT_tr = G_I*dPliq/dT_tr (or 1*dPliq/dT_tr if I_max<=0)
if I_max > 0
    dPe_dT_tr = G_I * dPliq_dT_tr;
else
    dPe_dT_tr = dPliq_dT_tr;
end
drnf_dT_tr = dPe_dT_tr * F_F;

j = 8;
Jth_f(1,j) = dPsnow_dT_tr - dM_dT_tr;   % Swe eq
Jth_f(2,j) = dPliq_dT_tr - dPe_dT_tr;   % Si eq (EvapI indpndnt of T_tr)
Jth_f(3,j) = dPe_dT_tr - drnf_dT_tr;    % Su eq (Ea/prc indpndnt of T_tr)
Jth_f(4,j) = drnf_dT_tr;                % Sf eq

% ---- theta(9) = f_dd  (degree-day factor)
if I_max > 0
    dPe_f_dd = G_I * dPliq_f_dd;
else
    dPe_f_dd = dPliq_f_dd;
end
drnf_f_dd = dPe_f_dd * F_F;

j = 9;
Jth_f(1,j) = -dM_f_dd;                  % Swe eq (Psnow indpndnt of f_dd)
Jth_f(2,j) = dPliq_f_dd - dPe_f_dd;     % Si eq
Jth_f(3,j) = dPe_f_dd - drnf_f_dd;      % Su eq
Jth_f(4,j) = drnf_f_dd;                 % Sf eq

% --------------------------------------------------
% 9) Sensitivity update: dS/dt = Jf_x * S + Jf_theta
% --------------------------------------------------
dSdt = Jx_f * S + Jth_f;

end


%% ------------------------------------------------------------------------
%%  Helper derivatives: expFlux wrt z and a
%% ------------------------------------------------------------------------
function dFdz = d_expFlux_dz(z, a)
% clamp z to [0,1] as in expFlux

if abs(a) <= 1e-5
    % linear limit: F ≈ z
    dFdz = 1.0;
else
    den = -expm1(-a);
    dFdz = a * exp(-a*z) / den;
end

end

function dFda = d_expFlux_da(z, a)
% derivative w.r.t. a of (1-exp(-a z))/(1-exp(-a))
if abs(a) < 1e-6
    % limit a -> 0: F ≈ z, no dependence on a
    dFda = 0.0;
    return
end
ea = exp(-a);
ea_z = exp(-a * z);
N = 1.0 - ea_z;
D = 1.0 - ea;
Np = z * ea_z;   % dN/da
Dp = ea;         % dD/da
dFda = (D * Np - N * Dp) / (D * D);

end


%% Exponential flux
function Q_r = expFlux(S_r,a)
% Relative flux from exponential storage-flux relation
if abs(a) <= 1e-5
    Q_r = S_r;    % approximately linear
else
    Q_r = (1.-exponen(-a*S_r))/(1.-exponen(-a));
end

end

%% Exponent
function f = exponen(x)
% Exponential function with protection against overflow
f = exp(min(300,x));

end

%% 4. HMODEL: RHS-only for J(x)_f and J(th)_f check
function [dxdt,dSdt,Jth_f,Jx_f] = rhs_only(t,z,th,data,T_sm, ...
    eps_m,m,d)

Smat = ones(m,d);
% Call RHS
[dxdt,dSdt,Jth_f,Jx_f] = hmodel_odefcn(t,z,th,Smat,data,T_sm, ...
    eps_m,m,d);

end