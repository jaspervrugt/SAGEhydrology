function varargout = gr4jA(par,mdl,data,ode,check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GR4JA: Runge Kutta implementation of Gr4jA conceptual watershed model
% SYNOPSIS: varargout = gr4jA(par,mdl,data,ode,check)
%  par          dx1 vector of parameter values
%   x1:par(1)    production store capacity (mm)
%   x2:par(2)    maximum groundwater exchange flux/time step (mm/T)
%   x3:par(3)    routing store capacity (mm)
%   x4:par(4)    routing time base (controls time constants, T)
%   x5:par(5)    routing flow partitioning factor (-)
%   f_p:par(6)   pan evaporation coefficient (-)
%   T_tr:par(7)  temperature threshold (°C)
%   f_dd:par(8)  degree-day factor (mm/°C/T)
%  mdl          structure with model state/parameter info
%   .mcode       scalar numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation ode_gr4ja C++
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
                        % 4: Runge Kutta implementation ode_gr4j C++
mem = ode.mem;          % state variable storage or not
if mcode == 2 && mem == 0
    warning(['gr4jA: ' ...
        'built-in ode45 ' ...
        'solver stores ' ...
        'states: mem = 1'])
    mem = 1;
end
T_sm = 1.0;                 % smoothing width (°C) partition/positive-part
eps_m = 1e-6;               % smoothing for min()
eps_s = 1e-12;              % smoothing for state variables
n1 = mdl.n1;                % UH1 cascade order
n2 = mdl.n2;                % UH2 cascade order
kappa = 4/9;                % percol. shape constant (Perrin et al., 2003)
                            % --> kappa = 0 means no soil water percolation
bg = 3.5;                   % groundwater exchange exponent (Mathevet, 2005)
bR = 5;                     % routing outflow exp. (Perrin et al., 2003)
mts = 1.0;                  % model time step [T]
m = n1 + n2 + 4;            % # state variables [+ snow: 1 state = SWE]
ns = mdl.tout + 1;          % # print times
d = numel(par);             % # parameters
n = mdl.idx(2)-mdl.idx(1);  % # elements discharge vector
nvar = m*(d+1);             % # number of variables
rho = 1e-2;                 % smoothing coefficient
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

switch mcode

    case 1 %% MATLAB: Runge Kutta implementation
        hin = ode.InitStep;     % Initial time step
        hmax_ = ode.MaxStep;    % Maximum time step
        hmin_ = ode.MinStep;    % minimum time step
        reltol = ode.RelTol;    % Relative tolerance
        abstol = ode.AbsTol;    % Absolute tolerance
        order = ode.Order;      % Order
        maxiter = ode.maxiter;  % Maximum iterations  
        iterCount = 0; flag = 0; fail = false;

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
                    data,n1,n2,kappa,bg,bR,mts, ...
                    T_sm,eps_m,eps_s,rho,m,d);
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
                h = min(hNext,t2-t);   % Don't step past t2

                if (iterCount >= maxiter)
                    if (flag == 0)
                        warning("WARNING: Max step limit " + ...
                            "reached at t = %.5f\n", t);
                        warning("WARNING: Parameter x1  %8.5f\n", th(1));
                        warning("WARNING: Parameter x2  %8.5f\n", th(2));
                        warning("WARNING: Parameter x3  %8.5f\n", th(3));
                        warning("WARNING: Parameter x4  %8.5f\n", th(4));
                        warning("WARNING: Parameter x5  %8.5f\n", th(5));
                        warning("WARNING: Parameter f_p  %8.5f\n", th(6));
                        warning("WARNING: Parameter T_tr  %8.5f\n", th(7));
                        warning("WARNING: Parameter f_dd %8.5f\n", th(8));
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
            % -----------------------------
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
        hin = ode.InitStep;     % Initial time step
        hmax_ = ode.MaxStep;    % Maximum time step
        reltol = ode.RelTol;    % Relative tolerance
        abstol = ode.AbsTol;    % Absolute tolerance
        mdl.tout = mdl.tout-1e-10;    % Loop ode goes to maxT and then still one try
        % Then, time index goes out of bound of forcing data
        ode_options = odeset('InitialStep',hin,...  % initial time-step (T)
            'MaxStep',hmax_, ...                    % maximum time-step (T)
            'RelTol',reltol, ...                    % relative tolerance
            'AbsTol',abstol);                       % absolute tol (mm)
        [~,Z] = ode45(@(t,z) gr4jA_aug_ode(t,z,th,data,n1,n2, ...
            kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d), ...
            0:mdl.tout,Z(1,1:nvar),ode_options);
        s = size(Z,1); if s < ns, fail = true; s = s+1; end

    case 3 %% MATLAB: Explicit Euler hymod_odefcn with int_steps steps
        int_steps = 1000;                       % # integration steps
        dt = 1/int_steps;                       % integration timestep
        for s = 2:ns                            % Start time loop
            if mem == 1
                z = Z(s-1,1:nvar);              % Init. state variables
            else
                z = Z(1,1:nvar);
            end
            for it = 1:int_steps                % integrate int_steps steps
                dzdt = gr4jA_aug_ode(s-2,...    % compute dzdt based on 
                    z,th,data,n1,n2, ...        % current state, par, P, Ep
                    kappa,bg,bR,mts,T_sm, ...
                    eps_m,eps_s,rho,m,d)';
                if s == 150 && check
                    [~,~,Jth_f,Jx_f] = rhs_only(s-2,z,th,data, ...
                        n1,n2,kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
                    % --- Jth_f via 2-sided finite differences ---
                    Jth_fn = zeros(m,d);
                    for j = 1:d
                        h = 1e-6*max(1,abs(th(j)));
                        thp = th; thm = th;
                        thp(j) = thp(j) + h; thm(j) = thm(j) - h;
                        fp = rhs_only(s-2,z,thp,data,n1,n2, ...
                            kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
                        fm = rhs_only(s-2,z,thm,data,n1,n2, ...
                            kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
                        Jth_fn(:,j) = (fp - fm)/(2*h);
                    end
                    % --- Jx_f via 2-sided finite differences ---
                    Jx_fn = zeros(m,m);
                    for i = 1:m
                        h = 1e-6*max(1,abs(z(i)));
                        zp = z; zm = z;
                        zp(i) = zp(i) + h; zm(i) = zm(i) - h;
                        fp = rhs_only(s-2,zp,th,data,n1,n2, ...
                            kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
                        fm = rhs_only(s-2,zm,th,data,n1,n2, ...
                            kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
                        Jx_fn(:,i) = (fp - fm)/(2*h);
                    end
                    err_x = Jx_f(1:m,1:m) - Jx_fn(1:m,1:m);
                    err_th = Jth_f(1:m,1:d) - Jth_fn(1:m,1:d);
                    disp(err_th); disp(err_x); pause
                end
                z = z + dzdt * dt;              % update states
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

    case 4 %% C++: Runge Kutta implementation
        data.x1 = th(1);    % production store capacity (mm)
        data.x2 = th(2);    % exchange coefficient (can be negative) (mm/T)
        data.x3 = th(3);    % routing store capacity (mm)
        data.x4 = th(4);    % routing time base (controls time constant, T)
        data.x5 = th(5);    % flow partioning factor (-)
        data.f_p = th(6);   % pan evaporation coefficient (-)
        data.T_tr = th(7);  % Temperature threshold (°C)
        data.f_dd = th(8);  % Degree-day factor (mm/°C/T)

        data.n1 = n1;       % UH1-Cascade order
        data.n2 = n2;       % UH2-Cascade order
        data.kappa = kappa; % percolation shape cnst (Perrin et al., 2003)
        data.bg = bg;       % groundwater exchange exp (Mathevet, 2005)
        data.bR = bR;       % routing store exponent (Perrin et al., 2003)
        data.mts = mts;       % model time step [day]
        data.T_sm = T_sm;   % smoothing width (°C) partition/positive-part
        data.eps_m = eps_m; % smoothing for min (°C)
        data.eps_s = eps_s; % smoothing for state variables      
        data.rho = rho;     % Smoothing parameter
        data.ipr = ipr;     % Time to print

        if nargout == 1
            if mem == 1
                [Z,~] = crr_gr4jA(mdl.tout,Z(1,1:nvar)',data,ode);
                q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
            else
                [~,q_n] = crr_gr4jA(mdl.tout,Z(1,1:nvar)',data,ode);
            end
            J = [];
        else
            [Z,q_n,J] = crr_gr4jA(mdl.tout,Z(1,1:nvar)',data,ode);
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
function [z,LTE] = rk2(t,z,h,th,data,n1,n2,kappa,bg,bR,mts,T_sm, ...
    eps_m,eps_s,rho,m,d)

dzdtE = gr4jA_aug_ode(t,z,th,data,n1,n2,kappa,bg,bR,mts,T_sm,eps_m, ...
    eps_s,rho,m,d)';                % Euler
zE = z + h*dzdtE;
dzdtH = gr4jA_aug_ode(t,zE,th,data,n1,n2,kappa,bg,bR,mts,T_sm,eps_m, ...
    eps_s,rho,m,d)';                % Heun
z = z + 0.5*h*(dzdtE + dzdtH);      % New z
LTE = abs(zE - z);                  % LTE

end

%% 2. gr4jA augmented ode with sensitivities as state variables
function dzdt = gr4jA_aug_ode(t,z,th,data,n1,n2,kappa,bg,bR,mts, ...
    T_sm,eps_m,eps_s,rho,m,d)

x = z(1:m);                                     % x = [Swe S R U1 U2 Q]
S = reshape(z(m+1:m*(d+1)),m,d);                % mxd sensitivity matrix
[dxdt,dSdt] = gr4jA_odefcn(t,x,th,S,data,n1,... % compute dxdt & dSdt
    n2,kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
dzdt = [dxdt; dSdt(:)];                         % Repack into a single vector

end

%% 3) GR4JA ODE core + analytic sensitivities (Nash-Cascade routing)
function [dudt,dSdt,Jth_f,Jx_f] = gr4jA_odefcn(t,u,th,S,data,n1,n2, ...
    kappa,bg,bR,mts,T_sm,eps_m,eps_s,~,m,d)

% input arguments: rho = ~

dudt = nan(m,1);
% -------
% helpers
% -------
smooth_pos = @(a,eps) 0.5*(a + sqrt(a.^2 + eps.^2));
d_smooth_pos_da = @(a,eps) 0.5*(1 + a./sqrt(a.^2 + eps.^2));

smooth_min = @(A,B,eps) 0.5*(A + B - sqrt((A-B).^2 + eps.^2));
d_smooth_min_dA = @(A,B,eps) 0.5*(1 - (A-B)./sqrt((A-B).^2 + eps.^2));
d_smooth_min_dB = @(A,B,eps) 0.5*(1 + (A-B)./sqrt((A-B).^2 + eps.^2));

% -------------
% 1) parameters
% -------------
x1 = th(1);
x2 = th(2);
x3 = th(3);
x4 = th(4);
x5 = th(5);
f_p = th(6);
T_tr = th(7);
f_dd = th(8);

% ---------------------------------------------------------------------
% 2) indices (MATLAB 1-based), u = [Swe Sp Rr U1(1..n1) U2(1..n2) Qinf]
% ---------------------------------------------------------------------
iSwe = 1;
iSp = 2;
iRr = 3;
iU1_1 = 4;
iU1_e = 3 + n1;
iU2_1 = 4 + n1;
iU2_e = 3 + n1 + n2;
iQinf = m;

% ----------------------------------------
% 3) forcing (piecewise constant per step)
% ----------------------------------------
id = floor(t) + 1;
P = data.P(id,1);
Ep = data.Ep(id,1);
T = data.T(id,1);

% ----------------------------------------------
% 4) smooth positivity on storages (states only)
% ----------------------------------------------
Swe_u = u(iSwe);
Sp_u = u(iSp);
Rr_u = u(iRr);

Swe = smooth_pos(Swe_u, eps_s);
Sp = smooth_pos(Sp_u,  eps_s);
Rr = smooth_pos(Rr_u,  eps_s);

dSwe_du = d_smooth_pos_da(Swe_u, eps_s);
dSp_du = d_smooth_pos_da(Sp_u,  eps_s);
dRr_du = d_smooth_pos_da(Rr_u,  eps_s);

U1_u = u(iU1_1:iU1_e);
U2_u = u(iU2_1:iU2_e);

U1 = smooth_pos(U1_u, eps_s);       % Fast n1 reservoirs
U2 = smooth_pos(U2_u, eps_s);       % Slow n2 reservoirs

dU1_du = d_smooth_pos_da(U1_u, eps_s);
dU2_du = d_smooth_pos_da(U2_u, eps_s);

% -----------------------------
% 5) UH analogue time constants
% -----------------------------
tau_1 = x4 / max(1,n1);
tau_2 = (2.0*x4) / max(1,n2);

inv_tau1 = 1.0 / tau_1;
inv_tau2 = 1.0 / tau_2;

% --------------------------------------------
% 6) snow module (smooth HBV-style degree-day)
% --------------------------------------------
T_smeps_m = max(T_sm,eps_m);
utemp = (T - T_tr) / T_smeps_m;
snow_fr = 0.5*(1 - tanh(utemp));
rain_fr = 1 - snow_fr;

p_snow = P * snow_fr;
p_rain = P * rain_fr;

aT = (T - T_tr);
posT = 0.5*(aT + sqrt(aT^2 + T_smeps_m^2));
p_pmelt = f_dd * posT;

p_amelt = smooth_min(Swe,p_pmelt,eps_m);
p_liq = p_rain + p_amelt;

% snow derivatives
sech2 = 1/(cosh(utemp)^2);
dsnowFrac_dT_tr = 0.5 * sech2 / T_smeps_m;          % d(snow_fr)/dT_tr

dposT_da = 0.5*(1 + aT/sqrt(aT^2 + T_smeps_m^2));
dposT_dT_tr = -dposT_da;

dp_pmelt_dT_tr = f_dd * dposT_dT_tr;
dp_pmelt_df_dd = posT;

dM_dSwe = d_smooth_min_dA(Swe, p_pmelt,eps_m);
dM_dp_pmelt = d_smooth_min_dB(Swe,p_pmelt,eps_m);

dp_liq_dSwe = dM_dSwe;

% -------------------------------------------------------------------
% 7) production block (paper-consistent scaling; smooth wet/dry gate)
% -------------------------------------------------------------------
Ep_eff = f_p * Ep;              % mm/T
aPE = p_liq - Ep_eff;           % mm/T
daPE_dfp = -Ep;                 % f_p derivative enters only through aPE

Pn = smooth_pos( aPE, eps_s);   % mm/T
En = smooth_pos(-aPE, eps_s);   % mm/T

dPn_daPE =  d_smooth_pos_da( aPE, eps_s);     % = dPn/dp_liq
dEn_daPE = -d_smooth_pos_da(-aPE, eps_s);     % = dEn/dp_liq

% smooth wetness weight: wWet = 0.5*(1 + tanh(aPE/x1))
uPE = aPE / x1;
wWet = 0.5*(1 + tanh(uPE));
sech2PE = 1/(cosh(uPE)^2);

dwWet_daPE = 0.5 * sech2PE / x1;
dwWet_dx1 = 0.5 * sech2PE * (-aPE/(x1^2));

Sp_n = Sp / x1;

% --- Percolation: g = 1 + (kappa*Sp_n)^4 ---
%g = 1 + (kappa*Sp_n)^4;
Sp2 = Sp_n*Sp_n; Sp4 = Sp2*Sp2;
g = 1 + (kappa^4)*Sp4;
g_m14 = g^(-0.25);
Pp = Sp * (1 - g_m14);

dSpn_dSp = 1/x1;
dSpn_dx1 = -Sp/(x1^2);

dg_dSpn = 4*(kappa^4) * (Sp_n^3);
dg_m14_dg = -0.25 * g^(-1.25);

dPp_dSp = (1 - g_m14) + Sp * ( -dg_m14_dg * dg_dSpn * dSpn_dSp );
dPp_dx1 = Sp * ( -dg_m14_dg * dg_dSpn * dSpn_dx1 );

dPerc_dS = dPp_dSp;
dPerc_dx1 = dPp_dx1;

% --- Wet branch Ps_wet (uses Pn) ---
zP = Pn / x1;
aP = tanh(zP);
daP_dPn = (1 - aP^2)/x1;

Nw = (1 - Sp_n^2) * aP;
Dw = (1 + Sp_n*aP);
Ps_wet = x1 * Nw / Dw;

dNw_dSpn = (-2*Sp_n)*aP;
dNw_dPn = (1 - Sp_n^2)*daP_dPn;

dDw_dSpn = aP;
dDw_dPn = Sp_n*daP_dPn;

dPsWet_dSpn = x1 * (Dw*dNw_dSpn - Nw*dDw_dSpn) / (Dw^2);
dPsWet_dPn = x1 * (Dw*dNw_dPn  - Nw*dDw_dPn ) / (Dw^2);

dPsWet_dS = dPsWet_dSpn * dSpn_dSp;

% x1 derivative for Ps_wet
dzP_dx1 = -Pn/(x1^2);
daP_dx1 = (1 - aP^2) * dzP_dx1;

dNw_dx1 = dNw_dSpn*dSpn_dx1 + (1 - Sp_n^2)*daP_dx1;
dDw_dx1 = dDw_dSpn*dSpn_dx1 + Sp_n*daP_dx1;

dNoverD_dx1 = (Dw*dNw_dx1 - Nw*dDw_dx1)/(Dw^2);
dPsWet_dx1 = (Nw/Dw) + x1 * dNoverD_dx1;

% --- Dry branch Es_dry (paper form: Nd = Sp*(2-Sp_n)*aE) ---
zE = En / x1;
aE = tanh(zE);
daE_dEn = (1 - aE^2)/x1;

Nd = Sp * (2 - Sp_n) * aE;
Dd = 1 + (1 - Sp_n) * aE;

invD = 1/Dd;
invD2 = invD^2;

Es_dry = Nd * invD;

dNd_dSp = (2 - Sp_n) * aE;
dNd_dSpn = -Sp * aE;
dNd_daE = Sp * (2 - Sp_n);

dDd_dSpn = -aE;
dDd_daE = (1 - Sp_n);

dEs_dSp = dNd_dSp * invD;
dEs_dSpn = (dNd_dSpn*Dd - Nd*dDd_dSpn) * invD2;
dEs_daE = (dNd_daE*Dd  - Nd*dDd_daE ) * invD2;

dEsDry_dS = dEs_dSp + dEs_dSpn * dSpn_dSp;
dEsDry_dEn = dEs_daE * daE_dEn;

dzE_dx1 = -En/(x1^2);
daE_dx1 = (1 - aE^2) * dzE_dx1;
dEsDry_dx1 = dEs_dSpn * dSpn_dx1 + dEs_daE * daE_dx1;

% --- Smooth blend (include dwWet/dx1!) ---
Ps = wWet * Ps_wet;
Es = (1 - wWet) * Es_dry;

dPs_dS = wWet * dPsWet_dS;
dEs_dS = (1 - wWet) * dEsDry_dS;

dPsWet_daPE = dPsWet_dPn * dPn_daPE;
dEsDry_daPE = dEsDry_dEn * dEn_daPE;

dPs_daPE = wWet*dPsWet_daPE + Ps_wet*dwWet_daPE;
dEs_daPE = (1 - wWet)*dEsDry_daPE - Es_dry*dwWet_daPE;

% keep your state-derivative names (via p_liq)
dPs_dp_liq = dPs_daPE;    % because daPE/dp_liq = 1
dEs_dp_liq = dEs_daPE;

dPs_dx1 = wWet*dPsWet_dx1 + Ps_wet*dwWet_dx1;
dEs_dx1 = (1 - wWet)*dEsDry_dx1 - Es_dry*dwWet_dx1;

% Remaining production fluxes (consistent with blended Ps)
P1 = Pn - Ps;

% Since P1 = Pn - Ps:
dP1_dS = -dPs_dS;
dP1_daPE = dPn_daPE - dPs_dp_liq;   % dPn/dp_liq = dPn/daPE; dPs_dp_liq already includes gate terms
dP1_dx1 = -dPs_dx1;

Pr = P1 + Pp;

dPr_dS = dP1_dS + dPerc_dS;
dPr_dp_liq = dP1_daPE;              % Pp depends on Sp only
dPr_dx1 = dP1_dx1 + dPerc_dx1;

% f_p derivatives (aPE-only)
dPs_dfp = dPs_daPE * daPE_dfp;
dEs_dfp = dEs_daPE * daPE_dfp;
dP1_dfp = dP1_daPE * daPE_dfp;     % since dP1_daPE is dP1/daPE
dPr_dfp = dP1_dfp;                  % Pp independent of f_p

% -------------------------------------------------------------------------
% 8) routing + cascades
% -------------------------------------------------------------------------
in1 = x5 * Pr;
in2 = (1 - x5) * Pr;

q_1 = U1(n1) * inv_tau1;
q_2 = U2(n2) * inv_tau2;

Rratio = Rr / x3;
q_g = x2 * (Rratio^bg);
% R2 = Rratio*Rratio;
% R3 = R2*Rratio;
% q_g = x2 * R3 * sqrt(Rratio);

q_r = x3/((bR-1)*mts) * (Rratio^bR);
% q_r = x3/((bR-1)*mts) * R3 * R2;   % Rratio^5 = Rratio^3 * Rratio^2

q_d = q_2 + q_g;
q_dpos = smooth_pos(q_d, eps_s);
q_out = q_r + q_dpos;

dqg_dRr = x2 * bg * (Rratio^(bg-1)) * (1/x3);
dqg_dx2 = Rratio^bg;
dqg_dx3 = x2 * bg * (Rratio^(bg-1)) * (-Rr/(x3^2));

c0 = 1/((bR-1)*mts);
dqr_dRr = c0 * bR * (Rr^(bR-1)) * (x3^(1-bR));
dqr_dx3 = c0 * (Rr^bR) * (1-bR) * (x3^(-bR));

% ----------
% 9) ODE RHS
% ----------
dudt(iSwe) = p_snow - p_amelt;
dudt(iSp) = Ps - Es - Pp;
dudt(iRr) = q_1 + q_g - q_r;

dudt(iU1_1) = in1 - U1(1)*inv_tau1;
for k = 1:n1-1
    dudt(iU1_1 + k) = U1(k)*inv_tau1 - U1(k+1)*inv_tau1;
end
% for k = 2:n1
%     idx = iU1_1 + (k-1);
%     dudt(idx) = U1(k-1)*inv_tau1 - U1(k)*inv_tau1;
% end

dudt(iU2_1) = in2 - U2(1)*inv_tau2;
for k = 1:n2-1
    dudt(iU2_1 + k) = U2(k)*inv_tau2 - U2(k+1)*inv_tau2;
end
% for k = 2:n2
%     idx = iU2_1 + (k-1);
%     dudt(idx) = U2(k-1)*inv_tau2 - U2(k)*inv_tau2;
% end

dudt(iQinf) = q_out;

% -----------------------------------------------------------
% 10) Jx_f = df/du (then chain rule for smooth_pos on states)
% -----------------------------------------------------------
Jx_f = zeros(m,m);

Jx_f(iSwe,iSwe) = -dM_dSwe;

Jx_f(iSp,iSwe) = (dPs_dp_liq - dEs_dp_liq) * dp_liq_dSwe;
Jx_f(iSp,iSp) = dPs_dS - dEs_dS - dPerc_dS;

din1_dSwe = x5 * dPr_dp_liq * dp_liq_dSwe;
din1_dS = x5 * dPr_dS;

din2_dSwe = (1 - x5) * dPr_dp_liq * dp_liq_dSwe;
din2_dS = (1 - x5) * dPr_dS;

Jx_f(iU1_1,iSwe) = din1_dSwe;
Jx_f(iU1_1,iSp) = din1_dS;
Jx_f(iU1_1,iU1_1) = -inv_tau1;
for k = 2:n1
    idx = iU1_1 + (k-1);
    Jx_f(idx,idx-1) = inv_tau1;
    Jx_f(idx,idx) = -inv_tau1;
end

Jx_f(iU2_1,iSwe) = din2_dSwe;
Jx_f(iU2_1,iSp) = din2_dS;
Jx_f(iU2_1,iU2_1) = -inv_tau2;
for k = 2:n2
    idx = iU2_1 + (k-1);
    Jx_f(idx,idx-1) = inv_tau2;
    Jx_f(idx,idx) = -inv_tau2;
end

Jx_f(iRr,iU1_e) = inv_tau1;
Jx_f(iRr,iRr) = dqg_dRr - dqr_dRr;

dQdpos_dQd = d_smooth_pos_da(q_d, eps_s);
Jx_f(iQinf,iRr) = dqr_dRr + dQdpos_dQd * dqg_dRr;
Jx_f(iQinf,iU2_e) = dQdpos_dQd * inv_tau2;

% chain rule for smoothed nonnegative states
Jx_f(:,iSwe) = Jx_f(:,iSwe) * dSwe_du;
Jx_f(:,iSp) = Jx_f(:,iSp)  * dSp_du;
Jx_f(:,iRr) = Jx_f(:,iRr)  * dRr_du;

for k = 1:n1
    Jx_f(:, iU1_1 + (k-1)) = Jx_f(:, iU1_1 + (k-1)) * dU1_du(k);
end
for k = 1:n2
    Jx_f(:, iU2_1 + (k-1)) = Jx_f(:, iU2_1 + (k-1)) * dU2_du(k);
end

% ---------------------------------------------------------------------
% 11) Jth_f = df/dtheta (m x d), theta = [x1 x2 x3 x4 x5 f_p T_tr f_dd]
% ---------------------------------------------------------------------
Jth_f = zeros(m,d);

% snow chain terms
dPsnow_dT_tr = P * dsnowFrac_dT_tr;
dPrain_dT_tr = -P * dsnowFrac_dT_tr;

dM_dT_tr = dM_dp_pmelt * dp_pmelt_dT_tr;
dM_df_dd = dM_dp_pmelt * dp_pmelt_df_dd;

dp_liq_dT_tr = dPrain_dT_tr + dM_dT_tr;
dp_liq_df_dd = dM_df_dd;

dPs_dT_tr = dPs_dp_liq * dp_liq_dT_tr;
dPs_df_dd = dPs_dp_liq * dp_liq_df_dd;

dEs_dT_tr = dEs_dp_liq * dp_liq_dT_tr;
dEs_df_dd = dEs_dp_liq * dp_liq_df_dd;

dP1_dT_tr = dP1_daPE * dp_liq_dT_tr;
dP1_df_dd = dP1_daPE * dp_liq_df_dd;

dPr_dT_tr = dP1_dT_tr;
dPr_df_dd = dP1_df_dd;

% column mapping (1-based): 1:x1 2:x2 3:x3 4:x4 5:x5 6:f_p 7:T_tr 8:f_dd

% T_tr, f_dd
Jth_f(iSwe,7) = dPsnow_dT_tr - dM_dT_tr;
Jth_f(iSwe,8) = -dM_df_dd;

Jth_f(iSp,7) = dPs_dT_tr - dEs_dT_tr;
Jth_f(iSp,8) = dPs_df_dd - dEs_df_dd;

% x1
Jth_f(iSp,1) = dPs_dx1 - dEs_dx1 - dPerc_dx1;

% f_p
Jth_f(iSp,6) = dPs_dfp - dEs_dfp;

% x2, x3
Jth_f(iRr,2) = dqg_dx2;
Jth_f(iQinf,2) = dQdpos_dQd * dqg_dx2;

Jth_f(iRr,3) = dqg_dx3 - dqr_dx3;
Jth_f(iQinf,3) = dqr_dx3 + dQdpos_dQd * dqg_dx3;

% x4 affects inv_tau1 and inv_tau2
dinv_tau1_dx4 = -(1/(tau_1*tau_1)) * (1/max(1,n1));
dinv_tau2_dx4 = -(1/(tau_2*tau_2)) * (2/max(1,n2));

Jth_f(iU1_1,4) = -(U1(1)) * dinv_tau1_dx4;
for k = 2:n1
    idx = iU1_1 + (k-1);
    Jth_f(idx,4) = (U1(k-1) - U1(k)) * dinv_tau1_dx4;
end

Jth_f(iU2_1,4) = -(U2(1)) * dinv_tau2_dx4;
for k = 2:n2
    idx = iU2_1 + (k-1);
    Jth_f(idx,4) = (U2(k-1) - U2(k)) * dinv_tau2_dx4;
end

Jth_f(iRr,4) = U1(end) * dinv_tau1_dx4;
Jth_f(iQinf,4) = dQdpos_dQd * (U2(end) * dinv_tau2_dx4);

% x5 split
Jth_f(iU1_1,5) = Pr;
Jth_f(iU2_1,5) = -Pr;

% injections: sensitivity via Pr to x1, f_p, T_tr, f_dd
Jth_f(iU1_1,1) = x5 * dPr_dx1;
Jth_f(iU2_1,1) = (1 - x5) * dPr_dx1;

Jth_f(iU1_1,6) = x5 * dPr_dfp;
Jth_f(iU2_1,6) = (1 - x5) * dPr_dfp;

Jth_f(iU1_1,7) = x5 * dPr_dT_tr;
Jth_f(iU2_1,7) = (1 - x5) * dPr_dT_tr;

Jth_f(iU1_1,8) = x5 * dPr_df_dd;
Jth_f(iU2_1,8) = (1 - x5) * dPr_df_dd;

% -------------------
% 12) sensitivity ODE
% -------------------
dSdt = Jx_f * S + Jth_f;

end

%% 4. GR4JA: RHS-only for J(x)_f and J(th)_f check
function [dxdt,dSdt,Jth_f,Jx_f] = rhs_only(t,z,th,data,n1,n2,kappa, ...
    bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d)

Smat = ones(m,d);

[dxdt,dSdt,Jth_f,Jx_f] = gr4jA_odefcn(t,z,th,Smat,data,n1,n2,kappa, ...
    bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);

end