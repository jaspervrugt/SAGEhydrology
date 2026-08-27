function varargout = gr4jB(par,mdl,data,ode,check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GR4JB: Runge Kutta implementation of Gr4j conceptual watershed model with
% analytic routing 
% SYNOPSIS: varargout = gr4jB(par,mdl,data,ode,check)
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
%                 4 Runge Kutta implementation ode_gr4jb C++
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
                        % 4: Runge Kutta implementation ode_gr4jB C++
mem = ode.mem;          % state variable storage or not
if mcode == 2 && mem == 0
    warning(['gr4jB: ' ...
        'built-in ode45 ' ...
        'solver stores ' ...
        'states: mem = 1'])
    mem = 1;
end
T_sm = 1.0;                 % smoothing width (°C) partition/positive-part
eps_m = 1e-6;               % smoothing for min()
eps_s = 1e-12;              % smoothing for state variables
kappa = 4/9;                % percolation shape constant (Perrin et al. 2003)
                            % --> kappa = 0 means no soil water percolation
bg = 3.5;                   % groundwater exchange exponent (Mathevet, 2005)
bR = 5;                     % routing outflow exponent (Perrin et al., 2003)
mts = 1.0;                  % model time step (T)
eta = 2.5; tau = 2.5;
ns = mdl.tout + 1;          % # print times
d = numel(par);             % # parameters
n = mdl.idx(2)-mdl.idx(1);  % # elements discharge vector
rho = 1e-2;                 % smoothing coefficient
fail = false;               % Default: model completes run

switch mdl.pspace
    case 0 % hydrologic parameter values
        th = par;
        if (any(th<mdl.th_min) || any(th>mdl.th_max))
            varargout = {nan(n,1),nan(n,d),nan(d,1),nan}; return
        end
        Jth = ones(d,1);                    % return dq_n/dth
    case 1 % normalized hydrologic parameter values
        nth = par;
        if (any(nth<0) || any(nth>1))
            varargout = {nan(n,1),nan(n,d),nan(d,1),nan}; return
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

x4 = th(4);                     % routing time base
[U,dUdx4] = gr4j_UH_ord( ...    % read UH ordinates
    x4,mts,eta); 
L = numel(U);                   % # of unit hydrograph ordinates
m = 4 + 2*L;                    % # state variables
nvar = m*(d+1);                 % # number of variables
id = m + (1:d)*m;               % Indices of sensitivity state variables
if mem == 0
    q_n = nan(n,1);
    J = nan(n,d); 
    ipr = mdl.idx(1);           % --> C++ code
else
    Z = nan(ns,nvar); 
    ipr = 0;
end
Z(1,1:m) = mdl.y0(1);               % Initialize state variables at time 0
Z(1,m+1:nvar) = 0;                  % Initialize sensitivity at time 0

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
                    data,U,dUdx4,L,eta,kappa,bg,bR,mts, ...
                    T_sm,eps_m,eps_s,rho,m,d);
                if any(~isfinite(ztmp)) || any(abs(ztmp) > 1e12)
                    fail = true; break;
                end
                w = 1 ./ (reltol*abs(ztmp) + abstol);           % Weights
                wrms = sqrt(sum((w.*LTE).^2)/nvar);     % CORRECTED, Nov. 2022
                if (wrms <= 1) || (h <= hmin_)          % Accept if error is small enough
                    z = ztmp; t = t + h;
                    iterCount = 0;
                else
                    iterCount = iterCount + 1;
                end
                hNext = h*max(0.2,min(5.0,0.9*wrms^(-1/order)));    % Compute new step
                hNext = max(hmin_,min(hNext,hmax_));
                hCarry = hNext;                         % retain before boundary clipping
                h = min(hNext,t2-t);   % Cannot exceed t2-t
                if (iterCount >= maxiter)
                    if (flag == 0)
                        warning("WARNING: Max step limit " + ...
                            "reached at t = %.5f\n", t);
                        warning("WARNING: Parameter x1  %8.5f\n", th(1));
                        warning("WARNING: Parameter x2  %8.5f\n", th(2));
                        warning("WARNING: Parameter x3  %8.5f\n", th(3));
                        warning("WARNING: Parameter x4  %8.5f\n", th(4));
                        warning("WARNING: Parameter x5  %8.5f\n", th(5));
                        warning("WARNING: Parameter f_p %8.5f\n", th(6));
                        warning("WARNING: Parameter T_tr %8.5f\n", th(7));
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
        [~,Z] = ode45(@(t,z) gr4jB_aug_ode(t,z,th,data, ...
            U,dUdx4,L,eta,kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d), ...
            0:mdl.tout,Z(1,1:nvar),ode_options);
        s = size(Z,1); if s < ns, fail = true; s = s+1; end

    case 3 %% MATLAB: Explicit Euler hymod_odefcn with int_steps steps
        int_steps = 1000;                   % # integration steps
        dt = 1/int_steps;                   % integration timestep
        for s = 2:ns                        % Start time loop
            if mem == 1
                z = Z(s-1,1:nvar);          % Initialize state variables
            else
                z = Z(1,1:nvar);
            end
            for it = 1:int_steps             % integration int_steps steps
                dzdt = gr4jB_aug_ode(s-2,... % compute dzdt based on 
                    z,th,data,U,dUdx4,L, ... % current state, par and P,Ep
                    eta,kappa,bg,bR,mts, ...
                    T_sm,eps_m,eps_s,rho,m,d)';
                if s == 150 && check
                    [~,~,Jth_f,Jx_f] = rhs_only(s-2,z,th,data, ...
                        U,dUdx4,L,eta,kappa,bg,bR,mts,T_sm,eps_m, ...
                        eps_s,rho,m,d);
                    % --- Jth_f via 2-sided finite differences ---
                    Jth_fn = zeros(m,d);
                    for j = 1:d
                        h = 1e-6*max(1,abs(th(j)));
                        thp = th; thm = th;
                        thp(j) = thp(j) + h; thm(j) = thm(j) - h;
                        fp = rhs_only(s-2,z,thp,data,U,dUdx4,L,eta, ...
                            kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
                        fm = rhs_only(s-2,z,thm,data,U,dUdx4,L,eta, ...
                            kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
                        Jth_fn(:,j) = (fp - fm)/(2*h);
                    end
                    % --- Jx_f via 2-sided finite differences ---
                    Jx_fn = zeros(m,m);
                    for i = 1:m
                        h = 1e-6*max(1,abs(z(i)));
                        zp = z; zm = z;
                        zp(i) = zp(i) + h; zm(i) = zm(i) - h;
                        fp = rhs_only(s-2,zp,th,data,U,dUdx4,L,eta, ...
                            kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
                        fm = rhs_only(s-2,zm,th,data,U,dUdx4,L,eta, ...
                            kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);
                        Jx_fn(:,i) = (fp - fm)/(2*h);
                    end
                    err_x = Jx_f(1:m,1:m) - Jx_fn(1:m,1:m);
                    err_th = Jth_f(1:m,1:d) - Jth_fn(1:m,1:d);
                    disp(err_th); disp(err_x); pause
                end
                z = z + dzdt * dt;          % update states
            end
            if any(~isfinite(z)) || any(abs(z) > 1e12)
                fail = true; break;
            else
                % -----------------------------
                if mem == 1
                    Z(s,1:nvar) = z;        % State at t
                else
                    if s >= ipr+1
                        q_n(s-ipr) = z(m) - Z(m);
                        J(s-ipr,1:d) = z(id) - Z(id);
                    end
                    Z(1,1:nvar) = z;
                end
            end
        end                                 % End of time loop

    case 4 %% C++: Runge Kutta implementation
        data.x1 = th(1);    % production store capacity (mm)
        data.x2 = th(2);    % exchange coefficient (can be negative) (mm/T)
        data.x3 = th(3);    % routing store capacity (mm)
        data.x4 = th(4);    % routing time base (controls time constant, T)
        data.x5 = th(5);    % flow partioning factor (-)
        data.f_p = th(6);   % pan evaporation coefficient (-)
        data.T_tr = th(7);  % Temperature threshold (°C)
        data.f_dd = th(8);  % Degree-day factor (mm/°C/T)

        data.U = U;         % unit hydrograph ordinates
        data.dUdx4 = dUdx4; % derivative UH ordinates with respect to x4
        data.L = L;         % # ordinates of the unit hydrograph
        data.eta = eta;     % eta value of GR4J
        data.tau = tau;     % ADD
        data.kappa = kappa; % percolation shape con. (Perrin et al., 2003)
        data.bg = bg;       % groundwater exchange exp. (Mathevet, 2005)
        data.bR = bR;       % routing store exponent (Perrin et al., 2003)
        data.mts = mts;     % model time step [T]
        data.T_sm = T_sm;   % smoothing width (°C) partition/positive-part
        data.eps_m = eps_m; % smoothing for min (°C)
        data.eps_s = eps_s; % smoothing for state svariables
        data.rho = rho;     % Smoothing parameter
        data.ipr = ipr;     % Time to print

        if nargout == 1
            if mem == 1
                [Z,~] = crr_gr4jB(mdl.tout,Z(1,1:nvar)',data,ode);
                q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
            else
                [~,q_n] = crr_gr4jB(mdl.tout,Z(1,1:nvar)',data,ode);
            end
            J = [];
        else
            [Z,q_n,J] = crr_gr4jB(mdl.tout,Z(1,1:nvar)',data,ode);
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
function [z,LTE] = rk2(t,z,h,th,data,U,dUdx4,L,eta,kappa,bg,bR,mts, ...
    T_sm,eps_m,eps_s,rho,m,d)

dzdtE = gr4jB_aug_ode(t,z,th,data,U,dUdx4,L,eta,kappa,bg,bR,mts, ...
    T_sm,eps_m,eps_s,rho,m,d)';                 % Euler
zE = z + h*dzdtE;
dzdtH = gr4jB_aug_ode(t,zE,th,data,U,dUdx4,L,eta,kappa,bg,bR,mts, ...
    T_sm,eps_m,eps_s,rho,m,d)';                 % Heun
z = z + 0.5*h*(dzdtE + dzdtH);                  % New z
LTE = abs(zE - z);                              % LTE

end

%% 2. GR4J augmented ode with sensitivities as state variables
function dzdt = gr4jB_aug_ode(t,z,th,data,U,dUdx4,L,eta,kappa,bg,bR, ...
    mts,T_sm,eps_m,eps_s,rho,m,d)

x = z(1:m);                                     % x = [Swe S R U1 U2 Q]
S = reshape(z(m+1:m*(d+1)),m,d);                % mxd sensitivity matrix
[dxdt,dSdt] = gr4jB_odefcn(t,x,th,S,data, ...   % compute dxdt & dSdt
    U,dUdx4,L,eta,kappa,bg,bR,mts,T_sm,eps_m, ...
    eps_s,rho,m,d);
dzdt = [dxdt; dSdt(:)];                         % Repack single vector

end

function [dudt,dSdt,Jth_f,Jx_f] = gr4jB_odefcn(t,u,th,S,data, ...
    U,dUdx4,L,~,kappa,bg,bR,mts,T_sm,eps_m,eps_s,~,m,d)
%GR4JB_ODEFCN  ODE stores + analytic routing via UH convolution
% (Mathias-style CIF)
%
% Notes
% -----
% This "B" variant follows the unit-hydrograph convolution approach used by
% Mathias (CIF/CKF routing kernel), i.e., the two routing branches are
% obtained via discrete convolution of the effective rainfall (after
% production/percolation) rather than via an explicit Nash-cascade state
% augmentation. 
% eta = ~ and rho = ~
%
% State layout:
% [Swe Sp Rr z1(1:L) z2(1:L) Qinf]
% Returns:
% dudt  : m x 1
% dSdt  : m x d
% Jx_f  : m x m
% Jth_f : m x d

u = u(:)';

% -------
% helpers
% -------
smooth_pos = @(a,eps) 0.5*(a + sqrt(a.^2 + eps.^2));
d_smooth_pos_da = @(a,eps) 0.5*(1 + a./sqrt(a.^2 + eps.^2));

smooth_min = @(A,B,eps) 0.5*(A + B - sqrt((A-B).^2 + eps.^2));
d_smooth_min_dA = @(A,B,eps) 0.5*(1 - (A-B)./sqrt((A-B).^2 + eps.^2));
d_smooth_min_dB = @(A,B,eps) 0.5*(1 + (A-B)./sqrt((A-B).^2 + eps.^2));

% -------
% indices
% -------
iSwe = 1;
iSp = 2;
iRr = 3;
iz1 = 4;
iz2 = iz1 + L;
iQ = iz2 + L;

% ----------
% parameters
% ----------
x1 = th(1);
x2 = th(2);
x3 = th(3);
%x4 = th(4); % used through U outside; gradient w.r.t x4 needs dU/dx4 if desired
x5 = th(5);
f_p = th(6);
T_tr = th(7);
f_dd = th(8);

% ----------------------------
% forcing (piecewise constant)
% ----------------------------
id = floor(t) + 1;
P = data.P(id,1);
Ep = data.Ep(id,1);
T = data.T(id,1);

% -----------------------------
% smooth positivity on storages
% -----------------------------
Swe_u = u(iSwe);
Sp_u = u(iSp);
Rr_u = u(iRr);

Swe = smooth_pos(Swe_u,eps_s);
Sp = smooth_pos(Sp_u,eps_s);
Rr = smooth_pos(Rr_u,eps_s);

dSwe_du = d_smooth_pos_da(Swe_u,eps_s);
dSp_du = d_smooth_pos_da(Sp_u,eps_s);
dRr_du = d_smooth_pos_da(Rr_u,eps_s);

z_1 = u(iz1:iz1+L-1);    % do NOT smooth: these are linear memory states
z_2 = u(iz2:iz2+L-1);

% ---------------------------------------------
% 6) snow module (from GR4JB, plus derivatives)
% ---------------------------------------------
T_smeps_m = max(T_sm, eps_m);
utemp = (T - T_tr) / T_smeps_m;
snow_fr = 0.5*(1 - tanh(utemp));
rain_fr = 1 - snow_fr;

p_snow = P * snow_fr;
p_rain = P * rain_fr;

aT = (T - T_tr);
posT = 0.5*(aT + sqrt(aT^2 + T_smeps_m^2));
p_pmelt = f_dd * posT;

p_amelt = smooth_min(Swe, p_pmelt, eps_m);
p_liq = p_rain + p_amelt;

% snow derivatives
sech2 = 1/(cosh(utemp)^2);
dsnowFrac_dT_tr = 0.5 * sech2 / T_smeps_m;

dposT_da = 0.5*(1 + aT/sqrt(aT^2 + T_smeps_m^2));
dposT_dTtr = -dposT_da;

dp_pmelt_dTtr = f_dd * dposT_dTtr;
dp_pmelt_dfdd = posT;

dM_dSwe = d_smooth_min_dA(Swe, p_pmelt, eps_m);
dM_dp_pmelt = d_smooth_min_dB(Swe, p_pmelt, eps_m);

dp_liq_dSwe = dM_dSwe;

% -------------------
% 7) production block
% -------------------
Ep_eff = f_p * Ep;              % mm/T
aPE = p_liq - Ep_eff;           % mm/T
daPE_dfp = -Ep;                 % f_p derivative enters only through aPE

Pn = smooth_pos( aPE, eps_s);
En = smooth_pos(-aPE, eps_s);

dPn_daPE =  d_smooth_pos_da( aPE, eps_s);
dEn_daPE = -d_smooth_pos_da(-aPE, eps_s);

% wetness gate
uPE = aPE / x1;
wWet = 0.5*(1 + tanh(uPE));
sech2PE = 1/(cosh(uPE)^2);

dwWet_daPE = 0.5 * sech2PE / x1;
dwWet_dx1 = 0.5 * sech2PE * (-aPE/(x1^2));

Sp_n = Sp / x1;

% percolation
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

% wet branch Ps_wet
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

dzP_dx1 = -Pn/(x1^2);
daP_dx1 = (1 - aP^2) * dzP_dx1;

dNw_dx1 = dNw_dSpn*dSpn_dx1 + (1 - Sp_n^2)*daP_dx1;
dDw_dx1 = dDw_dSpn*dSpn_dx1 + Sp_n*daP_dx1;

dNoverD_dx1 = (Dw*dNw_dx1 - Nw*dDw_dx1)/(Dw^2);
dPsWet_dx1 = (Nw/Dw) + x1 * dNoverD_dx1;

% dry branch Es_dry
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

% smooth blend
Ps = wWet * Ps_wet;
Es = (1 - wWet) * Es_dry;

dPs_dS = wWet * dPsWet_dS;
dEs_dS = (1 - wWet) * dEsDry_dS;

dPsWet_daPE = dPsWet_dPn * dPn_daPE;
dEsDry_daPE = dEsDry_dEn * dEn_daPE;

dPs_daPE = wWet*dPsWet_daPE + Ps_wet*dwWet_daPE;
dEs_daPE = (1 - wWet)*dEsDry_daPE - Es_dry*dwWet_daPE;

dPs_dp_liq = dPs_daPE;
dEs_dp_liq = dEs_daPE;

dPs_dx1 = wWet*dPsWet_dx1 + Ps_wet*dwWet_dx1;
dEs_dx1 = (1 - wWet)*dEsDry_dx1 - Es_dry*dwWet_dx1;

% remaining production fluxes
P1 = Pn - Ps;

dP1_dS = -dPs_dS;
dP1_daPE = dPn_daPE - dPs_dp_liq;
dP1_dx1 = -dPs_dx1;

Pr = P1 + Pp;

dPr_dS = dP1_dS   + dPerc_dS;
dPr_dp_liq = dP1_daPE;
dPr_dx1 = dP1_dx1  + dPerc_dx1;

% f_p derivatives (aPE-only)
dPs_dfp = dPs_daPE * daPE_dfp;
dEs_dfp = dEs_daPE * daPE_dfp;
dP1_dfp = dP1_daPE * daPE_dfp;
dPr_dfp = dP1_dfp;

% snow chain terms for parameters
dPsnow_dTtr = P * dsnowFrac_dT_tr;
dPrain_dTtr = -P * dsnowFrac_dT_tr;

dM_dTtr = dM_dp_pmelt * dp_pmelt_dTtr;
dM_dfdd = dM_dp_pmelt * dp_pmelt_dfdd;

dp_liq_dTtr = dPrain_dTtr + dM_dTtr;
dp_liq_dfdd = dM_dfdd;

dPs_dTtr = dPs_dp_liq * dp_liq_dTtr;
dPs_dfdd = dPs_dp_liq * dp_liq_dfdd;

dEs_dTtr = dEs_dp_liq * dp_liq_dTtr;
dEs_dfdd = dEs_dp_liq * dp_liq_dfdd;

dP1_dTtr = dP1_daPE * dp_liq_dTtr;
dP1_dfdd = dP1_daPE * dp_liq_dfdd;

dPr_dTtr = dP1_dTtr;   % Pp independent of Ttr, fdd
dPr_dfdd = dP1_dfdd;

% derivatives wrt Swe (through melt -> p_liq)
dPr_dSwe = dPr_dp_liq * dp_liq_dSwe;

% ---------------------------------
% routing inputs and UH memory ODEs
% ---------------------------------
q_in1 = x5 * Pr;
q_in2 = (1 - x5) * Pr;

dz_1 = zeros(L,1);
dz_2 = zeros(L,1);

dz_1(1) = q_in1 - z_1(1);
dz_2(1) = q_in2 - z_2(1);
for k=2:L
    dz_1(k) = z_1(k-1) - z_1(k);
    dz_2(k) = z_2(k-1) - z_2(k);
end

q_1 = z_1 * U;   % scalar
q_2 = z_2 * U;   % scalar

% routing store + discharge
Rratio = Rr / x3;
q_g = x2 * (Rratio^bg);
q_r = x3/((bR-1)*mts) * (Rratio^bR);

q_d = q_2 + q_g;
dQdpos_dQd = d_smooth_pos_da(q_d, eps_s);
q_out = q_r + smooth_pos(q_d, eps_s);

% derivatives for routing store
dqg_dRr = x2 * bg * (Rratio^(bg-1)) * (1/x3);
dqg_dx2 = (Rratio^bg);
dqg_dx3 = x2 * bg * (Rratio^(bg-1)) * (-Rr/(x3^2));

c0 = 1/((bR-1)*mts);
dqr_dRr = c0 * bR * (Rr^(bR-1)) * (x3^(1-bR));
dqr_dx3 = c0 * (Rr^bR) * (1-bR) * (x3^(-bR));

% explicit x4 dependence through U(x4)
dq1_dx4 = z_1 * dUdx4;   % scalar
dq2_dx4 = z_2 * dUdx4;   % scalar

% -------
% ODE RHS
% -------
dudt = zeros(m,1);
dudt(iSwe) = p_snow - p_amelt;
dudt(iSp) = Ps - Es - Pp;
dudt(iRr) = q_1 + q_g - q_r;
dudt(iz1:iz1+L-1) = dz_1;
dudt(iz2:iz2+L-1) = dz_2;
dudt(iQ) = q_out;

% ------------------------
% Jacobians Jx_f and Jth_f
% ------------------------
Jx_f = zeros(m,m);
Jth_f = zeros(m,d);

% --- physical-store Jacobian pieces (same logic as GR4JB, mapped) ---
% Swe equation: dudt(Swe)=p_snow - p_amelt
Jx_f(iSwe,iSwe) = -dM_dSwe;

% Sp equation: dudt(Sp)=Ps - Es - Pp
Jx_f(iSp,iSwe) = (dPs_dp_liq - dEs_dp_liq) * dp_liq_dSwe; % via melt
Jx_f(iSp,iSp) = (dPs_dS - dEs_dS - dPerc_dS);            % via Sp

% UH memory injections depend on Pr(Swe,Sp)
Jx_f(iz1, iSwe) = x5 * dPr_dSwe;
Jx_f(iz1, iSp) = x5 * dPr_dS;
Jx_f(iz1, iz1) = -1;

% dz2(1) = (1-x5)*Pr - z2(1)
Jx_f(iz2, iSwe) = (1-x5) * dPr_dSwe;
Jx_f(iz2, iSp) = (1-x5) * dPr_dS;
Jx_f(iz2, iz2) = -1;

% shift register dynamics
for k=2:L
    Jx_f(iz1+k-1, iz1+k-2) = 1;
    Jx_f(iz1+k-1, iz1+k-1) = -1;
    Jx_f(iz2+k-1, iz2+k-2) = 1;
    Jx_f(iz2+k-1, iz2+k-1) = -1;
end

% routed inflow to routing store: q1 depends on z1; q2 affects q_out (via qd)
Jx_f(iRr, iz1:iz1+L-1) = U';
Jx_f(iQ,  iz2:iz2+L-1) = dQdpos_dQd * U';

% routing store dependence on Rr
Jx_f(iRr, iRr) = dqg_dRr - dqr_dRr;
Jx_f(iQ,  iRr) = dqr_dRr + dQdpos_dQd * dqg_dRr;

% chain rule for smooth_pos'd physical states
Jx_f(:,iSwe) = Jx_f(:,iSwe) * dSwe_du;
Jx_f(:,iSp) = Jx_f(:,iSp)  * dSp_du;
Jx_f(:,iRr) = Jx_f(:,iRr)  * dRr_du;

% ------------------------
% Parameter Jacobian Jth_f
% ------------------------

% (7) T_tr and (8) f_dd affect Swe + Sp via snow/melt and p_liq
Jth_f(iSwe,7) = dPsnow_dTtr - dM_dTtr;
Jth_f(iSwe,8) = -dM_dfdd;

Jth_f(iSp,7) = dPs_dTtr - dEs_dTtr;
Jth_f(iSp,8) = dPs_dfdd - dEs_dfdd;

% x1 affects Sp equation (Ps,Es,Pp)
Jth_f(iSp,1) = dPs_dx1 - dEs_dx1 - dPerc_dx1;

% f_p affects Sp equation through aPE
Jth_f(iSp,6) = dPs_dfp - dEs_dfp;

% routing store parameter effects (x2,x3)
Jth_f(iRr,2) = dqg_dx2;
Jth_f(iRr,3) = dqg_dx3 - dqr_dx3;

Jth_f(iQ,2) = dQdpos_dQd * dqg_dx2;
Jth_f(iQ,3) = dqr_dx3 + dQdpos_dQd * dqg_dx3;

Jth_f(iRr,4) = Jth_f(iRr,4) + dq1_dx4;
Jth_f(iQ, 4) = Jth_f(iQ, 4) + dQdpos_dQd * dq2_dx4;

% UH memory injections: explicit parameter effects
% dz1(1) depends on x5 and on Pr's parameter sensitivities
Jth_f(iz1,5) = Pr;                         % d(x5*Pr)/dx5
Jth_f(iz2,5) = -Pr;                        % d((1-x5)*Pr)/dx5

Jth_f(iz1,1) = x5 * dPr_dx1;
Jth_f(iz2,1) = (1-x5) * dPr_dx1;

Jth_f(iz1,6) = x5 * dPr_dfp;
Jth_f(iz2,6) = (1-x5) * dPr_dfp;

Jth_f(iz1,7) = x5 * dPr_dTtr;
Jth_f(iz2,7) = (1-x5) * dPr_dTtr;

Jth_f(iz1,8) = x5 * dPr_dfdd;
Jth_f(iz2,8) = (1-x5) * dPr_dfdd;

% ---------------
% sensitivity ODE
% ---------------
dSdt = Jx_f * S + Jth_f;

end

%% 4. GR4J: RHS-only for J(x)_f and J(th)_f check
function [dxdt,dSdt,Jth_f,Jx_f] = rhs_only(t,z,th,data,U,L,eta, ...
    kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d)

Smat = ones(m,d);

[dxdt,dSdt,Jth_f,Jx_f] = gr4j_odefcn(t,z,th,Smat,data,U,L,eta, ...
    kappa,bg,bR,mts,T_sm,eps_m,eps_s,rho,m,d);

end

% GR4J unit hydrograph functions
function [U,dUdx4] = gr4j_UH_ord(x4,mts,eta)

p = eta + 1;
N = max(1, ceil((2*x4)/mts));

t_k = (1:N)' * mts;
t_km1 = t_k - mts;

Gk = gr4j_CIF(t_k,x4,eta);
Gkm1 = gr4j_CIF(t_km1,x4,eta);

dGk = gr4j_dCIF_dx4(t_k,x4,p);
dGkm1 = gr4j_dCIF_dx4(t_km1,x4,p);

Uraw = Gk - Gkm1;
dUrawdx = dGk - dGkm1;
eps_U = 1e-12 + 1e-6*max(1, max(abs(Uraw)));   % scale-aware

Uraw_pos = 0.5*(Uraw + sqrt(Uraw.^2 + eps_U^2));
dpos_dUraw = 0.5*(1 + Uraw ./ sqrt(Uraw.^2 + eps_U^2));

Uraw = Uraw_pos;
dUrawdx = dUrawdx .* dpos_dUraw;

s = sum(Uraw);
ds = sum(dUrawdx);

s_safe = max(s, 1e-30);
U = Uraw / s_safe;
dUdx4 = (dUrawdx*s_safe - Uraw*ds) / (s_safe*s_safe);
end

function dG = gr4j_dCIF_dx4(t, x4, p)
% p = eta+1
dG = zeros(size(t));
u = t ./ x4;

% 0 < t < x4: G = u^p  => dG/dx4 = -(p/x4)*u^p
i2 = (t > 0) & (t < x4);
dG(i2) = -(p/x4) * (u(i2).^p);

% x4 <= t < 2x4: G = 1-(2-u)^p => dG/dx4 = -p*(2-u)^(p-1) * (t/x4^2)
i3 = (t >= x4) & (t < 2*x4);
dG(i3) = -p * (2 - u(i3)).^(p-1) .* (t(i3) ./ (x4^2));
end

function G = gr4j_CIF(t, x4, eta)
%GR4J_CIF  Cumulative integral function G(t) for the GR4J CKF.
% Piecewise over t in (-inf,0), [0,x4], [x4,2*x4], [2*x4,inf).
%
% This corresponds to the symmetric GR4J CKF support [0,2*x4]
% used by Mathias for runoff attenuation.

td = t ./ x4;   % nondimensional time
G = zeros(size(t));

% region 1: t <= 0 -> 0
i1 = (t <= 0);
G(i1) = 0;

% region 2: 0 < t < x4  (0 < td < 1)
i2 = (t > 0) & (t < x4);
% integral of g(t) = ((eta+1)/x4) * (t/x4)^eta over [0,t]
% => G = (t/x4)^(eta+1)
G(i2) = td(i2).^(eta+1);

% region 3: x4 <= t < 2*x4 (1 <= td < 2)
i3 = (t >= x4) & (t < 2*x4);
% symmetric second branch:
% G = 1 - (2 - td)^(eta+1)
G(i3) = 1 - (2 - td(i3)).^(eta+1);

% region 4: t >= 2*x4 -> 1
i4 = (t >= 2*x4);
G(i4) = 1;
end