function varargout = hymod(par,mdl,data,ode,check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HYMOD: Runge Kutta implementation of Hymod conceptual watershed model
% SYNOPSIS: varargout = hymod(par,mdl,data,ode,check)
%  par          dx1 vector of parameter values
%   s_umax:par(1) maximum storage of surface reservoir (mm)
%   beta:par(2)  beta coefficient (-)
%   alfa:par(3)  flow partitioning factor (-)
%   K_s:par(4)   residence time slow reservoir (1/T)
%   K_f:par(5)   residence time quick reservoir (1/T)
%   T_tr:par(6)  temperature threshold (°C)
%   f_dd:par(7)  degree-day factor (mm/°C/T)
%  mdl          structure with model state/parameter info
%   .mcode       scalar numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation ode_hymod C++
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
                            % 4: Runge Kutta implementation ode_hymod C++
    mem = ode.mem;          % state variable storage or not
    if mcode == 2 ...
            && mem == 0
        warning(['hymod: ' ...
            'built-in ode45 ' ...
            'solver stores ' ...
            'states: mem = 1'])
        mem = 1;
    end
    T_sm = 1.0;                 % smoothing width (°C) partition/positive-part
    eps_m = 1e-6;               % smoothing for min()
    m = 7;                      % # state variables [+ snow: 1 state = SWE]
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
    Z(1,1:m) = mdl.y0;          % initial state
    Z(1,m+1:nvar) = 0;          % initial sensitivities
    
    switch mdl.pspace
        case 0 % hydrologic parameter values
            th = par;
            if (any(th<mdl.th_min) ...
                    || any(th>mdl.th_max))
                varargout = {nan(n,1),nan(n,d),nan(d,1),Z}; return
            end
            Jth = ones(d,1);                    % return dq_n/dth
        case 1 % normalized hydrologic parameter values
            nth = par;
            if (any(nth<0) ...
                    || any(nth>1))
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
                t = t1;                                 % Initial time
                % Integrate from t1 to t2
                while (t < t2)
                    [ztmp,LTE] = rk2(t,z,h,th,...       % Evaluate rk2
                        data,T_sm,eps_m,rho,m,d);
                    if any(~isfinite(ztmp)) ...
                            || any(abs(ztmp) > 1e12)
                        fail = true; break;
                    end
                    w = 1 ./ (reltol*abs(ztmp) + abstol);   % Weights
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
                            warning("WARNING: Parameter s_umax %8.5f\n", th(1));
                            warning("WARNING: Parameter beta   %8.5f\n", th(2));
                            warning("WARNING: Parameter alfa   %8.5f\n", th(3));
                            warning("WARNING: Parameter K_s    %8.5f\n", th(4));
                            warning("WARNING: Parameter K_f    %8.5f\n", th(5));
                            warning("WARNING: Parameter T_tr   %8.5f\n", th(6));
                            warning("WARNING: Parameter f_dd   %8.5f\n", th(7));
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
            hin = ode.InitStep;  % Initial time step
            hmax_ = ode.MaxStep;   % Maximum time step
            reltol = ode.RelTol;    % Relative tolerance
            abstol = ode.AbsTol;    % Absolute tolerance
            mdl.tout = mdl.tout-1e-10;    % Loop ode goes to maxT and then still one try
            % Then, time index goes out of bound of forcing data
            ode_options = odeset('InitialStep',hin,...  % initial time-step (T)
                'MaxStep',hmax_, ...                    % maximum time-step (T)
                'RelTol',reltol, ...                    % relative tolerance
                'AbsTol',abstol);                       % absolute tol (mm)
            [~,Z] = ode45(@(t,z) hymod_aug_ode(t,z,th,data,T_sm,eps_m, ...
                rho,m,d),0:mdl.tout,Z(1,1:nvar),ode_options);
            s = size(Z,1); if s < ns, fail = true; s = s+1; end
            
        case 3 %% MATLAB: Explicit Euler hymod_odefcn with int_steps steps
            int_steps = 50;                     % # integration steps
            dt = 1/int_steps;                   % integration timestep
            for s = 2:ns                        % Start time loop
                if mem == 1
                    z = Z(s-1,1:nvar);          % Initialize state variables
                else
                    z = Z(1,1:nvar);
                end
                for it = 1:int_steps            % integrate int_steps steps
                    dzdt = hymod_aug_ode(s-2, ...       % compute dzdt based on 
                        z,th,data,T_sm,eps_m,rho,m,d)'; % current state, par,P,Ep
                    if s == 150 && check
                        [~,~,Jth_f,Jx_f] = rhs_only(s-2,z,th,data, ...
                            T_sm,eps_m,rho,m,d);
                        % --- Jth_f via 2-sided finite differences ---
                        Jth_fn = zeros(m,d);
                        for j = 1:d
                            h = 1e-6*max(1,abs(th(j)));
                            thp = th; thm = th;
                            thp(j) = thp(j) + h; thm(j) = thm(j) - h;
                            fp = rhs_only(s-2,z,thp,data, ...
                                T_sm,eps_m,rho,m,d);
                            fm = rhs_only(s-2,z,thm,data, ...
                                T_sm,eps_m,rho,m,d);
                            Jth_fn(:,j) = (fp - fm)/(2*h);
                        end
                        % --- Jx_f via 2-sided finite differences ---
                        Jx_fn = zeros(m,m);
                        for i = 1:m
                            h = 1e-6*max(1,abs(z(i)));
                            zp = z; zm = z;
                            zp(i) = zp(i) + h; zm(i) = zm(i) - h;
                            fp = rhs_only(s-2,zp,th,data, ...
                                T_sm,eps_m,rho,m,d);
                            fm = rhs_only(s-2,zm,th,data, ...
                                T_sm,eps_m,rho,m,d);
                            Jx_fn(:,i) = (fp - fm)/(2*h);
                        end
                        err_x = Jx_f(1:m,1:m) - Jx_fn(1:m,1:m);
                        err_th = Jth_f(1:m,1:d) - Jth_fn(1:m,1:d);
                        disp(err_th); disp(err_x); pause
                    end
                    z = z + dzdt * dt;                  % update states
                end
                if any(~isfinite(z)) ...
                        || any(abs(z) > 1e12)
                    fail = true; break;
                else
                    % -----------------------------
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
            end                                 % End of time loop
    
        case 4 %% C++: Runge Kutta implementation
            data.S_umax = th(1);    % Maximum storage of unsaturated zone
            data.beta = th(2);      % Spatial variability soil moisture cap.
            data.alfa = th(3);      % Flow partitioning coefficient (-)
            data.K_s = th(4);       % Recession constant slow reservoir (1/T) 
            data.K_f = th(5);       % Recession constant fast reservoir (1/T) 
            data.T_tr = th(6);      % Temperature threshold (°C)
            data.f_dd = th(7);      % Degree-day factor (mm/°C/T)
    
            data.T_sm = T_sm;       % smoothing width (°C) partition/positive
            data.eps_m = eps_m;     % smoothing for min (°C)
            data.rho = rho;         % Smoothing parameter
            data.ipr = ipr;         % Time to print
    
            if nargout == 1
                if mem == 1
                    [Z,~] = crr_hymod(mdl.tout,Z(1,1:nvar)',data,ode);
                    q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
                else
                    [~,q_n] = crr_hymod(mdl.tout,Z(1,1:nvar)',data,ode);
                end
                J = [];
            else
                [Z,q_n,J] = crr_hymod(mdl.tout,Z(1,1:nvar)',data,ode);
            end
    
    end
    
    if fail == true ...
            && (mem == 1)
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
function [z,LTE] = rk2(t,z,h,th,data,T_sm,eps_m,rho,m,d)

    dzdtE = hymod_aug_ode(t,z,th,data,T_sm,eps_m,rho,m,d)';  % Euler
    zE = z + h*dzdtE;
    dzdtH = hymod_aug_ode(t,zE,th,data,T_sm,eps_m,rho,m,d)'; % Heun
    z = z + 0.5*h*(dzdtE + dzdtH);                           % New z
    LTE = abs(zE - z);                                       % LTE

end

%% 2. HYMOd augmented ode with sensitivities as state variables
function dzdt = hymod_aug_ode(t,z,th,data,T_sm,eps_m,rho,m,d)

    x = z(1:m);                                         % x = [Su Ss 
                                                        %      Sf1 Sf2
                                                        %      Sf3 Q]
    S = reshape(z(m+1:m*(d+1)),m,d);                    % mxd sensitivity matrix
    [dxdt,dSdt] = hymod_odefcn(t,x,th,S,data,T_sm, ...  % compute dxdt & dSdt
        eps_m,rho,m,d);    
    dzdt = [dxdt; dSdt(:)];                             % Repack into a single vector

end

%% 3. HYMOD: Secondary function, ODE solver
function [dxdt,dSdt,Jth_f,Jx_f] = hymod_odefcn(t,x,th,S,data,T_sm, ...
    eps_m,rho,m,d)
% Fully smoothed HYMOD ODE + sensitivities:
% - smooth positivity for SWE
% - smooth snow partition + melt
% - smooth clamp of s_ur = Su/s_umax to [0,1] (no branch flips)
% - (optional) smooth floor for def = max(1-s_ur, 1e-6)

    dxdt = nan(m,1);        % Initialize return argument
    s_umax = th(1);         % Maximum storage of surface reservoir (mm)
    beta = th(2);           % Beta coefficient (-)
    alfa = th(3);           % Flow partitioning factor (-)
    K_s = th(4);            % Residence time slow reservoir (1/T)
    K_f = th(5);            % Residence time quick reservoir (1/T)
    T_tr = th(6);           % Temperature threshold (°C)
    f_dd = th(7);           % Degree-day factor (mm/°C/T)
    
    Swe = x(1);             % Snow water equivalent (mm)
    Su = x(2);             % Surface storage (mm)
    Ss = x(3);             % Storage of slow reservoir (mm)
    Sf = x(4:6);           % Storage of fast/quick reservoirs (mm)
    
    id = floor(t) + 1;      % Truncate time to current time index
    P = data.P(id,1);      % Get current rainfall (mm/T)
    Ep = data.Ep(id,1);     % Get current Ep (mm/T)
    T = data.T(id,1);      % Get current temperature (°C)
    
    % -------------------------------------------------------------------------
    % Smooth primitives
    % -------------------------------------------------------------------------
    smooth_pos = @(a,ep) 0.5*(a + sqrt(a.^2 + ep.^2));
    dsmooth_pos_da = @(a,ep) 0.5*(1 + a./sqrt(a.^2 + ep.^2));
    
    smooth_neg = @(a,ep) 0.5*(-a + sqrt(a.^2 + ep.^2));     % ~ max(-a,0)
    dsmooth_neg_da = @(a,ep) 0.5*(-1 + a./sqrt(a.^2 + ep.^2));
    
    % smooth clamp z to [0,1]
    smooth_clamp01 = @(z,ep) z - smooth_pos(z-1,ep) + smooth_neg(z,ep);
    dsmooth_clamp01_dz = @(z,ep) 1 - dsmooth_pos_da(z-1,ep) + dsmooth_neg_da(z,ep);
    
    % smoothing for clamp and for "def floor"
    eps_c = 1e-6;      % dimensionless smoothing for Su/s_umax clamp
    def0 = 1e-6;      % deficit floor
    eps_def = 1e-10;    % smoothing for deficit floor (can be smaller)
    
    % Smooth max(a,b): 0.5*(a+b+sqrt((a-b)^2+eps^2))
    smooth_max2 = @(a,b,ep) 0.5*(a + b + sqrt((a-b).^2 + ep.^2));
    dsmooth_max2_da = @(a,b,ep) 0.5*(1 + (a-b)./sqrt((a-b).^2 + ep.^2));% d/da
    
    % ------------------------
    % 2) SWE smooth positivity
    % ------------------------
    Swe_u = Swe;
    Swe = smooth_pos(Swe_u, eps_m);
    dSwe_dx1 = dsmooth_pos_da(Swe_u, eps_m);
    
    % --------------------------------------------
    % 3) Snow module (smooth HBV-style degree-day)
    % --------------------------------------------
    T_smeps_m = max(T_sm,eps_m);
    uT = (T - T_tr)/T_smeps_m;
    
    snow_fr = 0.5*(1 - tanh(uT));
    rain_fr = 1 - snow_fr;
    P_snow = P * snow_fr;
    P_rain = P * rain_fr;
    
    aT = (T - T_tr);
    posT = 0.5*(aT + sqrt(aT^2 + T_smeps_m^2)); % smooth max(T-T_tr,0)
    M_pot = f_dd * posT;
    
    dxy = (Swe - M_pot);
    sqrtm = sqrt(dxy^2 + eps_m^2);
    M = 0.5*(Swe + M_pot - sqrtm);          % smooth min(Swe, M_pot)
    Pliq = P_rain + M;
    
    % ---------------------------------------------------------
    % 4) Production store with FULLY SMOOTH clamp s_ur in [0,1]
    % ---------------------------------------------------------
    Su_over = Su / s_umax;                      % raw ratio
    s_ur = smooth_clamp01(Su_over, eps_c);   % smooth clamp [0,1]
    dsur_dz = dsmooth_clamp01_dz(Su_over, eps_c);
    
    ds_ur_dSu = dsur_dz * (1/s_umax);
    ds_ur_ds_umax = dsur_dz * (-Su/(s_umax^2));
    
    % deficit: def = max(1 - s_ur, def0), smoothed
    def_raw = 1 - s_ur;
    def = smooth_max2(def_raw, def0, eps_def);
    ddef_ddefraw = dsmooth_max2_da(def_raw, def0, eps_def); % d(def)/d(def_raw)
    ddef_ds_ur = -ddef_ddefraw;
    
    logdef = log(def);
    def_beta = exp(beta * logdef);
    A = 1 - def_beta;
    qu = Pliq * A;
    
    Ea = Ep * s_ur * (1 + rho)/(s_ur + rho);
    
    % ----------
    % 5) ODE RHS
    % ----------
    dxdt(1) = P_snow - M;
    dxdt(2) = Pliq - qu - Ea;
    
    qs = (1-alfa) * qu;
    qs_o = K_s * Ss;
    dxdt(3) = qs - qs_o;
    
    qf = alfa * qu;
    qf_o = K_f * Sf;
    dxdt(4) = qf - qf_o(1);
    dxdt(5) = qf_o(1) - qf_o(2);
    dxdt(6) = qf_o(2) - qf_o(3);
    dxdt(7) = qf_o(3) + qs_o;
    
    % ---------------------
    % 6) Build Jx_f = df/dx
    % ---------------------
    Jx_f = zeros(m,m);
    
    % Snow derivatives needed for coupling
    sech2 = 1/(cosh(uT)^2);
    dsnowFrac_dT_tr = 0.5 * sech2 / T_smeps_m;
    
    dposT_da = 0.5*(1 + aT/sqrt(aT^2 + T_smeps_m^2));
    dposT_dT_tr = -dposT_da;
    
    dMpot_dT_tr = f_dd * dposT_dT_tr;
    
    dM_dSWE = 0.5*(1 - dxy/sqrtm);
    dM_dMpot = 0.5*(1 + dxy/sqrtm);
    
    dPliq_dSWE = dM_dSWE;
    
    % qu derivatives wrt Pliq and s_ur
    dqu_dPliq = A;
    
    % dA/d(def) and dA/d(s_ur)
    dA_ddef = -(def_beta) * (beta / def);
    dA_ds_ur = dA_ddef * ddef_ds_ur;
    
    dqu_ds_ur = Pliq * dA_ds_ur;
    
    % chain to Su
    dqu_dSu = dqu_ds_ur * ds_ur_dSu;
    
    % Ea derivatives wrt s_ur then Su
    dEa_ds_ur = Ep * (1 + rho) * rho / (s_ur + rho)^2;
    dEa_dSu = dEa_ds_ur * ds_ur_dSu;
    
    % also need dqu/dSWE via Pliq
    dqu_dSWE = dqu_dPliq * dPliq_dSWE;
    
    % SWE equation
    Jx_f(1,1) = -dM_dSWE;
    
    % Su equation: f2 = Pliq - qu - Ea
    Jx_f(2,1) = dPliq_dSWE - dqu_dSWE;
    Jx_f(2,2) = -dqu_dSu - dEa_dSu;
    
    % Ss equation: f3 = (1-alfa)qu - Ks Ss
    Jx_f(3,1) = (1 - alfa) * dqu_dSWE;
    Jx_f(3,2) = (1 - alfa) * dqu_dSu;
    Jx_f(3,3) = -K_s;
    
    % Sf1 equation: f4 = alfa qu - Kf Sf1
    Jx_f(4,1) = alfa * dqu_dSWE;
    Jx_f(4,2) = alfa * dqu_dSu;
    Jx_f(4,4) = -K_f;
    
    % Sf2
    Jx_f(5,4) = K_f;
    Jx_f(5,5) = -K_f;
    
    % Sf3
    Jx_f(6,5) = K_f;
    Jx_f(6,6) = -K_f;
    
    % Q equation
    Jx_f(7,3) = K_s;
    Jx_f(7,6) = K_f;
    
    % Map df/d(Swe_smooth) -> df/d(Swe_raw)
    Jx_f(:,1) = Jx_f(:,1) * dSwe_dx1;
    
    % --------------------------------
    % 7) Build Jth_f = df/dtheta (7x7)
    % --------------------------------
    Jth_f = zeros(m,d);
    
    % parameter 1: s_umax
    % qu depends on s_ur, s_ur depends on s_umax
    dqu_ds_umax = dqu_ds_ur * ds_ur_ds_umax;
    dEa_ds_umax = dEa_ds_ur * ds_ur_ds_umax;
    
    Jth_f(2,1) = -dqu_ds_umax - dEa_ds_umax;
    Jth_f(3,1) = (1 - alfa) * dqu_ds_umax;
    Jth_f(4,1) = alfa * dqu_ds_umax;
    
    % parameter 2: beta
    % qu = Pliq*(1 - exp(beta*log(def))) => dqu/dbeta = -Pliq*def_beta*log(def)
    dqu_dbeta = -Pliq * def_beta * logdef;
    
    Jth_f(2,2) = -dqu_dbeta;
    Jth_f(3,2) = (1 - alfa) * dqu_dbeta;
    Jth_f(4,2) = alfa * dqu_dbeta;
    
    % parameter 3: alfa
    Jth_f(3,3) = -qu;
    Jth_f(4,3) = qu;
    
    % parameter 4: Ks
    Jth_f(3,4) = -Ss;
    Jth_f(7,4) = Ss;
    
    % parameter 5: Kf
    Jth_f(4,5) = -Sf(1);
    Jth_f(5,5) = Sf(1) - Sf(2);
    Jth_f(6,5) = Sf(2) - Sf(3);
    Jth_f(7,5) = Sf(3);
    
    % parameter 6: T_tr
    dPsnow_dT_tr = P * dsnowFrac_dT_tr;
    dPrain_dT_tr = -P * dsnowFrac_dT_tr;
    
    dM_dT_tr = dM_dMpot * dMpot_dT_tr;
    dPliq_dT_tr = dPrain_dT_tr + dM_dT_tr;
    dqu_dT_tr = dqu_dPliq * dPliq_dT_tr;   % A * dPliq_dT_tr
    
    Jth_f(1,6) = dPsnow_dT_tr - dM_dT_tr;
    Jth_f(2,6) = dPliq_dT_tr - dqu_dT_tr;
    Jth_f(3,6) = (1 - alfa) * dqu_dT_tr;
    Jth_f(4,6) = alfa * dqu_dT_tr;
    
    % parameter 7: f_dd
    dMpot_df_dd = posT;
    dM_df_dd = dM_dMpot * dMpot_df_dd;
    
    dPliq_df_dd = dM_df_dd;
    dqu_df_dd = dqu_dPliq * dPliq_df_dd; % A * dPliq_df_dd
    
    Jth_f(1,7) = -dM_df_dd;
    Jth_f(2,7) = dPliq_df_dd - dqu_df_dd;
    Jth_f(3,7) = (1 - alfa) * dqu_df_dd;
    Jth_f(4,7) = alfa * dqu_df_dd;
    
    % ---------------------
    % 8) Sensitivity update
    % ---------------------
    dSdt = Jx_f * S + Jth_f;

end

%% 4. HYMOD: RHS-only for J(x)_f and J(th)_f check
function [dxdt,dSdt,Jth_f,Jx_f] = rhs_only(t,z,th,data,T_sm, ...
    eps_m,rho,m,d)

    Smat = ones(m,d);
    % Call RHS
    [dxdt,dSdt,Jth_f,Jx_f] = hymod_odefcn(t,z,th,Smat,data,T_sm, ...
        eps_m,rho,m,d);

end