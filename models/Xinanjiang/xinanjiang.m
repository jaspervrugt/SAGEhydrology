function varargout = xinanjiang(par,mdl,data,ode,check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%XINANJIANG: Runge Kutta implementation of Xinanjiang conceptual watershed
% model
% SYNOPSIS: varargout = xinanjiang(par,mdl,data,ode,check)
%  par          dx1 vector of parameter values
%   f_p:par(1)   ratio of potential evaporion to pan evaporation
%   A_im:par(2)  fraction impervious area (-)
%   a:par(3)     tension water distribution inflection parameter (-)
%   b:par(4)     tension water distribution shape parameter (-)
%   c:par(7)     fraction of LM for second evaporation change (-)
%   Ex:par(9)    free water distribution shape parameter (-)
%   k_i:par(10)  free water interflow parameter (1/T)
%   k_g:par(11)  free water groundwater parameter (1/T)
%   c_i:par(12)  interflow time coefficient (1/T)
%   c_g:par(13)  baseflow time coefficient (1/T)
%   k_f:par(14)  recession constant fast reservoir (1/T)
%   T_tr:par(15) temperature threshold (°C)
%   f_dd:par(16) degree-day factor (mm/°C/T)
%  mdl          structure with model state/parameter info
%   .mcode       scalar numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation ode_xinanjiang C++
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 5
    check = 0;          % no check of J(x)_f and J(x)_th
end
mcode = mdl.mcode;      % Formulation/language
                        % 1: Runge Kutta implementation MATLAB
                        % 2: ode45 implementation MATLAB
                        % 3: Explicit Euler int_steps MATLAB
                        % 4: Runge Kutta implementation ode_xinanjiang C++
mem = ode.mem;          % state variable storage or not
if mcode == 2 && mem == 0
    warning(['xinanjiang: ' ...
        'built-in ode45 ' ...
        'solver stores ' ...
        'states: mem = 1'])
    mem = 1;
end
T_sm = 1.0;                 % smoothing width (°C) partition/positive-part
eps_m = 1e-6;               % smoothing for min()
eps_r = 1e-12;              % smoothing for regimes
m = 9;                      % # state variables [+ snow: 1 state = SWE]
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

f_wm = th(5);           % Fraction of S_tot that is W_max
f_lm = th(6);           % Fraction of W_max that is LM
S_tot = th(8);          % Total soil moisture storage (W+S) (mm)
W_max = f_wm*S_tot;     % Maximum tension water depth (mm)
S_max = (1-f_wm)*S_tot; % Maximum free water depth (mm)
LM = f_lm*W_max;        % Tension water threshold evaporation change (mm)

% Initialization
switch mcode
    case {0,1,2,3}
        %% Unpack parameters
        f_p = th(1);    % Ratio of potential evaporion to pan evaporation
        A_im = th(2);   % 0-0.01: Fraction impervious area (-)
        a = th(3);      % Tension distribution inflection parameter (-)
        b = th(4);      % 0.1-0.4: Tension distribution shape parameter (-)
        c = th(7);      % 0.1-0.3: Fraction LM secnd evaporation change (-)
        Ex = th(9);     % 0.5 - 2.5: Free distribution shape parameter (-)
        k_i = th(10);   % 0.2-0.5: Free water interflow parameter (1/T)
        k_g = th(11);   % 0.1-0.4: Free water groundwater parameter (1/T)
                        % Some papers use: ki+kg = 0.7;
        c_i = th(12);   % 0.1-0.99: Interflow time coefficient (1/T)
        c_g = th(13);   % 0.7-0.99: Baseflow time coefficient (1/T)
        k_f = th(14);   % Recession constant fast reservoir (1/T)
                        % --> to compute channel inflow
        T_tr = th(15);  % Temperature threshold (°C)
        f_dd = th(16);  % Degree-day factor (mm/°C/T)

    case 4
        data.f_p = th(1);   % Ratio potential evaporion to pan evaporation
        data.A_im = th(2);  % Fraction impervious area (-)
        data.a = th(3);     % Tension water distr. inflection parameter (-)
        data.b = th(4);     % Tension water distr. shape parameter (-)
        data.W_max = W_max; % Maximum tension water storage (mm)
        data.LM = LM;       % Tension water threshold evaprtion change (mm)
        data.c = th(7);     % Fraction of LM second evaporation change (-)
        data.S_max = S_max; % Maximum free water depth (mm)
        data.Ex = th(9);    % Free water distribution shape parameter (-)
        data.k_i = th(10);  % Free water interflow parameter (1/T)
        data.k_g = th(11);  % Free water groundwater parameter (1/T)
        data.c_i = th(12);  % Interflow time coefficient (1/T)
        data.c_g = th(13);  % Baseflow time coefficient (1/T)
        data.k_f = th(14);  % Recession constant fast reservoir (1/T)
                            % --> to compute channel inflow
        data.S_tot = S_tot; % Total soil moisture storage (W+S) (mm)
        data.f_wm = f_wm;   % Fraction of Stot that is Wmax (-)
        data.f_lm = f_lm;   % Fraction of wmax that is LM (-)
        data.T_tr = th(15); % Temperature threshold (°C)
        data.f_dd = th(16); % Degree-day factor (mm/°C/T)

        data.T_sm = T_sm;   % smoothing width (°C) partition/positive-part
        data.eps_m = eps_m; % smoothing for min (°C)
        data.eps_r = eps_r; % smoothing for regimes
        data.eps = eps;     % Dimensionless smoothing coef. [inactive]
        data.rho = rho;     % Dimensionless smoothing coef. [inactive]
        data.ipr = ipr;     % Time to print

end

% Execute model
switch mcode

    case 0 %% TEST
        hin = ode.InitStep;     % Initial time step
        hmax_ = ode.MaxStep;    % Maximum time step
        hmin_ = ode.MinStep;    % minimum time step
        reltol = ode.RelTol;    % Relative tolerance
        abstol = ode.AbsTol;    % Absolute tolerance
        order = ode.Order;      % Order
        maxiter = ode.maxiter;  % Maximum iterations
        flag = 0;
        fail = false;

        % Tell ODE how to index forcing (needed so ode45 works)
        % 0 = print-interval forcing (it)
        % 1 = time-based forcing (floor(t)+1),
        data.ForceByPrIdx = 0;
        EPS = 2.220446049250313e-16;

        hCarry = max(hmin_,min(hin,hmax_)); % carry adaptive recommendation

        for s = 2:ns
            t1 = s-2;
            t2 = s-1;
            % C++ uses P[s-1],Ep[s-1],T[s-1] fixed over this print interval
            data.it = s-1;           % 1-based index into data.it
            h = min(hCarry,t2-t1);   % reuse prior recommendation
            if mem == 1
                Z(s,1:nvar) = Z(s-1,1:nvar);
                z = Z(s,1:nvar);
            else
                z = Z(1,1:nvar);
            end
            tcur = t1; iterCount = 0;
            while tcur < t2
                tleft = t2 - tcur;      % always limit to remaining time
                h = min(h, tleft);
                % C++ snap tolerance
                tTol = 10 * EPS * max(1, max(abs(t2), abs(tcur)));
                if tleft <= tTol
                    break;
                end
                if h < tTol
                    % "step underflow" fail
                    fail = true;
                    break;
                end
                ztmp = Z(s,1:nvar);
                % one RK2 step (updates ztmp, returns LTE)
                [ztmp,LTE] = rk2(tcur,ztmp,h,f_p,A_im,a,b, ...
                    W_max,LM,c,S_max,S_tot,f_wm,f_lm,Ex, ...
                    k_i,k_g,c_i,c_g,k_f,T_tr,f_dd,data, ...
                    T_sm,eps_m,eps_r,eps,rho,m,d);
                % sanity checks
                if any(~isfinite(ztmp)) || any(abs(ztmp) > 1e12)
                    fail = true;
                    break;
                end
                w = 1 ./ (reltol * abs(ztmp) + abstol);
                wl = w .* LTE;
                wrms = sqrt(sum(wl.^2)/nvar);
                if ~isfinite(wrms), wrms = Inf; end
                accepted = (wrms <= 1.0) || (h <= hmin_);
                if accepted
                    Z(s,1:nvar) = ztmp;
                    tcur = tcur + h;
                    iterCount = 0;          % reset on progress
                else
                    iterCount = iterCount + 1;
                end
                % fac update
                if wrms <= 0
                    fac = 5.0;
                else
                    fac = 0.9 * wrms^(-1/order);
                    fac = min(5.0, max(0.2, fac));
                end
                hNext = h * fac;
                hNext = max(hmin_, min(hNext, hmax_));
                hCarry = hNext;                         % retain before boundary clipping
                h = min(hNext, t2 - tcur);
                % max reject check
                if iterCount >= maxiter
                    if flag == 0
                        warning("WARNING: Max rejection limit " + ...
                            "reached at t = %.5f\n", tcur);
                        warning("WARNING: Parameter f_p   %8.5f\n", f_p);
                        warning("WARNING: Parameter A_im  %8.5f\n", A_im);
                        warning("WARNING: Parameter a     %8.5f\n", a);
                        warning("WARNING: Parameter b     %8.5f\n", b);
                        warning("WARNING: Parameter W_max %8.5f\n", W_max);
                        warning("WARNING: Parameter LM    %8.5f\n", LM);
                        warning("WARNING: Parameter c     %8.5f\n", c);
                        warning("WARNING: Parameter S_max %8.5f\n", S_max);
                        warning("WARNING: Parameter S_tot %8.5f\n", S_tot);
                        warning("WARNING: Parameter f_wm  %8.5f\n", f_wm);
                        warning("WARNING: Parameter f_lm  %8.5f\n", f_lm);
                        warning("WARNING: Parameter Ex    %8.5f\n", Ex);
                        warning("WARNING: Parameter k_i   %8.5f\n", k_i);
                        warning("WARNING: Parameter k_g   %8.5f\n", k_g);
                        warning("WARNING: Parameter c_i   %8.5f\n", c_i);
                        warning("WARNING: Parameter c_g   %8.5f\n", c_g);
                        warning("WARNING: Parameter k_f   %8.5f\n", k_f);
                        warning("WARNING: Parameter T_tr  %8.5f\n", T_tr);
                        warning("WARNING: Parameter f_dd  %8.5f\n", f_dd);
                        flag = 1;
                    end
                    fail = true;
                    break;
                end
            end
            if fail
                break;
            else
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
        end

        % Return to default forcing mode for safety
        data.ForceByPrIdx = 1;

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
            t1 = s-2; t2 = s-1;                     % Set start, end times
            h = min(hCarry,t2-t1);              % reuse prior recommendation
            h = min(hCarry,t2-t1);              % reuse prior recommendation
            if mem == 1
                Z(s,1:nvar) = Z(s-1,1:nvar);
                z = Z(s,1:nvar);            % row vector
            else
                z = Z(1,1:nvar);
            end
            t = t1;                         % Initial time
            % Integrate from t1 to t2
            while (t < t2)
                [ztmp,LTE] = rk2(t,z,h,f_p,A_im,a,b, ...     % Evaluate rk2
                    W_max,LM,c,S_max,S_tot,f_wm,f_lm,Ex, ...
                    k_i,k_g,c_i,c_g,k_f,T_tr,f_dd,data, ...
                    T_sm,eps_m,eps_r,eps,rho,m,d);
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
                        warning("WARNING: Parameter f_p   %8.5f\n", f_p);
                        warning("WARNING: Parameter A_im  %8.5f\n", A_im);
                        warning("WARNING: Parameter a     %8.5f\n", a);
                        warning("WARNING: Parameter b     %8.5f\n", b);
                        warning("WARNING: Parameter W_max %8.5f\n", W_max);
                        warning("WARNING: Parameter LM    %8.5f\n", LM);
                        warning("WARNING: Parameter c     %8.5f\n", c);
                        warning("WARNING: Parameter S_max %8.5f\n", S_max);
                        warning("WARNING: Parameter S_tot %8.5f\n", S_tot);
                        warning("WARNING: Parameter f_wm  %8.5f\n", f_wm);
                        warning("WARNING: Parameter f_lm  %8.5f\n", f_lm);
                        warning("WARNING: Parameter Ex    %8.5f\n", Ex);
                        warning("WARNING: Parameter k_i   %8.5f\n", k_i);
                        warning("WARNING: Parameter k_g   %8.5f\n", k_g);
                        warning("WARNING: Parameter c_i   %8.5f\n", c_i);
                        warning("WARNING: Parameter c_g   %8.5f\n", c_g);
                        warning("WARNING: Parameter k_f   %8.5f\n", k_f);
                        warning("WARNING: Parameter T_tr  %8.5f\n", T_tr);
                        warning("WARNING: Parameter f_dd  %8.5f\n", f_dd);
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

    case 2 %% MATLAB: ode45 implementation of xinanjiang_odefcn
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
        [~,Z] = ode45(@(t,z) xinanjiang_aug_ode(t,z,f_p, ...
            A_im,a,b,W_max,LM,c,S_max,S_tot,f_wm,f_lm,Ex, ...
            k_i,k_g,c_i,c_g,k_f,T_tr,f_dd,data,T_sm,eps_m, ...
            eps_r,eps,rho,m,d),0:mdl.tout,Z(1,1:nvar),ode_options);
        s = size(Z,1); if s < ns, fail = true; s = s+1; end

    case 3 %% MATLAB: Explicit Euler xinanjiang_odefcn with int_steps steps
        int_steps = 500;                        % number of integration steps
        dt = 1/int_steps;                       % time of each int. step
        for s = 2:ns                            % Start time loop
            if mem == 1
                z = Z(s-1,1:nvar);              % Initialize state variables
            else
                z = Z(1,1:nvar);
            end
            for it = 1:int_steps                % Do integration in int_steps steps
                dzdt = xinanjiang_aug_ode(s-2,z,f_p,A_im, ...
                    a,b,W_max,LM,c,S_max,S_tot,f_wm,f_lm, ...
                    Ex,k_i,k_g,c_i,c_g,k_f,T_tr,f_dd,data, ...
                    T_sm,eps_m,eps_r,eps,rho,m,d)';
                if s == 1540 && check
                    [~,~,Jth_f,Jx_f] = rhs_only(s-2,z,th,data, ...
                        T_sm,eps_m,eps_r,eps,rho,m,d);
                    % --- Jth_f via 2-sided finite differences ---
                    Jth_fn = zeros(m,d);
                    for j = 1:d
                        h = 1e-6*max(1,abs(th(j)));
                        thp = th; thm = th;
                        thp(j) = thp(j) + h; thm(j) = thm(j) - h;
                        fp = rhs_only(s-2,z,thp,data,T_sm, ...
                            eps_m,eps_r,eps,rho,m,d);
                        fm = rhs_only(s-2,z,thm,data,T_sm, ...
                            eps_m,eps_r,eps,rho,m,d);
                        Jth_fn(:,j) = (fp - fm)/(2*h);
                    end
                    % --- Jx_f via 2-sided finite differences ---
                    Jx_fn = zeros(m,m);
                    for i = 1:m
                        h = 1e-6*max(1,abs(z(i)));
                        zp = z; zm = z;
                        zp(i) = zp(i) + h; zm(i) = zm(i) - h;
                        fp = rhs_only(s-2,zp,th,data,T_sm,eps_m, ...
                            eps_r,eps,rho,m,d);
                        fm = rhs_only(s-2,zm,th,data,T_sm,eps_m, ...
                            eps_r,eps,rho,m,d);
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

    case 4 %% C++: Runge Kutta implementation ( = similar to xinanjiang_odefcn)
        if nargout == 1
            if mem == 1
                [Z,~] = crr_xinanjiang(mdl.tout,Z(1,1:nvar)',data,ode);
                q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
            else
                [~,q_n] = crr_xinanjiang(mdl.tout,Z(1,1:nvar)',data,ode);
            end
            J = [];
        else
            [Z,q_n,J] = crr_xinanjiang(mdl.tout,Z(1,1:nvar)',data,ode);
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
function [z,LTE] = rk2(t,z,h,f_p,A_im,a,b,W_max,LM,c,S_max,S_tot, ...
    f_wm,f_lm,Ex,k_i,k_g,c_i,c_g,k_f,T_tr,f_dd,data,T_sm,eps_m, ...
    eps_r,eps,rho,m,d)

dzdtE = xinanjiang_aug_ode(t,z,f_p,A_im,a,b,W_max,LM,c,S_max,S_tot, ...
    f_wm,f_lm,Ex,k_i,k_g,c_i,c_g,k_f,T_tr,f_dd,data,T_sm,eps_m,eps_r, ...
    eps,rho,m,d)';                                      % Euler
zE = z + h*dzdtE;
dzdtH = xinanjiang_aug_ode(t,zE,f_p,A_im,a,b,W_max,LM,c,S_max,S_tot, ...
    f_wm,f_lm,Ex,k_i,k_g,c_i,c_g,k_f,T_tr,f_dd,data,T_sm,eps_m,eps_r, ...
    eps,rho,m,d)';                                      % Heun
z = z + 0.5*h*(dzdtE + dzdtH);                          % New z
LTE = abs(zE - z);                                      % LTE

end

%% 2. XINANJIANG augmented ode with sensitivities as state variables
function dzdt = xinanjiang_aug_ode(t,z,f_p,A_im,a,b,W_max,LM,c, ...
    S_max,S_tot,f_wm,f_lm,Ex,k_i,k_g,c_i,c_g,k_f,T_tr,f_dd,data, ...
    T_sm,eps_m,eps_r,eps,rho,m,d)

x = z(1:m);                                         % x = [s_we w s_w s_i
                                                    %      s_g s_f1 s_f2
                                                    %      s_f3 r_q]
S = reshape(z(m+1:m*(d+1)),m,d);                    % sensitivity matrix
[dxdt,dSdt] = xinanjiang_odefcn(t,x,S,f_p,A_im, ... % compute dxdt & dSdt
    a,b,W_max,LM,c,S_max,S_tot,f_wm,f_lm,Ex,k_i,k_g, ...
    c_i,c_g,k_f,T_tr,f_dd,data,T_sm,eps_m,eps_r,eps,rho,m,d);
dzdt = [dxdt; dSdt(:)];                             % repack single vector

end

%% 3. XINANJIANG: Secondary function, ODE solver
function [dxdt,dSdt,Jth_f,Jx_f] = xinanjiang_odefcn(t,x,Smat,f_p,A_im, ...  % compute dxdt & dSdt
    a,b,W_max,LM,c,S_max,S_tot,f_wm,f_lm,Ex,k_i,k_g,c_i,c_g,k_f,T_tr,f_dd, ...
    data,T_sm,eps_m,eps_r,~,~,m,d)

%% Input arguments
% eps = ~ and rho = ~

%% Define return argument
dxdt = nan(m,1);        % Initialize return argument

%% Unpack states and forcing
Swe = x(1);             % Snow water equivalent (mm)
W = x(2);               % Tension water storage (mm)
Sw = x(3);              % Free water storage (mm)
Si = x(4);              % Storage interflow routing reservoir (mm)
Sg = x(5);              % Storage baseflow routing reservoir (mm)
Sf = x(6:8);            % Storage 1st, 2nd & 3rd fast reservoir (rout, mm)

it = floor(t) + 1;      % Truncate time to current time index
P = data.P(it,1);       % Instantaneous rainfall (mm/T)
Ep = data.Ep(it,1);     % Instantaneous Ep (mm/T)
T = data.T(it,1);       % Get current temperature (°C)

smooth_pos = @(a,eps) 0.5*(a + sqrt(a.^2 + eps.^2));
d_smooth_pos_da = @(a,eps) 0.5*(1 + a./sqrt(a.^2 + eps.^2));

Swe_u = Swe;  W_u = W;  Sw_u = Sw;

Swe = smooth_pos(Swe_u,eps_m);
W = smooth_pos(W_u,eps_m);
Sw = smooth_pos(Sw_u,eps_m);

dSwe_dx1 = d_smooth_pos_da(Swe_u,eps_m);
dW_dx2 = d_smooth_pos_da(W_u,eps_m);
dSw_dx3 = d_smooth_pos_da(Sw_u,eps_m);

% --------------------------------------------
% 3) Snow module (smooth HBV-style degree-day)
% --------------------------------------------
T_smeps_m = max(T_sm,eps_m);
u = (T - T_tr)/T_smeps_m;
% Smooth snow fraction: ~1 when T<T_tr, ~0 when T>T_tr
snow_fr = 0.5*(1 - tanh(u));
rain_fr = 1 - snow_fr;
P_snow = P * snow_fr;
P_rain = P * rain_fr;

% Smooth positive part: posT = max(T-T_tr,0) (smooth)
aT = (T - T_tr);
posT = 0.5*(aT + sqrt(aT^2 + T_smeps_m^2)); % >=0, smooth
M_pot = f_dd * posT;                        % potential melt (mm/T)

% Smooth min(SWE, Mpot)
dxy = (Swe - M_pot);
sqrtm = sqrt(dxy^2 + eps_m^2);
M = 0.5*(Swe + M_pot - sqrtm);          % actual melt (mm/T)
Pliq = P_rain + M;                      % Into unsaturated reservoir (mm/T)

% -----------------------------------------------------
% 2) Snow derivatives needed for coupling (x and theta)
% -----------------------------------------------------
% dsnow_fr/dT_tr
sech2 = 1/(cosh(u)^2);
dsnow_dT_tr = 0.5 * sech2 / T_smeps_m;            % d(snow_fr)/dT_tr

% posT derivatives wrt T_tr
da_dT_tr = -1;
dposT_da = 0.5*(1 + aT/sqrt(aT^2 + T_smeps_m^2));
dposT_dT_tr = dposT_da * da_dT_tr;

% Mpot derivatives
dMpot_dT_tr = f_dd * dposT_dT_tr;
dMpot_f_dd = posT;

% M derivatives (smooth min)
dM_dSwe = 0.5*(1 - dxy/sqrtm);
dM_dMpot = 0.5*(1 + dxy/sqrtm);

% Melt derivatives wrt parameters
dM_dT_tr = dM_dMpot * dMpot_dT_tr;
dM_f_dd = dM_dMpot * dMpot_f_dd;

% Liquid precip derivatives
dPliq_dSwe = dM_dSwe;
dPliq_dT_tr = (-P * dsnow_dT_tr) + dM_dT_tr;
dPliq_f_dd = dM_f_dd;

% ------------------------------------------------------------
% 3) Compute fluxes between layers (Xinanjiang), but P -> Pliq
% ------------------------------------------------------------
Pi = (1-A_im)*Pliq;             % Infiltration (mm/T)
Rb = A_im*Pliq;                 % Direct runoff impervious

dPi_dPliq = (1 - A_im);
dRb_dPliq = A_im;

Ea = f_p * Ep;                  % Pan evaporation (mm/T)

% --- Runoff generation R(Pi,W) ---  (SMOOTH BLEND, no hard clamps)
uW = W / W_max;                 % W already smoothed
duW_dW = 1 / W_max;
duW_dWmax = -W / W_max^2;

% Smooth switch weight between lower/upper regimes around uW = 0.5 - a
g = uW - (0.5 - a);            % = uW - 0.5 + a
gt = g / eps_r;

w = 0.5*(1 + tanh(gt));       % 0..1
sech2 = 1/(cosh(gt)^2);
dw_dg = 0.5 * sech2 / eps_r;

dg_dW = duW_dW;
dg_da = 1;
dg_dWmax = duW_dWmax;

% --------------------------------------------------------------
% Lower branch: Rlow = Pi*c1*uW^b with smooth_pos(uW) protection
% --------------------------------------------------------------
base1 = max(0.5 - a, 1e-12);
c1 = base1^(1 - b);
dc1_da = -(1 - b) * base1^(-b);
dc1_db = -c1 * log(base1);

uWe = max(uW, eps_m);           % match C++: fmax(uW,eps_m)
duWe_duW = 1.0;                % C++ treats clamp derivative as 1
uWb = uWe^b;
uWbm1 = uWe^(b-1);
loguW = log(uWe);

Rlow = Pi * c1 * uWb;
dRlow_dPi = c1 * uWb;
dRlow_dW = Pi * c1 * b * uWbm1 * duWe_duW * duW_dW;
dRlow_da = Pi * dc1_da * uWb;              % uW does not depend on a
dRlow_db = Pi * (dc1_db * uWb + c1 * uWb * loguW);
dRlow_dWmax = Pi * c1 * b * uWbm1 * duWe_duW * duW_dWmax;

% -------------------------------------------------------------
% Upper branch: Rup = Pi*(1 - c2*v^b), v = max(1-uW,0) smoothly
% -------------------------------------------------------------
base2 = max(0.5 + a, 1e-12);   % should not need as a in [-0.4999,0.4999]
c2 = base2^(1 - b);
dc2_da = (1 - b) * base2^(-b);
dc2_db = -c2 * log(base2);

v_raw = 1 - uW;
v = smooth_pos(v_raw, eps_m);
dv_dvraw = d_smooth_pos_da(v_raw, eps_m);

dvraw_dW = -duW_dW;
dvraw_dWmax = -duW_dWmax;

dv_dW = dv_dvraw * dvraw_dW;
dv_dWmax = dv_dvraw * dvraw_dWmax;

ve = max(v, eps_m);              % match C++: fmax(v,eps_m)
vb = ve^b;
vbm1 = ve^(b-1);
logv = log(ve);

Rup = Pi * (1 - c2 * vb);
dRup_dPi = 1 - c2 * vb;
dRup_dW = Pi * (-c2 * b * vbm1) * dv_dW;
dRup_da = Pi * (-dc2_da * vb);
dRup_db = Pi * (-(dc2_db * vb + c2 * vb * logv));
dRup_dWmax = Pi * (-c2 * b * vbm1) * dv_dWmax;

% ------------
% Smooth blend
% ------------
R = (1-w)*Rlow + w*Rup;

dR_dPi = (1-w)*dRlow_dPi + w*dRup_dPi;
dR_dW = (1-w)*dRlow_dW + w*dRup_dW + dw_dg*dg_dW*(Rup - Rlow);
dR_da = (1-w)*dRlow_da + w*dRup_da + dw_dg*dg_da*(Rup - Rlow);
dR_db = (1-w)*dRlow_db + w*dRup_db;
dR_dWmax = (1-w)*dRlow_dWmax + w*dRup_dWmax + dw_dg*dg_dWmax*(Rup - Rlow);

dR_dPliq = dR_dPi * dPi_dPliq;

% --- Evaporation piece ---
if W > LM
    E = Ea;
    dE_df_p = Ep;
    dE_dLM = 0.0;
    dE_dc = 0.0;
    dE_dW = 0.0;
elseif W >= c*LM
    E = (W/LM) * Ea;
    dE_df_p = (W/LM) * Ep;
    dE_dLM = -(W/LM^2) * Ea;
    dE_dc = 0.0;
    dE_dW = Ea / LM;
else % W < c*LM
    E0 = c * Ea;
    W_phi = 1;
    % Water-limiting factor near W = 0
    phiW = W / (W + W_phi);             % ~1 for big W, ~0 near W=0
    dphi_dW = W_phi / (W + W_phi)^2;
    % Scaled evaporation
    E = E0 * phiW;
    % Derivatives
    dE_dW = E0 * dphi_dW;               % since E0 independent of W here
    dE_dc = Ea * phiW;                  % scale by phiW
    dE_df_p = (c * Ep) * phiW;          % scale by phiW (Ea depends on f_p)
    dE_dLM = 0.0;
end

% --- Partition runoff into Rs/Ri/Rg (hard z clamp) ---
z = 1 - Sw / S_max;
if z < 0, z = 0; end

dum = 1 - z^Ex;

Rs = R * dum;
Ri = k_i * Sw * dum;
Rg = k_g * Sw * dum;

Qi = c_i * Si;
Qg = c_g * Sg;
Qs = Rs + Rb;

q_fout = k_f * Sf;

% ddum_dS = 0.0;
% ddum_dSmax = 0.0;
% ddum_dEx = 0.0;

% if (z > 0) && (Sw < S_max)
%     ddum_dS = Ex * z^(Ex-1) * (1.0 / S_max);
%     ddum_dSmax = -Ex * z^(Ex-1) * (Sw / (S_max * S_max));
%     ddum_dEx = - (z^Ex) * log(z);
% end
sig = @(x) 0.5*(1 + tanh(x));

g1 = sig(z/eps_r);                % ~1 if z>0, ~0 if z<0
g2 = sig((S_max - Sw)/eps_r);     % ~1 if Sw<S_max, ~0 if Sw>S_max
g12 = g1 .* g2;                    % combined smooth gate

% Base derivatives (no gate)
ddS_base = Ex * z^(Ex-1) * (1.0 / S_max);
ddSmax_base = -Ex * z^(Ex-1) * (Sw / (S_max * S_max));
ddEx_base = -(z^Ex) * log(max(z,eps_m));

% Apply smooth gate
ddum_dS = g12 * ddS_base;
ddum_dSmax = g12 * ddSmax_base;
ddum_dEx = g12 * ddEx_base;

% ----------------------------------
% 4) Net flux into/out of reservoirs
% ----------------------------------
dxdt(1) = P_snow - M;                   % dSwe/dt
dxdt(2) = Pi - E - R;                   % dW/dt
dxdt(3) = R - Rs - Ri - Rg;             % dSw/dt
dxdt(4) = Ri - Qi;                      % dSi/dt
dxdt(5) = Rg - Qg;                      % dSg/dt
dxdt(6) = Qs + Qi + Qg - q_fout(1);     % dSf1/dt
dxdt(7) = q_fout(1) - q_fout(2);        % dSf2/dt
dxdt(8) = q_fout(2) - q_fout(3);        % dSf3/dt
dxdt(9) = q_fout(3);                    % dQ/dt

% ------------------------------------------------------------
% 5) Flux derivatives wrt states (existing + new SWE coupling)
% ------------------------------------------------------------
dRs_dW = dR_dW * dum;
dRs_dS = R * ddum_dS;

dRi_dS = k_i * (dum + Sw * ddum_dS);
dRg_dS = k_g * (dum + Sw * ddum_dS);

dQi_dSi = c_i;
dQg_dSg = c_g;

dQs_dW = dRs_dW;
dQs_dS = dRs_dS;

% NEW: derivatives wrt Pliq (for SWE coupling)
dPi_dSwe = dPi_dPliq * dPliq_dSwe;
dPi_dT_tr = dPi_dPliq * dPliq_dT_tr;
dPi_f_dd = dPi_dPliq * dPliq_f_dd;

dRb_dSwe = dRb_dPliq * dPliq_dSwe;
dRb_dT_tr = dRb_dPliq * dPliq_dT_tr;
dRb_f_dd = dRb_dPliq * dPliq_f_dd;

dR_dSwe = dR_dPliq  * dPliq_dSwe;
dR_dT_tr = dR_dPliq  * dPliq_dT_tr;
dR_f_dd = dR_dPliq  * dPliq_f_dd;

dRs_dSwe = dR_dSwe * dum;
dQs_dSwe = dRs_dSwe + dRb_dSwe;

% ------------------------------
% 6) Build Jx_f = df/dx  (m x m)
% ------------------------------
Jx_f = zeros(m,m);

% Row 1: dSwe/dt = P_snow - M
Jx_f(1,1) = -dM_dSwe;

% Row 2: dW/dt = Pi - E - R
Jx_f(2,1) = dPi_dSwe - dR_dSwe;     % SWE coupling through Pliq
Jx_f(2,2) = -dE_dW - dR_dW;         % wrt W

% Row 3: dSw/dt = R - Rs - Ri - Rg
Jx_f(3,1) = dR_dSwe - dRs_dSwe;     % SWE coupling
Jx_f(3,2) = dR_dW - dRs_dW;         % wrt W
Jx_f(3,3) = -dRs_dS - dRi_dS - dRg_dS;

% Row 4: dSi/dt = Ri - Qi
Jx_f(4,3) = dRi_dS;
Jx_f(4,4) = -dQi_dSi;

% Row 5: dSg/dt = Rg - Qg
Jx_f(5,3) = dRg_dS;
Jx_f(5,5) = -dQg_dSg;

% Row 6: dSf1/dt = Qs + Qi + Qg - qf1
Jx_f(6,1) = dQs_dSwe;       % SWE coupling
Jx_f(6,2) = dQs_dW;
Jx_f(6,3) = dQs_dS;
Jx_f(6,4) = dQi_dSi;
Jx_f(6,5) = dQg_dSg;
Jx_f(6,6) = -k_f;

% Row 7: dSf2/dt = qf1 - qf2
Jx_f(7,6) = k_f;
Jx_f(7,7) = -k_f;

% Row 8: dSf3/dt = qf2 - qf3
Jx_f(8,7) = k_f;
Jx_f(8,8) = -k_f;

% Row 9: dQ/dt = qf3
Jx_f(9,8) =  k_f;

Jx_f(:,1) = Jx_f(:,1) * dSwe_dx1;
Jx_f(:,2) = Jx_f(:,2) * dW_dx2;
Jx_f(:,3) = Jx_f(:,3) * dSw_dx3;

% -----------------------------------
% 7) Build Jth_f = df/dtheta  (m x d)
% -----------------------------------
Jth_f = zeros(m,d);

% Precompute derivatives used repeatedly
dRs_dSmax = R * ddum_dSmax;
dRi_dSmax = k_i * Sw * ddum_dSmax;
dRg_dSmax = k_g * Sw * ddum_dSmax;
dRs_dEx = R * ddum_dEx;
dRi_dEx = k_i * Sw * ddum_dEx;
dRg_dEx = k_g * Sw * ddum_dEx;

% --- theta(1) = f_p: PET -> pan ET coefficient
j = 1;
Jth_f(2,j) = - dE_df_p;   % only dx1 depends on f_p

% --- theta(2) = A_im: impervious fraction
dPi_dAim = -Pliq;
dRb_dAim = Pliq;
dR_dAim = dR_dPi * dPi_dAim;
dRs_dAim = dR_dAim * dum;
dQs_dAim = dRs_dAim + dRb_dAim;
j = 2;
Jth_f(2,j) = dPi_dAim - dR_dAim;   % df1 = Pi - E - R
Jth_f(3,j) = dR_dAim - dRs_dAim;   % df2 = R - Rs - Ri - Rg
Jth_f(6,j) = dQs_dAim;             % df5 = Qs + Qi + Qg - qf1

% --- theta(3) = a
dRs_da = dR_da * dum;
j = 3;
Jth_f(2,j) = - dR_da;
Jth_f(3,j) = dR_da - dRs_da;
Jth_f(6,j) = dRs_da;

% --- theta(4) = b
dRs_db = dR_db * dum;
j = 4;
Jth_f(2,j) = - dR_db;
Jth_f(3,j) = dR_db - dRs_db;
Jth_f(6,j) = dRs_db;

% --- theta(5) but this is Wmax
J_Wmax = zeros(m,1);
dRs_dWmax = dR_dWmax * dum;
J_Wmax(2) = - dR_dWmax;
J_Wmax(3) = dR_dWmax - dRs_dWmax;
J_Wmax(6) = dRs_dWmax;

% --- theta(6) but this is LM
J_LM = zeros(m,1);
J_LM(2) = - dE_dLM;

% --- theta(7) = c (second evaporation threshold)
j = 7;
Jth_f(2,j) = - dE_dc;

% --- theta(8) but this is Smax
J_Smax = zeros(m,1);
J_Smax(3) = - dRs_dSmax - dRi_dSmax - dRg_dSmax;    % df2 = R - Rs - Ri - Rg
J_Smax(4) = dRi_dSmax;                              % df3 = Ri - Qi
J_Smax(5) = dRg_dSmax;                              % df4 = Rg - Qg
J_Smax(6) = dRs_dSmax;                              % df5 = Qs + Qi + Qg - qf1

% --- theta(9) = Ex
j = 9;
Jth_f(3,j) = - dRs_dEx - dRi_dEx - dRg_dEx;
Jth_f(4,j) = dRi_dEx;
Jth_f(5,j) = dRg_dEx;
Jth_f(6,j) = dRs_dEx;

% --- theta(10) = ki
j = 10;
dRi_dki = Sw * dum;
Jth_f(3,j) = - dRi_dki;     % in dx2
Jth_f(4,j) =  dRi_dki;      % in dx3

% --- theta(11) = kg
j = 11;
dRg_dkg = Sw * dum;
Jth_f(3,j) = - dRg_dkg;     % add to existing
Jth_f(5,j) = dRg_dkg;

% --- theta(12) = ci
j = 12;
Jth_f(4,j) = -Si;           % dx3 = Ri - Qi
Jth_f(6,j) =  Si;           % dx5 = ... + Qi ...

% --- theta(13) = cg
j = 13;
Jth_f(5,j) = -Sg;           % dx4 = Rg - Qg
Jth_f(6,j) = Sg;            % dx5 = ... + Qg ...

% --- theta(14) = kf
j = 14;
Jth_f(6,j) = - Sf(1);       % from -qf1
Jth_f(7,j) = Sf(1) - Sf(2);
Jth_f(8,j) = Sf(2) - Sf(3);
Jth_f(9,j) = Sf(3);

% --------------------------------------------------------
% 8) Add snow parameter columns (assumed last two columns)
% --------------------------------------------------------

% Row 1: dSwe/dt = P_snow - M
Jth_f(1,d-1) = (P * dsnow_dT_tr) - dM_dT_tr;
Jth_f(1,d) = -dM_f_dd;

% Row 2: dW/dt = Pi - E - R (snow affects Pi and R through Pliq)
Jth_f(2,d-1) = dPi_dT_tr - dR_dT_tr;
Jth_f(2,d) = dPi_f_dd - dR_f_dd;

% Row 3: dSw/dt = R - Rs - Ri - Rg  (snow affects R and Rs through Pliq)
% Rs depends on R, so use the same paT_trern you used for SWE:
dRs_dT_tr = dR_dT_tr * dum;
dRs_f_dd = dR_f_dd * dum;

Jth_f(3,d-1) = dR_dT_tr - dRs_dT_tr;
Jth_f(3,d) = dR_f_dd - dRs_f_dd;

% Row 6: dSf1/dt includes Qs = Rs + Rb (snow affects both via Pliq)
dQs_dT_tr = dRs_dT_tr  + dRb_dT_tr;
dQs_f_dd = dRs_f_dd  + dRb_f_dd;

Jth_f(6,d-1) = dQs_dT_tr;
Jth_f(6,d) = dQs_f_dd;

% ---------------------------------------
%  5. chain rule for f_wm, f_lm and S_tot
% ---------------------------------------

Jth_f(:,5) = S_tot * J_Wmax + (f_lm*S_tot) * J_LM - S_tot * J_Smax;
Jth_f(:,6) = W_max * J_LM;
Jth_f(:,8) = f_wm * J_Wmax + (f_lm*f_wm) * J_LM + (1 - f_wm) * J_Smax;

% -------------------------------------------------
%  6. Sensitivity update: dS/dt = Jf_x*S + Jf_theta
% -------------------------------------------------
dSdt = Jx_f * Smat + Jth_f;         % (mxm)*(mxd) + (mxd) = mxd

end

%% 4. XINANJIANG: RHS-only for J(x)_f and J(th)_f check
function [dxdt,dSdt,Jth_f,Jx_f] = rhs_only(t,z,th,data,T_sm,eps_m, ...
    eps_r,eps,rho,m,d)

% Unpack once
f_p = th(1);            % Ratio of potential evaporion to pan evaporation
A_im = th(2);           % 0-0.01: Fraction impervious area (-)
a = th(3);              % Tension distribution inflection parameter (-)
b = th(4);              % 0.1-0.4: Tension distribution shape parameter (-)
f_wm = th(5);           % Fraction of Stot that is Wmax (-)
f_lm = th(6);           % Fraction of wmax that is LM (-)
S_tot = th(8);          % Total soil moisture storage (W+S) (mm)
W_max = f_wm*S_tot;     % 100-200: Maximum tension water depth (mm)
S_max = (1-f_wm)*S_tot; % 10-50: Maximum free water depth (mm)
LM = f_lm*W_max;        % Tension threshold for evaporation change (mm)
c = th(7);              % 0.1-0.3: Fraction LM secnd evaporation change (-)
Ex = th(9);             % 0.5 - 2.5: Free distribution shape parameter (-)
k_i = th(10);           % 0.2-0.5: Free water interflow parameter (1/T)
k_g = th(11);           % 0.1-0.4: Free water groundwater parameter (1/T)
                        % Some papers use: ki+kg = 0.7;
c_i = th(12);           % 0.1-0.99: Interflow time coefficient (1/T)
c_g = th(13);           % 0.7-0.99: Baseflow time coefficient (1/T)
k_f = th(14);           % Recession constant fast reservoir (1/T)
                        % --> to compute channel inflow
T_tr = th(15);          % Temperature threshold (°C)
f_dd = th(16);          % Degree-day factor (mm/°C/T)

Smat = ones(m,d);
% Call your scalar-argument RHS
[dxdt,dSdt,Jth_f,Jx_f] = xinanjiang_odefcn(t,z,Smat,f_p,A_im,a,b,W_max, ...
    LM,c,S_max,S_tot,f_wm,f_lm,Ex,k_i,k_g,c_i,c_g,k_f,T_tr,f_dd, ...
    data,T_sm,eps_m,eps_r,eps,rho,m,d);

end
