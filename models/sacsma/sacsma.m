function varargout = sacsma(par,mdl,data,ode,check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SACSMA: Runge Kutta implementation of sacsma conceptual watershed model
% SYNOPSIS: varargout = sacsma(par,mdl,data,ode,check)
%  par          dx1 vector of parameter values
%   uzfwm:par(1) maximum free water storage of upper zone (mm)
%   uztwm:par(2) maximum tension water storage of upper zone (mm)
%   lzfpm:par(3) maximum free water storage of lower zone primary (mm)
%   lzfsm:par(4) maximum free water storage of lower zone secundary (mm)
%   lztwm:par(5) maximum tension water storage of lower zone (mm)
%   zperc:par(6) multiplier percolation function (-)
%   rexp:par(7)  power of percolation function (-)
%   uzk:par(8)   interflow rate (1/T)
%   pfree:par(9) fraction percolation to tension storage lower layer (-)
%   lzpk:par(10) base flow depletion rate for primary reservoir (1/T)
%   lzsk:par(11) base flow depletion rate for secondary reservoir (1/T)
%   acm:par(12)  maximum fraction of saturated area (-)
%   kf:par(13)   recession constant fast reservoir (1/T)
%                --> to compute channel inflow
%   T_tr:par(14) temperature threshold (°C)
%   f_dd:par(15) degree-day factor (mm/°C/T)
%  mdl          structure with model state/parameter info
%   .mcode       scalar numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation ode_sacsma C++
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
                        % 4: Runge Kutta implementation ode_sacsma C++
mem = ode.mem;          % state variable storage or not
if mcode == 2 && mem == 0
    warning(['sacsma: ' ...
        'built-in ode45 ' ...
        'solver stores ' ...
        'states: mem = 1'])
    mem = 1;
end
T_sm = 1.0;                 % smoothing width (°C) partition/positive-part
eps_m = 1e-6;               % smoothing for min()
eps_s = 1e-12;              % smoothing for state variables
m = 10;                     % # state variables [+ snow: 1 state = SWE]
ns = mdl.tout + 1;          % # print times
d = numel(par);             % # parameters
n = mdl.idx(2)-mdl.idx(1);  % # elements discharge vector
nvar = m*(d+1);             % # number of variables
eps = 5.0;                  % Dimensionless smoothing coefficient
rho = 0.01;                 % Dimensionless smoothing coefficient
fail = false;               % Default: model completes run
data.ForceByPrIdx = 1;      % Test for this model only
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

% Initialization
switch mcode
    case {0,1,2,3}
        %% Unpack parameters
        uzfwm = th(1);  % Maximum free storage upper zone (mm)
        uztwm = th(2);  % Maximum tension storage upper zone (mm)
        lzfpm = th(3);  % Maximum free storage lower zone primary (mm)
        lzfsm = th(4);  % Maximum free storage lower zone secundary (mm)
        lztwm = th(5);  % Maximum tension water storage lower zone (mm)
        zperc = th(6);  % Multiplier percolation function (-)
        rexp = th(7);   % Power of percolation function (-)
        uzk = th(8);    % Interflow rate (1/T)
        pfree = th(9);  % Fraction percolation to tension lower layer (-)
        lzpk = th(10);  % Base flow depletion rate primary reservoir (1/T)
        lzsk = th(11);  % Base flow depletion rate secndary reservoir (1/T)
        acm = th(12);   % Maximum fraction of saturated area (-)
        kf = th(13);    % Recession constant fast reservoir (1/T)
                        % --> to compute channel inflow
        T_tr = th(14);  % Temperature threshold (°C)
        f_dd = th(15);  % Degree-day factor (mm/°C/T)

    case 4
        data.uzfwm = th(1); % Maximum free storage upper zone (mm)
        data.uztwm = th(2); % Maximum tension storage upper zone (mm)
        data.lzfpm = th(3); % Maximum free storage lower zone primary (mm)
        data.lzfsm = th(4); % Maximum free storage lower zone secndry (mm)
        data.lztwm = th(5); % Maximum tension storage lower zone (mm)
        data.zperc = th(6); % Multiplier percolation function (-)
        data.rexp = th(7);  % Power of percolation function (-)
        data.uzk = th(8);   % Interflow rate (1/T)
        data.pfree = th(9); % Fraction percolation tension lower layer (-)
        data.lzpk = th(10); % Base flow depletion rate prim reservoir (1/T)
        data.lzsk = th(11); % Base flow depletion rate secd reservoir (1/T)
        data.acm = th(12);  % Maximum fraction of saturated area (-)
        data.kf = th(13);   % Recession constant fast reservoir (1/T)
                            % --> to compute channel inflow
        data.T_tr = th(14); % Temperature threshold (°C)
        data.f_dd = th(15); % Degree-day factor (mm/°C/T)

        data.T_sm = T_sm;   % smoothing width (°C) partition/positive-part
        data.eps_m = eps_m; % smoothing for min (°C)
        data.eps_s = eps_s; % smoothing for state variables
        data.eps = eps;     % Dimensionless smoothing coefficient
        data.rho = rho;     % Dimensionless smoothing coefficient
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

        % Tell ODE how to index forcing (needed so ode45 still works)
        % 0 = print-interval forcing (data.it)
        % 1 = time-based forcing (floor(t)+1),
        data.ForceByPrIdx = 0;
        EPS = 2.220446049250313e-16;

        hCarry = max(hmin_,min(hin,hmax_)); % carry adaptive recommendation

        for s = 2:ns
            t1 = s-2; t2 = s-1;
            % C++ uses P[s-1],Ep[s-1],T[s-1] fixed over this print interval
            data.it = s-1;                      % 1-based index into data.it
            h = min(hCarry,t2-t1);              % reuse prior recommendation
            if mem == 1
                Z(s,1:nvar) = Z(s-1,1:nvar);
                z = Z(s,1:nvar);
            else
                z = Z(1,1:nvar);
            end
            tcur = t1; iterCount = 0;
            while tcur < t2
                tleft = t2 - tcur;          % always limit to remaining time
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
                % one RK2 step (updates ztmp, returns LTE)
                [ztmp,LTE] = rk2(tcur, z, h, uzfwm, uztwm, ...
                    lzfpm, lzfsm, lztwm, zperc, rexp, uzk, ...
                    pfree, lzpk, lzsk, acm, kf, T_tr, f_dd, data, ...
                    T_sm, eps_m, eps_s, eps, rho, m, d);
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
                    z = ztmp;
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
                        warning("WARNING: Parameter uzfwm %8.5f\n", uzfwm);
                        warning("WARNING: Parameter uztwm %8.5f\n", uztwm);
                        warning("WARNING: Parameter lzfpm %8.5f\n", lzfpm);
                        warning("WARNING: Parameter lzfsm %8.5f\n", lzfsm);
                        warning("WARNING: Parameter lztwm %8.5f\n", lztwm);
                        warning("WARNING: Parameter zperc %8.5f\n", zperc);
                        warning("WARNING: Parameter rexp  %8.5f\n", rexp);
                        warning("WARNING: Parameter uzk   %8.5f\n", uzk);
                        warning("WARNING: Parameter pfree %8.5f\n", pfree);
                        warning("WARNING: Parameter lzpk  %8.5f\n", lzpk);
                        warning("WARNING: Parameter lzsk  %8.5f\n", lzsk);
                        warning("WARNING: Parameter acm   %8.5f\n", acm);
                        warning("WARNING: Parameter kf    %8.5f\n", kf);
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
        % Return to default forcing mode for safety (optional)
        data.ForceByPrIdx = 1;

    case 1 %% MATLAB: Runge Kutta implementation (= model similar to hymod_odefcn)
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
                z = Z(s,1:nvar);                    % row vector
            else
                z = Z(1,1:nvar);
            end
            t = t1;                                 % Initial time
            % Integrate from t1 to t2
            while (t < t2)
                [ztmp,LTE] = rk2(t,z,h,uzfwm,uztwm, ...  % Evaluate rk2
                    lzfpm,lzfsm,lztwm,zperc,rexp,uzk, ...
                    pfree,lzpk,lzsk,acm,kf,T_tr,f_dd,data, ...
                    T_sm,eps_m,eps_s,eps,rho,m,d);
                if any(~isfinite(ztmp)) || any(abs(ztmp) > 1e12)
                    fail = true; break;
                end
                w = 1 ./ (reltol*abs(ztmp(1:nvar)) + abstol);   % Weights
                wrms = sqrt(sum((w.*LTE).^2)/nvar);             % CORRECTED, Nov. 2022
                if (wrms <= 1) || (h <= hmin_)                  % Accept if error is small enough
                    z = ztmp; t = t + h;
                end
                hNext = h*max(0.2,min(5.0,0.9*wrms^(-1/order)));    % new step
                hNext = max(hmin_,min(hNext,hmax_));
                hCarry = hNext;                         % retain before boundary clipping
                h = min(hNext,t2-t);   % Another turn
                if (iterCount >= maxiter)
                    if (flag == 0)
                        warning("WARNING: Max step limit " + ...
                            "reached at t = %.5f\n", t);
                        warning("WARNING: Parameter uzfwm %8.5f\n", uzfwm);
                        warning("WARNING: Parameter uztwm %8.5f\n", uztwm);
                        warning("WARNING: Parameter lzfpm %8.5f\n", lzfpm);
                        warning("WARNING: Parameter lzfsm %8.5f\n", lzfsm);
                        warning("WARNING: Parameter lztwm %8.5f\n", lztwm);
                        warning("WARNING: Parameter zperc %8.5f\n", zperc);
                        warning("WARNING: Parameter rexp  %8.5f\n", rexp);
                        warning("WARNING: Parameter uzk   %8.5f\n", uzk);
                        warning("WARNING: Parameter pfree %8.5f\n", pfree);
                        warning("WARNING: Parameter lzpk  %8.5f\n", lzpk);
                        warning("WARNING: Parameter lzsk  %8.5f\n", lzsk);
                        warning("WARNING: Parameter acm   %8.5f\n", acm);
                        warning("WARNING: Parameter kf    %8.5f\n", kf);
                        warning("WARNING: Parameter T_tr %8.5f\n", T_tr);
                        warning("WARNING: Parameter f_dd %8.5f\n", f_dd);
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

    case 2 %% MATLAB: ode45 implementation of sacsma_odefcn
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
        [~,Z] = ode45(@(t,z) sacsma_aug_ode(t,z,uzfwm,uztwm, ...
            lzfpm,lzfsm,lztwm,zperc,rexp,uzk,pfree,lzpk,lzsk, ...
            acm,kf,T_tr,f_dd,data,T_sm,eps_m,eps_s,eps,rho,m,d), ...
            0:mdl.tout,Z(1,1:nvar),ode_options);
        s = size(Z,1); if s < ns, fail = true; s = s+1; end

    case 3 %% MATLAB: Explicit Euler sacsma_odefcn with int_steps steps
        int_steps = 1000;                       % number of integration steps
        dt = 1/int_steps;                       % time of each int. step
        for s = 2:ns                            % Start time loop
            if mem == 1
                z = Z(s-1,1:nvar);              % Initialize state variables
            else
                z = Z(1,1:nvar);
            end
            for it = 1:int_steps                % Integrate int_steps steps
                dzdt = sacsma_aug_ode(s-2,z,uzfwm, ...
                    uztwm,lzfpm,lzfsm,lztwm,zperc, ...
                    rexp,uzk,pfree,lzpk,lzsk,acm, ...
                    kf,T_tr,f_dd,data,T_sm,eps_m, ...
                    eps_s,eps,rho,m,d)';
                if s == 150 && check
                    [~,~,Jth_f,Jx_f] = rhs_only(s-2,z,th,data, ...
                        T_sm,eps_m,eps_s,eps,rho,m,d);
                    % --- Jth_f via 2-sided finite differences ---
                    Jth_fn = zeros(m,d);
                    for j = 1:d
                        h = 1e-6*max(1,abs(th(j)));
                        thp = th; thm = th;
                        thp(j) = thp(j) + h; thm(j) = thm(j) - h;
                        fp = rhs_only(s-2,z,thp,data,T_sm, ...
                            eps_m,eps_s,eps,rho,m,d);
                        fm = rhs_only(s-2,z,thm,data,T_sm, ...
                            eps_m,eps_s,eps,rho,m,d);
                        Jth_fn(:,j) = (fp - fm)/(2*h);
                    end
                    % --- Jx_f via 2-sided finite differences ---
                    Jx_fn = zeros(m,m);
                    for i = 1:m
                        h = 1e-6*max(1,abs(z(i)));
                        zp = z; zm = z;
                        zp(i) = zp(i) + h; zm(i) = zm(i) - h;
                        fp = rhs_only(s-2,zp,th,data,T_sm,eps_m, ...
                            eps_s,eps,rho,m,d);
                        fm = rhs_only(s-2,zm,th,data,T_sm,eps_m, ...
                            eps_s,eps,rho,m,d);
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
                if mem == 1
                    Z(s,1:nvar) = z;            % State at t
                else
                    if tpr(s) == 1
                        q_n(s-ds) = z(m) - Z(m);
                        J(s-ds,1:d) = z(id) - Z(id);
                    end
                    Z(1,1:nvar) = z;
                end
            end
        end                                     % End of time loop

    case 4 %% C++: Runge Kutta implementation ( = similar to sacsma_odefcn)
        if nargout == 1
            if mem == 1
                [Z,~] = crr_sacsma(mdl.tout,Z(1,1:nvar)',data,ode);
                q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
            else
                [~,q_n] = crr_sacsma(mdl.tout,Z(1,1:nvar)',data,ode);
            end
            J = [];
        else
            [Z,q_n,J] = crr_sacsma(mdl.tout,Z(1,1:nvar)',data,ode);
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
function [z,LTE] = rk2(t,z,h,uzfwm,uztwm,lzfpm,lzfsm,lztwm, ...
    zperc,rexp,uzk,pfree,lzpk,lzsk,acm,kf,T_tr,f_dd,data,T_sm, ...
    eps_m,eps_s,eps,rho,m,d)

dzdtE = sacsma_aug_ode(t,z,uzfwm,uztwm,lzfpm,lzfsm,lztwm,...
    zperc,rexp,uzk,pfree,lzpk,lzsk,acm,kf,T_tr,f_dd,data,T_sm, ...
    eps_m,eps_s,eps,rho,m,d)';                                  % Euler
zE = z + h*dzdtE;
dzdtH = sacsma_aug_ode(t,zE,uzfwm,uztwm,lzfpm,lzfsm,lztwm,...
    zperc,rexp,uzk,pfree,lzpk,lzsk,acm,kf,T_tr,f_dd,data,T_sm, ...
    eps_m,eps_s,eps,rho,m,d)';                                  % Heun
z = z + 0.5*h*(dzdtE + dzdtH);                                  % New z
LTE = abs(zE - z);                                              % LTE

end

%% 2. SACSMA augmented ode with sensitivities as state variables
function dzdt = sacsma_aug_ode(t,z,uzfwm,uztwm,lzfpm,lzfsm,lztwm,...
    zperc,rexp,uzk,pfree,lzpk,lzsk,acm,kf,T_tr,f_dd,data,T_sm, ...
    eps_m,eps_s,eps,rho,m,d)

x = z(1:m);                                         % x = [Swe Uztw Uzfw 
                                                    %      Lztw Lzps Lzfs 
                                                    %      Sf1 Sf2 Sf3 Q]
S = reshape(z(m+1:m*(d+1)),m,d);                    % sensitivity matrix
[dxdt,dSdt] = sacsma_odefcn(t,x,S,uzfwm,uztwm, ...  % compute dxdt & dSdt
    lzfpm,lzfsm,lztwm,zperc,rexp,uzk,pfree, ...
    lzpk,lzsk,acm,kf,T_tr,f_dd,data,T_sm,eps_m, ...
    eps_s,eps,rho,m,d);
dzdt = [dxdt; dSdt(:)];                             % repack single vector

end

%% 3. SACSMA: Secondary function, ODE solver + sensitivities
function [dxdt,dSdt,Jth_f,Jx_f] = sacsma_odefcn( ...
    t,x,Smat,uzfwm,uztwm,lzfpm,lzfsm,lztwm,zperc,rexp,uzk,pfree, ...
    lzpk,lzsk,acm,kf,T_tr,f_dd,data,T_sm,eps_m,eps_s,eps,rho,m,d)
%SACSMA_ODEFCN_SMOOTH Fully smoothed SAC-SMA(+snow) ODE + sensitivities
%
% States (m=10):
%   x = [Swe_raw, Uztw, Uzfw, Lztw, Lzps, Lzfs, Sf1, Sf2, Sf3, Sq]'
%
% Parameters (d=15):
%   1 uzfwm, 2 uztwm, 3 lzfpm, 4 lzfsm, 5 lztwm, 6 zperc, 7 rexp,
%   8 uzk, 9 pfree, 10 lzpk, 11 lzsk, 12 acm, 13 kf, 14 T_tr, 15 f_dd
%
% Smoothing:
%   - Swe = smooth_pos(Swe_raw)
%   - Uzfw_pos = smooth_pos(Uzfw), in B and all its derivatives
%   - Uztw_pos = smooth_pos(Uztw), for ET and ac (removes if/else kinks)
%   - Lztw_pos = smooth_pos(Lztw), for ET lower
%   - Aeff = smooth_pos(Araw) or softplus-like via smooth_pos(Araw) 
%     (keeps dlz smooth at Araw=0)
%   - All "min/max" style switches removed. Logistic overflow is smooth.

%% helpers (C1 smoothing)
smooth_pos = @(a,ee) 0.5*(a + sqrt(a.^2 + ee.^2));          % ~max(a,0)
dsmooth_pos_da = @(a,ee) 0.5*(1 + a./sqrt(a.^2 + ee.^2));       % d/d a

smooth_min = @(a,b,ee) 0.5*(a + b - sqrt((a-b).^2 + ee.^2));
dsmoothmin_da = @(a,b,ee) 0.5*(1 - (a-b)./sqrt((a-b).^2 + ee.^2));
dsmoothmin_db = @(a,b,ee) 0.5*(1 + (a-b)./sqrt((a-b).^2 + ee.^2));

% smooth clip to [0,1] using two smooth_pos operations
smooth_clip01 = @(y,ee) 1 - smooth_pos(1 - smooth_pos(y,ee), ee); % ~min(1,max(y,0))

% logistic overflow (original smoothing)
lf_func = @(S,Smax) 1.0./(1.0 + exp(- ((S - (Smax - rho * Smax * eps)) ./ (rho * Smax))));
dlf_dS = @(lf_val,Smax) (lf_val .* (1.0 - lf_val)) ./ (rho * Smax);

%% unpack
dxdt = nan(m,1);

Swe_raw = x(1);
Uztw = x(2);
Uzfw = x(3);
Lztw = x(4);
Lzps = x(5);
Lzfs = x(6);
Sf1 = x(7);
Sf2 = x(8);
Sf3 = x(9);

if data.ForceByPrintIndex == 1
    it = floor(t) + 1;
else
    it = data.it;
end
P = data.P(it,1);
Ep = data.Ep(it,1);
T = data.T(it,1);

%% fully smoothed "positive" states where needed
Swe = smooth_pos(Swe_raw, eps_s);
dSwe_dRaw = dsmooth_pos_da(Swe_raw, eps_s);

Uztw_pos = smooth_pos(Uztw, eps_s);
dUztwpos_dUztw = dsmooth_pos_da(Uztw, eps_s);

Uzfw_pos = smooth_pos(Uzfw, eps_s);
dUzfwpos_dUzfw = dsmooth_pos_da(Uzfw, eps_s);

Lztw_pos = smooth_pos(Lztw, eps_s);
dLztwpos_dLztw = dsmooth_pos_da(Lztw, eps_s);

%% aggregates
Lztot = Lztw + Lzps + Lzfs;
Lzm = lzfpm + lzfsm + lztwm;

% ------------------------------------------------------------------
% 1) Snow module (fully smoothed, as you had)
% ------------------------------------------------------------------
T_smeps_m = max(T_sm, eps_m);
uT = (T - T_tr) / T_smeps_m;

snow_fr = 0.5 * (1.0 - tanh(uT));
rain_fr = 1.0 - snow_fr;

P_snow = P * snow_fr;
P_rain = P * rain_fr;

aT = (T - T_tr);
posT = smooth_pos(aT, T_smeps_m);   % = 0.5*(aT + sqrt(aT^2+T_smeps_m^2))

P_pmelt = f_dd * posT;

P_amelt = smooth_min(Swe, P_pmelt, eps_m);
P_liq = P_rain + P_amelt;

% derivatives for snow
coshu = cosh(uT);
sech2 = 1.0 ./ (coshu .* coshu);
dsnow_dT_tr = 0.5 * sech2 / T_smeps_m;

dposT_da = dsmooth_pos_da(aT, T_smeps_m);
dposT_dT_tr = -dposT_da;

dP_pmelt_dT_tr = f_dd * dposT_dT_tr;
dP_pmelt_f_dd = posT;

dP_amelt_dSwe_s = dsmoothmin_da(Swe, P_pmelt, eps_m);   % dP_amelt/dSwe
dP_amelt_dSwe = dP_amelt_dSwe_s * dSwe_dRaw;          % chain to Swe_raw

dP_amelt_dP_pmelt = dsmoothmin_db(Swe, P_pmelt, eps_m);% dP_amelt/dP_pmelt

dP_amelt_dT_tr = dP_amelt_dP_pmelt * dP_pmelt_dT_tr;
dP_amelt_f_dd = dP_amelt_dP_pmelt * dP_pmelt_f_dd;

dP_liq_dSwe = dP_amelt_dSwe;
dP_liq_dT_tr = (-P * dsnow_dT_tr) + dP_amelt_dT_tr;
dP_liq_f_dd = dP_amelt_f_dd;

% ------------------------------------------------
% 2) Fluxes (FULLY SMOOTHED: remove if/else kinks)
% ------------------------------------------------

% ET upper: e1 = Ep * min(Uztw,uztwm)/uztwm, but smooth and nonnegative
ratioU = Uztw_pos / uztwm;                              % >=0
ratioU_clip = smooth_clip01(ratioU, eps_s);             % ~min(1,ratioU)
e_1 = Ep * ratioU_clip;
% derivative wrt Uztw: Ep * d(ratioU_clip)/dUztw
% chain: ratioU depends on Uztw_pos and dUztwpos/dUztw
% smooth_clip01(y) = 1 - smooth_pos(1 - smooth_pos(y))
% d/dy:  dsmooth_pos(y)*dsmooth_pos(1 - smooth_pos(y))  (with signs)
y = ratioU;
sp1 = smooth_pos(y, eps_s);
dsp1 = dsmooth_pos_da(y, eps_s);
%sp2 = smooth_pos(1 - sp1, eps_s);
dsp2 = dsmooth_pos_da(1 - sp1, eps_s);
dratio_clip_dy = dsp2 * dsp1;                           % der. clip01 wrt y
de1_dUztw = Ep * dratio_clip_dy * (dUztwpos_dUztw / uztwm);
de1_duztwm = Ep * dratio_clip_dy * (-Uztw_pos / (uztwm*uztwm));  % treat Uztw_pos const wrt uztwm

% ET lower: e2 = (Ep - e1) * min(Lztw,lztwm)/lztwm, smooth and nonnegative
ratioL = Lztw_pos / lztwm;
y = ratioL;
sp1 = smooth_pos(y, eps_s); dsp1 = dsmooth_pos_da(y, eps_s);
%sp2 = smooth_pos(1 - sp1, eps_s);
dsp2 = dsmooth_pos_da(1 - sp1, eps_s);
dratioL_clip_dy = dsp2 * dsp1;
cL = smooth_clip01(ratioL, eps_s);
e_2 = (Ep - e_1) * cL;

dcL_dLztw = dratioL_clip_dy * (dLztwpos_dLztw / lztwm);
dcL_dlztwm = dratioL_clip_dy * (-Lztw_pos / (lztwm*lztwm));

de2_dLztw = (Ep - e_1) * dcL_dLztw;
de2_dlztwm = (Ep - e_1) * dcL_dlztwm;
de2_dUztw = -cL * de1_dUztw;
de2_duztwm = -cL * de1_duztwm;

% baseflow at saturation
q_0 = lzpk * lzfpm + lzsk * lzfsm;

% dlz factor: Araw = 1 - Lztot/Lzm. Fully smooth A = max(Araw,0)
Araw = 1.0 - Lztot / Lzm;
Apos = smooth_pos(Araw, eps_s);                 % ~max(Araw,0)
dApos_dAraw = dsmooth_pos_da(Araw, eps_s);

Afloor = 1e-12;
Aeff = smooth_pos(Apos - Afloor, eps_s) + Afloor;
Apow = Aeff.^rexp;
dlz = 1.0 + zperc * Apow;
% derivatives of dlz
ddlz_dzperc = Apow;
ddlz_drexp = zperc .* Apow .* log(Aeff);

% dAraw/d states and Lzm-params
dAraw_dLztw = -(1/Lzm);
dAraw_dLzps = -(1/Lzm);
dAraw_dLzfs = -(1/Lzm);

dAraw_dlzfpm = (Lztot)/(Lzm*Lzm);
dAraw_dlzfsm = (Lztot)/(Lzm*Lzm);
dAraw_dlztwm = (Lztot)/(Lzm*Lzm);

% chain Apos
dApos_dLztw = dApos_dAraw * dAraw_dLztw;
dApos_dLzps = dApos_dAraw * dAraw_dLzps;
dApos_dLzfs = dApos_dAraw * dAraw_dLzfs;

dApos_dlzfpm = dApos_dAraw * dAraw_dlzfpm;
dApos_dlzfsm = dApos_dAraw * dAraw_dlzfsm;
dApos_dlztwm = dApos_dAraw * dAraw_dlztwm;

% derivative of Apow wrt Apos: rexp * Aeff^(rexp-1) * dAeff/dApos
dAeff_dApos = dsmooth_pos_da(Apos - Afloor, eps_s);
dApow_dApos = rexp .* (Aeff.^(rexp-1)) .* dAeff_dApos;

ddlz_dLztw = zperc * dApow_dApos * dApos_dLztw;
ddlz_dLzps = zperc * dApow_dApos * dApos_dLzps;
ddlz_dLzfs = zperc * dApow_dApos * dApos_dLzfs;

ddlz_dlzfpm = zperc * dApow_dApos * dApos_dlzfpm;
ddlz_dlzfsm = zperc * dApow_dApos * dApos_dlzfsm;
ddlz_dlztwm = zperc * dApow_dApos * dApos_dlztwm;

% B uses Uzfw_pos
B = Uzfw_pos / uzfwm;
dB_dUzfw = dUzfwpos_dUzfw / uzfwm;
dB_duzfwm = -Uzfw_pos / (uzfwm*uzfwm);

q_12 = q_0 * dlz * B;
dq12_dUzfw = q_0 * dlz * dB_dUzfw;

dq12_dLztw = q_0 * B * ddlz_dLztw;
dq12_dLzps = q_0 * B * ddlz_dLzps;
dq12_dLzfs = q_0 * B * ddlz_dLzfs;

dq12_dzperc = q_0 * B * ddlz_dzperc;
dq12_drexp = q_0 * B * ddlz_drexp;

q_if = uzk * B;
dqif_dUzfw = uzk * dB_dUzfw;

q_bp = lzpk * Lzps;
q_bs = lzsk * Lzfs;

dqbp_dLzps = lzpk;
dqbp_dlzpk = Lzps;
dqbs_dLzfs = lzsk;
dqbs_dlzsk = Lzfs;

% ac uses Uztw_pos (smooth)
ac = acm * (Uztw_pos / uztwm);
q_sx = ac * P_liq;

dac_dUztw = acm * (dUztwpos_dUztw / uztwm);
dac_duztwm = acm * (-Uztw_pos / (uztwm*uztwm));
dac_dacm = Uztw_pos / uztwm;

dqsx_dUztw = P_liq * dac_dUztw;
dqsx_duztwm = P_liq * dac_duztwm;
dqsx_dacm = P_liq * dac_dacm;

dqsx_dSwe = ac * dP_liq_dSwe;
dqsx_dT_tr = ac * dP_liq_dT_tr;
dqsx_f_dd = ac * dP_liq_f_dd;

% logistic factors
lf_Uztw = lf_func(Uztw,uztwm);
lf_Uzfw = lf_func(Uzfw,uzfwm);
lf_Lztw = lf_func(Lztw,lztwm);
lf_Lzps = lf_func(Lzps,lzfpm);
lf_Lzfs = lf_func(Lzfs,lzfsm);

dlf_Uztw = dlf_dS(lf_Uztw,uztwm);
dlf_Uzfw = dlf_dS(lf_Uzfw,uzfwm);
dlf_Lztw = dlf_dS(lf_Lztw,lztwm);
dlf_Lzps = dlf_dS(lf_Lzps,lzfpm);
dlf_Lzfs = dlf_dS(lf_Lzfs,lzfsm);

q_utof = (P_liq - q_sx) * lf_Uztw;
q_ufof = q_utof * lf_Uzfw;

q_stof = pfree * q_12 * lf_Lztw;
prc_s = 0.5*(1.0 - pfree)*q_12 + 0.5*q_stof;

q_sfofa = prc_s * lf_Lzps;
q_sfofb = prc_s * lf_Lzfs;

q_fout1 = kf * Sf1;
q_fout2 = kf * Sf2;
q_fout3 = kf * Sf3;

% ------
% 3) ODE
% ------
dxdt(1) = P_snow - P_amelt;
dxdt(2) = P_liq - q_sx - e_1 - q_utof;
dxdt(3) = q_utof - q_12 - q_if - q_ufof;
dxdt(4) = pfree * q_12 - e_2 - q_stof;
dxdt(5) = prc_s - q_bp - q_sfofa;
dxdt(6) = prc_s - q_bs - q_sfofb;
dxdt(7) = q_if + q_sx + q_ufof + q_sfofa + q_sfofb + q_bp + q_bs - q_fout1;
dxdt(8) = q_fout1 - q_fout2;
dxdt(9) = q_fout2 - q_fout3;
dxdt(10) = q_fout3;

% ----------------------------------
% 4) Jacobian pieces (states/params)
% ----------------------------------

% q_utof derivatives
dqutof_dSwe = (dP_liq_dSwe - dqsx_dSwe) * lf_Uztw;
dqutof_dT_tr = (dP_liq_dT_tr  - dqsx_dT_tr ) * lf_Uztw;
dqutof_f_dd = (dP_liq_f_dd  - dqsx_f_dd ) * lf_Uztw;

dqutof_dUztw = -(dqsx_dUztw) * lf_Uztw + (P_liq - q_sx) * dlf_Uztw;
dqutof_dacm = -(dqsx_dacm ) * lf_Uztw;

% dlf/dSmax terms (for param derivatives) stay as your identities
dlfU_dUztwm = - dlf_Uztw * Uztw / uztwm;
dqutof_duztwm = -(dqsx_duztwm) * lf_Uztw + (P_liq - q_sx) * dlfU_dUztwm;

% q_ufof derivatives
dqufof_dSwe = dqutof_dSwe  * lf_Uzfw;
dqufof_dT_tr = dqutof_dT_tr   * lf_Uzfw;
dqufof_f_dd = dqutof_f_dd   * lf_Uzfw;

dqufof_dUztw = dqutof_dUztw * lf_Uzfw;
dqufof_dUzfw = q_utof * dlf_Uzfw;
dqufof_dacm = dqutof_dacm  * lf_Uzfw;

% q_stof derivatives
dqstof_dUzfw = pfree * dq12_dUzfw * lf_Lztw;
dqstof_dLztw = pfree * dq12_dLztw * lf_Lztw + pfree * q_12 * dlf_Lztw;
dqstof_dLzps = pfree * dq12_dLzps * lf_Lztw;
dqstof_dLzfs = pfree * dq12_dLzfs * lf_Lztw;

dqstof_dzperc = pfree * dq12_dzperc * lf_Lztw;
dqstof_drexp = pfree * dq12_drexp  * lf_Lztw;
dqstof_dpfree = q_12 * lf_Lztw;

% prc_s derivatives
dprc_dUzfw = 0.5*(1.0-pfree)*dq12_dUzfw + 0.5*dqstof_dUzfw;
dprc_dLztw = 0.5*(1.0-pfree)*dq12_dLztw + 0.5*dqstof_dLztw;
dprc_dLzps = 0.5*(1.0-pfree)*dq12_dLzps + 0.5*dqstof_dLzps;
dprc_dLzfs = 0.5*(1.0-pfree)*dq12_dLzfs + 0.5*dqstof_dLzfs;

dprc_dzperc = 0.5*(1.0-pfree)*dq12_dzperc + 0.5*dqstof_dzperc;
dprc_drexp = 0.5*(1.0-pfree)*dq12_drexp  + 0.5*dqstof_drexp;
dprc_dpfree = -0.5*q_12 + 0.5*dqstof_dpfree;

% q_sfofa derivatives
dqsfpf_dUzfw = dprc_dUzfw * lf_Lzps;
dqsfpf_dLztw = dprc_dLztw * lf_Lzps;
dqsfpf_dLzps = dprc_dLzps * lf_Lzps + prc_s * dlf_Lzps;
dqsfpf_dLzfs = dprc_dLzfs * lf_Lzps;

dqsfpf_dzperc = dprc_dzperc * lf_Lzps;
dqsfpf_drexp = dprc_drexp  * lf_Lzps;
dqsfpf_dpfree = dprc_dpfree * lf_Lzps;

% q_sfofb derivatives
dqsfss_dUzfw = dprc_dUzfw * lf_Lzfs;
dqsfss_dLztw = dprc_dLztw * lf_Lzfs;
dqsfss_dLzps = dprc_dLzps * lf_Lzfs;
dqsfss_dLzfs = dprc_dLzfs * lf_Lzfs + prc_s * dlf_Lzfs;

dqsfss_dzperc = dprc_dzperc * lf_Lzfs;
dqsfss_drexp = dprc_drexp  * lf_Lzfs;
dqsfss_dpfree = dprc_dpfree * lf_Lzfs;

% fast reservoirs
dqfout1_dSf1 = kf;
dqfout2_dSf2 = kf;
dqfout3_dSf3 = kf;

dqfout1_dkf = Sf1;
dqfout2_dkf = Sf2;
dqfout3_dkf = Sf3;

% -------------------------
% Parameter-specific pieces
% -------------------------

% uzfwm
dq12_duzfwm = q_0 * dlz * dB_duzfwm;
dqif_duzfwm = uzk * dB_duzfwm;

dqstof_duzfwm = pfree * lf_Lztw * dq12_duzfwm;
dprc_duzfwm = 0.5*(1.0-pfree)*dq12_duzfwm + 0.5*dqstof_duzfwm;

dqsfpf_duzfwm = dprc_duzfwm * lf_Lzps;
dqsfss_duzfwm = dprc_duzfwm * lf_Lzfs;

dlfUzfw_duzfwm = -dlf_Uzfw * Uzfw / uzfwm;        % your identity
dqufof_duzfwm = q_utof * dlfUzfw_duzfwm;

% uztwm
dqufof_duztwm = dqutof_duztwm * lf_Uzfw;

% lzfpm
dq0_dlzfpm = lzpk;
dq12_dlzfpm = dq0_dlzfpm * dlz * B + q_0 * ddlz_dlzfpm * B;

dqstof_dlzfpm = pfree * lf_Lztw * dq12_dlzfpm;
dprc_dlzfpm = 0.5*(1.0-pfree)*dq12_dlzfpm + 0.5*dqstof_dlzfpm;

dlfLzps_dlzfpm = -dlf_Lzps * Lzps / lzfpm;
dqsfpf_dlzfpm = dprc_dlzfpm * lf_Lzps + prc_s * dlfLzps_dlzfpm;
dqsfss_dlzfpm = dprc_dlzfpm * lf_Lzfs;

% lzfsm
dq0_dlzfsm = lzsk;
dq12_dlzfsm = dq0_dlzfsm * dlz * B + q_0 * ddlz_dlzfsm * B;

dqstof_dlzfsm = pfree * lf_Lztw * dq12_dlzfsm;
dprc_dlzfsm = 0.5*(1.0-pfree)*dq12_dlzfsm + 0.5*dqstof_dlzfsm;

dlfLzfs_dlzfsm = -dlf_Lzfs * Lzfs / lzfsm;
dqsfpf_dlzfsm = dprc_dlzfsm * lf_Lzps;
dqsfss_dlzfsm = dprc_dlzfsm * lf_Lzfs + prc_s * dlfLzfs_dlzfsm;

% lztwm
dq12_dlztwm = q_0 * ddlz_dlztwm * B;

dlfLztw_dlztwm = -dlf_Lztw * Lztw / lztwm;
dqstof_dlztwm = pfree * (dq12_dlztwm * lf_Lztw + q_12 * dlfLztw_dlztwm);
dprc_dlztwm = 0.5*(1.0-pfree)*dq12_dlztwm + 0.5*dqstof_dlztwm;

dqsfpf_dlztwm = dprc_dlztwm * lf_Lzps;
dqsfss_dlztwm = dprc_dlztwm * lf_Lzfs;

% ---------------------
% 5) Build Jx_f (m x m)
% ---------------------
Jx_f = zeros(m,m);

Jx_f(1,1) = -dP_amelt_dSwe;

Jx_f(2,1) = dP_liq_dSwe - dqsx_dSwe - dqutof_dSwe;
Jx_f(2,2) = -dqsx_dUztw - de1_dUztw - dqutof_dUztw;

Jx_f(3,1) = dqutof_dSwe - dqufof_dSwe;
Jx_f(3,2) = dqutof_dUztw - dqufof_dUztw;
Jx_f(3,3) = -dqif_dUzfw - dq12_dUzfw - dqufof_dUzfw;
Jx_f(3,4) = -dq12_dLztw;
Jx_f(3,5) = -dq12_dLzps;
Jx_f(3,6) = -dq12_dLzfs;

Jx_f(4,2) = -de2_dUztw;
Jx_f(4,3) = pfree*dq12_dUzfw - dqstof_dUzfw;
Jx_f(4,4) = pfree*dq12_dLztw - de2_dLztw - dqstof_dLztw;
Jx_f(4,5) = pfree*dq12_dLzps - dqstof_dLzps;
Jx_f(4,6) = pfree*dq12_dLzfs - dqstof_dLzfs;

Jx_f(5,3) = dprc_dUzfw - dqsfpf_dUzfw;
Jx_f(5,4) = dprc_dLztw - dqsfpf_dLztw;
Jx_f(5,5) = dprc_dLzps - dqbp_dLzps - dqsfpf_dLzps;
Jx_f(5,6) = dprc_dLzfs - dqsfpf_dLzfs;

Jx_f(6,3) = dprc_dUzfw - dqsfss_dUzfw;
Jx_f(6,4) = dprc_dLztw - dqsfss_dLztw;
Jx_f(6,5) = dprc_dLzps - dqsfss_dLzps;
Jx_f(6,6) = dprc_dLzfs - dqbs_dLzfs - dqsfss_dLzfs;

Jx_f(7,1) = dqsx_dSwe + dqufof_dSwe;
Jx_f(7,2) = dqsx_dUztw + dqufof_dUztw;
Jx_f(7,3) = dqif_dUzfw + dqufof_dUzfw + dqsfpf_dUzfw + dqsfss_dUzfw;
Jx_f(7,4) = dqsfpf_dLztw + dqsfss_dLztw;
Jx_f(7,5) = dqsfpf_dLzps + dqbp_dLzps + dqsfss_dLzps;
Jx_f(7,6) = dqsfpf_dLzfs + dqsfss_dLzfs + dqbs_dLzfs;
Jx_f(7,7) = -dqfout1_dSf1;

Jx_f(8,7) = dqfout1_dSf1;
Jx_f(8,8) = -dqfout2_dSf2;

Jx_f(9,8) = dqfout2_dSf2;
Jx_f(9,9) = -dqfout3_dSf3;

Jx_f(10,9) = dqfout3_dSf3;

% ----------------------
% 6) Build Jth_f (m x d)
% ----------------------
Jth_f = zeros(m,d);

% T_tr/f_dd
Jth_f(1,14) = (P * dsnow_dT_tr) - dP_amelt_dT_tr;
Jth_f(1,15) = -dP_amelt_f_dd;

Jth_f(2,14) = dP_liq_dT_tr - dqsx_dT_tr - dqutof_dT_tr;
Jth_f(2,15) = dP_liq_f_dd - dqsx_f_dd - dqutof_f_dd;

Jth_f(3,14) = dqutof_dT_tr - dqufof_dT_tr;
Jth_f(3,15) = dqutof_f_dd - dqufof_f_dd;

Jth_f(7,14) = dqsx_dT_tr + dqufof_dT_tr;
Jth_f(7,15) = dqsx_f_dd + dqufof_f_dd;

% 1..13 blocks (same layout as your code)
% uzfwm (1)
j=1;
Jth_f(2,j) = 0.0;
Jth_f(3,j) = -dq12_duzfwm - dqif_duzfwm - dqufof_duzfwm;
Jth_f(4,j) =  pfree*dq12_duzfwm - dqstof_duzfwm;
Jth_f(5,j) =  dprc_duzfwm - dqsfpf_duzfwm;
Jth_f(6,j) =  dprc_duzfwm - dqsfss_duzfwm;
Jth_f(7,j) =  dqif_duzfwm + dqufof_duzfwm + dqsfpf_duzfwm + dqsfss_duzfwm;

% uztwm (2)
j=2;
Jth_f(2,j) = -dqsx_duztwm - de1_duztwm - dqutof_duztwm;
Jth_f(3,j) =  dqutof_duztwm - dqufof_duztwm;
Jth_f(4,j) = -de2_duztwm;
Jth_f(7,j) =  dqsx_duztwm + dqufof_duztwm;

% lzfpm (3)
j=3;
Jth_f(3,j) = -dq12_dlzfpm;
Jth_f(4,j) =  pfree*dq12_dlzfpm - dqstof_dlzfpm;
Jth_f(5,j) =  dprc_dlzfpm - dqsfpf_dlzfpm;
Jth_f(6,j) =  dprc_dlzfpm - dqsfss_dlzfpm;
Jth_f(7,j) =  dqsfpf_dlzfpm + dqsfss_dlzfpm;

% lzfsm (4)
j=4;
Jth_f(3,j) = -dq12_dlzfsm;
Jth_f(4,j) =  pfree*dq12_dlzfsm - dqstof_dlzfsm;
Jth_f(5,j) =  dprc_dlzfsm - dqsfpf_dlzfsm;
Jth_f(6,j) =  dprc_dlzfsm - dqsfss_dlzfsm;
Jth_f(7,j) =  dqsfpf_dlzfsm + dqsfss_dlzfsm;

% lztwm (5)
j=5;
Jth_f(3,j) = -dq12_dlztwm;
Jth_f(4,j) =  pfree*dq12_dlztwm - de2_dlztwm - dqstof_dlztwm;
Jth_f(5,j) =  dprc_dlztwm - dqsfpf_dlztwm;
Jth_f(6,j) =  dprc_dlztwm - dqsfss_dlztwm;
Jth_f(7,j) =  dqsfpf_dlztwm + dqsfss_dlztwm;

% zperc (6)
j=6;
Jth_f(3,j) = -dq12_dzperc;
Jth_f(4,j) =  pfree*dq12_dzperc - dqstof_dzperc;
Jth_f(5,j) =  dprc_dzperc - dqsfpf_dzperc;
Jth_f(6,j) =  dprc_dzperc - dqsfss_dzperc;
Jth_f(7,j) =  dqsfpf_dzperc + dqsfss_dzperc;

% rexp (7)
j=7;
Jth_f(3,j) = -dq12_drexp;
Jth_f(4,j) =  pfree*dq12_drexp - dqstof_drexp;
Jth_f(5,j) =  dprc_drexp - dqsfpf_drexp;
Jth_f(6,j) =  dprc_drexp - dqsfss_drexp;
Jth_f(7,j) =  dqsfpf_drexp + dqsfss_drexp;

% uzk (8)
j=8;
Jth_f(3,j) = -(B);     % dqif/duzk = B, appears with minus in Uzfw eq
Jth_f(7,j) =  (B);

% pfree (9)
j=9;
Jth_f(4,j) = q_12 - dqstof_dpfree;
Jth_f(5,j) = dprc_dpfree - dqsfpf_dpfree;
Jth_f(6,j) = dprc_dpfree - dqsfss_dpfree;
Jth_f(7,j) = dqsfpf_dpfree + dqsfss_dpfree;

% lzpk (10)
j=10;
q12_lzpk = lzfpm * dlz * B;
qstof_lzpk = pfree * q12_lzpk * lf_Lztw;
prc_lzpk = 0.5*(1.0-pfree)*q12_lzpk + 0.5*qstof_lzpk;
qsfp_lzpk = prc_lzpk * lf_Lzps;
qsfs_lzpk = prc_lzpk * lf_Lzfs;

Jth_f(3,j) = -q12_lzpk;
Jth_f(4,j) =  pfree*q12_lzpk - qstof_lzpk;
Jth_f(5,j) =  prc_lzpk - dqbp_dlzpk - qsfp_lzpk;
Jth_f(6,j) =  prc_lzpk - qsfs_lzpk;
Jth_f(7,j) =  qsfp_lzpk + qsfs_lzpk + dqbp_dlzpk;

% lzsk (11)
j=11;
q12_lzsk = lzfsm * dlz * B;
qstof_lzsk = pfree * q12_lzsk * lf_Lztw;
prc_lzsk = 0.5*(1.0-pfree)*q12_lzsk + 0.5*qstof_lzsk;
qsfp_lzsk = prc_lzsk * lf_Lzps;
qsfs_lzsk = prc_lzsk * lf_Lzfs;

Jth_f(3,j) = -q12_lzsk;
Jth_f(4,j) =  pfree*q12_lzsk - qstof_lzsk;
Jth_f(5,j) =  prc_lzsk - qsfp_lzsk;
Jth_f(6,j) =  prc_lzsk - dqbs_dlzsk - qsfs_lzsk;
Jth_f(7,j) =  qsfp_lzsk + qsfs_lzsk + dqbs_dlzsk;

% acm (12)
j=12;
Jth_f(2,j) = -dqsx_dacm - dqutof_dacm;
Jth_f(3,j) =  dqutof_dacm - dqufof_dacm;
Jth_f(7,j) =  dqsx_dacm + dqufof_dacm;

% kf (13)
j=13;
Jth_f(7,j) = -dqfout1_dkf;
Jth_f(8,j) =  dqfout1_dkf - dqfout2_dkf;
Jth_f(9,j) =  dqfout2_dkf - dqfout3_dkf;
Jth_f(10,j) =  dqfout3_dkf;

% ---------------------
% 7) Sensitivity update
% ---------------------
dSdt = Jx_f * Smat + Jth_f;

end

%% 4. SACSMA: RHS-only for J(x)_f and J(th)_f check
function [dxdt,dSdt,Jth_f,Jx_f] = rhs_only(t,z,th,data,T_sm,eps_m, ...
    eps_s,eps,rho,m,d)

% Unpack once
uzfwm = th(1);  % Maximum free water storage of upper zone (mm)
uztwm = th(2);  % Maximum tension water storage of upper zone (mm)
lzfpm = th(3);  % Maximum free water storage of lower zone primary (mm)
lzfsm = th(4);  % Maximum free water storage of lower zone secundary (mm)
lztwm = th(5);  % Maximum tension water storage of lower zone (mm)
zperc = th(6);  % Multiplier percolation function (-)
rexp = th(7);   % Power of percolation function (-)
uzk = th(8);    % Interflow rate (1/T)
pfree = th(9);  % Fraction percolation to tension storage lower layer (-)
lzpk = th(10);  % Base flow depletion rate for primary reservoir (1/T)
lzsk = th(11);  % Base flow depletion rate for secondary reservoir (1/T)
acm = th(12);   % Maximum fraction of saturated area (-)
kf = th(13);    % Recession constant fast reservoir (1/T)
                % --> to compute channel inflow
T_tr = th(14);  % Temperature threshold (°C)
f_dd = th(15);  % Degree-day factor (mm/°C/T)

Smat = ones(m,d);
% Call RHS
[dxdt,dSdt,Jth_f,Jx_f] = sacsma_odefcn(t,z,Smat,uzfwm,uztwm,lzfpm, ...
    lzfsm,lztwm,zperc,rexp,uzk,pfree,lzpk,lzsk,acm,kf,T_tr,f_dd, ...
    data,T_sm,eps_m,eps_s,eps,rho,m,d);

end

% QUESTIONS: 1. Do we have to account for evaporation of channel inflow?
%            2. Do we have to treat impervious fraction ( = state variable)
%               instead of saturation excess?
% NOTE:      1. Current formulation [no snow] has 13 parameters instead of
%               14 as impervious fraction is not modeled. This has a state
%               variable, evaporation, impervious area and direct runoff
%            2. q_if + q_sx + q_ufof + q_sfofa + q_sfofb + q_bp + q_bs - q_fout1
%               routed through three linear fast reservoirs with same
%               recession constant