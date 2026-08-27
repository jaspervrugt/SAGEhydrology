function varargout = hbv(par,mdl,data,ode,check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HBV: Runge Kutta implementation of Hbv conceptual watershed model
% SYNOPSIS: varargout = hbv(par,mdl,data,ode,check)
%  par          dx1 vector of parameter values
%   f_c:par(1)   field capacity (mm)
%   beta:par(2)  shape exponent (-)
%   lp:par(3)    evap limit fraction (-)
%   k_0:par(4)   near-surface recession (1/T)
%   uzl:par(5)   threshold for k0 (mm)
%   k_1:par(6)   upper zone recession (1/T)
%   k_2:par(7)   lower zone recession (1/T)
%   perc:par(8)  percolation max (mm/T)
%   T_tr:par(9)  temperature threshold (°C)
%   f_dd:par(10) melt factor (mm/°C/T)
%   sfcf:par(11) snowfall correction factor (-)
%   cfr:par(12)  refreezing factor (-)
%  mdl          structure with model state/parameter info
%   .mcode       scalar numerical solution of watershed model
%                 1 Runge Kutta implementation MATLAB
%                 2 ode45 implementation MATLAB
%                 3 Explicit Euler int_steps MATLAB
%                 4 Runge Kutta implementation ode_hbv C++
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
%   .mem        storage of state variables [0: no, 1: yes]
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
                            % 4: Runge Kutta implementation ode_hbv C++
    mem = ode.mem;          % state variable storage or not
    if mcode == 2 ...
            && mem == 0
        warning(['hbv: ' ...
            'built-in ode45 ' ...
            'solver stores ' ...
            'states: mem = 1'])
        mem = 1;
    end
    m = 5;                      % # state variables [+ snow: 1 state = SWE]
    ns = mdl.tout + 1;          % # print times
    d_par = numel(par);         % # parameters
    d = 12;                     % HBV ODE parameters; MAXBAS is external
    if ~ismember(d_par,[d,d+1])
        error('hbv:InvalidParameterCount', ...
            ['HBV requires 12 parameters without routing or ' ...
             '13 parameters with MAXBAS routing; received %d.'],d_par);
    end
    n = mdl.idx(2)-mdl.idx(1);  % # elements discharge vector
    nvar = m*(d+1);             % # number of variables
    eps_s = 1e-3;               % storage/flux smoothing (mm)
    eps_t = 0.1;                % temperature smoothing (°C)
    eps_x = 1e-12;              % smoothing coefficient (mm)
    rho = 0.01;                 % Dimensionless smoothing coefficient
    eps = 5;                    % Dimensionless smoothing coefficient
    fail = false;               % Default: model completes run
    id = m + (1:d)*m;           % Indices of sensitivity state variables
    if mem == 0
        q_n = nan(n,1);
        % Preallocate all requested parameter columns, including the optional
        % MAXBAS routing parameter. This avoids expanding the full Jacobian
        % after a long 15-minute simulation.
        J = nan(n,d_par);
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
            if (any(th<mdl.th_min) ...
                    || any(th>mdl.th_max))
                varargout = {nan(n,1),nan(n,d_par),nan(d_par,1),Z}; 
                return
            end
            Jth = ones(d_par,1);                 % return dq_n/dth
        case 1 % normalized hydrologic parameter values
            nth = par;
            if (any(nth<0) ...
                    || any(nth>1))
                varargout = {nan(n,1),nan(n,d_par),nan(d_par,1),Z}; return
            end
            dth_dnth = mdl.th_max - mdl.th_min; % dth/dnth
            th = mdl.th_min + nth.*dth_dnth;    % hydrologic parameter values
            Jth = dth_dnth;                     % return dq_n/dnth
        case 2 % unconstrained parameters (for training)
            varth = par;
            nth = 1./(1 + exp(-varth));         % normalized parameter values
            dth_dnth = mdl.th_max-mdl.th_min;   % dth/dnth
            th = mdl.th_min + nth.*dth_dnth;    % hydrologic parameter values
            dnth_dvarth = nth.*(1-nth);         % dnth/dvarth
            Jth = dth_dnth .* dnth_dvarth;      % return dq_n/dvarth
    end
    
    doRouting = (d_par == d+1);    % if caller provides extra parameter
    if doRouting
        B = th(d+1);            % MAXBAS is parameter 13 when d = 12
        [w_rt,dw_rtdb] = hbv_maxbas(B);
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
                t1 = s-2; t2 = s-1;             % Set start, end times
                h = min(hCarry,t2-t1);          % reuse prior recommendation
                if mem == 1
                    Z(s,1:nvar) = Z(s-1,1:nvar);
                    z = Z(s,1:nvar);                % row vector
                else
                    z = Z(1,1:nvar);
                end
                t = t1;                         % Initial time
                % Integrate from t1 to t2
                while (t < t2)
                    [ztmp,LTE] = rk2(t,z,h,th,data,m,... % Evaluate rk2
                        d,eps_s,eps_t,eps_x);
                    if any(~isfinite(ztmp)) || any(abs(ztmp) > 1e12)
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
                    h = min(hNext,t2-t);   % Another turn
                    if (iterCount >= maxiter)
                        if (flag == 0)
                            warning("WARNING: Max step limit " + ...
                                "reached at t = %.5f\n",t);
                            warning("WARNING: Parameter f_c    " + ...
                                "%8.5f\n",th(1));
                            warning("WARNING: Parameter beta   " + ...
                                "%8.5f\n",th(2));
                            warning("WARNING: Parameter lp     " + ...
                                "%8.5f\n",th(3));
                            warning("WARNING: Parameter k_0    " + ...
                                "%8.5f\n",th(4));
                            warning("WARNING: Parameter uzl    " + ...
                                "%8.5f\n",th(5));
                            warning("WARNING: Parameter k_1    " + ...
                                "%8.5f\n",th(6));
                            warning("WARNING: Parameter k_2    " + ...
                                "%8.5f\n",th(7));
                            warning("WARNING: Parameter perc   " + ...
                                "%8.5f\n",th(8));
                            warning("WARNING: Parameter T_tr   " + ...
                                "%8.5f\n",th(9));
                            warning("WARNING: Parameter f_dd   " + ...
                                "%8.5f\n",th(10));
                            warning("WARNING: Parameter sfcf   " + ...
                                "%8.5f\n",th(11));
                            warning("WARNING: Parameter cfr    " + ...
                                "%8.5f\n",th(12));
                            if doRouting
                                warning("WARNING: Parameter maxbas " + ...
                                    "%8.5f\n",th(13));
                            end
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
            [~,Z] = ode45(@(t,z) hbv_aug_ode(t,z,th,data, ...
                m,d,eps_s,eps_t,eps_x),0:mdl.tout,Z(1,1:nvar),ode_options);
            s = size(Z,1); if s < ns, fail = true; s = s+1; end
    
        case 3 %% MATLAB: Explicit Euler hbv_odefcn with int_steps steps
            int_steps = 500;                        % # integration steps
            dt = 1/int_steps;                       % integration timestep
            for s = 2:ns                            % Start time loop
                if mem == 1
                    z = Z(s-1,1:nvar);              % Initialize state variables
                else
                    z = Z(1,1:nvar);
                end
                for it = 1:int_steps                % integration in int_steps steps
                    dzdt = hbv_aug_ode(s-2,z, ...       % compute dzdt based on 
                        th,data,m,d,eps_s,eps_t,eps_x)';% current state, x,P,Ep
                    if s == 150 ...
                            && check
                        [~,~,Jth_f,Jx_f] = rhs_only(s-2,z,th,data,m,d, ...
                            eps_s,eps_t,eps_x);
                        % --- Jth_f via 2-sided finite differences ---
                        Jth_fn = zeros(m,d);
                        for j = 1:d
                            h = 1e-6*max(1,abs(th(j)));
                            thp = th; thm = th;
                            thp(j) = thp(j) + h; thm(j) = thm(j) - h;
                            fp = rhs_only(s-2,z,thp,data,m,d, ...
                                eps_s,eps_t,eps_x);
                            fm = rhs_only(s-2,z,thm,data,m,d, ...
                                eps_s,eps_t,eps_x);
                            Jth_fn(:,j) = (fp - fm)/(2*h);
                        end
                        % --- Jx_f via 2-sided finite differences ---
                        Jx_fn = zeros(m,m);
                        for i = 1:m
                            h = 1e-6*max(1,abs(z(i)));
                            zp = z; zm = z;
                            zp(i) = zp(i) + h; zm(i) = zm(i) - h;
                            fp = rhs_only(s-2,zp,th,data,m,d, ...
                                eps_s,eps_t,eps_x);
                            fm = rhs_only(s-2,zm,th,data,m,d, ...
                                eps_s,eps_t,eps_x);
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
            data.f_c = th(1);       % field capacity (mm)
            data.beta = th(2);      % shape exponent (-)
            data.lp = th(3);        % evap limit fraction (-)
            data.k_0 = th(4);       % near-surface recession (1/T)
            data.uzl = th(5);       % threshold for k0 (mm)
            data.k_1 = th(6);       % upper zone recession (1/T)
            data.k_2 = th(7);       % lower zone recession (1/T)
            data.perc = th(8);      % percolation max (mm/T)
            data.T_tr = th(9);      % temperature threshold (°C)
            data.f_dd = th(10);     % melt factor (mm/°C/T)
            data.sfcf = th(11);     % snowfall correction factor (-)
            data.cfr = th(12);      % refreezing factor (-)
    
            data.eps = eps;         % Dimensionless smoothing coef. [inactive]
            data.rho = rho;         % Dimensionless smoothing coef. [inactive]
            data.eps_s = eps_s;     % storage/flux smoothing (mm)
            data.eps_t = eps_t;     % temperature smoothing (°C)
            data.eps_x = eps_x;     % smoothing coefficient (mm)
            data.ipr = ipr;         % Time to print

            %[Z,q_n,J] = crr_hbv(mdl.tout,Z(1,1:nvar)',data,ode);
            if nargout == 1
                if mem == 1
                    [Z,~] = crr_hbv( ...
                        mdl.tout,Z(1,1:nvar)',data,ode);
                    q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
                else
                    [~,q_n] = crr_hbv( ...
                        mdl.tout,Z(1,1:nvar)',data,ode);
                end
                J = [];
            else
                [Z,q_n,J] = crr_hbv( ...
                    mdl.tout,Z(1,1:nvar)',data,ode);
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
    
                % Add the optional routing-parameter column before filtering.
                % Preallocation here prevents MATLAB from copying and growing
                % the complete long Jacobian at J(:,d+1).
                if doRouting
                    J(:,d+1) = 0;
                end
        end
    end
    
    if doRouting
        if nargout > 1
            % Compute dq/dB while q_n still contains unrouted discharge.
            % This removes the additional full-length q_norout copy.
            J(:,d+1) = filter(dw_rtdb,1,q_n);
        end    
        % Route discharge only after its derivative with respect to B has
        % been computed from the unrouted series.
        q_n = filter(w_rt,1,q_n);
    end
    if nargout == 1
        varargout = {q_n}; return
    else
        if doRouting
            % Route discharge derivatives wrt the ODE parameters.
            for j = 1:d
                J(:,j) = filter(w_rt,1,J(:,j));
            end
        end
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
function [z,LTE] = rk2(t,z,h,th,data,m,d,eps_s,eps_t,eps_x)

    dzdtE = hbv_aug_ode(t,z,th,data,m,d,eps_s,eps_t,eps_x)';      % Euler
    zE = z + h*dzdtE;
    dzdtH = hbv_aug_ode(t,zE,th,data,m,d,eps_s,eps_t,eps_x)';     % Heun
    z = z + 0.5*h*(dzdtE + dzdtH);                                % New z
    LTE = abs(zE - z);                                            % LTE

end

%% 2. HBV augmented ode with sensitivities as state variables
function dzdt = hbv_aug_ode(t,z,th,data,m,d,eps_s,eps_t,eps_x)

    x = z(1:m);                                     % x = [Swe Sm Uz Lz Q]
    S = reshape(z(m+1:m*(d+1)),m,d);                % mxd sensitivity matrix
    [dxdt,dSdt] = hbv_odefcn(t,x,th,S,data,m,d, ...
        eps_s,eps_t,eps_x);                         % compute dxdt & dSdt
    dzdt = [dxdt; dSdt(:)];                         % Repack single vector

end

function [dxdt,dSdt,Jth_f,Jx_f] = hbv_odefcn(t,x,th,Smat,data,m,d, ...
    eps_s,eps_t,eps_x)
% HBV ODE core + forward sensitivities (fully smooth incl. "safety" maxes)
% States x = [Swe Sm Uz Lz Qcum]  (m = 5)
% Parameters th(1:12) are used by the ODE. Parameter 13 (MAXBAS)
% is applied outside the ODE as differentiable routing.

    dxdt = nan(m,1);
    
    % -------------------------
    % Unpack parameters
    % -------------------------
    f_c = th(1);
    beta = th(2);
    lp = th(3);
    k_0 = th(4);
    uzl = th(5);
    k_1 = th(6);
    k_2 = th(7);
    perc = th(8);
    T_tr = th(9);
    f_dd = th(10);
    sfcf = th(11);
    cfr = th(12);
    
    % -------------------------
    % Unpack states (RAW)
    % -------------------------
    Swe_raw = x(1);
    Sm_raw = x(2);
    Uz_raw = x(3);
    Lz_raw = x(4);
    
    % -----------------------------
    % Forcing at current time index
    % -----------------------------
    it = floor(t) + 1;
    P = data.P(it,1);
    Ep = data.Ep(it,1);
    T = data.T(it,1);
    
    % -----------------
    % Smooth primitives
    % -----------------
    smooth_pos = @(a,ep) 0.5*(a + sqrt(a.^2 + ep.^2));
    dsmooth_pos_da = @(a,ep) 0.5*(1 + a./sqrt(a.^2 + ep.^2));
    
    smooth_min = @(A,B,ep) 0.5*(A + B - sqrt((A-B).^2 + ep.^2));
    dmin_dA = @(A,B,ep) 0.5*(1 - (A-B)./sqrt((A-B).^2 + ep.^2)); % d/dA
    dmin_dB = @(A,B,ep) 0.5*(1 + (A-B)./sqrt((A-B).^2 + ep.^2)); % d/dB
    
    smooth_max = @(A,B,ep) 0.5*(A + B + sqrt((A-B).^2 + ep.^2));
    dmax_dA = @(A,B,ep) 0.5*(1 + (A-B)./sqrt((A-B).^2 + ep.^2)); % d/dA
    
    % clamp01(z) = min(max(z,0),1) with consistent derivative
    clamp01 = @(z,ep) smooth_min( smooth_max(z,0,ep), 1, ep );
    
    dclamp01_dz = @(z,ep) ...
        dmin_dA( smooth_max(z,0,ep), 1, ep ) .* ...
        dmax_dA( z, 0, ep );
    
    % evap limiter: phiE = smooth_min(ratio,1)
    phiE_fun = @(ratio,ep) smooth_min(ratio,1,ep);
    dphiE_dratio = @(ratio,ep) dmin_dA(ratio,1,ep);
    
    % ------------------------------------
    % Smooth-clamp storages used in fluxes
    % ------------------------------------
    Swe = smooth_pos(Swe_raw,eps_s);
    Sm = smooth_pos(Sm_raw,eps_s);
    Uz = smooth_pos(Uz_raw,eps_s);
    Lz = smooth_pos(Lz_raw,eps_s);
    
    dSwe_dSweRaw = dsmooth_pos_da(Swe_raw,eps_s);
    dSm_dSmRaw = dsmooth_pos_da(Sm_raw,eps_s);
    dUz_dUzRaw = dsmooth_pos_da(Uz_raw,eps_s);
    dLz_dLzRaw = dsmooth_pos_da(Lz_raw,eps_s);
    
    % ------------------------------------------
    % 1) Precip partition (smooth snow fraction)
    % ------------------------------------------
    uT = (T - T_tr) / eps_t;
    snow_fr = 0.5*(1 - tanh(uT));
    rain_fr = 1 - snow_fr;
    
    sech2 = 1/(cosh(uT)^2);
    dsnow_dT_tr = 0.5 * sech2 / eps_t;  % as in your code
    
    Ps = sfcf * P * snow_fr;
    Pr = P * rain_fr;
    
    dPs_dsfcf = P * snow_fr;
    dPs_dT_tr = sfcf * P * dsnow_dT_tr;
    dPr_dT_tr = -P * dsnow_dT_tr;
    
    % ------------------
    % 2) Melt & refreeze
    % ------------------
    Tm = T - T_tr;
    
    posTm = smooth_pos(Tm,  eps_t);
    negTm = smooth_pos(-Tm, eps_t);
    
    dposTm_dTm = dsmooth_pos_da(Tm,  eps_t);
    dnegTm_dTm = -dsmooth_pos_da(-Tm, eps_t);
    
    Mpot = f_dd * posTm;
    Rpot = cfr * f_dd * negTm;
    
    dMpot_dT_tr = -f_dd * dposTm_dTm;
    dMpot_df_dd = posTm;
    
    dRpot_dT_tr = -(cfr*f_dd) * dnegTm_dTm;
    dRpot_df_dd = cfr * negTm;
    dRpot_dcfr = f_dd * negTm;
    
    M = smooth_min(Swe, Mpot, eps_s);
    R = smooth_min(Pr,  Rpot, eps_s);
    
    dM_dSwe = dmin_dA(Swe, Mpot, eps_s);
    dM_dMpot = dmin_dB(Swe, Mpot, eps_s);
    
    dR_dPr = dmin_dA(Pr, Rpot, eps_s);
    dR_dRpot = dmin_dB(Pr, Rpot, eps_s);
    
    dM_dT_tr = dM_dMpot * dMpot_dT_tr;
    dM_df_dd = dM_dMpot * dMpot_df_dd;
    
    dR_dT_tr = dR_dPr * dPr_dT_tr + dR_dRpot * dRpot_dT_tr;
    dR_df_dd = dR_dRpot * dRpot_df_dd;
    dR_dcfr = dR_dRpot * dRpot_dcfr;
    
    % ------------------------------------------------------
    % 3) Liquid reaching soil: Pliq = smooth_pos(Pr + M - R)
    % ------------------------------------------------------
    Win_raw = Pr + M - R;
    Pliq = smooth_pos(Win_raw, eps_s);
    dPliq_dWinraw = dsmooth_pos_da(Win_raw, eps_s);
    
    % wrt Swe (via M)
    dPliq_dSwe = dPliq_dWinraw * dM_dSwe;
    
    % wrt T_tr/f_dd/cfr
    dPliq_dT_tr = dPliq_dWinraw * (dPr_dT_tr + dM_dT_tr - dR_dT_tr);
    dPliq_df_dd = dPliq_dWinraw * (dM_df_dd - dR_df_dd);
    dPliq_dcfr = dPliq_dWinraw * (0 - dR_dcfr);
    
    % ------------------------------------------------------
    % 4) Soil evap + recharge (FIXED: den and up are smooth)
    % ------------------------------------------------------
    
    % --- Smooth "den = max(lp*f_c,eps_x)" with derivatives
    den0 = lp * f_c;
    den = smooth_max(den0,eps_x,eps_x);          % smooth max
    dDen_dDen0 = dmax_dA(den0,eps_x,eps_x);      % d den / d den0
    
    % ratio = Sm/den
    ratio = Sm / den;
    
    % derivatives of ratio (now consistent)
    dratio_dSm = 1/den;
    dratio_dDen = -Sm/(den^2);
    dratio_dDen0 = dratio_dDen * dDen_dDen0;
    
    dratio_dlp = dratio_dDen0 * (f_c);
    dratio_dfc = dratio_dDen0 * (lp);
    
    phiE = phiE_fun(ratio, eps_s);
    dphiE_drat = dphiE_dratio(ratio, eps_s);
    
    E = Ep * phiE;
    
    dE_dSm = Ep * dphiE_drat * dratio_dSm;
    dE_dlp = Ep * dphiE_drat * dratio_dlp;
    dE_dfc = Ep * dphiE_drat * dratio_dfc;
    
    % --- Recharge: Re = Pliq * clamp01(Sm/fc_safe)^beta
    % Smooth "fc_safe = max(f_c, eps_x)" with derivatives
    fc_safe = smooth_max(f_c, eps_x, eps_x);
    dfcSafe_dfc = dmax_dA(f_c, eps_x, eps_x);
    
    u_raw = Sm / fc_safe;
    
    % derivatives of u_raw
    duraw_dSm = 1/fc_safe;
    duraw_dfc = -Sm/(fc_safe^2) * dfcSafe_dfc;
    
    u = clamp01(u_raw, eps_s);
    du_duraw = dclamp01_dz(u_raw, eps_s);
    
    du_dSm = du_duraw * duraw_dSm;
    du_dfc = du_duraw * duraw_dfc;
    
    % Smooth lower bound used by both recharge and its derivatives.
    up = smooth_max(u,eps_x,eps_x);
    dup_du = dmax_dA(u,eps_x,eps_x);
    up_beta = up^beta;

    Re = Pliq * up_beta;
    dRe_dPliq = up_beta;
    dRe_dSm = Pliq * beta * up^(beta-1) * dup_du * du_dSm;
    dRe_dfc = Pliq * beta * up^(beta-1) * dup_du * du_dfc;
    dRe_dbeta = Pliq * up_beta * log(up);
    
    % ---------------------------
    % 5) Response routine (Uz/Lz)
    % ---------------------------
    h = smooth_pos(Uz - uzl, eps_s);
    dh_dUz = dsmooth_pos_da(Uz - uzl, eps_s);
    dh_duzl = -dh_dUz;
    
    q0 = k_0 * h;
    dq0_dUz = k_0 * dh_dUz;
    dq0_dk0 = h;
    dq0_duzl = k_0 * dh_duzl;
    
    q1 = k_1 * Uz;
    dq1_dUz = k_1;
    dq1_dk1 = Uz;
    
    perc_flux = smooth_min(perc, Uz, eps_s);
    dperc_dperc = dmin_dA(perc, Uz, eps_s);
    dperc_dUz = dmin_dB(perc, Uz, eps_s);
    
    q2 = k_2 * Lz;
    dq2_dLz = k_2;
    dq2_dk2 = Lz;
    
    Q = q0 + q1 + q2;
    
    % ----------------------------
    % State ODEs (RAW state space)
    % ----------------------------
    dxdt(1) = Ps + R - M;
    dxdt(2) = Pliq - Re - E;
    dxdt(3) = Re - q0 - q1 - perc_flux;
    dxdt(4) = perc_flux - q2;
    dxdt(5) = Q;
    
    % ---------------------------------------------------------------
    % Jacobian wrt states: build in SMOOTHED-space, then chain to RAW
    % ---------------------------------------------------------------
    Jx_f = zeros(m,m);
    
    % f1
    Jx_f(1,1) = -dM_dSwe;
    
    % f2
    Jx_f(2,1) = dPliq_dSwe - dRe_dPliq*dPliq_dSwe;
    Jx_f(2,2) = -dRe_dSm - dE_dSm;
    
    % f3
    Jx_f(3,1) = dRe_dPliq*dPliq_dSwe;
    Jx_f(3,2) = dRe_dSm;
    Jx_f(3,3) = -(dq0_dUz + dq1_dUz + dperc_dUz);
    
    % f4
    Jx_f(4,3) = dperc_dUz;
    Jx_f(4,4) = -dq2_dLz;
    
    % f5
    Jx_f(5,3) = dq0_dUz + dq1_dUz;
    Jx_f(5,4) = dq2_dLz;
    
    % Chain columns back to RAW states
    Jx_f(:,1) = Jx_f(:,1) * dSwe_dSweRaw;
    Jx_f(:,2) = Jx_f(:,2) * dSm_dSmRaw;
    Jx_f(:,3) = Jx_f(:,3) * dUz_dUzRaw;
    Jx_f(:,4) = Jx_f(:,4) * dLz_dLzRaw;
    
    % ------------------
    % Parameter Jacobian
    % ------------------
    Jth_f = zeros(m,d);
    
    % indices
    j_fc=1; j_beta=2; j_lp=3; j_k0=4; j_uzl=5; j_k1=6; j_k2=7; j_perc=8;
    j_T_tr=9; j_f_dd=10; j_sfcf=11; j_cfr=12;
    
    % f1: Ps + R - M
    Jth_f(1,j_T_tr) = dPs_dT_tr + dR_dT_tr - dM_dT_tr;
    Jth_f(1,j_f_dd) = dR_df_dd - dM_df_dd;
    Jth_f(1,j_sfcf) = dPs_dsfcf;
    Jth_f(1,j_cfr) = dR_dcfr;
    
    % f2: Pliq - Re - E
    Jth_f(2,j_T_tr) = dPliq_dT_tr    - dRe_dPliq*dPliq_dT_tr;
    Jth_f(2,j_f_dd) = dPliq_df_dd - dRe_dPliq*dPliq_df_dd;
    Jth_f(2,j_cfr) = dPliq_dcfr   - dRe_dPliq*dPliq_dcfr;
    
    Jth_f(2,j_fc) = -dRe_dfc - dE_dfc;
    Jth_f(2,j_beta) = -dRe_dbeta;
    Jth_f(2,j_lp) = -dE_dlp;
    
    % f3: Re - q0 - q1 - perc_flux
    Jth_f(3,j_T_tr) = dRe_dPliq*dPliq_dT_tr;
    Jth_f(3,j_f_dd) = dRe_dPliq*dPliq_df_dd;
    Jth_f(3,j_cfr) = dRe_dPliq*dPliq_dcfr;
    
    Jth_f(3,j_fc) = dRe_dfc;
    Jth_f(3,j_beta) = dRe_dbeta;
    Jth_f(3,j_k0) = -dq0_dk0;
    Jth_f(3,j_uzl) = -dq0_duzl;
    Jth_f(3,j_k1) = -dq1_dk1;
    Jth_f(3,j_perc) = -dperc_dperc;
    
    % f4: perc_flux - q2
    Jth_f(4,j_k2) = -dq2_dk2;
    Jth_f(4,j_perc) =  dperc_dperc;
    
    % f5: Q
    Jth_f(5,j_k0) = dq0_dk0;
    Jth_f(5,j_uzl) = dq0_duzl;
    Jth_f(5,j_k1) = dq1_dk1;
    Jth_f(5,j_k2) = dq2_dk2;
    
    % ---------------
    % Sensitivity ODE
    % ---------------
    dSdt = Jx_f * Smat + Jth_f;

end

% ----------------
% helper functions
% ----------------
function [w,dwdb] = hbv_maxbas(b,L_max,eps_w)
%HBV_MAXBAS_RELAXED  Differentiable (relaxed) maxbas routing.
%
% Qout(t) = sum_{k=0}^{L_max-1} w_k(b) * q_raw(t-k)
% b is continuous analogue of maxbas (in time steps).
% SYNOPSIS: [w,dwdb] = hbv_maxbas(b,L_max,eps_w)
%   b       scalar > 0, continuous width parameter (time steps)
%   L_max   integer >= 1, fixed kernel length [def: 60]
%   eps_w   small smoothing constant [def: 1e-8]
%   w       OUTPUT: L_max x 1 weights, sum(w) = 1
%   dwdb    OUTPUT: L_max x 1 derivative dw/db

    if nargin < 3 || isempty(eps_w), eps_w = 1e-8; end
    if nargin < 2 || isempty(L_max), L_max = 60; end
    
    % ---- smooth primitives ----
    smooth_pos = @(a,ep) 0.5*(a + sqrt(a.^2 + ep.^2));
    dsmooth_pos_da = @(a,ep) 0.5*(1 + a./sqrt(a.^2 + ep.^2));
    smooth_max = @(A,B,ep) 0.5*(A + B + sqrt((A-B).^2 + ep.^2));
    dmax_dA = @(A,B,ep) 0.5*(1 + (A-B)./sqrt((A-B).^2 + ep.^2)); % d/dA
    smooth_abs = @(z,ep) sqrt(z.^2 + ep.^2);
    dabs_dz = @(z,ep) z ./ sqrt(z.^2 + ep.^2);
    
    % ---- continuous triangle on k = 0..Lmax-1 ----
    k = (0:L_max-1)';               % causal lags
    b = max(b,eps_w);
    c = 0.5*b;                      % center
    half = 0.5*b;                   % half-width
    denom = smooth_max(half, ...
        eps_w,eps_w);
    ddenom_db = 0.5 * dmax_dA( ...  % d(denom)/d(half) * d(half)/db
        half,eps_w,eps_w);
    z = k - c;                      % distance from center
    az = smooth_abs(z,eps_w);
    tri = 1 - az./denom;            % unthresholded triangle
    % smooth cutoff at 0 (instead of hard max(0,tri))
    a = smooth_pos(tri,eps_w);
    dz_db = -0.5;                   % derivatives da/db
    daz_db = dabs_dz(z,eps_w) .* dz_db;
    dtri_db = -(daz_db)./denom + az./(denom.^2) .* ddenom_db;
    da_db = dsmooth_pos_da(tri,eps_w) .* dtri_db;
    
    sa = sum(a);                    % normalize: w = a / sum(a)
    sda = sum(da_db);
    w = a / sa;
    dwdb = (da_db*sa - a*sda) / (sa^2);

end

%% 4. HBV: RHS-only for J(x)_f and J(th)_f check
function [dxdt,dSdt,Jth_f,Jx_f] = rhs_only(t,z,th,data,m,d, ...
    eps_s,eps_t,eps_x)
    
    Smat = ones(m,d);
    % Call RHS
    [dxdt,dSdt,Jth_f,Jx_f] = hbv_odefcn(t,z,th,Smat,data,m,d, ...
        eps_s,eps_t,eps_x);
end

% function varargout = hbv(par,mdl,data,ode,check)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %HBV: Runge Kutta implementation of Hbv conceptual watershed model
% % SYNOPSIS: varargout = hbv(par,mdl,data,ode,check)
% %  par          dx1 vector of parameter values
% %   f_c:par(1)   field capacity (mm)
% %   beta:par(2)  shape exponent (-)
% %   lp:par(3)    evap limit fraction (-)
% %   k_0:par(4)   near-surface recession (1/T)
% %   uzl:par(5)   threshold for k0 (mm)
% %   k_1:par(6)   upper zone recession (1/T)
% %   k_2:par(7)   lower zone recession (1/T)
% %   perc:par(8)  percolation max (mm/T)
% %   T_tr:par(9)  temperature threshold (°C)
% %   f_dd:par(10) melt factor (mm/°C/T)
% %   sfcf:par(11) snowfall correction factor (-)
% %   cfr:par(12)  refreezing factor (-)
% %  mdl          structure with model state/parameter info
% %   .mcode       scalar numerical solution of watershed model
% %                 1 Runge Kutta implementation MATLAB
% %                 2 ode45 implementation MATLAB
% %                 3 Explicit Euler int_steps MATLAB
% %                 4 Runge Kutta implementation ode_hbv C++
% %   .y0          mx1 vector of initial states
% %   .pspace      0:hydrologic, 1:unit cube, 2:unconstrained parameters
% %   .th_min      dx1 vector of lower parameter values [= in pspace]
% %   .th_max      dx1 vector of upper parameter values [= in pspace]
% %   .par_names   1xd cell parameter names
% %   .id_train    1x2 vector of start and end index training period
% %   .id_eval     1x2 vector of start and end index evaluation period
% %   .eval_mode   evaluation design used during SAGE training
% %     'per'       training basins evaluated on evaluation period
% %     'bas'       evaluation basins evaluated on training period
% %     'basper'    evaluation basins evaluated on evaluation period
% %     'none'      no evaluation
% %   .tout        final model print time (scalar)
% %   .idx         1x2 vector of indices train&val periods
% %  data         structure with meteorological data and other info
% %   .P           (n+m)x1 record of precipitation (mm/T)
% %   .Ep          (n+m)x1 record of potential evapotranspiration (mm/T)
% %   .T           (n+m)x1 record of air temperature (°C)
% %  ode          structure with numerical settings ODE solver
% %   .InitStep    Initial time step
% %   .MaxStep     Maximum time step
% %   .MinStep     Minimum time step
% %   .RelTol      Relative tolerance
% %   .AbsTol      Absolute tolerance
% %   .Order       Order
% %   .maxiter     Maximum number of iterations
% %   .mem        storage of state variables [0: no, 1: yes]
% %  check        numerical check of J(x)_f and J(x)_th matrices (or not)
% %   0            do not check
% %   1            check Jacobian matrices of states and parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if nargin < 5
%     check = 0;          % no check of J(x)_f and J(x)_th
% end
% mcode = mdl.mcode;      % Formulation/language
%                         % 1: Runge Kutta implementation MATLAB
%                         % 2: ode45 implementation MATLAB
%                         % 3: Explicit Euler int_steps MATLAB
%                         % 4: Runge Kutta implementation ode_hbv C++
% mem = ode.mem;          % state variable storage or not
% if mcode == 2 && mem == 0
%     warning(['hbv: ' ...
%         'built-in ode45 ' ...
%         'solver stores ' ...
%         'states: mem = 1'])
%     mem = 1;
% end
% m = 5;                      % # state variables [+ snow: 1 state = SWE]
% ns = mdl.tout + 1;          % # print times
% d_par = numel(par);         % # parameters
% d = 12;                     % HBV only: maxbas [routing] outside odefnc
% n = mdl.idx(2)-mdl.idx(1);  % # elements discharge vector
% nvar = m*(d+1);             % # number of variables
% eps_s = 1e-3;               % storage/flux smoothing (mm)
% eps_t = 0.1;                % temperature smoothing (°C)
% eps_x = 1e-12;              % smoothing coefficient (mm)
% rho = 0.01;                 % Dimensionless smoothing coefficient
% eps = 5;                    % Dimensionless smoothing coefficient
% fail = false;               % Default: model completes run
% id = m + (1:d)*m;           % Indices of sensitivity state variables
% if mem == 0
%     q_n = nan(n,1);
%     J = nan(n,d); 
%     ipr = mdl.idx(1);       % --> C++ code
% else
%     Z = nan(ns,nvar); 
%     ipr = 0;
% end
% Z(1,1:m) = mdl.y0;          % Initialize state variables at time 0
% Z(1,m+1:nvar) = 0;          % Initialize sensitivity at time 0
% 
% switch mdl.pspace
%     case 0 % hydrologic parameter values
%         th = par;
%         if (any(th<mdl.th_min) || any(th>mdl.th_max))
%             varargout = {nan(n,1),nan(n,d_par),nan(d_par,1),Z}; return
%         end
%         Jth = ones(d_par,1);                 % return dq_n/dth
%     case 1 % normalized hydrologic parameter values
%         nth = par;
%         if (any(nth<0) || any(nth>1))
%             varargout = {nan(n,1),nan(n,d_par),nan(d_par,1),Z}; return
%         end
%         dth_dnth = mdl.th_max - mdl.th_min; % dth/dnth
%         th = mdl.th_min + nth.*dth_dnth;    % hydrologic parameter values
%         Jth = dth_dnth;                     % return dq_n/dnth
%     case 2 % unconstrained parameters (for training)
%         varth = par;
%         nth = 1./(1 + exp(-varth));         % normalized parameter values
%         dth_dnth = mdl.th_max-mdl.th_min;   % dth/dnth
%         th = mdl.th_min + nth.*dth_dnth;    % hydrologic parameter values
%         dnth_dvarth = nth.*(1-nth);         % dnth/dvarth
%         Jth = dth_dnth .* dnth_dvarth;      % return dq_n/dvarth
% end
% 
% doRouting = (d_par > d);    % if caller provides extra parameter
% if doRouting
%     B = th(d+1);            % expects B = par(15) when d=14
%     [w_rt,dw_rtdb] = hbv_maxbas(B);
% end
% 
% switch mcode
% 
%     case 1 %% MATLAB: Runge Kutta implementation
%         hin = ode.InitStep;     % Initial time step
%         hmax_ = ode.MaxStep;    % Maximum time step
%         hmin_ = ode.MinStep;    % minimum time step
%         reltol = ode.RelTol;    % Relative tolerance
%         abstol = ode.AbsTol;    % Absolute tolerance
%         order = ode.Order;      % Order
%         maxiter = ode.maxiter;  % Maximum iterations
%         iterCount = 0; flag = 0;
% 
%         for s = 2:ns            % Loop over elements
%             t1 = s-2; t2 = s-1;                     % Set start, end times
%             h = hin;                                % Set initial step
%             h = max(hmin_,min(h,hmax_));            % Make sure in range
%             h = min(h,t2-t1);                       % h based on t2-t1
%             if mem == 1
%                 Z(s,1:nvar) = Z(s-1,1:nvar);
%                 z = Z(s,1:nvar);                    % row vector
%             else
%                 z = Z(1,1:nvar);
%             end
%             t = t1;                                 % Initial time
%             % Integrate from t1 to t2
%             while (t < t2)
%                 [ztmp,LTE] = rk2(t,z,h,th,data,m,... % Evaluate rk2
%                     d,eps_s,eps_t,eps_x);
%                 if any(~isfinite(ztmp)) || any(abs(ztmp) > 1e12)
%                     fail = true; break;
%                 end
%                 w = 1 ./ (reltol*abs(ztmp) + abstol);   % Weights
%                 wrms = sqrt(sum((w.*LTE).^2)/nvar);     % CORRECTED, Nov. 2022
%                 if (wrms <= 1) || (h <= hmin_)          % Accept if error is small enough
%                     z = ztmp; t = t + h;
%                 end
%                 h = h*max(0.2,min(5.0,0.9*wrms^(-1/order)));    % Compute new step
%                 h = max(hmin_,min(h,hmax_)); h = min(h,t2-t);   % Another turn
%                 if (iterCount >= maxiter)
%                     if (flag == 0)
%                         warning("WARNING: Max step limit " + ...
%                             "reached at t = %.5f\n",t);
%                         warning("WARNING: Parameter f_c    %8.5f\n",th(1));
%                         warning("WARNING: Parameter beta   %8.5f\n",th(2));
%                         warning("WARNING: Parameter lp     %8.5f\n",th(3));
%                         warning("WARNING: Parameter k_0    %8.5f\n",th(4));
%                         warning("WARNING: Parameter uzl    %8.5f\n",th(5));
%                         warning("WARNING: Parameter k_1    %8.5f\n",th(6));
%                         warning("WARNING: Parameter k_2    %8.5f\n",th(7));
%                         warning("WARNING: Parameter perc   %8.5f\n",th(8));
%                         warning("WARNING: Parameter T_tr   %8.5f\n",th(9));
%                         warning("WARNING: Parameter f_dd   %8.5f\n",th(10));
%                         warning("WARNING: Parameter sfcf   %8.5f\n",th(11));
%                         warning("WARNING: Parameter cfr    %8.5f\n",th(12));
%                         warning("WARNING: Parameter maxbas %8.5f\n",th(13));
%                         %  warning("WARNING: Parameter cwh    %8.5f\n",th(14));
%                         flag = 1; fail = true; break;
%                     end
%                 else
%                     iterCount = iterCount + 1;
%                 end
%                 if flag == 1
%                     fprintf('Time: %7.5f  h: %7.5f  iterCount: %d', ...
%                         t,h,iterCount);
%                 end                
%             end
%             if fail
%                 break;
%             else
%                 iterCount = 0;
%             end
%             % -----------------------------
%             if mem == 1
%                 Z(s,1:nvar) = z;
%             else
%                 if s >= ipr+1
%                     q_n(s-ipr) = z(m) - Z(m);
%                     J(s-ipr,1:d) = z(id) - Z(id);
%                 end
%                 Z(1,1:nvar) = z;
%             end
%         end
% 
%     case 2 %% MATLAB: ode45 implementation + augmented sensitivities
%         hin = ode.InitStep;     % Initial time step
%         hmax_ = ode.MaxStep;    % Maximum time step
%         reltol = ode.RelTol;    % Relative tolerance
%         abstol = ode.AbsTol;    % Absolute tolerance
%         mdl.tout = mdl.tout-1e-10;    % Loop ode goes to maxT and then still one try
%         % Then, time index goes out of bound of forcing data
%         ode_options = odeset('InitialStep',hin,...  % initial time-step (T)
%             'MaxStep',hmax_, ...                    % maximum time-step (T)
%             'RelTol',reltol, ...                    % relative tolerance
%             'AbsTol',abstol);                       % absolute tol (mm)
%         [~,Z] = ode45(@(t,z) hbv_aug_ode(t,z,th,data, ...
%             m,d,eps_s,eps_t,eps_x),0:mdl.tout,Z(1,1:nvar),ode_options);
%         s = size(Z,1); if s < ns, fail = true; s = s+1; end
% 
%     case 3 %% MATLAB: Explicit Euler hbv_odefcn with int_steps steps
%         int_steps = 500;                        % # integration steps
%         dt = 1/int_steps;                       % integration timestep
%         for s = 2:ns                            % Start time loop
%             if mem == 1
%                 z = Z(s-1,1:nvar);              % Initialize state variables
%             else
%                 z = Z(1,1:nvar);
%             end
%             for it = 1:int_steps                % integration in int_steps steps
%                 dzdt = hbv_aug_ode(s-2,z, ...       % compute dzdt based on 
%                     th,data,m,d,eps_s,eps_t,eps_x)';% current state, x,P,Ep
%                 if s == 150 && check
%                     [~,~,Jth_f,Jx_f] = rhs_only(s-2,z,th,data,m,d, ...
%                         eps_s,eps_t,eps_x);
%                     % --- Jth_f via 2-sided finite differences ---
%                     Jth_fn = zeros(m,d);
%                     for j = 1:d
%                         h = 1e-6*max(1,abs(th(j)));
%                         thp = th; thm = th;
%                         thp(j) = thp(j) + h; thm(j) = thm(j) - h;
%                         fp = rhs_only(s-2,z,thp,data,m,d, ...
%                             eps_s,eps_t,eps_x);
%                         fm = rhs_only(s-2,z,thm,data,m,d, ...
%                             eps_s,eps_t,eps_x);
%                         Jth_fn(:,j) = (fp - fm)/(2*h);
%                     end
%                     % --- Jx_f via 2-sided finite differences ---
%                     Jx_fn = zeros(m,m);
%                     for i = 1:m
%                         h = 1e-6*max(1,abs(z(i)));
%                         zp = z; zm = z;
%                         zp(i) = zp(i) + h; zm(i) = zm(i) - h;
%                         fp = rhs_only(s-2,zp,th,data,m,d, ...
%                             eps_s,eps_t,eps_x);
%                         fm = rhs_only(s-2,zm,th,data,m,d, ...
%                             eps_s,eps_t,eps_x);
%                         Jx_fn(:,i) = (fp - fm)/(2*h);
%                     end
%                     err_x = Jx_f(1:m,1:m) - Jx_fn(1:m,1:m);
%                     err_th = Jth_f(1:m,1:d) - Jth_fn(1:m,1:d);
%                     disp(err_th); disp(err_x); pause
%                 end
%                 z = z + dzdt * dt;              % update states
%             end
%             if any(~isfinite(z)) || any(abs(z) > 1e12)
%                 fail = true; break;
%             else
%                 % -----------------------------
%                 if mem == 1
%                     Z(s,1:nvar) = z;            % State at t
%                 else
%                     if s >= ipr+1
%                         q_n(s-ipr) = z(m) - Z(m);
%                         J(s-ipr,1:d) = z(id) - Z(id);
%                     end
%                     Z(1,1:nvar) = z;
%                 end
%             end
%         end                                     % End of time loop
% 
%     case 4 %% C++: Runge Kutta implementation
%         data.f_c = th(1);       % field capacity (mm)
%         data.beta = th(2);      % shape exponent (-)
%         data.lp = th(3);        % evap limit fraction (-)
%         data.k_0 = th(4);       % near-surface recession (1/T)
%         data.uzl = th(5);       % threshold for k0 (mm)
%         data.k_1 = th(6);       % upper zone recession (1/T)
%         data.k_2 = th(7);       % lower zone recession (1/T)
%         data.perc = th(8);      % percolation max (mm/T)
%         data.T_tr = th(9);      % temperature threshold (°C)
%         data.f_dd = th(10);     % melt factor (mm/°C/T)
%         data.sfcf = th(11);     % snowfall correction factor (-)
%         data.cfr = th(12);      % refreezing factor (-)
% 
%         data.eps = eps;         % Dimensionless smoothing coef. [inactive]
%         data.rho = rho;         % Dimensionless smoothing coef. [inactive]
%         data.eps_s = eps_s;     % storage/flux smoothing (mm)
%         data.eps_t = eps_t;     % temperature smoothing (°C)
%         data.eps_x = eps_x;     % smoothing coefficient (mm)
%         data.ipr = ipr;         % Time to print
% 
%         [Z,q_n,J] = crr_hbv(mdl.tout,Z(1,1:nvar)',data,ode);
% end
% 
% if fail == true && (mem == 1)
%     Z(s:ns,:) = repmat(Z(s-1,:), ns-s+1, 1);
% end
% 
% if mem == 1
%     q_n = diff(Z(mdl.idx(1):mdl.idx(2),m));
%     switch nargout
%         case {2,3,4}
%             % diff appropriate elements of sensitivity state variables
%             J = diff(Z(mdl.idx(1):mdl.idx(2),id));
%     end
% end
% 
% if doRouting
%     q_norout = q_n;
%     q_n = filter(w_rt,1,q_n);
% end
% if nargout == 1
%     varargout = {q_n}; return
% else
%     if doRouting
%         % route discharge and derivatives wrt the ODE parameters
%         for j = 1:d
%             J(:,j) = filter(w_rt,1,J(:,j));
%         end
%         % derivative wrt routing parameter B
%         J(:,d+1) = filter(dw_rtdb,1,q_norout);
%     end
%     J = J .* reshape(Jth,1,[]);
%     if nargout == 2
%         varargout = {q_n,J}; return
%     elseif nargout == 3
%         varargout = {q_n,J,Jth}; return
%     else
%         varargout = {q_n,J,Jth,Z}; return
%     end
% end
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% %%                   Secondary functions listed below                    %%
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% 
% %% 1. Runge Kutta integration
% function [z,LTE] = rk2(t,z,h,th,data,m,d,eps_s,eps_t,eps_x)
% 
% dzdtE = hbv_aug_ode(t,z,th,data,m,d,eps_s,eps_t,eps_x)';      % Euler
% zE = z + h*dzdtE;
% dzdtH = hbv_aug_ode(t,zE,th,data,m,d,eps_s,eps_t,eps_x)';     % Heun
% z = z + 0.5*h*(dzdtE + dzdtH);                                % New z
% LTE = abs(zE - z);                                            % LTE
% 
% end
% 
% %% 2. HBV augmented ode with sensitivities as state variables
% function dzdt = hbv_aug_ode(t,z,th,data,m,d,eps_s,eps_t,eps_x)
% 
% x = z(1:m);                                         % x = [Swe Sm Uz Lz Q]
% S = reshape(z(m+1:m*(d+1)),m,d);                    % mxd sensitivity matrix
% [dxdt,dSdt] = hbv_odefcn(t,x,th,S,data,m,d, ...
%     eps_s,eps_t,eps_x);                             % compute dxdt & dSdt
% dzdt = [dxdt; dSdt(:)];                             % Repack single vector
% 
% end
% 
% function [dxdt,dSdt,Jth_f,Jx_f] = hbv_odefcn(t,x,th,Smat,data,m,d, ...
%     eps_s,eps_t,eps_x)
% % HBV ODE core + forward sensitivities (fully smooth incl. "safety" maxes)
% % States x = [Swe Sm Uz Lz Qcum]  (m = 5)
% % Params th(1:14) order:
% % 1 f_c, 2 beta, 3 lp, 4 k_0, 5 uzl, 6 k_1, 7 k_2, 8 perc, 9 T_tr,
% % 10 f_dd, 11 sfcf, 12 cfr, 13 maxbas, 14 cwh (unused)
% 
% dxdt = nan(m,1);
% 
% % -------------------------
% % Unpack parameters
% % -------------------------
% f_c = th(1);
% beta = th(2);
% lp = th(3);
% k_0 = th(4);
% uzl = th(5);
% k_1 = th(6);
% k_2 = th(7);
% perc = th(8);
% T_tr = th(9);
% f_dd = th(10);
% sfcf = th(11);
% cfr = th(12);
% 
% % -------------------------
% % Unpack states (RAW)
% % -------------------------
% Swe_raw = x(1);
% Sm_raw = x(2);
% Uz_raw = x(3);
% Lz_raw = x(4);
% 
% % -----------------------------
% % Forcing at current time index
% % -----------------------------
% it = floor(t) + 1;
% P = data.P(it,1);
% Ep = data.Ep(it,1);
% T = data.T(it,1);
% 
% % -----------------
% % Smooth primitives
% % -----------------
% smooth_pos = @(a,ep) 0.5*(a + sqrt(a.^2 + ep.^2));
% dsmooth_pos_da = @(a,ep) 0.5*(1 + a./sqrt(a.^2 + ep.^2));
% 
% smooth_min = @(A,B,ep) 0.5*(A + B - sqrt((A-B).^2 + ep.^2));
% dmin_dA = @(A,B,ep) 0.5*(1 - (A-B)./sqrt((A-B).^2 + ep.^2)); % d/dA
% dmin_dB = @(A,B,ep) 0.5*(1 + (A-B)./sqrt((A-B).^2 + ep.^2)); % d/dB
% 
% smooth_max = @(A,B,ep) 0.5*(A + B + sqrt((A-B).^2 + ep.^2));
% dmax_dA = @(A,B,ep) 0.5*(1 + (A-B)./sqrt((A-B).^2 + ep.^2)); % d/dA
% 
% % clamp01(z) = min(max(z,0),1) with consistent derivative
% clamp01 = @(z,ep) smooth_min( smooth_max(z,0,ep), 1, ep );
% 
% dclamp01_dz = @(z,ep) ...
%     dmin_dA( smooth_max(z,0,ep), 1, ep ) .* ...
%     dmax_dA( z, 0, ep );
% 
% % evap limiter: phiE = smooth_min(ratio,1)
% phiE_fun = @(ratio,ep) smooth_min(ratio,1,ep);
% dphiE_dratio = @(ratio,ep) dmin_dA(ratio,1,ep);
% 
% % ------------------------------------
% % Smooth-clamp storages used in fluxes
% % ------------------------------------
% Swe = smooth_pos(Swe_raw,eps_s);
% Sm = smooth_pos(Sm_raw,eps_s);
% Uz = smooth_pos(Uz_raw,eps_s);
% Lz = smooth_pos(Lz_raw,eps_s);
% 
% dSwe_dSweRaw = dsmooth_pos_da(Swe_raw,eps_s);
% dSm_dSmRaw = dsmooth_pos_da(Sm_raw,eps_s);
% dUz_dUzRaw = dsmooth_pos_da(Uz_raw,eps_s);
% dLz_dLzRaw = dsmooth_pos_da(Lz_raw,eps_s);
% 
% % ------------------------------------------
% % 1) Precip partition (smooth snow fraction)
% % ------------------------------------------
% uT = (T - T_tr) / eps_t;
% snow_fr = 0.5*(1 - tanh(uT));
% rain_fr = 1 - snow_fr;
% 
% sech2 = 1/(cosh(uT)^2);
% dsnow_dT_tr = 0.5 * sech2 / eps_t;  % as in your code
% 
% Ps = sfcf * P * snow_fr;
% Pr = P * rain_fr;
% 
% dPs_dsfcf = P * snow_fr;
% dPs_dT_tr = sfcf * P * dsnow_dT_tr;
% dPr_dT_tr = -P * dsnow_dT_tr;
% 
% % ------------------
% % 2) Melt & refreeze
% % ------------------
% Tm = T - T_tr;
% 
% posTm = smooth_pos(Tm,  eps_t);
% negTm = smooth_pos(-Tm, eps_t);
% 
% dposTm_dTm = dsmooth_pos_da(Tm,  eps_t);
% dnegTm_dTm = -dsmooth_pos_da(-Tm, eps_t);
% 
% Mpot = f_dd * posTm;
% Rpot = cfr * f_dd * negTm;
% 
% dMpot_dT_tr = -f_dd * dposTm_dTm;
% dMpot_df_dd = posTm;
% 
% dRpot_dT_tr = -(cfr*f_dd) * dnegTm_dTm;
% dRpot_df_dd = cfr * negTm;
% dRpot_dcfr = f_dd * negTm;
% 
% M = smooth_min(Swe, Mpot, eps_s);
% R = smooth_min(Pr,  Rpot, eps_s);
% 
% dM_dSwe = dmin_dA(Swe, Mpot, eps_s);
% dM_dMpot = dmin_dB(Swe, Mpot, eps_s);
% 
% dR_dPr = dmin_dA(Pr, Rpot, eps_s);
% dR_dRpot = dmin_dB(Pr, Rpot, eps_s);
% 
% dM_dT_tr = dM_dMpot * dMpot_dT_tr;
% dM_df_dd = dM_dMpot * dMpot_df_dd;
% 
% dR_dT_tr = dR_dPr * dPr_dT_tr + dR_dRpot * dRpot_dT_tr;
% dR_df_dd = dR_dRpot * dRpot_df_dd;
% dR_dcfr = dR_dRpot * dRpot_dcfr;
% 
% % ------------------------------------------------------
% % 3) Liquid reaching soil: Pliq = smooth_pos(Pr + M - R)
% % ------------------------------------------------------
% Win_raw = Pr + M - R;
% Pliq = smooth_pos(Win_raw, eps_s);
% dPliq_dWinraw = dsmooth_pos_da(Win_raw, eps_s);
% 
% % wrt Swe (via M)
% dPliq_dSwe = dPliq_dWinraw * dM_dSwe;
% 
% % wrt T_tr/f_dd/cfr
% dPliq_dT_tr = dPliq_dWinraw * (dPr_dT_tr + dM_dT_tr - dR_dT_tr);
% dPliq_df_dd = dPliq_dWinraw * (dM_df_dd - dR_df_dd);
% dPliq_dcfr = dPliq_dWinraw * (0 - dR_dcfr);
% 
% % ------------------------------------------------------
% % 4) Soil evap + recharge (FIXED: den and up are smooth)
% % ------------------------------------------------------
% 
% % --- Smooth "den = max(lp*f_c,eps_x)" with derivatives
% den0 = lp * f_c;
% den = smooth_max(den0,eps_x,eps_x);          % smooth max
% dDen_dDen0 = dmax_dA(den0,eps_x,eps_x);      % d den / d den0
% 
% % ratio = Sm/den
% ratio = Sm / den;
% 
% % derivatives of ratio (now consistent)
% dratio_dSm = 1/den;
% dratio_dDen = -Sm/(den^2);
% dratio_dDen0 = dratio_dDen * dDen_dDen0;
% 
% dratio_dlp = dratio_dDen0 * (f_c);
% dratio_dfc = dratio_dDen0 * (lp);
% 
% phiE = phiE_fun(ratio, eps_s);
% dphiE_drat = dphiE_dratio(ratio, eps_s);
% 
% E = Ep * phiE;
% 
% dE_dSm = Ep * dphiE_drat * dratio_dSm;
% dE_dlp = Ep * dphiE_drat * dratio_dlp;
% dE_dfc = Ep * dphiE_drat * dratio_dfc;
% 
% % --- Recharge: Re = Pliq * clamp01(Sm/fc_safe)^beta
% % Smooth "fc_safe = max(f_c, eps_x)" with derivatives
% fc_safe = smooth_max(f_c, eps_x, eps_x);
% dfcSafe_dfc = dmax_dA(f_c, eps_x, eps_x);
% 
% u_raw = Sm / fc_safe;
% 
% % derivatives of u_raw
% duraw_dSm = 1/fc_safe;
% duraw_dfc = -Sm/(fc_safe^2) * dfcSafe_dfc;
% 
% u = clamp01(u_raw, eps_s);
% du_duraw = dclamp01_dz(u_raw, eps_s);
% 
% du_dSm = du_duraw * duraw_dSm;
% du_dfc = du_duraw * duraw_dfc;
% 
% Re = Pliq * (u^beta);
% dRe_dPliq = u^beta;
% 
% % Smooth "up = max(u, eps_x)" with derivatives (for pow/log stability)
% up = smooth_max(u, eps_x, eps_x);
% dup_du = dmax_dA(u, eps_x, eps_x);
% 
% % Now derivatives of Re w.r.t Sm/fc/beta are consistent with up(u)
% dRe_dSm = Pliq * beta * up^(beta-1) * (dup_du * du_dSm);
% dRe_dfc = Pliq * beta * up^(beta-1) * (dup_du * du_dfc);
% dRe_dbeta = Pliq * up^beta * log(up);
% 
% % ---------------------------
% % 5) Response routine (Uz/Lz)
% % ---------------------------
% h = smooth_pos(Uz - uzl, eps_s);
% dh_dUz = dsmooth_pos_da(Uz - uzl, eps_s);
% dh_duzl = -dh_dUz;
% 
% q0 = k_0 * h;
% dq0_dUz = k_0 * dh_dUz;
% dq0_dk0 = h;
% dq0_duzl = k_0 * dh_duzl;
% 
% q1 = k_1 * Uz;
% dq1_dUz = k_1;
% dq1_dk1 = Uz;
% 
% perc_flux = smooth_min(perc, Uz, eps_s);
% dperc_dperc = dmin_dA(perc, Uz, eps_s);
% dperc_dUz = dmin_dB(perc, Uz, eps_s);
% 
% q2 = k_2 * Lz;
% dq2_dLz = k_2;
% dq2_dk2 = Lz;
% 
% Q = q0 + q1 + q2;
% 
% % ----------------------------
% % State ODEs (RAW state space)
% % ----------------------------
% dxdt(1) = Ps + R - M;
% dxdt(2) = Pliq - Re - E;
% dxdt(3) = Re - q0 - q1 - perc_flux;
% dxdt(4) = perc_flux - q2;
% dxdt(5) = Q;
% 
% % ---------------------------------------------------------------
% % Jacobian wrt states: build in SMOOTHED-space, then chain to RAW
% % ---------------------------------------------------------------
% Jx_f = zeros(m,m);
% 
% % f1
% Jx_f(1,1) = -dM_dSwe;
% 
% % f2
% Jx_f(2,1) = dPliq_dSwe - dRe_dPliq*dPliq_dSwe;
% Jx_f(2,2) = -dRe_dSm - dE_dSm;
% 
% % f3
% Jx_f(3,1) = dRe_dPliq*dPliq_dSwe;
% Jx_f(3,2) = dRe_dSm;
% Jx_f(3,3) = -(dq0_dUz + dq1_dUz + dperc_dUz);
% 
% % f4
% Jx_f(4,3) = dperc_dUz;
% Jx_f(4,4) = -dq2_dLz;
% 
% % f5
% Jx_f(5,3) = dq0_dUz + dq1_dUz;
% Jx_f(5,4) = dq2_dLz;
% 
% % Chain columns back to RAW states
% Jx_f(:,1) = Jx_f(:,1) * dSwe_dSweRaw;
% Jx_f(:,2) = Jx_f(:,2) * dSm_dSmRaw;
% Jx_f(:,3) = Jx_f(:,3) * dUz_dUzRaw;
% Jx_f(:,4) = Jx_f(:,4) * dLz_dLzRaw;
% 
% % ------------------
% % Parameter Jacobian
% % ------------------
% Jth_f = zeros(m,d);
% 
% % indices
% j_fc=1; j_beta=2; j_lp=3; j_k0=4; j_uzl=5; j_k1=6; j_k2=7; j_perc=8;
% j_T_tr=9; j_f_dd=10; j_sfcf=11; j_cfr=12;
% 
% % f1: Ps + R - M
% Jth_f(1,j_T_tr) = dPs_dT_tr + dR_dT_tr - dM_dT_tr;
% Jth_f(1,j_f_dd) = dR_df_dd - dM_df_dd;
% Jth_f(1,j_sfcf) = dPs_dsfcf;
% Jth_f(1,j_cfr) = dR_dcfr;
% 
% % f2: Pliq - Re - E
% Jth_f(2,j_T_tr) = dPliq_dT_tr    - dRe_dPliq*dPliq_dT_tr;
% Jth_f(2,j_f_dd) = dPliq_df_dd - dRe_dPliq*dPliq_df_dd;
% Jth_f(2,j_cfr) = dPliq_dcfr   - dRe_dPliq*dPliq_dcfr;
% 
% Jth_f(2,j_fc) = -dRe_dfc - dE_dfc;
% Jth_f(2,j_beta) = -dRe_dbeta;
% Jth_f(2,j_lp) = -dE_dlp;
% 
% % f3: Re - q0 - q1 - perc_flux
% Jth_f(3,j_T_tr) = dRe_dPliq*dPliq_dT_tr;
% Jth_f(3,j_f_dd) = dRe_dPliq*dPliq_df_dd;
% Jth_f(3,j_cfr) = dRe_dPliq*dPliq_dcfr;
% 
% Jth_f(3,j_fc) = dRe_dfc;
% Jth_f(3,j_beta) = dRe_dbeta;
% Jth_f(3,j_k0) = -dq0_dk0;
% Jth_f(3,j_uzl) = -dq0_duzl;
% Jth_f(3,j_k1) = -dq1_dk1;
% Jth_f(3,j_perc) = -dperc_dperc;
% 
% % f4: perc_flux - q2
% Jth_f(4,j_k2) = -dq2_dk2;
% Jth_f(4,j_perc) =  dperc_dperc;
% 
% % f5: Q
% Jth_f(5,j_k0) = dq0_dk0;
% Jth_f(5,j_uzl) = dq0_duzl;
% Jth_f(5,j_k1) = dq1_dk1;
% Jth_f(5,j_k2) = dq2_dk2;
% 
% % ---------------
% % Sensitivity ODE
% % ---------------
% dSdt = Jx_f * Smat + Jth_f;
% 
% end
% 
% % ----------------
% % helper functions
% % ----------------
% function [w,dwdb] = hbv_maxbas(b,L_max,eps_w)
% %HBV_MAXBAS_RELAXED  Differentiable (relaxed) maxbas routing.
% %
% % Qout(t) = sum_{k=0}^{L_max-1} w_k(b) * q_raw(t-k)
% % b is continuous analogue of maxbas (in time steps).
% % SYNOPSIS: [w,dwdb] = hbv_maxbas(b,L_max,eps_w)
% %   b       scalar > 0, continuous width parameter (time steps)
% %   L_max   integer >= 1, fixed kernel length [def: 60]
% %   eps_w   small smoothing constant [def: 1e-8]
% %   w       OUTPUT: L_max x 1 weights, sum(w) = 1
% %   dwdb    OUTPUT: L_max x 1 derivative dw/db
% 
% if nargin < 3 || isempty(eps_w), eps_w = 1e-8; end
% if nargin < 2 || isempty(L_max), L_max = 60; end
% 
% % ---- smooth primitives ----
% smooth_pos = @(a,ep) 0.5*(a + sqrt(a.^2 + ep.^2));
% dsmooth_pos_da = @(a,ep) 0.5*(1 + a./sqrt(a.^2 + ep.^2));
% smooth_max = @(A,B,ep) 0.5*(A + B + sqrt((A-B).^2 + ep.^2));
% dmax_dA = @(A,B,ep) 0.5*(1 + (A-B)./sqrt((A-B).^2 + ep.^2)); % d/dA
% smooth_abs = @(z,ep) sqrt(z.^2 + ep.^2);
% dabs_dz = @(z,ep) z ./ sqrt(z.^2 + ep.^2);
% 
% % ---- continuous triangle on k = 0..Lmax-1 ----
% k = (0:L_max-1)';               % causal lags
% b = max(b,eps_w);
% c = 0.5*b;                      % center
% half = 0.5*b;                   % half-width
% denom = smooth_max(half, ...
%     eps_w,eps_w);
% ddenom_db = 0.5 * dmax_dA( ...  % d(denom)/d(half) * d(half)/db
%     half,eps_w,eps_w);
% z = k - c;                      % distance from center
% az = smooth_abs(z,eps_w);
% tri = 1 - az./denom;            % unthresholded triangle
% % smooth cutoff at 0 (instead of hard max(0,tri))
% a = smooth_pos(tri,eps_w);
% dz_db = -0.5;                   % derivatives da/db
% daz_db = dabs_dz(z,eps_w) .* dz_db;
% dtri_db = -(daz_db)./denom + az./(denom.^2) .* ddenom_db;
% da_db = dsmooth_pos_da(tri,eps_w) .* dtri_db;
% 
% sa = sum(a);                    % normalize: w = a / sum(a)
% sda = sum(da_db);
% w = a / sa;
% dwdb = (da_db*sa - a*sda) / (sa^2);
% 
% end
% 
% %% 4. HBV: RHS-only for J(x)_f and J(th)_f check
% function [dxdt,dSdt,Jth_f,Jx_f] = rhs_only(t,z,th,data,m,d, ...
%     eps_s,eps_t,eps_x)
% 
% Smat = ones(m,d);
% % Call RHS
% [dxdt,dSdt,Jth_f,Jx_f] = hbv_odefcn(t,z,th,Smat,data,m,d, ...
%     eps_s,eps_t,eps_x);
% end
