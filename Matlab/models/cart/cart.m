% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cart/cart.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function [q,J,initial] = cart(p,mdl,data,initial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CART Evaluate the CaRT rainfall-runoff model and discharge Jacobian.
%
% SYNOPSIS:
% cart(...)
%
% INPUT ARGUMENTS:
% p : 11 CaRT calibration coordinates, including logarithmic entries.
% mdl : Model configuration, including solver/backend case.
% data : Basin forcing, elevation bands, and model settings.
% initial : Initial stores (mm), or [] to use the spin-up policy.
% data.forcing : N-by-4 [P,T,PET,training]; P/PET mm per interval, T deg C.
% data.elevation : Five band elevations (m).
% data.area : Five nonnegative band-area fractions summing to one.
% data.resolution : Intervals per day (1, 24, or 96).
% p order: ln(C), X4, ln(X3), ln(tau), X2, ln(rKp), Tsnow, ln(cm), fsc,
% ln(Gamma), ln(ks).
% X4 and tau use days; capacities use mm. P/PET are interval depths.
%
% OUTPUT ARGUMENTS:
% q : Interval discharge depths (mm).
% J : Discharge Jacobian with respect to the wrapper coordinates.
% initial : Detached baseline spin-up stores (mm).
%
% NOTES:
% CART Capacity-Residence-Time baseline after Xia (2026). P contains the 11
% parameters in the upstream calibration order. DATA.forcing is
% [P,T,PET,training]. DATA supplies five elevation bands. Groundwater
% recession uses midpoint RK2 (1), ode45 (2), Euler (3), or the native
% closed-form implementation (4). J differentiates discharge with the
% initial stores held fixed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 4
        initial = [];
    end
    f = double(data.forcing);
    z = double(data.elevation(:));
    a = double(data.area(:));
    p = double(p(:));
    res = data.resolution;
    assert(ismember(res,[1 24 96]),'Unsupported CaRT resolution.');
    assert(size(f,2)==4 && all(isfinite(f),'all') && any(f(:,4)), ...
        'CaRT requires finite forcing and a training mask.');
    assert(numel(p)==11 && numel(z)==5 && numel(a)==5, ...
        'CaRT requires 11 parameters and five elevation bands.');
    assert(isempty(initial) || numel(initial)==3, ...
        'CaRT requires three initial store values.');
    J = [];
    if mdl.mcode == 4
        run = @crr_cart_mex;
        if isfield(mdl,'config') && ...
                isfield(mdl.config,'cart_uniform_input') && ...
                mdl.config.cart_uniform_input
            run = @crr_cart_uniform_f986af2;
        end
        [q,J,initial] = run(f,z,a,p,initial,nargout>1,res);
        return
    end
    if isfield(mdl,'config') && ...
            isfield(mdl.config,'cart_uniform_input') && ...
            mdl.config.cart_uniform_input
        error('cart:scheme','Uniform input requires mcode 4.');
    end
    assert(ismember(mdl.mcode,1:3),'Unknown CaRT solver case.');
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'cartplus'));
    [q,initial] = daily(p,f,z,a,initial,res,mdl.mcode);
    if nargout > 1
        J = zeros(size(f,1),11);
        seed = eye(11);
        for j = 1:11
            y = daily(CartplusTangent(p,seed(:,j)),f,z,a, ...
                initial,res,mdl.mcode);
            J(:,j) = y.slope;
        end
    end
end

function [q,initial] = daily(p,f,z,a,initial,res,mcode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAILY Advance interval snow, production, and routing recursions.
%
% SYNOPSIS:
% daily(...)
%
% INPUT ARGUMENTS:
% p : 11 CaRT calibration coordinates, including logarithmic entries.
% f : Forcing [P, T, PET, training] in the caller layout.
% z : Five elevation-band heights (m).
% a : Five nonnegative band fractions, summing to one.
% initial : Initial stores (mm), or [] to use the spin-up policy.
% res : Forcing intervals per day: 1, 24, or 96.
% mcode : Backend/solver case selected by the model wrapper.
%
% OUTPUT ARGUMENTS:
% q : Interval discharge depths (mm).
% initial : Detached baseline spin-up stores (mm).
%
% NOTES:
% DAILY Apply the snow, production, hydrograph and routing recursions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = size(f,1);
    dt = 1/res;
    cap = exp(p(1));
    x4 = p(2);
    x3 = exp(p(3));
    tau = exp(p(4));
    kp = (4/9)^4*exp(p(6));
    offset = -exp(p(10))*(z-a.'*z)/1000;
    snow = zeros_like(p,5);
    liquid = zeros_like(p,n);
    for i = 1:n
        for j = 1:5
            warm = f(i,2)+offset(j)-p(7);
            cold = real(warm) < 0;
            if cold
                snow(j) = snow(j)+p(9)*f(i,1);
            end
            snow(j) = snow(j)-low(exp(p(11))*f(i,3),snow(j));
            melt = low(high(dt*exp(p(8))*warm,0),snow(j));
            snow(j) = snow(j)-melt;
            liquid(i) = liquid(i)+a(j)*(melt+~cold*f(i,1));
        end
    end
    if isempty(initial)
        [~,~,s0] = production(.7*cap,liquid,f,cap,kp*dt,min(1460*res,n));
    else
        s0 = initial(1);
    end
    [direct,recharge] = production(s0,liquid,f,cap,kp*dt,n);
    main = .9*hydrograph(direct,x4,false,res);
    bypass = .1*hydrograph(direct,x4,true,res);
    if isempty(initial)
        [~,r0] = quick(.5*x3,main,bypass,dt*p(5),x3,min(730*res,n),dt);
    else
        r0 = initial(2);
    end
    fast = quick(r0,main,bypass,dt*p(5),x3,n,dt);
    gbar = high(mean(recharge(f(:,4)~=0))/dt,1e-6);
    k = 1/(gbar*tau*tau);
    if isempty(initial)
        initial = [s0;r0;gbar*tau];
    end
    water = initial(3);
    q = zeros_like(p,n);
    for i = 1:n
        v = water+recharge(i);
        water = recession(v,k,dt,mcode);
        q(i) = fast(i)+v-water;
    end
end

function [d,g,s] = production(s,liquid,f,c,kp,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRODUCTION Evaluate production-store runoff, evaporation, and recharge.
%
% SYNOPSIS:
% production(...)
%
% INPUT ARGUMENTS:
% s : S used by this numerical operator.
% liquid : Liquid precipitation plus snowmelt on the selected time base.
% f : Forcing [P, T, PET, training] in the caller layout.
% c : Prepared model constants for the selected parameter vector.
% kp : Production-store percolation coefficient.
% n : Number of intervals or state entries in this helper.
%
% OUTPUT ARGUMENTS:
% d : Computed numerical result; tangent lanes are retained.
% g : Computed numerical result; tangent lanes are retained.
% s : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% PRODUCTION Update the GR4J production store and percolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    d = zeros_like(c,n);
    g = zeros_like(c,n);
    for i = 1:n
        pn = high(liquid(i)-f(i,3),0);
        en = high(f(i,3)-liquid(i),0);
        u = s/c;
        tp = tanh(pn/c);
        te = tanh(en/c);
        add = c*(1-u*u)*tp/(1+u*tp);
        evap = s*(2-u)*te/(1+(1-u)*te);
        v = s-evap+add;
        g(i) = v*(1-(1+kp*(v/c)^4)^(-.25));
        d(i) = pn-add;
        s = v-g(i);
    end
end

function y = hydrograph(x,x4,symmetric,res)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYDROGRAPH Route input depths through the selected unit hydrograph.
%
% SYNOPSIS:
% hydrograph(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
% x4 : First unit-hydrograph time base (days).
% symmetric : Whether to use the symmetric two-time-base hydrograph.
% res : Forcing intervals per day: 1, 24, or 96.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% HYDROGRAPH Convolve daily runoff with the original GR4J ordinates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = (6+7*symmetric)*res;
    s = zeros_like(x4,n+1);
    for i = 1:n
        u = (i/res)/x4;
        if ~symmetric
            s(i+1) = low(u,1)^2.5;
        elseif real(u) <= 1
            s(i+1) = .5*u^2.5;
        else
            s(i+1) = 1-.5*high(2-u,0)^2.5;
        end
    end
    y = conv(x,diff(s));
    y = y(1:length(real(x)));
end

function [q,r] = quick(r,main,bypass,x2,x3,n,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUICK Advance the quick store and direct routing branch.
%
% SYNOPSIS:
% quick(...)
%
% INPUT ARGUMENTS:
% r : R used by this numerical operator.
% main : Main used by this numerical operator.
% bypass : Bypass input after hydrograph routing.
% x2 : Exchange coefficient in the CaRT parameterization.
% x3 : Quick-store capacity scale (mm).
% n : Number of intervals or state entries in this helper.
% dt : Integration or observation interval duration (days).
%
% OUTPUT ARGUMENTS:
% q : Interval discharge depths (mm).
% r : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% QUICK Apply loss-only exchange and GR4J routing discharge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    q = zeros_like(x2,n);
    for i = 1:n
        loss = x2*(r/x3)^3.5;
        v = high(r+main(i)+.9*loss,0);
        r = v*(1+dt*(v/x3)^4)^(-.25);
        q(i) = v-r+high(bypass(i)+.1*loss,0);
    end
end

function w = recession(w,k,dt,mcode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECESSION Advance nonlinear groundwater recession for one interval.
%
% SYNOPSIS:
% recession(...)
%
% INPUT ARGUMENTS:
% w : Integration workspace or scalar operand in this helper.
% k : State index or release coefficient in this helper.
% dt : Integration or observation interval duration (days).
% mcode : Backend/solver case selected by the model wrapper.
%
% OUTPUT ARGUMENTS:
% w : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% RECESSION Advance the quadratic groundwater recession and tangents.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch mcode
        case 2
            opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
            rate = real(k);
            if isa(w,'CartplusTangent')
                dk = k.slope;
                rhs = @(~,x) [-rate*x(1)^2; ...
                    -2*rate*x(1)*x(2)-dk*x(1)^2];
                [~,out] = ode45(rhs,[0 dt],[real(w);w.slope],opts);
                w = CartplusTangent(out(end,1),out(end,2));
            else
                [~,out] = ode45(@(~,x) -rate*x^2,[0 dt],w,opts);
                w = out(end);
            end
        case {1,3}
            steps = max(8,ceil(8*dt*real(k)*real(w)));
            h = dt/steps;
            for j = 1:steps
                stage = w;
                if mcode == 1
                    stage = w-.5*h*k*w*w;
                end
                w = w-h*k*stage*stage;
            end
    end
end

function y = zeros_like(x,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZEROS_LIKE Allocate zeros preserving the prototype numerical type.
%
% SYNOPSIS:
% zeros_like(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
% n : Number of intervals or state entries in this helper.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% ZEROS_LIKE Allocate ordinary values or forward-mode sensitivities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y = zeros(n,1);
    if isa(x,'CartplusTangent')
        y = CartplusTangent(y,y);
    end
end

function y = high(a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HIGH Select the larger operand without changing its tangent type.
%
% SYNOPSIS:
% high(...)
%
% INPUT ARGUMENTS:
% a : First numerical operand or elevation-band fractions.
% b : Second numerical operand or rate coefficient.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% HIGH Select the greater primal value and carry its derivative.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if real(a) >= real(b)
        y = a;
    else
        y = b;
    end
end

function y = low(a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOW Select the smaller operand without changing its tangent type.
%
% SYNOPSIS:
% low(...)
%
% INPUT ARGUMENTS:
% a : First numerical operand or elevation-band fractions.
% b : Second numerical operand or rate coefficient.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% LOW Select the lesser primal value and carry its derivative.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if real(a) <= real(b)
        y = a;
    else
        y = b;
    end
end
