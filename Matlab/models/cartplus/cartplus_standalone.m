% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cartplus/cartplus_standalone.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function [out,states] = cartplus_standalone(data,p,mcode,fixed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CARTPLUS_STANDALONE Run sequential CaRT+ without the native MEX backend.
%
% SYNOPSIS:
% cartplus_standalone(...)
%
% INPUT ARGUMENTS:
% data : Basin forcing, elevation bands, and model settings.
% p : 42 physical CaRT+ parameters in native order.
% mcode : Backend/solver case selected by the model wrapper.
% fixed : Detached spin-up stores, or [] to determine them.
% data.forcing : N-by-4 [P,T,PET,training]; P/PET mm per interval, T deg C.
% data.elevation : Five band elevations (m).
% data.area : Five nonnegative band-area fractions summing to one.
% data.resolution : Intervals per day (1, 24, or 96).
%
% OUTPUT ARGUMENTS:
% out : Model flux array or numerical result in the caller contract.
% states : Captured spin-up states for a fixed-state derivative replay.
%
% NOTES:
% CARTPLUS_STANDALONE Advance the CaRT+ water balance without a MEX.
% data.forcing contains [P,T,PET,training]. P and PET are interval depths. p
% contains 42 physical parameters in the native CaRT+ parameter order. fixed
% optionally supplies seven baseline spin-up stores. out contains Q, ET,
% external flux, and storage change per interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f = data.forcing;
    dt = 1/data.resolution;
    n = size(f,1);
    spin = min(n,1460*data.resolution);
    train = logical(f(:,4));
    if ~any(train)
        error('cartplus:training','CaRT+ requires training intervals.');
    end
    P = f(:,1);
    T = f(:,2);
    PET = f(:,3);
    area = data.area(:);
    height = data.elevation(:);
    vegetation = logistic((T-p(23))/p(24));
    average = mean(vegetation(1:min(n,365*data.resolution)));
    filtered = memory(vegetation-average,exp(-dt/p(25)), ...
        2048*data.resolution);
    demand = PET.*(1-p(22)*(1-clip(average+filtered,0,1)));
    mean_rain = max(mean(P(train)),.1*dt);

    snow = zero_store(p,5,1);
    rain = zero_store(p,n,1);
    sublimation = zero_store(p,n,1);
    snow_extra = zero_store(p,n,1);
    snow_change = zero_store(p,n,1);
    offset = p(6)*(height-area.'*height)/1000;
    for i = 1:n
        before = area.'*snow;
        warmth = T(i)-offset-p(1);
        cold = real(warmth)<0;
        snow = snow+cold*p(3)*P(i);
        sub = lesser(snow,p(4)*PET(i));
        snow = snow-sub;
        melt = lesser(snow,(~cold).*(p(2)*dt*warmth+p(5)*PET(i)));
        snow = snow-melt;
        rain(i) = area.'*(melt+(~cold)*P(i));
        sublimation(i) = area.'*sub;
        snow_extra(i) = (area.'*cold)*P(i)*(p(3)-1);
        snow_change(i) = area.'*snow-before;
    end

    drainage = -expm1(-dt/p(8));
    [surface,pre,pre_pet] = surface_series(.5*p(7), ...
        rain(1:spin),demand(1:spin),p(7),drainage);
    if ~isempty(fixed)
        surface = fixed(1);
    end
    surface_start = surface;
    [~,rain,demand,surface_et,change] = surface_series( ...
        surface,rain,demand,p(7),drainage);
    pre = infiltrate(pre,pre_pet,mean_rain,p(12),p(13), ...
        data.floor,data.resolution);
    [rain,excess] = infiltrate(rain,demand,mean_rain,p(12),p(13), ...
        data.floor,data.resolution);

    soil = .7*p(9);
    near = .3*p(15);
    wet = .3*p(26);
    if isempty(fixed)
        for i = 1:spin
            soil = pareto(soil,pre(i),pre_pet(i),p(9),p(10),p(11), ...
                dt,true);
            near = zone_step(near,pre(i),pre_pet(i),p(15),p(16), ...
                dt,mcode);
            wet = zone_step(wet,pre(i),pre_pet(i),p(26),p(27), ...
                dt,mcode);
        end
    else
        soil = fixed(2);
        near = fixed(3);
        wet = fixed(4);
    end
    soil_start = soil;
    near_start = near;
    wet_start = wet;
    fast = zero_store(p,n,1);
    recharge = zero_store(p,n,1);
    et = zero_store(p,n,1);
    near_q = zero_store(p,n,1);
    saturation = zero_store(p,n,1);
    wetness = zero_store(p,n,1);
    upland = 1-p(14);
    for i = 1:n
        before = upland*soil+p(14)*near;
        [soil,runoff,percolation,evap] = pareto(soil,rain(i), ...
            demand(i),p(9),p(10),p(11),dt,true);
        [near,zone_runoff,zone_q,zone_et] = zone_step(near,rain(i), ...
            demand(i),p(15),p(16),dt,mcode);
        wet = zone_step(wet,rain(i),demand(i),p(26),p(27),dt,mcode);
        fast(i) = upland*runoff+p(14)*zone_runoff+excess(i);
        recharge(i) = upland*percolation;
        et(i) = upland*evap+p(14)*zone_et;
        near_q(i) = p(14)*zone_q;
        saturation(i) = clip(soil/p(9),0,1);
        wetness(i) = wet/p(26);
        change(i) = change(i)+upland*soil+p(14)*near-before;
    end

    depression = .5*p(17);
    if isempty(fixed)
        for i = 1:spin
            depression = pareto(depression,fast(i), ...
                greater(demand(i)-et(i),0),p(17),p(19),p(18), ...
                dt,false);
        end
    else
        depression = fixed(5);
    end
    depression_start = depression;
    for i = 1:n
        before = depression;
        [depression,fast(i),add,evap] = pareto(depression,fast(i), ...
            greater(demand(i)-et(i),0),p(17),p(19),p(18),dt,false);
        recharge(i) = recharge(i)+add;
        et(i) = et(i)+evap;
        change(i) = change(i)+depression-before;
    end

    equilibrium = greater(mean(recharge(train))/dt,1e-6)*p(20);
    groundwater_start = equilibrium;
    if ~isempty(fixed)
        groundwater_start = fixed(6);
    end
    groundwater = groundwater_start;
    recession = dt/(equilibrium*p(20));
    base_q = zero_store(p,n,1);
    connected = zero_store(p,n,1);
    for i = 1:n
        before = groundwater;
        groundwater = groundwater+recharge(i);
        x = p(21)*greater(demand(i)-et(i),0);
        x = x/(groundwater+groundwater_start);
        if real(x) < 1e-3
            uptake = groundwater*(x-.5*x*x+x*x*x/6);
        else
            uptake = -groundwater*expm1(-x);
        end
        groundwater = groundwater-uptake;
        next = groundwater/(1+recession*groundwater);
        base_q(i) = groundwater-next;
        groundwater = next;
        change(i) = change(i)+groundwater-before;
        et(i) = et(i)+uptake;
        signal = wetness(i)-p(29)+p(31)*saturation(i) ...
            +p(32)*(groundwater/groundwater_start-1);
        connected(i) = p(28)+(1-p(28))*logistic(signal/p(30));
    end

    direct = p(33)*connected.*fast;
    routed = hydrograph(fast-direct,p(34),dt,false);
    bypass = hydrograph(direct,p(34),dt,true);
    change = change+fast-routed-bypass;
    exterior = dt*p(37);
    routed = routed+.9*exterior;
    bypass = bypass+.1*exterior;
    kelvin = max(T,0)+273.15;
    viscosity = -log(10)*247.8*(1./(kelvin-140)-1/153.15);
    connection = log(connected)+viscosity ...
        -p(40)*logistic((p(38)-T)/p(39));
    quick = .5*p(35);
    if isempty(fixed)
        for i = 1:min(n,730*data.resolution)
            quick = quick_step(quick,routed(i),bypass(i), ...
                connection(i),p(35),p(36),p(41),dt);
        end
    else
        quick = fixed(7);
    end
    states = [surface_start,soil_start,near_start,wet_start, ...
        depression_start,groundwater_start,quick];
    out = zero_store(p,n,4);
    for i = 1:n
        before = quick;
        [quick,flow,exchange] = quick_step(quick,routed(i),bypass(i), ...
            connection(i),p(35),p(36),p(41),dt);
        flow = flow+base_q(i)+near_q(i);
        cap = dt*p(42);
        x = greater(flow,0)/cap;
        loss = cap*saturated(x);
        out(i,:) = [flow-loss,et(i)+surface_et(i)+sublimation(i), ...
            snow_extra(i)+exterior+exchange-loss, ...
            snow_change(i)+change(i)+quick-before];
    end
end

function y = logistic(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOGISTIC Evaluate the logistic response with forward derivatives.
%
% SYNOPSIS:
% logistic(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% LOGISTIC Evaluate the CaRT+ logistic response.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y = 1./(1+exp(-x));
end

function y = lesser(a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LESSER Select the smaller primal branch and retain its derivatives.
%
% SYNOPSIS:
% lesser(...)
%
% INPUT ARGUMENTS:
% a : First numerical operand or elevation-band fractions.
% b : Second numerical operand or rate coefficient.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% LESSER Select the smaller value using the unperturbed ordering.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mask = real(a)<real(b);
    y = a.*mask+b.*(~mask);
end

function y = greater(a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GREATER Select the larger primal branch and retain its derivatives.
%
% SYNOPSIS:
% greater(...)
%
% INPUT ARGUMENTS:
% a : First numerical operand or elevation-band fractions.
% b : Second numerical operand or rate coefficient.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% GREATER Select the larger value using the unperturbed ordering.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mask = real(a)>real(b);
    y = a.*mask+b.*(~mask);
end

function y = clip(x,lo,hi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLIP Bound the numerical operand to its admissible interval.
%
% SYNOPSIS:
% clip(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
% lo : Lower numerical bound.
% hi : Upper numerical bound.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% CLIP Restrict a model value to its admissible interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y = lesser(greater(x,lo),hi);
end

function out = memory(values,decay,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEMORY Apply the finite exponentially weighted memory kernel.
%
% SYNOPSIS:
% memory(...)
%
% INPUT ARGUMENTS:
% values : Numerical values or vector of state entries.
% decay : Memory-retention coefficient (-).
% window : Finite memory-kernel length in observation intervals.
%
% OUTPUT ARGUMENTS:
% out : Model flux array or numerical result in the caller contract.
%
% NOTES:
% MEMORY Apply the finite causal exponential filter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    out = zero_store(values,size(values));
    previous = 0;
    tail = (1-decay)*decay^window;
    for i = 1:numel(real(values))
        previous = decay*previous+(1-decay)*values(i);
        if i > window
            previous = previous-tail*values(i-window);
        end
        out(i) = previous;
    end
end

function [s,q,demand,et,change] = surface_series(s,rain,pet,cap,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SURFACE_SERIES Advance surface storage, overflow, drainage, and
% evaporation.
%
% SYNOPSIS:
% surface_series(...)
%
% INPUT ARGUMENTS:
% s : S used by this numerical operator.
% rain : Liquid input depth/rate on this operator time base.
% pet : Potential evaporation depth/rate on this operator time base.
% cap : Store capacity (mm).
% k : State index or release coefficient in this helper.
%
% OUTPUT ARGUMENTS:
% s : Computed numerical result; tangent lanes are retained.
% q : Interval discharge depths (mm).
% demand : Computed numerical result; tangent lanes are retained.
% et : Computed numerical result; tangent lanes are retained.
% change : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% SURFACE_SERIES Advance surface evaporation, overflow, and drainage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = numel(real(rain));
    q = zero_store(k,n,1);
    demand = zero_store(k,n,1);
    et = zero_store(k,n,1);
    change = zero_store(k,n,1);
    for i = 1:n
        before = s;
        s = s+rain(i);
        et(i) = lesser(pet(i),s);
        s = s-et(i);
        overflow = greater(s-cap,0);
        s = s-overflow;
        drain = k*s;
        s = s-drain;
        q(i) = overflow+drain;
        demand(i) = greater(pet(i)-et(i),0);
        change(i) = s-before;
    end
end

function [rain,excess] = infiltrate(rain,pet,average,cap,shape,floor,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFILTRATE Partition liquid input into infiltration and excess runoff.
%
% SYNOPSIS:
% infiltrate(...)
%
% INPUT ARGUMENTS:
% rain : Liquid input depth/rate on this operator time base.
% pet : Potential evaporation depth/rate on this operator time base.
% average : Mean training rainfall on the forcing time base.
% cap : Store capacity (mm).
% shape : Dimensionless storage/infiltration shape parameter.
% floor : Minimum infiltration capacity in daily-rate units (mm/day).
% r : R used by this numerical operator.
%
% OUTPUT ARGUMENTS:
% rain : Computed numerical result; tangent lanes are retained.
% excess : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% INFILTRATE Partition runoff using antecedent water availability.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    decay = exp(-1/(3*r));
    history = memory(rain,decay,36*r);
    previous = [0;history(1:end-1)];
    initial = decay.^(0:numel(real(rain))-1)';
    cap = greater(cap*exp(-shape*(previous/average+initial)),floor)/r;
    available = greater(rain-pet,0);
    excess = available-cap.*saturated(available./cap);
    rain = rain-excess;
end

function y = saturated(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SATURATED Limit the soil-saturation response to the admissible range.
%
% SYNOPSIS:
% saturated(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% SATURATED Evaluate one minus the exponential with a small-input branch.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y = -expm1(-x);
    small = real(x)<1e-3;
    y(small) = x(small)-.5*x(small).^2;
end

function [s,q,recharge,et] = pareto(s,rain,pet,cap,k,shape,dt,net)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARETO Advance the Pareto production-store water balance.
%
% SYNOPSIS:
% pareto(...)
%
% INPUT ARGUMENTS:
% s : S used by this numerical operator.
% rain : Liquid input depth/rate on this operator time base.
% pet : Potential evaporation depth/rate on this operator time base.
% cap : Store capacity (mm).
% k : State index or release coefficient in this helper.
% shape : Dimensionless storage/infiltration shape parameter.
% dt : Integration or observation interval duration (days).
% net : Whether to remove interception before the production update.
%
% OUTPUT ARGUMENTS:
% s : Computed numerical result; tangent lanes are retained.
% q : Interval discharge depths (mm).
% recharge : Computed numerical result; tangent lanes are retained.
% et : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% PARETO Advance a production store and return runoff, recharge, and ET.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    water = rain;
    potential = pet;
    interception = 0;
    if net
        water = greater(rain-pet,0);
        potential = greater(pet-rain,0);
        interception = lesser(rain,pet);
    end
    deficit = greater(1-clip(s/cap,0,1),1e-7);
    exponent = shape+1;
    remaining = deficit^(1/exponent)-water/(exponent*cap);
    retained = clip(cap*(deficit-greater(remaining,0)^exponent), ...
        0,water);
    response = tanh(potential/cap);
    dry = deficit^2;
    evap = lesser(cap*(1-dry)*response/(1+dry/deficit*response),s);
    s = s+retained-evap;
    recharge = -s*expm1(-.25*dt*log1p(k*(s/cap)^4));
    s = s-recharge;
    q = water-retained;
    et = evap+interception;
end

function [s,q,flow,et] = zone_step(s,rain,pet,cap,days,dt,mcode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZONE_STEP Advance near-stream storage and its interval fluxes.
%
% SYNOPSIS:
% zone_step(...)
%
% INPUT ARGUMENTS:
% s : S used by this numerical operator.
% rain : Liquid input depth/rate on this operator time base.
% pet : Potential evaporation depth/rate on this operator time base.
% cap : Store capacity (mm).
% days : Residence time (days).
% dt : Integration or observation interval duration (days).
% mcode : Backend/solver case selected by the model wrapper.
%
% OUTPUT ARGUMENTS:
% s : Computed numerical result; tangent lanes are retained.
% q : Interval discharge depths (mm).
% flow : Computed numerical result; tangent lanes are retained.
% et : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% ZONE_STEP Integrate a near-stream or wetness zone for one interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    net = (rain-pet)/cap;
    rate = dt/days;
    inputs = [s/cap,net,rate];
    if isa(inputs,'CartplusTangent')
        local = cartplus_zone_solver([inputs.value,dt,1],mcode);
        local = reshape(local,2,4);
        value = CartplusTangent(local(:,1),local(:,2:4)*inputs.slope.');
    elseif isreal(inputs)
        value = cartplus_zone_solver([inputs,dt,0],mcode);
    else
        tangent = cartplus_zone_solver([real(inputs),dt,1],mcode);
        tangent = reshape(tangent,2,4);
        value = tangent(:,1)+1i*tangent(:,2:4)*imag(inputs(:));
    end
    next = cap*value(1);
    flow = cap*rate*value(2);
    change = next-s;
    q = 0;
    if real(net) > 0
        q = rain-pet-flow-change;
    end
    s = next;
    et = rain-q-flow-change;
end

function q = hydrograph(input,days,dt,symmetric)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYDROGRAPH Route input depths through the selected unit hydrograph.
%
% SYNOPSIS:
% hydrograph(...)
%
% INPUT ARGUMENTS:
% input : Forcing/storage input for the current operator.
% days : Residence time (days).
% dt : Integration or observation interval duration (days).
% symmetric : Whether to use the symmetric two-time-base hydrograph.
%
% OUTPUT ARGUMENTS:
% q : Interval discharge depths (mm).
%
% NOTES:
% HYDROGRAPH Convolve interval depths with the discrete routing weights.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    scale = 1+symmetric;
    count = ceil(real(days)*scale/dt);
    if dt == 1
        count = 6+7*symmetric;
    end
    lag = (0:count)'*dt/days;
    if symmetric
        cumulative = .5*lag.^2.5;
        right = real(lag)>1;
        cumulative(right) = 1-.5*greater(2-lag(right),0).^2.5;
    else
        cumulative = lesser(lag,1).^2.5;
    end
    q = conv(input,diff(cumulative));
    q = q(1:numel(real(input)));
end

function [s,q,exchange] = quick_step(s,routed,bypass,signal,cap,k,b,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUICK_STEP Advance quick routing, exchange, and bypass discharge.
%
% SYNOPSIS:
% quick_step(...)
%
% INPUT ARGUMENTS:
% s : S used by this numerical operator.
% routed : Hydrograph-routed input to the quick store.
% bypass : Bypass input after hydrograph routing.
% signal : Signal used by this numerical operator.
% cap : Store capacity (mm).
% k : State index or release coefficient in this helper.
% b : Second numerical operand or rate coefficient.
% dt : Integration or observation interval duration (days).
%
% OUTPUT ARGUMENTS:
% s : Computed numerical result; tangent lanes are retained.
% q : Interval discharge depths (mm).
% exchange : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% QUICK_STEP Advance the quick store and account for external exchange.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    flux = dt*k*(greater(s,0)/cap)^3.5;
    available = greater(s+routed+.9*flux,0);
    direct = greater(bypass+.1*flux,0);
    factor = signal+b*log(greater(available,1e-3)/cap);
    rate = .25*dt*cap^-4*exp(6*tanh(factor/6));
    next = available*exp(-.25*log1p(4*rate*available^4));
    q = available-next+direct;
    exchange = available-s-routed+direct-bypass;
    s = next;
end

function out = zero_store(prototype,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZERO_STORE Allocate zeros preserving the prototype numerical type.
%
% SYNOPSIS:
% zero_store(...)
%
% INPUT ARGUMENTS:
% prototype : Numerical prototype preserving tangent/complex scalar type.
% varargin : Additional MATLAB operands or dimension arguments.
%
% OUTPUT ARGUMENTS:
% out : Model flux array or numerical result in the caller contract.
%
% NOTES:
% ZERO_STORE Allocate values and, when seeded, their tangent storage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    out = zeros(varargin{:});
    if isa(prototype,'CartplusTangent')
        out = CartplusTangent(out);
    end
end
