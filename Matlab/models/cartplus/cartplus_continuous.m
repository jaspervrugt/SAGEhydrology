% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cartplus/cartplus_continuous.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function [q,jac,balance] = cartplus_continuous(data,p,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CARTPLUS_CONTINUOUS Integrate fully coupled CaRT+ and its water-balance
% ledger.
%
% SYNOPSIS:
% cartplus_continuous(...)
%
% INPUT ARGUMENTS:
% data : Basin forcing, elevation bands, and model settings.
% p : 42 physical CaRT+ parameters in native order.
% options : Optional solver controls and regularization settings.
% data.forcing : N-by-4 [P,T,PET,training]; P/PET mm per interval, T deg C.
% data.elevation : Five band elevations (m).
% data.area : Five nonnegative band-area fractions summing to one.
% data.resolution : Intervals per day (1, 24, or 96).
%
% OUTPUT ARGUMENTS:
% q : Interval discharge depths (mm).
% jac : Discharge Jacobian in the selected parameter directions.
% balance : Water-balance residual per observation interval (mm).
%
% NOTES:
% CARTPLUS_CONTINUOUS Fully coupled CaRT+ with DOP853 and 42 tangents.
% Forcing columns are [P, T, PET, training]. Depths use mm, times days.
% Optional data.initial holds detached [N1..N5,H,S,N,B,A,J,V,W,R,Z].
% Optional options.seed scales physical sensitivities to other coordinates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    arguments
        data (1,1) struct
        p (42,1) double {mustBeFinite}
        options.refinement (1,1) double {mustBeInteger,mustBePositive} = 16
        options.rtol (1,1) double {mustBePositive,mustBeFinite} = 1e-8
        options.atol (1,1) double {mustBePositive,mustBeFinite} = 1e-10
        options.floor (1,1) double {mustBeNonnegative,mustBeFinite} = 0
        options.availability (1,1) double {mustBePositive} = .02
        options.transition (1,1) double {mustBePositive} = .01
        options.seed (42,1) double {mustBeFinite} = ones(42,1)
        options.columns (1,:) double {mustBeInteger,mustBePositive} = 1:42
        options.sensitivities (1,1) logical = true
    end
    f = double(data.forcing);
    a = double(data.area(:));
    z = double(data.elevation(:));
    n = size(f,1);
    if isfield(data,'interval_days')
        dt = double(data.interval_days(:)).*ones(n,1);
    else
        dt = ones(n,1)/data.resolution;
    end
    if size(f,2) ~= 4 || numel(dt) ~= n
        error('cartplus:input','Invalid forcing shape or interval count.');
    end
    train = f(:,4) == 1;
    if ~any(train) || ~all(ismember(f(:,4),[0,1])) || ...
            any(f(:,[1,3]) < 0,'all') || ...
            any(~isfinite(f),'all') || numel(a) ~= 5 || numel(z) ~= 5 || ...
            any(~isfinite([z;a;dt])) || any(a < 0) || any(dt <= 0) || ...
            abs(sum(a)-1) > 1e-10 || options.transition > 1 || ...
            ~isfinite(options.availability)
        error('cartplus:input','Invalid forcing, bands or solver options.');
    end
    positive = [2,6:12,15:20,24:27,30,34,35,39,42];
    fractions = p([14,22,28,33]);
    if any(p(positive) <= 0) || any(fractions < 0 | fractions > 1) || ...
            any(p([3,4,5,13,21,31,32,37,40,41]) < 0) || p(36) > 0
        error('cartplus:input','Invalid physical parameters.');
    end
    m = options.refinement;
    ns = 14+3*m;
    fixed = [];
    if isfield(data,'initial') && ~isempty(data.initial)
        fixed = double(data.initial(:));
        if numel(fixed) ~= ns || any(~isfinite(fixed)) || ...
                any(fixed < 0) || fixed(12) > 1 || ...
                any(fixed(6:10) > p([7,9,15,26,17]))
            error('cartplus:input','Invalid initial stores.');
        end
    end

    % The imaginary lanes carry simultaneous directional derivatives.
    % Branches use the first (real) lane, including at exact ties.
    epsilon = 1e-25;
    tangents = nargout >= 2 && options.sensitivities;
    columns = options.columns;
    if any(columns > 42) || numel(unique(columns)) ~= numel(columns)
        error('cartplus:input','Invalid sensitivity columns.');
    end
    if tangents
        seeds = diag(options.seed);
        p = p+1i*epsilon*[zeros(42,1),seeds(:,columns)];
    end
    lanes = size(p,2);
    [rkA,rkB,rkE3,rkE5] = cartplus_dop853();
    eta9 = exit_fractions(m,false);
    eta1 = exit_fractions(m,true);
    lambda = m./p(34,:);
    upland = 1-p(14,:);
    trans = options.transition;
    avail = options.availability;
    train_days = sum(dt(train));
    mean_rain = max(sum(f(train,1))/train_days,.1);
    offset = (z-a.'*z)/1000;
    vegetation = sigmoid((f(:,2)-p(23,:))./p(24,:));
    starts = [0;cumsum(dt(1:end-1))];
    year_weights = min(dt,max(365-starts,0));
    veg0 = (year_weights.'*vegetation)/sum(year_weights);
    pilot_k = zeros(13*lanes,13,'like',p);
    full_k = zeros((ns+7)*lanes,13,'like',p);
    y = initialize(true,0);
    recharge = zeros(1,lanes);
    pilot_end = find(train,1,'last');
    for i = 1:pilot_end
        before = y(13,:);
        y = advance(y,i,true,0);
        if train(i)
            recharge = recharge+y(13,:)-before;
        end
    end
    reference = upper(recharge/train_days,1e-6).*p(20,:);
    y = initialize(false,reference);
    q = zeros(n,1);
    jac = zeros(n,numel(columns)*tangents);
    balance = zeros(n,1);
    weights = [a;1;real(upland(1));real(p(14,1));0;1;0;0;1;1; ...
        ones(3*m,1)];
    for i = 1:n
        before = y(ns+1,:);
        old_balance = ledger(y);
        y = advance(y,i,false,reference);
        delta = y(ns+1,:)-before;
        q(i) = real(delta(1));
        if tangents
            jac(i,:) = imag(delta(2:end))/epsilon;
        end
        balance(i) = ledger(y)-old_balance;
    end

    function total = ledger(y)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LEDGER Evaluate the water-balance ledger from storage and integrals.
    %
    % SYNOPSIS:
    % ledger(...)
    %
    % INPUT ARGUMENTS:
    % y : Numerical state/operand; augmented with sensitivities when
    % enabled.
    %
    % OUTPUT ARGUMENTS:
    % total : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = real(y(:,1));
        total = weights.'*x(1:ns)-x(ns+4)-x(ns+3)+x(ns+2)+ ...
            x(ns+6)+x(ns+1);
    end

    function y = initialize(pilot,reference)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIALIZE Initialize water stores and detach the prescribed spin-up.
    %
    % SYNOPSIS:
    % initialize(...)
    %
    % INPUT ARGUMENTS:
    % pilot : Whether to run the reference-recharge initialization pass.
    % reference : Pilot mean recharge or reference storage used by the
    % caller.
    %
    % OUTPUT ARGUMENTS:
    % y : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y = zeros(13,lanes);
        count = 12;
        if ~pilot
            y = zeros(ns+7,lanes);
            count = ns;
        end
        if ~isempty(fixed)
            y(1:count,:) = repmat(fixed(1:count),1,lanes);
            return
        end
        y(6:12,:) = [.5*p(7,:);.7*p(9,:);.3*p(15,:);.3*p(26,:); ...
            .5*p(17,:);mean_rain*ones(1,lanes);veg0];
        if ~pilot
            y(13,:) = reference;
            y(14,:) = .5*p(35,:);
        end
    end

    function y = advance(y,i,pilot,reference)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADVANCE Advance routing memory and return its outlet depth.
    %
    % SYNOPSIS:
    % advance(...)
    %
    % INPUT ARGUMENTS:
    % y : Numerical state/operand; augmented with sensitivities when
    % enabled.
    % i : Observation/state index.
    % pilot : Whether to run the reference-recharge initialization pass.
    % reference : Pilot mean recharge or reference storage used by the
    % caller.
    %
    % OUTPUT ARGUMENTS:
    % y : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        precip = f(i,1)/dt(i);
        temp = f(i,2);
        pet = f(i,3)/dt(i);
        warm = temp-offset*p(6,:)-p(1,:);
        cold = real(warm(:,1)) < 0;
        snowfall = cold*(p(3,:)*precip);
        melt_demand = (warm.*p(2,:)+p(5,:)*pet).*(~cold);
        liquid_rain = precip*sum(a(~cold));
        effective = precip+(p(3,:)-1)*precip*sum(a(cold));
        tk = max(temp,0)+273.15;
        thermal = -log(10)*247.8*(1/(tk-140)-1/153.15)- ...
            p(40,:).*sigmoid((p(38,:)-temp)./p(39,:));
        sublimation = p(4,:)*pet;
        rhs = @(x) rates(x,pilot,reference,pet,snowfall,melt_demand, ...
            liquid_rain,effective,thermal,sublimation,vegetation(i,:));
        y = integrate(y,dt(i),rhs,pilot);
        if pilot
            count = 12;
        else
            count = ns;
        end
        if any(real(y(1:count,1)) < -max(1e-7,100*options.atol)) || ...
                any(~isfinite(y),'all')
            error('cartplus:state','Invalid coupled CaRT+ state.');
        end
    end

    function dy = rates(x,pilot,ref,pet,snowfall,melt_demand, ...
            liquid_rain,effective,thermal,sublimation,veg_target)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RATES Evaluate state and cumulative-flux derivatives.
    %
    % SYNOPSIS:
    % rates(...)
    %
    % INPUT ARGUMENTS:
    % x : Numerical operand; scalar, vector, or forward-mode scalar.
    % pilot : Whether to run the reference-recharge initialization pass.
    % ref : Ref used by this numerical operator.
    % pet : Potential evaporation depth/rate on this operator time base.
    % snowfall : Band snowfall depth/rate on the selected time base.
    % melt_demand : Potential snowmelt rate for each band.
    % liquid_rain : Liquid precipitation rate supplied to the snow
    % calculation.
    % effective : Effective used by this numerical operator.
    % thermal : Thermal used by this numerical operator.
    % sublimation : Snow sublimation depth/rate on the selected time base.
    % veg_target : Veg target used by this numerical operator.
    %
    % OUTPUT ARGUMENTS:
    % dy : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = real(x(:,1))+1i*imag(x);
        if ~tangents
            x = real(x);
        end
        water = upper(x(1:5,:),0)/avail;
        sub = lower(repmat(sublimation,5,1),water);
        melt = lower(upper(water-sub,0),melt_demand);
        liquid = liquid_rain+a.'*melt;
        demand = pet*(1-p(22,:).*(1-clip(x(12,:))));
        eh = lower(demand,upper(x(6,:),0)/avail);
        qh = upper(x(6,:),0)./p(8,:);
        fh = qh+upper(liquid-eh-qh,0).* ...
            clip((x(6,:)./p(7,:)-(1-trans))/trans);
        ep = upper(demand-eh,0);
        cap = upper(p(12,:).*exp(-p(13,:).* ...
            upper(x(11,:),0)/mean_rain),max(options.floor,1e-12));
        qx = excess(upper(fh-ep,0),cap);
        ps = fh-qx;
        [qs,es,ds] = pareto(x(7,:),ps,ep,p(9,:),p(11,:),p(10,:),true);
        [qn,en,fn] = zone(x(8,:),ps,ep,p(15,:),p(16,:));
        [qw,ew,fw] = zone(x(9,:),ps,ep,p(26,:),p(27,:));
        ez = upland.*es+p(14,:).*en;
        fast = qx+upland.*qs+p(14,:).*qn;
        [qa,ea,da] = pareto(x(10,:),fast,upper(ep-ez,0), ...
            p(17,:),p(18,:),p(19,:),false);
        dr = upland.*ds+da;
        upstream = [snowfall-sub-melt;liquid-eh-fh;ps-qs-es-ds; ...
            ps-qn-en-fn;ps-qw-ew-fw;fast-qa-ea-da;(fh-x(11,:))/3; ...
            (veg_target-x(12,:))./p(25,:)];
        if pilot
            dy = [upstream;dr];
            return
        end
        gw = upper(x(13,:),0);
        r = upper(x(14,:),0);
        qg = gw.^2./(ref.*p(20,:));
        eg = p(21,:).*gw.*upper(ep-ez-ea,0)./(gw+ref);
        signal = (x(9,:)./p(26,:)-p(29,:)+ ...
            p(31,:).*clip(x(7,:)./p(9,:))+p(32,:).*(gw./ref-1)) ...
            ./p(30,:);
        c = p(28,:)+(1-p(28,:)).*sigmoid(signal);
        bypass = p(33,:).*c.*qa;
        [dz9,u9] = cascade(x(15:14+m,:),qa-bypass,eta9);
        [dz1,u1] = cascade(x(15+m:ns,:),bypass,eta1);
        exchange = p(36,:).*(r./p(35,:)).^3.5;
        ir = u9+.9*(p(37,:)+exchange);
        if real(r(1)) <= 0
            ir = upper(ir,0);
        end
        qb = upper(u1+.1*(p(37,:)+exchange),0);
        modulation = log(upper(c,1e-12))+thermal+ ...
            p(41,:).*log(upper(r,1e-3)./p(35,:));
        qr = .25*r.*(r./p(35,:)).^4.*exp(6*tanh(modulation/6));
        preloss = qr+qb+p(14,:).*fn+qg;
        reported = excess(preloss,p(42,:));
        dy = [upstream;dr-qg-eg;ir-qr;dz9;dz1;reported; ...
            eh+ez+ea+eg+a.'*sub;ir+qb-u9-u1;effective;dr; ...
            preloss-reported;preloss];
    end

    function [dz,out] = cascade(store,in,eta)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CASCADE Evaluate routing-cascade derivatives and outlet release.
    %
    % SYNOPSIS:
    % cascade(...)
    %
    % INPUT ARGUMENTS:
    % store : Current water store (mm).
    % in : Input flux entering the first cascade stage.
    % eta : Routing-stage exit fractions (-).
    %
    % OUTPUT ARGUMENTS:
    % dz : Computed numerical result; tangent lanes are retained.
    % out : Model flux array or numerical result in the caller contract.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        release = store.*lambda;
        through = (1-eta).*release;
        dz = [in;through(1:end-1,:)]-release;
        out = eta.'*release;
    end

    function [q,e,d] = pareto(store,rain,pet,cap,shape,perc,net)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PARETO Advance the Pareto production-store water balance.
    %
    % SYNOPSIS:
    % pareto(...)
    %
    % INPUT ARGUMENTS:
    % store : Current water store (mm).
    % rain : Liquid input depth/rate on this operator time base.
    % pet : Potential evaporation depth/rate on this operator time base.
    % cap : Store capacity (mm).
    % shape : Dimensionless storage/infiltration shape parameter.
    % perc : Percolation coefficient in the model parameterization.
    % net : Whether to remove interception before the production update.
    %
    % OUTPUT ARGUMENTS:
    % q : Interval discharge depths (mm).
    % e : Computed numerical result; tangent lanes are retained.
    % d : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        wet = clip(store./cap);
        deficit = 1-wet;
        nu = shape./(1+shape);
        interception = zeros(1,lanes);
        if net
            interception = lower(rain,pet);
        end
        if real(deficit(1)) < trans
            retain = deficit/trans.*exp(nu*log(trans));
        else
            retain = deficit.^nu;
        end
        q = (rain-interception).*(1-retain);
        e = interception+(pet-interception).*(1-deficit.^2);
        d = .25*upper(store,0).*safe_log1p(perc.*wet.^4);
    end

    function y = integrate(y,dt,rhs,pilot)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INTEGRATE Advance the augmented state with the configured integrator.
    %
    % SYNOPSIS:
    % integrate(...)
    %
    % INPUT ARGUMENTS:
    % y : Numerical state/operand; augmented with sensitivities when
    % enabled.
    % dt : Integration or observation interval duration (days).
    % rhs : Callback evaluating state derivatives on the model time base.
    % pilot : Whether to run the reference-recharge initialization pass.
    %
    % OUTPUT ARGUMENTS:
    % y : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        shape = size(y);
        y = y(:);
        if pilot
            k = pilot_k;
            pilot_k = [];
        else
            k = full_k;
            full_k = [];
        end
        k(:,1) = reshape(rhs(reshape(y,shape)),[],1);
        t = 0;
        h = min(dt,.001);
        rejected = false;
        trials = 0;
        while t < dt
            trials = trials+1;
            if trials > 100000 || h < 10*eps(t)
                error('cartplus:step','Continuous CaRT+ step limit.');
            end
            stop = min(dt,t+h);
            h = stop-t;
            for s = 2:12
                stage = y+h*(k(:,1:s-1)*rkA(s,1:s-1).');
                k(:,s) = reshape(rhs(reshape(stage,shape)),[],1);
            end
            next = y+h*(k(:,1:12)*rkB.');
            k(:,13) = reshape(rhs(reshape(next,shape)),[],1);
            old = augmented(y,shape(1));
            proposed = augmented(next,shape(1));
            scale = options.atol+options.rtol*max(abs(old),abs(proposed));
            err5 = augmented(k*rkE5.',shape(1))./scale;
            err3 = augmented(k*rkE3.',shape(1))./scale;
            n5 = sum(err5.^2);
            n3 = sum(err3.^2);
            err = 0;
            if n5 ~= 0 || n3 ~= 0
                err = h*n5/sqrt((n5+.01*n3)*numel(scale));
            end
            factor = .2;
            if isfinite(err)
                factor = min(10,max(.2,.9*err^(-1/8)));
            end
            if err < 1
                y = next;
                k(:,1) = k(:,13);
                t = stop;
                if rejected
                    factor = min(1,factor);
                end
                rejected = false;
            else
                rejected = true;
            end
            h = h*factor;
        end
        y = reshape(y,shape);
        if pilot
            pilot_k = k;
        else
            full_k = k;
        end
    end

    function values = augmented(x,count)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUGMENTED Pack the primal state and its complex-step derivative lanes.
    %
    % SYNOPSIS:
    % augmented(...)
    %
    % INPUT ARGUMENTS:
    % x : Numerical operand; scalar, vector, or forward-mode scalar.
    % count : Number of intervals or state entries in this helper.
    %
    % OUTPUT ARGUMENTS:
    % values : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        values = real(x(1:count));
        if tangents
            values = [values;imag(x(count+1:end))/epsilon];
        end
    end
end

function [q,e,f] = zone(store,rain,pet,cap,days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZONE Evaluate near-stream runoff, evaporation, and drainage.
%
% SYNOPSIS:
% zone(...)
%
% INPUT ARGUMENTS:
% store : Current water store (mm).
% rain : Liquid input depth/rate on this operator time base.
% pet : Potential evaporation depth/rate on this operator time base.
% cap : Store capacity (mm).
% days : Residence time (days).
%
% OUTPUT ARGUMENTS:
% q : Interval discharge depths (mm).
% e : Computed numerical result; tangent lanes are retained.
% f : Computed numerical result; tangent lanes are retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wet = clip(store./cap);
    i = lower(rain,pet);
    q = (rain-i).*wet.^2;
    e = i+(pet-i).*wet.*(2-wet);
    f = upper(store,0)./days;
end

function eta = exit_fractions(m,symmetric)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXIT_FRACTIONS Construct the routing-cascade exit fractions.
%
% SYNOPSIS:
% exit_fractions(...)
%
% INPUT ARGUMENTS:
% m : Positive routing refinement: stages per X4 time base.
% symmetric : Whether to use the symmetric two-time-base hydrograph.
%
% OUTPUT ARGUMENTS:
% eta : Computed numerical result; tangent lanes are retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if symmetric
        t = (0:2*m)'/m;
        h = .5*min(t,1).^2.5;
        h(t>1) = 1-.5*(2-t(t>1)).^2.5;
    else
        h = ((0:m)'/m).^2.5;
    end
    w = diff(h);
    eta = w./flipud(cumsum(flipud(w)));
end

function q = excess(water,cap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCESS Evaluate smooth infiltration-excess runoff.
%
% SYNOPSIS:
% excess(...)
%
% INPUT ARGUMENTS:
% water : Available water or store depth (mm).
% cap : Store capacity (mm).
%
% OUTPUT ARGUMENTS:
% q : Interval discharge depths (mm).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = water./cap;
    if real(x(1)) < .001
        q = cap.*(x.^2/2-x.^3/6+x.^4/24-x.^5/120);
    else
        q = water+cap.*safe_expm1(-x);
    end
end

function y = safe_log1p(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAFE_LOG1P Evaluate log(1+x), including complex-step perturbations.
%
% SYNOPSIS:
% safe_log1p(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = real(x);
    y = log1p(r);
    if ~isreal(x)
        y = y+1i*imag(x)./(1+r);
    end
end

function y = safe_expm1(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAFE_EXPM1 Evaluate exp(x)-1, including complex-step perturbations.
%
% SYNOPSIS:
% safe_expm1(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = real(x);
    y = expm1(r);
    if ~isreal(x)
        y = y+1i*exp(r).*imag(x);
    end
end

function y = clip(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLIP Bound the numerical operand to its admissible interval.
%
% SYNOPSIS:
% clip(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = lower(upper(x,0),1);
end

function y = sigmoid(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGMOID Evaluate the dimensionless logistic response.
%
% SYNOPSIS:
% sigmoid(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    positive = real(x(:,1)) >= 0;
    y = zeros(size(x),'like',x);
    y(positive,:) = 1./(1+exp(-x(positive,:)));
    e = exp(x(~positive,:));
    y(~positive,:) = e./(1+e);
end

function y = lower(x,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOWER Apply a lower numerical bound while retaining derivatives.
%
% SYNOPSIS:
% lower(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
% b : Second numerical operand or rate coefficient.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = x;
    rows = real(b(:,1)) < real(x(:,1));
    if isscalar(b)
        y(rows,:) = b;
    else
        y(rows,:) = b(rows,:);
    end
end

function y = upper(x,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPPER Apply an upper numerical bound while retaining derivatives.
%
% SYNOPSIS:
% upper(...)
%
% INPUT ARGUMENTS:
% x : Numerical operand; scalar, vector, or forward-mode scalar.
% b : Second numerical operand or rate coefficient.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = x;
    rows = real(x(:,1)) < real(b(:,1));
    if isscalar(b)
        y(rows,:) = b;
    else
        y(rows,:) = b(rows,:);
    end
end
