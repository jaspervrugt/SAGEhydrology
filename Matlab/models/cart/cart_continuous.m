% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cart/cart_continuous.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function [q,jac] = cart_continuous(f,z,a,p,initial,rtol,atol,withJac)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CART_CONTINUOUS Integrate fully coupled CaRT with fixed initial stores.
%
% SYNOPSIS:
% cart_continuous(...)
%
% INPUT ARGUMENTS:
% f : Forcing [P, T, PET, training] in the caller layout.
% z : Five elevation-band heights (m).
% a : Five nonnegative band fractions, summing to one.
% p : 11 CaRT calibration coordinates, including logarithmic entries.
% initial : Initial stores (mm), or [] to use the spin-up policy.
% rtol : Positive relative integration-error tolerance (-).
% atol : Positive absolute integration-error tolerance for augmented states.
% withJac : Whether to evaluate discharge sensitivities.
%
% OUTPUT ARGUMENTS:
% q : Interval discharge depths (mm).
% jac : Discharge Jacobian in the selected parameter directions.
%
% NOTES:
% CART_CONTINUOUS Daily-boundary DOP853 for the coupled CaRT equations.
% Parameters use calibration coordinates. Initial stores are detached.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = p(:);
    z = z(:);
    a = a(:);
    m = (numel(initial)-8)/3;
    ns = numel(initial);
    epsCS = 1e-25;
    if withJac
        p = p + 1i*epsCS*[zeros(11,1),eye(11)];
    end
    lanes = size(p,2);
    [rkA,rkB,e3,e5] = cartplus_dop853();
    C = exp(p(1,:));
    x3 = exp(p(3,:));
    tau = exp(p(4,:));
    kp = (4/9)^4*exp(p(6,:));
    cf = exp(p(8,:));
    gamma = exp(p(10,:));
    ks = exp(p(11,:));
    lambda = m./p(2,:);
    offset = (z-a.'*z)/1000;
    t9 = (0:m)'/m;
    t1 = (0:2*m)'/m;
    H1 = .5*min(t1,1).^2.5;
    H1(t1>1) = 1-.5*(2-t1(t1>1)).^2.5;
    w9 = diff(t9.^2.5);
    w1 = diff(H1);
    eta9 = w9./flipud(cumsum(flipud(w9)));
    eta1 = w1./flipud(cumsum(flipud(w1)));
    train = f(:,4) ~= 0;
    n = size(f,1);
    gbar = zeros(1,lanes,'like',p);
    run(true);
    gbar = greater(gbar/sum(train),1e-6);
    hydrograph = run(false);
    q = real(hydrograph(:,1));
    jac = imag(hydrograph(:,2:end))/epsCS;

    function hydrograph = run(pilot)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RUN Advance the pilot or full model over the forcing sequence.
    %
    % SYNOPSIS:
    % run(...)
    %
    % INPUT ARGUMENTS:
    % pilot : Whether to run the reference-recharge initialization pass.
    %
    % OUTPUT ARGUMENTS:
    % hydrograph : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        count = ns+5;
        if pilot
            count = 7;
        end
        y = zeros(count,lanes,'like',p);
        if pilot
            y(1:6,:) = repmat(reshape(initial(1:6),[],1),1,lanes);
        else
            y(1:ns,:) = repmat(initial(:),1,lanes);
        end
        hydrograph = zeros(n,lanes,'like',p);
        for i = 1:n
            previous = y;
            rain = f(i,1);
            pet = f(i,3);
            warm = f(i,2)-offset.*gamma-p(7,:);
            cold = real(warm(:,1)) < 0;
            snowfall = cold.*p(9,:)*rain;
            potential = cf.*greater(warm,0);
            active = real(y(1:5,1)) > 0 | ...
                real(snowfall(:,1)-ks(1)*pet-potential(:,1)) > 0;
            t = 0;
            while t < 1
                [snow,net,demand,atmosphere] = forcing(active);
                hits = inf(5,1);
                draining = active & real(snow(:,1)) < 0;
                hits(draining) = max(0,real(y(draining,1)))./ ...
                    (-real(snow(draining,1)));
                remaining = 1-t;
                dt = min(remaining,min(hits));
                rhs = @(x) rates(x,pilot,snow,net,demand,atmosphere, ...
                    rain,cold);
                if dt > 0
                    y = integrate(y,dt,rhs);
                end
                t = t+dt;
                for hit = find(hits <= dt).'
                    if withJac
                        [sn,pn,en,atm] = forcing(active);
                        before = rates(y,pilot,sn,pn,en,atm,rain,cold);
                        event = -imag(y(hit,2:end))/real(before(hit,1));
                    end
                    active(hit) = false;
                    y(hit,:) = 0;
                    if withJac
                        [sn,pn,en,atm] = forcing(active);
                        after = rates(y,pilot,sn,pn,en,atm,rain,cold);
                        y(:,2:end) = y(:,2:end)+1i*(real(before(:, ...
                            1)-after(:,1))*event);
                        y(hit,:) = 0;
                    end
                end
                if 1-t < 1e-14
                    t = 1;
                end
            end
            if pilot && train(i)
                gbar = gbar+y(7,:)-previous(7,:);
            elseif ~pilot
                hydrograph(i,:) = y(ns+1,:)-previous(ns+1,:);
            end
        end

        function [snow,net,demand,atm] = forcing(active)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FORCING Partition band precipitation and snowmelt into liquid
        % supply.
        %
        % SYNOPSIS:
        % forcing(...)
        %
        % INPUT ARGUMENTS:
        % active : Flags for snow bands with an active snow store.
        %
        % OUTPUT ARGUMENTS:
        % snow : Computed numerical result; tangent lanes are retained.
        % net : Computed numerical result; tangent lanes are retained.
        % demand : Computed numerical result; tangent lanes are retained.
        % atm : Computed numerical result; tangent lanes are retained.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sub = repmat(ks*pet,5,1);
            melt = potential;
            sub(~active,:) = lesser(sub(~active,:),snowfall(~active,:));
            melt(~active,:) = lesser(melt(~active,:), ...
                greater(snowfall(~active,:)-sub(~active,:),0));
            snow = snowfall-sub-melt;
            liquid = a.'*(~cold*rain+melt);
            e0 = lesser(liquid,pet);
            net = liquid-e0;
            demand = pet-e0;
            atm = e0+a.'*sub;
        end
    end

    function dy = rates(y,pilot,snow,net,demand,atmosphere,rain,cold)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RATES Evaluate state and cumulative-flux derivatives.
    %
    % SYNOPSIS:
    % rates(...)
    %
    % INPUT ARGUMENTS:
    % y : Numerical state/operand; augmented with sensitivities when
    % enabled.
    % pilot : Whether to run the reference-recharge initialization pass.
    % snow : Snow-store depths for the five elevation bands (mm).
    % net : Whether to remove interception before the production update.
    % demand : Evaporation demand on this operator time base.
    % atmosphere : Evaporation/sublimation rate leaving the catchment.
    % rain : Liquid input depth/rate on this operator time base.
    % cold : Boolean flags for freezing band temperatures.
    %
    % OUTPUT ARGUMENTS:
    % dy : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dy = zeros(size(y),'like',p);
        dy(1:5,:) = snow;
        soil = y(6,:);
        s = soil./C;
        s2 = s.*s;
        infiltration = (1-s2).*net;
        evaporation = s.*(2-s).*demand;
        recharge = .25*kp.*soil.*s2.*s2;
        dy(6,:) = infiltration-evaporation-recharge;
        if pilot
            dy(7,:) = recharge;
            return
        end
        direct = net-infiltration;
        [dy(9:8+m,:),u9] = cascade(y(9:8+m,:),.9*direct,eta9);
        [dy(9+m:ns,:),u1] = cascade(y(9+m:ns,:),.1*direct,eta1);
        R = greater(y(7,:),0);
        W = greater(y(8,:),0);
        relative = R./x3;
        base = max(real(relative(1)),0);
        lossPower = base^3.5 + 1i*3.5*base^2.5*imag(relative);
        if ~withJac
            lossPower = real(lossPower);
        end
        loss = -p(5,:).*lossPower;
        relative2 = relative.*relative;
        qr = .25*R.*relative2.*relative2;
        qb = greater(u1-.1*loss,0);
        denominator = gbar.*tau.*tau;
        qg = (W.*W)./denominator;
        dy(7,:) = u9-qr-.9*loss;
        dy(8,:) = recharge-qg;
        dy(ns+1,:) = qr+qb+qg;
        dy(ns+2,:) = atmosphere+evaporation;
        dy(ns+3,:) = .9*loss+lesser(u1,.1*loss);
        dy(ns+4,:) = rain+(p(9,:)-1)*rain*sum(a(cold));
        dy(ns+5,:) = recharge;
    end

    function [dy,out] = cascade(y,in,eta)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CASCADE Evaluate routing-cascade derivatives and outlet release.
    %
    % SYNOPSIS:
    % cascade(...)
    %
    % INPUT ARGUMENTS:
    % y : Numerical state/operand; augmented with sensitivities when
    % enabled.
    % in : Input flux entering the first cascade stage.
    % eta : Routing-stage exit fractions (-).
    %
    % OUTPUT ARGUMENTS:
    % dy : Computed numerical result; tangent lanes are retained.
    % out : Model flux array or numerical result in the caller contract.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        release = lambda.*y;
        dy = [in;(1-eta(1:end-1)).*release(1:end-1,:)]-release;
        out = eta.'*release;
    end

    function y = integrate(y,dt,rhs)
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
    %
    % OUTPUT ARGUMENTS:
    % y : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        shape = size(y);
        y = y(:);
        k = zeros(numel(y),13,'like',p);
        k(:,1) = reshape(rhs(reshape(y,shape)),[],1);
        t = 0;
        h = min(dt,.001);
        rejected = false;
        while t < dt
            stop = min(dt,t+h);
            h = stop-t;
            if h <= 0
                error('cart:step', ...
                    'DOP853 step underflow at %g of segment %g.',t,dt);
            end
            for stage = 2:12
                next = y+h*weighted(k(:,1:stage-1),rkA(stage,1:stage-1).');
                k(:,stage) = reshape(rhs(reshape(next,shape)),[],1);
            end
            next = y+h*weighted(k(:,1:12),rkB.');
            k(:,13) = reshape(rhs(reshape(next,shape)),[],1);
            scale = atol+rtol*max(abs(augmented(y,shape(1))), ...
                abs(augmented(next,shape(1))));
            n5 = sum((augmented(weighted(k,e5.'),shape(1))./scale).^2);
            n3 = sum((augmented(weighted(k,e3.'),shape(1))./scale).^2);
            err = 0;
            if n5 ~= 0 || n3 ~= 0
                err = h*n5/sqrt((n5+.01*n3)*numel(scale));
            end
            if ~isfinite(err) && h < 1e-12
                [badRows,badLanes] = find(~isfinite( ...
                    reshape(k(:,1),shape)));
                error('cart:step', ...
                    'Nonfinite derivative at %g, rows %s, lanes %s.', ...
                    t,mat2str(badRows.'),mat2str(badLanes.'));
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
    end

    function product = weighted(values,weights)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WEIGHTED Sum elevation-band values using the supplied area weights.
    %
    % SYNOPSIS:
    % weighted(...)
    %
    % INPUT ARGUMENTS:
    % values : Numerical values or vector of state entries.
    % weights : Dimensionless quadrature or band weights.
    %
    % OUTPUT ARGUMENTS:
    % product : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        product = real(values)*weights;
        if withJac
            product = complex(product,imag(values)*weights);
        end
    end

    function values = augmented(y,count)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUGMENTED Pack the primal state and its complex-step derivative lanes.
    %
    % SYNOPSIS:
    % augmented(...)
    %
    % INPUT ARGUMENTS:
    % y : Numerical state/operand; augmented with sensitivities when
    % enabled.
    % count : Number of intervals or state entries in this helper.
    %
    % OUTPUT ARGUMENTS:
    % values : Computed numerical result; tangent lanes are retained.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        values = real(y(1:count));
        if withJac
            values = [values;imag(y(count+1:end))/epsCS];
        end
    end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = choose(a,b,@le);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = choose(a,b,@ge);
end

function y = choose(a,b,compare)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE Select a numerical branch using the supplied comparison.
%
% SYNOPSIS:
% choose(...)
%
% INPUT ARGUMENTS:
% a : First numerical operand or elevation-band fractions.
% b : Second numerical operand or rate coefficient.
% compare : Branch-comparison callback.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rows = max(size(a,1),size(b,1));
    cols = max(size(a,2),size(b,2));
    if isempty(a) || isempty(b)
        y = zeros(0,cols,'like',a);
        return
    end
    a = a+zeros(rows,cols);
    b = b+zeros(rows,cols);
    take = compare(real(a(:,1)),real(b(:,1)));
    y = a;
    y(~take,:) = b(~take,:);
end
