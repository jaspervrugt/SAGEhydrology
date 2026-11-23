% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cartplus/cartplus_zone_solver.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function y = cartplus_zone_solver(input,mcode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CARTPLUS_ZONE_SOLVER Integrate near-stream storage and its forward
% derivatives.
%
% SYNOPSIS:
% cartplus_zone_solver(...)
%
% INPUT ARGUMENTS:
% input : Forcing/storage input for the current operator.
% mcode : Backend/solver case selected by the model wrapper.
%
% OUTPUT ARGUMENTS:
% y : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% CARTPLUS_ZONE_SOLVER Integrate normalized zone storage over one interval.
% input = [initial storage, net forcing, recession, days, derivatives]. When
% requested, columns 2–4 carry tangents to the first three inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    initial = input(1);
    net = input(2);
    recession = input(3);
    dt = input(4);
    derivative = logical(input(5));
    state = zeros(2,1+3*derivative);
    state(1,1) = initial;
    if derivative
        state(1,2) = 1;
    end

    rhs = @(~,y) zone_rhs(y,net,recession);
    switch mcode
        case 2
            options = odeset('RelTol',1e-10,'AbsTol',1e-12);
            [~,trajectory] = ode45(rhs,[0 1],state(:),options);
            y = trajectory(end,:).';
        case {1,3}
            steps = max([1,ceil(16*dt),ceil(4*(abs(net)+recession))]);
            h = 1/steps;
            y = state(:);
            for k = 1:steps
                rate = rhs(0,y);
                if mcode == 1
                    rate = rhs(0,y+.5*h*rate);
                end
                y = y+h*rate;
            end
    end
end

function rate = zone_rhs(y,net,recession)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZONE_RHS Evaluate near-stream storage and tangent derivatives.
%
% SYNOPSIS:
% zone_rhs(...)
%
% INPUT ARGUMENTS:
% y : Numerical state/operand; augmented with sensitivities when enabled.
% net : Whether to remove interception before the production update.
% recession : Groundwater recession coefficient in the model time base.
%
% OUTPUT ARGUMENTS:
% rate : Computed numerical result; tangent lanes are retained.
%
% NOTES:
% ZONE_RHS Return storage, integrated storage and variational rates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    state = reshape(y,2,[]);
    s = state(1,1);
    if net >= 0
        shape = 1-s*s;
        slope = -2*s;
    else
        shape = s*(2-s);
        slope = 2-2*s;
    end
    rate = zeros(size(state));
    rate(1,:) = (net*slope-recession)*state(1,:);
    rate(1,1) = net*shape-recession*s;
    rate(2,:) = state(1,:);
    if size(state,2) > 1
        rate(1,3) = rate(1,3)+shape;
        rate(1,4) = rate(1,4)-s;
    end
    rate = rate(:);
end
