% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cartplus/cartplus.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function [out,J] = cartplus(par,mdl,data,ode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CARTPLUS Evaluate the CaRT+ rainfall-runoff model and discharge Jacobian.
%
% SYNOPSIS:
% cartplus(...)
%
% INPUT ARGUMENTS:
% par : 42 physical CaRT+ parameters in native order.
% mdl : Model configuration, including solver/backend case.
% data : Basin forcing, elevation bands, and model settings.
% ode : Numerical integration and derivative settings.
% data.forcing : N-by-4 [P,T,PET,training]; P/PET mm per interval, T deg C.
% data.elevation : Five band elevations (m).
% data.area : Five nonnegative band-area fractions summing to one.
% data.resolution : Intervals per day (1, 24, or 96).
%
% OUTPUT ARGUMENTS:
% out : Model flux array or numerical result in the caller contract.
% J : Discharge Jacobian with respect to the wrapper coordinates.
%
% NOTES:
% CARTPLUS Run the CaRT+ rainfall-runoff model. par             42-by-1
% vector of physical CaRT+ parameters. mdl.mcode       Zone ODE solver: 1
% RK2, 2 ode45, 3 Euler, 4 C++. Only case 4 loads the native C++ MEX.
% data.forcing    N-by-4 matrix [P, T, PET, training mask]. P and PET are
% depths per interval (mm), T is in deg C. data.elevation  Five band
% elevations (m). data.area       Five band area fractions. data.resolution
% Number of steps per day (1, 24, or 96). data.floor      Optional daily
% infiltration capacity floor. ode.jacobian    'analytic' (default),
% 'numerical' (centered difference), or 'complex-step' (MATLAB cases only).
% out             [Q, ET, external flux, storage change] in mm/interval. J
% N-by-42 discharge Jacobian at fixed spin-up stores.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 4
        ode = struct();
    end
    if ~isfield(data,'floor')
        data.floor = 0;
    end
    if ~ismember(data.resolution,[1 24 96])
        error('cartplus:interval','Expected 1, 24, or 96 steps per day.');
    end
    if size(data.forcing,2) ~= 4 || isempty(data.forcing) || ...
        numel(data.elevation) ~= 5 || numel(data.area) ~= 5 || ...
        numel(par) ~= 42
        error('cartplus:input','Invalid forcing, bands, or parameters.');
    end
    data.forcing = double(data.forcing);
    data.elevation = double(data.elevation);
    data.area = double(data.area);
    data.resolution = double(data.resolution);
    data.floor = double(data.floor);
    par = double(par);
    if ~isfield(ode,'jacobian')
        ode.jacobian = 'analytic';
    end

    native = @crr_cartplus_mex;
    if isfield(mdl,'config') && ...
            isfield(mdl.config,'cart_uniform_input') && ...
            mdl.config.cart_uniform_input
        if mdl.mcode ~= 4
            error('cartplus:scheme','Uniform input requires mcode 4.');
        end
        native = @crr_cartplus_uniform_f986af2;
    end
    switch mdl.mcode
        case {1,2,3}
            [out,states] = cartplus_standalone(data,par,mdl.mcode,[]);
            if nargout > 1
                run = @(p,s) cartplus_standalone(data,p,mdl.mcode,s);
                if strcmpi(ode.jacobian,'analytic')
                    J = analytic_jacobian(run,par,states,out);
                else
                    J = difference(run,par,states,out,ode.jacobian);
                end
            end
            return
        case 4
        otherwise
            error('cartplus:solver','CaRT+ requires mcode 1, 2, 3 or 4.');
    end
    args = {data.forcing,data.elevation,data.area,par, ...
        data.resolution,data.floor};
    if nargout < 2
        out = native(args{:});
    elseif strcmpi(ode.jacobian,'analytic')
        [out,J] = native(args{:});
    elseif strcmpi(ode.jacobian,'numerical')
        [out,~,states] = native(args{:});
        run = @(p,s) native(args{1:3},p,args{5:6},s);
        J = difference(run,par,states,out,ode.jacobian);
    else
        error('cartplus:jacobian','Unknown Jacobian method.');
    end
end

function J = analytic_jacobian(run,par,states,out)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTIC_JACOBIAN Propagate seeded derivatives with fixed spin-up stores.
%
% SYNOPSIS:
% analytic_jacobian(...)
%
% INPUT ARGUMENTS:
% run : Forward-model callback used for Jacobian evaluation.
% par : Model parameter vector in native model order.
% states : Hydrological state history or detached spin-up stores.
% out : Flux/state output; modified in place for mutating routines.
%
% OUTPUT ARGUMENTS:
% J : Discharge Jacobian with respect to the wrapper coordinates.
%
% NOTES:
% ANALYTIC_JACOBIAN Propagate seeded tangents at fixed spin-up.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    J = zeros(size(out,1),numel(par));
    for j = 1:numel(par)
        seed = zeros(size(par));
        seed(j) = 1;
        result = run(CartplusTangent(par,seed),states);
        J(:,j) = result.slope(:,1);
    end
    if any(~isfinite(J),'all')
        error('cartplus:jacobian','The analytic Jacobian is nonfinite.');
    end
end

function J = difference(run,par,states,out,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENCE Check sensitivities by complex-step or centered differences.
%
% SYNOPSIS:
% difference(...)
%
% INPUT ARGUMENTS:
% run : Forward-model callback used for Jacobian evaluation.
% par : Model parameter vector in native model order.
% states : Hydrological state history or detached spin-up stores.
% out : Flux/state output; modified in place for mutating routines.
% method : Selected derivative or solver method.
%
% OUTPUT ARGUMENTS:
% J : Discharge Jacobian with respect to the wrapper coordinates.
%
% NOTES:
% DIFFERENCE Complex-step or centered differences at fixed spin-up.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    J = zeros(size(out,1),numel(par));
    for j = 1:numel(par)
        trial = par;
        switch lower(method)
            case 'complex-step'
                step = 1e-24;
                trial(j) = par(j)+1i*step;
                output = run(trial,states);
                J(:,j) = imag(output(:,1))/step;
            case 'numerical'
                step = 1e-5*max(abs(par(j)),1e-3);
                trial(j) = par(j)+step;
                plus = run(trial,states);
                trial(j) = par(j)-step;
                minus = run(trial,states);
                J(:,j) = (plus(:,1)-minus(:,1))/(2*step);
            otherwise
                error('cartplus:jacobian','Unknown Jacobian method.');
        end
    end
end
