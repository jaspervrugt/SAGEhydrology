function varargout = user_model(par,mdl,data,ode,check,request)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER_MODEL: Compiled modular implementation of user-defined basin model
% SYNOPSIS: varargout = user_model(par,mdl,data,ode,check,request)
%  par          dx1 vector of parameter values
%  mdl          structure with model state/parameter info
%   .mcode       scalar numerical solution of watershed model
%                 1 MATLAB adaptive RK2; 2 MATLAB ode45;
%                 3 fixed-step Euler; 4 compiled C++ adaptive RK2
%   .y0          mx1 vector of initial states
%   .L           fixed routing-memory length
%                 20 daily; 100 hourly and 15-minute
%   .pspace      0:hydrologic, 1:unit cube, 2:unconstrained parameters
%   .th_min      dx1 vector of lower parameter values [= in pspace]
%   .th_max      dx1 vector of upper parameter values [= in pspace]
%   .par_names   1xd cell parameter names
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
%  check        retained for compatibility with built-in model interface
%   0            accepted; currently unused by compiled interface
%   1            accepted; currently unused by compiled interface
%  q_n          nx1 vector of simulated discharge
%  J            nxd discharge Jacobian in space selected by mdl.pspace
%                 0 dq_n/dth'
%                 1 dq_n/dnth'
%                 2 dq_n/dvth'
%  Jth          dx1 vector with diagonal of parameter transformation
%                 0 ones(d,1)
%                 1 dth/dnth
%                 2 dth/dvth
%  Z            complete augmented state trajectory when ode.mem = 1
%                final augmented state when ode.mem = 0
%  request      OPTIONAL named-output request for Q, SWE, SM, or states
%  varargout    requested result structure, or q_n, J, Jth, Z
%  notes        solver, parameter mapping, state preparation, sensitivity
%               propagation and output handling are implemented by
%               crr_user_model_mex.cpp, user_model_prepare.cpp and
%               user_model.cpp. The gateway requests only outputs needed
%               by the caller and avoids a complete Jacobian allocation
%               for discharge-only calls when ode.mem = 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, May 2026 / updated last Sept. 2026        %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    narginchk(4,6);
    nargoutchk(0,4);

    if nargin < 5
        check = 0;
    end
    structured = nargin >= 6 && ~isempty(request);
    if structured
        request = crr_request(request);
        if ~isempty(request.states)
            ode.mem = 1;
        else
            ode.mem = 0;
        end
    end

    if ~isstruct(mdl) ...
            || ~isfield(mdl,'mcode')
        error('user_model:InvalidModelStructure', ...
            ['mdl must be a structure ' ...
            'containing field mcode.']);
    end

    if ~ismember(mdl.mcode,1:4)
        error('user_model:InvalidSolver', ...
            'mcode must be 1, 2, 3, or 4.');
    end

    if ~isstruct(data)
        error('user_model:InvalidDataStructure', ...
            'data must be a structure.');
    end

    if ~isstruct(ode)
        error('user_model:InvalidODEStructure', ...
            'ode must be a structure.');
    end

    if ~isnumeric(par) ...
        || ~isreal(par) ...
        || isempty(par) ...
        || ~isvector(par)
        error('user_model:InvalidParameters', ...
            ['par must be a nonempty ' ...
            'real numeric vector.']);
    end

    if ~isscalar(check) ...
            || ~(islogical(check) ...
            || isnumeric(check))
        error('user_model:InvalidCheck', ...
            ['check must be a scalar ' ...
             'logical or numeric value.']);
    end

    mdl=apply_model_configuration(mdl,ode);
    if mdl.mcode~=4
        [q_n,J,Jth,Z,swe,Jswe,sm,Jsm,fail]=user_model_matlab(par,mdl,data,ode);
        if structured
            if ~isempty(request.states)
                out=crr_result(q_n,J,Jth,Z,mdl,request);
            else
                out=struct('failed',fail,'obs',struct(),'jac',struct());
                for observable=request.obs(:)'
                    switch observable,case "Q",out.obs.Q=q_n;case "SWE",out.obs.SWE=swe;case "SM",out.obs.SM=sm;end
                end
                for observable=request.jac(:)'
                    switch observable,case "Q",out.jac.Q=J;case "SWE",out.jac.SWE=Jswe;case "SM",out.jac.SM=Jsm;end
                end
            end
            varargout={out};
        else
            output={q_n,J,Jth,Z};varargout=output(1:nargout);
        end
        return
    end

    if exist('crr_user_model','file') ~= 3
        if isdeployed
            error('user_model:MissingDeployedMEX', ...
                ['The deployed SAGE application cannot compile a user ' ...
                 'model. Install a precompiled crr_user_model binary for ' ...
                 'this operating system and MATLAB Runtime release in ' ...
                 'the user_model folder.']);
        end
        error('user_model:MissingMEX', ...
            ['crr_user_model is not compiled. Run user_model_compile ' ...
             'from the user_model folder with a configured C++ compiler.']);
    end

    if structured
        if isempty(request.states)
            ode.want_q = any(request.obs=="Q");
            ode.want_jq = any(request.jac=="Q");
            ode.want_swe = any(request.obs=="SWE");
            ode.want_jswe = any(request.jac=="SWE");
            ode.want_sm = any(request.obs=="SM");
            ode.want_jsm = any(request.jac=="SM");
            [q_n,J,~,~,fail,swe,Jswe,sm,Jsm] = ...
                crr_user_model(par,mdl,data,ode,check);
            out = struct('failed',logical(fail));
            if ~isempty(request.obs), out.obs=struct(); end
            if any(request.obs=="Q"), out.obs.Q=q_n; end
            if any(request.obs=="SWE"), out.obs.SWE=swe; end
            if any(request.obs=="SM"), out.obs.SM=sm; end
            if ~isempty(request.jac), out.jac=struct(); end
            if any(request.jac=="Q"), out.jac.Q=J; end
            if any(request.jac=="SWE"), out.jac.SWE=Jswe; end
            if any(request.jac=="SM"), out.jac.SM=Jsm; end
        else
            [q_n,J,Jth,Z] = ...
                crr_user_model(par,mdl,data,ode,check);
            out = crr_result(q_n,J,Jth,Z,mdl,request);
        end
        varargout = {out};
        return
    end

    [q_n,J,Jth,Z] = ...
        crr_user_model(par,mdl,data,ode,check);

    output = {q_n,J,Jth,Z};
    varargout = output(1:nargout);

end

