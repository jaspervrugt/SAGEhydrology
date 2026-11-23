function [q,J] = cart_ode_discharge(data,parameters,model,config)
%CART_ODE_DISCHARGE Run the native CaRT or CaRT+ ODE kernel from MATLAB.
    arguments
        data (1,1) struct
        parameters double
        model (1,1) double
        config.cart_solver (1,1) string ...
            {mustBeMember(config.cart_solver,["dop853","rk2","rk23", ...
            "rk45","heun_fixed","euler"])} = "dop853"
        config.cart_refinement (1,1) double ...
            {mustBePositive,mustBeInteger} = 16
        config.cart_rtol double {mustBePositive,mustBeFinite} = []
        config.cart_atol double {mustBePositive,mustBeFinite} = []
        config.cart_max_step double {mustBePositive,mustBeFinite} = []
        config.cart_min_step (1,1) double ...
            {mustBeNonnegative,mustBeFinite} = 0
    end
    persistent gateways buildFolder
    if isempty(gateways), gateways = cell(1,2); end
    k = 1 + (model == 18);
    name = ["crr_cart_numerics_v1","crr_cartplus_numerics_v1"];
    if isempty(gateways{k})
        if exist(name(k),'file') ~= 3
            compiler = mex.getCompilerConfigurations('C++','Selected');
            if isempty(compiler)
                error('cart:ode:compiler', ...
                    ['Select a MATLAB C++ compiler with mex -setup C++ ' ...
                    'to build the CaRT ODE gateway.']);
            end
            if isempty(buildFolder)
                buildFolder = tempname;
                mkdir(buildFolder);
                addpath(buildFolder);
            end
            if contains(compiler.Name,'Microsoft')
                flags = {'COMPFLAGS=$COMPFLAGS /std:c++17'};
            else
                flags = {'CXXFLAGS=$CXXFLAGS -std=c++17'};
            end
            if k == 2, flags{end+1} = '-DMATLAB_CARTPLUS'; end
            source = fullfile(fileparts(mfilename('fullpath')), ...
                'cart_ode_mex.cpp');
            mex('-O','-R2017b',flags{:},'-outdir',buildFolder, ...
                '-output',char(name(k)),source);
        end
        gateways{k} = str2func(char(name(k)));
    end
    codes = [0 1 2 9 10 11];
    methods = ["dop853","rk2","rk23","rk45","heun_fixed","euler"];
    method = codes(find(config.cart_solver == methods,1));
    dt = repmat(1/data.resolution,size(data.forcing,1),1);
    floor = 0;
    if isfield(data,'floor'), floor = data.floor; end
    [rtol,atol] = cart_ode_tolerances(config.cart_solver);
    if ~isempty(config.cart_rtol), rtol = config.cart_rtol; end
    if ~isempty(config.cart_atol), atol = config.cart_atol; end
    if ~isscalar(rtol) || ~isscalar(atol)
        error('cart:ode:tolerance', ...
            'ODE tolerances must be positive finite scalars.');
    end
    maximum = config.cart_max_step;
    if isempty(maximum)
        maximum = 1e100;
        if method >= 9, maximum = 1/96; end
    end
    if ~isscalar(maximum) || config.cart_min_step > maximum
        error('cart:ode:step','Invalid maximum/minimum step in days.');
    end
    options = [rtol,atol,floor,.02,.01,maximum,config.cart_min_step];
    run = gateways{k};
    if nargout > 1
        [q,J] = run(data.forcing,data.elevation,data.area,parameters, ...
            dt,config.cart_refinement,options,method);
    else
        q = run(data.forcing,data.elevation,data.area,parameters, ...
            dt,config.cart_refinement,options,method);
    end
end
