function [rtol,atol] = cart_ode_tolerances(method)
%CART_ODE_TOLERANCES Default error tolerances for the CaRT ODE solvers.
    if strcmp(method,'rk2')
        rtol = 1e-3;
        atol = 1e-3;
    else
        rtol = 1e-8;
        atol = 1e-10;
    end
end
