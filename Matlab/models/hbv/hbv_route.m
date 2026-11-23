function [q,J] = hbv_route(qraw,Jraw,width,length)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HBV_ROUTE Shared native MAXBAS routing for discharge and sensitivities.
%
% qraw is a discharge vector; Jraw has the 12 ODE sensitivity columns (or
% is empty). The returned Jacobian has 13 columns, including dQ/dMAXBAS.
% Inputs already contain the selected post-spin-up interval. Routing starts
% with an empty filter state. The model equations are not evaluated here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026 Jasper A. Vrugt
% University of California, Irvine | jasper@uci.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist('hbv_route_mex','file') ~= 3
        compile_hbv_routing();
    end
    [q,J] = hbv_route_mex(qraw,Jraw,width,double(length));
end
