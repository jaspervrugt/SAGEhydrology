% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cart/test_cart_jacobian.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function test_cart_jacobian()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_CART_JACOBIAN Verify cart jacobian.
%
% SYNOPSIS:
% test_cart_jacobian(...)
%
% INPUT ARGUMENTS:
% None.
%
% OUTPUT ARGUMENTS:
% No return value; performs the documented build or validation.
%
% NOTES:
% TEST_CART_JACOBIAN Check all CaRT cases against centered differences.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rng(52,'twister');
    n = 36;
    f = [1+9*rand(n,1), -6+21*rand(n,1), ...
        .1+4*rand(n,1), (1:n)'<=24];
    p = [log(225);1.26;log(33);log(78);-.35;0;0; ...
        log(3);1;log(6.5);log(.001)];
    d.elevation = (0:4)'*200;
    d.area = ones(5,1)/5;
    step = 1e-5;
    for resolution = [1 24 96]
        d.resolution = resolution;
        d.forcing = f;
        d.forcing(:,[1 3]) = f(:,[1 3])/resolution;
        for code = 1:4
            mdl.mcode = code;
            [~,J,initial] = cart(p,mdl,d);
            err = 0;
            for j = 1:numel(p)
                delta = zeros(size(p));
                delta(j) = step;
                fd = (cart(p+delta,mdl,d,initial) ...
                    - cart(p-delta,mdl,d,initial))/(2*step);
                err = max(err,max(abs(fd-J(:,j))));
            end
            fprintf('CaRT %d steps/day, case %d: %.3g\n', ...
                resolution,code,err);
            assert(err < 1e-7,'CaRT finite-difference mismatch.');
        end
    end
end
