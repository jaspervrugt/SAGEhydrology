% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cartplus/test_cartplus_jacobian.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

function test_cartplus_jacobian(cases)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_CARTPLUS_JACOBIAN Verify cartplus jacobian.
%
% SYNOPSIS:
% test_cartplus_jacobian(...)
%
% INPUT ARGUMENTS:
% cases : Solver cases selected for derivative validation.
%
% OUTPUT ARGUMENTS:
% No return value; performs the documented build or validation.
%
% NOTES:
% TEST_CARTPLUS_JACOBIAN Compare derivatives with five-point differences.
% Exercises all four MATLAB solver cases at three forcing intervals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 1
        cases = 1:4;
    end
    par = [0 3 1 .001 .01 6.5 20 10 225 .039018442310623375 ...
        .7 1000 .3 .05 15 4 .2 1 .01 78 .05 .1 8 3 20 50 ...
        20 .75 .4 .15 .05 .05 .13 1.26 33 -.35 .001 0 2.4 ...
        1 .05 .003];
    data.elevation = [100 220 350 510 740];
    data.area = [.1 .2 .3 .25 .15];
    data.floor = 0;

    for r = [1 24 96]
        t = (0:9)';
        P = (2 + 1.5*sin(t) + 8*(sin(2*t) > .6)) / r;
        T = 4 - 14*cos(t/2);
        PET = (.7 + .5*cos(t)) / r;
        data.forcing = [P,T,PET,double(t >= 1 & t < 8)];
        data.resolution = r;
        for mcode = cases
            mdl.mcode = mcode;
            [out,analytic] = cartplus(par,mdl,data);
            [~,numerical] = cartplus(par,mdl,data, ...
                struct('jacobian','numerical'));
            if mcode < 4
                [~,states] = cartplus_standalone(data,par,mcode,[]);
                run = @(p) cartplus_standalone(data,p,mcode,states);
            else
                args = {data.forcing,data.elevation,data.area,par,r,0};
                [~,~,states] = crr_cartplus_mex(args{:});
                run = @(p) crr_cartplus_mex(args{1:3},p,args{5:6},states);
            end

            reference = zeros(size(analytic));
            for j = 1:numel(par)
                step = 3e-6 * max(abs(par(j)),1e-3);
                direction = zeros(size(par));
                direction(j) = step;
                plus2 = discharge(run,par+2*direction);
                plus1 = discharge(run,par+direction);
                minus1 = discharge(run,par-direction);
                minus2 = discharge(run,par-2*direction);
                reference(:,j) = ...
                    (-plus2+8*plus1-8*minus1+minus2)/(12*step);
            end

            tol = 2e-5 + .005 * abs(reference);
            assert(all(abs(analytic-reference) <= tol,'all'));
            assert(all(abs(numerical-reference) <= tol,'all'));
            balance = P + out(:,3) - out(:,1) - out(:,2) - out(:,4);
            assert(max(abs(balance)) < 1e-10);
            format = ['%d steps/day, case %d: ' ...
                'primary %.3g, centered %.3g\n'];
            fprintf(format,r,mcode, ...
                max(abs(analytic-reference),[],'all'), ...
                max(abs(numerical-reference),[],'all'));
        end
    end
end

function q = discharge(run,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCHARGE Extract interval discharge from a model callback.
%
% SYNOPSIS:
% discharge(...)
%
% INPUT ARGUMENTS:
% run : Forward-model callback used for Jacobian evaluation.
% par : Model parameter vector in native model order.
%
% OUTPUT ARGUMENTS:
% q : Interval discharge depths (mm).
%
% NOTES:
% DISCHARGE Run a parameter perturbation with the baseline stores fixed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    out = run(par);
    q = out(:,1);
end
