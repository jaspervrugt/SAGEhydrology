function [value,out] = zscore_q_result(out,dat,mdl,request)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZSCORE_Q_RESULT Score a legacy Q-only result without changing the backend.
%  The trajectory Jacobian is the physical model derivative. The objective
%  derivative includes the simulated mean and standard-deviation terms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    id = double(dat.stats.Q.train.indices(:));
    y = dat.obs.Q.value(id);
    [value,delta] = zscore_loss(y,out.q(id));
    J = out.jacobian;
    if size(J,1) ~= numel(id)
        error('SAGE:ZscoreJacobian','Training Jacobian length mismatch.');
    end
    gradient = J'*delta;
    if request.gradient
        out.gradient = gradient;
    end
    if request.attribution
        [At,An] = sage_attribution(J,delta,mdl,gradient);
        out.attribution = struct('total',At,'net',An);
    end
end
