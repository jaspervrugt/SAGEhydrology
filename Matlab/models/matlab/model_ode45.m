function Z = model_ode45(rhs,tout,z0,options)
%MODEL_ODE45 Integrate each piecewise-constant forcing interval separately.
% All seven public callers use time only to select a forcing row. Hold that
% row fixed through the segment endpoint and restart the adaptive solver at
% each change. Carry the entire augmented state, including sensitivities.
% Local time 0..1 avoids calendar-time rounding; InitialStep applies anew.
    if tout<0 || ...
        tout~=fix(tout)
        error('model_ode45:Intervals', ...
            'tout must be a nonnegative integer.');
    end
    Z=zeros(tout+1,numel(z0));Z(1,:)=z0(:)';
    for interval=1:tout
        % Three tspan entries make the single-interval endpoints explicit.
        [~,history]=ode45(@(~,z) rhs(interval-1,z), ...
            [0,.5,1],Z(interval,:)',options);
        if size(history,1)~=3
            error('model_ode45:Incomplete', ...
                'Interval %d did not complete.',interval);
        end
        Z(interval+1,:)=history(end,:);
    end
end
