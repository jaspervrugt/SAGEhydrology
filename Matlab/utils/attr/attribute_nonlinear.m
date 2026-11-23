function y = attribute_nonlinear(x,state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ATTRIBUTE_NONLINEAR Apply fitted nonlinear transforms using base MATLAB.
%  Copyright (c) 2026 Jasper A. Vrugt, UCI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = state.method; y = x;
    if m == 5
        if any(x < 0), error('attributes:Log','Log1p requires x >= 0.'); end
        y = log1p(x);
    elseif m == 6
        p = state.lambda; pos = x >= 0;
        if abs(p) < eps
            y(pos) = log1p(x(pos));
        else
            y(pos) = expm1(p*log1p(x(pos)))/p;
        end
        if abs(p-2) < eps
            y(~pos) = -log1p(-x(~pos));
        else
            y(~pos) = -expm1((2-p)*log1p(-x(~pos)))/(2-p);
        end
    elseif ismember(m,[7 8])
        if numel(state.knots) == 1
            y = repmat(0.5,size(x));
        else
            bounded = min(max(x,state.knots(1)),state.knots(end));
            y = interp1(state.knots,state.probabilities,bounded,'linear');
        end
        if m == 8
            y = -sqrt(2)*erfcinv(2*min(max(y,1e-7),1-1e-7));
        end
    end
end
