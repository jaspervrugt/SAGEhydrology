function length = hbv_routing_length(mdl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HBV_ROUTING_LENGTH Fixed 60-day MAXBAS window at the active resolution.
%
% mdl.dt is the number of timesteps per day: 1 daily, 24 hourly, or 96
% for 15-minute data. Legacy callers without dt retain the daily default.
% This window is independent of the fitted routing-width parameter B.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026 Jasper A. Vrugt
% University of California, Irvine | jasper@uci.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = 1;
    if isfield(mdl,'dt') && ...
        ~isempty(mdl.dt)
        dt = mdl.dt;
    end
    if ~(isnumeric(dt) && ...
        isreal(dt) && ...
        isscalar(dt) && ...
        isfinite(dt) && ...
        ismember(dt,[1,24,96]))
        error('hbv_routing_length:Resolution', ...
            'mdl.dt must be 1, 24, or 96 timesteps per day.');
    end
    length = 60 * double(dt);
end
