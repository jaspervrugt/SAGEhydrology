function varargout = parameter_scale(action,mdl,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETER_SCALE Map normalized/physical parameters with analytic slopes.
%  Each parameter independently uses Linear or Logarithmic. pspace retains
%  its meaning: 0 physical, 1 normalized, 2 unconstrained sigmoid input.
%  [theta,slope] = parameter_scale('physical',mdl,x)
%  x = parameter_scale('normalized',mdl,theta)
%  modes = parameter_scale('modes',mdl)
%  bounds = parameter_scale('bounds',mdl)
%  mdl = parameter_scale('configure',mdl,bounds)
%  For CaRT, bounds reports physical values while the solver receives
%  the same encoded coordinates as before.
%  Bounds and parameter dimensions are columns; basin columns may repeat.
%  Uses base MATLAB. Copyright (c) 2026 Jasper A. Vrugt, UCI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(action,'configure')
        varargout = {configure_bounds(mdl,x)};
        return
    end

    lo = double(mdl.th_min(:));
    hi = double(mdl.th_max(:));
    n = numel(lo);
    modes = repmat("linear",n,1);
    if isfield(mdl,'parameter_scales') && ...
        ~isempty(mdl.parameter_scales)
        modes = lower(string(mdl.parameter_scales(:)));
    end
    encoded = cart_log_mask(mdl,n);
    linearPhysical = modes == "physical_linear";
    if numel(hi) ~= n || ...
        numel(modes) ~= n || ...
        any(~isfinite(lo) | ...
        ~isfinite(hi) | ...
        hi <= lo) || ...
        any(~ismember(modes,["linear","log","physical_linear"])) || ...
        any(linearPhysical & ~encoded)
        error('parameter_scale:Invalid','Invalid bounds or scale modes.');
    end
    logged = modes == "log";
    if any(lo(logged) <= 0)
        error('parameter_scale:LogBounds', ...
            ['Log parameters require strictly positive lower ' ...
                'bounds.']);
    end
    physicalLo = exp(lo(linearPhysical));
    physicalHi = exp(hi(linearPhysical));
    if any(physicalLo <= 0 | ...
        ~isfinite(physicalHi) | ...
        physicalHi <= physicalLo)
        error('parameter_scale:PhysicalBounds', ...
            'CaRT physical bounds must be finite, positive and distinct.');
    end
    if strcmpi(action,'bounds')
        if nargout > 1
            error('parameter_scale:Outputs','Request one bounds struct.');
        end
        lo(encoded) = exp(lo(encoded));
        hi(encoded) = exp(hi(encoded));
        scales = modes;
        scales(logged | (encoded & modes == "linear")) = "logarithmic";
        scales(linearPhysical) = "linear";
        scales(encoded & logged) = "legacy_log";
        bounds = struct('lower',lo,'upper',hi,'scales',scales);
        if isfield(mdl,'par_names')
            bounds.names = string(mdl.par_names(:));
            bounds.names(encoded) = extractAfter(bounds.names(encoded),4);
        end
        if isfield(mdl,'par_units')
            bounds.units = string(mdl.par_units(:));
            strip = encoded & startsWith(bounds.units,"log of ");
            bounds.units(strip) = extractAfter(bounds.units(strip),7);
        end
        varargout = {bounds};
        return
    end
    if strcmpi(action,'modes')
        varargout = {modes}; return
    end
    if size(x,1) ~= n || ...
        any(~isfinite(x(:)))
        error('parameter_scale:Size','Invalid parameter matrix.');
    end
    span = hi-lo;
    if strcmpi(action,'normalized')
        if any(x(logged,:) <= 0,'all')
            error('parameter_scale:LogValue', ...
                'Physical Log parameters must be positive.');
        end
        y = (x-lo)./span;
        y(logged,:) = (log(x(logged,:))-log(lo(logged))) ...
            ./log(hi(logged)./lo(logged));
        y(linearPhysical,:) = (exp(x(linearPhysical,:))-physicalLo) ...
            ./(physicalHi-physicalLo);
        varargout = {y}; return
    end
    if ~strcmpi(action,'physical')
        error('parameter_scale:Action','Unknown mapping action.');
    end
    space = 1;
    if isfield(mdl,'pspace'), space = mdl.pspace; end
    if space == 0
        varargout = {x,ones(size(x))}; return
    elseif space == 2
        x = 1./(1+exp(-x));
        outer = x.*(1-x);
    elseif space == 1
        outer = ones(size(x));
    else
        error('parameter_scale:Space','Unknown parameter space.');
    end
    y = lo+x.*span;
    slope = repmat(span,1,size(x,2));
    y(logged,:) = lo(logged).*exp(x(logged,:) ...
        .*log(hi(logged)./lo(logged)));
    valid = x(logged,:) >= 0 & ...
        x(logged,:) <= 1;
    z = y(logged,:);
    bounded = min(max(z,lo(logged)),hi(logged));
    z(valid) = bounded(valid);
    a = x(logged,:) == 0; b = x(logged,:) == 1;
    lowerBound = repmat(lo(logged),1,size(x,2));
    upper = repmat(hi(logged),1,size(x,2));
    z(a) = lowerBound(a); z(b) = upper(b);
    y(logged,:) = z;
    slope(logged,:) = z.*log(hi(logged)./lo(logged));
    weights = x(linearPhysical,:);
    physical = (1-weights).*physicalLo+weights.*physicalHi;
    y(linearPhysical,:) = log(physical);
    y(linearPhysical,:) = replace_endpoints( ...
        y(linearPhysical,:),weights,lo(linearPhysical),hi(linearPhysical));
    slope(linearPhysical,:) = (physicalHi-physicalLo)./physical;
    varargout = {y,slope.*outer};
end

function encoded = cart_log_mask(mdl,n)
    encoded = false(n,1);
    if isfield(mdl,'model') && ...
        ismember(mdl.model,[18,19]) && ...
        isfield(mdl,'par_names')
        encoded = startsWith(string(mdl.par_names(:)),"\ln ");
    end
end

function mdl = configure_bounds(mdl,bounds)
    lo = double(bounds.lower(:));
    hi = double(bounds.upper(:));
    scales = lower(string(bounds.scales(:)));
    n = numel(mdl.th_min);
    encoded = cart_log_mask(mdl,n);
    if numel(lo) ~= n || ...
        numel(hi) ~= n || ...
        numel(scales) ~= n || ...
        any(~ismember(scales, ...
            ["linear","log","logarithmic","legacy_log"])) || ...
        any(scales(~encoded) == "legacy_log")
        error('parameter_scale:Invalid','Invalid bounds or scale modes.');
    end
    if any(lo(encoded) <= 0 | hi(encoded) <= 0)
        error('parameter_scale:PhysicalBounds', ...
            'CaRT log-encoded parameters require positive bounds.');
    end
    modes = scales;
    modes(scales == "logarithmic" | scales == "legacy_log") = "log";
    modes(encoded & modes == "log" & scales ~= "legacy_log") = "linear";
    modes(encoded & scales == "linear") = "physical_linear";
    originalLo = mdl.th_min(:);
    originalHi = mdl.th_max(:);
    lo(encoded) = log(lo(encoded));
    hi(encoded) = log(hi(encoded));
    % Preserve exact encoded defaults on an unchanged round trip.
    sameLo = encoded & bounds.lower(:) == exp(originalLo);
    sameHi = encoded & bounds.upper(:) == exp(originalHi);
    lo(sameLo) = originalLo(sameLo);
    hi(sameHi) = originalHi(sameHi);
    mdl.th_min = lo;
    mdl.th_max = hi;
    mdl.parameter_scales = modes;
    parameter_scale('modes',mdl);
end

function values = replace_endpoints(values,x,lo,hi)
    first = x == 0;
    last = x == 1;
    lo = repmat(lo,1,size(x,2));
    hi = repmat(hi,1,size(x,2));
    values(first) = lo(first);
    values(last) = hi(last);
end
