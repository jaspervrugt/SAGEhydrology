function [value,gradient,detail] = ...
    observation_loss(name,sim,J,dat,loss,period)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBSERVATION_LOSS Evaluate a named observation loss.
%
%  Applies the selected loss to Q, SWE, or SM and computes its gradient.
%
% SYNOPSIS:
%   [value,gradient,detail] = ...
%       observation_loss(name,sim,J,dat,loss,period)
%
% INPUT ARGUMENTS:
%   name            observation name: 'Q', 'SWE', or 'SM'
%   sim             complete simulated trajectory
%   J               complete trajectory Jacobian, n-by-d
%   dat             basin observations and prepared statistics
%   loss            selected loss-function settings
%   period          optional 'train' or 'eval'; default 'train'
%
% OUTPUT ARGUMENTS:
%   value           scalar loss for the selected period
%   gradient        d-by-1 loss gradient
%   detail          availability, indices, metrics, and FDC diagnostics
%
% NOTES:
%   Missing observations are excluded using indices from PREP_STATS. The
%   gradient is J'*(dL/dsim) in the same parameter space as J.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Sep. 2026                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6 ...
            || isempty(period)
        period = 'train';
    end
    period = lower(char(string(period)));
    if ~ismember(period,{'train','eval'})
        error('observation_loss:BadPeriod', ...
            'Period must be ''train'' or ''eval''.');
    end

    name = upper(string(name));
    if ~isscalar(name) ...
            || ~ismember(name,["Q","SWE","SM"])
        error('observation_loss:UnknownObservation', ...
            'Observation name must be Q, SWE, or SM.');
    end

    field = char(name);
    d = size(J,2);
    value = NaN;
    gradient = nan(d,1);
    detail = struct('name',name, ...
        'period',period, ...
        'indices',zeros(0,1), ...
        'n',0, ...
        'available',false, ...
        'delta',zeros(0,1), ...
        'NSE',NaN, ...
        'KGE',NaN, ...
        'KGE_r',NaN, ...
        'KGE_alpha',NaN, ...
        'KGE_beta',NaN, ...
        'D_fdc',NaN, ...
        'D_p',NaN, ...
        'D_logp',NaN, ...
        'JKGE',NaN, ...
        'JKGE_M',NaN, ...
        'JKGE_V',NaN, ...
        'JKGE_C',NaN);

    if ~isfield(dat,'obs') ...
            || ~isfield(dat.obs,field) ...
            || isempty(dat.obs.(field)) ...
            || ~isfield(dat.obs.(field),'value') ...
            || isempty(dat.obs.(field).value) ...
            || ~isfield(dat,'stats') ...
            || ~isfield(dat.stats,field) ...
            || ~isfield(dat.stats.(field),period)
        return
    end

    stats = dat.stats.(field).(period);
    if ~isfield(stats,'indices') ...
            || isempty(stats.indices)
        return
    end

    id = double(stats.indices(:));
    y = double(dat.obs.(field).value(:));
    sim = double(sim(:));
    J = double(J);

    if numel(sim) ~= numel(y) ...
            || size(J,1) ~= numel(sim)
        error('observation_loss:SizeMismatch', ...
            ['Simulation, observation, and ' ...
            'Jacobian lengths differ ' ...
             'for %s.'],field);
    end
    if any(id < 1) ...
            || any(id > numel(y)) ...
            || any(id ~= fix(id))
        error('observation_loss:BadIndex', ...
            ['Prepared %s training indices ' ...
            'are invalid.'],field);
    end

    good = isfinite(y(id)) ...
        & isfinite(sim(id)) ...
        & all(isfinite(J(id,:)),2);
    id = id(good);
    if isempty(id)
        return
    end

    yt = y(id);
    qt = sim(id);
    Jt = J(id,:);
    res = yt - qt;
    fnc = double(loss.fnc);
    delta = [];

    if isfinite(stats.TSS) ...
            && stats.TSS > 0
        detail.NSE = 1-sum(res.^2)/stats.TSS;
    end
    if isfinite(stats.mean) ...
            && isfinite(stats.std) ...
            && stats.std > 0
        [detail.KGE,detail.KGE_r, ...
            detail.KGE_alpha,detail.KGE_beta] = ...
            kge(yt,qt,stats.mean,stats.std);
    end
    fdcKey = period(1);
    if isfield(dat,'fdc') ...
            && isfield(dat.fdc,field) ...
            && isfield(dat.fdc.(field),fdcKey)
        fdc = dat.fdc.(field).(fdcKey);
        if numel(qt) == fdc.n
            [detail.D_fdc,detail.D_p,detail.D_logp] = ...
                fdc_metrics_cached(qt,fdc);
        end
    end

    switch fnc
        case 1
            value = sum(abs(res));
            delta = delta_n(1,yt,qt);

        case 2
            value = sum(res.^2);
            delta = delta_n(2,yt,qt,1);

        case 3
            if isfinite(stats.TSS) ...
                    && stats.TSS > 0
                value = sum(res.^2)/stats.TSS;
                delta = -2*res/stats.TSS;
            end

        case 4
            if isfinite(stats.mean) ...
                    && isfinite(stats.std) ...
                    && stats.std > 0
                score = detail.KGE;
                value = 1-score;
                delta = delta_n(4,yt,qt);
            end

        case 5
            if isfinite(stats.huber_scale) ...
                    && stats.huber_scale > 0
                [value,delta] = huber_loss( ...
                    res,stats.huber_scale);
            end

        case 6
            if ~isfield(dat,'fdc') ...
                    || ~isfield(dat.fdc,field) ...
                    || ~isfield(dat.fdc.(field),fdcKey)
                error('observation_loss:MissingFDC', ...
                    ['Prepared %s FDC statistics ' ...
                    'are unavailable.'],field);
            end
            fdc = dat.fdc.(field).(fdcKey);
            if numel(qt) ~= fdc.n
                error('observation_loss:FDCSize', ...
                    ['%s FDC cache length does ' ...
                    'not match observations.'],field);
            end
            values = [detail.D_fdc,detail.D_p,detail.D_logp];
            formulation = local_fdc_formulation(loss);
            value = values(formulation);
            delta = delta_n(6,yt,qt,formulation);

        case 7
            jkgeData = local_jkge_data(dat,field);
            if isempty(jkgeData)
                error('observation_loss:MissingJKGE', ...
                    ['Prepared JKGE information is unavailable ' ...
                     'for %s. Run PREP_STATS after selecting ' ...
                     'the training observations.'],field);
            end
            method = double(loss.method);
            if method < 3
                aux = loss.n_win;
            elseif method == 4
                aux = loss.meta.mo_all;
            else
                aux = [];
            end
            Mdef = 2;
            if isfield(loss,'M') ...
                    && ~isempty(loss.M)
                Mdef = loss.M;
            end
            [score,deltaFull,~,detail.JKGE_M, ...
                detail.JKGE_V,detail.JKGE_C] = ...
                jkge_grad(y,sim,jkgeData.m_y,id, ...
                method,aux,jkgeData.cache,Mdef);
            detail.JKGE = score;
            value = 1-score;
            delta = deltaFull;

        otherwise
            error('observation_loss:UnknownLoss', ...
                'Unknown loss function index: %g.',fnc);
    end

    if isfinite(value) ...
            && ~isempty(delta) ...
            && all(isfinite(delta))
        if fnc == 7 ...
                && numel(delta) == size(J,1)
            % JKGE's benchmark transformation uses the complete simulated
            % trajectory. Retain that full chain rule even though the
            % objective score is evaluated on the selected period.
            gradient = J'*delta;
        else
            gradient = Jt'*delta;
        end
        if numel(delta) == size(J,1)
            detail.delta = delta(:);
        else
            detail.delta = zeros(size(J,1),1);
            detail.delta(id) = delta(:);
        end
        detail.indices = id;
        detail.n = numel(id);
        detail.available = all(isfinite(gradient));
    end
end

function value = local_jkge_data(dat,field)
%LOCAL_JKGE_DATA Return prepared named JKGE data.

    value = [];
    if ~isfield(dat,'jkge') ...
            || ~isstruct(dat.jkge)
        return
    end
    if isfield(dat.jkge,field) ...
            && isstruct(dat.jkge.(field)) ...
            && isfield(dat.jkge.(field),'m_y') ...
            && isfield(dat.jkge.(field),'cache')
        value = dat.jkge.(field);
        return
    end
    if strcmp(field,'Q') ...
            && isfield(dat.jkge,'m_y') ...
            && isfield(dat.jkge,'cache')
        value = struct('m_y',dat.jkge.m_y, ...
            'cache',dat.jkge.cache);
    end
end

function formulation = local_fdc_formulation(loss)
%LOCAL_FDC_FORMULATION Return the selected FDC formulation.

    formulation = 1;
    if isfield(loss,'fdc') ...
            && isstruct(loss.fdc) ...
            && isfield(loss.fdc,'formulation') ...
            && ~isempty(loss.fdc.formulation)
        formulation = loss.fdc.formulation;
    elseif isfield(loss,'formulation') ...
            && ~isempty(loss.formulation)
        formulation = loss.formulation;
    end

    if isnumeric(formulation) ...
            && isscalar(formulation) ...
            && isfinite(formulation)
        formulation = double(formulation);
    else
        key = lower(regexprep(char(string(formulation)), ...
            '[^a-z0-9]',''));
        switch key
            case {'1','a','fdc', ...
                    'dfdc','physical','physicalcdf'}
                formulation = 1;
            case {'2','b','p','dp','quantile'}
                formulation = 2;
            case {'3','c','logp','dlogp','logquantile'}
                formulation = 3;
            otherwise
                formulation = NaN;
        end
    end

    if ~ismember(formulation,1:3)
        error('observation_loss:BadFDCFormulation', ...
            'FDC formulation must be 1, 2, or 3.');
    end
end
