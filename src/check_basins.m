function eligibility = check_basins(dat,mdl,bas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK_BASINS Verify forcing and discharge eligibility of selected basins
%
% SYNOPSIS: eligibility = check_basins(dat,mdl,bas)
%   dat         cell structure with basin-specific forcing and discharge
%    {k}.meteo.P   meteorological precipitation vector for basin k
%    {k}.meteo.Ep  potential evapotranspiration vector for basin k
%    {k}.meteo.T   temperature vector for basin k
%    {k}.y_n       observed discharge vector for basin k
%   mdl         structure with model and split information
%    .tout       number of meteorological time steps used by the model
%    .id_train   training-period indices into observed discharge
%    .id_eval    evaluation-period indices into observed discharge
%   bas         basin structure
%    .K          number of selected basins
%   eligibility  basin-level data-quality structure
%    .valid             true when forcing and discharge are complete
%    .forcing_complete  true for finite P, Ep, and T over integration
%    .discharge_complete true when scored discharge contains usable data
%    .coverage_q_train   usable discharge coverage in training period (%)
%    .coverage_q_eval    usable discharge coverage in evaluation period (%)
%    .discharge_variable_train true when finite training Q is nonconstant
%    .discharge_variable_eval  true when finite evaluation Q is nonconstant
%    .discharge_variable true when Q is nonconstant in both periods
%    .runoff_ratio_train long-term Q/P ratio in the training period
%    .runoff_ratio_eval  long-term Q/P ratio in the evaluation period
%    .hydrologic_alert   true when Q/P is outside the diagnostic range
%    .documented_exclusion true for a region-configured excluded gauge
%    .minimum_q_coverage minimum coverage required in each period (%)
%    .reason            exclusion reason (empty for valid basins)
%
% NOTES:
%   1. Meteorological inputs must have mdl.tout entries because they include
%      the full model integration window, including spin-up.
%   2. Observed discharge must be long enough to support all training and
%      evaluation indices used by the loss functions.
%   3. This function is intended to catch reader/indexing inconsistencies
%      immediately after read_Q and before prep_stats or model execution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026 / updated Jun. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K = bas.K;
    incompleteForcing = false(K,1);
    incompleteDischarge = false(K,1);
    coverageQTrain = nan(K,1);
    coverageQEval = nan(K,1);
    dischargeVariableTrain = false(K,1);
    dischargeVariableEval = false(K,1);
    runoffRatioTrain = nan(K,1);
    runoffRatioEval = nan(K,1);
    documentedExclusion = false(K,1);
    documentedReason = strings(K,1);
    minimumQCoverage = 5;
    maximumRunoffRatio = 1.25;

    [excludedGauge,excludedReason] = ...
        local_documented_exclusions(bas);
    
    for k = 1:K
        nP = numel(dat{k}.meteo.P);
        nT = numel(dat{k}.meteo.T);
        nE = numel(dat{k}.meteo.Ep);
        nQ = numel(dat{k}.y_n);
    
        if nP ~= mdl.tout ...
                || nT ~= mdl.tout ...
                || nE ~= mdl.tout
            error(['Length mismatch basin %d: ' ...
                'P=%d, T=%d, Ep=%d, ' ...
                'expected mdl.tout=%d'], ...
                k,nP,nT,nE,mdl.tout);
        end
    
        if nQ < max(mdl.id_train)
            error(['Length mismatch basin %d: ' ...
                'Q=%d, max(id_train)=%d'], ...
                k,nQ,max(mdl.id_train));
        end

        incompleteForcing(k) = ...
            any(~isfinite(dat{k}.meteo.P(:))) ...
            || any(~isfinite(dat{k}.meteo.Ep(:))) ...
            || any(~isfinite(dat{k}.meteo.T(:)));

        % Compute training/evaluation coverage separately. Basin-local
        % indices take precedence for rainfall-block splitting.
        idTrain = local_period_indices(dat{k},mdl,'train');
        idEval = local_period_indices(dat{k},mdl,'eval');
        coverageQTrain(k) = local_q_coverage(dat{k}.y_n,idTrain);
        coverageQEval(k) = local_q_coverage(dat{k}.y_n,idEval);

        % Use one paired basin population in every ECDF scenario. A basin
        % therefore needs sufficient, nonconstant discharge in both the
        % training and evaluation periods, irrespective of whether it is a
        % training or an evaluation basin spatially.
        dischargeVariableTrain(k) = local_q_is_variable( ...
            dat{k}.y_n,idTrain);
        dischargeVariableEval(k) = local_q_is_variable( ...
            dat{k}.y_n,idEval);
        runoffRatioTrain(k) = local_runoff_ratio( ...
            dat{k}.meteo.P,dat{k}.y_n,idTrain);
        runoffRatioEval(k) = local_runoff_ratio( ...
            dat{k}.meteo.P,dat{k}.y_n,idEval);
        if ~isempty(excludedGauge) && isfield(bas,'id_gauge') ...
                && numel(bas.id_gauge) >= k
            gauge = local_normalize_gauge(bas.id_gauge(k));
            match = find(excludedGauge == gauge,1);
            if ~isempty(match)
                documentedExclusion(k) = true;
                documentedReason(k) = excludedReason(match);
            end
        end
        incompleteDischarge(k) = ...
            ~isfinite(coverageQTrain(k)) ...
            || coverageQTrain(k) < minimumQCoverage ...
            || ~isfinite(coverageQEval(k)) ...
            || coverageQEval(k) < minimumQCoverage ...
            || ~dischargeVariableTrain(k) ...
            || ~dischargeVariableEval(k);
    end

    eligibility = struct();
    eligibility.forcing_complete = ~incompleteForcing;
    eligibility.discharge_complete = ~incompleteDischarge;
    eligibility.coverage_q_train = coverageQTrain;
    eligibility.coverage_q_eval = coverageQEval;
    eligibility.discharge_variable_train = dischargeVariableTrain;
    eligibility.discharge_variable_eval = dischargeVariableEval;
    eligibility.discharge_variable = ...
        dischargeVariableTrain & dischargeVariableEval;
    eligibility.minimum_q_coverage = minimumQCoverage;
    eligibility.runoff_ratio_train = runoffRatioTrain;
    eligibility.runoff_ratio_eval = runoffRatioEval;
    eligibility.maximum_runoff_ratio = maximumRunoffRatio;
    eligibility.hydrologic_alert = ...
        (isfinite(runoffRatioTrain) ...
        & (runoffRatioTrain < 0 ...
        | runoffRatioTrain > maximumRunoffRatio)) ...
        | (isfinite(runoffRatioEval) ...
        & (runoffRatioEval < 0 ...
        | runoffRatioEval > maximumRunoffRatio));
    eligibility.documented_exclusion = documentedExclusion;
    eligibility.valid = ~incompleteForcing & ~incompleteDischarge ...
        & ~documentedExclusion;
    eligibility.reason = strings(K,1);
    eligibility.reason(incompleteForcing) = "incomplete forcing";
    onlyQ = ~incompleteForcing & incompleteDischarge;
    lowCoverage = incompleteDischarge & ...
        (~isfinite(coverageQTrain) ...
        | coverageQTrain < minimumQCoverage ...
        | ~isfinite(coverageQEval) ...
        | coverageQEval < minimumQCoverage);
    constantQ = incompleteDischarge & ~lowCoverage ...
        & (~dischargeVariableTrain | ~dischargeVariableEval);
    eligibility.reason(onlyQ & lowCoverage) = ...
        "less than 5% discharge in one or both periods";
    eligibility.reason(onlyQ & constantQ) = ...
        "constant discharge in one or both periods";
    both = incompleteForcing & incompleteDischarge;
    eligibility.reason(both & lowCoverage) = ...
        "incomplete forcing and less than 5% discharge in one or both periods";
    eligibility.reason(both & constantQ) = ...
        "incomplete forcing and constant discharge in one or both periods";
    eligibility.reason(documentedExclusion) = documentedReason( ...
        documentedExclusion);

    if any(incompleteForcing)
        basinIndex = find(incompleteForcing);
        basinLabel = string(basinIndex);
        try
            if isfield(bas,'id_gauge') ...
                    && numel(bas.id_gauge) >= K
                basinLabel = string(bas.id_gauge(basinIndex));
            end
        catch
        end
        nShow = min(10,numel(basinLabel));
        example = strjoin(cellstr(basinLabel(1:nShow)),', ');
        if numel(basinLabel) > nShow
            example = sprintf('%s, ...',example);
        end
        fprintf(['      Data-quality screening: %d of %d basins excluded.\n' ...
                 '      Reason : incomplete P, Ep, or T over the ' ...
                 'integration window\n' ...
                 '      Gauges : %s\n'], ...
                nnz(incompleteForcing),K,example);
    end

    if any(incompleteDischarge)
        basinIndex = find(incompleteDischarge);
        basinLabel = string(basinIndex);
        try
            if isfield(bas,'id_gauge') && numel(bas.id_gauge) >= K
                basinLabel = string(bas.id_gauge(basinIndex));
            end
        catch
        end
        nShow = min(10,numel(basinLabel));
        example = strjoin(cellstr(basinLabel(1:nShow)),', ');
        if numel(basinLabel) > nShow
            example = sprintf('%s, ...',example);
        end
        fprintf(['      Data-quality screening: %d of %d basins excluded.\n' ...
                 '      Reason : < %.1f%% usable discharge or constant ' ...
                 'discharge in >= 1 operative period\n' ...
                 '      Gauges : %s\n'], ...
                nnz(incompleteDischarge),K,minimumQCoverage,example);
    end
end

function ratio = local_runoff_ratio(P,Q,index)
%LOCAL_RUNOFF_RATIO Return paired period-total discharge/precipitation.
    ratio = NaN;
    if isempty(index)
        return
    end
    P = double(P(:));
    Q = double(Q(:));
    index = index(index <= numel(P) & index <= numel(Q));
    if isempty(index)
        return
    end
    p = P(index);
    q = Q(index);
    paired = isfinite(p) & p >= 0 & isfinite(q) & q >= 0;
    if ~any(paired)
        return
    end
    totalP = sum(p(paired));
    if totalP > 0
        ratio = sum(q(paired))/totalP;
    end
end

function [gauge,reason] = local_documented_exclusions(bas)
%LOCAL_DOCUMENTED_EXCLUSIONS Read optional region-specific exclusions.
    gauge = strings(0,1);
    reason = strings(0,1);
    if ~isfield(bas,'excluded_gauges') ...
            || isempty(bas.excluded_gauges)
        return
    end
    rawGauge = string(bas.excluded_gauges(:));
    gauge = strings(size(rawGauge));
    for j = 1:numel(rawGauge)
        gauge(j) = local_normalize_gauge(rawGauge(j));
    end
    defaultReason = "documented regional data-quality exclusion";
    reason = repmat(defaultReason,numel(gauge),1);
    if isfield(bas,'exclusion_reason') ...
            && ~isempty(bas.exclusion_reason)
        supplied = string(bas.exclusion_reason(:));
        if isscalar(supplied)
            reason(:) = supplied;
        elseif numel(supplied) == numel(gauge)
            reason = supplied;
        else
            error('SAGE:dataQuality:ExclusionReasonSize', ...
                ['bas.exclusion_reason must contain one entry or one ' ...
                'entry per excluded gauge.']);
        end
    end
end

function gauge = local_normalize_gauge(value)
%LOCAL_NORMALIZE_GAUGE Compare numeric gauge IDs independent of padding.
    gauge = upper(strip(string(value)));
    if ~isempty(regexp(char(gauge),'^\d+$','once'))
        gauge = string(sprintf('%d',str2double(gauge)));
    end
end

function index = local_period_indices(entry,mdl,period)
%LOCAL_PERIOD_INDICES Return basin-local or model-wide scoring indices.
    field = ['id_' period];
    index = [];
    if isstruct(entry) && isfield(entry,field) && ~isempty(entry.(field))
        index = entry.(field)(:);
    elseif isstruct(mdl) && isfield(mdl,field) && ~isempty(mdl.(field))
        index = mdl.(field)(:);
        % Model-wide contiguous periods are stored compactly as
        % [first last]; expand them before computing coverage/statistics.
        if numel(index) == 2 && index(2) >= index(1)
            index = (index(1):index(2)).';
        end
    end
    index = unique(double(index(isfinite(index) & index >= 1 ...
        & index == floor(index))));
end

function tf = local_q_is_variable(Q,index)
%LOCAL_Q_IS_VARIABLE True for at least two unequal finite, nonnegative Q.
    tf = false;
    if isempty(index)
        return
    end
    Q = Q(:);
    index = index(index <= numel(Q));
    q = double(Q(index));
    q = q(isfinite(q) & q >= 0);
    if numel(q) < 2
        return
    end
    tf = max(q) > min(q);
end

function coverage = local_q_coverage(Q,index)
%LOCAL_Q_COVERAGE Percentage finite, nonnegative observations in a period.
    coverage = NaN;
    if isempty(index)
        return
    end
    Q = Q(:);
    usable = false(numel(index),1);
    inside = index <= numel(Q);
    if any(inside)
        q = Q(index(inside));
        usable(inside) = isfinite(q) & q >= 0;
    end
    coverage = 100*sum(usable)/numel(index);
end
