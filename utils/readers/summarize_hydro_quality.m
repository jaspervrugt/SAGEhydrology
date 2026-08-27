function quality = summarize_hydro_quality(dat,bas,split)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUMMARIZE_HYDRO_QUALITY Summarize meteorological and discharge quality.
%
% SYNOPSIS: quality = summarize_hydro_quality(dat,bas,split)
%   dat       Kx1 SAGE time-series cell array
%   bas       structure with selected basin information
%   split     optional structure with scored-period indices in .idx
%   quality   compact structure used by the SAGE data-quality dashboard
%
% DESCRIPTION:
%   Coverage is calculated directly from the canonical P, Ep, T, and y_n
%   vectors produced by the generic readers. Existing bad fields are
%   reported separately as consistency information. No source data are
%   changed and no complete time series are retained in the summary.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K = numel(dat);
    gauge = strings(K,1);
    use = strings(K,1);
    nMeteo = zeros(K,1);
    nQ = zeros(K,1);
    nPBad = zeros(K,1);
    nEpBad = zeros(K,1);
    nTBad = zeros(K,1);
    nQBad = zeros(K,1);
    nMeteoBad = zeros(K,1);
    nReaderQBad = zeros(K,1);
    nForcingBad = zeros(K,1);
    longestQGap = zeros(K,1);
    meanP = nan(K,1);
    meanEp = nan(K,1);
    meanQ = nan(K,1);
    runoffRatioAnnualMean = nan(K,1);
    runoffRatioAnnualStd = nan(K,1);
    runoffRatioAnnualYears = zeros(K,1);

    if nargin < 3
        split = struct();
    end

    for k = 1:K
        entry = dat{k};
        gauge(k) = local_gauge(entry,bas,k);
        use(k) = local_use(entry,bas,k);

        if isstruct(entry) ...
                && isfield(entry,'meteo') ...
                && isstruct(entry.meteo)
            P = local_vector(entry.meteo,'P');
            Ep = local_vector(entry.meteo,'Ep');
            T = local_vector(entry.meteo,'T');
            nMeteo(k) = max([numel(P),numel(Ep),numel(T)]);
            nPBad(k) = local_invalid_count(P,nMeteo(k),false);
            nEpBad(k) = local_invalid_count(Ep,nMeteo(k),false);
            nTBad(k) = local_invalid_count(T,nMeteo(k),false);
            nForcingBad(k) = local_combined_invalid_count( ...
                P,Ep,T,nMeteo(k));
            if isfield(entry.meteo,'bad')
                nMeteoBad(k) = local_bad_count( ...
                    entry.meteo.bad,nMeteo(k));
            end
        end

        Q = local_vector(entry,'y_n');
        nQ(k) = numel(Q);
        qInvalid = ~isfinite(Q) | Q < 0;
        nQBad(k) = sum(qInvalid);
        longestQGap(k) = local_longest_run(qInvalid);
        if isstruct(entry) ...
                && isfield(entry,'bad')
            nReaderQBad(k) = local_bad_count(entry.bad,nQ(k));
        end

        if exist('P','var')
            [Pscore,Epscore] = local_scored_forcing(P,Ep,Q,split);
            nPair = min(numel(Pscore),numel(Q));
            pair = isfinite(Pscore(1:nPair)) ...
                & isfinite(Q(1:nPair)) ...
                & Q(1:nPair) >= 0;
            meanP(k) = mean(Pscore(pair),'omitnan');
            meanQ(k) = mean(Q(pair),'omitnan');
            if numel(Epscore) >= nPair
                epPair = pair & isfinite(Epscore(1:nPair));
                meanEp(k) = mean(Epscore(epPair),'omitnan');
            end
            annual = local_annual_runoff_ratios( ...
                Pscore(1:nPair),Q(1:nPair),pair,split);
            runoffRatioAnnualMean(k) = mean(annual,'omitnan');
            runoffRatioAnnualStd(k) = std(annual,0,'omitnan');
            runoffRatioAnnualYears(k) = sum(isfinite(annual));
        else
            meanQ(k) = mean(Q(isfinite(Q) & Q >= 0),'omitnan');
        end
        clear P Ep T
    end

    quality = struct();
    quality.gauge = gauge;
    quality.use = use;
    quality.n_meteo = nMeteo;
    quality.n_q = nQ;
    quality.n_bad_p = nPBad;
    quality.n_bad_ep = nEpBad;
    quality.n_bad_t = nTBad;
    quality.n_bad_q = nQBad;
    quality.n_bad_meteo = nMeteoBad;
    quality.n_bad_q_reader = nReaderQBad;
    quality.n_bad_forcing = nForcingBad;
    quality.longest_q_gap = longestQGap;
    quality.coverage_p = local_coverage(nPBad,nMeteo);
    quality.coverage_ep = local_coverage(nEpBad,nMeteo);
    quality.coverage_t = local_coverage(nTBad,nMeteo);
    quality.coverage_q = local_coverage(nQBad,nQ);
    quality.coverage_meteo = local_coverage(nForcingBad,nMeteo);
    quality.meteo_bad_consistent = nMeteoBad == nForcingBad;
    quality.q_bad_consistent = nReaderQBad == nQBad;
    quality.mean_p = meanP;
    quality.mean_ep = meanEp;
    quality.mean_q = meanQ;
    quality.runoff_ratio = meanQ./meanP;
    quality.runoff_ratio_annual_mean = runoffRatioAnnualMean;
    quality.runoff_ratio_annual_std = runoffRatioAnnualStd;
    quality.runoff_ratio_annual_years = runoffRatioAnnualYears;
    quality.aridity_ratio = meanEp./meanP;
    quality.water_balance_remainder = meanP - meanQ;
    quality.hydrologic_flag = isfinite(quality.runoff_ratio) ...
        & (quality.runoff_ratio < 0 | quality.runoff_ratio > 1.25);
    quality.hydrologic_alert_uncertain = quality.hydrologic_flag ...
        & quality.coverage_q < 30;
    quality.created = datetime('now');

end

function ratio = local_annual_runoff_ratios(P,Q,pair,split)

    ratio = [];
    n = numel(pair);
    if n == 0 ...
            || ~isstruct(split) ...
            || ~isfield(split,'dt0') ...
            || ~isfield(split,'dt')
        return
    end
    if isfield(split,'idx') ...
            && ~isempty(split.idx)
        idx = double(split.idx(:));
    else
        idx = (1:n).';
    end
    if numel(idx) < n
        return
    end
    idx = idx(1:n);
    switch double(split.dt)
        case 1
            time = split.dt0 + days(idx - 1);
        case 24
            time = split.dt0 + hours(idx - 1);
        case 96
            time = split.dt0 + minutes(15*(idx - 1));
        otherwise
            return
    end
    waterYear = year(time);
    waterYear(month(time) >= 10) = waterYear(month(time) >= 10) + 1;
    years = unique(waterYear(:));
    ratio = nan(numel(years),1);
    for j = 1:numel(years)
        inYear = waterYear == years(j) & isfinite(P);
        valid = inYear & pair;
        nExpected = sum(inYear);
        if nExpected == 0 ...
                || 100*sum(valid)/nExpected < 30
            continue
        end
        pTotal = sum(P(valid));
        if isfinite(pTotal) ...
                && pTotal > 0
            ratio(j) = sum(Q(valid))/pTotal;
        end
    end

end

function [Pscore,Epscore] = local_scored_forcing(P,Ep,Q,split)

    Pscore = P;
    Epscore = Ep;
    if isstruct(split) ...
            && isfield(split,'idx') ...
            && ~isempty(split.idx)
        idx = double(split.idx(:));
        idx = idx(isfinite(idx) & idx >= 1 ...
            & idx == floor(idx));
        if ~isempty(idx) ...
                && max(idx) <= numel(P)
            Pscore = P(idx);
        end
        if ~isempty(idx) ...
                && max(idx) <= numel(Ep)
            Epscore = Ep(idx);
        end
    elseif ~isempty(Q) ...
            && numel(Q) == numel(P)
        Pscore = P;
        Epscore = Ep;
    end

end

function n = local_combined_invalid_count(P,Ep,T,nExpected)

    if nExpected == 0
        n = 0;
        return
    end
    invalid = true(nExpected,3);
    values = {P,Ep,T};
    for j = 1:3
        x = values{j};
        nx = min(numel(x),nExpected);
        if nx > 0
            invalid(1:nx,j) = ~isfinite(x(1:nx));
        end
    end
    n = sum(any(invalid,2));

end

function x = local_vector(s,name)

    x = [];
    if isstruct(s) ...
            && isfield(s,name) ...
            && isnumeric(s.(name))
        x = double(s.(name)(:));
    end

end

function n = local_invalid_count(x,nExpected,negativeInvalid)

    if nargin < 3
        negativeInvalid = false;
    end
    if nExpected == 0
        n = 0;
        return
    end
    n = nExpected - numel(x) + sum(~isfinite(x));
    if negativeInvalid
        n = n + sum(isfinite(x) & x < 0);
    end

end

function n = local_bad_count(bad,nExpected)

    if isempty(bad)
        n = 0;
    elseif islogical(bad)
        n = sum(bad(:));
    elseif isnumeric(bad)
        idx = double(bad(:));
        idx = idx(isfinite(idx) & idx >= 1 & idx <= nExpected);
        n = numel(unique(idx));
    else
        n = 0;
    end

end

function coverage = local_coverage(nBad,nTotal)

    coverage = nan(size(nTotal));
    ok = nTotal > 0;
    coverage(ok) = 100*(1 - nBad(ok)./nTotal(ok));
    coverage(ok) = min(100,max(0,coverage(ok)));

end

function n = local_longest_run(mask)

    mask = logical(mask(:));
    if isempty(mask) ...
            || ~any(mask)
        n = 0;
        return
    end
    edge = diff([false;mask;false]);
    first = find(edge == 1);
    last = find(edge == -1) - 1;
    n = max(last - first + 1);

end

function gauge = local_gauge(entry,bas,k)

    gauge = "";
    if isstruct(entry) ...
            && isfield(entry,'gauge') ...
            && ~isempty(entry.gauge)
        gauge = string(entry.gauge);
    elseif isstruct(bas) ...
            && isfield(bas,'id_gauge') ...
            && numel(bas.id_gauge) >= k
        ids = bas.id_gauge;
        if iscell(ids)
            gauge = string(ids{k});
        else
            gauge = string(ids(k));
        end
    end

end

function use = local_use(entry,bas,k)

    use = "";
    if isstruct(entry) ...
            && isfield(entry,'use') ...
            && ~isempty(entry.use)
        use = string(entry.use);
    elseif isstruct(bas) ...
            && isfield(bas,'K_t') ...
            && k <= double(bas.K_t)
        use = "training";
    elseif isstruct(bas) ...
            && isfield(bas,'K_e')
        use = "evaluation";
    end

end
