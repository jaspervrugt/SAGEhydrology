function [ax,frmt] = print_SAGE(mdl,ax,prf,i,dirres,loss,net)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRINT_SAGE Print and update live SAGE diagnostics.
%
%  Formats iteration metrics and refreshes the training diagnostic figures.
%
% SYNOPSIS:
%   [ax,frmt] = print_SAGE(mdl,ax,prf,i,dirres,loss,net)
%
% INPUT ARGUMENTS:
%   mdl             assessment design settings
%    .mode           basin/period design code
%   ax              optional graphics handles; empty on first call
%   prf             current basin metrics and iteration histories
%   i               training iteration
%   dirres          results directory
%   loss            optional selected loss settings
%    .observed       names of observations in the objective
%   net             network settings, including parameter count
%
% OUTPUT ARGUMENTS:
%   ax              updated graphics handles
%   frmt            format string used for iteration printing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6 ...
            || isempty(loss)
        loss = struct('fnc',4);
        loss_fnc = 4;
    else
        loss_fnc = loss.fnc;
    end
    if nargin < 7 ...
            || isempty(net) ...
            || ~isfield(net,'l') ...
            || isempty(net.l)
        l = nan;
    else
        l = net.l;
    end
    
    if ~isfield(mdl,'mode') ...
            || isempty(mdl.mode) ...
            || ~isscalar(mdl.mode) ...
            || ~isnumeric(mdl.mode) ...
            || ~ismember(double(mdl.mode),[1 2 3 4])
        error(['      Error: print_SAGE: ' ...
            'mdl.mode must be one of 1,2,3,4.']);
    end
    
    fdcFormulation = local_fdc_formulation(loss);
    if loss_fnc ~= 6
        fdcFormulation = 1;
    end
    [fdcSkill,fdcMetricName] = local_selected_fdc_skill( ...
        prf.curr,fdcFormulation);
    frmt = print_stats(mdl.mode,prf,i,l,loss_fnc,dirres, ...
        fdcFormulation,loss);
    ax = print_figs(mdl,prf,i,prf.curr.NSE,prf.curr.KGE, ...
        fdcSkill,prf.curr.JKGE,loss_fnc,ax, ...
        fdcFormulation,fdcMetricName,loss);

end

function formulation = local_fdc_formulation(loss)
%LOCAL_FDC_FORMULATION Backward-compatible 6a/6b/6c selector.

    formulation = 1;
    if isstruct(loss) && isfield(loss,'fdc') ...
            && isstruct(loss.fdc) ...
            && isfield(loss.fdc,'formulation') ...
            && ~isempty(loss.fdc.formulation)
        formulation = double(loss.fdc.formulation);
    end
    if ~isscalar(formulation) ...
            || ~isfinite(formulation) ...
            || ~ismember(formulation,1:3)
        formulation = 1;
    end
end

function [metric,name] = local_selected_fdc_skill(curr,formulation)
%LOCAL_SELECTED_FDC_SKILL Select the score paired with loss 6a, 6b, or 6c.

    [name,~,~,~] = local_fdc_print_names(formulation);
    if isstruct(curr) ...
            && isfield(curr,name)
        metric = curr.(name);
    elseif isstruct(curr) ...
            && isfield(curr,'S_fdc')
        metric = curr.S_fdc;
        name = 'S_fdc';
    else
        metric = struct('tt',[],'te',[],'et',[],'ee',[]);
    end
end

function [metric,medianField,sibField,lossName,lossTex] = ...
    local_fdc_print_names(formulation)
%LOCAL_FDC_PRINT_NAMES Storage names and compact printed notation.

    switch double(formulation)
        case 2
            metric = 'S_p';
            medianField = 'mS_p';
            sibField = 'Sib_S_p';
            lossName = 'd_p';
            lossTex = 'd_p';
        case 3
            metric = 'S_logp';
            medianField = 'mS_logp';
            sibField = 'Sib_S_logp';
            lossName = 'd_logp';
            lossTex = 'd_{\log p}';
        otherwise
            metric = 'S_fdc';
            medianField = 'mS_fdc';
            sibField = 'Sib_S_fdc';
            lossName = 'd_fdc';
            lossTex = 'd_{\mathrm{fdc}}';
    end
end

function tex = local_fdc_loss_scenario_tex(formulation,scenario)
%LOCAL_FDC_LOSS_SCENARIO_TEX Combine formulation and scenario subscripts.

    scn = lower(strtrim(char(string(scenario))));
    switch double(formulation)
        case 2
            tex = sprintf('d_{p,\\mathrm{%s}}',scn);
        case 3
            tex = sprintf('d_{\\log p,\\mathrm{%s}}',scn);
        otherwise
            tex = sprintf('d_{\\mathrm{fdc},\\mathrm{%s}}',scn);
    end
end

function tex = local_metric_scenario_tex(metricName,scenario)
%LOCAL_METRIC_SCENARIO_TEX Publication notation with one valid subscript.

    metric = lower(strtrim(char(string(metricName))));
    scn = lower(strtrim(char(string(scenario))));
    switch metric
        case 's_fdc'
            tex = sprintf('S_{\\mathrm{fdc},\\mathrm{%s}}',scn);
        case 's_p'
            tex = sprintf('S_{p,\\mathrm{%s}}',scn);
        case 's_logp'
            tex = sprintf('S_{\\log p,\\mathrm{%s}}',scn);
        case 'nse'
            tex = sprintf('\\mathrm{NSE}_{\\mathrm{%s}}',scn);
        case 'kge'
            tex = sprintf('\\mathrm{KGE}_{\\mathrm{%s}}',scn);
        case 'jkge'
            tex = sprintf('\\mathrm{JKGE}_{\\mathrm{%s}}',scn);
        otherwise
            tex = sprintf('\\mathrm{%s}_{\\mathrm{%s}}', ...
                char(string(metricName)),scn);
    end
end

function tex = local_median_metric_subscript(metricName)
%LOCAL_MEDIAN_METRIC_SUBSCRIPT Lowercase subscripts of the T estimator.

    metric = lower(strtrim(char(string(metricName))));
    switch metric
        case 'nse'
            tex = '\mathrm{nse}';
        case 'kge'
            tex = '\mathrm{kge}';
        case 'jkge'
            tex = '\mathrm{jkge}';
        case 's_fdc'
            tex = 'S_{\mathrm{fdc}}';
        case 's_p'
            tex = 'S_{p}';
        case 's_logp'
            tex = 'S_{\log p}';
        otherwise
            tex = sprintf('\\mathrm{%s}',metric);
    end
end

function update_ecdf_metric_labels(ax,useJKGE,fdcMetricName, ...
    rightField,rightAvailable,fontSize)
%UPDATE_ECDF_METRIC_LABELS Keep ECDF labels synchronized with loss 6a/b/c.

    scenarioTrain = 'tt';
    if rightAvailable
        scenarioRight = lower(strtrim(char(string(rightField))));
    else
        scenarioRight = 'na';
    end
    if useJKGE
        bottomMetric = 'JKGE';
    else
        bottomMetric = fdcMetricName;
    end

    local_set_ecdf_ylabel(ax.nse_t,'NSE',scenarioTrain,fontSize);
    local_set_ecdf_ylabel(ax.kge_t,'KGE',scenarioTrain,fontSize);
    local_set_ecdf_ylabel(ax.jkge_t,bottomMetric,scenarioTrain,fontSize);
    local_set_ecdf_xlabel(ax.jkge_t,bottomMetric,scenarioTrain,fontSize);
    xlabel(ax.nse_t,'');
    xlabel(ax.kge_t,'');

    local_set_ecdf_ylabel(ax.nse_r,'NSE',scenarioRight,fontSize);
    local_set_ecdf_ylabel(ax.kge_r,'KGE',scenarioRight,fontSize);
    local_set_ecdf_ylabel(ax.jkge_r,bottomMetric,scenarioRight,fontSize);
    local_set_ecdf_xlabel(ax.jkge_r,bottomMetric,scenarioRight,fontSize);
    xlabel(ax.nse_r,'');
    xlabel(ax.kge_r,'');
end

function local_set_ecdf_ylabel(axh,metricName,scenario,fontSize)
%LOCAL_SET_ECDF_YLABEL Set the scenario-aware ECDF y-axis label.

    symbol = local_metric_scenario_tex(metricName,scenario);
    ylabel(axh,sprintf('$F(%s)$',symbol), ...
        'Interpreter','latex','FontName', ...
        get(groot,'DefaultAxesFontName'), ...
        'FontSize',fontSize);
end

function local_set_ecdf_xlabel(axh,metricName,scenario,fontSize)
%LOCAL_SET_ECDF_XLABEL Set the scenario-aware ECDF x-axis label.

    symbol = local_metric_scenario_tex(metricName,scenario);
    xlabel(axh,sprintf('$%s$',symbol), ...
        'Interpreter','latex','FontName', ...
        get(groot,'DefaultAxesFontName'), ...
        'FontSize',fontSize);
end

function frmt = print_stats(mode,prf,i,l,loss_fnc,dirres, ...
    fdcFormulation,loss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRINT_STATS Prints iteration statistics to screen and log file
%
% SYNOPSIS:
%  frmt = print_stats(mode,prf,i,l,loss_fnc,dirres, ...
%      fdcFormulation,loss)
%
%   mode        assessment design
%                1 = training basins only | training period only
%                2 = training basins only | training and evaluation
%                    period/mask
%                3 = training and evaluation basins | training period only
%                4 = training and evaluation basins | training and
%                    evaluation period/mask
%   prf         structure with performance histories
%   i           iteration number
%   l           number of network weights and biases
%   loss_fnc    loss function (scalar, optional if stored in prf)
%   dirres      directory with SAGE results
%   fdcFormulation OPTIONAL duration-curve formulation, 1, 2, or 3
%   loss        OPTIONAL complete loss settings for joint-loss reporting
%
% OUTPUT:
%   frmt        print format string retained for backward compatibility
%
% NOTES:
%   - Statistics are printed in compact scenario-row table format.
%   - Depending on mode, scenarios tt, te, et, and/or ee are reported.
%   - Reported quantities include selected loss, RSS, median NSE,
%     median KGE or JKGE, and S_ib.
%   - The header is reprinted periodically and RSS is scaled adaptively.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Apr. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6 ...
            || isempty(dirres)
        dirres = pwd;
    end
    if ~exist(dirres,'dir')
        mkdir(dirres);
    end
    
    if nargin < 5 ...
            || isempty(loss_fnc)
        if isfield(prf,'loss_fnc') ...
                && ~isempty(prf.loss_fnc)
            loss_fnc = prf.loss_fnc;
        else
            error(['      Error: print_stats: ' ...
                'loss_fnc must be provided, ' ...
                   'or stored in prf.loss_fnc.']);
        end
    end
    
    log_file = fullfile(dirres, ...
        'demo_SAGE_iterations.txt');
    
    if i == 1
        fid = fopen(log_file,'w');
    else
        fid = fopen(log_file,'a');
    end
    
    if fid == -1
        error(['      Error: print_stats: ' ...
            'Could not open log file.']);
    end

    if nargin >= 8 ...
            && local_is_joint_loss(loss) ...
            && isfield(prf.iter,'joint')
        frmt = local_print_joint_stats( ...
            mode,prf,i,l,loss,fdcFormulation,fid);
        fclose(fid);
        return
    end
    
    % ------------
    % User options
    % ------------
    headerEvery = 25;   % Reprint header every N iterations
    rssDigits = 3;      % Digits after decimal for scaled RSS
    lossDigits = 3;
    if nargin < 7
        fdcFormulation = 1;
    end
    [fdcMetricName,fdcMedianField,fdcSibField,fdcLossName] = ...
        local_fdc_print_names(fdcFormulation);
    
    % -------------------------------------------------
    % Selected loss function label and associated field
    % -------------------------------------------------
    switch loss_fnc
        case 1
            lossName = 'ΣSAR';
            lossField = 'SAR';
            useRSS = true;
    
        case 2
            lossName = 'ΣGLS';
            lossField = 'GLS';
            useRSS = false;   % GLS already is weighted RSS-like quantity
    
        case 3
            lossName = 'Σ1-NSE';
            lossField = 'L';
            useRSS = true;
    
        case 4
            lossName = 'Σ1-KGE';
            lossField = 'L';
            useRSS = true;
    
        case 5
            lossName = 'ΣHuber';
            lossField = 'Huber';
            useRSS = true;
    
        case 6
            lossName = ['Σ' fdcLossName];
            lossField = 'L';
            useRSS = true;
    
        case 7
            lossName = 'Σ1-JKGE';
            lossField = 'L';
            useRSS = true;
    
        otherwise
            fclose(fid);
            error(['      Error: print_stats: ' ...
                'unknown loss_fnc = %g.'],loss_fnc);
    end
    
    K_t = prf.K_t;
    K_e = prf.K_e;
    
    % ------------------------------------------
    % Determine scenarios to print based on mode
    % ------------------------------------------
    switch mode
        case 1
            scens = {'tt'};
        case 2
            scens = {'tt','te'};
        case 3
            scens = {'tt','et'};
        case 4
            scens = {'tt','te','et','ee'};
        otherwise
            fclose(fid);
            error(['      Error: print_stats: ' ...
                'mode must be one of 1,2,3,4.']);
    end
    
    nScen = numel(scens);
    
    % --------------------------------
    % Collect values for all scenarios
    % --------------------------------
    cpuT = prf.iter.cpuT(i);
    
    Loss = nan(1,nScen);
    RSS = nan(1,nScen);
    mNSE = nan(1,nScen);
    mKGE = nan(1,nScen);
    mS_fdc = nan(1,nScen);
    mJKGE = nan(1,nScen);
    SibNSE = nan(1,nScen);
    SibKGE = nan(1,nScen);
    SibJKGE = nan(1,nScen);
    SibFDC = nan(1,nScen);
    
    for j = 1:nScen
        scn = scens{j};
    
        switch scn
            case {'tt','te'}
                Kscale = K_t;
            case {'et','ee'}
                Kscale = K_e;
            otherwise
                Kscale = 1;
        end
    
        if isfield(prf.iter.(lossField),scn) ...
                && i <= numel(prf.iter.(lossField).(scn))
            Loss(j) = Kscale * prf.iter.(lossField).(scn)(i);
        end
    
        if isfield(prf.iter.RSS,scn) ...
                && i <= numel(prf.iter.RSS.(scn))
            RSS(j) = Kscale * prf.iter.RSS.(scn)(i);
        end
    
        if isfield(prf.iter.mNSE,scn) ...
                && i <= numel(prf.iter.mNSE.(scn))
            mNSE(j) = prf.iter.mNSE.(scn)(i);
        end
    
        if isfield(prf.iter,'mKGE') ...
                && isfield(prf.iter.mKGE,scn) ...
                && i <= numel(prf.iter.mKGE.(scn))
            mKGE(j) = prf.iter.mKGE.(scn)(i);
        end

        if isfield(prf.iter,fdcMedianField) ...
                && isfield(prf.iter.(fdcMedianField),scn) ...
                && i <= numel(prf.iter.(fdcMedianField).(scn))
            mS_fdc(j) = prf.iter.(fdcMedianField).(scn)(i);
        end
        
        if isfield(prf.iter,'mJKGE') ...
                && isfield(prf.iter.mJKGE,scn) ...
                && i <= numel(prf.iter.mJKGE.(scn))
            mJKGE(j) = prf.iter.mJKGE.(scn)(i);
        end
    
        if isfield(prf.iter,'Sib_NSE') ...
                && isfield(prf.iter.Sib_NSE,scn) ...
                && i <= numel(prf.iter.Sib_NSE.(scn))
            SibNSE(j) = prf.iter.Sib_NSE.(scn)(i);
        end

        if isfield(prf.iter,'Sib_KGE') ...
                && isfield(prf.iter.Sib_KGE,scn) ...
                && i <= numel(prf.iter.Sib_KGE.(scn))
            SibKGE(j) = prf.iter.Sib_KGE.(scn)(i);
        end

        if isfield(prf.iter,'Sib_JKGE') ...
                && isfield(prf.iter.Sib_JKGE,scn) ...
                && i <= numel(prf.iter.Sib_JKGE.(scn))
            SibJKGE(j) = prf.iter.Sib_JKGE.(scn)(i);
        end

        if isfield(prf.iter,fdcSibField) ...
                && isfield(prf.iter.(fdcSibField),scn) ...
                && i <= numel(prf.iter.(fdcSibField).(scn))
            SibFDC(j) = prf.iter.(fdcSibField).(scn)(i);
        end
    end
    
    % JKGE is an on-demand metric. Show both JKGE columns only when it is
    % the selected loss and a finite score is actually available.
    useJKGE = loss_fnc == 7 ...
        && any(isfinite(mJKGE));

    % Use the ordinary four-displayed-digit convention for T and compact
    % notation for each integrated score. Large integrated scores switch
    % to one-decimal scientific notation so they cannot overflow a field.
    SibNSETxt = arrayfun(@local_format_sib,SibNSE, ...
        'UniformOutput',false);
    SibKGETxt = arrayfun(@local_format_sib,SibKGE, ...
        'UniformOutput',false);
    SibJKGETxt = arrayfun(@local_format_sib,SibJKGE, ...
        'UniformOutput',false);
    SibFDCTxt = arrayfun(@local_format_sib,SibFDC, ...
        'UniformOutput',false);
    
    % ---------------------------------
    % Adaptive RSS scale
    % Recompute when header is printed.
    % ---------------------------------
    printHeader = (i == 1) ...
        || (mod(i-1,headerEvery) == 0);
    
    persistent rssPow_cached rssScale_cached
    if isempty(rssPow_cached)
        rssPow_cached = 0;
        rssScale_cached = 1;
    end
    
    if printHeader
        rssVals = RSS(isfinite(RSS) & RSS ~= 0);
    
        if isempty(rssVals)
            rssPow_cached = 0;
            rssScale_cached = 1;
        else
            [rssPow_cached,rssScale_cached] = ...
                nice_eng_scale(rssVals);
        end
    end
    
    RSSs = RSS ./ rssScale_cached;
    
    persistent lossPow_cached lossScale_cached
    if isempty(lossPow_cached)
        lossPow_cached = 0;
        lossScale_cached = 1;
    end
    
    if printHeader
        lossVals = Loss(isfinite(Loss) & Loss ~= 0);
    
        if isempty(lossVals)
            lossPow_cached = 0;
            lossScale_cached = 1;
        else
            [lossPow_cached,lossScale_cached] = ...
                nice_eng_scale(lossVals);
        end
    end
    
    Losss = Loss ./ lossScale_cached;
    
    % ------------------------------------------------------------
    % Column definitions. Metric names use a compact two-level header:
    % scale factors sit above the loss/RSS columns, while T and S_ib are
    % centered over their respective metric groups.
    % ------------------------------------------------------------
    headers = {'scen',lossName};
    widths = [6 7];
    if useRSS
        % The diagnostic is the unweighted GLS/RSS sum. Name it by the
        % quantity displayed in SAGE rather than its storage field.
        headers{end + 1} = 'ΣGLS';
        widths(end + 1) = 7;
    end
    headers = [headers {'nse','kge'}];
    widths = [widths 7 7];
    if useJKGE
        headers{end + 1} = 'jkge';
        widths(end + 1) = 7;
    end
    headers{end + 1} = fdcMetricName;
    widths(end + 1) = 7;
    headers = [headers {'nse','kge'}];
    widths = [widths 7 7];
    if useJKGE
        headers{end + 1} = 'jkge';
        widths(end + 1) = 7;
    end
    headers{end + 1} = fdcMetricName;
    % Six characters are sufficient for the final integrated-score field
    % and keep the complete table at the requested 68-character width.
    widths(end + 1) = 6;
    % A single vertical rule after the optimized loss keeps the objective
    % visually distinct without boxing every individual metric.
    nMetric = 3 + double(useJKGE);
    sGroupStart = 2 + double(useRSS) + nMetric + 1;
    hdr1 = local_iteration_line(headers,widths,2,sGroupStart);
    hdrGroup = local_iteration_group_lines(widths,2 + double(useRSS), ...
        nMetric,nMetric,lossPow_cached,rssPow_cached,useRSS);

    % ------------------
    % Print header block
    % ------------------
    if printHeader

        % Match the separator to the exact rendered header width. hdr1
        % contains one trailing newline, which is excluded from the count.
        headerWidth = numel(hdr1) - 1;
        sep = repmat('-',1,headerWidth);
        hdr0 = sprintf(['\n%s\nit %4d   ' ...
            'l %5d   CPU %5.1f s\n'],sep,i,l,cpuT);
    
        fprintf('%s',hdr0);
        fprintf('%s',hdrGroup);
        fprintf('%s',hdr1);
        fprintf('%s\n',sep);
    
        fprintf(fid,'%s',hdr0);
        fprintf(fid,'%s',hdrGroup);
        fprintf(fid,'%s',hdr1);
        fprintf(fid,'%s\n',sep);
    else
        % Print compact iteration leader
        lead = sprintf(['\nit %4d   l %5d   ' ...
            'CPU %5.1f s\n'],i,l,cpuT);
        fprintf('%s',lead);
        fprintf(fid,'%s',lead);
    end
    
    % -------------------
    % Print scenario rows
    % -------------------
    for j = 1:nScen
        scn = scens{j};
        values = cell(1,numel(headers));
        c = 1;
        values{c} = scn;
        c = c + 1;
        values{c} = sprintf('%.*f',lossDigits,Losss(j));
        if useRSS
            c = c + 1;
            values{c} = sprintf('%.*f',rssDigits,RSSs(j));
        end
        c = c + 1;
        values{c} = local_format_metric(mNSE(j));
        c = c + 1;
        values{c} = local_format_metric(mKGE(j));
        if useJKGE
            c = c + 1;
            values{c} = local_format_metric(mJKGE(j));
        end
        c = c + 1;
        values{c} = local_format_metric(mS_fdc(j));
        c = c + 1;
        values{c} = SibNSETxt{j};
        c = c + 1;
        values{c} = SibKGETxt{j};
        if useJKGE
            c = c + 1;
            values{c} = SibJKGETxt{j};
        end
        c = c + 1;
        values{c} = SibFDCTxt{j};
        row = local_iteration_line(values,widths,2,sGroupStart);
    
        fprintf('%s',row);
        fprintf(fid,'%s',row);
    end
    
    frmt = '';   % retained only for backward compatibility
    fclose(fid);
end

function tf = local_is_joint_loss(loss)
%LOCAL_IS_JOINT_LOSS True for any nonlegacy named observation objective.

    tf = false;
    if isstruct(loss) ...
            && isfield(loss,'observed') ...
            && ~isempty(loss.observed)
        names = unique(upper(strtrim( ...
            string(loss.observed(:)))),'stable');
        names = names(strlength(names) > 0);
        tf = ~(isscalar(names) && names(1)=="Q");
    end
end

function frmt = local_print_joint_stats( ...
    mode,prf,i,l,loss,fdcFormulation,fid)
%LOCAL_PRINT_JOINT_STATS Print joint loss and per-observation diagnostics.

    switch mode
        case 1
            scens = {'tt'};
        case 2
            scens = {'tt','te'};
        case 3
            scens = {'tt','et'};
        case 4
            scens = {'tt','te','et','ee'};
        otherwise
            error('print_SAGE:BadMode', ...
                'Unknown assessment mode: %g.',mode);
    end

    names = unique(upper(strtrim( ...
        string(loss.observed(:)))),'stable');
    names = names(strlength(names) > 0);
    if ~isequal(names(:),["Q";"SWE"])
        frmt = local_print_named_stats( ...
            scens,prf,i,l,loss,fdcFormulation,fid,names);
        return
    end

    durationQ = 'S_fdc';
    durationSWE = 'S_sdc';
    if fdcFormulation == 2
        durationQ = 'S_p';
        durationSWE = 'S_p';
    elseif fdcFormulation == 3
        durationQ = 'S_logp';
        durationSWE = 'S_logp';
    end

    cpuT = prf.iter.cpuT(i);
    printHeader = i == 1 || mod(i-1,25) == 0;
    sep = repmat('-',1,105);

    J = prf.iter.joint;
    values = nan(numel(scens),9);
    for j = 1:numel(scens)
        sc = scens{j};
        values(j,:) = [J.loss.total.(sc)(i), ...
            J.loss.Q.(sc)(i), ...
            J.loss.SWE.(sc)(i), ...
            J.Q.NSE.(sc)(i),J.Q.KGE.(sc)(i), ...
            J.Q.duration.(sc)(i), ...
            J.SWE.NSE.(sc)(i),J.SWE.KGE.(sc)(i), ...
            J.SWE.duration.(sc)(i)];
    end

    persistent qPow qScale swePow sweScale
    if isempty(qPow)
        [qPow,qScale,swePow,sweScale] = deal(0,1,0,1);
    end
    if printHeader
        [qPow,qScale] = local_joint_scale(values(:,2));
        [swePow,sweScale] = local_joint_scale(values(:,3));
    end
    scaledQ = values(:,2)/qScale;
    scaledSWE = values(:,3)/sweScale;

    if printHeader
        lead = sprintf(['\n%s\nit %4d   l %5d   ' ...
            'CPU %5.1f s\n'],sep,i,l,cpuT);
        [cQ,cSWE,wQ,wSWE,sQ,sSWE] = ...
            local_joint_coefficients(loss);
        lossName = local_joint_loss_name(loss.fnc,fdcFormulation);
        names = upper(strtrim(string(loss.observed(:))));
        names = names(strlength(names) > 0);
        dataLine = sprintf(['data types: {%s}, loss functions: ' ...
            '%s (loss = %d)\n'], ...
            strjoin(cellstr(names),','),lossName,loss.fnc);
        formula = sprintf(['total loss: L_tot = c_Q x %s^Q + ' ...
            'c_SWE x %s^SWE\n'],lossName,lossName);
        whereQ = sprintf(['where: c_Q = w_Q/s_Q = %.4g/%.7g ' ...
            '= %.4g\n'],wQ,sQ,cQ);
        whereSWE = sprintf(['       c_SWE = w_SWE/s_SWE = ' ...
            '%.4g/%.4g = %.4g\n\n'],wSWE,sSWE,cSWE);
        group = sprintf(['            L_tot  |  %s^Q  %s^SWE  ' ...
            '|---------- Q --------|  |-------- SWE --------|\n'], ...
            lossName,lossName);
        header = sprintf(['scen               |   x10^%d   x10^%d     ' ...
            'nse     kge   %6s     nse     kge   %6s\n'], ...
            qPow,swePow,durationQ,durationSWE);
        tableSep = repmat('-',1,numel(header));
        block = [lead dataLine formula whereQ whereSWE ...
            group header tableSep newline];
        fprintf('%s',block);
        fprintf(fid,'%s',block);
    else
        lead = sprintf(['\nit %4d   l %5d   ' ...
            'CPU %5.1f s\n'],i,l,cpuT);
        fprintf('%s',lead);
        fprintf(fid,'%s',lead);
    end

    for j = 1:numel(scens)
        sc = scens{j};
        totalText = local_joint_total(values(j,1));
        row = sprintf(['%-6s %10s  | %6.3f   %6.3f   ' ...
            '%7.3f %7.3f %7.3f  ' ...
            '%7.3f %7.3f %7.3f\n'], ...
            sc,totalText,scaledQ(j),scaledSWE(j), ...
            values(j,4:9));
        fprintf('%s',row);
        fprintf(fid,'%s',row);
    end
    frmt = '';
end

function frmt = local_print_named_stats( ...
    scens,prf,i,l,loss,fdcFormulation,fid,names)
%LOCAL_PRINT_NAMED_STATS Print an arbitrary named-observation objective.

    J = prf.iter.joint;
    cpuT = prf.iter.cpuT(i);
    printHeader = i == 1 || mod(i-1,25) == 0;
    lossName = local_joint_loss_name(loss.fnc,fdcFormulation);
    if printHeader
        lead = sprintf(['\n%s\nit %4d   l %5d   ' ...
            'CPU %5.1f s\n'],repmat('-',1,105),i,l,cpuT);
        dataLine = sprintf(['data types: {%s}, loss function: ' ...
            '%s (loss = %d)\n'],strjoin(cellstr(names),','), ...
            lossName,loss.fnc);
        formula = 'total loss: L_tot = ';
        for k = 1:numel(names)
            name = char(names(k));
            [coefficient,weight,scale] = ...
                local_named_coefficient(loss,name,numel(names));
            if k > 1
                formula = [formula ' + ']; %#ok<AGROW>
            end
            formula = [formula sprintf('c_%s x L_%s', ...
                name,name)]; %#ok<AGROW>
            detail = sprintf(['  c_%s = %.4g ' ...
                '(w = %.4g, scale = %.7g)\n'], ...
                name,coefficient,weight,scale);
            dataLine = [dataLine detail]; %#ok<AGROW>
        end
        formula = [formula newline];
        header = sprintf('%-6s %11s','scen','L_tot');
        for k = 1:numel(names)
            name = char(names(k));
            header = [header sprintf([' | %11s %8s ' ...
                '%8s %8s'],['L_' name],['NSE_' name], ...
                ['KGE_' name],['JKGE_' name])]; %#ok<AGROW>
        end
        block = [lead dataLine formula header newline ...
            repmat('-',1,numel(header)) newline];
        fprintf('%s',block);
        fprintf(fid,'%s',block);
    else
        lead = sprintf(['\nit %4d   l %5d   ' ...
            'CPU %5.1f s\n'],i,l,cpuT);
        fprintf('%s',lead);
        fprintf(fid,'%s',lead);
    end

    for j = 1:numel(scens)
        sc = scens{j};
        row = sprintf('%-6s %11s',sc,local_joint_total( ...
            J.loss.total.(sc)(i)));
        for k = 1:numel(names)
            name = char(names(k));
            row = [row sprintf(' | %11.4g %8.3f %8.3f %8.3f', ...
                J.loss.(name).(sc)(i), ...
                J.(name).NSE.(sc)(i), ...
                J.(name).KGE.(sc)(i), ...
                J.(name).JKGE.(sc)(i))]; %#ok<AGROW>
        end
        row = [row newline]; %#ok<AGROW>
        fprintf('%s',row);
        fprintf(fid,'%s',row);
    end
    frmt = '';
end

function [coefficient,weight,scale] = ...
    local_named_coefficient(loss,name,namesCount)
%LOCAL_NAMED_COEFFICIENT Return one effective objective coefficient.

    weight = local_named_value(loss,'weight',name,1/namesCount);
    if isfield(loss,'normalization') ...
            && isstruct(loss.normalization)
        norm = loss.normalization;
    else
        norm = struct();
    end
    scale = local_named_value(norm,'scale',name,1);
    coefficient = weight/scale;
    if isfield(norm,'n_total') && isfield(norm,'n_basin')
        nTotal = double(norm.n_total);
        nBasin = local_named_value(norm,'n_basin',name,nTotal);
        coefficient = coefficient*nTotal/nBasin;
    end
end

function text = local_joint_total(value)
%LOCAL_JOINT_TOTAL Print five significant digits, including trailing zeros.

    if ~isfinite(value)
        text = char(string(value));
        return
    end
    if value == 0
        decimals = 4;
    else
        decimals = max(0,4-floor(log10(abs(value))));
    end
    text = sprintf('%.*f',decimals,value);
end

function [power,scale] = local_joint_scale(values)
%LOCAL_JOINT_SCALE Select a compact decimal scale for one loss column.

    values = values(isfinite(values) & values ~= 0);
    if isempty(values)
        power = 0;
        scale = 1;
    else
        power = floor(log10(max(abs(values))));
        scale = 10^power;
    end
end

function name = local_joint_loss_name(lossFnc,fdcFormulation)
%LOCAL_JOINT_LOSS_NAME Compact name of the shared component loss.

    switch double(lossFnc)
        case 1
            name = 'ΣSAR';
        case 2
            name = 'ΣGLS';
        case 3
            name = 'Σ1-NSE';
        case 4
            name = 'Σ1-KGE';
        case 5
            name = 'ΣHuber';
        case 6
            [~,~,~,fdcName] = ...
                local_fdc_print_names(fdcFormulation);
            name = ['Σ' fdcName];
        case 7
            name = 'sum(1-JKGE)';
        otherwise
            name = sprintf('loss%d',lossFnc);
    end
end

function [cQ,cSWE,wQ,wSWE,sQ,sSWE] = ...
    local_joint_coefficients(loss)
%LOCAL_JOINT_COEFFICIENTS Return effective printed Q and SWE coefficients.

    wQ = local_named_value(loss,'weight','Q',0.5);
    wSWE = local_named_value(loss,'weight','SWE',0.5);
    if isfield(loss,'normalization') ...
            && isstruct(loss.normalization)
        norm = loss.normalization;
    else
        norm = struct();
    end
    sQ = local_named_value(norm,'scale','Q',1);
    sSWE = local_named_value(norm,'scale','SWE',1);
    cQ = wQ/sQ;
    cSWE = wSWE/sSWE;
    if isfield(norm,'n_total') ...
            && isfield(norm,'n_basin')
        nTotal = double(norm.n_total);
        nQ = local_named_value(norm,'n_basin','Q',nTotal);
        nSWE = local_named_value(norm,'n_basin','SWE',nTotal);
        cQ = cQ*nTotal/nQ;
        cSWE = cSWE*nTotal/nSWE;
    end
end

function value = local_named_value(S,field,name,defaultValue)
%LOCAL_NAMED_VALUE Read one scalar from a named or numeric setting.

    value = defaultValue;
    if ~isstruct(S) ...
            || ~isfield(S,field) ...
            || isempty(S.(field))
        return
    end
    item = S.(field);
    if isstruct(item) ...
            && isfield(item,name)
        value = double(item.(name));
    elseif isnumeric(item) ...
            && isscalar(item)
        value = double(item);
    end
end

function ax = print_joint_figs(mdl,prf,i,loss,ax)
%PRINT_JOINT_FIGS Plot total and two observation-specific loss histories.

    names = unique(upper(strtrim( ...
        string(loss.observed(:)))),'stable');
    names = names(strlength(names) > 0);
    if numel(names) ~= 2
        ax = print_multi_loss_figs(mdl,prf,i,names,ax);
        return
    end

    fnt_ax = 17;
    fnt_lab = 18;
    fnt_tit = 19;
    fnt_med = 16;
    fnt_leg = 18;
    cT = [0 0 1];
    cE = [0.00 0.55 0.00];
    cAx = [0.15 0.15 0.15];
    c = local_sage_model_color(mdl.model);

    [rightField,nameE,~,rightAvailable] = ...
        choose_joint_compare_scenario(prf,names,'period');
    needInit = isempty(ax) || ~isstruct(ax) ...
        || ~isfield(ax,'fig') || ~isgraphics(ax.fig) ...
        || ~isfield(ax,'layout') ...
        || ~strcmp(ax.layout,'joint');
    if needInit && isstruct(ax) ...
            && isfield(ax,'fig') && isgraphics(ax.fig)
        delete(ax.fig);
    end

    if needInit
        ax = struct('layout','joint');
        scr = get(0,'ScreenSize');
        figW = 0.95*scr(3);
        figH = 0.95*scr(4);
        figX = scr(1)+0.5*(scr(3)-figW);
        figY = scr(2)+0.5*(scr(4)-figH);
        modelName = upper(local_print_model_name(mdl));
        ax.fig = figure('Units','pixels','Color','w', ...
            'Name',sprintf('%s: SAGE diagnostics',modelName), ...
            'NumberTitle','off', ...
            'Position',[figX figY figW figH], ...
            'SizeChangedFcn', ...
            @(src,evt) refresh_history_top_frames(src));
        setappdata(ax.fig,'SAGEPreserveCaptureSize',true);
        ax.flag = add_country_flag(ax.fig,mdl);

        marginL = 0.05;
        marginR = 0.05;
        marginT = 0.04;
        marginB = 0.115;
        gapX1 = 0.08;
        gapX2 = 0.08;
        gapY = 0.04;
        totalH = 1-marginT-marginB;
        rowH = (totalH-2*gapY)/3;
        colW_ecdf = 0.175;
        xLcol = marginL;
        xRcol = 1-marginR-colW_ecdf;
        xMcol = xLcol+colW_ecdf+gapX1;
        colW_mid = 0.92*(xRcol-gapX2-xMcol);
        y3 = marginB;
        y2 = y3+rowH+gapY;
        y1 = y2+rowH+gapY;

        ax.blank_t = axes('Parent',ax.fig, ...
            'Position',[xLcol y1 colW_ecdf rowH], ...
            'Visible','off');
        ax.blank_r = axes('Parent',ax.fig, ...
            'Position',[xRcol y1 colW_ecdf rowH], ...
            'Visible','off');
        ax.totalAx = axes('Parent',ax.fig, ...
            'Position',[xMcol y1 colW_mid rowH]);
        ax.obs1_t = axes('Parent',ax.fig, ...
            'Position',[xLcol y2 colW_ecdf rowH]);
        ax.obs1Ax = axes('Parent',ax.fig, ...
            'Position',[xMcol y2 colW_mid rowH]);
        ax.obs1_r = axes('Parent',ax.fig, ...
            'Position',[xRcol y2 colW_ecdf rowH]);
        ax.obs2_t = axes('Parent',ax.fig, ...
            'Position',[xLcol y3 colW_ecdf rowH]);
        ax.obs2Ax = axes('Parent',ax.fig, ...
            'Position',[xMcol y3 colW_mid rowH]);
        ax.obs2_r = axes('Parent',ax.fig, ...
            'Position',[xRcol y3 colW_ecdf rowH]);

        historyAxes = [ax.totalAx ax.obs1Ax ax.obs2Ax];
        ecdfAxes = [ax.obs1_t ax.obs1_r ax.obs2_t ax.obs2_r];
        for ah = [historyAxes ecdfAxes]
            hold(ah,'on');
        end
        set(ecdfAxes,'FontSize',fnt_ax,'LineWidth',1, ...
            'TickDir','out','Box','off','Layer','top', ...
            'XLim',[-1 1],'YLim',[0 1]);
        set([ax.obs1_t ax.obs2_t],'XColor',cT,'YColor',cT);
        set([ax.obs1_r ax.obs2_r],'XColor',cE,'YColor',cE);
        for ah = [ax.obs1_t ax.obs2_t]
            add_top_right_frame(ah,-1,1,0,1,cT);
        end
        for ah = [ax.obs1_r ax.obs2_r]
            add_top_right_frame(ah,-1,1,0,1,cE);
        end
        set([ax.obs1_t ax.obs1_r],'XTickLabel',[]);
        for ah = historyAxes
            init_history_axis(ah,fnt_ax);
        end
        setappdata(ax.fig,'HistoryAxesHandles',historyAxes);
        setappdata(ax.fig,'HistoryTopFrameColor',cAx);

        ax.totalLineT = local_joint_history_line( ...
            ax.totalAx,'left',cT,'-o',true);
        ax.totalLineR = local_joint_history_line( ...
            ax.totalAx,'right',cE,'-s',true);
        ax.obs1LineT = local_joint_history_line( ...
            ax.obs1Ax,'left',cT,'-o',false);
        ax.obs1LineR = local_joint_history_line( ...
            ax.obs1Ax,'right',cE,'-s',false);
        ax.obs2LineT = local_joint_history_line( ...
            ax.obs2Ax,'left',cT,'-o',false);
        ax.obs2LineR = local_joint_history_line( ...
            ax.obs2Ax,'right',cE,'-s',false);
        legend(ax.totalAx,'off');

        tags = {'obs1_t','obs1_r','obs2_t','obs2_r'};
        for k = 1:numel(tags)
            tag = tags{k};
            ax.ecdf.(tag) = struct('patch',gobjects(1), ...
                'line',gobjects(1),'med_v',gobjects(1), ...
                'med_h',gobjects(1),'med_s',gobjects(1), ...
                'med_txt',gobjects(1));
        end
    end

    title(ax.totalAx,'$\mathrm{History\ of\ total\ loss}$', ...
        'Interpreter','latex','FontWeight','normal', ...
        'FontSize',fnt_tit);
    title(ax.obs1Ax,sprintf( ...
        '$\\mathrm{History\\ of\\ %s\\ loss}$',names(1)), ...
        'Interpreter','latex','FontWeight','normal', ...
        'FontSize',fnt_tit);
    title(ax.obs2Ax,sprintf( ...
        '$\\mathrm{History\\ of\\ %s\\ loss}$',names(2)), ...
        'Interpreter','latex','FontWeight','normal', ...
        'FontSize',fnt_tit);
    set(ax.totalLineT,'DisplayName', ...
        'Train basins | train period');
    if rightAvailable
        set(ax.totalLineR,'DisplayName',nameE);
    else
        set(ax.totalLineR,'DisplayName','Unavailable');
    end
    xlabel(ax.totalAx,'');
    xlabel(ax.obs1Ax,'');
    xlabel(ax.obs2Ax,'SAGE iteration, $i$', ...
        'Interpreter','latex','FontSize',fnt_lab);

    for j = 1:2
        name = char(names(j));
        leftTag = sprintf('obs%d_t',j);
        rightTag = sprintf('obs%d_r',j);
        leftAx = ax.(leftTag);
        rightAx = ax.(rightTag);
        zT = local_joint_current_nse(prf,name,'tt');
        zR = [];
        if rightAvailable
            zR = local_joint_current_nse( ...
                prf,name,char(rightField));
        end
        [zT,~] = finite_vec(zT);
        [zR,~] = finite_vec(zR);
        title(leftAx,local_joint_ecdf_title( ...
            name,'Train basins | train period'), ...
            'Interpreter','latex','FontWeight','normal', ...
            'FontSize',fnt_tit);
        title(rightAx,local_joint_ecdf_title(name,nameE), ...
            'Interpreter','latex','FontWeight','normal', ...
            'FontSize',fnt_tit);
        ax = update_manual_ecdf(ax,leftTag,leftAx,zT, ...
            c,-1,1,0.20,1.4,'NSE','tt',fnt_med,name);
        if rightAvailable
            ax = update_manual_ecdf(ax,rightTag,rightAx,zR, ...
                c,-1,1,0.20,1.4,'NSE', ...
                char(rightField),fnt_med,name);
        else
            clear_manual_ecdf(ax,rightTag);
            set_unavailable_label(rightAx,true);
        end
        ylabel(leftAx,sprintf( ...
            '$F(\\mathrm{NSE}^{\\mathrm{%s}}_{\\mathrm{tt}})$',name), ...
            'Interpreter','latex','FontSize',fnt_lab);
        ylabel(rightAx,sprintf( ...
            '$F(\\mathrm{NSE}^{\\mathrm{%s}}_{\\mathrm{%s}})$',name, ...
            local_joint_scenario(rightField,rightAvailable)), ...
            'Interpreter','latex','FontSize',fnt_lab);
        if j == 2
            xlabel(leftAx,sprintf( ...
                '$\\mathrm{NSE}^{\\mathrm{%s}}_{\\mathrm{tt}}$',name), ...
                'Interpreter','latex','FontSize',fnt_lab);
            xlabel(rightAx,sprintf( ...
                '$\\mathrm{NSE}^{\\mathrm{%s}}_{\\mathrm{%s}}$',name, ...
                local_joint_scenario(rightField,rightAvailable)), ...
                'Interpreter','latex','FontSize',fnt_lab);
        end
    end

    [totalT,totalR] = local_joint_history_pair( ...
        prf,'total',i,rightField,rightAvailable);
    [obs1T,obs1R] = local_joint_history_pair( ...
        prf,char(names(1)),i,rightField,rightAvailable);
    [obs2T,obs2R] = local_joint_history_pair( ...
        prf,char(names(2)),i,rightField,rightAvailable);
    rightScn = local_joint_scenario(rightField,rightAvailable);
    update_history_side(ax.totalAx,'left',ax.totalLineT, ...
        totalT,cT,'$L_{\mathrm{tot}}$', ...
        fnt_lab,false,false);
    update_history_side(ax.totalAx,'right',ax.totalLineR, ...
        totalR,cE,'$L_{\mathrm{tot}}$', ...
        fnt_lab,false,true);
    obs1LabT = local_joint_loss_label(loss,names(1),'tt');
    obs1LabR = local_joint_loss_label(loss,names(1),rightScn);
    obs2LabT = local_joint_loss_label(loss,names(2),'tt');
    obs2LabR = local_joint_loss_label(loss,names(2),rightScn);
    update_dual_history_loss(ax.obs1Ax, ...
        ax.obs1LineT,ax.obs1LineR,obs1T,obs1R,cT,cE, ...
        obs1LabT,obs1LabR,fnt_lab);
    update_dual_history_loss(ax.obs2Ax, ...
        ax.obs2LineT,ax.obs2LineR,obs2T,obs2R,cT,cE, ...
        obs2LabT,obs2LabR,fnt_lab);

    xmax = max([2 numel(totalT) numel(totalR) ...
        numel(obs1T) numel(obs1R) numel(obs2T) numel(obs2R)]);
    xt = local_history_xticks(xmax);
    for ah = [ax.totalAx ax.obs1Ax ax.obs2Ax]
        xlim(ah,[1 xmax]);
    end
    set_history_xticks(ax.totalAx,xt,false);
    set_history_xticks(ax.obs1Ax,xt,false);
    set_history_xticks(ax.obs2Ax,xt,true);
    if rightAvailable
        legendRight = nameE;
    else
        legendRight = 'Unavailable';
    end
    local_top_history_legend(ax.totalAx, ...
        ax.totalLineT,ax.totalLineR,cT,cE, ...
        'Train basins | train period',legendRight,fnt_leg);
    if i == 1 || needInit
        for ah = [ax.totalAx ax.obs1Ax ax.obs2Ax]
            add_history_top_frame(ah,cAx);
        end
    end
    if ~isdeployed
        drawnow expose
    else
        drawnow limitrate nocallbacks
    end
end

function ax = print_multi_loss_figs(mdl,prf,i,names,ax)
%PRINT_MULTI_LOSS_FIGS Plot total and component losses for 3+ variables.

    cT = [0 0 1];
    cE = [0.00 0.55 0.00];
    [rightField,nameE,~,rightAvailable] = ...
        choose_joint_compare_scenario(prf,names,'period');
    needInit = isempty(ax) || ~isstruct(ax) ...
        || ~isfield(ax,'fig') || ~isgraphics(ax.fig) ...
        || ~isfield(ax,'layout') || ~strcmp(ax.layout,'multi');
    if needInit && isstruct(ax) ...
            && isfield(ax,'fig') && isgraphics(ax.fig)
        delete(ax.fig);
    end
    if needInit
        ax = struct('layout','multi');
        scr = get(0,'ScreenSize');
        ax.fig = figure('Units','pixels','Color','w', ...
            'Name',sprintf('%s: SAGE diagnostics', ...
            upper(local_print_model_name(mdl))), ...
            'NumberTitle','off', ...
            'Position',[scr(1)+0.05*scr(3) scr(2)+0.05*scr(4) ...
            0.9*scr(3) 0.9*scr(4)]);
        ax.totalAx = subplot(2,1,1,'Parent',ax.fig);
        ax.componentAx = subplot(2,1,2,'Parent',ax.fig);
        hold(ax.totalAx,'on');
        hold(ax.componentAx,'on');
        ax.totalLineT = local_joint_history_line( ...
            ax.totalAx,'left',cT,'-o',true);
        ax.totalLineR = local_joint_history_line( ...
            ax.totalAx,'right',cE,'-s',true);
    end

    [totalT,totalR] = local_joint_history_pair( ...
        prf,'total',i,rightField,rightAvailable);
    rightScn = local_joint_scenario(rightField,rightAvailable);
    update_history_side(ax.totalAx,'left',ax.totalLineT, ...
        totalT,cT,'$L_{\mathrm{tot},\mathrm{tt}}$',18,false,false);
    update_history_side(ax.totalAx,'right',ax.totalLineR, ...
        totalR,cE,sprintf('$L_{\mathrm{tot},\mathrm{%s}}$', ...
        rightScn),18,false,true);
    if rightAvailable
        legendRight = nameE;
    else
        legendRight = 'Unavailable';
    end
    local_top_history_legend(ax.totalAx, ...
        ax.totalLineT,ax.totalLineR,cT,cE, ...
        'Train basins | train period',legendRight,16);
    title(ax.totalAx,'History of total loss', ...
        'Interpreter','none','FontSize',19);

    cla(ax.componentAx);
    hold(ax.componentAx,'on');
    colors = lines(numel(names));
    for k = 1:numel(names)
        name = char(names(k));
        [trainLoss,evalLoss] = local_joint_history_pair( ...
            prf,name,i,rightField,rightAvailable);
        plot(ax.componentAx,1:numel(trainLoss),trainLoss, ...
            '-o','Color',colors(k,:),'LineWidth',1.3, ...
            'DisplayName',sprintf('%s tt',name));
        if rightAvailable
            plot(ax.componentAx,1:numel(evalLoss),evalLoss, ...
                '--s','Color',colors(k,:),'LineWidth',1.3, ...
                'DisplayName',sprintf('%s %s',name,nameE));
        end
    end
    xlabel(ax.componentAx,'SAGE iteration, $i$', ...
        'Interpreter','latex','FontSize',18);
    ylabel(ax.componentAx,'Individual loss', ...
        'Interpreter','none','FontSize',18);
    title(ax.componentAx,'Observation-specific losses', ...
        'Interpreter','none','FontSize',19);
    set(ax.componentAx,'FontSize',17,'LineWidth',1, ...
        'TickDir','out','Box','off');
    legend(ax.componentAx,'Location','best','Box','off');
    grid(ax.componentAx,'on');
    drawnow limitrate nocallbacks
end

function h = local_joint_history_line(axh,side,color,style,showLegend)
%LOCAL_JOINT_HISTORY_LINE Create one joint-objective history line.

    yyaxis(axh,side);
    if strcmp(side,'left')
        label = 'Train basins | train period';
    else
        label = 'Eval basins | eval period';
    end
    if showLegend
        visibility = 'on';
    else
        visibility = 'off';
    end
    h = plot(axh,nan,nan,style,'Color',color, ...
        'MarkerFaceColor',color,'LineWidth',1.3, ...
        'DisplayName',label,'HandleVisibility',visibility);
end

function local_top_history_legend( ...
    axh,hT,hE,cT,cE,labelT,labelE,fontSize)
%LOCAL_TOP_HISTORY_LEGEND Match the single-data top-panel legend.

    labelT = local_colored_legend_text(cT,labelT);
    labelE = local_colored_legend_text(cE,labelE);
    set([hT hE],'HandleVisibility','off');

    key = 'TopHistoryLegendHandles';
    if isappdata(axh,key)
        old = getappdata(axh,key);
        delete(old(isgraphics(old)));
    end

    yyaxis(axh,'left');
    pT = plot(axh,nan,nan,'-o','Color',cT, ...
        'MarkerFaceColor',cT,'LineWidth',2.3, ...
        'MarkerSize',8,'DisplayName',labelT);
    yyaxis(axh,'right');
    pE = plot(axh,nan,nan,'-s','Color',cE, ...
        'MarkerFaceColor',cE,'LineWidth',2.3, ...
        'MarkerSize',8,'DisplayName',labelE);
    setappdata(axh,key,[pT pE]);

    lgd = legend(axh,[pT pE], ...
        'Location','northeast','Interpreter','tex','Box','off');
    lgd.FontSize = fontSize;
end

function label = local_colored_legend_text(color,textValue)
%LOCAL_COLORED_LEGEND_TEXT Return a TeX-colored legend label.

    label = sprintf('\\color[rgb]{%.4f,%.4f,%.4f} %s', ...
        color(1),color(2),color(3),char(string(textValue)));
end

function z = local_joint_current_nse(prf,name,scenario)
%LOCAL_JOINT_CURRENT_NSE Return the current named NSE value.

    z = [];
    if isfield(prf,'curr') && isfield(prf.curr,'joint') ...
            && isfield(prf.curr.joint,name) ...
            && isfield(prf.curr.joint.(name),'NSE') ...
            && isfield(prf.curr.joint.(name).NSE,scenario)
        z = prf.curr.joint.(name).NSE.(scenario);
    end
end

function [yT,yR] = local_joint_history_pair( ...
    prf,name,i,rightField,rightAvailable)
%LOCAL_JOINT_HISTORY_PAIR Return left- and right-scenario histories.

    yT = [];
    yR = [];
    if ~isfield(prf,'iter') || ~isfield(prf.iter,'joint') ...
            || ~isfield(prf.iter.joint,'loss') ...
            || ~isfield(prf.iter.joint.loss,name)
        return
    end
    hist = prf.iter.joint.loss.(name);
    yT = local_joint_history_values(hist,'tt',i);
    if rightAvailable
        scn = char(rightField);
        yR = local_joint_history_values(hist,scn,i);
    end
end

function y = local_joint_history_values(hist,scenario,i)
%LOCAL_JOINT_HISTORY_VALUES Return finite history values for a scenario.

    y = [];
    if ~isfield(hist,scenario)
        return
    end
    n = min(i,numel(hist.(scenario)));
    id = 1:n;
    id = id(isfinite(hist.(scenario)(id)));
    y = hist.(scenario)(id);
end

function value = local_joint_ecdf_title(observable,label)
%LOCAL_JOINT_ECDF_TITLE Return an upright LaTeX ECDF-panel title.

    observable = upper(strtrim(char(string(observable))));
    parts = strtrim(strsplit(char(string(label)),'|'));
    parts = cellfun(@(x)strrep(x,' ','\ '), ...
        parts,'UniformOutput',false);
    if numel(parts) == 2
        value = sprintf('$\\mathrm{%s:\\ %s}\\mid\\mathrm{%s}$', ...
            observable,parts{1},parts{2});
    else
        value = sprintf('$\\mathrm{%s:\\ %s}$', ...
            observable,parts{1});
    end
end

function scenario = local_joint_scenario(rightField,rightAvailable)
%LOCAL_JOINT_SCENARIO Return the available right-side scenario tag.

    if rightAvailable
        scenario = char(rightField);
    else
        scenario = 'na';
    end
end

function ax = print_figs(mdl,prf,i,NSE,KGE,S_fdc,JKGE,loss_fnc,ax, ...
    fdcFormulation,fdcMetricName,loss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRINT_FIGS Live ECDFs of NSE, KGE and JKGE for training scenario and the
% mode-dependent comparison scenario, plus RSS and S_ib history plots.
%
% SYNOPSIS:
%   ax = print_figs(mdl,prf,i,NSE,KGE,S_fdc,JKGE,loss_fnc,ax, ...
%       fdcFormulation,fdcMetricName,loss)
%
%   mdl         structure with model state/parameter and assessment info
%    .mode       assessment design
%                 1 = train basins | train period/mask
%                 2 = train basins | train and eval period/mask
%                 3 = train and eval basins | train period/mask
%                 4 = train and eval basins | train and eval period/mask
%    .names      list of model names
%    .model      selected model index
%    .sp_method  split design used to label period/mask in titles
%   prf         structure with performance histories for explicit scenarios
%                tt = training basins | training period
%                te = training basins | evaluation period/mask
%                et = evaluation basins | training period
%                ee = evaluation basins | evaluation period/mask
%   i           iteration number
%   NSE         structure with scenario-wise NSE vectors
%                .tt, .te, .et, .ee
%   KGE         structure with scenario-wise KGE vectors
%                .tt, .te, .et, .ee
%   S_fdc       structure with scenario-wise S_fdc skill score vectors
%                .tt, .te, .et, .ee
%   JKGE        structure with scenario-wise JKGE vectors
%                .tt, .te, .et, .ee
%   loss_fnc    loss function (scalar)
%   ax          OPTIONAL: structure with graphics handles
%   fdcFormulation selected duration-curve formulation, 1, 2, or 3
%   fdcMetricName selected duration-curve skill field
%   loss        loss settings, including .observed for joint objectives
%
% NOTES:
%   Left ECDF panels always correspond to scenario tt.
%   Right ECDF panels correspond to:
%     mode 1 -> unavailable
%     mode 2 -> te
%     mode 3 -> et
%     mode 4 -> ee
%   For one observation, RSS and S_ib histories plot tt on the left axis
%   and the mode-dependent comparison scenario on the right axis. For two
%   observations, the top history shows total loss and the next two rows
%   show NSE distributions and individual losses for each observation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Mar. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 9 ...
            || isempty(ax)
        ax = struct();
        ax.layout = 'single';
    end
    if nargin >= 12 ...
            && local_is_joint_loss(loss)
        ax = print_joint_figs(mdl,prf,i,loss,ax);
        return
    end
    
    % ----------------
    % Basic formatting
    % ----------------
    fnt_ax = 17;    % tick labels
    fnt_lab = 18;   % axis labels
    fnt_tit = 19;
    fnt_med = 16;   % median annotations
    fnt_leg = 18;   % legends
    lw = 1.4;
    fa = 0.20;
    xL = -1; xR = 1;
    yL = 0;  yR = 1;

    cT = [0 0 1];
    %cE = [1 0.5 0];
    %cE = [0.75 0.00 0.00]; % crimson: right y-axis
    cE = [0.00 0.55 0.00]; 	% green
    cAx = [0.15 0.15 0.15]; % 2nd x-axis middle history panels
    
    % ----------------------
    % Model names and colors
    % ----------------------
    if ~isfield(mdl,'names') ...
            || isempty(mdl.names)
        error(['      Error: print_figs_new: ' ...
            'mdl.names is missing or empty.']);
    end
    if ~isfield(mdl,'model') ...
            || isempty(mdl.model)
        error(['      Error: print_figs_new: ' ...
            'mdl.model is missing or empty.']);
    end
    
    model = mdl.model;
    mdl_name_raw = local_print_model_name(mdl);
    
    if strcmpi(mdl_name_raw,'xinanjiang')
        mdl_name_disp = 'Xinanjiang';
    else
        mdl_name_disp = upper(mdl_name_raw);
    end
    mdl_tex = strrep(mdl_name_disp,'_','\_');
    
    colors = [ ...
        0.0000 0.4470 0.7410;
        0.9500 0.4500 0.0000;
        0.9290 0.6940 0.1250;
        0.6000 0.2000 0.8000;
        0.0000 0.5500 0.0000;
        0.6500 0.3500 0.0000;
        0.0000 0.7000 0.7000;
        0.9000 0.0000 0.9000;
        0.6500 0.6500 0.6500];
    
    if model >= 1 ...
            && model <= size(colors,1)
        c = colors(model,:);
    else
        c = [0 0 0];
    end
    
    % ---------------------------
    % Split method / period label
    % ---------------------------
    if isfield(mdl,'sp_method') ...
            && ~isempty(mdl.sp_method)
        sp_method = lower(string(mdl.sp_method));
    else
        sp_method = "manual";
    end
    
    if any(strcmp(sp_method, ...
            ["manual","traditional_block","block", ...
            "deterministic_block","random_block", ...
            "deterministic_kfold"]))
        samp_word = 'period';
    else
        samp_word = 'mask';
    end
    
    nameT = sprintf(['Train basins ' ...
        '| train %s'],samp_word);
    
    % ------------------------------------------
    % Comparison scenario priority: ee > et > te
    % ------------------------------------------
    [rightField,nameE,klabel_right,rightAvailable] = ...
        choose_compare_scenario(prf,samp_word);
    
    % ----------------
    % Pull metric data
    % ----------------
    NSEtt = get_metric_vector(NSE,'tt');
    KGEtt = get_metric_vector(KGE,'tt');
    if nargin < 10
        fdcFormulation = 1;
    end
    if nargin < 11 || isempty(fdcMetricName)
        [fdcMetricName,~,~,~] = local_fdc_print_names(fdcFormulation);
    end
    Sfdctt = get_metric_vector(S_fdc,'tt');
    JKGEtt = get_metric_vector(JKGE,'tt');
    
    if rightAvailable
        NSEr = get_metric_vector(NSE,char(rightField));
        KGEr = get_metric_vector(KGE,char(rightField));
        Sfdcr = get_metric_vector(S_fdc,char(rightField));
        JKGEr = get_metric_vector(JKGE,char(rightField));
    else
        NSEr = [];
        KGEr = [];
        Sfdcr = [];
        JKGEr = [];
    end
    
    [NSEtt_use,nt_nse] = finite_vec(NSEtt);
    [KGEtt_use,~] = finite_vec(KGEtt);
    [Sfdctt_use,~] = finite_vec(Sfdctt);
    [JKGEtt_use,~] = finite_vec(JKGEtt);
    
    [NSEr_use,nr_nse] = finite_vec(NSEr);
    [KGEr_use,~] = finite_vec(KGEr);
    [Sfdcr_use,~] = finite_vec(Sfdcr);
    [JKGEr_use,~] = finite_vec(JKGEr);
    
    % -------------------------------------------------
    % Bottom ECDF: JKGE when optimized, otherwise S_FDC
    % -------------------------------------------------
    useJKGEpanel = (loss_fnc == 7);

    if useJKGEpanel
        bottomT = JKGEtt_use;
        bottomR = JKGEr_use;
        bottomName = 'JKGE';
    else
        bottomT = Sfdctt_use;
        bottomR = Sfdcr_use;
        bottomName = fdcMetricName;
    end

    % ----------------------------------------------------
    % Horizontal limits for bottom JKGE / S_fdc ECDF panel
    % ----------------------------------------------------
    xLb = -1;
    xRb = 1;

    % -----------------------
    % Loss label / loss field
    % -----------------------
    [~,lossLeftLabel,lossField] = ...
        get_loss_strings(loss_fnc,fdcFormulation);
    
    K_t = prf.K_t;
    K_e = prf.K_e;
    K_right = K_e;
    if rightAvailable && strcmp(string(rightField),"te")
        K_right = K_t;
    end
    
    % -------------------
    % Need initialization
    % -------------------
    needInit = isempty(ax) ...
        || ~isstruct(ax) ...
        || ~isfield(ax,'fig') ...
        || ~isgraphics(ax.fig);
    
    if ~needInit
        req = {'nse_t','kge_t', ...
            'jkge_t','lossAx', ...
            'rssAx','sibAx', ...
            'nse_r','kge_r','jkge_r'};
        for k = 1:numel(req)
            if ~isfield(ax,req{k}) ...
                    || ~isgraphics(ax.(req{k}))
                needInit = true;
                break
            end
        end
    end
    
    if needInit
        ax = struct();
        scr = get(0,'ScreenSize');
        
        figW = 0.95 * scr(3);
        figH = 0.95 * scr(4);
        
        figX = scr(1) + 0.5 * (scr(3) - figW);
        figY = scr(2) + 0.5 * (scr(4) - figH);
        
        ax.fig = figure( ...
            'Units','pixels', ...
            'Color','w', ...
            'Name',sprintf('%s: SAGE diagnostics',mdl_name_disp), ...
            'NumberTitle','off', ...
            'Position',[figX figY figW figH], ...
            'SizeChangedFcn', ...
            @(src,evt) refresh_history_top_frames(src));
        setappdata(ax.fig,'SAGEPreserveCaptureSize',true);
        ax.flag = add_country_flag(ax.fig,mdl);
        % --------------------------------------------------------
        % Manual layout with more whitespace around history panels
        % --------------------------------------------------------
        marginL = 0.05;
        marginR = 0.05;
        marginT = 0.04;
        % Keep the bottom history tick labels and xlabel inside the figure
        % canvas during PNG/PPTX export.
        marginB = 0.115;
        gapX1 = 0.08;   % left gap: left ECDF -> history
        gapX2 = 0.08;   % right gap: history -> right ECDF
        gapY = 0.04;
    
        totalH = 1 - marginT - marginB;
        rowH = (totalH - 2*gapY) / 3;
    
        colW_ecdf = 0.175;
    
        % Fixed left and right ECDF columns
        xLcol = marginL;
        xRcol = 1 - marginR - colW_ecdf;
    
        % Fixed history start
        xMcol = xLcol + colW_ecdf + gapX1;
    
        % History width fills the space up to the fixed right ECDF column
        colW_mid0 = xRcol - gapX2 - xMcol;
        colW_mid = 0.92 * colW_mid0;   % reduce width by 15%
    
        y3 = marginB;
        y2 = y3 + rowH + gapY;
        y1 = y2 + rowH + gapY;
    
        % Left ECDF column
        ax.nse_t = axes('Parent',ax.fig, ...
            'Position',[xLcol y1 colW_ecdf rowH]); 
        hold(ax.nse_t,'on');
        ax.kge_t = axes('Parent',ax.fig, ...
            'Position',[xLcol y2 colW_ecdf rowH]); 
        hold(ax.kge_t,'on');
        ax.jkge_t = axes('Parent',ax.fig, ...
            'Position',[xLcol y3 colW_ecdf rowH]); 
        hold(ax.jkge_t,'on');
    
        % Middle history column
        ax.lossAx = axes('Parent',ax.fig, ...
            'Position',[xMcol y1 colW_mid rowH]); 
        hold(ax.lossAx,'on');
        ax.rssAx = axes('Parent',ax.fig, ...
            'Position',[xMcol y2 colW_mid rowH]); 
        hold(ax.rssAx,'on');
        ax.sibAx = axes('Parent',ax.fig, ...
            'Position',[xMcol y3 colW_mid rowH]); 
        hold(ax.sibAx,'on');
    
        setappdata(ax.fig, ...
            'HistoryAxesHandles', ...
            [ax.lossAx ax.rssAx ax.sibAx]);
        setappdata(ax.fig, ...
            'HistoryTopFrameColor',cAx);
    
        % Right ECDF column
        ax.nse_r = axes('Parent',ax.fig, ...
            'Position',[xRcol y1 colW_ecdf rowH]); 
        hold(ax.nse_r,'on');
        ax.kge_r = axes('Parent',ax.fig, ...
            'Position',[xRcol y2 colW_ecdf rowH]); 
        hold(ax.kge_r,'on');
        ax.jkge_r = axes('Parent',ax.fig, ...
            'Position',[xRcol y3 colW_ecdf rowH]); 
        hold(ax.jkge_r,'on');
    
        % ECDF axes style
        ecdf_axes = [ax.nse_t ...
            ax.kge_t ax.jkge_t ...
            ax.nse_r ax.kge_r ...
            ax.jkge_r];
        set(ecdf_axes, ...
            'FontName',get(groot,'DefaultAxesFontName'), ...
            'FontSize',fnt_ax, ...
            'LineWidth',1, ...
            'TickDir','out', ...
            'Box','off', ...
            'Layer','top', ...
            'XMinorTick','off', ...
            'YMinorTick','off', ...
            'TickLength',[0.030 0.030]);   % doubled tick length
    
        % Left ECDF column uses left-axis color
        set([ax.nse_t ax.kge_t ax.jkge_t], ...
            'XColor',cT, ...
            'YColor',cT);
        
        % Right ECDF column uses right-axis color
        set([ax.nse_r ax.kge_r ax.jkge_r], ...
            'XColor',cE, ...
            'YColor',cE);
    
        for ah = ecdf_axes
            xlim(ah,[xL xR]);
            ylim(ah,[yL yR]);
        end
        % Overwrite if other than x in [-1,1]
        xlim(ax.jkge_t,[xLb xRb]);
        xlim(ax.jkge_r,[xLb xRb]);

        % Halve the font-relative gap on the bottom ECDF tick labels.
        % This runs only when creating the axes, not on training updates.
        for ah = [ax.jkge_t ax.jkge_r]
            if isprop(ah.XRuler,'TickLabelGapMultiplier')
                ah.XRuler.TickLabelGapMultiplier = ...
                    0.5*ah.XRuler.TickLabelGapMultiplier;
            end
        end

        % No xtick labels on top 2 rows
        set(ax.nse_t,'XTickLabel',[]);
        set(ax.kge_t,'XTickLabel',[]);
        set(ax.nse_r,'XTickLabel',[]);
        set(ax.kge_r,'XTickLabel',[]);
    
        % Labels
        ylabel(ax.nse_t ,'F(NSE_{tt})', ...
            'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
            'FontSize',fnt_lab);
        ylabel(ax.kge_t ,'F(KGE_{tt})', ...
            'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
            'FontSize',fnt_lab);

        if useJKGEpanel
            ylabel(ax.jkge_t, ...
                'F(JKGE_{tt})', ...
                'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                'FontSize',fnt_lab);
            xlabel(ax.jkge_t, ...
                'JKGE_{tt}', ...
                'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                'FontSize',fnt_lab);
        else
            ylabel(ax.jkge_t, ...
                sprintf('$F(%s)$',local_metric_scenario_tex( ...
                fdcMetricName,'tt')), ...
                'Interpreter','latex','FontName',get(groot,'DefaultAxesFontName'), ...
                'FontSize',fnt_lab);
            xlabel(ax.jkge_t, ...
                sprintf('$%s$',local_metric_scenario_tex( ...
                fdcMetricName,'tt')), ...
                'Interpreter','latex','FontName',get(groot,'DefaultAxesFontName'), ...
                'FontSize',fnt_lab);
        end
    
        if rightAvailable
            ylabel(ax.nse_r, ...
                sprintf('F(NSE_{%s})',char(rightField)), ...
                'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                'FontSize',fnt_lab);
            ylabel(ax.kge_r, ...
                sprintf('F(KGE_{%s})',char(rightField)), ...
                'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                'FontSize',fnt_lab);
            if useJKGEpanel
                ylabel(ax.jkge_r, ...
                    sprintf('F(JKGE_{%s})', ...
                    char(rightField)), ...
                    'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                    'FontSize',fnt_lab);
                xlabel(ax.jkge_r, ...
                    sprintf('JKGE_{%s}', ...
                    char(rightField)), ...
                    'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                    'FontSize',fnt_lab);
            else
                ylabel(ax.jkge_r, ...
                    sprintf('$F(%s)$',local_metric_scenario_tex( ...
                    fdcMetricName,char(rightField))), ...
                    'Interpreter','latex','FontName',get(groot,'DefaultAxesFontName'), ...
                    'FontSize',fnt_lab);
                xlabel(ax.jkge_r, ...
                    sprintf('$%s$',local_metric_scenario_tex( ...
                    fdcMetricName,char(rightField))), ...
                    'Interpreter','latex','FontName',get(groot,'DefaultAxesFontName'), ...
                    'FontSize',fnt_lab);
            end
            xlabel(ax.nse_r ,'');
            xlabel(ax.kge_r ,'');
        else
            ylabel(ax.nse_r,'F(NSE_{na})', ...
                'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                'FontSize',fnt_lab);
            ylabel(ax.kge_r,'F(KGE_{na})', ...
                'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                'FontSize',fnt_lab);

            if useJKGEpanel
                ylabel(ax.jkge_r, ...
                    'F(JKGE_{na})', ...
                    'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                    'FontSize',fnt_lab);
                xlabel(ax.jkge_r, ...
                    'JKGE_{na}', ...
                    'Interpreter','tex','FontName',get(groot,'DefaultAxesFontName'), ...
                    'FontSize',fnt_lab);
            else
                ylabel(ax.jkge_r, ...
                    sprintf('$F(%s)$',local_metric_scenario_tex( ...
                    fdcMetricName,'na')), ...
                    'Interpreter','latex','FontName',get(groot,'DefaultAxesFontName'), ...
                    'FontSize',fnt_lab);
                xlabel(ax.jkge_r, ...
                    sprintf('$%s$',local_metric_scenario_tex( ...
                    fdcMetricName,'na')), ...
                    'Interpreter','latex','FontName',get(groot,'DefaultAxesFontName'), ...
                    'FontSize',fnt_lab);
            end
        end
    
        xlabel(ax.nse_t ,'');
        xlabel(ax.kge_t ,'');
    
        % Top/right frame lines
        add_top_right_frame(ax.nse_t, ...
            xL,xR,yL,yR,cT);
        add_top_right_frame(ax.kge_t, ...
            xL,xR,yL,yR,cT);
        add_top_right_frame(ax.jkge_t, ...
            xLb,xRb,yL,yR,cT);
        
        add_top_right_frame(ax.nse_r, ...
            xL,xR,yL,yR,cE);
        add_top_right_frame(ax.kge_r, ...
            xL,xR,yL,yR,cE);
        add_top_right_frame(ax.jkge_r, ...
            xLb,xRb,yL,yR,cE);
    
        % Titles: only top row ECDF
        set_panel_title(ax.nse_t,mdl_tex, ...
            nameT,K_t,nt_nse,fnt_tit, ...
            'K_{\rm t}',true,false);
        set_panel_title(ax.nse_r,mdl_tex, ...
            nameE,K_right,nr_nse,fnt_tit, ...
            klabel_right,rightAvailable,false);
        ax.nse_t.Title.Color = cT;
        ax.nse_r.Title.Color = cE;
    
        title(ax.kge_t ,'');
        title(ax.jkge_t,'');
        title(ax.kge_r ,'');
        title(ax.jkge_r,'');
    
        % History axes
        init_history_axis(ax.lossAx, ...
            fnt_ax);
        init_history_axis(ax.rssAx, ...
            fnt_ax);
        init_history_axis(ax.sibAx, ...
            fnt_ax);
    
        % Only top history title
        title(ax.lossAx, ...
            sprintf(['\\texttt{%s}: ' ...
            'History of loss functions'],mdl_tex), ...
            'Interpreter','latex','FontSize',fnt_tit);
        title(ax.rssAx,'');
        title(ax.sibAx,'');
    
        xlabel(ax.lossAx,'');
        xlabel(ax.rssAx ,'');
        xlabel(ax.sibAx,'SAGE iteration, $i$', ...
            'Interpreter','latex', ...
            'FontSize',fnt_lab);
        set(ax.lossAx,'XTickLabel',[]);
        set(ax.rssAx,'XTickLabel',[]);
    
        % History lines
        yyaxis(ax.lossAx,'left');
        ax.lossLineT = plot(ax.lossAx, ...
            nan,nan,'-o', ...
            'Color',cT, ...
            'MarkerFaceColor',cT, ...
            'LineWidth',1.3, ...
            'DisplayName',nameT);
        yyaxis(ax.lossAx,'right');
        ax.lossLineR = plot(ax.lossAx, ...
            nan,nan,'-s', ...
            'Color',cE, ...
            'MarkerFaceColor',cE, ...
            'LineWidth',1.3, ...
            'DisplayName',nameE);
        
        yyaxis(ax.rssAx,'left');
        ax.rssLineT = plot(ax.rssAx, ...
            nan,nan,'-o', ...
            'Color',cT, ...
            'MarkerFaceColor',cT, ...
            'LineWidth',1.3, ...
            'HandleVisibility','off');
        yyaxis(ax.rssAx,'right');
        ax.rssLineR = plot(ax.rssAx, ...
            nan,nan,'-s', ...
            'Color',cE, ...
            'MarkerFaceColor',cE, ...
            'LineWidth',1.3, ...
            'HandleVisibility','off');
    
        yyaxis(ax.sibAx,'left');
        ax.sibLineT = plot(ax.sibAx, ...
            nan,nan,'-o', ...
            'Color',cT, ...
            'MarkerFaceColor',cT, ...
            'LineWidth',1.3, ...
            'HandleVisibility','off');
        yyaxis(ax.sibAx,'right');
        ax.sibLineR = plot(ax.sibAx, ...
            nan,nan,'-s', ...
            'Color',cE, ...
            'MarkerFaceColor',cE, ...
            'LineWidth',1.3, ...
            'HandleVisibility','off');
    
        % Only top history panel has legend, inside the axes
        local_top_history_legend(ax.lossAx, ...
            ax.lossLineT,ax.lossLineR,cT,cE, ...
            nameT,nameE,fnt_leg);
    
        % ECDF handles
        tags = {'nse_t','kge_t', ...
            'jkge_t','nse_r','kge_r','jkge_r'};
        for k = 1:numel(tags)
            tg = tags{k};
            ax.ecdf.(tg).patch = ...
                gobjects(1);
            ax.ecdf.(tg).line = ...
                gobjects(1);
            ax.ecdf.(tg).med_v = ...
                gobjects(1);
            ax.ecdf.(tg).med_h = ...
                gobjects(1);
            ax.ecdf.(tg).med_s = ...
                gobjects(1);
            ax.ecdf.(tg).med_txt = ...
                gobjects(1);
        end
    end

    % Refresh these labels on every call so an existing dashboard follows
    % the selected 6a/6b/6c FDC formulation instead of retaining S_fdc.
    update_ecdf_metric_labels(ax,useJKGEpanel,fdcMetricName, ...
        rightField,rightAvailable,fnt_lab);
    
    % Top-row titles update only
    set_panel_title(ax.nse_t, ...
        mdl_tex,nameT,K_t,nt_nse, ...
        fnt_tit,'K_{\rm t}',true,false);
    set_panel_title(ax.nse_r, ...
        mdl_tex,nameE,K_right,nr_nse, ...
        fnt_tit,klabel_right, ...
        rightAvailable,false);
    
    ax.nse_t.Title.Color = cT;
    ax.nse_r.Title.Color = cE;
    
    % ------------
    % Update ECDFs
    % ------------
    ax = update_manual_ecdf(ax, ...
        'nse_t',ax.nse_t,NSEtt_use, ...
        c,xL,xR,fa,lw,'NSE','tt', ...
        fnt_med);
    ax = update_manual_ecdf(ax, ...
        'kge_t',ax.kge_t,KGEtt_use, ...
        c,xL,xR,fa,lw,'KGE','tt', ...
        fnt_med);
    ax = update_manual_ecdf(ax, ...
        'jkge_t',ax.jkge_t,bottomT, ...
        c,xLb,xRb,fa,lw,bottomName,'tt', ...
        fnt_med);
    
    if rightAvailable
        ax = update_manual_ecdf(ax, ...
            'nse_r',ax.nse_r ,NSEr_use, ...
            c,xL,xR,fa,lw,'NSE', ...
            char(rightField),fnt_med);
        ax = update_manual_ecdf(ax, ...
            'kge_r',ax.kge_r ,KGEr_use, ...
            c,xL,xR,fa,lw,'KGE', ...
            char(rightField),fnt_med);
        ax = update_manual_ecdf(ax, ...
            'jkge_r',ax.jkge_r,bottomR, ...
            c,xLb,xRb,fa,lw,bottomName, ...
            char(rightField),fnt_med);
    else
        clear_manual_ecdf(ax,'nse_r');  
        set_unavailable_label(ax.nse_r,true);
        clear_manual_ecdf(ax,'kge_r');  
        set_unavailable_label(ax.kge_r,true);
        clear_manual_ecdf(ax,'jkge_r'); 
        set_unavailable_label(ax.jkge_r,true);
    end

    % ----------------
    % Update histories
    % ----------------
    [LossT,LossR] = get_history_pair( ...
        prf,lossField,i,rightField, ...
        K_t,K_e,rightAvailable);
    [RSST,RSSR] = get_history_pair( ...
        prf,'RSS',i,rightField, ...
        K_t,K_e,rightAvailable);
    [SibT,SibR] = get_history_pair( ...
        prf,'Sib',i,rightField, ...
        1,1,rightAvailable);
    
    update_dual_history_loss( ...
        ax.lossAx,ax.lossLineT, ...
        ax.lossLineR, ...
        LossT,LossR,cT,cE, ...
        lossLeftLabel,right_loss_label( ...
        loss_fnc,rightField, ...
        rightAvailable,fdcFormulation),fnt_lab);
    
    update_dual_history(ax.rssAx, ...
        ax.rssLineT,ax.rssLineR, ...
        RSST,RSSR,cT,cE, ...
        ['$\Sigma \mathrm{RSS}_' ...
        '{\rm tt}$'],right_label( ...
        'RSS',rightField, ...
        rightAvailable),fnt_lab);
    
    update_dual_history(ax.sibAx, ...
        ax.sibLineT,ax.sibLineR, ...
        SibT,SibR,cT,cE, ...
        ['$\widehat{\mathcal{S}}_' ...
        '{\mathrm{ib}_{\rm tt}}$'], ...
        right_label('Sib', ...
        rightField,rightAvailable), ...
        fnt_lab);
    
    xmax = max([2 numel(LossT) ...
        numel(LossR) numel(RSST) ...
        numel(RSSR) numel(SibT) ...
        numel(SibR)]); 
    
    xt = local_history_xticks(xmax);
    
    xlim(ax.lossAx,[1 xmax]);
    xlim(ax.rssAx,[1 xmax]);
    xlim(ax.sibAx,[1 xmax]);
    
    set_history_xticks(ax.lossAx,xt,false);
    set_history_xticks(ax.rssAx,xt,false);
    set_history_xticks(ax.sibAx,xt,true);
    
    % Add horizontal axis on top
    if i == 1 ...
            || needInit
        add_history_top_frame(ax.lossAx,cAx);
        add_history_top_frame(ax.rssAx,cAx);
        add_history_top_frame(ax.sibAx,cAx);
    end
    
    if ~isdeployed
        drawnow expose
        pause(0.001)
    else
        drawnow limitrate nocallbacks
    end

end

% =================
% Helpers functions
% =================
function [rightField,nameE,klabel_right, ...
    rightAvailable] = choose_compare_scenario(prf,samp_word)
    % ECDFs use the current basin-wise metric vectors. Since the prf
    % refactor, these live in prf.curr; iteration histories such as RSS
    % live separately in prf.iter and must not determine ECDF availability.
    if isfield(prf,'curr') && isfield(prf.curr,'NSE')
        metric = prf.curr.NSE;
    else
        metric = struct();
    end

    if has_nonempty_scenario(metric,'ee')
        rightField = "ee";
        nameE = sprintf(['Eval basins |' ...
            ' eval %s'],samp_word);
        klabel_right = 'K_{\rm e}';
        rightAvailable = true;
    elseif has_nonempty_scenario(metric,'et')
        rightField = "et";
        nameE = sprintf(['Eval basins |' ...
            ' train %s'],samp_word);
        klabel_right = 'K_{\rm e}';
        rightAvailable = true;
    elseif has_nonempty_scenario(metric,'te')
        rightField = "te";
        nameE = sprintf(['Train basins |' ...
            ' eval %s'],samp_word);
        klabel_right = 'K_{\rm t}';
        rightAvailable = true;
    else
        rightField = "";
        nameE = 'Unavailable';
        klabel_right = '';
        rightAvailable = false;
    end
end

function [rightField,nameE,klabel_right, ...
    rightAvailable] = choose_joint_compare_scenario( ...
    prf,names,samp_word)
%CHOOSE_JOINT_COMPARE_SCENARIO Select a stored named-observation scenario.

    candidates = ["ee" "et" "te"];
    labels = { ...
        sprintf('Eval basins | eval %s',samp_word), ...
        sprintf('Eval basins | train %s',samp_word), ...
        sprintf('Train basins | eval %s',samp_word)};
    basinLabels = {'K_{\rm e}','K_{\rm e}','K_{\rm t}'};
    rightField = "";
    nameE = 'Unavailable';
    klabel_right = '';
    rightAvailable = false;
    for s = 1:numel(candidates)
        if has_joint_scenario(prf,names,candidates(s))
            rightField = candidates(s);
            nameE = labels{s};
            klabel_right = basinLabels{s};
            rightAvailable = true;
            return
        end
    end
end

function tf = has_joint_scenario(prf,names,scenario)
%HAS_JOINT_SCENARIO Test named NSE vectors for one basin/period scenario.

    tf = false;
    if ~isfield(prf,'curr') ...
            || ~isfield(prf.curr,'joint')
        return
    end
    for j = 1:numel(names)
        name = char(names(j));
        if isfield(prf.curr.joint,name) ...
                && isfield(prf.curr.joint.(name),'NSE') ...
                && has_nonempty_scenario( ...
                prf.curr.joint.(name).NSE,char(scenario))
            tf = true;
            return
        end
    end
end

function c = local_sage_model_color(model)
%LOCAL_SAGE_MODEL_COLOR Return the canonical SAGE model color.

    colors = [ ...
        0.0000 0.4470 0.7410;
        0.9500 0.4500 0.0000;
        0.9290 0.6940 0.1250;
        0.6000 0.2000 0.8000;
        0.0000 0.5500 0.0000;
        0.6500 0.3500 0.0000;
        0.0000 0.7000 0.7000;
        0.9000 0.0000 0.9000;
        0.6500 0.6500 0.6500];
    if isnumeric(model) && isscalar(model) ...
            && model >= 1 && model <= size(colors,1)
        c = colors(model,:);
    else
        c = [0 0 0];
    end
end

function tf = has_nonempty_scenario(S,scn)
%HAS_NONEMPTY_SCENARIO Test whether a scenario contains finite values.

    tf = false;
    if isstruct(S) && isfield(S,scn)
        x = S.(scn);
        tf = ~isempty(x) ...
            && any(isfinite(x));
    end
end

function x = get_metric_vector(S,scn)
%GET_METRIC_VECTOR Return one scenario metric as a column vector.

    x = [];
    if isstruct(S) && isfield(S,scn)
        x = S.(scn)(:);
    end
end

function [x,n] = finite_vec(v)
%FINITE_VEC Return finite entries and their count.

    if isempty(v) 
        x = []; n = 0; 
        return; 
    end
    v = v(:);
    I = isfinite(v);
    x = v(I);
    n = nnz(I);
end

function [lossTitle,leftLab,fld] = ...
    get_loss_strings(loss_fnc,fdcFormulation)
%GET_LOSS_STRINGS Return labels and history field for a loss function.

    if nargin < 2
        fdcFormulation = 1;
    end
    switch loss_fnc
        case 1
            lossTitle = '\Sigma SAR';
            leftLab = ['$\Sigma {\rm SAR}_' ...
                '{\rm tt}$'];
            fld = 'SAR';
        case 2
            lossTitle = '\Sigma GLS';
            leftLab = ['$\Sigma {\rm GLS}_' ...
                '{\rm tt}$'];
            fld = 'GLS';
        case 3
            lossTitle = '\Sigma(1-NSE)';
            leftLab = ['$\Sigma(1-{\rm NSE})_' ...
                '{\rm tt}$'];
            fld = 'L';
        case 4
            lossTitle = '\Sigma(1-KGE)';
            leftLab = ['$\Sigma(1-{\rm KGE})_' ...
                '{\rm tt}$'];
            fld = 'L';
        case 5
            lossTitle = '\Sigma Huber';
            leftLab = ['$\Sigma {\rm Huber}_' ...
                '{\rm tt}$'];
            fld = 'Huber';
        case 6
            [~,~,~,~,fdcLossTex] = ...
                local_fdc_print_names(fdcFormulation);
            lossTitle = ['\Sigma ' fdcLossTex];
            leftLab = ['$\Sigma ' local_fdc_loss_scenario_tex( ...
                fdcFormulation,'tt') '$'];
            fld = 'L';
        case 7
            lossTitle = '\Sigma(1-JKGE)';
            leftLab = ['$\Sigma(1-{\rm JKGE})_' ...
                '{\rm tt}$'];
            fld = 'L';
        otherwise
            error('Unknown loss function.');
    end
end

function lab = local_joint_loss_label(loss,name,scenario)
%LOCAL_JOINT_LOSS_LABEL Add an upright observation to a loss label.

    obs = upper(char(string(name)));
    scn = lower(char(string(scenario)));
    switch double(loss.fnc)
        case 1
            core = '{\rm SAR}';
        case 2
            core = '{\rm GLS}';
        case 3
            core = '(1-{\rm NSE})';
        case 4
            core = '(1-{\rm KGE})';
        case 5
            core = '{\rm Huber}';
        case 6
            formulation = 1;
            if isfield(loss,'fdc') ...
                    && isfield(loss.fdc,'formulation')
                formulation = loss.fdc.formulation;
            elseif isfield(loss,'formulation')
                formulation = loss.formulation;
            end
            [~,~,~,~,core] = ...
                local_fdc_print_names(formulation);
        case 7
            core = '(1-{\rm JKGE})';
        otherwise
            core = '\mathcal{L}';
    end
    lab = sprintf('$\\Sigma %s^{\\rm %s}_{\\rm %s}$', ...
        core,obs,scn);
end

function init_history_axis(axh,fs)
%INIT_HISTORY_AXIS Apply common formatting to a history axis.

    set(axh,'FontSize',fs, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off');
    
    xtickformat(axh,'%d');
    
    try
        axh.YAxis(1).Exponent = 0;
        axh.YAxis(2).Exponent = 0;
    catch
    end

end

function set_panel_title(axh, ...
    mdl_tex,str,K,nFinite,fnt, ...
    klabel,isAvailable,showModel)
%SET_PANEL_TITLE Set the title for one monitoring panel.

    if nargin < 8 ...
            || isempty(showModel)
        showModel = true;
    end
    
    if ~isAvailable
        if showModel
            axh.Title.String = ...
                sprintf('%s: %s', ...
                mdl_tex,str);
        else
            axh.Title.String = str;
        end
    elseif isempty(klabel)
        if showModel
            axh.Title.String = ...
                sprintf('%s: %s', ...
                mdl_tex,str);
        else
            axh.Title.String = str;
        end
    else
        if showModel
            axh.Title.String = ...
                sprintf(['%s: %s ' ...
                '($%s = %d$)'],mdl_tex, ...
                str,klabel,K);
        else
            axh.Title.String = ...
                sprintf('%s ($%s = %d$)', ...
                str,klabel,K);
        end
        if nFinite ~= K
            axh.Title.String = sprintf( ...
                '%s; $n_{\rm NSE} = %d$', ...
                axh.Title.String,nFinite);
        end
    end
    axh.Title.Interpreter = 'latex';
    axh.Title.FontName = get(groot,'DefaultAxesFontName');
    axh.Title.FontSize = fnt;
end

function h = add_top_right_frame(axh, ...
    xL,xR,yL,yR,clr)
%ADD_TOP_RIGHT_FRAME Draw the upper and right frame segments.

    hold(axh,'on');
    if isappdata(axh, ...
            'TopRightFrameHandles')
        hh = getappdata(axh, ...
            'TopRightFrameHandles');
        try
            if isfield(hh,'top') ...
                    && isgraphics(hh.top)
                delete(hh.top); 
            end
            if isfield(hh,'right') ...
                    && isgraphics(hh.right)
                delete(hh.right); 
            end
        catch
        end
    end
    h.top = line(axh,[xL xR],[yR yR], ...
        'Color',clr, ...
        'LineWidth',1, ...
        'HandleVisibility','off', ...
        'Clipping','off');
    h.right = line(axh,[xR xR],[yL yR], ...
        'Color',clr, ...
        'LineWidth',1, ...
        'HandleVisibility','off', ...
        'Clipping','off');
    setappdata(axh,'TopRightFrameHandles',h);
end

function ax = update_manual_ecdf(ax, ...
    tag,axh,z,c,xL,xR,fa,lw, ...
    metName,scn,fnt_med,observable)
%UPDATE_MANUAL_ECDF Refresh one manually managed ECDF panel.

    if nargin < 13
        observable = '';
    end

    if isempty(z)
        set_unavailable_label(axh,true);
        clear_manual_ecdf(ax,tag);
        return
    else
        set_unavailable_label(axh,false);
    end
    
    [f,x] = sage_ecdf(z);
    [xs,fs] = ecdf_to_stairs_fixed(x,f,xL,xR);
    [xp,yp] = stairs_fill_poly(xs,fs);
    
    if ~isgraphics(ax.ecdf.(tag).patch)
        ax.ecdf.(tag).patch = patch(axh, ...
            xp,yp,c, ...
            'FaceAlpha',fa, ...
            'EdgeAlpha',0);
        ax.ecdf.(tag).line = plot(axh, ...
            xs,fs, ...
            'Color',c, ...
            'LineWidth',lw);
    else
        set(ax.ecdf.(tag).patch, ...
            'XData',xp, ...
            'YData',yp, ...
            'FaceColor',c);
        set(ax.ecdf.(tag).line, ...
            'XData',xs, ...
            'YData',fs, ...
            'Color',c);
    end
    
    medX = median(z);
    medF = ecdf_value_from_stairs( ...
        xs,fs,medX);
    [xMark,yMark,x1m,x2m,xt,ha,mode] = ...
        median_marker_geom(medX,medF,xL,xR);
    
    if ~isgraphics(ax.ecdf.(tag).med_v)
        ax.ecdf.(tag).med_v = line(axh, ...
            nan,nan, ...
            'Color','k', ...
            'LineWidth',1.0, ...
            'HandleVisibility','off');
        ax.ecdf.(tag).med_h = line(axh, ...
            nan,nan, ...
            'Color','k', ...
            'LineWidth',1.0, ...
            'HandleVisibility','off');
        ax.ecdf.(tag).med_s = line(axh, ...
            nan,nan, ...
            'Marker','s', ...
            'MarkerFaceColor','k', ...
            'MarkerEdgeColor','k', ...
            'LineStyle','none', ...
            'MarkerSize',5, ...
            'HandleVisibility','off');
        ax.ecdf.(tag).med_txt = ...
            text(axh, ...
            nan,nan,'', ...
            'Interpreter','latex', ...
            'FontWeight','bold', ...
            'HorizontalAlignment',ha, ...
            'VerticalAlignment','middle', ...
            'Color','k', ...
            'HandleVisibility','off', ...
            'FontSize',fnt_med);
    end
    
    if ~mode.specialHalfAtLeft
        set(ax.ecdf.(tag).med_v, ...
            'XData',[xMark xMark], ...
            'YData',[0 yMark], ...
            'Visible','on');
    else
        set(ax.ecdf.(tag).med_v, ...
            'XData',[], ...
            'YData',[], ...
            'Visible','off');
    end
    set(ax.ecdf.(tag).med_h, ...
        'XData',[x1m x2m], ...
        'YData',[yMark yMark], ...
        'Visible','on');
    set(ax.ecdf.(tag).med_s, ...
        'XData',xMark, ...
        'YData',yMark, ...
        'Visible','on');
    scnTex = lower(strtrim(char(string(scn))));
    if isempty(observable)
        metricTex = local_median_metric_subscript(metName);
        medLabel = sprintf( ...
            '$\\widehat{T}_{%s,\\mathrm{%s}} = %.3f$', ...
            metricTex,scnTex,medX);
    else
        obsTex = upper(strtrim(char(string(observable))));
        medLabel = sprintf( ...
            '$\\widehat{T}^{\\mathrm{%s}}_{\\mathrm{%s}} = %.3f$', ...
            obsTex,scnTex,medX);
    end
    set(ax.ecdf.(tag).med_txt, ...
        'Position',[xt yMark 0], ...
        'String',medLabel, ...
        'HorizontalAlignment',ha, ...
        'Visible','on', ...
        'FontSize',fnt_med);
    
    xlim(axh,[xL xR]);
    ylim(axh,[0 1]);
    end
    
    function clear_manual_ecdf(ax,tag)
    %CLEAR_MANUAL_ECDF Remove graphics from one managed ECDF panel.

    flds = {'patch','line', ...
        'med_v','med_h','med_s','med_txt'};
    for k = 1:numel(flds)
        h = ax.ecdf.(tag).(flds{k});
        if isgraphics(h)
            if strcmp(flds{k},'med_txt')
                set(h,'String','', ...
                    'Visible','off');
            else
                set(h,'XData',[], ...
                    'YData',[], ...
                    'Visible','off');
            end
        end
    end
end

function [yT,yR] = get_history_pair(prf, ...
    fld,i,rightField,K_t,K_e,rightAvailable)

    % Iteration histories are stored under prf.iter. The generic Sib
    % diagnostic plotted here is the NSE-based integrated basin score.
    if strcmp(fld,'Sib')
        histFld = 'Sib_NSE';
        sT = 1;
    else
        histFld = fld;
        sT = K_t;
    end

    if ~isfield(prf,'iter') ...
            || ~isfield(prf.iter,histFld)
        error('print_SAGE:MissingHistory', ...
            'Missing iteration history prf.iter.%s.',histFld);
    end
    hist = prf.iter.(histFld);
    
    it = min(i,numel(hist.tt));
    iiT = 1:it;
    iiT = iiT(isfinite(hist.tt(iiT)));
    yT = sT * hist.tt(iiT);
    
    if ~rightAvailable
        yR = [];
        return
    end
    
    switch char(rightField)
        case 'te'
            sR = strcmp(fld,'Sib') * 1 ...
                + ~strcmp(fld,'Sib') * K_t;
        case {'et','ee'}
            sR = strcmp(fld,'Sib') * 1 ...
                + ~strcmp(fld,'Sib') * K_e;
        otherwise
            sR = 1;
    end
    
    rightField = char(rightField);
    iv = min(i,numel(hist.(rightField)));
    iiR = 1:iv;
    iiR = iiR(isfinite(hist.(rightField)(iiR)));
    yR = sR * hist.(rightField)(iiR);
end

function lab = right_label(fld, ...
    rightField,rightAvailable)
%RIGHT_LABEL Return the right-axis metric label.

    if ~rightAvailable
        lab = 'Unavailable';
        return
    end
    
    switch fld
        case 'L'
            core = '\mathcal{L}';
        case 'RSS'
            core = '\Sigma \mathrm{RSS}';
        case 'Sib'
            lab = ...
                sprintf(['$\\widehat{\\mathcal{S}}_' ...
                '{\\mathrm{ib}_{\\rm %s}}$'], ...
                char(rightField));
            return
        case 'SAR'
            core = '\Sigma {\rm SAR}';
        case 'GLS'
            core = '\Sigma {\rm GLS}';
        case 'Huber'
            core = '\Sigma {\rm Huber}';
        otherwise
            core = fld;
    end
    lab = sprintf('$%s_{\\rm %s}$', ...
        core,char(rightField));
end

function lab = right_loss_label(loss_fnc, ...
    rightField,rightAvailable,fdcFormulation)
%RIGHT_LOSS_LABEL Return the right-axis loss label.

    if nargin < 4
        fdcFormulation = 1;
    end

    if ~rightAvailable
        lab = 'Unavailable';
        return
    end
    
    scn = char(rightField);
    
    switch loss_fnc
        case 1
            lab = sprintf(['$\\Sigma ' ...
                '{\\rm SAR}_{\\rm %s}$'],scn);
        case 2
            lab = sprintf(['$\\Sigma ' ...
                '{\\rm GLS}_{\\rm %s}$'],scn);
        case 3
            lab = sprintf(['$\\Sigma' ...
                '(1-{\\rm NSE})_{\\rm %s}$'],scn);
        case 4
            lab = sprintf(['$\\Sigma' ...
                '(1-{\\rm KGE})_{\\rm %s}$'],scn);
        case 5
            lab = sprintf(['$\\Sigma ' ...
                '{\\rm Huber}_{\\rm %s}$'],scn);
        case 6
            fdcLossTex = local_fdc_loss_scenario_tex( ...
                fdcFormulation,scn);
            lab = sprintf('$\\Sigma %s$',fdcLossTex);
        case 7
            lab = sprintf(['$\\Sigma' ...
                '(1-{\\rm JKGE})_{\\rm %s}$'],scn);
        otherwise
            lab = sprintf(['$\\mathcal{L}_' ...
                '{\\rm %s}$'],scn);
    end

end

function update_dual_history(axh, ...
    hT,hR,yT,yR,cT,cE,leftLab, ...
    rightLab,fsLab)
%UPDATE_DUAL_HISTORY Refresh a dual-axis history panel.

    useLogT = should_use_log_history_scale(yT(:));
    useLogR = should_use_log_history_scale(yR(:));

    update_history_side(axh,'left',hT,yT,cT, ...
        leftLab,fsLab,useLogT,false);
    update_history_side(axh,'right',hR,yR,cE, ...
        rightLab,fsLab,useLogR,true);
end

function update_history_side(axh,side,hLine,y,clr, ...
    axisLabel,fsLab,useLog,allowUnavailable)
%UPDATE_HISTORY_SIDE Refresh one side of a dual-axis history panel.

    yyaxis(axh,side);
    axh.YColor = clr;
    try
        if strcmp(side,'left')
            axh.YAxis(1).Exponent = 0;
        else
            axh.YAxis(2).Exponent = 0;
        end
    catch
    end

    if isempty(y)
        set(hLine,'XData',nan,'YData',nan);
        axh.YScale = 'linear';
        if allowUnavailable
            ylabel(axh,'Unavailable', ...
                'Interpreter','none', ...
                'FontSize',fsLab);
        else
            ylabel(axh,axisLabel, ...
                'Interpreter','latex', ...
                'FontSize',fsLab);
        end
        return
    end

    if useLog
        axh.YScale = 'log';
        yPlot = y;
        yPlot(~isfinite(yPlot) ...
            | yPlot <= 0) = nan;
        set(hLine, ...
            'XData',1:numel(yPlot), ...
            'YData',yPlot);
        local_set_dynamic_ylim_log(axh,side,yPlot);
        local_no_axis_exponent(axh);
        ylabel(axh,axisLabel, ...
            'Interpreter','latex', ...
            'FontSize',fsLab);
        return
    end

    axh.YScale = 'linear';
    [scalePow,scaleVal] = nice_eng_scale(y(:));
    yPlot = y./scaleVal;
    set(hLine, ...
        'XData',1:numel(yPlot), ...
        'YData',yPlot);
    local_set_dynamic_ylim(axh,side,yPlot);
    local_no_axis_exponent(axh);

    if scalePow ~= 0
        ylabel(axh, ...
            sprintf('%s $\\times 10^{%d}$', ...
            axisLabel,scalePow), ...
            'Interpreter','latex', ...
            'FontSize',fsLab);
    else
        ylabel(axh,axisLabel, ...
            'Interpreter','latex', ...
            'FontSize',fsLab);
    end
end
function [scalePow,scaleVal] = nice_eng_scale(y)
%NICE_ENG_SCALE Scale so plotted values stay visually compact.
% Example: 6e5 -> x10^4 gives 60 instead of x10^3 gives 600.

    y = y(isfinite(y));
    
    if isempty(y) ...
            || all(y == 0)
        scalePow = 0;
        scaleVal = 1;
        return
    end
    
    m = max(abs(y));
    
    if m >= 100
        scalePow = floor(log10(m));
    else
        scalePow = 0;
    end
    
    scaleVal = 10^scalePow;
end

function set_unavailable_label(axh,tf)
%SET_UNAVAILABLE_LABEL Show or hide the unavailable-data annotation.

    key = 'SAGE_UnavailableText';
    if isappdata(axh,key)
        h = getappdata(axh,key);
    else
        h = text(axh,0.5,0.5,'Unavailable', ...
            'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontWeight','bold', ...
            'FontSize',18, ...
            'Color',[0.35 0.35 0.35], ...
            'Interpreter','none', ...
            'Visible','off');
        setappdata(axh,key,h);
    end
    if isgraphics(h)
        if tf
            set(h,'Visible','on');
        else
            set(h,'Visible','off');
        end
    end
end

function update_dual_history_loss(axh, ...
    hT,hR,yT,yR,cT,cE,leftLab, ...
    rightLab,fsLab)

    % useLog = should_use_log_loss_scale([yT(:); yR(:)], ...
    %     numel(yT));
    useLogT = should_use_log_loss_scale(yT(:),numel(yT));
    useLogR = should_use_log_loss_scale(yR(:),numel(yR));
    useLog = useLogT || useLogR;

    if useLog
        yyaxis(axh,'left');
        try
            axh.YAxis(1).Exponent = 0;
        catch
        end
        axh.YColor = cT;
        axh.YScale = 'log';
        if isempty(yT)
            set(hT,'XData',nan, ...
                'YData',nan);
        else
            yTplot = yT;
            yTplot(~isfinite(yTplot) ...
                | yTplot <= 0) = nan;
            set(hT,'XData',1:numel(yTplot), ...
                'YData',yTplot);
            local_set_dynamic_ylim_log(axh, ...
                'left',yTplot);
            local_no_axis_exponent(axh);
        end
        ylabel(axh,leftLab, ...
            'Interpreter','latex', ...
            'FontSize',fsLab);
    
        yyaxis(axh,'right');
        try
            axh.YAxis(2).Exponent = 0;
        catch
        end
        axh.YColor = cE;
        axh.YScale = 'log';
        if isempty(yR)
            set(hR,'XData',nan, ...
                'YData',nan);
            ylabel(axh,'Unavailable', ...
                'Interpreter','none', ...
                'FontSize',fsLab);
        else
            yRplot = yR;
            yRplot(~isfinite(yRplot) ...
                | yRplot <= 0) = nan;
            set(hR,'XData',1:numel(yRplot), ...
                'YData',yRplot);
            local_set_dynamic_ylim_log(axh, ...
                'right',yRplot);
            local_no_axis_exponent(axh);
            ylabel(axh,rightLab, ...
                'Interpreter','latex', ...
                'FontSize',fsLab);
        end
    else
        %[pow,sc] = nice_eng_scale([yT(:); yR(:)]);
        [powT,scT] = nice_eng_scale(yT(:));
        [powR,scR] = nice_eng_scale(yR(:));
        yyaxis(axh,'left');
        try
            axh.YAxis(1).Exponent = 0;
        catch
        end
        axh.YColor = cT;
        axh.YScale = 'linear';
        if isempty(yT)
            set(hT,'XData',nan, ...
                'YData',nan);
        else
            yTplot = yT./scT;
            set(hT,'XData',1:numel(yT), ...
                'YData',yTplot);
            local_set_dynamic_ylim(axh, ...
                'left',yTplot);
            local_no_axis_exponent(axh);
        end
        if powT ~= 0
            ylabel(axh, ...
                sprintf('%s $\\times 10^{%d}$', ...
                leftLab,powT), ...
                'Interpreter','latex', ...
                'FontSize',fsLab);
        else
            ylabel(axh,leftLab, ...
                'Interpreter','latex', ...
                'FontSize',fsLab);
        end
    
        yyaxis(axh,'right');
        try
            axh.YAxis(2).Exponent = 0;
        catch
        end
        axh.YColor = cE;
        axh.YScale = 'linear';
        if isempty(yR)
            set(hR,'XData',nan, ...
                'YData',nan);
            ylabel(axh,'Unavailable', ...
                'Interpreter','none', ...
                'FontSize',fsLab);
        else
            yRplot = yR./scR;
            set(hR,'XData',1:numel(yR), ...
                'YData',yRplot);
            local_set_dynamic_ylim(axh, ...
                'right',yRplot);
            local_no_axis_exponent(axh);
            if powR ~= 0
                ylabel(axh, ...
                    sprintf('%s $\\times 10^{%d}$', ...
                    rightLab,powR), ...
                    'Interpreter','latex', ...
                    'FontSize',fsLab);
            else
                ylabel(axh,rightLab, ...
                    'Interpreter','latex', ...
                    'FontSize',fsLab);
            end
        end
    end
end

function tf = should_use_log_loss_scale(y,nIter)
%SHOULD_USE_LOG_LOSS_SCALE Select logarithmic scaling for loss history.

    y = y(isfinite(y) & y > 0);
    if nIter < 5 ...
            || numel(y) < 2
        tf = false;
        return
    end
    tf = (max(y) / min(y)) >= 100;
end

function h = add_history_top_frame(axh,clr)
%ADD_HISTORY_TOP_FRAME Draw solid top frame line for one history axes
% in normalized figure coordinates so it behaves well with yyaxis.

    if ~isgraphics(axh)
        h = gobjects(1);
        return
    end
    fig = ancestor(axh,'figure');
    
    % delete old line if present
    if isappdata(axh,'HistoryTopFrameHandle')
        hOld = getappdata(axh,'HistoryTopFrameHandle');
        try
            if isgraphics(hOld)
                delete(hOld);
            end
        catch
        end
    end
    
    % drawnow limitrate nocallbacks;
    
    pos = axh.Position;   % normalized figure coordinates
    x1 = pos(1);
    x2 = pos(1) + pos(3);
    y  = pos(2) + pos(4);
    
    h = annotation(fig,'line',[x1 x2],[y y], ...
        'Color',clr, ...
        'LineWidth',1, ...
        'LineStyle','-', ...
        'HitTest','off');
    
    setappdata(axh,'HistoryTopFrameHandle',h);

end

function refresh_history_top_frames(fig)
%REFRESH_HISTORY_TOP_FRAMES Redraw top frame lines after figure resize.

    if ~isgraphics(fig)
        return
    end
    
    if ~isappdata(fig,'HistoryAxesHandles') ...
            || ~isappdata(fig,'HistoryTopFrameColor')
        return
    end
    
    axs = getappdata(fig,'HistoryAxesHandles');
    clr = getappdata(fig,'HistoryTopFrameColor');
    
    axs = axs(isgraphics(axs));
    
    for k = 1:numel(axs)
        add_history_top_frame(axs(k),clr);
    end

end

function set_history_xticks(axh,xt,showLabels)
%SET_HISTORY_XTICKS Force unique integer ticks on yyaxis history panels.

    if isempty(xt)
        xt = 1;
    end
    
    xt = unique(round(xt(:).'),'stable');
    
    % Important: x-axis is shared between yyaxis left/right.
    % Do not clear XTickLabel again after switching yyaxis.
    set(axh,'XTick',xt);
    
    if showLabels
        set(axh,'XTickLabel',compose('%d',xt), ...
            'XTickLabelMode','manual', ...
            'XTickLabelRotation',0);
    else
        set(axh,'XTickLabel',[], ...
            'XTickLabelMode','manual', ...
            'XTickLabelRotation',0);
    end
    
    set(axh,'TickDir','out', ...
        'TickLength',[0.012 0.012]);
end

function xt = local_history_xticks(xmax)
%LOCAL_HISTORY_XTICKS Clean iteration ticks for history plots.

    xmax = max(1,round(xmax));
    
    if xmax <= 5
        xt = 1:xmax;
        return
    end

    % At most about six horizontal labels. Select a conventional 1/2/5
    % interval so values remain easy to read (for example, a 250-iteration
    % run uses 1, 50, 100, 150, 200, 250).
    rawStep = xmax/5;
    decade = 10^floor(log10(rawStep));
    scaled = rawStep/decade;
    if scaled <= 1
        multiplier = 1;
    elseif scaled <= 2
        multiplier = 2;
    elseif scaled <= 5
        multiplier = 5;
    else
        multiplier = 10;
    end
    step = max(1,round(multiplier*decade));
    xt = step:step:xmax;

    if ~isempty(xt) && xmax-xt(end) < 0.45*step
        xt(end) = xmax;
    end
    xt = unique([1 xt xmax]);
end

function local_set_dynamic_ylim(axh,side,y)
%LOCAL_SET_DYNAMIC_YLIM Set robust linear limits on one y-axis.

    yyaxis(axh,side);
    
    y = y(isfinite(y));
    
    if isempty(y)
        set(axh,'YLimMode','auto');
        return
    end
    
    ymin = min(y);
    ymax = max(y);
    
    if ymin == ymax
        pad = max(0.05*abs(ymin),0.05);
    else
        pad = 0.08*(ymax - ymin);
    end
    
    yl = [ymin - pad, ymax + pad];
    
    % Keep nonnegative histories anchored at zero when appropriate
    if ymin >= 0
        yl(1) = max(0,yl(1));
    end
    
    set(axh,'YLim',yl);

end

function local_set_dynamic_ylim_log(axh,side,y)
%LOCAL_SET_DYNAMIC_YLIM_LOG Set robust logarithmic limits on one y-axis.

    yyaxis(axh,side);
    
    y = y(isfinite(y) & y > 0);
    
    if isempty(y)
        set(axh,'YLimMode','auto');
        return
    end
    
    ymin = min(y);
    ymax = max(y);
    
    if ymin == ymax
        fac = 1.25;
        yl = [ymin/fac, ymax*fac];
    else
        yl = [ymin/1.25, ymax*1.25];
    end
    
    set(axh,'YLim',yl);

end

function tf = should_use_log_history_scale(y)
%SHOULD_USE_LOG_HISTORY_SCALE Test whether history spans a log-scale range.

    y = y(isfinite(y) ...
        & y > 0);
    if numel(y) < 2
        tf = false;
        return
    end
    tf = max(y)/min(y) >= 100;
end

function local_no_axis_exponent(axh)
%LOCAL_NO_AXIS_EXPONENT Suppress automatic y-axis exponent labels.

    try
        axh.YAxis(1).Exponent = 0;
    catch
    end
    try
        axh.YAxis(2).Exponent = 0;
    catch
    end
end

function txt = local_format_sib(value)
%LOCAL_FORMAT_SIB Compact integrated score for the iteration table.
% Values below 1e5 retain the ordinary four-digit display convention.
% Larger values use one-decimal scientific notation without a plus sign or
% a leading zero in the exponent (for example, 6048029 -> 6.0e6).
    if isfinite(value) ...
            && abs(value) >= 1e5
        exponent = floor(log10(abs(value)));
        mantissa = round(value/10^exponent,1);
        if abs(mantissa) >= 10
            mantissa = mantissa/10;
            exponent = exponent + 1;
        end
        txt = sprintf('%.1fe%d',mantissa,exponent);
    else
        txt = local_format_metric(value);
    end
end

function txt = local_format_metric(value)
%LOCAL_FORMAT_METRIC Four displayed digits, excluding sign/decimal point.
    if isnan(value)
        txt = 'NaN';
    elseif isinf(value)
        if value > 0
            txt = 'Inf';
        else
            txt = '-Inf';
        end
    else
        magnitude = abs(value);
        if magnitude < 10
            nDecimals = 3;
        else
            nDecimals = max(0,3-floor(log10(magnitude)));
        end
        roundedMagnitude = abs(round(value,nDecimals));
        if roundedMagnitude >= 10
            nDecimals = max(0,3-floor(log10(roundedMagnitude)));
        end
        txt = sprintf('%.*f',nDecimals,value);
    end
end

function line = local_iteration_line(values,widths,dividerAfter,groupBefore)
%LOCAL_ITERATION_LINE Right-align fields with compact inter-column spacing.
    if nargin < 3
        dividerAfter = [];
    end
    if nargin < 4
        groupBefore = [];
    end
    n = numel(values);
    parts = cell(1,n);
    parts{1} = sprintf('%-*s',widths(1),char(values{1}));
    for k = 2:n
        value = char(values{k});
        if ~isempty(dividerAfter) && k > dividerAfter
            % Short metric names are centred over the five-character
            % positive-value footprint; numerical strings are left aligned.
            if any(strcmp(value, ...
                    {'nse','kge','jkge','S_fdc','S_p','S_logp'}))
                if numel(value) > 5
                    % Long FDC names use the complete field.  Reserving a
                    % sign column here clipped S_logp to S_log.
                    parts{k} = sprintf('%-*s',widths(k),value);
                else
                    % Reserve the first character for the sign used by the
                    % numerical rows, then center the heading over the usual
                    % five-character positive-value footprint.
                    parts{k} = [' ' local_center_field( ...
                        value,widths(k)-1,5)];
                end
            elseif strcmp(value,'ΣGLS')
                % Centre the four-character diagnostic heading over its
                % five-character numeric value (for example, 1.820).
                parts{k} = ['  ' local_center_field( ...
                    value,widths(k)-2,5)];
            else
                numericValue = str2double(value);
                if ~isnan(numericValue) && ~startsWith(value,{'-','+'})
                    % Reserve the first character for the sign so positive
                    % and negative values share the same decimal position.
                    value = sprintf(' %s',value);
                end
                parts{k} = sprintf('%-*s',widths(k),value);
            end
        else
            parts{k} = sprintf('%*s',widths(k),value);
        end
    end
    % Five-character metric values are right-aligned in seven-character
    % fields. Their two leading blanks plus this one structural separator
    % produce three visible spaces between adjacent positive entries.
    line = parts{1};
    for k = 2:n
        if ismember(k-1,dividerAfter)
            separator = '  | ';
        elseif ismember(k,groupBefore)
            separator = '  ';
        else
            separator = ' ';
        end
        line = [line separator parts{k}]; %#ok<AGROW>
    end
    line = [line newline];
end

function field = local_center_field(value,width,footprint)
%LOCAL_CENTER_FIELD Centre a header over the usual positive metric value.
    footprint = min(width,footprint);
    left = max(0,floor((footprint-numel(value))/2));
    field = [repmat(' ',1,left) value];
    field = [field repmat(' ',1,max(0,width-numel(field)))];
    field = field(1:width);
end

function line = local_iteration_group_lines( ...
    widths,nPrefix,nT,nS,lossPow,rssPow,useRSS)
%LOCAL_ITERATION_GROUP_LINES Two-level grouped header for iteration table.
    ordinaryGap = 1;
    objectiveGap = 4; % two blanks, vertical rule, one blank
    metricGroupGap = 2;
    tGroupEnd = nPrefix + nT;
    starts = zeros(size(widths));
    stops = zeros(size(widths));
    pos = 1;
    for k = 1:numel(widths)
        starts(k) = pos;
        stops(k) = pos + widths(k) - 1;
        if k == 2
            gap = objectiveGap;
        elseif k == tGroupEnd
            gap = metricGroupGap;
        else
            gap = ordinaryGap;
        end
        pos = stops(k) + gap + 1;
    end

    lineWidth = stops(end);
    top = repmat(' ',1,lineWidth);
    lower = repmat(' ',1,lineWidth);

    top(stops(2)-numel('loss')+1:stops(2)) = 'loss';

    % Major boundary between the optimized objective and diagnostics.
    objectiveRule = stops(2) + 3;
    top(objectiveRule) = '|';
    lower(objectiveRule) = '|';

    if lossPow ~= 0
        scaleLabel = sprintf('x10^%d',lossPow);
        lower(stops(2)-numel(scaleLabel)+1:stops(2)) = scaleLabel;
    end
    if useRSS && rssPow ~= 0
        lower = local_place_group_label(lower,starts,stops,3, ...
            sprintf('x10^%d',rssPow));
    end
    idxT = nPrefix + (1:nT);
    idxS = nPrefix + nT + (1:nS);
    % Place the group rules over the five-character numerical footprints,
    % rather than over the surrounding sign/alignment padding.
    leftT = starts(idxT(1))+1;
    rightT = stops(idxT(end))-1;
    leftS = starts(idxS(1))+1;
    % End the outer performance/S_ib rules one character earlier so the
    % closing bar aligns with the final character of S_fdc below it.
    rightS = stops(idxS(end));
    lower = local_draw_group_box(lower,leftT,rightT,'T_');
    lower = local_draw_group_box(lower,leftS,rightS,'S_ib');

    % The upper box shares the exact outer boundaries of the two lower
    % boxes, so all three vertical rules form one coherent header.
    metricLeft = leftT;
    metricRight = rightS;
    top(metricLeft:metricRight) = '-';
    top(metricLeft) = '|';
    top(metricRight) = '|';
    label = ' performance metrics ';
    first = ceil((metricLeft+metricRight-numel(label)+1)/2);
    top(first:first+numel(label)-1) = label;
    line = [deblank(top) newline deblank(lower) newline];
end

function chars = local_draw_group_box(chars,left,right,label)
%LOCAL_DRAW_GROUP_BOX Draw a group box with a centred, padded label.
    chars(left:right) = '-';
    chars(left) = '|';
    chars(right) = '|';
    paddedLabel = [' ' label ' '];
    first = floor((left+right-numel(paddedLabel)+1)/2);
    chars(first:first+numel(paddedLabel)-1) = paddedLabel;
end

function chars = local_place_group_label(chars,starts,stops,idx,label)
%LOCAL_PLACE_GROUP_LABEL Center one label over a contiguous column group.
    left = starts(idx(1));
    right = stops(idx(end));
    first = floor((left + right - numel(label) + 1) / 2);
    chars(first:first + numel(label) - 1) = label;
end

function flagAx = add_country_flag(fig,mdl)
%ADD_COUNTRY_FLAG Add the active region's packaged flag to diagnostics.
    flagAx = gobjects(0);
    if ~isfield(mdl,'region') ...
        || isempty(mdl.region)
        return
    end
    try
        regionCode = char(region_helpers('code',mdl.region));
        short = char(region_helpers('short',regionCode));
        if any(strcmpi(regionCode, ...
                {'CAMELS_US','CAMELSH_US','MACH_US'}))
            short = 'US';
        elseif any(strcmpi(regionCode, ...
                {'CAMELS_KR','CAMELSH_KR'}))
            short = 'KR';
        elseif strcmpi(regionCode,'CAMELS_DEH')
            short = 'DE';
        elseif strcmpi(regionCode,'BULL_ES')
            short = 'ES';
        end
        projectRoot = fileparts(fileparts(mfilename('fullpath')));
        flagFile = fullfile(projectRoot,'flags', ...
            [lower(short) '.png']);
        if ~isfile(flagFile)
            return
        end
        [rgb,~,alpha] = imread(flagFile);
        flagAx = axes('Parent',fig, ...
            'Units','normalized', ...
            'Position',[0.958 0.954 0.030 0.030], ...
            'Visible','off','Color','none', ...
            'HandleVisibility','off','HitTest','off');
        h = image(flagAx,rgb);
        if ~isempty(alpha)
            h.AlphaData = alpha;
        end
        h.HitTest = 'off';
        axis(flagAx,'image');
        axis(flagAx,'off');
    catch
        if ~isempty(flagAx) ...
                && isgraphics(flagAx)
            delete(flagAx);
        end
        flagAx = gobjects(0);
    end
end

function name = local_print_model_name(mdl)
    if isfield(mdl,'variant') ...
            && strcmpi(string(mdl.variant),'gchm_ode')
        name = 'gchm_ode';
    else
        name = sage_model_name(mdl.model);
    end
end
