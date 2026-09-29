function zoneTbl = sage_zone_summary(bas,curr,mdl,prd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZONE_SUMMARY_SAGE Summarize basin skill by zone and scenario.
%
%  Aggregate current per-basin metrics without figure postprocessing. 
%  Return one row per nonempty zone/scenario group.
%
% SYNOPSIS:
%   zoneTbl = sage_zone_summary(bas,curr,mdl,prd)
%
% INPUT ARGUMENTS:
%   bas             basin and zone information
%    .zone           zone metadata
%     .num            zone number for each selected basin
%    .K_t            number of training basins
%    .K_e            number of evaluation basins
%   curr            current basin-wise NSE, KGE, JKGE, and FDC scores
%   mdl             model identification and assessment settings
%   prd             period and temporal-resolution information
%
% OUTPUT ARGUMENTS:
%   zoneTbl         table of zone/scenario summaries
%    .n              basins with a finite score
%    .T_*            median skill scores by metric
%    .Sib_*          integrated basin loss scores by metric
%
% NOTES:
%   The table is empty when required zone or metric data are absent.
%   Scenario tags are tt, te, et, and ee.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    names = {'Zone','ZoneName','Scenario','ScenarioName','Model','n', ...
        'T_NSE','T_KGE','T_JKGE','T_S_fdc', ...
        'Sib_NSE','Sib_KGE','Sib_JKGE','Sib_S_fdc'};
    types = {'double','string','string','string','string','double', ...
        'double','double','double','double', ...
        'double','double','double','double'};
    zoneTbl = table('Size',[0 numel(names)], ...
        'VariableTypes',types,'VariableNames',names);

    if ~isstruct(bas) ...
            || ~isfield(bas,'zone') ...
            || ~isfield(bas.zone,'num') ...
            || isempty(bas.zone.num) ...
            || ~isfield(bas,'K_t') ...
            || ~isfield(bas,'K_e') ...
            || ~isstruct(curr) ...
            || ~isfield(curr,'NSE')
        return
    end

    zoneNo = double(bas.zone.num(:));
    zoneIDs = unique(zoneNo(isfinite(zoneNo) & zoneNo > 0));
    if isempty(zoneIDs)
        return
    end
    zoneIDs = (1:max(zoneIDs)).';
    zoneNames = localZoneNames(bas);
    modelName = localModelName(mdl);
    [resolution,periodName] = localPeriodWords(mdl,prd);

    tags = {'tt','te','et','ee'};
    groups = {'training','training', ...
        'evaluation','evaluation'};
    periods = {'training','evaluation', ...
        'training','evaluation'};
    groupIDs = {1:bas.K_t,1:bas.K_t, ...
        bas.K_t+(1:bas.K_e),bas.K_t + (1:bas.K_e)};
    rows = cell(0,numel(names));
    for iz = 1:numel(zoneIDs)
        zid = zoneIDs(iz);
        if iz <= numel(zoneNames)
            zoneName = zoneNames(iz);
        else
            zoneName = "Z" + string(zid);
        end
        for is = 1:numel(tags)
            ids = groupIDs{is};
            ids = ids(ids >= 1 ...
                & ids <= numel(zoneNo));
            if isempty(ids)
                continue
            end
            inZone = zoneNo(ids) == zid;
            if ~any(inZone)
                continue
            end
            tag = tags{is};
            nse = localMetric(curr,'NSE',tag,inZone);
            kge = localMetric(curr,'KGE',tag,inZone);
            jkge = localMetric(curr,'JKGE',tag,inZone);
            sfdc = localMetric(curr,'S_fdc',tag,inZone);
            nValid = max([nnz(isfinite(nse)),nnz(isfinite(kge)), ...
                nnz(isfinite(jkge)),nnz(isfinite(sfdc))]);
            scenarioName = sprintf('%s %s basin | %s %s', ...
                resolution,groups{is},periods{is},periodName);
            rows(end+1,:) = {double(zid),zoneName,string(tag), ...
                string(scenarioName),modelName,double(nValid), ...
                median(nse,'omitnan'),median(kge,'omitnan'), ...
                median(jkge,'omitnan'),median(sfdc,'omitnan'), ...
                localIntegratedLoss(nse),localIntegratedLoss(kge), ...
                localIntegratedLoss(jkge), ...
                localIntegratedLoss(sfdc)}; %#ok<AGROW>
        end
    end
    if ~isempty(rows)
        zoneTbl = cell2table(rows,'VariableNames',names);
    end
end

function values = localMetric(curr,name,tag,inZone)

    values = nan(0,1);
    if ~isfield(curr,name) ...
            || ~isstruct(curr.(name)) ...
            || ~isfield(curr.(name),tag)
        return
    end
    metric = curr.(name).(tag);
    if isempty(metric)
        return
    end
    metric = double(metric(:));
    if numel(metric) ~= numel(inZone)
        return
    end
    values = metric(inZone);
end

function score = localIntegratedLoss(values)

    values = values(isfinite(values));
    if isempty(values)
        score = NaN;
    else
        score = mean(1-values);
    end
end

function names = localZoneNames(bas)

    names = strings(0,1);
    if isfield(bas.zone,'names_short') ...
            && ~isempty(bas.zone.names_short)
        names = string(bas.zone.names_short(:));
    elseif isfield(bas.zone,'names') ...
            && ~isempty(bas.zone.names)
        names = string(bas.zone.names(:));
        names = strtrim(strrep(names,'_',' '));
    end
end

function name = localModelName(mdl)

    name = "SAGE";
    if ~isstruct(mdl) ...
            || ~isfield(mdl,'model') ...
            || isempty(mdl.model)
        return
    end
    modelID = double(mdl.model(1));
    try
        if isfield(mdl,'names') ...
                && ~isempty(mdl.names) ...
                && modelID >= 1 ...
                && modelID <= numel(mdl.names) ...
                && ~ismember(modelID,[11 99])
            if iscell(mdl.names)
                name = string(mdl.names{modelID});
            else
                name = string(mdl.names(modelID));
            end
        else
            name = string(sage_model_name(modelID));
        end
    catch
        name = "model_" + string(modelID);
    end
    if strcmpi(name,"xinanjiang")
        name = "Xinanjiang";
    end
end

function [resolution,periodName] = localPeriodWords(mdl,prd)

    resolution = 'daily';
    if isstruct(prd) ...
            && isfield(prd,'dt')
        if prd.dt == 24
            resolution = 'hourly';
        elseif prd.dt == 96
            resolution = '15-minute';
        end
    end
    
    periodName = 'period';
    if isstruct(mdl) ...
            && isfield(mdl,'sp_method') ...
            && ~isempty(mdl.sp_method)
        method = lower(string(mdl.sp_method));
        if method == "rainfall_block"
            periodName = 'split';
        elseif ~ismember(method, ...
                ["manual","traditional_block","block", ...
                 "deterministic_block"])
            periodName = 'mask';
        end
    end
end
