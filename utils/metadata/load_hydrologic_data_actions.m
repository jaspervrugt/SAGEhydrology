function [actions,note] = load_hydrologic_data_actions(region,resolution)
%LOAD_HYDROLOGIC_DATA_ACTIONS Read curated data corrections/exclusions.
%
% Actions are deliberately kept outside generic quality-screening code.
% The registry can therefore document exceptional source-data decisions
% for any region without making the common diagnostics region-specific.

    file = fullfile(fileparts(mfilename('fullpath')), ...
        'hydrologic_data_actions.csv');
    actions = local_empty_actions();
    note = "";
    if ~isfile(file)
        return
    end

    opts = detectImportOptions(file,'VariableNamingRule','preserve');
    textFields = {'region','resolution','gauge_id','variable', ...
        'action','date','reason','evidence'};
    for i = 1:numel(textFields)
        opts = setvartype(opts,textFields{i},'string');
    end
    opts = setvartype(opts,'factor','double');
    opts = setvartype(opts,'enabled','logical');
    T = readtable(file,opts);
    required = [textFields {'factor','enabled'}];
    missing = setdiff(required,T.Properties.VariableNames);
    if ~isempty(missing)
        error('SAGE:dataActions:BadRegistry', ...
            'Data-action registry is missing column(s): %s.', ...
            strjoin(missing,', '));
    end

    use = T.enabled ...
        & strcmpi(strip(T.region),strip(string(region))) ...
        & (strcmpi(strip(T.resolution),strip(string(resolution))) ...
        | strip(T.resolution) == "*");
    T = T(use,:);
    if isempty(T)
        return
    end

    actions = repmat(local_action_template(),height(T),1);
    for i = 1:height(T)
        action = lower(strip(T.action(i)));
        variable = strip(T.variable(i));
        gauge = local_normalize_gauge(T.gauge_id(i));
        if ~ismember(action,["exclude","scale","mask"])
            error('SAGE:dataActions:BadAction', ...
                'Unsupported action "%s" for %s gauge %s.', ...
                action,region,gauge);
        end
        factor = T.factor(i);
        if action == "scale" ...
                && (~isfinite(factor) ...
                || factor <= 0)
            error('SAGE:dataActions:BadFactor', ...
                ['Scale action for %s gauge ' ...
                '%s requires factor > 0.'], ...
                region,gauge);
        end
        if action ~= "scale"
            factor = NaN;
        end
        actionDate = strip(T.date(i));
        if ismissing(actionDate)
            actionDate = "";
        end
        if action == "mask" ...
                && (strlength(actionDate) == 0 ...
                || isnat(datetime(actionDate, ...
                'InputFormat','yyyy-MM-dd')))
            error('SAGE:dataActions:BadDate', ...
                'Mask action for %s gauge %s requires a valid date.', ...
                region,gauge);
        end
        actions(i).region = char(strip(T.region(i)));
        actions(i).resolution = char(strip(T.resolution(i)));
        actions(i).gauge_id = char(gauge);
        actions(i).variable = char(variable);
        actions(i).action = char(action);
        actions(i).factor = factor;
        actions(i).date = char(actionDate);
        actions(i).reason = char(strip(T.reason(i)));
        actions(i).evidence = char(strip(T.evidence(i)));
        actions(i).registry_file = file;
    end
    note = local_note(actions);
end

function actions = local_empty_actions()
    actions = repmat(local_action_template(),0,1);
end

function action = local_action_template()
    action = struct('region','','resolution','','gauge_id','', ...
        'variable','','action','','date','','factor',NaN,'reason','', ...
        'evidence','','registry_file','');
end

function gauge = local_normalize_gauge(value)
    gauge = upper(strip(string(value)));
    if ~isempty(regexp(char(gauge),'^\d+(?:\.0+)?$','once'))
        gauge = string(sprintf('%d',str2double(gauge)));
    end
end

function note = local_note(actions)
    lines = strings(numel(actions),1);
    for i = 1:numel(actions)
        a = actions(i);
        if strcmpi(a.action,'exclude')
            lines(i) = sprintf('Excluded gauge %s: %s.', ...
                a.gauge_id,a.reason);
        elseif strcmpi(a.action,'mask')
            lines(i) = sprintf(['Masked %s for gauge %s on %s: ' ...
                '%s.'],a.variable,a.gauge_id,a.date,a.reason);
        else
            lines(i) = sprintf(['Scaled %s for gauge %s by %.12g: ' ...
                '%s.'],a.variable,a.gauge_id,a.factor,a.reason);
        end
    end
    note = strjoin(lines,newline);
end
