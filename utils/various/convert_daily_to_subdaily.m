function [th_min,th_max,par_units] = ...
    convert_daily_to_subdaily(th_min,th_max,par_units,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONVERT_DAILY_TO_SUBDAILY Convert daily parameter bounds to model-step
%units for a subdaily simulation.
%
% SYNOPSIS:
%   [th_min,th_max,par_units] = ...
%       convert_daily_to_subdaily(th_min,th_max,par_units,dt)
%
% INPUTS:
%   th_min      d x 1 lower parameter bounds in daily units
%   th_max      d x 1 upper parameter bounds in daily units
%   par_units   1 x d or d x 1 cell/string array of parameter units
%   dt          number of model time steps per day
%               1  = daily
%               24 = hourly
%               96 = 15-minute
%
% OUTPUTS:
%   th_min      converted lower bounds
%   th_max      converted upper bounds
%   par_units   converted unit labels
%
% CONVERSION RULES:
%   Rates per day are divided by dt:
%       1/d, mm/d, mm/d/degC, mm/degC/d
%
%   Durations in days are multiplied by dt:
%       d -> h for dt = 24
%       d -> 15 min for dt = 96
%       d -> time step for other subdaily dt values
%
%   Parameters without an explicit daily time unit remain unchanged.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    validateattributes(dt,{'numeric'}, ...
        {'scalar','real','finite','positive'}, ...
        mfilename,'dt',4);

    if abs(dt - round(dt)) > 1e-10
        error('convert_daily_to_subdaily:nonIntegerDt', ...
            ['dt must be the integer number of model ' ...
            'time steps per day. Received dt = %g.'],dt);
    end
    dt = double(round(dt));

    if dt < 1
        error('convert_daily_to_subdaily:badDt', ...
            'dt must be at least 1.');
    end

    if ~isnumeric(th_min) ...
            || ~isnumeric(th_max) ...
            || ~isvector(th_min) ...
            || ~isvector(th_max) ...
            || numel(th_min) ~= numel(th_max)
        error('convert_daily_to_subdaily:badBounds', ...
            ['th_min and th_max must be numeric vectors ' ...
            'with equal length.']);
    end

    wasString = isstring(par_units);
    wasRow = isrow(par_units);
    units = cellstr(string(par_units(:)));

    if numel(units) ~= numel(th_min)
        error('convert_daily_to_subdaily:unitCountMismatch', ...
            ['The number of parameter units (%d) must equal ' ...
            'the number of parameter bounds (%d).'], ...
            numel(units),numel(th_min));
    end

    th_min = double(th_min);
    th_max = double(th_max);

    % Daily data require no conversion.
    if dt == 1
        par_units = local_restore_unit_shape(units,wasString,wasRow);
        return
    end

    for j = 1:numel(units)
        uOriginal = strtrim(units{j});
        u = local_normalize_unit(uOriginal);

        if local_is_inverse_day_unit(u)
            th_min(j) = th_min(j) / dt;
            th_max(j) = th_max(j) / dt;
            units{j} = local_replace_inverse_day_unit(uOriginal,dt);

        elseif local_is_day_rate_unit(u)
            th_min(j) = th_min(j) / dt;
            th_max(j) = th_max(j) / dt;
            units{j} = local_replace_day_rate_unit(uOriginal,dt);

        elseif local_is_day_duration_unit(u)
            th_min(j) = th_min(j) * dt;
            th_max(j) = th_max(j) * dt;
            units{j} = local_subdaily_duration_label(dt);
        end
    end

    par_units = local_restore_unit_shape(units,wasString,wasRow);
end

function tf = local_is_inverse_day_unit(u)
    tf = any(strcmp(u,{ ...
        '1/d','d^-1','d-1','day^-1','day-1','1/day'}));
end

function tf = local_is_day_rate_unit(u)
    % Match explicit '/d' or '/day' factors anywhere in the unit.
    tf = contains(u,'/d') ...
        || contains(u,'/day');
end

function tf = local_is_day_duration_unit(u)
    tf = any(strcmp(u,{'d','day','days'}));
end

function u = local_normalize_unit(u)
    u = lower(strtrim(char(u)));
    u = strrep(u,' ', '');
    u = strrep(u,'°c','degc');
    u = strrep(u,'degrees c','degc');
    u = strrep(u,'degreec','degc');
end

function out = local_replace_inverse_day_unit(~,dt)
    switch dt
        case 24
            out = '1/h';
        case 96
            out = '1/15 min';
        otherwise
            out = '1/time step';
    end
end

function out = local_replace_day_rate_unit(unitIn,dt)
    out = char(unitIn);

    switch dt
        case 24
            token = 'h';
        case 96
            token = '15 min';
        otherwise
            token = 'time step';
    end

    % Replace only the daily denominator. Preserve the remaining unit order.
    out = regexprep(out,'/days?(?=($|/))',['/' token], ...
        'ignorecase');
    out = regexprep(out,'/d(?=($|/))',['/' token], ...
        'ignorecase');

    % Handle inverse-day notation that may appear inside a compound unit.
    out = regexprep(out,'day\^?-?1',token,'ignorecase');
    out = regexprep(out,'d\^?-?1',token,'ignorecase');
end

function out = local_subdaily_duration_label(dt)
    switch dt
        case 24
            out = 'h';
        case 96
            out = '15 min';
        otherwise
            out = 'time steps';
    end
end

function unitsOut = local_restore_unit_shape(units,wasString,wasRow)
    if wasString
        unitsOut = string(units);
    else
        unitsOut = units;
    end

    if wasRow
        unitsOut = reshape(unitsOut,1,[]);
    else
        unitsOut = reshape(unitsOut,[],1);
    end
end
