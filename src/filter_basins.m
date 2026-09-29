function [dat,A,latlon,bas] = filter_basins( ...
    dat,A,latlon,bas,eligibility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FILTER_BASINS Remove data-ineligible basins.
%
%  Retain eligible training basins before eligible evaluation basins and
%  compact basin-indexed data, attributes, coordinates, and metadata.
%
% SYNOPSIS:
%   [dat,A,latlon,bas] = filter_basins( ...
%       dat,A,latlon,bas,eligibility)
%
% INPUT ARGUMENTS:
%   dat             K-by-1 cell array of basin time-series structures
%   A               r-by-K static attributes; may be empty
%   latlon          K-by-2 basin coordinates; may be empty
%   bas             selected-basin metadata
%    .K              total requested basins
%    .K_t            training basins, listed first
%    .K_e            evaluation basins
%   eligibility     basin eligibility returned by check_basins
%    .valid          K-by-1 logical retention mask
%
% OUTPUT ARGUMENTS:
%   dat             time-series entries for eligible basins
%   A               attributes for eligible basins
%   latlon          coordinates for eligible basins
%   bas             compacted metadata and exclusion audit
%    .K              retained basin count
%    .K_t            retained training count
%    .K_e            retained evaluation count
%    .data_quality   requested/active/excluded basin audit
%
% NOTES:
%   Basin ordering is preserved within each training/evaluation group.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oldK = double(bas.K);
    oldKt = min(double(bas.K_t),oldK);
    valid = logical(eligibility.valid(:));
    if numel(valid) ~= oldK
        error('SAGE:dataQuality:EligibilitySize', ...
            ['Eligibility vector length ' ...
            'does not match bas.K.']);
    end
    if numel(dat) ~= oldK
        error('SAGE:dataQuality:DataSize', ...
            ['Number of dat entries ' ...
            'does not match bas.K.']);
    end

    keepT = find(valid(1:oldKt));
    keepE = oldKt + find(valid(oldKt+1:oldK));
    keep = [keepT(:); keepE(:)];
    drop = find(~valid);

    requestedId = string((1:oldK).');
    if isfield(bas,'id_gauge') ...
            && numel(bas.id_gauge) >= oldK
        requestedId = string(bas.id_gauge(1:oldK));
    end
    bas.data_quality = struct( ...
        'requested_K',oldK, ...
        'requested_K_t',oldKt, ...
        'requested_K_e',oldK-oldKt, ...
        'requested_id',requestedId, ...
        'active_id',requestedId(keep), ...
        'excluded_id',requestedId(drop), ...
        'excluded_reason',eligibility.reason(drop), ...
        'forcing_complete',eligibility.forcing_complete, ...
        'discharge_complete',eligibility.discharge_complete, ...
        'coverage_q_train',eligibility.coverage_q_train, ...
        'coverage_q_eval',eligibility.coverage_q_eval, ...
        'runoff_ratio_train',eligibility.runoff_ratio_train, ...
        'runoff_ratio_eval',eligibility.runoff_ratio_eval, ...
        'hydrologic_alert',eligibility.hydrologic_alert, ...
        'documented_exclusion',eligibility.documented_exclusion, ...
        'minimum_q_coverage',eligibility.minimum_q_coverage);

    dat = dat(keep);
    if ~isempty(A)
        if size(A,2) ~= oldK
            error('SAGE:dataQuality:AttributeSize', ...
                ['Number of attribute columns ' ...
                'does not match bas.K.']);
        end
        A = A(:,keep);
    end
    if ~isempty(latlon)
        if size(latlon,1) ~= oldK
            error('SAGE:dataQuality:CoordinateSize', ...
                ['Number of coordinate rows ' ...
                'does not match bas.K.']);
        end
        latlon = latlon(keep,:);
    end

    for name = {'id_gauge','gname'}
        field = name{1};
        if isfield(bas,field) ...
                && numel(bas.(field)) == oldK
            value = bas.(field);
            bas.(field) = value(keep);
        end
    end
    if isfield(bas,'zone')
        bas.zone = local_subset_basin_value(bas.zone, ...
            keep,oldK);
    end

    bas.K_t = numel(keepT);
    bas.K_e = numel(keepE);
    bas.K = numel(keep);
    bas.id_t = (1:bas.K_t).';
    bas.id_e = (bas.K_t+1:bas.K).';
    bas.id_plot = (1:bas.K).';
end

function value = local_subset_basin_value(value,keep,K)
%LOCAL_SUBSET_BASIN_VALUE Subset vector, table, or structure zone formats.
    if istable(value) ...
            && height(value) == K
        value = value(keep,:);
    elseif (isvector(value) ...
            || isstring(value) ...
            || iscell(value)) ...
            && numel(value) == K
        value = value(keep);
    elseif isnumeric(value) ...
            && size(value,1) == K
        value = value(keep,:);
    elseif isstruct(value)
        fields = fieldnames(value);
        for j = 1:numel(fields)
            field = fields{j};
            item = value.(field);
            if istable(item) ...
                    && height(item) == K
                value.(field) = item(keep,:);
            elseif (isvector(item) ...
                    || isstring(item) ...
                    || iscell(item)) ...
                    && numel(item) == K
                value.(field) = item(keep);
            elseif isnumeric(item) ...
                    && size(item,1) == K
                value.(field) = item(keep,:);
            end
        end
    end
end
