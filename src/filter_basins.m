function [dat,A,latlon,bas] = filter_basins( ...
    dat,A,latlon,bas,eligibility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FILTER_BASINS Remove data-ineligible basins from a SAGE population
%
% SYNOPSIS:
%   [dat,A,latlon,bas] = filter_basins( ...
%       dat,A,latlon,bas,eligibility)
%
% INPUTS:
%   dat          K-by-1 cell array with basin time-series structures
%   A            r-by-K matrix of static basin attributes
%   latlon       K-by-2 basin-coordinate matrix (may be empty)
%   bas          basin structure; training basins must precede evaluation
%                basins and .K, .K_t, and .K_e must be defined
%   eligibility output structure returned by check_basins
%
% OUTPUTS:
%   dat,A,latlon basin-indexed inputs compacted to eligible basins
%   bas          compacted basin metadata and updated K, K_t, K_e, id_t,
%                id_e, and id_plot. The requested population and exclusion
%                audit are retained in bas.data_quality.
%
% DESCRIPTION:
%   Eligible training basins are retained first, followed by eligible
%   evaluation basins. This preserves the ordering assumed throughout
%   SAGE. Known basin-indexed metadata fields (id_gauge, gname, and zone)
%   are subset consistently. Scalar basin configuration remains unchanged.
%
% SEE ALSO: check_basins, prep_stats, init_args
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
