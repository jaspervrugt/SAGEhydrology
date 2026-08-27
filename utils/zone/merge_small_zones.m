function zone = merge_small_zones(zone,min_per_zone)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MERGE_SMALL_ZONES Merge hydroclimatic zones with too few basins.
%
% SYNOPSIS:
%   zone = merge_small_zones(zone,min_per_zone)
%
% INPUT:
%   zone          Structure with hydroclimatic zone information
%    .id            K x 1 string array with zone labels for each basin
%    .num           K x 1 numeric array with integer zone numbers
%    .names         m x 1 string array with unique zone names
%
%   min_per_zone  Minimum number of basins required for a zone to be kept
%
% OUTPUT:
%   zone          Updated zone structure in which zones with fewer than
%                 min_per_zone basins have been merged into larger zones.
%                 Zone numbers are compacted and zone.names is updated.
%
% DESCRIPTION:
%   This helper enforces a minimum sample size per hydroclimatic zone.
%   Small zones are reassigned to compatible larger zones so that
%   zone-stratified basin sampling remains stable and avoids creating
%   sparsely populated classes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, June 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 || isempty(min_per_zone)
        min_per_zone = 5;
    end

    % The GUI can occasionally pass a non-scalar numeric value here.
    % Zone merging requires one scalar minimum count.
    min_per_zone = double(min_per_zone);
    min_per_zone = min_per_zone(1);
    if ~isfinite(min_per_zone) || min_per_zone < 1
        min_per_zone = 5;
    end
    min_per_zone = max(1,round(min_per_zone));

    id = string(zone.id(:));
    K = numel(id);
    maxZones = max(1,floor(K ./ min_per_zone));

    X = [ ...
        zone.aridity(:), ...
        zone.elev_mean(:), ...
        zone.frac_snow(:), ...
        abs(zone.lat(:)) ];

    for j = 1:size(X,2)
        x = X(:,j);
        mu = mean(x,'omitnan');
        sg = std(x,'omitnan');

        if ~isfinite(sg) ...
                || sg == 0
            sg = 1;
        end

        X(:,j) = (x - mu) ./ sg;
    end

    while true

        names = unique(id,'stable');
        counts = zeros(numel(names),1);

        for j = 1:numel(names)
            counts(j) = sum(id == names(j));
        end

        if numel(names) <= 1
            break
        end

        if all(counts >= min_per_zone) ...
                && numel(names) <= maxZones
            break
        end

        % Merge the smallest class
        [~,js] = min(counts);

        smallName = names(js);
        ixSmall = id == smallName;

        cand = setdiff(1:numel(names),js,'stable');

        % Prefer valid target zones, but only if available
        bigCand = cand(counts(cand) >= min_per_zone);
        if ~isempty(bigCand)
            cand = bigCand;
        end

        xs = mean(X(ixSmall,:),1,'omitnan');

        d = inf(numel(cand),1);
        for q = 1:numel(cand)
            ixC = id == names(cand(q));
            xc = mean(X(ixC,:),1,'omitnan');

            ok = isfinite(xs) & isfinite(xc);
            if any(ok)
                d(q) = sqrt(sum((xs(ok) - xc(ok)).^2));
            end
        end

        [~,best] = min(d);
        parentName = names(cand(best));

        id(ixSmall) = parentName;

    end

    names = unique(id,'stable');
    newNum = nan(K,1);

    for j = 1:numel(names)
        newNum(id == names(j)) = j;
    end

    zone.id = id(:);
    zone.num = newNum(:);
    zone.names = names(:);
    zone.codes = names(:);
    zone.long = names(:);

    % fprintf('\nFinal merged zones:\n');
    % u = unique(zone.id,'stable');
    % for j = 1:numel(u)
    %     fprintf('%2d  %-25s  %3d\n', ...
    %         j,u(j),sum(zone.id == u(j)));
    % end
    % fprintf('Total zones = %d; max allowed = %d\n\n', ...
    %     numel(u),floor(numel(zone.id)/min_per_zone));

end
