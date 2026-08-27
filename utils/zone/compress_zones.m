function zone = compress_zones(zone,max_zones)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPRESS_ZONES Merge small hydroclimatic zones into a compact zone set.
%
% SYNOPSIS:
%   zone = compress_zones(zone,max_zones)
%
% INPUT:
%   zone        Structure with hydroclimatic zone information
%    .id          K x 1 string array with zone labels for each basin
%    .num         K x 1 numeric array with integer zone numbers
%    .names       m x 1 string array with unique zone names
%
%   max_zones   Maximum number of zones retained after compression
%
% OUTPUT:
%   zone        Updated zone structure with compressed zone labels,
%               renumbered zone IDs, and updated zone names
%
% DESCRIPTION:
%   This function reduces the number of hydroclimatic zones to at most
%   max_zones by merging low-priority or small zones into larger existing
%   zones. Output zone structure preserves one zone assignment per basin
%   while ensuring that zone.num contains compact integer identifiers.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, June 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    id = string(zone.id(:));   % FORCE SAFE TYPE
    
    % ensure no categorical/cell leakage
    id = string(cellstr(id));
    
    while numel(unique(id)) > max_zones
    
        u = unique(id,'stable');
    
        % counts
        c = zeros(numel(u),1);
        for i = 1:numel(u)
            c(i) = sum(id == u(i));
        end
    
        % smallest zone
        [~,iSmall] = min(c);
        smallName = u(iSmall);
    
        ixSmall = id == smallName;
    
        % candidates (exclude small itself)
        cand = u;
        cand(iSmall) = [];
    
        % FIX: avoid max() on strings
        % choose largest zone by count instead
        c2 = zeros(numel(cand),1);
        for i = 1:numel(cand)
            c2(i) = sum(id == cand(i));
        end
    
        [~,iBig] = max(c2);
        target = cand(iBig);
    
        % merge
        id(ixSmall) = target;
    end
    
    % re-index
    u = unique(id,'stable');
    newNum = nan(size(id));
    
    for i = 1:numel(u)
        newNum(id == u(i)) = i;
    end
    
    zone.id = id;
    zone.num = newNum;
    
    zone.codes = u;
    zone.names = u;
    zone.long  = u;

end
