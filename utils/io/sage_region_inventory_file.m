function file = sage_region_inventory_file(name,root,bundleRoot)
% Resolve packaged basin inventories independently of compiled code paths.
    if nargin<3
        bundleRoot='';
        if isdeployed, bundleRoot=ctfroot; end
    end
    roots={fullfile(root,'regions')};
    if ~isempty(bundleRoot)
        roots=[{bundleRoot},roots];
    end
    for k=1:numel(roots)
        matches=dir(fullfile(roots{k},'**',name));
        matches=matches(~[matches.isdir]);
        if ~isempty(matches)
            file=fullfile(matches(1).folder,matches(1).name);
            return
        end
    end
    error('SAGE:MissingRegionInventory', ...
        'Cannot find the packaged basin inventory %s.',name);
end
