function dat = attach_cart_bands(dat,bas,region,dirQ)
%ATTACH_CART_BANDS Add basin hypsometry to the shared MATLAB inputs.

    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'projects'));
    folder = char(dirQ);
    while true
        [parent,leaf] = fileparts(folder);
        if strcmpi(leaf,region) || ...
            isfile(fullfile(folder,'cartplus_terrain_bands.csv')) || ...
            isfile(fullfile(folder,'elev_bands671.npz'))
            break
        end
        if isempty(parent) || strcmp(parent,folder)
            error('CaRT:DataRoot','Cannot find the data root for %s.',region);
        end
        folder = parent;
    end
    C = struct('data_root',folder,'region',region);
    C.gauges = bas.id_gauge;
    bands = read_cartplus_project(C,true);
    for i = 1:numel(dat)
        dat{i}.cart_bands.elevation = bands.elevation(:,i);
        dat{i}.cart_bands.area = bands.area(:,i);
        dat{i}.cart_bands.floor = 1e-4*strcmp(region,'CAMELS_AU');
    end
end
