function folder = sage_native_build_dir()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAGE_NATIVE_BUILD_DIR Platform-specific native output outside source
% folders.
%
% Copyright (c) 2026 Jasper A. Vrugt
% University of California, Irvine | jasper@uci.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    here = fileparts(mfilename('fullpath'));
    folder = fullfile(fileparts(fileparts(here)), ...
        'release_artifacts','native','Matlab',mexext);
    if ~isfolder(folder), mkdir(folder); end
    addpath(folder,'-begin');
end
