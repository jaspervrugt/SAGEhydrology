function root = sage_native_root()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_NATIVE_ROOT Locate the single shared C++ source tree.
%
% MATLAB gateways and Python bindings compile the same Native/models files.
% Keep Native, Matlab and Python together when building from source.
% Packaged applications use prebuilt native binaries instead.
%
% Copyright (c) 2026 Jasper A. Vrugt
% University of California, Irvine | jasper@uci.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    here = fileparts(mfilename('fullpath'));
    root = fullfile(fileparts(fileparts(here)),'Native');
    assert(isfile(fullfile(root,'models','hbv','hbv_routing.hpp')), ...
        'SAGE:MissingNativeSources', ...
        'Shared Native sources are missing: %s',root);
end
