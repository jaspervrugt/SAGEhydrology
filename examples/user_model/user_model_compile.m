function flag = user_model_compile(verbose)
%USER_MODEL_COMPILE Compile the standalone SAGE user-model MEX plug-in.
%
% SYNOPSIS:
%   flag = user_model_compile()
%   flag = user_model_compile(verbose)
%
% INPUT:
%   verbose     Optional logical/scalar flag.
%               0 = normal MEX output [default]
%               1 = verbose MEX output
%
% OUTPUT:
%   flag        Compilation status:
%                1 = successfully compiled user_model.cpp
%               -1 = compilation failed or no C++ compiler configured
%               -2 = user_model.cpp not found
%
% DESCRIPTION:
%   The MEX is built from crr_user_model_mex.cpp, user_model_prepare.cpp,
%   and user_model.cpp. It remains independent of crr_model_mex so a
%   precompiled user model can be installed with a deployed SAGE GUI.
%
%   Platform-specific optimization flags are used for Windows, macOS, and
%   Linux/Unix. If custom flags fail, the function can be extended to retry
%   with MATLAB default MEX flags.
%
% REQUIREMENTS:
%   A configured MATLAB-supported C++ MEX compiler:
%
%       mex -setup C++
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, May 2026                                  %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    verbose = 0;
end

flag = 1;

cc = mex.getCompilerConfigurations('C++','Selected');

if isempty(cc)
    flag = -1;
    warning('user_model_compile:noCompiler', ...
        'No C++ MEX compiler configured. Run: mex -setup C++');
    return
end

this_dir = fileparts(mfilename('fullpath'));
cppFiles = { ...
    fullfile(this_dir,'crr_user_model_mex.cpp'), ...
    fullfile(this_dir,'user_model_prepare.cpp'), ...
    fullfile(this_dir,'user_model.cpp')};

for k = 1:numel(cppFiles)
    if ~isfile(cppFiles{k})
        flag = -2;
        warning('user_model_compile:missingCppFile', ...
            'Could not find cpp file: %s',cppFiles{k});
        return
    end
end

mexName = 'crr_user_model';
mexOut = fullfile(this_dir,[mexName '.' mexext]);

oldFiles = {mexOut};
if ispc
    oldFiles = [oldFiles,{ ...
        fullfile(this_dir,[mexName '.lib']), ...
        fullfile(this_dir,[mexName '.exp']), ...
        fullfile(this_dir,[mexName '.obj'])}];
else
    oldFiles = [oldFiles,{fullfile(this_dir,[mexName '.o'])}];
end
for k = 1:numel(oldFiles)
    if isfile(oldFiles{k})
        delete(oldFiles{k});
    end
end

try
    clear mex
    rehash toolboxcache

    mexArgs = [{"-R2018a","-O", ...
        "-outdir",this_dir, ...
        "-output",mexName}, cppFiles];

    if verbose
        mexArgs = [{"-v"}, mexArgs];
    end

    if ispc
        if contains(lower(cc.Name),'mingw')
            mexArgs = [mexArgs, ...
                {"CXXFLAGS=$CXXFLAGS -O2 -DNDEBUG"}];
        else
            mexArgs = [mexArgs, ...
                {"COMPFLAGS=$COMPFLAGS /O2 " + ...
                "/DNDEBUG /fp:fast /GL /Gw /Gy", ...
                 "CXXFLAGS=$CXXFLAGS /O2 /DNDEBUG /EHsc", ...
                 "LINKFLAGS=$LINKFLAGS /LTCG /OPT:REF /OPT:ICF"}];
        end
    elseif ismac
        mexArgs = [mexArgs, ...
            {"CXXFLAGS=$CXXFLAGS -O3 -DNDEBUG -std=c++11", ...
             "LDFLAGS=$LDFLAGS -Wl,-dead_strip"}];
    elseif isunix
        mexArgs = [mexArgs, ...
            {"CXXFLAGS=$CXXFLAGS -O3 -DNDEBUG -std=c++11"}];
    end
    mex(mexArgs{:});

    if ~isfile(mexOut)
        flag = -1;
        warning('user_model_compile:noMexOutput', ...
            'Compilation did not produce %s',mexOut);
    end

catch ME
    flag = -1;
    warning('user_model_compile:compileFailed', ...
        'Compilation failed for crr_user_model:\n%s',ME.message);
end
end
