function mexFile = compile_crr_model_mex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPILE_CRR_MODEL_MEX Build the unified CRR_MODEL MEX gateway.
%
%  Links the native cores for HYMOD, HMODEL, SAC-SMA, Xinanjiang,
%  GR4J-A+, HBV, and CFE-NWM into one MEX file.
%
% SYNOPSIS:
%   mexFile = compile_crr_model_mex
%
% INPUT ARGUMENTS:
%   (none)          no input arguments
%
% OUTPUT ARGUMENTS:
%   mexFile         path to the compiled unified MEX file
%
% NOTES:
%   Model 99 (user_model) remains a separately compiled plug-in.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    here = fileparts(mfilename('fullpath'));
    native = sage_native_root();
    core = fullfile(native,'models');
    output = sage_native_build_dir();

    old = pwd;
    cleanupObj = onCleanup(@() cd(old));
    cd(here);

    sources = {
        fullfile(native,'adapters','matlab','crr_model_mex.cpp')
        fullfile(core,'hymod','hymod.cpp')
        fullfile(core,'hmodel','hmodel.cpp')
        fullfile(core,'sacsma','sacsma.cpp')
        fullfile(core,'Xinanjiang','xinanjiang.cpp')
        fullfile(core,'gr4jA','gr4jA.cpp')

        fullfile(core,'hbv','hbv.cpp')
        fullfile(core,'cfe_nwm','cfe_nwm.cpp')

        };

    % Check compiler
    cfg = mex.getCompilerConfigurations('C++','Selected');

    if isempty(cfg)
        error('compile_crr_model_mex:NoCompiler', ...
            'No selected C++ MEX compiler.');
    end

    fprintf('Build folder: %s\n',here);
    fprintf('Compiler: %s\n',cfg.Name);

    % Check sources
    for j = 1:numel(sources)
        if ~isfile(sources{j})
            error('compile_crr_model_mex:MissingSource', ...
                'Missing source file:\n%s',sources{j});
        end
    end

    fprintf(['Building direct-core CRR_MODEL MEX: ' ...
        'HYMOD, HMODEL, SAC-SMA, Xinanjiang, ' ...
        'GR4J-A, HBV, CFE-NWM...\n']);

    mex( ...
        '-R2018a', ...
        '-O', ...
        ['-I' core], ...
        sources{:}, ...
        '-outdir',output, ...
        '-output','crr_model_mex');

    mexFile = fullfile(output, ...
        ['crr_model_mex.' mexext]);

    if ~isfile(mexFile)
        error('compile_crr_model_mex:BuildFailed', ...
            'MEX output was not created:\n%s',mexFile);
    end

    fprintf('Created: %s\n',mexFile);

end
