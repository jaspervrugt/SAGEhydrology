function mexFile = compile_crr_model_mex
%COMPILE_CRR_MODEL_MEX Build direct-core central CRR_MODEL MEX.
%
% Models 1--7 are linked directly against their extracted native cores:
%   1 HYMOD
%   2 HMODEL
%   3 SAC-SMA
%   4 Xinanjiang
%   5 GR4J-A+
%   6 HBV
%   7 CFE-NWM
%
% GR4J-B (model 21) is linked as a validated built-in core.
% user_model (model 8) remains an independently compiled plug-in.


    here = fileparts(mfilename('fullpath'));

    old = pwd;
    cleanupObj = onCleanup(@() cd(old));
    cd(here);

    sources = {
        fullfile(here,'crr_model_mex.cpp')
        fullfile(here,'hymod','hymod.cpp')
        fullfile(here,'hmodel','hmodel.cpp')
        fullfile(here,'sacsma','sacsma.cpp')
        fullfile(here,'Xinanjiang','xinanjiang.cpp')
        fullfile(here,'gr4jA','gr4jA.cpp')
        fullfile(here,'gr4jB','gr4jB.cpp')
        fullfile(here,'hbv','hbv.cpp')
        fullfile(here,'cfe_nwm','cfe_nwm.cpp')
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
        'GR4J-A, GR4J-B, HBV, CFE-NWM...\n']);

    mex( ...
        '-R2018a', ...
        '-O', ...
        sources{:}, ...
        '-outdir',here, ...
        '-output','crr_model_mex');

    mexFile = fullfile(here, ...
        ['crr_model_mex.' mexext]);

    if ~isfile(mexFile)
        error('compile_crr_model_mex:BuildFailed', ...
            'MEX output was not created:\n%s',mexFile);
    end

    fprintf('Created: %s\n',mexFile);

end
