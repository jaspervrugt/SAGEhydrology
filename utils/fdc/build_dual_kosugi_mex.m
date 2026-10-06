function mexFile = build_dual_kosugi_mex()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD_DUAL_KOSUGI_MEX Compile the dual-Kosugi objective for this platform.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sourceDir = fileparts(mfilename('fullpath'));
    sourceFile = fullfile(sourceDir,'dual_kosugi_objective_mex.cpp');
    if ~isfile(sourceFile)
        error('build_dual_kosugi_mex:MissingSource', ...
            'Cannot find %s.',sourceFile);
    end
    oldDir = pwd;
    cleanup = onCleanup(@()cd(oldDir));
    cd(sourceDir);
    mex('-R2018a',sourceFile);
    mexFile = fullfile(sourceDir, ...
        ['dual_kosugi_objective_mex.' mexext]);
    clear cleanup
end
