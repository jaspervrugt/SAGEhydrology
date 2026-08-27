function mexFile = compile_information(force)
%COMPILE_INFORMATION Compile the native SAGE Information diagnostic kernel.
%
%   mexFile = compile_information()
%   mexFile = compile_information(force)
%
% Source MATLAB only. Deployed SAGE never compiles C++ code.
%
% This helper tries to compile information_bottleneck_SAGE_mex with OpenMP
% first. If OpenMP flags are unsupported on the current platform/compiler,
% it falls back automatically to a serial build.
%
% Usage notes
% -----------
% Serial (safest with MATLAB parfeval/parpool workers):
%   opts.info_threads = 1;
%
% Parallel native computation:
%   opts.info_threads = 4;   % explicit number of threads
%   opts.info_threads = 0;   % automatic maximum available OpenMP threads
%
% Master switch:
%   opts.info_parallel = false;   % forces serial mode

    if nargin < 1 ...
            || isempty(force)
        force = false;
    end
    force = logical(force);

    if isdeployed
        error('SAGE:Information:DeployedCompilation', ...
            ['Information MEX compilation ' ...
            'is unavailable in deployed mode. ' ...
             'Compile it on the target ' ...
             'platform before packaging SAGE.']);
    end

    here = fileparts(mfilename('fullpath'));
    sourceFile = fullfile(here, ...
        'information_bottleneck_SAGE_mex.cpp');
    mexFile = fullfile(here, ...
        ['information_bottleneck_SAGE_mex.' mexext]);

    assert(isfile(sourceFile), ...
        'SAGE:Information:MissingSource', ...
        'Information C++ source not found: %s',sourceFile);

    if ~force && isfile(mexFile)
        s = dir(sourceFile);
        m = dir(mexFile);
        if ~isempty(s) ...
                && ~isempty(m) ...
                && m.datenum >= s.datenum
            fprintf(['Information MEX is ' ...
                'already up to date:\n  %s\n'],mexFile);
            return
        end
    end

    cc = mex.getCompilerConfigurations('C++','Selected');
    if isempty(cc)
        error('SAGE:Information:NoCompiler', ...
            ['No C++ MEX compiler is configured. Run:\n\n' ...
             '    mex -setup C++\n\nand retry.']);
    end

    clear information_bottleneck_SAGE_mex
    fprintf('\nCompiling SAGE Information MEX\n');
    fprintf('  Platform : %s\n',computer);
    fprintf('  MEX ext  : .%s\n',mexext);
    fprintf('  Source   : %s\n',sourceFile);
    fprintf('  Output   : %s\n',mexFile);

    compilerName = '';
    compilerMfr  = '';
    try
        compilerName = string(cc.Name);
        compilerMfr  = string(cc.Manufacturer);
    catch
    end

    baseArgs = {'-R2018a','-O','-outdir',here,sourceFile};

    % openmpArgs = baseArgs;
    % openmpNote = "serial";
    % tryOpenMP = true;

    if contains(lower(compilerMfr),'gnu') ...
            || contains(lower(compilerName),'mingw') ...
            || contains(lower(compilerName),'g++') ...
            || contains(lower(compilerName),'gcc')
        openmpArgs = [baseArgs, ...
            {'CXXFLAGS=$CXXFLAGS -fopenmp', ...
             'LDFLAGS=$LDFLAGS -fopenmp'}];
        openmpNote = "GNU/MinGW OpenMP";
    elseif contains(lower(compilerMfr),'microsoft') ...
            || contains(lower(compilerName),'visual')
        openmpArgs = [baseArgs, ...
            {'COMPFLAGS=$COMPFLAGS /openmp'}];
        openmpNote = "MSVC OpenMP";
    else
        % Unknown toolchain: try GCC-style OpenMP and fall back if it fails.
        openmpArgs = [baseArgs, ...
            {'CXXFLAGS=$CXXFLAGS -fopenmp', ...
             'LDFLAGS=$LDFLAGS -fopenmp'}];
        openmpNote = "generic OpenMP attempt";
    end

    %compiledWithOpenMP = false;

    try
        fprintf('  Trying %s build ...\n',openmpNote);
        mex(openmpArgs{:});
        compiledWithOpenMP = true;
    catch ME1
        warning('SAGE:Information:OpenMPCompileFailed', ...
            ['OpenMP build failed; ' ...
            'falling back to serial build.\n' ...
             'Reason: %s'],ME1.message);

        try
            mex(baseArgs{:});
            compiledWithOpenMP = false;
        catch ME2
            error('SAGE:Information:CompileFailed', ...
                ['Could not compile ' ...
                'information_bottleneck_SAGE_mex:\n%s'], ...
                ME2.message);
        end
    end

    assert(isfile(mexFile), ...
        'SAGE:Information:MissingOutput', ...
        'Compilation did not create: %s',mexFile);

    rehash;
    if compiledWithOpenMP
        fprintf(['Information MEX ready (OpenMP ' ...
            'enabled when requested):\n  %s\n\n'],mexFile);
    else
        fprintf(['Information MEX ready ' ...
            '(serial build):\n  %s\n\n'],mexFile);
    end
end
