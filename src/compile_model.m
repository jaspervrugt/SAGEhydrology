function [flag,status] = compile_model(mdl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPILE_MODEL Locate or compile C++/MEX model core
%
% SYNOPSIS:
%   flag = compile_model(mdl)
%   [flag,status] = compile_model(mdl)
%
% INPUT:
%   mdl         structure with model selection/settings
%    .model      choice of model
%                 1 hymod
%                 2 hmodel
%                 3 sacsma
%                 4 xinanjiang
%                 5 gr4jA
%                 6 hbv
%                 7 gr4jB
%    .mcode      numerical implementation
%                 1 Runge Kutta MATLAB
%                 2 ode45 MATLAB
%                 3 Euler MATLAB
%                 4 Runge Kutta C++/MEX
%    .mode       evaluation mode
%      'seq'      sequential
%      'par'      parallel
%    .names      list/cell array of model names
%
% OUTPUT:
%   flag        execution flag
%                0 no compilation needed; existing MEX found
%                1 successfully compiled
%               -1 compilation failed
%               -2 required file/folder not found
%               -3 deployed mode: MEX missing and compilation not allowed
%
%   status      short status message suitable for GUI logging
%
% NOTES:
%   When called with one output:
%
%       flag = compile_model(mdl);
%
%   the status is printed directly in the MATLAB command window.
%
%   When called with two outputs:
%
%       [flag,status] = compile_model(mdl);
%
%   no status is printed. The caller can then log status itself.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025
%   Revised for robust deployed/non-deployed and GUI behavior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Store the output behavior of the parent function call.
    %
    % This must be evaluated here. Calling nargout inside a nested helper
    % would return the output count of that helper rather than that of
    % compile_model.
    print_status = (nargout < 2);

    model       = mdl.model;
    model_names = mdl.names;
    mcode       = mdl.mcode;

    % Resolve selected model name
    if iscell(model_names)
        mname = char(model_names{model});

    elseif isstring(model_names)
        mname = char(model_names(model));

    elseif ischar(model_names)
        if size(model_names,1) > 1
            mname = model_names(model,:);
        else
            mname = model_names;
        end

    else
        error('compile_model:invalidModelNames', ...
            'mdl.names must be a cell array, string array, or character array.');
    end

    mname = strtrim(mname);

    flag   = 0;
    status = '';

    % -------------------------------------------------------------
    % Compilation is relevant only for the C++/MEX implementation
    % -------------------------------------------------------------
    if mcode ~= 4
        flag   = 0;
        status = 'MATLAB implementation selected';

        printStatus();
        return
    end

    % ----------------------------------------------------------
    % Locate SAGE root robustly from this function's location
    %
    % Assumption: compile_model.m is somewhere inside SAGE tree
    % ----------------------------------------------------------
    this_file = mfilename('fullpath');
    this_dir  = fileparts(this_file);

    sage_root  = this_dir;
    dir_models = fullfile(sage_root,'models');

    while ~isfolder(dir_models)

        parent_dir = fileparts(sage_root);

        % Stop when the file-system root has been reached
        if strcmp(parent_dir,sage_root)
            flag   = -2;
            status = 'SAGE models folder not found';

            warning('compile_model:modelsFolderNotFound', ...
                ['      Warning: compile_model: Could not locate the ' ...
                 'SAGE "models" directory from:' newline ...
                 '    %s'], ...
                 this_file);

            printStatus();
            return
        end

        sage_root  = parent_dir;
        dir_models = fullfile(sage_root,'models');
    end

    % -----------------------------------------
    % Resolve model folder name and source stem
    % -----------------------------------------
    is_user_model = strcmpi(mname,'user_model');
    switch lower(mname)

        case 'hymod'
            folder_name = 'hymod';
            cpp_stem    = 'crr_hymod';

        case 'hmodel'
            folder_name = 'hmodel';
            cpp_stem    = 'crr_hmodel';

        case 'sacsma'
            folder_name = 'sacsma';
            cpp_stem    = 'crr_sacsma';

        case 'xinanjiang'
            folder_name = 'Xinanjiang';
            cpp_stem    = 'crr_xinanjiang';

        case 'gr4ja'
            folder_name = 'gr4jA';
            cpp_stem    = 'crr_gr4jA';

        case 'gr4jb'
            folder_name = 'gr4jB';
            cpp_stem    = 'crr_gr4jB';

        case 'hbv'
            folder_name = 'hbv';
            cpp_stem    = 'crr_hbv';

        otherwise
            folder_name = mname;
            cpp_stem    = ['crr_',mname];
    end

    if is_user_model
        dir_model = local_user_model_dir(mdl,sage_root);
    else
        dir_model = fullfile(dir_models,folder_name);
    end

    if ~isfolder(dir_model)
        flag   = -2;
        status = 'model folder not found';

        warning('compile_model:modelFolderNotFound', ...
            ['      Warning: compile_model: Could not find model folder:' ...
             newline ...
             '    %s'], ...
             dir_model);

        printStatus();
        return
    end

    split_models = {'hymod','hmodel','sacsma','xinanjiang', ...
        'gr4ja','gr4jb','hbv','cfe_nwm'};
    is_split_model = ismember(lower(mname),split_models);

    if is_user_model
        source_names = {'crr_user_model_mex.cpp', ...
            'user_model_prepare.cpp','user_model.cpp'};
    elseif is_split_model
        if strcmpi(mname,'gr4ja')
            native_stem = 'gr4jA';
        else
            native_stem = lower(mname);
        end
        source_names = {[cpp_stem,'_mex.cpp'],[native_stem,'.cpp']};
    else
        source_names = {[cpp_stem,'.cpp']};
    end
    cpp_file = fullfile(dir_model,source_names{1});
    mex_file = fullfile(dir_model,[cpp_stem,'.',mexext]);

    % ------------------------------------------------------------
    % 1. Best case: compiled binary already exists on hard disk
    % ------------------------------------------------------------
    if isfile(mex_file)
        flag   = 0;
        status = 'MEX found on hard disk';

        printStatus();
        return
    end

    % -------------------------------------------------------------------
    % 2. Deployed mode: do not compile new MEX binaries during runtime
    % -------------------------------------------------------------------
    if isdeployed
        flag   = -3;
        status = 'MEX missing in deployed mode';

        warning('compile_model:mexMissingDeployed', ...
            ['      Warning: compile_model: Required MEX file not found ' ...
             'in deployed mode:' newline ...
             '    %s' newline ...
             '      Runtime recompilation is disabled. Ship the ' ...
             'precompiled platform-specific MEX file with the app.'], ...
             mex_file);

        printStatus();
        return
    end

    % ------------------------------------------------------------
    % 3. MATLAB desktop mode: source is required for compilation
    % ------------------------------------------------------------
    if ~isfile(cpp_file)
        flag   = -2;
        status = 'C++ source file not found';

        warning('compile_model:cppFileNotFound', ...
            ['      Warning: compile_model: Could not find C++ file:' ...
             newline ...
             '    %s'], ...
             cpp_file);

        printStatus();
        return
    end

    % ------------------------------------------------------------
    % Confirm that a C++ MEX compiler has been configured
    % ------------------------------------------------------------
    cc = mex.getCompilerConfigurations('C++','Selected');

    if isempty(cc)
        flag   = -1;
        status = 'C++ compiler not configured';

        warning('compile_model:noCompiler', ...
            ['      Warning: compile_model: No C++ MEX compiler is ' ...
             'configured.' newline ...
             '      Run: mex -setup C++']);

        printStatus();
        return
    end

    % ------------------------------------------------------------
    % Unload previously loaded MEX functions and refresh MATLAB
    % ------------------------------------------------------------
    try
        clear mex %#ok<CLMEX>
        rehash toolbox
        rehash toolboxcache
    catch
        % Failure to refresh is not necessarily fatal
    end

    % ------------------------------------------------------------
    % Remove output files left by previous compilation attempts
    % ------------------------------------------------------------
    tmp_files = {mex_file};

    if ispc
        tmp_files = [tmp_files, ...
            {fullfile(dir_model,[cpp_stem,'.lib'])}, ...
            {fullfile(dir_model,[cpp_stem,'.exp'])}, ...
            {fullfile(dir_model,[cpp_stem,'.obj'])}];

    elseif ismac || isunix
        tmp_files = [tmp_files, ...
            {fullfile(dir_model,[cpp_stem,'.o'])}];
    end

    for k = 1:numel(tmp_files)
        try
            if isfile(tmp_files{k})
                delete(tmp_files{k});
            end

        catch MEdelete
            warning('compile_model:temporaryFileDeleteFailed', ...
                ['      Warning: compile_model: Could not delete ' ...
                 'temporary compilation file:' newline ...
                 '    %s' newline ...
                 '      %s'], ...
                 tmp_files{k},MEdelete.message);
        end
    end

    % ------------------------------------------------------------
    % Print traditional compilation message for direct MATLAB use
    % ------------------------------------------------------------
    if print_status
        fprintf('    Compiling %s\n',cpp_file);
    end

    % ------------------------------------------------------------
    % Compile from inside the model directory
    % ------------------------------------------------------------
    old_dir    = pwd;
    cleanupObj = onCleanup(@() cd(old_dir)); %#ok<NASGU>

    try
        cd(dir_model);
        rehash toolboxcache

        if is_user_model
            user_compile = fullfile(dir_model,'user_model_compile.m');
            if ~isfile(user_compile)
                error('compile_model:userModelCompilerMissing', ...
                    'Cannot find user-model compiler: %s',user_compile);
            end
            user_model_compile(1);
            rehash toolbox
            rehash toolboxcache
            if isfile(mex_file)
                flag = 1;
                status = 'MEX compiled successfully';
            else
                flag = -1;
                status = 'MEX compilation failed';
            end
            printStatus();
            return
        end

        if is_split_model
            mex_args = [ ...
                {'-v','-R2018a','-O','-output',cpp_stem},source_names];
        else
            mex_args = [{'-v','-O','-largeArrayDims'},source_names];
        end

        % --------------------------------------------------------
        % Add platform/compiler-specific optimization flags
        % --------------------------------------------------------
        if is_split_model
            % The split gateway/core build uses MATLAB's selected compiler
            % defaults; force_compile uses the same validated configuration.
        elseif ispc
            mex_args = [mex_args, ...
                {'COMPFLAGS=$COMPFLAGS /O2 /DNDEBUG /fp:fast /GL /Gw /Gy', ...
                 'CXXFLAGS=$CXXFLAGS /O2 /DNDEBUG /EHsc', ...
                 'LINKFLAGS=$LINKFLAGS /LTCG /OPT:REF /OPT:ICF'}];

        elseif ismac
            mex_args = [mex_args, ...
                {'CXXFLAGS=$CXXFLAGS -O3 -DNDEBUG -std=c++11', ...
                 'LDFLAGS=$LDFLAGS'}];

        elseif isunix
            mex_args = [mex_args, ...
                {'CXXFLAGS=$CXXFLAGS -O3 -DNDEBUG -std=c++11', ...
                 'LDFLAGS=$LDFLAGS'}];
        end

        % --------------------------------------------------------
        % First attempt: custom optimized compilation settings
        % --------------------------------------------------------
        try
            mex(mex_args{:});

        catch MEcustom
            warning('compile_model:customFlagsFailed', ...
                ['      Warning: compile_model: Custom MEX flags failed; ' ...
                 'retrying with MATLAB default optimization settings.' ...
                 newline ...
                 '%s'], ...
                 MEcustom.message);

            % Remove any incomplete output before retrying
            try
                if isfile(mex_file)
                    delete(mex_file);
                end
            catch
            end

            % Second attempt: MATLAB default optimization settings
            if is_split_model
                mex('-v','-R2018a','-O','-output',cpp_stem, ...
                    source_names{:});
            else
                mex('-v','-O','-largeArrayDims',source_names{:});
            end
        end

        % Refresh MATLAB's view of the newly compiled binary
        rehash toolbox
        rehash toolboxcache

        % --------------------------------------------------------
        % Confirm compilation produced the expected MEX file
        % --------------------------------------------------------
        if isfile(mex_file)
            flag   = 1;
            status = 'MEX compiled successfully';

        else
            flag   = -1;
            status = 'MEX compilation failed';

            warning('compile_model:mexNotProduced', ...
                ['      Warning: compile_model: Compilation did not ' ...
                 'produce the expected MEX file:' newline ...
                 '    %s'], ...
                 mex_file);
        end

    catch ME
        flag   = -1;
        status = 'MEX compilation failed';

        warning('compile_model:compilationFailed', ...
            ['      Warning: compile_model: Compilation failed for %s:' ...
             newline ...
             '%s'], ...
             cpp_stem,ME.message);
    end

    printStatus();

    % =====================================================================
    % Nested utility: print only when status was not requested by caller
    % =====================================================================
    function printStatus()

        if ~print_status || isempty(status)
            return
        end

        switch flag
            case 0
                if mcode == 4
                    fprintf(['    Using compiled MEX file ' ...
                        'on hard disk:' newline ...
                        '      %s' newline], ...
                        mex_file);
                else
                    fprintf('    %s\n',status);
                end

            case 1
                fprintf(['    compile_model: ' ...
                    'Successfully compiled:' newline ...
                    '      %s' newline], ...
                    mex_file);

            otherwise
                % Detailed warning has already been issued above.
                fprintf('    compile_model: %s\n',status);
        end

    end

    function folder = local_user_model_dir(modelInfo,sageRoot)

        candidates = strings(0,1);
        located = which(['crr_user_model.' mexext]);
        if ~isempty(located)
            candidates(end+1) = string(fileparts(located));
        end
        if isstruct(modelInfo) ...
                && isfield(modelInfo,'root') ...
                && ~isempty(modelInfo.root)
            candidates(end+1) = string(fullfile( ...
                char(string(modelInfo.root)),'user_model'));
        end
        candidates(end+1) = string(fullfile( ...
            fileparts(sageRoot),'user_model'));
        candidates(end+1) = string(fullfile(sageRoot,'user_model'));
        candidates(end+1) = string(fullfile(pwd,'user_model'));

        folder = '';
        for ii = 1:numel(candidates)
            candidate = char(candidates(ii));
            if isfolder(candidate) ...
                    && (isfile(fullfile(candidate, ...
                    ['crr_user_model.' mexext])) ...
                    || isfile(fullfile(candidate, ...
                    'user_model_compile.m')))
                folder = candidate;
                return
            end
        end

        % Return the conventional external location so the caller's
        % existing warning includes the most useful expected path.
        folder = fullfile(fileparts(sageRoot),'user_model');

    end

end
