function flag = force_compile(mdl,verbose)
%FORCE_COMPILE Compile standalone MEX files for validated native cores.
%
% SYNOPSIS:
%   flag = force_compile(mdl,verbose)
%
% INPUT:
%   mdl         OPTIONAL model structure
%    .mcode      numerical solution code; compilation requires mcode = 4
%    .names      model names to compile
%                default:
%                ["hymod","hmodel","sacsma","xinanjiang", ...
%                 "gr4jA","gr4jB","hbv","cfe_nwm"]
%
%   verbose     OPTIONAL, 0/1. Add -v to the MEX command when nonzero.
%
% OUTPUT:
%   flag         1  all requested models compiled successfully
%                0  skipped because mdl.mcode ~= 4
%               -1  one or more compilations failed
%               -2  one or more required source files were missing
%
% Each standalone model MEX is built from TWO source files:
%
%   crr_<model>_mex.cpp          MATLAB MEX gateway
%   <model>.cpp                  MATLAB-independent native model
%
% The build intentionally uses the same MEX options as the individually
% validated compile_crr_<model>_core_mex scripts:
%
%   mex('-R2018a','-O',...)
%
% Do NOT combine -R2018a with -largeArrayDims: recent MATLAB releases reject
% that combination.

    if nargin < 2 ...
            || isempty(verbose)
        verbose = 1;
    end

    if nargin < 1 ...
            || isempty(mdl)
        mdl = struct();
        mdl.mcode = 4;
        mdl.names = ["hymod","hmodel","sacsma","xinanjiang", ...
                     "gr4jA","gr4jB","hbv","cfe_nwm"];
    end

    if ~isfield(mdl,'mcode') ...
            || isempty(mdl.mcode)
        mdl.mcode = 4;
    end

    if ~isfield(mdl,'names') ...
            || isempty(mdl.names)
        mdl.names = ["hymod","hmodel","sacsma","xinanjiang", ...
                     "gr4jA","gr4jB","hbv","cfe_nwm"];
    end

    if mdl.mcode ~= 4
        warning('force_compile:mcodeNotCpp', ...
            ['mdl.mcode = %d; native-core ' ...
            'MEX compilation skipped.'], ...
            mdl.mcode);
        flag = 0;
        return
    end

    % -------------------------------------------------------------
    % Locate models directory
    % -------------------------------------------------------------
    this_dir = fileparts(mfilename('fullpath'));
    [~,this_name] = fileparts(this_dir);

    if strcmpi(this_name,'models')
        dir_models = this_dir;
    else
        dir_models = fullfile(this_dir,'models');
    end

    if ~isfolder(dir_models)
        warning('force_compile:missingModelsDir', ...
            'Could not find models directory: %s',dir_models);
        flag = -2;
        return
    end

    % -------------------------------------------------------------
    % Check compiler
    % -------------------------------------------------------------
    cc = mex.getCompilerConfigurations('C++','Selected');

    if isempty(cc)
        warning('force_compile:noCompiler', ...
            ['No C++ MEX compiler is configured. ' ...
             'Run: mex -setup C++']);
        flag = -1;
        return
    end

    fprintf('      force_compile: Compiler: %s\n',cc.Name);

    % requested alias | folder hint | core stem | gateway/MEX stem
    model_map = {
        'hymod',      'hymod',       'hymod',       'crr_hymod'
        'hmodel',     'hmodel',      'hmodel',      'crr_hmodel'
        'sacsma',     'sacsma',      'sacsma',      'crr_sacsma'
        'xinanjiang', 'xinanjiang',  'xinanjiang',  'crr_xinanjiang'
        'gr4ja',      'gr4jA',       'gr4jA',       'crr_gr4jA'
        'gr4jb',      'gr4jB',       'gr4jB',       'crr_gr4jB'
        'hbv',        'hbv',         'hbv',         'crr_hbv'
        'cfe_nwm',    'cfe_nwm',     'cfe_nwm',     'crr_cfe_nwm'
        };

    names = string(mdl.names(:));
    flag = 1;

    for im = 1:numel(names)

        req_name = lower(strtrim(names(im)));

        % Common aliases
        if req_name == "gr4j" || req_name == "gr4j-a"
            req_name = "gr4ja";
        elseif req_name == "gr4j-b"
            req_name = "gr4jb";
        end

        row = find(strcmpi(model_map(:,1),char(req_name)),1,'first');

        if isempty(row)
            warning('force_compile:unknownModel', ...
                'No validated native-core definition for "%s".',names(im));
            flag = local_merge_flag(flag,-2);
            continue
        end

        folder_hint = model_map{row,2};
        core_stem   = model_map{row,3};
        mex_name    = model_map{row,4};

        model_dir = local_find_model_dir(dir_models,folder_hint);

        if isempty(model_dir)
            warning('force_compile:missingModelDir', ...
                'Could not find model folder "%s" in %s.', ...
                folder_hint,dir_models);
            flag = local_merge_flag(flag,-2);
            continue
        end

        core_cpp = fullfile(model_dir,[core_stem,'.cpp']);
        gate_cpp = fullfile(model_dir,[mex_name,'_mex.cpp']);
        mex_out  = fullfile(model_dir,[mex_name,'.',mexext]);

        fprintf('\n');
        fprintf('      force_compile: Preparing %s\n', ...
            upper(char(req_name)));
        fprintf('      folder : %s\n',model_dir);
        fprintf('      core   : %s\n',core_cpp);
        fprintf('      gateway: %s\n',gate_cpp);

        missing = false;

        if ~isfile(core_cpp)
            warning('force_compile:missingCore', ...
                'Missing native core source: %s',core_cpp);
            missing = true;
        end

        if ~isfile(gate_cpp)
            warning('force_compile:missingGateway', ...
                'Missing standalone MEX gateway: %s',gate_cpp);
            missing = true;
        end

        if missing
            flag = local_merge_flag(flag,-2);
            continue
        end

        % ---------------------------------------------------------
        % Clear loaded MEX files and remove old output
        % ---------------------------------------------------------
        clear mex  %#ok
        rehash toolboxcache

        old_files = {
            mex_out
            fullfile(model_dir,[mex_name,'.lib'])
            fullfile(model_dir,[mex_name,'.exp'])
            fullfile(model_dir,[mex_name,'.obj'])
            fullfile(model_dir,[mex_name,'.o'])
            };

        for jf = 1:numel(old_files)
            try
                if isfile(old_files{jf})
                    delete(old_files{jf});
                end
            catch MEdelete
                warning('force_compile:deleteFailed', ...
                    'Could not delete %s:\n%s', ...
                    old_files{jf},MEdelete.message);
            end
        end

        % ---------------------------------------------------------
        % Compile using exactly the same options as the individually
        % validated core build scripts.
        % ---------------------------------------------------------
        try

            args = { ...
                '-R2018a', ...
                '-O', ...
                '-outdir',model_dir, ...
                '-output',mex_name, ...
                gate_cpp, ...
                core_cpp};

            if verbose
                args = [{'-v'},args];  %#ok
            end

            mex(args{:});

            if isfile(mex_out)
                fprintf(['      force_compile: Successfully ' ...
                         'compiled %s\n'],mex_out);
            else
                warning('force_compile:noMexOutput', ...
                    ['MEX command completed but ' ...
                    'output was not found: %s'], ...
                    mex_out);
                flag = local_merge_flag(flag,-1);
            end

        catch ME

            warning('force_compile:compileFailed', ...
                'Compilation failed for %s:\n%s', ...
                mex_name,ME.message);

            flag = local_merge_flag(flag,-1);
        end
    end
end


function model_dir = local_find_model_dir(dir_models,mname)
%LOCAL_FIND_MODEL_DIR Return model directory, case-insensitive.

    model_dir = '';

    candidate = fullfile(dir_models,mname);

    if isfolder(candidate)
        model_dir = candidate;
        return
    end

    dd = dir(dir_models);
    dd = dd([dd.isdir]);

    names = {dd.name};
    dd = dd(~ismember(names,{'.','..'}));

    if isempty(dd)
        return
    end

    folder_names = {dd.name};
    idx = find(strcmpi(folder_names,mname),1,'first');

    if ~isempty(idx)
        model_dir = fullfile(dir_models,folder_names{idx});
    end
end


function flag = local_merge_flag(flag,new_flag)
%LOCAL_MERGE_FLAG Preserve the most informative failure state.
%
% -2 indicates a missing source/folder and takes precedence over -1.

    if flag == -2
        return
    end

    if new_flag == -2
        flag = -2;
    elseif new_flag == -1 && flag == 1
        flag = -1;
    end
end
