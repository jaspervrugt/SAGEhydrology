function folder=sage_runtime_asset_folder(name,root,bundleRoot)
%SAGE_RUNTIME_ASSET_FOLDER Find external or packaged public runtime assets.
    if nargin<3
        bundleRoot='';
        if isdeployed, bundleRoot=ctfroot; end
    end
    candidates={fullfile(root,name),fullfile(root,'SAGEhydrology',name), ...
        fullfile(pwd,name)};
    if isdeployed && ispc
        % The internal getExecutablePath function is unavailable in Runtime.
        % Windows process metadata is optional; discovery must never fail
        % merely because executable-location inspection is unavailable.
        try
            process=System.Diagnostics.Process.GetCurrentProcess();
            executable=fileparts(char(process.MainModule.FileName));
            candidates=[{fullfile(executable,name)},candidates];
        catch
            % Configured-root, current-folder, and CTF fallbacks remain valid.
        end
    end
    if ~isempty(bundleRoot)
        candidates{end+1}=fullfile(bundleRoot,name);
        candidates{end+1}=fullfile(bundleRoot,'SAGEhydrology',name);
        marker='user_model_info.mat';
        if strcmp(name,'projects'), marker='project.json'; end
        matches=dir(fullfile(bundleRoot,'**',marker));
        for k=1:numel(matches)
            candidate=matches(k).folder;
            if strcmp(name,'projects'), candidate=fileparts(candidate); end
            candidates{end+1}=candidate; %#ok<AGROW>
        end
    end
    folder=fullfile(root,name);
    for k=1:numel(candidates)
        candidate=candidates{k};
        if strcmp(name,'projects')
            available=~isempty(dir(fullfile(candidate,'*','project.json')));
        else
            available=isfile(fullfile(candidate,'user_model_info.mat'));
        end
        if available, folder=candidate; return, end
    end
end
