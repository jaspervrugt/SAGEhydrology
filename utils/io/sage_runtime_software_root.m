function root = sage_runtime_software_root(fallback,executableFolder)
% Prefer the user's installation folder over the extracted Runtime cache.
    if nargin<2
        executableFolder='';
        if ispc
            try
                process=System.Diagnostics.Process.GetCurrentProcess();
                executableFolder=fileparts(char(process.MainModule.FileName));
            catch
            end
        end
    end
    candidates={executableFolder,pwd};
    root=fallback;
    for k=1:numel(candidates)
        folder=candidates{k};
        if isempty(folder), continue, end
        if isfolder(fullfile(folder,'Data')) ...
                || isfolder(fullfile(folder,'SAGEhydrology'))
            root=folder;
            return
        end
    end
end
