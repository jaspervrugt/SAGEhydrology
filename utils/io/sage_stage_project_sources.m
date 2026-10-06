function count = sage_stage_project_sources(projectRoot,stageRoot)
% Package readable public project exports as data, not compiled MATLAB code.
% Only staged descriptors change; development project files remain intact.
    projects=sage_project_catalog('discover',projectRoot);
    count=0;
    for k=1:numel(projects)
        p=projects(k);
        source=p.source_export;
        if isempty(source), continue, end
        assert(isfile(source),'build_SAGE:MissingProjectSource', ...
            'Cannot package project source: %s',source);
        [~,folderName]=fileparts(p.folder);
        folder=fullfile(stageRoot,folderName);
        [~,name,extension]=fileparts(source);
        asset=[name extension '.txt'];
        copyfile(source,fullfile(folder,asset),'f');
        descriptor=fullfile(folder,'project.json');
        meta=jsondecode(fileread(descriptor));
        meta.source_export=asset;
        fid=fopen(descriptor,'w','n','UTF-8');
        assert(fid>=0,'build_SAGE:ProjectDescriptor', ...
            'Cannot write staged project descriptor: %s',descriptor);
        cleanup=onCleanup(@()fclose(fid));
        fprintf(fid,'%s\n',jsonencode(meta,PrettyPrint=true));
        clear cleanup
        count=count+1;
    end
    fprintf('Including %d readable project MATLAB source exports.\n',count);
end
