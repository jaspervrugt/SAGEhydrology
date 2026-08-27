function write_SAGE_live_script(templateFile,outputFile,codeBlocks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WRITE_SAGE_LIVE_SCRIPT Populate code cells in a SAGE Live Script template.
%
% The formatted text, equations, hyperlinks, image, and styles in the MLX
% template are retained. Only the gray MATLAB code paragraphs are replaced.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isfile(templateFile)
        error('write_SAGE_live_script:MissingTemplate', ...
            'Live Script template not found: %s',templateFile);
    end
    if ~iscell(codeBlocks)
        error('write_SAGE_live_script:InvalidBlocks', ...
            'codeBlocks must be a cell array of character vectors.');
    end

    workDir = tempname;
    mkdir(workDir);
    cleanupObj = onCleanup(@() local_cleanup(workDir));
    unzip(templateFile,workDir);

    documentFile = fullfile(workDir,'matlab','document.xml');
    document = fileread(documentFile);
    % [\s\S] matches every character, including line breaks. Each code
    % cell in the current template contains multiline MATLAB code.
    pattern = ['(<w:p><w:pPr><w:pStyle w:val="code"/>' ...
        '</w:pPr><w:r><w:t><!\[CDATA\[)([\s\S]*?)(\]\]>' ...
        '</w:t></w:r></w:p>)'];
    [starts,ends,tokens] = regexp(document,pattern, ...
        'start','end','tokens');
    if numel(starts) ~= numel(codeBlocks)
        error('write_SAGE_live_script:CellCount', ...
            ['Template contains %d code cells, but %d replacement ' ...
             'blocks were supplied.'],numel(starts),numel(codeBlocks));
    end

    for k = numel(starts):-1:1
        code = char(string(codeBlocks{k}));
        if contains(code,']]>')
            error('write_SAGE_live_script:InvalidCode', ...
                'Code block %d contains the XML CDATA terminator.',k);
        end
        replacement = [tokens{k}{1} code tokens{k}{3}];
        document = [document(1:starts(k)-1) replacement ...
            document(ends(k)+1:end)];
    end

    fid = fopen(documentFile,'w','n','UTF-8');
    if fid < 0
        error('write_SAGE_live_script:DocumentWrite', ...
            'Could not write temporary Live Script document.');
    end
    fileCleanup = onCleanup(@() fclose(fid));
    fwrite(fid,unicode2native(document,'UTF-8'),'uint8');
    clear fileCleanup

    outputFile = char(string(outputFile));
    [~,~,extension] = fileparts(outputFile);
    if ~strcmpi(extension,'.mlx')
        outputFile = [outputFile '.mlx'];
    end
    [outputDir,~,~] = fileparts(outputFile);
    if isempty(outputDir)
        outputDir = pwd;
        outputFile = fullfile(outputDir,outputFile);
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    % Assemble the MLX completely before replacing an existing export.
    % Writing the staging file beside the destination also avoids a
    % cross-volume move from the MATLAB Runtime cache on deployed systems.
    temporaryZip = fullfile(workDir,'populated_template.zip');
    topLevel = {'[Content_Types].xml','_rels','matlab','mathml', ...
        'media','metadata'};
    zip(temporaryZip,topLevel,workDir);

    stagedFile = [tempname(outputDir) '.mlx'];
    stagedCleanup = onCleanup(@() local_delete_file(stagedFile));
    [copied,message] = copyfile(temporaryZip,stagedFile,'f');
    if ~copied
        error('write_SAGE_live_script:StageWrite', ...
            'Could not stage Live Script in %s: %s',outputDir,message);
    end

    local_validate_mlx(stagedFile,codeBlocks);
    [moved,message] = movefile(stagedFile,outputFile,'f');
    if ~moved
        error('write_SAGE_live_script:OutputWrite', ...
            'Could not write Live Script %s: %s',outputFile,message);
    end
    clear stagedCleanup

    if ~isfile(outputFile)
        error('write_SAGE_live_script:MissingOutput', ...
            'Live Script was not created: %s',outputFile);
    end
end

function local_validate_mlx(mlxFile,codeBlocks)
% Confirm that the staged archive is readable and contains populated code.

    checkDir = tempname;
    mkdir(checkDir);
    cleanupObj = onCleanup(@() local_cleanup(checkDir));
    unzip(mlxFile,checkDir);
    documentFile = fullfile(checkDir,'matlab','document.xml');
    if ~isfile(documentFile)
        error('write_SAGE_live_script:InvalidOutput', ...
            'Generated Live Script does not contain matlab/document.xml.');
    end
    document = fileread(documentFile);
    for k = 1:numel(codeBlocks)
        code = char(string(codeBlocks{k}));
        if ~isempty(code) && ~contains(document,code)
            error('write_SAGE_live_script:InvalidOutput', ...
                'Generated Live Script is missing populated code block %d.',k);
        end
    end
end

function local_delete_file(fileName)
    try
        if isfile(fileName)
            delete(fileName);
        end
    catch
    end
end

function local_cleanup(workDir)
    try
        if isfolder(workDir)
            rmdir(workDir,'s');
        end
    catch
    end
end
