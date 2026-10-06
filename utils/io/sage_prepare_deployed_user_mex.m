function sage_prepare_deployed_user_mex(mdl)
% Runtime binds compiled MEX calls to the extracted packaged binary.
% Install the user's compatible binary in that slot before it is loaded.
    if ~isdeployed, return, end
    root=pwd;
    if isfield(mdl,'root') && ~isempty(mdl.root), root=char(mdl.root); end
    folder=sage_runtime_asset_folder('user_model',root);
    external=fullfile(folder,['crr_user_model.' mexext]);
    assert(isfile(external),'SAGE:UserModelMexMissing', ...
        'The precompiled user-model MEX is missing: %s',external);
    packaged=which('crr_user_model');
    assert(~isempty(packaged) && isfile(packaged), ...
        'SAGE:UserModelMexNotPackaged', ...
        'This build has no packaged user-model MEX interface.');
    if strcmpi(external,packaged), return, end
    persistent checkedExternal checkedPackaged checkedStamp
    info=dir(external);
    stamp=[info.bytes info.datenum];
    if isequal(external,checkedExternal) && isequal(packaged,checkedPackaged) ...
            && isequal(stamp,checkedStamp)
        return
    end
    if ~isequal(local_bytes(external),local_bytes(packaged))
        try
            clear crr_user_model
            copyfile(external,packaged,'f');
        catch ME
            error('SAGE:UserModelMexInstall', ...
                ['Cannot load the replacement user-model MEX. Close SAGE ' ...
                 'and its parallel workers, then reopen it. %s'],ME.message);
        end
    end
    checkedExternal=external;
    checkedPackaged=packaged;
    checkedStamp=stamp;
end

function bytes=local_bytes(file)
    fid=fopen(file,'rb');
    assert(fid>=0,'SAGE:UserModelMexRead','Cannot read MEX: %s',file);
    cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    bytes=fread(fid,Inf,'*uint8');
end
