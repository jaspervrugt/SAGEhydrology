function out=dataset_installation_manifest(action,varargin)
%DATASET_INSTALLATION_MANIFEST Installation provenance shared by SAGE/SITE.
% Unknown provenance is never inferred from current software defaults.
persistent pending
out=[];
switch action
    case 'matches'
        m=local_read(varargin{1});
        out=~isempty(m)&& ...
            strcmp(m.region,varargin{2})&& ...
            strcmp( m.installedVersion,varargin{3});
    case 'installed'
        region=varargin{1};root=varargin{2};version=varargin{3};m= ...
            local_read(root);
        if isempty(m),m=local_observe(region,root);end
        assert(strcmp(m.region,region),'Manifest region differs');
        m.installedVersion=version;m.versionEvidence= ...
            'completed regional installer';
        m.lastInstallationCompletedAt=local_time();m.recordedAt= ...
            local_time();
        if numel(varargin)>3,m.sourceRecord=varargin{4};end
        local_write(root,m);out=m;
    case 'begin'
        pending=struct('region',varargin{1},'root',varargin{2}, ...
            'resolution',varargin{3},'sources',{{}},'archives',struct( ...
            'file',{},'algorithm',{},'checksum',{},'verified',{}));
    case 'clear'
        pending=[];
    case 'source'
        if isempty(pending),return,end
        url=char(varargin{1});url=regexprep(url,'[?#].*$','');
        pending.sources=unique([pending.sources {url}],'stable');
    case 'checksum'
        if isempty(pending),return,end
        [~,name,ext]=fileparts(varargin{1});
        r=struct('file',[name ext],'algorithm','MD5','checksum',lower( ...
            char(varargin{2})),'verified',true);
        if isempty(pending.archives)|| ...
            ~any(strcmp({ pending.archives.file},r.file))
            pending.archives(end+1)=r;
        end
    case 'finish'
        if isempty(pending),return,end
        out=local_read(pending.root);
        if isempty(out),out=local_observe(pending.region,pending.root);end
        observed=local_observe(pending.region,pending.root);
        if ~strcmp(observed.installedVersion,'unknown')
            out.installedVersion=observed.installedVersion; ...
                out.versionEvidence=observed.versionEvidence;
            out.sourceRecord=observed.sourceRecord;
        end
        out.legacyMarkers=observed.legacyMarkers;out.markerContents= ...
            observed.markerContents;
        event=struct('resolution',pending.resolution,'recordedAt', ...
            local_time(),'downloadCompleted',~isempty(pending.sources), ...
            'sourceURLs',{pending.sources},'archives',pending.archives);
        if ~isfield(out,'events')|| ...
            isempty(out.events)
            out.events=event;
        else
            out.events(end+1)=event;
        end
        if event.downloadCompleted,out.lastInstallationCompletedAt= ...
            event.recordedAt;end
        out.recordedAt=local_time();local_write(pending.root,out); ...
            pending=[];
    case 'migrate'
        region=varargin{1};root=varargin{2};out=local_read(root);
        if isempty(out),out=local_observe(region,root);local_write(root, ...
            out);end
    case 'read'
        out=local_read(varargin{1});
    otherwise
        error('dataset_installation_manifest:Action', ...
            'Unknown manifest action');
end
end
function m=local_observe(region,root)
m=struct('schema',1,'region',region,'installedVersion','unknown', ...
    'versionEvidence','none','recordedAt',local_time(), ...
    'lastInstallationCompletedAt','','sourceRecord','','sourceURLs',{{ ...
    }},'archives',struct('file',{},'algorithm',{},'checksum',{}, ...
    'verified',{}),'events',[]);
markers=dir(fullfile(root,'*.installed'));
m.legacyMarkers={markers.name};
if strcmp(region,'CAMELS_JP')&& ...
    isfile(fullfile(root, 'MERV_Jp_v2_1.installed'))
    m.installedVersion='2.1';m.versionEvidence='MERV_Jp_v2_1.installed'; ...
        m.sourceRecord='https://zenodo.org/records/22037429';
elseif strcmp(region,'CAMELS_AU')&& ...
    isfile(fullfile(root, 'CAMELS_AUS_v2.03.installed'))
    m.installedVersion='2.03';m.versionEvidence= ...
        'CAMELS_AUS_v2.03.installed';
end
% Preserve original marker evidence; it is not a new checksum verification.
m.markerContents=cell(size(markers));
for i=1:numel(markers),try,m.markerContents{i}=fileread(fullfile(root, ...
    markers(i).name));catch,m.markerContents{i}='';end,end
end
function m=local_read(root)
file=fullfile(root,'dataset_installation.json');m=[];
if isfile(file),m=jsondecode(fileread(file));assert(m.schema==1, ...
    'Unsupported installation manifest schema');end
end
function local_write(root,m)
assert(isfolder(root),'Dataset directory does not exist');
file=fullfile(root,'dataset_installation.json');temporary=[tempname( ...
    root) '.json'];
fid=fopen(temporary,'w','n','UTF-8');assert(fid>=0, ...
    'Cannot write installation manifest');
cleanup=onCleanup(@()local_remove(temporary));
try,fwrite(fid,jsonencode(m,'PrettyPrint',true),'char');fclose(fid); ...
    catch ME,fclose(fid);rethrow(ME);end
movefile(temporary,file,'f');clear cleanup
end
function local_remove(file),if isfile(file),delete(file);end,end
function value=local_time,value=char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));end
