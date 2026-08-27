function out = install_SAGEhydrology(rootDir,what,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INSTALL_SAGEHYDROLOGY  Installer for SAGEhydrology + CAMELS data layout
%
% SYNOPSIS:
%   out = install_SAGEhydrology(rootDir,what,opts)
%
% INPUT:
%   rootDir   root folder that will contain:
%              rootDir/SAGEhydrology/   (code; assumed already downloaded)
%              rootDir/Data/            (regional data folders)
%   what      data product to install:
%              'daily'   daily data for the active/default region
%              'hourly'  hourly data for CAMELS-US only
%              'all'     all supported data for the active/default region
%   opts      OPTIONAL structure with fields:
%    .region        region code:
%                     'CAMELS_AU'  Australia
%                     'CAMELS_BR'  Brazil
%                     'CAMELS_CL'  Chile
%                     'CAMELS_GB'  Great Britain
%                     'CAMELS_US'  United States
%                   default: 'CAMELS_US'
%    .force         0/1 re-download + re-extract even if target exists
%    .tmpDir        temp folder for downloads
%                   default: user Downloads folder
%    .keepArchives  0/1 keep downloaded archives
%                   default: 0
%    .verbose       0/1 print progress
%                   default: 1
%    .logFcn        function handle: logFcn(str)
%                   default: @disp
%    .timeoutSec    web timeout in seconds
%                   default: 60
%
% OUTPUT:
%   out       structure with fields:
%              .rootDir
%              .what
%              .ok
%              .messages
%              .paths
%
% NOTES:
%   - This function does not redistribute CAMELS data. It only automates
%     downloading from the official URLs provided in camels_urls().
%
%   - Supported products:
%       CAMELS_AU: daily
%       CAMELS_BR: daily
%       CAMELS_CL: daily
%       CAMELS_GB: daily
%       CAMELS_US: daily and hourly
%
%   - CAMELS-US daily:
%       Zenodo record: https://zenodo.org/records/15529996
%       ZIP archive:
%          basin_timeseries_v1p2_metForcing_obsFlow.zip
%
%   - CAMELS-US hourly:
%       Zenodo record: https://zenodo.org/records/4072701
%       TAR.GZ archives:
%          nldas_hourly_csv.tar.gz
%          usgs_streamflow_csv.tar.gz
%
%   - Folder conventions created:
%
%       rootDir/Data/CAMELS_US/daily/v1p2/forcing/{daymet,maurer,nldas}/...
%       rootDir/Data/CAMELS_US/daily/v1p2/streamflow/...
%       rootDir/Data/CAMELS_US/daily/v1p2/metadata/...
%       rootDir/Data/CAMELS_US/hourly/forcing/...
%       rootDir/Data/CAMELS_US/hourly/streamflow/...
%
%       rootDir/Data/CAMELS_GB/daily/...
%
%       rootDir/Data/CAMELS_BR/daily/...
%
%   - Shared CAMELS-US metadata files, when installed, are placed in:
%
%       rootDir/Data/CAMELS_US/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Mar. 2026                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defaults + bookkeeping
if nargin < 2 ...
        || isempty(what)
    what = 'all';
end
if nargin < 3, opts = struct(); end
opts = local_set_default_opts(opts);

out = struct();
out.rootDir = char(rootDir);
out.what = char(lower(string(what)));
out.ok = false;
out.messages = {};
out.paths = struct();

logf = opts.logFcn;

try
    % validate + make dirs
    if ~isfolder(out.rootDir)
        error('Root folder does not exist: %s', ...
            out.rootDir);
    end
    if ~isfield(opts,'region') ...
            || isempty(opts.region)
        opts.region = 'CAMELS_US';
    end
    dirDataRoot = fullfile(out.rootDir,'Data');
    dirData = fullfile(dirDataRoot,opts.region);
    
    if ~exist(dirDataRoot,'dir')
        mkdir(dirDataRoot);
    end
    if ~exist(dirData,'dir')
        mkdir(dirData);
    end
    % dirData = fullfile(out.rootDir, ...
    %     'Data');
    % if ~exist(dirData,'dir')
    %     mkdir(dirData);
    % end

    dirTmp = opts.tmpDir;
    % Downloads to user Downloads directory
    if isempty(dirTmp)
        dirTmp = local_default_download_dir();
    end
    if ~exist(dirTmp,'dir')
        mkdir(dirTmp);
    end

    out.paths.dirData = dirData;
    out.paths.dirTmp = dirTmp;

    U = camels_urls(); % central URLs
   
    % route by "what"
    switch out.what
        case 'daily'
            local_install_daily(dirData,dirTmp,U,opts,logf);
            local_copy_gauge_information_to_data(dirData,logf);
    
        case 'hourly'
            local_install_hourly(dirData,dirTmp,U,opts,logf);
    
        case 'all'
            local_install_daily(dirData,dirTmp,U,opts,logf);
            local_install_hourly(dirData,dirTmp,U,opts,logf);
            local_copy_gauge_information_to_data(dirData,logf);
        otherwise
            error(['Unknown "what" = "%s". ' ...
                'Use daily|hourly|all.'],out.what);
    end

    % install shared CAMELS metadata text files into Data/
    local_install_shared_metadata(dirData,dirTmp,U,opts,logf);  

    out.ok = true;
    out.messages = [out.messages; {'Install finished.'}];
    if opts.verbose
        logf('install_SAGEhydrology: done.');
    end

catch ME
    out.ok = false;
    out.messages = [out.messages; {ME.message}];
    if opts.verbose
        logf('install_SAGEhydrology: ERROR');
        logf(ME.message);
        for ii = 1:numel(ME.stack)
            logf(sprintf('  at %s (line %d)', ...
                ME.stack(ii).name, ...
                ME.stack(ii).line));
        end
    end
end

end

% ================
% Local installers
% ================
function local_install_daily(dirData,dirTmp,U,opts,logf)

logf('--- Installing DAILY CAMELS (v1p2) ---');

dirDaily = fullfile(dirData,'daily');
if ~exist(dirDaily,'dir')
    mkdir(dirDaily); 
end

zipName = U.daily.zip_name;
zipUrl = U.daily.zip_url;

% where we expect final layout
dirV1p2 = fullfile(dirDaily,'v1p2');

% Robust "already installed?" check
dailyOK = local_tree_has_files(fullfile( ...
    dirV1p2,'forcing'),'*.txt',1) && ...
    local_tree_has_files(fullfile(dirV1p2, ...
    'streamflow'),'*.txt',1);

if dailyOK && ~opts.force
    logf(['Daily target exists ' ...
        'and looks complete: ' dirV1p2]);
    logf('Skipping daily install.');
    return
elseif exist(dirV1p2,'dir')
    logf(['Daily install appears incomplete. ' ...
        'Cleaning and reinstalling ...']);
    try
        rmdir(dirV1p2,'s');
    catch
    end
end

% download
zipPath = fullfile(dirTmp,zipName);

% if exist(dirV1p2,'dir') && ~opts.force
%     logf(['Daily target exists: ' dirV1p2]);
%     logf(['Skipping daily install ' ...
%         '(set opts.force=1 to re-install).']);
%     return
% end

% download
%zipPath = fullfile(dirTmp,zipName);
local_download_if_needed(zipUrl,zipPath,opts,logf, ...
    U.daily.size_bytes,'daily',U.daily.md5,1e6);

logf('Daily download complete.');

% validate archive before extraction
local_assert_file_nontrivial(zipPath,1e6,'daily ZIP');   % at least 1 MB
local_assert_valid_zip(zipPath);
logf('Daily archive is valid.');

% extract (into dirDaily)
logf('Unzipping daily archive. This can take several minutes ...');
drawnow;
local_safe_unzip(zipPath,dirDaily,logf);
logf('Finished unzipping daily archive.');
drawnow;

logf('Daily extraction complete.');

% detect extracted top folder(s)
% expected from CAMELS daily zip:
%   basin_dataset_public_v1p2
%   basin_dataset_public (maybe)
dirA = fullfile(dirDaily, ...
    'basin_dataset_public_v1p2');
dirB = fullfile(dirDaily, ...
    'basin_dataset_public');

if exist(dirB,'dir')
    logf('Deleting folder basin_dataset_public ...');
    try
        rmdir(dirB,'s');
    catch
        % if delete fails, keep going
    end
end

if ~exist(dirA,'dir')
    % sometimes it may unzip into a nested folder; search one level
    dirA = local_find_subdir(dirDaily, ...
        'basin_dataset_public_v1p2');
end
if isempty(dirA) || ~exist(dirA,'dir')
    error(['Daily extraction did not ' ...
        'produce basin_dataset_public_v1p2. ' ...
        'Check archive contents in %s'], ...
        dirDaily);
end

% rename basin_dataset_public_v1p2 -> v1p2
if exist(dirV1p2,'dir')
    rmdir(dirV1p2,'s');
end
logf('Renaming basin_dataset_public_v1p2 -> v1p2');
movefile(dirA,dirV1p2);

% inside v1p2: rename folders to match SAGEhydrology conventions
ren = { ...
    {'basin_mean_forcing',   'forcing'}, ...
    {'usgs_streamflow',      'streamflow'}, ...
    {'basin_metadata',       'metadata'}, ...
    {'hru_forcing',          'forcing_hru'}, ...
    {'elev_bands_forcing',   'forcing_elev_bands'} ...
    };

% some zips have a space typo in the folder name; handle defensively
ren_alias = { ...
    {'elev_ bands_forcing',  'forcing_elev_bands'} ...
    };

for i = 1:numel(ren)
    src = fullfile(dirV1p2,ren{i}{1});
    dst = fullfile(dirV1p2,ren{i}{2});
    if exist(src,'dir') ...
            && ~exist(dst,'dir')
        logf(sprintf('Renaming %s -> %s', ...
            ren{i}{1},ren{i}{2}));
        movefile(src,dst);
    end
end

for i = 1:numel(ren_alias)
    src = fullfile(dirV1p2,ren_alias{i}{1});
    dst = fullfile(dirV1p2,ren_alias{i}{2});
    if exist(src,'dir') ...
            && ~exist(dst,'dir')
        logf(sprintf('Renaming %s -> %s', ...
            ren_alias{i}{1},ren_alias{i}{2}));
        movefile(src,dst);
    end
end

% quick sanity
mustExist = { ...
    fullfile(dirV1p2,'forcing'), ...
    fullfile(dirV1p2,'streamflow') ...
    };
for i = 1:numel(mustExist)
    if ~exist(mustExist{i},'dir')
        logf(['WARNING: missing expected ' ...
            'folder: ' mustExist{i}]);
    end
end

% cleanup archive
if ~opts.keepArchives
    local_try_delete(zipPath);
end

logf('Daily install complete.');
end

% ----------------------
% local helper functions
% ----------------------
function local_install_hourly(dirData,dirTmp,U,opts,logf)
%LOCAL_INSTALL_HOURLY

logf('--- Installing HOURLY CAMELS (NLDAS + USGS) ---');

dirHourly = fullfile(dirData,'hourly');
if ~exist(dirHourly,'dir') 
    mkdir(dirHourly); 
end

dirForcing = fullfile(dirHourly, ...
    'forcing');
dirFlow = fullfile(dirHourly, ...
    'streamflow');

% ---------------------------------------------------------
% Robust "already installed?" check (must contain real files)
% ---------------------------------------------------------
% forcing:    *_hourly_nldas.csv
% streamflow: *-usgs-hourly.csv
forcingOK = local_dir_has_files(dirForcing, ...
    '*_hourly_nldas.csv',1);
flowOK = local_dir_has_files(dirFlow,   ...
    '*-usgs-hourly.csv',1);
hourlyOK = forcingOK && flowOK;

% If fully installed and not forcing, skip
if hourlyOK && ~opts.force
    logf(['Hourly targets exist + look complete: ' ...
        dirForcing ' and ' dirFlow]);
    logf('Skipping hourly install (set opts.force=1 to re-install).');
    return
end

% ---------------------------------------------------------
% Half-broken install detected -> clean + FORCE RE-DOWNLOAD
% ---------------------------------------------------------
% If folders exist but are incomplete, we:
%  1) remove them (so we don't mix old/new)
%  2) force download even if archives exist (avoid unzip of bad tgz)
optsDL = opts;
if ~hourlyOK && (exist(dirForcing,'dir') ...
        || exist(dirFlow,'dir'))
    logf(['Hourly install appears incomplete. ' ...
        'Cleaning and re-installing ...']);
    logf(['Forcing fresh downloads ' ...
        '(ignoring any existing archives).']);
    optsDL.force = 1;  % <-- key: re-download first, then extract

    try
        if exist(dirForcing,'dir')
            rmdir(dirForcing,'s');
        end
    catch
    end
    try
        if exist(dirFlow,'dir')
            rmdir(dirFlow,'s');
        end
    catch
    end
end

% Download archives
tgzF_name = U.hourly.forcing_tgz_name;
tgzF_url = U.hourly.forcing_tgz_url;
tgzF_path = fullfile(dirTmp,tgzF_name);
local_download_if_needed(tgzF_url, ...
    tgzF_path,optsDL,logf, ...
    U.hourly.forcing_size_bytes, ...
    'hourly forcing','',1e6);

local_assert_file_nontrivial(tgzF_path, ...
    1e6,'hourly forcing tar.gz');
local_assert_valid_gzip(tgzF_path);

tgzQ_name = U.hourly.flow_tgz_name;
tgzQ_url = U.hourly.flow_tgz_url;
tgzQ_path = fullfile(dirTmp,tgzQ_name);
local_download_if_needed(tgzQ_url, ...
    tgzQ_path,optsDL,logf, ...
    U.hourly.flow_size_bytes, ...
    'hourly streamflow','',1e6);

local_assert_file_nontrivial(tgzQ_path, ...
    1e6,'hourly streamflow tar.gz');
local_assert_valid_gzip(tgzQ_path);

logf('Hourly downloads complete.');

% Extract forcing + streamflow
logf('Extracting hourly forcing (tar.gz) ...');
dirExtractF = local_extract_targz(tgzF_path, ...
    dirHourly,logf);

logf('Extracting hourly streamflow (tar.gz) ...');
dirExtractQ = local_extract_targz(tgzQ_path, ...
    dirHourly,logf);

logf('Hourly extraction complete.');

% Rename to SAGE conventions
candF = { ...
    fullfile(dirHourly,'nldas_hourly'), ...
    fullfile(dirHourly,'nldas_hourly_csv'), ...
    fullfile(dirHourly,'nldas_hourly_csv', ...
    'nldas_hourly'), ...
    dirExtractF};
candQ = { ...
    fullfile(dirHourly,'usgs_streamflow'), ...
    fullfile(dirHourly,'usgs-streamflow_csv'), ...
    fullfile(dirHourly,'usgs_streamflow_csv'), ...
    dirExtractQ};

srcF = local_first_existing_dir(candF);
srcQ = local_first_existing_dir(candQ);

if isempty(srcF)
    logf(['WARNING: could not locate ' ...
        'extracted forcing folder; ' ...
        'leaving as-is.']);
else
    if exist(dirForcing,'dir') && opts.force
        rmdir(dirForcing,'s');
    end
    if ~exist(dirForcing,'dir')
        logf(['Renaming extracted forcing ' ...
            'folder -> hourly/forcing']);
        movefile(srcF,dirForcing);
    end
end

if isempty(srcQ)
    logf(['WARNING: could not locate ' ...
        'extracted streamflow folder; ' ...
        'leaving as-is.']);
else
    if exist(dirFlow,'dir') && opts.force
        rmdir(dirFlow,'s');
    end
    if ~exist(dirFlow,'dir')
        logf(['Renaming extracted streamflow folder ' ...
            '-> hourly/streamflow']);
        movefile(srcQ,dirFlow);
    end
end

% Final sanity + cleanup
forcingOK2 = local_dir_has_files(dirForcing, ...
    '*_hourly_nldas.csv',1);
flowOK2 = local_dir_has_files(dirFlow,   ...
    '*-usgs-hourly.csv',1);

if ~forcingOK2
    logf(['WARNING: hourly forcing ' ...
        'folder present ' ...
        'but no *_hourly_nldas.csv ' ...
        'files found: ' dirForcing]);
end
if ~flowOK2
    logf(['WARNING: hourly streamflow ' ...
        'folder present ' ...
        'but no *-usgs-hourly.csv ' ...
        'files found: ' dirFlow]);
end

% --------------------------------------
% Cleanup leftover extracted directories
% --------------------------------------
local_cleanup_hourly_leftovers(dirHourly, ...
    dirForcing,dirFlow,logf);

if ~opts.keepArchives
    local_try_delete(tgzF_path);
    local_try_delete(tgzQ_path);
end

logf('Hourly install complete.');

end

function ok = local_dir_has_files(patDir, ...
    pat,minCount)
ok = false;
try
    if nargin < 3, minCount = 1; end
    if isempty(patDir) ...
            || ~exist(patDir,'dir')
        return
    end
    dd = dir(fullfile(patDir,pat));
    ok = numel(dd) >= minCount;
catch
    ok = false;
end
end

function ok = local_tree_has_files(patDir,pat,minCount)
ok = false;
try
    if nargin < 3
        minCount = 1;
    end
    if isempty(patDir) ...
            || ~exist(patDir,'dir')
        return
    end
    dd = dir(fullfile(patDir,'**',pat));
    dd = dd(~[dd.isdir]);
    ok = numel(dd) >= minCount;
catch
    ok = false;
end
end

% ---------------
% Local utilities
% ---------------
function opts = local_set_default_opts(opts)
if ~isfield(opts,'force'), opts.force = 0; end
if ~isfield(opts,'tmpDir'), opts.tmpDir = ''; end
if ~isfield(opts,'keepArchives'), opts.keepArchives = 0; end
if ~isfield(opts,'verbose'), opts.verbose = 1; end
if ~isfield(opts,'timeoutSec'), opts.timeoutSec = 60; end
if ~isfield(opts,'logFcn') || isempty(opts.logFcn)
    opts.logFcn = @(s) disp(s);
end
if ~isfield(opts,'progressFcn')
    opts.progressFcn = [];   % @(info) ...
end
end

function local_download_if_needed(url, ...
    filepath,opts,logf,expectedBytes, ...
    label,expectedMD5,minBytes)
%LOCAL_DOWNLOAD_IF_NEEDED Download file with safe temp-file staging
%
% expectedBytes OPTIONAL: numeric bytes, fallback if Content-Length missing
% label OPTIONAL: string for progress messages
if nargin < 5 ...
        || isempty(expectedBytes)
    expectedBytes = NaN;
end
if nargin < 6 ...
        || isempty(label)
    label = 'download';
end
if nargin < 7 ...
        || isempty(expectedMD5)
    expectedMD5 = '';
end
if nargin < 8 ...
        || isempty(minBytes)
    minBytes = 1e6;
end
% skip if exists and not forcing
if exist(filepath,'file') ...
        && ~opts.force
    reuseOK = true;

    if isfinite(expectedBytes) ...
            && expectedBytes > 0
        s = dir(filepath);
        if isempty(s) || ~isfinite(s.bytes) || ...
                abs(double(s.bytes) - ...
                double(expectedBytes)) > 0.01*double(expectedBytes)
            reuseOK = false;
        end
    end

    if reuseOK && ~isempty(expectedMD5)
        try
            logf(['Computing MD5 for ' ...
                'existing file: ' filepath]);
            md5_here = local_md5_file(filepath);
            logf(['Finished MD5 for ' ...
                'existing file: ' filepath]);
            if ~strcmpi(md5_here,expectedMD5)
                reuseOK = false;
                logf(['Archive MD5 mismatch ' ...
                    '-> re-download: ' filepath]);
            end
        catch ME
            reuseOK = false;
            logf(['Could not compute MD5 ' ...
                '-> re-download: ' ME.message]);
        end
    end

    if reuseOK
        logf(['Archive exists and passed ' ...
            'validation (skip download): ' filepath]);
        return
    else
        opts.force = 1;
    end
end

% delete old final file if forcing
if exist(filepath,'file') ...
        && opts.force
    try
        delete(filepath);
    catch
    end
end

% temp partial filename
tmpPath = [filepath '.part'];

% delete old temp file if present
if exist(tmpPath,'file')
    try
        delete(tmpPath);
    catch
    end
end

logf(['Downloading: ' url]);
logf(['      -> ' filepath]);

% ensure folder
[fp,~,~] = fileparts(filepath);
if ~exist(fp,'dir')
    mkdir(fp);
end

% download to temp file, then rename atomically
local_download_stream(url, ...
    tmpPath,opts,logf,expectedBytes,label);

% validate downloaded temp file before promoting to final name
local_assert_file_nontrivial(tmpPath, ...
    minBytes,[label ' archive']);

if contains(lower(filepath),'.zip')
    local_assert_valid_zip(tmpPath);
elseif endsWith(lower(filepath),'.gz')
    local_assert_valid_gzip(tmpPath);
end

% promote temp file to final filename
try
    if exist(filepath,'file')
        delete(filepath);
    end
catch
end

ok = movefile(tmpPath,filepath,'f');
if ~ok
    error(['Could not rename ' ...
        'temp download to final file: %s'], ...
        filepath);
end
end

function local_download_stream(url, ...
    filepath,opts,logf,expectedBytes,label)
%LOCAL_DOWNLOAD_STREAM Robust download; prefer curl on Windows

if ispc
    try
        local_download_stream_curl(url, ...
            filepath,opts,logf,expectedBytes,label);
        return
    catch ME
        logf(['curl download failed; ' ...
            'falling back to Java: ' ME.message]);
    end
end

% Fallback to legacy Java downloader
local_download_stream_java(url, ...
    filepath,opts,logf,expectedBytes,label);
end

function local_download_stream_java(url, ...
    filepath,opts,logf,expectedBytes,label)
%LOCAL_DOWNLOAD_STREAM_JAVA Stream download using Java to update progress

% open connection
u = java.net.URL(url);
c = u.openConnection();
try
    c.setConnectTimeout(1000*opts.timeoutSec);
    c.setReadTimeout(1000*opts.timeoutSec);
catch
end

cancelFcn = [];
if isfield(opts,'cancelFcn') ...
        && ~isempty(opts.cancelFcn)
    cancelFcn = opts.cancelFcn;
end

% Content-Length (bytes); may be -1
try
    total = double(c.getContentLengthLong());
catch
    total = NaN;
end
if ~(isfinite(total) && total > 0) ...
        && isfinite(expectedBytes) ...
        && expectedBytes > 0
    total = double(expectedBytes);
end

% streams
in = c.getInputStream();
out = java.io.FileOutputStream(filepath);

% download loop
bufSize = 1024*1024;                % 1 MB chunks
buf = int8(zeros(bufSize,1));
%nread = 0;

t0 = tic;
tLast = tic;
bytes = 0;

% initial ping
local_progress_ping(opts, ...
    label,bytes,total,t0,logf,true);

try
    while true

        % user cancel?
        if ~isempty(cancelFcn)
            try
                if logical(cancelFcn())
                    error(['Download canceled ' ...
                        'by user.']);
                end
            catch MEc
                error(MEc.message);
            end
        end
        n = in.read(buf,0,bufSize);
        if n < 0
            break
        end
        out.write(buf,0,n);

        bytes = bytes + double(n);

        % throttle UI/log updates
        if toc(tLast) >= 0.35
            local_progress_ping(opts, ...
                label,bytes,total, ...
                t0,logf,false);
            tLast = tic;
        end
    end
catch ME
    try out.close(); catch, end
    try in.close(); catch, end
    % delete partial file
    try
        if exist(filepath,'file')
            delete(filepath);
        end
    catch
    end
    error(['Download failed for: ' url newline ...
           'Message: ' ME.message]);
end

% close
try out.close(); catch, end
try in.close(); catch, end

% final ping
local_progress_ping(opts, ...
    label,bytes,total,t0,logf,true);

% post-download size sanity
if isfinite(expectedBytes) ...
        && expectedBytes > 0
    d = dir(filepath);
    if isempty(d) || ~isfinite(d.bytes)
        error(['Downloaded file missing ' ...
            'after transfer: %s'], ...
            filepath);
    end
    if abs(double(d.bytes) ...
            - double(expectedBytes)) > 0.01*double(expectedBytes)
        error(['Downloaded file size mismatch ' ...
            'for %s: got %.0f bytes, ' ...
            'expected about %.0f bytes'], ...
            filepath,double(d.bytes), ...
            double(expectedBytes));
    end
end

logf([upper(label) ' download complete.']);
end

function local_download_stream_curl(url, ...
    filepath,opts,logf,expectedBytes,label)
%LOCAL_DOWNLOAD_STREAM_CURL Download using Windows curl.exe

if exist(filepath,'file')
    try
        delete(filepath);
    catch
    end
end

[fp,~,~] = fileparts(filepath);
if ~exist(fp,'dir')
    mkdir(fp);
end

cmd = sprintf(['curl.exe -L --fail --retry 3 --retry-delay 2 ' ...
    '--retry-all-errors --silent --show-error --ssl-no-revoke ' ...
    '--output "%s" "%s"'],filepath,url);
logf(['Downloading with curl: ' url]);

if isfield(opts,'progressFcn') ...
        && ~isempty(opts.progressFcn)
    try
        info = struct( ...
            'label',label, ...
            'bytes',NaN, ...
            'totalBytes',NaN, ...
            'frac',NaN, ...
            'speedMBs',NaN, ...
            'etaSec',NaN, ...
            'indeterminate',true);
        opts.progressFcn(info);
    catch
    end
end

try
    [status,txt] = system(cmd);

    if status ~= 0
        error('curl failed (%d): %s', ...
            status,strtrim(txt));
    end

    % post-download size sanity
    d = dir(filepath);
    if isempty(d) || ~isfinite(d.bytes)
        error(['Downloaded file ' ...
            'missing after curl transfer: %s'], ...
            filepath);
    end

    if isfinite(expectedBytes) ...
            && expectedBytes > 0
        if abs(double(d.bytes) - ...
                double(expectedBytes)) > 0.01*double(expectedBytes)
            error(['Downloaded file size ' ...
                'mismatch for %s: got %.0f bytes, ' ...
                'expected about %.0f bytes'], ...
                filepath,double(d.bytes), ...
                double(expectedBytes));
        end
    end

    % archive-type validation
    if contains(lower(filepath),'.zip')
        local_assert_valid_zip(filepath);
    elseif endsWith(lower(filepath),'.gz')
        local_assert_valid_gzip(filepath);
    end

catch ME
    pause(0.5);
    try
        if exist(filepath,'file')
            delete(filepath);
        end
    catch
        logf(['WARNING: could not ' ...
            'delete partial file: ' filepath]);
    end
    rethrow(ME);
end

logf([upper(label) ' download complete.']);
end

function local_progress_ping(opts, ...
    label,bytes,total,t0,logf,forcePrint)
%LOCAL_PROGRESS_PING Send progress to opts.progressFcn and/or logf

dt = toc(t0);
if dt <= 0
    dt = eps;
end
speedMBs = (bytes/1024/1024) / dt;

frac = NaN;
etaSec = NaN;
if isfinite(total) && total > 0
    frac = min(max(bytes/total,0),1);
    if speedMBs > 0
        remMB = (total-bytes)/1024/1024;
        etaSec = max(remMB/speedMBs,0);
    end
end

% callback (GUI)
if isfield(opts,'progressFcn') ...
        && ~isempty(opts.progressFcn)
    try
        info = struct('label',label, ...
            'bytes',bytes, ...
            'totalBytes',total, ...
            'frac',frac, ...
            'speedMBs',speedMBs, ...
            'etaSec',etaSec);
        opts.progressFcn(info);
    catch
    end
end

% optional console/log printing
if opts.verbose

    % initialize / keep state
    persistent lastPct %#ok
    if isempty(lastPct)
        lastPct = -1;
    end

    % print at start/end always
    if forcePrint ...
            || (isfinite(frac) ...
            && (frac==0 || frac==1))

        lastPct = -1;   % reset for next download / next file

        if isfinite(frac)
            logf(sprintf('%s: %5.1f%%%%  (%.2f MB/s)', ...
                label,100*frac,speedMBs));
        else
            logf(sprintf('%s: %s  (%.2f MB/s)', ...
                label,'... ',speedMBs));
        end

        return
    end

    % otherwise: print at 10%% increments
    if isfinite(frac)
        pct10 = floor((100*frac)/10);     % 0..10 (each step = 10%)
        if pct10 ~= lastPct
            lastPct = pct10;
            logf(sprintf('%s: %5.1f%%%%  (%.2f MB/s)', ...
                label,100*frac,speedMBs));
        end
    end

end

end

function local_safe_unzip(zipPath, ...
    dstDir,logf)
if ~exist(zipPath,'file')
    error('Zip file not found: %s', ...
        zipPath);
end
if ~exist(dstDir,'dir')
    mkdir(dstDir);
end

local_assert_file_nontrivial(zipPath, ...
    1e6,'ZIP archive');
local_assert_valid_zip(zipPath);

try
    logf(['Starting unzip of: ' zipPath]);
    unzip(zipPath,dstDir);
    logf('Unzip done.');
catch ME
    error(['unzip failed: ' zipPath ...
        newline ME.message]);
end
logf('Unzip done.');
end

function outDir = local_extract_targz(tgzPath, ...
    dstDir,logf)
% returns the extraction folder we can guess (may be dstDir)
outDir = '';

if ~exist(tgzPath,'file')
    error('Archive not found: %s',tgzPath);
end
local_assert_file_nontrivial(tgzPath, ...
    1e6,'tar.gz archive');
local_assert_valid_gzip(tgzPath);

if ~exist(dstDir,'dir')
    mkdir(dstDir);
end

% gunzip -> .tar
[~,~,ext] = fileparts(tgzPath);
if strcmpi(ext,'.gz')
    logf(['  gunzip: ' tgzPath]);
    try
        tarFiles = gunzip(tgzPath,dstDir);
    catch ME
        error(['gunzip failed: ' tgzPath ...
            newline ME.message]);
    end
    if isempty(tarFiles) ...
            || ~exist(tarFiles{1},'file')
        error(['gunzip did not produce ' ...
            'a .tar file for: %s'],tgzPath);
    end
    tarPath = tarFiles{1};
else
    tarPath = tgzPath;
end

% untar
logf(['  untar: ' tarPath]);
try
    listing = untar(tarPath,dstDir);
catch ME
    error(['untar failed: ' tarPath ...
        newline ME.message]);
end

% attempt to infer top folder from listing
if ~isempty(listing)
    % listing can be files; infer common root
    top = local_infer_common_root(listing);
    if ~isempty(top) ...
            && exist(top,'dir')
        outDir = top;
    end
end

% delete tar (optional)
try
    if exist(tarPath,'file') ...
            && ~strcmpi(tarPath,tgzPath)
        delete(tarPath);
    end
catch
end

logf('  extract done.');
end

function top = local_infer_common_root(listing)
% listing is cellstr of extracted paths (absolute or relative)
top = '';
try
    p = listing{1};
    if iscell(p), p = p{1}; end
    p = char(p);
    % take first folder in path
    [d1,~,~] = fileparts(p);
    if isempty(d1)
        return
    end
    % if this is a file path, its parent might be the top folder
    top = d1;
catch
    top = '';
end
end

function dirFound = local_find_subdir(parentDir, ...
    folderName)
dirFound = '';
try
    dd = dir(parentDir);
    for i = 1:numel(dd)
        if dd(i).isdir ...
                && ~strcmp(dd(i).name,'.') ...
                && ~strcmp(dd(i).name,'..')
            if strcmp(dd(i).name,folderName)
                dirFound = fullfile(parentDir, ...
                    dd(i).name);
                return
            end
        end
    end
catch
end
end

function p = local_first_existing_dir(cands)
p = '';
for i = 1:numel(cands)
    if ~isempty(cands{i}) ...
            && exist(cands{i},'dir')
        p = cands{i};
        return
    end
end
end

function local_try_delete(fp)
try
    if exist(fp,'file')
        delete(fp);
    end
catch
end
end

function local_assert_file_nontrivial(fp, ...
    minBytes,label)
%LOCAL_ASSERT_FILE_NONTRIVIAL Basic existence + size sanity check

if nargin < 2 ...
        || isempty(minBytes)
    minBytes = 1;
end
if nargin < 3 ...
        || isempty(label)
    label = 'archive';
end

if ~exist(fp,'file')
    error('%s not found: %s', ...
        label,fp);
end

d = dir(fp);
if isempty(d) ...
        || ~isfinite(d.bytes)
    error(['Could not determine size ' ...
        'of %s: %s'], ...
        label,fp);
end

if d.bytes < minBytes
    error(['%s looks too small to be ' ...
        'valid: %s ' ...
        '(%.0f bytes)'], ...
        label,fp,double(d.bytes));
end
end

function local_assert_valid_zip(zipPath)
%LOCAL_ASSERT_VALID_ZIP Check ZIP signature before unzip

fid = fopen(zipPath,'r');
if fid < 0
    error('Could not open ZIP file: %s', ...
        zipPath);
end

cleanup = onCleanup(@() fclose(fid));

sig = fread(fid,4,'uint8=>uint8')';
if numel(sig) < 4
    error(['ZIP file is too ' ...
        'short / truncated: %s'], ...
        zipPath);
end

% Common ZIP signatures:
%   50 4B 03 04  local file header
%   50 4B 05 06  empty archive end-of-central-directory
%   50 4B 07 08  spanned/special
if ~(sig(1)==80 && sig(2)==75)
    headTxt = char(sig);
    headTxt(~isstrprop(headTxt,'print')) = '.';
    error(['Downloaded file is not a ' ...
        'valid ZIP archive: %s ' ...
        '(header="%s")'], ...
        zipPath,headTxt);
end
end

function local_assert_valid_gzip(gzPath)
%LOCAL_ASSERT_VALID_GZIP Check GZIP signature before gunzip

fid = fopen(gzPath,'r');
if fid < 0
    error('Could not open gzip file: %s', ...
        gzPath);
end

cleanup = onCleanup(@() fclose(fid));

sig = fread(fid,3,'uint8=>uint8')';
if numel(sig) < 3
    error(['Gzip file is too ' ...
        'short / truncated: %s'], ...
        gzPath);
end

% GZIP magic bytes: 1F 8B 08
if ~(sig(1)==31 ...
        && sig(2)==139 ...
        && sig(3)==8)
    headTxt = char(sig);
    headTxt(~isstrprop(headTxt,'print')) = '.';
    error(['Downloaded file is not a ' ...
        'valid GZIP archive: %s ' ...
        '(header="%s")'],gzPath,headTxt);
end
end

% ------------
% URL registry
% ------------
function U = camels_urls()
%CAMELS_URLS  Central place for download URLs + filenames.
%
% You can update links here without touching installer logic.

U = struct();

% -----------------------------
% Daily CAMELS-US (v1p2) zip
% -----------------------------
U.daily = struct();
U.daily.record = 'https://zenodo.org/records/15529996';
U.daily.zip_name = 'basin_timeseries_v1p2_metForcing_obsFlow.zip';
% NOTE:
% Zenodo supports stable "records/<id>/files/<filename>" style.
% If Zenodo changes, update this one line.
U.daily.zip_url = [ ...
    'https://zenodo.org/records/15529996/files/' ...
    U.daily.zip_name '?download=1' ];

U.daily.size_bytes = 3326784 * 1024;  % basin_timeseries... (KB -> bytes)
U.daily.md5 = '8e9a466710e8270b58f01d332a87184f';

% -----------------------------------------
% Hourly CAMELS-US (NeuralHydrology) tar.gz
% -----------------------------------------
U.hourly = struct();
U.hourly.record = 'https://zenodo.org/records/4072701';

U.hourly.forcing_tgz_name = 'nldas_hourly_csv.tar.gz';
U.hourly.flow_tgz_name = 'usgs_streamflow_csv.tar.gz';

U.hourly.forcing_tgz_url = [ ...
    'https://zenodo.org/records/4072701/files/' ...
    U.hourly.forcing_tgz_name '?download=1' ];

U.hourly.flow_tgz_url = [ ...
    'https://zenodo.org/records/4072701/files/' ...
    U.hourly.flow_tgz_name '?download=1' ];

U.hourly.forcing_size_bytes = 16578615 * 1024; % nldas_hourly_csv.tar.gz
U.hourly.flow_size_bytes = 1011629 * 1024; % usgs_streamflow_csv.tar.gz

% ----------------------------------
% Shared CAMELS attribute text files
% Installed into root/Data/
% ----------------------------------
U.shared = struct();

U.shared.files = [ ...
    struct('name','camels_clim.txt', ...
           'url',['https://zenodo.org/records/15529996/' ...
           'files/camels_clim.txt?download=1'], ...
           'size_bytes',NaN), ...
    struct('name','camels_geol.txt', ...
           'url',['https://zenodo.org/records/15529996/' ...
           'files/camels_geol.txt?download=1'], ...
           'size_bytes',NaN), ...
    struct('name','camels_hydro.txt', ...
           'url',['https://zenodo.org/records/15529996/' ...
           'files/camels_hydro.txt?download=1'], ...
           'size_bytes',NaN), ...
    struct('name','camels_name.txt', ...
           'url',['https://zenodo.org/records/15529996/' ...
           'files/camels_name.txt?download=1'], ...
           'size_bytes',NaN), ...
    struct('name','camels_soil.txt', ...
           'url',['https://zenodo.org/records/15529996/' ...
           'files/camels_soil.txt?download=1'], ...
           'size_bytes',NaN), ...
    struct('name','camels_topo.txt', ...
           'url',['https://zenodo.org/records/15529996/' ...
           'files/camels_topo.txt?download=1'], ...
           'size_bytes',NaN), ...
    struct('name','camels_vege.txt', ...
           'url',['https://zenodo.org/records/15529996/' ...
           'files/camels_vege.txt?download=1'], ...
           'size_bytes',NaN)];
end

function local_cleanup_hourly_leftovers(dirHourly, ...
    dirForcing,dirFlow,logf)
%LOCAL_CLEANUP_HOURLY_LEFTOVERS Remove empty/leftover extraction folders

cands = { ...
    fullfile(dirHourly,'nldas_hourly'), ...
    fullfile(dirHourly,'nldas_hourly_csv'), ...
    fullfile(dirHourly,'usgs_streamflow'), ...
    fullfile(dirHourly,'usgs-streamflow_csv'), ...
    fullfile(dirHourly,'usgs_streamflow_csv') ...
    };

for i = 1:numel(cands)
    p = cands{i};

    % never delete the final targets
    if strcmpi(p,dirForcing) || strcmpi(p,dirFlow)
        continue
    end

    if exist(p,'dir')
        try
            % if folder is empty, remove; otherwise remove only if forced
            if local_is_empty_dir(p)
                logf(['Deleting leftover ' ...
                    'empty folder: ' p]);
                rmdir(p,'s');
            end
        catch
        end
    end
end

end

function tf = local_is_empty_dir(p)
%LOCAL_IS_EMPTY_DIR True if directory contains no files/folders

try
    dd = dir(p);
    names = {dd.name};
    names = names(~ismember(names,{'.','..'}));
    tf = isempty(names);
catch
    tf = false;
end
end

function dirTmp = local_default_download_dir()
%LOCAL_DEFAULT_DOWNLOAD_DIR Return the user's Downloads folder

try
    homeDir = char(java.lang.System.getProperty('user.home'));
    cand = fullfile(homeDir, ...
        'Downloads');
    if exist(cand,'dir')
        dirTmp = cand;
        return
    end
catch
end

try
    userprof = getenv('USERPROFILE');
    cand = fullfile(userprof, ...
        'Downloads');
    if ~isempty(userprof) ...
            && exist(cand,'dir')
        dirTmp = cand;
        return
    end
catch
end

try
    dirTmp = tempdir;
catch
    dirTmp = pwd;
end
end

function local_install_shared_metadata(dirData,dirTmp,U,opts,logf)
%LOCAL_INSTALL_SHARED_METADATA Download CAMELS attribute txt files to Data/

logf(['--- Installing shared CAMELS ' ...
    'metadata text files ---']);

files = U.shared.files;

for i = 1:numel(files)
    name = files(i).name;
    url  = files(i).url;
    fp   = fullfile(dirData,name);

    % store downloaded archive/file first in Downloads/tmpDir
    tmpFp = fullfile(dirTmp,name);

    local_download_if_needed(url,tmpFp,opts,logf, ...
        files(i).size_bytes,name,'',100);
    local_assert_file_nontrivial(tmpFp, ...
        100,[name ' metadata file']);

    % promote/copy into Data/
    if ~exist(dirData,'dir')
        mkdir(dirData);
    end

    try
        copyfile(tmpFp,fp,'f');
        logf(['Installed shared metadata file: ' fp]);
    catch ME
        error(['Could not copy metadata file to Data: ' fp ...
            newline ME.message]);
    end
end

logf('Shared CAMELS metadata text files installed.');
end

function local_copy_gauge_information_to_data(dirData,logf)
%LOCAL_COPY_GAUGE_INFORMATION_TO_DATA Copy gauge_information.txt to Data/

src = fullfile(dirData,'daily','v1p2', ...
    'metadata','gauge_information.txt');
dst = fullfile(dirData,'gauge_information.txt');

if exist(src,'file')
    try
        copyfile(src,dst,'f');
        logf(['Copied gauge_information.txt to: ' dst]);
    catch ME
        logf(['WARNING: could not copy ' ...
            'gauge_information.txt to Data: ' ...
            ME.message]);
    end
else
    logf(['WARNING: gauge_information.txt ' ...
        'not found at expected location: ' ...
        src]);
end
end

function h = local_md5_file(fp)
%LOCAL_MD5_FILE Compute MD5 hash of a file, lowercase hex string

import java.security.*
import java.io.*

md = MessageDigest.getInstance('MD5');
fis = FileInputStream(java.io.File(fp));
dis = DigestInputStream(fis,md);

buf = zeros(1024*1024,1,'int8');   % 1 MB buffer

try
    while dis.read(buf,0,numel(buf)) ~= -1
    end
    dis.close();
catch ME
    try dis.close(); catch, end
    rethrow(ME);
end

raw = typecast(md.digest(),'uint8');
h = lower(reshape(dec2hex(raw)',1,[]));
end