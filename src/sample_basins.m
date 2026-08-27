function [A,bas,latlon] = sample_basins(A_camels,id_gauge,bas, ...
    prd,gname,zone,dirD,file_u,file_s)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAMPLE_BASINS Divides the CONUS watersheds into training and evaluation
% basins
%  SYNOPSIS: [A,bas,latlon] = sample_basins(A_camels,id_gauge,bas,prd, ...
%   gname,zone,dirD,file_u,file_s)
%   A_camels    rxK_all matrix catchment attributes all CONUS watersheds
%                [] -> skip attribute handling (SITE training)
%   id_gauge    K_allx1 vector gauge catchment codes
%   bas         structure basin information
%    .K          total number of allowed basins
%                [671 = all daily CAMELS / 531 = restricted daily /
%                499 = hourly]
%    .K_t        number of training watersheds
%    .K_e        number of evaluation watersheds (OPTIONAL; default = 0)
%    .sample     basin sampling mode
%                'random' -> random train/evaluation selection
%                'file'   -> user-defined file order
%   prd         structure about training/evaluation/spin-up period
%    .dt         Temporal data resolution [1:daily, 24:hourly]
%    .dts        First day of training [d/m/y]
%    .dte        Last day of training [d/m/y]
%    .des        First day of evaluation [d/m/y]
%    .dee        Last day of evaluation [d/m/y]
%    .spinup     Spin up period in days
%   gname       K_allx1 vector gauge basin names, state
%   zone        INPUT: structure with hydroclimatic basin classification
%               for all basins in id_gauge before sampling
%    .id         (K_t+K_e)x1 string array with zone labels per basin
%                e.g., "humid_rain", "subhumid_snow", "dry_rain"
%    .num        (K_t+K_e)x1 numeric vector with integer zone identifiers
%    .names      mx1 string array with unique zone names
%    .aridity    (K_t+K_e)x1 vector of aridity index values
%    .frac_snow  (K_t+K_e)x1 vector of snow-fraction values
%   dirD        Main data directory with gauge_information.txt metadat file
%   file_u      file with gauge codes of all available basins
%               (the basin universe; e.g., 531_basins.txt,
%               499_basins.txt, 671_basins.txt)
%   file_s      split file with gauge codes of training and
%               evaluation basins
%               first K_t entries     -> training basins
%               next  K_e entries     -> evaluation basins
%               used only if bas.sample = 'file'
%   A           OUTPUT: rxK matrix attributes K = K_t+K_e train/eval basins
%   bas         OUTPUT: updated structure basin information
%    .K          number of training and evaluation basins
%    .r          number of basin attributes
%    .id_t       K_tx1 vector integers training basins (in final list)
%    .id_e       K_ex1 vector integers evaluation basins (in final list)
%    .id_gauge    revised list gauge catchment codes (order: train; eval)
%    .gname      corresponding gauge names (order: train; eval)
%    .id_plot    display order [train sorted by gauge; eval sorted by gauge]
%    .mode       assessment mode returned by get_assess_mode
%    .zone       hydroclimatic zone structure for selected basins
%   latlon      OUTPUT: (K_t+K_e)x2 matrix lat/lon (°) of K basins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Mar. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 9
        file_s = [];
    end
    if nargin < 7 ...
            || isempty(dirD)
        fname_gauge = 'gauge_information.txt';
    else
        fname_gauge = fullfile(dirD, ...
            'gauge_information.txt');
    end
    
    fprintf('... Dividing basins into training/evaluation');
    
    % ---------------------
    % Unpack basin settings
    % ---------------------
    K_t = bas.K_t;
    
    if isfield(bas,'K_e') ...
            && ~isempty(bas.K_e)
    else
        bas.K_e = 0;
    end
    K_e = bas.K_e;
    
    if isfield(bas,'sample') ...
            && ~isempty(bas.sample)
        sp_mode = lower(strtrim(char(bas.sample)));
    else
        sp_mode = 'random';
    end
    
    K = K_t + K_e;
    
    % ------------------
    % Basic input checks
    % ------------------
    if K_t < 1
        fprintf('\n');
        error(['      Error: sample_basins: ' ...
            'bas.K_t = %d, ' ...
            'but K_t must be >= 1.'],K_t);
    end
    
    if K_e < 0
        fprintf('\n');
        error(['      Error: sample_basins: ' ...
            'bas.K_e = %d, but K_e must be >= 0.'],K_e);
    end
    
    switch sp_mode
        case {'random','zone','file'}
            % ok
        otherwise
            fprintf('\n');
            error(['      Error: sample_basins: ' ...
                'Unknown bas.sample = %s. ' ...
                'Use ''random'', ''zone'', or ''file''.'], ...
                sp_mode);
    end
    
    % -------------------------------
    % Standardize incoming dimensions
    % -------------------------------    
    if istable(gname)
        gname = gname{:,1};
    end
    if iscell(gname)
        gname = string(gname);
    elseif ischar(gname)
        gname = string(cellstr(gname));
    else
        gname = string(gname);
    end
    gname = gname(:);
    
    if local_has_prepared_basins(bas) ...
        && (nargin < 9 ...
        || isempty(file_s))
        [A,bas,latlon] = local_use_prepared_basins( ...
            A_camels,id_gauge,bas,prd,gname,zone,dirD);
        return
    end
    % --------------------------------------------
    % Determine allowed basin universe from file_u
    % --------------------------------------------
    if nargin < 8 ...
            || isempty(file_u)
        fprintf('\n');
        error(['      Error: sample_basins: ' ...
            'No basin-universe file was provided.']);
    end
    
    file_univ = local_resolve_basin_file( ...
        file_u,dirD);

    isINDfile = contains( ...
        upper(string(file_univ)), ...
        'IND_');
    
    id_gauge = local_id_key( ...
        id_gauge,isINDfile);
    id_univ  = local_load_vec( ...
        file_univ,isINDfile);
    
    % Harmonize regional basin-ID conventions.  CAMELS-DE basin-list
    % files use IDs such as DE110000, DEA10000, DEB10000, ... whereas
    % some attribute readers may still return 110000, A10000, B10000, ... .
    % If needed, add the leading DE prefix to the in-memory id_gauge vector
    % so matching and downstream file names remain consistent with the
    % basin/split files.
    [id_gauge,id_univ] = local_harmonize_de_ids(id_gauge,id_univ);
    
    % Keep only basins that are both:
    %   1) present in the requested universe file
    %   2) present in the provided id_gauge vector
    %
    % Preserve the order of the universe file
    [tf_univ,loc_gauge] = ismember(id_univ,id_gauge);
    id_keep = loc_gauge(tf_univ);
    
    if isempty(id_keep)
        fprintf('\n');
        error(['      Error: sample_basins: ' ...
            'No overlap between id_gauge and %s.'], ...
            file_univ);
    end
    
    id_gauge = id_gauge(id_keep);
    gname = gname(id_keep);
    
    if nargin < 6 ...
            || isempty(zone)
        zone_keep = [];
    else
        zone_keep = ...
            local_subset_zone(zone,id_keep);
    end
    
    if ~isempty(A_camels)
        A = A_camels(:,id_keep);
    else
        A = [];
    end
    
    K_avail = numel(id_gauge);
    
    if K > K_avail
        fprintf('\n');
        error(['      Error: sample_basins: ' ...
            'Requested K_t + K_e = %d basins, ' ...
            'but only %d are available in %s.'], ...
            K,K_avail,file_univ);
    end
    
    if isfield(bas,'K') ...
            && ~isempty(bas.K)
        if numel(id_univ) ~= bas.K
            warning(['sample_basins: ' ...
                'basin universe file contains %d IDs, ' ...
                'but bas.K = %d.'], ...
                numel(id_univ),bas.K);
        end
    end
    
    % ---------------------------------------------
    % Select training/evaluation within allowed set
    % ---------------------------------------------
    switch sp_mode
    
        case 'random'
            id_sel = randperm(K_avail,K);
    
            id_t = id_sel(1:K_t);
            id_e = id_sel(K_t+1:K);
    
            id_t = id_t(:);
            id_e = id_e(:);
    
        case 'zone'
            if isempty(zone_keep) ...
                    || ~isfield(zone_keep,'num') ...
                    || numel(zone_keep.num) ~= K_avail
                fprintf('\n');
                error(['      Error: sample_basins: ' ...
                    'bas.sample = ''zone'' requires ' ...
                    'a valid zone structure ' ...
                    'with zone.num for all ' ...
                    'available basins.']);
            end
            [id_t,id_e] = local_stratified_sample( ...
                zone_keep.num,K_t,K_e);
            id_t = id_t(:);
            id_e = id_e(:);
    
        case 'file'
            if isempty(file_s)
                fprintf('\n');
                error(['      Error: sample_basins: ' ...
                    'bas.sample = ''file'' but no ' ...
                    'split file was provided.']);
            end
            file_s = local_resolve_basin_file(file_s,dirD);
            id_file_raw = local_load_vec(file_s,isINDfile);
            id_file = local_normalize_split_ids(id_file_raw, id_gauge);
            n_file = numel(id_file);
    
            if n_file < K
                fprintf('\n');
                error(['      Error: sample_basins: ' ...
                    'Split file %s contains %d basin ' ...
                    'IDs, but K_t + K_e = %d ' ...
                    'are required.'],file_s,n_file,K);
            end
    
            id_req_t = id_file(1:K_t);
            id_req_e = id_file(K_t+1:K);
    
            id_req = [id_req_t(:); id_req_e(:)];
            if numel(unique(id_req)) ~= numel(id_req)
                fprintf('\n');
                error(['      Error: sample_basins: ' ...
                    'Split file %s contains duplicate ' ...
                    'basin IDs in the requested ' ...
                    'training/evaluation selection.'], ...
                    file_s);
            end
    
            [tf_t,loc_t] = ismember(id_req_t,id_gauge);
            [tf_e,loc_e] = ismember(id_req_e,id_gauge);
    
            if ~all(tf_t)
                nmiss = sum(~tf_t);
                fprintf('\n');
                error(['      Error: sample_basins: ' ...
                    '%d requested training basin(s) ' ...
                    'from file %s were not found ' ...
                    'in the selected universe %s.'], ...
                    nmiss,file_s,file_univ);
            end
    
            if ~all(tf_e)
                nmiss = sum(~tf_e);
                fprintf('\n');
                error(['      Error: sample_basins: ' ...
                    '%d requested evaluation basin(s) ' ...
                    'from file %s were not found ' ...
                    'in the selected universe %s.'], ...
                    nmiss,file_s,file_univ);
            end
    
            % Preserve exact file order
            id_t = loc_t(:);
            id_e = loc_e(:);
    end
    
    % ------------------------------------------
    % Final output order is always [train; eval]
    % ------------------------------------------
    id_out = [id_t; id_e];
    
    if ~isempty(A)
        A = A(:,id_out);
        bad_attr = any(~isfinite(A),2);
        bad_basin = any(~isfinite(A),1);
        if any(bad_attr)
            bad_idx = find(bad_attr);

            if isfield(bas,'id_attr') ...
                    && numel(bas.id_attr) == numel(bad_attr)
                id_bad_attr = bas.id_attr(bad_attr);
            else
                id_bad_attr = bad_idx;
            end
            fprintf('\n');
            error(['      Error: sample_basins: ' ...
                'Selected attributes contain ' ...
                'NaN/Inf values ' ...
                'after combining with the ' ...
                'selected basin file. ' ...
                'Remove attribute ID(s): %s.'], ...
                mat2str(id_bad_attr));
        end
        if any(bad_basin)
            id_bad_basin = id_gauge( ...
                id_out(bad_basin));
            fprintf('\n');
            error(['      Error: sample_basins: ' ...
                'Selected basin(s) contain ' ...
                'NaN/Inf attribute values. ' ...
                'Example basin ID: %s.'], ...
                local_id_str(id_bad_basin(1)));
        end
    else
        A = [];
    end
    id_gauge = id_gauge(id_out);
    % Store CAMELS-DE basin IDs without leading DE
    if contains(upper(string(dirD)),"CAMELS_DE") ...
            || contains(upper(string(file_univ)),"\DE\") ...
            || contains(upper(string(file_univ)),"/DE/")
        id_gauge = local_strip_de_prefix(id_gauge);
    end

    gname = gname(id_out);
    if ~isempty(zone_keep)
        zone_out = ...
            local_subset_zone(zone_keep,id_out);
    else
        zone_out = [];
    end
    
    % ------------------------------------------
    % Read latitude/longitude of selected basins
    % ------------------------------------------
    if ~isempty(zone_out) ...
            && isfield(zone_out,'lat') ...
            && isfield(zone_out,'lon') ...
            && numel(zone_out.lat) == K ...
            && numel(zone_out.lon) == K
    
        latlon = [double(zone_out.lat(:)) ...
                  double(zone_out.lon(:))];
    
    else
        if ~isfile(fname_gauge)
            fprintf('\n');
            error(['      Error: sample_basins: ' ...
                'Cannot find gauge_information.txt ' ...
                'at:\n      %s'],fname_gauge);
        end
    
        T = readtable(fname_gauge, ...
            'Delimiter','\t', ...
        'VariableNamingRule','preserve');
    
        vn = lower(strtrim(string( ...
            T.Properties.VariableNames)));
        
        iGauge = find(ismember(vn, ...
            ["gauge_id","gauge","gaugeid", ...
            "id","station_id"]),1);
        
        if isempty(iGauge)
            if width(T) >= 2
                iGauge = 2;
            else
                iGauge = 1;
            end
        end
        
        id_meta = local_id_key(T{:,iGauge});
        id_meta = local_harmonize_de_reference( ...
            id_gauge,id_meta);
        
        [tf_meta,loc_meta] = ...
            ismember(id_gauge,id_meta);
        
        if ~all(tf_meta)
            nmiss = sum(~tf_meta);
            fprintf('\n');
            error(['      Error: sample_basins: ' ...
                '%d selected basin(s) not found ' ...
                'in gauge_information.txt.'],nmiss);
        end
        
        iLat = find(ismember(vn, ...
            ["lat","latitude","lat_wgs84", ...
            "latitude_wgs84"]),1);
        iLon = find(ismember(vn, ...
            ["lon","long","longitude", ...
            "lon_wgs84","longitude_wgs84"]),1);
        
        if ~isempty(iLat) ...
                && ~isempty(iLon)
            latlon = [double(T{loc_meta,iLat}) ...
                      double(T{loc_meta,iLon})];
        elseif width(T) >= 5
            latlon = double(T{loc_meta,[4 5]});
        elseif width(T) >= 4
            latlon = double(T{loc_meta,[3 4]});
        else
            latlon = nan(numel(id_gauge),2);
        end
    end
    
    % ------------------------------
    % Final train/eval index vectors
    % ------------------------------
    id_t = (1:K_t).';
    id_e = (K_t+1:K).';
    
    % ---------------------------------------------------------------------
    % Display order: first training basins sorted by gauge, then evaluation
    % basins sorted by gauge
    % ---------------------------------------------------------------------
    id_plot_t = (1:K_t).';
    if ~isempty(id_plot_t)
        [~,ord_t] = sort(id_gauge(id_plot_t));
        id_plot_t = id_plot_t(ord_t);
    end
    
    id_plot_e = (K_t+1:K).';
    if ~isempty(id_plot_e)
        [~,ord_e] = sort(id_gauge(id_plot_e));
        id_plot_e = id_plot_e(ord_e);
    end
    
    id_plot = [id_plot_t; id_plot_e];
    
    % -------------------------
    % Determine assessment mode
    % -------------------------
    try
        mode = get_assess_mode(prd,bas);
    catch
        mode = 'none';
    end
    
    % --------------------
    % Number of attributes
    % --------------------
    if isfield(bas,'id_attr') ...
            && ~isempty(bas.id_attr)
        r = numel(bas.id_attr);
    else
        r = 0;
    end
    
    % ----------------------
    % Update basin structure
    % ----------------------
    add = struct('K',K, ...
        'r',r, ...
        'id_t',id_t, ...
        'id_e',id_e, ...
        'id_gauge',id_gauge(:), ...
        'gname',gname(:), ...
        'zone',zone_out, ...
        'id_plot',id_plot(:), ...
        'mode',mode);
    
    fn = fieldnames(add);
    for i = 1:numel(fn)
        bas.(fn{i}) = add.(fn{i});
    end
    
    fprintf(' ... Done\n');

end

% ------------------------------------------------------
% helper: robust load of a numeric vector from .txt/.mat
% ------------------------------------------------------
function v = local_load_vec(fname,pad5)

    if nargin < 2
        pad5 = false;
    end

    if ~isfile(fname)
        error(['      Error: sample_basins: ' ...
            'Cannot find file:' ...
            '\n      %s'],fname);
    end
    
    [~,~,ext] = fileparts(fname);
    ext = lower(ext);
    
    switch ext
        case '.mat'
            S = load(fname);
            f = fieldnames(S);
            x = S.(f{1});
            v = local_id_key(x,pad5);

        otherwise
            % Some basin identifiers must be read literally. CAMELS-PE
            % IDs can contain E, while CAMELSH-US IDs commonly begin
            % with zero. READCELL may convert either form to a number,
            % changing the identifier. Read these basin files
            % literally as text lines. Keep the established loader for
            % all other regions so their behavior is unchanged.
            if local_requires_literal_basin_ids(fname)
                C = readlines(fname);
                v = local_id_key(C(:),pad5);
            else
                C = readcell(fname, ...
                    'FileType','text');
                v = local_id_key(C(:),pad5);
            end
    end

    v = v(v ~= "" & ~ismissing(v));
    
    if isempty(v)
        error(['      Error: sample_basins: ' ...
            'File %s does not ' ...
            'contain any basin IDs.'],fname);
    end
end

% -------------------------------------------------------
% helper: identify basin files requiring literal text IDs
% -------------------------------------------------------
function tf = local_requires_literal_basin_ids(fname)
%LOCAL_REQUIRES_LITERAL_BASIN_IDS Preserve significant ID characters.
%
% Peru IDs such as 47E01126 and zero-prefixed CAMELSH-US and MACH-US
% IDs must be read literally. Limit this handling to their regional
% directories or filename prefixes so behavior elsewhere is unchanged.

    f = upper(strrep(char(string(fname)),'/','\'));
    [~,base,~] = fileparts(f);

    tf = contains(f,[filesep 'PE' filesep]) ...
        || startsWith(upper(string(base)),"PE_") ...
        || contains(f,[filesep 'USH' filesep]) ...
        || startsWith(upper(string(base)),"USH_") ...
        || startsWith(upper(string(base)),"SPLIT_USH_") ...
        || contains(f,[filesep 'MACH' filesep]) ...
        || startsWith(upper(string(base)),"MACH_") ...
        || startsWith(upper(string(base)),"SPLIT_MACH_");

end

% -------------------------------------------------------
% helper: resolve basin file for GUI, scripts, and deploy
% -------------------------------------------------------
function fname = local_resolve_basin_file(fname,dirD)

    % If already a valid full/relative path, keep it
    if isfile(fname)
        return
    end
    
    cands = {};
    
    % --------------------------
    % 1. Relative to dirD itself
    % --------------------------
    if nargin >= 2 && ~isempty(dirD)
        dirD = char(dirD);
        cands{end+1} = fullfile(dirD,fname);
    
        % parent of data directory
        root1 = fileparts(dirD);
        cands{end+1} = fullfile(root1, ...
            fname);
        cands{end+1} = fullfile(root1, ...
            'basins',fname);
        cands{end+1} = fullfile(root1, ...
            'SAGEhydrology','basins',fname);
    end
    
    % ---------------------------------------
    % 2. Relative to this function's location
    % ---------------------------------------
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        thisDir = fileparts(thisFile);
        cands{end+1} = fullfile(thisDir, ...
            fname);
        cands{end+1} = fullfile(thisDir, ...
            '..','..','basins',fname);
        cands{end+1} = fullfile(thisDir, ...
            '..','..','..','basins',fname);
    end
    
    % -------------------------------
    % 3. Relative to deployed archive
    % -------------------------------
    if isdeployed
        cands{end+1} = fullfile(ctfroot, ...
            fname);
        cands{end+1} = fullfile(ctfroot, ...
            'basins',fname);
        cands{end+1} = fullfile(ctfroot, ...
            'SAGEhydrology','basins',fname);
    end
    
    % -----------------
    % 4. On MATLAB path
    % -----------------
    p = which(fname);
    if ~isempty(p)
        fname = p;
        return
    end
    
    % ---------------
    % Test candidates
    % ---------------
    for i = 1:numel(cands)
        ftry = cands{i};
        try
            ftry = char(string(ftry));
        catch
        end
        if isfile(ftry)
            fname = ftry;
            return
        end
    end
    
    % ---------------------------
    % Fail with clear diagnostics
    % ---------------------------
    msg = sprintf(['      Error: sample_basins: ' ...
        'Cannot find file: \n' ...
        '%s\n'],fname);
    
    if ~isempty(cands)
        msg = sprintf('%s      Tried:\n', ...
            msg);
        for i = 1:numel(cands)
            msg = sprintf('%s      %s\n', ...
                msg,cands{i});
        end
    end
    
    error(msg);
end

function zout = local_subset_zone(zone,idx)
%LOCAL_SUBSET_ZONE Subset hydroclimatic zone structure by basin indices

    zout = zone;
    
    flds = fieldnames(zone);
    
    for k = 1:numel(flds)
        f = flds{k};
        x = zone.(f);
    
        if isvector(x) ...
                && numel(x) == numel(zone.id)
            zout.(f) = x(idx);
        end
    end
    
    % Recompute unique names/numbers for selected basins
    if isfield(zout,'id') ...
            && ~isempty(zout.id)
        [names,~,ic] = unique(string(zout.id(:)),'stable');
        zout.names = names;
        zout.num = ic;
        if isfield(zout,'id_short')
            [namesShort,~,~] = ...
                unique(string(zout.id_short(:)),'stable');
            zout.names_short = namesShort;
        end
    end

end

function [id_t,id_e] = local_stratified_sample(zone_num,K_t,K_e)
%LOCAL_STRATIFIED_SAMPLE Stratified train/eval basin sampling by zone

    zone_num = zone_num(:);
    K_avail = numel(zone_num);
    K = K_t + K_e;
    
    if K > K_avail
        error(['      Error:local_stratified_sample: ' ...
            'Requested %d basins, ' ...
            'but only %d available.'], ...
            K,K_avail);
    end
    
    zones = unique(zone_num(isfinite(zone_num)), ...
        'stable');
    nz = numel(zones);
    
    % First choose K basins while preserving zone proportions
    id_sel = [];
    
    for j = 1:nz
        z = zones(j);
        iz = find(zone_num == z);
        nz_avail = numel(iz);
    
        kz = round(K * nz_avail / K_avail);
        kz = min(kz,nz_avail);
    
        if kz > 0
            iz = iz(randperm(nz_avail,kz));
            id_sel = [id_sel; iz(:)];          %#ok
        end
    end
    
    % Correct possible rounding mismatch
    if numel(id_sel) < K
        pool = setdiff((1:K_avail)',id_sel,'stable');
        add = pool(randperm(numel(pool),K-numel(id_sel)));
        id_sel = [id_sel; add(:)];
    elseif numel(id_sel) > K
        id_sel = id_sel(randperm(numel(id_sel),K));
    end
    
    % Now split selected basins into train/eval, again stratified by zone
    id_t = [];
    id_e = [];
    
    for j = 1:nz
        z = zones(j);
        iz = id_sel(zone_num(id_sel) == z);
        iz = iz(randperm(numel(iz)));
    
        ktz = round(K_t * numel(iz) / K);
    
        ktz = min(ktz,numel(iz));
        kez = numel(iz) - ktz;
    
        id_t = [id_t; iz(1:ktz)];              %#ok
        id_e = [id_e; iz(ktz+1:ktz+kez)];      %#ok
    end
    
    % Correct train/eval rounding mismatch
    if numel(id_t) < K_t
        need = K_t - numel(id_t);
        move = id_e(randperm(numel(id_e),need));
        id_e = setdiff(id_e,move,'stable');
        id_t = [id_t; move(:)];
    elseif numel(id_t) > K_t
        move = id_t(randperm(numel(id_t),numel(id_t)-K_t));
        id_t = setdiff(id_t,move,'stable');
        id_e = [id_e; move(:)];
    end
    
    if numel(id_e) < K_e
        need = K_e - numel(id_e);
        pool = setdiff(id_sel,[id_t; id_e],'stable');
        add = pool(randperm(numel(pool),need));
        id_e = [id_e; add(:)];
    elseif numel(id_e) > K_e
        id_e = id_e(randperm(numel(id_e),K_e));
    end
    
    % Shuffle within train/eval so zones are not block ordered
    id_t = id_t(randperm(numel(id_t)));
    id_e = id_e(randperm(numel(id_e)));

end

function [id_gauge,id_ref] = local_harmonize_de_ids(id_gauge,id_ref)
%LOCAL_HARMONIZE_DE_IDS Harmonize CAMELS-DE prefixed/unprefixed IDs.
%
% CAMELS-DE file names and basin-list files use IDs like DE110000,
% DEA10000, DEB10000, ... .  Some source tables may store the same basins
% as 110000, A10000, B10000, ... .  When the reference list is clearly
% CAMELS-DE and there is no direct overlap, prefix DE to the in-memory
% gauge IDs for matching.

    id_gauge = local_id_key(id_gauge);
    id_ref = local_id_key(id_ref);
    
    [tf0,~] = ismember(id_ref,id_gauge);
    if any(tf0)
        return
    end
    
    isRefDE = any(startsWith(id_ref,"DE"));
    isGaugeDE = any(startsWith(id_gauge,"DE"));
    
    if isRefDE ...
            && ~isGaugeDE
        id_try = "DE" + id_gauge;
        [tf1,~] = ismember(id_ref,id_try);
        if any(tf1)
            id_gauge = id_try;
            return
        end
    end
    
    if ~isRefDE ...
            && isGaugeDE
        id_ref_try = regexprep(id_ref,'^DE','');
        [tf2,~] = ismember(id_ref_try,regexprep(id_gauge,'^DE',''));
        if any(tf2)
            id_ref = "DE" + id_ref;
        end
    end

end

function id_ref = local_harmonize_de_reference(id_target,id_ref)
%LOCAL_HARMONIZE_DE_REFERENCE Align metadata IDs to selected gauge IDs.

    id_target = local_id_key(id_target);
    id_ref = local_id_key(id_ref);
    
    [tf0,~] = ismember(id_target,id_ref);
    if all(tf0)
        return
    end
    
    isTargetDE = any(startsWith(id_target,"DE"));
    isRefDE = any(startsWith(id_ref,"DE"));
    
    if isTargetDE ...
            && ~isRefDE
        id_try = "DE" + id_ref;
        [tf1,~] = ismember(id_target,id_try);
        if sum(tf1) > sum(tf0)
            id_ref = id_try;
        end
    elseif ~isTargetDE ...
            && isRefDE
        id_try = regexprep(id_ref,'^DE','');
        [tf1,~] = ismember(id_target,id_try);
        if sum(tf1) > sum(tf0)
            id_ref = id_try;
        end
    end

end

function id = local_id_key(x,pad5)

    if nargin < 2
        pad5 = false;
    end
    
    if iscell(x)
        x = string(x);
    elseif isnumeric(x)
        x = string(round(double(x(:))));
    else
        x = string(x);
    end
    
    id = upper(strtrim(x(:)));
    id = regexprep(id,'\.0+$','');
    
    if pad5
        isNum = ~isnan(str2double(id));
        id(isNum) = compose('%05.0f', ...
            str2double(id(isNum)));
    end

end

function tf = local_has_prepared_basins(bas)

    tf = isstruct(bas) ...
        && isfield(bas,'K') ...
        && isfield(bas,'id_t') ...
        && isfield(bas,'id_e') ...
        && isfield(bas,'id_gauge') ...
        && ~isempty(bas.id_gauge) ...
        && numel(bas.id_gauge) == bas.K;

end

function [A,bas,latlon] = local_use_prepared_basins( ...
    A_cns,ID,bas,prd,gname,zone,dirD)

    id_pre = local_clean_basin_id(bas.id_gauge);
    id_all = local_clean_basin_id(ID);
    
    [tf,loc] = ismember(id_pre,id_all);
    
    if ~all(tf)
        miss = id_pre(~tf);
        error(['      Error:sample_basins: prepared basin ' ...
            'selection contains basin(s) not found in ' ...
            'read_attr output, for example %s.'], ...
            char(miss(1)));
    end
    
    if ~isempty(A_cns)
        A = A_cns(:,loc);
    else
        A = [];
    end
    
    % Preserve/refresh basin metadata in prepared order
    bas.id_gauge = id_pre(:);
    bas.K = numel(id_pre);
    bas.r = numel(bas.id_attr);
    bas.id_plot = (1:bas.K).';
    try
        bas.mode = get_assess_mode(prd,bas);
    catch
        bas.mode = 'none';
    end

    if isfield(bas,'id_t') && isfield(bas,'id_e')
        bas.K_t = numel(bas.id_t);
        bas.K_e = numel(bas.id_e);
    else
        bas.K_t = bas.K;
        bas.K_e = 0;
    end
    
    try
        gname_all = string(gname(:));
        bas.gname = gname_all(loc);
    catch
    end
    
    try
        bas.zone = local_subset_zone(zone,loc);
    catch
    end
    
    try
        latlon = get_latlon(dirD,bas);
    catch
        latlon = [];
    end

end

function id = local_clean_basin_id(id)

    id = string(id(:));
    id = strtrim(id);
    id = erase(id,'"');
    id = erase(id,"'");
    id = regexprep(id,'\.0+$','');

end

function id = local_strip_de_prefix(id)

    id = string(id(:));
    id = regexprep(id,'^DE','');

end

function id_out = local_normalize_split_ids(id_in, id_univ)

    id_in   = string(id_in(:));
    id_univ = string(id_univ(:));
    
    % -----------------------------
    % Case 1: already valid IDs
    % -----------------------------
    [tf,~] = ismember(id_in, id_univ);
    if all(tf)
        id_out = id_in;
        return
    end
    
    % -----------------------------
    % Case 2: numeric / index file
    % -----------------------------
    id_num = str2double(id_in);

    if all(~isnan(id_num))

        % ===================================
        % CAMELS-STYLE DECODE (BLOCK + BASIN)
        % ===================================

        id_num = round(id_num);

        % split into components
        block = floor(id_num / 1000);
        basin = mod(id_num, 1000);

        id_out = strings(size(id_in));

        for i = 1:numel(id_in)

            if isnan(block(i)) ...
                    || isnan(basin(i))
                continue
            end

            % reconstruct pattern search key
            key_block = string(block(i));

            % find matching IDs in universe
            cand = id_univ(contains(id_univ, ...
                key_block));

            if isempty(cand)
                error('sample_basins:splitMismatch', ...
                    'Block %s not found in universe.', ...
                    key_block);
            end

            % pick basin index within block group
            if basin(i) <= numel(cand)
                id_out(i) = cand(basin(i));
            else
                error('sample_basins:splitMismatch', ...
                    ['Invalid basin index in ' ...
                    'split file: %d'], id_num(i));
            end
        end

        return
    end
    % ------------------------------
    % Case 3: partial match fallback
    % ------------------------------
    % try stripping leading zeros
    id_clean = regexprep(id_in,'^0+','');
    id_univ_clean = regexprep(id_univ,'^0+','');
    
    [tf,loc] = ismember(id_clean, ...
        id_univ_clean);
    
    if all(tf)
        id_out = id_univ(loc);
        return
    end
    
    error('sample_basins:splitMismatch', ...
        ['Split file IDs do not ' ...
        'match universe IDs.']);

end
