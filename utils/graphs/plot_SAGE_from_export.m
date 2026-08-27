function zoneTbl = plot_SAGE_from_export(matFile,varargin)
%PLOT_SAGE_FROM_EXPORT Recreate SAGE postprocessing figures from an export 
% file.
%
% This utility reconstructs the data structures required by PLOT_SAGE from
% the final SAGE export structure E and regenerates the standard suite of
% postprocessing figures without rerunning model simulations or training.
%
% The function supports replay of ECDF diagnostics, hydroclimatic-zone
% summaries, parameter maps, parameter traces, variograms, flow-duration
% curves, and basin-scale hydrograph diagnostics directly from exported
% SAGE results.
%
% USAGE
% zoneTbl = plot_SAGE_from_export('SAGE_export_....mat');
% zoneTbl = plot_SAGE_from_export(E);
%
% OPTIONAL NAME-VALUE PAIRS
% 'LastOutFile' '' Optional SAGEnew_lastOut.mat file containing
% auxiliary postprocessing information such as
% basin coordinates, parameter traces, parameter
% maps, flow-duration curves, and GUI metadata.
%
% 'Scenarios' {'tt','te','et','ee'}
% Benchmark scenarios to include in the replay.
%
% 'PlotFDC' true Plot flow-duration-curve diagnostics when
% available.
%
% 'Gaugescen' [] Override the default basin selections used for
% basin-scale diagnostics and hydrograph plots.
%
% INPUT
% matFile Path to a SAGE export MAT-file or an export structure E.
%
% OUTPUT
% zoneTbl Hydroclimatic-zone performance summary table returned by
% PLOT_SAGE.
%
% NOTES
% This function does not rerun SAGE and does not perform any additional
% optimization. All figures are generated from information stored in the
% export file and, optionally, the accompanying SAGEnew_lastOut.mat file.
%
% EXAMPLE
% zoneTbl = plot_SAGE_from_export( ...
%   'SAGE_export_CAMELS_DE_4_hourly_20260616_081730.mat', ...
%   'LastOutFile','SAGEnew_lastOut.mat');
%
% See also: PLOT_SAGE

    p = inputParser;
    p.addRequired('matFile');
    p.addParameter('Scenarios',{'tt','te','et','ee'}, ...
        @(x)iscellstr(x) || isstring(x));
    p.addParameter('PlotFDC',true,@(x)islogical(x) ...
        && isscalar(x));
    p.addParameter('Gaugescen',[],@(x)isstruct(x) ...
        || isempty(x));
    p.addParameter('LastOutFile','', ...
        @(x)ischar(x) ...
        || isstring(x));
    p.parse(matFile,varargin{:});

    scen = cellstr(string(p.Results.Scenarios));

    % -----------------
    % Load final export
    % -----------------
    if isstruct(matFile)
        E = matFile;
    else
        S = load(matFile,'E');
        if ~isfield(S,'E')
            error(['The MAT file must ' ...
                'contain a structure named E.']);
        end
        E = S.E;
    end
    out = struct();
    if ~isempty(p.Results.LastOutFile)
        S2 = load(p.Results.LastOutFile,'out');
        if isfield(S2,'out')
            out = S2.out;
        end
    end
    
    % -------------------------
    % Rebuild required objects
    % -------------------------
    mdl = get_first_existing(E,{'mdl', ...
        'model','codes.mdl','metadata.mdl'});
    bas = get_first_existing(E,{'bas', ...
        'basins_struct','metadata.bas','metadata.basins'});
    prd = get_first_existing(E,{'prd', ...
        'period','metadata.prd','metadata.period'});

    if isempty(mdl)
        mdl = struct(); 
    end
    if isempty(bas)
        bas = struct(); 
    end
    if isempty(prd)
        prd = struct(); 
    end

    % Prefer richer information from SAGEnew_lastOut.mat

    if isfield(out,'bas') ...
            && ~isempty(out.bas)
        bas = out.bas;
    end
    if isfield(out,'prd') ...
            && ~isempty(out.prd)
        prd = out.prd;
    end

    % Fill basin structure from E.basins table when possible.
    if isfield(E,'basins') ...
            && istable(E.basins)
        B = sort_by_row_id(E.basins);
        K = height(B);
        if ~isfield(bas,'K')
            bas.K = K; 
        end
        if ~isfield(bas,'K_t')
            bas.K_t = nrows_scenario(E,'tt'); 
        end
        if ~isfield(bas,'K_e')
            bas.K_e = max(0,K - bas.K_t); 
        end
        if ~isfield(bas,'gname')
            bas.gname = table_column_or_default(B, ...
            {'name','gname','gauge_name','basin_name'}, ...
            repmat({''},K,1));
        end
        if ~isfield(bas,'id_gauge')
            bas.id_gauge = table_column_or_default(B, ...
            {'gauge','id_gauge','gauge_id', ...
            'basin_id','camels_id'},(1:K)');
        end
        if ~isfield(bas,'zone')
            bas.zone = build_zone_from_table(B,K);
        end
        latlon = [];
        if isfield(out,'latlon') ...
                && ~isempty(out.latlon)
            latlon = out.latlon;
        elseif isfield(E,'basins') ...
                && istable(E.basins)
            latlon = build_latlon_from_table(B);
        end
    else
        latlon = [];
    end

    % Make model bookkeeping robust enough for plot_SAGE.
    mdl = complete_mdl_defaults(mdl,bas,prd);
    prd = complete_prd_defaults(prd);

    % Metrics: plot_SAGE wants structures NSE, KGE, JKGE with fields tt/te/et/ee.
    [NSE,KGE,S_fdc,JKGE] = performance_to_metric_structs(E,scen);
    curr = struct('NSE',NSE,'KGE',KGE, ...
        'S_fdc',S_fdc,'JKGE',JKGE);

    % Discharge: plot_SAGE wants Q plus dat{k}.y_n/bad/gauge for time series.
    [Q,dat] = discharge_to_Q_dat(E,bas);
    
    % Rainfall-block exports need local train/evaluation masks in dat{k}.
    % Reconstruct them from E.discharge.split_code when available:
    %   1 = training time step
    %   2 = evaluation time step
    [dat,mdl] = complete_rainfall_block_masks(E,dat,mdl);

    % % Rainfall-block exports need local train/evaluation masks in dat{k}.
    % % If those masks are not present in the export, fall back to a
    % % non-local postprocessing label so ECDF/FDC figures can still be made.
    % [dat,mdl] = complete_rainfall_block_masks(E,dat,mdl,bas);

    % Optional final products.
    Qfdc = [];
    if p.Results.PlotFDC
        if isfield(out,'Qfdc') ...
                && ~isempty(out.Qfdc)
            Qfdc = out.Qfdc;
        else
            Qfdc = get_first_existing(E, ...
                {'Qfdc','discharge.Qfdc', ...
                'discharge.fdc'});
        end
    end
    region = get_first_existing(E,{'region', ...
        'metadata.region','codes.region'});
    if isempty(region)
        region = 'US'; 
    end
    if isfield(out,'tTheta') ...
            && ~isempty(out.tTheta)
        tTheta = out.tTheta;
    else
        tTheta = get_first_existing(E, ...
            {'tTheta','theta_trace', ...
            'metadata.tTheta', ...
            'performance.tTheta'});
    end
    if isfield(out,'nTheta_last') ...
            && ~isempty(out.nTheta_last)    
        nTheta = out.nTheta_last;    
    else    
        nTheta = theta_to_nTheta(E,bas);   
    end
    At = get_first_existing(E,{'At', ...
        'attribution.At','metadata.At'});
    An = get_first_existing(E,{'An', ...
        'attribution.An','metadata.An'});

    gaugescen = p.Results.Gaugescen;
    if isempty(gaugescen)
        gaugescen = get_first_existing(E, ...
        {'gaugescen','metadata.gaugescen', ...
        'discharge.gaugescen'});
    end
    if isempty(gaugescen)
        gaugescen = struct(); 
    end

    % -------------------------
    % Create figures
    % -------------------------
    zoneTbl = plot_SAGE('sage',mdl,dat,bas, ...
        prd,Q,curr,Qfdc,region,tTheta,latlon, ...
        nTheta,At,An,gaugescen);
end

% ================
% Helper functions
% ================
function v = get_first_existing(S,paths)
    v = [];
    for i = 1:numel(paths)
        parts = strsplit(paths{i},'.');
        x = S; ok = true;
        for j = 1:numel(parts)
            if isstruct(x) ...
                && isfield(x,parts{j})
                x = x.(parts{j});
            else
                ok = false; 
                break
            end
        end
        if ok ...
                && ~isempty(x)
            v = x; return
        end
    end
end

function T = sort_by_row_id(T)
    if istable(T) && any(strcmp( ...
        T.Properties.VariableNames,'row_id'))
        T = sortrows(T,'row_id');
    end
end

function n = nrows_scenario(E,field)
    n = 0;
    if isfield(E,'performance') ...
            && isfield(E.performance,field) ...
            && istable(E.performance.(field))
        n = height(E.performance.(field));
    end
end

function c = table_column_or_default(T,names,defaultVal)
    c = defaultVal;
    for i = 1:numel(names)
        hit = strcmpi(T.Properties.VariableNames,names{i});
        if any(hit)
            c = T{:,find(hit,1)};
            return
        end
    end
end

function zone = build_zone_from_table(T,K)
    zone = struct();
    zid = table_column_or_default(T,{'zone', ...
        'zone_id','hydroclimatic_zone'}, ...
        repmat("all",K,1));
    zid = string(zid);
    [names,~,num] = unique(zid,'stable');
    zone.id = zid;
    zone.num = num(:);
    zone.names = names(:);
    zone.aridity = table_column_or_default(T, ...
        {'aridity','aridity_index'},nan(K,1));
    zone.frac_snow = table_column_or_default(T, ...
        {'frac_snow','snow_fraction', ...
        'snow_frac'},nan(K,1));
end

function latlon = build_latlon_from_table(T)
    lat = table_column_or_default(T,{'lat', ...
        'latitude','gauge_lat'},[]);
    lon = table_column_or_default(T,{'lon', ...
        'long','longitude','gauge_lon'},[]);
    if isempty(lat) ...
            || isempty(lon)
        latlon = [];
    else
        latlon = [double(lat(:)) double(lon(:))];
    end
end

function mdl = complete_mdl_defaults(mdl,bas,prd)
    if ~isfield(mdl,'names') ...
            || isempty(mdl.names)
        mdl.names = {'hymod','hmodel','sacsma', ...
            'xinanjiang','gr4jA','hbv', ...
            'gr4jB','cfe_nwm'};
    end
    if ~isfield(mdl,'model') ...
            || isempty(mdl.model)
        mdl.model = infer_model_index(mdl);
    end
    if ~isfield(mdl,'mcode')
        mdl.mcode = 4; end
    if ~isfield(mdl,'calc')
        mdl.calc = 'seq'; end
    if ~isfield(mdl,'mode')        
        mdl.mode = infer_mode_from_bas(bas); 
    end
    if ~isfield(mdl,'sp_method')
        mdl.sp_method = 'manual'; 
    end
    if ~isfield(mdl,'pspace')
        mdl.pspace = 1; 
    end
    if ~isfield(mdl,'id_train') ...
            || isempty(mdl.id_train)
        if isfield(prd,'id_train')
            mdl.id_train = prd.id_train;
        else
            mdl.id_train = 1; 
        end
    end
    if ~isfield(mdl,'id_eval')
        if isfield(prd,'id_eval')
            mdl.id_eval = prd.id_eval;
        else
            mdl.id_eval = []; 
        end
    end
end

function idx = infer_model_index(mdl)
    idx = 1;
    if isfield(mdl,'model_name')
        hit = find(strcmpi(mdl.names, ...
        mdl.model_name),1);
        if ~isempty(hit)
            idx = hit; 
        end
    end
end

function mode = infer_mode_from_bas(bas)
    if isfield(bas,'K_e') ...
            && bas.K_e > 0
        mode = 4;
    else
        mode = 2;
    end
end

function prd = complete_prd_defaults(prd)
    if ~isfield(prd,'dt') ...
            || isempty(prd.dt) 
        prd.dt = 1; 
    end
    if ~isfield(prd,'spinup') ...
            || isempty(prd.spinup)
        prd.spinup = 0; 
    end
end

function [NSE,KGE,S_fdc,JKGE] = performance_to_metric_structs(E,scen)
    NSE = struct(); 
    KGE = struct(); 
    S_fdc = struct();
    JKGE = struct();
    if ~isfield(E,'performance')
        error(['E.performance is required ' ...
            'for ECDF/CDF figures.']);
    end
    for i = 1:numel(scen)
        s = scen{i};
        if ~isfield(E.performance,s) ...
                || ~istable(E.performance.(s))
            continue; 
        end
        T = sort_by_row_id(E.performance.(s));
        NSE.(s) = col_as_vector(T,'NSE');
        KGE.(s) = col_as_vector(T,'KGE');
        S_fdc.(s) = col_as_vector(T,'S_fdc');
        JKGE.(s) = col_as_vector(T,'JKGE');
    end
end

function x = col_as_vector(T,name)
    hit = strcmpi(T.Properties.VariableNames,name);
    if any(hit)
        x = T{:,find(hit,1)};
        x = double(x(:));
    else
        x = nan(height(T),1);
    end
end

function [Q,dat] = discharge_to_Q_dat(E,bas)

    Q = struct();

    K = bas.K;
    dat = cell(K,1);

    hasObs = isfield(E,'discharge') ...
            && isfield(E.discharge,'Q_measured') ...
            && ~isempty(E.discharge.Q_measured);

    if ~hasObs
        error(['E.discharge.Q_measured is ' ...
            'required for hydrograph replay.']);
    end

    Y = double(E.discharge.Q_measured);

    if size(Y,1) == K
        Y = Y.';
    end

    for k = 1:K
        dat{k} = struct();

        try
            dat{k}.gauge = bas.id_gauge(k);
        catch
            dat{k}.gauge = k;
        end

        dat{k}.y_n = Y(:,k);
        dat{k}.bad = ~isfinite(dat{k}.y_n);
    end

    if isfield(E.discharge,'Q_scenario') ...
            && isstruct(E.discharge.Q_scenario)

        Q = E.discharge.Q_scenario;

    else

        hasSim = isfield(E.discharge,'Q_simulated') ...
                && ~isempty(E.discharge.Q_simulated);

        if ~hasSim
            error(['E.discharge.Q_scenario or ' ...
                'E.discharge.Q_simulated is required.']);
        end

        Qs = double(E.discharge.Q_simulated);

        if size(Qs,1) == K
            Qs = Qs.';
        end

        Q.tt = Qs(:,1:bas.K_t);
        Q.te = Qs(:,1:bas.K_t);

        if bas.K_e > 0
            Q.et = Qs(:,bas.K_t+1:bas.K_t+bas.K_e);
            Q.ee = Qs(:,bas.K_t+1:bas.K_t+bas.K_e);
        else
            Q.et = [];
            Q.ee = [];
        end
    end
end

function [dat,mdl] = complete_rainfall_block_masks(E,dat,mdl)

    if ~isfield(mdl,'sp_method') ...
            || ~strcmpi(string(mdl.sp_method), ...
        'rainfall_block')
        return
    end

    if isfield(E,'discharge') ...
            && isfield(E.discharge,'split_code') ...
            && ~isempty(E.discharge.split_code)

        split_code = E.discharge.split_code;

        if size(split_code,1) ~= numel(dat)
            error(['E.discharge.split_code has ' ...
                    '%d rows, but dat has %d basins. ' ...
                    'These must match to reconstruct ' ...
                    'rainfall-block masks.'], ...
                   size(split_code,1),numel(dat));
        end

        for k = 1:numel(dat)
            dat{k}.id_train = find(split_code(k,:) == 1);
            dat{k}.id_eval  = find(split_code(k,:) == 2);
        end

        mdl.sp_method = 'rainfall_block';
        mdl.local = 1;
        return
    end

    warning(['The export does not contain ' ...
             'E.discharge.split_code, so local ' ...
             'rainfall-block train/evaluation ' ...
             'masks cannot be reconstructed. ' ...
             'Switching mdl.sp_method to manual ' ...
             'for limited postprocessing replay.']);

    mdl.sp_method = 'manual';
    mdl.local = 0;

end

function nTheta = theta_to_nTheta(E,bas)
    nTheta = get_first_existing(E, ...
        {'nTheta','theta_n_array','metadata.nTheta'});
    if ~isempty(nTheta), return; end
    nTheta = [];
    if isfield(E,'theta_n') && istable(E.theta_n)
        T = sort_by_row_id(E.theta_n);
        vars = T.Properties.VariableNames;
        keep = ~strcmpi(vars,'row_id');
        X = T{:,keep};       % K x d
        nTheta = double(X');    % d x K, as expected by plot_SAGE
        if nargin > 1 && isfield(bas,'K') && size(nTheta,2) ~= bas.K
            warning(['theta_n rows do not match ' ...
                    'bas.K; parameter maps may be skipped.']);
            nTheta = [];
        end
    end
end
