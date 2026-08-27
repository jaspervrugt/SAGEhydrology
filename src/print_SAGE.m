function [ax,frmt] = print_SAGE(mdl,ax,prf,i,dirres,loss,net)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRINT_SAGE Prints iteration statistics and updates live diagnostic figures
%
% SYNOPSIS:
%  [ax,frmt] = print_SAGE(mdl,ax,prf,i,dirres,loss,net)
%
%   mdl         structure with model state/parameter info
%    .mode       assessment design
%                 1 = train basins | train period/mask
%                 2 = train basins | train and eval period/mask
%                 3 = train and eval basins | train period/mask
%                 4 = train and eval basins | train and eval period/mask
%   ax          optional struct with axes handles; pass [] on first call
%   prf         structure with current basin metrics in .curr and
%               iteration histories in .iter
%   i           iteration number
%   dirres      directory with SAGE results
%   loss        loss-function settings; loss.fnc selects the objective
%   net         network settings; net.l is the number of FFN parameters
%
% OUTPUT:
%   ax          updated graphics-handle structure
%   frmt        print format string used by print_stats
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6 ...
            || isempty(loss)
        loss_fnc = 4;
    else
        loss_fnc = loss.fnc;
    end
    if nargin < 7 ...
            || isempty(net) ...
            || ~isfield(net,'l') ...
            || isempty(net.l)
        l = nan;
    else
        l = net.l;
    end
    
    if ~isfield(mdl,'mode') ...
            || isempty(mdl.mode) ...
            || ~isscalar(mdl.mode) ...
            || ~isnumeric(mdl.mode) ...
            || ~ismember(double(mdl.mode),[1 2 3 4])
        error(['      Error: print_SAGE: ' ...
            'mdl.mode must be one of 1,2,3,4.']);
    end
    
    frmt = print_stats(mdl.mode,prf,i,l,loss_fnc,dirres);
    ax = print_figs(mdl,prf,i,prf.curr.NSE,prf.curr.KGE, ...
        prf.curr.S_fdc,prf.curr.JKGE,loss_fnc,ax);

end

function frmt = print_stats(mode,prf,i,l,loss_fnc,dirres)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRINT_STATS Prints iteration statistics to screen and log file
%
% SYNOPSIS:
%  frmt = print_stats(mode,prf,i,l,loss_fnc,dirres)
%
%   mode        assessment design
%                1 = training basins only | training period only
%                2 = training basins only | training and evaluation
%                    period/mask
%                3 = training and evaluation basins | training period only
%                4 = training and evaluation basins | training and
%                    evaluation period/mask
%   prf         structure with performance histories
%   i           iteration number
%   l           number of network weights and biases
%   loss_fnc    loss function (scalar, optional if stored in prf)
%   dirres      directory with SAGE results
%
% OUTPUT:
%   frmt        print format string retained for backward compatibility
%
% NOTES:
%   - Statistics are printed in compact scenario-row table format.
%   - Depending on mode, scenarios tt, te, et, and/or ee are reported.
%   - Reported quantities include selected loss, RSS, median NSE,
%     median KGE or JKGE, and S_IB.
%   - The header is reprinted periodically and RSS is scaled adaptively.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Apr. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 6 ...
            || isempty(dirres)
        dirres = pwd;
    end
    if ~exist(dirres,'dir')
        mkdir(dirres);
    end
    
    if nargin < 5 ...
            || isempty(loss_fnc)
        if isfield(prf,'loss_fnc') ...
                && ~isempty(prf.loss_fnc)
            loss_fnc = prf.loss_fnc;
        else
            error(['      Error: print_stats: ' ...
                'loss_fnc must be provided, ' ...
                   'or stored in prf.loss_fnc.']);
        end
    end
    
    log_file = fullfile(dirres, ...
        'demo_SAGE_iterations.txt');
    
    if i == 1
        fid = fopen(log_file,'w');
    else
        fid = fopen(log_file,'a');
    end
    
    if fid == -1
        error(['      Error: print_stats: ' ...
            'Could not open log file.']);
    end
    
    % ------------
    % User options
    % ------------
    headerEvery = 25;   % Reprint header every N iterations
    rssDigits = 3;      % Digits after decimal for scaled RSS
    lossDigits = 3;
    metDigits = 3;
    
    % -------------------------------------------------
    % Selected loss function label and associated field
    % -------------------------------------------------
    switch loss_fnc
        case 1
            lossName = 'ΣSAR';
            lossField = 'SAR';
            useRSS = true;
    
        case 2
            lossName = 'ΣGLS';
            lossField = 'GLS';
            useRSS = false;   % GLS already is weighted RSS-like quantity
    
        case 3
            lossName = 'Σ1-NSE';
            lossField = 'L';
            useRSS = true;
    
        case 4
            lossName = 'Σ1-KGE';
            lossField = 'L';
            useRSS = true;
    
        case 5
            lossName = 'ΣHuber';
            lossField = 'Huber';
            useRSS = true;
    
        case 6
            lossName = 'Σd_FDC';
            lossField = 'L';
            useRSS = true;
    
        case 7
            lossName = 'Σ1-JKGE';
            lossField = 'L';
            useRSS = true;
    
        otherwise
            fclose(fid);
            error(['      Error: print_stats: ' ...
                'unknown loss_fnc = %g.'],loss_fnc);
    end
    
    K_t = prf.K_t;
    K_e = prf.K_e;
    
    % ------------------------------------------
    % Determine scenarios to print based on mode
    % ------------------------------------------
    switch mode
        case 1
            scens = {'tt'};
        case 2
            scens = {'tt','te'};
        case 3
            scens = {'tt','et'};
        case 4
            scens = {'tt','te','et','ee'};
        otherwise
            fclose(fid);
            error(['      Error: print_stats: ' ...
                'mode must be one of 1,2,3,4.']);
    end
    
    nScen = numel(scens);
    
    % --------------------------------
    % Collect values for all scenarios
    % --------------------------------
    cpuT = prf.iter.cpuT(i);
    
    Loss = nan(1,nScen);
    RSS = nan(1,nScen);
    mNSE = nan(1,nScen);
    mKGE = nan(1,nScen);
    mS_fdc = nan(1,nScen);
    mJKGE = nan(1,nScen);
    SibNSE = nan(1,nScen);
    SibKGE = nan(1,nScen);
    SibFDC = nan(1,nScen);
    
    for j = 1:nScen
        scn = scens{j};
    
        switch scn
            case {'tt','te'}
                Kscale = K_t;
            case {'et','ee'}
                Kscale = K_e;
            otherwise
                Kscale = 1;
        end
    
        if isfield(prf.iter.(lossField),scn) ...
                && i <= numel(prf.iter.(lossField).(scn))
            Loss(j) = Kscale * prf.iter.(lossField).(scn)(i);
        end
    
        if isfield(prf.iter.RSS,scn) ...
                && i <= numel(prf.iter.RSS.(scn))
            RSS(j) = Kscale * prf.iter.RSS.(scn)(i);
        end
    
        if isfield(prf.iter.mNSE,scn) ...
                && i <= numel(prf.iter.mNSE.(scn))
            mNSE(j) = prf.iter.mNSE.(scn)(i);
        end
    
        if isfield(prf.iter,'mKGE') ...
                && isfield(prf.iter.mKGE,scn) ...
                && i <= numel(prf.iter.mKGE.(scn))
            mKGE(j) = prf.iter.mKGE.(scn)(i);
        end

        if isfield(prf.iter,'mS_fdc') ...
                && isfield(prf.iter.mS_fdc,scn) ...
                && i <= numel(prf.iter.mS_fdc.(scn))
            mS_fdc(j) = prf.iter.mS_fdc.(scn)(i);
        end
        
        if isfield(prf.iter,'mJKGE') ...
                && isfield(prf.iter.mJKGE,scn) ...
                && i <= numel(prf.iter.mJKGE.(scn))
            mJKGE(j) = prf.iter.mJKGE.(scn)(i);
        end
    
        if isfield(prf.iter,'Sib_NSE') ...
                && isfield(prf.iter.Sib_NSE,scn) ...
                && i <= numel(prf.iter.Sib_NSE.(scn))
            SibNSE(j) = prf.iter.Sib_NSE.(scn)(i);
        end

        if isfield(prf.iter,'Sib_KGE') ...
                && isfield(prf.iter.Sib_KGE,scn) ...
                && i <= numel(prf.iter.Sib_KGE.(scn))
            SibKGE(j) = prf.iter.Sib_KGE.(scn)(i);
        end

        if isfield(prf.iter,'Sib_S_fdc') ...
                && isfield(prf.iter.Sib_S_fdc,scn) ...
                && i <= numel(prf.iter.Sib_S_fdc.(scn))
            SibFDC(j) = prf.iter.Sib_S_fdc.(scn)(i);
        end
    end
    
    useJKGE = any(isfinite(mJKGE));

    % Use scientific notation for the complete S_IB column when values
    % are large or span more than four orders of magnitude. Keeping one
    % format per iteration preserves vertical alignment across scenarios.
    allSib = [SibNSE SibKGE SibFDC];
    sibVals = abs(allSib(isfinite(allSib) & allSib ~= 0));
    useSciSib = ~isempty(sibVals) && (max(sibVals) >= 1e4 ...
        || max(sibVals) / min(sibVals) >= 1e4);
    if useSciSib
        sibType = 'e';
        sibDigits = 3;
    else
        sibType = 'f';
        sibDigits = metDigits;
    end
    
    % ---------------------------------
    % Adaptive RSS scale
    % Recompute when header is printed.
    % ---------------------------------
    printHeader = (i == 1) ...
        || (mod(i-1,headerEvery) == 0);
    
    persistent rssPow_cached rssScale_cached
    if isempty(rssPow_cached)
        rssPow_cached = 0;
        rssScale_cached = 1;
    end
    
    if printHeader
        rssVals = RSS(isfinite(RSS) & RSS ~= 0);
    
        if isempty(rssVals)
            rssPow_cached = 0;
            rssScale_cached = 1;
        else
            [rssPow_cached,rssScale_cached] = ...
                nice_eng_scale(rssVals);
        end
    end
    
    RSSs = RSS ./ rssScale_cached;
    
    persistent lossPow_cached lossScale_cached
    if isempty(lossPow_cached)
        lossPow_cached = 0;
        lossScale_cached = 1;
    end
    
    if printHeader
        lossVals = Loss(isfinite(Loss) & Loss ~= 0);
    
        if isempty(lossVals)
            lossPow_cached = 0;
            lossScale_cached = 1;
        else
            [lossPow_cached,lossScale_cached] = ...
                nice_eng_scale(lossVals);
        end
    end
    
    Losss = Loss ./ lossScale_cached;
    
    % -----------
    % Loss header
    % -----------
    if lossPow_cached ~= 0
        lossHdr = sprintf('%s (x10^%d)', ...
            lossName,lossPow_cached);
    else
        lossHdr = lossName;
    end
    
    % -------------
    % Column widths
    % -------------
    wScen = 6;
    wLoss = max(12,numel(lossHdr) + 2);
    wRSS = 12;
    wNSE = 8;
    wKGE = 8;
    wFDC = 8;
    wJKGE = 8;
    wSIB = 9;
    
    % ------------------
    % Print header block
    % ------------------
    if printHeader
        if useRSS
            if useJKGE
                hdr1 = sprintf(['%-', ...
                    num2str(wScen),'s  %', ...
                    num2str(wLoss),'s  %', ...
                    num2str(wRSS),'s  %', ...
                    num2str(wNSE),'s  %', ...
                    num2str(wKGE),'s  %', ...
                    num2str(wFDC),'s  %', ...
                    num2str(wJKGE),'s  %', ...
                    num2str(wSIB),'s  %', ...
                    num2str(wSIB),'s  %', ...
                    num2str(wSIB),'s\n'], ...
                    'scen',lossHdr, ...
                    ['RSS (x10^', ...
                    num2str(rssPow_cached),')'], ...
                    'T_NSE','T_KGE','T_S_fdc','T_JKGE', ...
                    'Sib_NSE','Sib_KGE','Sib_Sfdc');
            else
                hdr1 = sprintf(['%-', ...
                    num2str(wScen),'s  %', ...
                    num2str(wLoss),'s  %', ...
                    num2str(wRSS),'s  %', ...
                    num2str(wNSE),'s  %', ...
                    num2str(wKGE),'s  %', ...
                    num2str(wFDC),'s  %', ...
                    num2str(wSIB),'s  %', ...
                    num2str(wSIB),'s  %', ...
                    num2str(wSIB),'s\n'], ...
                    'scen',lossHdr, ...
                    ['RSS (x10^', ...
                    num2str(rssPow_cached),')'], ...
                    'T_NSE','T_KGE','T_S_fdc', ...
                    'Sib_NSE','Sib_KGE','Sib_Sfdc');
            end
        else
            if useJKGE
                hdr1 = sprintf(['%-', ...
                    num2str(wScen),'s  %', ...
                    num2str(wLoss),'s  %', ...
                    num2str(wNSE),'s  %', ...
                    num2str(wKGE),'s  %', ...
                    num2str(wFDC),'s  %', ...
                    num2str(wJKGE),'s  %', ...
                    num2str(wSIB),'s  %', ...
                    num2str(wSIB),'s  %', ...
                    num2str(wSIB),'s\n'], ...
                    'scen',lossHdr, ...
                    'T_NSE','T_KGE','T_S_fdc','T_JKGE', ...
                    'Sib_NSE','Sib_KGE','Sib_Sfdc');
            else
                hdr1 = sprintf(['%-', ...
                    num2str(wScen),'s  %', ...
                    num2str(wLoss),'s  %', ...
                    num2str(wNSE),'s  %', ...
                    num2str(wKGE),'s  %', ...
                    num2str(wFDC),'s  %', ...
                    num2str(wSIB),'s  %', ...
                    num2str(wSIB),'s  %', ...
                    num2str(wSIB),'s\n'], ...
                    'scen',lossHdr, ...
                    'T_NSE','T_KGE','T_S_fdc', ...
                    'Sib_NSE','Sib_KGE','Sib_Sfdc');
            end
        end

        % Match the separator to the exact rendered header width. hdr1
        % contains one trailing newline, which is excluded from the count.
        headerWidth = numel(hdr1) - 1;
        sep = repmat('-',1,headerWidth);
        hdr0 = sprintf(['\n%s\nit %4d   ' ...
            'l %5d   CPU %5.1f s\n'],sep,i,l,cpuT);
    
        fprintf('%s',hdr0);
        fprintf('%s',hdr1);
        fprintf('%s\n',sep);
    
        fprintf(fid,'%s',hdr0);
        fprintf(fid,'%s',hdr1);
        fprintf(fid,'%s\n',sep);
    else
        % Print compact iteration leader
        lead = sprintf(['\nit %4d   l %5d   ' ...
            'CPU %5.1f s\n'],i,l,cpuT);
        fprintf('%s',lead);
        fprintf(fid,'%s',lead);
    end
    
    % -------------------
    % Print scenario rows
    % -------------------
    for j = 1:nScen
        scn = scens{j};
    
        if useRSS
            if useJKGE
                row = sprintf(['%-', ...
                    num2str(wScen),'s  %', ...
                    num2str(wLoss),'.', ...
                    num2str(lossDigits),'f  %', ...
                    num2str(wRSS),'.', ...
                    num2str(rssDigits),'f  %', ...
                    num2str(wNSE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wKGE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wFDC),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wJKGE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'\n'], ...
                    scn,Losss(j),RSSs(j), ...
                    mNSE(j),mKGE(j),mS_fdc(j),mJKGE(j), ...
                    SibNSE(j),SibKGE(j),SibFDC(j));
            else
                row = sprintf(['%-', ...
                    num2str(wScen),'s  %', ...
                    num2str(wLoss),'.', ...
                    num2str(lossDigits),'f  %', ...
                    num2str(wRSS),'.', ...
                    num2str(rssDigits),'f  %', ...
                    num2str(wNSE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wKGE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wFDC),'.',num2str(metDigits),'f  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'\n'], ...
                    scn,Losss(j),RSSs(j), ...
                    mNSE(j),mKGE(j),mS_fdc(j), ...
                    SibNSE(j),SibKGE(j),SibFDC(j));
            end
        else
            if useJKGE
                row = sprintf(['%-', ...
                    num2str(wScen),'s  %', ...
                    num2str(wLoss),'.', ...
                    num2str(lossDigits),'f  %', ...
                    num2str(wNSE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wKGE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wFDC),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wJKGE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'\n'], ...
                    scn,Losss(j),mNSE(j), ...
                    mKGE(j),mS_fdc(j),mJKGE(j), ...
                    SibNSE(j),SibKGE(j),SibFDC(j));
            else
                row = sprintf(['%-', ...
                    num2str(wScen),'s  %', ...
                    num2str(wLoss),'.', ...
                    num2str(lossDigits),'f  %', ...
                    num2str(wNSE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wKGE),'.', ...
                    num2str(metDigits),'f  %', ...
                    num2str(wFDC),'.',num2str(metDigits),'f  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'  %', ...
                    num2str(wSIB),'.',num2str(sibDigits),sibType,'\n'], ...
                    scn,Losss(j),mNSE(j), ...
                    mKGE(j),mS_fdc(j), ...
                    SibNSE(j),SibKGE(j),SibFDC(j));
            end
        end
    
        fprintf('%s',row);
        fprintf(fid,'%s',row);
    end
    
    frmt = '';   % retained only for backward compatibility
    fclose(fid);
end

function ax = print_figs(mdl,prf,i,NSE,KGE,S_fdc,JKGE,loss_fnc,ax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRINT_FIGS Live ECDFs of NSE, KGE and JKGE for training scenario and the
% mode-dependent comparison scenario, plus RSS and S_IB history plots.
%
% SYNOPSIS:
%   ax = print_figs(mdl,prf,i,NSE,KGE,S_fdc,JKGE,loss_fnc,ax)
%
%   mdl         structure with model state/parameter and assessment info
%    .mode       assessment design
%                 1 = train basins | train period/mask
%                 2 = train basins | train and eval period/mask
%                 3 = train and eval basins | train period/mask
%                 4 = train and eval basins | train and eval period/mask
%    .names      list of model names
%    .model      selected model index
%    .sp_method  split design used to label period/mask in titles
%   prf         structure with performance histories for explicit scenarios
%                tt = training basins | training period
%                te = training basins | evaluation period/mask
%                et = evaluation basins | training period
%                ee = evaluation basins | evaluation period/mask
%   i           iteration number
%   NSE         structure with scenario-wise NSE vectors
%                .tt, .te, .et, .ee
%   KGE         structure with scenario-wise KGE vectors
%                .tt, .te, .et, .ee
%   S_fdc       structure with scenario-wise S_fdc skill score vectors
%                .tt, .te, .et, .ee
%   JKGE        structure with scenario-wise JKGE vectors
%                .tt, .te, .et, .ee
%   loss_fnc    loss function (scalar)
%   ax          OPTIONAL: structure with graphics handles
%
% NOTES:
%   Left ECDF panels always correspond to scenario tt.
%   Right ECDF panels correspond to:
%     mode 1 -> unavailable
%     mode 2 -> te
%     mode 3 -> et
%     mode 4 -> ee
%   RSS and S_IB histories always plot tt on the left axis and the
%   mode-dependent comparison scenario on the right axis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Mar. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 9 ...
            || isempty(ax)
        ax = struct();
    end
    
    % ----------------
    % Basic formatting
    % ----------------
    fnt_ax = 15;    % tick labels +2
    fnt_lab = 16;   % axis labels +2
    fnt_tit = 17;
    fnt_med = 14;   % median text +2
    fnt_leg = 16;   % legend +2
    lw = 1.4;
    fa = 0.20;
    xL = -1; xR = 1;
    yL = 0;  yR = 1;

    cT = [0 0 1];
    %cE = [1 0.5 0];
    %cE = [0.75 0.00 0.00]; % crimson: right y-axis
    cE = [0.00 0.55 0.00]; 	% green
    cAx = [0.15 0.15 0.15]; % 2nd x-axis middle history panels
    
    % ----------------------
    % Model names and colors
    % ----------------------
    if ~isfield(mdl,'names') ...
            || isempty(mdl.names)
        error(['      Error: print_figs_new: ' ...
            'mdl.names is missing or empty.']);
    end
    if ~isfield(mdl,'model') ...
            || isempty(mdl.model)
        error(['      Error: print_figs_new: ' ...
            'mdl.model is missing or empty.']);
    end
    
    model_names = mdl.names;
    model = mdl.model;
    
    if iscell(model_names)
        mdl_name_raw = char( ...
            model_names{model});
    else
        mdl_name_raw = char( ...
            string(model_names(model)));
    end
    
    if strcmpi(mdl_name_raw,'xinanjiang')
        mdl_name_disp = 'Xinanjiang';
    else
        mdl_name_disp = upper(mdl_name_raw);
    end
    mdl_tex = strrep(mdl_name_disp,'_','\_');
    
    colors = [ ...
        0.0000 0.4470 0.7410;
        0.9500 0.4500 0.0000;
        0.9290 0.6940 0.1250;
        0.6000 0.2000 0.8000;
        0.0000 0.5500 0.0000;
        0.6500 0.3500 0.0000;
        0.0000 0.7000 0.7000;
        0.9000 0.0000 0.9000;
        0.6500 0.6500 0.6500];
    
    if model >= 1 ...
            && model <= size(colors,1)
        c = colors(model,:);
    else
        c = [0 0 0];
    end
    
    % ---------------------------
    % Split method / period label
    % ---------------------------
    if isfield(mdl,'sp_method') ...
            && ~isempty(mdl.sp_method)
        sp_method = lower(string(mdl.sp_method));
    else
        sp_method = "manual";
    end
    
    if any(strcmp(sp_method, ...
            ["manual","traditional_block","block", ...
            "deterministic_block","random_block", ...
            "deterministic_kfold"]))
        samp_word = 'period';
    else
        samp_word = 'mask';
    end
    
    nameT = sprintf(['Train basins ' ...
        '| train %s'],samp_word);
    
    % ------------------------------------------
    % Comparison scenario priority: ee > et > te
    % ------------------------------------------
    [rightField,nameE,klabel_right,rightAvailable] = ...
        choose_compare_scenario(prf,samp_word);
    
    % ----------------
    % Pull metric data
    % ----------------
    NSEtt = get_metric_vector(NSE,'tt');
    KGEtt = get_metric_vector(KGE,'tt');
    Sfdctt = get_metric_vector(S_fdc,'tt');
    JKGEtt = get_metric_vector(JKGE,'tt');
    
    if rightAvailable
        NSEr = get_metric_vector(NSE,char(rightField));
        KGEr = get_metric_vector(KGE,char(rightField));
        Sfdcr = get_metric_vector(S_fdc,char(rightField));
        JKGEr = get_metric_vector(JKGE,char(rightField));
    else
        NSEr = [];
        KGEr = [];
        Sfdcr = [];
        JKGEr = [];
    end
    
    [NSEtt_use,nt_nse] = finite_vec(NSEtt);
    [KGEtt_use,~] = finite_vec(KGEtt);
    [Sfdctt_use,~] = finite_vec(Sfdctt);
    [JKGEtt_use,~] = finite_vec(JKGEtt);
    
    [NSEr_use,nr_nse] = finite_vec(NSEr);
    [KGEr_use,~] = finite_vec(KGEr);
    [Sfdcr_use,~] = finite_vec(Sfdcr);
    [JKGEr_use,~] = finite_vec(JKGEr);
    
    % -------------------------------------------------
    % Bottom ECDF: JKGE when optimized, otherwise S_FDC
    % -------------------------------------------------
    useJKGEpanel = (loss_fnc == 7);

    if useJKGEpanel
        bottomT = JKGEtt_use;
        bottomR = JKGEr_use;
        bottomName = 'JKGE';
    else
        bottomT = Sfdctt_use;
        bottomR = Sfdcr_use;
        bottomName = 'S_fdc';
    end

    % ----------------------------------------------------
    % Horizontal limits for bottom JKGE / S_fdc ECDF panel
    % ----------------------------------------------------
    % if useJKGEpanel
    %     xLb = -1;
    %     xRb = 1;
    % else
    %     z = [bottomT(:); bottomR(:)];
    %     z = z(isfinite(z));
    %     xLb = 0;
    %     if isempty(z)
    %         xRb = 1;
    %     else
    %         xmaxD = max(z);
    %         if xmaxD <= 0
    %             xRb = 1;
    %         else
    %             xRb = 1.15*xmaxD;
    %         end
    %     end
    % end
    % JKGE and S_fdc both use [-1,1]
    xLb = -1;
    xRb = 1;

    % -----------------------
    % Loss label / loss field
    % -----------------------
    [~,lossLeftLabel,lossField] = ...
        get_loss_strings(loss_fnc);
    
    K_t = prf.K_t;
    K_e = prf.K_e;
    K_right = K_e;
    if rightAvailable && strcmp(string(rightField),"te")
        K_right = K_t;
    end
    
    % -------------------
    % Need initialization
    % -------------------
    needInit = isempty(ax) ...
        || ~isstruct(ax) ...
        || ~isfield(ax,'fig') ...
        || ~isgraphics(ax.fig);
    
    if ~needInit
        req = {'nse_t','kge_t', ...
            'jkge_t','lossAx', ...
            'rssAx','sibAx', ...
            'nse_r','kge_r','jkge_r'};
        for k = 1:numel(req)
            if ~isfield(ax,req{k}) ...
                    || ~isgraphics(ax.(req{k}))
                needInit = true;
                break
            end
        end
    end
    
    if needInit
        ax = struct();
        scr = get(0,'ScreenSize');
        
        figW = 0.95 * scr(3);
        figH = 0.95 * scr(4);
        
        figX = scr(1) + 0.5 * (scr(3) - figW);
        figY = scr(2) + 0.5 * (scr(4) - figH);
        
        ax.fig = figure( ...
            'Units','pixels', ...
            'Color','w', ...
            'Name',sprintf('%s: SAGE diagnostics',mdl_name_disp), ...
            'NumberTitle','off', ...
            'Position',[figX figY figW figH], ...
            'SizeChangedFcn', ...
            @(src,evt) refresh_history_top_frames(src));
        % --------------------------------------------------------
        % Manual layout with more whitespace around history panels
        % --------------------------------------------------------
        marginL = 0.05;
        marginR = 0.05;
        marginT = 0.04;
        % Keep the bottom history tick labels and xlabel inside the figure
        % canvas during PNG/PPTX export.
        marginB = 0.115;
        gapX1 = 0.08;   % left gap: left ECDF -> history
        gapX2 = 0.08;   % right gap: history -> right ECDF
        gapY = 0.04;
    
        totalH = 1 - marginT - marginB;
        rowH = (totalH - 2*gapY) / 3;
    
        colW_ecdf = 0.175;
    
        % Fixed left and right ECDF columns
        xLcol = marginL;
        xRcol = 1 - marginR - colW_ecdf;
    
        % Fixed history start
        xMcol = xLcol + colW_ecdf + gapX1;
    
        % History width fills the space up to the fixed right ECDF column
        colW_mid0 = xRcol - gapX2 - xMcol;
        colW_mid = 0.92 * colW_mid0;   % reduce width by 15%
    
        y3 = marginB;
        y2 = y3 + rowH + gapY;
        y1 = y2 + rowH + gapY;
    
        % Left ECDF column
        ax.nse_t = axes('Parent',ax.fig, ...
            'Position',[xLcol y1 colW_ecdf rowH]); 
        hold(ax.nse_t,'on');
        ax.kge_t = axes('Parent',ax.fig, ...
            'Position',[xLcol y2 colW_ecdf rowH]); 
        hold(ax.kge_t,'on');
        ax.jkge_t = axes('Parent',ax.fig, ...
            'Position',[xLcol y3 colW_ecdf rowH]); 
        hold(ax.jkge_t,'on');
    
        % Middle history column
        ax.lossAx = axes('Parent',ax.fig, ...
            'Position',[xMcol y1 colW_mid rowH]); 
        hold(ax.lossAx,'on');
        ax.rssAx = axes('Parent',ax.fig, ...
            'Position',[xMcol y2 colW_mid rowH]); 
        hold(ax.rssAx,'on');
        ax.sibAx = axes('Parent',ax.fig, ...
            'Position',[xMcol y3 colW_mid rowH]); 
        hold(ax.sibAx,'on');
    
        setappdata(ax.fig, ...
            'HistoryAxesHandles', ...
            [ax.lossAx ax.rssAx ax.sibAx]);
        setappdata(ax.fig, ...
            'HistoryTopFrameColor',cAx);
    
        % Right ECDF column
        ax.nse_r = axes('Parent',ax.fig, ...
            'Position',[xRcol y1 colW_ecdf rowH]); 
        hold(ax.nse_r,'on');
        ax.kge_r = axes('Parent',ax.fig, ...
            'Position',[xRcol y2 colW_ecdf rowH]); 
        hold(ax.kge_r,'on');
        ax.jkge_r = axes('Parent',ax.fig, ...
            'Position',[xRcol y3 colW_ecdf rowH]); 
        hold(ax.jkge_r,'on');
    
        % ECDF axes style
        ecdf_axes = [ax.nse_t ...
            ax.kge_t ax.jkge_t ...
            ax.nse_r ax.kge_r ...
            ax.jkge_r];
        set(ecdf_axes, ...
            'FontSize',fnt_ax, ...
            'LineWidth',1, ...
            'TickDir','out', ...
            'Box','off', ...
            'Layer','top', ...
            'XMinorTick','off', ...
            'YMinorTick','off', ...
            'TickLength',[0.030 0.030]);   % doubled tick length
    
        % Left ECDF column uses left-axis color
        set([ax.nse_t ax.kge_t ax.jkge_t], ...
            'XColor',cT, ...
            'YColor',cT);
        
        % Right ECDF column uses right-axis color
        set([ax.nse_r ax.kge_r ax.jkge_r], ...
            'XColor',cE, ...
            'YColor',cE);
    
        for ah = ecdf_axes
            xlim(ah,[xL xR]);
            ylim(ah,[yL yR]);
        end
        % Overwrite if other than x in [-1,1]
        xlim(ax.jkge_t,[xLb xRb]);
        xlim(ax.jkge_r,[xLb xRb]);

        % No xtick labels on top 2 rows
        set(ax.nse_t,'XTickLabel',[]);
        set(ax.kge_t,'XTickLabel',[]);
        set(ax.nse_r,'XTickLabel',[]);
        set(ax.kge_r,'XTickLabel',[]);
    
        % Labels
        ylabel(ax.nse_t ,'$F({\rm NSE}_{\rm tt})$', ...
            'Interpreter','latex', ...
            'FontSize',fnt_lab);
        ylabel(ax.kge_t ,'$F({\rm KGE}_{\rm tt})$', ...
            'Interpreter','latex', ...
            'FontSize',fnt_lab);

        if useJKGEpanel
            ylabel(ax.jkge_t, ...
                '$F({\rm JKGE}_{\rm tt})$', ...
                'Interpreter','latex', ...
                'FontSize',fnt_lab);
            xlabel(ax.jkge_t, ...
                '${\rm JKGE}_{\rm tt}$', ...
                'Interpreter','latex', ...
                'FontSize',fnt_lab);
        else
            ylabel(ax.jkge_t, ...
                '$F(S_{{\rm FDC}_{\rm tt}})$', ...
                'Interpreter','latex', ...
                'FontSize',fnt_lab);
            xlabel(ax.jkge_t, ...
                '$S_{{\rm FDC}_{\rm tt}}$', ...
                'Interpreter','latex', ...
                'FontSize',fnt_lab);
        end
    
        if rightAvailable
            ylabel(ax.nse_r, ...
                sprintf(['$F({\\rm NSE}_' ...
                '{\\rm %s})$'],char(rightField)), ...
                'Interpreter','latex', ...
                'FontSize',fnt_lab);
            ylabel(ax.kge_r, ...
                sprintf(['$F({\\rm KGE}_' ...
                '{\\rm %s})$'],char(rightField)), ...
                'Interpreter','latex', ...
                'FontSize',fnt_lab);
            % ylabel(ax.jkge_r, ...
            %     sprintf(['$F({\\rm JKGE}_' ...
            %     '{\\rm %s})$'],char(rightField)), ...
            %     'Interpreter','latex', ...
            %     'FontSize',fnt_lab);
            % xlabel(ax.jkge_r, ...
            %     sprintf(['${\\rm JKGE}_' ...
            %     '{\\rm %s}$'],char(rightField)), ...
            %     'Interpreter','latex', ...
            %     'FontSize',fnt_lab);
            if useJKGEpanel
                ylabel(ax.jkge_r, ...
                    sprintf('$F({\\rm JKGE}_{\\rm %s})$', ...
                    char(rightField)), ...
                    'Interpreter','latex', ...
                    'FontSize',fnt_lab);
                xlabel(ax.jkge_r, ...
                    sprintf('${\\rm JKGE}_{\\rm %s}$', ...
                    char(rightField)), ...
                    'Interpreter','latex', ...
                    'FontSize',fnt_lab);
            else
                ylabel(ax.jkge_r, ...
                    sprintf('$F(S_{{\\rm FDC}_{\\rm %s}})$', ...
                    char(rightField)), ...
                    'Interpreter','latex', ...
                    'FontSize',fnt_lab);
                xlabel(ax.jkge_r, ...
                    sprintf('$S_{{\\rm FDC}_{\\rm %s}}$', ...
                    char(rightField)), ...
                    'Interpreter','latex', ...
                    'FontSize',fnt_lab);
            end
            xlabel(ax.nse_r ,'');
            xlabel(ax.kge_r ,'');
        else
            ylabel(ax.nse_r,['$F({\rm NSE}_' ...
                '{\rm na})$'] , ...
                'Interpreter','latex', ...
                'FontSize',fnt_lab);
            ylabel(ax.kge_r,['$F({\rm KGE}_' ...
                '{\rm na})$'] , ...
                'Interpreter','latex', ...
                'FontSize',fnt_lab);

            if useJKGEpanel
                ylabel(ax.jkge_r, ...
                    '$F({\rm JKGE}_{\rm na})$', ...
                    'Interpreter','latex', ...
                    'FontSize',fnt_lab);
                xlabel(ax.jkge_r, ...
                    '${\rm JKGE}_{\rm na}$', ...
                    'Interpreter','latex', ...
                    'FontSize',fnt_lab);
            else
                ylabel(ax.jkge_r, ...
                    '$F(S_{{\rm FDC},{\rm na}})$', ...
                    'Interpreter','latex', ...
                    'FontSize',fnt_lab);
                xlabel(ax.jkge_r, ...
                    '$S_{{\rm FDC},{\rm na}}$', ...
                    'Interpreter','latex', ...
                    'FontSize',fnt_lab);
            end
            % ylabel(ax.jkge_r,['$F({\rm JKGE}_' ...
            %     '{\rm na})$'], ...
            %     'Interpreter','latex', ...
            %     'FontSize',fnt_lab);
            % xlabel(ax.jkge_r,['${\rm JKGE}_' ...
            %     '{\rm na}$'], ...
            %     'Interpreter','latex', ...
            %     'FontSize',fnt_lab);
        end
    
        % xlabel(ax.jkge_t,['${\rm JKGE}_' ...
        %     '{\rm tt}$'], ...
        %     'Interpreter','latex', ...
        %     'FontSize',fnt_lab);
        xlabel(ax.nse_t ,'');
        xlabel(ax.kge_t ,'');
    
        % Top/right frame lines
        add_top_right_frame(ax.nse_t, ...
            xL,xR,yL,yR,cT);
        add_top_right_frame(ax.kge_t, ...
            xL,xR,yL,yR,cT);
        add_top_right_frame(ax.jkge_t, ...
            xLb,xRb,yL,yR,cT);
        
        add_top_right_frame(ax.nse_r, ...
            xL,xR,yL,yR,cE);
        add_top_right_frame(ax.kge_r, ...
            xL,xR,yL,yR,cE);
        add_top_right_frame(ax.jkge_r, ...
            xLb,xRb,yL,yR,cE);
    
        % Titles: only top row ECDF
        set_panel_title(ax.nse_t,mdl_tex, ...
            nameT,K_t,nt_nse,fnt_tit, ...
            'K_{\rm t}',true,false);
        set_panel_title(ax.nse_r,mdl_tex, ...
            nameE,K_right,nr_nse,fnt_tit, ...
            klabel_right,rightAvailable,false);
        ax.nse_t.Title.Color = cT;
        ax.nse_r.Title.Color = cE;
    
        title(ax.kge_t ,'');
        title(ax.jkge_t,'');
        title(ax.kge_r ,'');
        title(ax.jkge_r,'');
    
        % History axes
        init_history_axis(ax.lossAx, ...
            fnt_ax);
        init_history_axis(ax.rssAx, ...
            fnt_ax);
        init_history_axis(ax.sibAx, ...
            fnt_ax);
    
        % Only top history title
        title(ax.lossAx, ...
            sprintf(['\\texttt{%s}: ' ...
            'History of loss functions'],mdl_tex), ...
            'Interpreter','latex','FontSize',fnt_tit);
        title(ax.rssAx,'');
        title(ax.sibAx,'');
    
        xlabel(ax.lossAx,'');
        xlabel(ax.rssAx ,'');
        xlabel(ax.sibAx,'SAGE iteration, $i$', ...
            'Interpreter','latex', ...
            'FontSize',fnt_lab);
        set(ax.lossAx,'XTickLabel',[]);
        set(ax.rssAx,'XTickLabel',[]);
    
        % History lines
        yyaxis(ax.lossAx,'left');
        clrlegT = char(strcat('\color[rgb]{', ...
            num2str(cT(1)),',',num2str(cT(2)), ...
            ',',num2str(cT(3)),'}',{' '},nameT));
        ax.lossLineT = plot(ax.lossAx, ...
            nan,nan,'-o', ...
            'Color',cT, ...
            'MarkerFaceColor',cT, ...
            'LineWidth',1.3, ...
            'DisplayName',clrlegT);
           % 'DisplayName',['\color[rgb]{0,0,1} ' nameT]);
           % 'DisplayName',nameT);
        yyaxis(ax.lossAx,'right');
        clrlegE = char(strcat('\color[rgb]{', ...
            num2str(cE(1)),',',num2str(cE(2)), ...
            ',',num2str(cE(3)),'}',{' '},nameE));
        ax.lossLineR = plot(ax.lossAx, ...
            nan,nan,'-s', ...
            'Color',cE, ...
            'MarkerFaceColor',cE, ...
            'LineWidth',1.3, ...
            'DisplayName',clrlegE);
           % 'DisplayName',['\color[rgb]{0.75,0,0} ' nameE]);
           % 'DisplayName',nameE);
        
        yyaxis(ax.rssAx,'left');
        ax.rssLineT = plot(ax.rssAx, ...
            nan,nan,'-o', ...
            'Color',cT, ...
            'MarkerFaceColor',cT, ...
            'LineWidth',1.3, ...
            'HandleVisibility','off');
        yyaxis(ax.rssAx,'right');
        ax.rssLineR = plot(ax.rssAx, ...
            nan,nan,'-s', ...
            'Color',cE, ...
            'MarkerFaceColor',cE, ...
            'LineWidth',1.3, ...
            'HandleVisibility','off');
    
        yyaxis(ax.sibAx,'left');
        ax.sibLineT = plot(ax.sibAx, ...
            nan,nan,'-o', ...
            'Color',cT, ...
            'MarkerFaceColor',cT, ...
            'LineWidth',1.3, ...
            'HandleVisibility','off');
        yyaxis(ax.sibAx,'right');
        ax.sibLineR = plot(ax.sibAx, ...
            nan,nan,'-s', ...
            'Color',cE, ...
            'MarkerFaceColor',cE, ...
            'LineWidth',1.3, ...
            'HandleVisibility','off');
    
        % Only top history panel has legend, inside the axes
        lgd = legend(ax.lossAx, ...
            'Location','northeast', ...
            'Interpreter','tex', ...
            'Box','off');
        lgd.FontSize = fnt_leg;
    
        % ECDF handles
        tags = {'nse_t','kge_t', ...
            'jkge_t','nse_r','kge_r','jkge_r'};
        for k = 1:numel(tags)
            tg = tags{k};
            ax.ecdf.(tg).patch = ...
                gobjects(1);
            ax.ecdf.(tg).line = ...
                gobjects(1);
            ax.ecdf.(tg).med_v = ...
                gobjects(1);
            ax.ecdf.(tg).med_h = ...
                gobjects(1);
            ax.ecdf.(tg).med_s = ...
                gobjects(1);
            ax.ecdf.(tg).med_txt = ...
                gobjects(1);
        end
    end
    
    % Top-row titles update only
    set_panel_title(ax.nse_t, ...
        mdl_tex,nameT,K_t,nt_nse, ...
        fnt_tit,'K_{\rm t}',true,false);
    set_panel_title(ax.nse_r, ...
        mdl_tex,nameE,K_right,nr_nse, ...
        fnt_tit,klabel_right, ...
        rightAvailable,false);
    
    ax.nse_t.Title.Color = cT;
    ax.nse_r.Title.Color = cE;
    
    % ------------
    % Update ECDFs
    % ------------
    ax = update_manual_ecdf(ax, ...
        'nse_t',ax.nse_t,NSEtt_use, ...
        c,xL,xR,fa,lw,'NSE','tt', ...
        fnt_med);
    ax = update_manual_ecdf(ax, ...
        'kge_t',ax.kge_t,KGEtt_use, ...
        c,xL,xR,fa,lw,'KGE','tt', ...
        fnt_med);
    % ax = update_manual_ecdf(ax, ...
    %     'jkge_t',ax.jkge_t,JKGEtt_use, ...
    %     c,xL,xR,fa,lw,'JKGE','tt', ...
    %     fnt_med);
    ax = update_manual_ecdf(ax, ...
        'jkge_t',ax.jkge_t,bottomT, ...
        c,xLb,xRb,fa,lw,bottomName,'tt', ...
        fnt_med);
    
    if rightAvailable
        ax = update_manual_ecdf(ax, ...
            'nse_r',ax.nse_r ,NSEr_use, ...
            c,xL,xR,fa,lw,'NSE', ...
            char(rightField),fnt_med);
        ax = update_manual_ecdf(ax, ...
            'kge_r',ax.kge_r ,KGEr_use, ...
            c,xL,xR,fa,lw,'KGE', ...
            char(rightField),fnt_med);
        % ax = update_manual_ecdf(ax, ...
        %     'jkge_r',ax.jkge_r,JKGEr_use, ...
        %     c,xL,xR,fa,lw,'JKGE', ...
        %     char(rightField),fnt_med);
        ax = update_manual_ecdf(ax, ...
            'jkge_r',ax.jkge_r,bottomR, ...
            c,xLb,xRb,fa,lw,bottomName, ...
            char(rightField),fnt_med);
    else
        clear_manual_ecdf(ax,'nse_r');  
        set_unavailable_label(ax.nse_r,true);
        clear_manual_ecdf(ax,'kge_r');  
        set_unavailable_label(ax.kge_r,true);
        clear_manual_ecdf(ax,'jkge_r'); 
        set_unavailable_label(ax.jkge_r,true);
    end
    
    % % Refresh the manually drawn top/right frames accordingly.
    % add_top_right_frame(ax.jkge_t, ...
    %     xLb,xRb,yL,yR,cT);
    % add_top_right_frame(ax.jkge_r, ...
    %     xLb,xRb,yL,yR,cE);

    % ----------------
    % Update histories
    % ----------------
    [LossT,LossR] = get_history_pair( ...
        prf,lossField,i,rightField, ...
        K_t,K_e,rightAvailable);
    [RSST,RSSR] = get_history_pair( ...
        prf,'RSS',i,rightField, ...
        K_t,K_e,rightAvailable);
    [SibT,SibR] = get_history_pair( ...
        prf,'Sib',i,rightField, ...
        1,1,rightAvailable);
    
    update_dual_history_loss( ...
        ax.lossAx,ax.lossLineT, ...
        ax.lossLineR, ...
        LossT,LossR,cT,cE, ...
        lossLeftLabel,right_loss_label( ...
        loss_fnc,rightField, ...
        rightAvailable),fnt_lab);
    
    update_dual_history(ax.rssAx, ...
        ax.rssLineT,ax.rssLineR, ...
        RSST,RSSR,cT,cE, ...
        ['$\Sigma \mathrm{RSS}_' ...
        '{\rm tt}$'],right_label( ...
        'RSS',rightField, ...
        rightAvailable),fnt_lab);
    
    update_dual_history(ax.sibAx, ...
        ax.sibLineT,ax.sibLineR, ...
        SibT,SibR,cT,cE, ...
        ['$\widehat{\mathcal{S}}_' ...
        '{\mathrm{IB_{\rm tt}}}$'], ...
        right_label('Sib', ...
        rightField,rightAvailable), ...
        fnt_lab);
    
    xmax = max([2 numel(LossT) ...
        numel(LossR) numel(RSST) ...
        numel(RSSR) numel(SibT) ...
        numel(SibR)]); 
    
    xt = local_history_xticks(xmax);
    
    xlim(ax.lossAx,[1 xmax]);
    xlim(ax.rssAx,[1 xmax]);
    xlim(ax.sibAx,[1 xmax]);
    
    set_history_xticks(ax.lossAx,xt,false);
    set_history_xticks(ax.rssAx,xt,false);
    set_history_xticks(ax.sibAx,xt,true);
    
    % Add horizontal axis on top
    if i == 1 ...
            || needInit
        add_history_top_frame(ax.lossAx,cAx);
        add_history_top_frame(ax.rssAx,cAx);
        add_history_top_frame(ax.sibAx,cAx);
    end
    
    if ~isdeployed
        drawnow expose
        pause(0.001)
    else
        drawnow limitrate nocallbacks
    end

end

% =================
% Helpers functions
% =================
function [rightField,nameE,klabel_right, ...
    rightAvailable] = choose_compare_scenario(prf,samp_word)
    % ECDFs use the current basin-wise metric vectors. Since the prf
    % refactor, these live in prf.curr; iteration histories such as RSS
    % live separately in prf.iter and must not determine ECDF availability.
    if isfield(prf,'curr') && isfield(prf.curr,'NSE')
        metric = prf.curr.NSE;
    else
        metric = struct();
    end

    if has_nonempty_scenario(metric,'ee')
        rightField = "ee";
        nameE = sprintf(['Eval basins |' ...
            ' eval %s'],samp_word);
        klabel_right = 'K_{\rm e}';
        rightAvailable = true;
    elseif has_nonempty_scenario(metric,'et')
        rightField = "et";
        nameE = sprintf(['Eval basins |' ...
            ' train %s'],samp_word);
        klabel_right = 'K_{\rm e}';
        rightAvailable = true;
    elseif has_nonempty_scenario(metric,'te')
        rightField = "te";
        nameE = sprintf(['Train basins |' ...
            ' eval %s'],samp_word);
        klabel_right = 'K_{\rm t}';
        rightAvailable = true;
    else
        rightField = "";
        nameE = 'Unavailable';
        klabel_right = '';
        rightAvailable = false;
    end
end

function tf = has_nonempty_scenario(S,scn)
    tf = false;
    if isstruct(S) && isfield(S,scn)
        x = S.(scn);
        tf = ~isempty(x) ...
            && any(isfinite(x));
    end
end

function x = get_metric_vector(S,scn)
    x = [];
    if isstruct(S) && isfield(S,scn)
        x = S.(scn)(:);
    end
end

function [x,n] = finite_vec(v)
    if isempty(v) 
        x = []; n = 0; 
        return; 
    end
    v = v(:);
    I = isfinite(v);
    x = v(I);
    n = nnz(I);
end

function [lossTitle,leftLab,fld] = ...
    get_loss_strings(loss_fnc)
    switch loss_fnc
        case 1
            lossTitle = '\Sigma SAR';
            leftLab = ['$\Sigma {\rm SAR}_' ...
                '{\rm tt}$'];
            fld = 'SAR';
        case 2
            lossTitle = '\Sigma GLS';
            leftLab = ['$\Sigma {\rm GLS}_' ...
                '{\rm tt}$'];
            fld = 'GLS';
        case 3
            lossTitle = '\Sigma(1-NSE)';
            leftLab = ['$\Sigma(1-{\rm NSE})_' ...
                '{\rm tt}$'];
            fld = 'L';
        case 4
            lossTitle = '\Sigma(1-KGE)';
            leftLab = ['$\Sigma(1-{\rm KGE})_' ...
                '{\rm tt}$'];
            fld = 'L';
        case 5
            lossTitle = '\Sigma Huber';
            leftLab = ['$\Sigma {\rm Huber}_' ...
                '{\rm tt}$'];
            fld = 'Huber';
        case 6
            lossTitle = '\Sigma d_{FDC}';
            leftLab = ['$\Sigma d_{{\rm FDC},' ...
                '\rm tt}$'];
            fld = 'L';
        case 7
            lossTitle = '\Sigma(1-JKGE)';
            leftLab = ['$\Sigma(1-{\rm JKGE})_' ...
                '{\rm tt}$'];
            fld = 'L';
        otherwise
            error('Unknown loss function.');
    end
end

function init_history_axis(axh,fs)
    set(axh,'FontSize',fs, ...
        'LineWidth',1, ...
        'TickDir','out', ...
        'Box','off');
    
    xtickformat(axh,'%d');
    
    try
        axh.YAxis(1).Exponent = 0;
        axh.YAxis(2).Exponent = 0;
    catch
    end

end

function set_panel_title(axh, ...
    mdl_tex,str,K,nFinite,fnt, ...
    klabel,isAvailable,showModel)

    if nargin < 8 ...
            || isempty(showModel)
        showModel = true;
    end
    
    if ~isAvailable
        if showModel
            axh.Title.String = ...
                sprintf('\\texttt{%s}: %s', ...
                mdl_tex,str);
        else
            axh.Title.String = str;
        end
    elseif isempty(klabel)
        if showModel
            axh.Title.String = ...
                sprintf('\\texttt{%s}: %s', ...
                mdl_tex,str);
        else
            axh.Title.String = str;
        end
    else
        if showModel
            axh.Title.String = ...
                sprintf(['\\texttt{%s}: %s ' ...
                '($%s = %d$)'],mdl_tex, ...
                str,klabel,K);
        else
            axh.Title.String = ...
                sprintf('%s ($%s = %d$)', ...
                str,klabel,K);
        end
        if nFinite ~= K
            axh.Title.String = sprintf( ...
                '%s; $n_{\rm NSE} = %d$', ...
                axh.Title.String,nFinite);
        end
    end
    axh.Title.Interpreter = 'latex';
    axh.Title.FontSize = fnt;
end

function h = add_top_right_frame(axh, ...
    xL,xR,yL,yR,clr)
    hold(axh,'on');
    if isappdata(axh, ...
            'TopRightFrameHandles')
        hh = getappdata(axh, ...
            'TopRightFrameHandles');
        try
            if isfield(hh,'top') ...
                    && isgraphics(hh.top)
                delete(hh.top); 
            end
            if isfield(hh,'right') ...
                    && isgraphics(hh.right)
                delete(hh.right); 
            end
        catch
        end
    end
    h.top = line(axh,[xL xR],[yR yR], ...
        'Color',clr, ...
        'LineWidth',1, ...
        'HandleVisibility','off', ...
        'Clipping','off');
    h.right = line(axh,[xR xR],[yL yR], ...
        'Color',clr, ...
        'LineWidth',1, ...
        'HandleVisibility','off', ...
        'Clipping','off');
    setappdata(axh,'TopRightFrameHandles',h);
end

function ax = update_manual_ecdf(ax, ...
    tag,axh,z,c,xL,xR,fa,lw, ...
    metName,scn,fnt_med)
    if isempty(z)
        set_unavailable_label(axh,true);
        clear_manual_ecdf(ax,tag);
        return
    else
        set_unavailable_label(axh,false);
    end
    
    [f,x] = sage_ecdf(z);
    [xs,fs] = ecdf_to_stairs_fixed(x,f,xL,xR);
    [xp,yp] = stairs_fill_poly(xs,fs);
    
    if ~isgraphics(ax.ecdf.(tag).patch)
        ax.ecdf.(tag).patch = patch(axh, ...
            xp,yp,c, ...
            'FaceAlpha',fa, ...
            'EdgeAlpha',0);
        ax.ecdf.(tag).line = plot(axh, ...
            xs,fs, ...
            'Color',c, ...
            'LineWidth',lw);
    else
        set(ax.ecdf.(tag).patch, ...
            'XData',xp, ...
            'YData',yp, ...
            'FaceColor',c);
        set(ax.ecdf.(tag).line, ...
            'XData',xs, ...
            'YData',fs, ...
            'Color',c);
    end
    
    medX = median(z);
    medF = ecdf_value_from_stairs( ...
        xs,fs,medX);
    [xMark,yMark,x1m,x2m,xt,ha,mode] = ...
        median_marker_geom(medX,medF,xL,xR);
    
    if ~isgraphics(ax.ecdf.(tag).med_v)
        ax.ecdf.(tag).med_v = line(axh, ...
            nan,nan, ...
            'Color','k', ...
            'LineWidth',1.0, ...
            'HandleVisibility','off');
        ax.ecdf.(tag).med_h = line(axh, ...
            nan,nan, ...
            'Color','k', ...
            'LineWidth',1.0, ...
            'HandleVisibility','off');
        ax.ecdf.(tag).med_s = line(axh, ...
            nan,nan, ...
            'Marker','s', ...
            'MarkerFaceColor','k', ...
            'MarkerEdgeColor','k', ...
            'LineStyle','none', ...
            'MarkerSize',5, ...
            'HandleVisibility','off');
        ax.ecdf.(tag).med_txt = ...
            text(axh, ...
            nan,nan,'', ...
            'Interpreter','latex', ...
            'FontWeight','bold', ...
            'HorizontalAlignment',ha, ...
            'VerticalAlignment','middle', ...
            'Color','k', ...
            'HandleVisibility','off', ...
            'FontSize',fnt_med);
    end
    
    if ~mode.specialHalfAtLeft
        set(ax.ecdf.(tag).med_v, ...
            'XData',[xMark xMark], ...
            'YData',[0 yMark], ...
            'Visible','on');
    else
        set(ax.ecdf.(tag).med_v, ...
            'XData',[], ...
            'YData',[], ...
            'Visible','off');
    end
    set(ax.ecdf.(tag).med_h, ...
        'XData',[x1m x2m], ...
        'YData',[yMark yMark], ...
        'Visible','on');
    set(ax.ecdf.(tag).med_s, ...
        'XData',xMark, ...
        'YData',yMark, ...
        'Visible','on');
    if strcmpi(metName,'S_fdc')
        medLabel = sprintf( ...
            '$\\widehat{T}_{S_{{{\\rm FDC}_{\\rm %s}}}} = %.3f$', ...
            scn,medX);
    else
        medLabel = sprintf( ...
            ['$\\widehat{T}_' ...
            '{\\mathrm{%s}_{\\rm %s}} = %.3f$'], ...
            metName,scn,medX);
    end
    set(ax.ecdf.(tag).med_txt, ...
        'Position',[xt yMark 0], ...
        'String',medLabel, ...
        'HorizontalAlignment',ha, ...
        'Visible','on', ...
        'FontSize',fnt_med);
    
    xlim(axh,[xL xR]);
    ylim(axh,[0 1]);
    end
    
    function clear_manual_ecdf(ax,tag)
    flds = {'patch','line', ...
        'med_v','med_h','med_s','med_txt'};
    for k = 1:numel(flds)
        h = ax.ecdf.(tag).(flds{k});
        if isgraphics(h)
            if strcmp(flds{k},'med_txt')
                set(h,'String','', ...
                    'Visible','off');
            else
                set(h,'XData',[], ...
                    'YData',[], ...
                    'Visible','off');
            end
        end
    end
end

function [yT,yR] = get_history_pair(prf, ...
    fld,i,rightField,K_t,K_e,rightAvailable)

    % Iteration histories are stored under prf.iter. The generic Sib
    % diagnostic plotted here is the NSE-based integrated basin score.
    if strcmp(fld,'Sib')
        histFld = 'Sib_NSE';
        sT = 1;
    else
        histFld = fld;
        sT = K_t;
    end

    if ~isfield(prf,'iter') ...
            || ~isfield(prf.iter,histFld)
        error('print_SAGE:MissingHistory', ...
            'Missing iteration history prf.iter.%s.',histFld);
    end
    hist = prf.iter.(histFld);
    
    it = min(i,numel(hist.tt));
    iiT = 1:it;
    iiT = iiT(isfinite(hist.tt(iiT)));
    yT = sT * hist.tt(iiT);
    
    if ~rightAvailable
        yR = [];
        return
    end
    
    switch char(rightField)
        case 'te'
            sR = strcmp(fld,'Sib') * 1 ...
                + ~strcmp(fld,'Sib') * K_t;
        case {'et','ee'}
            sR = strcmp(fld,'Sib') * 1 ...
                + ~strcmp(fld,'Sib') * K_e;
        otherwise
            sR = 1;
    end
    
    rightField = char(rightField);
    iv = min(i,numel(hist.(rightField)));
    iiR = 1:iv;
    iiR = iiR(isfinite(hist.(rightField)(iiR)));
    yR = sR * hist.(rightField)(iiR);
end

function lab = right_label(fld, ...
    rightField,rightAvailable)
    if ~rightAvailable
        lab = 'Unavailable';
        return
    end
    
    switch fld
        case 'L'
            core = '\mathcal{L}';
        case 'RSS'
            core = '\Sigma \mathrm{RSS}';
        case 'Sib'
            lab = ...
                sprintf(['$\\widehat{\\mathcal{S}}_' ...
                '{\\mathrm{IB_{\\rm %s}}}$'], ...
                char(rightField));
            return
        case 'SAR'
            core = '\Sigma {\rm SAR}';
        case 'GLS'
            core = '\Sigma {\rm GLS}';
        case 'Huber'
            core = '\Sigma {\rm Huber}';
        otherwise
            core = fld;
    end
    lab = sprintf('$%s_{\\rm %s}$', ...
        core,char(rightField));
end

function lab = right_loss_label(loss_fnc, ...
    rightField,rightAvailable)

    if ~rightAvailable
        lab = 'Unavailable';
        return
    end
    
    scn = char(rightField);
    
    switch loss_fnc
        case 1
            lab = sprintf(['$\\Sigma ' ...
                '{\\rm SAR}_{\\rm %s}$'],scn);
        case 2
            lab = sprintf(['$\\Sigma ' ...
                '{\\rm GLS}_{\\rm %s}$'],scn);
        case 3
            lab = sprintf(['$\\Sigma' ...
                '(1-{\\rm NSE})_{\\rm %s}$'],scn);
        case 4
            lab = sprintf(['$\\Sigma' ...
                '(1-{\\rm KGE})_{\\rm %s}$'],scn);
        case 5
            lab = sprintf(['$\\Sigma ' ...
                '{\\rm Huber}_{\\rm %s}$'],scn);
        case 6
            lab = sprintf(['$\\Sigma ' ...
                'd_{\\mathrm{FDC},\\rm %s}$'],scn);
        case 7
            lab = sprintf(['$\\Sigma' ...
                '(1-{\\rm JKGE})_{\\rm %s}$'],scn);
        otherwise
            lab = sprintf(['$\\mathcal{L}_' ...
                '{\\rm %s}$'],scn);
    end

end

function update_dual_history(axh, ...
    hT,hR,yT,yR,cT,cE,leftLab, ...
    rightLab,fsLab)

    useLogT = should_use_log_history_scale(yT(:));
    useLogR = should_use_log_history_scale(yR(:));

    update_history_side(axh,'left',hT,yT,cT, ...
        leftLab,fsLab,useLogT,false);
    update_history_side(axh,'right',hR,yR,cE, ...
        rightLab,fsLab,useLogR,true);
end

function update_history_side(axh,side,hLine,y,clr, ...
    axisLabel,fsLab,useLog,allowUnavailable)

    yyaxis(axh,side);
    axh.YColor = clr;
    try
        if strcmp(side,'left')
            axh.YAxis(1).Exponent = 0;
        else
            axh.YAxis(2).Exponent = 0;
        end
    catch
    end

    if isempty(y)
        set(hLine,'XData',nan,'YData',nan);
        axh.YScale = 'linear';
        if allowUnavailable
            ylabel(axh,'Unavailable', ...
                'Interpreter','none', ...
                'FontSize',fsLab);
        else
            ylabel(axh,axisLabel, ...
                'Interpreter','latex', ...
                'FontSize',fsLab);
        end
        return
    end

    if useLog
        axh.YScale = 'log';
        yPlot = y;
        yPlot(~isfinite(yPlot) ...
            | yPlot <= 0) = nan;
        set(hLine, ...
            'XData',1:numel(yPlot), ...
            'YData',yPlot);
        local_set_dynamic_ylim_log(axh,side,yPlot);
        local_no_axis_exponent(axh);
        ylabel(axh,axisLabel, ...
            'Interpreter','latex', ...
            'FontSize',fsLab);
        return
    end

    axh.YScale = 'linear';
    [scalePow,scaleVal] = nice_eng_scale(y(:));
    yPlot = y./scaleVal;
    set(hLine, ...
        'XData',1:numel(yPlot), ...
        'YData',yPlot);
    local_set_dynamic_ylim(axh,side,yPlot);
    local_no_axis_exponent(axh);

    if scalePow ~= 0
        ylabel(axh, ...
            sprintf('%s $\\times 10^{%d}$', ...
            axisLabel,scalePow), ...
            'Interpreter','latex', ...
            'FontSize',fsLab);
    else
        ylabel(axh,axisLabel, ...
            'Interpreter','latex', ...
            'FontSize',fsLab);
    end
end
function [scalePow,scaleVal] = nice_eng_scale(y)
%NICE_ENG_SCALE Scale so plotted values stay visually compact.
% Example: 6e5 -> x10^4 gives 60 instead of x10^3 gives 600.

    y = y(isfinite(y));
    
    if isempty(y) ...
            || all(y == 0)
        scalePow = 0;
        scaleVal = 1;
        return
    end
    
    m = max(abs(y));
    
    if m >= 100
        scalePow = floor(log10(m));
    else
        scalePow = 0;
    end
    
    scaleVal = 10^scalePow;
end

function set_unavailable_label(axh,tf)
    key = 'SAGE_UnavailableText';
    if isappdata(axh,key)
        h = getappdata(axh,key);
    else
        h = text(axh,0.5,0.5,'Unavailable', ...
            'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontWeight','bold', ...
            'FontSize',18, ...
            'Color',[0.35 0.35 0.35], ...
            'Interpreter','none', ...
            'Visible','off');
        setappdata(axh,key,h);
    end
    if isgraphics(h)
        if tf
            set(h,'Visible','on');
        else
            set(h,'Visible','off');
        end
    end
end

function update_dual_history_loss(axh, ...
    hT,hR,yT,yR,cT,cE,leftLab, ...
    rightLab,fsLab)

    % useLog = should_use_log_loss_scale([yT(:); yR(:)], ...
    %     numel(yT));
    useLogT = should_use_log_loss_scale(yT(:),numel(yT));
    useLogR = should_use_log_loss_scale(yR(:),numel(yR));
    useLog = useLogT || useLogR;

    if useLog
        yyaxis(axh,'left');
        try
            axh.YAxis(1).Exponent = 0;
        catch
        end
        axh.YColor = cT;
        axh.YScale = 'log';
        if isempty(yT)
            set(hT,'XData',nan, ...
                'YData',nan);
        else
            yTplot = yT;
            yTplot(~isfinite(yTplot) ...
                | yTplot <= 0) = nan;
            set(hT,'XData',1:numel(yTplot), ...
                'YData',yTplot);
            local_set_dynamic_ylim_log(axh, ...
                'left',yTplot);
            local_no_axis_exponent(axh);
        end
        ylabel(axh,leftLab, ...
            'Interpreter','latex', ...
            'FontSize',fsLab);
    
        yyaxis(axh,'right');
        try
            axh.YAxis(2).Exponent = 0;
        catch
        end
        axh.YColor = cE;
        axh.YScale = 'log';
        if isempty(yR)
            set(hR,'XData',nan, ...
                'YData',nan);
            ylabel(axh,'Unavailable', ...
                'Interpreter','none', ...
                'FontSize',fsLab);
        else
            yRplot = yR;
            yRplot(~isfinite(yRplot) ...
                | yRplot <= 0) = nan;
            set(hR,'XData',1:numel(yRplot), ...
                'YData',yRplot);
            local_set_dynamic_ylim_log(axh, ...
                'right',yRplot);
            local_no_axis_exponent(axh);
            ylabel(axh,rightLab, ...
                'Interpreter','latex', ...
                'FontSize',fsLab);
        end
    else
        %[pow,sc] = nice_eng_scale([yT(:); yR(:)]);
        [powT,scT] = nice_eng_scale(yT(:));
        [powR,scR] = nice_eng_scale(yR(:));
        yyaxis(axh,'left');
        try
            axh.YAxis(1).Exponent = 0;
        catch
        end
        axh.YColor = cT;
        axh.YScale = 'linear';
        if isempty(yT)
            set(hT,'XData',nan, ...
                'YData',nan);
        else
            yTplot = yT./scT;
            set(hT,'XData',1:numel(yT), ...
                'YData',yTplot);
            local_set_dynamic_ylim(axh, ...
                'left',yTplot);
            local_no_axis_exponent(axh);
        end
        if powT ~= 0
            ylabel(axh, ...
                sprintf('%s $\\times 10^{%d}$', ...
                leftLab,powT), ...
                'Interpreter','latex', ...
                'FontSize',fsLab);
        else
            ylabel(axh,leftLab, ...
                'Interpreter','latex', ...
                'FontSize',fsLab);
        end
    
        yyaxis(axh,'right');
        try
            axh.YAxis(2).Exponent = 0;
        catch
        end
        axh.YColor = cE;
        axh.YScale = 'linear';
        if isempty(yR)
            set(hR,'XData',nan, ...
                'YData',nan);
            ylabel(axh,'Unavailable', ...
                'Interpreter','none', ...
                'FontSize',fsLab);
        else
            yRplot = yR./scR;
            set(hR,'XData',1:numel(yR), ...
                'YData',yRplot);
            local_set_dynamic_ylim(axh, ...
                'right',yRplot);
            local_no_axis_exponent(axh);
            if powR ~= 0
                ylabel(axh, ...
                    sprintf('%s $\\times 10^{%d}$', ...
                    rightLab,powR), ...
                    'Interpreter','latex', ...
                    'FontSize',fsLab);
            else
                ylabel(axh,rightLab, ...
                    'Interpreter','latex', ...
                    'FontSize',fsLab);
            end
        end
    end
end

function tf = should_use_log_loss_scale(y,nIter)
    y = y(isfinite(y) & y > 0);
    if nIter < 5 ...
            || numel(y) < 2
        tf = false;
        return
    end
    tf = (max(y) / min(y)) >= 100;
end

function h = add_history_top_frame(axh,clr)
%ADD_HISTORY_TOP_FRAME Draw solid top frame line for one history axes
% in normalized figure coordinates so it behaves well with yyaxis.

    if ~isgraphics(axh)
        h = gobjects(1);
        return
    end
    fig = ancestor(axh,'figure');
    
    % delete old line if present
    if isappdata(axh,'HistoryTopFrameHandle')
        hOld = getappdata(axh,'HistoryTopFrameHandle');
        try
            if isgraphics(hOld)
                delete(hOld);
            end
        catch
        end
    end
    
    % drawnow limitrate nocallbacks;
    
    pos = axh.Position;   % normalized figure coordinates
    x1 = pos(1);
    x2 = pos(1) + pos(3);
    y  = pos(2) + pos(4);
    
    h = annotation(fig,'line',[x1 x2],[y y], ...
        'Color',clr, ...
        'LineWidth',1, ...
        'LineStyle','-', ...
        'HitTest','off');
    
    setappdata(axh,'HistoryTopFrameHandle',h);

end

function refresh_history_top_frames(fig)
%REFRESH_HISTORY_TOP_FRAMES Redraw top frame lines after figure resize.

    if ~isgraphics(fig)
        return
    end
    
    if ~isappdata(fig,'HistoryAxesHandles') ...
            || ~isappdata(fig,'HistoryTopFrameColor')
        return
    end
    
    axs = getappdata(fig,'HistoryAxesHandles');
    clr = getappdata(fig,'HistoryTopFrameColor');
    
    axs = axs(isgraphics(axs));
    
    for k = 1:numel(axs)
        add_history_top_frame(axs(k),clr);
    end

end

function set_history_xticks(axh,xt,showLabels)
%SET_HISTORY_XTICKS Force unique integer ticks on yyaxis history panels.

    if isempty(xt)
        xt = 1;
    end
    
    xt = unique(round(xt(:).'),'stable');
    
    % Important: x-axis is shared between yyaxis left/right.
    % Do not clear XTickLabel again after switching yyaxis.
    set(axh,'XTick',xt);
    
    if showLabels
        set(axh,'XTickLabel',compose('%d',xt), ...
            'XTickLabelMode','manual', ...
            'XTickLabelRotation',0);
    else
        set(axh,'XTickLabel',[], ...
            'XTickLabelMode','manual', ...
            'XTickLabelRotation',0);
    end
    
    set(axh,'TickDir','out', ...
        'TickLength',[0.012 0.012]);
end

function xt = local_history_xticks(xmax)
%LOCAL_HISTORY_XTICKS Clean iteration ticks for history plots.

    xmax = max(1,round(xmax));
    
    if xmax <= 5
        xt = 1:xmax;
        return
    end

    % At most about six horizontal labels. Select a conventional 1/2/5
    % interval so values remain easy to read (for example, a 250-iteration
    % run uses 1, 50, 100, 150, 200, 250).
    rawStep = xmax/5;
    decade = 10^floor(log10(rawStep));
    scaled = rawStep/decade;
    if scaled <= 1
        multiplier = 1;
    elseif scaled <= 2
        multiplier = 2;
    elseif scaled <= 5
        multiplier = 5;
    else
        multiplier = 10;
    end
    step = max(1,round(multiplier*decade));
    xt = step:step:xmax;

    if ~isempty(xt) && xmax-xt(end) < 0.45*step
        xt(end) = xmax;
    end
    xt = unique([1 xt xmax]);
end

function local_set_dynamic_ylim(axh,side,y)

    yyaxis(axh,side);
    
    y = y(isfinite(y));
    
    if isempty(y)
        set(axh,'YLimMode','auto');
        return
    end
    
    ymin = min(y);
    ymax = max(y);
    
    if ymin == ymax
        pad = max(0.05*abs(ymin),0.05);
    else
        pad = 0.08*(ymax - ymin);
    end
    
    yl = [ymin - pad, ymax + pad];
    
    % Keep nonnegative histories anchored at zero when appropriate
    if ymin >= 0
        yl(1) = max(0,yl(1));
    end
    
    set(axh,'YLim',yl);

end

function local_set_dynamic_ylim_log(axh,side,y)

    yyaxis(axh,side);
    
    y = y(isfinite(y) & y > 0);
    
    if isempty(y)
        set(axh,'YLimMode','auto');
        return
    end
    
    ymin = min(y);
    ymax = max(y);
    
    if ymin == ymax
        fac = 1.25;
        yl = [ymin/fac, ymax*fac];
    else
        yl = [ymin/1.25, ymax*1.25];
    end
    
    set(axh,'YLim',yl);

end

function tf = should_use_log_history_scale(y)
    y = y(isfinite(y) ...
        & y > 0);
    if numel(y) < 2
        tf = false;
        return
    end
    tf = max(y)/min(y) >= 100;
end

function local_no_axis_exponent(axh)
    try
        axh.YAxis(1).Exponent = 0;
    catch
    end
    try
        axh.YAxis(2).Exponent = 0;
    catch
    end
end
