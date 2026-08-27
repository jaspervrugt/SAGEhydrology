function dat = rainfall_rank_split(dat,split)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RAINFALL_RANK_SPLIT Basin-specific rainfall-ranked train/eval partition
%
% SYNOPSIS: dat = rainfall_rank_split(dat,split)
%
% DESCRIPTION:
%   Creates basin-specific training and evaluation indices using annual
%   water-year precipitation totals. For each basin independently, water
%   years are ranked from wettest to driest and alternately assigned to
%   training and evaluation:
%
%       Wettest WY   -> training
%       2nd wettest  -> evaluation
%       3rd wettest  -> training
%       etc.
%
%   This approach creates rainfall-balanced train/eval partitions while
%   preserving full temporal continuity within each water year.
%
% WATER YEAR:
%   Water year XX+1 spans:
%
%       01-Oct-XX  through  30-Sep-XX+1
%
% INPUT:
%   dat{k}.meteo.P   precipitation time series
%   split            SAGE split structure
%
% OUTPUT:
%   dat{k}.id_train      training indices
%   dat{k}.id_eval       evaluation indices
%   dat{k}.mask_train    logical training mask
%   dat{k}.mask_eval     logical evaluation mask
%   dat{k}.train_WY      training water years
%   dat{k}.eval_WY       evaluation water years
%   dat{k}.P_WY          table with annual precipitation totals
%
% NOTES:
%   - Splitting is performed independently for each basin.
%   - Spin-up periods are excluded through split.idx.
%   - Supports daily, hourly, and 15-minute time steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Jasper A. Vrugt, Mar. 2026 / updated Jul. 2026               %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K = numel(dat);
    dt = double(split.dt);
    
    if ~isfield(split,'idx') ...
            || isempty(split.idx)
        error(['      Error:read_meteo_US: ' ...
            'split.idx is required for ' ...
            'local rainfall-rank split.']);
    end
    
    idx = double(split.idx(:));
    
    % local scored indices inside dat{k}.meteo.P
    % Usually split.idx = spinup+1 : end, e.g. 366:7306
    i0 = idx(1);
    i1 = idx(end);
    
    % Build datetime vector for the scored period
    if isfield(split,'ds') ...
            && ~isempty(split.ds)
        ymd0 = split.ds;
    elseif isfield(split,'dts') ...
            && ~isempty(split.dts)
        ymd0 = split.dts;
    else
        error(['      Error:read_meteo_US: ' ...
            'split.ds or split.dts ' ...
            'is required for local split.']);
    end
    
    t0 = datetime(ymd0(1),ymd0(2),ymd0(3),0,0,0);
    
    for k = 1:K
    
        if ~isfield(dat{k},'meteo') ...
                || ~isfield(dat{k}.meteo,'P') ...
                || isempty(dat{k}.meteo.P)
            error(['      Error:read_meteo_US: ' ...
                'dat{%d}.meteo.P is missing.'],k);
        end
    
        Pfull = double(dat{k}.meteo.P(:));
        
        i1k = min(i1,numel(Pfull));
    
        if i0 > numel(Pfull)
            error(['      Error:read_meteo_US: ' ...
                'Cannot locate scored period ' ...
                'in dat{%d}.meteo.P. ' ...
                'idx(1) = %d, length(P) = %d.'], ...
                k,i0,numel(Pfull));
        end
        
        P = Pfull(i0:i1k);
        n_scored_k = numel(P);
        
        switch dt
            case 1
                tt = t0 + days(0:n_scored_k-1);
            case 24
                tt = t0 + hours(0:n_scored_k-1);
            case 96
                tt = t0 + minutes(15*(0:n_scored_k-1));
            otherwise
                error('rainfall_rank_split:unsupportedDt', ...
                    ['      Error:rainfall_rank_split: ' ...
                    'Unsupported split.dt = %g. Expected 1 ' ...
                    '(daily), 24 (hourly), or 96 (15-minute).'],dt);
        end
        
        tt = tt(:);
    
        wy = year(tt);
        wy(month(tt) >= 10) ...
            = wy(month(tt) >= 10) + 1;
        waterYears = unique(wy,'stable');
        ny = numel(waterYears);
    
        Pyear = nan(ny,1);
        for j = 1:ny
            I = wy == waterYears(j);
            Pyear(j) = sum(P(I),'omitnan');
        end
        
        [~,ord] = sort(Pyear,'descend', ...
            'MissingPlacement','last');
        
        trainYears = waterYears(ord(1:2:end));
        evalYears = waterYears(ord(2:2:end));
        
        id_train = find(ismember(wy,trainYears));
        id_eval = find(ismember(wy,evalYears));
    
        dat{k}.id_train = id_train(:).';
        dat{k}.id_eval = id_eval(:).';
    
        dat{k}.mask_train = false(1,n_scored_k);
        dat{k}.mask_eval = false(1,n_scored_k);
    
        dat{k}.mask_train(dat{k}.id_train) = true;
        dat{k}.mask_eval(dat{k}.id_eval) = true;
    
        dat{k}.train_WY = trainYears(:);
        dat{k}.eval_WY = evalYears(:);
    
        dat{k}.P_WY = table(waterYears(:),Pyear(:), ...
            'VariableNames',{'WaterYear','P_total'});
    end

end

% function dat = rainfall_rank_split(dat,split)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %RAINFALL_RANK_SPLIT Basin-specific rainfall-ranked train/eval partition
% %
% % SYNOPSIS: dat = rainfall_rank_split(dat,split)
% %
% % DESCRIPTION:
% %   Creates basin-specific training and evaluation indices using annual
% %   water-year precipitation totals. For each basin independently, water
% %   years are ranked from wettest to driest and alternately assigned to
% %   training and evaluation:
% %
% %       Wettest WY   -> training
% %       2nd wettest  -> evaluation
% %       3rd wettest  -> training
% %       etc.
% %
% %   This approach creates rainfall-balanced train/eval partitions while
% %   preserving full temporal continuity within each water year.
% %
% % WATER YEAR:
% %   Water year XX+1 spans:
% %
% %       01-Oct-XX  through  30-Sep-XX+1
% %
% % INPUT:
% %   dat{k}.meteo.P   precipitation time series
% %   split            SAGE split structure
% %
% % OUTPUT:
% %   dat{k}.id_train      training indices
% %   dat{k}.id_eval       evaluation indices
% %   dat{k}.mask_train    logical training mask
% %   dat{k}.mask_eval     logical evaluation mask
% %   dat{k}.train_WY      training water years
% %   dat{k}.eval_WY       evaluation water years
% %   dat{k}.P_WY          table with annual precipitation totals
% %
% % NOTES:
% %   - Splitting is performed independently for each basin.
% %   - Spin-up periods are excluded through split.idx.
% %   - Supports daily and hourly time steps.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Written by Jasper A. Vrugt, Mar. 2026                                   %
% % University of California, Irvine                                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% K = numel(dat);
% dt = double(split.dt);
% 
% if ~isfield(split,'idx') ...
%         || isempty(split.idx)
%     error(['      Error:read_meteo_US: ' ...
%         'split.idx is required for ' ...
%         'local rainfall-rank split.']);
% end
% 
% idx = double(split.idx(:));
% 
% % local scored indices inside dat{k}.meteo.P
% % Usually split.idx = spinup+1 : end, e.g. 366:7306
% i0 = idx(1);
% i1 = idx(end);
% 
% % Build datetime vector for the scored period
% if isfield(split,'ds') ...
%         && ~isempty(split.ds)
%     ymd0 = split.ds;
% elseif isfield(split,'dts') ...
%         && ~isempty(split.dts)
%     ymd0 = split.dts;
% else
%     error(['      Error:read_meteo_US: ' ...
%         'split.ds or split.dts ' ...
%         'is required for local split.']);
% end
% 
% t0 = datetime(ymd0(1),ymd0(2),ymd0(3),0,0,0);
% 
% for k = 1:K
% 
%     if ~isfield(dat{k},'meteo') ...
%             || ~isfield(dat{k}.meteo,'P') ...
%             || isempty(dat{k}.meteo.P)
%         error(['      Error:read_meteo_US: ' ...
%             'dat{%d}.meteo.P is missing.'],k);
%     end
% 
%     Pfull = double(dat{k}.meteo.P(:));
% 
%     i1k = min(i1,numel(Pfull));
% 
%     if i0 > numel(Pfull)
%         error(['      Error:read_meteo_US: ' ...
%             'Cannot locate scored period ' ...
%             'in dat{%d}.meteo.P. ' ...
%             'idx(1) = %d, length(P) = %d.'], ...
%             k,i0,numel(Pfull));
%     end
% 
%     P = Pfull(i0:i1k);
%     n_scored_k = numel(P);
% 
%     switch dt
%         case 1
%             tt = t0 + days(0:n_scored_k-1);
%         case 24
%             tt = t0 + hours(0:n_scored_k-1);
%     end
% 
%     tt = tt(:);
% 
%     wy = year(tt);
%     wy(month(tt) >= 10) ...
%         = wy(month(tt) >= 10) + 1;
%     waterYears = unique(wy,'stable');
%     ny = numel(waterYears);
% 
%     Pyear = nan(ny,1);
%     for j = 1:ny
%         I = wy == waterYears(j);
%         Pyear(j) = sum(P(I),'omitnan');
%     end
% 
%     [~,ord] = sort(Pyear,'descend', ...
%         'MissingPlacement','last');
% 
%     trainYears = waterYears(ord(1:2:end));
%     evalYears = waterYears(ord(2:2:end));
% 
%     id_train = find(ismember(wy,trainYears));
%     id_eval = find(ismember(wy,evalYears));
% 
%     dat{k}.id_train = id_train(:).';
%     dat{k}.id_eval = id_eval(:).';
% 
%     dat{k}.mask_train = false(1,n_scored_k);
%     dat{k}.mask_eval = false(1,n_scored_k);
% 
%     dat{k}.mask_train(dat{k}.id_train) = true;
%     dat{k}.mask_eval(dat{k}.id_eval) = true;
% 
%     dat{k}.train_WY = trainYears(:);
%     dat{k}.eval_WY = evalYears(:);
% 
%     dat{k}.P_WY = table(waterYears(:),Pyear(:), ...
%         'VariableNames',{'WaterYear','P_total'});
% end
% 
% end