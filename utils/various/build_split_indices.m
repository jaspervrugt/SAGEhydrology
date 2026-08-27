function split = build_split_indices(split)
%BUILD_SPLIT_INDICES Build global model/data indices from split only
%
% Notes:
%   split.tout      model print times over full simulation window
%   split.idx       state/print indices of scoring window within tout
%   split.id_train  data indices of training period within scored arrays
%   split.id_eval   data indices of evaluation period within scored arrays

    dt0 = split.dt0;
    dtTrain0 = split.dt_train0;
    dtTrain1 = split.dt_train1;
    dtEval0 = split.dt_eval0;
    dtEval1 = split.dt_eval1;
    dtEnd = split.dt_end;
    dt = double(split.dt);
    
    switch dt
        case 1
            n0 = int64(days(dateshift( ...
                dtTrain0,'start','day') ...
                 - dateshift(dt0,'start','day')));
            n1 = int64(days(dateshift( ...
                dtTrain1,'start','day') ...
                 - dateshift(dt0,'start','day')));
            m0 = int64(days(dateshift( ...
                dtEval0,'start','day') ...
                 - dateshift(dt0,'start','day')));
            m1 = int64(days(dateshift( ...
                dtEval1,'start','day') ...
                 - dateshift(dt0,'start','day')));
            nEnd = int64(days(dateshift( ...
                dtEnd,'start','day') ...
                   - dateshift(dt0,'start','day')));
    
        case 24
            t0 = datetime(year(dt0), ...
                month(dt0),day(dt0),0,0,0);
            tTrain0 = datetime(year(dtTrain0), ...
                month(dtTrain0), ...
                day(dtTrain0),0,0,0);
            tTrain1 = datetime(year(dtTrain1), ...
                month(dtTrain1), ...
                day(dtTrain1),23,0,0);
            tEval0 = datetime(year(dtEval0), ...
                month(dtEval0), ...
                day(dtEval0),0,0,0);
            tEval1 = datetime(year(dtEval1), ...
                month(dtEval1), ...
                day(dtEval1),23,0,0);
            tEnd = datetime(year(dtEnd), ...
                month(dtEnd), ...
                day(dtEnd),23,0,0);
    
            n0 = int64(round(hours(tTrain0 - t0)));
            n1 = int64(round(hours(tTrain1 - t0)));
            m0 = int64(round(hours(tEval0 - t0)));
            m1 = int64(round(hours(tEval1 - t0)));
            nEnd = int64(round(hours(tEnd - t0)));
    
        case 96
            % Fifteen-minute records use complete calendar days:
            % 00:00 through 23:45.
            t0 = datetime(year(dt0), ...
                month(dt0),day(dt0),0,0,0);
            tTrain0 = datetime(year(dtTrain0), ...
                month(dtTrain0), ...
                day(dtTrain0),0,0,0);
            tTrain1 = datetime(year(dtTrain1), ...
                month(dtTrain1), ...
                day(dtTrain1),23,45,0);
            tEval0 = datetime(year(dtEval0), ...
                month(dtEval0), ...
                day(dtEval0),0,0,0);
            tEval1 = datetime(year(dtEval1), ...
                month(dtEval1), ...
                day(dtEval1),23,45,0);
            tEnd = datetime(year(dtEnd), ...
                month(dtEnd), ...
                day(dtEnd),23,45,0);
    
            n0 = int64(round(minutes(tTrain0 - t0)/15));
            n1 = int64(round(minutes(tTrain1 - t0)/15));
            m0 = int64(round(minutes(tEval0 - t0)/15));
            m1 = int64(round(minutes(tEval1 - t0)/15));
            nEnd = int64(round(minutes(tEnd - t0)/15));
    
        otherwise
            error(['      Error:build_split_indices: ' ...
                'Unsupported dt = %g. Expected 1, 24, or 96.'],dt);
    end
    
    % ------------------------------------------------------
    % Model print times over full simulation window
    %   If there are N data steps, there are N+1 print times
    % ------------------------------------------------------
    split.tout = 0:double(nEnd + 1);
    
    % ------------------------------------------------
    % Scoring window in model/state coordinates
    %   Starts at earliest of train/eval
    %   Ends at final state after latest of train/eval
    % ------------------------------------------------
    i0 = min(n0,m0);
    i1 = max(n1,m1);
    
    split.idx = (double(i0) + 1):(double(i1) + 2);
    
    % ------------------------------------------------------
    % Training/evaluation indices in DATA coordinates
    %   These index Q, P, Ep, T, etc. over the scored window
    % ------------------------------------------------------
    split.id_train = (double(n0 - i0) + 1):(double(n1 - i0) + 1);
    split.id_eval = (double(m0 - i0) + 1):(double(m1 - i0) + 1);

end