function timeline_helpers(ax,T)
%TIMELINE_HELPERS Draw the SAGE Period-tab split timeline.
%
% This function is intentionally external to SAGE_ui.m. Keeping the
% timeline drawing code outside the large GUI file avoids MATLAB parser
% complexity problems caused by many nested helper functions.

    if nargin < 2 || isempty(T)
        T = struct();
    end

    try
        cla(ax,'reset');
    catch
        cla(ax);
    end
    hold(ax,'on');

    method = lower(string(getfield_safe(T,'method','manual')));
    spinup = getfield_safe(T,'spinup',0);
    spinup = max(0,round(double(spinup)));

    switch method
        case "manual"
            draw_manual(ax,T,spinup);

        case "deterministic_block"
            ds = getfield_safe(T,'dblock_ds',[]);
            de = getfield_safe(T,'dblock_de',[]);
            blockDays = getfield_safe(T,'dblock_size',365);
            ftr = clip01(getfield_safe(T,'dblock_frac',0.8));
            draw_block(ax,ds,de,spinup,blockDays,ftr,false,1, ...
                sprintf(['Deterministic block split: block widths and spin-up ' ...
                'are proportional to days; train fraction = %.2f.'],ftr));

        case "random_block"
            ds = getfield_safe(T,'rblock_ds',[]);
            de = getfield_safe(T,'rblock_de',[]);
            blockDays = getfield_safe(T,'rblock_size',365);
            ftr = clip01(getfield_safe(T,'rblock_frac',0.8));
            seed = round(double(getfield_safe(T,'rblock_seed',1)));
            draw_block(ax,ds,de,spinup,blockDays,ftr,true,seed, ...
                sprintf(['Random block split: block widths and spin-up are ' ...
                'proportional to days; block size = %d d, train fraction = %.2f.'], ...
                max(1,round(blockDays)),ftr));

        case "random"
            ds = getfield_safe(T,'random_ds',[]);
            de = getfield_safe(T,'random_de',[]);
            ftr = clip01(getfield_safe(T,'random_frac',0.8));
            seed = round(double(getfield_safe(T,'random_seed',1)));
            draw_random_points(ax,ds,de,spinup,ftr,seed);

        case "deterministic_kfold"
            ds = getfield_safe(T,'dkfold_ds',[]);
            de = getfield_safe(T,'dkfold_de',[]);
            nFolds = max(2,round(double(getfield_safe(T,'dkfold_n',5))));
            held = max(1,min(nFolds,round(double(getfield_safe(T,'dkfold_held',1)))));
            draw_kfold(ax,ds,de,spinup,nFolds,held,false);

        case "random_kfold"
            ds = getfield_safe(T,'rkfold_ds',[]);
            de = getfield_safe(T,'rkfold_de',[]);
            nFolds = max(2,round(double(getfield_safe(T,'rkfold_n',5))));
            held = max(1,min(nFolds,round(double(getfield_safe(T,'rkfold_held',1)))));
            draw_kfold(ax,ds,de,spinup,nFolds,held,true);

        case "rainfall_block"
            ds = getfield_safe(T,'rain_ds',[]);
            de = getfield_safe(T,'rain_de',[]);
            draw_rainfall_block(ax,ds,de,spinup);

        otherwise
            draw_blocks_days(ax,[0 100],"period",{[0.86 0.86 0.86]}, ...
                'Select a split method.');
    end

    hold(ax,'off');
end

function draw_manual(ax,T,spinupDays)

    dts = getfield_safe(T,'dts',[]);
    dte = getfield_safe(T,'dte',[]);
    des = getfield_safe(T,'des',[]);
    dee = getfield_safe(T,'dee',[]);

    try
        if isempty(dts) ...
                || isempty(dte) ...
                || isempty(des) ...
                || isempty(dee)
            fallback(ax,spinupDays,['Manual dates: ' ...
                'enter valid training/evaluation dates.']);
            return
        end

        tTrain1 = datetime(dts(3),dts(2),dts(1));
        tTrain2 = datetime(dte(3),dte(2),dte(1));
        tEval1  = datetime(des(3),des(2),des(1));
        tEval2  = datetime(dee(3),dee(2),dee(1));

        if tTrain2 < tTrain1 || tEval2 < tEval1
            fallback(ax,spinupDays,['Manual dates: ' ...
                'end date must follow start date.']);
            return
        end

        draw_manual_days(ax,tTrain1, ...
            tTrain2,tEval1,tEval2,spinupDays);
    catch
        fallback(ax,spinupDays,['Manual dates: ' ...
            'enter valid training/evaluation dates.']);
    end
end

function draw_manual_days(ax,tTrain1, ...
    tTrain2,tEval1,tEval2,spinupDays)

    cSpin  = [0.86 0.86 0.86];
    cTrain = [0.55 0.80 0.90];
    cEval  = [1.00 0.72 0.10];
    cGap   = [0.94 0.94 0.94];

    T0 = min([tTrain1 tEval1]);
    T1 = max([tTrain2 tEval2]) + days(1); %#ok
    spinupDays = max(0,round(spinupDays));

    intervals = struct( ...
        't1',{tTrain1,tEval1}, ...
        't2',{tTrain2 + days(1),tEval2 + days(1)}, ...
        'lab',{"training","evaluation"}, ...
        'col',{cTrain,cEval});

    [~,ord] = sort([intervals.t1]);
    intervals = intervals(ord);

    edges = [];
    labels = strings(0);
    colors = {};

    if spinupDays > 0
        edges = [edges 0 spinupDays];
        labels(end+1) = "spin-up"; 
        colors{end+1} = cSpin;
    end

    lastT = T0;
    for ii = 1:numel(intervals)
        if intervals(ii).t1 > lastT
            edges = [edges, ...
                spinupDays + days(lastT - T0), ...
                spinupDays + days(intervals(ii).t1 - T0)]; %#ok<AGROW>
            labels(end+1) = "gap"; %#ok<AGROW>
            colors{end+1} = cGap; %#ok<AGROW>
        end

        edges = [edges, ...
            spinupDays + days(intervals(ii).t1 - T0), ...
            spinupDays + days(intervals(ii).t2 - T0)]; %#ok<AGROW>
        labels(end+1) = intervals(ii).lab; %#ok<AGROW>
        colors{end+1} = intervals(ii).col; %#ok<AGROW>
        lastT = max(lastT,intervals(ii).t2);
    end

    draw_blocks_days(ax, ...
        pair_edges_to_sequence(edges,labels), ...
        labels,colors, ...
        ['Manual dates: blocks are chronological ' ...
        'and proportional to calendar days.']);
end

function draw_block(ax,ds,de, ...
    spinupDays,blockDays,ftr,isRandom,seed,note)

    [ok,~,~,nDays] = timeline_dates(ds,de);
    if ~ok
        fallback(ax,spinupDays,['Block split: ' ...
            'enter valid start/end dates.']);
        return
    end

    spinupDays = max(0,round(spinupDays));
    blockDays = max(1,round(double(blockDays)));
    nBlocks = ceil(nDays/blockDays);
    dur = repmat(blockDays,1,nBlocks);
    dur(end) = nDays - blockDays*(nBlocks-1);

    nTrain = max(0,min(nBlocks, ...
        round(ftr*nBlocks)));
    idTrain = false(1,nBlocks);
    idTrain(1:nTrain) = true;

    if isRandom
        rng(seed);
        idTrain = idTrain(randperm(nBlocks));
    end

    edges = [0 spinupDays spinupDays + cumsum(dur)];
    labels = strings(1,nBlocks+1);
    colors = cell(1,nBlocks+1);
    labels(1) = "spin-up";
    colors{1} = [0.86 0.86 0.86];

    for ii = 1:nBlocks
        if idTrain(ii)
            labels(ii+1) = "train";
            colors{ii+1} = [0.55 0.80 0.90];
        else
            labels(ii+1) = "eval";
            colors{ii+1} = [1.00 0.72 0.10];
        end
    end

    draw_blocks_days(ax, ...
        edges,labels,colors,note);
end

function draw_random_points(ax, ...
    ds,de,spinupDays,ftr,seed)

    [ok,~,~,nDays] = timeline_dates(ds,de);
    if ~ok
        fallback(ax,spinupDays,['Random split: ' ...
            'enter valid start/end dates.']);
        return
    end

    spinupDays = max(0,round(spinupDays));
    n = min(80,max(20,round(nDays/90)));
    x = spinupDays + linspace(0,nDays,n);

    rng(seed);
    nTrain = max(0,min(n,round(ftr*n)));
    idTrain = false(1,n);
    idTrain(1:nTrain) = true;
    idTrain = idTrain(randperm(n));

    draw_blocks_days(ax,[0 spinupDays spinupDays+nDays], ...
        ["spin-up","record"], ...
        {[0.86 0.86 0.86],[0.94 0.94 0.94]}, ...
        sprintf(['Random day split: ' ...
        'dots reflect train fraction = %.2f; ' ...
        'axis widths show days.'],ftr));

    plot(ax,x(idTrain),0.50*ones(1,sum(idTrain)), ...
        'o','MarkerSize',6, ...
        'MarkerFaceColor',[0.55 0.80 0.90], ...
        'MarkerEdgeColor',[0.2 0.2 0.2]);

    plot(ax,x(~idTrain),0.50*ones(1,sum(~idTrain)), ...
        'o','MarkerSize',6, ...
        'MarkerFaceColor',[1.00 0.72 0.10], ...
        'MarkerEdgeColor',[0.2 0.2 0.2]);

end

function draw_kfold(ax,ds, ...
    de,spinupDays,nFolds,held,isRandom)

    [ok,~,~,nDays] = timeline_dates(ds,de);
    if ~ok
        fallback(ax,spinupDays,['K-fold split: ' ...
            'enter valid start/end dates.']);
        return
    end

    spinupDays = max(0,round(spinupDays));
    nFolds = max(2,round(double(nFolds)));
    held = max(1,min(nFolds, ...
        round(double(held))));

    base = floor(nDays/nFolds);
    remd = mod(nDays,nFolds);
    dur = base*ones(1,nFolds);
    dur(1:remd) = dur(1:remd) + 1;

    edges = [0 spinupDays spinupDays ...
        + cumsum(dur)];
    labels = strings(1,nFolds+1);
    colors = cell(1,nFolds+1);
    labels(1) = "spin-up";
    colors{1} = [0.86 0.86 0.86];

    for ii = 1:nFolds
        labels(ii+1) = "fold " + string(ii);
        if ii == held
            colors{ii+1} = [1.00 0.72 0.10];
        else
            colors{ii+1} = [0.55 0.80 0.90];
        end
    end

    if isRandom
        note = sprintf(['Random k-fold split: ' ...
            '%d folds; fold %d is evaluation; ' ...
            'fold widths show days.'], ...
            nFolds,held);
    else
        note = sprintf(['Deterministic k-fold split: ' ...
            '%d folds; fold %d is evaluation; ' ...
            'fold widths show days.'], ...
            nFolds,held);
    end

    draw_blocks_days(ax, ...
        edges,labels,colors,note);
end

function draw_rainfall_block(ax, ...
    ds,de,spinupDays)

    [ok,dsT,deT,nDays] = timeline_dates(ds,de);
    if ~ok
        fallback(ax,spinupDays,['Rainfall-block split: ' ...
            'enter valid start/end dates.']);
        return
    end

    spinupDays = max(0,round(spinupDays));
    cSpin  = [0.86 0.86 0.86];
    cTrain = [0.55 0.80 0.90];
    cEval  = [1.00 0.72 0.10];

    edges = [0 spinupDays];
    labels = "spin-up";
    colors = {cSpin};

    wyStart = dsT;
    k = 1;
    while wyStart <= deT
        wyEndExcl = min(wyStart + ...
            calyears(1),deT + days(1));
        edges(end+1) = spinupDays + ...
            days(wyEndExcl - dsT); %#ok<AGROW>
        labels(end+1) = "WY" + string(k); %#ok<AGROW>
        if mod(k,2) == 1
            colors{end+1} = cTrain; %#ok<AGROW>
        else
            colors{end+1} = cEval; %#ok<AGROW>
        end
        wyStart = wyEndExcl;
        k = k + 1;
    end

    draw_blocks_days(ax,edges,labels,colors, ...
        sprintf(['Rainfall-block split: ' ...
        'spin-up and water-year widths show ' ...
        'calendar days; record = %d d.'], ...
        round(nDays)));
end

function fallback(ax,spinupDays,note)
    spinupDays = max(0,round(spinupDays));
    mainDays = 5*365;
    draw_blocks_days(ax,[0 spinupDays ...
        spinupDays+mainDays], ...
        ["spin-up","period"], ...
        {[0.86 0.86 0.86], ...
        [0.55 0.80 0.90]},note);
end

function draw_blocks_days(ax, ...
    edges,labels,colors,note)

    h  = 0.45;
    y0 = 0.50 - h/2;
    fsBlock = 12;
    fsNote  = 12;

    edges = double(edges(:).');
    edges = edges(isfinite(edges));

    if numel(edges) < 2
        edges = [0 100];
        labels = "period";
        colors = {[0.86 0.86 0.86]};
    end

    xMin = min(edges);
    xMax = max(edges);
    if xMax <= xMin
        xMax = xMin + 1;
    end

    for ii = 1:min(numel(labels),numel(edges)-1)
        x = edges(ii);
        w = edges(ii+1) - edges(ii);
        if w <= 0
            continue
        end

        rectangle(ax, ...
            'Position',[x y0 w h], ...
            'FaceColor',colors{ii}, ...
            'EdgeColor',[0.25 0.25 0.25], ...
            'LineWidth',0.75);

        if w >= 0.045*(xMax-xMin)
            text(ax,x+w/2,y0+h/2,labels(ii), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize',fsBlock, ...
                'Interpreter','none');
        end
    end

    text(ax,xMin+0.01*(xMax-xMin),0.90,note, ...
        'FontSize',fsNote, ...
        'Interpreter','none', ...
        'HorizontalAlignment','left');

    xlim(ax,[xMin xMax]);
    ylim(ax,[0 1]);
    ax.YTick = [];
    ax.FontSize = 12;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.XAxisLocation = 'bottom';
    ax.XColor = [0.15 0.15 0.15];
    ax.YColor = 'none';
    set_day_ticks(ax,xMin,xMax);

end

function set_day_ticks(ax,xMin,xMax)
    span = xMax - xMin;
    if span <= 0
        ax.XTick = [xMin xMax];
        ax.XTickLabel = {'start','end'};
        return
    end

    ticks = unique(round( ...
        linspace(xMin,xMax,5)),'stable');
    labels = strings(size(ticks));
    for ii = 1:numel(ticks)
        if ii == 1
            labels(ii) = "start";
        elseif ii == numel(ticks)
            labels(ii) = "end";
        else
            labels(ii) = sprintf('%d d', ...
                round(ticks(ii)-xMin));
        end
    end

    ax.XTick = ticks;
    ax.XTickLabel = cellstr(labels);
end

function [ok,dsT,deT,nDays] = timeline_dates(ds,de)
    ok = false;
    dsT = NaT;
    deT = NaT;
    nDays = NaN;
    try
        if isempty(ds) || isempty(de)
            return
        end
        dsT = datetime(ds(3),ds(2),ds(1));
        deT = datetime(de(3),de(2),de(1));
        nDays = days(deT - dsT) + 1;
        ok = isfinite(nDays) ...
            && nDays > 0;
    catch
        ok = false;
    end
end

function e = pair_edges_to_sequence(pairEdges,labels)
    e = zeros(1,numel(labels)+1);
    for ii = 1:numel(labels)
        e(ii) = pairEdges(2*ii-1);
        e(ii+1) = pairEdges(2*ii);
    end
end

function f = clip01(f)
    f = double(f);
    if ~isfinite(f)
        f = 0.8;
    end
    f = max(0,min(1,f));
end

function v = getfield_safe(S,name,defaultValue)
    if isstruct(S) ...
            && isfield(S,name) ...
            && ~isempty(S.(name))
        v = S.(name);
    else
        v = defaultValue;
    end
end
