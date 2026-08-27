function draw_FFN_schematic(ax,nIn,nHidden,nOut,tf,selectedOutputLabel)
%DRAW_FFN_SCHEMATIC Draw the shared SAGE FFN schematic
%   ax        target axes
%   nIn       number of input attributes (r)
%   nHidden   scalar/vector with hidden-layer sizes, e.g. 8 or [12 6]
%   nOut      number of outputs (d)
%   tf        cell transfer functions, e.g. {'tanh'} or {'tanh','relu'}
%   selectedOutputLabel optional LaTeX parameter name. When supplied for a
%             single output, that parameter replaces the generic theta_1.
%
% Notes:
%   - shows up to 15 positions in each of as many as five hidden layers
%   - if a layer has more neurons, the plot uses a compact layout with
%     top nodes, a central ellipsis, and bottom nodes
%   - inputs/outputs are also shown compactly to keep diagram readable

    if nargin < 5 ...
            || isempty(tf)
        tf = {'tanh'};
    end
    selectedMode = false;
    selectedOutputIndex = 1;
    if nargin < 6
        selectedOutputLabel = '';
    elseif isstruct(selectedOutputLabel)
        selectedMode = true;
        selectedOutputIndex = double(selectedOutputLabel.index);
        selectedOutputLabel = char(string(selectedOutputLabel.label));
    else
        selectedOutputLabel = char(string(selectedOutputLabel));
        selectedMode = ~isempty(strtrim(selectedOutputLabel)) && nOut == 1;
    end

    nHidden = double(nHidden(:).');
    nLayers = max(1,min(5,numel(nHidden)));
    nHidden = nHidden(1:nLayers);
    nIn = max(1,round(double(nIn)));
    nOut = max(1,round(double(nOut)));

    if numel(tf) < nLayers
        tf(end+1:nLayers) = {'tanh'};
    end

    cla(ax,'reset');
    hold(ax,'on');
    axis(ax,[0 1 0 1]);
    axis(ax,'off');
    ax.Clipping = 'off';

    compactScale = max(0.72,1-0.07*(nLayers-1));
    fsNodeSmall = round(16*compactScale) + 2;
    fsNode = round(16*compactScale) + 2;
    fsEllipsis = round(17*compactScale) + 2;
    fsLayer = round(16*compactScale) + 2;

    % Evenly distribute input, hidden, and output columns. The positions
    % remain vertically fixed when the number of layers changes.
    xLayer = linspace(0.10,0.90,nLayers+2);
    xIn = xLayer(1);
    xHidden = xLayer(2:end-1);
    xOut = xLayer(end);

    % ------------------------
    % display limits per layer
    % ------------------------
    maxShowIn = 8;
    maxShowH = 15;
    maxShowOut = 8;

    % ----------
    % node sizes
    % ----------
    sIn = 0.0105;
    sOut = 0.0105;
    rH = 0.0150*compactScale;

    lineCol = [0.62 0.62 0.62];
    %inputCol = [0.15 0.35 0.85];
    inputCol = [0.70 0.70 0.70];
    hiddenColors = [ ...
        0.90 0.55 0.20; ...
        0.55 0.78 0.28; ...
        0.35 0.70 0.75; ...
        0.56 0.55 0.82; ...
        0.78 0.48 0.72];
    outCol = [0.20 0.20 0.20];

    % % ----------------------------
    % % build displayed node layouts
    % % ----------------------------
    % Lin = local_compact_layer_layout( ...
    %     nIn,maxShowIn);
    % Lout = local_compact_layer_layout( ...
    %     nOut,maxShowOut);
    % LH = cell(1,nLayers);
    % for layer = 1:nLayers
    %     LH{layer} = local_compact_layer_layout( ...
    %         nHidden(layer),maxShowH);
    % end

    % ----------------------------
    % build displayed node layouts
    % ----------------------------

    % Scale the vertical extent of each layer according to its actual width.
    % This makes narrowing architectures taper visually and allows the output
    % layer to widen again when nOut > nHidden(end).
    allLayerSizes = [nIn, nHidden(:).', nOut];
    maxLayerSize = max(allLayerSizes);

    fullSpan = 0.76;     % same overall height as old 0.12 -> 0.88 range
    minSpan  = 0.10;     % prevents very small layers from collapsing visually

    spanIn = minSpan + ...
        (fullSpan-minSpan)*(nIn/maxLayerSize);

    spanOut = minSpan + ...
        (fullSpan-minSpan)*(nOut/maxLayerSize);

    Lin = local_compact_layer_layout( ...
        nIn,maxShowIn,spanIn);

    Lout = local_compact_layer_layout( ...
        nOut,maxShowOut,spanOut);

    if selectedMode
        selectedOutputIndex = min(nOut,max(1,round(selectedOutputIndex)));
        if nOut == 1
            selectedY = 0.5;
        else
            selectedY = 0.5 + spanOut/2 ...
                - (selectedOutputIndex-1)*spanOut/(nOut-1);
        end
        Lout = struct('y',selectedY,'labels',selectedOutputIndex, ...
            'hasGap',false,'gapIndex',NaN,'n',nOut);
    end

    LH = cell(1,nLayers);
    for layer = 1:nLayers

        spanH = minSpan + ...
            (fullSpan-minSpan)*(nHidden(layer)/maxLayerSize);

        LH{layer} = local_compact_layer_layout( ...
            nHidden(layer),maxShowH,spanH);
    end

    drawnow limitrate nocallbacks
    xyScale = local_ffn_axis_xy_scale(ax);
    sInX = sIn*xyScale;
    sOutX = sOut*xyScale;
    rHX = rH*xyScale;

    % ----------------------
    % draw connections first
    % ----------------------
    local_draw_layer_connections(ax, ...
        xIn,local_real_node_y(Lin), ...
        xHidden(1),local_real_node_y(LH{1}), ...
        sInX,rHX,lineCol);
    for layer = 1:nLayers-1
        local_draw_layer_connections(ax, ...
            xHidden(layer),local_real_node_y(LH{layer}), ...
            xHidden(layer+1),local_real_node_y(LH{layer+1}), ...
            rHX,rHX,lineCol);
    end
    local_draw_layer_connections(ax, ...
        xHidden(end),local_real_node_y(LH{end}), ...
        xOut,local_real_node_y(Lout),rHX,sOutX,lineCol, ...
        'FFNSelectedOutputConnection');

    % -------------------------
    % draw input nodes + labels
    % -------------------------
    for k = 1:numel(Lin.y)
        if ~isfinite(Lin.labels(k))
            continue
        end
        local_draw_square(ax, ...
            xIn,Lin.y(k),sIn,inputCol,lineCol);
    end
    local_draw_compact_labels(ax,xIn,Lin,'input', ...
        fsNode,fsNodeSmall,fsEllipsis);

    % ------------------
    % draw hidden layers
    % ------------------
    for layer = 1:nLayers
        for k = 1:numel(LH{layer}.y)
            if ~isfinite(LH{layer}.labels(k))
                continue
            end
            local_draw_circle(ax, ...
                xHidden(layer),LH{layer}.y(k),rH, ...
                hiddenColors(layer,:),lineCol);
        end
        local_draw_compact_labels(ax, ...
            xHidden(layer),LH{layer}, ...
            sprintf('hidden%d',layer),fsNode, ...
            fsNodeSmall,fsEllipsis);
        text(ax,xHidden(layer),0.94, ...
            sprintf('$%s$',local_tf_latex(tf{layer})), ...
            'interpreter','latex', ...
            'horizontalalignment','center', ...
            'verticalalignment','middle', ...
            'fontsize',fsLayer);
    end

    % --------------------------
    % draw output nodes + labels
    % --------------------------
    for k = 1:numel(Lout.y)
        if ~isfinite(Lout.labels(k))
            continue
        end
        hOut = local_draw_square(ax, ...
            xOut,Lout.y(k),sOut,outCol,lineCol);
        if selectedMode
            hOut.Tag = 'FFNSelectedOutputNode';
        end
    end
    if selectedMode
        hOutputLabel = text(ax,xOut+0.0175,Lout.y(1), ...
            ['$' local_selected_output_latex(selectedOutputLabel) '$'], ...
            'interpreter','latex', ...
            'horizontalalignment','left', ...
            'verticalalignment','middle', ...
            'fontsize',fsNode);
        hOutputLabel.Tag = 'FFNSelectedOutputLabel';
    else
        local_draw_compact_labels(ax, ...
            xOut,Lout,'output', ...
            fsNode,fsNodeSmall,fsEllipsis);
    end

    % --------------
    % layer captions
    % --------------
    text(ax,xIn,0.045,sprintf(['inputs: ' ...
        '$r=%d$'],nIn), ...
        'interpreter','latex', ...
        'horizontalalignment','center', ...
        'verticalalignment','middle', ...
        'fontsize',fsLayer);
    for layer = 1:nLayers
        text(ax,xHidden(layer),0.045, ...
            sprintf('$h_{%d}=%d$',layer,nHidden(layer)), ...
            'interpreter','latex', ...
            'horizontalalignment','center', ...
            'verticalalignment','middle', ...
            'fontsize',fsLayer);
    end
    if selectedMode
        outputCaption = 'selected output';
    else
        outputCaption = sprintf('outputs: $d=%d$',nOut);
    end
    text(ax,xOut,0.045,outputCaption, ...
        'interpreter','latex', ...
        'horizontalalignment','center', ...
        'verticalalignment','middle', ...
        'fontsize',fsLayer);

    hold(ax,'off');
end

function s = local_selected_output_latex(name)
%LOCAL_SELECTED_OUTPUT_LATEX Underline the parameter symbol, not its subscript.
    s = strtrim(char(string(name)));
    s = regexprep(s,'^\$|\$$','');
    if startsWith(s,'\underline{')
        return
    end
    parts = regexp(s,'^([^_]+)(.*)$','tokens','once');
    if isempty(parts)
        s = ['\underline{' s '}'];
    else
        s = ['\underline{' parts{1} '}' parts{2}];
    end
end

function L = local_compact_layer_layout(n,maxShow,span)
%function L = local_compact_layer_layout(n,maxShow)
%LOCAL_COMPACT_LAYER_LAYOUT Return y-positions and displayed ids
%
% Output fields:
%   .y          y positions of displayed markers
%   .labels     integer ids for displayed nodes
%   .hasGap     true if compressed with ellipsis
%   .gapIndex   index where ellipsis should be drawn

    n = max(1,round(double(n)));
    maxShow = max(3,round(double(maxShow)));

    % added with new function
    if nargin < 3 ...
            || isempty(span)
        span = 0.76;
    end
    span = max(0.08,min(0.76,double(span)));
    yMin = 0.5 - span/2;
    yMax = 0.5 + span/2;
    % end added with new function

    % if n <= maxShow
    %     y = linspace(0.12,0.88,n);
    %     labels = 1:n;
    %     hasGap = false;
    %     gapIndex = NaN;
    % else
    %     nTop = floor((maxShow-1)/2);
    %     nBot = maxShow - nTop - 1;
    %     yTop = linspace(0.82,0.56,nTop);
    %     yMid = 0.50;
    %     yBot = linspace(0.44,0.18,nBot);
    %     y = [yTop, yMid, yBot];
    %     labels = [1:nTop, NaN, (n-nBot+1):n];
    %     hasGap = true;
    %     gapIndex = nTop + 1;
    % end

    if n <= maxShow

        if n == 1
            y = 0.5;
        else
            y = linspace(yMax,yMin,n);
        end

        labels = 1:n;
        hasGap = false;
        gapIndex = NaN;

    else

        nTop = floor((maxShow-1)/2);
        nBot = maxShow - nTop - 1;

        % Keep the ellipsis centered while shrinking/expanding the entire
        % displayed layer according to its true number of neurons.
        gapHalf = 0.08*span/0.76;

        upperInner = 0.5 + gapHalf;
        lowerInner = 0.5 - gapHalf;

        yTop = linspace(yMax,upperInner,nTop);
        yMid = 0.50;
        yBot = linspace(lowerInner,yMin,nBot);

        y = [yTop, yMid, yBot];
        labels = [1:nTop, NaN, (n-nBot+1):n];
        hasGap = true;
        gapIndex = nTop + 1;
    end

    L = struct();
    L.y = y(:).';
    L.labels = labels(:).';
    L.hasGap = hasGap;
    L.gapIndex = gapIndex;
    L.n = n;
end

function local_draw_layer_connections(ax, ...
        x1,y1,x2,y2,r1,r2,col,varargin)  %#ok

    tag = '';
    if ~isempty(varargin), tag = char(string(varargin{1})); end

    for i = 1:numel(y1)
        for j = 1:numel(y2)
            %            plot(ax,[x1+r1 x2-r2],[y1(i) y2(j)],'-', ...
            h = plot(ax,[x1 x2],[y1(i) y2(j)],'-', ...
                'color',col, ...
                'linewidth',0.30);
            if ~isempty(tag), h.Tag = tag; end
        end
    end
end

function y = local_real_node_y(layout)
    y = layout.y(isfinite(layout.labels));
end

function scale = local_ffn_axis_xy_scale(ax)
    scale = 1;
    try
        position = getpixelposition(ax,true);
        if numel(position) >= 4 ...
                && position(3) > 0 ...
                && position(4) > 0
            scale = position(4)/position(3);
        end
    catch
    end
    scale = min(10,max(0.1,scale));
end

function local_draw_compact_labels(ax, ...
        x,L,kind,fsNode,fsNodeSmall,fsEllipsis)

    switch lower(kind)
        case 'input'
            xTxt = x - 0.035;
            ha = 'right';
            if L.n == 1
                text(ax,xTxt,L.y(1), ...
                    '$a_{1}$', ...
                    'interpreter','latex', ...
                    'horizontalalignment',ha, ...
                    'verticalalignment','middle', ...
                    'fontsize',fsNode);
                return
            end

            if ~L.hasGap
                for k = 1:numel(L.labels)
                    lab = sprintf('$a_{%d}$', ...
                        L.labels(k));
                    text(ax,xTxt,L.y(k),lab, ...
                        'interpreter','latex', ...
                        'horizontalalignment',ha, ...
                        'verticalalignment','middle', ...
                        'fontsize',fsNodeSmall);
                end
            else
                nShow = numel(L.labels);
                kTopLast = L.gapIndex - 1;
                kBotFirst = L.gapIndex + 1;

                text(ax,xTxt,L.y(1),'$a_{1}$', ...
                    'interpreter','latex', ...
                    'horizontalalignment',ha, ...
                    'verticalalignment','middle', ...
                    'fontsize',fsNode);

                if kTopLast >= 2
                    text(ax,xTxt,L.y(2), ...
                        sprintf('$a_{%d}$', ...
                        L.labels(2)), ...
                        'interpreter','latex', ...
                        'horizontalalignment',ha, ...
                        'verticalalignment','middle', ...
                        'fontsize',fsNodeSmall);
                end

                text(ax,xTxt,L.y(L.gapIndex), ...
                    sprintf('%c',hex2dec('22EE')), ...
                    'interpreter','none', ...
                    'horizontalalignment',ha, ...
                    'verticalalignment','middle', ...
                    'fontsize',fsEllipsis);

                if kBotFirst <= nShow-1
                    text(ax,xTxt,L.y(kBotFirst), ...
                        sprintf('$a_{%d}$', ...
                        L.labels(kBotFirst)), ...
                        'interpreter','latex', ...
                        'horizontalalignment',ha, ...
                        'verticalalignment','middle', ...
                        'fontsize',fsNodeSmall);
                end

                text(ax,xTxt,L.y(end), ...
                    sprintf('$a_{%d}$', ...
                    L.labels(end)), ...
                    'interpreter','latex', ...
                    'horizontalalignment',ha, ...
                    'verticalalignment','middle', ...
                    'fontsize',fsNode);
            end

        case 'output'
            xTxt = x + 0.035;
            ha = 'left';

            if L.n == 1
                text(ax,xTxt,L.y(1), ...
                    '$\underline{\theta}_{1}$', ...
                    'interpreter','latex', ...
                    'horizontalalignment',ha, ...
                    'verticalalignment','middle', ...
                    'fontsize',fsNode);
                return
            end

            if ~L.hasGap
                for k = 1:numel(L.labels)
                    lab = ...
                        sprintf(['$\\underline{\\theta}_' ...
                        '{%d}$'],L.labels(k));
                    text(ax,xTxt,L.y(k),lab, ...
                        'interpreter','latex', ...
                        'horizontalalignment',ha, ...
                        'verticalalignment','middle', ...
                        'fontsize',fsNodeSmall);
                end
            else
                nShow = numel(L.labels);
                kTopLast = L.gapIndex - 1;
                kBotFirst = L.gapIndex + 1;

                text(ax,xTxt,L.y(1), ...
                    '$\underline{\theta}_{1}$', ...
                    'interpreter','latex', ...
                    'horizontalalignment',ha, ...
                    'verticalalignment','middle', ...
                    'fontsize',fsNode);

                if kTopLast >= 2
                    text(ax,xTxt,L.y(2), ...
                        sprintf(['$\\underline{\\theta}_' ...
                        '{%d}$'],L.labels(2)), ...
                        'interpreter','latex', ...
                        'horizontalalignment',ha, ...
                        'verticalalignment','middle', ...
                        'fontsize',fsNodeSmall);
                end

                text(ax,xTxt,L.y(L.gapIndex), ...
                    sprintf('%c',hex2dec('22EE')), ...
                    'interpreter','none', ...
                    'horizontalalignment',ha, ...
                    'verticalalignment','middle', ...
                    'fontsize',fsEllipsis);

                if kBotFirst <= nShow-1
                    text(ax,xTxt,L.y(kBotFirst), ...
                        sprintf(['$\\underline{\\theta}' ...
                        '_{%d}$'],L.labels(kBotFirst)), ...
                        'interpreter','latex', ...
                        'horizontalalignment',ha, ...
                        'verticalalignment','middle', ...
                        'fontsize',fsNodeSmall);
                end

                text(ax,xTxt,L.y(end), ...
                    sprintf(['$\\underline{\\theta}' ...
                    '_{%d}$'],L.labels(end)), ...
                    'interpreter','latex', ...
                    'horizontalalignment',ha, ...
                    'verticalalignment','middle', ...
                    'fontsize',fsNode);
            end

        case {'hidden1','hidden2','hidden3','hidden4','hidden5'}
            if L.hasGap
                text(ax,x, ...
                    L.y(L.gapIndex), ...
                    sprintf('%c',hex2dec('22EE')), ...
                    'interpreter','none', ...
                    'horizontalalignment','center', ...
                    'verticalalignment','middle', ...
                    'fontsize',fsEllipsis, ...
                    'BackgroundColor','w', ...
                    'Margin',0.5);
            end
    end
end

function s = local_tf_latex(tf)

    tf = lower(strtrim(char(string(tf))));

    switch tf
        case 'tanh'
            s = '{\rm tanh}';
        case 'relu'
            s = '{\rm ReLU}';
        otherwise
            s = ['{\rm ' tf '}'];
    end
end

function hp = local_draw_square(ax,xc,yc,h,fc,ec)
    hx = h*local_ffn_axis_xy_scale(ax);
    hp = patch(ax,[xc-hx xc+hx xc+hx xc-hx], ...
        [yc-h yc-h yc+h yc+h],fc, ...
        'edgecolor',ec, ...
        'linewidth',0.7);
end

function local_draw_circle(ax,xc,yc,r,fc,ec)
    th = linspace(0,2*pi,80);
    rx = r*local_ffn_axis_xy_scale(ax);
    patch(ax,xc + rx*cos(th),yc + r*sin(th),fc, ...
        'edgecolor',ec, ...
        'linewidth',0.7);
end
