function sage_apply_figure_theme(fig,theme,profile)
%SAGE_APPLY_FIGURE_THEME Apply canonical styling to a complete figure.
if nargin < 2 || isempty(theme), theme = sage_visual_theme(); end
if nargin < 3 || isempty(profile), profile = 'standalone'; end
fig = fig(isgraphics(fig));
for k = 1:numel(fig)
    f = fig(k);
    if isappdata(f,'SAGEPreservePublicationLayout') ...
            && getappdata(f,'SAGEPreservePublicationLayout')
        continue
    end
    name = char(lower(string(profile)));
    if isfield(theme.profiles,name)
        backgroundName = theme.profiles.(name).background;
    else
        backgroundName = 'standalone';
    end
    if isfield(theme.backgrounds,backgroundName)
        try
            f.Color = theme.backgrounds.(backgroundName);
        catch
        end
    end
    axesObjects = findall(f,'Type','axes');
    axesObjects = axesObjects(~arrayfun(@local_is_legend_axes,axesObjects));
    sage_apply_axes_theme(axesObjects,theme,profile);
    legends = findall(f,'Type','legend');
    sage_apply_legend_theme(legends,theme,profile);
    bars = findall(f,'Type','colorbar');
    sage_apply_colorbar_theme(bars,theme,profile);
    textObjects = findall(f,'Type','text');
    excluded = gobjects(0);
    for j = 1:numel(axesObjects)
        excluded(end+1:end+4) = [axesObjects(j).Title; ...
            axesObjects(j).XLabel; axesObjects(j).YLabel; ...
            axesObjects(j).ZLabel];
    end
    textObjects = textObjects(~ismember(textObjects,excluded));
    % Preserve semantic size differences such as ECDF median labels while
    % removing platform-dependent and legacy font-family assignments.
    for j = 1:numel(textObjects)
        try
            textObjects(j).FontName = theme.plots.fontName;
        catch
        end
    end
    local_publication_typography(f,axesObjects,theme);
end
end

function local_publication_typography(fig,axesObjects,theme)
% Preserve publication sizing and short ticks in standalone and PPTX output.
if ~isappdata(fig,'SAGEFigureFamily'), return; end
family = char(string(getappdata(fig,'SAGEFigureFamily')));
switch family
    case 'variogram'
        tickSize=11; labelSize=10; titleSize=12; tickLength=[0.012 0.012];
    case 'fdc'
        tickSize=13; labelSize=15; titleSize=16; tickLength=[0.010 0.010];
    case 'timeseries'
        tickSize=13; labelSize=14; titleSize=16; tickLength=[0.0015 0.0015];
    case 'zone_ecdf'
        tickSize=10; labelSize=10; titleSize=11; tickLength=[0.012 0.012];
    case 'scenario_ecdf'
        tickSize=12; labelSize=13; titleSize=14; tickLength=[0.030 0.030];
    otherwise
        return
end
for j=1:numel(axesObjects)
    ax=axesObjects(j);
    set(ax,'FontName',theme.plots.fontName,'FontSize',tickSize, ...
        'TickLength',tickLength);
    ax.XLabel.FontSize=labelSize;
    ax.YLabel.FontSize=labelSize;
    ax.Title.FontSize=titleSize;
    if strcmp(family,'fdc')
        ax.TickLabelInterpreter='latex';
    end
    if strcmp(family,'zone_ecdf')
        % These labels are constructed as plain text and TeX, not LaTeX.
        ax.XLabel.Interpreter='tex';
        ax.YLabel.Interpreter='none';
        ax.Title.Interpreter='none';
    end
end
if strcmp(family,'timeseries')
    set(findall(fig,'Type','text','Tag','SAGEFixedLeftYLabel'),'FontSize',14);
end
if strcmp(family,'fdc')
    set(findall(fig,'Type','text','Tag','SAGEFDCBasinLabel'), ...
        'FontSize',14,'Interpreter','none');
end
end

function tf = local_is_legend_axes(ax)
tf = false;
try
    tf = strcmpi(ax.Tag,'legend') || strcmpi(ax.Tag,'Colorbar');
catch
end
end
