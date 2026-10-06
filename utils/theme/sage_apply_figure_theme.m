function sage_apply_figure_theme(fig,theme,profile)
%SAGE_APPLY_FIGURE_THEME Apply canonical styling to a complete figure.
if nargin < 2 || isempty(theme), theme = sage_visual_theme(); end
if nargin < 3 || isempty(profile), profile = 'standalone'; end
fig = fig(isgraphics(fig));
for k = 1:numel(fig)
    f = fig(k);
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
end
end

function tf = local_is_legend_axes(ax)
tf = false;
try
    tf = strcmpi(ax.Tag,'legend') || strcmpi(ax.Tag,'Colorbar');
catch
end
end
