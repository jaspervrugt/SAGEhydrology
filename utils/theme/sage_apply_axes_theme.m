function sage_apply_axes_theme(ax,theme,profile)
%SAGE_APPLY_AXES_THEME Apply the canonical scientific axes contract.

if nargin < 2 || isempty(theme), theme = sage_visual_theme(); end
if nargin < 3 || isempty(profile), profile = 'gui'; end
p = local_profile(theme,profile);
if isempty(ax), return; end
ax = ax(isgraphics(ax));
for k = 1:numel(ax)
    a = ax(k);
    local_set(a,'FontName',theme.plots.fontName);
    local_set(a,'FontSize',round(theme.plots.tickSize*p.scale));
    local_set(a,'FontWeight','normal');
    local_set(a,'TickDir',theme.axes.tickDirection);
    local_set(a,'TickLength',theme.axes.tickLength);
    local_set(a,'LineWidth',theme.plots.axisLineWidth);
    local_set(a,'XColor',theme.plots.axisColor);
    local_set(a,'YColor',theme.plots.axisColor);
    local_set(a,'GridColor',theme.plots.gridColor);
    local_set(a,'GridAlpha',theme.plots.gridAlpha);
    local_set(a,'MinorGridColor',theme.plots.gridColor);
    local_set(a,'MinorGridAlpha',theme.plots.gridAlpha/2);
    local_set(a,'Box',p.box);
    local_set(a,'Layer',p.layer);
    local_set(a,'TickLabelInterpreter',theme.interpreters.tick);
    local_set(a,'Color',local_background(theme,p.background));
    local_text(a.Title,theme,'title',p.scale,[]);
    local_text(a.XLabel,theme,'label',p.scale,[]);
    local_text(a.YLabel,theme,'label',p.scale,[]);
    if isprop(a,'ZLabel')
        local_text(a.ZLabel,theme,'label',p.scale,[]);
    end
end
end

function p = local_profile(theme,profile)
if isstruct(profile), p = profile; return; end
name = char(lower(string(profile)));
if isfield(theme.profiles,name), p = theme.profiles.(name);
else, p = theme.profiles.gui;
end
end

function color = local_background(theme,name)
if isnumeric(name), color = name; return; end
name = char(string(name));
if isfield(theme.backgrounds,name), color = theme.backgrounds.(name);
else, color = theme.axes.background;
end
end

function local_text(h,theme,role,scale,interpreter)
if isempty(h) || ~isgraphics(h), return; end
if strcmp(role,'title')
    sizeValue = theme.plots.titleSize;
    weight = theme.plots.titleWeight;
    defaultInterpreter = theme.interpreters.title;
else
    sizeValue = theme.plots.labelSize;
    weight = theme.plots.labelWeight;
    defaultInterpreter = theme.interpreters.label;
end
local_set(h,'FontName',theme.plots.fontName);
local_set(h,'FontSize',round(sizeValue*scale));
local_set(h,'FontWeight',weight);
if isempty(interpreter)
    % An explicitly selected LaTeX interpreter is part of the semantic
    % content contract and must survive later figure-wide theming. Use the
    % role default only when the object does not already expose a setting.
    try
        interpreter = h.Interpreter;
    catch
        interpreter = defaultInterpreter;
    end
    if isempty(interpreter)
        interpreter = defaultInterpreter;
    end
end
local_set(h,'Interpreter',interpreter);
end

function local_set(object,name,value)
try
    if isprop(object,name), object.(name) = value; end
catch
end
end
