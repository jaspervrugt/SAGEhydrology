function sage_apply_text_theme(object,theme,role,interpreter)
%SAGE_APPLY_TEXT_THEME Apply semantic typography to text-like objects.

if nargin < 2 || isempty(theme), theme = sage_visual_theme(); end
if nargin < 3 || isempty(role), role = 'annotation'; end
if nargin < 4, interpreter = ''; end
objects = object(isgraphics(object));
for k = 1:numel(objects)
    h = objects(k);
    [family,sizeValue,weight,color,defaultInterpreter] = ...
        local_role(theme,lower(char(string(role))));
    local_set(h,'FontName',family);
    local_set(h,'FontSize',sizeValue);
    local_set(h,'FontWeight',weight);
    local_set(h,'Color',color);
    if isempty(interpreter), interpreterValue = defaultInterpreter;
    else, interpreterValue = interpreter;
    end
    local_set(h,'Interpreter',interpreterValue);
end
end

function [family,sizeValue,weight,color,interpreter] = local_role(t,role)
family = t.plots.fontName;
sizeValue = t.plots.annotationSize;
weight = 'normal';
color = t.plots.annotationColor;
interpreter = t.interpreters.annotation;
switch role
    case 'title'
        sizeValue=t.plots.titleSize; weight=t.plots.titleWeight;
        color=t.colors.neutral.dark; interpreter=t.interpreters.title;
    case 'label'
        sizeValue=t.plots.labelSize; weight=t.plots.labelWeight;
        color=t.plots.axisColor; interpreter=t.interpreters.label;
    case 'legend'
        sizeValue=t.plots.legendSize; color=t.legend.textColor;
    case {'placeholder','unavailable'}
        family=t.placeholders.fontName; sizeValue=t.placeholders.fontSize;
        weight=t.placeholders.fontWeight; color=t.placeholders.color;
    case {'gui','control'}
        family=t.fonts.family.gui; sizeValue=t.fonts.gui.gui;
        color=t.colors.neutral.dark; interpreter=t.interpreters.ui;
    case {'table','tableheader'}
        family=t.tables.fontName; sizeValue=t.tables.fontSize;
        color=t.colors.neutral.dark; interpreter=t.interpreters.ui;
        if strcmp(role,'tableheader'), weight=t.tables.headerWeight; end
    case {'log','code','monospace'}
        family=t.fonts.family.monospace; sizeValue=t.fonts.gui.notes;
        color=t.colors.neutral.dark; interpreter='none';
    case 'equation'
        interpreter=t.interpreters.equation;
end
end

function local_set(object,name,value)
try
    if isprop(object,name), object.(name) = value; end
catch
end
end
