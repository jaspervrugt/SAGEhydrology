function sage_apply_ui_theme(container,theme)
%SAGE_APPLY_UI_THEME Apply canonical fonts to GUI controls and tables.
if nargin < 2 || isempty(theme), theme = sage_visual_theme(); end
if isempty(container) || ~isgraphics(container), return; end
objects = findall(container);
for k = 1:numel(objects)
    h = objects(k);
    className = class(h);
    if contains(className,'UITable')
        local_set(h,'FontName',theme.tables.fontName);
        local_set(h,'FontSize',theme.tables.fontSize);
    elseif contains(className,'UIControl') || contains(className,'ui.control')
        currentFont = '';
        try
            currentFont = char(string(h.FontName));
        catch
        end
        if strcmpi(currentFont,theme.fonts.family.monospace)
            continue
        end
        local_set(h,'FontName',theme.fonts.family.gui);
    end
end
end

function local_set(object,name,value)
try
    if isprop(object,name), object.(name) = value; end
catch
end
end
