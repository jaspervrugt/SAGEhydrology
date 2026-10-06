function sage_apply_colorbar_theme(object,theme,profile)
%SAGE_APPLY_COLORBAR_THEME Apply canonical colorbar typography.
if nargin < 2 || isempty(theme), theme = sage_visual_theme(); end
if nargin < 3 || isempty(profile), profile = 'gui'; end
name = char(lower(string(profile))); scale = 1;
if isfield(theme.profiles,name), scale = theme.profiles.(name).scale; end
objects = object(isgraphics(object));
for k = 1:numel(objects)
    h = objects(k);
    local_set(h,'FontName',theme.colorbar.fontName);
    local_set(h,'FontSize',round(theme.colorbar.fontSize*scale));
    local_set(h,'Color',theme.colorbar.color);
    local_set(h,'LineWidth',theme.colorbar.lineWidth);
    local_set(h,'TickLabelInterpreter',theme.interpreters.tick);
    if isprop(h,'Label')
        sage_apply_text_theme(h.Label,theme,'label',theme.interpreters.label);
    end
end
end

function local_set(object,name,value)
try
    if isprop(object,name), object.(name) = value; end
catch
end
end
