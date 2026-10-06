function sage_apply_legend_theme(object,theme,profile)
%SAGE_APPLY_LEGEND_THEME Apply canonical legend typography and framing.
if nargin < 2 || isempty(theme), theme = sage_visual_theme(); end
if nargin < 3 || isempty(profile), profile = 'gui'; end
scale = 1;
name = char(lower(string(profile)));
if isfield(theme.profiles,name), scale = theme.profiles.(name).scale; end
objects = object(isgraphics(object));
for k = 1:numel(objects)
    h = objects(k);
    local_set(h,'FontName',theme.legend.fontName);
    local_set(h,'FontSize',round(theme.legend.fontSize*scale));
    local_set(h,'TextColor',theme.legend.textColor);
    local_set(h,'Box',theme.legend.box);
    % Preserve explicit LaTeX legends; only supply the canonical default
    % when no interpreter has been assigned.
    try
        currentInterpreter = h.Interpreter;
    catch
        currentInterpreter = '';
    end
    if isempty(currentInterpreter)
        local_set(h,'Interpreter',theme.interpreters.annotation);
    end
end
end

function local_set(object,name,value)
try
    if isprop(object,name), object.(name) = value; end
catch
end
end
