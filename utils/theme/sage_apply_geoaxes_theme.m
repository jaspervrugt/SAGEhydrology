function sage_apply_geoaxes_theme(ax,theme,profile)
%SAGE_APPLY_GEOAXES_THEME Apply canonical typography to geographic axes.
if nargin < 2 || isempty(theme), theme = sage_visual_theme(); end
if nargin < 3 || isempty(profile), profile = 'gui'; end
name = char(lower(string(profile))); scale = 1;
if isfield(theme.profiles,name), scale = theme.profiles.(name).scale; end
ax = ax(isgraphics(ax));
for k = 1:numel(ax)
    a = ax(k);
    local_set(a,'FontName',theme.plots.fontName);
    local_set(a,'FontSize',round(theme.plots.tickSize*scale));
    local_set(a,'FontWeight','normal');
    local_set(a,'GridColor',theme.plots.gridColor);
    local_set(a,'GridAlpha',theme.plots.gridAlpha);
    local_set(a,'LineWidth',theme.plots.axisLineWidth);
    if isprop(a,'Title')
        sage_apply_text_theme(a.Title,theme,'title',theme.interpreters.title);
    end
end
end

function local_set(object,name,value)
try
    if isprop(object,name), object.(name) = value; end
catch
end
end
