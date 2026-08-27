function plot_metric_ecdf_figure(figNo,part, ...
    mdl,dt_str,samp_word,scenarios,n_m,id_model, ...
    model_names,colors,model_name_fig, ...
    metric_tag,xlab_str)
%PLOT_METRIC_ECDF_FIGURE

figure(figNo); clf;
set(gcf,'Name', ...
    sprintf(['%s: Empirical Cumulative ' ...
    'Distribution Functions: %s'], ...
    model_name_fig,metric_tag), ...
    'NumberTitle','off', ...
    'color','w', ...
    'Units','inches', ...
    'Position',[0.5 0.5 24 7.5]);

ax = gobjects(4,1);
marginLeft = 0.055;
marginRight = 0.020;
panelBottom = 0.180;
panelHeight = 0.680;
panelWidth = panelHeight*(7.5/24);
gap = (1-marginLeft-marginRight-4*panelWidth)/3;
for ii = 1:4
    panelLeft = marginLeft + (ii-1)*(panelWidth+gap);
    ax(ii) = axes(gcf,'Units','normalized', ...
        'Position',[panelLeft panelBottom panelWidth panelHeight]);
    box(ax(ii),'on');
    hold(ax(ii),'on');
end

xL = -1;
xR = 1;
fill_alpha = 0.20;
line_width = 1.8;
med_lw = 1.1;
med_ms = 6;
% Keep the standalone helper synchronized with the ECDF typography used by
% plot_SAGE. These are the former S_FDC sizes, now shared by every metric.
fnt_axis = 18;
fnt_label = 18;
fnt_title = 16;
fnt_med = 16;

panelTitles = ...
    {'train basins | train period', ...
     'train basins | eval period', ...
     'eval basins | train period', ...
     'eval basins | eval period'};

for ii = 1:4
    set(ax(ii),'TickDir','out', ...
        'fontsize',fnt_axis, ...
        'linewidth',1);

    xlim(ax(ii),[xL xR]);
    ylim(ax(ii),[0 1]);
    axis(ax(ii),'square');

    xlabel(ax(ii),xlab_str, ...
        'interpreter','latex', ...
        'fontsize',fnt_label);

    ht = title(ax(ii),panelTitles{ii});
    set(ht,'interpreter','none', ...
        'fontsize',fnt_title);

    if ii == 1
        ylabel(ax(ii),'ECDF', ...
            'interpreter','latex', ...
            'fontsize',fnt_label);
    else
        set(ax(ii),'YTickLabel',[]);
        ylabel(ax(ii),'');
    end
end

for is = 1:4
    M = scenarios{is}.metric;

    if ~local_has_data(M)
        cla(ax(is));
        box(ax(is),'on');
        xlim(ax(is),[xL xR]);
        ylim(ax(is),[0 1]);
        axis(ax(is),'square');

        if is == 1
            ylabel(ax(is),'ECDF', ...
                'interpreter','latex', ...
                'fontsize',fnt_label);
        else
            set(ax(is),'YTickLabel',[]);
        end

        xlabel(ax(is),xlab_str, ...
            'interpreter','latex', ...
            'fontsize',fnt_label);
        ht = title(ax(is),panelTitles{is});
        set(ht,'interpreter','none', ...
            'fontsize',fnt_title);

        text(ax(is),0,0.5,'Not available', ...
            'horizontalalignment','center', ...
            'verticalalignment','middle', ...
            'interpreter','latex', ...
            'fontsize',20);
        drawnow;
        continue
    end

    for mdl_i = 1:n_m
        c = colors(id_model( ...
            min(mdl_i,numel(id_model))),:);
        im = id_model(min(mdl_i,numel(id_model)));
        nam = local_model_display_name(model_names{im});
        nFinite = nnz(isfinite(M(:,mdl_i)));
        lbl = sprintf('\\texttt{%s}\\; ($n = %d$)',nam,nFinite);

        plot_ecdf_panel(ax(is), ...
            M(:,mdl_i),c,lbl, ...
            part,metric_tag, ...
            scenarios{is}.tag, ...
            xL,xR,fill_alpha, ...
            line_width,med_lw, ...
            med_ms,fnt_med);
    end

    legend(ax(is),'show', ...
        'interpreter','latex', ...
        'fontsize',16, ...
        'location','northwest', ...
        'box','off');
    drawnow;
end

part_name = local_part_name(part,mdl);
part_name_tex = strrep(part_name,'_','\_');
if any(strcmpi(metric_tag,{'NSE','KGE','JKGE'}))
    metric_title = sprintf('$\\mathrm{%s}$',metric_tag);
else
    metric_title = sprintf('$%s$',metric_tag);
end

figureTitleSize = 21;
annotation(gcf,'textbox',[0.02 0.935 0.96 0.055], ...
    'String',sprintf(['$\\texttt{%s}$: ' ...
    '%s %s performance summary for basin/period scenarios'], ...
    part_name_tex,dt_str,metric_title), ...
    'Interpreter','latex', ...
    'FontSize',figureTitleSize, ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle');
end

% ================
% helper functions
% ================
function tf = local_has_data(A)
% LOCAL_HAS_DATA

tf = ~isempty(A);

if tf
    try
        tf = any(isfinite(A(:)));
    catch
        tf = true;
    end
end
end

function s = local_model_display_name(nameIn)
%LOCAL_MODEL_DISPLAY_NAME Return model name for display only

s = char(string(nameIn));

if strcmpi(strtrim(s),'xinanjiang')
    s = 'Xinanjiang';
end
end

function plot_ecdf_panel(axh,z,c, ...
    lbl,part,metric_tag, ...
    scenario_tag, ...
    xL,xR,fill_alpha,line_width, ...
    med_lw,med_ms,fnt_med)
%PLOT_ECDF_PANEL

z = z(:);
z = z(isfinite(z));

if isempty(z)
    return
end

[f,x] = sage_ecdf(z);
[xs,fs] = ecdf_to_stairs_fixed(x,f,xL,xR);

% -------------------------------------
% SITE: line only
% SAGE: line + fill + median annotation
% -------------------------------------
if strcmpi(part,'sage')
    [xp,yp] = stairs_fill_poly(xs,fs);
    p = patch(axh,xp,yp,c, ...
        'facealpha',fill_alpha, ...
        'edgealpha',0);
    set(p,'handlevisibility','off');
end

hLine = plot(axh,xs,fs, ...
    'color',c, ...
    'linewidth',line_width);
set(hLine,'displayname',lbl);

if strcmpi(part,'site')
    return
end

medX = median(z);
medF = ecdf_value_from_stairs(xs,fs,medX);

[xMark,yMark,x1m,x2m,xt,ha,mode] = ...
    median_marker_geom(medX,medF,xL,xR);

if ~mode.specialHalfAtLeft
    line(axh,[xMark xMark], ...
        [0 yMark], ...
        'color','k', ...
        'linewidth',med_lw, ...
        'handlevisibility','off');
end

line(axh,[x1m x2m], ...
    [yMark yMark], ...
    'color','k', ...
    'linewidth',med_lw, ...
    'handlevisibility','off');

line(axh,xMark,yMark, ...
    'marker','s', ...
    'markerfacecolor','k', ...
    'markeredgecolor','k', ...
    'linestyle','none', ...
    'markersize',med_ms, ...
    'handlevisibility','off');

txt = sprintf(['$\\widehat{T}_' ...
    '{\\mathrm{%s}_{\\rm %s}} = %.3f$'], ...
    metric_tag,scenario_tag,medX);

text(axh,xt,yMark,txt, ...
    'interpreter','latex', ...
    'fontweight','bold', ...
    'horizontalalignment',ha, ...
    'verticalalignment','middle', ...
    'color','k', ...
    'fontsize',fnt_med, ...
    'handlevisibility','off');
end

function name = local_part_name(part,mdl)
%LOCAL_PART_NAME

if strcmpi(part,'sage')
    if isempty(mdl.model)
        name = 'SAGE';
    else
        name = char(string( ...
            mdl.names(mdl.model(1))));
    end
else
    name = 'SITE';
end
end
