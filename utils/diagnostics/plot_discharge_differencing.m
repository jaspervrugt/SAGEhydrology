function handles = plot_discharge_differencing(result,axesHandle,response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT_DISCHARGE_DIFFERENCING Plot hourly discharge-error diagnostics.
%
% SYNOPSIS:
%   handles = plot_discharge_differencing(result)
%   handles = plot_discharge_differencing(result,axesHandle,response)
%   response = 'sigma' (default) or 'variance'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(axesHandle)
        figureHandle = figure('Color','w');
        axesHandle = axes(figureHandle);
    end
    if nargin < 3 ...
            || isempty(response)
        response = 'sigma';
    end
    response = lower(strtrim(char(string(response))));
    if ~isgraphics(axesHandle,'axes')
        error('plot_discharge_differencing:InvalidAxes', ...
            'axesHandle must be a valid axes object.');
    end

    cla(axesHandle);
    hold(axesHandle,'on');
    rawColor = [0.76 0.76 0.76];
    smoothColor = [0.05 0.35 0.75];
    fitColor = [0.05 0.05 0.05];
    paperColor = [0.15 0.15 0.15];

    switch response
        case 'sigma'
            yRaw = result.raw.sigma;
            ySmooth = result.smooth.sigma;
            yFit = result.fit_sigma.predicted;
            ylabelText = 'Estimated hourly error standard deviation';
            fit = result.fit_sigma;
            equation = sprintf( ...
                '\sigma = %.4g Q + %.4g', ...
                fit.slope,fit.intercept);
        case 'variance'
            yRaw = result.raw.sigma2;
            ySmooth = result.smooth.sigma2;
            yFit = result.fit_variance.predicted;
            ylabelText = 'Estimated hourly error variance';
            fit = result.fit_variance;
            equation = sprintf( ...
                '\sigma^2 = %.4g Q + %.4g', ...
                fit.slope,fit.intercept);
        otherwise
            error('plot_discharge_differencing:InvalidResponse', ...
                'response must be ''sigma'' or ''variance''.');
    end

    handles = struct();
    handles.raw = scatter(axesHandle, ...
        result.raw.discharge,yRaw,8,rawColor,'filled', ...
        'MarkerEdgeColor','none');
    handles.smooth = plot(axesHandle, ...
        result.smooth.discharge,ySmooth,'s', ...
        'Color',smoothColor, ...
        'MarkerFaceColor',smoothColor, ...
        'MarkerSize',3);
    handles.fit = plot(axesHandle, ...
        result.smooth.discharge,yFit,'-', ...
        'Color',fitColor,'LineWidth',2);
    handles.paper = gobjects(0);
    if strcmp(response,'variance')
        handles.paper = plot(axesHandle, ...
            result.smooth.discharge, ...
            result.fit_variance.paper_curve,'--', ...
            'Color',paperColor,'LineWidth',1.5);
    end

    xlabel(axesHandle,'Hourly discharge');
    ylabel(axesHandle,ylabelText);
    title(axesHandle,sprintf( ...
        'Order %d | %d valid pairs | window %d', ...
        result.options.order,result.n_pairs,result.smooth.window));
    grid(axesHandle,'on');
    box(axesHandle,'on');
    text(axesHandle,0.02,0.98,sprintf('%s\nR^2 = %.3f', ...
        equation,fit.r2), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'BackgroundColor','w', ...
        'Margin',4);

    if strcmp(response,'variance')
        legend(axesHandle, ...
            {'Raw estimate','Moving average', ...
            'Linear variance fit','Squared paper model'}, ...
            'Location','best');
    else
        legend(axesHandle, ...
            {'Raw estimate','Moving average','Linear fit'}, ...
            'Location','best');
    end
    hold(axesHandle,'off');

end
