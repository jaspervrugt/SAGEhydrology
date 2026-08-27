function outFile = save_to_pptx(C,outFile)
%SAVE_TO_PPTX Export SAGE plot figures to PPTX, with PDF fallback.
%
% SYNOPSIS:
%   save_to_pptx(C)
%   save_to_pptx(C,outFile)
%   outFile = save_to_pptx(...)
%
% DESCRIPTION:
%   Exports visible MATLAB plot figures to a PowerPoint file when
%   Microsoft PowerPoint is available on Windows. If PowerPoint is not
%   available, or if the requested output extension is .pdf, the function
%   writes a multi-page PDF instead. Figures are first rasterized to PNG
%   files, then inserted into PPTX/PDF. This is usually more stable in
%   deployed MATLAB applications than exporting live graphics objects
%   directly to PowerPoint or PDF.
%
% INPUT:
%   C        SAGE configuration structure (optional)
%   outFile  output filename, .pptx or .pdf (optional)
%
% OUTPUT:
%   outFile  path to the file that was actually written
%
% NOTES:
%   - Exports visible MATLAB plot figures only
%   - Skips the main SAGE GUI/uifigure window
%   - PPTX writing through ActiveX requires Microsoft PowerPoint on Windows
%   - If PowerPoint is unavailable, a PDF file is written instead
%   - PDF files are assembled from PNG images for robustness
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 1
        C = [];
    end
    
    if nargin < 2 ...
            || isempty(outFile)
        outFile = local_default_pptx_name(C);
    end
    
    outFile = char(string(outFile));
    [outDir,outBase,outExt] = fileparts(outFile);
    
    if isempty(outDir)
        outDir = pwd;
    end
    if isempty(outExt)
        outExt = '.pptx';
    end
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    
    figs = local_get_export_figures();
    
    if isempty(figs)
        error(['No exportable ' ...
            'plot figures found.']);
    end
    
    tmpDir = tempname;
    mkdir(tmpDir);
    
    try
        imgFiles = local_export_figures_to_png( ...
            figs,tmpDir);
    
        if isempty(imgFiles)
            error(['No figures could ' ...
                'be exported to PNG.']);
        end
    
        titlePng = '';
        schematicFile = local_find_sage_schematic(C,true);
        schematicPng = local_find_sage_schematic(C,false);
        modelDiagramFile = local_find_model_diagram(C,true);
        modelDiagramPng = local_find_model_diagram(C,false);
    
        if ~isempty(C)
            titlePng = local_make_title_png( ...
                C,tmpDir);
        end

        previewPng = local_make_preview_png(C,tmpDir);
        mapsPng = local_make_region_zones_png(C,tmpDir);
        networkPng = local_make_network_png(C,tmpDir);
        reportFiles = {networkPng; previewPng; mapsPng};
        reportFiles = reportFiles(~cellfun(@isempty,reportFiles));
        imgFiles = [reportFiles; imgFiles(:)];
    
        wantPdf = strcmpi(outExt,'.pdf');
    
        if ~wantPdf && ispc
            pptFile = fullfile(outDir, ...
                [outBase '.pptx']);
            try
                local_write_pptx_from_png( ...
                    C,pptFile,schematicFile, ...
                    titlePng,modelDiagramFile, ...
                    modelDiagramPng,imgFiles);
                outFile = pptFile;
            catch MEppt
                fprintf(['PPTX export failed: ' ...
                    '%s\n'],MEppt.message);
                fprintf(['Writing PDF fallback ' ...
                    'instead.\n']);
        
                outFile = fullfile(outDir, ...
                    [outBase '.pdf']);
                local_write_pdf_from_png( ...
                    outFile,schematicPng, ...
                    titlePng,modelDiagramPng,imgFiles);
            end
        else
            outFile = fullfile(outDir, ...
                [outBase '.pdf']);
            local_write_pdf_from_png( ...
                outFile,schematicPng, ...
                titlePng,modelDiagramPng,imgFiles);
        end
    
        fprintf(['Saved figure ' ...
            'report to %s\n'],outFile);
        local_close_export_figures();
        
        if ~isdeployed
            local_open_file_safely(outFile);
        end
    
    catch ME
        try
            rmdir(tmpDir,'s');
        catch
        end
        rethrow(ME);
    end
    
    try
        rmdir(tmpDir,'s');
    catch
    end
end

function networkPng = local_make_network_png(C,tmpDir)
%LOCAL_MAKE_NETWORK_PNG Export the configured FFN schematic.

% Reconstruct the schematic in an independent figure using the same
% drawing function as the GUI. This keeps the slide independent of GUI
% tab visibility and gives GUI, script, and deployed runs identical output.

    networkPng = fullfile(tmpDir,'ffn_network.png');

    fig = [];
    try
        fig = figure('Visible','off', ...
            'Color','white', ...
            'Units','pixels', ...
            'Position',[100 100 1400 760]);
        axExport = axes(fig, ...
            'Units','normalized', ...
            'Position',[0.02 0.02 0.96 0.96], ...
            'Color','white');

        nIn = local_input_size(C);
        nHidden = local_hidden_vector(C.net);
        nOut = local_output_size(C);
        tf = local_ann_tf(C.net);
        draw_FFN_schematic( ...
            axExport,nIn,nHidden,nOut,tf);

        exportgraphics(fig,networkPng, ...
            'ContentType','image', ...
            'Resolution',220, ...
            'BackgroundColor','white');
        close(fig);
    catch ME
        if ~isempty(fig) && isgraphics(fig)
            close(fig);
        end
        networkPng = '';
        fprintf(['Skipping FFN network slide: ' ...
            '%s\n'],ME.message);
    end

end

function previewPng = local_make_preview_png(C,tmpDir)
%LOCAL_MAKE_PREVIEW_PNG Render the GUI Preview text as a report slide.

    previewPng = '';
    try
        if isempty(C) || ~isfield(C,'runtime') ...
                || ~isfield(C.runtime,'previewText') ...
                || isempty(C.runtime.previewText)
            return
        end

        value = string(C.runtime.previewText(:));
        value = value(strlength(value) > 0);
        if isempty(value)
            return
        end

        previewPng = fullfile(tmpDir,'preview_text.png');
        fig = figure('Visible','off','Color','w','Units','pixels', ...
            'Position',[100 100 1600 900],'PaperPositionMode','auto');
        ax = axes(fig,'Position',[0 0 1 1]);
        axis(ax,'off');
        text(ax,0.04,0.95,'SAGE configuration preview', ...
            'Units','normalized','FontSize',24,'FontWeight','bold', ...
            'Interpreter','none','VerticalAlignment','top');
        text(ax,0.04,0.89,char(join(value,newline)), ...
            'Units','normalized','FontName','Consolas','FontSize',10, ...
            'Interpreter','none','VerticalAlignment','top');
        exportgraphics(fig,previewPng,'Resolution',150, ...
            'BackgroundColor','white','ContentType','image');
        close(fig);
    catch
        previewPng = '';
    end
end

function mapsPng = local_make_region_zones_png(C,tmpDir)
%LOCAL_MAKE_REGION_ZONES_PNG Combine the GUI Region and Zones maps.

    mapsPng = '';
    try
        if isempty(C) || ~isfield(C,'runtime') ...
                || ~isfield(C.runtime,'regionAxes') ...
                || ~isfield(C.runtime,'zonesAxes')
            return
        end
        regionAxes = C.runtime.regionAxes;
        zonesAxes = C.runtime.zonesAxes;
        if ~isgraphics(regionAxes) || ~isgraphics(zonesAxes)
            return
        end

        regionFile = fullfile(tmpDir,'region_map.png');
        zonesFile = fullfile(tmpDir,'zones_map.png');
        exportgraphics(regionAxes,regionFile,'Resolution',150);
        exportgraphics(zonesAxes,zonesFile,'Resolution',150);
        regionImage = imread(regionFile);
        zonesImage = imread(zonesFile);

        mapsPng = fullfile(tmpDir,'region_zones_maps.png');
        fig = figure('Visible','off','Color','w','Units','pixels', ...
            'Position',[100 100 1600 900],'PaperPositionMode','auto');
        layout = tiledlayout(fig,1,2,'Padding','compact', ...
            'TileSpacing','compact');
        axRegion = nexttile(layout,1);
        image(axRegion,regionImage);
        axis(axRegion,'image');
        axis(axRegion,'off');
        title(axRegion,'Region','FontSize',18,'FontWeight','normal');
        axZones = nexttile(layout,2);
        image(axZones,zonesImage);
        axis(axZones,'image');
        axis(axZones,'off');
        title(axZones,'Zones','FontSize',18,'FontWeight','normal');
        exportgraphics(fig,mapsPng,'Resolution',150, ...
            'BackgroundColor','white','ContentType','image');
        close(fig);
    catch
        mapsPng = '';
    end
end

function figs = local_get_export_figures()
%LOCAL_GET_EXPORT_FIGURES Return visible plot figures safe for export.

    figs0 = findall(groot,'Type','figure');
    figs = matlab.ui.Figure.empty;
    
    for i = 1:numel(figs0)
    
        f = figs0(i);
    
        try
            if ~isvalid(f)
                continue
            end
    
            if isprop(f,'Visible') ...
                    && strcmpi(f.Visible,'off')
                continue
            end
    
            % Skip App Designer/uifigure GUI windows.
            try
                kids = findall(f);
                cls = arrayfun(@(h) string(class(h)),kids, ...
                    'UniformOutput',false);
                cls = string(cls);
    
                if any(contains(cls, ...
                        ["matlab.ui.control.HTML", ...
                         "matlab.ui.container.GridLayout", ...
                         "matlab.ui.container.TabGroup", ...
                         "matlab.ui.container.Tab", ...
                         "matlab.ui.control.Table", ...
                         "matlab.ui.control.Button", ...
                         "matlab.ui.control.DropDown", ...
                         "matlab.ui.control.EditField"]))
                    continue
                end
            catch
            end
    
            nm = "";
            if isprop(f,'Name')
                nm = string(f.Name);
            end

            nml = lower(nm);
            % Skip main SAGE GUI/app window, but keep SAGE result figures.
            if contains(nml,"sage") ...
                    && ~contains(nml, ...
                    ["flow","duration","ecdf", ...
                    "trace","performance","time","series","fdc", ...
                    "conus","zones","variogram","attribution", ...
                    "summary","preview","diagnostic","diagnostics", ...
                    "history","parameter"])
                continue
            end
    
            ax = findall(f,'Type','axes');
            ax = ax(arrayfun(@(a) isvalid(a),ax));
    
            if isempty(ax)
                continue
            end
    
            figs(end+1,1) = f; %#ok<AGROW>
    
        catch
        end
    end
    
    try
        [~,idx] = sort(double(figs));
        figs = figs(idx);
    catch
    end
end

function imgFiles = local_export_figures_to_png(figs,tmpDir)
%LOCAL_EXPORT_FIGURES_TO_PNG Rasterize MATLAB figures to PNG files.

    imgFiles = cell(numel(figs),1);
    
    for k = 1:numel(figs)
    
        f = figs(k);
        figureName = lower(local_fig_name(f));
        isFDC = contains(figureName, ...
            ["fdc","flow duration"]);
        isMetricECDF = contains(figureName, ...
            "empirical cumulative distribution functions");
        isTimeSeries = contains(figureName,"time") ...
            || contains(figureName,"series") ...
            || contains(figureName,"training basin") ...
            || contains(figureName,"evaluation basin");
        imgFiles{k} = fullfile(tmpDir, ...
            sprintf('fig_%03d.png',k));
    
        try
            fprintf('Exporting figure %d/%d: %s\n', ...
                k,numel(figs), ...
                local_fig_name(f));
    
            drawnow limitrate nocallbacks
            pause(0.2);
            % Avoid enormous bitmap exports in deployed mode.
            try
                oldUnits = f.Units;
                cleanupUnits = onCleanup(@() ...
                    set(f,'Units',oldUnits));
                f.Units = 'pixels';
                pos = f.Position;
                % All metric ECDF figures must be rasterized at the same
                % dimensions. Previously S_FDC matched the generic FDC
                % branch (1600 px), while NSE/KGE used 1400 px. PowerPoint
                % then fitted both images to the same slide area, making
                % the NSE/KGE typography appear about 14% larger.
                if isMetricECDF || isFDC
                    maxWidth = 1600;
                    maxHeight = 1100;
                elseif isTimeSeries
                    maxWidth = 1800;
                    maxHeight = 1150;
                else
                    maxWidth = 1400;
                    maxHeight = 900;
                end
                scale = min([1, ...
                    maxWidth/max(pos(3),1), ...
                    maxHeight/max(pos(4),1)]);
                pos(3:4) = max(1,round(pos(3:4)*scale));
                f.Position = pos;
                drawnow;
                clear cleanupUnits
            catch
            end
    
        % Robust export path. FDC/time-series figures can contain log axes,
        % many markers, legends, and annotations. For those figures, MATLAB's
        % painters renderer/exportgraphics path can hang or crash in deployed GUI
        % workflows. Prefer an painters screen-capture raster export.
        local_export_one_figure_png(f,imgFiles{k},isFDC);
    
        catch ME1
	        fprintf(['Primary export failed ' ...
                'for figure %d: %s\n'],k,ME1.message);
            try
                exportgraphics(f, ...
                    imgFiles{k},'Resolution',120);
            catch ME2
                fprintf(['Skipping figure ' ...
                    '%d: %s\n'],k,ME2.message);
                imgFiles{k} = '';
            end
        end
        drawnow;
        pause(0.02);
        close(f)
    end
    
    imgFiles = imgFiles(~cellfun(@isempty,imgFiles));
end

function local_export_one_figure_png(f,pngFile,isFDC)
%LOCAL_EXPORT_ONE_FIGURE_PNG Robust figure-to-PNG export for PPT.
%
% For FDC and time-series figures, avoid the painters renderer. These
% figures often contain log axes, many line objects, legend-only handles,
% and annotations; in compiled/AppDesigner workflows those are much more
% stable as raster screen captures.

    if nargin < 3
        isFDC = false;
    end
    
    try
        set(f,'Color','w');
    catch
    end
    try
        set(f,'InvertHardcopy','off');
    catch
    end
    try
        set(f,'Renderer','painters');
    catch
    end
    
    drawnow;
    pause(0.05);
    
    % Safest path for FDC/time-series figures: capture the rendered figure.
    nm = lower(local_fig_name(f));
    isTimeSeries = contains(nm,"time") ...
        || contains(nm,"series") ...
        || contains(nm,"training basin") ...
        || contains(nm,"evaluation basin");
    
    nm = lower(local_fig_name(f));
    if contains(nm,"parameter maps") ...
            || contains(nm,"parameter")
        drawnow;
        pause(2.0);
        drawnow;
    end

    if isFDC || isTimeSeries
        figure(f);
        drawnow;
        pause(0.1);
        fr = getframe(f);
        if isempty(fr.cdata)
            error('save_to_pptx:EmptyFigureCapture', ...
                'The rendered figure capture is empty.');
        end
        imwrite(fr.cdata,pngFile);
        return
    end

    fr = getframe(f);
    imwrite(fr.cdata,pngFile);

end

function local_write_pptx_from_png(C, ...
    outFile,schematicPng,titlePng, ...
    modelDiagramFile,modelDiagramPng,imgFiles)
%LOCAL_WRITE_PPTX_FROM_PNG Create PowerPoint file from image files.
%
% The SAGE schematic and MATLAB figures are normally inserted as raster
% images. The selected hydrologic-model diagram is inserted as SVG when
% available so that it remains a vector graphic in PowerPoint. If the
% installed PowerPoint version cannot insert SVG, the matching PNG file is
% used automatically.
%
% This function is safe in deployed mode provided that the docs directory
% is included in compiler.build AdditionalFiles and the file lookup helper
% searches ctfroot.

    ppt = [];
    pres = [];

    try
        % -----------------------------------
        % Start Microsoft PowerPoint via COM
        % -----------------------------------
        ppt = actxserver('PowerPoint.Application');

        try
            set(ppt,'DisplayAlerts',0);
        catch
        end

        try
            set(ppt,'Visible',0);
        catch
        end

        pres = invoke(ppt.Presentations,'Add');
        ppLayoutBlank = 12;

        % ------------------------
        % SAGE schematic slide
        % ------------------------
        if ~isempty(schematicPng) ...
                && isfile(schematicPng)
            try
                local_add_image_slide( ...
                    pres,ppLayoutBlank, ...
                    schematicPng,5);
                pause(0.2);
            catch MEschematic
                fprintf(['Skipping SAGE schematic ' ...
                    'slide: %s\n'], ...
                    MEschematic.message);
            end
        end

        % ------------------------
        % Run-information slide
        % ------------------------
        if ~isempty(C)
            local_add_title_slide( ...
                pres,ppLayoutBlank,C);
            pause(0.2);

            % Place the screening audit immediately after the run settings
            % and before model diagrams or performance results.
            local_add_data_quality_slide( ...
                pres,ppLayoutBlank,C);
            pause(0.2);

        elseif ~isempty(titlePng) ...
                && isfile(titlePng)
            local_add_image_slide( ...
                pres,ppLayoutBlank,titlePng);
            pause(0.2);
        end

        % ---------------------------------
        % Hydrologic-model diagram slide
        % ---------------------------------
        % Prefer SVG because it remains a vector graphic in PowerPoint.
        % Some older PowerPoint versions cannot insert SVG through ActiveX;
        % in that case retry automatically using the matching PNG.
        modelSlideAdded = false;

        if ~isempty(modelDiagramFile) ...
                && isfile(modelDiagramFile)

            [~,~,modelExt] = fileparts( ...
                modelDiagramFile);

            try
                fprintf(['Adding model diagram: ' ...
                    '%s\n'],modelDiagramFile);

                local_add_image_slide( ...
                    pres,ppLayoutBlank, ...
                    modelDiagramFile,20);

                pause(0.2);
                modelSlideAdded = true;

            catch MEmodel
                fprintf(['Could not insert model diagram ' ...
                    '%s: %s\n'], ...
                    upper(modelExt), ...
                    MEmodel.message);
            end
        end

        % Use PNG if:
        %   1. no SVG/model file was found, or
        %   2. PowerPoint rejected the SVG.
        if ~modelSlideAdded ...
                && ~isempty(modelDiagramPng) ...
                && isfile(modelDiagramPng)

            try
                fprintf(['Adding PNG fallback for ' ...
                    'model diagram: %s\n'], ...
                    modelDiagramPng);

                local_add_image_slide( ...
                    pres,ppLayoutBlank, ...
                    modelDiagramPng,20);

                pause(0.2);
                modelSlideAdded = true;

            catch MEmodelPng
                fprintf(['Skipping model diagram ' ...
                    'slide: %s\n'], ...
                    MEmodelPng.message);
            end
        end

        if ~modelSlideAdded
            fprintf(['No model diagram slide was ' ...
                'added.\n']);
        end

        % ------------------------
        % MATLAB results figures
        % ------------------------
        for k = 1:numel(imgFiles)
            try
                fprintf(['Adding PPT slide ' ...
                    '%d/%d: %s\n'], ...
                    k,numel(imgFiles), ...
                    imgFiles{k});

                if ~isfile(imgFiles{k})
                    fprintf(['Skipping missing ' ...
                        'PNG: %s\n'], ...
                        imgFiles{k});
                    continue
                end

                local_add_image_slide( ...
                    pres,ppLayoutBlank, ...
                    imgFiles{k});

                pause(0.05);

                if mod(k,10) == 0
                    pause(0.5);
                end

            catch MEslide
                fprintf(['Skipping PPT slide ' ...
                    '%d: %s\n'], ...
                    k,MEslide.message);
            end
        end

        % ----------------
        % Save presentation
        % ----------------
        if exist(outFile,'file')
            delete(outFile);
        end

        fprintf('Saving PPTX: %s\n',outFile);
        pause(0.5);

        invoke(pres,'SaveAs',outFile);

        fprintf('Closing PPTX\n');
        invoke(pres,'Close');
        pres = [];

        fprintf('Quitting PowerPoint\n');
        invoke(ppt,'Quit');
        %ppt = [];

    catch ME
        % Ensure that PowerPoint is not left running following an error.
        try
            if ~isempty(pres)
                invoke(pres,'Close');
            end
        catch
        end

        try
            if ~isempty(ppt)
                invoke(ppt,'Quit');
            end
        catch
        end

        rethrow(ME);
    end
end

function local_write_pdf_from_png(pdfFile, ...
    schematicPng,titlePng,modelDiagramPng,imgFiles)
%LOCAL_WRITE_PDF_FROM_PNG Create a multi-page PDF from PNG images.

    if exist(pdfFile,'file')
        delete(pdfFile);
    end
    
    allFiles = {};

    if ~isempty(schematicPng)
        allFiles = [allFiles; {schematicPng}];
    end

    if ~isempty(titlePng)
        allFiles = [allFiles; {titlePng}];
    end

    % PDF fallback cannot reliably rasterize SVG files here. Use the
    % matching PNG diagram when it is available.
    if ~isempty(modelDiagramPng)
        allFiles = [allFiles; {modelDiagramPng}];
    end

    allFiles = [allFiles; imgFiles(:)];

    page = 0;

    for k = 1:numel(allFiles)

        imageFile = char(string(allFiles{k}));
        [~,~,imageExt] = fileparts(imageFile);

        if ~isfile(imageFile) ...
                || strcmpi(imageExt,'.svg')
            fprintf(['Skipping non-raster PDF ' ...
                'asset: %s\n'],imageFile);
            continue
        end

        try
            img = imread(imageFile);
        catch MEread
            fprintf(['Skipping unreadable PDF ' ...
                'asset %s: %s\n'], ...
                imageFile,MEread.message);
            continue
        end
        [h,w,~] = size(img);
    
        figW = min(w,2200);
        figH = min(h,1400);
    
        fig = figure('Visible','off', ...
            'Color','w', ...
            'Units','pixels', ...
            'Position',[100 100 figW figH], ...
            'PaperPositionMode','auto');
    
        ax = axes(fig, ...
            'Units','normalized', ...
            'Position',[0 0 1 1]);
    
        image(ax,img);
        axis(ax,'image');
        axis(ax,'off');
    
        drawnow limitrate nocallbacks
    
        page = page + 1;

        if page == 1
            exportgraphics(fig,pdfFile, ...
                'ContentType','image', ...
                'Resolution',150);
        else
            exportgraphics(fig,pdfFile, ...
                'ContentType','image', ...
                'Resolution',150, ...
                'Append',true);
        end
    
        close(fig);
    end

    if page == 0
        error('save_to_pptx:NoPdfPages', ...
            'No raster images were available for PDF export.');
    end
end

function titlePng = local_make_title_png(C,tmpDir)
%LOCAL_MAKE_TITLE_PNG Create a raster title page for PDF fallback.

    titlePng = fullfile(tmpDir,'title_page.png');
    
    fig = figure('Visible','off', ...
        'Color','w', ...
        'Units','pixels', ...
        'Position',[100 100 1200 675], ...
        'PaperPositionMode','auto');
    
    ax = axes(fig,'Position',[0 0 1 1]);
    axis(ax,'off');
    
    modelName = local_model_name(C);
    dtName = local_dt_name(C);
    lossStr = local_loss_string(C);
    regionName = local_region_string(C);
    
    nTrain = local_getfield_or_nan(C.bas,'K_t');
    nEval = local_getfield_or_nan(C.bas,'K_e');
    
    createdStr = char(string(datetime('now'), ...
        'dd-MMM-uuuu HH:mm'));
    
    trainPeriodStr = local_training_period_string(C);
    evalPeriodStr = local_evaluation_period_string(C);
    spinupStr = local_spinup_string(C);

    splitStr = local_split_method_string(C);
    sampleStr = local_basin_sample_string(C);

    txt = local_title_text(C,modelName,dtName,lossStr, ...
        regionName,nTrain,nEval,createdStr, ...
        trainPeriodStr,evalPeriodStr,spinupStr, ...
        splitStr,sampleStr);

    text(0.05,0.92,'SAGE results', ...
        'FontSize',28, ...
        'FontWeight','bold', ...
        'Interpreter','none');
    
    text(0.05,0.82,txt, ...
        'FontSize',14, ...
        'Interpreter','none', ...
        'VerticalAlignment','top');
    
    note = local_run_note_string(C);
    if strlength(note) > 0
        if strlength(note) > 1000
            note = extractBefore(note,1000) + " ...";
        end
    
        text(0.05,0.08,['Run note:' newline char(note)], ...
            'FontSize',11, ...
            'FontAngle','italic', ...
            'Interpreter','none', ...
            'VerticalAlignment','top');
    end
    
    % drawnow;
    % print(fig,titlePng,'-dpng','-r150');
    drawnow limitrate nocallbacks
    exportgraphics(fig,titlePng, ...
        'Resolution',150, ...
        'BackgroundColor','white', ...
        'ContentType','image');
    
    close(fig);
end

function local_open_file_safely(outFile)
%LOCAL_OPEN_FILE_SAFELY Try to open output file after writing.

    try
        if ispc
            winopen(outFile);
        elseif ismac
            system(sprintf('open "%s"',outFile));
        else
            system(sprintf('xdg-open "%s" &',outFile));
        end
    catch
        fprintf(['Could not automatically ' ...
            'open file: %s\n'],outFile);
    end
end

function pptFile = local_default_pptx_name(C)
%LOCAL_DEFAULT_PPTX_NAME Build default PPTX filename from SAGE settings.

    if isempty(C)
        pptFile = fullfile(pwd,'all_figures.pptx');
        return
    end
    
    modelName = local_model_name(C);
    modelName = regexprep(char( ...
        string(modelName)),'[^A-Za-z0-9_]+','');
    
    dtName = local_dt_name(C);
    stamp = char(string(datetime('now'), ...
        'yyyyMMdd_HHmm'));
    
    if isfield(C,'dirres') ...
            && ~isempty(C.dirres)
        outDir = C.dirres;
    else
        outDir = pwd;
    end
    
    regionCode = local_region_code(C);

    pptFile = fullfile(outDir,sprintf( ...
        'SAGE_results_%s_%s_%s_%s_%s.pptx', ...
        regionCode,modelName,dtName, ...
        local_split_method(C),stamp));
end

function local_add_title_slide(pres,ppLayoutBlank,C)
%LOCAL_ADD_TITLE_SLIDE Add first slide with SAGE run settings.

    slides = get(pres,'Slides');
    slideCount = get(slides,'Count');
    slide = invoke(slides, ...
        'Add',slideCount+1,ppLayoutBlank);
    
    shapes = get(slide,'Shapes');
    
    slideW = get(pres.PageSetup, ...
        'SlideWidth');
    
    titleBox = invoke(shapes, ...
        'AddTextbox',1,40,35,slideW-80,45);
    titleBox.TextFrame.TextRange.Text = 'SAGE results';
    titleBox.TextFrame.TextRange.Font.Size = 28;
    titleBox.TextFrame.TextRange.Font.Bold = 1;
    
    modelName = local_model_name(C);
    dtName = local_dt_name(C);
    lossStr = local_loss_string(C);
    regionName = local_region_string(C);
    
    nTrain = local_getfield_or_nan( ...
        C.bas,'K_t');
    nEval = local_getfield_or_nan( ...
        C.bas,'K_e');
    
    createdStr = char(string( ...
        datetime('now'),'dd-MMM-uuuu HH:mm'));
    trainPeriodStr = ...
        local_training_period_string(C);
    evalPeriodStr = ...
        local_evaluation_period_string(C);
    spinupStr = ...
        local_spinup_string(C);
    
    splitStr = local_split_method_string(C);
    sampleStr = local_basin_sample_string(C);

    txt = local_title_text(C,modelName,dtName,lossStr, ...
        regionName,nTrain,nEval,createdStr, ...
        trainPeriodStr,evalPeriodStr,spinupStr, ...
        splitStr,sampleStr);
    
    textBox = invoke(shapes,'AddTextbox', ...
        1,60,90,slideW-120,335);
    textBox.TextFrame.TextRange.Text = txt;
    textBox.TextFrame.TextRange.Font.Size = 11;
    
    note = local_run_note_string(C);
    if strlength(note) > 0
        if strlength(note) > 1000
            note = extractBefore(note,1000) + " ...";
        end
    
        noteBox = invoke(shapes, ...
            'AddTextbox',1,60,430,slideW-120,90);
        noteBox.TextFrame.TextRange.Text = ...
            ['Run note:' newline char(note)];
        noteBox.TextFrame.TextRange.Font.Size = 11;
        noteBox.TextFrame.TextRange.Font.Italic = 1;
    end
end

function local_add_data_quality_slide(pres,ppLayoutBlank,C)
%LOCAL_ADD_DATA_QUALITY_SLIDE Summarize requested and active populations.

    if ~isfield(C,'bas') || ~isstruct(C.bas) ...
            || ~isfield(C.bas,'data_quality') ...
            || isempty(C.bas.data_quality)
        return
    end
    dq = C.bas.data_quality;
    required = {'requested_K','requested_K_t','requested_K_e', ...
        'forcing_complete','discharge_complete'};
    if ~all(isfield(dq,required))
        return
    end

    forcingComplete = logical(dq.forcing_complete(:));
    dischargeComplete = logical(dq.discharge_complete(:));
    if numel(forcingComplete) ~= numel(dischargeComplete)
        return
    end
    badForcing = ~forcingComplete;
    badDischarge = ~dischargeComplete;
    requestedK = double(dq.requested_K);
    requestedKt = double(dq.requested_K_t);
    requestedKe = double(dq.requested_K_e);
    isTrain = false(requestedK,1);
    isTrain(1:min(requestedKt,requestedK)) = true;
    isEval = ~isTrain;
    forcingMask = badForcing & ~badDischarge;
    dischargeMask = ~badForcing & badDischarge;
    bothMask = badForcing & badDischarge;
    forcingCount = [nnz(forcingMask),nnz(forcingMask & isTrain), ...
        nnz(forcingMask & isEval)];
    dischargeCount = [nnz(dischargeMask),nnz(dischargeMask & isTrain), ...
        nnz(dischargeMask & isEval)];
    bothCount = [nnz(bothMask),nnz(bothMask & isTrain), ...
        nnz(bothMask & isEval)];
    activeKt = double(C.bas.K_t);
    activeKe = double(C.bas.K_e);
    activeK = activeKt + activeKe;
    excludedK = requestedK - activeK;

    slides = get(pres,'Slides');
    slideCount = get(slides,'Count');
    slide = invoke(slides,'Add',slideCount+1,ppLayoutBlank);
    shapes = get(slide,'Shapes');
    slideW = get(pres.PageSetup,'SlideWidth');

    titleBox = invoke(shapes,'AddTextbox',1,40,30,slideW-80,45);
    titleBox.TextFrame.TextRange.Text = 'Data-quality screening';
    titleBox.TextFrame.TextRange.Font.Size = 28;
    titleBox.TextFrame.TextRange.Font.Bold = 1;

    introBox = invoke(shapes,'AddTextbox',1,55,85,slideW-110,42);
    introBox.TextFrame.TextRange.Text = [ ...
        'The universal basin inventory is screened before SAGE training.'];
    introBox.TextFrame.TextRange.Font.Size = 16;

    requestedBox = invoke(shapes,'AddTextbox', ...
        1,70,135,slideW-140,42);
    requestedBox.TextFrame.TextRange.Text = sprintf( ...
        'Requested:   K = %d      K_t = %d      K_e = %d', ...
        requestedK,requestedKt,requestedKe);
    requestedBox.TextFrame.TextRange.Font.Size = 20;
    requestedBox.TextFrame.TextRange.Font.Bold = 1;

    activeBox = invoke(shapes,'AddTextbox', ...
        1,70,185,slideW-140,42);
    activeBox.TextFrame.TextRange.Text = sprintf( ...
        'Active:          K = %d      K_t = %d      K_e = %d', ...
        activeK,activeKt,activeKe);
    activeBox.TextFrame.TextRange.Font.Size = 20;
    activeBox.TextFrame.TextRange.Font.Bold = 1;
    activeBox.TextFrame.TextRange.Font.Color.RGB = 32768;

    excludedBox = invoke(shapes,'AddTextbox', ...
        1,70,235,slideW-140,42);
    excludedBox.TextFrame.TextRange.Text = sprintf( ...
        'Excluded:      K = %d      training = %d      evaluation = %d', ...
        excludedK,requestedKt-activeKt,requestedKe-activeKe);
    excludedBox.TextFrame.TextRange.Font.Size = 20;
    excludedBox.TextFrame.TextRange.Font.Bold = 1;
    excludedBox.TextFrame.TextRange.Font.Color.RGB = 192;

    categoryTitle = invoke(shapes,'AddTextbox', ...
        1,70,300,slideW-140,32);
    categoryTitle.TextFrame.TextRange.Text = ...
        'Mutually exclusive exclusion categories';
    categoryTitle.TextFrame.TextRange.Font.Size = 18;
    categoryTitle.TextFrame.TextRange.Font.Bold = 1;

    % Use a native PowerPoint table rather than space-padded text. Font
    % metrics vary between PowerPoint installations and otherwise cause the
    % Total, Train, and Eval headings to drift away from their values.
    tableShape = invoke(shapes,'AddTable',4,4,70,335,slideW-140,112);
    pptTable = get(tableShape,'Table');
    labels = { ...
        'Criterion','Total','Train','Eval'; ...
        'Incomplete meteorological forcing only', ...
            forcingCount(1),forcingCount(2),forcingCount(3); ...
        'Insufficient or constant discharge only', ...
            dischargeCount(1),dischargeCount(2),dischargeCount(3); ...
        'Both forcing and discharge criteria', ...
            bothCount(1),bothCount(2),bothCount(3)};
    try
        pptColumns = get(pptTable,'Columns');
        firstColumn = invoke(pptColumns,'Item',1);
        set(firstColumn,'Width',slideW-380);
        for column = 2:4
            numberColumn = invoke(pptColumns,'Item',column);
            set(numberColumn,'Width',80);
        end
    catch
    end
    for row = 1:4
        for column = 1:4
            tableCell = invoke(pptTable,'Cell',row,column);
            cellShape = get(tableCell,'Shape');
            if isnumeric(labels{row,column})
                cellText = sprintf('%d',labels{row,column});
            else
                cellText = labels{row,column};
            end
            cellShape.TextFrame.TextRange.Text = cellText;
            cellShape.TextFrame.TextRange.Font.Size = 14;
            cellShape.TextFrame.VerticalAnchor = 3;
            if column == 1
                cellShape.TextFrame.TextRange.ParagraphFormat.Alignment = 1;
            else
                cellShape.TextFrame.TextRange.ParagraphFormat.Alignment = 2;
            end
            if row == 1
                cellShape.TextFrame.TextRange.Font.Bold = 1;
                cellShape.Fill.ForeColor.RGB = 15132390;
            end
        end
    end

    noteBox = invoke(shapes,'AddTextbox', ...
        1,55,455,slideW-110,45);
    noteBox.TextFrame.TextRange.Text = [ ...
        'Eligibility requires complete P, E_p, and T, plus at least 5% ' ...
        'finite, nonconstant discharge in both periods.'];
    noteBox.TextFrame.TextRange.Font.Size = 13;
    noteBox.TextFrame.TextRange.Font.Italic = 1;
end

function [imgW,imgH] = local_graphic_size(imgFile)
%LOCAL_GRAPHIC_SIZE Return raster or SVG dimensions for slide placement.

    imgW = 1000;
    imgH = 700;

    [~,~,ext] = fileparts(char(string(imgFile)));

    if strcmpi(ext,'.svg')
        try
            txt = fileread(imgFile);
            tok = regexp(txt, ...
                'viewBox\s*=\s*["'']\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*["'']', ...
                'tokens','once');
            if ~isempty(tok)
                imgW = str2double(tok{1});
                imgH = str2double(tok{2});
                if isfinite(imgW) && isfinite(imgH) ...
                        && imgW > 0 && imgH > 0
                    return
                end
            end
        catch
        end

        try
            txt = fileread(imgFile);
            wtok = regexp(txt, ...
                '<svg[^>]*\swidth\s*=\s*["'']\s*([-+0-9.eE]+)', ...
                'tokens','once');
            htok = regexp(txt, ...
                '<svg[^>]*\sheight\s*=\s*["'']\s*([-+0-9.eE]+)', ...
                'tokens','once');
            if ~isempty(wtok) && ~isempty(htok)
                imgW = str2double(wtok{1});
                imgH = str2double(htok{1});
                if isfinite(imgW) && isfinite(imgH) ...
                        && imgW > 0 && imgH > 0
                    return
                end
            end
        catch
        end
        return
    end

    try
        info = imfinfo(imgFile);
        imgW = info.Width;
        imgH = info.Height;
    catch
    end
end

function local_add_image_slide(pres,ppLayoutBlank,imgFile,margin)
%LOCAL_ADD_IMAGE_SLIDE Add one figure image as a centered PPT slide.

    if nargin < 4 ...
            || isempty(margin)
        margin = 20;
    end
    imgFile = char(string(imgFile));
    
    if ~isfile(imgFile)
        error(['PNG file not ' ...
            'found: %s'],imgFile);
    end
    
    try
        slideW = get(pres.PageSetup, ...
            'SlideWidth');
        slideH = get(pres.PageSetup, ...
            'SlideHeight');
    catch ME
        error(['PowerPoint connection ' ...
            'lost while reading ' ...
            'slide size: %s'], ...
            ME.message);
    end
    
    slides = get(pres,'Slides');
    slideCount = get(slides,'Count');
    slide = invoke(slides,'Add', ...
        slideCount+1,ppLayoutBlank);
    
    % slideW = get(pres.PageSetup, ...
    %     'SlideWidth');
    % slideH = get(pres.PageSetup, ...
    %     'SlideHeight');
    
    [imgW,imgH] = local_graphic_size(imgFile);

    scale = min((slideW-2*margin)/imgW, ...
        (slideH-2*margin)/imgH);
    
    w = imgW * scale;
    h = imgH * scale;
    left = (slideW - w)/2;
    top = (slideH - h)/2;
    
    shapes = get(slide,'Shapes');
    invoke(shapes,'AddPicture', ...
        imgFile,0,1,left,top,w,h);
end

function dtName = local_dt_name(C)
%LOCAL_DT_NAME Return text label for model time step.

    try
        switch C.prd.dt
            case 1
                dtName = 'daily';
            case 24
                dtName = 'hourly';
            case 96
                dtName = '15min';
            otherwise
                dtName = sprintf('dt%g',C.prd.dt);
        end
    catch
        dtName = 'n/a';
    end
end

function lossStr = local_loss_string(C)

    if ~isfield(C,'loss') || ~isfield(C.loss,'fnc')
        lossStr = 'n/a';
        return
    end

    lossStr = local_loss_name(double(C.loss.fnc));

end

function imgFile = local_find_model_diagram(C,preferSvg)
%LOCAL_FIND_MODEL_DIAGRAM Find the SVG/PNG diagram for the selected model.

    imgFile = '';

    if nargin < 2
        preferSvg = true;
    end

    modelName = lower(strtrim(char(string(local_model_name(C)))));
    modelName = regexprep(modelName,'[^a-z0-9_]+','');

    switch modelName
        case {'hymod'}
            svgName = 'hymod_snow.svg';
            pngName = 'hymod.png';
        case {'hmodel'}
            svgName = 'hmodel_snow.svg';
            pngName = 'hmodel.png';
        case {'sacsma','sac_sma','sacramento'}
            svgName = 'sacsma_snow.svg';
            pngName = 'sacsma.png';
        case {'xinanjiang'}
            svgName = 'Xinanjiang_snow.svg';
            pngName = 'Xinanjiang.png';
        case {'gr4ja','gr4j'}
            svgName = 'gr4j_snow.svg';
            pngName = 'gr4j.png';
        case {'hbv'}
            svgName = 'hbv_snow.svg';
            pngName = 'hbv.png';
        case {'cfe_nwm','cfenwm'}
            svgName = 'cfe_nwm.svg';
            pngName = 'cfe_nwm.png';
        case {'user_model','usermodel'}
            svgName = 'user_model.svg';
            pngName = 'user_model.png';
        otherwise
            if strcmp(modelName,'na') || isempty(modelName)
                return
            end
            svgName = [modelName '.svg'];
            pngName = [modelName '.png'];
    end

    docs = local_doc_dir_candidates(C);

    candidates = strings(0,1);
    for i = 1:numel(docs)
        if preferSvg
            candidates(end+1) = fullfile(docs(i),svgName); %#ok<AGROW>
            candidates(end+1) = fullfile(docs(i),pngName); %#ok<AGROW>
        else
            candidates(end+1) = fullfile(docs(i),pngName); %#ok<AGROW>
        end
    end

    for i = 1:numel(candidates)
        f = char(candidates(i));
        if isfile(f)
            imgFile = f;
            return
        end
    end
end

function docs = local_doc_dir_candidates(C)
%LOCAL_DOC_DIR_CANDIDATES Return possible packaged/development docs folders.
%
% In a normal MATLAB session, documentation images generally reside in:
%
%   SAGEhydrology/docs
%
% In a deployed MATLAB Compiler application, folders passed through
% AdditionalFiles are extracted beneath ctfroot. Depending on how the
% source hierarchy was packaged, the directory may appear as either:
%
%   ctfroot/docs
%
% or:
%
%   ctfroot/SAGEhydrology/docs

    docs = strings(0,1);

    % --------------------------------
    % MATLAB Compiler deployed paths
    % --------------------------------
    if isdeployed
        try
            docs(end+1) = string( ...
                fullfile(ctfroot,'docs'));
        catch
        end

        try
            docs(end+1) = string(fullfile( ...
                ctfroot,'SAGEhydrology','docs'));
        catch
        end
    end

    % --------------------------------
    % Paths stored in configuration C
    % --------------------------------
    try
        if isstruct(C) ...
                && isfield(C,'SAGEhydro') ...
                && ~isempty(C.SAGEhydro)

            docs(end+1) = string(fullfile( ...
                char(string(C.SAGEhydro)), ...
                'docs'));
        end
    catch
    end

    try
        if isstruct(C) ...
                && isfield(C,'cfg') ...
                && isstruct(C.cfg) ...
                && isfield(C.cfg,'SAGEhydro') ...
                && ~isempty(C.cfg.SAGEhydro)

            docs(end+1) = string(fullfile( ...
                char(string(C.cfg.SAGEhydro)), ...
                'docs'));
        end
    catch
    end

    try
        if isstruct(C) ...
                && isfield(C,'root') ...
                && ~isempty(C.root)

            docs(end+1) = string(fullfile( ...
                char(string(C.root)), ...
                'SAGEhydrology','docs'));
        end
    catch
    end

    try
        if isstruct(C) ...
                && isfield(C,'cfg') ...
                && isstruct(C.cfg) ...
                && isfield(C.cfg,'root') ...
                && ~isempty(C.cfg.root)

            docs(end+1) = string(fullfile( ...
                char(string(C.cfg.root)), ...
                'SAGEhydrology','docs'));
        end
    catch
    end

    % --------------------------------
    % Paths relative to this function
    % --------------------------------
    try
        here = fileparts(mfilename('fullpath'));

        docs(end+1) = string(fullfile( ...
            here,'docs'));

        docs(end+1) = string(fullfile( ...
            here,'..','docs'));

        docs(end+1) = string(fullfile( ...
            here,'..','..','docs'));
    catch
    end

    % Normalize paths and remove duplicates. Force column orientation
    % because unique/logical indexing can otherwise return differently
    % oriented string arrays in deployed MATLAB.
    docs = docs(strlength(docs) > 0);
    docs = docs(:);
    docs = unique(docs,'stable');
    docs = docs(:);

    % Prefer directories that actually exist, while retaining nonexistent
    % candidates at the end for diagnostic consistency.
    exists = false(numel(docs),1);

    for i = 1:numel(docs)
        try
            exists(i) = isfolder(char(docs(i)));
        catch
            exists(i) = false;
        end
    end

    docsExisting = docs(exists);
    docsMissing = docs(~exists);

    docs = [docsExisting(:); docsMissing(:)];
end

function imgFile = local_find_sage_schematic(C,preferSvg)

    imgFile = '';

    if nargin < 2
        preferSvg = true;
    end

    docs = local_doc_dir_candidates(C);
    if preferSvg
        names = {'SAGE_front_software.svg', ...
            'SAGE_front_software.png', ...
            'SAGE_schematic.png'};
    else
        names = {'SAGE_front_software.png', ...
            'SAGE_schematic.png'};
    end

    candidates = strings(0,1);
    for i = 1:numel(docs)
        for j = 1:numel(names)
            candidates(end+1) = fullfile(docs(i),names{j}); %#ok<AGROW>
        end
    end

    for i = 1:numel(candidates)
        f = char(candidates(i));
        if isfile(f)
            imgFile = f;
            return
        end
    end

end

function name = local_loss_name(fnc)
%LOCAL_LOSS_NAME Convert numeric loss-function code to name.
    
    switch double(fnc)
        case 1
            name = 'SAR';
        case 2
            name = 'GLS';
        case 3
            name = 'NSE';
        case 4
            name = 'KGE';
        case 5
            name = 'Huber';
        case 6
            name = 'FDC';
        case 7
            name = 'JKGE';
        otherwise
            name = sprintf('loss %g',fnc);
    end
end

function s = local_fmtDMY(d)
%LOCAL_FMTDMY Format [day month year] vector as dd/mm/yyyy.

    if isempty(d) ...
            || numel(d) < 3
        s = 'n/a';
    else
        d = double(d(:));
        s = sprintf('%02d/%02d/%04d', ...
            d(1),d(2),d(3));
    end
end

function v = local_getfield_or_nan(S,field)
%LOCAL_GETFIELD_OR_NAN Return structure field value or NaN.

    if isstruct(S) ...
            && isfield(S,field) ...
            && ~isempty(S.(field))
        v = S.(field);
    else
        v = NaN;
    end
end

function outDir = local_results_dir(C)
%LOCAL_RESULTS_DIR Return SAGE results directory or current folder.

    if isfield(C,'dirres') ...
            && ~isempty(C.dirres)
        outDir = C.dirres;
    else
        outDir = pwd;
    end
end

function s = local_int_string(x)
%LOCAL_INT_STRING Print integer-like values without decimals.
    
    if isempty(x) ...
            || any(isnan(x))
        s = 'n/a';
    else
        s = sprintf('%d',round(double(x(1))));
    end
end

function s = local_attribute_string(C)
%LOCAL_ATTRIBUTE_STRING Summarize attributes used.

    s = 'n/a';
    
    try
        if isfield(C,'bas') ...
                && isfield(C.bas,'id_attr') ...
                && ~isempty(C.bas.id_attr)
            id_attr = C.bas.id_attr(:).';
            s = sprintf('%d attributes: %s', ...
                numel(id_attr), ...
                char(strjoin(string(id_attr),', ')));
            return
        end
    
    catch
        s = 'n/a';
    end
end

function s = local_meteo_pet_string(C)
%LOCAL_METEO_PET_STRING Return GUI-provided forcing summary for export.
%
% The GUI is the single source of truth for meteorological display text.
% It stores the selected product, precipitation, temperature, PET, and a
% ready-to-use summary in C.misc.meteo.display when building the run
% configuration. save_to_pptx therefore does not maintain region-specific
% meteorological lookup tables.

    s = 'n/a';

    try
        if ~isfield(C,'misc') ...
                || ~isfield(C.misc,'meteo') ...
                || ~isfield(C.misc.meteo,'display') ...
                || ~isstruct(C.misc.meteo.display)
            return
        end

        D = C.misc.meteo.display;

        if isfield(D,'summary') ...
                && strlength(strtrim(string(D.summary))) > 0
            s = char(strtrim(string(D.summary)));
            return
        end

        product = local_display_field(D,'product');
        precip = local_display_field(D,'precip');
        temp = local_display_field(D,'temp');
        pet = local_display_field(D,'pet');

        s = sprintf( ...
            'product=%s; precip=%s; temp=%s; PET=%s', ...
            product,precip,temp,pet);

    catch
        s = 'n/a';
    end
end

function s = local_region_code(C)
%LOCAL_REGION_CODE Return short CAMELS region code.

    s = 'REG';
    try
        s = region_helpers('short',C);
    catch
        try
            s = regexprep(upper(char(string( ...
                region_helpers('code',C)))), ...
                '^CAMELS[_-]','');
        catch
        end
    end
end

function v = local_get_first_existing(S,fields)
%LOCAL_GET_FIRST_EXISTING Return first existing nonempty field.

    v = [];    
    for i = 1:numel(fields)
        f = fields{i};
        if isstruct(S) ...
                && isfield(S,f) ...
                && ~isempty(S.(f))
            v = S.(f);
            return
        end
    end
end

function s = local_ffn_string(C)
%LOCAL_FFN_STRING Summarize FFN architecture.

    try
        ann = C.net;
    
        h = local_hidden_vector(ann);
        nout = local_output_size(C);
    
        if isempty(h)
            s = sprintf('output layer: %d',nout);
        elseif isscalar(h)
            s = sprintf(['hidden layer: %d, ' ...
                'output layer: %d'],h(1),nout);
        else
            parts = strings(1,numel(h));
            for i = 1:numel(h)
                parts(i) = sprintf(['hidden ' ...
                    'layer %d: %d'],i,h(i));
            end
            s = sprintf('%s, output layer: %d', ...
                char(strjoin(parts,', ')),nout);
        end
    catch
        s = 'n/a';
    end
end

function s = local_hidden_layers_string(C)
%LOCAL_HIDDEN_LAYERS_STRING Number of hidden layers.

    try
        tf = local_ann_tf(C.net);
        if ~isempty(tf)
            s = sprintf('%d',numel(tf));
        else
            h = local_hidden_vector(C.net);
            s = sprintf('%d',numel(h));
        end
    catch
        s = 'n/a';
    end
end

function s = local_transfer_string(C)
%LOCAL_TRANSFER_STRING Summarize transfer functions.

    try
        ann = C.net;
    
        tf = local_get_first_existing(ann, ...
            {'transfer','tf','act','activation'});
    
        if isempty(tf)
            s = 'n/a';
            return
        end
    
        tf = string(tf(:)).';
    
        if isscalar(tf)
            s = sprintf('transfer function: %s',tf);
        else
            parts = strings(1,numel(tf));
            for i = 1:numel(tf)
                parts(i) = sprintf(['transfer ' ...
                    'function %d: %s'],i,tf(i));
            end
            s = char(strjoin(parts,', '));
        end
    catch
        s = 'n/a';
    end
end

function s = local_nffn_parameters_string(C)
%LOCAL_NFFN_PARAMETERS_STRING Number of FFN weights and biases.
    
    try
        if isfield(C.net,'l') && ~isempty(C.net.l)
            s = sprintf('%d',round( ...
                double(C.net.l)));
            return
        end
    
        nin = local_input_size(C);
        h = local_hidden_vector(C.net);
        nout = local_output_size(C);
    
        dims = [nin h nout];
    
        npar = 0;
        for j = 1:numel(dims)-1
            npar = npar + dims(j)*dims(j+1) + dims(j+1);
        end
    
        s = sprintf('%d',round(npar));
    catch
        s = 'n/a';
    end
end

function s = local_optimizer_string(C)

    s = 'n/a';
    
    try
        alg = C.alg;
    
        if isfield(alg,'method') ...
                && ~isempty(alg.method)
            switch double(alg.method)
                case 1
                    s = 'Gradient descent';
                case 2
                    s = 'Adam/AdamW';
                otherwise
                    s = sprintf('algorithm %g', ...
                        double(alg.method));
            end
        end
    catch
        s = 'n/a';
    end
end

function s = local_alg_field_string(C,field)
%LOCAL_ALG_FIELD_STRING Return optimizer scalar field as text.

    s = 'n/a';
    
    try
        if isfield(C,'alg') && isfield(C.alg,field) ...
                && ~isempty(C.alg.(field))
            v = C.alg.(field);
            if isnumeric(v)
                s = sprintf('%.6g',double(v(1)));
            else
                s = char(string(v));
            end
        end
    catch
        s = 'n/a';
    end
end

function s = local_iterations_done_string(C)

    s = 'n/a';
    
    try
        if isfield(C,'runtime') ...
                && isfield(C.runtime,'lastIter') ...
                && ~isempty(C.runtime.lastIter)
            s = sprintf('%d',round(double(C.runtime.lastIter)));
        elseif isfield(C,'lastIter') ...
                && ~isempty(C.lastIter)
            s = sprintf('%d',round(double(C.lastIter)));
        elseif isfield(C,'prf') ...
                && isfield(C.prf,'iter') ...
                && isfield(C.prf.iter,'cpuT') ...
                && ~isempty(C.prf.iter.cpuT)
            % A finite timing entry is written only for a completed
            % iteration; trailing preallocated entries remain NaN.
            nDone = find(isfinite(double(C.prf.iter.cpuT(:))),1,'last');
            if ~isempty(nDone)
                s = sprintf('%d',nDone);
            end
        end

        % An exported script calls save_to_pptx immediately after the
        % optimization loop, so i is the completed iteration count.
        if strcmp(s,'n/a') ...
                && evalin('base','exist(''i'',''var'') == 1')
            nDone = evalin('base','i');
            if isnumeric(nDone) && isscalar(nDone) ...
                    && isfinite(double(nDone)) && double(nDone) >= 0
                s = sprintf('%d',round(double(nDone)));
            end
        end

        % A timing profile is a safer fallback if the loop variable was
        % cleared or renamed in a hand-edited script.
        if strcmp(s,'n/a') ...
                && evalin('base','exist(''prf'',''var'') == 1')
            profile = evalin('base','prf');
            if isstruct(profile) && isfield(profile,'iter') ...
                    && isfield(profile.iter,'cpuT') ...
                    && ~isempty(profile.iter.cpuT)
                nDone = find(isfinite(double(profile.iter.cpuT(:))),1,'last');
                if ~isempty(nDone)
                    s = sprintf('%d',nDone);
                end
            end
        end
    catch
        s = 'n/a';
    end
end

function s = local_total_runtime_string(C)
%LOCAL_TOTAL_RUNTIME_STRING Return elapsed run time as hh:mm:ss when available.

    s = 'n/a';
    
    try
        t0 = [];
    
        if isfield(C,'runtime') ...
                && isfield(C.runtime,'startTime') ...
                && ~isempty(C.runtime.startTime)
            t0 = C.runtime.startTime;
        elseif isfield(C,'startTime') ...
                && ~isempty(C.startTime)
            t0 = C.startTime;
        end
    
        if isempty(t0) ...
                && evalin('base', ...
                ['exist(''SAGEnew_lastOut'',''var'') ...' ...
                '&& isfield(SAGEnew_lastOut,''runtime'')'])
            rt = evalin('base','SAGEnew_lastOut.runtime');
            if isstruct(rt) ...
                    && isfield(rt,'startTime')
                t0 = rt.startTime;
            end
        end
    
        if isempty(t0)
            return
        end
    
        if isa(t0,'datetime')
            dt = datetime('now') - t0;
        elseif isnumeric(t0)
            dt = seconds(toc(t0));
        else
            try
                dt = datetime('now') - datetime(t0);
            catch
                return
            end
        end
    
        try
            dt.Format = 'hh:mm:ss';
            s = char(dt);
        catch
            s = char(string(dt));
        end
    catch
        s = 'n/a';
    end
end

function note = local_run_note_string(C)
%LOCAL_RUN_NOTE_STRING Return run note from C.runNote or C.misc.runNote.

    note = "";
    
    try
        if isfield(C,'runNote') ...
                && ~isempty(C.runNote)
            note = string(C.runNote);
        elseif isfield(C,'misc') ...
                && isfield(C.misc,'runNote') ...
                && ~isempty(C.misc.runNote)
            note = string(C.misc.runNote);
        end
    
        note = strjoin(note(:),newline);
    catch
        note = "";
    end
end

function nin = local_input_size(C)
%LOCAL_INPUT_SIZE Determine FFN input size.

    nin = [];
    
    try
        if isfield(C,'net') && isfield(C.net,'r') ...
                && ~isempty(C.net.r)
            nin = double(C.net.r(1));
            return
        end
    
        if isfield(C,'bas') ...
                && isfield(C.bas,'id_attr') ...
                && ~isempty(C.bas.id_attr)
            nin = numel(C.bas.id_attr);
            return
        end
    catch
        nin = [];
    end
end

function v = local_get_numeric(S,field)
%LOCAL_GET_NUMERIC Safely return numeric scalar field.

    v = [];
    
    try
        if isstruct(S) ...
                && isfield(S,field) ...
                && ~isempty(S.(field))
            v = double(S.(field));
            v = v(1);
        end
    catch
        v = [];
    end
end

function h = local_hidden_vector(ann)
%LOCAL_HIDDEN_VECTOR Extract hidden layer widths from net.h.

    h = [];
    
    try
        if isfield(ann,'h') ...
                && ~isempty(ann.h)
            h = ann.h;
    
            while iscell(h) ...
                    && isscalar(h)
                h = h{1};
            end
    
            if isnumeric(h)
                h = double(h(:).');
            elseif isstring(h)
                if isscalar(h)
                    h = str2num(char(h)); %#ok<ST2NM>
                else
                    h = str2double(h(:).');
                end
            elseif ischar(h)
                h = str2num(h); %#ok<ST2NM>
            elseif iscell(h)
                hh = nan(1,numel(h));
                for i = 1:numel(h)
                    hh(i) = str2double(string(h{i}));
                end
                h = hh;
            end
    
            h = double(h(:).');
        end
    catch
        h = [];
    end
end

function nout = local_output_size(C)
%LOCAL_OUTPUT_SIZE Determine output layer size.

    nout = [];
    
    try
        if isfield(C.mdl,'npar') ...
                && ~isempty(C.mdl.npar)
            nout = double(C.mdl.npar);
        elseif isfield(C.mdl,'par_names') ...
                && ~isempty(C.mdl.par_names)
            nout = numel(C.mdl.par_names);
        elseif isfield(C,'net') && isfield(C.net,'nout') ...
                && ~isempty(C.net.nout)
            nout = double(C.net.nout);
        end
    
        if isempty(nout)
            nout = NaN;
        else
            nout = nout(1);
        end
    catch
        nout = NaN;
    end
end

function s = local_fig_name(f)
%LOCAL_FIG_NAME Return figure name for logging.

    s = "unnamed figure";
    
    try
        if isprop(f,'Name') ...
                && strlength(string(f.Name)) > 0
            s = string(f.Name);
        elseif isprop(f,'Number')
            s = "Figure " + string(f.Number);
        end
    catch
    end
    
    s = char(s);
end

function s = local_model_name(C)
%LOCAL_MODEL_NAME Return current model name.

    s = 'n/a';
    
    try
        if isfield(C,'mdl') ...
                && isfield(C.mdl,'names') ...
                && isfield(C.mdl,'model') ...
                && ~isempty(C.mdl.names)
            s = C.mdl.names{C.mdl.model};
        end
    catch
        try
            s = C.mdl.names(C.mdl.model);
            s = char(string(s));
        catch
            s = 'n/a';
        end
    end
end

function s = local_region_string(C)

    try
        s = region_helpers('name',C);
        return
    catch
    end

    s = 'unknown';

end

function s = local_split_method(C)
%LOCAL_SPLIT_METHOD Return split method string.

    s = 'n/a';
    
    try
        if isfield(C,'prd') ...
                && isfield(C.prd,'method') ...
                && ~isempty(C.prd.method)
            s = char(string(C.prd.method));
        elseif isfield(C,'mdl') ...
                && isfield(C.mdl,'sp_method') ...
                && ~isempty(C.mdl.sp_method)
            s = char(string(C.mdl.sp_method));
        end
    catch
        s = 'n/a';
    end
end

function tf = local_ann_tf(ann)
%LOCAL_ANN_TF Return transfer functions as row string array.

    tf = strings(1,0);
    
    try
        if ~isfield(ann,'tf') ...
                || isempty(ann.tf)
            return
        end
    
        x = ann.tf;
        while iscell(x) ...
                && isscalar(x)
            x = x{1};
        end
    
        if ischar(x) ...
                || isstring(x)
            tf = string(x);
            tf = tf(:).';
        elseif iscell(x)
            tf = string(x(:)).';
        end
    catch
        tf = strings(1,0);
    end
end

function s = local_training_period_string(C)

    s = '-';
    
    try
        P = C.prd;
        M = C.mdl;
    
        if isfield(P,'dts') ...
                && isfield(P,'dte') ...
                && ~isempty(P.dts) ...
                && ~isempty(P.dte)
            s = [local_fmtDMY(P.dts) ' - ' local_fmtDMY(P.dte)];
        elseif isfield(M,'dts') ...
                && isfield(M,'dte') ...
                && ~isempty(M.dts) ...
                && ~isempty(M.dte)
            s = [local_fmtDMY(M.dts) ' - ' local_fmtDMY(M.dte)];
        elseif isfield(P,'ds') ...
                && isfield(P,'de') ...
                && ~isempty(P.ds) ...
                && ~isempty(P.de)
            s = [local_fmtDMY(P.ds) ' - ' local_fmtDMY(P.de)];
        elseif isfield(M,'ds') ...
                && isfield(M,'de') ...
                && ~isempty(M.ds) ...
                && ~isempty(M.de)
            s = [local_fmtDMY(M.ds) ' - ' local_fmtDMY(M.de)];
        end
    catch
    end
end

function s = local_evaluation_period_string(C)

    s = '-';
    
    try
        P = C.prd;
        M = C.mdl;
    
        if isfield(P,'des') ...
                && isfield(P,'dee') ...
                && ~isempty(P.des) ...
                && ~isempty(P.dee)
            s = [local_fmtDMY(P.des) ' - ' local_fmtDMY(P.dee)];
        elseif isfield(M,'des') ...
                && isfield(M,'dee') ...
                && ~isempty(M.des) ...
                && ~isempty(M.dee)
            s = [local_fmtDMY(M.des) ' - ' local_fmtDMY(M.dee)];
        elseif isfield(P,'ds') ...
                && isfield(P,'de') ...
                && ~isempty(P.ds) ...
                && ~isempty(P.de)
            s = [local_fmtDMY(P.ds) ' - ' local_fmtDMY(P.de)];
        elseif isfield(M,'ds') ...
                && isfield(M,'de') ...
                && ~isempty(M.ds) ...
                && ~isempty(M.de)
            s = [local_fmtDMY(M.ds) ' - ' local_fmtDMY(M.de)];
        end
    catch
    end
end

function s = local_spinup_string(C)

    s = '-';
    
    try
        if isfield(C,'prd') ...
                && isfield(C.prd,'spinup') ...
                && ~isempty(C.prd.spinup)
            s = sprintf('%g days',double(C.prd.spinup));
        end
    catch
    end
end

function local_close_export_figures()
%LOCAL_CLOSE_EXPORT_FIGURES Close remaining MATLAB plot figures, not GUI.

    figs = findall(groot,'Type','figure');
    
    for i = 1:numel(figs)
        f = figs(i);
    
        try
            if ~isvalid(f)
                continue
            end
    
            % Do not close App Designer/uifigure GUI windows.
            kids = findall(f);
            cls = string(arrayfun(@(h) class(h),kids, ...
                'UniformOutput',false));
    
            isGUI = any(contains(cls, ...
                ["matlab.ui.control.HTML", ...
                "matlab.ui.container.GridLayout", ...
                "matlab.ui.container.TabGroup", ...
                "matlab.ui.container.Tab", ...
                "matlab.ui.control.Table", ...
                "matlab.ui.control.Button", ...
                "matlab.ui.control.DropDown", ...
                "matlab.ui.control.EditField"]));
    
            if isGUI
                continue
            end
    
            close(f)
    
        catch
        end
    end
end

function s = local_split_method_string(C)

    method = local_split_method(C);

    switch lower(method)
        case 'manual'
            s = 'manual dates';
        case 'deterministic_block'
            s = 'deterministic block';
        case 'random_block'
            s = 'random block';
        case 'random'
            s = 'random';
        case 'deterministic_kfold'
            s = 'deterministic k-fold';
        case 'random_kfold'
            s = 'random k-fold';
        case 'rainfall_block'
            s = 'rainfall block';
        otherwise
            s = method;
    end

end

function s = local_basin_sample_string(C)
    s = 'n/a';
    
    try
        sample = char(string(C.bas.sample));
    
        switch lower(sample)
            case 'file'
                splitFile = local_basin_split_file_string(C);
                s = sprintf('user file; split file = %s',splitFile);

            case 'random'
                seed = local_basin_numeric_field(C,'seed',1);
                generator = local_basin_text_field( ...
                    C,'rng_type','twister');
                splitFile = local_basin_split_file_string(C);
                s = sprintf([ ...
                    'random; seed = %d; generator = %s; ' ...
                    'split file = %s'], ...
                    round(seed),generator,splitFile);

            case 'zone'
                minZone = local_basin_numeric_field( ...
                    C,'min_per_zone',NaN);
                maxZones = local_basin_numeric_field( ...
                    C,'max_zones',NaN);
                s = sprintf([ ...
                    'zone-stratified; min/zone = %s; ' ...
                    'max zones = %s'], ...
                    local_optional_integer(minZone), ...
                    local_optional_integer(maxZones));

            otherwise
                s = sample;
        end
    catch
    end
end

function v = local_basin_numeric_field(C,name,defaultValue)
    v = defaultValue;
    try
        if isfield(C,'bas') ...
                && isfield(C.bas,name) ...
                && ~isempty(C.bas.(name))
            x = double(C.bas.(name));
            if isfinite(x(1))
                v = x(1);
            end
        end
    catch
    end
end

function s = local_basin_text_field(C,name,defaultValue)
    s = defaultValue;
    try
        if isfield(C,'bas') ...
                && isfield(C.bas,name) ...
                && ~isempty(C.bas.(name))
            s = char(string(C.bas.(name)));
        end
    catch
    end
end

function s = local_basin_split_file_string(C)
    s = 'n/a';
    try
        if isfield(C,'file_split') && ~isempty(C.file_split)
            s = char(string(C.file_split));
        elseif isfield(C,'bas') ...
                && isfield(C.bas,'file_split') ...
                && ~isempty(C.bas.file_split)
            s = char(string(C.bas.file_split));
        end
    catch
    end
end

function s = local_optional_integer(x)
    if isempty(x) || ~isfinite(double(x(1)))
        s = 'n/a';
    else
        s = sprintf('%d',round(double(x(1))));
    end
end

function s = local_jkge_settings_string(C)
    s = 'n/a';
    
    try
        if ~isfield(C,'loss') || C.loss.fnc ~= 7
            return
        end
    
        op = local_jkge_operator_name(C.loss.method);
    
        if isfield(C.loss,'n_win') && ~isempty(C.loss.n_win)
            win = sprintf('%g days',C.loss.n_win);
        else
            win = 'n/a';
        end
    
        if isfield(C.loss,'M') && ~isempty(C.loss.M)
            M = round(double(C.loss.M));
            switch M
                case 1
                    mdef = '1 paper';
                case 2
                    mdef = '2 norm-based';
                otherwise
                    mdef = sprintf('%g',M);
            end
        else
            mdef = 'n/a';
        end
    
        s = sprintf('%s; window = %s; M definition = %s', ...
            op,win,mdef);
    catch
        s = 'n/a';
    end
end

function s = local_jkge_operator_name(method)
    switch round(double(method))
        case 1
            s = '1 moving average-mean';
        case 2
            s = '2 section-wise mean';
        case 3
            s = '3 long-term mean';
        case 4
            s = '4 climatological mean';
        otherwise
            s = 'n/a';
    end
end

function txt = local_title_text(C,modelName,dtName,lossStr, ...
    regionName,nTrain,nEval,createdStr, ...
    trainPeriodStr,evalPeriodStr,spinupStr, ...
    splitStr,sampleStr)

    jkgeLine = '';
    jkgeArg = {};
    if isfield(C,'loss') && isfield(C.loss,'fnc') ...
            && C.loss.fnc == 7
        jkgeLine = 'JKGE settings: %s\n';
        jkgeArg = {local_jkge_settings_string(C)};
    end

    txt = sprintf([ ...
        'Model: %s\n' ...
        'Region: %s\n' ...
        'Resolution: %s\n' ...
        'Training period: %s\n' ...
        'Evaluation period: %s\n' ...
        'Spin-up: %s\n' ...
        'Split method: %s\n' ...
        'Training basins: %s\n' ...
        'Evaluation basins: %s\n' ...
        'Basin sampling: %s\n' ...
        'Loss function: %s\n' ...
        jkgeLine ...
        'Attributes: %s\n' ...
        'Meteo: %s\n' ...
        'FFN architecture: %s\n' ...
        '# hidden layers: %s\n' ...
        'Transfer function(s): %s\n' ...
        '# FFN parameters: %s\n' ...
        'Optimization method: %s\n' ...
        'Max iterations: %s\n' ...
        'Iterations done: %s\n' ...
        'Gradient clipping: %s\n' ...
        'Weight decay: %s\n' ...
        'Learning rate: %s\n' ...
        'Data loading: %s\n' ...
        'Initialization: %s\n' ...
        'SAGE training: %s\n' ...
        'Information theory: %s\n' ...
        'Postprocessing: %s\n' ...
        'GUI overhead: %s\n' ...
        'Total elapsed: %s\n' ...
        'Results directory: %s\n' ...
        'Created: %s'], ...
        modelName,regionName,dtName, ...
        trainPeriodStr, ...
        evalPeriodStr, ...
        spinupStr, ...
        splitStr, ...
        local_int_string(nTrain), ...
        local_int_string(nEval), ...
        sampleStr, ...
        lossStr, ...
        jkgeArg{:}, ...
        local_attribute_string(C), ...
        local_meteo_pet_string(C), ...
        local_ffn_string(C), ...
        local_hidden_layers_string(C), ...
        local_transfer_string(C), ...
        local_nffn_parameters_string(C), ...
        local_optimizer_string(C), ...
        local_alg_field_string(C,'i_max'), ...
        local_iterations_done_string(C), ...
        local_alg_field_string(C,'clipn'), ...
        local_alg_field_string(C,'wdecay'), ...
        local_alg_field_string(C,'lr'), ...
        local_phase_runtime_string(C,'data'), ...
        local_phase_runtime_string(C,'initialization'), ...
        local_phase_runtime_string(C,'training'), ...
        local_phase_runtime_string(C,'information'), ...
        local_phase_runtime_string(C,'postprocessing'), ...
        local_phase_runtime_string(C,'gui'), ...
        local_phase_runtime_string(C,'total'), ...
        local_results_dir(C),createdStr);
end

function s = local_phase_runtime_string(C,name)
%LOCAL_PHASE_RUNTIME_STRING Format one phase wall time stored in seconds.
    s = 'n/a';
    try
        if isfield(C,'runtime') && isfield(C.runtime,'timing') ...
                && isfield(C.runtime.timing,name)
            value = double(C.runtime.timing.(name));
            if isscalar(value) && isfinite(value) && value >= 0
                s = sprintf('%.1f s',value);
            end
        end
    catch
    end
end

function s = local_display_field(D,name)
%LOCAL_DISPLAY_FIELD Read one presentation metadata field safely.

    s = 'n/a';
    try
        if isfield(D,name) && ~isempty(D.(name))
            s = char(string(D.(name)));
        end
    catch
        s = 'n/a';
    end
end
