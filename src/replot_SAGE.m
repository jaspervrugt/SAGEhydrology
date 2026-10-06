function output = replot_SAGE(matFile,outputFile)
%REPLOT_SAGE Redraw saved results without running training or the model.
% replot_SAGE('results/SAGE_export_....mat','results/review.pptx')
% Older exports use the companion SAGEnew_lastOut.mat for ECDF/time series.
    root=fileparts(fileparts(mfilename('fullpath')));
    if ~isdeployed, addpath(genpath(root)); end
    if nargin<2, [folder,name]=fileparts(matFile); outputFile=fullfile(folder,[name '_replot.pptx']); end
    outputFile=char(outputFile);
    if isempty(regexp(outputFile,'^([A-Za-z]:[\\/]|[\\/])','once'))
        outputFile=fullfile(pwd,outputFile);
    end
    saved=load(matFile); E=saved.E;
    prf=struct('curr',E.performance.current,'iter',E.performance.iteration);
    if isfield(E,'figure_inputs')
        P=E.figure_inputs; C=P.config;
        gauges=struct('ecdf_prf',prf,'loss',C.misc.loss, ...
            'observation_out',P.observation_out,'observation_cfg',C);
        theta=E.theta_n{:,2:end}';
        plot_SAGE('sage',E.metadata.mdl,P.dat,E.metadata.bas,E.metadata.prd, ...
            E.discharge.Q_scenario,prf.curr,P.Qfdc,C.region, ...
            P.tTheta,P.latlon,theta,P.At,P.An,gauges);
    else
        companion=fullfile(fileparts(matFile),'SAGEnew_lastOut.mat');
        assert(isfile(companion),'This older export needs its companion SAGEnew_lastOut.mat.');
        last=load(companion); C=last.cfg_save; out=last.out;
        assert(isequal(out.bas.K_t,E.metadata.bas.K_t) && ...
            isequaln(out.prf.curr.NSE,E.performance.current.NSE), ...
            'The companion checkpoint belongs to a different run.');
        warning('replot_SAGE:LegacyExport', ...
            'Older MAT lacks attribution arrays: redrawing available figures from the matching checkpoint.');
        dat=cell(height(E.basins),1);
        for k=1:numel(dat)
            y=E.discharge.Q_measured(k,:)';
            dat{k}=struct('gauge',char(E.basins.basin_id(k)), ...
                'obs',struct('Q',struct('value',y,'bad',~isfinite(y))), ...
                'id_train',find(E.discharge.split_code(k,:)==1), ...
                'id_eval',find(E.discharge.split_code(k,:)==2));
        end
        for pos=1:numel(out.Qfdc.id)
            k=out.Qfdc.id(pos); R=out.Qfdc.observations{pos};
            if isfield(R,'meteo'), dat{k}.meteo=R.meteo; end
        end
        gauges=struct('ecdf_prf',prf,'loss',C.misc.loss, ...
            'observation_out',out,'observation_cfg',C);
        plot_SAGE('sage',E.metadata.mdl,dat,E.metadata.bas,E.metadata.prd, ...
            E.discharge.Q_scenario,prf.curr,out.Qfdc,C.region, ...
            out.tTheta,out.latlon,E.theta_n{:,2:end}',[],[],gauges);
    end
    C.net=C.misc.net; C.alg=C.misc.alg; C.loss=C.misc.loss;
    output=save_to_pptx(C,outputFile);
end
