% RUN_HYMOD_MOEAD_ADAM One-catchment population-based calibration.
% Run this file from MATLAB. Edit settings below; no SITE GUI changes.
% Each of 21 parameter vectors has a different NSE/S_logp trade-off.
% ONE HYMOD trajectory and sensitivity solve supplies BOTH objective gradients.
% Archive and physical parameters are saved separately from paper/SITE results.
% Gradient correction to the supplied description: multiply softmax by lambda.
% Smooth max uses losses minus the best observed ideal point; no absolute value
% is needed because ideal is updated before comparisons and is componentwise
% no greater than every current or proposed loss.
% References: Zhang & Li (2007), doi:10.1109/TEVC.2007.892759;
% Kingma & Ba (2015), https://arxiv.org/abs/1412.6980.

sageRoot = fileparts(fileparts(mfilename('fullpath')));
root = fileparts(fileparts(sageRoot));
gaugeID = "02472000"; % = Leaf River near Collins, MS
        % "01022500";                 % CAMELS-US daily catchment
populationSize = 31; 
maxIterations = 100; 
neighborhoodSize = 5;
learningRate = 0.01; 
smoothingTau = 0.02; 
randomSeed = 10;
useParallel = false;                  % true parallelizes candidates, not basins
checkGradient = true;                % finite-difference check before optimization
makePlots = true;
resultDir = fullfile(fileparts(sageRoot), ...
    'results','Matlab','moead_hymod');
addpath(fullfile(sageRoot,'utils'));
region = bootstrap_SAGE(root, ...
    'CAMELS_US',resultDir);
dirM = fullfile(root,'Data', ...
    'CAMELS_US','daily', ...
    'v1p2','forcing');
dirQ = fullfile(root,'Data', ...
    'CAMELS_US','daily', ...
    'v1p2','streamflow');
assert(isfolder(dirM) ...
    && isfolder(dirQ), ...
    'Edit dirM/dirQ to your CAMELS-US daily folders.');
prd = struct('dt',1, ...
    'method','manual', ...
    'dts',[1 10 1999], ...
    'dte',[30 9 2008], ...
    'des',[1 10 2008], ...
    'dee',[30 9 2013], ...
    'spinup',365);
bas = struct('K',1, ...
    'K_t',1, ...
    'K_e',0, ...
    'id_t',1, ...
    'id_e',[], ...
    'id_gauge',gaugeID, ...
    'dt',1, ...
    'resolution','Daily', ...
    'stream','daily');
mdl = struct('model',1, ...
    'mcode',4, ...
    'calc','seq', ...
    'pspace',1);
ode = read_numsettings(struct());
% Production tolerances: keep the calibration evaluations fast.
ode.RelTol = 1e-3; 
ode.AbsTol = 1e-3;
[mdl,d] = read_model(mdl,prd);
[split,mdl] = build_split(mdl,prd,bas);
meteo = struct('data',3, ...
    'precip',NaN, ...
    'temp',NaN, ...
    'pet',1);
[dat,aux] = read_meteo(region,dirM,bas,split,meteo);
dat=read_Q(region,dirQ,mdl,dat,bas,split,aux);
[eligibility,dat] = check_basins(dat,mdl,bas);
assert(eligibility.valid(1), ...
    'The requested catchment failed data-quality screening.');
loss = struct('fnc',3,'observed',{{'Q'}}, ...
    'fdc',struct('formulation',3));
[dat,loss] = prep_stats(dat,mdl,split,loss);
basinData = dat{1};
objective = @(x)hymodObjectives(x,mdl,basinData,ode,loss);
if checkGradient
    % Verify the derivative implementation using a tightly integrated model.
    % This does not alter the production objective or ode tolerances.
    verificationODE = ode;
    verificationODE.RelTol = 1e-8;
    verificationODE.AbsTol = 1e-8;
    verificationObjective = @(x)hymodObjectives( ...
        x,mdl,basinData,verificationODE,loss);
    x = .5*ones(d,1); 
    [~,J] = verificationObjective(x); 
    fd = zeros(d,2); 
    h = 1e-4;
    for j=1:d
        xp = x; xm = x; xp(j) = xp(j) + h; xm(j) = xm(j) - h;
        fp = verificationObjective(xp);
        fm = verificationObjective(xm);
        fd(j,:) = (fp - fm)'/(2*h);
    end
    relativeError = norm(J-fd,'fro')/max(1,norm(fd,'fro'));
    fprintf(['HYMOD objective ' ...
        'gradient relative error: %.3g\n'],relativeError);
    assert(relativeError<1e-1, ...
        'Gradient check failed; inspect model/solver settings.');
end
runClock = tic;
moead = moead_adam(objective,zeros(d,1),ones(d,1), ...
    population = populationSize, ...
    iterations = maxIterations, ...
    neighbors = neighborhoodSize, ...
    learningRate = learningRate, ...
    tau = smoothingTau, ...
    seed = randomSeed, ...
    parallel = useParallel);
moead.theta = mdl.th_min(:) + (mdl.th_max(:) - mdl.th_min(:)) .* moead.x;
moead.archiveTheta = mdl.th_min(:) + ...
    (mdl.th_max(:) - mdl.th_min(:)) .* moead.archiveX;
moead.gaugeID = gaugeID; 
moead.parameterNames = mdl.par_names;
moead.totalElapsedSeconds = toc(runClock);
paretoTable = array2table([moead.archiveLosses' moead.archiveTheta'], ...
    'VariableNames',cellstr(matlab.lang.makeUniqueStrings( ...
    matlab.lang.makeValidName(['oneMinusNSE';'oneMinusSlogp'; ...
    string(mdl.par_names(:))]))));
stamp = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
outputFile = fullfile(resultDir, ...
    ['HYMOD_MOEAD_' char(gaugeID) '_' stamp '.mat']);
save(outputFile,'moead','paretoTable', ...
    'mdl','ode','prd','bas','loss');
writetable(paretoTable,replace(outputFile,'.mat','.csv'));
fprintf('Saved %d nondominated solutions: %s\n', ...
    height(paretoTable),outputFile);
if makePlots
    figure('Name','HYMOD MOEA/D-Adam Pareto approximation');
    scatter(moead.losses(1,:), ...
        moead.losses(2,:),35,[.7 .7 .7],'filled'); hold on
    scatter(moead.archiveLosses(1,:), ...
        moead.archiveLosses(2,:),45,[0 .45 .74],'filled');
    xlabel('1 - NSE (minimize)'); 
    ylabel('1 - S_{log p} (minimize)'); grid on
    legend('Final population','Nondominated archive', ...
        'Location','best');
    title(sprintf('HYMOD: %s, %d preference vectors', ...
        gaugeID,populationSize));
end

% Plot of nondominated solutions
%% Rainfall, PET, and Pareto discharge envelope
nPareto = size(moead.archiveX,2);
Qpareto = nan(numel(basinData.obs.Q.value),nPareto);
request = crr_request(struct('obs',{{'Q'}}));

for k = 1:nPareto
    [~,sim] = crr_model(moead.archiveX(:,k), ...
        mdl,basinData,ode,loss,request);
    Qpareto(:,k) = sim.obs.Q(:);
end
assert(all(isfinite(Qpareto),'all'), ...
    'A Pareto solution returned invalid discharge.');

Qlo = min(Qpareto,[],2);
Qhi = max(Qpareto,[],2);
Qobs = basinData.obs.Q.value(:);
Qobs(logical(basinData.obs.Q.bad(:))) = NaN;

%% Compare equal-weight hydrograph aggregates with training-selected members
[~,bestNSE] = min(moead.archiveLosses(1,:));
[~,bestSlogp] = min(moead.archiveLosses(2,:));
Qmean = mean(Qpareto,2);
Qmedian = median(Qpareto,2);
comparisonQ = [Qpareto(:,bestNSE),Qpareto(:,bestSlogp),Qmean,Qmedian];
comparisonNames = ["Best archive NSE";"Best archive S_logp"; ...
    "Pareto mean";"Pareto median"];
% RMSSE denominator uses only adjacent, valid TRAINING observations.
% RMSSE=sqrt(MSE / mean((Qobs(t)-Qobs(t-1))^2)); lag=1 day.
trainingIDs=double(basinData.stats.Q.train.indices(:));
trainingMask=false(size(Qobs)); trainingMask(trainingIDs)=true;
adjacent=trainingMask(2:end) & trainingMask(1:end-1) ...
    & isfinite(Qobs(2:end)) & isfinite(Qobs(1:end-1));
differences=diff(double(Qobs));
rmsseDenominator=mean(differences(adjacent).^2);
assert(isfinite(rmsseDenominator) && rmsseDenominator>0, ...
    'RMSSE requires varying, adjacent training observations.');
comparisonTable=table();
for periodName=["train","eval"]
    period=char(periodName);
    if isempty(basinData.stats.Q.(period).indices), continue; end
    for k=1:4
        metrics=hydrographMetrics(comparisonQ(:,k),basinData,loss, ...
            period,rmsseDenominator);
        row=table(string(period),comparisonNames(k),metrics(1), ...
            metrics(2),metrics(3),metrics(4),metrics(5),metrics(6), ...
            'VariableNames',{'Period','Solution','RMSE','RMSSE', ...
            'NSE','KGE','S_logp','PBIAS_percent'});
        comparisonTable=[comparisonTable;row]; %#ok<AGROW>
    end
end
disp(comparisonTable);
% Positive improvement means the aggregate outperforms BOTH selected members
% for that metric in that period. Evaluation members stay training-selected.
improvementTable=table();
for periodName=unique(comparisonTable.Period,'stable')'
    rows=comparisonTable(comparisonTable.Period==periodName,:);
    for k=3:4
        row=table(periodName,rows.Solution(k), ...
            min(rows.RMSE(1:2))-rows.RMSE(k), ...
            min(rows.RMSSE(1:2))-rows.RMSSE(k), ...
            rows.NSE(k)-max(rows.NSE(1:2)), ...
            rows.KGE(k)-max(rows.KGE(1:2)), ...
            rows.S_logp(k)-max(rows.S_logp(1:2)), ...
            'VariableNames',{'Period','Aggregate','RMSE_improvement', ...
            'RMSSE_improvement','NSE_improvement','KGE_improvement', ...
            'S_logp_improvement'});
        improvementTable=[improvementTable;row]; %#ok<AGROW>
    end
end
fprintf('Positive values below indicate improvement over both objective-best members.\n');
disp(improvementTable);
moead.objectiveNames={'1-NSE','1-S_logp'};
save(outputFile,'moead','Qpareto','Qmean','Qmedian', ...
    'comparisonTable','improvementTable','bestNSE','bestSlogp','-append');
writetable(comparisonTable,replace(outputFile,'.mat','_comparison.csv'));
writetable(improvementTable,replace(outputFile,'.mat','_improvement.csv'));


P   = basinData.meteo.P(:);
PET = basinData.meteo.Ep(:);

% Plot the calibration period; these are indices into the scored series.
if isfield(mdl,'local') && mdl.local == 1
    idx = expand_index(basinData.id_train);
else
    idx = expand_index(mdl.id_train);
end
idx = idx(:);
% Meteorology includes spin-up; observations contain only scored days.
% mdl.idx brackets cumulative discharge states: the last endpoint is
% a state boundary, not an additional daily forcing/discharge sample.
forcingIdx = (double(mdl.idx(1)):double(mdl.idx(2))-1)';
assert(numel(forcingIdx)==numel(Qobs) && ...
    all(forcingIdx>=1 & forcingIdx<=numel(P)) && ...
    all(forcingIdx<=numel(PET)), ...
    'The scored forcing indices do not match the discharge series.');
P = P(forcingIdx);
PET = PET(forcingIdx);

% Daily calendar dates for the configured manual training period.
t0 = datetime(prd.dts(3),prd.dts(2),prd.dts(1));
t = t0 + days(idx-idx(1));
x = datenum(t);  % Numeric dates allow a shaded patch.

figure('Name','HYMOD Pareto discharge envelope', ...
    'Color','w','Position',[100 100 1150 750]);
tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

ax(1) = nexttile;
bar(x,P(idx),1,'FaceColor',[0.2 0.45 0.75],'EdgeColor','none');
ylabel('Rainfall (mm/day)');
set(gca,'YDir','reverse');
grid on

ax(2) = nexttile;
plot(x,PET(idx),'Color',[0.85 0.4 0.1],'LineWidth',1);
ylabel('PET (mm/day)');
grid on

ax(3) = nexttile([2 1]);
hold on
hBand = fill([x; flipud(x)], ...
    [Qlo(idx); flipud(Qhi(idx))],[0.7 0.7 0.7], ...
    'EdgeColor','none','FaceAlpha',0.6);
hObs = plot(x,Qobs(idx),'k','LineWidth',0.9);
hMean = plot(x,Qmean(idx),'Color',[0 0.45 0.74],'LineWidth',1);
hMedian = plot(x,Qmedian(idx),'Color',[0.85 0.3 0.1],'LineWidth',1);
ylabel('Discharge (mm/day)');
xlabel('Date');
legend([hBand hObs hMean hMedian], ...
    {'Nondominated solution range','Observed discharge','Pareto mean','Pareto median'}, ...
    'Location','best');
grid on

linkaxes(ax,'x');
for j = 1:numel(ax)
    xlim(ax(j),[x(1) x(end)]);
    datetick(ax(j),'x','mmm yyyy','keeplimits');
end
ax(1).XTickLabel = [];
ax(2).XTickLabel = [];
title(tl,sprintf('HYMOD — gauge %s — %d Pareto solutions', ...
    char(gaugeID),nPareto));


function [L,J] = hymodObjectives(x,mdl,dat,ode,loss)
    request = crr_request(struct('obs',{{'Q'}},'jac',{{'Q'}}));
    [~,out] = crr_model(x,mdl,dat,ode,loss,request);
    L = zeros(2,1); 
    J = zeros(numel(x),2);
    for j = 1:2
        selectedLoss = loss; 
        selectedLoss.fnc = 3;
        if j==2, selectedLoss.fnc=6; end;
        [L(j),J(:,j)] = observation_loss('Q', ...
            out.obs.Q,out.jac.Q,dat,selectedLoss,'train');
        if j==2
            reference=dat.fdc.Q.t.D0_logp;
            assert(isfinite(reference) && reference>0, ...
                'S_logp requires a positive log-FDC reference divergence.');
            L(j)=L(j)/reference; J(:,j)=J(:,j)/reference;
        end
    end
end

function values=hydrographMetrics(q,dat,loss,period,rmsseDenominator)
    id=double(dat.stats.Q.(period).indices(:));
    y=double(dat.obs.Q.value(:)); q=double(q(:));
    % obs.Q.bad can be row or column; enforce column masks.
    bad=logical(dat.obs.Q.bad(:));
    valid=isfinite(y(id)) & isfinite(q(id)) & ~bad(id);
    id=id(valid); y=y(id); prediction=q(id);
    mse=mean((prediction-y).^2);
    selectedLoss=loss; selectedLoss.fnc=3;
    [nseLoss,~,detail]=observation_loss('Q',q,zeros(numel(q),0), ...
        dat,selectedLoss,period);
    cache=dat.fdc.Q.(period(1));
    [~,~,Dlogp]=fdc_metrics_cached(prediction,cache);
    slogp=1-Dlogp/cache.D0_logp;
    bias=100*sum(prediction-y)/sum(y);
    values=[sqrt(mse),sqrt(mse/rmsseDenominator),1-nseLoss, ...
        detail.KGE,slogp,bias];
end
