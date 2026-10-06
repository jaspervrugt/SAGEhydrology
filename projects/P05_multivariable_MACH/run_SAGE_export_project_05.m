% ========================================================================
%   S A G E
%   Sensitivity-Aware Learning for Rapid Training of Process-Based
%   Hydrologic Models across Large Basin Samples
% ========================================================================
%
% Project 5      : Switzerland multivariable
% Title          : Beyond Streamflow: Dual- and Triple-Response Hydrologic
%                  Learning across Switzerland with SAGE
% Summary        : Daily HYD-RESPONSES experiment fitting discharge,
%                  snow-water equivalent, and soil moisture jointly with
%                  HBV.
%
% Model          : hbv
% Region         : Switzerland
% Resolution     : Daily
% Basins         : 184 (training: 150; evaluation: 34)
% Period         : 01 Jan 2007 - 31 Dec 2022
% Temporal split : Manual dates
% Training       : 01 Jan 2007 - 31 Dec 2014
% Evaluation     : 01 Jan 2015 - 31 Dec 2022
% Loss function  : Generalized Least Squares
% Data type      : Discharge, snow water equivalent and soil moisture
%                  content
% Optimizer      : Adam/AdamW
% Exported       : 04-Oct-2026 22:42:09
%
% Generated from the canonical SAGE GUI run-export pipeline.

clear
clc

wallTimerTotal = tic; wallTData = 0;
%% Installation and data paths
projectFolder = fileparts(mfilename('fullpath'));
SAGEhydro = fileparts(fileparts(projectFolder));
root = fileparts(SAGEhydro);
SAGEdir = SAGEhydro;
region = 'HYD_RESPONSES';
regionName = 'Switzerland';
file_univ = 'HYD_184_basins.txt';
file_split = 'split_HYD_184_150_34_1.txt';
dirDroot = fullfile(root,'Data');
dirD = fullfile(dirDroot,'HYD_RESPONSES');
dirM = fullfile(dirD,'daily','timeseries');
dirQ = fullfile(dirD,'daily','timeseries');
dirres = fullfile(SAGEhydro,'results');

%% Model configuration
mdl.model = 6;
mdl.mcode = 4;
mdl.calc = 'parfeval';
mdl.mode = 4;

%% Basin selection and attributes
bas.K = 184;
bas.K_t = 150;
bas.K_e = 34;
bas.sample = 'file';
bas.id_attr = [1 2 4 5 6 7 9 10 16 18 32 33 34 35 36 37 38 39 ...
    45 46 47 49 50 51 52 53];
bas.pr_attr = 0;
bas.dt = 1;
bas.mode = 4;

%% Training and evaluation periods
prd.dt = 1;
prd.method = 'manual';
prd.dts = [1 1 2007];
prd.dte = [31 12 2014];
prd.des = [1 1 2015];
prd.dee = [31 12 2022];
prd.spinup = 365;
prd.eval_was_empty = false;
prd.eval_partial = false;

%% Feedforward neural-network settings
net.h = 32;
net.tf = {'tanh'};
net.seed = 0;

%% Loss-function settings
loss.fnc = 2;
loss.method = 1;
loss.n_win = 7;
loss.M = 2;
loss.observed = {'Q', 'SWE', 'SM'};
loss.weight.Q = 1/3;
loss.weight.SWE = 1/3;
loss.weight.SM = 1/3;
loss.normalization.method = 'auto';

%% Optimization settings
alg.method = 2;
alg.i_max = 250;
alg.lr = 0.01;
alg.beta_1 = 0.9;
alg.beta_2 = 0.999;
alg.vareps = 1e-8;
alg.wdecay = 1e-4;
alg.clipn = 1;
alg.lr_min = 1e-4;
alg.lr_scheme = 'constant';
alg.lr_hold = 0.2;

%% Numerical solver settings
ode.InitStep = 0.01;
ode.MaxStep = 1;
ode.MinStep = 1e-4;
ode.RelTol = 1e-3;
ode.AbsTol = 1e-3;
ode.Order = 2;
ode.maxiter = 10000;
ode.mem = 0;

%% Miscellaneous forcing, output, and plotting settings
misc.meteo.data = 1;
misc.meteo.precip = 1;
misc.meteo.temp = 1;
misc.meteo.pet = 1;
misc.io.prt = 1;
misc.io.file = 1;
misc.attr = 0;
misc.attr_view = 'abs';
misc.plot.Btr = 5;
misc.plot.Bev = 5;
misc.plot.gaugescen.tt = {};
misc.plot.gaugescen.te = {};
misc.plot.gaugescen.et = {};
misc.plot.gaugescen.ee = {};
misc.crr_backend = 'cpp';

parOverride = struct();

%% Initialize SAGEhydrology
% Make the bootstrap helper visible in a fresh session.
sageUtils = fullfile(char(root),'SAGEhydrology','utils');
if ~isfolder(sageUtils)
    error('SAGE utility folder not found: %s',sageUtils);
end
addpath(sageUtils);
region = bootstrap_SAGE(root,region,dirres);
mdl.region = region; % Region metadata for diagnostics and country flag
clear sageUtils
 
%% Prepare data and hydrologic model
% Write SAGEhydrology header 
projectInfo = jsondecode(fileread( ...
    fullfile(projectFolder,'project.json')));
projectHeader = {'Project 5',projectInfo.name; ...
    'Title',projectInfo.paper_title; ...
    'Summary',projectInfo.summary};
write_header(mdl,'sage',projectHeader) 
% Numerical solver 
ode = read_numsettings(ode);
% Catchment attributes, hydroclimatic zones and basin names 
wallTimerData = tic;
[A_reg,ID,gname,zone] = read_attr(region,dirD,bas);
wallTData = wallTData + toc(wallTimerData);
% Select publication basins and reproduce their stored split 
[A,bas,latlon] = sample_basins(A_reg,ID,bas,prd, ...
    gname,zone,dirD,file_univ,file_split);
% Select and prepare the fastest available CRR backend 
[mdl,misc,crrStatus] = crr_prepare_backend(mdl,misc);
mdl.crr_backend = misc.crr_backend; % Keep training and
    % postprocessing synchronized
fprintf('... Preparing CRR backend ... Done\n');
fprintf('      %s.\n',crrStatus);
% Read model parameter information 
[mdl,d] = read_model(mdl,prd);
% Parameter override by user in GUI 
mdl = apply_parameter_override( ... 
    mdl,parOverride,mdl.model,prd.dt);
% Derive network dimensions from data and model.
net.r = size(A,1); net.d = d;
% Determine train/eval split of discharge data 
[split,mdl] = build_split(mdl,prd,bas);
% Read meteorological data 
wallTimerData = tic;
[dat,aux] = read_meteo(region,dirM,bas,split,misc.meteo);
wallTData = wallTData + toc(wallTimerData);
% Read discharge data 
wallTimerData = tic;
dat = read_Q(region,dirQ,mdl,dat,bas,split,aux);
wallTData = wallTData + toc(wallTimerData);
% Verify data and retain active basins 
[eligibility,dat] = check_basins(dat,mdl,bas);
[dat,A,latlon,bas] = filter_basins(dat,A,latlon,bas,eligibility);
% Compute statistics measured data 
[dat,loss] = prep_stats(dat,mdl,split,loss);
% Fire up parallel environment 
setup_parpool(mdl,bas);
% Initialize remaining variables 
[prf,ax,tTheta,At,An,nTheta] = init_args(bas,mdl,alg,misc.attr); 
prf.wallT_data = wallTData;
prf.wallT_init = max(0,toc(wallTimerTotal)-wallTData);
 
%% Train neural network
for i = 1:alg.i_max
    T = tic;
    if i == 1
        % Initialize the FFN and optimizer state
        [phi,net] = ffn_theta('init',net);
        opts = descent('init',alg,phi);
    else
        % Descent step: new parameter values 
        [phi,opts] = descent('dyn',alg,phi,dLdphi,i,opts);
        prf.iter.lr(i) = opts.lr_current;
    end
    % Compute normalized hydrologic parameters 
    [nTheta,cache] = ffn_theta('eval',phi,A);
    % Evaluate basins and return loss and gradients
    [L,G,met,At(:,:,i),An(:,:,i),Qfdc,loss] = ... 
        camels(nTheta,mdl,dat,bas,ode,loss,misc,d,i,dirres);
    % Compute performance metrics
    prf = pmetrics(bas,loss,L,met,i,prf);
    % Update attribution metrics
    tTheta = trace_theta(tTheta,nTheta,bas,i);
    % Compute gradient for the next evaluated state
    if i < alg.i_max
        dLdphi = ffn_theta('back',phi,cache,alg,G);
    end
    prf.iter.cpuT(i) = toc(T);
    % Display performance metrics and cdfs
    ax = print_SAGE(mdl,ax,prf,i,dirres,loss,net);
end
prf.wallT_train = sum(prf.iter.cpuT,'omitnan');
prf.wallT_info = 0; prf.wallT_gui = 0;
%% Postprocessing and export
wallTimerPost = tic;
% Postprocess SAGE results 
Q = postproc_SAGE('sage',nTheta,mdl, ...
    dat,bas,prd,ode,loss,dirres,met);
% Plot SAGE results 
plot_SAGE('sage',mdl,dat,bas,prd,Q,prf.curr, ...
    Qfdc,region,tTheta,latlon,nTheta, ...
    At,An,misc.plot.gaugescen);
% Write SAGE results to binary file 
export_SAGE(mdl,dat,bas,prd,dirres, ... 
    Q,prf,region,nTheta);
prf.wallT_post = toc(wallTimerPost);
prf.wallT_total = toc(wallTimerTotal);
% Write SAGE figures to Powerpoint file 
C = struct('root',root,'SAGEhydro',SAGEhydro, ...
    'SAGEdir',SAGEdir,'region',region,'regionName',regionName, ...
    'dirDroot',dirDroot,'dirD',dirD,'dirM',dirM,'dirQ',dirQ, ...
    'dirres',dirres,'file_univ',file_univ,'file_split',file_split, ...
    'mdl',mdl,'bas',bas,'prd',prd,'ode',ode,'net',net, ...
    'alg',alg,'loss',loss,'misc',misc);
C.runtime = struct('lastIter',i,'timing',struct( ...
    'data',prf.wallT_data, ...
    'initialization',prf.wallT_init, ...
    'training',prf.wallT_train, ...
    'information',prf.wallT_info, ...
    'postprocessing',prf.wallT_post, ...
    'gui',prf.wallT_gui,'total',prf.wallT_total));
save_to_pptx(C);

%% References
% [1] Vrugt, Frame, and Bollman (2026), Reclaiming First Principles: An
%     Analytic Differentiable Framework for Conceptual Hydrologic Models,
%     Water Resources Research. https://arxiv.org/pdf/2602.06429
% [2] Vrugt and Frame (2026), SAGE: Sensitivity-Aware Learning for Rapid 
%     Training of Process-Based Hydrologic Models across Large 
%     Basin Samples, Hydrology and Earth System Sciences Discussions.
%     https://doi.org/10.5194/egusphere-2026-693
% [3] Vrugt, Frame, and Gao (2026), Continental-Scale Hydrologic Model
%     Training at Hourly Resolution with Sensitivity-Aware Gradient
%     Estimation, ARC Geophysical Research.
% [4] Vrugt and Frame (2026), SAGE GUI: Fast and Flexible
%     Differentiable Hydrologic Modeling with Analytic Gradients,
%     Environmental Modelling & Software.
