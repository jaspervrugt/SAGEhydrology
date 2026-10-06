% =========================================================================
%   S A G E
%   Sensitivity-Aware Learning for Rapid Training of Process-Based
%   Hydrologic Models across Large Basin Samples
% =========================================================================
%
% Project 1 : Analytic gradients and Jacobians
% Title     : Reclaiming First Principles: An Analytic Differentiable
%             Framework for Conceptual Hydrologic Models
% Summary   : Verification of analytic Jacobians and loss-function gradients
%             against fixed and adaptive numerical differentiation.
%
% Project 1 verifies both Jacobians and loss-function gradients for the
% four study catchments used in the WRR paper. Its data provenance is:
%
%   1. French Broad River at Asheville, North Carolina (USGS 03451500).
%      The bundled MOPEX record spans 01-Jan-1960 through 31-Dec-1966;
%      the default verification window is 01-Oct-1961 through 30-Sep-1964.
%
%   2. Leaf River near Collins, Mississippi (USGS 02472000). The record
%      is read from the installed CAMELS-US daily product using Maurer
%      meteorological forcing. Install CAMELS-US on the Region tab before
%      running the four-catchment verification. The source record spans
%      01-Oct-1980 through 30-Sep-2008; the default verification window is
%      01-Oct-1981 through 30-Sep-1984.
%
%   3. Wye at Cefn Brwyn (CAMELS-GB 55008) and Severn at Plynlimon flume
%      (CAMELS-GB 54022). Their precipitation and discharge are read from
%      the installed CAMELS-GB hourly product; daily CAMELS-GB PET is
%      divided over the corresponding hours. Install both CAMELS-GB daily
%      and hourly on the Region tab before running the four-catchment
%      verification. The public
%      CAMELS-GB record spans 01-Oct-1990 through 01-Oct-2022; the default
%      verification window is 01-Oct-1993 through 30-Sep-1994.
%
%      The CAMELS-GB discharge series closely match the 1992-1996 workbook
%      supplied by Prof. Jim Kirchner, but the precipitation series differ.
%      The workbook likely uses locally measured experimental-catchment
%      rainfall, whereas CAMELS-GB uses gridded CEH-GEAR precipitation.
%      CAMELS-GB hourly also has no hourly PET field; Project 1 derives its
%      hourly PET reproducibly from the public CAMELS-GB daily product.
%
% Model identifiers
%   1 = HYMOD       2 = HMODEL       3 = SAC-SMA
%   4 = Xinanjiang  5 = GR4J         6 = HBV
%   7 = CFE-NWM
%
% By default, the verification window is the final water year of the
% resolution-specific CAMELS-US training period. Set verification_start
% and verification_end (dd/MM/yyyy) to use another interval.
% =========================================================================

clear
clc

%% Project 1 settings -- edit these values
study_data = 'wrr_paper';             % 'wrr_paper' or 'camels_us'
model_id = 4;                         % integer from 1 through 7
time_resolution = 'daily';            % 'daily' or 'hourly'
basin_ids = ["01022500";"01031500"];  % selected CAMELS-US gauge IDs
verification_start = "";              % dd/MM/yyyy; empty = last water year
verification_end = "";                % dd/MM/yyyy; empty = last water year
parameter_values = [];                % normalized values; [] = midpoint 0.5
h = 1e-4;                             % numerical-differentiation step
fixed_step = [];                      % [] = adaptive WRR ODE settings below
N = 1;                                % WRR paper uses 10 parameter vectors
random_seed = 1;                      % reproducible trials 2,...,N
parameter_margin = 0.15;              % keep generated values off bounds
plot_trial = 1;                       % detailed figures; 0 disables them
run_derivest = true;                  % WRR adaptive Richardson benchmark

%% Initialize SAGEhydrology and resolve CAMELS-US paths
assert(isscalar(model_id) ...
    && isnumeric(model_id) ...
    && isfinite(model_id) ...
    && model_id == fix(model_id) ...
    && ismember(model_id,1:7), ...
    'model_id must be an integer from 1 through 7.');

resolution = validatestring(lower(string(time_resolution)), ...
    {'daily','hourly'});
if strcmp(resolution,'daily')
    resolution_label = 'Daily';
else
    resolution_label = 'Hourly';
end

basin_ids = strip(string(basin_ids(:)));
basin_ids = basin_ids(strlength(basin_ids) > 0);
assert(~isempty(basin_ids), ...
    'Provide at least one CAMELS-US gauge ID in basin_ids.');
assert(numel(unique(basin_ids)) == numel(basin_ids), ...
    'basin_ids contains duplicate gauge IDs.');
assert(isscalar(h) ...
    && isnumeric(h) ...
    && isfinite(h) ...
    && h > 0, ...
    'h must be a positive finite scalar.');

assert(isscalar(N) ...
    && isnumeric(N) ...
    && isfinite(N) ...
    && N == fix(N) ...
    && N >= 1, ...
    'N must be a positive integer.');
assert(isscalar(parameter_margin) ...
    && isfinite(parameter_margin) ...
    && parameter_margin > 0 ...
    && parameter_margin < 0.5, ...
    'parameter_margin must lie strictly between 0 and 0.5.');

script_file = mfilename('fullpath');
SAGEhydro = fileparts(fileparts(script_file));
software_root = fileparts(SAGEhydro);
data_root = fullfile(software_root,'Data');
results_dir = fullfile(SAGEhydro,'results','project_1');
addpath(fullfile(SAGEhydro,'utils'));
addpath(fullfile(SAGEhydro,'utils','pet'));
derivest_dir = fullfile(SAGEhydro,'projects','P01_analytic_gradients', ...
    'third_party','DERIVESTsuite');
if run_derivest
    assert(isfile(fullfile(derivest_dir,'jacobianest.m')) ...
        && isfile(fullfile(derivest_dir,'gradest.m')) ...
        && isfile(fullfile(derivest_dir,'derivest.m')), ...
        ['The vendored DERIVESTsuite subset is incomplete. Restore the ' ...
        'Project 1 third_party/DERIVESTsuite folder.']);
    addpath(derivest_dir);
end

study_data = validatestring(lower(study_data), ...
    {'wrr_paper','camels_us'});
if strcmp(study_data,'wrr_paper')
    [jac_check,benchmark_runs,benchmark_summary] = ...
        verify_wrr_paper(model_id,parameter_values,h,fixed_step,N, ...
        random_seed,plot_trial,run_derivest, ...
        SAGEhydro,software_root,results_dir);
    return
end

region = bootstrap_SAGE(software_root,'CAMELS_US',results_dir);

region_config = region_helpers('config',SAGEhydro,region);
resolution_config = region_helpers( ...
    'resolution',SAGEhydro,region,resolution_label);
paths = region_helpers('resolvepaths', ...
    region_config,resolution_config,data_root);

assert(isfolder(paths.dirM), ...
    'CAMELS-US %s NLDAS forcing folder not found: %s', ...
    resolution,paths.dirM);
assert(isfolder(paths.dirQ), ...
    'CAMELS-US %s streamflow folder not found: %s', ...
    resolution,paths.dirQ);

%% Model, basin, period, and NLDAS configuration
model_names = ["hymod","hmodel","sacsma","Xinanjiang", ...
    "gr4j","hbv","cfe_nwm"];
mdl = struct();
mdl.model = model_id;
mdl.names = model_names;
mdl.mcode = 4;
mdl.calc = 'sequential';
mdl.mode = 1;
mdl.region = region;

bas = struct();
bas.K = numel(basin_ids);
bas.K_t = bas.K;
bas.K_e = 0;
bas.id_t = (1:bas.K).';
bas.id_e = zeros(0,1);
bas.id_gauge = basin_ids;
bas.dt = resolution_config.dt;
bas.resolution = resolution_label;
bas.stream = resolution;
bas.data_actions = struct();

[date_start,date_end] = verification_dates( ...
    resolution_config,verification_start,verification_end);
prd = struct();
prd.dt = resolution_config.dt;
prd.method = 'manual';
prd.dts = date_vector(date_start);
prd.dte = date_vector(date_end);
prd.des = prd.dts;
prd.dee = prd.dte;
prd.spinup = min(30,double(resolution_config.period.spinup_days));
prd.eval_was_empty = true;
prd.eval_partial = false;

misc = struct();
misc.crr_backend = 'cpp';
misc.meteo.data = 3;                  % NLDAS
misc.meteo.precip = 3;                % NLDAS precipitation
misc.meteo.temp = 3;                  % NLDAS temperature
misc.meteo.pet = 1;                   % Penman-Monteith PET

header_rows = project1_header_rows();
header_rows = [header_rows; { ...
    'Resolution',resolution_label; ...
    'Forcing','NLDAS'; ...
    'Basins',char(strjoin(basin_ids,', ')); ...
    'Trials, N',sprintf('%d parameter vector(s)',N); ...
    'FD step, h',sprintf('%.3g (normalized parameter space)',h); ...
    'Period',sprintf('%s to %s (+ %d-day spin-up)', ...
        char(string(date_start,'dd-MMM-yyyy')), ...
        char(string(date_end,'dd-MMM-yyyy')),prd.spinup)}];
if run_derivest
    header_rows(end+1,:) = {'Adaptive', ...
        'D''Errico DERIVESTsuite (bounded Richardson)'};
end
write_header(mdl,'sage',header_rows);

model_message = evalc('[mdl,d] = read_model(mdl,prd);');
print_project1_setup_message(model_message);
[mdl,misc,backend_status] = crr_prepare_backend(mdl,misc);
fprintf('... CRR backend: %s\n',backend_status);
[split,mdl] = build_split(mdl,prd,bas);
[dat,aux] = read_meteo(region,paths.dirM,bas,split,misc.meteo);
dat = read_Q(region,paths.dirQ,mdl,dat,bas,split,aux);

if isempty(parameter_values)
    old_rng = rng;
    cleanup_rng = onCleanup(@() rng(old_rng));
    rng(random_seed,'twister');
    parameter_values = parameter_margin + (1-2*parameter_margin)*rand(d,N);
    parameter_values(:,1) = 0.5;
elseif isvector(parameter_values)
    parameter_values = parameter_values(:);
end
assert(size(parameter_values,1) == d ...
    && size(parameter_values,2) == N, ...
    ['parameter_values must be empty or a %d-by-N matrix containing one ' ...
    'normalized parameter vector per trial.'],d);
assert(all(isfinite(parameter_values)) ...
    && all(parameter_values >= 0 ...
    & parameter_values <= 1,'all'), ...
    'parameter_values must lie in the normalized interval [0,1].');
ode = project1_ode_settings();
if isempty(fixed_step)
    fixed_step = 0.01/prd.dt;
end
assert(isscalar(fixed_step) ...
    && isfinite(fixed_step) ...
    && fixed_step > 0, ...
    'fixed_step must be a positive finite scalar.');
ode.InitStep = fixed_step;
ode.MinStep = fixed_step;
ode.MaxStep = fixed_step;

%% Compare analytic and finite-difference Jacobians for every basin
requested_observables = ["Q";"SWE";"SM"];
observed = supported_observables( ...
    parameter_values(:,1),mdl,dat{1},ode,requested_observables);
assert(~isempty(observed), ...
    'The selected model did not expose any verifiable observables.');
fprintf('  Jacobians  : %s\n',strjoin(observed,', '));
jac_check = struct([]);
for k = 1:bas.K
    for trial = 1:N
        fprintf('\nVerifying basin %s (%d of %d), trial %d of %d\n', ...
            basin_ids(k),k,bas.K,trial,N);
        result = verify_basin(parameter_values(:,trial),mdl,dat{k},ode, ...
            observed,h,basin_ids(k),"CAMELS-US",run_derivest, ...
            trial == plot_trial,trial);
        if k == 1 && trial == 1
            jac_check = result;
        else
            jac_check(k,trial) = result;
        end
    end
end

benchmark_runs = benchmark_run_table(jac_check);
benchmark_summary = benchmark_summary_table(benchmark_runs);
fprintf('\nPer-run benchmark results\n');
disp(benchmark_runs);
fprintf('\nAggregate benchmark (all basins and parameter vectors)\n');
disp(benchmark_summary);
fprintf('\nProject 1 verification completed: %d basin(s) x %d trial(s).\n', ...
    bas.K,N);

function rows = project1_header_rows()
%PROJECT1_HEADER_ROWS Publication identity shared by both verification modes.

    rows = { ...
        'Project','Project 1: Analytic gradients and Jacobians'; ...
        'Title',['Reclaiming First Principles: An Analytic Differentiable ' ...
            'Framework for Conceptual Hydrologic Models']; ...
        'Summary',['Verification of analytic Jacobians and loss-function ' ...
            'gradients against fixed and adaptive numerical differentiation']};
end

function [jac_check,benchmark_runs,benchmark_summary] = verify_wrr_paper( ...
        model_id,parameter_values,h,fixed_step,N,random_seed, ...
        plot_trial,run_derivest,SAGEhydro, ...
        software_root,results_dir)
%VERIFY_WRR_PAPER Reproduce the four-catchment WRR verification design.

    model_names = ["hymod","hmodel","sacsma","Xinanjiang", ...
        "gr4j","hbv","cfe_nwm"];
    paper_data = fullfile(SAGEhydro,'projects', ...
        'P01_analytic_gradients','data');
    data_root = fullfile(software_root,'Data');
    studies = load_wrr_studies(paper_data,data_root);
    bootstrap_SAGE(software_root,'CAMELS_US',results_dir);

    mdl_base = struct('model',model_id,'names',model_names, ...
        'mcode',4,'calc','sequential','mode',1, ...
        'region','WRR_PAPER');
    header_rows = project1_header_rows();
    header_rows = [header_rows; { ...
        'Catchments',['Leaf River near Collins, Mississippi; French Broad ' ...
            'River at Asheville, North Carolina; Wye; Severn']; ...
        'Data',['Bundled French Broad MOPEX; CAMELS-US Maurer; ' ...
            'CAMELS-GB hourly']; ...
        'Trials, N',sprintf('%d parameter vector(s)',N); ...
        'FD step, h',sprintf('%.3g (normalized parameter space)',h); ...
        'ODE solver',['Heun order 2; initial/minimum/maximum step = ' ...
            '1e-3/1e-4/1 day; RelTol = AbsTol = 1e-5']}];
    if run_derivest
        header_rows(end+1,:) = {'Adaptive', ...
            'D''Errico DERIVESTsuite (bounded Richardson)'};
    end
    write_header(mdl_base,'sage',header_rows);
    misc = struct('crr_backend','cpp');
    ode_base = project1_ode_settings();
    jac_check = struct([]); parameter_matrix = [];

    for k = 1:numel(studies)
        study = studies(k);
        prd = struct('dt',study.dt,'spinup',0); %#ok<NASGU>
        model_message = evalc( ...
            '[mdl,d] = read_model(mdl_base,prd,k == 1);');
        print_project1_setup_message(model_message);
        [mdl,misc,backend_status] = crr_prepare_backend(mdl,misc);
        if k == 1
            fprintf('%-11s: %s\n','CRR backend',backend_status);
        end
        n = numel(study.P);
        mdl.tout = n;
        mdl.idx = [1 n+1];
        mdl.id_train = [1 n];
        mdl.id_eval = zeros(1,0);
        mdl.local = false;

        dat = struct();
        dat.meteo = struct('P',study.P(:),'Ep',study.Ep(:), ...
            'T',10*ones(n,1));
        dat.obs.Q.value = study.Q(:);
        dat.obs.Q.bad = ~isfinite(study.Q(:));

        if isempty(parameter_matrix)
            parameter_matrix = publication_parameters(parameter_values, ...
                d,N,random_seed,run_derivest);
        else
            assert(size(parameter_matrix,1) == d, ...
                'All WRR study runs must use the same selected model.');
        end

        ode = ode_base;
        if ~isempty(fixed_step)
            assert(isscalar(fixed_step) ...
                && isfinite(fixed_step) ...
                && fixed_step > 0, ...
                'fixed_step must be empty or a positive finite scalar.');
            ode.InitStep = fixed_step;
            ode.MinStep = fixed_step;
            ode.MaxStep = fixed_step;
        end

        fprintf('\n%s: %s, %s to %s, %d observations\n', ...
            study.name,study.resolution, ...
            char(string(study.start_date,'dd-MMM-yyyy')), ...
            char(string(study.end_date,'dd-MMM-yyyy')),n);
        for trial = 1:N
            fprintf('  Trial %d of %d\n',trial,N);
            result = verify_basin(parameter_matrix(:,trial),mdl,dat,ode, ...
                "Q",h,study.name,study.source,run_derivest, ...
                trial == plot_trial,trial);
            if k == 1 && trial == 1
                jac_check = result;
            else
                jac_check(k,trial) = result;
            end
        end
    end

    benchmark_runs = benchmark_run_table(jac_check);
    benchmark_summary = benchmark_summary_table(benchmark_runs);
    fprintf('\nPer-run benchmark results\n');
    disp(benchmark_runs);
    fprintf('\nAggregate WRR benchmark\n');
    disp(benchmark_summary);
    fprintf('\nProject 1 WRR verification completed: 4 catchments x %d trials.\n',N);
end

function values = publication_parameters(values,d,N,seed,~)
%PUBLICATION_PARAMETERS Resolve reproducible normalized parameter vectors.

    if isempty(values)
        old_rng = rng;
        cleanup_rng = onCleanup(@() rng(old_rng));
        rng(seed,'twister');
        % The WRR experiments sampled unconstrained values on [-3,3] and
        % mapped them through the logistic function to the unit hypercube.
        values = 1./(1+exp(-(-3+6*rand(d,N))));
    elseif isvector(values)
        values = values(:);
    end
    assert(size(values,1) == d ...
        && size(values,2) == N, ...
        'parameter_values must be empty or a %d-by-N matrix.',d);
    assert(all(isfinite(values),'all') ...
        && all(values >= 0 ...
        & values <= 1,'all'), ...
        'parameter_values must lie in [0,1].');
end

function ode = project1_ode_settings()
%PROJECT1_ODE_SETTINGS Numerical settings used in the WRR experiments.
%
% Keep these values explicit and local to Project 1. The general SAGE
% defaults are intentionally looser and may evolve independently.

    settings = struct();
    settings.InitStep = 1e-3;
    settings.MaxStep = 1.0;
    settings.MinStep = 1e-4;
    settings.RelTol = 1e-5;
    settings.AbsTol = 1e-5;
    settings.Order = 2;
    settings.maxiter = 1e4;
    settings.mem = 0; %#ok<STRNU> Used by captured command.
    ode = struct(); % Assigned by the captured command below.
    setup_message = evalc('ode = read_numsettings(settings);');
    print_project1_setup_message(setup_message);
end

function print_project1_setup_message(message)
%PRINT_PROJECT1_SETUP_MESSAGE Indent setup output for this script only.

    lines = splitlines(string(message));
    lines = strip(lines);
    lines = lines(strlength(lines) > 0);
    for k = 1:numel(lines)
        fprintf('             %s\n',lines(k));
    end
end

function studies = load_wrr_studies(folder,data_root)
%LOAD_WRR_STUDIES Load the public products used by Project 1.

    french_file = fullfile(folder,'03451500.dly');
    assert(isfile(french_file), ...
        ['The bundled Project 1 French Broad MOPEX file is missing: %s. ' ...
        'Repair or reinstall SAGEhydrology.'],french_file);
    french = load(fullfile(folder,'03451500.dly'));
    french_dates = datetime(french(:,1),french(:,2),french(:,3));
    [leaf_dates,leaf_P,leaf_Ep,leaf_Q] = ...
        load_leaf_camels_us(data_root);
    [wye_dates,wye_P,wye_Ep,wye_Q] = ...
        load_camels_gb_hourly(data_root,'55008');
    [severn_dates,severn_P,severn_Ep,severn_Q] = ...
        load_camels_gb_hourly(data_root,'54022');

    studies = repmat(struct('name',"",'source',"",'resolution',"",'dt',1, ...
        'start_date',NaT,'end_date',NaT,'P',[],'Ep',[],'Q',[]),4,1);
    studies(1) = make_wrr_study( ...
        "Leaf River near Collins, Mississippi","CAMELS-US", ...
        "daily",1,leaf_dates, ...
        leaf_P,leaf_Ep,leaf_Q,datetime(1981,10,1),datetime(1984,9,30));
    studies(2) = make_wrr_study( ...
        "French Broad River at Asheville, North Carolina", ...
        "MOPEX","daily",1,french_dates, ...
        french(:,4),french(:,5),french(:,6),datetime(1961,10,1),datetime(1964,9,30));
    studies(3) = make_wrr_study("Wye","CAMELS-GB v2", ...
        "hourly",24,wye_dates, ...
        wye_P,wye_Ep,wye_Q, ...
        datetime(1993,10,1),datetime(1994,9,30));
    studies(4) = make_wrr_study("Severn","CAMELS-GB v2", ...
        "hourly",24,severn_dates, ...
        severn_P,severn_Ep,severn_Q, ...
        datetime(1993,10,1),datetime(1994,9,30));
end

function [dates,P,Ep,Q] = load_leaf_camels_us(data_root)
%LOAD_LEAF_CAMELS_US Read Leaf River from CAMELS-US daily/Maurer.

    forcing_root = fullfile(data_root,'CAMELS_US','daily','v1p2', ...
        'forcing','maurer');
    discharge_root = fullfile(data_root,'CAMELS_US','daily','v1p2', ...
        'streamflow');
    forcing = dir(fullfile(forcing_root,'**', ...
        '02472000_lump_maurer_forcing_leap.txt'));
    discharge = dir(fullfile(discharge_root,'**', ...
        '02472000_streamflow_qc.txt'));
    assert(~isempty(forcing) ...
        && ~isempty(discharge), ...
        ['Project 1 requires CAMELS-US daily data with Maurer forcing ' ...
        'for Leaf River (02472000). Install CAMELS-US Daily on the ' ...
        'Region tab, then run this script again.']);

    meteo_file = fullfile(forcing(1).folder,forcing(1).name);
    M = readmatrix(meteo_file,'FileType','text','NumHeaderLines',4);
    dates_m = datetime(M(:,1),M(:,2),M(:,3));
    P_m = M(:,6);
    radiation_24h = M(:,7).*M(:,5)/86400;
    Ep_m = et0_fao56_daily(M(:,9),M(:,10),radiation_24h,M(:,11), ...
        day(dates_m,'dayofyear'),32.13,137,2,1);

    discharge_file = fullfile(discharge(1).folder,discharge(1).name);
    fid = fopen(discharge_file,'r');
    cleanup = onCleanup(@() fclose(fid));
    C = textscan(fid,'%s %f %f %f %f %s');
    dates_q = datetime(C{2},C{3},C{4});
    % The original Project 1 conversion used 1,927.13 km^2.
    Q_q = C{5}*0.028316846592*86400/1927130000*1000;

    [dates,im,iq] = intersect(dates_m,dates_q,'stable');
    P = P_m(im);
    Ep = Ep_m(im);
    Q = Q_q(iq);
end

function [dates,P,Ep,Q] = load_camels_gb_hourly(data_root,gauge)
%LOAD_CAMELS_GB_HOURLY Read hourly P/Q and derive hourly PET from daily PET.

    hourly_root = fullfile(data_root,'CAMELS_GB','hourly','timeseries');
    daily_root = fullfile(data_root,'CAMELS_GB','daily','timeseries');
    hourly = dir(fullfile(hourly_root, ...
        ['*hydromet_hourly_timeseries_' gauge '_*.csv']));
    daily = dir(fullfile(daily_root, ...
        ['*hydromet_timeseries_' gauge '_*.csv']));
    assert(~isempty(hourly) ...
        && ~isempty(daily), ...
        ['Project 1 requires both CAMELS-GB Hourly and Daily data for ' ...
        'Wye (55008) and Severn (54022). Install both resolutions on ' ...
        'the Region tab, then run this script again.']);

    H = readtable(fullfile(hourly(1).folder,hourly(1).name), ...
        'VariableNamingRule','preserve');
    D = readtable(fullfile(daily(1).folder,daily(1).name), ...
        'VariableNamingRule','preserve');
    dates = H.date;
    if ~isdatetime(dates)
        dates = datetime(dates,'InputFormat','yyyy-MM-dd HH:mm:ss');
    end
    daily_dates = D.date;
    if ~isdatetime(daily_dates)
        daily_dates = datetime(daily_dates,'InputFormat','yyyy-MM-dd');
    end
    [found,id] = ismember(dateshift(dates,'start','day'),daily_dates);
    P = H.precipitation_cehgear;
    Q = H.discharge_spec;
    Ep = nan(height(H),1);
    Ep(found) = D.pet(id(found))/24;
end

function study = make_wrr_study(name,source,resolution,dt,dates,P,Ep,Q,d0,d1)
%MAKE_WRR_STUDY Select the publication window and validate its forcing.

    use = dates >= d0 ...
        & dates < d1 + days(1);
    study = struct('name',name,'source',source, ...
        'resolution',resolution,'dt',dt, ...
        'start_date',d0,'end_date',d1,'P',double(P(use)), ...
        'Ep',double(Ep(use)),'Q',double(Q(use)));
    assert(~isempty(study.P) ...
        && numel(study.P) == numel(study.Ep) ...
        && numel(study.P) == numel(study.Q), ...
        'Incomplete WRR data for %s.',name);
    assert(all(isfinite(study.P)) ...
        && all(isfinite(study.Ep)), ...
        'Nonfinite meteorological forcing in WRR data for %s.',name);
end

function observed = supported_observables(x,mdl,dat,ode,candidates)
%SUPPORTED_OBSERVABLES Retain model outputs exposed by the selected model.

    observed = strings(0,1);
    skipped = strings(0,1);
    for v = 1:numel(candidates)
        name = candidates(v);
        try
            request = crr_request(struct('obs',name,'jac',name));
            crr_model(x,mdl,dat,ode,struct(),request);
            observed(end+1,1) = name; %#ok<AGROW>
        catch ME
            if contains(ME.message,'unavailable','IgnoreCase',true) ...
                    || contains(ME.message,'unsupported','IgnoreCase',true)
                skipped(end+1,1) = name; %#ok<AGROW>
            else
                rethrow(ME)
            end
        end
    end
    if ~isempty(skipped)
        fprintf('  Not exposed : %s (skipped for model %s)\n', ...
            strjoin(skipped,', '),string(mdl.names(mdl.model)));
    end
end

function result = verify_basin(x,mdl,dat,ode,observed,h,gauge,source, ...
        run_derivest,do_plot,trial)
%VERIFY_BASIN Compare analytic and numerical Jacobians for one basin.

    request = crr_request(struct('obs',observed,'jac',observed));
    request_fd = crr_request(struct('obs',observed));

    % Warm up both execution paths before measuring wall-clock time. This
    % avoids charging one method for first-call MEX/JIT initialization.
    crr_model(x,mdl,dat,ode,struct(),request);
    crr_model(x,mdl,dat,ode,struct(),request_fd);

    analytic_timer = tic;
    [~,ana] = crr_model(x,mdl,dat,ode,struct(),request);
    analytic_seconds = toc(analytic_timer);

    value = struct();
    J = struct();
    Jfd = struct();
    for v = 1:numel(observed)
        name = char(observed(v));
        value.(name) = ana.obs.(name)(:);
        J.(name) = ana.jac.(name);
        Jfd.(name) = nan(size(J.(name)));
        assert(any(abs(value.(name)) > 10*eps) ...
            || any(abs(J.(name)) > 10*eps,'all'), ...
            ['Degenerate %s verification for gauge %s: both the ' ...
            'simulated trajectory and its analytic Jacobian are zero. ' ...
            'Choose another parameter vector or reduce fixed_step.'], ...
            name,gauge);
    end

    d = numel(x);
    scheme = strings(d,1);
    step = nan(d,1);
    [lower,upper] = parameter_bounds(mdl,d);
    numerical_evaluations = 0;
    numerical_timer = tic;

    for j = 1:d
        h_j = h*max(1,abs(x(j)));
        h_j = max(h_j,10*eps(max(1,abs(x(j)))));
        room_lo = x(j)-lower(j);
        room_hi = upper(j)-x(j);

        if room_lo >= h_j && room_hi >= h_j
            xp = x; xm = x;
            xp(j) = x(j)+h_j; xm(j) = x(j)-h_j;
            [~,out_p] = crr_model(xp,mdl,dat,ode,struct(),request_fd);
            [~,out_m] = crr_model(xm,mdl,dat,ode,struct(),request_fd);
            numerical_evaluations = numerical_evaluations+2;
            for v = 1:numel(observed)
                name = char(observed(v));
                Jfd.(name)(:,j) = (out_p.obs.(name)(:) ...
                    - out_m.obs.(name)(:))/(2*h_j);
            end
            scheme(j) = "central";
        elseif room_hi > 0
            h_j = min(h_j,room_hi);
            xp = x; xp(j) = x(j)+h_j;
            [~,out_p] = crr_model(xp,mdl,dat,ode,struct(),request_fd);
            numerical_evaluations = numerical_evaluations+1;
            for v = 1:numel(observed)
                name = char(observed(v));
                Jfd.(name)(:,j) = (out_p.obs.(name)(:) ...
                    - value.(name))/h_j;
            end
            scheme(j) = "forward";
        elseif room_lo > 0
            h_j = min(h_j,room_lo);
            xm = x; xm(j) = x(j)-h_j;
            [~,out_m] = crr_model(xm,mdl,dat,ode,struct(),request_fd);
            numerical_evaluations = numerical_evaluations+1;
            for v = 1:numel(observed)
                name = char(observed(v));
                Jfd.(name)(:,j) = (value.(name) ...
                    - out_m.obs.(name)(:))/h_j;
            end
            scheme(j) = "backward";
        else
            error('Parameter %d has no admissible perturbation.',j);
        end
        step(j) = h_j;
    end
    numerical_seconds = toc(numerical_timer);

    % D'Errico's estimators use trial steps as large as 100 times the
    % current coordinate.  The legacy WRR model accepted those trials even
    % when they left the normalized parameter cube; current SAGE models
    % intentionally reject them.  Evaluate Richardson extrapolation in a
    % smooth unconstrained coordinate and transform the resulting local
    % derivatives back to the selected parameter space.  Thus the reported
    % Jacobian and gradient are still derivatives with respect to x.
    Jr = struct(); Er = struct();
    richardson_seconds = NaN; gradest_seconds = NaN;
    richardson_evaluations = 0;
    gradient = struct('analytic',[],'fixed',[],'adaptive',[], ...
        'adaptive_error_estimate',[],'adaptive_final_step',[], ...
        'relative_error',NaN);
    [y0,Ja_all,part] = concatenate_observables(value,J,observed);
    gradient.analytic = Ja_all.'*y0/max(numel(y0),1);
    [~,Jfd_all] = concatenate_observables(value,Jfd,observed);
    gradient.fixed = Jfd_all.'*y0/max(numel(y0),1);
    if run_derivest
        [z,dxdz] = adaptive_coordinates(x,mdl);
        richardson_timer = tic;
        [Jz_all,Ez_all] = jacobianest(@(coordinate) observations_at( ...
            parameters_from_adaptive(coordinate,mdl),mdl,dat,ode, ...
            request_fd,observed),z);
        richardson_seconds = toc(richardson_timer);
        richardson_evaluations = 1+52*d;
        Jr_all = Jz_all./reshape(dxdz,1,[]);
        Er_all = Ez_all./reshape(abs(dxdz),1,[]);
        for v = 1:numel(observed)
            name = char(observed(v));
            rows = part{v};
            Jr.(name) = Jr_all(rows,:);
            Er.(name) = Er_all(rows,:);
        end

        gradest_timer = tic;
        [gradient.adaptive,gradient.adaptive_error_estimate, ...
            gradient.adaptive_final_step] = gradest(@(coordinate) ...
            objective_at(parameters_from_adaptive(coordinate,mdl), ...
            mdl,dat,ode,request_fd,observed),z);
        gradest_seconds = toc(gradest_timer);
        gradient.adaptive = gradient.adaptive(:)./dxdz;
        gradient.adaptive_error_estimate = ...
            abs(gradient.adaptive_error_estimate(:)./dxdz);
        gradient.adaptive_final_step = ...
            abs(gradient.adaptive_final_step(:).*dxdz);
        gradient.relative_error = vector_relative_error( ...
            gradient.analytic,gradient.adaptive);
    end

    names = parameter_names(mdl,d);
    runtime = struct('analytic_seconds',analytic_seconds, ...
        'numerical_seconds',numerical_seconds, ...
        'richardson_seconds',richardson_seconds, ...
        'gradest_seconds',gradest_seconds, ...
        'analytic_evaluations',1, ...
        'numerical_evaluations',numerical_evaluations, ...
        'richardson_evaluations',richardson_evaluations, ...
        'fixed_over_analytic',numerical_seconds/max(analytic_seconds,eps), ...
        'richardson_over_analytic', ...
            richardson_seconds/max(analytic_seconds,eps));
    result = struct('gauge',gauge,'trial',trial,'x',x,'h',h,'step',step, ...
        'scheme',scheme,'obs',value,'analytic',J,'numerical',Jfd, ...
        'richardson',Jr,'richardson_error_estimate',Er, ...
        'gradient',gradient,'error',struct(),'runtime',runtime);
    for v = 1:numel(observed)
        name = char(observed(v));
        [fixed_relative,fixed_maximum] = ...
            jacobian_error(J.(name),Jfd.(name));
        adaptive_relative = nan(d,1); adaptive_maximum = nan(d,1);
        adaptive_estimated_maximum = nan(d,1);
        if run_derivest
            [adaptive_relative,adaptive_maximum] = ...
                jacobian_error(J.(name),Jr.(name));
            adaptive_estimated_maximum = max(Er.(name),[],1).';
        end
        result.error.(name) = table(names,scheme,step, ...
            fixed_relative,fixed_maximum,adaptive_relative, ...
            adaptive_maximum,adaptive_estimated_maximum, ...
            'VariableNames',{'parameter','scheme','step', ...
            'fixed_relative_error','fixed_maximum_absolute_error', ...
            'adaptive_relative_error','adaptive_maximum_absolute_error', ...
            'adaptive_estimated_maximum_error'});
        fprintf('\n%s Jacobian, gauge %s, trial %d\n',name,gauge,trial);
        disp(result.error.(name));
        if do_plot
            plot_jacobian(J.(name),Jfd.(name),Jr.(name),names, ...
                fixed_relative,adaptive_relative,name,gauge,source, ...
                run_derivest);
        end
    end
    fprintf(['\nRuntime, gauge %s\n' ...
        '  Analytic Jacobian : %.4f s (%d model evaluation)\n' ...
        '  Fixed-step central : %.4f s (%d model evaluations)\n'], ...
        gauge,runtime.analytic_seconds,runtime.analytic_evaluations, ...
        runtime.numerical_seconds,runtime.numerical_evaluations);
    if run_derivest
        fprintf(['  Adaptive Richardson: %.4f s (%d model evaluations)\n' ...
            '  Adaptive gradest   : %.4f s\n' ...
            '  Fixed/analytic     : %.2fx\n' ...
            '  Adaptive/analytic  : %.2fx\n' ...
            '  gradest rel. error : %.3e\n'], ...
            runtime.richardson_seconds,runtime.richardson_evaluations, ...
            runtime.gradest_seconds,runtime.fixed_over_analytic, ...
            runtime.richardson_over_analytic,gradient.relative_error);
    else
        fprintf('  Fixed/analytic     : %.2fx\n',runtime.fixed_over_analytic);
    end
end

function [z,dxdz] = adaptive_coordinates(x,mdl)
%ADAPTIVE_COORDINATES Map bounded parameters to unconstrained coordinates.

    x = x(:);
    [lower,upper] = parameter_bounds(mdl,numel(x));
    bounded = isfinite(lower) & isfinite(upper);
    z = x;
    dxdz = ones(size(x));
    if any(bounded)
        width = upper(bounded)-lower(bounded);
        unit = (x(bounded)-lower(bounded))./width;
        margin = sqrt(eps(class(x)));
        unit = min(max(unit,margin),1-margin);
        z(bounded) = log(unit./(1-unit));
        dxdz(bounded) = width.*unit.*(1-unit);
    end
end

function x = parameters_from_adaptive(z,mdl)
%PARAMETERS_FROM_ADAPTIVE Map unconstrained trials to valid parameters.

    z = z(:);
    [lower,upper] = parameter_bounds(mdl,numel(z));
    bounded = isfinite(lower) & isfinite(upper);
    x = z;
    if any(bounded)
        positive = z(bounded) >= 0;
        unit = zeros(sum(bounded),1);
        unit(positive) = 1./(1+exp(-z(bounded & z >= 0)));
        negative_values = z(bounded & z < 0);
        exponential = exp(negative_values);
        unit(~positive) = exponential./(1+exponential);
        x(bounded) = lower(bounded) + ...
            (upper(bounded)-lower(bounded)).*unit;
    end
end

function [date_start,date_end] = verification_dates(X,start_text,end_text)
%VERIFICATION_DATES Resolve an explicit or default one-water-year window.

    has_start = strlength(strip(string(start_text))) > 0;
    has_end = strlength(strip(string(end_text))) > 0;
    assert(has_start == has_end, ...
        'Set both verification_start and verification_end, or neither.');
    if has_start
        date_start = datetime(start_text,'InputFormat','dd/MM/yyyy');
        date_end = datetime(end_text,'InputFormat','dd/MM/yyyy');
    else
        date_end = datetime(X.period.manual.train_end, ...
            'InputFormat','dd/MM/yyyy');
        date_start = date_end-calyears(1)+days(1);
    end
    assert(date_end >= date_start, ...
        'verification_end must be on or after verification_start.');
end

function ymd = date_vector(value)
%DATE_VECTOR Convert datetime to the SAGE [day month year] convention.

    ymd = [day(value) month(value) year(value)];
end

function [lower,upper] = parameter_bounds(mdl,d)
%PARAMETER_BOUNDS Return bounds in the selected parameter space.

    switch mdl.pspace
        case 0
            lower = mdl.th_min(:); upper = mdl.th_max(:);
        case 1
            lower = zeros(d,1); upper = ones(d,1);
        case 2
            lower = -inf(d,1); upper = inf(d,1);
        otherwise
            error('Unknown mdl.pspace value: %g.',mdl.pspace);
    end
end

function names = parameter_names(mdl,d)
%PARAMETER_NAMES Return printable parameter names.

    if isfield(mdl,'par_names') && numel(mdl.par_names) == d
        names = string(mdl.par_names(:));
    else
        names = "parameter " + string((1:d)');
    end
end

function [relative,maximum] = jacobian_error(Ja,Jn)
%JACOBIAN_ERROR Return columnwise relative and maximum absolute errors.

    d = size(Ja,2);
    relative = nan(d,1); maximum = nan(d,1);
    for j = 1:d
        good = isfinite(Ja(:,j)) & isfinite(Jn(:,j));
        if ~any(good), continue, end
        delta = Ja(good,j)-Jn(good,j);
        scale = max([norm(Ja(good,j)),norm(Jn(good,j)),eps]);
        relative(j) = norm(delta)/scale;
        maximum(j) = max(abs(delta));
        if norm(Ja(good,j)) <= 10*eps ...
                && norm(Jn(good,j)) <= 10*eps
            relative(j) = NaN;
            maximum(j) = NaN;
        end
    end
end

function plot_jacobian(Ja,Jn,Jr,names,fixed_error,adaptive_error, ...
        observable,gauge,source,has_adaptive)
%PLOT_JACOBIAN Plot analytic, fixed-step, and adaptive Jacobian columns.

    d = size(Ja,2);
    ncol = min(4,ceil(sqrt(d))); nrow = ceil(d/ncol);
    figure('Color','w','Name',sprintf('%s Jacobian: %s %s', ...
        observable,source,gauge));
    tl = tiledlayout(nrow,ncol,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('%s Jacobian: %s %s',observable,source,gauge), ...
        'Interpreter','none');
    xlabel(tl,'Model output index'); ylabel(tl,['d' observable '/dx']);
    for j = 1:d
        ax = nexttile(tl);
        plot(ax,Ja(:,j),'b-','LineWidth',1.1); hold(ax,'on');
        plot(ax,Jn(:,j),'r--','LineWidth',1.0);
        if has_adaptive
            plot(ax,Jr(:,j),'Color',[0.00 0.55 0.25], ...
                'LineStyle',':','LineWidth',1.2);
        end
        grid(ax,'on');
        title(ax,sprintf('%s; fixed %.2e; adaptive %.2e', ...
            names(j),fixed_error(j),adaptive_error(j)), ...
            'Interpreter','none');
        if j == 1
            if has_adaptive
                labels = {'analytic','fixed central','adaptive Richardson'};
            else
                labels = {'analytic','fixed central'};
            end
            legend(ax,labels,'Location','best');
        end
    end
end

function y = observations_at(theta,mdl,dat,ode,request,observed)
%OBSERVATIONS_AT Evaluate concatenated outputs in the selected parameter space.

    [~,out] = crr_model(theta(:),mdl,dat,ode,struct(),request);
    y = zeros(0,1);
    for v = 1:numel(observed)
        y = [y; out.obs.(char(observed(v)))(:)]; %#ok<AGROW>
    end
end

function value = objective_at(theta,mdl,dat,ode,request,observed)
%OBJECTIVE_AT Scalar least-squares functional used by GRADEST.

    y = observations_at(theta,mdl,dat,ode,request,observed);
    value = 0.5*mean(y.^2);
end

function [y,Jall,parts] = concatenate_observables(value,J,observed)
%CONCATENATE_OBSERVABLES Stack output vectors and matching Jacobian rows.

    y = zeros(0,1); Jall = zeros(0,size(J.(char(observed(1))),2));
    parts = cell(numel(observed),1);
    cursor = 0;
    for v = 1:numel(observed)
        name = char(observed(v));
        n = numel(value.(name));
        parts{v} = cursor+(1:n);
        y = [y; value.(name)(:)]; %#ok<AGROW>
        Jall = [Jall; J.(name)]; %#ok<AGROW>
        cursor = cursor+n;
    end
end

function value = vector_relative_error(reference,estimate)
%VECTOR_RELATIVE_ERROR Symmetric normwise relative error.

    good = isfinite(reference) & isfinite(estimate);
    if ~any(good)
        value = NaN;
        return
    end
    value = norm(reference(good)-estimate(good))/max( ...
        [norm(reference(good)),norm(estimate(good)),eps]);
end

function runs = benchmark_run_table(results)
%BENCHMARK_RUN_TABLE Flatten basin/trial results into one table.

    count = numel(results);
    gauge = strings(count,1); trial = zeros(count,1);
    analytic_seconds = nan(count,1); fixed_seconds = nan(count,1);
    adaptive_seconds = nan(count,1); gradest_seconds = nan(count,1);
    fixed_relative_error = nan(count,1);
    adaptive_relative_error = nan(count,1);
    gradest_relative_error = nan(count,1);
    for q = 1:count
        item = results(q);
        gauge(q) = string(item.gauge); trial(q) = item.trial;
        analytic_seconds(q) = item.runtime.analytic_seconds;
        fixed_seconds(q) = item.runtime.numerical_seconds;
        adaptive_seconds(q) = item.runtime.richardson_seconds;
        gradest_seconds(q) = item.runtime.gradest_seconds;
        fixed_relative_error(q) = pooled_jacobian_error(item,'fixed');
        adaptive_relative_error(q) = pooled_jacobian_error(item,'adaptive');
        gradest_relative_error(q) = item.gradient.relative_error;
    end
    runs = table(gauge,trial,analytic_seconds,fixed_seconds, ...
        adaptive_seconds,gradest_seconds,fixed_relative_error, ...
        adaptive_relative_error,gradest_relative_error);
end

function value = pooled_jacobian_error(item,method)
%POOLED_JACOBIAN_ERROR Median finite columnwise error over all outputs.

    fields = fieldnames(item.error); values = zeros(0,1);
    variable = [method '_relative_error'];
    for q = 1:numel(fields)
        values = [values; item.error.(fields{q}).(variable)]; %#ok<AGROW>
    end
    values = values(isfinite(values));
    if isempty(values), value = NaN; else, value = median(values); end
end

function summary = benchmark_summary_table(runs)
%BENCHMARK_SUMMARY_TABLE Aggregate timing and accuracy across all trials.

    method = ["analytic";"fixed central";"adaptive Richardson";"gradest"];
    timing = {runs.analytic_seconds,runs.fixed_seconds, ...
        runs.adaptive_seconds,runs.gradest_seconds};
    errors = {nan(height(runs),1),runs.fixed_relative_error, ...
        runs.adaptive_relative_error,runs.gradest_relative_error};
    median_seconds = nan(4,1); q25_seconds = nan(4,1);
    q75_seconds = nan(4,1); minimum_seconds = nan(4,1);
    maximum_seconds = nan(4,1); median_relative_error = nan(4,1);
    for q = 1:4
        stats = five_stats(timing{q});
        median_seconds(q) = stats(3); q25_seconds(q) = stats(2);
        q75_seconds(q) = stats(4); minimum_seconds(q) = stats(1);
        maximum_seconds(q) = stats(5);
        error_stats = five_stats(errors{q});
        median_relative_error(q) = error_stats(3);
    end
    summary = table(method,median_seconds,q25_seconds,q75_seconds, ...
        minimum_seconds,maximum_seconds,median_relative_error);
end

function stats = five_stats(values)
%FIVE_STATS Return min, quartile, median, quartile, and max without toolboxes.

    values = sort(values(isfinite(values)));
    values = values(:).';
    if isempty(values)
        stats = nan(1,5); return
    end
    p = [0 0.25 0.5 0.75 1]; n = numel(values);
    position = 1+(n-1)*p;
    lo = floor(position); hi = ceil(position); weight = position-lo;
    stats = values(lo).*(1-weight) + values(hi).*weight;
end
