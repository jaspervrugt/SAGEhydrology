function varargout = descent(stage,alg,x,dLdphi,i,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DESCENT Initialize or update FFN optimization state.
%
%  Implements normalized gradient descent or Adam/AdamW with an optional
%  learning-rate schedule.
%
% SYNOPSIS:
%   opts = descent('init',alg,phi)
%   [phi,opts] = descent('dyn',alg,phi,dLdphi,i,opts)
%
% INPUT ARGUMENTS:
%   stage           'init' or 'dyn'
%   alg             optimizer settings
%    .method         1 normalized descent; 2 Adam/AdamW
%    .i_max          number of evaluated SAGE iterations
%    .lr             initial learning rate
%    .lr_min         optional minimum/final learning rate
%    .lr_scheme      constant, cosine, or hold_cosine
%    .lr_hold        optional initial hold fraction
%    .beta_1         Adam first-moment decay
%    .beta_2         Adam second-moment decay
%    .vareps         Adam stabilizer
%    .wdecay         optional AdamW weight decay
%   phi             network weights, biases, layers, and transfers
%   dLdphi          network-parameter loss gradient for 'dyn'
%   i               outer SAGE iteration number
%   opts            optimizer state from 'init' or prior update
%
% OUTPUT ARGUMENTS:
%   phi             updated network parameters for 'dyn'
%   opts            optimizer state and learning-rate diagnostics
%    .t              completed optimizer updates
%    .lr_current     learning rate of latest proposal
%    .n_phi          number of network weights and biases
%
% NOTES:
%   Global gradient clipping is applied upstream by FFN_THETA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3
        x = [];
    end
    if nargin < 4
        dLdphi = [];
    end
    if nargin < 5
        i = [];
    end
    if nargin < 6
        opts = [];
    end

    stage = lower(strtrim(char(string(stage))));

    switch stage
        case 'init'
            phi = local_validate_phi_struct(x);
            opts = local_initialize_opts(alg,phi);
            varargout = {opts};

        case 'dyn'
            phi = local_validate_phi_struct(x);

            if isempty(dLdphi)
                error(['      Error: descent(''dyn''): ' ...
                    'dLdphi is empty. Use ' ...
                    'ffn_theta(''init'',net) to create ' ...
                    'the initial phi before optimizer setup.']);
            end
            if isempty(opts) ...
                    || ~isstruct(opts)
                error(['      Error: descent(''dyn''): ' ...
                    'opts must be the state returned ' ...
                    'by descent(''init'',alg,phi).']);
            end
            if ~isempty(i) ...
                    && (~isnumeric(i) ...
                    || ~isscalar(i) ...
                    || ~isfinite(i) ...
                    || i < 1)
                error(['      Error: descent(''dyn''): ' ...
                    'i must be a positive finite scalar.']);
            end

            [phi,opts] = local_dynamic_update(phi,dLdphi,opts);
            varargout = {phi,opts};

        otherwise
            error(['      Error: descent: ' ...
                'unknown definition "%s". ' ...
                'Use ''init'' or ''dyn''.'],stage);
    end
end

% =============
% local helpers
% =============
function opts = local_initialize_opts(alg,phi)
%LOCAL_INITIALIZE_OPTS Validate settings and initialize optimizer state.

    if ~isstruct(alg)
        error(['      Error: descent(''init''): ' ...
            'alg must be a structure.']);
    end

    if ~isfield(alg,'method') ...
            || isempty(alg.method)
        alg.method = 2;
    end
    if ~isfield(alg,'i_max') ...
            || isempty(alg.i_max)
        alg.i_max = 500;
    end
    if ~isfield(alg,'clipn') ...
            || isempty(alg.clipn)
        alg.clipn = 0;
    end
    if ~isfield(alg,'wdecay') ...
            || isempty(alg.wdecay)
        alg.wdecay = 0;
    end
    if ~isfield(alg,'lr') ...
            || isempty(alg.lr)
        alg.lr = (alg.method == 1) * 1e-4 + (alg.method == 2) * 1e-2;
    end
    if ~isfield(alg,'lr_scheme') ...
            || isempty(alg.lr_scheme)
        alg.lr_scheme = 'constant';
    end
    if ~isfield(alg,'lr_min') ...
            || isempty(alg.lr_min)
        alg.lr_min = alg.lr;
    end
    if ~isfield(alg,'lr_hold') ...
            || isempty(alg.lr_hold)
        alg.lr_hold = 0.20;
    end

    if ~isnumeric(alg.method) ...
            || ~isscalar(alg.method) ...
            || ~ismember(alg.method,[1 2])
        error(['      Error: descent: alg.method must be ' ...
            '1 or 2.']);
    end
    local_positive_integer(alg.i_max,'alg.i_max');

    if alg.method == 2
        required = {'beta_1','beta_2','vareps'};
        for k = 1:numel(required)
            if ~isfield(alg,required{k}) ...
                    || isempty(alg.(required{k}))
                error(['      Error: descent: alg.%s is required ' ...
                    'for Adam/AdamW.'],required{k});
            end
        end
    end

    if ~isnumeric(alg.lr) ...
            || ~isscalar(alg.lr) ...
            || ~isfinite(alg.lr) ...
            || alg.lr <= 0
        error(['      Error: descent: alg.lr must be a positive ' ...
            'finite scalar.']);
    end
    if alg.method == 2
        if ~isnumeric(alg.beta_1) ...
                || ~isscalar(alg.beta_1) ...
                || ~isfinite(alg.beta_1) ...
                || alg.beta_1 <= 0 ...
                || alg.beta_1 >= 1
            error(['      Error: descent: alg.beta_1 must be a finite ' ...
                'scalar in the open interval (0,1).']);
        end
        if ~isnumeric(alg.beta_2) ...
                || ~isscalar(alg.beta_2) ...
                || ~isfinite(alg.beta_2) ...
                || alg.beta_2 <= 0 ...
                || alg.beta_2 >= 1
            error(['      Error: descent: alg.beta_2 must be a finite ' ...
                'scalar in the open interval (0,1).']);
        end
        if ~isnumeric(alg.vareps) ...
                || ~isscalar(alg.vareps) ...
                || ~isfinite(alg.vareps) ...
                || alg.vareps <= 0
            error(['      Error: descent: alg.vareps must be a positive ' ...
                'finite scalar.']);
        end
    end
    if ~isnumeric(alg.lr_min) ...
            || ~isscalar(alg.lr_min) ...
            || ~isfinite(alg.lr_min) ...
            || alg.lr_min <= 0 ...
            || alg.lr_min > alg.lr
        error(['      Error: descent: alg.lr_min must be positive ' ...
            'and cannot exceed alg.lr.']);
    end
    if ~isnumeric(alg.lr_hold) ...
            || ~isscalar(alg.lr_hold) ...
            || ~isfinite(alg.lr_hold) ...
            || alg.lr_hold < 0 ...
            || alg.lr_hold >= 1
        error(['      Error: descent: alg.lr_hold must be in ' ...
            'the interval [0,1).']);
    end
    if ~isnumeric(alg.wdecay) ... 
            || ~isscalar(alg.wdecay) ...
            || ~isfinite(alg.wdecay) ...
            || alg.wdecay < 0
        error(['      Error: descent: alg.wdecay must be a ' ...
            'nonnegative finite scalar.']);
    end

    scheme = lower(strtrim(char(string(alg.lr_scheme))));
    if ~ismember(scheme,{'constant','cosine','hold_cosine'})
        error(['      Error: descent: alg.lr_scheme must be ' ...
            '''constant'', ''cosine'', or ''hold_cosine''.']);
    end

    opts = struct( ...
        'method',alg.method, ...
        'i_max',double(alg.i_max), ...
        'lr',double(alg.lr), ...
        'lr_min',double(alg.lr_min), ...
        'lr_scheme',scheme, ...
        'lr_hold',double(alg.lr_hold), ...
        'lr_current',NaN, ...
        'clipn',double(alg.clipn), ...
        'wdecay',double(alg.wdecay), ...
        't',0, ...
        'n_phi',local_count_phi(phi));

    if opts.method == 2
        opts.beta_1 = double(alg.beta_1);
        opts.beta_2 = double(alg.beta_2);
        opts.vareps = double(alg.vareps);

        nL = numel(phi.W);
        opts.mW = cell(1,nL);
        opts.vW = cell(1,nL);
        opts.mB = cell(1,nL);
        opts.vB = cell(1,nL);

        for li = 1:nL
            opts.mW{li} = zeros(size(phi.W{li}));
            opts.vW{li} = zeros(size(phi.W{li}));
            opts.mB{li} = zeros(size(phi.b{li}));
            opts.vB{li} = zeros(size(phi.b{li}));
        end
    end
end

function [phi,opts] = local_dynamic_update(phi,dLdphi,opts)
%LOCAL_DYNAMIC_UPDATE Generate one new network proposal.

    local_validate_gradient(dLdphi,phi);

    opts.t = opts.t + 1;
    t = opts.t;

    lr = local_learning_rate(opts,t);
    opts.lr_current = lr;
    phiOld = phi;
    opts.weight_decay_norm = 0;

    nL = numel(phi.W);

    switch opts.method
        case 1
            for li = 1:nL
                gW = dLdphi.W{li};
                gB = dLdphi.b{li};

                nW = norm(gW(:));
                nB = norm(gB(:));
                if nW > 0
                    gW = gW./nW;
                end
                if nB > 0
                    gB = gB./nB;
                end

                phi.W{li} = phi.W{li} - lr*single(gW);
                phi.b{li} = phi.b{li} - lr*single(gB);
            end

        case 2
            beta_1 = opts.beta_1;
            beta_2 = opts.beta_2;
            vareps = opts.vareps;
            wd = opts.wdecay;

            for li = 1:nL
                gW = dLdphi.W{li};
                gB = dLdphi.b{li};

                opts.mW{li} = beta_1 * opts.mW{li} ...
                    + (1 - beta_1)*gW;
                opts.vW{li} = beta_2 * opts.vW{li} ...
                    + (1 - beta_2)*(gW.^2);
                mW_hat = opts.mW{li} / (1 - beta_1^t);
                vW_hat = opts.vW{li} / (1 - beta_2^t);
                stepW = mW_hat ./ (sqrt(vW_hat) + vareps);

                if wd > 0
                    phi.W{li} = phi.W{li} ...
                        - lr * (stepW + wd*phi.W{li});
                else
                    phi.W{li} = phi.W{li} - lr * stepW;
                end

                opts.mB{li} = beta_1 * opts.mB{li} ...
                    + (1 - beta_1) * gB;
                opts.vB{li} = beta_2 * opts.vB{li} ...
                    + (1 - beta_2) * (gB.^2);
                mB_hat = opts.mB{li} / (1 - beta_1^t);
                vB_hat = opts.vB{li} / (1 - beta_2^t);
                stepB = mB_hat ./ (sqrt(vB_hat) + vareps);

                phi.b{li} = phi.b{li} - lr * stepB;
            end

        otherwise
            error(['      Error: descent: ' ...
                'unknown opts.method. ' ...
                'Use 1 or 2.']);
    end
    opts.update_norm = local_phi_difference_norm(phi,phiOld);
    opts.relative_update_norm = opts.update_norm / ...
        max(local_phi_norm(phiOld),eps);
    if opts.method == 2
        opts.moment1_norm = local_cell_pair_norm(opts.mW,opts.mB);
        opts.moment2_norm = local_cell_pair_norm(opts.vW,opts.vB);
        opts.weight_decay_norm = lr * opts.wdecay ...
            * local_phi_weight_norm(phiOld);
    else
        opts.moment1_norm = NaN;
        opts.moment2_norm = NaN;
    end
end

function n = local_phi_norm(phi)
%LOCAL_PHI_NORM 
    n = local_cell_pair_norm(phi.W,phi.b);
end

function n = local_phi_weight_norm(phi)
%LOCAL_PHI_WEIGHT_NORM 
    n2 = 0; 
    for k = 1:numel(phi.W)
        n2 = n2 + sum(double(phi.W{k}(:)).^2); 
    end
    n = sqrt(n2);
end

function n = local_phi_difference_norm(a,b)
%LOCAL_PHI_DIFFERENCE_NORM 
    n2 = 0; 
    for k = 1:numel(a.W) 
        n2 = n2 + sum(double(a.W{k}(:)-b.W{k}(:)).^2) ...
            + sum(double(a.b{k}(:)-b.b{k}(:)).^2); 
    end
    n = sqrt(n2);
end

function n = local_cell_pair_norm(A,B)
%LOCAL_CELL_PAIR_NORM 
    n2 = 0; 
    for k = 1:numel(A)
        n2 = n2+sum(double(A{k}(:)).^2) + ...
            sum(double(B{k}(:)).^2);
    end
    n = sqrt(n2);
end

function lr = local_learning_rate(opts,t)
%LOCAL_LEARNING_RATE Scheduled rate for optimizer update t.

    if ~isnumeric(t) ...
            || ~isscalar(t) ...
            || ~isfinite(t) ...
            || t < 1
        error(['      Error: descent: ' ...
            'optimizer update t must be a ' ...
            'positive finite scalar.']);
    end

    maxUpdates = max(double(opts.i_max) - 1,1);
    progress = (double(t) - 1)/max(maxUpdates - 1,1);
    progress = min(max(progress,0),1);

    switch lower(opts.lr_scheme)
        case 'constant'
            lr = opts.lr;

        case 'cosine'
            lr = opts.lr_min ...
                + 0.5 * (opts.lr - opts.lr_min) * (1 + cos(pi*progress));

        case 'hold_cosine'
            if progress <= opts.lr_hold
                lr = opts.lr;
            else
                decayProgress = (progress - opts.lr_hold) ...
                    / max(1 - opts.lr_hold,eps);
                decayProgress = min(max(decayProgress,0),1);
                lr = opts.lr_min ...
                    + 0.5 * (opts.lr - opts.lr_min) ...
                    * (1 + cos(pi * decayProgress));
            end

        otherwise
            error(['      Error: descent: ' ...
                'unknown learning-rate ' ...
                'scheme "%s".'],opts.lr_scheme);
    end
end

function n = local_count_phi(phi)
%LOCAL_COUNT_PHI Count all network weights and biases.

    n = 0;
    for li = 1:numel(phi.W)
        n = n + numel(phi.W{li}) + numel(phi.b{li});
    end
end

function local_validate_gradient(g,phi)
%LOCAL_VALIDATE_GRADIENT Validate dL/dphi against phi.

    if ~isstruct(g) ...
            || ~isfield(g,'W') ...
            || ~isfield(g,'b') ...
            || ~iscell(g.W) ...
            || ~iscell(g.b)
        error(['      Error: descent(''dyn''): ' ...
            'dLdphi must contain ' ...
            'cell arrays W and b.']);
    end
    if numel(g.W) ~= numel(phi.W) ...
            || numel(g.b) ~= numel(phi.b)
        error(['      Error: descent(''dyn''): ' ...
            'dLdphi and phi have ' ...
            'different numbers of layers.']);
    end

    for li = 1:numel(phi.W)
        if ~isequal(size(g.W{li}),size(phi.W{li})) ...
                || ~isequal(size(g.b{li}),size(phi.b{li}))
            error(['      Error: descent(''dyn''): ' ...
                'gradient size mismatch in layer %d.'],li);
        end
        if any(~isfinite(g.W{li}),'all') ...
                || any(~isfinite(g.b{li}),'all')
            error(['      Error: descent(''dyn''): ' ...
                'non-finite FFN gradient ' ...
                'received in layer %d.'],li);
        end
    end
end

function phi = local_validate_phi_struct(phi)
%LOCAL_VALIDATE_PHI_STRUCT Validate canonical network-variable structure.

    if ~isstruct(phi)
        error(['      Error: descent: ' ...
            'phi must be a structure.']);
    end

    req = {'W','b','layers','tf'};
    for k = 1:numel(req)
        if ~isfield(phi,req{k})
            error(['      Error: descent: ' ...
                'phi.%s is missing.'],req{k});
        end
    end

    if ~iscell(phi.W) ...
            || ~iscell(phi.b)
        error(['      Error: descent: ' ...
            'phi.W and phi.b must be cell arrays.']);
    end

    phi.layers = double(phi.layers(:).');
    nL = numel(phi.layers) - 1;
    nH = nL - 1;

    if nH < 1 || nH > 5 ...
            || numel(phi.W) ~= nL ...
            || numel(phi.b) ~= nL
        error(['      Error: descent: ' ...
            'inconsistent phi architecture.']);
    end

    phi.tf = local_parse_transfer_functions(phi.tf,nH);

    for li = 1:nL
        ni = phi.layers(li);
        no = phi.layers(li+1);
        if ~isequal(size(phi.W{li}),[no ni]) ...
                || ~isequal(size(phi.b{li}),[no 1])
            error(['      Error: descent: ' ...
                'inconsistent phi dimensions ' ...
                'in layer %d.'],li);
        end
    end
end

function value = local_positive_integer(value,name)
%LOCAL_POSITIVE_INTEGER Validate and return a positive integer scalar.

    if ~isnumeric(value) ...
            || ~isscalar(value) ...
            || ~isfinite(value) ...
            || value < 1 ...
            || mod(value,1) ~= 0
        error(['      Error: descent: ' ...
            '%s must be a positive integer.'],name);
    end
    value = double(value);
end

function tf = local_parse_transfer_functions(tfin,nH)
%LOCAL_PARSE_TRANSFER_FUNCTIONS Canonicalize hidden transfer functions.

    while iscell(tfin) ...
            && isscalar(tfin)
        tfin = tfin{1};
    end

    if ischar(tfin)
        s = strtrim(lower(tfin));
        if ~ismember(s,{'tanh','relu'})
            error(['      Error: descent: ' ...
                'transfer functions must be ' ...
                '''tanh'' or ''relu''.']);
        end
        tf = repmat({s},1,nH);

    elseif isstring(tfin)
        tf = cell(1,numel(tfin));
        for k = 1:numel(tfin)
            s = strtrim(lower(char(tfin(k))));
            if ~ismember(s,{'tanh','relu'})
                error(['      Error: descent: ' ...
                    'transfer functions must ' ...
                    'be ''tanh'' or ''relu''.']);
            end
            tf{k} = s;
        end

    elseif iscell(tfin)
        tf = cell(1,numel(tfin));
        for k = 1:numel(tfin)
            x = tfin{k};
            while iscell(x) ...
                    && isscalar(x)
                x = x{1};
            end
            if ~(ischar(x) ...
                    || (isstring(x) ...
                    && isscalar(x)))
                error(['      Error: descent: ' ...
                    'transfer-function entries ' ...
                    'must be text scalars.']);
            end
            s = strtrim(lower(char(x)));
            if ~ismember(s,{'tanh','relu'})
                error(['      Error: descent: ' ...
                    'transfer functions must ' ...
                    'be ''tanh'' or ''relu''.']);
            end
            tf{k} = s;
        end
    else
        error(['      Error: descent: ' ...
            'net.tf must be text or a ' ...
            'cell/string array.']);
    end

    if isscalar(tf) ...
            && nH > 1
        tf = repmat(tf,1,nH);
    end
    if numel(tf) ~= nH
        error(['      Error: descent: ' ...
            'number of transfer functions must ' ...
            'match the number of hidden layers.']);
    end
end
