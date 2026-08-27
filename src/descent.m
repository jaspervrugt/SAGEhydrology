function [phi,opts] = descent(stage,alg,x,dLdphi,i,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DESCENT Initializes and updates feedforward-network decision variables.
%
% SYNOPSIS:
%   [phi,opts] = descent('init',alg,net)
%   [phi,opts] = descent('dyn',alg,phi,dLdphi,i,opts)
%
%   stage       operation
%                'init' initialize phi and optimizer state
%                'dyn'  generate a new phi proposal
%
%   alg         descent-algorithm structure
%    .method      1 normalized gradient descent, 2 Adam/AdamW
%    .i_max       number of evaluated SAGE iterations
%    .clipn       global gradient-norm clipping threshold; clipping is
%                 applied upstream by ffn_theta (default 0 = off)
%    .wdecay      AdamW weight decay (default 0 = off)
%    .lr          initial learning rate (> 0)
%    .lr_min      minimum/final learning rate (default = lr)
%    .lr_scheme   'constant' (default), 'cosine', or 'hold_cosine'
%    .lr_hold     fraction of optimizer updates held at lr for
%                 'hold_cosine' (default 0.20; 0 <= lr_hold < 1)
%
%   net         network specification used with stage = 'init'
%    .r           number of input attributes
%    .d           number of output hydrologic parameters
%    .h           hidden-layer widths
%    .tf          hidden-layer transfer functions
%    .seed        optional nonnegative integer RNG seed (default 0)
%
%   phi         network weights and biases used with stage = 'dyn'
%    .W{1:nL}    weight matrices
%    .b{1:nL}    bias vectors
%    .layers      layer sizes
%    .tf          hidden-layer transfer functions
%
%   dLdphi     gradient of loss with respect to phi
%   i           outer SAGE iteration number; retained for diagnostics
%   opts        internal optimizer state
%
% OUTPUT:
%   phi         initialized or updated network variables
%   opts        optimizer settings/state
%    .t           number of completed optimizer updates
%    .lr_current  learning rate used for the latest proposal
%    .n_phi       total number of network weights and biases
%
% NOTES:
%   1. descent owns creation and updating of phi.
%   2. ffn_theta only evaluates phi -> nTheta and computes dL/dphi.
%   3. The first dynamic call is optimizer update t = 1, even when it
%      occurs at outer SAGE iteration i = 2.
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
            net = x;
            phi = local_initialize_phi(net);
            opts = local_initialize_opts(alg,phi);

        case 'dyn'
            phi = local_validate_phi_struct(x);

            if isempty(dLdphi)
                error(['      Error: descent(''dyn''): ' ...
                    'dLdphi is empty. Use ' ...
                    'descent(''init'',alg,net) to create ' ...
                    'the initial phi.']);
            end
            if isempty(opts) ...
                    || ~isstruct(opts)
                error(['      Error: descent(''dyn''): ' ...
                    'opts must be the state returned ' ...
                    'by descent(''init'',alg,net).']);
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

        otherwise
            error(['      Error: descent: ' ...
                'unknown definition "%s". ' ...
                'Use ''init'' or ''dyn''.'],stage);
    end
end

function phi = local_initialize_phi(net)
%LOCAL_INITIALIZE_PHI Initialize ANN weights and biases.

    if ~isstruct(net)
        error(['      Error: descent(''init''): ' ...
            'net must be a structure.']);
    end

    if ~isfield(net,'h') ...
            && isfield(net,'ann') ...
            && isstruct(net.ann) ...
            && isfield(net.ann,'h')
        net.h = net.ann.h;
    end
    if ~isfield(net,'tf') ...
            && isfield(net,'ann') ...
            && isstruct(net.ann) ...
            && isfield(net.ann,'tf')
        net.tf = net.ann.tf;
    end

    req = {'r','d','h','tf'};
    for k = 1:numel(req)
        if ~isfield(net,req{k}) ...
                || isempty(net.(req{k}))
            error(['      Error: descent(''init''): ' ...
                'net.%s is missing or empty.'],req{k});
        end
    end

    r = local_positive_integer(net.r,'net.r');
    d = local_positive_integer(net.d,'net.d');
    h = local_parse_hidden_widths(net.h);

    if numel(h) > 5
        error(['      Error: descent(''init''): ' ...
            'at most five hidden layers are supported.']);
    end

    tf = local_parse_transfer_functions(net.tf,numel(h));
    layers = [double(r),h,double(d)];
    nL = numel(layers)-1;
    nH = numel(h);

    if isfield(net,'seed') ...
            && ~isempty(net.seed)
        seed = net.seed;
    else
        seed = 0;
    end
    if ~isnumeric(seed) ...
            || ~isscalar(seed) ...
            || ~isfinite(seed) ...
            || seed < 0 ...
            || mod(seed,1) ~= 0
        error(['      Error: descent(''init''): ' ...
            'net.seed must be a ' ...
            'nonnegative integer scalar.']);
    end

    fprintf(['... Initializing network ' ...
        'weights and biases']);
    rng(double(seed),'twister');

    phi = struct();
    phi.W = cell(1,nL);
    phi.b = cell(1,nL);
    phi.layers = layers;
    phi.tf = tf;

    for li = 1:nL
        ni = layers(li);
        no = layers(li+1);

        lim = sqrt(6/(ni+no));
        W = (2*rand(no,ni)-1)*lim;

        if li <= nH ...
                && strcmpi(tf{li},'relu')
            W = randn(no,ni)*sqrt(2/ni);
        end

        phi.W{li} = W;
        phi.b{li} = zeros(no,1);
    end

    fprintf(' ... Done\n');
end

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
        alg.lr = (alg.method==1)*1e-4 + (alg.method==2)*1e-2;
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

    if ~isnumeric(alg.lr) ...
            || ~isscalar(alg.lr) ...
            || ~isfinite(alg.lr) ...
            || alg.lr <= 0
        error(['      Error: descent: alg.lr must be a positive ' ...
            'finite scalar.']);
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
        opts.beta_1 = 0.9;
        opts.beta_2 = 0.999;
        opts.vareps = 1e-8;

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

                opts.mW{li} = beta_1*opts.mW{li} ...
                    + (1-beta_1)*gW;
                opts.vW{li} = beta_2*opts.vW{li} ...
                    + (1-beta_2)*(gW.^2);
                mW_hat = opts.mW{li}/(1-beta_1^t);
                vW_hat = opts.vW{li}/(1-beta_2^t);
                stepW = mW_hat./(sqrt(vW_hat)+vareps);

                if wd > 0
                    phi.W{li} = phi.W{li} ...
                        - lr*(stepW + wd*phi.W{li});
                else
                    phi.W{li} = phi.W{li} - lr*stepW;
                end

                opts.mB{li} = beta_1*opts.mB{li} ...
                    + (1-beta_1)*gB;
                opts.vB{li} = beta_2*opts.vB{li} ...
                    + (1-beta_2)*(gB.^2);
                mB_hat = opts.mB{li}/(1-beta_1^t);
                vB_hat = opts.vB{li}/(1-beta_2^t);
                stepB = mB_hat./(sqrt(vB_hat)+vareps);

                phi.b{li} = phi.b{li} - lr*stepB;
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
    n = local_cell_pair_norm(phi.W,phi.b);
end
function n = local_phi_weight_norm(phi)
    n2 = 0; 
    for k = 1:numel(phi.W)
        n2 = n2 + sum(double(phi.W{k}(:)).^2); 
    end
    n = sqrt(n2);
end
function n = local_phi_difference_norm(a,b)
    n2 = 0; 
    for k = 1:numel(a.W) 
        n2 = n2 + sum(double(a.W{k}(:)-b.W{k}(:)).^2) ...
            + sum(double(a.b{k}(:)-b.b{k}(:)).^2); 
    end
    n = sqrt(n2);
end
function n = local_cell_pair_norm(A,B)
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

    maxUpdates = max(double(opts.i_max)-1,1);
    progress = (double(t)-1)/max(maxUpdates-1,1);
    progress = min(max(progress,0),1);

    switch lower(opts.lr_scheme)
        case 'constant'
            lr = opts.lr;

        case 'cosine'
            lr = opts.lr_min ...
                + 0.5*(opts.lr-opts.lr_min)*(1+cos(pi*progress));

        case 'hold_cosine'
            if progress <= opts.lr_hold
                lr = opts.lr;
            else
                decayProgress = (progress-opts.lr_hold) ...
                    / max(1-opts.lr_hold,eps);
                decayProgress = min(max(decayProgress,0),1);
                lr = opts.lr_min ...
                    + 0.5*(opts.lr-opts.lr_min) ...
                    *(1+cos(pi*decayProgress));
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
                'non-finite ANN gradient ' ...
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
    nL = numel(phi.layers)-1;
    nH = nL-1;

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

function h = local_parse_hidden_widths(hin)
%LOCAL_PARSE_HIDDEN_WIDTHS Convert hidden widths to a numeric row vector.

    while iscell(hin) ...
            && isscalar(hin)
        hin = hin{1};
    end

    if isnumeric(hin)
        h = double(hin(:).');
    elseif isstring(hin)
        if isscalar(hin)
            h = str2num(char(hin)); %#ok<ST2NM>
        else
            h = str2double(hin(:).');
        end
    elseif ischar(hin)
        h = str2num(hin); %#ok<ST2NM>
    elseif iscell(hin)
        h = nan(1,numel(hin));
        for k = 1:numel(hin)
            x = hin{k};
            while iscell(x) ...
                    && isscalar(x)
                x = x{1};
            end
            if isnumeric(x) ...
                    && isscalar(x)
                h(k) = double(x);
            elseif (isstring(x) ...
                    && isscalar(x)) || ischar(x)
                h(k) = str2double(x);
            else
                error(['      Error: descent: ' ...
                    'hidden-layer widths ' ...
                    'must be numeric or text scalars.']);
            end
        end
    else
        error(['      Error: descent: ' ...
            'net.h must be numeric, text, ' ...
            'or a cell array.']);
    end

    h = double(h(:).');
    if isempty(h) ...
            || any(~isfinite(h)) ...
            || any(h < 1) ...
            || any(mod(h,1) ~= 0)
        error(['      Error: descent: ' ...
            'net.h must contain positive ' ...
            'integers.']);
    end
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
            if ~(ischar(x) || (isstring(x) ...
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
