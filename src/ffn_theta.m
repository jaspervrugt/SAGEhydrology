function varargout = ffn_theta(stage,x1,x2,x3,x4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FFN_THETA Map basin attributes to normalized parameters.
%
%  Initializes, evaluates, and differentiates the feedforward network.
%
% SYNOPSIS:
%   [phi,net] = ffn_theta('init',net)
%   [nTheta,cache] = ffn_theta('eval',phi,A)
%   [nTheta,H,cache] = ffn_theta('info',phi,A)
%   dLdphi = ffn_theta('back',phi,cache,alg,G)
%
% INPUT ARGUMENTS:
%   stage           'init', 'eval', 'info', or 'back'
%   x1              stage-specific net or canonical phi structure
%   x2              attributes A, or forward cache for 'back'
%   x3              optimizer settings alg for 'back'
%   x4              model-parameter gradient G for 'back'
%
% OUTPUT ARGUMENTS:
%   varargout       stage-specific network results
%    'init'          initialized phi and network settings
%    'eval'          normalized parameters and forward cache
%    'info'          parameters, hidden activations, and cache
%    'back'          gradient with respect to phi
%
% NOTES:
%   The sigmoid output layer maps hydrologic parameters to (0,1).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 1
        error(['      ffn_theta: ' ...
            'stage input is missing.']);
    end
    if nargin < 5
        x4 = [];
    end
    
    switch lower(stage)

        case 'init'
            if nargin < 2
                error(['      Error: ffn_theta(''init''): ' ...
                    'requires the network specification net.']);
            end
            [phi,net] = local_initialize_phi(x1);
            varargout = {phi,net};

        case 'eval'
            if nargin < 3
                error(['      Error: ffn_theta(''eval''): ' ...
                    'requires inputs phi and A.']);
            end
            phi = local_validate_phi_struct(x1);
            A = x2;
            if nargout > 1
                [nTheta,cache] = forward_cache(phi,A);
                varargout = {nTheta,cache};
            else
                nTheta = forward_struct(phi,A);
                varargout = {nTheta};
            end

        case 'back'

            if nargin < 5
                error(['      Error: ffn_theta(''back''): ' ...
                    'requires inputs phi, cache, alg, and G.']);
            end
            phi = local_validate_phi_struct(x1);
            clipn = 0;
            if isstruct(x3) ...
                    && isfield(x3,'clipn') ...
                    && ~isempty(x3.clipn)
                clipn = x3.clipn;
            end
            cache = x2;
            G = x4;    
            if ~isstruct(cache) ...
                    || ~isfield(cache,'A') ...
                    || ~isfield(cache,'Z') ...
                    || ~isfield(cache,'Y')
                error(['      Error: ffn_theta(''back''): ' ...
                    'cache must be returned by ffn_theta(''eval'') ' ...
                    'or ffn_theta(''info'').']);
            end
            if ~isnumeric(G) ...
                    || ~ismatrix(G)
                error(['      Error: ffn_theta(''back''): ' ...
                    'G must be a numeric matrix.']);
            end
            nL = numel(phi.W);
            if ~iscell(cache.A) || ~iscell(cache.Z) ...
                    || numel(cache.A) ~= nL || numel(cache.Z) ~= nL
                error(['      Error: ffn_theta(''back''): ' ...
                    'cache is incompatible with phi.']);
            end
            K_t = size(G,2);
            if size(cache.Y,1) ~= phi.layers(end) ...
                    || size(cache.Y,2) < K_t
                error(['      Error: ffn_theta(''back''): ' ...
                    'cache does not contain the requested training basins.']);
            end
            if size(G,1) ~= phi.layers(end)
                error(['      Error: ffn_theta(''back''): ' ...
                    'size(G,1) = %d, but expected %d.'], ...
                    size(G,1),phi.layers(end));
            end
            good = all(isfinite(G),1) ...
                & all(isfinite(cache.Y(:,1:K_t)),1);
            for li = 1:nL
                if size(cache.A{li},2) < K_t ...
                        || size(cache.Z{li},2) < K_t
                    error(['      Error: ffn_theta(''back''): ' ...
                        'cache layer %d has too few basin columns.'],li);
                end
                good = good ...
                    & all(isfinite(cache.A{li}(:,1:K_t)),1) ...
                    & all(isfinite(cache.Z{li}(:,1:K_t)),1);
            end
            K_t = nnz(good);
            if K_t == 0
                dLdphi = local_zero_like_phi(phi);
                gradInfo = struct('normPre',0, ...
                    'normPost',0,'clipped',false);
                varargout = {dLdphi,gradInfo};
                return
            end    
            % backward pass
            nH = numel(phi.tf);    
            dLdphi.W = cell(1,nL);
            dLdphi.b = cell(1,nL);    
            % output layer: sigmoid
            Y = cache.Y(:,1:size(G,2));
            Y = Y(:,good);
            G = G(:,good);
            dsg = Y .* (1 - Y);             % d x K_t
            dZ = G .* dsg;                  % d x K_t
    
            for li = nL:-1:1
                Aprev = cache.A{li}(:,1:size(good,2));
                Aprev = Aprev(:,good);      % input to layer li
                dLdphi.W{li} = dZ * Aprev.';
                dLdphi.b{li} = sum(dZ,2);    
                if li > 1
                    dAprev = phi.W{li}.' * dZ;    
                    % apply derivative only if previous layer is hidden
                    if (li-1) <= nH
                        Zprev = cache.Z{li-1}(:,1:size(good,2));
                        Zprev = Zprev(:,good);
                        tfi = phi.tf{li-1};
                        dAct = activate_deriv(Zprev,tfi);
                        dZ = dAprev .* dAct;
                    else
                        dZ = dAprev;
                    end
                end
            end
    
            % average over training watersheds
            dLdphi = struct_divide(dLdphi,K_t);    
            % sanitize NaN/Inf entries and optionally clip global norm
            [dLdphi,flag,gradInfo] = ...
                local_postprocess_grad(dLdphi,clipn); %#ok    
            varargout = {dLdphi,gradInfo};
    
        case 'info'            
            if nargin < 3
                error(['      Error: ffn_theta(''info''): ' ...
                    'requires inputs phi and A.']);
            end
            phi = local_validate_phi_struct(x1);
            A = x2;
            if ~isnumeric(A) ...
                    || ~ismatrix(A)
                error(['      Error: ffn_theta(''info''): ' ...
                    'A must be a numeric matrix.']);
            end
            if size(A,1) ~= phi.layers(1)
                error(['      Error: ffn_theta(''info''): ' ...
                    'size(A,1) = %d, but expected ' ...
                    '%d from phi.layers(1).'], ...
                    size(A,1),phi.layers(1));
            end
            [nTheta,H,cache] = forward_info(phi,A);
            varargout = {nTheta,H,cache};

        otherwise
            error(['      Error: ffn_theta: unknown stage "%s". ' ...
                'Use ''init'', ''eval'', ''info'', or ''back''.'],stage);
    end

end

function [phi,net] = local_initialize_phi(net)
%LOCAL_INITIALIZE_PHI Initialize FFN weights and biases.

    if ~isstruct(net)
        error(['      Error: ffn_theta(''init''): ' ...
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
            error(['      Error: ffn_theta(''init''): ' ...
                'net.%s is missing or empty.'],req{k});
        end
    end

    r = local_positive_integer(net.r,'net.r');
    d = local_positive_integer(net.d,'net.d');
    h = local_parse_hidden_widths(net.h);

    if numel(h) > 5
        error(['      Error: ffn_theta(''init''): ' ...
            'at most five hidden layers are supported.']);
    end

    tf = local_parse_transfer_functions(net.tf,numel(h));
    layers = [r,h,d];
    nL = numel(layers)-1;
    nH = numel(h);

    if isfield(net,'seed') && ~isempty(net.seed)
        seed = net.seed;
    else
        seed = 0;
    end
    if ~isnumeric(seed) ...
            || ~isscalar(seed) ...
            || ~isfinite(seed) ...
            || seed < 0 ...
            || mod(seed,1) ~= 0
        error(['      Error: ffn_theta(''init''): ' ...
            'net.seed must be a nonnegative integer scalar.']);
    end

    fprintf('... Initializing FFN weights and biases');
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

        if li <= nH && strcmpi(tf{li},'relu')
            W = randn(no,ni)*sqrt(2/ni);
        end

        phi.W{li} = W;
        phi.b{li} = zeros(no,1);
    end

    net.r = r;
    net.d = d;
    net.h = h;
    net.tf = tf;
    net.seed = double(seed);
    net.l = sum(layers(1:end-1).*layers(2:end) + layers(2:end));

    fprintf(' ... Done\n');
end

function value = local_positive_integer(value,name)
%LOCAL_POSITIVE_INTEGER Validate and return a positive integer scalar.

    if ~isnumeric(value) ...
            || ~isscalar(value) ...
            || ~isfinite(value) ...
            || value < 1 ...
            || mod(value,1) ~= 0
        error(['      Error: ffn_theta(''init''): ' ...
            '%s must be a positive integer.'],name);
    end
    value = double(value);
end

function h = local_parse_hidden_widths(hin)
%LOCAL_PARSE_HIDDEN_WIDTHS Convert hidden widths to a numeric row vector.

    while iscell(hin) && isscalar(hin)
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
            while iscell(x) && isscalar(x)
                x = x{1};
            end
            if isnumeric(x) && isscalar(x)
                h(k) = double(x);
            elseif (isstring(x) && isscalar(x)) || ischar(x)
                h(k) = str2double(x);
            else
                error(['      Error: ffn_theta(''init''): ' ...
                    'hidden-layer widths must be numeric or text scalars.']);
            end
        end
    else
        error(['      Error: ffn_theta(''init''): ' ...
            'net.h must be numeric, text, or a cell array.']);
    end

    h = double(h(:).');
    if isempty(h) ...
            || any(~isfinite(h)) ...
            || any(h < 1) ...
            || any(mod(h,1) ~= 0)
        error(['      Error: ffn_theta(''init''): ' ...
            'net.h must contain positive integers.']);
    end
end

function tf = local_parse_transfer_functions(tfin,nH)
%LOCAL_PARSE_TRANSFER_FUNCTIONS Convert tfin to canonical 1 x nH cell
% array of chars

    % unwrap nested scalar cells, e.g. {{'tanh'}} or {{{'tanh','relu'}}}
    while iscell(tfin) ...
            && isscalar(tfin)
        tfin = tfin{1};
    end
    
    if ischar(tfin)
    
        % Accept single transfer function only, e.g. 'tanh'
        s = strtrim(lower(tfin));
    
        if ismember(s,{'tanh','relu'})
            tf = repmat({s},1,nH);
        else
            error(['      Error: ffn_theta: ' ...
                'net.tf = ''%s'' is not valid. ' ...
                'Use ''tanh'', ''relu'', ' ...
                'or a cell array such as ' ...
                '{''tanh'',''relu''}.'],s);
        end
    
    elseif isstring(tfin)
    
        if isscalar(tfin)
            s = strtrim(lower(char(tfin)));
            if ismember(s,{'tanh','relu'})
                tf = repmat({s},1,nH);
            else
                error(['      Error: ffn_theta: ' ...
                    'net.tf = "%s" is not valid. ' ...
                    'Use "tanh", "relu", or a ' ...
                    'string array such as ' ...
                    '["tanh","relu"].'],s);
            end
        else
            tf = cell(1,numel(tfin));
            for i = 1:numel(tfin)
                s = strtrim(lower(char(tfin(i))));
                if ~ismember(s,{'tanh','relu'})
                    error(['      Error: ffn_theta: ' ...
                        'net.tf entry %d must be ' ...
                        '''tanh'' or ''relu''.'],i);
                end
                tf{i} = s;
            end
        end
    
    elseif iscell(tfin)
    
        tf = cell(1,numel(tfin));
    
        for i = 1:numel(tfin)
            xi = tfin{i};
    
            while iscell(xi) ...
                && isscalar(xi)
                xi = xi{1};
            end
    
            if ischar(xi) ...
                    || (isstring(xi) ...
                    && isscalar(xi))
                s = strtrim(lower(char(xi)));
            else
                error(['      Error: ffn_theta: ' ...
                    'net.tf cell entries must be ' ...
                    'text scalars.']);
            end
    
            if ~ismember(s,{'tanh','relu'})
                error(['      Error: ffn_theta: ' ...
                    'net.tf entry %d must be ' ...
                    '''tanh'' or ''relu''.'],i);
            end
    
            tf{i} = s;
        end
    
    else
        error(['      Error: ffn_theta: ' ...
            'net.tf must be char, ' ...
            'string, or cell ' ...
            'array of text.']);
    end
    
    % expand scalar specification to all hidden layers
    if isscalar(tf) ...
            && nH > 1
        tf = repmat(tf,1,nH);
    end
    
    if numel(tf) ~= nH
        error(['      Error: ffn_theta: number of ' ...
            'transfer functions in net.tf ' ...
            'must match the number of ' ...
            'hidden layers in net.h.']);
    end

end

function phi = local_validate_phi_struct(phi)
%LOCAL_VALIDATE_PHI_STRUCT Ensure phi has a consistent canonical format

    if ~isstruct(phi)
        error(['      Error: ffn_theta: ' ...
            'phi must be a structure.']);
    end
    
    req = {'W','b','layers','tf'};
    for i = 1:numel(req)
        if ~isfield(phi,req{i})
            error(['      Error: ffn_theta: ' ...
                'phi.%s is missing.'], ...
                req{i});
        end
    end
    
    if ~iscell(phi.W) ...
            || ~iscell(phi.b)
        error(['      Error: ffn_theta: ' ...
            'phi.W and phi.b must be ' ...
            'cell arrays.']);
    end
    
    if ~isnumeric(phi.layers) ...
            || isempty(phi.layers)
        error(['      Error: ffn_theta: ' ...
            'phi.layers must be numeric.']);
    end
    phi.layers = double(phi.layers(:).');
    
    nL = numel(phi.layers) - 1;
    if nL < 1
        error(['      Error: ffn_theta: ' ...
            'phi.layers must define at ' ...
            'least 1 layer.']);
    end
    
    if numel(phi.W) ~= nL || numel(phi.b) ~= nL
        error(['      Error: ffn_theta: ' ...
            'inconsistent phi: numel(W), ' ...
            'numel(b), and phi.layers do ' ...
            'not match.']);
    end
    
    nH = nL - 1;
    if nH < 1 || nH > 5
        error(['      Error: ffn_theta: ' ...
            'phi implies %d hidden layers; ' ...
            'only 1 to 5 hidden layers are ' ...
            'supported.'],nH);
    end
    
    phi.tf = local_parse_transfer_functions(phi.tf,nH);
    
    for li = 1:nL
        Wi = phi.W{li};
        bi = phi.b{li};
    
        if ~isnumeric(Wi) || ~isnumeric(bi)
            error(['      Error: ffn_theta: ' ...
                'phi.W{%d} and phi.b{%d} ' ...
                'must be numeric.'],li,li);
        end
    
        ni = phi.layers(li);
        no = phi.layers(li+1);
    
        if ~isequal(size(Wi),[no ni])
            error(['      Error: ffn_theta: ' ...
                'phi.W{%d} has size [%d %d], ' ...
                'expected [%d %d].'], ...
                li,size(Wi,1),size(Wi,2),no,ni);
        end
    
        if ~isequal(size(bi),[no 1])
            error(['      Error: ffn_theta: ' ...
                'phi.b{%d} has size [%d %d], ' ...
                'expected [%d 1].'], ...
                li,size(bi,1),size(bi,2),no);
        end
    end

end


function [nTheta,H,cache] = forward_info(phi,A)
%FORWARD_INFO Forward pass with hidden-layer activations for diagnostics.
%
%   [nTheta,H,cache] = forward_info(phi,A)
%
%   H{ell} is K x h_ell (basins x neurons), matching the orientation
%   expected by sage_information_bottleneck. The ordinary SAGE network
%   orientation remains neurons x basins internally.

    nH = numel(phi.tf);
    H = cell(1,nH);
    [nTheta,cache] = forward_cache(phi,A);
    for li = 1:nH
        H{li} = activate(cache.Z{li},phi.tf{li}).';
    end
end


function [Y,cache] = forward_cache(phi,A)
%FORWARD_CACHE Forward pass with storage of intermediate results
%  SYNOPSIS: [Y,cache] = forward_cache(phi,A)
    
    nL = numel(phi.W);
    nH = numel(phi.tf);
    
    cache.A = cell(1,nL);
    cache.Z = cell(1,nL);
    
    a = A;
    for li = 1:nL
        cache.A{li} = a;
        z = phi.W{li}*a + phi.b{li};
        cache.Z{li} = z;
    
        if li <= nH
            a = activate(z,phi.tf{li});
        else
            a = sigmoid_dl(z);
        end
    end
    
    Y = a;
    cache.Y = Y;

end

function nTheta = forward_struct(phi,A)
%FORWARD_STRUCT Forward pass using phi struct
%  SYNOPSIS: nTheta = forward_struct(phi,A)

    if ~isnumeric(A) ...
            || ~ismatrix(A)
        error(['      Error: ffn_theta(''eval''): ' ...
            'A must be a numeric matrix.']);
    end
    if size(A,1) ~= phi.layers(1)
        error(['      Error: ffn_theta(''eval''): ' ...
            'size(A,1) = %d, but expected ' ...
            '%d from phi.layers(1).'], ...
            size(A,1),phi.layers(1));
    end
    
    nL = numel(phi.W);
    nH = numel(phi.tf);
    
    a = A;
    for li = 1:nL
        a = phi.W{li}*a + phi.b{li};
        if li <= nH
            a = activate(a,phi.tf{li});
        else
            a = sigmoid_dl(a);
        end
    end
    
    nTheta = a;

end

function y = activate(x,tf)
%ACTIVATE Activation function

    switch lower(tf)
        case 'tanh'
            y = tanh(x);
        case 'relu'
            y = max(x,0);
        otherwise
            error(['      Error: ffn_theta: ' ...
                'unknown activation "%s".'], ...
                tf);
    end

end

function y = activate_deriv(x,tf)
%ACTIVATE_DERIV Derivative of activation function

    switch lower(tf)
        case 'tanh'
            y = 1 - tanh(x).^2;
        case 'relu'
            y = (x > 0);
        otherwise
            error(['     Error: ffn_theta: ' ...
                'unknown activation "%s".'], ...
                tf);
    end

end

function y = sigmoid_dl(x)
%SIGMOID_DL Sigmoidal (logistic) transfer function

    y = 1 ./ (1 + exp(-x));

end

function S = struct_divide(S,K)
%STRUCT_DIVIDE Divides each field of structure S by K

    fn = fieldnames(S);
    for i = 1:numel(fn)
        f = fn{i};
        v = S.(f);
        if isstruct(v)
            S.(f) = struct_divide(v,K);
        elseif iscell(v)
            for j = 1:numel(v)
                if isnumeric(v{j})
                    v{j} = v{j} ./ K;
                end
            end
            S.(f) = v;
        elseif isnumeric(v)
            S.(f) = v ./ K;
        end
    end

end

function [dLdphi,flag,info] = local_postprocess_grad(dLdphi,clipn)
%LOCAL_POSTPROCESS_GRAD Sanitize FFN gradients and optionally clip
%their global Euclidean norm

    vareps = 1e-8;
    flag = 0;
    info = struct('normPre',NaN, ...
        'normPost',NaN, ...
        'clipped',false);
    
    % --------------------
    % 1. Replace NaN / Inf
    % --------------------
    for pi = ["W","b"]
        for l = 1:numel(dLdphi.(pi))
            g = dLdphi.(pi){l};
    
            id_nan = isnan(g);
            id_inf = isinf(g);
    
            if any(id_nan,'all')
                fprintf(['      Warning: ffn_theta: ' ...
                    'NaN detected in dLdphi.%s{%d}; ' ...
                    'replacing with %g\n'],pi,l,vareps);
                g(id_nan) = vareps;
                flag = 1;
            end
    
            if any(id_inf,'all')
                fprintf(['      Warning: ffn_theta: ' ...
                    'Inf detected in dLdphi.%s{%d}; ' ...
                    'replacing with %g\n'],pi,l,vareps);
                g(id_inf) = vareps;
                flag = max(flag,2);
            end
    
            dLdphi.(pi){l} = g;
        end
    end
    
    % --------------------------------
    % 2. Optional global norm clipping
    % --------------------------------
    info.normPre = local_gradient_norm(dLdphi);
    if ~isempty(clipn) ...
            && isnumeric(clipn) ...
            && isscalar(clipn) ...
            && isfinite(clipn) ...
            && clipn > 0
        gn2 = 0;
        for li = 1:numel(dLdphi.W)
            gW = dLdphi.W{li};
            gB = dLdphi.b{li};
            gn2 = gn2 + sum(gW(:).^2) + sum(gB(:).^2);
        end
    
        gn = sqrt(gn2);
    
        if ~isfinite(gn)
            warning(['      Warning: ffn_theta: ' ...
                'non-finite FFN gradient ' ...
                'norm after sanitation.']);
            flag = max(flag,3);
        elseif gn > clipn
            info.clipped = true;
            scale = clipn / (gn + vareps);
            for li = 1:numel(dLdphi.W)
                dLdphi.W{li} = dLdphi.W{li} * scale;
                dLdphi.b{li} = dLdphi.b{li} * scale;
            end
            fprintf(['      Warning: ffn_theta: ' ...
                '||dLdphi|| = %.3e exceeds ' ...
            'clipn = %.3e; scaling FFN ' ...
                'gradient by %.3e\n'], ...
                gn,clipn,scale);
        end
    end
    info.normPost = local_gradient_norm(dLdphi);

end

function gn = local_gradient_norm(g)
%LOCAL_GRADIENT_NORM Return the Euclidean norm of network gradients.

    gn2 = 0;
    for li = 1:numel(g.W)
        gn2 = gn2 + sum(double(g.W{li}(:)).^2) ...
            + sum(double(g.b{li}(:)).^2);
    end
    gn = sqrt(gn2);
end

function dLdphi = local_zero_like_phi(phi)
%LOCAL_ZERO_LIKE_PHI Return zero FFN-gradient structure.

    nL = numel(phi.W);
    
    dLdphi.W = cell(1,nL);
    dLdphi.b = cell(1,nL);
    
    for li = 1:nL
        dLdphi.W{li} = zeros(size(phi.W{li}));
        dLdphi.b{li} = zeros(size(phi.b{li}));
    end
end
