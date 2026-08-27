function varargout = ffn_theta(stage,x1,x2,x3,x4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FFN_THETA  MLP mapping catchment attributes -> normalized parameters (0,1)
% using a feedforward neural network
% 
% SYNOPSIS:
%   nTheta = ffn_theta('eval',phi,A)
%   [nTheta,H] = ffn_theta('eval_info',phi,A)
%   dLdphi = ffn_theta('grad',phi,A,alg,G)
%
%   stage   task that should be completed
%    'eval' to evaluate network and return normalized parameter values
%       x1: phi = struct with fields:
%           phi.W{1:L} = weight matrices
%           phi.b{1:L} = bias vectors
%           phi.layers = layer sizes
%           phi.tf{1:nH} = transfer function(s) hidden layers
%       x2: A = r x K matrix of catchment attributes
%    'grad' to compute gradient of network
%       x1: phi = struct with fields:
%           phi.W{1:L} = weight matrices
%           phi.b{1:L} = bias vectors
%           phi.layers = layer sizes
%           phi.tf{1:nH} = transfer function(s) hidden layers
%       x2: A = r x K matrix of catchment attributes
%       x3: alg = optimizer settings; alg.clipn controls gradient clipping
%       x4: G = d x K matrix of gradients of hydrologic model
%    'eval_info' to evaluate network and also return hidden-layer
%       activations for information-bottleneck diagnostics
%       x1: phi = canonical network structure
%       x2: A = r x K matrix of catchment attributes
%
% OUTPUT ARGUMENTS
%    'eval': nTheta = d x K matrix of normalized parameter values
%    'grad': dLdphi = structure with mean gradients w.r.t. all W and b
%    'eval_info':
%       nTheta = d x K matrix of normalized parameter values
%       H      = 1 x nH cell array; H{ell} is K x h_ell and contains
%                post-activation hidden-layer values for all basins
%
% NOTES
%   1. Output layer uses sigmoid to enforce unit-cube parameterization.
%   2. Supports one to five hidden layers.
%   3. Network-variable initialization is owned by descent('init').
%   4. ffn_theta only turns phi -> nTheta and backpropagates G -> dLdphi.
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
        case 'eval'
            if nargin < 3
                error(['      Error: ffn_theta(''eval''): ' ...
                    'requires inputs phi and A.']);
            end
    
            phi = local_validate_phi_struct(x1);
            A = x2;
    
            nTheta = forward_struct(phi,A);
            varargout = {nTheta};

        case 'grad'
            if nargin < 5
                error(['      Error: ffn_theta(''grad''): ' ...
                    'requires inputs phi, A, alg, and G.']);
            end
    
            phi = local_validate_phi_struct(x1);
            clipn = 0;
            if isstruct(x3) && isfield(x3,'clipn') ...
                    && ~isempty(x3.clipn)
                clipn = x3.clipn;
            end
            A = x2;
            G = x4;
    
            if ~isnumeric(A) ...
                    || ~ismatrix(A)
                error(['      Error: ffn_theta(''grad''): ' ...
                    'A must be a numeric matrix.']);
            end
            if ~isnumeric(G) ...
                    || ~ismatrix(G)
                error(['      Error: ffn_theta(''grad''): ' ...
                    'G must be a numeric matrix.']);
            end
    
            % Remove basins with invalid attributes or invalid hydrologic gradients.
            % These basins have no valid training contribution and should not enter
            % the ANN-gradient average.
            good = all(isfinite(A),1) ...
                & all(isfinite(G),1);

            if ~all(good)
                A = A(:,good);
                G = G(:,good);
            end

            % How many watersheds for training
            K_t = size(A,2);

            if K_t == 0
                dLdphi = local_zero_like_phi(phi);
                gradInfo = struct('normPre',0, ...
                    'normPost',0,'clipped',false);
                varargout = {dLdphi,gradInfo};
                return
            end
    
            if size(A,1) ~= phi.layers(1)
                error(['      Error: ffn_theta(''grad''): ' ...
                    'size(A,1) = %d, but ' ...
                    'expected %d from phi.layers(1).'], ...
                    size(A,1),phi.layers(1));
            end
            if size(G,1) ~= phi.layers(end)
                error(['      Error: ffn_theta(''grad''): ' ...
                    'size(G,1) = %d, but ' ...
                    'expected %d from phi.layers(end).'], ...
                    size(G,1),phi.layers(end));
            end
            if size(G,2) ~= K_t
                error(['      Error: ffn_theta(''grad''): ' ...
                    'size(G,2) must match ' ...
                    'size(A,2).']);
            end
    
            % forward pass cache
            [Y,cache] = forward_cache(phi,A);
    
            % backward pass
            nL = numel(phi.W);
            nH = numel(phi.tf);
    
            dLdphi.W = cell(1,nL);
            dLdphi.b = cell(1,nL);
    
            % output layer: sigmoid
            dsg = Y .* (1 - Y);             % d x K_t
            dZ = G .* dsg;                  % d x K_t
    
            for li = nL:-1:1
                Aprev = cache.A{li};        % input to layer li
                dLdphi.W{li} = dZ * Aprev.';
                dLdphi.b{li} = sum(dZ,2);
    
                if li > 1
                    dAprev = phi.W{li}.' * dZ;
    
                    % apply derivative only if previous layer is hidden
                    if (li-1) <= nH
                        Zprev = cache.Z{li-1};
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
    
        case 'eval_info'            
            if nargin < 3
                error(['      Error: ffn_theta(''eval_info''): ' ...
                    'requires inputs phi and A.']);
            end

            phi = local_validate_phi_struct(x1);
            A = x2;

            if ~isnumeric(A) ...
                    || ~ismatrix(A)
                error(['      Error: ffn_theta(''eval_info''): ' ...
                    'A must be a numeric matrix.']);
            end
            if size(A,1) ~= phi.layers(1)
                error(['      Error: ffn_theta(''eval_info''): ' ...
                    'size(A,1) = %d, but expected ' ...
                    '%d from phi.layers(1).'], ...
                    size(A,1),phi.layers(1));
            end

            [nTheta,H] = forward_info(phi,A);
            varargout = {nTheta,H};

        otherwise
            error(['      Error: ffn_theta: unknown stage "%s". ' ...
                'Use ''eval'' or ''grad''.'],stage);
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


function [nTheta,H] = forward_info(phi,A)
%FORWARD_INFO Forward pass with hidden-layer activations for diagnostics.
%
%   [nTheta,H] = forward_info(phi,A)
%
%   H{ell} is K x h_ell (basins x neurons), matching the orientation
%   expected by information_bottleneck_SAGE. The ordinary SAGE network
%   orientation remains neurons x basins internally.

    nL = numel(phi.W);
    nH = numel(phi.tf);

    H = cell(1,nH);

    a = A;
    for li = 1:nL
        z = phi.W{li}*a + phi.b{li};

        if li <= nH
            a = activate(z,phi.tf{li});
            H{li} = a.';
        else
            a = sigmoid_dl(z);
        end
    end

    nTheta = a;
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
%LOCAL_POSTPROCESS_GRAD Sanitize ANN gradients and optionally clip
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
                'non-finite ANN gradient ' ...
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
                'clipn = %.3e; scaling ANN ' ...
                'gradient by %.3e\n'], ...
                gn,clipn,scale);
        end
    end
    info.normPost = local_gradient_norm(dLdphi);

end

function gn = local_gradient_norm(g)
    gn2 = 0;
    for li = 1:numel(g.W)
        gn2 = gn2 + sum(double(g.W{li}(:)).^2) ...
            + sum(double(g.b{li}(:)).^2);
    end
    gn = sqrt(gn2);
end

function dLdphi = local_zero_like_phi(phi)
%LOCAL_ZERO_LIKE_PHI Return zero ANN-gradient structure.

    nL = numel(phi.W);
    
    dLdphi.W = cell(1,nL);
    dLdphi.b = cell(1,nL);
    
    for li = 1:nL
        dLdphi.W{li} = zeros(size(phi.W{li}));
        dLdphi.b{li} = zeros(size(phi.b{li}));
    end
end
