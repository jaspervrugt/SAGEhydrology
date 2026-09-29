function info = sage_information_bottleneck(A,H,Y,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INFORMATION_BOTTLENECK_SAGE Estimate information quantities.
%
%  Call the native information estimator, or its MATLAB reference
%  implementation when the native MEX is unavailable.
%
% SYNOPSIS:
%   info = sage_information_bottleneck(A,H,Y,opts)
%   info = sage_information_bottleneck(A,H,Y)
%
% INPUT ARGUMENTS:
%   A               basin-attribute matrix
%   H               hidden-layer representation accepted by the estimator
%   Y               target matrix
%   opts            optional settings; omitted fields use defaults
%    .info_threads   OpenMP thread count; default 1
%    .info_parallel  false forces serial evaluation
%    .nbins          attribute discretization bins; default 10
%    .target_nbins   target discretization bins; default 8
%
% OUTPUT ARGUMENTS:
%   info            information quantities returned by the selected solver
%
% NOTES:
%   The default is one native thread, including inside MATLAB workers.
%   This wrapper creates no figures; visualization belongs to the GUI.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 4 ...
            || isempty(opts)
        opts = struct();
    end

    opts = local_defaults(opts);

    % Prevent nested parallelism by default inside MATLAB workers.
    if ~isfield(opts,'info_threads') ...
            || isempty(opts.info_threads)
        opts.info_threads = 1;   % force serial inside worker
    end

    if exist('sage_information_bottleneck_mex','file') ~= 3
        if exist('sage_information_bottleneck_matlab','file') == 2
            warning('SAGE:Information:MexUnavailable', ...
                ['sage_information_bottleneck_mex ' ...
                'is unavailable; using the ' ...
                 'MATLAB reference implementation.']);
            info = sage_information_bottleneck_matlab(A,H,Y,opts);
            return
        end

        error('SAGE:Information:MexUnavailable', ...
            ['sage_information_bottleneck_mex is unavailable. Run ' ...
             'compile_information in source ' ...
             'MATLAB or include the native ' ...
             'MEX in the deployed SAGE build.']);
    end

    info = sage_information_bottleneck_mex( ...
        double(A),H,double(Y),opts);
end


function opts = local_defaults(opts)
%LOCAL_DEFAULTS Apply defaults for information-bottleneck settings.

    D = struct( ...
        'nbins',10, ...
        'target_nbins',8, ...
        'quantization','quantile', ...
        'ksg',true, ...
        'ksg_k',5, ...
        'ksg_target_indices',1, ...
        'ksg_standardize',true, ...
        'ksg_jitter',1e-10, ...
        'ksg_permutation_n',0, ...
        'ksg_uncertainty_neurons',false, ...
        'ksg_uncertainty_target_indices',1, ...
        'random_seed',1729, ...
        'active_entropy_threshold',0.10, ...
        'info_threads',1, ...
        'info_parallel',true);

    fn = fieldnames(D);
    for i = 1:numel(fn)
        if ~isfield(opts,fn{i}) ...
                || isempty(opts.(fn{i}))
            opts.(fn{i}) = D.(fn{i});
        end
    end
end
