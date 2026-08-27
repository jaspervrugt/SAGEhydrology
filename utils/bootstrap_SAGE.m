function region = bootstrap_SAGE(root,region,dirres)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BOOTSTRAP_SAGE Initialize a standalone SAGEhydrology run.
%
% SYNOPSIS:
%   region = bootstrap_SAGE(root,region,dirres)
%
% INPUT:
%   root       Root directory containing the SAGEhydrology installation
%   region     Region name, alias, short code, or canonical region code
%   dirres     Directory in which SAGE output and results are written
%
% OUTPUT:
%   region     Canonical region code returned by region_helpers
%
% DESCRIPTION:
%   This helper initializes a SAGE run launched outside SAGE-GUI. It adds
%   the required SAGEhydrology source directories to the MATLAB path,
%   normalizes the selected region, and creates the results directory when
%   necessary. The folder containing this function must be added to the
%   MATLAB path before bootstrap_SAGE is called.
%
% EXAMPLE:
%   sageUtils = fullfile(root,'SAGEhydrology','utils');
%   addpath(sageUtils);
%   region = bootstrap_SAGE(root,region,dirres);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Aug. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin ~= 3
        error('bootstrap_SAGE:InvalidInputCount', ...
            ['Expected root, region, and dirres ' ...
            'as three input arguments.']);
    end

    root = local_text_scalar(root,'root');
    region = local_text_scalar(region,'region');
    dirres = local_text_scalar(dirres,'dirres');

    ensureSAGEpath(root);
    region = region_helpers('code',region);

    if ~isfolder(dirres)
        [status,message] = mkdir(dirres);
        if ~status
            error('bootstrap_SAGE:CannotCreateResultsDirectory', ...
                'Cannot create results directory %s: %s', ...
                dirres,message);
        end
    end
end

% ============================================
function value = local_text_scalar(value,name)
% ============================================
% Validate and normalize a character vector or string scalar.

    if ~(ischar(value) ...
            || (isstring(value) ...
            && isscalar(value)))
        error('bootstrap_SAGE:InvalidTextInput', ...
            ['%s must be a character vector ' ...
            'or string scalar.'],name);
    end

    value = char(string(value));
    if isempty(strtrim(value))
        error('bootstrap_SAGE:EmptyTextInput', ...
            '%s cannot be empty.',name);
    end
end
