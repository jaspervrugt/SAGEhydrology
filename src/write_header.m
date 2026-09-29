function write_header(mdl,part)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WRITE_HEADER Print the selected model and training banner.
%
%  Writes the command-window heading for SITE or SAGE training.
%
% SYNOPSIS:
%   write_header(mdl)
%   write_header(mdl,part)
%
% INPUT ARGUMENTS:
%   mdl             model-selection settings
%    .model          hydrologic model identifier
%    .variant        optional model variant
%   part            optional 'site' or 'sage'; default 'sage'
%
% OUTPUT ARGUMENTS:
%   none            text is written to the command window
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(part)
        part = 'sage';
    end
    
    model = mdl.model;
    mname = sage_model_name(model);
    if isfield(mdl,'variant') ...
            && strcmpi(string(mdl.variant),'gchm_ode')
        mname = 'gchm_ode';
    end
    
    w = 43;   % inside width (characters)
    
    if strcmpi(part,'sage')
        msg = sprintf('Hydrologic model: %s', ...
            mname);
    elseif strcmpi(part,'site')
        msg = sprintf('Individual training: %s', ...
            mname);
    else
        error(['      Error: write_header: ' ...
            'unknown entry of variable part']);
    end
    
    if numel(msg) > w
        msg = msg(1:w);
    end
    
    padL = floor((w - numel(msg))/2);
    padR = w - numel(msg) - padL;
    line = [repmat(' ',1,padL) msg repmat(' ',1,padR)];
    
    % -----------------------
    % SAGE header
    % -----------------------
    disp('');
    disp('           ┌─────────────────────────────────────────────┐ ');
    disp('           |                                             | ');
    disp('           ├─────────────── SAGEhydrology ───────────────┤ ');
    disp('           |                                             | ');
    disp('           |    Sensitivity-Aware Gradient Estimation    | ');
    disp('           |                                             | ');
    disp('           ├─────────────────────────────────────────────┤ ');
    fprintf('           ├─%s─┤ \n', line);
    disp('           ├─────────────────────────────────────────────┤ ');
    disp('           |                                             | ');
    disp('           |  © 2026 Jasper A. Vrugt                     | ');
    disp('           |    All rights reserved                      | ');
    disp('           |                                             | ');
    disp('           | ✉ jasper@uci.edu                           | ');
    disp('           └─────────────────────────────────────────────┘ ');
    disp('');
    fprintf('\n');

end
