function write_header(mdl,part)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WRITE_HEADER Write SAGE header
%
% SYNOPSIS:
%   flag = write_header(mdl,part)
%
% INPUT:
%   mdl         structure with model selection/settings
%    .model      choice of model
%                 1 hymod
%                 2 hmodel
%                 3 sacsma
%                 4 xinanjiang
%                 5 gr4jA
%                 6 hbv
%                 7 gr4jB
%    .mcode      numerical implementation
%                 1 Runge Kutta MATLAB
%                 2 ode45 MATLAB
%                 3 Euler MATLAB
%                 4 Runge Kutta C++/MEX
%    .mode       evaluation mode
%      'seq'      sequential
%      'par'      parallel
%    .names      list/cell array of model names
%
%   part        OPTIONAL: training approach
%    'site'      single-site training
%    'sage'      SAGE training
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(part)
        part = 'sage';
    end
    
    model = mdl.model;
    model_names = mdl.names;
    
    if iscell(model_names)
        mname = char(model_names{model});
    else
        mname = char(model_names(model));
    end
    
    w = 43;   % inside width (characters)
    
    if strcmpi(part,'sage')
        msg = sprintf('Hydrologic model: %s', ...
            mname);
    elseif strcmpi(part,'site')
        msg = sprintf('Individual training: %s', ...
            mname);
    else
        error(['      Error: compile_model: ' ...
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