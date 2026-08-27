function infoFile = make_user_model_info(outDir)
%MAKE_USER_MODEL_INFO Create user_model_info.mat for SAGE user_model.
%
% SYNOPSIS:
%   infoFile = make_user_model_info()
%   infoFile = make_user_model_info(outDir)
%
% OUTPUT:
%   user_model_info.mat with fields:
%       par_info
%       m
%       y0
%       pspace
%
% NOTES:
%   par_info columns are:
%       {index, name, description, daily_unit, daily_min, daily_max}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, May 2026                                  %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 ...
        || isempty(outDir)
    outDir = fileparts( ...
        mfilename('fullpath'));
end

if ~isfolder(outDir)
    mkdir(outDir);
end

%  =======================================================
%  USER-EDIT BLOCK 1:
%  Replace with parameters of your hydrologic model
%  Make sure you use latex language for parameter names
%  "f_{\rm p}" means f is italic, subscript "p" is upright
%  for subscripts that are not variables but abbreviations
%  ======================================================= 
par_info = {
    1, 'x_{1}',      'production store capacity',       ...
    'mm',       25, 4000; % 1,4000 --> avoid enormous sensitivities
    2, 'x_{2}',      'groundwater exchange flux',       ...
    'mm/d',    -20,   15;
    3, 'x_{3}',      'routing store capacity',          ...
    'mm',       10, 2000; % 1,2000 --> avoid enormous sensitivities
    4, 'x_{4}',      'routing time base',               ...
    'd',       0.5,   10; % 0.1,10 --> avoid enormous sensitivities
    5, 'x_{5}',      'routing flow partition factor',   ...
    '-',      1e-3,  1.0;
    6, 'f_{\rm p}',  'pan evaporation coefficient',     ...
    '-',       0.3,  1.7;
    7, 'T_{\rm tr}', 'snow/rain temperature threshold', ...
    '°C',       -3,    3;
    8, 'f_{\rm dd}', 'degree-day melt factor',          ...
    'mm/d/°C', 0.1,   10};

%Small \(x_3\): routing discharge scales approximately as \(R^5/x_3^4\). Values near zero create enormous discharge sensitivities.
%Small \(x_1\): several production-store equations contain \(S/x_1\), \(P/x_1\), and inverse powers of \(x_1\), producing steep gradients.
%\(x_4\leq0.5\) day: dynamic routing gives \(L=1\). After normalization \(U_1=1\), so \(x_4\) no longer affects discharge and its derivative is zero.

% Number of physical states used for initial state vector.
% For this GR4J-B-style user model: 1: Swe, 2: production store, 
% 3: routing store, 4: UH memories, 5: Q bookkeeping.
m = 5;

%  =======================================================
%  END USER-EDIT BLOCK 1
%  ======================================================= 

% Initial physical states
y0 = 1e-5 * ones(m,1);

% Parameter space:
%   0 = hydrologic space
%   1 = normalized [0,1] space
%   2 = transformed/unconstrained space
pspace = 1;

infoFile = fullfile(outDir, ...
    'user_model_info.mat');
save(infoFile,'par_info', ...
    'm','y0','pspace');
fprintf(['Wrote user model ' ...
    'metadata:\n  %s\n'], ...
    infoFile);

end