function name = sage_model_name(model,rootFolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_MODEL_NAME Return a canonical SAGE model name.
%
%  Map supported numeric model identifiers to names. Unknown codes are
%  returned as "model N".
%
% SYNOPSIS:
%   name = sage_model_name(model)
%
% INPUT ARGUMENTS:
%   model           scalar or array of numeric model identifiers
%
% OUTPUT ARGUMENTS:
%   name            char for scalar input; string array otherwise
%
% NOTES:
%   Codes: 1 HYMOD; 2 HMODEL; 3 SAC-SMA; 4 Xinanjiang; 5 GR4J; 6 HBV; 7
%   CFE-NWM; 11 GCHM; 12 private MN(DS)_SALO(5); 99 user model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2
        rootFolder = '';
    end
    ids = [1 2 3 4 5 6 7 11 99 12];
    names = ["hymod","hmodel","sacsma","Xinanjiang", ...
        "gr4jA","hbv","cfe_nwm","gchm","user_model","mcp_salo"];

    model = double(model);
    name = strings(size(model));
    for i = 1:numel(model)
        j = find(ids == model(i),1,'first');
        if isempty(j)
            name(i) = "model " + string(model(i));
        else
            name(i) = names(j);
        end
        if model(i) == 99
            authoredName = local_ai_model_name(rootFolder);
            if strlength(authoredName) > 0
                name(i) = authoredName;
            end
        end
    end

    if isscalar(name)
        name = char(name);
    end
end

function name = local_ai_model_name(rootFolder)
    name = "";
    rootFolder = char(string(rootFolder));
    here = fileparts(mfilename('fullpath'));
    sageRoot = fileparts(fileparts(here));
    softwareRoot = fileparts(sageRoot);

    candidateDirs = {};
    if ~isempty(rootFolder)
        candidateDirs{end+1} = fullfile(rootFolder, ...
            'AI_assisted_model','active');
        candidateDirs{end+1} = fullfile(fileparts(rootFolder), ...
            'AI_assisted_model','active');
        candidateDirs{end+1} = fullfile(rootFolder,'user_model');
        candidateDirs{end+1} = fullfile(fileparts(rootFolder), ...
            'user_model');
    end
    candidateDirs{end+1} = fullfile(softwareRoot, ...
        'AI_assisted_model','active');
    candidateDirs{end+1} = fullfile(softwareRoot,'user_model');

    try
        mexPath = which(['crr_user_model.' mexext]);
        if ~isempty(mexPath)
            candidateDirs{end+1} = fileparts(mexPath);
        end
    catch
    end
    if isdeployed
        candidateDirs{end+1} = fullfile(ctfroot,'user_model');
        candidateDirs{end+1} = fullfile(ctfroot, ...
            'AI_assisted_model','active');
    end

    candidateDirs = unique(candidateDirs,'stable');
    if isdeployed
        candidateDirs=candidateDirs(~contains(string(candidateDirs), ...
            [filesep 'AI_assisted_model' filesep]));
    end
    for k = 1:numel(candidateDirs)
        manifestFile = fullfile(candidateDirs{k}, ...
            'ai_model_manifest.json');
        if ~isfile(manifestFile)
            continue
        end
        try
            manifest = jsondecode(fileread(manifestFile));
            candidate = string(manifest.model_name);
        catch
            continue
        end
        if isscalar(candidate) && ~isempty(regexp(char(candidate), ...
                '^[A-Za-z][A-Za-z0-9_]*$', 'once'))
            name = candidate;
            return
        end
    end
end
