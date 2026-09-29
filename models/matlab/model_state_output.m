function [swe,Jswe,sm,Jsm] = ...
    model_state_output(z,m,d,swe_ind,sm_ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL_STATE_OUTPUT Extract named state values and sensitivities.
%
%  Sums selected states and their parameter-sensitivity blocks.
%
% SYNOPSIS:
%   [swe,Jswe,sm,Jsm] = model_state_output(z,m,d,swe_ind,sm_ind)
%
% INPUT ARGUMENTS:
%   z               current augmented state vector
%   m               number of physical model states
%   d               number of model parameters
%   swe_ind         indices of SWE component states
%   sm_ind          indices of SM component states
%
% OUTPUT ARGUMENTS:
%   swe             sum of SWE component states
%   Jswe            1-by-d SWE sensitivities
%   sm              sum of SM component states
%   Jsm             1-by-d SM sensitivities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    swe = sum(z(swe_ind));
    sm = sum(z(sm_ind));
    Jswe = zeros(1,d);
    Jsm = zeros(1,d);
    for j = 1:d
        offset = j*m;
        Jswe(j) = sum(z(offset+swe_ind));
        Jsm(j) = sum(z(offset+sm_ind));
    end
end
