function [swe,Jswe,sm,Jsm] = model_capture_states( ...
    z,row,m,d,swe_ind,sm_ind,needSwe,needJswe, ...
    needSm,needJsm,swe,Jswe,sm,Jsm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL_CAPTURE_STATES Store named states at one output time.
%
%  Stores requested snow and soil-moisture state outputs.
%
% SYNOPSIS:
%   [swe,Jswe,sm,Jsm] = model_capture_states(...)
%
% INPUT ARGUMENTS:
%   z               current augmented model state
%   row             output row index
%   m               number of physical model states
%   d               number of model parameters
%   swe_ind         state indices contributing to SWE
%   sm_ind          state indices contributing to SM
%   needSwe         whether SWE was requested
%   needJswe        whether SWE Jacobian was requested
%   needSm          whether SM was requested
%   needJsm         whether SM Jacobian was requested
%   swe             preallocated SWE output
%   Jswe            preallocated SWE Jacobian
%   sm              preallocated SM output
%   Jsm             preallocated SM Jacobian
%
% OUTPUT ARGUMENTS:
%   swe             SWE output with current row populated
%   Jswe            SWE Jacobian with current row populated
%   sm              SM output with current row populated
%   Jsm             SM Jacobian with current row populated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~(needSwe ...
            || needJswe ...
            || needSm || ...
            needJsm)
        return
    end
    [swe_i,Jswe_i,sm_i,Jsm_i] = ...
        model_state_output(z,m,d,swe_ind,sm_ind);
    if needSwe
        swe(row) = swe_i;
    end
    if needJswe
        Jswe(row,1:d) = Jswe_i;
    end
    if needSm
        sm(row) = sm_i;
    end
    if needJsm
        Jsm(row,1:d) = Jsm_i;
    end
end
