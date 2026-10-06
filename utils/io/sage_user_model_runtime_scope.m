function cleanup = sage_user_model_runtime_scope(mdl)
% Ensure each Runtime worker uses the selected external precompiled MEX.
% Compiled calls bind to the packaged slot; no path/folder change is needed.
    cleanup=[];
    if ~isdeployed, return, end
    sage_prepare_deployed_user_mex(mdl);
end
