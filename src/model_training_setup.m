function [net,alg,status]=model_training_setup(mdl,net,alg)
%MODEL_TRAINING_SETUP Apply installed model-specific source training settings.
% Built-in models retain their settings. Private hooks are resolved dynamically.
    status='';
    if isfield(net,'initial_output_model') && net.initial_output_model==12
        old=intersect(fieldnames(net),{'initial_output','initial_output_weight_scale','initial_output_model'});
        net=rmfield(net,old);
    end
    if isequal(double(mdl.model),12)
        configure=str2func('mcp_salo_training_setup');
        [net,alg,status]=configure(mdl,net,alg);
    end
end
