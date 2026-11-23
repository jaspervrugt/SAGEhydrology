function [A,bas] = prepare_attribute_inputs(A,bas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE_ATTRIBUTE_INPUTS Fit selected active basins after quality
% screening.
%  A uses MATLAB's attributes-by-basins orientation. Older configurations
%  without attribute_preprocessing retain the historical reader behavior.
%  Copyright (c) 2026 Jasper A. Vrugt, UCI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(bas,'attribute_preprocessing'), return; end
    options = bas.attribute_preprocessing;
    if strcmp(options.mode,'regional'), return; end
    [x,state] = preprocess_attributes(A.',1:bas.K_t,options);
    A = x.';
    state.attribute_ids = bas.id_attr;
    bas.attribute_state = state;
end
