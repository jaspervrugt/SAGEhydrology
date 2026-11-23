function loss = configure_multi_loss( ...
        loss,observed,weight,transform,normalization)
%CONFIGURE_MULTI_LOSS Configure response-specific objectives for SAGE.
%  Untransformed responses retain loss.fnc. Zscore responses use the
%  standardized MSE and analytic derivatives in observation_loss and
%  zscore_loss. Model equations and sensitivities are unchanged.
%
%  loss = configure_multi_loss(loss,observed,weight,transform,normalization)
%
% Copyright (c) 2026 Jasper A. Vrugt, University of California, Irvine.
    loss.observed = observed;
    loss.weight = weight;
    loss.transform = transform;
    loss.normalization = normalization;
    for name = string(observed(:)).'
        observation_transform(loss,name);
    end
end
