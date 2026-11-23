function A = apply_attribute_preprocessing(raw,state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%APPLY_ATTRIBUTE_PREPROCESSING Reuse fitted transforms on new basins.
%  RAW is basins-by-attributes; STATE is from PREPROCESS_ATTRIBUTES.
%  Copyright (c) 2026 Jasper A. Vrugt, UCI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = double(raw);
    if size(x,2) ~= numel(state.columns)
        error('attributes:Shape', ...
            'Saved attribute count differs.');
    end
    if strcmp(state.missing,'reject') && ...
        any(~isfinite(x),'all')
        error('attributes:Missing', ...
            'Saved preprocessing rejects missing.');
    end
    A = zeros(size(x));
    for j = 1:size(x,2)
        v = x(:,j); v(~isfinite(v)) = state.median(j);
        if iscell(state.columns)
            c = state.columns{j};
        else
            c = state.columns(j);
        end
        z = (attribute_nonlinear(v,c)-c.center)/c.scale;
        if c.method == 3
            z = 2*z-1;
        end
        A(:,j) = z;
    end
    if any(~isfinite(A),'all')
        error('attributes:Nonfinite', ...
            'Transform produced nonfinite values.');
    end
end
