function loss = fdc_loss_cached(q_n,fdc,formulation)
%FDC_LOSS_CACHED Selected FDC divergence using cached observations.

    if nargin < 3 || isempty(formulation)
        formulation = 1;
    end
    [D_fdc,D_p,D_logp] = fdc_metrics_cached(q_n,fdc);
    switch local_formulation(formulation)
        case 1
            loss = D_fdc;
        case 2
            loss = D_p;
        case 3
            loss = D_logp;
    end
end

function formulation = local_formulation(value)
    if isnumeric(value) && isscalar(value) && isfinite(value)
        formulation = double(value);
    else
        key = lower(regexprep(char(string(value)),'[^a-z0-9]',''));
        switch key
            case {'1','a','fdc','dfdc','physical','physicalcdf'}
                formulation = 1;
            case {'2','b','p','dp','quantile'}
                formulation = 2;
            case {'3','c','logp','dlogp','logquantile'}
                formulation = 3;
            otherwise
                formulation = NaN;
        end
    end
    if ~ismember(formulation,1:3)
        error('fdc_loss_cached:BadFormulation', ...
            'FDC formulation must be 1 (d_fdc), 2 (d_p), or 3 (d_logp).');
    end
end
