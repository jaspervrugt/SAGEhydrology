function validate_initial_storage(mdl,th)
% VALIDATE_INITIAL_STORAGE Check defined capacities against actual basin
% parameters.
% No artificial capacity is imposed on unbounded routing or groundwater
% stores.
    y=mdl.y0(:); limit=inf(size(y));
    switch mdl.model
        case 1, limit(2)=th(1);
        case 2, limit(2:3)=th(1:2);
        case 3, limit(2:6)=th([2 1 5 3 4]);
        case 4, limit(2:3)=[th(5)*th(8);(1-th(5))*th(8)];
        case 5, limit(2:3)=th([1 3]);
        case 6, limit(2)=th(1);
        case 7, limit(2)=th(1);
        case 99
            if isfield(mdl,'model_template') && ...
                strcmp( mdl.model_template,'gr4jb')
                limit(2:3)=th([1 3]);
            end
    end
    bad=find(y>limit+1e-12*max(1,limit),1);
    if ~isempty(bad)
        error('validate_initial_storage:Capacity', ...
            ['Initial state %d (%g) exceeds this basin''s storage ' ...
                'capacity (%g).'],bad,y(bad),limit(bad));
    end
end
