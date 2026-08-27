function state = reader_progress(action,arg,k)
%READER_PROGRESS One-line console progress for regional data readers.
%
%   state = reader_progress('start',label,K)
%   state = reader_progress('update',state,k)
%           reader_progress('finish',state,K)
%
% The transient [k/K, percent done] suffix is rewritten in place. FINISH
% always terminates the console line, so output from the following reader
% cannot be appended to an unfinished progress message.

    action = lower(char(string(action)));

    switch action
        case 'start'
            label = char(string(arg));
            K = max(0,round(double(k)));
            state = struct('label',label,'K',K,'suffix','', ...
                'last',-1);
            fprintf('%s',label);
            state = local_update(state,0);

        case 'update'
            state = arg;
            state = local_update(state,k);

        case 'finish'
            state = arg;
            state = local_update(state,k);
            fprintf(repmat('\b',1,numel(state.suffix)));
            fprintf(' ... Done\n');
            state.suffix = '';

        otherwise
            error('reader_progress:unknownAction', ...
                'Unknown action: %s.',action);
    end
end

function state = local_update(state,k)
    k = max(0,min(state.K,round(double(k))));
    if k == state.last
        return
    end

    suffix = sprintf(' [%d/%d, %5.1f%% done]', ...
        k,state.K,100*k/max(1,state.K));
    if ~isempty(state.suffix)
        fprintf(repmat('\b',1,numel(state.suffix)));
    end
    fprintf('%s',suffix);
    state.suffix = suffix;
    state.last = k;
    drawnow limitrate
end
