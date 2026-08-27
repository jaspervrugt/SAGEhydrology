function args = jkge_args(loss,block)
%JKGE_ARGS Build JKGE varargin cell array for JKGE benchmark methods
%
% SYNOPSIS: args = jkge_args(loss,block)
%   loss        structure with JKGE settings
%     .method   scalar JKGE benchmark method
%                1 = moving-average mean
%                2 = section-wise mean
%                3 = long-term mean
%                4 = monthly climatology
%     .n_win    scalar window/block size for methods 1 and 2
%     .meta     structure with auxiliary data for method 4
%                .mo_t   month labels for training record
%                .mo_e   month labels for evaluation record
%                .mo_all month labels for full record
%   block       string specifying record subset
%                't'   = training record
%                'e'   = evaluation record
%                'all' = full record
%   args        OUTPUT: 1xN cell array additional inputs JKGE routines
%               such that:
%                 method 1: args = {1, n_win}
%                 method 2: args = {2, n_win}
%                 method 3: args = {3}
%                 method 4: args = {4, mo_block}
%
% DESCRIPTION:
%   Constructs the variable input argument list required by the JKGE
%   benchmark and efficiency functions. The returned cell array can be
%   expanded using args{:} in calls to jkge and related routines.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Optional M-definition selector for the benchmark-matching term
    %   M = 1: original ratio form from the paper
    %   M = 2: norm-based formulation
    if isfield(loss,'M') ...
            && ~isempty(loss.M)
        Mdef = loss.M;
    else
        Mdef = 2;
    end
    
    if ~(isscalar(Mdef) ...
            && isnumeric(Mdef) ...
            && isfinite(Mdef) ...
            && any(Mdef == [1 2]))
        error(['      Error:jkge_args: ' ...
            'loss.M must be 1 or 2.']);
    end
    
    % switch loss.method
    %     case 1 % moving-average mean
    %         args = {1, loss.n_win};
    % 
    %     case 2 % section-wise mean
    %         args = {2, loss.n_win};
    % 
    %     case 3 % long-term mean
    %         args = {3, []};
    % 
    %     case 4 % monthly climatology
    %         switch lower(char(block))
    %             case 't'
    %                 args = {4, loss.meta.mo_t};
    %             case 'e'
    %                 args = {4, loss.meta.mo_e};
    %             case 'all'
    %                 args = {4, loss.meta.mo_all};
    %             otherwise
    %                 error(['      Error:jkge_args: ' ...
    %                     'unknown block = %s.'],block);
    %         end
    % 
    %     otherwise
    %         error(['      Error:jkge_args: ' ...
    %             'unknown JKGE method = %g.'], ...
    %             loss.method);
    % end
    switch loss.method
        case 1
            args = {1, loss.n_win, Mdef};

        case 2
            args = {2, loss.n_win, Mdef};

        case 3
            args = {3, Mdef};

        case 4
            switch lower(char(block))
                case 't'
                    args = {4, loss.meta.mo_t, Mdef};
                case 'e'
                    args = {4, loss.meta.mo_e, Mdef};
                case 'all'
                    args = {4, loss.meta.mo_all, Mdef};
                otherwise
                    error(['      Error:jkge_args: ' ...
                        'unknown block = %s.'],block);
            end
    end

end