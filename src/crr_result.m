function out = crr_result(q,J,Jth,Z,mdl,request,named)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRR_RESULT Package requested MATLAB model outputs.
%
%  Collects selected observables, Jacobians, and state trajectories from a
%  model run.
%
% SYNOPSIS:
%   out = crr_result(q,J,Jth,Z,mdl,request,named)
%
% INPUT ARGUMENTS:
%   q               simulated discharge vector
%   J               discharge Jacobian in the selected parameter space
%   Jth             physical-to-selected parameter derivative vector
%   Z               augmented state and sensitivity history
%   mdl             model and output-time settings
%    .y0             initial state vector
%    .idx            output indices for complete trajectories
%    .state_name     optional model-state names
%    .swe_ind        optional SWE state indices
%   request         selection of named observations, Jacobians, and states
%   named           optional named results that take precedence over Z
%
% OUTPUT ARGUMENTS:
%   out             structure with only requested model results
%    .failed         true if simulated discharge is nonfinite
%    .obs            requested Q, SWE, and SM trajectories
%    .jac            requested Q, SWE, and SM Jacobians
%    .states         requested model-state trajectories
%
% NOTES:
%   J is already in the selected parameter space. SWE sensitivities
%   extracted from Z are scaled columnwise by Jth.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Sept. 2026                                %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    request = crr_request(request);
    if nargin < 7 ...
            || isempty(named)
        named = struct();
    end
    failed = any(~isfinite(q(:)));
    if isfield(named,'obs')
        fields = fieldnames(named.obs);
        for k = 1:numel(fields)
            failed = failed || any(~isfinite( ...
                named.obs.(fields{k})(:)));
        end
    end
    out = struct('failed',failed);
    
    if ~isempty(request.obs)
        out.obs = struct();
        if any(request.obs == "Q")
            out.obs.Q = q; 
        end
    end

    if ~isempty(request.jac)
        out.jac = struct();
        if any(request.jac == "Q")
            out.jac.Q = J; 
        end
    end
    
    needSwe = any(request.obs == "SWE") ...
        || any(request.jac == "SWE");
    
    if needSwe
        swe_ind = 1;
        if isfield(mdl,'swe_ind') && ~isempty(mdl.swe_ind)
            swe_ind = mdl.swe_ind(:)';
        elseif mdl.model == 11
            swe_ind = [1 3];
        end
        directObs = isfield(named,'obs') ...
            && isfield(named.obs,'SWE');
        directJac = isfield(named,'jac') ...
            && isfield(named.jac,'SWE');
        rows = mdl.idx(1)+1:mdl.idx(2);
        
        if any(request.obs == "SWE")
            if ~isfield(out,'obs')
                out.obs = struct(); 
            end
            if directObs
                out.obs.SWE = named.obs.SWE;
            else
                out.obs.SWE = sum(Z(rows,swe_ind),2);
            end
        end
        
        if any(request.jac == "SWE")
            if directJac
                out.jac.SWE = named.jac.SWE;
            else
                d = numel(Jth);
                if mdl.model == 6
                    m = numel(mdl.y0);
                else
                    m = size(Z,2)/(d+1);
                end
                dz = size(Z,2)/m - 1;
            
                if dz ~= fix(dz)
                    error('crr_result:History', ...
                        'Invalid augmented history.');
                end
                Jswe = zeros(numel(rows),d);
            
                for j = 1:min(d,dz)
                    cols = j*m+swe_ind;
                    Jswe(:,j) = sum(Z(rows,cols),2).*Jth(j);
                end
            
            if ~isfield(out,'jac')
                out.jac = struct(); 
            end
                out.jac.SWE = Jswe;
            end
        end
    end

    needSm = any(request.obs == "SM") ...
        || any(request.jac == "SM");
    if needSm
        sm_ind = local_sm_indices(mdl);
        if isempty(sm_ind)
            error('crr_result:UnavailableObservable', ...
                ['SM is not a supported training observable ' ...
                'for model %d.'],mdl.model);
        end
        rows = mdl.idx(1)+1:mdl.idx(2);
        directObs = isfield(named,'obs') ...
            && isfield(named.obs,'SM');
        directJac = isfield(named,'jac') ...
            && isfield(named.jac,'SM');
        if any(request.obs == "SM")
            if ~isfield(out,'obs')
                out.obs = struct(); 
            end
            if directObs
                out.obs.SM = named.obs.SM;
            else
                out.obs.SM = sum(Z(rows,sm_ind),2);
            end
        end
        if any(request.jac == "SM")
            if ~isfield(out,'jac')
                out.jac = struct(); 
            end
            if directJac
                out.jac.SM = named.jac.SM;
            else
                d = numel(Jth);
                if mdl.model == 6
                    m = numel(mdl.y0);
                else
                    m = size(Z,2)/(d+1);
                end
                dz = size(Z,2)/m - 1;
                if dz ~= fix(dz)
                    error('crr_result:History', ...
                        'Invalid augmented history.');
                end
                Jsm = zeros(numel(rows),d);
                for j = 1:min(d,dz)
                    cols = j*m + sm_ind;
                    Jsm(:,j) = sum(Z(rows,cols),2).*Jth(j);
                end
                out.jac.SM = Jsm;
            end
        end
    end
    
    if ~isempty(request.states)
        out.states = model_states(mdl,Z,request.states,numel(q)); 
    end
end

function ind = local_sm_indices(mdl)
%LOCAL_SM_INDICES Return configured soil-moisture state indices.

    if isfield(mdl,'sm_ind') ...
            && ~isempty(mdl.sm_ind)
        ind = mdl.sm_ind(:)';
        return
    end
    switch mdl.model
        case 2
            ind = 3;
        case 3
            ind = 2:6;
        case 4
            ind = 2:3;
        case {5,6,7}
            ind = 2;
        otherwise
            ind = [];
    end
end
