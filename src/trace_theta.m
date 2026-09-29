function tTheta = trace_theta(tTheta,nTheta,bas,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRACE_THETA Store parameter percentiles for one iteration.
%
%  Compute seven percentiles across training basins for each normalized
%  model parameter and update the current trace slice.
%
% SYNOPSIS:
%   tTheta = trace_theta(tTheta,nTheta,bas,i)
%
% INPUT ARGUMENTS:
%   tTheta          preallocated i_max-by-d-by-7 percentile trace
%   nTheta          d-by-K normalized-parameter matrix
%   bas             basin partition information
%    .K_t            number of training basins in leading columns
%   i               current iteration index
%
% OUTPUT ARGUMENTS:
%   tTheta          trace updated at tTheta(i,:,:)
%
% NOTES:
%   Percentiles are 5, 15, 25, 50, 75, 85, and 95.
%   Nonfinite parameter values are omitted; empty rows become NaN.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Mar. 2026                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(nTheta)
        return
    end

    if ~ismatrix(nTheta)
        error(['      Error: trace_theta: ' ...
            'nTheta must be a 2-D array ' ...
            'of size [d x K].']);
    end

    if bas.K_t < 1
        return
    end

    if size(nTheta,2) < bas.K_t
        error(['      Error: trace_theta: ' ...
            'size(nTheta,2) = %d but ' ...
            'bas.K_t = %d.'],size(nTheta,2), ...
            bas.K_t);
    end

    d = size(nTheta,1);

    if isempty(tTheta)
        error(['      Error: trace_theta: ' ...
            'tTheta is empty. ' ...
            'Initialize tTheta first as [i_max x d x 7].']);
    end

    if ndims(tTheta) ~= 3 ...
            || size(tTheta,2) ~= d ...
            || size(tTheta,3) ~= 7
        error(['      Error: trace_theta: ' ...
            'tTheta must have size ' ...
            '[i_max x d x 7], with d = size(nTheta,1).']);
    end

    p = [5 15 25 50 75 85 95];
    X = nTheta(:,1:bas.K_t);
    has_prctile = exist('prctile','file') ~= 0;

    for j = 1:d
        x = X(j,:);
        x = x(isfinite(x));

        if isempty(x)
            tTheta(i,j,:) = NaN(1,1,7);
        else
            if has_prctile
                pct = prctile(x,p);
            else
                pct = sage_prctile(x,p);
            end
            tTheta(i,j,:) = reshape(pct,1,1,7);
        end
    end

end
