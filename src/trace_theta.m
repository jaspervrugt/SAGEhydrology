function tTheta = trace_theta(tTheta,nTheta,bas,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRACE_THETA Store normalized-parameter percentile traces for iteration i
%
% SYNOPSIS:
%  tTheta = trace_theta(tTheta,nTheta,bas,i)
%   tTheta    percentile traces [i_max x d x 7]
%   nTheta    normalized parameter matrix [d x K]
%   bas       structure with basin information
%    .K_t     number of training basins
%   i         current iteration
%   tTheta    OUTPUT: updated percentile traces with
%             tTheta(i,j,:) = [5 15 25 50 75 85 95] percentiles
%             of normalized parameter j across training basins
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

for j = 1:d
    x = X(j,:);
    x = x(isfinite(x));

    if isempty(x)
        tTheta(i,j,:) = NaN(1,1,7);
    else
        tTheta(i,j,:) = reshape(prctile(x,p),1,1,7);
    end
end

end