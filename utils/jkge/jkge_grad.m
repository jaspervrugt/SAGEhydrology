function [JKGE,delta,m_q,M_qy,V_qy,C_qy] = ...
    jkge_grad(y,q,m_y,idx,method,n_win_or_mo,cache,Mdef)
%JKGE_GRAD Compute JKGE and its loss-sensitivity vector
%
% SYNOPSIS:
%  [JKGE,delta,m_q,M_qy,V_qy,C_qy] = ...
%      jkge_grad(y,q,m_y,idx,method,n_win_or_mo,cache,Mdef)
%
% INPUT:
%   y             nx1 observed discharge vector on full model time axis
%   q             nx1 simulated discharge vector on full model time axis
%   m_y           nx1 observed benchmark vector on full model time axis
%   idx           kx1 vector or nx1 logical mask with entries used in the
%                 JKGE objective
%   method        scalar JKGE benchmark method
%                   1 = moving-average mean
%                   2 = section-wise mean
%                   3 = long-term mean
%                   4 = monthly climatology
%   n_win_or_mo   scalar benchmark window length for methods 1 and 2,
%                 or nx1 month labels in {1,...,12} for method 4
%   cache         optional precomputed benchmark/cache information
%   Mdef          optional M-component definition
%                   1 = paper definition, mean squared pointwise ratio err
%                   2 = revised norm-based benchmark mismatch
%                 If empty or missing, Mdef = 2 is used.
%
% OUTPUT:
%   JKGE     scalar JKGE efficiency
%   delta    nx1 loss-sensitivity vector dL/dq on full model time axis
%   m_q      nx1 simulated benchmark vector on full model time axis
%   M_qy     scalar benchmark-matching component
%   V_qy     scalar variability component
%   C_qy     scalar correlation component
%
% NOTES:
%   - This routine combines jkge and delta_n (loss_fnc = 7) to avoid
%     redundant computations of anomalies, variances, and correlation.
%   - The simulated benchmark m_q is computed once on the full record
%     and reused for both objective evaluation and gradient calculation.
%   - The loss-sensitivity vector delta is returned on the full model
%     time axis; reduction to training indices must be done externally.
%   - Mdef controls only the M component of JKGE; the variability and
%     correlation components are unchanged.
%   - Mdef = 1 uses the paper formulation based on pointwise benchmark
%     ratios m_q ./ m_y, excluding near-zero m_y values.
%   - Mdef = 2 uses the revised norm-based benchmark mismatch.
%
%   Reference:
%     Jawad, M., Gupta, H.V., Wang, Y.H., Farmani, M.A., Behrangi, A., and
%     Niu, G.Y. (2026), Improving model performance by adapting the
%     KGE metric to account for system non-stationarity,
%     https://arxiv.org/pdf/2604.03906.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method < 3
    n_win = n_win_or_mo;
    mo = [];
elseif method == 4
    n_win = [];
    mo = n_win_or_mo;
else
    n_win = [];
    mo = [];
end

if nargin < 7
    cache = [];
end
if nargin < 8 ...
        || isempty(Mdef)
    Mdef = 2;
end
if ~(isscalar(Mdef) ...
        && isnumeric(Mdef) ...
        && isfinite(Mdef) ...
        && any(Mdef == [1 2]))
    error(['      Error:jkge_grad: ' ...
        'M must be 1 or 2.']);
end
Mdef = double(Mdef);

y = double(y(:));
q = double(q(:));
m_y = double(m_y(:));

n_full = numel(q);
delta = zeros(n_full,1);

[JKGE,M_qy,V_qy,C_qy] = deal(NaN);

if numel(y) ~= n_full ...
        || numel(m_y) ~= n_full
    error(['      Error:jkge_grad: ' ...
        'y, q, and m_y must have ' ...
        'the same length.']);
end

% Simulated benchmark, computed once
m_q = jkge_benchmark(q,method,n_win,mo,cache);
m_q = double(m_q(:));

% Valid objective support
idx = idx(:);
idx = idx(isfinite(idx) ...
    & idx >= 1 ...
    & idx <= n_full ...
    & mod(idx,1) == 0);

good = false(n_full,1);
good(idx) = true;

good = good ...
    & isfinite(y) ...
    & isfinite(q) ...
    & isfinite(m_y) ...
    & isfinite(m_q);

n = nnz(good);

if n < 2
    return
end

y_red = y(good);
q_red = q(good);
m_y_red = m_y(good);
m_q_red = m_q(good);

a_y = y_red - m_y_red;
a_q = q_red - m_q_red;

s_y = sqrt(mean(a_y.^2));
s_q = sqrt(mean(a_q.^2));

if ~(isfinite(s_y) ...
        && s_y > 0 ...
        && isfinite(s_q) ...
        && s_q > 0)
    return
end

z_y = a_y / s_y;
z_q = a_q / s_q;

rho_qy = mean(z_y .* z_q);
rho_qy = max(-1,min(1,rho_qy));

alpha_qy = s_q / s_y;

switch Mdef
    case 1
        tol_m = max(1e-8,1e-6 * ...
            mean(abs(y_red),'omitnan'));
        goodM_red = abs(m_y_red) > tol_m;

        if ~any(goodM_red)
            return
        end

        ratio_qy = m_q_red(goodM_red) ./ ...
            m_y_red(goodM_red);
        errM = 1 - ratio_qy;
        M_qy = mean(errM.^2);

    case 2
        goodM_red = true(size(m_y_red));
        ybar = mean(y_red,'omitnan');
        denM = norm(m_y_red - ybar);

        if ~(isfinite(denM) ...
                && denM > 0)
            return
        end

        M_qy = (norm(m_q_red - m_y_red) / ...
            denM)^2;
end

V_qy = (1 - alpha_qy)^2;
C_qy = (1 - rho_qy)^2;

L_jkge = sqrt(M_qy + V_qy + C_qy);

if ~(isfinite(L_jkge) && L_jkge > 0)
    return
end

JKGE = 1 - L_jkge;

% M term
dM_dm_q_red = zeros(n,1);
dM_dq_red = zeros(n,1);

switch Mdef
    case 1
        jj = find(goodM_red);
        ngM = numel(jj);
        dM_dm_q_red(jj) = -2 * (1 - ...
            m_q_red(jj) ./ m_y_red(jj)) ./ ...
            (ngM * m_y_red(jj));

    case 2
        dM_dm_q_red = 2 * (m_q_red ...
            - m_y_red) / denM^2;
end

% V term
dV_da_q = -2 * (1 - alpha_qy) * ...
    a_q / (n * s_q * s_y);

% C term
S = sum(z_y .* a_q);
drho_da_q = z_y / (n * s_q) ...
    - S * a_q / (n^2 * s_q^3);

dC_da_q = -2 * (1 ...
    - rho_qy) * drho_da_q;

% a_q = q_red - m_q_red
dPhi_dq_red = dM_dq_red ...
    + dV_da_q + dC_da_q;

dPhi_dm_q_red = dM_dm_q_red ...
    - dV_da_q - dC_da_q;

% L = sqrt(Phi)
dL_dq_red = dPhi_dq_red / ...
    (2 * L_jkge);
dL_dm_q_red = dPhi_dm_q_red / ...
    (2 * L_jkge);

% Direct path
delta(good) = delta(good) ...
    + dL_dq_red;

% Benchmark path through B'
v_full = zeros(n_full,1);
v_full(good) = dL_dm_q_red;

delta = delta + jkge_transpose( ...
    v_full,q,method,n_win,mo,cache);

end