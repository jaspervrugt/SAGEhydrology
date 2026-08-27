function [JKGE,m_q,M_qy,V_qy,C_qy,rho_qy,a_y,a_q, ...
    s_y,s_q,good,goodM,eM,alpha_qy,z_y,z_q,m_y_red] = ...
    jkge(y,q,m_y,varargin)
%JKGE Jawad-Kling-Gupta efficiency
%
% SYNOPSIS:
%  [JKGE,m_q,M_qy,V_qy,C_qy,rho_qy,a_y,a_q, ...
%   s_y,s_q,good,goodM,eM,alpha_qy,z_y,z_q,m_y_red] = ...
%      jkge(y,q,m_y,varargin)
%
%   y        nx1 observed discharge vector on full model time axis
%   q        nx1 simulated discharge vector on full model time axis
%   m_y      nx1 observed benchmark vector on full model time axis
%   varargin additional inputs defining optional valid indices and the
%            JKGE benchmark
%
% CALLING CONVENTION:
%   Full-record calculation:
%     method 1: jkge(y,q,m_y,1,n_win,M)
%     method 2: jkge(y,q,m_y,2,n_win,M)
%     method 3: jkge(y,q,m_y,3,M)
%     method 4: jkge(y,q,m_y,4,mo,M)
%
%   Indexed calculation, recommended for hourly data with missing discharge
%     method 1: jkge(y,q,m_y,idx,1,n_win,M)
%     method 2: jkge(y,q,m_y,idx,2,n_win,M)
%     method 3: jkge(y,q,m_y,idx,3,M)
%     method 4: jkge(y,q,m_y,idx,4,mo,M)
%
%   idx      optional kx1 vector or nx1 logical mask with entries used in
%            the JKGE calculation. The simulated benchmark m_q is always
%            computed first on the full q record and only then reduced to
%            idx and finite entries.
%   method   scalar JKGE benchmark method
%              1 = moving-average mean
%              2 = section-wise mean
%              3 = long-term mean
%              4 = monthly climatology
%   n_win    scalar benchmark scale in samples required for methods 1 and 2
%   mo       nx1 month labels in {1,...,12} required for method 4
%
% OUTPUT:
%   JKGE     scalar JKGE efficiency
%   m_q      nx1 simulated benchmark vector on full model time axis
%   M_qy     scalar benchmark-matching component
%   V_qy     scalar variability component
%   C_qy     scalar correlation component
%   rho_qy   scalar anomaly correlation observed/simulated anomalies
%   a_y      kx1 observed anomaly vector after idx/finite reduction
%   a_q      kx1 simulated anomaly vector after idx/finite reduction
%   s_y      scalar root-mean-square observed anomaly
%   s_q      scalar root-mean-square simulated anomaly
%   good     nx1 logical index of entries used in JKGE calculation
%   goodM    nx1 logical index of entries used in benchmark mismatch term
%   eM       scalar mean benchmark error, used for def = 1; NaN for def = 2
%   alpha_qy scalar variability ratio s_q / s_y
%   z_y      kx1 standardized observed anomalies
%   z_q      kx1 standardized simulated anomalies
%   m_y_red  kx1 reduced observed benchmark vector
%
%   Reference:
%     Jawad, M., Gupta, H.V., Wang, Y.H., Farmani, M.A., Behrangi, A., and
%     Niu, G.Y. (2026), Improving model performance by adapting the
%     KGE metric to account for system non-stationarity,
%     https://arxiv.org/pdf/2604.03906.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, April 2026                                %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize scalars
[JKGE,M_qy,V_qy,C_qy,rho_qy, ...
    s_y,s_q,eM,alpha_qy] = deal(nan);

% initialize vectors
m_q = [];
a_y = nan(0,1);
a_q = nan(0,1);
z_y = nan(0,1);
z_q = nan(0,1);
m_y_red = nan(0,1);
good = false(0,1);
goodM = false(0,1);

if isempty(y) ...
        || isempty(q) ...
        || isempty(m_y)
    return
end

y = double(y(:));
q = double(q(:));
m_y = double(m_y(:));

n = numel(q);
goodM = false(n,1);

if numel(y) ~= n ...
        || numel(m_y) ~= n
    error(['      Error:jkge: ' ...
        'y, q, and m_y must ' ...
        'have the same length.']);
end

% ----------------------------------------
% parse optional idx + benchmark arguments
% ----------------------------------------
[idx,method,n_win,mo,Mdef] = ...
    jkge_parse_inputs(n,varargin{:});

if isempty(idx)
    return
end

% -----------------------------------------------------
% compute simulated benchmark on the full record first
% -----------------------------------------------------
m_q = jkge_benchmark(q,method,n_win,mo);

if numel(m_q) ~= n
    error(['      Error:jkge: ' ...
        'jkge_benchmark returned ' ...
        'm_q with wrong length.']);
end

m_q = double(m_q(:));

% -----------------------------------------------------
% now reduce to requested valid entries + finite values
% -----------------------------------------------------
idx_mask = false(n,1);
idx_mask(idx) = true;

good = idx_mask ...
    & isfinite(y) ...
    & isfinite(q) ...
    & isfinite(m_y) ...
    & isfinite(m_q);

if nnz(good) < 2
    return
end

y_red = y(good);
q_red = q(good);
m_y_red = m_y(good);
m_q_red = m_q(good);

% --------------
% anomaly series
% --------------
a_y = y_red - m_y_red;
a_q = q_red - m_q_red;

s_y = sqrt(mean(a_y.^2));
s_q = sqrt(mean(a_q.^2));

% ---------------------
% variability term V_qy
% ---------------------
if ~(isfinite(s_y) ...
        && s_y > 0 ...
        && isfinite(s_q) ...
        && s_q > 0)
    return
end

alpha_qy = s_q / s_y;
V_qy = (1 - alpha_qy)^2;

% ---------------------
% correlation term C_qy
% ---------------------
z_y = a_y / s_y;
z_q = a_q / s_q;

rho_qy = mean(z_q .* z_y);
rho_qy = max(-1,min(1,rho_qy));

C_qy = (1 - rho_qy)^2;

% ----------------------------
% benchmark-matching term M_qy
% ----------------------------
switch Mdef

    case 1
        tol_m = max(1e-8,1e-6 * ...
            mean(abs(y_red),'omitnan'));

        goodM_red = abs(m_y_red) > tol_m;
        id_good = find(good);
        goodM(id_good(goodM_red)) = true;

        if ~any(goodM_red)
            return
        end

        ratio_qy = m_q_red(goodM_red) ./ ...
            m_y_red(goodM_red);
        errM = 1 - ratio_qy;

        eM = mean(errM);
        M_qy = mean(errM.^2);

    case 2
        goodM_red = true(size(m_y_red));
        id_good = find(good);
        goodM(id_good(goodM_red)) = true;

        ybar = mean(y_red,'omitnan');
        denM = norm(m_y_red - ybar);

        if ~(isfinite(denM) ...
                && denM > 0)
            return
        end

        Mraw = norm(m_q_red - m_y_red) ...
            / denM;
        M_qy = Mraw^2;

    otherwise
        error(['      Error:jkge: ' ...
            'unknown definition ' ...
            'selector M.']);
end

% ----
% JKGE
% ----
JKGE = 1 - sqrt(M_qy + V_qy + C_qy);

end