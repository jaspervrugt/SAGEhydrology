function [JKGE,M_qy,V_qy,C_qy] = ...
    jkge_score_given_mq(y,q,m_y,m_q,idx,loss)
%JKGE_SCORE_GIVEN_MQ Compute JKGE for a selected index set.
%
% Computes the Jawad Kling-Gupta Efficiency (JKGE) using observed
% discharge y, simulated discharge q, the observed benchmark m_y, and
% the simulated benchmark m_q. The calculation is restricted to the
% entries listed in idx and ignores nonfinite paired values.
%
% INPUT
%   y      nx1 observed discharge vector
%   q      nx1 simulated discharge vector
%   m_y    nx1 benchmark/mean-flow vector for observed discharge
%   m_q    nx1 benchmark/mean-flow vector for simulated discharge
%   idx    kx1 indices used for the JKGE calculation
%   loss   structure with JKGE options
%          .M = 1 : paper definition of the M component
%          .M = 2 : revised norm-based M component
%
% OUTPUT
%   JKGE   scalar Jawad Kling-Gupta Efficiency value
%   M_qy   scalar benchmark-matching component
%   V_qy   scalar variability component
%   C_qy   scalar correlation component
%
% NOTES
%   - If loss.M is missing/empty, the revised norm-based definition is used
%   - Returns NaN when fewer than two valid paired values are available
%     or when required normalization terms are zero/nonfinite.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, April 2026                                %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6 ...
        || ~isfield(loss,'M') ...
        || isempty(loss.M)
    def = 2;   % revised norm-based version of M
else
    def = loss.M;
end

[JKGE,M_qy,V_qy,C_qy] = deal(NaN);

y = double(y(:));
q = double(q(:));
m_y = double(m_y(:));
m_q = double(m_q(:));

n = numel(q);

idx = idx(:);
idx = idx(idx >= 1 & idx <= n);

good = false(n,1);
good(idx) = true;

good = good ...
    & isfinite(y) ...
    & isfinite(q) ...
    & isfinite(m_y) ...
    & isfinite(m_q);

if nnz(good) < 2
    return
end

y_red = y(good);
m_y_red = m_y(good);
m_q_red = m_q(good);

a_y = y_red - m_y_red;
a_q = q(good) - m_q_red;

s_y = sqrt(mean(a_y.^2));
s_q = sqrt(mean(a_q.^2));

if ~(isfinite(s_y) ...
        && s_y > 0 ...
        && isfinite(s_q) ...
        && s_q > 0)
    return
end

alpha_qy = s_q / s_y;
V_qy = (1 - alpha_qy)^2;

z_y = a_y / s_y;
z_q = a_q / s_q;

rho_qy = mean(z_q .* z_y);
rho_qy = max(-1,min(1,rho_qy));
C_qy = (1 - rho_qy)^2;

switch def
    case 1
        Mraw = mean(m_q_red) ...
            / mean(m_y_red) - 1;
        M_qy = Mraw^2;

    case 2
        ybar = mean(y_red,'omitnan');
        denM = norm(m_y_red - ybar);

        if ~(isfinite(denM) ...
                && denM > 0)
            return
        end

        Mraw = norm(m_q_red - ...
            m_y_red) / denM;
        M_qy = Mraw^2;
    otherwise
        error(['Unsupported JKGE ' ...
            'M definition: %g'],def);
end

JKGE = 1 - sqrt(M_qy + V_qy + C_qy);

end
