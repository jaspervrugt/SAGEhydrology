function [KGE,r_qy,nu_qy,zeta_qy] = kge(y,q,mu_y,std_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%KGE Compute Kling-Gupta efficiency and its components
%
% SYNOPSIS: [KGE,r_qy,nu_qy,zeta_qy] = kge(y,q,mu_y,std_y)
%   y       nx1 vector of observed discharge
%   q       nx1 vector of simulated discharge
%   mu_y    scalar mean of observed discharge
%   std_y   scalar standard deviation of observed discharge
%   KGE     OUTPUT: Kling-Gupta efficiency
%   r_qy    OUTPUT: linear correlation coefficient
%   nu_qy   OUTPUT: variability ratio (std_q / std_y)
%   zeta_qy OUTPUT: bias ratio (mean_q / mean_y)
%
% DESCRIPTION:
%   This function evaluates the Kling-Gupta efficiency (KGE), a composite
%   performance metric that balances correlation, variability, and bias
%   between simulated and observed discharge.
%
%   The KGE is defined as:
%       KGE = 1 - sqrt((r_qy - 1)^2 + (nu_qy - 1)^2 + (zeta_qy - 1)^2)
%
%   where r is the Pearson correlation coefficient between y and q, nu is
%   the ratio of standard deviations, and zeta is the ratio of means.
%
%   Missing or non-finite values are removed prior to computation.
%
%   Reference:
%     Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009),
%     Decomposition of the mean squared error and NSE performance criteria:
%     Implications for improving hydrological modelling,
%     Journal of Hydrology, 377(1–2), 80–91.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Apr. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [r_qy,nu_qy,zeta_qy,KGE] = deal(NaN);
    
    if isempty(y) ...
            || isempty(q)
        return
    end
    
    y = y(:);
    q = q(:);
    
    good = isfinite(y) ...
        & isfinite(q);
    y = y(good);
    q = q(good);

    n = numel(y);

    if n < 2
        return
    end
    
    % mu_q = sum(q) / n;
    % dq = q - mu_q;
    % std_q = sqrt(dot(dq,dq) / (n - 1));
    % 
    % nu_qy = std_q / std_y;
    % zeta_qy = mu_q / mu_y;
    % r_qy = dot(y - mu_y,dq) ...
    %     / ((n - 1) * std_y * std_q);
    % 
    % r_qy = max(-1,min(1,r_qy));
    % 
    % KGE = 1 - sqrt((r_qy - 1)^2 ...
    %         + (nu_qy - 1)^2 ...
    %         + (zeta_qy - 1)^2);

    mu_q = mean(q);
    std_q = std(q);

    if isfinite(std_y) ...
            && std_y > 0
        nu_qy = std_q / std_y;
    end

    if isfinite(mu_y) ...
            && mu_y ~= 0
        zeta_qy = mu_q / mu_y;
    end

    if isfinite(std_q) ...
            && std_q > 0 ...
            && isfinite(std_y) ...
            && std_y > 0
        r_qy = dot(y - mu_y, q - mu_q) ...
            / ((n - 1) * std_y * std_q);
        % Protect against tiny floating-point excursions.
        r_qy = max(-1, min(1, r_qy));        
    end

    if isfinite(r_qy) ...
            && isfinite(nu_qy) ...
            && isfinite(zeta_qy)
        KGE = 1 - sqrt((r_qy - 1)^2 + ...
            (nu_qy - 1)^2 + (zeta_qy - 1)^2);
    end
end