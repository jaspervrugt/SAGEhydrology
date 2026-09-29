function delta = delta_n(loss_fnc,y_n,q_n,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DELTA_N Compute loss sensitivity to simulated discharge.
%
%  Return the derivative or subgradient of a selected discharge loss with
%  respect to each simulated value.
%
% SYNOPSIS:
%   delta = delta_n(loss_fnc,y_n,q_n,varargin)
%
% INPUT ARGUMENTS:
%   loss_fnc        loss code: 1 SAR; 2 GLS; 3 NSE; 4 KGE; 6 FDC
%   y_n             measured-discharge vector
%   q_n             simulated-discharge vector
%   varargin        loss-specific additional argument
%    {1}             measurement-error covariance for GLS (loss 2)
%    {1}             FDC formulation 1, 2, or 3 for loss 6
%
% OUTPUT ARGUMENTS:
%   delta           sensitivity vector d(loss)/d(q_n)
%
% NOTES:
%   FDC formulations 1-3 use physical discharge, probability, and
%   log-probability, respectively; omitted formulation defaults to 1.
%   Huber (5) and JKGE (7) use dedicated gradient routines.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025 / updated Apr. 2026             %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    e_n = y_n - q_n;                    % n x 1 vector of residuals
    
    switch loss_fnc
        
        case 1 % Sum of absolute residuals
            s_n = sign(e_n);            % sample mean of y_n
            delta = - s_n;              % δ(theta)
        
        case 2 % Generalized least squares
            Sigma_eps = varargin{1};
            delta = - 2*(Sigma_eps\e_n);% δ(theta)
        
        case 3 % Nash-Sutcliffe efficiency
            m_y = mean(y_n);            % sample mean of y_n
            SSt = sum((y_n - m_y).^2);  % total sum of squares
            delta = -2/SSt * e_n;       % δ(theta)
        
        case 4 % Kling-Gupta efficiency
            n = numel(q_n);             % # elements of q_n = (q_1,...,q_n)'
            m_y = mean(y_n);            % sample mean of y_n
            m_q = mean(q_n);            % sample mean of q_n
            C = cov(y_n,q_n);           % sample covariance matrix y_n and q_n
            s_y = sqrt(C(1));           % sample standard deviation of y_n
            s_q = sqrt(C(4));           % sample standard deviation of q_n
            r = C(2)/(s_y*s_q);         % sample correlation coefficient
            v = s_q/s_y;                % variability ratio
            z = m_q/m_y;                % bias ratio
            L_jkge = sqrt((r-1)^2 + ...
                (v-1)^2+(z-1)^2);
            dm_qdqt = ones(n,1)/n;
            ds_qdqt = (q_n-m_q)/((n-1)*s_q);
            drdqt = (y_n-m_y)/((n-1)*s_q*s_y) ...
                - r*(q_n-m_q)/((n-1)*s_q^2);
            dvdqt = ds_qdqt/s_y;
            dzdqt = dm_qdqt/m_y;
            delta = ((r-1)*drdqt+(v-1)*dvdqt+(z-1)*dzdqt)/L_jkge; % δ(theta)
    
        case 5 % Huber loss
            % Without need for statistics toolbox
            error(['      Error:delta_n: ' ...
                'Use huber_loss function ' ...
                'instead']);
        
        case 6 % Flow-duration-curve divergence
            % O(n log n) derivative of the energy-distance/FDC objective.
            %
            % For each simulated discharge q_i,
            %
            %   dL/dq_i = (1/n^2) [sum_j sign(q_i-y_j)
            %                      - sum_j sign(q_i-q_j)].
            %
            % Instead of evaluating the two sign sums explicitly in O(n^2),
            % obtain them from sorted ranks/counts:
            %
            %   sum_j sign(x-v_j) = N_< - N_>
            %                     = 2*N_< + N_= - n.
            %
            % Sorting y_n and q_n dominates the computational cost, giving
            % O(n log n) complexity. Ties handled exactly via lower/upper
            % rank counts, so this is algebraically equivalent to original
            % pairwise implementation.    
            y = double(y_n(:));
            q = double(q_n(:));
            n = numel(q);
    
            if numel(y) ~= n
                error(['      Error:delta_n: ' ...
                    'y_n and q_n must have the same length.']);
            end
    
            if n == 0
                delta = zeros(0,1);
                return
            end
    
            ys = sort(y);
            [qs,ord] = sort(q);
            if isempty(varargin)
                formulation = 1;
            else
                formulation = local_fdc_formulation(varargin{1});
            end

            % Probability-space formulations pair equal-probability
            % order statistics. Mapping through ord returns dL/dq_n in
            % the original time order. At exact ties this is a valid
            % stable-sort subgradient.
            if formulation == 2
                delta = zeros(n,1);
                delta(ord) = (2/n)*(qs - ys);
                return
            elseif formulation == 3
                yp = ys(ys > 0);
                if isempty(yp)
                    delta = nan(n,1);
                    return
                end
                q0 = 0.01*median(yp);
                if any(qs + q0 <= 0)
                    delta = nan(n,1);
                    return
                end
                delta = zeros(n,1);
                delta(ord) = (2/n)*( ...
                    log(qs + q0) - log(ys + q0))./(qs + q0);
                return
            end

            % -------------------------------------
            % First sign sum: sum_j sign(q_i - y_j)
            % -------------------------------------
            %
            % Sweep through sorted q and sorted y. For each group equal q,
            % count the observations strictly below q and equal to q.
            s_qy_sorted = zeros(n,1);
            jy = 1;    
            iq = 1;
            while iq <= n
                x = qs(iq);    
                % Find the complete tie block in q.
                iq2 = iq;
                while iq2 < n ...
                        && qs(iq2+1) == x
                    iq2 = iq2 + 1;
                end    
                % Advance y pointer over entries strictly below x.
                while jy <= n ...
                        && ys(jy) < x
                    jy = jy + 1;
                end
                nLess = jy - 1;    
                % Count entries equal to x.
                jy2 = jy;
                while jy2 <= n ...
                        && ys(jy2) == x
                    jy2 = jy2 + 1;
                end
                nEqual = jy2 - jy;    
                s = 2*nLess + nEqual - n;
                s_qy_sorted(iq:iq2) = s;    
                iq = iq2 + 1;
            end    
            % --------------------------------------
            % Second sign sum: sum_j sign(q_i - q_j)
            % --------------------------------------
            %
            % For a sorted tie block iq:iq2,
            %   N_< = iq - 1
            %   N_= = iq2 - iq + 1
            % exactly.
            s_qq_sorted = zeros(n,1);
            iq = 1;
            while iq <= n
                x = qs(iq);
                iq2 = iq;
                while iq2 < n ...
                        && qs(iq2+1) == x
                    iq2 = iq2 + 1;
                end    
                nLess = iq - 1;
                nEqual = iq2 - iq + 1;
                s = 2*nLess + nEqual - n;
                s_qq_sorted(iq:iq2) = s;    
                iq = iq2 + 1;
            end    
            % Map from sorted-q order back to the original time order.
            delta = zeros(n,1);
            delta(ord) = (s_qy_sorted - s_qq_sorted) / n^2;
    
        case 7 % Jawad-Kling-Gupta efficiency
    
            error(['      Error:delta_n: ' ...
                'Use jkge_grad function ' ...
                'instead']);
    
        otherwise
            error(['      Error:delta_n: ' ...
                'Unknown loss function ' ...
                'index: %d.'],loss_fnc);
    end

end

function formulation = local_fdc_formulation(value)
%LOCAL_FDC_FORMULATION Normalize the public 6a/6b/6c selection.

    if isnumeric(value) && isscalar(value) && isfinite(value)
        formulation = double(value);
    else
        key = lower(regexprep(char(string(value)),'[^a-z0-9]',''));
        switch key
            case {'1','a','fdc','dfdc','physical','physicalcdf'}
                formulation = 1;
            case {'2','b','p','dp','quantile'}
                formulation = 2;
            case {'3','c','logp','dlogp','logquantile'}
                formulation = 3;
            otherwise
                formulation = NaN;
        end
    end
    if ~ismember(formulation,1:3)
        error('delta_n:BadFDCFormulation', ...
            'FDC formulation must be 1 (d_fdc), 2 (d_p), or 3 (d_logp).');
    end
end
