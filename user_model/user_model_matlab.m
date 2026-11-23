function [q,J,Jth,Z,swe,Jswe,sm,Jsm,failed]=user_model_matlab(par,mdl,data,ode)
%USER_MODEL_MATLAB GR4JB template, analytic sensitivities and MATLAB solvers.
% The fixed routing-memory layout and equations match user_model_prepare.cpp.
    par=double(par(:));d=numel(par);N=mdl.tout;first=mdl.idx(1);last=mdl.idx(2)-1;
    switch mdl.pspace
        case 0,th=par;Jth=ones(d,1);valid=all(th>=mdl.th_min & th<=mdl.th_max);
        case 1,th=mdl.th_min+par.*(mdl.th_max-mdl.th_min);Jth=mdl.th_max-mdl.th_min;valid=all(par>=0 & par<=1);
        case 2,h=.5*(1+tanh(par/2));Jth=(mdl.th_max-mdl.th_min).*h.*(1-h);th=mdl.th_min+h.*(mdl.th_max-mdl.th_min);valid=true;
        otherwise,error('user_model:ParameterSpace','Unknown parameter space.');
    end
    validate_initial_storage(mdl,th);
    L=mdl.L;m=4+2*L;C=mdl.config;
    if numel(mdl.y0)~=m,error('user_model:InitialStates','Expected %d initial states.',m);end
    validateattributes(ode.InitStep,{'numeric'},{'scalar','positive','finite'});
    [U,dUdx4]=gr4j_UH_ord(th(4),C.mts,C.eta,L);
    z=[mdl.y0(:);zeros(m*d,1)];mem=isfield(ode,'mem') && logical(ode.mem);
    if mem,Z=nan(N+1,m*(d+1));Z(1,:)=z';else,Z=z';end
    q=nan(N,1);J=nan(N,d);swe=q;sm=q;Jswe=J;Jsm=J;failed=~valid;
    for s=1:N
        if failed,break,end
        before=z;rhs=@(~,state)user_model_aug_ode(s-1,state,th,data,U,dUdx4,L,C.eta,C.kappa,C.bg,C.bR,C.mts,C.T_sm,C.eps_m,C.eps_s,C.rho,m,d);
        switch mdl.mcode
            case 1
                t=0;step=max(ode.MinStep,min(ode.InitStep,ode.MaxStep));count=0;
                while t<1
                    step=min(step,1-t);a=rhs(0,z);e=z+step*a;b=rhs(0,e);candidate=z+.5*step*(a+b);
                    errorNorm=sqrt(mean(((e-candidate)./(ode.RelTol*abs(candidate)+ode.AbsTol)).^2));
                    if any(~isfinite(candidate)) || any(abs(candidate)>1e12),failed=true;break,end
                    if errorNorm<=1 || step<=ode.MinStep,z=candidate;t=t+step;end
                    count=count+1;if count>=ode.maxiter,failed=true;break,end
                    factor=max(.2,min(5,.9*errorNorm^(-1/ode.Order)));step=max(ode.MinStep,min(ode.MaxStep,step*factor));
                end
            case 2
                settings=odeset('InitialStep',min(1,ode.InitStep),'MaxStep',ode.MaxStep,'RelTol',ode.RelTol,'AbsTol',ode.AbsTol);
                [~,trajectory]=ode45(rhs,[0 1],z,settings);z=trajectory(end,:)';
            case 3
                steps=1000;if isfield(ode,'int_steps'),steps=ode.int_steps;end
                validateattributes(steps,{'numeric'},{'scalar','integer','positive','finite'});
                for k=1:steps,z=z+rhs(0,z)/steps;end
        end
        if failed || any(~isfinite(z)) || any(abs(z)>1e12),failed=true;break,end
        S=reshape(z(m+1:end),m,d);S0=reshape(before(m+1:end),m,d);
        q(s)=z(m)-before(m);J(s,:)=(S(m,:)-S0(m,:)).*Jth';
        swe(s)=z(1);sm(s)=z(2);Jswe(s,:)=S(1,:).*Jth';Jsm(s,:)=S(2,:).*Jth';
        if mem,Z(s+1,:)=z';else,Z=z';end
    end
    take=first:last;q=q(take);J=J(take,:);swe=swe(take);sm=sm(take);Jswe=Jswe(take,:);Jsm=Jsm(take,:);
end
function dzdt = user_model_aug_ode(t,z,th,data,U,dUdx4,L,eta,kappa,bg,bR, ...
    mts,T_sm,eps_m,eps_s,rho,m,d)

x = z(1:m);                                     % x = [Swe S R U1 U2 Q]
S = reshape(z(m+1:m*(d+1)),m,d);                % mxd sensitivity matrix
[dxdt,dSdt] = user_model_odefcn(t,x,th,S,data, ...   % compute dxdt & dSdt
    U,dUdx4,L,eta,kappa,bg,bR,mts,T_sm,eps_m, ...
    eps_s,rho,m,d);
dzdt = [dxdt; dSdt(:)];                         % Repack single vector

end

function [dudt,dSdt,Jth_f,Jx_f] = user_model_odefcn(t,u,th,S,data, ...
    U,dUdx4,L,~,kappa,bg,bR,mts,T_sm,eps_m,eps_s,~,m,d)
%user_model_ODEFCN  ODE stores + analytic routing via UH convolution
% (Mathias-style CIF)
%
% Notes
% -----
% This "B" variant follows the unit-hydrograph convolution approach used by
% Mathias (CIF/CKF routing kernel), i.e., the two routing branches are
% obtained via discrete convolution of the effective rainfall (after
% production/percolation) rather than via an explicit Nash-cascade state
% augmentation. 
% eta = ~ and rho = ~
%
% State layout:
% [Swe Sp Rr z1(1:L) z2(1:L) Qinf]
% Returns:
% dudt  : m x 1
% dSdt  : m x d
% Jx_f  : m x m
% Jth_f : m x d

u = u(:)';

% -------
% helpers
% -------
smooth_pos = @(a,eps) 0.5*(a + sqrt(a.^2 + eps.^2));
d_smooth_pos_da = @(a,eps) 0.5*(1 + a./sqrt(a.^2 + eps.^2));

smooth_min = @(A,B,eps) 0.5*(A + B - sqrt((A-B).^2 + eps.^2));
d_smooth_min_dA = @(A,B,eps) 0.5*(1 - (A-B)./sqrt((A-B).^2 + eps.^2));
d_smooth_min_dB = @(A,B,eps) 0.5*(1 + (A-B)./sqrt((A-B).^2 + eps.^2));

% -------
% indices
% -------
iSwe = 1;
iSp = 2;
iRr = 3;
iz1 = 4;
iz2 = iz1 + L;
iQ = iz2 + L;

% ----------
% parameters
% ----------
x1 = th(1);
x2 = th(2);
x3 = th(3);
%x4 = th(4); % used through U outside; gradient w.r.t x4 needs dU/dx4 if desired
x5 = th(5);
f_p = th(6);
T_tr = th(7);
f_dd = th(8);

% ----------------------------
% forcing (piecewise constant)
% ----------------------------
id = floor(t) + 1;
P = data.P(id,1);
Ep = data.Ep(id,1);
T = data.T(id,1);

% -----------------------------
% smooth positivity on storages
% -----------------------------
Swe_u = u(iSwe);
Sp_u = u(iSp);
Rr_u = u(iRr);

Swe = smooth_pos(Swe_u,eps_s);
Sp = smooth_pos(Sp_u,eps_s);
Rr = smooth_pos(Rr_u,eps_s);

dSwe_du = d_smooth_pos_da(Swe_u,eps_s);
dSp_du = d_smooth_pos_da(Sp_u,eps_s);
dRr_du = d_smooth_pos_da(Rr_u,eps_s);

z_1 = u(iz1:iz1+L-1);    % do NOT smooth: these are linear memory states
z_2 = u(iz2:iz2+L-1);

% ---------------------------------------------
% 6) snow module (from user_model, plus derivatives)
% ---------------------------------------------
T_smeps_m = max(T_sm, eps_m);
utemp = (T - T_tr) / T_smeps_m;
snow_fr = 0.5*(1 - tanh(utemp));
rain_fr = 1 - snow_fr;

p_snow = P * snow_fr;
p_rain = P * rain_fr;

aT = (T - T_tr);
posT = 0.5*(aT + sqrt(aT^2 + T_smeps_m^2));
p_pmelt = f_dd * posT;

p_amelt = smooth_min(Swe, p_pmelt, eps_m);
p_liq = p_rain + p_amelt;

% snow derivatives
sech2 = 1/(cosh(utemp)^2);
dsnowFrac_dT_tr = 0.5 * sech2 / T_smeps_m;

dposT_da = 0.5*(1 + aT/sqrt(aT^2 + T_smeps_m^2));
dposT_dTtr = -dposT_da;

dp_pmelt_dTtr = f_dd * dposT_dTtr;
dp_pmelt_dfdd = posT;

dM_dSwe = d_smooth_min_dA(Swe, p_pmelt, eps_m);
dM_dp_pmelt = d_smooth_min_dB(Swe, p_pmelt, eps_m);

dp_liq_dSwe = dM_dSwe;

% -------------------
% 7) production block
% -------------------
Ep_eff = f_p * Ep;              % mm/T
aPE = p_liq - Ep_eff;           % mm/T
daPE_dfp = -Ep;                 % f_p derivative enters only through aPE

Pn = smooth_pos( aPE, eps_s);
En = smooth_pos(-aPE, eps_s);

dPn_daPE =  d_smooth_pos_da( aPE, eps_s);
dEn_daPE = -d_smooth_pos_da(-aPE, eps_s);

% wetness gate
uPE = aPE / x1;
wWet = 0.5*(1 + tanh(uPE));
sech2PE = 1/(cosh(uPE)^2);

dwWet_daPE = 0.5 * sech2PE / x1;
dwWet_dx1 = 0.5 * sech2PE * (-aPE/(x1^2));

Sp_n = Sp / x1;

% percolation
Sp2 = Sp_n*Sp_n; Sp4 = Sp2*Sp2;
g = 1 + (kappa^4)*Sp4;
g_m14 = g^(-0.25);
Pp = Sp * (1 - g_m14);

dSpn_dSp = 1/x1;
dSpn_dx1 = -Sp/(x1^2);

dg_dSpn = 4*(kappa^4) * (Sp_n^3);
dg_m14_dg = -0.25 * g^(-1.25);

dPp_dSp = (1 - g_m14) + Sp * ( -dg_m14_dg * dg_dSpn * dSpn_dSp );
dPp_dx1 = Sp * ( -dg_m14_dg * dg_dSpn * dSpn_dx1 );

dPerc_dS = dPp_dSp;
dPerc_dx1 = dPp_dx1;

% wet branch Ps_wet
zP = Pn / x1;
aP = tanh(zP);
daP_dPn = (1 - aP^2)/x1;

Nw = (1 - Sp_n^2) * aP;
Dw = (1 + Sp_n*aP);
Ps_wet = x1 * Nw / Dw;

dNw_dSpn = (-2*Sp_n)*aP;
dNw_dPn = (1 - Sp_n^2)*daP_dPn;

dDw_dSpn = aP;
dDw_dPn = Sp_n*daP_dPn;

dPsWet_dSpn = x1 * (Dw*dNw_dSpn - Nw*dDw_dSpn) / (Dw^2);
dPsWet_dPn = x1 * (Dw*dNw_dPn  - Nw*dDw_dPn ) / (Dw^2);

dPsWet_dS = dPsWet_dSpn * dSpn_dSp;

dzP_dx1 = -Pn/(x1^2);
daP_dx1 = (1 - aP^2) * dzP_dx1;

dNw_dx1 = dNw_dSpn*dSpn_dx1 + (1 - Sp_n^2)*daP_dx1;
dDw_dx1 = dDw_dSpn*dSpn_dx1 + Sp_n*daP_dx1;

dNoverD_dx1 = (Dw*dNw_dx1 - Nw*dDw_dx1)/(Dw^2);
dPsWet_dx1 = (Nw/Dw) + x1 * dNoverD_dx1;

% dry branch Es_dry
zE = En / x1;
aE = tanh(zE);
daE_dEn = (1 - aE^2)/x1;

Nd = Sp * (2 - Sp_n) * aE;
Dd = 1 + (1 - Sp_n) * aE;

invD = 1/Dd;
invD2 = invD^2;

Es_dry = Nd * invD;

dNd_dSp = (2 - Sp_n) * aE;
dNd_dSpn = -Sp * aE;
dNd_daE = Sp * (2 - Sp_n);

dDd_dSpn = -aE;
dDd_daE = (1 - Sp_n);

dEs_dSp = dNd_dSp * invD;
dEs_dSpn = (dNd_dSpn*Dd - Nd*dDd_dSpn) * invD2;
dEs_daE = (dNd_daE*Dd  - Nd*dDd_daE ) * invD2;

dEsDry_dS = dEs_dSp + dEs_dSpn * dSpn_dSp;
dEsDry_dEn = dEs_daE * daE_dEn;

dzE_dx1 = -En/(x1^2);
daE_dx1 = (1 - aE^2) * dzE_dx1;
dEsDry_dx1 = dEs_dSpn * dSpn_dx1 + dEs_daE * daE_dx1;

% smooth blend
Ps = wWet * Ps_wet;
Es = (1 - wWet) * Es_dry;

dPs_dS = wWet * dPsWet_dS;
dEs_dS = (1 - wWet) * dEsDry_dS;

dPsWet_daPE = dPsWet_dPn * dPn_daPE;
dEsDry_daPE = dEsDry_dEn * dEn_daPE;

dPs_daPE = wWet*dPsWet_daPE + Ps_wet*dwWet_daPE;
dEs_daPE = (1 - wWet)*dEsDry_daPE - Es_dry*dwWet_daPE;

dPs_dp_liq = dPs_daPE;
dEs_dp_liq = dEs_daPE;

dPs_dx1 = wWet*dPsWet_dx1 + Ps_wet*dwWet_dx1;
dEs_dx1 = (1 - wWet)*dEsDry_dx1 - Es_dry*dwWet_dx1;

% remaining production fluxes
P1 = Pn - Ps;

dP1_dS = -dPs_dS;
dP1_daPE = dPn_daPE - dPs_dp_liq;
dP1_dx1 = -dPs_dx1;

Pr = P1 + Pp;

dPr_dS = dP1_dS   + dPerc_dS;
dPr_dp_liq = dP1_daPE;
dPr_dx1 = dP1_dx1  + dPerc_dx1;

% f_p derivatives (aPE-only)
dPs_dfp = dPs_daPE * daPE_dfp;
dEs_dfp = dEs_daPE * daPE_dfp;
dP1_dfp = dP1_daPE * daPE_dfp;
dPr_dfp = dP1_dfp;

% snow chain terms for parameters
dPsnow_dTtr = P * dsnowFrac_dT_tr;
dPrain_dTtr = -P * dsnowFrac_dT_tr;

dM_dTtr = dM_dp_pmelt * dp_pmelt_dTtr;
dM_dfdd = dM_dp_pmelt * dp_pmelt_dfdd;

dp_liq_dTtr = dPrain_dTtr + dM_dTtr;
dp_liq_dfdd = dM_dfdd;

dPs_dTtr = dPs_dp_liq * dp_liq_dTtr;
dPs_dfdd = dPs_dp_liq * dp_liq_dfdd;

dEs_dTtr = dEs_dp_liq * dp_liq_dTtr;
dEs_dfdd = dEs_dp_liq * dp_liq_dfdd;

dP1_dTtr = dP1_daPE * dp_liq_dTtr;
dP1_dfdd = dP1_daPE * dp_liq_dfdd;

dPr_dTtr = dP1_dTtr;   % Pp independent of Ttr, fdd
dPr_dfdd = dP1_dfdd;

% derivatives wrt Swe (through melt -> p_liq)
dPr_dSwe = dPr_dp_liq * dp_liq_dSwe;

% ---------------------------------
% routing inputs and UH memory ODEs
% ---------------------------------
q_in1 = x5 * Pr;
q_in2 = (1 - x5) * Pr;

dz_1 = zeros(L,1);
dz_2 = zeros(L,1);

dz_1(1) = q_in1 - z_1(1);
dz_2(1) = q_in2 - z_2(1);
for k=2:L
    dz_1(k) = z_1(k-1) - z_1(k);
    dz_2(k) = z_2(k-1) - z_2(k);
end

q_1 = z_1 * U;   % scalar
q_2 = z_2 * U;   % scalar

% routing store + discharge
Rratio = Rr / x3;
q_g = x2 * (Rratio^bg);
q_r = x3/((bR-1)*mts) * (Rratio^bR);

q_d = q_2 + q_g;
dQdpos_dQd = d_smooth_pos_da(q_d, eps_s);
q_out = q_r + smooth_pos(q_d, eps_s);

% derivatives for routing store
dqg_dRr = x2 * bg * (Rratio^(bg-1)) * (1/x3);
dqg_dx2 = (Rratio^bg);
dqg_dx3 = x2 * bg * (Rratio^(bg-1)) * (-Rr/(x3^2));

c0 = 1/((bR-1)*mts);
dqr_dRr = c0 * bR * (Rr^(bR-1)) * (x3^(1-bR));
dqr_dx3 = c0 * (Rr^bR) * (1-bR) * (x3^(-bR));

% explicit x4 dependence through U(x4)
dq1_dx4 = z_1 * dUdx4;   % scalar
dq2_dx4 = z_2 * dUdx4;   % scalar

% -------
% ODE RHS
% -------
dudt = zeros(m,1);
dudt(iSwe) = p_snow - p_amelt;
dudt(iSp) = Ps - Es - Pp;
dudt(iRr) = q_1 + q_g - q_r;
dudt(iz1:iz1+L-1) = dz_1;
dudt(iz2:iz2+L-1) = dz_2;
dudt(iQ) = q_out;

% ------------------------
% Jacobians Jx_f and Jth_f
% ------------------------
Jx_f = zeros(m,m);
Jth_f = zeros(m,d);

% --- physical-store Jacobian pieces (same logic as user_model, mapped) ---
% Swe equation: dudt(Swe)=p_snow - p_amelt
Jx_f(iSwe,iSwe) = -dM_dSwe;

% Sp equation: dudt(Sp)=Ps - Es - Pp
Jx_f(iSp,iSwe) = (dPs_dp_liq - dEs_dp_liq) * dp_liq_dSwe; % via melt
Jx_f(iSp,iSp) = (dPs_dS - dEs_dS - dPerc_dS);            % via Sp

% UH memory injections depend on Pr(Swe,Sp)
Jx_f(iz1, iSwe) = x5 * dPr_dSwe;
Jx_f(iz1, iSp) = x5 * dPr_dS;
Jx_f(iz1, iz1) = -1;

% dz2(1) = (1-x5)*Pr - z2(1)
Jx_f(iz2, iSwe) = (1-x5) * dPr_dSwe;
Jx_f(iz2, iSp) = (1-x5) * dPr_dS;
Jx_f(iz2, iz2) = -1;

% shift register dynamics
for k=2:L
    Jx_f(iz1+k-1, iz1+k-2) = 1;
    Jx_f(iz1+k-1, iz1+k-1) = -1;
    Jx_f(iz2+k-1, iz2+k-2) = 1;
    Jx_f(iz2+k-1, iz2+k-1) = -1;
end

% routed inflow to routing store: q1 depends on z1; q2 affects q_out (via qd)
Jx_f(iRr, iz1:iz1+L-1) = U';
Jx_f(iQ,  iz2:iz2+L-1) = dQdpos_dQd * U';

% routing store dependence on Rr
Jx_f(iRr, iRr) = dqg_dRr - dqr_dRr;
Jx_f(iQ,  iRr) = dqr_dRr + dQdpos_dQd * dqg_dRr;

% chain rule for smooth_pos'd physical states
Jx_f(:,iSwe) = Jx_f(:,iSwe) * dSwe_du;
Jx_f(:,iSp) = Jx_f(:,iSp)  * dSp_du;
Jx_f(:,iRr) = Jx_f(:,iRr)  * dRr_du;

% ------------------------
% Parameter Jacobian Jth_f
% ------------------------

% (7) T_tr and (8) f_dd affect Swe + Sp via snow/melt and p_liq
Jth_f(iSwe,7) = dPsnow_dTtr - dM_dTtr;
Jth_f(iSwe,8) = -dM_dfdd;

Jth_f(iSp,7) = dPs_dTtr - dEs_dTtr;
Jth_f(iSp,8) = dPs_dfdd - dEs_dfdd;

% x1 affects Sp equation (Ps,Es,Pp)
Jth_f(iSp,1) = dPs_dx1 - dEs_dx1 - dPerc_dx1;

% f_p affects Sp equation through aPE
Jth_f(iSp,6) = dPs_dfp - dEs_dfp;

% routing store parameter effects (x2,x3)
Jth_f(iRr,2) = dqg_dx2;
Jth_f(iRr,3) = dqg_dx3 - dqr_dx3;

Jth_f(iQ,2) = dQdpos_dQd * dqg_dx2;
Jth_f(iQ,3) = dqr_dx3 + dQdpos_dQd * dqg_dx3;

Jth_f(iRr,4) = Jth_f(iRr,4) + dq1_dx4;
Jth_f(iQ, 4) = Jth_f(iQ, 4) + dQdpos_dQd * dq2_dx4;

% UH memory injections: explicit parameter effects
% dz1(1) depends on x5 and on Pr's parameter sensitivities
Jth_f(iz1,5) = Pr;                         % d(x5*Pr)/dx5
Jth_f(iz2,5) = -Pr;                        % d((1-x5)*Pr)/dx5

Jth_f(iz1,1) = x5 * dPr_dx1;
Jth_f(iz2,1) = (1-x5) * dPr_dx1;

Jth_f(iz1,6) = x5 * dPr_dfp;
Jth_f(iz2,6) = (1-x5) * dPr_dfp;

Jth_f(iz1,7) = x5 * dPr_dTtr;
Jth_f(iz2,7) = (1-x5) * dPr_dTtr;

Jth_f(iz1,8) = x5 * dPr_dfdd;
Jth_f(iz2,8) = (1-x5) * dPr_dfdd;

% ---------------
% sensitivity ODE
% ---------------
dSdt = Jx_f * S + Jth_f;

end

%% 4. GR4J: RHS-only for J(x)_f and J(th)_f check
% GR4J unit hydrograph functions
function [U,dUdx4] = gr4j_UH_ord(x4,mts,eta,N)

p = eta + 1;
N = max(1,N);

t_k = (1:N)' * mts;
t_km1 = t_k - mts;

Gk = gr4j_CIF(t_k,x4,eta);
Gkm1 = gr4j_CIF(t_km1,x4,eta);

dGk = gr4j_dCIF_dx4(t_k,x4,p);
dGkm1 = gr4j_dCIF_dx4(t_km1,x4,p);

Uraw = Gk - Gkm1;
dUrawdx = dGk - dGkm1;
eps_U = 1e-12 + 1e-6*max(1, max(abs(Uraw)));   % scale-aware

Uraw_pos = 0.5*(Uraw + sqrt(Uraw.^2 + eps_U^2));
dpos_dUraw = 0.5*(1 + Uraw ./ sqrt(Uraw.^2 + eps_U^2));

Uraw = Uraw_pos;
dUrawdx = dUrawdx .* dpos_dUraw;
% Match the C++ preparation layer's exact zero tail.
tail=((0:N-1)'*mts)>=2*x4;Uraw(tail)=0;dUrawdx(tail)=0;

s = sum(Uraw);
ds = sum(dUrawdx);

s_safe = max(s, 1e-30);
U = Uraw / s_safe;
dUdx4 = (dUrawdx*s_safe - Uraw*ds) / (s_safe*s_safe);
end

function dG = gr4j_dCIF_dx4(t, x4, p)
% p = eta+1
dG = zeros(size(t));
u = t ./ x4;

% 0 < t < x4: G = u^p  => dG/dx4 = -(p/x4)*u^p
i2 = (t > 0) & (t < x4);
dG(i2) = -(p/x4) * (u(i2).^p);

% x4 <= t < 2x4: G = 1-(2-u)^p => dG/dx4 = -p*(2-u)^(p-1) * (t/x4^2)
i3 = (t >= x4) & (t < 2*x4);
dG(i3) = -p * (2 - u(i3)).^(p-1) .* (t(i3) ./ (x4^2));
end

function G = gr4j_CIF(t, x4, eta)
%GR4J_CIF  Cumulative integral function G(t) for the GR4J CKF.
% Piecewise over t in (-inf,0), [0,x4], [x4,2*x4], [2*x4,inf).
%
% This corresponds to the symmetric GR4J CKF support [0,2*x4]
% used by Mathias for runoff attenuation.

td = t ./ x4;   % nondimensional time
G = zeros(size(t));

% region 1: t <= 0 -> 0
i1 = (t <= 0);
G(i1) = 0;

% region 2: 0 < t < x4  (0 < td < 1)
i2 = (t > 0) & (t < x4);
% integral of g(t) = ((eta+1)/x4) * (t/x4)^eta over [0,t]
% => G = (t/x4)^(eta+1)
G(i2) = td(i2).^(eta+1);

% region 3: x4 <= t < 2*x4 (1 <= td < 2)
i3 = (t >= x4) & (t < 2*x4);
% symmetric second branch:
% G = 1 - (2 - td)^(eta+1)
G(i3) = 1 - (2 - td(i3)).^(eta+1);

% region 4: t >= 2*x4 -> 1
i4 = (t >= 2*x4);
G(i4) = 1;
end