# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_cfe_nwm_equations.py
#
# Python equations and analytic sensitivities matching the shared
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Python equations and analytic sensitivities matching the shared
cfe_nwm core.
Equations by Jasper A. Vrugt; Python translation, October 2026.
No C++ or MATLAB runtime is used by these equations.
"""

import numpy as np
from numba import njit


@njit(cache=True)
def sigmoid(x):
    """Evaluate the smooth tanh-based switching function."""
    return 0.5 * (1 + np.tanh(x))


@njit(cache=True, error_model="numpy")
def smooth_pos(a, ep):
    """Evaluate the differentiable positive-part approximation."""
    return 0.5 * (a + np.sqrt(((a * a) + (ep * ep))))


@njit(cache=True, error_model="numpy")
def dsmooth_pos_da(a, ep):
    """Evaluate the analytic derivative of the smooth_pos response."""
    return 0.5 * (1.0 + (a / np.sqrt(((a * a) + (ep * ep)))))


@njit(cache=True, error_model="numpy")
def smooth_min2(a, b, ep):
    """Evaluate the smooth minimum of two values."""
    d = a - b
    return 0.5 * ((a + b) - np.sqrt(((d * d) + (ep * ep))))


@njit(cache=True, error_model="numpy")
def dmin2_da(a, b, ep):
    """Evaluate the analytic derivative of the min2 response."""
    d = a - b
    r = np.sqrt(((d * d) + (ep * ep)))
    return 0.5 * (1.0 - (d / r))


@njit(cache=True, error_model="numpy")
def cfe_nwm_odefcn(
    u,
    udot,
    Smat,
    dSdt,
    P,
    Ep,
    T,
    s_max,
    s_fc,
    s_wp,
    k_sch,
    a1,
    k_perc,
    lf_thr,
    a2,
    k_lf,
    g_max,
    c_gw,
    mm,
    k_nsh,
    giuh,
    T_tr,
    f_dd,
    K,
    L,
    dT,
    T_sm,
    eps_m,
    eps,
    rho,
    m,
    d,
):
    """Evaluate the cfe_nwm state equations and their derivatives."""
    Swe_u = u[0]
    S = u[1]
    G = u[2]
    Q = u[3:]
    F = u[(3 + L) :]
    Swe = smooth_pos(Swe_u, eps_m)
    dSwe_dSweu = dsmooth_pos_da(Swe_u, eps_m)
    T_smeps_m = max(T_sm, eps_m)
    uT = (T - T_tr) / T_smeps_m
    tanhu = np.tanh(uT)
    coshu = np.cosh(uT)
    sech2 = 1.0 / (coshu * coshu)
    snow_fr = 0.5 * (1.0 - tanhu)
    rain_fr = 1.0 - snow_fr
    P_snow = P * snow_fr
    P_rain = P * rain_fr
    aT = T - T_tr
    denom_pos = np.sqrt(((aT * aT) + (T_smeps_m * T_smeps_m)))
    posT = 0.5 * (aT + denom_pos)
    dpos_daT = 0.5 * (1.0 + (aT / denom_pos))
    dpos_dT_tr = -dpos_daT
    M_pot = f_dd * posT
    dxy = Swe - M_pot
    sqrtm = np.sqrt(((dxy * dxy) + (eps_m * eps_m)))
    M = 0.5 * ((Swe + M_pot) - sqrtm)
    dM_dSwe = 0.5 * (1.0 - (dxy / sqrtm))
    dM_dMpot = 0.5 * (1.0 + (dxy / sqrtm))
    dM_dSweu = dM_dSwe * dSwe_dSweu
    dsnow_dT_tr = (0.5 * sech2) * (1.0 / T_smeps_m)
    drain_dT_tr = -dsnow_dT_tr
    dP_snow_dT_tr = P * dsnow_dT_tr
    dP_rain_dT_tr = P * drain_dT_tr
    dMpot_dT_tr = f_dd * dpos_dT_tr
    dMpot_df_dd = posT
    dM_dT_tr = dM_dMpot * dMpot_dT_tr
    dM_df_dd = dM_dMpot * dMpot_df_dd
    Pliq = P_rain + M
    dPliq_dSweu = dM_dSweu
    dPliq_dT_tr = dP_rain_dT_tr + dM_dT_tr
    dPliq_df_dd = dM_df_dd
    E_r = smooth_min2(Pliq, Ep, eps_m)
    dEr_dPliq = dmin2_da(Pliq, Ep, eps_m)
    P_e = Pliq - E_r
    Ep_star = Ep - E_r
    dPe_dPliq = 1.0 - dEr_dPliq
    dEpstar_dPliq = -dEr_dPliq
    dPe_dSweu = dPe_dPliq * dPliq_dSweu
    dPe_dT_tr = dPe_dPliq * dPliq_dT_tr
    dPe_df_dd = dPe_dPliq * dPliq_df_dd
    dEpstar_dSweu = dEpstar_dPliq * dPliq_dSweu
    dEpstar_dT_tr = dEpstar_dPliq * dPliq_dT_tr
    dEpstar_df_dd = dEpstar_dPliq * dPliq_df_dd
    I = 0.0
    R_s = 0.0
    dI_dS = 0.0
    dI_dSmax = 0.0
    dI_dKsch = 0.0
    dI_dPe = 0.0
    dRs_dS = 0.0
    dRs_dSmax = 0.0
    dRs_dKsch = 0.0
    dRs_dPe = 0.0
    D_soil = s_max - S
    eTerm = np.exp(((-k_sch) * dT))
    eps_pe = 1e-6
    eps_D = 1e-6
    Pe_pos = smooth_pos(P_e, eps_pe)
    dPePos_dPe = dsmooth_pos_da(P_e, eps_pe)
    D_pos = smooth_pos(D_soil, eps_D)
    dDpos_dD = dsmooth_pos_da(D_soil, eps_D)
    Ic = D_pos * (1.0 - eTerm)
    Px = Pe_pos * dT
    denPI = Px + Ic
    I_amt = Px * (Ic / denPI)
    R_amt = Px - I_amt
    I = I_amt / dT
    R_s = R_amt / dT
    dI_dPx_amt = (Ic * Ic) / (denPI * denPI)
    dI_dIc_amt = (Px * Px) / (denPI * denPI)
    dIc_dDpos = 1.0 - eTerm
    dIc_dS = (dIc_dDpos * dDpos_dD) * (-1.0)
    dIc_dSmax = (dIc_dDpos * dDpos_dD) * (+1.0)
    dIc_dKsch = D_pos * (dT * eTerm)
    dI_dS = (dI_dIc_amt * dIc_dS) / dT
    dI_dSmax = (dI_dIc_amt * dIc_dSmax) / dT
    dI_dKsch = (dI_dIc_amt * dIc_dKsch) / dT
    dRs_dS = -dI_dS
    dRs_dSmax = -dI_dSmax
    dRs_dKsch = -dI_dKsch
    dI_dPe = dI_dPx_amt * dPePos_dPe
    dRs_dPe = (1.0 - dI_dPx_amt) * dPePos_dPe
    dI_dSweu = dI_dPe * dPe_dSweu
    dI_dT_tr = dI_dPe * dPe_dT_tr
    dI_df_dd = dI_dPe * dPe_df_dd
    dRs_dSweu = dRs_dPe * dPe_dSweu
    dRs_dT_tr = dRs_dPe * dPe_dT_tr
    dRs_df_dd = dRs_dPe * dPe_df_dd
    dQdt = np.zeros(L)
    for i in range(0, L):
        Qi = Q[i]
        Qshift = Q[(i + 1)] if (i < (L - 1)) else 0.0
        dQdt[i] = ((-Qi) + Qshift) + (giuh[i] * R_s)
    Q_giuh = Q[0]
    Q_perc = 0.0
    dQperc_dS = 0.0
    dQperc_dSmax = 0.0
    dQperc_dSfc = 0.0
    dQperc_da1 = 0.0
    dQperc_dkperc = 0.0
    den1 = s_max - s_fc
    above_raw = S - s_fc
    eps_ab = 1e-6
    above_pos = smooth_pos(above_raw, eps_ab)
    dAbove_dS = dsmooth_pos_da(above_raw, eps_ab)
    dAbove_dSfc = -dAbove_dS
    r1 = above_pos / den1
    eps_r1 = 1e-12
    r1e = r1 + eps_r1
    Qunc = k_perc * np.power(r1e, a1)
    dQunc_dkperc = np.power(r1e, a1)
    dQunc_da1 = Qunc * np.log(r1e)
    dQunc_dr1 = (k_perc * a1) * np.power(r1e, (a1 - 1.0))
    dr1_dS = dAbove_dS / den1
    dr1_dSmax = (-above_pos) / (den1 * den1)
    dr1_dSfc = ((dAbove_dSfc * den1) + above_pos) / (den1 * den1)
    dQunc_dS1 = dQunc_dr1 * dr1_dS
    dQunc_dSmax1 = dQunc_dr1 * dr1_dSmax
    dQunc_dSfc = dQunc_dr1 * dr1_dSfc
    eps_cap_perc = 1e-6
    Q_perc = smooth_min2(Qunc, above_pos, eps_cap_perc)
    w_unc_perc = dmin2_da(Qunc, above_pos, eps_cap_perc)
    w_cap_perc = 1.0 - w_unc_perc
    dQperc_dkperc = w_unc_perc * dQunc_dkperc
    dQperc_da1 = w_unc_perc * dQunc_da1
    dQperc_dS = (w_unc_perc * dQunc_dS1) + (w_cap_perc * dAbove_dS)
    dQperc_dSmax = w_unc_perc * dQunc_dSmax1
    dQperc_dSfc = (w_unc_perc * dQunc_dSfc) + (w_cap_perc * dAbove_dSfc)
    Q_lf = 0.0
    dQlf_dS = 0.0
    dQlf_dSmax = 0.0
    dQlf_dlfthr = 0.0
    dQlf_da2 = 0.0
    dQlf_dklf = 0.0
    eps_on = 1e-6
    eps_cap = 1e-6
    eps_r = 1e-12
    a2raw = S - lf_thr
    R_on = np.sqrt(((a2raw * a2raw) + (eps_on * eps_on)))
    above2 = 0.5 * (a2raw + R_on)
    dabove2_dS = 0.5 * (1.0 + (a2raw / R_on))
    dabove2_dlfthr = -dabove2_dS
    den2 = s_max - lf_thr
    r2 = above2 / den2
    dr2_dS = dabove2_dS / den2
    dr2_dSmax = (-above2) / (den2 * den2)
    dr2_dlfthr = ((dabove2_dlfthr * den2) + above2) / (den2 * den2)
    r2e = r2 + eps_r
    r2e_a2 = np.power(r2e, a2)
    Q_unc = k_lf * r2e_a2
    dQunc_dklf = r2e_a2
    dQunc_da2 = Q_unc * np.log(r2e)
    dQunc_dr2 = (k_lf * a2) * np.power(r2e, (a2 - 1.0))
    dQunc_dS2 = dQunc_dr2 * dr2_dS
    dQunc_dSmax2 = dQunc_dr2 * dr2_dSmax
    dQunc_dlfthr = dQunc_dr2 * dr2_dlfthr
    cap2 = above2 - Q_perc
    dcap2_dS = dabove2_dS - dQperc_dS
    dcap2_dSmax = -dQperc_dSmax
    dcap2_dlfthr = dabove2_dlfthr
    Dcap = Q_unc - cap2
    Rcap = np.sqrt(((Dcap * Dcap) + (eps_cap * eps_cap)))
    Q_lf = 0.5 * ((Q_unc + cap2) - Rcap)
    w_unc = 0.5 * (1.0 - (Dcap / Rcap))
    w_cap = 1.0 - w_unc
    dQlf_dklf = w_unc * dQunc_dklf
    dQlf_da2 = w_unc * dQunc_da2
    dQlf_dS = (w_unc * dQunc_dS2) + (w_cap * dcap2_dS)
    dQlf_dSmax = (w_unc * dQunc_dSmax2) + (w_cap * dcap2_dSmax)
    dQlf_dlfthr = (w_unc * dQunc_dlfthr) + (w_cap * dcap2_dlfthr)
    E_s = 0.0
    dEs_dS = 0.0
    dEs_dSfc = 0.0
    dEs_dSwp = 0.0
    dEs_dSweu = 0.0
    dEs_dT_tr = 0.0
    dEs_df_dd = 0.0
    dEs_dEpstar = 0.0
    if Ep_star > 0.0:
        if S <= s_wp:
            E_s = 0.0
            dEs_dEpstar = 0.0
        else:
            if S >= s_fc:
                E_s = Ep_star
                dEs_dEpstar = 1.0
            else:
                den = s_fc - s_wp
                phi = (S - s_wp) / den
                E_s = Ep_star * phi
                dEs_dS = Ep_star / den
                dEs_dSfc = Ep_star * ((-(S - s_wp)) / (den * den))
                dEs_dSwp = Ep_star * ((S - s_fc) / (den * den))
                dEs_dEpstar = phi
    dEs_dSweu = dEs_dEpstar * dEpstar_dSweu
    dEs_dT_tr = dEs_dEpstar * dEpstar_dT_tr
    dEs_df_dd = dEs_dEpstar * dEpstar_df_dd
    expterm = np.exp(((mm * G) / g_max))
    flux_exp = expterm - 1.0
    Q_gw = c_gw * flux_exp
    dQgw_dG = (c_gw * expterm) * (mm / g_max)
    dQgw_dC = flux_exp
    dQgw_dmm = (c_gw * expterm) * (G / g_max)
    dQgw_dGmax = (c_gw * expterm) * (((-mm) * G) / (g_max * g_max))
    dF = np.zeros(K)
    for k in range(0, K):
        Fk = F[k]
        if k == 0:
            dF[k] = Q_lf - (k_nsh * Fk)
        else:
            dF[k] = (k_nsh * F[(k - 1)]) - (k_nsh * Fk)
    Q_fast = k_nsh * F[(K - 1)]
    for i in range(0, m):
        udot[i] = 0.0
    udot[0] = P_snow - M
    udot[1] = ((I - Q_perc) - Q_lf) - E_s
    udot[2] = Q_perc - Q_gw
    for i in range(0, L):
        udot[(3 + i)] = dQdt[i]
    for k in range(0, K):
        udot[((3 + L) + k)] = dF[k]
    udot[(m - 1)] = (Q_fast + Q_giuh) + Q_gw
    Jx_f = np.zeros((m * m))
    Jp_f = np.zeros((m * d))
    Jth_f = np.zeros((m * d))
    iSwe = 0
    iSoil = 1
    iGW = 2
    iq1 = 3
    iF1 = 3 + L
    iQcum = m - 1
    Jx_f[((iSwe * m) + iSwe)] = -dM_dSweu
    Jx_f[((iSoil * m) + iSwe)] = dI_dSweu - dEs_dSweu
    Jx_f[((iSoil * m) + iSoil)] = ((dI_dS - dQperc_dS) - dQlf_dS) - dEs_dS
    Jx_f[((iGW * m) + iSoil)] = dQperc_dS
    Jx_f[((iGW * m) + iGW)] = -dQgw_dG
    for i in range(0, L):
        row = iq1 + i
        col = iq1 + i
        Jx_f[((row * m) + col)] = -1.0
        if i < (L - 1):
            Jx_f[((row * m) + (col + 1))] = 1.0
        Jx_f[((row * m) + iSoil)] = giuh[i] * dRs_dS
        Jx_f[((row * m) + iSwe)] = giuh[i] * dRs_dSweu
    Jx_f[((iF1 * m) + iSoil)] = dQlf_dS
    Jx_f[((iF1 * m) + iF1)] = -k_nsh
    for k in range(1, K):
        row = iF1 + k
        Jx_f[((row * m) + (row - 1))] = k_nsh
        Jx_f[((row * m) + row)] = -k_nsh
    Jx_f[((iQcum * m) + iq1)] = 1.0
    Jx_f[((iQcum * m) + ((iF1 + K) - 1))] = k_nsh
    Jx_f[((iQcum * m) + iGW)] = dQgw_dG
    pSmax = 0
    pSfc = 1
    pSwp = 2
    pKsch = 3
    pa1 = 4
    pkperc = 5
    pLfthr = 6
    pa2 = 7
    pkLf = 8
    pGmax = 9
    pCgw = 10
    pmm = 11
    pkNsh = 12
    pT_tr = 13
    pf_dd = 14
    Jp_f[(iSwe + (m * pT_tr))] = dP_snow_dT_tr - dM_dT_tr
    Jp_f[(iSwe + (m * pf_dd))] = -dM_df_dd
    Jp_f[(iSoil + (m * pSmax))] += dI_dSmax
    Jp_f[(iSoil + (m * pKsch))] += dI_dKsch
    Jp_f[(iSoil + (m * pa1))] -= dQperc_da1
    Jp_f[(iSoil + (m * pkperc))] -= dQperc_dkperc
    Jp_f[(iSoil + (m * pSmax))] -= dQperc_dSmax
    Jp_f[(iSoil + (m * pSfc))] -= dQperc_dSfc
    Jp_f[(iSoil + (m * pSmax))] -= dQlf_dSmax
    Jp_f[(iSoil + (m * pLfthr))] -= dQlf_dlfthr
    Jp_f[(iSoil + (m * pa2))] -= dQlf_da2
    Jp_f[(iSoil + (m * pkLf))] -= dQlf_dklf
    Jp_f[(iSoil + (m * pSfc))] -= dEs_dSfc
    Jp_f[(iSoil + (m * pSwp))] -= dEs_dSwp
    Jp_f[(iSoil + (m * pT_tr))] += dI_dT_tr - dEs_dT_tr
    Jp_f[(iSoil + (m * pf_dd))] += dI_df_dd - dEs_df_dd
    Jp_f[(iGW + (m * pa1))] += dQperc_da1
    Jp_f[(iGW + (m * pkperc))] += dQperc_dkperc
    Jp_f[(iGW + (m * pSmax))] += dQperc_dSmax
    Jp_f[(iGW + (m * pSfc))] += dQperc_dSfc
    Jp_f[(iGW + (m * pGmax))] -= dQgw_dGmax
    Jp_f[(iGW + (m * pCgw))] -= dQgw_dC
    Jp_f[(iGW + (m * pmm))] -= dQgw_dmm
    for i in range(0, L):
        row = iq1 + i
        Jp_f[(row + (m * pSmax))] = giuh[i] * dRs_dSmax
        Jp_f[(row + (m * pKsch))] = giuh[i] * dRs_dKsch
        Jp_f[(row + (m * pT_tr))] = giuh[i] * dRs_dT_tr
        Jp_f[(row + (m * pf_dd))] = giuh[i] * dRs_df_dd
    Jp_f[(iF1 + (m * pSmax))] += dQlf_dSmax
    Jp_f[(iF1 + (m * pLfthr))] += dQlf_dlfthr
    Jp_f[(iF1 + (m * pa2))] += dQlf_da2
    Jp_f[(iF1 + (m * pkLf))] += dQlf_dklf
    Jp_f[(iF1 + (m * pkNsh))] += -F[0]
    for k in range(1, K):
        Jp_f[((iF1 + k) + (m * pkNsh))] += F[(k - 1)] - F[k]
    Jp_f[(iQcum + (m * pGmax))] += dQgw_dGmax
    Jp_f[(iQcum + (m * pCgw))] += dQgw_dC
    Jp_f[(iQcum + (m * pmm))] += dQgw_dmm
    Jp_f[(iQcum + (m * pkNsh))] += F[(K - 1)]
    th2 = (s_fc / s_max) if (s_max != 0.0) else 0.0
    th3 = (s_wp / s_max) if (s_max != 0.0) else 0.0
    denom = s_fc - s_wp
    th7 = ((lf_thr - s_wp) / denom) if (abs(denom) > 0.0) else 0.0
    dsmax_dth1 = 1.0
    dsfc_dth1 = th2
    dsfc_dth2 = s_max
    dswp_dth1 = th3
    dswp_dth3 = s_max
    dlf_dth7 = s_fc - s_wp
    dlf_dth2 = th7 * s_max
    dlf_dth3 = (1.0 - th7) * s_max
    dlf_dth1 = ((1.0 - th7) * th3) + (th7 * th2)
    for i in range(0, m):
        Jth_f[(i + (m * 0))] = (
            (
                (Jp_f[(i + (m * pSmax))] * dsmax_dth1)
                + (Jp_f[(i + (m * pSfc))] * dsfc_dth1)
            )
            + (Jp_f[(i + (m * pSwp))] * dswp_dth1)
        ) + (Jp_f[(i + (m * pLfthr))] * dlf_dth1)
        Jth_f[(i + (m * 1))] = (Jp_f[(i + (m * pSfc))] * dsfc_dth2) + (
            Jp_f[(i + (m * pLfthr))] * dlf_dth2
        )
        Jth_f[(i + (m * 2))] = (Jp_f[(i + (m * pSwp))] * dswp_dth3) + (
            Jp_f[(i + (m * pLfthr))] * dlf_dth3
        )
        Jth_f[(i + (m * 3))] = Jp_f[(i + (m * pKsch))]
        Jth_f[(i + (m * 4))] = Jp_f[(i + (m * pa1))]
        Jth_f[(i + (m * 5))] = Jp_f[(i + (m * pkperc))]
        Jth_f[(i + (m * 6))] = Jp_f[(i + (m * pLfthr))] * dlf_dth7
        Jth_f[(i + (m * 7))] = Jp_f[(i + (m * pa2))]
        Jth_f[(i + (m * 8))] = Jp_f[(i + (m * pkLf))]
        Jth_f[(i + (m * 9))] = Jp_f[(i + (m * pGmax))]
        Jth_f[(i + (m * 10))] = Jp_f[(i + (m * pCgw))]
        Jth_f[(i + (m * 11))] = Jp_f[(i + (m * pmm))]
        Jth_f[(i + (m * 12))] = Jp_f[(i + (m * pkNsh))]
        Jth_f[(i + (m * 13))] = Jp_f[(i + (m * pT_tr))]
        Jth_f[(i + (m * 14))] = Jp_f[(i + (m * pf_dd))]
    for j in range(0, d):
        for i in range(0, m):
            dSdt[(i + (m * j))] = Jth_f[(i + (m * j))]
    for i in range(0, m):
        Jrow = Jx_f[(i * m) :]
        for k in range(0, m):
            coefficient = Jrow[k]
            if coefficient == 0.0:
                continue
            for j in range(0, d):
                dSdt[(i + (m * j))] += coefficient * Smat[(k + (m * j))]


FIELDS = [
    "s_max",
    "s_fc",
    "s_wp",
    "k_sch",
    "a1",
    "k_perc",
    "lf_thr",
    "a2",
    "k_lf",
    "g_max",
    "c_gw",
    "mm",
    "k_nsh",
    "giuh",
    "T_tr",
    "f_dd",
    "K",
    "L",
    "dT",
    "T_sm",
    "eps_m",
    "eps",
    "rho",
]


@njit(cache=True, error_model="numpy")
def augmented(z, forcing, values):
    """Evaluate state derivatives and analytic parameter sensitivities."""
    m = 14
    d = 15
    m = 4 + int(values[16]) + int(values[17])
    out = np.empty_like(z)
    cfe_nwm_odefcn(
        z[:m],
        out[:m],
        z[m:],
        out[m:],
        forcing[0],
        forcing[1],
        forcing[2],
        values[0],
        values[1],
        values[2],
        values[3],
        values[4],
        values[5],
        values[6],
        values[7],
        values[8],
        values[9],
        values[10],
        values[11],
        values[12],
        values[13],
        values[14],
        values[15],
        int(values[16]),
        int(values[17]),
        values[18],
        values[19],
        values[20],
        values[21],
        values[22],
        m,
        d,
    )
    return out
