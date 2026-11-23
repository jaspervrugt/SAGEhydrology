# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_gr4ja_equations.py
#
# Python equations and analytic sensitivities matching the shared
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Python equations and analytic sensitivities matching the shared
gr4jA core.
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
def smooth_pos(a, eps):
    """Evaluate the differentiable positive-part approximation."""
    return 0.5 * (a + np.sqrt(((a * a) + (eps * eps))))


@njit(cache=True, error_model="numpy")
def d_smooth_pos_da(a, eps):
    """Evaluate the analytic derivative of the smooth_pos response."""
    return 0.5 * (1.0 + (a / np.sqrt(((a * a) + (eps * eps)))))


@njit(cache=True, error_model="numpy")
def smooth_min(A, B, eps):
    """Evaluate the differentiable minimum approximation."""
    d = A - B
    return 0.5 * ((A + B) - np.sqrt(((d * d) + (eps * eps))))


@njit(cache=True, error_model="numpy")
def d_smooth_min_dA(A, B, eps):
    """Evaluate the analytic derivative of the smooth_min response."""
    d = A - B
    s = np.sqrt(((d * d) + (eps * eps)))
    return 0.5 * (1.0 - (d / s))


@njit(cache=True, error_model="numpy")
def d_smooth_min_dB(A, B, eps):
    """Evaluate the analytic derivative of the smooth_min response."""
    d = A - B
    s = np.sqrt(((d * d) + (eps * eps)))
    return 0.5 * (1.0 + (d / s))


@njit(cache=True, error_model="numpy")
def gr4jA_odefcn(
    u,
    udot,
    Smat,
    dSdt,
    P,
    Ep,
    T,
    x1,
    x2,
    x3,
    x4,
    x5,
    f_p,
    T_tr,
    f_dd,
    n1,
    n2,
    kappa,
    bg,
    bR,
    mts,
    T_sm,
    eps_m,
    eps_s,
    rho,
    nvar,
    m,
    d,
):
    """Evaluate the gr4jA state equations and their derivatives."""
    iSwe = 0
    iSp = 1
    iRr = 2
    iU1_1 = 3
    iU1_e = (3 + n1) - 1
    iU2_1 = 3 + n1
    iU2_e = ((3 + n1) + n2) - 1
    iQinf = m - 1
    Swe_u = u[iSwe]
    Sp_u = u[iSp]
    Rr_u = u[iRr]
    Swe = smooth_pos(Swe_u, eps_s)
    Sp = smooth_pos(Sp_u, eps_s)
    Rr = smooth_pos(Rr_u, eps_s)
    dSwe_du = d_smooth_pos_da(Swe_u, eps_s)
    dSp_du = d_smooth_pos_da(Sp_u, eps_s)
    dRr_du = d_smooth_pos_da(Rr_u, eps_s)
    U1 = np.zeros(0)
    U2 = np.zeros(0)
    dU1_du = np.zeros(0)
    dU2_du = np.zeros(0)
    if n1 > 0:
        U1 = np.zeros(n1)
        dU1_du = np.zeros(n1)
        for k in range(0, n1):
            uk = u[(iU1_1 + k)]
            U1[k] = smooth_pos(uk, eps_s)
            dU1_du[k] = d_smooth_pos_da(uk, eps_s)
    if n2 > 0:
        U2 = np.zeros(n2)
        dU2_du = np.zeros(n2)
        for k in range(0, n2):
            uk = u[(iU2_1 + k)]
            U2[k] = smooth_pos(uk, eps_s)
            dU2_du[k] = d_smooth_pos_da(uk, eps_s)
    n1_eff = n1
    n2_eff = n2
    tau_1 = x4 / n1_eff
    tau_2 = (2.0 * x4) / n2_eff
    inv_tau1 = 1.0 / tau_1
    inv_tau2 = 1.0 / tau_2
    T_smeps_m = max(T_sm, eps_m)
    utemp = (T - T_tr) / T_smeps_m
    snow_fr = 0.5 * (1.0 - np.tanh(utemp))
    rain_fr = 1.0 - snow_fr
    p_snow = P * snow_fr
    p_rain = P * rain_fr
    aT = T - T_tr
    posT = 0.5 * (aT + np.sqrt(((aT * aT) + (T_smeps_m * T_smeps_m))))
    p_pmelt = f_dd * posT
    p_amelt = smooth_min(Swe, p_pmelt, eps_m)
    p_liq = p_rain + p_amelt
    sech2 = 1.0 / (np.cosh(utemp) * np.cosh(utemp))
    dsnowFrac_dT_tr = (0.5 * sech2) / T_smeps_m
    dposT_da = 0.5 * (
        1.0 + (aT / np.sqrt(((aT * aT) + (T_smeps_m * T_smeps_m))))
    )
    dposT_dT_tr = -dposT_da
    dp_pmelt_dT_tr = f_dd * dposT_dT_tr
    dp_pmelt_df_dd = posT
    dM_dSwe = d_smooth_min_dA(Swe, p_pmelt, eps_m)
    dM_dp_pmelt = d_smooth_min_dB(Swe, p_pmelt, eps_m)
    dp_liq_dSwe = dM_dSwe
    dPsnow_dT_tr = P * dsnowFrac_dT_tr
    dPrain_dT_tr = -dPsnow_dT_tr
    dM_dT_tr = dM_dp_pmelt * dp_pmelt_dT_tr
    dM_df_dd = dM_dp_pmelt * dp_pmelt_df_dd
    dp_liq_dT_tr = dPrain_dT_tr + dM_dT_tr
    dp_liq_df_dd = dM_df_dd
    Ep_eff = f_p * Ep
    aPE = p_liq - Ep_eff
    daPE_dfp = -Ep
    Pn = smooth_pos(aPE, eps_s)
    En = smooth_pos((-aPE), eps_s)
    dPn_dp_liq = d_smooth_pos_da(aPE, eps_s)
    dEn_dp_liq = -d_smooth_pos_da((-aPE), eps_s)
    uPE = aPE / x1
    wWet = 0.5 * (1.0 + np.tanh(uPE))
    sech2PE = 1.0 / (np.cosh(uPE) * np.cosh(uPE))
    dwWet_daPE = (0.5 * sech2PE) / x1
    dwWet_dp_liq = dwWet_daPE
    dwWet_dx1 = (0.5 * sech2PE) * ((-aPE) / (x1 * x1))
    Sp_n = Sp / x1
    kappa4 = ((kappa * kappa) * kappa) * kappa
    Sp_n3 = (Sp_n * Sp_n) * Sp_n
    g = 1.0 + ((kappa4 * Sp_n3) * Sp_n)
    g_m14 = np.power(g, (-0.25))
    Pp = Sp * (1.0 - g_m14)
    dSpn_dSp = 1.0 / x1
    dSpn_dx1 = (-Sp) / (x1 * x1)
    dg_dSpn = (4.0 * kappa4) * Sp_n3
    dg_m14_dg = (-0.25) * np.power(g, (-1.25))
    dPp_dSp = (1.0 - g_m14) + (Sp * (((-dg_m14_dg) * dg_dSpn) * dSpn_dSp))
    dPp_dx1 = Sp * (((-dg_m14_dg) * dg_dSpn) * dSpn_dx1)
    dPerc_dS = dPp_dSp
    dPerc_dx1 = dPp_dx1
    Ps_wet = 0.0
    dPsWet_dS = 0.0
    dPsWet_dPn = 0.0
    dPsWet_dx1 = 0.0
    zP = Pn / x1
    aP = np.tanh(zP)
    daP_dPn = (1.0 - (aP * aP)) / x1
    Nw = (1.0 - (Sp_n * Sp_n)) * aP
    Dw = 1.0 + (Sp_n * aP)
    Ps_wet = (x1 * Nw) / Dw
    dNw_dSpn = ((-2.0) * Sp_n) * aP
    dNw_dPn = (1.0 - (Sp_n * Sp_n)) * daP_dPn
    dDw_dSpn = aP
    dDw_dPn = Sp_n * daP_dPn
    dPsWet_dSpn_loc = (x1 * ((Dw * dNw_dSpn) - (Nw * dDw_dSpn))) / (Dw * Dw)
    dPsWet_dPn = (x1 * ((Dw * dNw_dPn) - (Nw * dDw_dPn))) / (Dw * Dw)
    dPsWet_dS = dPsWet_dSpn_loc * dSpn_dSp
    dzP_dx1 = (-Pn) / (x1 * x1)
    daP_dx1 = (1.0 - (aP * aP)) * dzP_dx1
    dNw_dx1 = (dNw_dSpn * dSpn_dx1) + ((1.0 - (Sp_n * Sp_n)) * daP_dx1)
    dDw_dx1 = (dDw_dSpn * dSpn_dx1) + (Sp_n * daP_dx1)
    dNoverD_dx1 = ((Dw * dNw_dx1) - (Nw * dDw_dx1)) / (Dw * Dw)
    dPsWet_dx1 = (Nw / Dw) + (x1 * dNoverD_dx1)
    Es_dry = 0.0
    dEsDry_dS = 0.0
    dEsDry_dEn = 0.0
    dEsDry_dx1 = 0.0
    zE = En / x1
    aE = np.tanh(zE)
    daE_dEn = (1.0 - (aE * aE)) / x1
    Nd = (Sp * (2.0 - Sp_n)) * aE
    Dd = 1.0 + ((1.0 - Sp_n) * aE)
    invD = 1.0 / Dd
    invD2 = invD * invD
    Es_dry = Nd * invD
    dNd_dSp = (2.0 - Sp_n) * aE
    dNd_dSpn = (-Sp) * aE
    dNd_daE = Sp * (2.0 - Sp_n)
    dDd_dSpn = -aE
    dDd_daE = 1.0 - Sp_n
    dEs_dSp = dNd_dSp * invD
    dEs_dSpn = ((dNd_dSpn * Dd) - (Nd * dDd_dSpn)) * invD2
    dEs_daE = ((dNd_daE * Dd) - (Nd * dDd_daE)) * invD2
    dEsDry_dS = dEs_dSp + (dEs_dSpn * dSpn_dSp)
    dEsDry_dEn = dEs_daE * daE_dEn
    dzE_dx1 = (-En) / (x1 * x1)
    daE_dx1 = (1.0 - (aE * aE)) * dzE_dx1
    dEsDry_dx1 = (dEs_dSpn * dSpn_dx1) + (dEs_daE * daE_dx1)
    Ps = wWet * Ps_wet
    Es = (1.0 - wWet) * Es_dry
    dPs_dS = wWet * dPsWet_dS
    dEs_dS = (1.0 - wWet) * dEsDry_dS
    dPsWet_dp_liq = dPsWet_dPn * dPn_dp_liq
    dEsDry_dp_liq = dEsDry_dEn * dEn_dp_liq
    dPs_dp_liq = (wWet * dPsWet_dp_liq) + (Ps_wet * dwWet_dp_liq)
    dEs_dp_liq = ((1.0 - wWet) * dEsDry_dp_liq) - (Es_dry * dwWet_dp_liq)
    dPs_dx1 = (wWet * dPsWet_dx1) + (Ps_wet * dwWet_dx1)
    dEs_dx1 = ((1.0 - wWet) * dEsDry_dx1) - (Es_dry * dwWet_dx1)
    P1 = Pn - Ps
    dP1_dS = -dPs_dS
    dP1_daPE = dPn_dp_liq - dPs_dp_liq
    dP1_dx1 = -dPs_dx1
    dP1_dfp = (dPn_dp_liq - dPs_dp_liq) * daPE_dfp
    dP1_dT_tr = dP1_daPE * dp_liq_dT_tr
    dP1_df_dd = dP1_daPE * dp_liq_df_dd
    Pr = P1 + Pp
    dPr_dS = dP1_dS + dPerc_dS
    dPr_dp_liq = dP1_daPE
    dPr_dx1 = dP1_dx1 + dPerc_dx1
    dPr_dfp = dP1_dfp
    dPr_dT_tr = dP1_dT_tr
    dPr_df_dd = dP1_df_dd
    in1 = x5 * Pr
    in2 = (1.0 - x5) * Pr
    q_1 = (U1[(n1 - 1)] * inv_tau1) if (n1 > 0) else 0.0
    q_2 = (U2[(n2 - 1)] * inv_tau2) if (n2 > 0) else 0.0
    Rratio = Rr / x3
    Rratio2 = Rratio * Rratio
    Rratio3 = Rratio2 * Rratio
    sqrtRratio = np.sqrt(Rratio)
    q_g = (x2 * Rratio3) * sqrtRratio
    q_r = ((x3 / ((bR - 1.0) * mts)) * Rratio3) * Rratio2
    q_d = q_2 + q_g
    q_dpos = smooth_pos(q_d, eps_s)
    q_out = q_r + q_dpos
    dqg_dRr = (((x2 * bg) * Rratio2) * sqrtRratio) * (1.0 / x3)
    dqg_dx2 = Rratio3 * sqrtRratio
    dqg_dx3 = (((x2 * bg) * Rratio2) * sqrtRratio) * ((-Rr) / (x3 * x3))
    c0 = 1.0 / ((bR - 1.0) * mts)
    Rr4 = ((Rr * Rr) * Rr) * Rr
    x3_4 = ((x3 * x3) * x3) * x3
    dqr_dRr = (((c0 * bR) * Rr4) * 1) / x3_4
    dqr_dx3 = ((((c0 * Rr) * Rr4) * (1.0 - bR)) * 1) / (x3_4 * x3)
    udot[iSwe] = p_snow - p_amelt
    udot[iSp] = (Ps - Es) - Pp
    udot[iRr] = (q_1 + q_g) - q_r
    if n1 > 0:
        udot[iU1_1] = in1 - (U1[0] * inv_tau1)
        for k in range(1, n1):
            idx = iU1_1 + k
            udot[idx] = (U1[(k - 1)] * inv_tau1) - (U1[k] * inv_tau1)
    if n2 > 0:
        udot[iU2_1] = in2 - (U2[0] * inv_tau2)
        for k in range(1, n2):
            idx = iU2_1 + k
            udot[idx] = (U2[(k - 1)] * inv_tau2) - (U2[k] * inv_tau2)
    udot[iQinf] = q_out
    Jx_f_ = np.zeros((m * m))
    Jth_f_ = np.zeros((m * d))
    Jx_f_[(iSwe) * m + (iSwe)] = -dM_dSwe
    Jx_f_[(iSp) * m + (iSwe)] = (dPs_dp_liq - dEs_dp_liq) * dp_liq_dSwe
    Jx_f_[(iSp) * m + (iSp)] = (dPs_dS - dEs_dS) - dPerc_dS
    din1_dSwe = (x5 * dPr_dp_liq) * dp_liq_dSwe
    din1_dS = x5 * dPr_dS
    din2_dSwe = ((1.0 - x5) * dPr_dp_liq) * dp_liq_dSwe
    din2_dS = (1.0 - x5) * dPr_dS
    if n1 > 0:
        Jx_f_[(iU1_1) * m + (iSwe)] = din1_dSwe
        Jx_f_[(iU1_1) * m + (iSp)] = din1_dS
        Jx_f_[(iU1_1) * m + (iU1_1)] = -inv_tau1
        for k in range(1, n1):
            idx = iU1_1 + k
            Jx_f_[(idx) * m + ((idx - 1))] = inv_tau1
            Jx_f_[(idx) * m + (idx)] = -inv_tau1
        Jx_f_[(iRr) * m + (iU1_e)] = inv_tau1
    if n2 > 0:
        Jx_f_[(iU2_1) * m + (iSwe)] = din2_dSwe
        Jx_f_[(iU2_1) * m + (iSp)] = din2_dS
        Jx_f_[(iU2_1) * m + (iU2_1)] = -inv_tau2
        for k in range(1, n2):
            idx = iU2_1 + k
            Jx_f_[(idx) * m + ((idx - 1))] = inv_tau2
            Jx_f_[(idx) * m + (idx)] = -inv_tau2
    Jx_f_[(iRr) * m + (iRr)] = dqg_dRr - dqr_dRr
    dQdpos_dQd = d_smooth_pos_da(q_d, eps_s)
    Jx_f_[(iQinf) * m + (iRr)] = dqr_dRr + (dQdpos_dQd * dqg_dRr)
    if n2 > 0:
        Jx_f_[(iQinf) * m + (iU2_e)] = dQdpos_dQd * inv_tau2
    for ii in range(0, m):
        Jx_f_[(ii) * m + (iSwe)] *= dSwe_du
        Jx_f_[(ii) * m + (iSp)] *= dSp_du
        Jx_f_[(ii) * m + (iRr)] *= dRr_du
    if n1 > 0:
        for k in range(0, n1):
            col = iU1_1 + k
            for ii in range(0, m):
                Jx_f_[(ii) * m + (col)] *= dU1_du[k]
    if n2 > 0:
        for k in range(0, n2):
            col = iU2_1 + k
            for ii in range(0, m):
                Jx_f_[(ii) * m + (col)] *= dU2_du[k]
    dPs_dT_tr = dPs_dp_liq * dp_liq_dT_tr
    dPs_df_dd = dPs_dp_liq * dp_liq_df_dd
    dEs_dT_tr = dEs_dp_liq * dp_liq_dT_tr
    dEs_df_dd = dEs_dp_liq * dp_liq_df_dd
    Jth_f_[(iSwe) + m * (6)] = dPsnow_dT_tr - dM_dT_tr
    Jth_f_[(iSwe) + m * (7)] = -dM_df_dd
    Jth_f_[(iSp) + m * (6)] = dPs_dT_tr - dEs_dT_tr
    Jth_f_[(iSp) + m * (7)] = dPs_df_dd - dEs_df_dd
    Jth_f_[(iSp) + m * (0)] = (dPs_dx1 - dEs_dx1) - dPerc_dx1
    dPs_dfp = dPs_dp_liq * daPE_dfp
    dEs_dfp = dEs_dp_liq * daPE_dfp
    Jth_f_[(iSp) + m * (5)] = dPs_dfp - dEs_dfp
    if n1 > 0:
        Jth_f_[(iU1_1) + m * (5)] = x5 * dPr_dfp
    if n2 > 0:
        Jth_f_[(iU2_1) + m * (5)] = (1.0 - x5) * dPr_dfp
    Jth_f_[(iRr) + m * (1)] = dqg_dx2
    Jth_f_[(iQinf) + m * (1)] = dQdpos_dQd * dqg_dx2
    Jth_f_[(iRr) + m * (2)] = dqg_dx3 - dqr_dx3
    Jth_f_[(iQinf) + m * (2)] = dqr_dx3 + (dQdpos_dQd * dqg_dx3)
    dinv_tau1_dx4 = (-(1.0 / (tau_1 * tau_1))) * (1.0 / n1_eff)
    dinv_tau2_dx4 = (-(1.0 / (tau_2 * tau_2))) * (2.0 / n2_eff)
    if n1 > 0:
        Jth_f_[(iU1_1) + m * (3)] = (-U1[0]) * dinv_tau1_dx4
        for k in range(1, n1):
            idx = iU1_1 + k
            Jth_f_[(idx) + m * (3)] = (U1[(k - 1)] - U1[k]) * dinv_tau1_dx4
        Jth_f_[(iRr) + m * (3)] = U1[(n1 - 1)] * dinv_tau1_dx4
    if n2 > 0:
        Jth_f_[(iU2_1) + m * (3)] = (-U2[0]) * dinv_tau2_dx4
        for k in range(1, n2):
            idx = iU2_1 + k
            Jth_f_[(idx) + m * (3)] = (U2[(k - 1)] - U2[k]) * dinv_tau2_dx4
        Jth_f_[(iQinf) + m * (3)] = dQdpos_dQd * (
            U2[(n2 - 1)] * dinv_tau2_dx4
        )
    if n1 > 0:
        Jth_f_[(iU1_1) + m * (4)] = Pr
    if n2 > 0:
        Jth_f_[(iU2_1) + m * (4)] = -Pr
    if n1 > 0:
        Jth_f_[(iU1_1) + m * (0)] = x5 * dPr_dx1
        Jth_f_[(iU1_1) + m * (6)] = x5 * dPr_dT_tr
        Jth_f_[(iU1_1) + m * (7)] = x5 * dPr_df_dd
    if n2 > 0:
        Jth_f_[(iU2_1) + m * (0)] = (1.0 - x5) * dPr_dx1
        Jth_f_[(iU2_1) + m * (6)] = (1.0 - x5) * dPr_dT_tr
        Jth_f_[(iU2_1) + m * (7)] = (1.0 - x5) * dPr_df_dd
    for j in range(0, d):
        for i in range(0, m):
            dSdt[(i + (m * j))] = Jth_f_[(i) + m * (j)]
    for i in range(0, m):
        for k in range(0, m):
            coefficient = Jx_f_[(i) * m + (k)]
            if coefficient == 0.0:
                continue
            for j in range(0, d):
                ij = i + (m * j)
                kj = k + (m * j)
                dSdt[ij] += coefficient * Smat[kj]


FIELDS = [
    "x1",
    "x2",
    "x3",
    "x4",
    "x5",
    "f_p",
    "T_tr",
    "f_dd",
    "n1",
    "n2",
    "kappa",
    "bg",
    "bR",
    "mts",
    "T_sm",
    "eps_m",
    "eps_s",
    "rho",
]


@njit(cache=True, error_model="numpy")
def augmented(z, forcing, values):
    """Evaluate state derivatives and analytic parameter sensitivities."""
    m = 10
    d = 8
    out = np.empty_like(z)
    gr4jA_odefcn(
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
        int(values[8]),
        int(values[9]),
        values[10],
        values[11],
        values[12],
        values[13],
        values[14],
        values[15],
        values[16],
        values[17],
        len(z),
        m,
        d,
    )
    return out
