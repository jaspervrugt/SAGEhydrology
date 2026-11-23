# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_xinanjiang_equations.py
#
# Python equations and analytic sensitivities matching the shared
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Python equations and analytic sensitivities matching the shared
Xinanjiang core.
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
def lf_func(S, Smax, rho, eps):
    """Evaluate the nonlinear storage-dependent flux."""
    num = S - (Smax - ((rho * Smax) * eps))
    den = rho * Smax
    z = -(num / den)
    ez = np.exp(z)
    return 1.0 / (1.0 + ez)


@njit(cache=True, error_model="numpy")
def dlf_dS(lf_val, Smax, rho):
    """Differentiate the nonlinear flux with respect to storage."""
    return (lf_val * (1.0 - lf_val)) / (rho * Smax)


@njit(cache=True, error_model="numpy")
def smooth_pos(a, eps):
    """Evaluate the differentiable positive-part approximation."""
    return 0.5 * (a + np.sqrt(((a * a) + (eps * eps))))


@njit(cache=True, error_model="numpy")
def d_smooth_pos_da(a, eps):
    """Evaluate the analytic derivative of the smooth_pos response."""
    denom = np.sqrt(((a * a) + (eps * eps)))
    return 0.5 * (1.0 + (a / denom))


@njit(cache=True, error_model="numpy")
def xinanjiang_odefcn(
    u,
    udot,
    Smat,
    dSdt,
    P,
    Ep,
    T,
    f_p,
    A_im,
    a,
    b,
    W_max,
    LM,
    c,
    S_max,
    S_tot,
    f_wm,
    f_lm,
    Ex,
    k_i,
    k_g,
    c_i,
    c_g,
    k_f,
    T_tr,
    f_dd,
    T_sm,
    eps_m,
    eps_r,
    eps,
    rho,
    m,
    d,
):
    """Evaluate the xinanjiang state equations and their derivatives."""
    Swe_u = u[0]
    W_u = u[1]
    Sw_u = u[2]
    Si = u[3]
    Sg = u[4]
    Sf1 = u[5]
    Sf2 = u[6]
    Sf3 = u[7]
    Swe = smooth_pos(Swe_u, eps_m)
    W = smooth_pos(W_u, eps_m)
    Sw = smooth_pos(Sw_u, eps_m)
    dSwe_dx1 = d_smooth_pos_da(Swe_u, eps_m)
    dW_dx2 = d_smooth_pos_da(W_u, eps_m)
    dSw_dx3 = d_smooth_pos_da(Sw_u, eps_m)
    T_smeps_m = max(T_sm, eps_m)
    utemp = (T - T_tr) / T_smeps_m
    snow_fr = 0.5 * (1.0 - np.tanh(utemp))
    rain_fr = 1.0 - snow_fr
    P_snow = P * snow_fr
    P_rain = P * rain_fr
    aT = T - T_tr
    posT = 0.5 * (aT + np.sqrt(((aT * aT) + (T_smeps_m * T_smeps_m))))
    M_pot = f_dd * posT
    dxy = Swe - M_pot
    sqrtm = np.sqrt(((dxy * dxy) + (eps_m * eps_m)))
    M = 0.5 * ((Swe + M_pot) - sqrtm)
    Pliq = P_rain + M
    coshu = np.cosh(utemp)
    sech2 = 1.0 / (coshu * coshu)
    dsnow_dT_tr = (0.5 * sech2) / T_smeps_m
    denom_pos = np.sqrt(((aT * aT) + (T_smeps_m * T_smeps_m)))
    dposT_da = 0.5 * (1.0 + (aT / denom_pos))
    dposT_dT_tr = -dposT_da
    dMpot_dT_tr = f_dd * dposT_dT_tr
    dMpot_f_dd = posT
    dM_dSwe = 0.5 * (1.0 - (dxy / sqrtm))
    dM_dMpot = 0.5 * (1.0 + (dxy / sqrtm))
    dM_dT_tr = dM_dMpot * dMpot_dT_tr
    dM_f_dd = dM_dMpot * dMpot_f_dd
    dPliq_dSwe = dM_dSwe
    dPliq_dT_tr = ((-P) * dsnow_dT_tr) + dM_dT_tr
    dPliq_f_dd = dM_f_dd
    Pi = (1.0 - A_im) * Pliq
    Rb = A_im * Pliq
    Ea = f_p * Ep
    dPi_dPliq = 1.0 - A_im
    dRb_dPliq = A_im
    W_max2 = W_max * W_max
    uW = W / W_max
    duW_dW = 1.0 / W_max
    duW_dWmax = (-W) / W_max2
    g = uW - (0.5 - a)
    gt = g / eps_r
    w = 0.5 * (1.0 + np.tanh(gt))
    sech22 = 1.0 / (np.cosh(gt) * np.cosh(gt))
    dw_dg = (0.5 * sech22) / eps_r
    dg_dW = duW_dW
    dg_da = 1.0
    dg_dWmax = duW_dWmax
    base1 = 0.5 - a
    base1e = max(base1, 1e-12)
    c1 = np.power(base1e, (1.0 - b))
    dc1_da = (-(1.0 - b)) * np.power(base1e, (-b))
    dc1_db = (-c1) * np.log(base1e)
    uWe = smooth_pos(uW, eps_m)
    duWe_duW = d_smooth_pos_da(uW, eps_m)
    uWb = np.power(uWe, b)
    uWbm1 = np.power(uWe, (b - 1.0))
    loguW = np.log(uWe)
    Rlow = (Pi * c1) * uWb
    dRlow_dPi = c1 * uWb
    dRlow_dW = ((((Pi * c1) * b) * uWbm1) * duWe_duW) * duW_dW
    dRlow_da = (Pi * dc1_da) * uWb
    dRlow_db = Pi * ((dc1_db * uWb) + ((c1 * uWb) * loguW))
    dRlow_dWmax = ((((Pi * c1) * b) * uWbm1) * duWe_duW) * duW_dWmax
    base2 = 0.5 + a
    base2e = max(base2, 1e-12)
    c2 = np.power(base2e, (1.0 - b))
    dc2_da = (1.0 - b) * np.power(base2e, (-b))
    dc2_db = (-c2) * np.log(base2e)
    v_raw = 1.0 - uW
    v = smooth_pos(v_raw, eps_m)
    dv_dvraw = d_smooth_pos_da(v_raw, eps_m)
    dvraw_dW = -duW_dW
    dvraw_dWmax = -duW_dWmax
    dv_dW = dv_dvraw * dvraw_dW
    dv_dWmax = dv_dvraw * dvraw_dWmax
    vb = np.power(v, b)
    vbm1 = np.power(v, (b - 1.0))
    logv = np.log(v)
    Rup = Pi * (1.0 - (c2 * vb))
    dRup_dPi = 1.0 - (c2 * vb)
    dRup_dW = (Pi * (((-c2) * b) * vbm1)) * dv_dW
    dRup_da = Pi * ((-dc2_da) * vb)
    dRup_db = Pi * (-((dc2_db * vb) + ((c2 * vb) * logv)))
    dRup_dWmax = (Pi * (((-c2) * b) * vbm1)) * dv_dWmax
    R = ((1.0 - w) * Rlow) + (w * Rup)
    dR_dPi = ((1.0 - w) * dRlow_dPi) + (w * dRup_dPi)
    dR_dW = (((1.0 - w) * dRlow_dW) + (w * dRup_dW)) + (
        (dw_dg * dg_dW) * (Rup - Rlow)
    )
    dR_da = (((1.0 - w) * dRlow_da) + (w * dRup_da)) + (
        (dw_dg * dg_da) * (Rup - Rlow)
    )
    dR_db = ((1.0 - w) * dRlow_db) + (w * dRup_db)
    dR_dWmax = (((1.0 - w) * dRlow_dWmax) + (w * dRup_dWmax)) + (
        (dw_dg * dg_dWmax) * (Rup - Rlow)
    )
    dR_dPliq = dR_dPi * dPi_dPliq
    E = 0.0
    dE_df_p = 0.0
    dE_dLM = 0.0
    dE_dc = 0.0
    dE_dW = 0.0
    # Evaporation uses actual available water, not its positive smoothing
    # offset. Normalize the former 1-mm limiter to match c*Ea at c*LM.
    We = max(W_u, 0.0)
    if We > LM:
        E = Ea
        dE_df_p = Ep
        dE_dW = 0.0
    else:
        if We >= (c * LM):
            E = (We / LM) * Ea
            dE_df_p = (We / LM) * Ep
            dE_dLM = (-(We / (LM * LM))) * Ea
            dE_dc = 0.0
            dE_dW = Ea / LM
        else:
            availability = We / (We + 1.0)
            factor = c + 1.0 / LM
            E = Ea * availability * factor
            dE_dW = Ea * factor / ((We + 1.0) ** 2) if W_u > 0 else 0.0
            dE_dc = Ea * availability
            dE_df_p = Ep * availability * factor
            dE_dLM = -Ea * availability / (LM * LM)
    z = 1.0 - (Sw / S_max)
    if z < 0.0:
        z = 0.0
    dum = 1.0 - np.power(z, Ex)
    Rs = R * dum
    Ri = (k_i * Sw) * dum
    Rg = (k_g * Sw) * dum
    Qi = c_i * Si
    Qg = c_g * Sg
    Qs = Rs + Rb
    qf1 = k_f * Sf1
    qf2 = k_f * Sf2
    qf3 = k_f * Sf3
    ddum_dS = 0.0
    ddum_dSmax = 0.0
    ddum_dEx = 0.0
    g1 = sigmoid((z / eps_r))
    g2 = sigmoid(((S_max - Sw) / eps_r))
    g12 = g1 * g2
    ddS_base = 0.0
    ddSmax_base = 0.0
    ddEx_base = 0.0
    zlog = max(z, eps_m)
    z_pow_Exm1 = np.power(zlog, (Ex - 1.0))
    z_pow_Ex = z_pow_Exm1 * zlog
    ddS_base = (Ex * z_pow_Exm1) * (1.0 / S_max)
    ddSmax_base = ((-Ex) * z_pow_Exm1) * (Sw / (S_max * S_max))
    ddEx_base = (-z_pow_Ex) * np.log(zlog)
    ddum_dS = g12 * ddS_base
    ddum_dSmax = g12 * ddSmax_base
    ddum_dEx = g12 * ddEx_base
    udot[0] = P_snow - M
    udot[1] = (Pi - E) - R
    udot[2] = ((R - Rs) - Ri) - Rg
    udot[3] = Ri - Qi
    udot[4] = Rg - Qg
    udot[5] = ((Qs + Qi) + Qg) - qf1
    udot[6] = qf1 - qf2
    udot[7] = qf2 - qf3
    udot[8] = qf3
    dRs_dW = dR_dW * dum
    dRs_dS = R * ddum_dS
    dRi_dS = k_i * (dum + (Sw * ddum_dS))
    dRg_dS = k_g * (dum + (Sw * ddum_dS))
    dQi_dSi = c_i
    dQg_dSg = c_g
    dQs_dW = dRs_dW
    dQs_dS = dRs_dS
    dPi_dSwe = dPi_dPliq * dPliq_dSwe
    dRb_dSwe = dRb_dPliq * dPliq_dSwe
    dR_dSwe = dR_dPliq * dPliq_dSwe
    dRs_dSwe = dR_dSwe * dum
    dQs_dSwe = dRs_dSwe + dRb_dSwe
    Jx_f = np.zeros(
        (
            9,
            9,
        )
    )
    Jth_f = np.zeros(
        (
            9,
            16,
        )
    )
    J_Wmax = np.zeros((9,))
    J_LM = np.zeros((9,))
    J_Smax = np.zeros((9,))
    Jx_f[0][0] = -dM_dSwe
    Jx_f[1][0] = dPi_dSwe - dR_dSwe
    Jx_f[1][1] = -dR_dW
    Jx_f[2][0] = dR_dSwe - dRs_dSwe
    Jx_f[2][1] = dR_dW - dRs_dW
    Jx_f[2][2] = ((-dRs_dS) - dRi_dS) - dRg_dS
    Jx_f[3][2] = dRi_dS
    Jx_f[3][3] = -dQi_dSi
    Jx_f[4][2] = dRg_dS
    Jx_f[4][4] = -dQg_dSg
    Jx_f[5][0] = dQs_dSwe
    Jx_f[5][1] = dQs_dW
    Jx_f[5][2] = dQs_dS
    Jx_f[5][3] = dQi_dSi
    Jx_f[5][4] = dQg_dSg
    Jx_f[5][5] = -k_f
    Jx_f[6][5] = k_f
    Jx_f[6][6] = -k_f
    Jx_f[7][6] = k_f
    Jx_f[7][7] = -k_f
    Jx_f[8][7] = k_f
    for ii in range(0, m):
        Jx_f[ii][0] *= dSwe_dx1
        Jx_f[ii][1] *= dW_dx2
        Jx_f[ii][2] *= dSw_dx3
    # dE_dW already differentiates the raw water state.
    Jx_f[1][1] -= dE_dW
    dRs_dSmax = R * ddum_dSmax
    dRi_dSmax = (k_i * Sw) * ddum_dSmax
    dRg_dSmax = (k_g * Sw) * ddum_dSmax
    dRs_dEx = R * ddum_dEx
    dRi_dEx = (k_i * Sw) * ddum_dEx
    dRg_dEx = (k_g * Sw) * ddum_dEx
    j = 0
    Jth_f[1][j] = -dE_df_p
    j = 1
    dPi_dAim = -Pliq
    dRb_dAim = Pliq
    dR_dAim = dR_dPi * dPi_dAim
    dRs_dAim = dR_dAim * dum
    dQs_dAim = dRs_dAim + dRb_dAim
    Jth_f[1][j] = dPi_dAim - dR_dAim
    Jth_f[2][j] = dR_dAim - dRs_dAim
    Jth_f[5][j] = dQs_dAim
    j = 2
    dRs_da = dR_da * dum
    Jth_f[1][j] = -dR_da
    Jth_f[2][j] = dR_da - dRs_da
    Jth_f[5][j] = dRs_da
    j = 3
    dRs_db = dR_db * dum
    Jth_f[1][j] = -dR_db
    Jth_f[2][j] = dR_db - dRs_db
    Jth_f[5][j] = dRs_db
    dRs_dWmax = dR_dWmax * dum
    J_Wmax[1] += -dR_dWmax
    J_Wmax[2] += dR_dWmax - dRs_dWmax
    J_Wmax[5] += dRs_dWmax
    J_LM[1] += -dE_dLM
    j = 6
    Jth_f[1][j] = -dE_dc
    J_Smax[2] += ((-dRs_dSmax) - dRi_dSmax) - dRg_dSmax
    J_Smax[3] += dRi_dSmax
    J_Smax[4] += dRg_dSmax
    J_Smax[5] = dRs_dSmax
    j = 8
    Jth_f[2][j] = ((-dRs_dEx) - dRi_dEx) - dRg_dEx
    Jth_f[3][j] = dRi_dEx
    Jth_f[4][j] = dRg_dEx
    Jth_f[5][j] = dRs_dEx
    j = 9
    dRi_dki = Sw * dum
    Jth_f[2][j] = -dRi_dki
    Jth_f[3][j] = dRi_dki
    j = 10
    dRg_dkg = Sw * dum
    Jth_f[2][j] = -dRg_dkg
    Jth_f[4][j] = dRg_dkg
    j = 11
    Jth_f[3][j] = -Si
    Jth_f[5][j] = Si
    j = 12
    Jth_f[4][j] = -Sg
    Jth_f[5][j] = Sg
    j = 13
    Jth_f[5][j] = -Sf1
    Jth_f[6][j] = Sf1 - Sf2
    Jth_f[7][j] = Sf2 - Sf3
    Jth_f[8][j] = Sf3
    j_T_tr = d - 2
    j_f_dd = d - 1
    dPi_dT_tr = dPi_dPliq * dPliq_dT_tr
    dPi_f_dd = dPi_dPliq * dPliq_f_dd
    dRb_dT_tr = dRb_dPliq * dPliq_dT_tr
    dRb_f_dd = dRb_dPliq * dPliq_f_dd
    dR_dT_tr = dR_dPliq * dPliq_dT_tr
    dR_f_dd = dR_dPliq * dPliq_f_dd
    dRs_dT_tr = dR_dT_tr * dum
    dRs_f_dd = dR_f_dd * dum
    dQs_dT_tr = dRs_dT_tr + dRb_dT_tr
    dQs_f_dd = dRs_f_dd + dRb_f_dd
    Jth_f[0][j_T_tr] = (P * dsnow_dT_tr) - dM_dT_tr
    Jth_f[0][j_f_dd] = -dM_f_dd
    Jth_f[1][j_T_tr] = dPi_dT_tr - dR_dT_tr
    Jth_f[1][j_f_dd] = dPi_f_dd - dR_f_dd
    Jth_f[2][j_T_tr] = dR_dT_tr - dRs_dT_tr
    Jth_f[2][j_f_dd] = dR_f_dd - dRs_f_dd
    Jth_f[5][j_T_tr] = dQs_dT_tr
    Jth_f[5][j_f_dd] = dQs_f_dd
    id_f_wm = 4
    id_f_lm = 5
    id_S_tot = 7
    for i in range(0, m):
        Jth_f[i][id_f_wm] += (
            (S_tot * J_Wmax[i]) + ((f_lm * S_tot) * J_LM[i])
        ) - (S_tot * J_Smax[i])
        Jth_f[i][id_f_lm] += W_max * J_LM[i]
        Jth_f[i][id_S_tot] += (
            (f_wm * J_Wmax[i]) + ((f_lm * f_wm) * J_LM[i])
        ) + ((1.0 - f_wm) * J_Smax[i])
    for j in range(0, d):
        for i in range(0, m):
            dSdt[(i + (m * j))] = Jth_f[i][j]
    for i in range(0, m):
        for k in range(0, m):
            coefficient = Jx_f[i][k]
            if coefficient == 0.0:
                continue
            for j in range(0, d):
                dSdt[(i + (m * j))] += coefficient * Smat[(k + (m * j))]


FIELDS = [
    "f_p",
    "A_im",
    "a",
    "b",
    "W_max",
    "LM",
    "c",
    "S_max",
    "S_tot",
    "f_wm",
    "f_lm",
    "Ex",
    "k_i",
    "k_g",
    "c_i",
    "c_g",
    "k_f",
    "T_tr",
    "f_dd",
    "T_sm",
    "eps_m",
    "eps_r",
    "eps",
    "rho",
]


@njit(cache=True, error_model="numpy")
def augmented(z, forcing, values):
    """Evaluate state derivatives and analytic parameter sensitivities."""
    m = 9
    d = 16
    out = np.empty_like(z)
    xinanjiang_odefcn(
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
        values[16],
        values[17],
        values[18],
        values[19],
        values[20],
        values[21],
        values[22],
        values[23],
        m,
        d,
    )
    return out
