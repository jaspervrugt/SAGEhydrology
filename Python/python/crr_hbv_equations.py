# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_hbv_equations.py
#
# Python equations and analytic sensitivities matching the shared hbv core.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Python equations and analytic sensitivities matching the shared hbv core.
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
    return 0.5 * (1.0 - (d / np.sqrt(((d * d) + (eps * eps)))))


@njit(cache=True, error_model="numpy")
def d_smooth_min_dB(A, B, eps):
    """Evaluate the analytic derivative of the smooth_min response."""
    d = A - B
    return 0.5 * (1.0 + (d / np.sqrt(((d * d) + (eps * eps)))))


@njit(cache=True, error_model="numpy")
def smooth_max(A, B, eps):
    """Evaluate the differentiable maximum approximation."""
    d = A - B
    return 0.5 * ((A + B) + np.sqrt(((d * d) + (eps * eps))))


@njit(cache=True, error_model="numpy")
def d_smooth_max_dA(A, B, eps):
    """Evaluate the analytic derivative of the smooth_max response."""
    d = A - B
    return 0.5 * (1.0 + (d / np.sqrt(((d * d) + (eps * eps)))))


@njit(cache=True, error_model="numpy")
def clamp01(z, eps):
    """Apply the smooth unit-interval clipping function."""
    return smooth_min(smooth_max(z, 0.0, eps), 1.0, eps)


@njit(cache=True, error_model="numpy")
def dclamp01_dz(z, eps):
    """Evaluate the analytic derivative of the clamp01 response."""
    A = smooth_max(z, 0.0, eps)
    dmin_dA = d_smooth_min_dA(A, 1.0, eps)
    dmax_dz = d_smooth_max_dA(z, 0.0, eps)
    return dmin_dA * dmax_dz


@njit(cache=True, error_model="numpy")
def phiE_fun(ratio, eps):
    """Evaluate the storage-dependent evaporation response."""
    return smooth_min(ratio, 1.0, eps)


@njit(cache=True, error_model="numpy")
def dphiE_dratio(ratio, eps):
    """Differentiate the evaporation response by its storage ratio."""
    return d_smooth_min_dA(ratio, 1.0, eps)


@njit(cache=True, error_model="numpy")
def sech2(x):
    """Evaluate the squared hyperbolic secant."""
    ax = abs(x)
    if ax > 350.0:
        return 0.0
    e = np.exp(((-2.0) * ax))
    return (4.0 * e) / ((1.0 + e) * (1.0 + e))


@njit(cache=True, error_model="numpy")
def hbv_odefcn(
    u,
    udot,
    Smat,
    dSdt,
    P,
    Ep,
    T,
    f_c,
    beta,
    lp,
    k_0,
    uzl,
    k_1,
    k_2,
    perc,
    T_tr,
    f_dd,
    sfcf,
    cfr,
    eps_s,
    eps_t,
    eps_x,
    m,
    d,
):
    """Evaluate the hbv state equations and their derivatives."""
    Swe_raw = u[0]
    Sm_raw = u[1]
    Uz_raw = u[2]
    Lz_raw = u[3]
    Swe = smooth_pos(Swe_raw, eps_s)
    Sm = smooth_pos(Sm_raw, eps_s)
    Uz = smooth_pos(Uz_raw, eps_s)
    Lz = smooth_pos(Lz_raw, eps_s)
    dSwe_dSweRaw = d_smooth_pos_da(Swe_raw, eps_s)
    dSm_dSmRaw = d_smooth_pos_da(Sm_raw, eps_s)
    dUz_dUzRaw = d_smooth_pos_da(Uz_raw, eps_s)
    dLz_dLzRaw = d_smooth_pos_da(Lz_raw, eps_s)
    uT = (T - T_tr) / eps_t
    snow_fr = 0.5 * (1.0 - np.tanh(uT))
    rain_fr = 1.0 - snow_fr
    dsnow_dT_tr = (0.5 * sech2(uT)) / eps_t
    Ps = (sfcf * P) * snow_fr
    Pr = P * rain_fr
    dPs_dsfcf = P * snow_fr
    dPs_dT_tr = (sfcf * P) * dsnow_dT_tr
    dPr_dT_tr = (-P) * dsnow_dT_tr
    Tm = T - T_tr
    posTm = smooth_pos(Tm, eps_t)
    negTm = smooth_pos((-Tm), eps_t)
    dposTm_dTm = d_smooth_pos_da(Tm, eps_t)
    dnegTm_dTm = -d_smooth_pos_da((-Tm), eps_t)
    Mpot = f_dd * posTm
    Rpot = (cfr * f_dd) * negTm
    dMpot_dT_tr = (-f_dd) * dposTm_dTm
    dMpot_df_dd = posTm
    dRpot_dT_tr = (-(cfr * f_dd)) * dnegTm_dTm
    dRpot_df_dd = cfr * negTm
    dRpot_dcfr = f_dd * negTm
    M = smooth_min(Swe, Mpot, eps_s)
    R = smooth_min(Pr, Rpot, eps_s)
    dM_dSwe = d_smooth_min_dA(Swe, Mpot, eps_s)
    dM_dMpot = d_smooth_min_dB(Swe, Mpot, eps_s)
    dR_dPr = d_smooth_min_dA(Pr, Rpot, eps_s)
    dR_dRpot = d_smooth_min_dB(Pr, Rpot, eps_s)
    dM_dT_tr = dM_dMpot * dMpot_dT_tr
    dM_df_dd = dM_dMpot * dMpot_df_dd
    dR_dT_tr = (dR_dPr * dPr_dT_tr) + (dR_dRpot * dRpot_dT_tr)
    dR_df_dd = dR_dRpot * dRpot_df_dd
    dR_dcfr = dR_dRpot * dRpot_dcfr
    Win_raw = (Pr + M) - R
    Win = smooth_pos(Win_raw, eps_s)
    dWin_dWinraw = d_smooth_pos_da(Win_raw, eps_s)
    Pliq = Win
    dPliq_dSwe = dWin_dWinraw * dM_dSwe
    dPliq_dT_tr = dWin_dWinraw * ((dPr_dT_tr + dM_dT_tr) - dR_dT_tr)
    dPliq_df_dd = dWin_dWinraw * (dM_df_dd - dR_df_dd)
    dPliq_dcfr = dWin_dWinraw * (0.0 - dR_dcfr)
    den0 = lp * f_c
    den = smooth_max(den0, eps_x, eps_x)
    dDen_dDen0 = d_smooth_max_dA(den0, eps_x, eps_x)
    ratio = Sm / den
    phiE = phiE_fun(ratio, eps_s)
    dphiE_drat = dphiE_dratio(ratio, eps_s)
    E = Ep * phiE
    dratio_dSm = 1.0 / den
    dratio_dDen = -(Sm / (den * den))
    dratio_dDen0 = dratio_dDen * dDen_dDen0
    dratio_dlp = dratio_dDen0 * f_c
    dratio_dfc = dratio_dDen0 * lp
    dE_dSm = (Ep * dphiE_drat) * dratio_dSm
    dE_dlp = (Ep * dphiE_drat) * dratio_dlp
    dE_dfc = (Ep * dphiE_drat) * dratio_dfc
    fc_safe = smooth_max(f_c, eps_x, eps_x)
    dfcSafe_dfc = d_smooth_max_dA(f_c, eps_x, eps_x)
    u_raw = Sm / fc_safe
    ucap = clamp01(u_raw, eps_s)
    du_duraw = dclamp01_dz(u_raw, eps_s)
    duraw_dSm = 1.0 / fc_safe
    duraw_dfc = (-(Sm / (fc_safe * fc_safe))) * dfcSafe_dfc
    du_dSm = du_duraw * duraw_dSm
    du_dfc = du_duraw * duraw_dfc
    up = smooth_max(ucap, eps_x, eps_x)
    dup_du = d_smooth_max_dA(ucap, eps_x, eps_x)
    up_beta = np.power(up, beta)
    Re = Pliq * up_beta
    dRe_dPliq = up_beta
    up_bm1 = np.power(up, (beta - 1.0))
    dRe_dSm = (((Pliq * beta) * up_bm1) * dup_du) * du_dSm
    dRe_dfc = (((Pliq * beta) * up_bm1) * dup_du) * du_dfc
    dRe_dbeta = (Pliq * up_beta) * np.log(up)
    h = smooth_pos((Uz - uzl), eps_s)
    dh_dUz = d_smooth_pos_da((Uz - uzl), eps_s)
    dh_duzl = -dh_dUz
    q0 = k_0 * h
    dq0_dUz = k_0 * dh_dUz
    dq0_dk0 = h
    dq0_duzl = k_0 * dh_duzl
    q1 = k_1 * Uz
    dq1_dUz = k_1
    dq1_dk1 = Uz
    perc_flux = smooth_min(perc, Uz, eps_s)
    dperc_dperc = d_smooth_min_dA(perc, Uz, eps_s)
    dperc_dUz = d_smooth_min_dB(perc, Uz, eps_s)
    q2 = k_2 * Lz
    dq2_dLz = k_2
    dq2_dk2 = Lz
    Q = (q0 + q1) + q2
    udot[0] = (Ps + R) - M
    udot[1] = (Pliq - Re) - E
    udot[2] = ((Re - q0) - q1) - perc_flux
    udot[3] = perc_flux - q2
    udot[4] = Q
    Jx_f = np.zeros(
        (
            5,
            5,
        )
    )
    Jth_f = np.zeros(
        (
            5,
            12,
        )
    )
    Jx_f[0][0] = -dM_dSwe
    Jx_f[1][0] = dPliq_dSwe - (dRe_dPliq * dPliq_dSwe)
    Jx_f[1][1] = (-dRe_dSm) - dE_dSm
    Jx_f[2][0] = dRe_dPliq * dPliq_dSwe
    Jx_f[2][1] = dRe_dSm
    Jx_f[2][2] = -((dq0_dUz + dq1_dUz) + dperc_dUz)
    Jx_f[3][2] = dperc_dUz
    Jx_f[3][3] = -dq2_dLz
    Jx_f[4][2] = dq0_dUz + dq1_dUz
    Jx_f[4][3] = dq2_dLz
    for i in range(0, m):
        Jx_f[i][0] *= dSwe_dSweRaw
        Jx_f[i][1] *= dSm_dSmRaw
        Jx_f[i][2] *= dUz_dUzRaw
        Jx_f[i][3] *= dLz_dLzRaw
    j_fc = 0
    j_beta = 1
    j_lp = 2
    j_k0 = 3
    j_uzl_ = 4
    j_k1_ = 5
    j_k2_ = 6
    j_perc_ = 7
    j_T_tr = 8
    j_f_dd = 9
    j_sfcf = 10
    j_cfr_ = 11
    Jth_f[0][j_T_tr] = (dPs_dT_tr + dR_dT_tr) - dM_dT_tr
    Jth_f[0][j_f_dd] = (0.0 + dR_df_dd) - dM_df_dd
    Jth_f[0][j_sfcf] = dPs_dsfcf
    Jth_f[0][j_cfr_] = dR_dcfr
    Jth_f[1][j_T_tr] = dPliq_dT_tr - (dRe_dPliq * dPliq_dT_tr)
    Jth_f[1][j_f_dd] = dPliq_df_dd - (dRe_dPliq * dPliq_df_dd)
    Jth_f[1][j_cfr_] = dPliq_dcfr - (dRe_dPliq * dPliq_dcfr)
    Jth_f[1][j_fc] = (-dRe_dfc) - dE_dfc
    Jth_f[1][j_beta] = -dRe_dbeta
    Jth_f[1][j_lp] = -dE_dlp
    Jth_f[2][j_T_tr] = dRe_dPliq * dPliq_dT_tr
    Jth_f[2][j_f_dd] = dRe_dPliq * dPliq_df_dd
    Jth_f[2][j_cfr_] = dRe_dPliq * dPliq_dcfr
    Jth_f[2][j_fc] = dRe_dfc
    Jth_f[2][j_beta] = dRe_dbeta
    Jth_f[2][j_k0] = -dq0_dk0
    Jth_f[2][j_uzl_] = -dq0_duzl
    Jth_f[2][j_k1_] = -dq1_dk1
    Jth_f[2][j_perc_] = -dperc_dperc
    Jth_f[3][j_k2_] = -dq2_dk2
    Jth_f[3][j_perc_] = dperc_dperc
    Jth_f[4][j_k0] = dq0_dk0
    Jth_f[4][j_uzl_] = dq0_duzl
    Jth_f[4][j_k1_] = dq1_dk1
    Jth_f[4][j_k2_] = dq2_dk2
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
    "f_c",
    "beta",
    "lp",
    "k_0",
    "uzl",
    "k_1",
    "k_2",
    "perc",
    "T_tr",
    "f_dd",
    "sfcf",
    "cfr",
    "eps_s",
    "eps_t",
    "eps_x",
]


@njit(cache=True, error_model="numpy")
def augmented(z, forcing, values):
    """Evaluate state derivatives and analytic parameter sensitivities."""
    m = 5
    d = 12
    out = np.empty_like(z)
    hbv_odefcn(
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
        m,
        d,
    )
    return out
