# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_hmodel_equations.py
#
# Python equations and analytic sensitivities matching the shared
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Python equations and analytic sensitivities matching the shared
hmodel core.
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
    denom = np.sqrt(((a * a) + (eps * eps)))
    return 0.5 * (1.0 + (a / denom))


@njit(cache=True, error_model="numpy")
def smooth_min(A, B, eps):
    """Evaluate the differentiable minimum approximation."""
    d = A - B
    return 0.5 * ((A + B) - np.sqrt(((d * d) + (eps * eps))))


@njit(cache=True, error_model="numpy")
def smooth_max(A, B, eps):
    """Evaluate the differentiable maximum approximation."""
    d = A - B
    return 0.5 * ((A + B) + np.sqrt(((d * d) + (eps * eps))))


@njit(cache=True, error_model="numpy")
def d_smooth_min_dA(A, B, eps):
    """Evaluate the analytic derivative of the smooth_min response."""
    d = A - B
    return 0.5 * (1.0 - (d / np.sqrt(((d * d) + (eps * eps)))))


@njit(cache=True, error_model="numpy")
def d_smooth_max_dA(A, B, eps):
    """Evaluate the analytic derivative of the smooth_max response."""
    d = A - B
    return 0.5 * (1.0 + (d / np.sqrt(((d * d) + (eps * eps)))))


@njit(cache=True, error_model="numpy")
def clamp01(z, eps):
    """Apply the smooth unit-interval clipping function."""
    zmax = smooth_max(z, 0.0, eps)
    return smooth_min(zmax, 1.0, eps)


@njit(cache=True, error_model="numpy")
def d_clamp01_dz(z, eps):
    """Evaluate the analytic derivative of the clamp01 response."""
    A = smooth_max(z, 0.0, eps)
    return d_smooth_min_dA(A, 1.0, eps) * d_smooth_max_dA(z, 0.0, eps)


@njit(cache=True, error_model="numpy")
def exponen(x):
    """Evaluate the smoothed exponential response."""
    return np.exp(min(300.0, x))


@njit(cache=True, error_model="numpy")
def expFlux(Sr, a):
    """Evaluate the exponential storage-flux relationship."""
    if abs(a) <= 1e-5:
        return Sr
    else:
        num = 1.0 - exponen(((-a) * Sr))
        den = 1.0 - exponen((-a))
        return num / den


@njit(cache=True, error_model="numpy")
def d_expFlux_dz(Sr, a):
    """Evaluate the analytic derivative of the expFlux response."""
    if abs(a) <= 1e-5:
        return 1.0
    den = -np.expm1((-a))
    ez = exponen(((-a) * Sr))
    return (a * ez) / den


@njit(cache=True, error_model="numpy")
def d_expFlux_da(Sr, a):
    """Evaluate the analytic derivative of the expFlux response."""
    if abs(a) < 1e-6:
        return 0.0
    ea = exponen((-a))
    ea_z = exponen(((-a) * Sr))
    N = 1.0 - ea_z
    D = 1.0 - ea
    Np = Sr * ea_z
    Dp = ea
    return ((D * Np) - (N * Dp)) / (D * D)


@njit(cache=True, error_model="numpy")
def hmodel_odefcn(
    u,
    udot,
    Smat,
    dSdt,
    P,
    Ep,
    T,
    I_max,
    Su_max,
    Q_max,
    a_E,
    a_F,
    a_S,
    r_f,
    r_s,
    T_tr,
    f_dd,
    T_sm,
    eps_m,
    m,
    d,
):
    """Evaluate the hmodel state equations and their derivatives."""
    a_I = 50.0
    a_P = -50.0
    Swe_raw = u[0]
    Si = u[1]
    Su = u[2]
    Sf = u[3]
    Ss = u[4]
    Swe = smooth_pos(Swe_raw, eps_m)
    dSwe_dx1 = d_smooth_pos_da(Swe_raw, eps_m)
    T_smeps_m = max(T_sm, eps_m)
    uT = (T - T_tr) / T_smeps_m
    snow_fr = 0.5 * (1.0 - np.tanh(uT))
    rain_fr = 1.0 - snow_fr
    P_snow = P * snow_fr
    P_rain = P * rain_fr
    a = T - T_tr
    posT = 0.5 * (a + np.sqrt(((a * a) + (T_smeps_m * T_smeps_m))))
    M_pot = f_dd * posT
    dxy = Swe - M_pot
    sqrtm = np.sqrt(((dxy * dxy) + (eps_m * eps_m)))
    M = 0.5 * ((Swe + M_pot) - sqrtm)
    Pliq = P_rain + M
    coshu = np.cosh(uT)
    sech2 = 1.0 / (coshu * coshu)
    dsnow_dT_tr = (0.5 * sech2) / T_smeps_m
    dposT_da = 0.5 * (
        1.0 + (a / np.sqrt(((a * a) + (T_smeps_m * T_smeps_m))))
    )
    dposT_dT_tr = -dposT_da
    dMpot_dT_tr = f_dd * dposT_dT_tr
    dMpot_f_dd = posT
    dM_dSwe = 0.5 * (1.0 - (dxy / sqrtm))
    dM_dMpot = 0.5 * (1.0 + (dxy / sqrtm))
    dPliq_dSwe = dM_dSwe
    dPliq_dT_tr = ((-P) * dsnow_dT_tr) + (dM_dMpot * dMpot_dT_tr)
    dPliq_f_dd = dM_dMpot * dMpot_f_dd
    dPsnow_dT_tr = P * dsnow_dT_tr
    dM_dT_tr = dM_dMpot * dMpot_dT_tr
    dM_f_dd = dM_dMpot * dMpot_f_dd
    EvapI = 0.0
    P_e = Pliq
    Ep_e = Ep
    dEvapI_dSi = 0.0
    dPe_dSi = 0.0
    dEp_e_dSi = 0.0
    dEvapI_dImax = 0.0
    dPe_dImax = 0.0
    dEp_e_dImax = 0.0
    dPe_dSwe = dPliq_dSwe
    G_I = 1.0
    if I_max > 0.0:
        Sir_raw = Si / I_max
        Sir = clamp01(Sir_raw, eps_m)
        dSir_dSirraw = d_clamp01_dz(Sir_raw, eps_m)
        dSir_dSi = dSir_dSirraw * (1.0 / I_max)
        dSir_dImax = dSir_dSirraw * ((-Si) / (I_max * I_max))
        F_I = expFlux(Sir, a_I)
        G_I = expFlux(Sir, a_P)
        EvapI = Ep * F_I
        P_e = Pliq * G_I
        dF_I_dz = d_expFlux_dz(Sir, a_I)
        dG_I_dz = d_expFlux_dz(Sir, a_P)
        dEvapI_dSi = (Ep * dF_I_dz) * dSir_dSi
        dEvapI_dImax = (Ep * dF_I_dz) * dSir_dImax
        dPe_dSi = (Pliq * dG_I_dz) * dSir_dSi
        dPe_dImax = (Pliq * dG_I_dz) * dSir_dImax
        arg = Ep - EvapI
        dEp_e_darg = d_smooth_pos_da(arg, eps_m)
        Ep_e = smooth_pos(arg, eps_m)
        dEp_e_dSi = dEp_e_darg * (-dEvapI_dSi)
        dEp_e_dImax = dEp_e_darg * (-dEvapI_dImax)
        dPe_dSwe = G_I * dPliq_dSwe
    else:
        dPe_dSwe = dPliq_dSwe
    Sur = 0.0
    dEa_dSu = 0.0
    dprc_dSu = 0.0
    drnf_dSu = 0.0
    dEa_dSumax = 0.0
    dprc_dSumax = 0.0
    drnf_dSumax = 0.0
    if Su_max > 0.0:
        Sur_raw = Su / Su_max
        Sur = clamp01(Sur_raw, eps_m)
        dSur_dSurraw = d_clamp01_dz(Sur_raw, eps_m)
        dSur_dSu = dSur_dSurraw * (1.0 / Su_max)
        dSur_dSumax = dSur_dSurraw * ((-Su) / (Su_max * Su_max))
        dF_E_dz = d_expFlux_dz(Sur, a_E)
        dF_S_dz = d_expFlux_dz(Sur, a_S)
        dF_F_dz = d_expFlux_dz(Sur, a_F)
        dEa_dSu = (Ep_e * dF_E_dz) * dSur_dSu
        dprc_dSu = (Q_max * dF_S_dz) * dSur_dSu
        drnf_dSu = (P_e * dF_F_dz) * dSur_dSu
        dEa_dSumax = (Ep_e * dF_E_dz) * dSur_dSumax
        dprc_dSumax = (Q_max * dF_S_dz) * dSur_dSumax
        drnf_dSumax = (P_e * dF_F_dz) * dSur_dSumax
    else:
        Sur = 0.0
    F_E = expFlux(Sur, a_E)
    F_S = expFlux(Sur, a_S)
    F_F = expFlux(Sur, a_F)
    Ea = Ep_e * F_E
    prc = Q_max * F_S
    rnf = P_e * F_F
    qf = Sf / r_f
    qs = Ss / r_s
    dEa_dSi = dEp_e_dSi * F_E
    drnf_dSi = dPe_dSi * F_F
    drnf_dSwe = dPe_dSwe * F_F
    dqf_dSf = 1.0 / r_f
    dqs_dSs = 1.0 / r_s
    udot[0] = P_snow - M
    udot[1] = (Pliq - EvapI) - P_e
    udot[2] = ((P_e - Ea) - prc) - rnf
    udot[3] = rnf - qf
    udot[4] = prc - qs
    udot[5] = qf + qs
    Jx_f = np.zeros(
        (
            6,
            6,
        )
    )
    Jth_f = np.zeros(
        (
            6,
            9,
        )
    )
    Jx_f[0][0] = -dM_dSwe
    Jx_f[1][0] = dPliq_dSwe - dPe_dSwe
    Jx_f[1][1] = (-dEvapI_dSi) - dPe_dSi
    Jx_f[2][0] = dPe_dSwe - drnf_dSwe
    Jx_f[2][1] = (dPe_dSi * (1.0 - F_F)) - (F_E * dEp_e_dSi)
    Jx_f[2][2] = ((-dEa_dSu) - dprc_dSu) - drnf_dSu
    Jx_f[3][0] = drnf_dSwe
    Jx_f[3][1] = drnf_dSi
    Jx_f[3][2] = drnf_dSu
    Jx_f[3][3] = -dqf_dSf
    Jx_f[4][2] = dprc_dSu
    Jx_f[4][4] = -dqs_dSs
    Jx_f[5][3] = dqf_dSf
    Jx_f[5][4] = dqs_dSs
    for ii in range(0, 6):
        Jx_f[ii][0] *= dSwe_dx1
    dEa_dImax = dEp_e_dImax * F_E
    drnf_dImax = dPe_dImax * F_F
    j = 0
    Jth_f[1][j] = (-dEvapI_dImax) - dPe_dImax
    Jth_f[2][j] = (dPe_dImax - dEa_dImax) - drnf_dImax
    Jth_f[3][j] = drnf_dImax
    j = 1
    Jth_f[2][j] = ((-dEa_dSumax) - dprc_dSumax) - drnf_dSumax
    Jth_f[3][j] = drnf_dSumax
    Jth_f[4][j] = dprc_dSumax
    j = 2
    dprc_dQmax = F_S
    Jth_f[2][j] = -dprc_dQmax
    Jth_f[4][j] = dprc_dQmax
    j = 3
    dF_E_da = d_expFlux_da(Sur, a_E)
    dEa_daE = Ep_e * dF_E_da
    Jth_f[2][j] = -dEa_daE
    j = 4
    dF_F_da = d_expFlux_da(Sur, a_F)
    drnf_daF = P_e * dF_F_da
    Jth_f[2][j] = -drnf_daF
    Jth_f[3][j] = drnf_daF
    j = 5
    Jth_f[3][j] = Sf / (r_f * r_f)
    Jth_f[5][j] = (-Sf) / (r_f * r_f)
    j = 6
    Jth_f[4][j] = Ss / (r_s * r_s)
    Jth_f[5][j] = (-Ss) / (r_s * r_s)
    j = 7
    dPe_dT_tr = (G_I * dPliq_dT_tr) if (I_max > 0.0) else dPliq_dT_tr
    drnf_dT_tr = dPe_dT_tr * F_F
    Jth_f[0][j] = dPsnow_dT_tr - dM_dT_tr
    Jth_f[1][j] = dPliq_dT_tr - dPe_dT_tr
    Jth_f[2][j] = dPe_dT_tr - drnf_dT_tr
    Jth_f[3][j] = drnf_dT_tr
    j = 8
    dPe_f_dd = (G_I * dPliq_f_dd) if (I_max > 0.0) else dPliq_f_dd
    drnf_f_dd = dPe_f_dd * F_F
    Jth_f[0][j] = -dM_f_dd
    Jth_f[1][j] = dPliq_f_dd - dPe_f_dd
    Jth_f[2][j] = dPe_f_dd - drnf_f_dd
    Jth_f[3][j] = drnf_f_dd
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
    "I_max",
    "Su_max",
    "Q_max",
    "a_E",
    "a_F",
    "a_S",
    "r_f",
    "r_s",
    "T_tr",
    "f_dd",
    "T_sm",
    "eps_m",
]


@njit(cache=True, error_model="numpy")
def augmented(z, forcing, values):
    """Evaluate state derivatives and analytic parameter sensitivities."""
    m = 6
    d = 9
    out = np.empty_like(z)
    hmodel_odefcn(
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
        m,
        d,
    )
    return out
