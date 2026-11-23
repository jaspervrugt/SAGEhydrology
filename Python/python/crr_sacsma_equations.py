# ==========================================================================
# SAGEhydrology | Python source
# File: python/crr_sacsma_equations.py
#
# Python equations and analytic sensitivities matching the shared
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Python equations and analytic sensitivities matching the shared
sacsma core.
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
def smooth_pos(a, ee):
    """Evaluate the differentiable positive-part approximation."""
    return 0.5 * (a + np.sqrt(((a * a) + (ee * ee))))


@njit(cache=True, error_model="numpy")
def dsmooth_pos_da(a, ee):
    """Evaluate the analytic derivative of the smooth_pos response."""
    return 0.5 * (1.0 + (a / np.sqrt(((a * a) + (ee * ee)))))


@njit(cache=True, error_model="numpy")
def smooth_min(a, b, ee):
    """Evaluate the differentiable minimum approximation."""
    d = a - b
    return 0.5 * ((a + b) - np.sqrt(((d * d) + (ee * ee))))


@njit(cache=True, error_model="numpy")
def dsmoothmin_da(a, b, ee):
    """Evaluate the analytic derivative of the smoothmin response."""
    d = a - b
    return 0.5 * (1.0 - (d / np.sqrt(((d * d) + (ee * ee)))))


@njit(cache=True, error_model="numpy")
def dsmoothmin_db(a, b, ee):
    """Evaluate the analytic derivative of the smoothmin response."""
    d = a - b
    return 0.5 * (1.0 + (d / np.sqrt(((d * d) + (ee * ee)))))


@njit(cache=True, error_model="numpy")
def smooth_clip01(y, ee):
    """Apply the smooth unit-interval clipping function."""
    sp1 = smooth_pos(y, ee)
    sp2 = smooth_pos((1.0 - sp1), ee)
    return 1.0 - sp2


@njit(cache=True, error_model="numpy")
def dsmooth_clip01_dy(y, ee):
    """Evaluate the analytic derivative of the smooth_clip01 response."""
    sp1 = smooth_pos(y, ee)
    dsp1 = dsmooth_pos_da(y, ee)
    dsp2 = dsmooth_pos_da((1.0 - sp1), ee)
    return dsp2 * dsp1


@njit(cache=True, error_model="numpy")
def lf_func(S, Smax, eps, rho):
    """Evaluate the nonlinear storage-dependent flux."""
    denom = rho * Smax
    z = (S - (Smax - ((rho * Smax) * eps))) / denom
    return 1.0 / (1.0 + np.exp((-z)))


@njit(cache=True, error_model="numpy")
def dlf_dS(lf_val, Smax, rho):
    """Differentiate the nonlinear flux with respect to storage."""
    return (lf_val * (1.0 - lf_val)) / (rho * Smax)


@njit(cache=True, error_model="numpy")
def sacsma_odefcn(
    u,
    udot,
    Smat,
    dSdt,
    P,
    Ep,
    T,
    uzfwm,
    uztwm,
    lzfpm,
    lzfsm,
    lztwm,
    zperc,
    rexp,
    uzk,
    pfree,
    lzpk,
    lzsk,
    acm,
    kf,
    T_tr,
    f_dd,
    T_sm,
    eps_m,
    eps_s,
    eps,
    rho,
    m,
    d,
):
    """Evaluate the sacsma state equations and their derivatives."""
    Swe_raw = u[0]
    Uztw = u[1]
    Uzfw = u[2]
    Lztw = u[3]
    Lzps = u[4]
    Lzfs = u[5]
    Sf1 = u[6]
    Sf2 = u[7]
    Sf3 = u[8]
    Swe = smooth_pos(Swe_raw, eps_s)
    dSwe_dRaw = dsmooth_pos_da(Swe_raw, eps_s)
    Uztw_pos = smooth_pos(Uztw, eps_s)
    dUztwpos_dUztw = dsmooth_pos_da(Uztw, eps_s)
    Uzfw_pos = smooth_pos(Uzfw, eps_s)
    dUzfwpos_dUzfw = dsmooth_pos_da(Uzfw, eps_s)
    Lztw_pos = smooth_pos(Lztw, eps_s)
    dLztwpos_dLztw = dsmooth_pos_da(Lztw, eps_s)
    Lztot = (Lztw + Lzps) + Lzfs
    Lzm = (lzfpm + lzfsm) + lztwm
    T_smeps_m = max(T_sm, eps_m)
    uT = (T - T_tr) / T_smeps_m
    snow_fr = 0.5 * (1.0 - np.tanh(uT))
    rain_fr = 1.0 - snow_fr
    P_snow = P * snow_fr
    P_rain = P * rain_fr
    aT = T - T_tr
    posT = smooth_pos(aT, T_smeps_m)
    M_pot = f_dd * posT
    M = smooth_min(Swe, M_pot, eps_m)
    Pliq = P_rain + M
    c = np.cosh(uT)
    sech2 = 1.0 / (c * c)
    dsnow_dT_tr = (0.5 * sech2) / T_smeps_m
    dposT_da = dsmooth_pos_da(aT, T_smeps_m)
    dposT_dT_tr = -dposT_da
    dMpot_dT_tr = f_dd * dposT_dT_tr
    dMpot_f_dd = posT
    dM_dSwe_s = dsmoothmin_da(Swe, M_pot, eps_m)
    dM_dSwe = dM_dSwe_s * dSwe_dRaw
    dM_dMpot = dsmoothmin_db(Swe, M_pot, eps_m)
    dM_dT_tr = dM_dMpot * dMpot_dT_tr
    dM_f_dd = dM_dMpot * dMpot_f_dd
    dPliq_dSwe = dM_dSwe
    dPliq_dT_tr = ((-P) * dsnow_dT_tr) + dM_dT_tr
    dPliq_f_dd = dM_f_dd
    ratioU = Uztw_pos / uztwm
    ratioU_clip = smooth_clip01(ratioU, eps_s)
    dratioUclip_dy = dsmooth_clip01_dy(ratioU, eps_s)
    e1 = Ep * ratioU_clip
    de1_dUztw = (Ep * dratioUclip_dy) * (dUztwpos_dUztw / uztwm)
    de1_duztwm = (Ep * dratioUclip_dy) * ((-Uztw_pos) / (uztwm * uztwm))
    ratioL = Lztw_pos / lztwm
    cL = smooth_clip01(ratioL, eps_s)
    dratioLclip_dy = dsmooth_clip01_dy(ratioL, eps_s)
    dcL_dLztw = dratioLclip_dy * (dLztwpos_dLztw / lztwm)
    dcL_dlztwm = dratioLclip_dy * ((-Lztw_pos) / (lztwm * lztwm))
    e2 = (Ep - e1) * cL
    de2_dLztw = (Ep - e1) * dcL_dLztw
    de2_dlztwm = (Ep - e1) * dcL_dlztwm
    de2_dUztw = (-cL) * de1_dUztw
    de2_duztwm = (-cL) * de1_duztwm
    q0 = (lzpk * lzfpm) + (lzsk * lzfsm)
    Araw = 1.0 - (Lztot / Lzm)
    Apos = smooth_pos(Araw, eps_s)
    dApos_dAraw = dsmooth_pos_da(Araw, eps_s)
    Afloor = 1e-12
    Aeff = smooth_pos((Apos - Afloor), eps_s) + Afloor
    dAeff_dApos = dsmooth_pos_da((Apos - Afloor), eps_s)
    logAeff = np.log(Aeff)
    Apow = np.exp((rexp * logAeff))
    dlz = 1.0 + (zperc * Apow)
    ddlz_dzperc = Apow
    ddlz_drexp = (zperc * Apow) * logAeff
    dAraw_dLztw = -(1.0 / Lzm)
    dAraw_dLzps = -(1.0 / Lzm)
    dAraw_dLzfs = -(1.0 / Lzm)
    dAraw_dLzm = Lztot / (Lzm * Lzm)
    dAraw_dlzfpm = dAraw_dLzm
    dAraw_dlzfsm = dAraw_dLzm
    dAraw_dlztwm = dAraw_dLzm
    dApos_dLztw = dApos_dAraw * dAraw_dLztw
    dApos_dLzps = dApos_dAraw * dAraw_dLzps
    dApos_dLzfs = dApos_dAraw * dAraw_dLzfs
    dApos_dlzfpm = dApos_dAraw * dAraw_dlzfpm
    dApos_dlzfsm = dApos_dAraw * dAraw_dlzfsm
    dApos_dlztwm = dApos_dAraw * dAraw_dlztwm
    dApow_dApos = (rexp * np.exp(((rexp - 1.0) * logAeff))) * dAeff_dApos
    ddlz_dLztw = (zperc * dApow_dApos) * dApos_dLztw
    ddlz_dLzps = (zperc * dApow_dApos) * dApos_dLzps
    ddlz_dLzfs = (zperc * dApow_dApos) * dApos_dLzfs
    ddlz_dlzfpm = (zperc * dApow_dApos) * dApos_dlzfpm
    ddlz_dlzfsm = (zperc * dApow_dApos) * dApos_dlzfsm
    ddlz_dlztwm = (zperc * dApow_dApos) * dApos_dlztwm
    B = Uzfw_pos / uzfwm
    dB_dUzfw = dUzfwpos_dUzfw / uzfwm
    dB_duzfwm = (-Uzfw_pos) / (uzfwm * uzfwm)
    q12 = (q0 * dlz) * B
    dq12_dUzfw = (q0 * dlz) * dB_dUzfw
    dq12_dLztw = (q0 * B) * ddlz_dLztw
    dq12_dLzps = (q0 * B) * ddlz_dLzps
    dq12_dLzfs = (q0 * B) * ddlz_dLzfs
    dq12_dzperc = (q0 * B) * ddlz_dzperc
    dq12_drexp = (q0 * B) * ddlz_drexp
    qif = uzk * B
    dqif_dUzfw = uzk * dB_dUzfw
    qbp = lzpk * Lzps
    qbs = lzsk * Lzfs
    dqbp_dLzps = lzpk
    dqbs_dLzfs = lzsk
    dqbp_dlzpk = Lzps
    dqbs_dlzsk = Lzfs
    ac = acm * (Uztw_pos / uztwm)
    qsx = ac * Pliq
    dac_dUztw = acm * (dUztwpos_dUztw / uztwm)
    dac_duztwm = acm * ((-Uztw_pos) / (uztwm * uztwm))
    dac_dacm = Uztw_pos / uztwm
    dqsx_dUztw = Pliq * dac_dUztw
    dqsx_duztwm = Pliq * dac_duztwm
    dqsx_dacm = Pliq * dac_dacm
    dqsx_dSwe = ac * dPliq_dSwe
    dqsx_dT_tr = ac * dPliq_dT_tr
    dqsx_f_dd = ac * dPliq_f_dd
    lfUztw = lf_func(Uztw, uztwm, eps, rho)
    lfUzfw = lf_func(Uzfw, uzfwm, eps, rho)
    lfLztw = lf_func(Lztw, lztwm, eps, rho)
    lfLzps = lf_func(Lzps, lzfpm, eps, rho)
    lfLzfs = lf_func(Lzfs, lzfsm, eps, rho)
    dlfUztw = dlf_dS(lfUztw, uztwm, rho)
    dlfUzfw = dlf_dS(lfUzfw, uzfwm, rho)
    dlfLztw = dlf_dS(lfLztw, lztwm, rho)
    dlfLzps = dlf_dS(lfLzps, lzfpm, rho)
    dlfLzfs = dlf_dS(lfLzfs, lzfsm, rho)
    qutof = (Pliq - qsx) * lfUztw
    qufof = qutof * lfUzfw
    qstof = (pfree * q12) * lfLztw
    prcs = ((0.5 * (1.0 - pfree)) * q12) + (0.5 * qstof)
    qsfofa = prcs * lfLzps
    qsfofb = prcs * lfLzfs
    qf1 = kf * Sf1
    qf2 = kf * Sf2
    qf3 = kf * Sf3
    udot[0] = P_snow - M
    udot[1] = ((Pliq - qsx) - e1) - qutof
    udot[2] = ((qutof - q12) - qif) - qufof
    udot[3] = ((pfree * q12) - e2) - qstof
    udot[4] = (prcs - qbp) - qsfofa
    udot[5] = (prcs - qbs) - qsfofb
    udot[6] = (
        (((((qif + qsx) + qufof) + qsfofa) + qsfofb) + qbp) + qbs
    ) - qf1
    udot[7] = qf1 - qf2
    udot[8] = qf2 - qf3
    udot[9] = qf3
    Jx_f = np.zeros(
        (
            10,
            10,
        )
    )
    Jth_f = np.zeros(
        (
            10,
            15,
        )
    )
    dlfU_dUztwm = (-dlfUztw) * (Uztw / uztwm)
    dqutof_duztwm = ((-dqsx_duztwm) * lfUztw) + ((Pliq - qsx) * dlfU_dUztwm)
    dqutof_dSwe = (dPliq_dSwe - dqsx_dSwe) * lfUztw
    dqutof_dT_tr = (dPliq_dT_tr - dqsx_dT_tr) * lfUztw
    dqutof_f_dd = (dPliq_f_dd - dqsx_f_dd) * lfUztw
    dqutof_dUztw = ((-dqsx_dUztw) * lfUztw) + ((Pliq - qsx) * dlfUztw)
    dqutof_dacm = (-dqsx_dacm) * lfUztw
    dqufof_dSwe = dqutof_dSwe * lfUzfw
    dqufof_dT_tr = dqutof_dT_tr * lfUzfw
    dqufof_f_dd = dqutof_f_dd * lfUzfw
    dqufof_dUztw = dqutof_dUztw * lfUzfw
    dqufof_dUzfw = qutof * dlfUzfw
    dqufof_dacm = dqutof_dacm * lfUzfw
    dqufof_duztwm = dqutof_duztwm * lfUzfw
    dqstof_dUzfw = (pfree * dq12_dUzfw) * lfLztw
    dqstof_dLztw = ((pfree * dq12_dLztw) * lfLztw) + (
        (pfree * q12) * dlfLztw
    )
    dqstof_dLzps = (pfree * dq12_dLzps) * lfLztw
    dqstof_dLzfs = (pfree * dq12_dLzfs) * lfLztw
    dqstof_dzperc = (pfree * dq12_dzperc) * lfLztw
    dqstof_drexp = (pfree * dq12_drexp) * lfLztw
    dqstof_dpfree = q12 * lfLztw
    dprc_dUzfw = ((0.5 * (1.0 - pfree)) * dq12_dUzfw) + (0.5 * dqstof_dUzfw)
    dprc_dLztw = ((0.5 * (1.0 - pfree)) * dq12_dLztw) + (0.5 * dqstof_dLztw)
    dprc_dLzps = ((0.5 * (1.0 - pfree)) * dq12_dLzps) + (0.5 * dqstof_dLzps)
    dprc_dLzfs = ((0.5 * (1.0 - pfree)) * dq12_dLzfs) + (0.5 * dqstof_dLzfs)
    dprc_dzperc = ((0.5 * (1.0 - pfree)) * dq12_dzperc) + (
        0.5 * dqstof_dzperc
    )
    dprc_drexp = ((0.5 * (1.0 - pfree)) * dq12_drexp) + (0.5 * dqstof_drexp)
    dprc_dpfree = ((-0.5) * q12) + (0.5 * dqstof_dpfree)
    dqsfpf_dUzfw = dprc_dUzfw * lfLzps
    dqsfpf_dLztw = dprc_dLztw * lfLzps
    dqsfpf_dLzps = (dprc_dLzps * lfLzps) + (prcs * dlfLzps)
    dqsfpf_dLzfs = dprc_dLzfs * lfLzps
    dqsfpf_dzperc = dprc_dzperc * lfLzps
    dqsfpf_drexp = dprc_drexp * lfLzps
    dqsfpf_dpfree = dprc_dpfree * lfLzps
    dqsfss_dUzfw = dprc_dUzfw * lfLzfs
    dqsfss_dLztw = dprc_dLztw * lfLzfs
    dqsfss_dLzps = dprc_dLzps * lfLzfs
    dqsfss_dLzfs = (dprc_dLzfs * lfLzfs) + (prcs * dlfLzfs)
    dqsfss_dzperc = dprc_dzperc * lfLzfs
    dqsfss_drexp = dprc_drexp * lfLzfs
    dqsfss_dpfree = dprc_dpfree * lfLzfs
    Jx_f[0][0] = -dM_dSwe
    Jx_f[1][0] = (dPliq_dSwe - dqsx_dSwe) - dqutof_dSwe
    Jx_f[1][1] = ((-dqsx_dUztw) - de1_dUztw) - dqutof_dUztw
    Jx_f[2][0] = dqutof_dSwe - dqufof_dSwe
    Jx_f[2][1] = dqutof_dUztw - dqufof_dUztw
    Jx_f[2][2] = ((-dqif_dUzfw) - dq12_dUzfw) - dqufof_dUzfw
    Jx_f[2][3] = -dq12_dLztw
    Jx_f[2][4] = -dq12_dLzps
    Jx_f[2][5] = -dq12_dLzfs
    Jx_f[3][1] = -de2_dUztw
    Jx_f[3][2] = (pfree * dq12_dUzfw) - dqstof_dUzfw
    Jx_f[3][3] = ((pfree * dq12_dLztw) - de2_dLztw) - dqstof_dLztw
    Jx_f[3][4] = (pfree * dq12_dLzps) - dqstof_dLzps
    Jx_f[3][5] = (pfree * dq12_dLzfs) - dqstof_dLzfs
    Jx_f[4][2] = dprc_dUzfw - dqsfpf_dUzfw
    Jx_f[4][3] = dprc_dLztw - dqsfpf_dLztw
    Jx_f[4][4] = (dprc_dLzps - dqbp_dLzps) - dqsfpf_dLzps
    Jx_f[4][5] = dprc_dLzfs - dqsfpf_dLzfs
    Jx_f[5][2] = dprc_dUzfw - dqsfss_dUzfw
    Jx_f[5][3] = dprc_dLztw - dqsfss_dLztw
    Jx_f[5][4] = dprc_dLzps - dqsfss_dLzps
    Jx_f[5][5] = (dprc_dLzfs - dqbs_dLzfs) - dqsfss_dLzfs
    Jx_f[6][0] = dqsx_dSwe + dqufof_dSwe
    Jx_f[6][1] = dqsx_dUztw + dqufof_dUztw
    Jx_f[6][2] = ((dqif_dUzfw + dqufof_dUzfw) + dqsfpf_dUzfw) + dqsfss_dUzfw
    Jx_f[6][3] = dqsfpf_dLztw + dqsfss_dLztw
    Jx_f[6][4] = (dqsfpf_dLzps + dqbp_dLzps) + dqsfss_dLzps
    Jx_f[6][5] = (dqsfpf_dLzfs + dqsfss_dLzfs) + dqbs_dLzfs
    Jx_f[6][6] = -kf
    Jx_f[7][6] = kf
    Jx_f[7][7] = -kf
    Jx_f[8][7] = kf
    Jx_f[8][8] = -kf
    Jx_f[9][8] = kf
    Jth_f[0][13] = (P * dsnow_dT_tr) - dM_dT_tr
    Jth_f[0][14] = -dM_f_dd
    Jth_f[1][13] = (dPliq_dT_tr - dqsx_dT_tr) - dqutof_dT_tr
    Jth_f[1][14] = (dPliq_f_dd - dqsx_f_dd) - dqutof_f_dd
    Jth_f[2][13] = dqutof_dT_tr - dqufof_dT_tr
    Jth_f[2][14] = dqutof_f_dd - dqufof_f_dd
    Jth_f[6][13] = dqsx_dT_tr + dqufof_dT_tr
    Jth_f[6][14] = dqsx_f_dd + dqufof_f_dd
    dq12_duzfwm = (q0 * dlz) * dB_duzfwm
    dqif_duzfwm = uzk * dB_duzfwm
    dqstof_duzfwm = (pfree * lfLztw) * dq12_duzfwm
    dprc_duzfwm = ((0.5 * (1.0 - pfree)) * dq12_duzfwm) + (
        0.5 * dqstof_duzfwm
    )
    dqsfpf_duzfwm = dprc_duzfwm * lfLzps
    dqsfss_duzfwm = dprc_duzfwm * lfLzfs
    dlfUzfw_duzfwm = (-dlfUzfw) * (Uzfw / uzfwm)
    dqufof_duzfwm = qutof * dlfUzfw_duzfwm
    Jth_f[1][0] = 0.0
    Jth_f[2][0] = ((-dq12_duzfwm) - dqif_duzfwm) - dqufof_duzfwm
    Jth_f[3][0] = (pfree * dq12_duzfwm) - dqstof_duzfwm
    Jth_f[4][0] = dprc_duzfwm - dqsfpf_duzfwm
    Jth_f[5][0] = dprc_duzfwm - dqsfss_duzfwm
    Jth_f[6][0] = (
        (dqif_duzfwm + dqufof_duzfwm) + dqsfpf_duzfwm
    ) + dqsfss_duzfwm
    Jth_f[1][1] = ((-dqsx_duztwm) - de1_duztwm) - dqutof_duztwm
    Jth_f[2][1] = dqutof_duztwm - dqufof_duztwm
    Jth_f[3][1] = -de2_duztwm
    Jth_f[6][1] = dqsx_duztwm + dqufof_duztwm
    dq0_dlzfpm = lzpk
    dq12_dlzfpm = ((dq0_dlzfpm * dlz) * B) + ((q0 * ddlz_dlzfpm) * B)
    dqstof_dlzfpm = (pfree * lfLztw) * dq12_dlzfpm
    dprc_dlzfpm = ((0.5 * (1.0 - pfree)) * dq12_dlzfpm) + (
        0.5 * dqstof_dlzfpm
    )
    dlfLzps_dlzfpm = (-dlfLzps) * (Lzps / lzfpm)
    dqsfpf_dlzfpm = (dprc_dlzfpm * lfLzps) + (prcs * dlfLzps_dlzfpm)
    dqsfss_dlzfpm = dprc_dlzfpm * lfLzfs
    Jth_f[2][2] = -dq12_dlzfpm
    Jth_f[3][2] = (pfree * dq12_dlzfpm) - dqstof_dlzfpm
    Jth_f[4][2] = dprc_dlzfpm - dqsfpf_dlzfpm
    Jth_f[5][2] = dprc_dlzfpm - dqsfss_dlzfpm
    Jth_f[6][2] = dqsfpf_dlzfpm + dqsfss_dlzfpm
    dq0_dlzfsm = lzsk
    dq12_dlzfsm = ((dq0_dlzfsm * dlz) * B) + ((q0 * ddlz_dlzfsm) * B)
    dqstof_dlzfsm = (pfree * lfLztw) * dq12_dlzfsm
    dprc_dlzfsm = ((0.5 * (1.0 - pfree)) * dq12_dlzfsm) + (
        0.5 * dqstof_dlzfsm
    )
    dlfLzfs_dlzfsm = (-dlfLzfs) * (Lzfs / lzfsm)
    dqsfpf_dlzfsm = dprc_dlzfsm * lfLzps
    dqsfss_dlzfsm = (dprc_dlzfsm * lfLzfs) + (prcs * dlfLzfs_dlzfsm)
    Jth_f[2][3] = -dq12_dlzfsm
    Jth_f[3][3] = (pfree * dq12_dlzfsm) - dqstof_dlzfsm
    Jth_f[4][3] = dprc_dlzfsm - dqsfpf_dlzfsm
    Jth_f[5][3] = dprc_dlzfsm - dqsfss_dlzfsm
    Jth_f[6][3] = dqsfpf_dlzfsm + dqsfss_dlzfsm
    dq12_dlztwm = (q0 * ddlz_dlztwm) * B
    dlfLztw_dlztwm = (-dlfLztw) * (Lztw / lztwm)
    dqstof_dlztwm = pfree * (
        (dq12_dlztwm * lfLztw) + (q12 * dlfLztw_dlztwm)
    )
    dprc_dlztwm = ((0.5 * (1.0 - pfree)) * dq12_dlztwm) + (
        0.5 * dqstof_dlztwm
    )
    dqsfpf_dlztwm = dprc_dlztwm * lfLzps
    dqsfss_dlztwm = dprc_dlztwm * lfLzfs
    Jth_f[2][4] = -dq12_dlztwm
    Jth_f[3][4] = ((pfree * dq12_dlztwm) - de2_dlztwm) - dqstof_dlztwm
    Jth_f[4][4] = dprc_dlztwm - dqsfpf_dlztwm
    Jth_f[5][4] = dprc_dlztwm - dqsfss_dlztwm
    Jth_f[6][4] = dqsfpf_dlztwm + dqsfss_dlztwm
    Jth_f[2][5] = -dq12_dzperc
    Jth_f[3][5] = (pfree * dq12_dzperc) - dqstof_dzperc
    Jth_f[4][5] = dprc_dzperc - dqsfpf_dzperc
    Jth_f[5][5] = dprc_dzperc - dqsfss_dzperc
    Jth_f[6][5] = dqsfpf_dzperc + dqsfss_dzperc
    Jth_f[2][6] = -dq12_drexp
    Jth_f[3][6] = (pfree * dq12_drexp) - dqstof_drexp
    Jth_f[4][6] = dprc_drexp - dqsfpf_drexp
    Jth_f[5][6] = dprc_drexp - dqsfss_drexp
    Jth_f[6][6] = dqsfpf_drexp + dqsfss_drexp
    Jth_f[2][7] = -B
    Jth_f[6][7] = B
    Jth_f[3][8] = q12 - dqstof_dpfree
    Jth_f[4][8] = dprc_dpfree - dqsfpf_dpfree
    Jth_f[5][8] = dprc_dpfree - dqsfss_dpfree
    Jth_f[6][8] = dqsfpf_dpfree + dqsfss_dpfree
    q12_lzpk = (lzfpm * dlz) * B
    qstof_lzpk = (pfree * q12_lzpk) * lfLztw
    prc_lzpk = ((0.5 * (1.0 - pfree)) * q12_lzpk) + (0.5 * qstof_lzpk)
    qsfp_lzpk = prc_lzpk * lfLzps
    qsfs_lzpk = prc_lzpk * lfLzfs
    Jth_f[2][9] = -q12_lzpk
    Jth_f[3][9] = (pfree * q12_lzpk) - qstof_lzpk
    Jth_f[4][9] = (prc_lzpk - dqbp_dlzpk) - qsfp_lzpk
    Jth_f[5][9] = prc_lzpk - qsfs_lzpk
    Jth_f[6][9] = (qsfp_lzpk + qsfs_lzpk) + dqbp_dlzpk
    q12_lzsk = (lzfsm * dlz) * B
    qstof_lzsk = (pfree * q12_lzsk) * lfLztw
    prc_lzsk = ((0.5 * (1.0 - pfree)) * q12_lzsk) + (0.5 * qstof_lzsk)
    qsfp_lzsk = prc_lzsk * lfLzps
    qsfs_lzsk = prc_lzsk * lfLzfs
    Jth_f[2][10] = -q12_lzsk
    Jth_f[3][10] = (pfree * q12_lzsk) - qstof_lzsk
    Jth_f[4][10] = prc_lzsk - qsfp_lzsk
    Jth_f[5][10] = (prc_lzsk - dqbs_dlzsk) - qsfs_lzsk
    Jth_f[6][10] = (qsfp_lzsk + qsfs_lzsk) + dqbs_dlzsk
    Jth_f[1][11] = (-dqsx_dacm) - dqutof_dacm
    Jth_f[2][11] = dqutof_dacm - dqufof_dacm
    Jth_f[6][11] = dqsx_dacm + dqufof_dacm
    Jth_f[6][12] = -Sf1
    Jth_f[7][12] = Sf1 - Sf2
    Jth_f[8][12] = Sf2 - Sf3
    Jth_f[9][12] = Sf3
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
    "uzfwm",
    "uztwm",
    "lzfpm",
    "lzfsm",
    "lztwm",
    "zperc",
    "rexp",
    "uzk",
    "pfree",
    "lzpk",
    "lzsk",
    "acm",
    "kf",
    "T_tr",
    "f_dd",
    "T_sm",
    "eps_m",
    "eps_s",
    "eps",
    "rho",
]


@njit(cache=True, error_model="numpy")
def augmented(z, forcing, values):
    """Evaluate state derivatives and analytic parameter sensitivities."""
    m = 10
    d = 15
    out = np.empty_like(z)
    sacsma_odefcn(
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
        m,
        d,
    )
    return out
