/*
 * sacsma.cpp
 *
 * Conceptual SAC-SMA rainfall-runoff model with an adaptive-step
 * explicit Runge-Kutta integrator. This file contains the MATLAB-independent
 * native numerical core shared by crr_sacsma and crr_model_mex.
 *
 * Written by Jasper A. Vrugt.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#include "sacsma.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <cstdio>

static void rk2(int nvar,
                int m,
                int d,
                double h,
                double* z,
                double* LTE,
                double* zdotE,
                double* zE,
                double* zdot,
                double P,
                double Ep,
                double T,
                double uzfwm,
                double uztwm,
                double lzfpm,
                double lzfsm,
                double lztwm,
                double zperc,
                double rexp,
                double uzk,
                double pfree,
                double lzpk,
                double lzsk,
                double acm,
                double kf,
                double T_tr,
                double f_dd,
                double T_sm,
                double eps_m,
                double eps_s,
                double eps,
                double rho);

static void sacsma_aug_ode(double* z,
                           double* zdot,
                           double P,
                           double Ep,
                           double T,
                           double uzfwm,
                           double uztwm,
                           double lzfpm,
                           double lzfsm,
                           double lztwm,
                           double zperc,
                           double rexp,
                           double uzk,
                           double pfree,
                           double lzpk,
                           double lzsk,
                           double acm,
                           double kf,
                           double T_tr,
                           double f_dd,
                           double T_sm,
                           double eps_m,
                           double eps_s,
                           double eps,
                           double rho,
                           int m,
                           int d);

static void sacsma_odefcn(double* u,
                          double* udot,
                          const double* Smat,
                          double* dSdt,
                          double P,
                          double Ep,
                          double T,
                          double uzfwm,
                          double uztwm,
                          double lzfpm,
                          double lzfsm,
                          double lztwm,
                          double zperc,
                          double rexp,
                          double uzk,
                          double pfree,
                          double lzpk,
                          double lzsk,
                          double acm,
                          double kf,
                          double T_tr,
                          double f_dd,
                          double T_sm,
                          double eps_m,
                          double eps_s,
                          double eps,
                          double rho,
                          int m,
                          int d);

static inline double smooth_pos(double a, double ee)
{
    return 0.5 * (a + sqrt(a * a + ee * ee));
}

static inline double dsmooth_pos_da(double a, double ee)
{
    return 0.5 * (1.0 + a / sqrt(a * a + ee * ee));
}

static inline double smooth_min(double a, double b, double ee)
{
    const double d = a - b;
    return 0.5 * (a + b - sqrt(d * d + ee * ee));
}

static inline double dsmoothmin_da(double a, double b, double ee)
{
    const double d = a - b;
    return 0.5 * (1.0 - d / sqrt(d * d + ee * ee));
}

static inline double dsmoothmin_db(double a, double b, double ee)
{
    const double d = a - b;
    return 0.5 * (1.0 + d / sqrt(d * d + ee * ee));
}

static inline double smooth_clip01(double y, double ee)
{
    const double sp1 = smooth_pos(y, ee);
    const double sp2 = smooth_pos(1.0 - sp1, ee);
    return 1.0 - sp2;
}

static inline double dsmooth_clip01_dy(double y, double ee)
{
    const double sp1 = smooth_pos(y, ee);
    const double dsp1 = dsmooth_pos_da(y, ee);
    const double dsp2 = dsmooth_pos_da(1.0 - sp1, ee);
    return dsp2 * dsp1;
}

static inline double lf_func(double S, double Smax, double eps, double rho)
{
    const double denom = rho * Smax;
    const double z = (S - (Smax - rho * Smax * eps)) / denom;
    return 1.0 / (1.0 + exp(-z));
}

static inline double dlf_dS(double lf_val, double Smax, double rho)
{
    return (lf_val * (1.0 - lf_val)) / (rho * Smax);
}

static inline bool isFinite(double x)
{
    return std::isfinite(x);
}

static void rk2(int nvar,
                int m,
                int d,
                double h,
                double* z,
                double* LTE,
                double* zdotE,
                double* zE,
                double* zdot,
                double P,
                double Ep,
                double T,
                double uzfwm,
                double uztwm,
                double lzfpm,
                double lzfsm,
                double lztwm,
                double zperc,
                double rexp,
                double uzk,
                double pfree,
                double lzpk,
                double lzsk,
                double acm,
                double kf,
                double T_tr,
                double f_dd,
                double T_sm,
                double eps_m,
                double eps_s,
                double eps,
                double rho)
{
    /* Euler */
    sacsma_aug_ode(z,
                   zdotE,
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
                   d);
    for (int i = 0; i < nvar; ++i) {
        zE[i] = z[i] + h * zdotE[i];
    }

    /* Heun slope */
    sacsma_aug_ode(zE,
                   zdot,
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
                   d);

    /* Final update */
    for (int i = 0; i < nvar; ++i) {
        z[i] = z[i] + 0.5 * h * (zdotE[i] + zdot[i]);
    }

    /* LTE */
    for (int i = 0; i < nvar; ++i) {
        LTE[i] = fabs(zE[i] - z[i]);
    }
}

static void sacsma_aug_ode(double* z,
                           double* zdot,
                           double P,
                           double Ep,
                           double T,
                           double uzfwm,
                           double uztwm,
                           double lzfpm,
                           double lzfsm,
                           double lztwm,
                           double zperc,
                           double rexp,
                           double uzk,
                           double pfree,
                           double lzpk,
                           double lzsk,
                           double acm,
                           double kf,
                           double T_tr,
                           double f_dd,
                           double T_sm,
                           double eps_m,
                           double eps_s,
                           double eps,
                           double rho,
                           int m,
                           int d)
{
    double* u = z;        /* m */
    double* Smat = z + m; /* m*d */

    double* udot = zdot;     /* m */
    double* dSdt = zdot + m; /* m*d */

    sacsma_odefcn(u,
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
                  d);
}

static void sacsma_odefcn(double* u,
                          double* udot,
                          const double* Smat,
                          double* dSdt,
                          double P,
                          double Ep,
                          double T,
                          double uzfwm,
                          double uztwm,
                          double lzfpm,
                          double lzfsm,
                          double lztwm,
                          double zperc,
                          double rexp,
                          double uzk,
                          double pfree,
                          double lzpk,
                          double lzsk,
                          double acm,
                          double kf,
                          double T_tr,
                          double f_dd,
                          double T_sm,
                          double eps_m,
                          double eps_s,
                          double eps,
                          double rho,
                          int m,
                          int d)
{
    /* ---------- 1) unpack ---------- */
    const double Swe_raw = u[0];
    const double Uztw = u[1];
    const double Uzfw = u[2];
    const double Lztw = u[3];
    const double Lzps = u[4];
    const double Lzfs = u[5];
    const double Sf1 = u[6];
    const double Sf2 = u[7];
    const double Sf3 = u[8];
    /* u[9] = Sq */

    /* ---------- 2) smoothed positives used in fluxes ---------- */
    const double Swe = smooth_pos(Swe_raw, eps_s);
    const double dSwe_dRaw = dsmooth_pos_da(Swe_raw, eps_s);

    const double Uztw_pos = smooth_pos(Uztw, eps_s);
    const double dUztwpos_dUztw = dsmooth_pos_da(Uztw, eps_s);

    const double Uzfw_pos = smooth_pos(Uzfw, eps_s);
    const double dUzfwpos_dUzfw = dsmooth_pos_da(Uzfw, eps_s);

    const double Lztw_pos = smooth_pos(Lztw, eps_s);
    const double dLztwpos_dLztw = dsmooth_pos_da(Lztw, eps_s);

    const double Lztot = Lztw + Lzps + Lzfs;
    const double Lzm = lzfpm + lzfsm + lztwm;

    /* ---------- 3) snow module ---------- */
    const double T_smeps_m = fmax(T_sm, eps_m);
    const double uT = (T - T_tr) / T_smeps_m;

    const double snow_fr = 0.5 * (1.0 - tanh(uT));
    const double rain_fr = 1.0 - snow_fr;

    const double P_snow = P * snow_fr;
    const double P_rain = P * rain_fr;

    const double aT = (T - T_tr);
    const double posT = smooth_pos(aT, T_smeps_m);
    const double M_pot = f_dd * posT;

    const double M = smooth_min(Swe, M_pot, eps_m);
    const double Pliq = P_rain + M;

    /* snow derivatives */
    const double c = cosh(uT);
    const double sech2 = 1.0 / (c * c);
    const double dsnow_dT_tr = 0.5 * sech2 / T_smeps_m;

    const double dposT_da = dsmooth_pos_da(aT, T_smeps_m);
    const double dposT_dT_tr = -dposT_da;

    const double dMpot_dT_tr = f_dd * dposT_dT_tr;
    const double dMpot_f_dd = posT;

    const double dM_dSwe_s = dsmoothmin_da(Swe, M_pot, eps_m);
    const double dM_dSwe = dM_dSwe_s * dSwe_dRaw; /* chain to raw */
    const double dM_dMpot = dsmoothmin_db(Swe, M_pot, eps_m);

    const double dM_dT_tr = dM_dMpot * dMpot_dT_tr;
    const double dM_f_dd = dM_dMpot * dMpot_f_dd;

    const double dPliq_dSwe = dM_dSwe;
    const double dPliq_dT_tr = (-P * dsnow_dT_tr) + dM_dT_tr;
    const double dPliq_f_dd = dM_f_dd;

    /* ---------- 4) ET (smooth clip01) ---------- */
    const double ratioU = Uztw_pos / uztwm;
    const double ratioU_clip = smooth_clip01(ratioU, eps_s);
    const double dratioUclip_dy = dsmooth_clip01_dy(ratioU, eps_s);

    const double e1 = Ep * ratioU_clip;
    const double de1_dUztw = Ep * dratioUclip_dy * (dUztwpos_dUztw / uztwm);
    const double de1_duztwm = Ep * dratioUclip_dy * (-Uztw_pos / (uztwm * uztwm));

    const double ratioL = Lztw_pos / lztwm;
    const double cL = smooth_clip01(ratioL, eps_s);
    const double dratioLclip_dy = dsmooth_clip01_dy(ratioL, eps_s);

    const double dcL_dLztw = dratioLclip_dy * (dLztwpos_dLztw / lztwm);
    const double dcL_dlztwm = dratioLclip_dy * (-Lztw_pos / (lztwm * lztwm));

    const double e2 = (Ep - e1) * cL;

    const double de2_dLztw = (Ep - e1) * dcL_dLztw;
    const double de2_dlztwm = (Ep - e1) * dcL_dlztwm;
    const double de2_dUztw = -cL * de1_dUztw;
    const double de2_duztwm = -cL * de1_duztwm;

    /* ---------- 5) q12 / dlz block ---------- */
    const double q0 = lzpk * lzfpm + lzsk * lzfsm;

    const double Araw = 1.0 - Lztot / Lzm;
    const double Apos = smooth_pos(Araw, eps_s);
    const double dApos_dAraw = dsmooth_pos_da(Araw, eps_s);

    const double Afloor = 1e-12;
    const double Aeff = smooth_pos(Apos - Afloor, eps_s) + Afloor;
    const double dAeff_dApos = dsmooth_pos_da(Apos - Afloor, eps_s);

    const double logAeff = log(Aeff);
    const double Apow = exp(rexp * logAeff); /* often faster than pow */
    const double dlz = 1.0 + zperc * Apow;

    const double ddlz_dzperc = Apow;
    const double ddlz_drexp = zperc * Apow * logAeff;

    const double dAraw_dLztw = -(1.0 / Lzm);
    const double dAraw_dLzps = -(1.0 / Lzm);
    const double dAraw_dLzfs = -(1.0 / Lzm);

    const double dAraw_dLzm = (Lztot) / (Lzm * Lzm);
    const double dAraw_dlzfpm = dAraw_dLzm;
    const double dAraw_dlzfsm = dAraw_dLzm;
    const double dAraw_dlztwm = dAraw_dLzm;

    const double dApos_dLztw = dApos_dAraw * dAraw_dLztw;
    const double dApos_dLzps = dApos_dAraw * dAraw_dLzps;
    const double dApos_dLzfs = dApos_dAraw * dAraw_dLzfs;

    const double dApos_dlzfpm = dApos_dAraw * dAraw_dlzfpm;
    const double dApos_dlzfsm = dApos_dAraw * dAraw_dlzfsm;
    const double dApos_dlztwm = dApos_dAraw * dAraw_dlztwm;

    const double dApow_dApos = rexp * exp((rexp - 1.0) * logAeff) * dAeff_dApos;

    const double ddlz_dLztw = zperc * dApow_dApos * dApos_dLztw;
    const double ddlz_dLzps = zperc * dApow_dApos * dApos_dLzps;
    const double ddlz_dLzfs = zperc * dApow_dApos * dApos_dLzfs;

    const double ddlz_dlzfpm = zperc * dApow_dApos * dApos_dlzfpm;
    const double ddlz_dlzfsm = zperc * dApow_dApos * dApos_dlzfsm;
    const double ddlz_dlztwm = zperc * dApow_dApos * dApos_dlztwm;

    const double B = Uzfw_pos / uzfwm;
    const double dB_dUzfw = dUzfwpos_dUzfw / uzfwm;
    const double dB_duzfwm = -Uzfw_pos / (uzfwm * uzfwm);

    const double q12 = q0 * dlz * B;
    const double dq12_dUzfw = q0 * dlz * dB_dUzfw;

    const double dq12_dLztw = q0 * B * ddlz_dLztw;
    const double dq12_dLzps = q0 * B * ddlz_dLzps;
    const double dq12_dLzfs = q0 * B * ddlz_dLzfs;

    const double dq12_dzperc = q0 * B * ddlz_dzperc;
    const double dq12_drexp = q0 * B * ddlz_drexp;

    const double qif = uzk * B;
    const double dqif_dUzfw = uzk * dB_dUzfw;

    /* baseflows + ac */
    const double qbp = lzpk * Lzps;
    const double qbs = lzsk * Lzfs;

    const double dqbp_dLzps = lzpk;
    const double dqbs_dLzfs = lzsk;

    const double dqbp_dlzpk = Lzps;
    const double dqbs_dlzsk = Lzfs;

    const double ac = acm * (Uztw_pos / uztwm);
    const double qsx = ac * Pliq;

    const double dac_dUztw = acm * (dUztwpos_dUztw / uztwm);
    const double dac_duztwm = acm * (-Uztw_pos / (uztwm * uztwm));
    const double dac_dacm = Uztw_pos / uztwm;

    const double dqsx_dUztw = Pliq * dac_dUztw;
    const double dqsx_duztwm = Pliq * dac_duztwm;
    const double dqsx_dacm = Pliq * dac_dacm;

    const double dqsx_dSwe = ac * dPliq_dSwe;
    const double dqsx_dT_tr = ac * dPliq_dT_tr;
    const double dqsx_f_dd = ac * dPliq_f_dd;

    /* logistics */
    const double lfUztw = lf_func(Uztw, uztwm, eps, rho);
    const double lfUzfw = lf_func(Uzfw, uzfwm, eps, rho);
    const double lfLztw = lf_func(Lztw, lztwm, eps, rho);
    const double lfLzps = lf_func(Lzps, lzfpm, eps, rho);
    const double lfLzfs = lf_func(Lzfs, lzfsm, eps, rho);

    const double dlfUztw = dlf_dS(lfUztw, uztwm, rho);
    const double dlfUzfw = dlf_dS(lfUzfw, uzfwm, rho);
    const double dlfLztw = dlf_dS(lfLztw, lztwm, rho);
    const double dlfLzps = dlf_dS(lfLzps, lzfpm, rho);
    const double dlfLzfs = dlf_dS(lfLzfs, lzfsm, rho);

    const double qutof = (Pliq - qsx) * lfUztw;
    const double qufof = qutof * lfUzfw;

    const double qstof = pfree * q12 * lfLztw;
    const double prcs = 0.5 * (1.0 - pfree) * q12 + 0.5 * qstof;

    const double qsfofa = prcs * lfLzps;
    const double qsfofb = prcs * lfLzfs;

    const double qf1 = kf * Sf1;
    const double qf2 = kf * Sf2;
    const double qf3 = kf * Sf3;

    /* ---------- primal ODE ---------- */
    udot[0] = P_snow - M;
    udot[1] = Pliq - qsx - e1 - qutof;
    udot[2] = qutof - q12 - qif - qufof;
    udot[3] = pfree * q12 - e2 - qstof;
    udot[4] = prcs - qbp - qsfofa;
    udot[5] = prcs - qbs - qsfofb;
    udot[6] = qif + qsx + qufof + qsfofa + qsfofb + qbp + qbs - qf1;
    udot[7] = qf1 - qf2;
    udot[8] = qf2 - qf3;
    udot[9] = qf3;

    /* --------------------------------------------------------------- */
    /* 6) Jacobians: Jx_f (m x m) and Jth_f (m x d)                    */
    /* --------------------------------------------------------------- */
    double Jx_f[10][10] = {{0}};
    double Jth_f[10][15] = {{0}};

    /* --- key partials used below --- */
    const double dlfU_dUztwm = -dlfUztw * (Uztw / uztwm);
    const double dqutof_duztwm = -(dqsx_duztwm)*lfUztw + (Pliq - qsx) * dlfU_dUztwm;

    const double dqutof_dSwe = (dPliq_dSwe - dqsx_dSwe) * lfUztw;
    const double dqutof_dT_tr = (dPliq_dT_tr - dqsx_dT_tr) * lfUztw;
    const double dqutof_f_dd = (dPliq_f_dd - dqsx_f_dd) * lfUztw;

    const double dqutof_dUztw = -(dqsx_dUztw)*lfUztw + (Pliq - qsx) * dlfUztw;
    const double dqutof_dacm = -(dqsx_dacm)*lfUztw;

    const double dqufof_dSwe = dqutof_dSwe * lfUzfw;
    const double dqufof_dT_tr = dqutof_dT_tr * lfUzfw;
    const double dqufof_f_dd = dqutof_f_dd * lfUzfw;

    const double dqufof_dUztw = dqutof_dUztw * lfUzfw;
    const double dqufof_dUzfw = qutof * dlfUzfw;
    const double dqufof_dacm = dqutof_dacm * lfUzfw;

    const double dqufof_duztwm = dqutof_duztwm * lfUzfw;

    const double dqstof_dUzfw = pfree * dq12_dUzfw * lfLztw;
    const double dqstof_dLztw = pfree * dq12_dLztw * lfLztw + pfree * q12 * dlfLztw;
    const double dqstof_dLzps = pfree * dq12_dLzps * lfLztw;
    const double dqstof_dLzfs = pfree * dq12_dLzfs * lfLztw;

    const double dqstof_dzperc = pfree * dq12_dzperc * lfLztw;
    const double dqstof_drexp = pfree * dq12_drexp * lfLztw;
    const double dqstof_dpfree = q12 * lfLztw;

    const double dprc_dUzfw = 0.5 * (1.0 - pfree) * dq12_dUzfw + 0.5 * dqstof_dUzfw;
    const double dprc_dLztw = 0.5 * (1.0 - pfree) * dq12_dLztw + 0.5 * dqstof_dLztw;
    const double dprc_dLzps = 0.5 * (1.0 - pfree) * dq12_dLzps + 0.5 * dqstof_dLzps;
    const double dprc_dLzfs = 0.5 * (1.0 - pfree) * dq12_dLzfs + 0.5 * dqstof_dLzfs;

    const double dprc_dzperc = 0.5 * (1.0 - pfree) * dq12_dzperc + 0.5 * dqstof_dzperc;
    const double dprc_drexp = 0.5 * (1.0 - pfree) * dq12_drexp + 0.5 * dqstof_drexp;
    const double dprc_dpfree = -0.5 * q12 + 0.5 * dqstof_dpfree;

    const double dqsfpf_dUzfw = dprc_dUzfw * lfLzps;
    const double dqsfpf_dLztw = dprc_dLztw * lfLzps;
    const double dqsfpf_dLzps = dprc_dLzps * lfLzps + prcs * dlfLzps;
    const double dqsfpf_dLzfs = dprc_dLzfs * lfLzps;

    const double dqsfpf_dzperc = dprc_dzperc * lfLzps;
    const double dqsfpf_drexp = dprc_drexp * lfLzps;
    const double dqsfpf_dpfree = dprc_dpfree * lfLzps;

    const double dqsfss_dUzfw = dprc_dUzfw * lfLzfs;
    const double dqsfss_dLztw = dprc_dLztw * lfLzfs;
    const double dqsfss_dLzps = dprc_dLzps * lfLzfs;
    const double dqsfss_dLzfs = dprc_dLzfs * lfLzfs + prcs * dlfLzfs;

    const double dqsfss_dzperc = dprc_dzperc * lfLzfs;
    const double dqsfss_drexp = dprc_drexp * lfLzfs;
    const double dqsfss_dpfree = dprc_dpfree * lfLzfs;

    /* ---------- Jx_f ---------- */
    Jx_f[0][0] = -dM_dSwe;

    Jx_f[1][0] = dPliq_dSwe - dqsx_dSwe - dqutof_dSwe;
    Jx_f[1][1] = -dqsx_dUztw - de1_dUztw - dqutof_dUztw;

    Jx_f[2][0] = dqutof_dSwe - dqufof_dSwe;
    Jx_f[2][1] = dqutof_dUztw - dqufof_dUztw;
    Jx_f[2][2] = -dqif_dUzfw - dq12_dUzfw - dqufof_dUzfw;
    Jx_f[2][3] = -dq12_dLztw;
    Jx_f[2][4] = -dq12_dLzps;
    Jx_f[2][5] = -dq12_dLzfs;

    Jx_f[3][1] = -de2_dUztw;
    Jx_f[3][2] = pfree * dq12_dUzfw - dqstof_dUzfw;
    Jx_f[3][3] = pfree * dq12_dLztw - de2_dLztw - dqstof_dLztw;
    Jx_f[3][4] = pfree * dq12_dLzps - dqstof_dLzps;
    Jx_f[3][5] = pfree * dq12_dLzfs - dqstof_dLzfs;

    Jx_f[4][2] = dprc_dUzfw - dqsfpf_dUzfw;
    Jx_f[4][3] = dprc_dLztw - dqsfpf_dLztw;
    Jx_f[4][4] = dprc_dLzps - dqbp_dLzps - dqsfpf_dLzps;
    Jx_f[4][5] = dprc_dLzfs - dqsfpf_dLzfs;

    Jx_f[5][2] = dprc_dUzfw - dqsfss_dUzfw;
    Jx_f[5][3] = dprc_dLztw - dqsfss_dLztw;
    Jx_f[5][4] = dprc_dLzps - dqsfss_dLzps;
    Jx_f[5][5] = dprc_dLzfs - dqbs_dLzfs - dqsfss_dLzfs;

    Jx_f[6][0] = dqsx_dSwe + dqufof_dSwe;
    Jx_f[6][1] = dqsx_dUztw + dqufof_dUztw;
    Jx_f[6][2] = dqif_dUzfw + dqufof_dUzfw + dqsfpf_dUzfw + dqsfss_dUzfw;
    Jx_f[6][3] = dqsfpf_dLztw + dqsfss_dLztw;
    Jx_f[6][4] = dqsfpf_dLzps + dqbp_dLzps + dqsfss_dLzps;
    Jx_f[6][5] = dqsfpf_dLzfs + dqsfss_dLzfs + dqbs_dLzfs;
    Jx_f[6][6] = -kf;

    Jx_f[7][6] = kf;
    Jx_f[7][7] = -kf;

    Jx_f[8][7] = kf;
    Jx_f[8][8] = -kf;

    Jx_f[9][8] = kf;

    /* ---------- Jth_f ---------- */
    /* T_tr (14), f_dd (15) => 0-based 13,14 */
    Jth_f[0][13] = (P * dsnow_dT_tr) - dM_dT_tr;
    Jth_f[0][14] = -dM_f_dd;

    Jth_f[1][13] = dPliq_dT_tr - dqsx_dT_tr - dqutof_dT_tr;
    Jth_f[1][14] = dPliq_f_dd - dqsx_f_dd - dqutof_f_dd;

    Jth_f[2][13] = dqutof_dT_tr - dqufof_dT_tr;
    Jth_f[2][14] = dqutof_f_dd - dqufof_f_dd;

    Jth_f[6][13] = dqsx_dT_tr + dqufof_dT_tr;
    Jth_f[6][14] = dqsx_f_dd + dqufof_f_dd;

    /* 1) uzfwm */
    {
        const double dq12_duzfwm = q0 * dlz * dB_duzfwm;
        const double dqif_duzfwm = uzk * dB_duzfwm;
        const double dqstof_duzfwm = pfree * lfLztw * dq12_duzfwm;
        const double dprc_duzfwm = 0.5 * (1.0 - pfree) * dq12_duzfwm + 0.5 * dqstof_duzfwm;
        const double dqsfpf_duzfwm = dprc_duzfwm * lfLzps;
        const double dqsfss_duzfwm = dprc_duzfwm * lfLzfs;
        const double dlfUzfw_duzfwm = -dlfUzfw * (Uzfw / uzfwm);
        const double dqufof_duzfwm = qutof * dlfUzfw_duzfwm;

        Jth_f[1][0] = 0.0;
        Jth_f[2][0] = -dq12_duzfwm - dqif_duzfwm - dqufof_duzfwm;
        Jth_f[3][0] = pfree * dq12_duzfwm - dqstof_duzfwm;
        Jth_f[4][0] = dprc_duzfwm - dqsfpf_duzfwm;
        Jth_f[5][0] = dprc_duzfwm - dqsfss_duzfwm;
        Jth_f[6][0] = dqif_duzfwm + dqufof_duzfwm + dqsfpf_duzfwm + dqsfss_duzfwm;
    }

    /* 2) uztwm */
    Jth_f[1][1] = -dqsx_duztwm - de1_duztwm - dqutof_duztwm;
    Jth_f[2][1] = dqutof_duztwm - dqufof_duztwm;
    Jth_f[3][1] = -de2_duztwm;
    Jth_f[6][1] = dqsx_duztwm + dqufof_duztwm;

    /* 3) lzfpm */
    {
        const double dq0_dlzfpm = lzpk;
        const double dq12_dlzfpm = dq0_dlzfpm * dlz * B + q0 * ddlz_dlzfpm * B;
        const double dqstof_dlzfpm = pfree * lfLztw * dq12_dlzfpm;
        const double dprc_dlzfpm = 0.5 * (1.0 - pfree) * dq12_dlzfpm + 0.5 * dqstof_dlzfpm;
        const double dlfLzps_dlzfpm = -dlfLzps * (Lzps / lzfpm);
        const double dqsfpf_dlzfpm = dprc_dlzfpm * lfLzps + prcs * dlfLzps_dlzfpm;
        const double dqsfss_dlzfpm = dprc_dlzfpm * lfLzfs;

        Jth_f[2][2] = -dq12_dlzfpm;
        Jth_f[3][2] = pfree * dq12_dlzfpm - dqstof_dlzfpm;
        Jth_f[4][2] = dprc_dlzfpm - dqsfpf_dlzfpm;
        Jth_f[5][2] = dprc_dlzfpm - dqsfss_dlzfpm;
        Jth_f[6][2] = dqsfpf_dlzfpm + dqsfss_dlzfpm;
    }

    /* 4) lzfsm */
    {
        const double dq0_dlzfsm = lzsk;
        const double dq12_dlzfsm = dq0_dlzfsm * dlz * B + q0 * ddlz_dlzfsm * B;
        const double dqstof_dlzfsm = pfree * lfLztw * dq12_dlzfsm;
        const double dprc_dlzfsm = 0.5 * (1.0 - pfree) * dq12_dlzfsm + 0.5 * dqstof_dlzfsm;
        const double dlfLzfs_dlzfsm = -dlfLzfs * (Lzfs / lzfsm);
        const double dqsfpf_dlzfsm = dprc_dlzfsm * lfLzps;
        const double dqsfss_dlzfsm = dprc_dlzfsm * lfLzfs + prcs * dlfLzfs_dlzfsm;

        Jth_f[2][3] = -dq12_dlzfsm;
        Jth_f[3][3] = pfree * dq12_dlzfsm - dqstof_dlzfsm;
        Jth_f[4][3] = dprc_dlzfsm - dqsfpf_dlzfsm;
        Jth_f[5][3] = dprc_dlzfsm - dqsfss_dlzfsm;
        Jth_f[6][3] = dqsfpf_dlzfsm + dqsfss_dlzfsm;
    }

    /* 5) lztwm */
    {
        const double dq12_dlztwm = q0 * ddlz_dlztwm * B;
        const double dlfLztw_dlztwm = -dlfLztw * (Lztw / lztwm);
        const double dqstof_dlztwm = pfree * (dq12_dlztwm * lfLztw + q12 * dlfLztw_dlztwm);
        const double dprc_dlztwm = 0.5 * (1.0 - pfree) * dq12_dlztwm + 0.5 * dqstof_dlztwm;
        const double dqsfpf_dlztwm = dprc_dlztwm * lfLzps;
        const double dqsfss_dlztwm = dprc_dlztwm * lfLzfs;

        Jth_f[2][4] = -dq12_dlztwm;
        Jth_f[3][4] = pfree * dq12_dlztwm - de2_dlztwm - dqstof_dlztwm;
        Jth_f[4][4] = dprc_dlztwm - dqsfpf_dlztwm;
        Jth_f[5][4] = dprc_dlztwm - dqsfss_dlztwm;
        Jth_f[6][4] = dqsfpf_dlztwm + dqsfss_dlztwm;
    }

    /* 6) zperc */
    Jth_f[2][5] = -dq12_dzperc;
    Jth_f[3][5] = pfree * dq12_dzperc - dqstof_dzperc;
    Jth_f[4][5] = dprc_dzperc - dqsfpf_dzperc;
    Jth_f[5][5] = dprc_dzperc - dqsfss_dzperc;
    Jth_f[6][5] = dqsfpf_dzperc + dqsfss_dzperc;

    /* 7) rexp */
    Jth_f[2][6] = -dq12_drexp;
    Jth_f[3][6] = pfree * dq12_drexp - dqstof_drexp;
    Jth_f[4][6] = dprc_drexp - dqsfpf_drexp;
    Jth_f[5][6] = dprc_drexp - dqsfss_drexp;
    Jth_f[6][6] = dqsfpf_drexp + dqsfss_drexp;

    /* 8) uzk */
    Jth_f[2][7] = -B;
    Jth_f[6][7] = B;

    /* 9) pfree */
    Jth_f[3][8] = q12 - dqstof_dpfree;
    Jth_f[4][8] = dprc_dpfree - dqsfpf_dpfree;
    Jth_f[5][8] = dprc_dpfree - dqsfss_dpfree;
    Jth_f[6][8] = dqsfpf_dpfree + dqsfss_dpfree;

    /* 10) lzpk */
    {
        const double q12_lzpk = lzfpm * dlz * B;
        const double qstof_lzpk = pfree * q12_lzpk * lfLztw;
        const double prc_lzpk = 0.5 * (1.0 - pfree) * q12_lzpk + 0.5 * qstof_lzpk;
        const double qsfp_lzpk = prc_lzpk * lfLzps;
        const double qsfs_lzpk = prc_lzpk * lfLzfs;

        Jth_f[2][9] = -q12_lzpk;
        Jth_f[3][9] = pfree * q12_lzpk - qstof_lzpk;
        Jth_f[4][9] = prc_lzpk - dqbp_dlzpk - qsfp_lzpk;
        Jth_f[5][9] = prc_lzpk - qsfs_lzpk;
        Jth_f[6][9] = qsfp_lzpk + qsfs_lzpk + dqbp_dlzpk;
    }

    /* 11) lzsk */
    {
        const double q12_lzsk = lzfsm * dlz * B;
        const double qstof_lzsk = pfree * q12_lzsk * lfLztw;
        const double prc_lzsk = 0.5 * (1.0 - pfree) * q12_lzsk + 0.5 * qstof_lzsk;
        const double qsfp_lzsk = prc_lzsk * lfLzps;
        const double qsfs_lzsk = prc_lzsk * lfLzfs;

        Jth_f[2][10] = -q12_lzsk;
        Jth_f[3][10] = pfree * q12_lzsk - qstof_lzsk;
        Jth_f[4][10] = prc_lzsk - qsfp_lzsk;
        Jth_f[5][10] = prc_lzsk - dqbs_dlzsk - qsfs_lzsk;
        Jth_f[6][10] = qsfp_lzsk + qsfs_lzsk + dqbs_dlzsk;
    }

    /* 12) acm */
    Jth_f[1][11] = -dqsx_dacm - dqutof_dacm;
    Jth_f[2][11] = dqutof_dacm - dqufof_dacm;
    Jth_f[6][11] = dqsx_dacm + dqufof_dacm;

    /* 13) kf */
    Jth_f[6][12] = -Sf1;
    Jth_f[7][12] = Sf1 - Sf2;
    Jth_f[8][12] = Sf2 - Sf3;
    Jth_f[9][12] = Sf3;

    /* ---------- sensitivity update (column-major Smat) ---------- */
    /* dSdt(:,j) = Jx_f * S(:,j) + Jth_f(:,j) */
    for (int j = 0; j < d; ++j) {
        for (int i = 0; i < m; ++i) {
            dSdt[i + m * j] = Jth_f[i][j];
        }
    }

    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < m; ++k) {
            const double coefficient = Jx_f[i][k];
            if (coefficient == 0.0) {
                continue;
            }
            for (int j = 0; j < d; ++j) {
                dSdt[i + m * j] += coefficient * Smat[k + m * j];
            }
        }
    }
}

namespace sage_sacsma {
bool run_into(int ns,
              const double* z0,
              std::size_t nz0,
              const Forcing& forcing,
              const Params& p,
              const Options& opt,
              bool mem,
              int ipr,
              bool needJ,
              const OutputView& out)
{
    const int d = 15;
    const int m = 10;
    const int nvar = m * (d + 1);
    const int nt = ns + 1;
    if (ns < 1 || !z0 || nz0 != (std::size_t)nvar) {
        return true;
    }
    if (!forcing.P || !forcing.Ep || !forcing.T || forcing.n < (std::size_t)ns) {
        return true;
    }
    if (ipr < 1) {
        ipr = 1;
    }
    if (ipr > ns + 1) {
        ipr = ns + 1;
    }
    const int n_q = (!mem && ipr <= ns) ? (ns - ipr + 1) : 0;
    const std::size_t zr = mem ? (std::size_t)nt : 1u;
    if (!out.Z || out.zrows != zr || out.zcols != (std::size_t)nvar) {
        return true;
    }
    if (!mem && n_q > 0) {
        if (!out.q || out.nq != (std::size_t)n_q) {
            return true;
        }
        if (needJ && (!out.J || out.nj != (std::size_t)d)) {
            return true;
        }
    }
    const double* P = forcing.P;
    const double* Ep = forcing.Ep;
    const double* T = forcing.T;
    const double hin = opt.InitStep, hmax_ = opt.MaxStep, hmin_ = opt.MinStep,
                 reltol = opt.RelTol, abstol = opt.AbsTol, order = opt.Order;
    const int maxiter = opt.maxiter;
    const double uzfwm = p.uzfwm;
    const double uztwm = p.uztwm;
    const double lzfpm = p.lzfpm;
    const double lzfsm = p.lzfsm;
    const double lztwm = p.lztwm;
    const double zperc = p.zperc;
    const double rexp = p.rexp;
    const double uzk = p.uzk;
    const double pfree = p.pfree;
    const double lzpk = p.lzpk;
    const double lzsk = p.lzsk;
    const double acm = p.acm;
    const double kf = p.kf;
    const double T_tr = p.T_tr;
    const double f_dd = p.f_dd;
    const double T_sm = p.T_sm;
    const double eps_m = p.eps_m;
    const double eps_s = p.eps_s;
    const double eps = p.eps;
    const double rho = p.rho;

    std::vector<double> LTE(nvar), ztmp(nvar), w(nvar), zdotE(nvar), zE(nvar), zdot(nvar),
        zcur(z0, z0 + nvar), prev;
    std::fill(out.Z, out.Z + out.zrows * out.zcols, 0.0);
    if (!mem && n_q > 0) {
        std::fill(out.q, out.q + out.nq, 0.0);
        if (needJ) {
            std::fill(out.J, out.J + out.nq * out.nj, 0.0);
        }
    }
    int last_stored_row = 0;
    if (mem) {
        for (int j = 0; j < nvar; ++j) {
            out.Z[(std::size_t)nt * j] = zcur[j];
        }
    } else {
        prev.assign(z0, z0 + nvar);
    }
    bool fail = false;
    int iterCount = 0, k_out = 0;
    double hCarry = std::max(hmin_, std::min(hin, hmax_));
    for (int ss = 1; ss <= ns; ++ss) {
        const double t1 = double(ss - 1), t2 = double(ss);
        double tcur = t1, h = std::min(hCarry, t2 - t1);
        while (tcur < t2) {
            std::copy(zcur.begin(), zcur.end(), ztmp.begin());
            rk2(nvar,
                m,
                d,
                h,
                ztmp.data(),
                LTE.data(),
                zdotE.data(),
                zE.data(),
                zdot.data(),
                P[ss - 1],
                Ep[ss - 1],
                T[ss - 1],
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
                rho);
            for (int ii = 0; ii < nvar; ++ii) {
                if (!std::isfinite(ztmp[ii]) || std::abs(ztmp[ii]) > 1e12) {
                    fail = true;
                    break;
                }
            }
            if (fail) {
                break;
            }
            double sumWL = 0.0;
            for (int ii = 0; ii < nvar; ++ii) {
                w[ii] = 1.0 / (reltol * std::abs(ztmp[ii]) + abstol);
                const double wl = w[ii] * LTE[ii];
                sumWL += wl * wl;
            }
            const double wrms = std::sqrt(sumWL / double(nvar));
            const bool accepted = (wrms <= 1.0) || (h <= hmin_);
            if (accepted) {
                zcur.swap(ztmp);
                tcur += h;
                iterCount = 0;
            } else {
                ++iterCount;
            }
            double fac =
                (wrms <= 0.0)
                    ? 5.0
                    : std::min(5.0, std::max(0.2, 0.9 * std::pow(wrms, -1.0 / order)));
            double hNext = std::max(hmin_, std::min(h * fac, hmax_));
            hCarry = hNext;
            const double tleft = t2 - tcur;
            const double tTol = 10.0 * std::numeric_limits<double>::epsilon() *
                                std::max(1.0, std::max(std::abs(t2), std::abs(tcur)));
            if (tleft <= tTol) {
                tcur = t2;
                break;
            }
            h = std::min(hNext, tleft);
            if (h < tTol) {
                fail = true;
                break;
            }
            if (iterCount >= maxiter) {
                fail = true;
                break;
            }
        }
        if (fail) {
            break;
        }
        if (mem) {
            for (int j = 0; j < nvar; ++j) {
                out.Z[(std::size_t)ss + (std::size_t)nt * j] = zcur[j];
            }
            last_stored_row = ss;
        } else {
            for (int j = 0; j < nvar; ++j) {
                out.Z[(std::size_t)j] = zcur[j];
            }
            if (ss >= ipr && out.q) {
                out.q[k_out] = zcur[m - 1] - prev[m - 1];
                if (needJ && out.J) {
                    for (int j = 0; j < d; ++j) {
                        const int idx = (j + 2) * m - 1;
                        out.J[(std::size_t)k_out + (std::size_t)n_q * j] =
                            zcur[idx] - prev[idx];
                    }
                }
                ++k_out;
            }
            std::copy(zcur.begin(), zcur.end(), prev.begin());
        }
    }
    if (fail && mem) {
        for (int row = last_stored_row + 1; row < nt; ++row) {
            for (int j = 0; j < nvar; ++j) {
                out.Z[(std::size_t)row + (std::size_t)nt * j] =
                    out.Z[(std::size_t)last_stored_row + (std::size_t)nt * j];
            }
        }
    }
    return fail;
}
} // namespace sage_sacsma
