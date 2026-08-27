/*
 * cfe_nwm.cpp
 *
 * Conceptual CFE-NWM rainfall-runoff model with an adaptive-step
 * explicit Runge-Kutta integrator. This file contains the MATLAB-independent
 * native numerical core shared by crr_cfe_nwm and crr_model_mex.
 *
 * Written by Jasper A. Vrugt, January 2025.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#include "cfe_nwm.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <cstdio>
#include <cstdlib>

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
                double s_max,
                double s_fc,
                double s_wp,
                double k_sch,
                double a1,
                double k_perc,
                double lf_thr,
                double a2,
                double k_lf,
                double g_max,
                double c_gw,
                double mm,
                double k_nsh,
                double* giuh,
                double T_tr,
                double f_dd,
                int K,
                int L,
                double dT,
                double T_sm,
                double eps_m,
                double eps,
                double rho);

static void cfe_nwm_aug_ode(double* z,
                            double* zdot,
                            double P,
                            double Ep,
                            double T,
                            double s_max,
                            double s_fc,
                            double s_wp,
                            double k_sch,
                            double a1,
                            double k_perc,
                            double lf_thr,
                            double a2,
                            double k_lf,
                            double g_max,
                            double c_gw,
                            double mm,
                            double k_nsh,
                            double* giuh,
                            double T_tr,
                            double f_dd,
                            int K,
                            int L,
                            double dT,
                            double T_sm,
                            double eps_m,
                            double eps,
                            double rho,
                            int m,
                            int d);

static void cfe_nwm_odefcn(const double* u,
                           double* udot,
                           const double* S,
                           double* dSdt,
                           double P,
                           double Ep,
                           double T,
                           double s_max,
                           double s_fc,
                           double s_wp,
                           double k_sch,
                           double a1,
                           double k_perc,
                           double lf_thr,
                           double a2,
                           double k_lf,
                           double g_max,
                           double c_gw,
                           double mm,
                           double k_nsh,
                           double* giuh,
                           double T_tr,
                           double f_dd,
                           int K,
                           int L,
                           double dT,
                           double T_sm,
                           double eps_m,
                           double eps,
                           double rho,
                           int m,
                           int d);

static inline double smooth_pos(double a, double ep)
{
    return 0.5 * (a + sqrt(a * a + ep * ep));
}

static inline double dsmooth_pos_da(double a, double ep)
{
    return 0.5 * (1.0 + a / sqrt(a * a + ep * ep));
}

static inline double smooth_min2(double a, double b, double ep)
{
    const double d = a - b;
    return 0.5 * (a + b - sqrt(d * d + ep * ep));
}

static inline double dmin2_da(double a, double b, double ep)
{
    const double d = a - b;
    const double r = sqrt(d * d + ep * ep);
    return 0.5 * (1.0 - d / r);
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
                double s_max,
                double s_fc,
                double s_wp,
                double k_sch,
                double a1,
                double k_perc,
                double lf_thr,
                double a2,
                double k_lf,
                double g_max,
                double c_gw,
                double mm,
                double k_nsh,
                double* giuh,
                double T_tr,
                double f_dd,
                int K,
                int L,
                double dT,
                double T_sm,
                double eps_m,
                double eps,
                double rho)
{
    /* Euler solution */
    cfe_nwm_aug_ode(z,
                    zdotE,
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
                    d);
    for (int i = 0; i < nvar; ++i) {
        zE[i] = z[i] + h * zdotE[i];
    }
    /* Heun solution */
    cfe_nwm_aug_ode(zE,
                    zdot,
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

static void cfe_nwm_aug_ode(double* z,
                            double* zdot,
                            double P,
                            double Ep,
                            double T,
                            double s_max,
                            double s_fc,
                            double s_wp,
                            double k_sch,
                            double a1,
                            double k_perc,
                            double lf_thr,
                            double a2,
                            double k_lf,
                            double g_max,
                            double c_gw,
                            double mm,
                            double k_nsh,
                            double* giuh,
                            double T_tr,
                            double f_dd,
                            int K,
                            int L,
                            double dT,
                            double T_sm,
                            double eps_m,
                            double eps,
                            double rho,
                            int m,
                            int d)
{
    double* u = z;        /* m */
    double* Smat = z + m; /* m*d */

    double* udot = zdot;     /* m */
    double* dSdt = zdot + m; /* m*d */

    /* Call your CFE_NWM ODE + sensitivity routine */
    cfe_nwm_odefcn(u,
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
                   d);
}

static void cfe_nwm_odefcn(const double* u,
                           double* udot,
                           const double* Smat,
                           double* dSdt,
                           double P,
                           double Ep,
                           double T,
                           double s_max,
                           double s_fc,
                           double s_wp,
                           double k_sch,
                           double a1,
                           double k_perc,
                           double lf_thr,
                           double a2,
                           double k_lf,
                           double g_max,
                           double c_gw,
                           double mm,
                           double k_nsh,
                           double* giuh,
                           double T_tr,
                           double f_dd,
                           int K,
                           int L,
                           double dT,
                           double T_sm,
                           double eps_m,
                           double eps,
                           double rho,
                           int m,
                           int d)
{
    (void)eps;
    (void)rho;

    /* ------------------------------------------------------------------
     * State order:
     *   u[0]=Swe, u[1]=S, u[2]=G,
     *   u[3..3+L-1]=GIUH queue Q (mm/T),
     *   u[3+L..3+L+K-1]=Nash storages F (mm),
     *   u[m-1]=Qcum
     *
     * Sensitivities: Smat and dSdt are m×d column-major: Smat[i + m*j]
     * ------------------------------------------------------------------ */

    /* --- unpack states --- */
    const double Swe_u = u[0];
    const double S = u[1];
    const double G = u[2];
    const double* Q = u + 3;     /* length L */
    const double* F = u + 3 + L; /* length K */

    /* ------------------------------------------------------------------
     * 0) SWE smooth positivity (REPLACES fmax(u[0],0))
     * ------------------------------------------------------------------ */
    const double Swe = smooth_pos(Swe_u, eps_m);
    const double dSwe_dSweu = dsmooth_pos_da(Swe_u, eps_m);

    /* ------------------------------------------------------------------
     * 1) Snow module (smooth partition + smooth melt min)
     * ------------------------------------------------------------------ */
    const double T_smeps_m = fmax(T_sm, eps_m);
    const double uT = (T - T_tr) / T_smeps_m;

    const double tanhu = tanh(uT);
    const double coshu = cosh(uT);
    const double sech2 = 1.0 / (coshu * coshu);

    const double snow_fr = 0.5 * (1.0 - tanhu);
    const double rain_fr = 1.0 - snow_fr;

    const double P_snow = P * snow_fr;
    const double P_rain = P * rain_fr;

    const double aT = (T - T_tr);
    const double denom_pos = sqrt(aT * aT + T_smeps_m * T_smeps_m);
    const double posT = 0.5 * (aT + denom_pos);
    const double dpos_daT = 0.5 * (1.0 + aT / denom_pos);
    const double dpos_dT_tr = -dpos_daT;

    const double M_pot = f_dd * posT;

    /* smooth min(Swe, M_pot) */
    const double dxy = (Swe - M_pot);
    const double sqrtm = sqrt(dxy * dxy + eps_m * eps_m);
    const double M = 0.5 * (Swe + M_pot - sqrtm);

    /* derivatives wrt SWE (smoothed) and M_pot */
    const double dM_dSwe = 0.5 * (1.0 - dxy / sqrtm);
    const double dM_dMpot = 0.5 * (1.0 + dxy / sqrtm);

    /* chain back to RAW Swe_u */
    const double dM_dSweu = dM_dSwe * dSwe_dSweu;

    const double dsnow_dT_tr = 0.5 * sech2 * (1.0 / T_smeps_m);
    const double drain_dT_tr = -dsnow_dT_tr;

    const double dP_snow_dT_tr = P * dsnow_dT_tr;
    const double dP_rain_dT_tr = P * drain_dT_tr;

    const double dMpot_dT_tr = f_dd * dpos_dT_tr;
    const double dMpot_df_dd = posT;

    const double dM_dT_tr = dM_dMpot * dMpot_dT_tr;
    const double dM_df_dd = dM_dMpot * dMpot_df_dd;

    const double Pliq = P_rain + M;

    /* chain derivatives wrt RAW Swe_u */
    const double dPliq_dSweu = dM_dSweu;
    const double dPliq_dT_tr = dP_rain_dT_tr + dM_dT_tr;
    const double dPliq_df_dd = dM_df_dd;

    /* ------------------------------------------------------------------
     * 2) Rainfall-first ET (SMOOTH MIN)
     * ------------------------------------------------------------------ */
    const double E_r = smooth_min2(Pliq, Ep, eps_m);
    const double dEr_dPliq = dmin2_da(Pliq, Ep, eps_m);

    const double P_e = Pliq - E_r;
    const double Ep_star = Ep - E_r;

    const double dPe_dPliq = 1.0 - dEr_dPliq;
    const double dEpstar_dPliq = -dEr_dPliq;

    const double dPe_dSweu = dPe_dPliq * dPliq_dSweu;
    const double dPe_dT_tr = dPe_dPliq * dPliq_dT_tr;
    const double dPe_df_dd = dPe_dPliq * dPliq_df_dd;

    const double dEpstar_dSweu = dEpstar_dPliq * dPliq_dSweu;
    const double dEpstar_dT_tr = dEpstar_dPliq * dPliq_dT_tr;
    const double dEpstar_df_dd = dEpstar_dPliq * dPliq_df_dd;

    /* ------------------------------------------------------------------
     * 3) Schaake partitioning (FULLY SMOOTHED)
     * ------------------------------------------------------------------ */
    double I = 0.0, R_s = 0.0;
    double dI_dS = 0.0, dI_dSmax = 0.0, dI_dKsch = 0.0, dI_dPe = 0.0;
    double dRs_dS = 0.0, dRs_dSmax = 0.0, dRs_dKsch = 0.0, dRs_dPe = 0.0;

    const double D_soil = s_max - S;
    const double eTerm = exp(-k_sch * dT);

    /* smooth activation for Pe and deficit */
    const double eps_pe = 1e-6;
    const double eps_D = 1e-6;

    const double Pe_pos = smooth_pos(P_e, eps_pe);
    const double dPePos_dPe = dsmooth_pos_da(P_e, eps_pe);

    const double D_pos = smooth_pos(D_soil, eps_D);
    const double dDpos_dD = dsmooth_pos_da(D_soil, eps_D);

    const double Ic = D_pos * (1.0 - eTerm); /* mm */
    const double Px = Pe_pos * dT;           /* mm */

    const double denPI = (Px + Ic);
    const double I_amt = Px * (Ic / denPI);
    const double R_amt = Px - I_amt;

    I = I_amt / dT;
    R_s = R_amt / dT;

    /* dI_amt/dPx and dI_amt/dIc */
    const double dI_dPx_amt = (Ic * Ic) / (denPI * denPI);
    const double dI_dIc_amt = (Px * Px) / (denPI * denPI);

    /* Ic derivatives */
    const double dIc_dDpos = (1.0 - eTerm);
    const double dIc_dS = dIc_dDpos * dDpos_dD * (-1.0);
    const double dIc_dSmax = dIc_dDpos * dDpos_dD * (+1.0);
    const double dIc_dKsch = D_pos * (dT * eTerm);

    dI_dS = (dI_dIc_amt * dIc_dS) / dT;
    dI_dSmax = (dI_dIc_amt * dIc_dSmax) / dT;
    dI_dKsch = (dI_dIc_amt * dIc_dKsch) / dT;

    dRs_dS = -dI_dS;
    dRs_dSmax = -dI_dSmax;
    dRs_dKsch = -dI_dKsch;

    /* chain through Pe_pos(Pe): dI/dPe and dRs/dPe */
    dI_dPe = (dI_dPx_amt)*dPePos_dPe;
    dRs_dPe = (1.0 - dI_dPx_amt) * dPePos_dPe;

    /* chain from snow variables to I and Rs through Pe */
    const double dI_dSweu = dI_dPe * dPe_dSweu;
    const double dI_dT_tr = dI_dPe * dPe_dT_tr;
    const double dI_df_dd = dI_dPe * dPe_df_dd;

    const double dRs_dSweu = dRs_dPe * dPe_dSweu;
    const double dRs_dT_tr = dRs_dPe * dPe_dT_tr;
    const double dRs_df_dd = dRs_dPe * dPe_df_dd;

    /* ------------------------------------------------------------------
     * 4) GIUH queue (avoid per-call malloc if you want; keep here for now)
     * ------------------------------------------------------------------ */
    double* dQdt =
        static_cast<double*>(std::calloc(static_cast<std::size_t>(L), sizeof(double)));
    for (int i = 0; i < L; ++i) {
        const double Qi = Q[i];
        const double Qshift = (i < L - 1) ? Q[i + 1] : 0.0;
        dQdt[i] = -Qi + Qshift + giuh[i] * R_s;
    }
    const double Q_giuh = Q[0];

    /* ------------------------------------------------------------------
     * 5) Soil outlets: percolation (SMOOTHED) + lateral flow (your smoothed)
     * ------------------------------------------------------------------ */

    /* ---- Percolation (SMOOTHED activation + cap) ---- */
    double Q_perc = 0.0, dQperc_dS = 0.0, dQperc_dSmax = 0.0, dQperc_dSfc = 0.0,
           dQperc_da1 = 0.0, dQperc_dkperc = 0.0;

    const double den1 = (s_max - s_fc); /* must be >0 by bounds */
    const double above_raw = S - s_fc;
    const double eps_ab = 1e-6;

    const double above_pos = smooth_pos(above_raw, eps_ab);
    const double dAbove_dS = dsmooth_pos_da(above_raw, eps_ab);
    const double dAbove_dSfc = -dAbove_dS;

    const double r1 = above_pos / den1;
    const double eps_r1 = 1e-12;
    const double r1e = r1 + eps_r1;

    const double Qunc = k_perc * pow(r1e, a1);

    /* Qunc derivatives */
    const double dQunc_dkperc = pow(r1e, a1);
    const double dQunc_da1 = Qunc * log(r1e);
    const double dQunc_dr1 = k_perc * a1 * pow(r1e, a1 - 1.0);

    const double dr1_dS = dAbove_dS / den1;
    const double dr1_dSmax = -above_pos / (den1 * den1);
    const double dr1_dSfc = (dAbove_dSfc * den1 + above_pos) / (den1 * den1); /* quotient */

    const double dQunc_dS1 = dQunc_dr1 * dr1_dS;
    const double dQunc_dSmax1 = dQunc_dr1 * dr1_dSmax;
    const double dQunc_dSfc = dQunc_dr1 * dr1_dSfc;

    /* cap at above_pos via smooth min */
    const double eps_cap_perc = 1e-6;
    Q_perc = smooth_min2(Qunc, above_pos, eps_cap_perc);

    const double w_unc_perc = dmin2_da(Qunc, above_pos, eps_cap_perc); /* dQperc/dQunc */
    const double w_cap_perc = 1.0 - w_unc_perc;                        /* dQperc/dAbove */

    dQperc_dkperc = w_unc_perc * dQunc_dkperc;
    dQperc_da1 = w_unc_perc * dQunc_da1;

    dQperc_dS = w_unc_perc * dQunc_dS1 + w_cap_perc * dAbove_dS;
    dQperc_dSmax = w_unc_perc * dQunc_dSmax1; /* above_pos has no Smax */
    dQperc_dSfc = w_unc_perc * dQunc_dSfc + w_cap_perc * dAbove_dSfc;

    /* ---- Lateral flow (your existing SMOOTH block) ---- */
    double Q_lf = 0.0, dQlf_dS = 0.0, dQlf_dSmax = 0.0, dQlf_dlfthr = 0.0, dQlf_da2 = 0.0,
           dQlf_dklf = 0.0;

    const double eps_on = 1e-6;  /* mm */
    const double eps_cap = 1e-6; /* mm/T */
    const double eps_r = 1e-12;  /* dimensionless */

    const double a2raw = S - lf_thr; /* mm */
    const double R_on = sqrt(a2raw * a2raw + eps_on * eps_on);
    const double above2 = 0.5 * (a2raw + R_on); /* mm */

    const double dabove2_dS = 0.5 * (1.0 + a2raw / R_on);
    const double dabove2_dlfthr = -dabove2_dS;

    const double den2 = (s_max - lf_thr); /* mm */
    const double r2 = above2 / den2;

    const double dr2_dS = dabove2_dS / den2;
    const double dr2_dSmax = -above2 / (den2 * den2);
    const double dr2_dlfthr = (dabove2_dlfthr * den2 + above2) / (den2 * den2);

    const double r2e = r2 + eps_r;
    const double r2e_a2 = pow(r2e, a2);
    const double Q_unc = k_lf * r2e_a2;

    const double dQunc_dklf = r2e_a2;
    const double dQunc_da2 = Q_unc * log(r2e);
    const double dQunc_dr2 = k_lf * a2 * pow(r2e, a2 - 1.0);

    const double dQunc_dS2 = dQunc_dr2 * dr2_dS;
    const double dQunc_dSmax2 = dQunc_dr2 * dr2_dSmax;
    const double dQunc_dlfthr = dQunc_dr2 * dr2_dlfthr;

    const double cap2 = above2 - Q_perc; /* mm/T */
    const double dcap2_dS = dabove2_dS - dQperc_dS;
    const double dcap2_dSmax = -dQperc_dSmax; /* above2 has no Smax */
    const double dcap2_dlfthr = dabove2_dlfthr;

    const double Dcap = Q_unc - cap2;
    const double Rcap = sqrt(Dcap * Dcap + eps_cap * eps_cap);
    Q_lf = 0.5 * (Q_unc + cap2 - Rcap);

    const double w_unc = 0.5 * (1.0 - Dcap / Rcap); /* dQlf/dQ_unc */
    const double w_cap = 1.0 - w_unc;               /* dQlf/dcap2 */

    dQlf_dklf = w_unc * dQunc_dklf;
    dQlf_da2 = w_unc * dQunc_da2;
    dQlf_dS = w_unc * dQunc_dS2 + w_cap * dcap2_dS;
    dQlf_dSmax = w_unc * dQunc_dSmax2 + w_cap * dcap2_dSmax;
    dQlf_dlfthr = w_unc * dQunc_dlfthr + w_cap * dcap2_dlfthr;

    /* ------------------------------------------------------------------
     * 6) Soil ET extraction (UNCHANGED piecewise; can be smoothed later)
     * ------------------------------------------------------------------ */
    double E_s = 0.0, dEs_dS = 0.0, dEs_dSfc = 0.0, dEs_dSwp = 0.0;
    double dEs_dSweu = 0.0, dEs_dT_tr = 0.0, dEs_df_dd = 0.0;
    double dEs_dEpstar = 0.0;

    if (Ep_star > 0.0) {
        if (S <= s_wp) {
            E_s = 0.0;
            dEs_dEpstar = 0.0;
        } else if (S >= s_fc) {
            E_s = Ep_star;
            dEs_dEpstar = 1.0;
        } else {
            const double den = (s_fc - s_wp);
            const double phi = (S - s_wp) / den;
            E_s = Ep_star * phi;

            dEs_dS = Ep_star / den;
            dEs_dSfc = Ep_star * (-(S - s_wp) / (den * den));
            dEs_dSwp = Ep_star * ((S - s_fc) / (den * den));
            dEs_dEpstar = phi;
        }
    }

    dEs_dSweu = dEs_dEpstar * dEpstar_dSweu;
    dEs_dT_tr = dEs_dEpstar * dEpstar_dT_tr;
    dEs_df_dd = dEs_dEpstar * dEpstar_df_dd;

    /* ------------------------------------------------------------------
     * 7) Groundwater flux (unchanged)
     * ------------------------------------------------------------------ */
    const double expterm = exp(mm * G / g_max);
    const double flux_exp = expterm - 1.0;
    const double Q_gw = c_gw * flux_exp;

    const double dQgw_dG = c_gw * expterm * (mm / g_max);
    const double dQgw_dC = flux_exp;
    const double dQgw_dmm = c_gw * expterm * (G / g_max);
    const double dQgw_dGmax = c_gw * expterm * (-mm * G / (g_max * g_max));

    /* ------------------------------------------------------------------
     * 8) Nash routing
     * ------------------------------------------------------------------ */
    double* dF =
        static_cast<double*>(std::calloc(static_cast<std::size_t>(K), sizeof(double)));
    for (int k = 0; k < K; ++k) {
        const double Fk = F[k];
        if (k == 0) {
            dF[k] = Q_lf - k_nsh * Fk;
        } else {
            dF[k] = k_nsh * F[k - 1] - k_nsh * Fk;
        }
    }
    const double Q_fast = k_nsh * F[K - 1];

    /* ------------------------------------------------------------------
     * 9) Assemble udot
     * ------------------------------------------------------------------ */
    for (int i = 0; i < m; ++i) {
        udot[i] = 0.0;
    }

    udot[0] = P_snow - M;
    udot[1] = I - Q_perc - Q_lf - E_s;
    udot[2] = Q_perc - Q_gw;

    for (int i = 0; i < L; ++i) {
        udot[3 + i] = dQdt[i];
    }
    for (int k = 0; k < K; ++k) {
        udot[3 + L + k] = dF[k];
    }

    udot[m - 1] = Q_fast + Q_giuh + Q_gw;

    /* ------------------------------------------------------------------
     * 10) Build Jx_f (m×m), Jp_f (m×15), Jth_f (m×d), then dSdt = Jx*S + Jth
     * ------------------------------------------------------------------ */
    //    std::vector<double> Jx_f_(m*m, 0.0);
    //    std::vector<double> Jth_f_(m*d, 0.0);
    //    std::vector<double> Jp_f_(m*d, 0.0);

    //    #define Jx_f(i,k)   (Jx_f_[(i)*(m) + (k)])     /* row-major */
    //    #define Jth_f(i,j)  (Jth_f_[(i) + (m)*(j)])    /* col-major */
    //    #define Jp_f(i,j)  (Jp_f_[(i) + (m)*(j)])      /* col-major */

    double* Jx_f = static_cast<double*>(std::calloc(
        static_cast<std::size_t>(m) * static_cast<std::size_t>(m), sizeof(double)));
    double* Jp_f = static_cast<double*>(std::calloc(
        static_cast<std::size_t>(m) * static_cast<std::size_t>(d), sizeof(double)));
    double* Jth_f = static_cast<double*>(std::calloc(
        static_cast<std::size_t>(m) * static_cast<std::size_t>(d), sizeof(double)));

    const int iSwe = 0, iSoil = 1, iGW = 2, iq1 = 3, iF1 = 3 + L, iQcum = m - 1;

    /* SWE' = P_snow(T_tr) - M(Swe_u,T_tr,f_dd)  -> MUST use dM_dSweu */
    Jx_f[iSwe * m + iSwe] = -dM_dSweu;

    /* Soil: S' = I - Qperc - Qlf - Es */
    Jx_f[iSoil * m + iSwe] = dI_dSweu - dEs_dSweu;
    Jx_f[iSoil * m + iSoil] = dI_dS - dQperc_dS - dQlf_dS - dEs_dS;

    /* GW: G' = Qperc - Qgw */
    Jx_f[iGW * m + iSoil] = dQperc_dS;
    Jx_f[iGW * m + iGW] = -dQgw_dG;

    /* GIUH block */
    for (int i = 0; i < L; ++i) {
        const int row = iq1 + i;
        const int col = iq1 + i;
        Jx_f[row * m + col] = -1.0;
        if (i < L - 1) {
            Jx_f[row * m + (col + 1)] = 1.0;
        }
        Jx_f[row * m + iSoil] = giuh[i] * dRs_dS;
        Jx_f[row * m + iSwe] = giuh[i] * dRs_dSweu;
    }

    /* Nash cascade */
    Jx_f[iF1 * m + iSoil] = dQlf_dS;
    Jx_f[iF1 * m + iF1] = -k_nsh;
    for (int k = 1; k < K; ++k) {
        const int row = iF1 + k;
        Jx_f[row * m + (row - 1)] = k_nsh;
        Jx_f[row * m + row] = -k_nsh;
    }

    /* Qcum' = Q_giuh + Q_fast + Q_gw */
    Jx_f[iQcum * m + iq1] = 1.0;
    Jx_f[iQcum * m + (iF1 + K - 1)] = k_nsh;
    Jx_f[iQcum * m + iGW] = dQgw_dG;

    /* ---- Jp_f (physical params p=15, 0-based) ---- */
    const int pSmax = 0, pSfc = 1, pSwp = 2, pKsch = 3, pa1 = 4, pkperc = 5, pLfthr = 6;
    const int pa2 = 7, pkLf = 8, pGmax = 9, pCgw = 10, pmm = 11, pkNsh = 12, pT_tr = 13,
              pf_dd = 14;

    /* SWE row */
    Jp_f[iSwe + m * pT_tr] = dP_snow_dT_tr - dM_dT_tr;
    Jp_f[iSwe + m * pf_dd] = -dM_df_dd;

    /* Soil row */
    Jp_f[iSoil + m * pSmax] += dI_dSmax;
    Jp_f[iSoil + m * pKsch] += dI_dKsch;

    Jp_f[iSoil + m * pa1] -= dQperc_da1;
    Jp_f[iSoil + m * pkperc] -= dQperc_dkperc;
    Jp_f[iSoil + m * pSmax] -= dQperc_dSmax;
    Jp_f[iSoil + m * pSfc] -= dQperc_dSfc;

    Jp_f[iSoil + m * pSmax] -= dQlf_dSmax;
    Jp_f[iSoil + m * pLfthr] -= dQlf_dlfthr;
    Jp_f[iSoil + m * pa2] -= dQlf_da2;
    Jp_f[iSoil + m * pkLf] -= dQlf_dklf;

    Jp_f[iSoil + m * pSfc] -= dEs_dSfc;
    Jp_f[iSoil + m * pSwp] -= dEs_dSwp;

    Jp_f[iSoil + m * pT_tr] += dI_dT_tr - dEs_dT_tr;
    Jp_f[iSoil + m * pf_dd] += dI_df_dd - dEs_df_dd;

    /* GW row */
    Jp_f[iGW + m * pa1] += dQperc_da1;
    Jp_f[iGW + m * pkperc] += dQperc_dkperc;
    Jp_f[iGW + m * pSmax] += dQperc_dSmax;
    Jp_f[iGW + m * pSfc] += dQperc_dSfc;

    Jp_f[iGW + m * pGmax] -= dQgw_dGmax;
    Jp_f[iGW + m * pCgw] -= dQgw_dC;
    Jp_f[iGW + m * pmm] -= dQgw_dmm;

    /* GIUH params via R_s */
    for (int i = 0; i < L; ++i) {
        const int row = iq1 + i;
        Jp_f[row + m * pSmax] = giuh[i] * dRs_dSmax;
        Jp_f[row + m * pKsch] = giuh[i] * dRs_dKsch;
        Jp_f[row + m * pT_tr] = giuh[i] * dRs_dT_tr;
        Jp_f[row + m * pf_dd] = giuh[i] * dRs_df_dd;
    }

    /* Nash: F1 depends on Q_lf params */
    Jp_f[iF1 + m * pSmax] += dQlf_dSmax;
    Jp_f[iF1 + m * pLfthr] += dQlf_dlfthr;
    Jp_f[iF1 + m * pa2] += dQlf_da2;
    Jp_f[iF1 + m * pkLf] += dQlf_dklf;

    /* Nash k_nsh */
    Jp_f[iF1 + m * pkNsh] += -F[0];
    for (int k = 1; k < K; ++k) {
        Jp_f[(iF1 + k) + m * pkNsh] += (F[k - 1] - F[k]);
    }

    /* Qcum row */
    Jp_f[iQcum + m * pGmax] += dQgw_dGmax;
    Jp_f[iQcum + m * pCgw] += dQgw_dC;
    Jp_f[iQcum + m * pmm] += dQgw_dmm;
    Jp_f[iQcum + m * pkNsh] += F[K - 1];

    /* ---- Chain rule to normalized theta (same as your old mapping) ---- */
    const double th2 = (s_max != 0.0) ? (s_fc / s_max) : 0.0;
    const double th3 = (s_max != 0.0) ? (s_wp / s_max) : 0.0;
    const double denom = (s_fc - s_wp);
    const double th7 = (fabs(denom) > 0.0) ? ((lf_thr - s_wp) / denom) : 0.0;

    const double dsmax_dth1 = 1.0;
    const double dsfc_dth1 = th2;
    const double dsfc_dth2 = s_max;
    const double dswp_dth1 = th3;
    const double dswp_dth3 = s_max;

    const double dlf_dth7 = (s_fc - s_wp);
    const double dlf_dth2 = th7 * s_max;
    const double dlf_dth3 = (1.0 - th7) * s_max;
    const double dlf_dth1 = (1.0 - th7) * th3 + th7 * th2;

    for (int i = 0; i < m; ++i) {
        /* th1 */
        Jth_f[i + m * 0] = Jp_f[i + m * pSmax] * dsmax_dth1 +
                           Jp_f[i + m * pSfc] * dsfc_dth1 + Jp_f[i + m * pSwp] * dswp_dth1 +
                           Jp_f[i + m * pLfthr] * dlf_dth1;

        /* th2 */
        Jth_f[i + m * 1] = Jp_f[i + m * pSfc] * dsfc_dth2 + Jp_f[i + m * pLfthr] * dlf_dth2;

        /* th3 */
        Jth_f[i + m * 2] = Jp_f[i + m * pSwp] * dswp_dth3 + Jp_f[i + m * pLfthr] * dlf_dth3;

        /* direct */
        Jth_f[i + m * 3] = Jp_f[i + m * pKsch];
        Jth_f[i + m * 4] = Jp_f[i + m * pa1];
        Jth_f[i + m * 5] = Jp_f[i + m * pkperc];
        Jth_f[i + m * 6] = Jp_f[i + m * pLfthr] * dlf_dth7;
        Jth_f[i + m * 7] = Jp_f[i + m * pa2];
        Jth_f[i + m * 8] = Jp_f[i + m * pkLf];
        Jth_f[i + m * 9] = Jp_f[i + m * pGmax];
        Jth_f[i + m * 10] = Jp_f[i + m * pCgw];
        Jth_f[i + m * 11] = Jp_f[i + m * pmm];
        Jth_f[i + m * 12] = Jp_f[i + m * pkNsh];
        Jth_f[i + m * 13] = Jp_f[i + m * pT_tr];
        Jth_f[i + m * 14] = Jp_f[i + m * pf_dd];
    }

    /* dSdt = Jx*S + Jth (m×d column-major) */
    for (int j = 0; j < d; ++j) {
        for (int i = 0; i < m; ++i) {
            dSdt[i + m * j] = Jth_f[i + m * j];
        }
    }

    for (int i = 0; i < m; ++i) {
        const double* Jrow = Jx_f + i * m;
        for (int k = 0; k < m; ++k) {
            const double coefficient = Jrow[k];
            if (coefficient == 0.0) {
                continue;
            }
            for (int j = 0; j < d; ++j) {
                dSdt[i + m * j] += coefficient * Smat[k + m * j];
            }
        }
    }

    std::free(dQdt);
    std::free(dF);
    std::free(Jx_f);
    std::free(Jp_f);
    std::free(Jth_f);
}

namespace sage_cfe_nwm {
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
    const int m = 3 + p.L + p.K + 1;
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
    const double s_max = p.s_max;
    const double s_fc = p.s_fc;
    const double s_wp = p.s_wp;
    const double k_sch = p.k_sch;
    const double a1 = p.a1;
    const double k_perc = p.k_perc;
    const double lf_thr = p.lf_thr;
    const double a2 = p.a2;
    const double k_lf = p.k_lf;
    const double g_max = p.g_max;
    const double c_gw = p.c_gw;
    const double mm = p.mm;
    const double k_nsh = p.k_nsh;
    const double T_tr = p.T_tr;
    const double f_dd = p.f_dd;
    const double dT = p.dT;
    const double T_sm = p.T_sm;
    const double eps_m = p.eps_m;
    const double eps = p.eps;
    const double rho = p.rho;
    const int K = p.K;
    const int L = p.L;
    const double* giuh = p.giuh;
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
                const_cast<double*>(giuh),
                T_tr,
                f_dd,
                K,
                L,
                dT,
                T_sm,
                eps_m,
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
} // namespace sage_cfe_nwm
