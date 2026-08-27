/*
 * gr4jA.cpp
 *
 * Conceptual GR4J-A rainfall-runoff model with an adaptive-step
 * explicit Runge-Kutta integrator. This file contains the MATLAB-independent
 * native numerical core shared by crr_gr4jA and crr_model_mex.
 *
 * Written by Jasper A. Vrugt.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#include "gr4jA.hpp"
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
                double x1,
                double x2,
                double x3,
                double x4,
                double x5,
                double f_p,
                double T_tr,
                double f_dd,
                int n1,
                int n2,
                double kappa,
                double bg,
                double bR,
                double mts,
                double T_sm,
                double eps_m,
                double eps_s,
                double rho);

static void gr4jA_aug_ode(double* z,
                          double* zdot,
                          double P,
                          double Ep,
                          double T,
                          double x1,
                          double x2,
                          double x3,
                          double x4,
                          double x5,
                          double f_p,
                          double T_tr,
                          double f_dd,
                          int n1,
                          int n2,
                          double kappa,
                          double bg,
                          double bR,
                          double mts,
                          double T_sm,
                          double eps_m,
                          double eps_s,
                          double rho,
                          int nvar,
                          int m,
                          int d);

static void gr4jA_odefcn(double* u,
                         double* udot,
                         const double* Smat,
                         double* dSdt,
                         double P,
                         double Ep,
                         double T,
                         double x1,
                         double x2,
                         double x3,
                         double x4,
                         double x5,
                         double f_p,
                         double T_tr,
                         double f_dd,
                         int n1,
                         int n2,
                         double kappa,
                         double bg,
                         double bR,
                         double mts,
                         double T_sm,
                         double eps_m,
                         double eps_s,
                         double rho,
                         int nvar,
                         int m,
                         int d);

static inline bool isFinite(double x)
{
    return std::isfinite(x);
}

static inline double smooth_pos(double a, double eps)
{
    return 0.5 * (a + sqrt(a * a + eps * eps));
}

static inline double d_smooth_pos_da(double a, double eps)
{
    return 0.5 * (1.0 + a / sqrt(a * a + eps * eps));
}

static inline double smooth_min(double A, double B, double eps)
{
    double d = A - B;
    return 0.5 * (A + B - sqrt(d * d + eps * eps));
}

static inline double d_smooth_min_dA(double A, double B, double eps)
{
    double d = A - B;
    double s = sqrt(d * d + eps * eps);
    return 0.5 * (1.0 - d / s);
}

static inline double d_smooth_min_dB(double A, double B, double eps)
{
    double d = A - B;
    double s = sqrt(d * d + eps * eps);
    return 0.5 * (1.0 + d / s);
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
                double x1,
                double x2,
                double x3,
                double x4,
                double x5,
                double f_p,
                double T_tr,
                double f_dd,
                int n1,
                int n2,
                double kappa,
                double bg,
                double bR,
                double mts,
                double T_sm,
                double eps_m,
                double eps_s,
                double rho)
{
    /* Euler */
    gr4jA_aug_ode(z,
                  zdotE,
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
                  d);
    for (int i = 0; i < nvar; ++i) {
        zE[i] = z[i] + h * zdotE[i];
    }

    /* Heun slope */
    gr4jA_aug_ode(zE,
                  zdot,
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

static void gr4jA_aug_ode(double* z,
                          double* zdot,
                          double P,
                          double Ep,
                          double T,
                          double x1,
                          double x2,
                          double x3,
                          double x4,
                          double x5,
                          double f_p,
                          double T_tr,
                          double f_dd,
                          int n1,
                          int n2,
                          double kappa,
                          double bg,
                          double bR,
                          double mts,
                          double T_sm,
                          double eps_m,
                          double eps_s,
                          double rho,
                          int nvar,
                          int m,
                          int d)
{
    double* u = z;        /* m */
    double* Smat = z + m; /* m*d */

    double* udot = zdot;     /* m */
    double* dSdt = zdot + m; /* m*d */

    gr4jA_odefcn(u,
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
                 d);
}

static void gr4jA_odefcn(double* u,
                         double* udot,
                         const double* Smat,
                         double* dSdt,
                         double P,
                         double Ep,
                         double T,
                         double x1,
                         double x2,
                         double x3,
                         double x4,
                         double x5,
                         double f_p,
                         double T_tr,
                         double f_dd,
                         int n1,
                         int n2,
                         double kappa,
                         double bg,
                         double bR,
                         double mts,
                         double T_sm,
                         double eps_m,
                         double eps_s,
                         double rho,
                         int nvar,
                         int m,
                         int d)
{
    (void)rho;
    (void)nvar;

    /* ----------------------------------------------------------------- */
    /* Indices (0-based)                                                 */
    /* ----------------------------------------------------------------- */
    const int iSwe = 0;
    const int iSp = 1;
    const int iRr = 2;
    const int iU1_1 = 3;
    const int iU1_e = 3 + n1 - 1;
    const int iU2_1 = 3 + n1;
    const int iU2_e = 3 + n1 + n2 - 1;
    const int iQinf = m - 1; /* should equal 3+n1+n2 */

    /* ----------------------------------------------------------------- */
    /* 4) Smooth positivity on storages (states only)                     */
    /* ----------------------------------------------------------------- */
    double Swe_u = u[iSwe];
    double Sp_u = u[iSp];
    double Rr_u = u[iRr];

    double Swe = smooth_pos(Swe_u, eps_s);
    double Sp = smooth_pos(Sp_u, eps_s);
    double Rr = smooth_pos(Rr_u, eps_s);

    double dSwe_du = d_smooth_pos_da(Swe_u, eps_s);
    double dSp_du = d_smooth_pos_da(Sp_u, eps_s);
    double dRr_du = d_smooth_pos_da(Rr_u, eps_s);

    /* smooth cascades */
    double *U1 = NULL, *U2 = NULL, *dU1_du = NULL, *dU2_du = NULL;
    if (n1 > 0) {
        U1 =
            static_cast<double*>(std::calloc(static_cast<std::size_t>(n1), sizeof(double)));
        dU1_du =
            static_cast<double*>(std::calloc(static_cast<std::size_t>(n1), sizeof(double)));
        for (int k = 0; k < n1; ++k) {
            double uk = u[iU1_1 + k];
            U1[k] = smooth_pos(uk, eps_s);
            dU1_du[k] = d_smooth_pos_da(uk, eps_s);
        }
    }
    if (n2 > 0) {
        U2 =
            static_cast<double*>(std::calloc(static_cast<std::size_t>(n2), sizeof(double)));
        dU2_du =
            static_cast<double*>(std::calloc(static_cast<std::size_t>(n2), sizeof(double)));
        for (int k = 0; k < n2; ++k) {
            double uk = u[iU2_1 + k];
            U2[k] = smooth_pos(uk, eps_s);
            dU2_du[k] = d_smooth_pos_da(uk, eps_s);
        }
    }

    /* ----------------------------------------------------------------- */
    /* 5) UH analogue time constants                                     */
    /* tau_1 = x4/max(1,n1), tau_2 = 2*x4/max(1,n2)                      */
    /* ----------------------------------------------------------------- */
    const double n1_eff = n1; // fmax(1.0, (double)n1);
    const double n2_eff = n2; // fmax(1.0, (double)n2);

    double tau_1 = x4 / n1_eff;
    double tau_2 = (2.0 * x4) / n2_eff;

    double inv_tau1 = 1.0 / tau_1;
    double inv_tau2 = 1.0 / tau_2;

    /* ----------------------------------------------------------------- */
    /* 6) Snow module (smooth degree-day)                                 */
    /* ----------------------------------------------------------------- */
    double T_smeps_m = fmax(T_sm, eps_m);

    double utemp = (T - T_tr) / T_smeps_m;
    double snow_fr = 0.5 * (1.0 - tanh(utemp));
    double rain_fr = 1.0 - snow_fr;

    double p_snow = P * snow_fr;
    double p_rain = P * rain_fr;

    double aT = (T - T_tr);
    double posT = 0.5 * (aT + sqrt(aT * aT + T_smeps_m * T_smeps_m));
    double p_pmelt = f_dd * posT;

    double p_amelt = smooth_min(Swe, p_pmelt, eps_m);
    double p_liq = p_rain + p_amelt;

    /* snow derivatives */
    double sech2 = 1.0 / (cosh(utemp) * cosh(utemp));
    double dsnowFrac_dT_tr = 0.5 * sech2 / T_smeps_m;

    double dposT_da = 0.5 * (1.0 + aT / sqrt(aT * aT + T_smeps_m * T_smeps_m));
    double dposT_dT_tr = -dposT_da;

    double dp_pmelt_dT_tr = f_dd * dposT_dT_tr;
    double dp_pmelt_df_dd = posT;

    double dM_dSwe = d_smooth_min_dA(Swe, p_pmelt, eps_m);
    double dM_dp_pmelt = d_smooth_min_dB(Swe, p_pmelt, eps_m);

    double dp_liq_dSwe = dM_dSwe;

    double dPsnow_dT_tr = P * dsnowFrac_dT_tr;
    double dPrain_dT_tr = -dPsnow_dT_tr;

    double dM_dT_tr = dM_dp_pmelt * dp_pmelt_dT_tr;
    double dM_df_dd = dM_dp_pmelt * dp_pmelt_df_dd;

    double dp_liq_dT_tr = dPrain_dT_tr + dM_dT_tr;
    double dp_liq_df_dd = dM_df_dd;

    /* ----------------------------------------------------------------- */
    /* 7) Production block                                               */
    /* ----------------------------------------------------------------- */
    double Ep_eff = f_p * Ep;
    double aPE = (p_liq - Ep_eff);
    double daPE_dfp = -Ep;

    double Pn = smooth_pos(aPE, eps_s);
    double En = smooth_pos(-aPE, eps_s);

    double dPn_dp_liq = d_smooth_pos_da(aPE, eps_s);   /* = dPn/daPE */
    double dEn_dp_liq = -d_smooth_pos_da(-aPE, eps_s); /* = dEn/daPE */

    /* wetness weight: wWet = 0.5*(1 + tanh(aPE/x1)) */
    double uPE = aPE / x1;
    double wWet = 0.5 * (1.0 + tanh(uPE));
    double sech2PE = 1.0 / (cosh(uPE) * cosh(uPE));

    double dwWet_daPE = 0.5 * sech2PE / x1;
    double dwWet_dp_liq = dwWet_daPE; /* daPE/dp_liq = 1 */
    double dwWet_dx1 = 0.5 * sech2PE * (-aPE / (x1 * x1));

    /* Percolation Pp: g = 1 + (kappa*(Sp/x1))^4, kappa=4/9 */
    double Sp_n = Sp / x1;

    double kappa4 = kappa * kappa * kappa * kappa;
    double Sp_n3 = Sp_n * Sp_n * Sp_n;
    // double g = 1.0 + pow(kappa*Sp_n, 4.0);
    double g = 1.0 + kappa4 * Sp_n3 * Sp_n;
    double g_m14 = pow(g, -0.25);
    double Pp = Sp * (1.0 - g_m14);

    double dSpn_dSp = 1.0 / x1;
    double dSpn_dx1 = -Sp / (x1 * x1);

    // double dg_dSpn = 4.0 * pow(kappa,4.0) * pow(Sp_n,3.0);
    double dg_dSpn = 4.0 * kappa4 * Sp_n3;
    double dg_m14_dg = -0.25 * pow(g, -1.25);

    double dPp_dSp = (1.0 - g_m14) + Sp * (-dg_m14_dg * dg_dSpn * dSpn_dSp);
    double dPp_dx1 = Sp * (-dg_m14_dg * dg_dSpn * dSpn_dx1);

    double dPerc_dS = dPp_dSp;
    double dPerc_dx1 = dPp_dx1;

    /* Wet branch Ps_wet (uses Pn) */
    double Ps_wet = 0.0, dPsWet_dS = 0.0, dPsWet_dPn = 0.0, dPsWet_dx1 = 0.0;

    {
        double zP = Pn / x1;
        double aP = tanh(zP);
        double daP_dPn = (1.0 - aP * aP) / x1;

        double Nw = (1.0 - Sp_n * Sp_n) * aP;
        double Dw = (1.0 + Sp_n * aP);

        Ps_wet = x1 * Nw / Dw;

        double dNw_dSpn = (-2.0 * Sp_n) * aP;
        double dNw_dPn = (1.0 - Sp_n * Sp_n) * daP_dPn;

        double dDw_dSpn = aP;
        double dDw_dPn = Sp_n * daP_dPn;

        double dPsWet_dSpn_loc = x1 * (Dw * dNw_dSpn - Nw * dDw_dSpn) / (Dw * Dw);
        dPsWet_dPn = x1 * (Dw * dNw_dPn - Nw * dDw_dPn) / (Dw * Dw);

        dPsWet_dS = dPsWet_dSpn_loc * dSpn_dSp;

        /* x1 derivative */
        double dzP_dx1 = -Pn / (x1 * x1);
        double daP_dx1 = (1.0 - aP * aP) * dzP_dx1;

        double dNw_dx1 = dNw_dSpn * dSpn_dx1 + (1.0 - Sp_n * Sp_n) * daP_dx1;
        double dDw_dx1 = dDw_dSpn * dSpn_dx1 + Sp_n * daP_dx1;

        double dNoverD_dx1 = (Dw * dNw_dx1 - Nw * dDw_dx1) / (Dw * Dw);
        dPsWet_dx1 = (Nw / Dw) + x1 * dNoverD_dx1;
    }

    /* Dry branch Es_dry (paper form: Nd = Sp*(2-Sp_n)*aE) */
    double Es_dry = 0.0, dEsDry_dS = 0.0, dEsDry_dEn = 0.0, dEsDry_dx1 = 0.0;

    {
        double zE = En / x1;
        double aE = tanh(zE);
        double daE_dEn = (1.0 - aE * aE) / x1;

        double Nd = Sp * (2.0 - Sp_n) * aE;
        double Dd = 1.0 + (1.0 - Sp_n) * aE;

        double invD = 1.0 / Dd;
        double invD2 = invD * invD;

        Es_dry = Nd * invD;

        double dNd_dSp = (2.0 - Sp_n) * aE;
        double dNd_dSpn = -Sp * aE;
        double dNd_daE = Sp * (2.0 - Sp_n);

        double dDd_dSpn = -aE;
        double dDd_daE = (1.0 - Sp_n);

        double dEs_dSp = dNd_dSp * invD;
        double dEs_dSpn = (dNd_dSpn * Dd - Nd * dDd_dSpn) * invD2;
        double dEs_daE = (dNd_daE * Dd - Nd * dDd_daE) * invD2;

        dEsDry_dS = dEs_dSp + dEs_dSpn * dSpn_dSp;
        dEsDry_dEn = dEs_daE * daE_dEn;

        double dzE_dx1 = -En / (x1 * x1);
        double daE_dx1 = (1.0 - aE * aE) * dzE_dx1;

        dEsDry_dx1 = dEs_dSpn * dSpn_dx1 + dEs_daE * daE_dx1;
    }

    /* Blend + derivatives */
    double Ps = wWet * Ps_wet;
    double Es = (1.0 - wWet) * Es_dry;

    double dPs_dS = wWet * dPsWet_dS;
    double dEs_dS = (1.0 - wWet) * dEsDry_dS;

    double dPsWet_dp_liq = dPsWet_dPn * dPn_dp_liq;
    double dEsDry_dp_liq = dEsDry_dEn * dEn_dp_liq;

    double dPs_dp_liq = wWet * dPsWet_dp_liq + Ps_wet * dwWet_dp_liq;
    double dEs_dp_liq = (1.0 - wWet) * dEsDry_dp_liq - Es_dry * dwWet_dp_liq;

    double dPs_dx1 = wWet * dPsWet_dx1 + Ps_wet * dwWet_dx1;
    double dEs_dx1 = (1.0 - wWet) * dEsDry_dx1 - Es_dry * dwWet_dx1;

    // double P1 = wWet * (Pn - Ps_wet);
    double P1 = Pn - Ps;
    // This single change restores the identity Ps+P1=Pn (mass conservation for the wet
    // side). If you want to keep everything as "raw branch + blending", the rule of thumb
    // is blend the fluxes, but don't multiply the "residual" by the same gate again.

    double dP1_dS = -dPs_dS;
    double dP1_daPE = dPn_dp_liq - dPs_dp_liq;
    double dP1_dx1 = -dPs_dx1;

    double dP1_dfp = (dPn_dp_liq - dPs_dp_liq) * daPE_dfp; // daPE_dfp = -Ep
    double dP1_dT_tr = dP1_daPE * dp_liq_dT_tr;
    double dP1_df_dd = dP1_daPE * dp_liq_df_dd;

    double Pr = P1 + Pp;

    double dPr_dS = dP1_dS + dPerc_dS;
    double dPr_dp_liq = dP1_daPE; // Pp depends on Sp only
    double dPr_dx1 = dP1_dx1 + dPerc_dx1;
    double dPr_dfp = dP1_dfp;
    double dPr_dT_tr = dP1_dT_tr;
    double dPr_df_dd = dP1_df_dd;

    /* ----------------------------------------------------------------- */
    /* 8) Routing + cascades                                             */
    /* ----------------------------------------------------------------- */
    double in1 = x5 * Pr;
    double in2 = (1.0 - x5) * Pr;

    double q_1 = (n1 > 0) ? (U1[n1 - 1] * inv_tau1) : 0.0;
    double q_2 = (n2 > 0) ? (U2[n2 - 1] * inv_tau2) : 0.0;

    double Rratio = Rr / x3;
    // double q_g = x2 * pow(Rratio, bg);
    const double Rratio2 = Rratio * Rratio;
    const double Rratio3 = Rratio2 * Rratio;
    const double sqrtRratio = sqrt(Rratio);
    double q_g = x2 * Rratio3 * sqrtRratio;

    // double q_r = x3 / ((bR - 1.0)*mts) * pow(Rratio, bR);
    double q_r = x3 / ((bR - 1.0) * mts) * Rratio3 * Rratio2;

    double q_d = q_2 + q_g;
    double q_dpos = smooth_pos(q_d, eps_s);
    double q_out = q_r + q_dpos;

    /* derivatives routing terms */
    // double dqg_dRr = x2 * bg * pow(Rratio, bg-1.0) * (1.0/x3);
    double dqg_dRr = x2 * bg * Rratio2 * sqrtRratio * (1.0 / x3);
    // double dqg_dx2 = pow(Rratio, bg);
    double dqg_dx2 = Rratio3 * sqrtRratio;
    // double dqg_dx3 = x2 * bg * pow(Rratio, bg-1.0) * (-Rr/(x3*x3));
    double dqg_dx3 = x2 * bg * Rratio2 * sqrtRratio * (-Rr / (x3 * x3));

    double c0 = 1.0 / ((bR - 1.0) * mts);
    // double dqr_dRr = c0 * bR * pow(Rr, bR-1.0) * pow(x3, 1.0-bR);
    const double Rr4 = Rr * Rr * Rr * Rr;
    const double x3_4 = x3 * x3 * x3 * x3;
    double dqr_dRr = c0 * bR * Rr4 * 1 / (x3_4);
    // double dqr_dx3 = c0 * pow(Rr, bR) * (1.0-bR) * pow(x3, -bR);
    double dqr_dx3 = c0 * Rr * Rr4 * (1.0 - bR) * 1 / (x3_4 * x3);

    /* ----------------------------------------------------------------- */
    /* 9) ODE RHS                                                        */
    /* ----------------------------------------------------------------- */
    udot[iSwe] = p_snow - p_amelt;
    udot[iSp] = Ps - Es - Pp;
    udot[iRr] = q_1 + q_g - q_r;

    if (n1 > 0) {
        udot[iU1_1] = in1 - U1[0] * inv_tau1;
        for (int k = 1; k < n1; ++k) {
            int idx = iU1_1 + k;
            udot[idx] = U1[k - 1] * inv_tau1 - U1[k] * inv_tau1;
        }
    }
    if (n2 > 0) {
        udot[iU2_1] = in2 - U2[0] * inv_tau2;
        for (int k = 1; k < n2; ++k) {
            int idx = iU2_1 + k;
            udot[idx] = U2[k - 1] * inv_tau2 - U2[k] * inv_tau2;
        }
    }
    udot[iQinf] = q_out;

    /* ----------------------------------------------------------------- */
    /* 10) Jx_f = df/dx (then chain rule for smooth_pos on states)       */
    /* ----------------------------------------------------------------- */
    std::vector<double> Jx_f_((size_t)m * (size_t)m, 0.0);
    std::vector<double> Jth_f_((size_t)m * (size_t)d, 0.0);

#define Jx_f(i, k) (Jx_f_[(size_t)(i) * (size_t)m + (size_t)(k)])   /* row-major */
#define Jth_f(i, j) (Jth_f_[(size_t)(i) + (size_t)m * (size_t)(j)]) /* col-major */

    /* Swe */
    Jx_f(iSwe, iSwe) = -dM_dSwe;

    /* Sp */
    Jx_f(iSp, iSwe) = (dPs_dp_liq - dEs_dp_liq) * dp_liq_dSwe;
    Jx_f(iSp, iSp) = dPs_dS - dEs_dS - dPerc_dS;

    /* cascade injections */
    double din1_dSwe = x5 * dPr_dp_liq * dp_liq_dSwe;
    double din1_dS = x5 * dPr_dS;

    double din2_dSwe = (1.0 - x5) * dPr_dp_liq * dp_liq_dSwe;
    double din2_dS = (1.0 - x5) * dPr_dS;

    /* U1 Jacobian */
    if (n1 > 0) {
        Jx_f(iU1_1, iSwe) = din1_dSwe;
        Jx_f(iU1_1, iSp) = din1_dS;
        Jx_f(iU1_1, iU1_1) = -inv_tau1;
        for (int k = 1; k < n1; ++k) {
            int idx = iU1_1 + k;
            Jx_f(idx, idx - 1) = inv_tau1;
            Jx_f(idx, idx) = -inv_tau1;
        }
        Jx_f(iRr, iU1_e) = inv_tau1; /* q_1 = U1(end)*inv_tau1 */
    }

    /* U2 Jacobian */
    if (n2 > 0) {
        Jx_f(iU2_1, iSwe) = din2_dSwe;
        Jx_f(iU2_1, iSp) = din2_dS;
        Jx_f(iU2_1, iU2_1) = -inv_tau2;
        for (int k = 1; k < n2; ++k) {
            int idx = iU2_1 + k;
            Jx_f(idx, idx - 1) = inv_tau2;
            Jx_f(idx, idx) = -inv_tau2;
        }
    }

    /* routing store + Qinf */
    Jx_f(iRr, iRr) = dqg_dRr - dqr_dRr;

    double dQdpos_dQd = d_smooth_pos_da(q_d, eps_s);
    Jx_f(iQinf, iRr) = dqr_dRr + dQdpos_dQd * dqg_dRr;
    if (n2 > 0) {
        Jx_f(iQinf, iU2_e) = dQdpos_dQd * inv_tau2;
    }

    /* chain rule scaling for smoothed states */
    for (int ii = 0; ii < m; ++ii) {
        Jx_f(ii, iSwe) *= dSwe_du;
        Jx_f(ii, iSp) *= dSp_du;
        Jx_f(ii, iRr) *= dRr_du;
    }
    if (n1 > 0) {
        for (int k = 0; k < n1; ++k) {
            int col = iU1_1 + k;
            for (int ii = 0; ii < m; ++ii) {
                Jx_f(ii, col) *= dU1_du[k];
            }
        }
    }
    if (n2 > 0) {
        for (int k = 0; k < n2; ++k) {
            int col = iU2_1 + k;
            for (int ii = 0; ii < m; ++ii) {
                Jx_f(ii, col) *= dU2_du[k];
            }
        }
    }

    /* ----------------------------------------------------------------- */
    /* 11) Jth_f = df/dtheta (mxd), theta=[x1 x2 x3 x4 x5 f_p T_tr f_dd] */
    /* column mapping: 0:x1 1:x2 2:x3 3:x4 4:x5 5:f_p 6:T_tr 7:f_dd      */
    /* ----------------------------------------------------------------- */

    /* --- snow chain terms --- */
    double dPs_dT_tr = dPs_dp_liq * dp_liq_dT_tr;
    double dPs_df_dd = dPs_dp_liq * dp_liq_df_dd;

    double dEs_dT_tr = dEs_dp_liq * dp_liq_dT_tr;
    double dEs_df_dd = dEs_dp_liq * dp_liq_df_dd;

    /* T_tr, f_dd */
    Jth_f(iSwe, 6) = dPsnow_dT_tr - dM_dT_tr;
    Jth_f(iSwe, 7) = -dM_df_dd;

    Jth_f(iSp, 6) = dPs_dT_tr - dEs_dT_tr;
    Jth_f(iSp, 7) = dPs_df_dd - dEs_df_dd;

    /* x1 */
    Jth_f(iSp, 0) = dPs_dx1 - dEs_dx1 - dPerc_dx1;

    /* f_p affects aPE only: d/d(f_p) = d/d(aPE) * (-Ep) */
    {
        double dPs_dfp = dPs_dp_liq * daPE_dfp; /* daPE_dfp = -Ep */
        double dEs_dfp = dEs_dp_liq * daPE_dfp;
        Jth_f(iSp, 5) = dPs_dfp - dEs_dfp;

        if (n1 > 0) {
            Jth_f(iU1_1, 5) = x5 * dPr_dfp;
        }
        if (n2 > 0) {
            Jth_f(iU2_1, 5) = (1.0 - x5) * dPr_dfp;
        }
    }

    /* x2, x3 */
    Jth_f(iRr, 1) = dqg_dx2;
    Jth_f(iQinf, 1) = dQdpos_dQd * dqg_dx2;

    Jth_f(iRr, 2) = dqg_dx3 - dqr_dx3;
    Jth_f(iQinf, 2) = dqr_dx3 + dQdpos_dQd * dqg_dx3;

    /* x4 affects inv_tau1, inv_tau2 */
    {
        double dinv_tau1_dx4 = -(1.0 / (tau_1 * tau_1)) * (1.0 / n1_eff);
        double dinv_tau2_dx4 = -(1.0 / (tau_2 * tau_2)) * (2.0 / n2_eff);

        if (n1 > 0) {
            Jth_f(iU1_1, 3) = -(U1[0]) * dinv_tau1_dx4;
            for (int k = 1; k < n1; ++k) {
                int idx = iU1_1 + k;
                Jth_f(idx, 3) = (U1[k - 1] - U1[k]) * dinv_tau1_dx4;
            }
            Jth_f(iRr, 3) = U1[n1 - 1] * dinv_tau1_dx4;
        }
        if (n2 > 0) {
            Jth_f(iU2_1, 3) = -(U2[0]) * dinv_tau2_dx4;
            for (int k = 1; k < n2; ++k) {
                int idx = iU2_1 + k;
                Jth_f(idx, 3) = (U2[k - 1] - U2[k]) * dinv_tau2_dx4;
            }
            Jth_f(iQinf, 3) = dQdpos_dQd * (U2[n2 - 1] * dinv_tau2_dx4);
        }
    }

    /* x5 split between in1 and in2 */
    if (n1 > 0) {
        Jth_f(iU1_1, 4) = Pr;
    }
    if (n2 > 0) {
        Jth_f(iU2_1, 4) = -Pr;
    }

    /* injection sensitivities to x1, T_tr, f_dd through Pr */
    if (n1 > 0) {
        Jth_f(iU1_1, 0) = x5 * dPr_dx1;
        Jth_f(iU1_1, 6) = x5 * dPr_dT_tr;
        Jth_f(iU1_1, 7) = x5 * dPr_df_dd;
    }
    if (n2 > 0) {
        Jth_f(iU2_1, 0) = (1.0 - x5) * dPr_dx1;
        Jth_f(iU2_1, 6) = (1.0 - x5) * dPr_dT_tr;
        Jth_f(iU2_1, 7) = (1.0 - x5) * dPr_df_dd;
    }

    /* ----------------------------------------------------------------- */
    /* 12) Sensitivity ODE: dS/dt = Jx*S + Jth                           */
    /* ----------------------------------------------------------------- */
    for (int j = 0; j < d; ++j) {
        for (int i = 0; i < m; ++i) {
            dSdt[(size_t)i + (size_t)m * (size_t)j] = Jth_f(i, j);
        }
    }

    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < m; ++k) {
            const double coefficient = Jx_f(i, k);
            if (coefficient == 0.0) {
                continue;
            }
            for (int j = 0; j < d; ++j) {
                const size_t ij = (size_t)i + (size_t)m * (size_t)j;
                const size_t kj = (size_t)k + (size_t)m * (size_t)j;
                dSdt[ij] += coefficient * Smat[kj];
            }
        }
    }

    /* free */
    if (U1) {
        std::free(U1);
    }
    if (U2) {
        std::free(U2);
    }
    if (dU1_du) {
        std::free(dU1_du);
    }
    if (dU2_du) {
        std::free(dU2_du);
    }

#undef Jx_f
#undef Jth_f
}

namespace sage_gr4ja {
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
    const int d = 8;
    const int m = p.n1 + p.n2 + 4;
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
    const double x1 = p.x1;
    const double x2 = p.x2;
    const double x3 = p.x3;
    const double x4 = p.x4;
    const double x5 = p.x5;
    const double f_p = p.f_p;
    const double T_tr = p.T_tr;
    const double f_dd = p.f_dd;
    const double kappa = p.kappa;
    const double bg = p.bg;
    const double bR = p.bR;
    const double mts = p.mts;
    const double T_sm = p.T_sm;
    const double eps_m = p.eps_m;
    const double eps_s = p.eps_s;
    const double rho = p.rho;
    const int n1 = p.n1;
    const int n2 = p.n2;

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
} // namespace sage_gr4ja
