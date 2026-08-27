#include "gr4jB.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

static void gr4jB_aug_ode(double* z,
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
                          const double* U,
                          const double* dUdx4,
                          int L,
                          double eta,
                          double tau,
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

static void gr4jB_odefcn(double* u,
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
                         const double* U,
                         const double* dUdx4,
                         int L,
                         double eta,
                         double tau,
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

/* ---------------- Smooth helpers ---------------- */
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
                const double* U,
                const double* dUdx4,
                int L,
                double eta,
                double tau,
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
    gr4jB_aug_ode(z,
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
                  U,
                  dUdx4,
                  L,
                  eta,
                  tau,
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
    gr4jB_aug_ode(zE,
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
                  U,
                  dUdx4,
                  L,
                  eta,
                  tau,
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

/* ------------------------------------------------------------------ */
/* Augmented ODE: z = [u; Smat]                                       */
/* ------------------------------------------------------------------ */
static void gr4jB_aug_ode(double* z,
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
                          const double* U,
                          const double* dUdx4,
                          int L,
                          double eta,
                          double tau,
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

    gr4jB_odefcn(u,
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
                 U,
                 dUdx4,
                 L,
                 eta,
                 tau,
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

/* ------------------------------------------------------------------ */
/* gr4jB ODE + sensitivities  (MATLAB-consistent version)              */
/* u = [Swe, Sp, Rr, U1(1..n1), U2(1..n2), Qinf]                  */
/* theta = [x1 x2 x3 x4 x5 f_p T_tr f_dd]                             */
/* ------------------------------------------------------------------ */
static void gr4jB_odefcn(double* u,
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
                         const double* U,
                         const double* dUdx4,
                         int L,
                         double eta,
                         double tau,
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
    (void)eta;
    (void)tau; /* keep if unused in this snippet */

    /* parameter column indices: theta=[x1 x2 x3 x4 x5 f_p T_tr f_dd] */
    enum {
        c_x1 = 0,
        c_x2 = 1,
        c_x3 = 2,
        c_x4 = 3,
        c_x5 = 4,
        c_fp = 5,
        c_Ttr = 6,
        c_fdd = 7
    };

    /* ----------------------------------------------------------------- */
    /* Indices (0-based) and layout check                                */
    /* ----------------------------------------------------------------- */
    const int iSwe = 0;
    const int iSp = 1;
    const int iRr = 2;
    const int iz1 = 3;
    const int iz2 = iz1 + L;
    const int iQ = iz2 + L;

    /* m must match [3 + L + L + 1] */
    if (m != (4 + 2 * L)) {
        throw std::invalid_argument(
            "GR4J-B state layout must satisfy m = 4 + 2*L.");
    }
    if (d != 8) {
        throw std::invalid_argument(
            "GR4J-B requires eight hydrologic parameters.");
    }

    /* pointers to UH states */
    double* z1 = &u[iz1];
    double* z2 = &u[iz2];

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

    /* ----------------------------------------------------------------- */
    /* 6) Snow module (smooth degree-day)                                */
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

    double dPn_dp_liq = d_smooth_pos_da(aPE, eps_s);   /* dPn/daPE */
    double dEn_dp_liq = -d_smooth_pos_da(-aPE, eps_s); /* dEn/daPE */

    /* wetness weight: wWet = 0.5*(1 + tanh(aPE/x1)) */
    double uPE = aPE / x1;
    double wWet = 0.5 * (1.0 + tanh(uPE));
    double sech2PE = 1.0 / (cosh(uPE) * cosh(uPE));

    double dwWet_daPE = 0.5 * sech2PE / x1;
    double dwWet_dp_liq = dwWet_daPE; /* daPE/dp_liq = 1 */
    double dwWet_dx1 = 0.5 * sech2PE * (-aPE / (x1 * x1));

    /* Percolation Pp: g = 1 + (kappa*(Sp/x1))^4 */
    double Sp_n = Sp / x1;

    double kappa4 = kappa * kappa * kappa * kappa;
    double Sp_n3 = Sp_n * Sp_n * Sp_n;
    double g = 1.0 + kappa4 * Sp_n3 * Sp_n; /* = 1 + (kappa*Sp_n)^4 */
    double g_m14 = pow(g, -0.25);
    double Pp = Sp * (1.0 - g_m14);

    double dSpn_dSp = 1.0 / x1;
    double dSpn_dx1 = -Sp / (x1 * x1);

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

    /* Dry branch Es_dry */
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

    /* KEY FIX (must match MATLAB): */
    double P1 = Pn - Ps;

    double dP1_dS = -dPs_dS;
    double dP1_daPE = dPn_dp_liq - dPs_dp_liq;
    double dP1_dx1 = -dPs_dx1;

    double dP1_dfp = (dPn_dp_liq - dPs_dp_liq) * daPE_dfp;
    double dP1_dT_tr = dP1_daPE * dp_liq_dT_tr;
    double dP1_df_dd = dP1_daPE * dp_liq_df_dd;

    double Pr = P1 + Pp;

    double dPr_dS = dP1_dS + dPerc_dS;
    double dPr_dp_liq = dP1_daPE; /* Pp depends on Sp only */
    double dPr_dx1 = dP1_dx1 + dPerc_dx1;
    double dPr_dfp = dP1_dfp;
    double dPr_dT_tr = dP1_dT_tr;
    double dPr_df_dd = dP1_df_dd;

    /* ----------------------------------------------------------------- */
    /* 2) UH convolution routing states                                  */
    /* ----------------------------------------------------------------- */
    const double q_in1 = x5 * Pr;
    const double q_in2 = (1.0 - x5) * Pr;

    udot[iz1 + 0] = q_in1 - z1[0];
    udot[iz2 + 0] = q_in2 - z2[0];
    for (int k = 1; k < L; ++k) {
        udot[iz1 + k] = z1[k - 1] - z1[k];
        udot[iz2 + k] = z2[k - 1] - z2[k];
    }

    /* q1 = z1'*U, q2 = z2'*U */
    double q_1 = 0.0, q_2 = 0.0;
    double dq1_dx4 = 0.0, dq2_dx4 = 0.0;
    for (int k = 0; k < L; ++k) {
        q_1 += z1[k] * U[k];
        q_2 += z2[k] * U[k];
        dq1_dx4 += z1[k] * dUdx4[k];
        dq2_dx4 += z2[k] * dUdx4[k];
    }

    /* ----------------------------------------------------------------- */
    /* 8) Routing (consistent pow-form everywhere)                        */
    /* ----------------------------------------------------------------- */
    double Rratio = Rr / x3;

    /* groundwater exchange (or q_g): x2*(Rr/x3)^bg */
    double q_g = x2 * pow(Rratio, bg);

    /* routing outflow: c0*x3*(Rr/x3)^bR  with c0=1/((bR-1)*mts) */
    double c0 = 1.0 / ((bR - 1.0) * mts);
    double q_r = c0 * x3 * pow(Rratio, bR);

    /* derivatives */
    double dqg_dRr = x2 * bg * pow(Rratio, bg - 1.0) * (1.0 / x3);
    double dqg_dx2 = pow(Rratio, bg);
    double dqg_dx3 = x2 * bg * pow(Rratio, bg - 1.0) * (-Rr / (x3 * x3));

    /* q_r = c0*x3*Rratio^bR */
    double dqr_dRr =
        c0 * bR * pow(Rratio, bR - 1.0); /* since dRratio/dRr = 1/x3, cancels x3 */
    double dqr_dx3 =
        c0 * (pow(Rratio, bR) + x3 * bR * pow(Rratio, bR - 1.0) * (-Rr / (x3 * x3)));

    double q_d = q_2 + q_g;
    double q_dpos = smooth_pos(q_d, eps_s);
    double q_out = q_r + q_dpos;

    /* ----------------------------------------------------------------- */
    /* 9) ODE RHS                                                        */
    /* ----------------------------------------------------------------- */
    udot[iSwe] = p_snow - p_amelt;
    udot[iSp] = Ps - Es - Pp;
    udot[iRr] = q_1 + q_g - q_r;
    udot[iQ] = q_out;

    /* ----------------------------------------------------------------- */
    /* 10) Build the parameter-source matrix Jth_f                        */
    /* ----------------------------------------------------------------- */
    std::vector<double> Jth_f_((size_t)m * (size_t)d, 0.0);

#define Jth_f(i, j) (Jth_f_[(size_t)(i) + (size_t)m * (size_t)(j)]) /* col-major */

    const double dfSwe_dSwe = -dM_dSwe * dSwe_du;
    const double dfSp_dSwe = (dPs_dp_liq - dEs_dp_liq) * dp_liq_dSwe * dSwe_du;
    const double dfSp_dSp = (dPs_dS - dEs_dS - dPerc_dS) * dSp_du;

    const double dqin1_dSwe = x5 * dPr_dp_liq * dp_liq_dSwe;
    const double dqin1_dSp = x5 * dPr_dS;
    const double dqin2_dSwe = (1.0 - x5) * dPr_dp_liq * dp_liq_dSwe;
    const double dqin2_dSp = (1.0 - x5) * dPr_dS;

    /* Qinf: qout = qr + smooth_pos(q2+qg) */
    const double dpos_dqd = d_smooth_pos_da(q_d, eps_s);
    const double dz1_dSwe = dqin1_dSwe * dSwe_du;
    const double dz1_dSp = dqin1_dSp * dSp_du;
    const double dz2_dSwe = dqin2_dSwe * dSwe_du;
    const double dz2_dSp = dqin2_dSp * dSp_du;
    const double dfRr_dRr = (dqg_dRr - dqr_dRr) * dRr_du;
    const double dfQ_dRr = (dqr_dRr + dpos_dqd * dqg_dRr) * dRr_du;

    /* ---- Jth_f: snow/prod terms ---- */
    Jth_f(iSwe, c_Ttr) = dPsnow_dT_tr - dM_dT_tr;
    Jth_f(iSwe, c_fdd) = -dM_df_dd;

    const double dPs_dTtr = dPs_dp_liq * dp_liq_dT_tr;
    const double dEs_dTtr = dEs_dp_liq * dp_liq_dT_tr;
    const double dPs_dfdd = dPs_dp_liq * dp_liq_df_dd;
    const double dEs_dfdd = dEs_dp_liq * dp_liq_df_dd;

    Jth_f(iSp, c_Ttr) = dPs_dTtr - dEs_dTtr;
    Jth_f(iSp, c_fdd) = dPs_dfdd - dEs_dfdd;

    /* x1 in production store equation (Ps,Es,Pp all depend on x1) */
    Jth_f(iSp, c_x1) = dPs_dx1 - dEs_dx1 - dPerc_dx1;

    /* f_p enters via aPE = p_liq - f_p*Ep: d/df_p = d/daPE * (-Ep) */
    /* If MATLAB includes fp in Sp equation explicitly, keep this: */
    Jth_f(iSp, c_fp) = (dPs_dp_liq - dEs_dp_liq) * daPE_dfp; /* daPE_dfp=-Ep */

    /* x2,x3 in routing store and Qinf */
    Jth_f(iRr, c_x2) = dqg_dx2;
    Jth_f(iQ, c_x2) = dpos_dqd * dqg_dx2;

    Jth_f(iRr, c_x3) = dqg_dx3 - dqr_dx3;
    Jth_f(iQ, c_x3) = dqr_dx3 + dpos_dqd * dqg_dx3;

    /* x5 split affects injections */
    Jth_f(iz1 + 0, c_x5) = Pr;
    Jth_f(iz2 + 0, c_x5) = -Pr;

    /* injection sensitivity via Pr */
    Jth_f(iz1 + 0, c_x1) += x5 * dPr_dx1;
    Jth_f(iz2 + 0, c_x1) += (1.0 - x5) * dPr_dx1;

    Jth_f(iz1 + 0, c_fp) += x5 * dPr_dfp;
    Jth_f(iz2 + 0, c_fp) += (1.0 - x5) * dPr_dfp;

    Jth_f(iz1 + 0, c_Ttr) += x5 * dPr_dT_tr;
    Jth_f(iz2 + 0, c_Ttr) += (1.0 - x5) * dPr_dT_tr;

    Jth_f(iz1 + 0, c_fdd) += x5 * dPr_df_dd;
    Jth_f(iz2 + 0, c_fdd) += (1.0 - x5) * dPr_df_dd;

    /* x4 enters through U(x4): q1=z1'*U, q2=z2'*U */
    Jth_f(iRr, c_x4) += dq1_dx4;
    Jth_f(iQ, c_x4) += dpos_dqd * dq2_dx4;

    /* ----------------------------------------------------------------- */
    /* 11) Sparse sensitivity ODE: dS/dt = Jx_f*S + Jth_f                */
    /* Smat is m x d column-major. Only nonzero state dependencies are   */
    /* evaluated; the dense m-by-m state Jacobian is not materialized.   */
    /* ----------------------------------------------------------------- */
    for (int j = 0; j < d; ++j) {
        const double* Sj = Smat + (size_t)m * (size_t)j;
        double* dSj = dSdt + (size_t)m * (size_t)j;

        dSj[iSwe] = dfSwe_dSwe * Sj[iSwe] + Jth_f(iSwe, j);
        dSj[iSp] = dfSp_dSwe * Sj[iSwe] + dfSp_dSp * Sj[iSp] + Jth_f(iSp, j);

        dSj[iz1] = dz1_dSwe * Sj[iSwe] + dz1_dSp * Sj[iSp] - Sj[iz1] + Jth_f(iz1, j);
        dSj[iz2] = dz2_dSwe * Sj[iSwe] + dz2_dSp * Sj[iSp] - Sj[iz2] + Jth_f(iz2, j);

        for (int k = 1; k < L; ++k) {
            dSj[iz1 + k] = Sj[iz1 + k - 1] - Sj[iz1 + k] + Jth_f(iz1 + k, j);
            dSj[iz2 + k] = Sj[iz2 + k - 1] - Sj[iz2 + k] + Jth_f(iz2 + k, j);
        }

        double routed1 = 0.0;
        double routed2 = 0.0;
        for (int k = 0; k < L; ++k) {
            routed1 += U[k] * Sj[iz1 + k];
            routed2 += U[k] * Sj[iz2 + k];
        }

        dSj[iRr] = routed1 + dfRr_dRr * Sj[iRr] + Jth_f(iRr, j);
        dSj[iQ] = dpos_dqd * routed2 + dfQ_dRr * Sj[iRr] + Jth_f(iQ, j);
    }

#undef Jth_f
}

} // namespace

namespace sage_gr4jb {
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
    const int m = 4 + 2 * p.L;
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
    const double eta = p.eta;
    const double tau = p.tau;
    const double kappa = p.kappa;
    const double bg = p.bg;
    const double bR = p.bR;
    const double mts = p.mts;
    const double T_sm = p.T_sm;
    const double eps_m = p.eps_m;
    const double eps_s = p.eps_s;
    const double rho = p.rho;
    const int L = p.L;
    const double* U = p.U;
    const double* dUdx4 = p.dUdx4;
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
                U,
                dUdx4,
                L,
                eta,
                tau,
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
} // namespace sage_gr4jb
