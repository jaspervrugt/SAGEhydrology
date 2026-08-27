/*
 * xinanjiang.cpp
 *
 * Conceptual Xinanjiang rainfall-runoff model with an adaptive-step
 * explicit Runge-Kutta integrator. This file contains the MATLAB-independent
 * native numerical core shared by crr_xinanjiang and crr_model_mex.
 *
 * Written by Jasper A. Vrugt.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#include "xinanjiang.hpp"
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
                double f_p,
                double A_im,
                double a,
                double b,
                double W_max,
                double LM,
                double c,
                double S_max,
                double S_tot,
                double f_wm,
                double f_lm,
                double Ex,
                double k_i,
                double k_g,
                double c_i,
                double c_g,
                double k_f,
                double T_tr,
                double f_dd,
                double T_sm,
                double eps_m,
                double eps_r,
                double eps,
                double rho);

static void xinanjiang_aug_ode(double* z,
                               double* zdot,
                               double P,
                               double Ep,
                               double T,
                               double f_p,
                               double A_im,
                               double a,
                               double b,
                               double W_max,
                               double LM,
                               double c,
                               double S_max,
                               double S_tot,
                               double f_wm,
                               double f_lm,
                               double Ex,
                               double k_i,
                               double k_g,
                               double c_i,
                               double c_g,
                               double k_f,
                               double T_tr,
                               double f_dd,
                               double T_sm,
                               double eps_m,
                               double eps_r,
                               double eps,
                               double rho,
                               int nvar,
                               int m,
                               int d);

static void xinanjiang_odefcn(const double* u,
                              double* udot,
                              const double* S,
                              double* dSdt,
                              double P,
                              double Ep,
                              double T,
                              double f_p,
                              double A_im,
                              double a,
                              double b,
                              double W_max,
                              double LM,
                              double c,
                              double S_max,
                              double S_tot,
                              double f_wm,
                              double f_lm,
                              double Ex,
                              double k_i,
                              double k_g,
                              double c_i,
                              double c_g,
                              double k_f,
                              double T_tr,
                              double f_dd,
                              double T_sm,
                              double eps_m,
                              double eps_r,
                              double eps,
                              double rho,
                              int m,
                              int d);

inline double lf_func(double S, double Smax, double rho, double eps)
{
    double num = S - (Smax - rho * Smax * eps);
    double den = rho * Smax;
    double z = -(num / den);
    double ez = exp(z);
    return 1.0 / (1.0 + ez);
}

inline double dlf_dS(double lf_val, double Smax, double rho)
{
    /* d lf / dS = lf * (1-lf) / (rho*Smax) */
    return (lf_val * (1.0 - lf_val)) / (rho * Smax);
}

static inline double smooth_pos(double a, double eps)
{
    return 0.5 * (a + sqrt(a * a + eps * eps));
}

static inline double d_smooth_pos_da(double a, double eps)
{
    const double denom = sqrt(a * a + eps * eps);
    return 0.5 * (1.0 + a / denom);
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
                double f_p,
                double A_im,
                double a,
                double b,
                double W_max,
                double LM,
                double c,
                double S_max,
                double S_tot,
                double f_wm,
                double f_lm,
                double Ex,
                double k_i,
                double k_g,
                double c_i,
                double c_g,
                double k_f,
                double T_tr,
                double f_dd,
                double T_sm,
                double eps_m,
                double eps_r,
                double eps,
                double rho)
{
    /* Euler */
    xinanjiang_aug_ode(z,
                       zdotE,
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
                       nvar,
                       m,
                       d);

    for (int i = 0; i < nvar; ++i) {
        zE[i] = z[i] + h * zdotE[i];
    }

    /* Heun */
    xinanjiang_aug_ode(zE,
                       zdot,
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
                       nvar,
                       m,
                       d);

    /* Update */
    for (int i = 0; i < nvar; ++i) {
        z[i] = z[i] + 0.5 * h * (zdotE[i] + zdot[i]);
    }

    /* LTE */
    for (int i = 0; i < nvar; ++i) {
        LTE[i] = fabs(zE[i] - z[i]);
    }
}

static void xinanjiang_aug_ode(double* z,
                               double* zdot,
                               double P,
                               double Ep,
                               double T,
                               double f_p,
                               double A_im,
                               double a,
                               double b,
                               double W_max,
                               double LM,
                               double c,
                               double S_max,
                               double S_tot,
                               double f_wm,
                               double f_lm,
                               double Ex,
                               double k_i,
                               double k_g,
                               double c_i,
                               double c_g,
                               double k_f,
                               double T_tr,
                               double f_dd,
                               double T_sm,
                               double eps_m,
                               double eps_r,
                               double eps,
                               double rho,
                               int nvar,
                               int m,
                               int d)
{
    (void)nvar; /* not used, but kept for symmetry */

    double* u = z;        /* m */
    double* Smat = z + m; /* m*d */

    double* udot = zdot;     /* m */
    double* dSdt = zdot + m; /* m*d */

    xinanjiang_odefcn(u,
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
                      d);
}

static void xinanjiang_odefcn(const double* u,
                              double* udot,
                              const double* Smat,
                              double* dSdt,
                              double P,
                              double Ep,
                              double T,
                              double f_p,
                              double A_im,
                              double a,
                              double b,
                              double W_max,
                              double LM,
                              double c,
                              double S_max,
                              double S_tot,
                              double f_wm,
                              double f_lm,
                              double Ex,
                              double k_i,
                              double k_g,
                              double c_i,
                              double c_g,
                              double k_f,
                              double T_tr,
                              double f_dd,
                              double T_sm,
                              double eps_m,
                              double eps_r,
                              double eps,
                              double rho,
                              int m,
                              int d)
{
    (void)eps;
    (void)rho; /* not used in the MATLAB snippet */

    /* -----------------------------
     * 1. Unpack states
     * ----------------------------- */
    const double Swe_u = u[0]; /* raw snow water equivalent */
    const double W_u = u[1];   /* raw tension water storage */
    const double Sw_u = u[2];  /* raw free water storage    */

    /* Routing states are NOT smoothed in your MATLAB code */
    const double Si = u[3];
    const double Sg = u[4];
    const double Sf1 = u[5];
    const double Sf2 = u[6];
    const double Sf3 = u[7];

    /* MATLAB: smooth_pos on Swe, W, Sw */
    const double Swe = smooth_pos(Swe_u, eps_m);
    const double W = smooth_pos(W_u, eps_m);
    const double Sw = smooth_pos(Sw_u, eps_m);

    const double dSwe_dx1 = d_smooth_pos_da(Swe_u, eps_m);
    const double dW_dx2 = d_smooth_pos_da(W_u, eps_m);
    const double dSw_dx3 = d_smooth_pos_da(Sw_u, eps_m);

    /* ------------------------------------------------------------------ */
    /* 2) Snow module (smooth HBV-style degree-day)                       */
    /* ------------------------------------------------------------------ */
    const double T_smeps_m = fmax(T_sm, eps_m);
    const double utemp = (T - T_tr) / T_smeps_m;

    const double snow_fr = 0.5 * (1.0 - tanh(utemp));
    const double rain_fr = 1.0 - snow_fr;

    const double P_snow = P * snow_fr;
    const double P_rain = P * rain_fr;

    const double aT = (T - T_tr);
    const double posT =
        0.5 * (aT + sqrt(aT * aT + T_smeps_m * T_smeps_m)); /* smooth max(T-T_tr,0) */
    const double M_pot = f_dd * posT;

    /* smooth min(Swe, M_pot) */
    const double dxy = (Swe - M_pot);
    const double sqrtm = sqrt(dxy * dxy + eps_m * eps_m);
    const double M = 0.5 * (Swe + M_pot - sqrtm);

    const double Pliq = P_rain + M;

    /* Snow derivatives (as MATLAB) */
    const double coshu = cosh(utemp);
    const double sech2 = 1.0 / (coshu * coshu);
    const double dsnow_dT_tr = 0.5 * sech2 / T_smeps_m;

    const double denom_pos = sqrt(aT * aT + T_smeps_m * T_smeps_m);
    const double dposT_da = 0.5 * (1.0 + aT / denom_pos);
    const double dposT_dT_tr = -dposT_da; /* da/dT_tr = -1 */

    const double dMpot_dT_tr = f_dd * dposT_dT_tr;
    const double dMpot_f_dd = posT;

    const double dM_dSwe = 0.5 * (1.0 - dxy / sqrtm);
    const double dM_dMpot = 0.5 * (1.0 + dxy / sqrtm);

    const double dM_dT_tr = dM_dMpot * dMpot_dT_tr;
    const double dM_f_dd = dM_dMpot * dMpot_f_dd;

    const double dPliq_dSwe = dM_dSwe;
    const double dPliq_dT_tr = (-P * dsnow_dT_tr) + dM_dT_tr;
    const double dPliq_f_dd = dM_f_dd;

    /* -------------------------------------
     * 3. Fluxes (Xinanjiang), but P -> Pliq
     * -------------------------------------- */
    const double Pi = (1.0 - A_im) * Pliq;
    const double Rb = A_im * Pliq;
    const double Ea = f_p * Ep;

    const double dPi_dPliq = (1.0 - A_im);
    const double dRb_dPliq = A_im;

    /* --- Runoff generation R(Pi,W) --- */
    const double W_max2 = W_max * W_max;

    /* uW = W/Wmax (W already smoothed) */
    const double uW = W / W_max;
    const double duW_dW = 1.0 / W_max;
    const double duW_dWmax = -W / W_max2;

    /* Smooth switch between lower/upper regimes around uW = 0.5 - a */
    const double g = uW - (0.5 - a); /* = uW - 0.5 + a */
    const double gt = g / eps_r;
    const double w = 0.5 * (1.0 + tanh(gt)); /* 0..1 */
    const double sech22 = 1.0 / (cosh(gt) * cosh(gt));
    const double dw_dg = 0.5 * sech22 / eps_r;

    const double dg_dW = duW_dW;
    const double dg_da = 1.0;
    const double dg_dWmax = duW_dWmax;

    /* --------------------------------------------------------
       Lower branch: uW <= 0.5-a
       R_low = Pi * c1 * (uW)^b, with smooth_pos(uW) protection
       -------------------------------------------------------- */
    const double base1 = 0.5 - a;
    const double base1e = fmax(base1, 1e-12); /* assume a is in valid range */
    const double c1 = pow(base1e, 1.0 - b);
    const double dc1_da = -(1.0 - b) * pow(base1e, -b); /* d/d a of (0.5-a)^(1-b) */
    const double dc1_db = -c1 * log(base1e);

    const double uWe = smooth_pos(uW, eps_m);
    const double duWe_duW = d_smooth_pos_da(uW, eps_m);

    const double uWb = pow(uWe, b);
    const double uWbm1 = pow(uWe, b - 1.0);
    const double loguW = log(uWe);

    const double Rlow = Pi * c1 * uWb;
    const double dRlow_dPi = c1 * uWb;
    const double dRlow_dW = Pi * c1 * b * uWbm1 * duWe_duW * duW_dW;
    const double dRlow_da = Pi * dc1_da * uWb; /* uW does not depend on a */
    const double dRlow_db = Pi * (dc1_db * uWb + c1 * uWb * loguW);
    const double dRlow_dWmax = Pi * c1 * b * uWbm1 * duWe_duW * duW_dWmax;

    /* -----------------------------------------------------------------
       Upper branch: uW > 0.5-a
       R_up = Pi * (1 - c2 * v^b), v = (1-uW) with smooth_pos protection
       ----------------------------------------------------------------- */
    const double base2 = 0.5 + a;
    const double base2e = fmax(base2, 1e-12);
    const double c2 = pow(base2e, 1.0 - b);
    const double dc2_da = (1.0 - b) * pow(base2e, -b);
    const double dc2_db = -c2 * log(base2e);

    const double v_raw = 1.0 - uW;
    const double v = smooth_pos(v_raw, eps_m);
    const double dv_dvraw = d_smooth_pos_da(v_raw, eps_m);

    const double dvraw_dW = -duW_dW;
    const double dvraw_dWmax = -duW_dWmax;

    const double dv_dW = dv_dvraw * dvraw_dW;
    const double dv_dWmax = dv_dvraw * dvraw_dWmax;

    const double vb = pow(v, b);
    const double vbm1 = pow(v, b - 1.0);
    const double logv = log(v);

    const double Rup = Pi * (1.0 - c2 * vb);
    const double dRup_dPi = 1.0 - c2 * vb;
    const double dRup_dW = Pi * (-c2 * b * vbm1) * dv_dW;
    const double dRup_da = Pi * (-dc2_da * vb); /* v does not depend on a */
    const double dRup_db = Pi * (-(dc2_db * vb + c2 * vb * logv));
    const double dRup_dWmax = Pi * (-c2 * b * vbm1) * dv_dWmax;

    /* -----------------------------
       Smooth blend: R = (1-w)*Rlow + w*Rup
       ----------------------------- */
    const double R = (1.0 - w) * Rlow + w * Rup;
    const double dR_dPi = (1.0 - w) * dRlow_dPi + w * dRup_dPi;
    const double dR_dW = (1.0 - w) * dRlow_dW + w * dRup_dW + dw_dg * dg_dW * (Rup - Rlow);
    const double dR_da = (1.0 - w) * dRlow_da + w * dRup_da + dw_dg * dg_da * (Rup - Rlow);
    const double dR_db = (1.0 - w) * dRlow_db + w * dRup_db;
    const double dR_dWmax =
        (1.0 - w) * dRlow_dWmax + w * dRup_dWmax + dw_dg * dg_dWmax * (Rup - Rlow);
    const double dR_dPliq = dR_dPi * dPi_dPliq; /* through Pi */

    /* --- Evaporation E --- */
    double E;
    double dE_df_p = 0.0;
    double dE_dLM = 0.0;
    double dE_dc = 0.0;
    double dE_dW = 0.0; /* dE/dW (W is smoothed already) */

    if (W > LM) {
        E = Ea;
        dE_df_p = Ep;
        dE_dW = 0.0;
    } else if (W >= c * LM) {
        E = (W / LM) * Ea;
        dE_df_p = (W / LM) * Ep;
        dE_dLM = -(W / (LM * LM)) * Ea;
        dE_dc = 0.0;
        dE_dW = Ea / LM;
    } else {
        const double E0 = c * Ea;
        // Water-limiting factor near W=0
        const double W_phi = 1;
        const double W_W_phi = W + W_phi;
        const double phiW = W / W_W_phi; // ~1 for big W, ~0 near W=0
        const double dphi_dW = W_phi / (W_W_phi * W_W_phi);
        // Scaled evaporation
        E = E0 * phiW;
        // Derivatives
        dE_dW = E0 * dphi_dW;      // since E0 independent of W here
        dE_dc = Ea * phiW;         // scale by phiW
        dE_df_p = (c * Ep) * phiW; // scale by phiW (Ea depends on f_p)
        dE_dLM = 0.0;
    }

    /* --- Partition runoff into Rs/Ri/Rg --- */
    double z = 1.0 - Sw / S_max;
    if (z < 0.0) {
        z = 0.0;
    }

    const double dum = 1.0 - pow(z, Ex);

    const double Rs = R * dum;
    const double Ri = k_i * Sw * dum;
    const double Rg = k_g * Sw * dum;

    const double Qi = c_i * Si;
    const double Qg = c_g * Sg;
    const double Qs = Rs + Rb;

    const double qf1 = k_f * Sf1;
    const double qf2 = k_f * Sf2;
    const double qf3 = k_f * Sf3;

    double ddum_dS = 0.0;
    double ddum_dSmax = 0.0;
    double ddum_dEx = 0.0;
    //
    // if (z > 0.0 && Sw < S_max) {
    //     ddum_dS = Ex * pow(z, Ex - 1.0) * (1.0 / S_max);
    //     ddum_dSmax = -Ex * pow(z, Ex - 1.0) * (Sw / (S_max * S_max));
    //     ddum_dEx = -pow(z, Ex) * log(z);
    // }

    /* Smooth gate for ddum_* derivatives */
    auto sig = [](double x) -> double { return 0.5 * (1.0 + tanh(x)); };

    /* smooth gates */
    const double g1 = sig(z / eps_r);            /* ~1 if z>0, ~0 if z<0 */
    const double g2 = sig((S_max - Sw) / eps_r); /* ~1 if Sw<S_max, ~0 if Sw>S_max */
    const double g12 = g1 * g2;

    /* Base derivatives (no gate) */
    double ddS_base = 0.0;
    double ddSmax_base = 0.0;
    double ddEx_base = 0.0;

    {
        /* protect log and pow at z ~ 0 */
        const double zlog = fmax(z, eps_m);

        /* if Ex can be < 1, z^(Ex-1) can blow up as z->0; zlog helps */
        const double z_pow_Exm1 = pow(zlog, Ex - 1.0);
        const double z_pow_Ex = z_pow_Exm1 * zlog; /* = zlog^Ex */

        ddS_base = Ex * z_pow_Exm1 * (1.0 / S_max);
        ddSmax_base = -Ex * z_pow_Exm1 * (Sw / (S_max * S_max));
        ddEx_base = -z_pow_Ex * log(zlog);
    }

    /* Apply smooth gate */
    ddum_dS = g12 * ddS_base;
    ddum_dSmax = g12 * ddSmax_base;
    ddum_dEx = g12 * ddEx_base;

    /* ------------------------------------------------------------------
     * 4) State derivatives udot
     * ------------------------------------------------------------------ */
    udot[0] = P_snow - M;
    udot[1] = Pi - E - R;
    udot[2] = R - Rs - Ri - Rg;
    udot[3] = Ri - Qi;
    udot[4] = Rg - Qg;
    udot[5] = Qs + Qi + Qg - qf1;
    udot[6] = qf1 - qf2;
    udot[7] = qf2 - qf3;
    udot[8] = qf3;

    /* ------------------------------------------------------------------
     * 5) Flux derivatives wrt states (+ SWE coupling)
     * ------------------------------------------------------------------ */
    const double dRs_dW = dR_dW * dum;
    const double dRs_dS = R * ddum_dS;

    const double dRi_dS = k_i * (dum + Sw * ddum_dS);
    const double dRg_dS = k_g * (dum + Sw * ddum_dS);

    const double dQi_dSi = c_i;
    const double dQg_dSg = c_g;

    const double dQs_dW = dRs_dW;
    const double dQs_dS = dRs_dS;

    /* SWE coupling through Pliq -> Pi, R, Rs, Qs */
    const double dPi_dSwe = dPi_dPliq * dPliq_dSwe;
    const double dRb_dSwe = dRb_dPliq * dPliq_dSwe;
    const double dR_dSwe = dR_dPliq * dPliq_dSwe;

    const double dRs_dSwe = dR_dSwe * dum;
    const double dQs_dSwe = dRs_dSwe + dRb_dSwe;

    /* ------------------------------------------------------------------
     * 6) Build Jx_f = df/dx  (m x m), m=9
     * ------------------------------------------------------------------ */
    double Jx_f[9][9] = {{0}};
    double Jth_f[9][16] = {{0}};
    double J_Wmax[9] = {0.0};
    double J_LM[9] = {0.0};
    double J_Smax[9] = {0.0};

    /* Row 0: Swe */
    Jx_f[0][0] = -dM_dSwe;

    /* Row 1: W */
    Jx_f[1][0] = dPi_dSwe - dR_dSwe; /* SWE coupling */
    Jx_f[1][1] = -dE_dW - dR_dW;     /* wrt W (smoothed) */

    /* Row 2: Sw */
    Jx_f[2][0] = dR_dSwe - dRs_dSwe;        /* SWE coupling */
    Jx_f[2][1] = dR_dW - dRs_dW;            /* wrt W (smoothed) */
    Jx_f[2][2] = -dRs_dS - dRi_dS - dRg_dS; /* wrt Sw (smoothed) */

    /* Row 3: Si */
    Jx_f[3][2] = dRi_dS;
    Jx_f[3][3] = -dQi_dSi;

    /* Row 4: Sg */
    Jx_f[4][2] = dRg_dS;
    Jx_f[4][4] = -dQg_dSg;

    /* Row 5: Sf1 */
    Jx_f[5][0] = dQs_dSwe;
    Jx_f[5][1] = dQs_dW;
    Jx_f[5][2] = dQs_dS;
    Jx_f[5][3] = dQi_dSi;
    Jx_f[5][4] = dQg_dSg;
    Jx_f[5][5] = -k_f;

    /* Row 6: Sf2 */
    Jx_f[6][5] = k_f;
    Jx_f[6][6] = -k_f;

    /* Row 7: Sf3 */
    Jx_f[7][6] = k_f;
    Jx_f[7][7] = -k_f;

    /* Row 8: Q */
    Jx_f[8][7] = k_f;

    /* map df/d(smoothed) back to df/d(raw) for Swe, W, Sw */
    for (int ii = 0; ii < m; ++ii) {
        Jx_f[ii][0] *= dSwe_dx1; /* column 1 */
        Jx_f[ii][1] *= dW_dx2;   /* column 2 */
        Jx_f[ii][2] *= dSw_dx3;  /* column 3 */
    }

    /* ------------------------------------------------------------------
     * 7) Build Jth_f = df/dtheta (m x d)
     * ------------------------------------------------------------------ */

    /* Precompute repeated */
    const double dRs_dSmax = R * ddum_dSmax;
    const double dRi_dSmax = k_i * Sw * ddum_dSmax;
    const double dRg_dSmax = k_g * Sw * ddum_dSmax;

    const double dRs_dEx = R * ddum_dEx;
    const double dRi_dEx = k_i * Sw * ddum_dEx;
    const double dRg_dEx = k_g * Sw * ddum_dEx;

    /* theta(1)=f_p */
    {
        int j = 0;
        Jth_f[1][j] = -dE_df_p;
    }

    /* theta(2)=A_im (NOTE: uses Pliq) */
    {
        int j = 1;
        const double dPi_dAim = -Pliq;
        const double dRb_dAim = Pliq;

        const double dR_dAim = dR_dPi * dPi_dAim;
        const double dRs_dAim = dR_dAim * dum;
        const double dQs_dAim = dRs_dAim + dRb_dAim;

        Jth_f[1][j] = dPi_dAim - dR_dAim;
        Jth_f[2][j] = dR_dAim - dRs_dAim;
        Jth_f[5][j] = dQs_dAim;
    }

    /* theta(3)=a */
    {
        int j = 2;
        const double dRs_da = dR_da * dum;
        Jth_f[1][j] = -dR_da;
        Jth_f[2][j] = dR_da - dRs_da;
        Jth_f[5][j] = dRs_da;
    }

    /* theta(4)=b */
    {
        int j = 3;
        const double dRs_db = dR_db * dum;
        Jth_f[1][j] = -dR_db;
        Jth_f[2][j] = dR_db - dRs_db;
        Jth_f[5][j] = dRs_db;
    }

    /* theta(5)= not W_max but f_wm */
    {
        const double dRs_dWmax = dR_dWmax * dum;

        J_Wmax[1] += -dR_dWmax;
        J_Wmax[2] += dR_dWmax - dRs_dWmax;
        J_Wmax[5] += dRs_dWmax;
    }

    /* theta(6)= not LM but f_lm */
    {
        J_LM[1] += -dE_dLM;
    }

    /* theta(7)=c */
    {
        int j = 6;
        Jth_f[1][j] = -dE_dc;
    }

    /* theta(8)= not S_max but S_tot */
    {
        J_Smax[2] += -dRs_dSmax - dRi_dSmax - dRg_dSmax;
        J_Smax[3] += dRi_dSmax;
        J_Smax[4] += dRg_dSmax;
        J_Smax[5] = dRs_dSmax;
    }

    /* theta(9)=Ex */
    {
        int j = 8;
        Jth_f[2][j] = -dRs_dEx - dRi_dEx - dRg_dEx;
        Jth_f[3][j] = dRi_dEx;
        Jth_f[4][j] = dRg_dEx;
        Jth_f[5][j] = dRs_dEx;
    }

    /* theta(10)=k_i */
    {
        int j = 9;
        const double dRi_dki = Sw * dum;
        Jth_f[2][j] = -dRi_dki;
        Jth_f[3][j] = dRi_dki;
    }

    /* theta(11)=k_g */
    {
        int j = 10;
        const double dRg_dkg = Sw * dum;
        Jth_f[2][j] = -dRg_dkg;
        Jth_f[4][j] = dRg_dkg;
    }

    /* theta(12)=c_i */
    {
        int j = 11;
        Jth_f[3][j] = -Si;
        Jth_f[5][j] = Si;
    }

    /* theta(13)=c_g */
    {
        int j = 12;
        Jth_f[4][j] = -Sg;
        Jth_f[5][j] = Sg;
    }

    /* theta(14)=k_f */
    {
        int j = 13;
        Jth_f[5][j] = -Sf1;
        Jth_f[6][j] = Sf1 - Sf2;
        Jth_f[7][j] = Sf2 - Sf3;
        Jth_f[8][j] = Sf3;
    }

    /* ------------------------------------------------------------------
     * 8) Snow parameter columns (last two cols): T_tr = d-2, f_dd = d-1
     * ------------------------------------------------------------------ */
    {
        const int j_T_tr = d - 2;
        const int j_f_dd = d - 1;

        /* derivatives wrt T_tr/f_dd through Pliq */
        const double dPi_dT_tr = dPi_dPliq * dPliq_dT_tr;
        const double dPi_f_dd = dPi_dPliq * dPliq_f_dd;

        const double dRb_dT_tr = dRb_dPliq * dPliq_dT_tr;
        const double dRb_f_dd = dRb_dPliq * dPliq_f_dd;

        const double dR_dT_tr = dR_dPliq * dPliq_dT_tr;
        const double dR_f_dd = dR_dPliq * dPliq_f_dd;

        const double dRs_dT_tr = dR_dT_tr * dum;
        const double dRs_f_dd = dR_f_dd * dum;

        const double dQs_dT_tr = dRs_dT_tr + dRb_dT_tr;
        const double dQs_f_dd = dRs_f_dd + dRb_f_dd;

        /* SWE row */
        Jth_f[0][j_T_tr] = (P * dsnow_dT_tr) - dM_dT_tr;
        Jth_f[0][j_f_dd] = -dM_f_dd;

        /* W row */
        Jth_f[1][j_T_tr] = dPi_dT_tr - dR_dT_tr;
        Jth_f[1][j_f_dd] = dPi_f_dd - dR_f_dd;

        /* Sw row */
        Jth_f[2][j_T_tr] = dR_dT_tr - dRs_dT_tr;
        Jth_f[2][j_f_dd] = dR_f_dd - dRs_f_dd;

        /* Sf1 row (via Qs = Rs + Rb) */
        Jth_f[5][j_T_tr] = dQs_dT_tr;
        Jth_f[5][j_f_dd] = dQs_f_dd;
    }

    /* ------------------------------------------------------------------
     * 9) Chain rule for f_wm, f_lm, S_tot (matches MATLAB indices)
     *     MATLAB (1-based): f_wm=5, f_lm=6, S_tot=8
     *     C (0-based):         4,      5,       7
     * ------------------------------------------------------------------ */
    {
        const int id_f_wm = 4;
        const int id_f_lm = 5;
        const int id_S_tot = 7;

        for (int i = 0; i < m; ++i) {
            // df/df_wm = (df/dWmax)*dWmax/df_wm + (df/dLM)*dLM/df_wm +
            // (df/dSmax)*dSmax/df_wm = (df/dWmax)*S_tot + (df/dLM)*(f_lm*S_tot) +
            // (df/dSmax)*(-S_tot)
            Jth_f[i][id_f_wm] +=
                S_tot * J_Wmax[i] + (f_lm * S_tot) * J_LM[i] - S_tot * J_Smax[i];

            // df/df_lm = (df/dLM) * dLM/df_lm = (df/dLM) * W_max
            Jth_f[i][id_f_lm] += W_max * J_LM[i];

            // df/dS_tot = (df/dWmax)*dWmax/dS_tot + (df/dLM)*dLM/dS_tot +
            // (df/dSmax)*dSmax/dS_tot = (df/dWmax)*f_wm + (df/dLM)*(f_lm*f_wm) +
            // (df/dSmax)*(1 - f_wm)
            Jth_f[i][id_S_tot] +=
                f_wm * J_Wmax[i] + (f_lm * f_wm) * J_LM[i] + (1.0 - f_wm) * J_Smax[i];
        }
    }

    /* ------------------------------------------------------------------
     * 10) Sensitivity update: dSdt = Jx_f*S + Jth_f
     * ------------------------------------------------------------------ */
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

namespace sage_xinanjiang {
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
    const int d = 16;
    const int m = 9;
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
    const double f_p = p.f_p;
    const double A_im = p.A_im;
    const double a = p.a;
    const double b = p.b;
    const double W_max = p.W_max;
    const double LM = p.LM;
    const double c = p.c;
    const double S_max = p.S_max;
    const double S_tot = p.S_tot;
    const double f_wm = p.f_wm;
    const double f_lm = p.f_lm;
    const double Ex = p.Ex;
    const double k_i = p.k_i;
    const double k_g = p.k_g;
    const double c_i = p.c_i;
    const double c_g = p.c_g;
    const double k_f = p.k_f;
    const double T_tr = p.T_tr;
    const double f_dd = p.f_dd;
    const double T_sm = p.T_sm;
    const double eps_m = p.eps_m;
    const double eps_r = p.eps_r;
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
} // namespace sage_xinanjiang
