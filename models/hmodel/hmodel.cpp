/*
 * hmodel.cpp
 *
 * Conceptual HMODEL rainfall-runoff model with an adaptive-step
 * explicit Runge-Kutta integrator. This file contains the MATLAB-independent
 * native numerical core shared by crr_hmodel and crr_model_mex.
 *
 * Written by Jasper A. Vrugt.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#include "hmodel.hpp"
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
                double I_max,
                double Su_max,
                double Q_max,
                double a_E,
                double a_F,
                double a_S,
                double r_f,
                double r_s,
                double T_tr,
                double f_dd,
                double T_sm,
                double eps_m,
                double rho);

static void hmodel_aug_ode(double* z,
                           double* zdot,
                           double P,
                           double Ep,
                           double T,
                           double I_max,
                           double Su_max,
                           double Q_max,
                           double a_E,
                           double a_F,
                           double a_S,
                           double r_f,
                           double r_s,
                           double T_tr,
                           double f_dd,
                           double T_sm,
                           double eps_m,
                           double rho,
                           int nvar,
                           int m,
                           int d);

static void hmodel_odefcn(double* u,
                          double* udot,
                          const double* Smat,
                          double* dSdt,
                          double P,
                          double Ep,
                          double T,
                          double I_max,
                          double Su_max,
                          double Q_max,
                          double a_E,
                          double a_F,
                          double a_S,
                          double r_f,
                          double r_s,
                          double T_tr,
                          double f_dd,
                          double T_sm,
                          double eps_m,
                          int m,
                          int d);

static inline double smooth_pos(double a, double eps)
{
    return 0.5 * (a + sqrt(a * a + eps * eps));
}

static inline double d_smooth_pos_da(double a, double eps)
{
    const double denom = sqrt(a * a + eps * eps);
    return 0.5 * (1.0 + a / denom);
}

static inline double smooth_min(double A, double B, double eps)
{
    const double d = (A - B);
    return 0.5 * (A + B - sqrt(d * d + eps * eps));
}

static inline double smooth_max(double A, double B, double eps)
{
    const double d = (A - B);
    return 0.5 * (A + B + sqrt(d * d + eps * eps));
}

static inline double d_smooth_min_dA(double A, double B, double eps)
{
    const double d = (A - B);
    return 0.5 * (1.0 - d / sqrt(d * d + eps * eps));
}

static inline double d_smooth_max_dA(double A, double B, double eps)
{
    const double d = (A - B);
    return 0.5 * (1.0 + d / sqrt(d * d + eps * eps));
}

static inline double clamp01(double z, double eps)
{
    const double zmax = smooth_max(z, 0.0, eps);
    return smooth_min(zmax, 1.0, eps);
}

static inline double d_clamp01_dz(double z, double eps)
{
    const double A = smooth_max(z, 0.0, eps);
    return d_smooth_min_dA(A, 1.0, eps) * d_smooth_max_dA(z, 0.0, eps);
}

static inline double exponen(double x)
{
    return exp(fmin(300.0, x));
}

static inline double expFlux(double Sr, double a)
{
    if (fabs(a) <= 1e-5) {
        return Sr;
    } else {
        const double num = 1.0 - exponen(-a * Sr);
        const double den = 1.0 - exponen(-a);
        return num / den;
    }
}

static inline double d_expFlux_dz(double Sr, double a)
{
    if (fabs(a) <= 1e-5) {
        return 1.0;
    }

    const double den = -expm1(-a);
    const double ez = exponen(-a * Sr);
    return a * ez / den;
}

static inline double d_expFlux_da(double Sr, double a)
{
    if (fabs(a) < 1e-6) {
        return 0.0;
    }

    const double ea = exponen(-a);
    const double ea_z = exponen(-a * Sr);

    const double N = 1.0 - ea_z;
    const double D = 1.0 - ea;
    const double Np = Sr * ea_z;
    const double Dp = ea;

    return (D * Np - N * Dp) / (D * D);
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
                double I_max,
                double Su_max,
                double Q_max,
                double a_E,
                double a_F,
                double a_S,
                double r_f,
                double r_s,
                double T_tr,
                double f_dd,
                double T_sm,
                double eps_m,
                double rho)
{
    /* Euler slope */
    hmodel_aug_ode(z,
                   zdotE,
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
                   rho,
                   nvar,
                   m,
                   d);

    for (int i = 0; i < nvar; ++i) {
        zE[i] = z[i] + h * zdotE[i];
    }

    /* Heun slope */
    hmodel_aug_ode(zE,
                   zdot,
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

static void hmodel_aug_ode(double* z,
                           double* zdot,
                           double P,
                           double Ep,
                           double T,
                           double I_max,
                           double Su_max,
                           double Q_max,
                           double a_E,
                           double a_F,
                           double a_S,
                           double r_f,
                           double r_s,
                           double T_tr,
                           double f_dd,
                           double T_sm,
                           double eps_m,
                           double rho,
                           int nvar,
                           int m,
                           int d)
{
    (void)nvar;

    double* u = z;        /* m */
    double* Smat = z + m; /* m*d */

    double* udot = zdot;     /* m */
    double* dSdt = zdot + m; /* m*d */

    hmodel_odefcn(u,
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
                  d);
}

static void hmodel_odefcn(double* u,
                          double* udot,
                          const double* Smat,
                          double* dSdt,
                          double P,
                          double Ep,
                          double T,
                          double I_max,
                          double Su_max,
                          double Q_max,
                          double a_E,
                          double a_F,
                          double a_S,
                          double r_f,
                          double r_s,
                          double T_tr,
                          double f_dd,
                          double T_sm,
                          double eps_m,
                          int m,
                          int d)
{
    /* Fixed constants (as in MATLAB) */
    const double a_I = 50.0;
    const double a_P = -50.0;

    /* --------------------------------------------------------------- */
    /* 1) Unpack states (u = [Swe, Si, Su, Sf, Ss, Q])                 */
    /* --------------------------------------------------------------- */
    const double Swe_raw = u[0];
    const double Si = u[1];
    const double Su = u[2];
    const double Sf = u[3];
    const double Ss = u[4];
    /* u[5] = Q */

    /* Only Swe is smoothed (matches MATLAB) */
    const double Swe = smooth_pos(Swe_raw, eps_m);
    const double dSwe_dx1 = d_smooth_pos_da(Swe_raw, eps_m);

    /* --------------------------------------------------------------- */
    /* 2) Snow module (smooth HBV-style degree-day)                    */
    /* --------------------------------------------------------------- */
    const double T_smeps_m = fmax(T_sm, eps_m);
    const double uT = (T - T_tr) / T_smeps_m;

    const double snow_fr = 0.5 * (1.0 - tanh(uT));
    const double rain_fr = 1.0 - snow_fr;

    const double P_snow = P * snow_fr;
    const double P_rain = P * rain_fr;

    const double a = (T - T_tr);
    const double posT = 0.5 * (a + sqrt(a * a + T_smeps_m * T_smeps_m));
    const double M_pot = f_dd * posT;

    /* Smooth min(Swe, M_pot) implemented explicitly as in MATLAB */
    const double dxy = (Swe - M_pot);
    const double sqrtm = sqrt(dxy * dxy + eps_m * eps_m);
    const double M = 0.5 * (Swe + M_pot - sqrtm);

    const double Pliq = P_rain + M;

    /* Snow derivatives for coupling (x and theta) */
    const double coshu = cosh(uT);
    const double sech2 = 1.0 / (coshu * coshu);
    const double dsnow_dT_tr = 0.5 * sech2 / T_smeps_m;

    const double dposT_da = 0.5 * (1.0 + a / sqrt(a * a + T_smeps_m * T_smeps_m));
    const double dposT_dT_tr = -dposT_da;

    const double dMpot_dT_tr = f_dd * dposT_dT_tr;
    const double dMpot_f_dd = posT;

    const double dM_dSwe = 0.5 * (1.0 - dxy / sqrtm);
    const double dM_dMpot = 0.5 * (1.0 + dxy / sqrtm);

    const double dPliq_dSwe = dM_dSwe;
    const double dPliq_dT_tr = (-P * dsnow_dT_tr) + dM_dMpot * dMpot_dT_tr;
    const double dPliq_f_dd = dM_dMpot * dMpot_f_dd;

    const double dPsnow_dT_tr = P * dsnow_dT_tr;
    const double dM_dT_tr = dM_dMpot * dMpot_dT_tr;
    const double dM_f_dd = dM_dMpot * dMpot_f_dd;

    /* --------------------------------------------------------------- */
    /* 3) Interception module                                          */
    /* --------------------------------------------------------------- */
    double EvapI = 0.0, P_e = Pliq, Ep_e = Ep;

    double dEvapI_dSi = 0.0, dPe_dSi = 0.0, dEp_e_dSi = 0.0;
    double dEvapI_dImax = 0.0, dPe_dImax = 0.0, dEp_e_dImax = 0.0;

    double dPe_dSwe = dPliq_dSwe; /* overridden when I_max>0 */
    double G_I = 1.0;             /* needed later for dPe/dT_tr, dPe/f_dd */

    if (I_max > 0.0) {
        const double Sir_raw = Si / I_max;
        const double Sir = clamp01(Sir_raw, eps_m);

        const double dSir_dSirraw = d_clamp01_dz(Sir_raw, eps_m);
        const double dSir_dSi = dSir_dSirraw * (1.0 / I_max);
        const double dSir_dImax = dSir_dSirraw * (-Si / (I_max * I_max));

        const double F_I = expFlux(Sir, a_I);
        G_I = expFlux(Sir, a_P);

        EvapI = Ep * F_I;
        P_e = Pliq * G_I;

        const double dF_I_dz = d_expFlux_dz(Sir, a_I);
        const double dG_I_dz = d_expFlux_dz(Sir, a_P);

        dEvapI_dSi = Ep * dF_I_dz * dSir_dSi;
        dEvapI_dImax = Ep * dF_I_dz * dSir_dImax;

        dPe_dSi = Pliq * dG_I_dz * dSir_dSi;
        dPe_dImax = Pliq * dG_I_dz * dSir_dImax;

        /* Ep_e = smooth_pos(Ep - EvapI) */
        {
            const double arg = Ep - EvapI;
            const double dEp_e_darg = d_smooth_pos_da(arg, eps_m);
            Ep_e = smooth_pos(arg, eps_m);

            dEp_e_dSi = dEp_e_darg * (-dEvapI_dSi);
            dEp_e_dImax = dEp_e_darg * (-dEvapI_dImax);
        }

        /* snow coupling */
        dPe_dSwe = G_I * dPliq_dSwe;
    } else {
        /* EvapI=0, P_e=Pliq, Ep_e=Ep already set */
        dPe_dSwe = dPliq_dSwe;
    }

    /* --------------------------------------------------------------- */
    /* 4) Unsaturated zone module (uses P_e and Ep_e)                  */
    /* --------------------------------------------------------------- */
    double Sur = 0.0;

    double dEa_dSu = 0.0, dprc_dSu = 0.0, drnf_dSu = 0.0;
    double dEa_dSumax = 0.0, dprc_dSumax = 0.0, drnf_dSumax = 0.0;

    if (Su_max > 0.0) {
        const double Sur_raw = Su / Su_max;
        Sur = clamp01(Sur_raw, eps_m);

        const double dSur_dSurraw = d_clamp01_dz(Sur_raw, eps_m);
        const double dSur_dSu = dSur_dSurraw * (1.0 / Su_max);
        const double dSur_dSumax = dSur_dSurraw * (-Su / (Su_max * Su_max));

        const double dF_E_dz = d_expFlux_dz(Sur, a_E);
        const double dF_S_dz = d_expFlux_dz(Sur, a_S);
        const double dF_F_dz = d_expFlux_dz(Sur, a_F);

        dEa_dSu = Ep_e * dF_E_dz * dSur_dSu;
        dprc_dSu = Q_max * dF_S_dz * dSur_dSu;
        drnf_dSu = P_e * dF_F_dz * dSur_dSu;

        dEa_dSumax = Ep_e * dF_E_dz * dSur_dSumax;
        dprc_dSumax = Q_max * dF_S_dz * dSur_dSumax;
        drnf_dSumax = P_e * dF_F_dz * dSur_dSumax;
    } else {
        Sur = 0.0;
    }

    const double F_E = expFlux(Sur, a_E);
    const double F_S = expFlux(Sur, a_S);
    const double F_F = expFlux(Sur, a_F);

    const double Ea = Ep_e * F_E;
    const double prc = Q_max * F_S;
    const double rnf = P_e * F_F;

    const double qf = Sf / r_f;
    const double qs = Ss / r_s;

    /* Helpful couplings used in Jx_f */
    const double dEa_dSi = dEp_e_dSi * F_E;
    const double drnf_dSi = dPe_dSi * F_F;
    const double drnf_dSwe = dPe_dSwe * F_F;

    const double dqf_dSf = 1.0 / r_f;
    const double dqs_dSs = 1.0 / r_s;

    /* --------------------------------------------------------------- */
    /* 5) ODE RHS (m=6)                                                */
    /* --------------------------------------------------------------- */
    udot[0] = P_snow - M;
    udot[1] = Pliq - EvapI - P_e;
    udot[2] = P_e - Ea - prc - rnf;
    udot[3] = rnf - qf;
    udot[4] = prc - qs;
    udot[5] = qf + qs;

    /* --------------------------------------------------------------- */
    /* 6) Jacobians: Jx_f (m x m) and Jth_f (m x d)                    */
    /* --------------------------------------------------------------- */
    double Jx_f[6][6] = {{0}};
    double Jth_f[6][9] = {{0}};

    /* Jx_f */
    Jx_f[0][0] = -dM_dSwe;

    Jx_f[1][0] = dPliq_dSwe - dPe_dSwe;
    Jx_f[1][1] = -dEvapI_dSi - dPe_dSi;

    Jx_f[2][0] = dPe_dSwe - drnf_dSwe;
    /* IMPORTANT: match your final MATLAB form */
    Jx_f[2][1] = dPe_dSi * (1.0 - F_F) - F_E * dEp_e_dSi;
    Jx_f[2][2] = -dEa_dSu - dprc_dSu - drnf_dSu;

    Jx_f[3][0] = drnf_dSwe;
    Jx_f[3][1] = drnf_dSi;
    Jx_f[3][2] = drnf_dSu;
    Jx_f[3][3] = -dqf_dSf;

    Jx_f[4][2] = dprc_dSu;
    Jx_f[4][4] = -dqs_dSs;

    Jx_f[5][3] = dqf_dSf;
    Jx_f[5][4] = dqs_dSs;

    /* PATCH: map df/d(Swe_smooth) back to df/d(Swe_raw) */
    for (int ii = 0; ii < 6; ++ii) {
        Jx_f[ii][0] *= dSwe_dx1;
    }

    /* Jth_f: theta=[Imax, Sumax, Qmax, aE, aF, rf, rs, T_tr, f_dd] */
    /* 1) I_max */
    {
        const double dEa_dImax = dEp_e_dImax * F_E;
        const double drnf_dImax = dPe_dImax * F_F;
        const int j = 0;
        Jth_f[1][j] = -dEvapI_dImax - dPe_dImax;
        Jth_f[2][j] = dPe_dImax - dEa_dImax - drnf_dImax;
        Jth_f[3][j] = drnf_dImax;
    }

    /* 2) Su_max */
    {
        const int j = 1;
        Jth_f[2][j] = -dEa_dSumax - dprc_dSumax - drnf_dSumax;
        Jth_f[3][j] = drnf_dSumax;
        Jth_f[4][j] = dprc_dSumax;
    }

    /* 3) Q_max */
    {
        const int j = 2;
        const double dprc_dQmax = F_S;
        Jth_f[2][j] = -dprc_dQmax;
        Jth_f[4][j] = dprc_dQmax;
    }

    /* 4) a_E */
    {
        const int j = 3;
        const double dF_E_da = d_expFlux_da(Sur, a_E);
        const double dEa_daE = Ep_e * dF_E_da;
        Jth_f[2][j] = -dEa_daE;
    }

    /* 5) a_F */
    {
        const int j = 4;
        const double dF_F_da = d_expFlux_da(Sur, a_F);
        const double drnf_daF = P_e * dF_F_da;
        Jth_f[2][j] = -drnf_daF;
        Jth_f[3][j] = drnf_daF;
    }

    /* 6) r_f */
    {
        const int j = 5;
        Jth_f[3][j] = Sf / (r_f * r_f);
        Jth_f[5][j] = -Sf / (r_f * r_f);
    }

    /* 7) r_s */
    {
        const int j = 6;
        Jth_f[4][j] = Ss / (r_s * r_s);
        Jth_f[5][j] = -Ss / (r_s * r_s);
    }

    /* 8) T_tr */
    {
        const int j = 7;
        const double dPe_dT_tr = (I_max > 0.0) ? (G_I * dPliq_dT_tr) : dPliq_dT_tr;
        const double drnf_dT_tr = dPe_dT_tr * F_F;

        Jth_f[0][j] = dPsnow_dT_tr - dM_dT_tr;
        Jth_f[1][j] = dPliq_dT_tr - dPe_dT_tr;
        Jth_f[2][j] = dPe_dT_tr - drnf_dT_tr;
        Jth_f[3][j] = drnf_dT_tr;
    }

    /* 9) f_dd */
    {
        const int j = 8;
        const double dPe_f_dd = (I_max > 0.0) ? (G_I * dPliq_f_dd) : dPliq_f_dd;
        const double drnf_f_dd = dPe_f_dd * F_F;

        Jth_f[0][j] = -dM_f_dd;
        Jth_f[1][j] = dPliq_f_dd - dPe_f_dd;
        Jth_f[2][j] = dPe_f_dd - drnf_f_dd;
        Jth_f[3][j] = drnf_f_dd;
    }

    /* --------------------------------------------------------------- */
    /* 7) Sensitivity update: dS/dt = Jx_f*S + Jth_f (Smat col-major)  */
    /* --------------------------------------------------------------- */
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

namespace sage_hmodel {
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
    const int d = 9;
    const int m = 6;
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
    const double I_max = p.I_max;
    const double Su_max = p.Su_max;
    const double Q_max = p.Q_max;
    const double a_E = p.a_E;
    const double a_F = p.a_F;
    const double a_S = p.a_S;
    const double r_f = p.r_f;
    const double r_s = p.r_s;
    const double T_tr = p.T_tr;
    const double f_dd = p.f_dd;
    const double T_sm = p.T_sm;
    const double eps_m = p.eps_m;
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
} // namespace sage_hmodel
