/*
 * hymod.cpp
 *
 * Conceptual HYMOD rainfall-runoff model with an adaptive-step
 * explicit Runge-Kutta integrator. This file contains the MATLAB-independent
 * native numerical core shared by crr_hymod and crr_model_mex.
 *
 * Written by Jasper A. Vrugt.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#include "hymod.hpp"
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
                double S_umax,
                double beta,
                double alfa,
                double K_s,
                double K_f,
                double T_tr,
                double f_dd,
                double T_sm,
                double eps_m,
                double rho);

static void hymod_aug_ode(double* z,
                          double* zdot,
                          double P,
                          double Ep,
                          double T,
                          double S_umax,
                          double beta,
                          double alfa,
                          double K_s,
                          double K_f,
                          double T_tr,
                          double f_dd,
                          double T_sm,
                          double eps_m,
                          double rho,
                          int nvar,
                          int m,
                          int d);

static void hymod_odefcn(double* u,
                         double* udot,
                         const double* Smat,
                         double* dSdt,
                         double P,
                         double Ep,
                         double T,
                         double S_umax,
                         double beta,
                         double alfa,
                         double K_s,
                         double K_f,
                         double T_tr,
                         double f_dd,
                         double T_sm,
                         double eps_m,
                         double rho,
                         int nvar,
                         int m,
                         int d);

static inline double smooth_pos(const double a, const double ep)
{
    return 0.5 * (a + sqrt(a * a + ep * ep));
}

static inline double dsmooth_pos_da(const double a, const double ep)
{
    return 0.5 * (1.0 + a / sqrt(a * a + ep * ep));
}

static inline double smooth_neg(const double a, const double ep)
{
    return 0.5 * (-a + sqrt(a * a + ep * ep));
}

static inline double dsmooth_neg_da(const double a, const double ep)
{
    return 0.5 * (-1.0 + a / sqrt(a * a + ep * ep));
}

static inline double smooth_clamp01(const double z, const double ep, double* d_dz)
{
    const double sp = smooth_pos(z - 1.0, ep);
    const double sn = smooth_neg(z, ep);
    if (d_dz) {
        const double dsp = dsmooth_pos_da(z - 1.0, ep);
        const double dsn = dsmooth_neg_da(z, ep);
        *d_dz = 1.0 - dsp + dsn;
    }
    return z - sp + sn;
}

static inline double
smooth_max2(const double a, const double b, const double ep, double* d_da /* d/da */)
{
    const double diff = a - b;
    const double s = sqrt(diff * diff + ep * ep);
    if (d_da) {
        *d_da = 0.5 * (1.0 + diff / s);
    }
    return 0.5 * (a + b + s);
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
                double S_umax,
                double beta,
                double alfa,
                double K_s,
                double K_f,
                double T_tr,
                double f_dd,
                double T_sm,
                double eps_m,
                double rho)
{
    /* Euler */
    hymod_aug_ode(z,
                  zdotE,
                  P,
                  Ep,
                  T,
                  S_umax,
                  beta,
                  alfa,
                  K_s,
                  K_f,
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
    hymod_aug_ode(zE,
                  zdot,
                  P,
                  Ep,
                  T,
                  S_umax,
                  beta,
                  alfa,
                  K_s,
                  K_f,
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

static void hymod_aug_ode(double* z,
                          double* zdot,
                          double P,
                          double Ep,
                          double T,
                          double S_umax,
                          double beta,
                          double alfa,
                          double K_s,
                          double K_f,
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

    hymod_odefcn(u,
                 udot,
                 Smat,
                 dSdt,
                 P,
                 Ep,
                 T,
                 S_umax,
                 beta,
                 alfa,
                 K_s,
                 K_f,
                 T_tr,
                 f_dd,
                 T_sm,
                 eps_m,
                 rho,
                 nvar,
                 m,
                 d);
}

static void hymod_odefcn(double* u,
                         double* udot,
                         const double* Smat,
                         double* dSdt,
                         double P,
                         double Ep,
                         double T,
                         double s_umax,
                         double beta,
                         double alfa,
                         double K_s,
                         double K_f,
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
    (void)m;
    (void)d;

    /* MATLAB constants inside the function */
    const double eps_c = 1e-6;    /* clamp smoothing for Su/s_umax */
    const double def0 = 1e-6;     /* deficit floor */
    const double eps_def = 1e-10; /* smoothing for deficit floor */

    /* ------------------------------------------------------------------ */
    /* 1) Unpack raw states                                                */
    /* ------------------------------------------------------------------ */
    const double Swe_u = u[0];
    const double Su = u[1];
    const double Ss = u[2];
    const double Sf1 = u[3];
    const double Sf2 = u[4];
    const double Sf3 = u[5];
    /* u[6] = Qinf */

    /* ------------------------------------------------------------------ */
    /* 2) SWE smooth positivity                                           */
    /* ------------------------------------------------------------------ */
    const double Swe = smooth_pos(Swe_u, eps_m);
    const double dSwe_dx = dsmooth_pos_da(Swe_u, eps_m);

    /* ------------------------------------------------------------------ */
    /* 3) Snow module                                                     */
    /* ------------------------------------------------------------------ */
    const double T_smeps_m = fmax(T_sm, eps_m);
    const double uT = (T - T_tr) / T_smeps_m;

    const double snow_fr = 0.5 * (1.0 - tanh(uT));
    const double rain_fr = 1.0 - snow_fr;
    const double P_snow = P * snow_fr;
    const double P_rain = P * rain_fr;

    const double aT = (T - T_tr);
    const double posT = 0.5 * (aT + sqrt(aT * aT + T_smeps_m * T_smeps_m));
    const double M_pot = f_dd * posT;

    const double dxy = (Swe - M_pot);
    const double sqrtm = sqrt(dxy * dxy + eps_m * eps_m);
    const double M = 0.5 * (Swe + M_pot - sqrtm);
    const double Pliq = P_rain + M;

    /* ------------------------------------------------------------------ */
    /* 4) Production store                                                */
    /* ------------------------------------------------------------------ */
    const double Su_over = Su / s_umax;

    double dsur_dz = 0.0;
    const double s_ur = smooth_clamp01(Su_over, eps_c, &dsur_dz);

    const double ds_ur_dSu = dsur_dz * (1.0 / s_umax);
    const double ds_ur_ds_umax = dsur_dz * (-Su / (s_umax * s_umax));

    const double def_raw = 1.0 - s_ur;

    double ddef_ddefraw = 0.0; /* d(def)/d(def_raw) */
    const double def = smooth_max2(def_raw, def0, eps_def, &ddef_ddefraw);

    const double ddef_ds_ur = -ddef_ddefraw;

    const double logdef = log(def);
    const double def_beta = exp(beta * logdef);
    const double A = 1.0 - def_beta;
    const double qu = Pliq * A;

    const double Ea = Ep * s_ur * (1.0 + rho) / (s_ur + rho);

    /* routing */
    const double qs = (1.0 - alfa) * qu;
    const double qs_o = K_s * Ss;

    const double qf = alfa * qu;
    const double qf_o1 = K_f * Sf1;
    const double qf_o2 = K_f * Sf2;
    const double qf_o3 = K_f * Sf3;

    /* ------------------------------------------------------------------ */
    /* 5) ODE RHS                                                         */
    /* ------------------------------------------------------------------ */
    udot[0] = P_snow - M;
    udot[1] = Pliq - qu - Ea;
    udot[2] = qs - qs_o;
    udot[3] = qf - qf_o1;
    udot[4] = qf_o1 - qf_o2;
    udot[5] = qf_o2 - qf_o3;
    udot[6] = qf_o3 + qs_o;

    /* ------------------------------------------------------------------ */
    /* 6) Jacobians                                                       */
    /* ------------------------------------------------------------------ */
    double Jx_f[7][7] = {{0}};
    double Jth_f[7][7] = {{0}};

    /* --- snow derivatives needed for coupling (MATCH MATLAB) --- */
    const double sech2 = 1.0 / (cosh(uT) * cosh(uT));
    const double dsnowFrac_dT_tr = 0.5 * sech2 / T_smeps_m;

    const double denom_pos = sqrt(aT * aT + T_smeps_m * T_smeps_m);
    const double dposT_da = 0.5 * (1.0 + aT / denom_pos);
    const double dposT_dT_tr = -dposT_da;

    const double dMpot_dT_tr = f_dd * dposT_dT_tr;

    const double dM_dSWE = 0.5 * (1.0 - dxy / sqrtm);
    const double dM_dMpot = 0.5 * (1.0 + dxy / sqrtm);

    const double dPliq_dSWE = dM_dSWE;

    /* qu derivatives */
    const double dqu_dPliq = A;

    const double dA_ddef = -(def_beta) * (beta / def);
    const double dA_ds_ur = dA_ddef * ddef_ds_ur;

    const double dqu_ds_ur = Pliq * dA_ds_ur;
    const double dqu_dSu = dqu_ds_ur * ds_ur_dSu;

    const double dEa_ds_ur = Ep * (1.0 + rho) * rho / ((s_ur + rho) * (s_ur + rho));
    const double dEa_dSu = dEa_ds_ur * ds_ur_dSu;

    const double dqu_dSWE = dqu_dPliq * dPliq_dSWE;

    /* --- Jx_f = df/dx, initially w.r.t. SWE (smooth) --- */
    Jx_f[0][0] = -dM_dSWE;

    Jx_f[1][0] = dPliq_dSWE - dqu_dSWE;
    Jx_f[1][1] = -dqu_dSu - dEa_dSu;

    Jx_f[2][0] = (1.0 - alfa) * dqu_dSWE;
    Jx_f[2][1] = (1.0 - alfa) * dqu_dSu;
    Jx_f[2][2] = -K_s;

    Jx_f[3][0] = alfa * dqu_dSWE;
    Jx_f[3][1] = alfa * dqu_dSu;
    Jx_f[3][3] = -K_f;

    Jx_f[4][3] = K_f;
    Jx_f[4][4] = -K_f;

    Jx_f[5][4] = K_f;
    Jx_f[5][5] = -K_f;

    Jx_f[6][2] = K_s;
    Jx_f[6][5] = K_f;

    /* Map df/d(Swe_smooth) -> df/d(Swe_raw) (MATCH MATLAB: Jx(:,1)*=dSwe_dx1) */
    for (int i = 0; i < 7; ++i) {
        Jx_f[i][0] *= dSwe_dx;
    }

    /* --- Jth_f = df/dtheta (MATCH MATLAB ordering) --- */

    /* param 1: s_umax */
    const double dqu_ds_umax = dqu_ds_ur * ds_ur_ds_umax;
    const double dEa_ds_umax = dEa_ds_ur * ds_ur_ds_umax;

    Jth_f[1][0] = -dqu_ds_umax - dEa_ds_umax;
    Jth_f[2][0] = (1.0 - alfa) * dqu_ds_umax;
    Jth_f[3][0] = alfa * dqu_ds_umax;

    /* param 2: beta */
    const double dqu_dbeta = -Pliq * def_beta * logdef;

    Jth_f[1][1] = -dqu_dbeta;
    Jth_f[2][1] = (1.0 - alfa) * dqu_dbeta;
    Jth_f[3][1] = alfa * dqu_dbeta;

    /* param 3: alfa */
    Jth_f[2][2] = -qu;
    Jth_f[3][2] = qu;

    /* param 4: Ks */
    Jth_f[2][3] = -Ss;
    Jth_f[6][3] = Ss;

    /* param 5: Kf */
    Jth_f[3][4] = -Sf1;
    Jth_f[4][4] = Sf1 - Sf2;
    Jth_f[5][4] = Sf2 - Sf3;
    Jth_f[6][4] = Sf3;

    /* param 6: T_tr */
    const double dPsnow_dT_tr = P * dsnowFrac_dT_tr;
    const double dPrain_dT_tr = -P * dsnowFrac_dT_tr;

    const double dM_dT_tr = dM_dMpot * dMpot_dT_tr;
    const double dPliq_dT_tr = dPrain_dT_tr + dM_dT_tr;
    const double dqu_dT_tr = dqu_dPliq * dPliq_dT_tr; /* A * dPliq_dT_tr */

    Jth_f[0][5] = dPsnow_dT_tr - dM_dT_tr;
    Jth_f[1][5] = dPliq_dT_tr - dqu_dT_tr;
    Jth_f[2][5] = (1.0 - alfa) * dqu_dT_tr;
    Jth_f[3][5] = alfa * dqu_dT_tr;

    /* param 7: f_dd */
    const double dMpot_df_dd = posT;
    const double dM_df_dd = dM_dMpot * dMpot_df_dd;

    const double dPliq_df_dd = dM_df_dd;
    const double dqu_df_dd = dqu_dPliq * dPliq_df_dd; /* A * dPliq_df_dd */

    Jth_f[0][6] = -dM_df_dd;
    Jth_f[1][6] = dPliq_df_dd - dqu_df_dd;
    Jth_f[2][6] = (1.0 - alfa) * dqu_df_dd;
    Jth_f[3][6] = alfa * dqu_df_dd;

    /* ------------------------------------------------------------------ */
    /* 7) Sensitivity update: dS/dt = Jx*S + Jth                          */
    /* Smat, dSdt are 7x7 column-major                                    */
    /* ------------------------------------------------------------------ */
    for (int j = 0; j < 7; ++j) {
        for (int i = 0; i < 7; ++i) {
            dSdt[i + 7 * j] = Jth_f[i][j];
        }
    }

    for (int i = 0; i < 7; ++i) {
        for (int k = 0; k < 7; ++k) {
            const double coefficient = Jx_f[i][k];
            if (coefficient == 0.0) {
                continue;
            }
            for (int j = 0; j < 7; ++j) {
                dSdt[i + 7 * j] += coefficient * Smat[k + 7 * j];
            }
        }
    }
}

namespace sage_hymod {
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
    const int d = 7;
    const int m = 7;
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
    const double S_umax = p.S_umax;
    const double beta = p.beta;
    const double alfa = p.alfa;
    const double K_s = p.K_s;
    const double K_f = p.K_f;
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
                S_umax,
                beta,
                alfa,
                K_s,
                K_f,
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
} // namespace sage_hymod
