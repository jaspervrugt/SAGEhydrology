/*
 * hbv.cpp
 *
 * Conceptual HBV rainfall-runoff model with an adaptive-step explicit
 * Runge-Kutta integrator. This file contains the MATLAB-independent native
 * numerical core shared by crr_hbv and crr_model_mex.
 *
 * Written by Jasper A. Vrugt, February 2022.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#include "hbv.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace sage_hbv {

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
                double f_c,
                double beta,
                double lp,
                double k_0,
                double uzl,
                double k_1,
                double k_2,
                double perc,
                double T_tr,
                double f_dd,
                double sfcf,
                double cfr,
                double eps_s,
                double eps_t,
                double eps_x);

static void hbv_aug_ode(double* z,
                        double* zdot,
                        double P,
                        double Ep,
                        double T,
                        double f_c,
                        double beta,
                        double lp,
                        double k_0,
                        double uzl,
                        double k_1,
                        double k_2,
                        double perc,
                        double T_tr,
                        double f_dd,
                        double sfcf,
                        double cfr,
                        double eps_s,
                        double eps_t,
                        double eps_x,
                        int m,
                        int d);

static void hbv_odefcn(double* u,
                       double* udot,
                       const double* Smat,
                       double* dSdt,
                       double P,
                       double Ep,
                       double T,
                       double f_c,
                       double beta,
                       double lp,
                       double k_0,
                       double uzl,
                       double k_1,
                       double k_2,
                       double perc,
                       double T_tr,
                       double f_dd,
                       double sfcf,
                       double cfr,
                       double eps_s,
                       double eps_t,
                       double eps_x,
                       int m,
                       int d);

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
    return 0.5 * (1.0 - d / sqrt(d * d + eps * eps));
}

static inline double d_smooth_min_dB(double A, double B, double eps)
{
    double d = A - B;
    return 0.5 * (1.0 + d / sqrt(d * d + eps * eps));
}

static inline double smooth_max(double A, double B, double eps)
{
    double d = A - B;
    return 0.5 * (A + B + sqrt(d * d + eps * eps));
}

static inline double d_smooth_max_dA(double A, double B, double eps)
{
    double d = A - B;
    return 0.5 * (1.0 + d / sqrt(d * d + eps * eps));
}

/* clamp01(z) = min(max(z,0),1) with consistent derivative */
static inline double clamp01(double z, double eps)
{
    return smooth_min(smooth_max(z, 0.0, eps), 1.0, eps);
}

static inline double dclamp01_dz(double z, double eps)
{
    double A = smooth_max(z, 0.0, eps);
    double dmin_dA = d_smooth_min_dA(A, 1.0, eps);
    double dmax_dz = d_smooth_max_dA(z, 0.0, eps);
    return dmin_dA * dmax_dz;
}

/* evap limiter: phiE = smooth_min(ratio, 1) */
static inline double phiE_fun(double ratio, double eps)
{
    return smooth_min(ratio, 1.0, eps);
}

static inline double dphiE_dratio(double ratio, double eps)
{
    return d_smooth_min_dA(ratio, 1.0, eps);
}

/* safe inverse cosh^2 (sech^2) */
static inline double sech2(double x)
{
    const double ax = fabs(x);
    if (ax > 350.0) {
        return 0.0;
    }
    const double e = exp(-2.0 * ax);
    return 4.0 * e / ((1.0 + e) * (1.0 + e));
}

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
    constexpr int m = 5;
    constexpr int d = 12;
    constexpr int nvar = m + m * d;

    if (ns < 1 || !z0 || nz0 != static_cast<std::size_t>(nvar)) {
        return true;
    }

    if (!forcing.P || !forcing.Ep || !forcing.T ||
        forcing.n < static_cast<std::size_t>(ns)) {
        return true;
    }

    if (ipr < 1) {
        ipr = 1;
    }
    if (ipr > ns + 1) {
        ipr = ns + 1;
    }

    const int nt = ns + 1;
    const int n_q = (!mem && ipr <= ns) ? (ns - ipr + 1) : 0;

    const std::size_t expected_zrows = mem ? static_cast<std::size_t>(nt) : 1u;
    if (!out.Z || out.zrows != expected_zrows ||
        out.zcols != static_cast<std::size_t>(nvar)) {
        return true;
    }

    if (!mem && n_q > 0) {
        if (!out.q || out.nq != static_cast<std::size_t>(n_q)) {
            return true;
        }
        if (needJ && (!out.J || out.nj != static_cast<std::size_t>(d))) {
            return true;
        }
    }

    std::fill(out.Z, out.Z + out.zrows * out.zcols, 0.0);
    if (!mem && n_q > 0) {
        std::fill(out.q, out.q + out.nq, 0.0);
        if (needJ) {
            std::fill(out.J, out.J + out.nq * out.nj, 0.0);
        }
    }

    std::vector<double> LTE(nvar), ztmp(nvar), w(nvar), zdotE(nvar), zE(nvar), zdot(nvar);
    std::vector<double> zcur(z0, z0 + nvar);
    std::vector<double> prev;

    int last_stored_row = 0;

    if (mem) {
        for (int j = 0; j < nvar; ++j) {
            out.Z[static_cast<std::size_t>(nt) * j] = zcur[j];
        }
    } else {
        prev.assign(z0, z0 + nvar);
    }

    bool fail = false;
    int warningPrinted = 0;
    int iterCount = 0;
    int k_out = 0;

    double hCarry = std::max(opt.MinStep, std::min(opt.InitStep, opt.MaxStep));

    for (int s = 1; s <= ns; ++s) {
        const double t1 = static_cast<double>(s - 1);
        const double t2 = static_cast<double>(s);
        double tcur = t1;
        double h = std::min(hCarry, t2 - t1);

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
                forcing.P[s - 1],
                forcing.Ep[s - 1],
                forcing.T[s - 1],
                p.f_c,
                p.beta,
                p.lp,
                p.k_0,
                p.uzl,
                p.k_1,
                p.k_2,
                p.perc,
                p.T_tr,
                p.f_dd,
                p.sfcf,
                p.cfr,
                p.eps_s,
                p.eps_t,
                p.eps_x);

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
                w[ii] = 1.0 / (opt.RelTol * std::abs(ztmp[ii]) + opt.AbsTol);
                const double wl = w[ii] * LTE[ii];
                sumWL += wl * wl;
            }

            const double wrms = std::sqrt(sumWL / static_cast<double>(nvar));
            const bool accepted = (wrms <= 1.0) || (h <= opt.MinStep);

            if (accepted) {
                zcur.swap(ztmp);
                tcur += h;
                iterCount = 0;
            } else {
                ++iterCount;
            }

            double fac;
            if (wrms <= 0.0) {
                fac = 5.0;
            } else {
                fac = 0.9 * std::pow(wrms, -1.0 / opt.Order);
                fac = std::min(5.0, std::max(0.2, fac));
            }

            double hNext = h * fac;
            hNext = std::max(opt.MinStep, std::min(hNext, opt.MaxStep));
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

            if (iterCount >= opt.maxiter) {
                if (!warningPrinted) {
                    std::fprintf(stderr,
                                 "WARNING: HBV max rejection limit reached at t = %.5f\\n",
                                 tcur);
                    warningPrinted = 1;
                }
                fail = true;
                break;
            }
        }

        if (fail) {
            break;
        }

        if (mem) {
            for (int j = 0; j < nvar; ++j) {
                out.Z[static_cast<std::size_t>(s) + static_cast<std::size_t>(nt) * j] =
                    zcur[j];
            }
            last_stored_row = s;
        } else {
            for (int j = 0; j < nvar; ++j) {
                out.Z[static_cast<std::size_t>(j)] = zcur[j];
            }

            if (s >= ipr && k_out < n_q) {
                out.q[static_cast<std::size_t>(k_out)] = zcur[m - 1] - prev[m - 1];

                if (needJ) {
                    for (int j = 0; j < d; ++j) {
                        const int idx = (j + 2) * m - 1;
                        out.J[static_cast<std::size_t>(k_out) +
                              static_cast<std::size_t>(n_q) * j] = zcur[idx] - prev[idx];
                    }
                }
                ++k_out;
            }

            prev = zcur;
        }
    }

    if (fail && mem) {
        for (int row = last_stored_row + 1; row < nt; ++row) {
            for (int j = 0; j < nvar; ++j) {
                out.Z[static_cast<std::size_t>(row) + static_cast<std::size_t>(nt) * j] =
                    out.Z[static_cast<std::size_t>(last_stored_row) +
                          static_cast<std::size_t>(nt) * j];
            }
        }
    }

    return fail;
}

Result run(int ns,
           const std::vector<double>& z0,
           const Forcing& forcing,
           const Params& p,
           const Options& opt,
           bool mem,
           int ipr,
           bool needJ)
{
    constexpr int d = 12;
    constexpr int nvar = 65;

    Result r;
    if (ns < 1 || z0.size() != static_cast<std::size_t>(nvar)) {
        r.fail = true;
        return r;
    }

    if (ipr < 1) {
        ipr = 1;
    }
    if (ipr > ns + 1) {
        ipr = ns + 1;
    }

    const int nt = ns + 1;
    const int nq = (!mem && ipr <= ns) ? (ns - ipr + 1) : 0;

    r.zrows = mem ? static_cast<std::size_t>(nt) : 1u;
    r.zcols = nvar;
    r.Z.resize(r.zrows * r.zcols);

    if (!mem && nq > 0) {
        r.q.resize(static_cast<std::size_t>(nq));
        if (needJ) {
            r.J.resize(static_cast<std::size_t>(nq) * d);
        }
    }

    OutputView ov;
    ov.Z = r.Z.data();
    ov.q = r.q.empty() ? nullptr : r.q.data();
    ov.J = r.J.empty() ? nullptr : r.J.data();
    ov.zrows = r.zrows;
    ov.zcols = r.zcols;
    ov.nq = r.q.size();
    ov.nj = needJ ? d : 0;

    r.fail = run_into(ns, z0.data(), z0.size(), forcing, p, opt, mem, ipr, needJ, ov);

    return r;
}

/* ------------------------------------------------------------------ */
/* RK2 step (Heun) + LTE estimate                                     */
/* ------------------------------------------------------------------ */
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
                double f_c,
                double beta,
                double lp,
                double k_0,
                double uzl,
                double k_1,
                double k_2,
                double perc,
                double T_tr,
                double f_dd,
                double sfcf,
                double cfr,
                double eps_s,
                double eps_t,
                double eps_x)
{
    /* Euler solution */
    hbv_aug_ode(z,
                zdotE,
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
                d);

    for (int i = 0; i < nvar; ++i) {
        zE[i] = z[i] + h * zdotE[i];
    }
    /* Heun solution */
    hbv_aug_ode(zE,
                zdot,
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

/* 4. hbv augmented ode with sensitivities as state variables */
static void hbv_aug_ode(double* z,
                        double* zdot,
                        double P,
                        double Ep,
                        double T,
                        double f_c,
                        double beta,
                        double lp,
                        double k_0,
                        double uzl,
                        double k_1,
                        double k_2,
                        double perc,
                        double T_tr,
                        double f_dd,
                        double sfcf,
                        double cfr,
                        double eps_s,
                        double eps_t,
                        double eps_x,
                        int m,
                        int d)
{
    double* u = z;        /* m */
    double* Smat = z + m; /* m*d */

    double* udot = zdot;     /* m */
    double* dSdt = zdot + m; /* m*d */

    /* Call your hbv ODE + sensitivity routine */
    hbv_odefcn(u,
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
               d);
}

/* ------------------------------------------------------------------ */
/* HBV ODE + sensitivities (combined): mirrors MATLAB hbv_odefcn      */
/* PATCH: smooth_pos clamps + column chain rule for Jx_f              */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
static void hbv_odefcn(double* u,
                       double* udot,
                       const double* Smat,
                       double* dSdt,
                       double P,
                       double Ep,
                       double T,
                       double f_c,
                       double beta,
                       double lp,
                       double k_0,
                       double uzl,
                       double k_1,
                       double k_2,
                       double perc,
                       double T_tr,
                       double f_dd,
                       double sfcf,
                       double cfr,
                       double eps_s,
                       double eps_t,
                       double eps_x,
                       int m,
                       int d)
{
    /* --------------------------------------------------------------- */
    /* 1) Unpack RAW states                                            */
    /* --------------------------------------------------------------- */
    const double Swe_raw = u[0];
    const double Sm_raw = u[1];
    const double Uz_raw = u[2];
    const double Lz_raw = u[3];

    /* --------------------------------------------------------------- */
    /* 2) Smooth-clamp storages for flux computations                   */
    /* --------------------------------------------------------------- */
    const double Swe = smooth_pos(Swe_raw, eps_s);
    const double Sm = smooth_pos(Sm_raw, eps_s);
    const double Uz = smooth_pos(Uz_raw, eps_s);
    const double Lz = smooth_pos(Lz_raw, eps_s);

    const double dSwe_dSweRaw = d_smooth_pos_da(Swe_raw, eps_s);
    const double dSm_dSmRaw = d_smooth_pos_da(Sm_raw, eps_s);
    const double dUz_dUzRaw = d_smooth_pos_da(Uz_raw, eps_s);
    const double dLz_dLzRaw = d_smooth_pos_da(Lz_raw, eps_s);

    /* =============================================================== */
    /* 1) Precip partition (smooth snow fraction)                       */
    /* =============================================================== */
    const double uT = (T - T_tr) / eps_t;
    const double snow_fr = 0.5 * (1.0 - tanh(uT));
    const double rain_fr = 1.0 - snow_fr;

    const double dsnow_dT_tr = 0.5 * sech2(uT) / eps_t;

    const double Ps = sfcf * P * snow_fr;
    const double Pr = P * rain_fr;

    const double dPs_dsfcf = P * snow_fr;
    const double dPs_dT_tr = sfcf * P * dsnow_dT_tr;

    const double dPr_dT_tr = -P * dsnow_dT_tr;

    /* =============================================================== */
    /* 2) Melt & refreeze (degree-day, smooth)                          */
    /* =============================================================== */
    const double Tm = T - T_tr;

    const double posTm = smooth_pos(Tm, eps_t);  /* ~ max(Tm,0) */
    const double negTm = smooth_pos(-Tm, eps_t); /* ~ max(-Tm,0) */

    const double dposTm_dTm = d_smooth_pos_da(Tm, eps_t);
    const double dnegTm_dTm = -d_smooth_pos_da(-Tm, eps_t);

    /* potentials */
    const double Mpot = f_dd * posTm;
    const double Rpot = cfr * f_dd * negTm;

    const double dMpot_dT_tr = -f_dd * dposTm_dTm;
    const double dMpot_df_dd = posTm;

    const double dRpot_dT_tr = -(cfr * f_dd) * dnegTm_dTm;
    const double dRpot_df_dd = cfr * negTm;
    const double dRpot_dcfr = f_dd * negTm;

    /* actual melt/refreeze (smooth min) */
    const double M = smooth_min(Swe, Mpot, eps_s);
    const double R = smooth_min(Pr, Rpot, eps_s);

    const double dM_dSwe = d_smooth_min_dA(Swe, Mpot, eps_s);
    const double dM_dMpot = d_smooth_min_dB(Swe, Mpot, eps_s);

    const double dR_dPr = d_smooth_min_dA(Pr, Rpot, eps_s);
    const double dR_dRpot = d_smooth_min_dB(Pr, Rpot, eps_s);

    /* chain to parameters */
    const double dM_dT_tr = dM_dMpot * dMpot_dT_tr;
    const double dM_df_dd = dM_dMpot * dMpot_df_dd;

    /* IMPORTANT: R depends on T_tr both via Pr and via Rpot */
    const double dR_dT_tr = dR_dPr * dPr_dT_tr + dR_dRpot * dRpot_dT_tr;
    const double dR_df_dd = dR_dRpot * dRpot_df_dd;
    const double dR_dcfr = dR_dRpot * dRpot_dcfr;

    /* =============================================================== */
    /* 3) Liquid reaching soil: Win = smooth_pos(Pr + M - R)             */
    /* =============================================================== */
    const double Win_raw = Pr + M - R;
    const double Win = smooth_pos(Win_raw, eps_s);
    const double dWin_dWinraw = d_smooth_pos_da(Win_raw, eps_s);

    const double Pliq = Win;

    /* derivatives wrt Swe (via M only) */
    const double dPliq_dSwe = dWin_dWinraw * dM_dSwe;

    /* wrt T_tr/f_dd/cfr */
    const double dPliq_dT_tr = dWin_dWinraw * (dPr_dT_tr + dM_dT_tr - dR_dT_tr);
    const double dPliq_df_dd = dWin_dWinraw * (dM_df_dd - dR_df_dd);
    const double dPliq_dcfr = dWin_dWinraw * (0.0 - dR_dcfr);

    /* =============================================================== */
    /* 4) Soil evap + recharge                                          */
    /* =============================================================== */

    /* ---- Smooth den = max(lp*f_c, eps_x) with consistent derivatives */
    const double den0 = lp * f_c;
    const double den = smooth_max(den0, eps_x, eps_x);             /* smooth max */
    const double dDen_dDen0 = d_smooth_max_dA(den0, eps_x, eps_x); /* d den / d den0 */

    const double ratio = Sm / den;

    /* phiE = smooth_min(ratio,1) */
    const double phiE = phiE_fun(ratio, eps_s);
    const double dphiE_drat = dphiE_dratio(ratio, eps_s);

    const double E = Ep * phiE;

    /* ratio derivatives (consistent with smooth max) */
    const double dratio_dSm = 1.0 / den;
    const double dratio_dDen = -(Sm / (den * den));
    const double dratio_dDen0 = dratio_dDen * dDen_dDen0;

    const double dratio_dlp = dratio_dDen0 * f_c; /* d(den0)/dlp = f_c */
    const double dratio_dfc = dratio_dDen0 * lp;  /* d(den0)/dfc = lp  */

    const double dE_dSm = Ep * dphiE_drat * dratio_dSm;
    const double dE_dlp = Ep * dphiE_drat * dratio_dlp;
    const double dE_dfc = Ep * dphiE_drat * dratio_dfc;

    /* ---- Smooth fc_safe = max(f_c, eps_x) with consistent derivatives */
    const double fc_safe = smooth_max(f_c, eps_x, eps_x);
    const double dfcSafe_dfc = d_smooth_max_dA(f_c, eps_x, eps_x);

    const double u_raw = Sm / fc_safe;
    const double ucap = clamp01(u_raw, eps_s);
    const double du_duraw = dclamp01_dz(u_raw, eps_s);

    /* u_raw derivatives (now include dfcSafe_dfc) */
    const double duraw_dSm = 1.0 / fc_safe;
    const double duraw_dfc = -(Sm / (fc_safe * fc_safe)) * dfcSafe_dfc;

    /* ucap derivatives */
    const double du_dSm = du_duraw * duraw_dSm;
    const double du_dfc = du_duraw * duraw_dfc;

    /* Recharge: primal and derivatives use the same smooth argument. */
    const double up = smooth_max(ucap, eps_x, eps_x);
    const double dup_du = d_smooth_max_dA(ucap, eps_x, eps_x);
    const double up_beta = pow(up, beta);
    const double Re = Pliq * up_beta;
    const double dRe_dPliq = up_beta;
    const double up_bm1 = pow(up, beta - 1.0);
    const double dRe_dSm = Pliq * beta * up_bm1 * dup_du * du_dSm;
    const double dRe_dfc = Pliq * beta * up_bm1 * dup_du * du_dfc;
    const double dRe_dbeta = Pliq * up_beta * log(up);

    /* =============================================================== */
    /* 5) Response routine (Uz/Lz) smooth                               */
    /* =============================================================== */
    const double h = smooth_pos(Uz - uzl, eps_s);
    const double dh_dUz = d_smooth_pos_da(Uz - uzl, eps_s);
    const double dh_duzl = -dh_dUz;

    const double q0 = k_0 * h;
    const double dq0_dUz = k_0 * dh_dUz;
    const double dq0_dk0 = h;
    const double dq0_duzl = k_0 * dh_duzl;

    const double q1 = k_1 * Uz;
    const double dq1_dUz = k_1;
    const double dq1_dk1 = Uz;

    const double perc_flux = smooth_min(perc, Uz, eps_s);
    const double dperc_dperc = d_smooth_min_dA(perc, Uz, eps_s);
    const double dperc_dUz = d_smooth_min_dB(perc, Uz, eps_s);

    const double q2 = k_2 * Lz;
    const double dq2_dLz = k_2;
    const double dq2_dk2 = Lz;

    const double Q = q0 + q1 + q2;

    /* =============================================================== */
    /* State ODEs (RAW state space)                                     */
    /* =============================================================== */
    udot[0] = Ps + R - M;
    udot[1] = Pliq - Re - E;
    udot[2] = Re - q0 - q1 - perc_flux;
    udot[3] = perc_flux - q2;
    udot[4] = Q;

    /* =============================================================== */
    /* Jacobians for sensitivity ODE: dS/dt = Jx_f*S + Jth_f           */
    /* Build Jx_f in SMOOTHED-state space first, then chain to RAW     */
    /* =============================================================== */
    double Jx_f[5][5] = {{0}};
    double Jth_f[5][12] = {{0}};

    /* --- Jx_f built w.r.t [Swe, Sm, Uz, Lz] smoothed --- */
    /* f1: Ps + R - M */
    Jx_f[0][0] = -dM_dSwe;

    /* f2: Pliq - Re - E */
    Jx_f[1][0] = dPliq_dSwe - dRe_dPliq * dPliq_dSwe;
    Jx_f[1][1] = -dRe_dSm - dE_dSm;

    /* f3: Re - q0 - q1 - perc_flux */
    Jx_f[2][0] = dRe_dPliq * dPliq_dSwe;
    Jx_f[2][1] = dRe_dSm;
    Jx_f[2][2] = -(dq0_dUz + dq1_dUz + dperc_dUz);

    /* f4: perc_flux - q2 */
    Jx_f[3][2] = dperc_dUz;
    Jx_f[3][3] = -dq2_dLz;

    /* f5: Q = q0 + q1 + q2 */
    Jx_f[4][2] = dq0_dUz + dq1_dUz;
    Jx_f[4][3] = dq2_dLz;

    /* --- chain columns back to RAW states --- */
    for (int i = 0; i < m; ++i) {
        Jx_f[i][0] *= dSwe_dSweRaw;
        Jx_f[i][1] *= dSm_dSmRaw;
        Jx_f[i][2] *= dUz_dUzRaw;
        Jx_f[i][3] *= dLz_dLzRaw;
        /* col 4 (Qcum) remains 0 */
    }

    /* Parameter indices (0-based) */
    const int j_fc = 0, j_beta = 1, j_lp = 2, j_k0 = 3, j_uzl_ = 4, j_k1_ = 5, j_k2_ = 6,
              j_perc_ = 7;
    const int j_T_tr = 8, j_f_dd = 9, j_sfcf = 10, j_cfr_ = 11;

    /* ---- f1: Ps + R - M ---- */
    Jth_f[0][j_T_tr] = dPs_dT_tr + dR_dT_tr - dM_dT_tr;
    Jth_f[0][j_f_dd] = 0.0 + dR_df_dd - dM_df_dd;
    Jth_f[0][j_sfcf] = dPs_dsfcf;
    Jth_f[0][j_cfr_] = dR_dcfr;

    /* ---- f2: Pliq - Re - E ---- */
    Jth_f[1][j_T_tr] = dPliq_dT_tr - dRe_dPliq * dPliq_dT_tr;
    Jth_f[1][j_f_dd] = dPliq_df_dd - dRe_dPliq * dPliq_df_dd;
    Jth_f[1][j_cfr_] = dPliq_dcfr - dRe_dPliq * dPliq_dcfr;

    Jth_f[1][j_fc] = -dRe_dfc - dE_dfc;
    Jth_f[1][j_beta] = -dRe_dbeta;
    Jth_f[1][j_lp] = -dE_dlp;

    /* ---- f3: Re - q0 - q1 - perc_flux ---- */
    Jth_f[2][j_T_tr] = dRe_dPliq * dPliq_dT_tr;
    Jth_f[2][j_f_dd] = dRe_dPliq * dPliq_df_dd;
    Jth_f[2][j_cfr_] = dRe_dPliq * dPliq_dcfr;

    Jth_f[2][j_fc] = dRe_dfc;
    Jth_f[2][j_beta] = dRe_dbeta;
    Jth_f[2][j_k0] = -dq0_dk0;
    Jth_f[2][j_uzl_] = -dq0_duzl;
    Jth_f[2][j_k1_] = -dq1_dk1;
    Jth_f[2][j_perc_] = -dperc_dperc;

    /* ---- f4: perc_flux - q2 ---- */
    Jth_f[3][j_k2_] = -dq2_dk2;
    Jth_f[3][j_perc_] = dperc_dperc;

    /* ---- f5: Q ---- */
    Jth_f[4][j_k0] = dq0_dk0;
    Jth_f[4][j_uzl_] = dq0_duzl;
    Jth_f[4][j_k1_] = dq1_dk1;
    Jth_f[4][j_k2_] = dq2_dk2;

    /* --------------------------------------------------------------- */
    /* Sensitivity update: dSdt = Jx_f*Smat + Jth_f                    */
    /* Smat and dSdt are m x d, column-major                           */
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

} // namespace sage_hbv
