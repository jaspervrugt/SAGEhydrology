/* ==================================================================
 * user_model.cpp
 *
 * Example hydrologic model template for user-defined conceptual
 * rainfall–runoff models within the SAGE framework.
 *
 * This file demonstrates how a custom hydrologic model can be coupled
     * to SAGE (Sensitivity-Aware Gradient Estimation) using:
     *   - explicit adaptive-step RK2 (Heun) time integration,
 *   - analytic forward sensitivities,
 *   - augmented-state sensitivity propagation,
 *   - direct MEX-based execution from MATLAB,
 *   - optional state-memory output,
 *   - efficient Jacobian evaluation for gradient-based calibration.
 * The present implementation uses the GR4JB conceptual rainfall–runoff
 * model as a worked example. The file illustrates:
 *   1) definition of physical state equations,
 *   2) routing and unit-hydrograph implementation,
 *   3) construction of analytic Jacobians,
 *   4) propagation of parameter sensitivities,
 *   5) extraction of discharge and dQ/dtheta,
 *   6) integration with SAGE parameter-space transformations.
 *
 * The model-specific portions of this file are clearly marked using
 *      USER-EDIT BLOCK
 * comments. Users adapting this template to another hydrologic model
 * typically only need to modify:
 *   - model parameters,
 *   - state equations,
 *   - routing structure,
 *   - Jacobian definitions,
 *   - discharge extraction logic.
 *
 * The generic SAGE/MEX interface, adaptive RK2 integration, augmented
 * sensitivity-state layout, and memory/output handling can normally
 * remain unchanged.
 *
 * Augmented-state layout:
 *      z = [u ; vec(S)]
 * where
 *      u = physical model states
 *      S = du/dtheta sensitivity matrix
 * with dimensions
 *      m = number of physical states
 *      d = number of parameters
 *      nvar = m*(d+1)
 *
 * Written by JA Vrugt, May 2026
 * University of California, Irvine
 * ================================================================== */

#include "mex.h"
#include "user_model.hpp"
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "log_fail_params.hpp" // to log bad parameter values

static void runge_kutta(int nvar,
                        int d,
                        int nt,
                        const double* z0,
                        const mxArray* data,
                        const mxArray* options,
                        double* Z,
                        int mem,
                        int ipr,
                        int ns,
                        double* q_out,
                        double* J_out,
                        int n_q,
                        bool* fail_out);

/* ================================================================
   USER-EDIT BLOCK 0:
   Update calls with those of your model
   RK2: update "double x1, ... , double rho" to what _odefcn uses
   aug_ode: update "double x1, ... , double rho" but end with
            "int nvar, int m, int d"
   _odefcn: update "double x1, ... , double rho" but end with
            "int nvar, int m, int d"
   ================================================================ */
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
                double rho);

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

/* ================================================================
   END USER-EDIT BLOCK 0
   ================================================================ */

static inline bool isFinite(double x)
{
    return mxIsFinite(x);
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

/* =====================================================================
   GR4JB mex wrapper using scalar final print time t_last:
   - accepts z0 as single or double
   - options.mem (default 1): mem==1 stores full trajectory (nt x nvar);
     mem==0 stores only final (1 x nvar)
   - optional 4th output: fail flag

   Internal print times are assumed to be:
       0, 1, 2, ..., t_last
   so
       ns = t_last
       nt = ns + 1
   ===================================================================== */

static const double*
getFieldAsDoublePtr(const mxArray* S, const char* fname, std::vector<double>& buf)
{
    const mxArray* A = mxGetField(S, 0, fname);
    if (!A) {
        mexErrMsgIdAndTxt("gr4jB_mex:data", "Missing field '%s' in data struct.", fname);
    }
    if (mxIsComplex(A)) {
        mexErrMsgIdAndTxt("gr4jB_mex:type", "Field '%s' must be real.", fname);
    }
    if (mxIsDouble(A)) {
        return mxGetPr(A);
    }
    if (mxIsSingle(A)) {
        const mwSize n = mxGetNumberOfElements(A);
        const float* p = (const float*)mxGetData(A);
        buf.resize((size_t)n);
        for (mwSize i = 0; i < n; i++) {
            buf[(size_t)i] = (double)p[i];
        }
        return buf.data();
    }

    mexErrMsgIdAndTxt("gr4jB_mex:type", "Field '%s' must be single or double.", fname);

    return nullptr; // never reached
}

static const mxArray* getOptionalField(const mxArray* S, const char* fname)
{
    return mxGetField(S, 0, fname); // nullptr if missing
}

inline const mxArray*
getRequiredField(const mxArray* s, const char* field, const char* parent)
{
    if (!s) {
        mexErrMsgIdAndTxt("mex:missingStruct", "Null struct pointer for '%s'.", parent);
    }

    const mxArray* out = mxGetField(s, 0, field);
    if (!out) {
        mexErrMsgIdAndTxt(
            "mex:missingField", "Missing required field '%s' in '%s'.", field, parent);
    }
    return out;
}

static double getRequiredScalar(const mxArray* S, const char* fname, const char* where)
{
    const mxArray* A = mxGetField(S, 0, fname);
    if (!A) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:missingField", "Missing field '%s' in %s.", fname, where);
    }
    if (mxGetNumberOfElements(A) != 1 || mxIsComplex(A)) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:badField", "Field '%s' in %s must be a real scalar.", fname, where);
    }
    if (mxIsDouble(A) || mxIsSingle(A)) {
        return mxGetScalar(A);
    }

    mexErrMsgIdAndTxt("crr_gr4jB:badField",
                      "Field '%s' in %s must be a real scalar single or double.",
                      fname,
                      where);

    return mxGetNaN(); // never reached
}

static int getRequiredInt(const mxArray* S, const char* fname, const char* where)
{
    const double v = getRequiredScalar(S, fname, where);
    if (std::floor(v) != v) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:badField", "Field '%s' in %s must be an integer.", fname, where);
    }
    return (int)v;
}

static const double* getRequiredPr(const mxArray* S, const char* fname, const char* where)
{
    const mxArray* A = mxGetField(S, 0, fname);
    if (!A) {
        mexErrMsgIdAndTxt(
            "crr_gr4jB:missingField", "Missing field '%s' in %s.", fname, where);
    }
    if (mxIsComplex(A) || !mxIsDouble(A)) {
        mexErrMsgIdAndTxt("crr_gr4jB:badField",
                          "Field '%s' in %s must be a real double array.",
                          fname,
                          where);
    }
    return mxGetPr(A);
}

/* ------------------------------------------------------------------ */
/* gr4jB_mex: Updated                                                 */
/*  - input 1 is scalar t_last                                        */
/*  - ipr is 1-based index in [1..ns+1]                               */
/*  - Before ipr: spin-up, do NOT store q/J                           */
/*  - From ipr..ns: store q and J (if mem==0)                         */
/*  - mem==1 unchanged (store full Z, return q=[], J=[])              */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* sage_user_model::run                                               */
/*                                                                    */
/* Core solver entry point used by the modular user_model MEX design. */
/* This function receives the already prepared inputs directly from   */
/* crr_user_model_mex.cpp / user_model_prepare.cpp:                   */
/*                                                                    */
/*     t_last  scalar final print time                                */
/*     z0_mx   initial augmented state vector                         */
/*     data    prepared data struct parameters, UH ordinates, etc.    */
/*     options ODE/settings struct                                    */
/*                                                                    */
/* Outputs follow the crr_user_model convention:                      */
/*     plhs[0] = Z                                                    */
/*     plhs[1] = q_n        if requested and mem == 0                 */
/*     plhs[2] = J_raw      if requested and mem == 0                 */
/*     plhs[3] = fail flag  if requested                              */
/* ------------------------------------------------------------------ */
void sage_user_model::run(int nlhs,
                          mxArray* plhs[],
                          double t_last_in,
                          const mxArray* z0_mx,
                          const mxArray* data,
                          const mxArray* options)
{
    if (nlhs < 1 || nlhs > 4) {
        mexErrMsgIdAndTxt("crr_user_model:nlhs",
                          "Expected 1 to 4 outputs: Z, q_n, J_raw, fail.");
    }

    if (!z0_mx || !data || !options) {
        mexErrMsgIdAndTxt("crr_user_model:nullInput",
                          "z0_mx, data, and options must be non-null mxArray pointers.");
    }

    if (!mxIsFinite(t_last_in) || t_last_in < 1.0) {
        mexErrMsgIdAndTxt("crr_user_model:tlast",
                          "Final print time must be finite and >= 1.");
    }

    const int ns = (int)std::llround(t_last_in);
    if (fabs(t_last_in - (double)ns) > 1e-12) {
        mexErrMsgIdAndTxt("crr_user_model:tlast", "Final print time must be an integer.");
    }
    const int nt = ns + 1;

    /* ================================================================
       USER-EDIT BLOCK 1:
       Replace with number of parameters of your model.
       ================================================================ */
    const int d = 8; /* number of parameters for this GR4JB user model */
    /* ================================================================
       END USER-EDIT BLOCK 1
       ================================================================ */

    const mwSize nvar_mw = mxGetNumberOfElements(z0_mx);
    if (nvar_mw < 1) {
        mexErrMsgIdAndTxt("crr_user_model:badSizes", "z0 must have >= 1 element.");
    }

    /* --- z0 pointer: accept single or double --- */
    std::vector<double> z0buf;
    const double* z0 = nullptr;
    if (mxIsDouble(z0_mx)) {
        z0 = mxGetPr(z0_mx);
    } else if (mxIsSingle(z0_mx)) {
        const mwSize n = mxGetNumberOfElements(z0_mx);
        const float* p = (const float*)mxGetData(z0_mx);
        z0buf.resize((size_t)n);
        for (mwSize i = 0; i < n; ++i) {
            z0buf[(size_t)i] = (double)p[i];
        }
        z0 = z0buf.data();
    } else {
        mexErrMsgIdAndTxt("crr_user_model:type", "z0 must be single or double.");
    }

    /* --- memory flag (default 1) --- */
    int mem = 1;
    {
        const mxArray* A = getOptionalField(options, "mem");
        if (A) {
            if (mxGetNumberOfElements(A) != 1 || mxIsComplex(A)) {
                mexErrMsgIdAndTxt("crr_user_model:mem",
                                  "options.mem must be a real scalar.");
            }
            mem = (int)mxGetScalar(A);
        }
    }
    if (mem != 0 && mem != 1) {
        mexErrMsgIdAndTxt("crr_user_model:mem", "options.mem must be 0 or 1.");
    }

    mwSize Z_rows = (mem == 1) ? (mwSize)nt : (mwSize)1;
    mwSize Z_cols = (mwSize)nvar_mw;

    /* Output #1: Z */
    plhs[0] = mxCreateDoubleMatrix(Z_rows, Z_cols, mxREAL);
    double* Z = mxGetPr(plhs[0]);

    /* --- ipr: 1-based start index for printing q/J --- */
    int ipr = 1;
    {
        const mxArray* A = getOptionalField(data, "ipr");
        if (A) {
            if (mxGetNumberOfElements(A) != 1) {
                mexErrMsgIdAndTxt("crr_user_model:iprNotScalar",
                                  "data.ipr must be a scalar (1-based).");
            }
            ipr = (int)mxGetScalar(A);
        }
    }
    if (ipr < 1) {
        ipr = 1;
    }
    if (ipr > ns + 1) {
        ipr = ns + 1;
    }

    /* Outputs #2 and #3: q and raw J, only when mem == 0 */
    double* q_out = nullptr;
    double* J_out = nullptr;

    int n_q = 0;
    if (mem == 0) {
        n_q = (ipr <= ns) ? (ns - ipr + 1) : 0;
    }

    if (nlhs >= 2) {
        if (mem == 0 && n_q > 0) {
            plhs[1] = mxCreateDoubleMatrix((mwSize)n_q, 1, mxREAL);
            q_out = mxGetPr(plhs[1]);

            if (nlhs >= 3) {
                plhs[2] = mxCreateDoubleMatrix((mwSize)n_q, (mwSize)d, mxREAL);
                J_out = mxGetPr(plhs[2]);
            }
        } else {
            plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
            if (nlhs >= 3) {
                plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
            }
        }
    }

    bool fail = false;

    runge_kutta(
        (int)nvar_mw, d, nt, z0, data, options, Z, mem, ipr, ns, q_out, J_out, n_q, &fail);

    if (nlhs >= 4) {
        plhs[3] = mxCreateLogicalScalar(fail);
    }
}

/* ------------------------------------------------------------------ */
/* Adaptive RK2 (Heun) integrator                                     */
/*  - updated: internal print times are 0,1,2,...,ns                  */
/* ------------------------------------------------------------------ */
static void runge_kutta(int nvar,
                        int d,
                        int nt,
                        const double* z0,
                        const mxArray* data,
                        const mxArray* options,
                        double* Z,
                        int mem,
                        int ipr,
                        int ns,
                        double* q_out,
                        double* J_out,
                        int n_q,
                        bool* fail_out)
{
    // ---- Integration options ----
    const double hin = getRequiredScalar(options, "InitStep", "options");
    const double hmax_ = getRequiredScalar(options, "MaxStep", "options");
    const double hmin_ = getRequiredScalar(options, "MinStep", "options");
    const double reltol = getRequiredScalar(options, "RelTol", "options");
    const double abstol = getRequiredScalar(options, "AbsTol", "options");
    const double order = getRequiredScalar(options, "Order", "options");
    const double maxiter = getRequiredScalar(options, "maxiter", "options");

    if (!mxIsFinite(hin) || !mxIsFinite(hmax_) || !mxIsFinite(hmin_) ||
        !mxIsFinite(reltol) || !mxIsFinite(abstol) || !mxIsFinite(order) ||
        !mxIsFinite(maxiter) || hin <= 0.0 || hmax_ <= 0.0 || hmin_ <= 0.0 ||
        hmin_ > hmax_ || reltol < 0.0 || abstol <= 0.0 || order <= 0.0 || maxiter < 1.0) {
        mexErrMsgIdAndTxt("crr_user_model:solverOptions",
                          "Invalid adaptive solver options.");
    }

    // ---- Forcings ----
    std::vector<double> Pbuf, Epbuf, Tbuf;
    const mxArray* Pmx = mxGetField(data, 0, "P");
    const mxArray* Epmx = mxGetField(data, 0, "Ep");
    const mxArray* Tmx = mxGetField(data, 0, "T");

    const double* P = getFieldAsDoublePtr(data, "P", Pbuf);
    const double* Ep = getFieldAsDoublePtr(data, "Ep", Epbuf);
    const double* T = getFieldAsDoublePtr(data, "T", Tbuf);

    if (mxGetNumberOfElements(Pmx) < (mwSize)ns ||
        mxGetNumberOfElements(Epmx) < (mwSize)ns ||
        mxGetNumberOfElements(Tmx) < (mwSize)ns) {
        mexErrMsgIdAndTxt(
            "gr4jB_mex:dataSize", "P/Ep/T must have at least ns=%d elements.", ns);
    }

    /* ================================================================
       USER-EDIT BLOCK 2:
       Replace with parameters & other inputs of your model.
       ================================================================ */
    // ---- Parameters (GR4JB-specific) ----
    const double x1 = getRequiredScalar(data, "x1", "data");
    const double x2 = getRequiredScalar(data, "x2", "data");
    const double x3 = getRequiredScalar(data, "x3", "data");
    const double x4 = getRequiredScalar(data, "x4", "data");
    const double x5 = getRequiredScalar(data, "x5", "data");
    const double f_p = getRequiredScalar(data, "f_p", "data");
    const double T_tr = getRequiredScalar(data, "T_tr", "data");
    const double f_dd = getRequiredScalar(data, "f_dd", "data");

    const double* U = getRequiredPr(data, "U", "data");
    const double* dUdx4 = getRequiredPr(data, "dUdx4", "data");
    const int L = getRequiredInt(data, "L", "data");
    if (L < 1) {
        mexErrMsgIdAndTxt("your_mex:data", "L must be a positive integer.");
    }
    const double eta = getRequiredScalar(data, "eta", "data");
    const double tau = getRequiredScalar(data, "tau", "data");
    const double kappa = getRequiredScalar(data, "kappa", "data");
    const double bg = getRequiredScalar(data, "bg", "data");
    const double bR = getRequiredScalar(data, "bR", "data");
    const double mts = getRequiredScalar(data, "mts", "data");

    const double T_sm = getRequiredScalar(data, "T_sm", "data");
    const double eps_m = getRequiredScalar(data, "eps_m", "data");
    const double eps_s = getRequiredScalar(data, "eps_s", "data");
    const double rho = getRequiredScalar(data, "rho", "data");

    const int m = 4 + 2 * L; // physical states
    /* ================================================================
       END USER-EDIT BLOCK 2
       ================================================================ */

    // ---- Layout validation ----
    const int expected_nvar = m + m * d;
    if (nvar != expected_nvar) {
        mexErrMsgIdAndTxt("gr4jB_mex:nvarMismatch",
                          "nvar=%d, expected %d = m + m*d (m=%d, d=%d). Check z0 layout.",
                          nvar,
                          expected_nvar,
                          m,
                          d);
    }

    // ---- Allocate temporaries ----
    double* LTE = (double*)mxMalloc((size_t)nvar * sizeof(double));
    double* ztmp = (double*)mxMalloc((size_t)nvar * sizeof(double));
    double* w = (double*)mxMalloc((size_t)nvar * sizeof(double));
    double* zdotE = (double*)mxMalloc((size_t)nvar * sizeof(double));
    double* zE = (double*)mxMalloc((size_t)nvar * sizeof(double));
    double* zdot = (double*)mxMalloc((size_t)nvar * sizeof(double));

    // current augmented state
    double* zcur = (double*)mxMalloc((size_t)nvar * sizeof(double));
    for (int i = 0; i < nvar; ++i) {
        zcur[i] = z0[i];
    }

    // store initial
    int last_stored_row = 0;
    if (mem == 1) {
        for (int j = 0; j < nvar; ++j) {
            Z[(size_t)0 + (size_t)nt * (size_t)j] = zcur[j];
        }
    } else {
        // Keep output #1 valid even when the first interval fails.
        for (int j = 0; j < nvar; ++j) {
            Z[(size_t)j] = zcur[j];
        }
    }

    // streaming previous (for q/J diffs)
    double* prev = nullptr;
    if (mem == 0) {
        prev = (double*)mxMalloc((size_t)nvar * sizeof(double));
        for (int i = 0; i < nvar; ++i) {
            prev[i] = z0[i];
        }
    }

    bool fail = false;
    int flag = 0;
    int iterCount = 0; // consecutive rejections
    int k_out = 0;

    /*
     * Carry the adaptive controller's recommended next step between
     * forcing intervals. Keep this value separate from the step h that
     * may be shortened only to land exactly on the next forcing boundary.
     */
    double hCarry = fmax(hmin_, fmin(hin, hmax_));

    for (int s = 1; s <= ns; ++s) {
        const double t1 = (double)(s - 1);
        const double t2 = (double)s;

        double tcur = t1;

        double h = fmin(hCarry, t2 - t1);

        while (tcur < t2) {
            for (int i = 0; i < nvar; ++i) {
                ztmp[i] = zcur[i];
            }

            /* =========================================================
               USER-EDIT BLOCK 3:
               Replace with rk2 call of your model: do not touch
               nvar, m, d, h, ztmp, LTE, zdotE, zE, zdot,
               P[s-1], Ep[s-1], T[s-1], but change x1, x2, ... , rho
               ========================================================= */
            rk2(nvar,
                m,
                d,
                h,
                ztmp,
                LTE,
                zdotE,
                zE,
                zdot,
                P[s - 1],
                Ep[s - 1],
                T[s - 1],
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
            /* =========================================================
               END USER-EDIT BLOCK 3
               ========================================================= */

            // sanity check
            for (int i = 0; i < nvar; ++i) {
                if (!isFinite(ztmp[i])) {
                    mexPrintf("FAIL nonfinite: s=%d t=%.10g h=%.3g i=%d z=%.3g\n",
                              s,
                              tcur,
                              h,
                              i,
                              ztmp[i]);
                    fail = true;
                    break;
                }
                /* Sensitivity states can legitimately exceed the physical
                   state scale over long records. Match the built-in model
                   safeguard by applying the magnitude limit only to the
                   m physical/cumulative states. Nonfinite values remain
                   fatal for every augmented state above. */
                if (fabs(ztmp[i]) > 1e12) {
                    mexPrintf("FAIL blowup:   s=%d t=%.10g h=%.3g i=%d z=%.3g\n",
                              s,
                              tcur,
                              h,
                              i,
                              ztmp[i]);
                    fail = true;
                    break;
                }
            }
            /* =========================================================
               USER-EDIT BLOCK 4:
               Replace with parameters and "name" of your model
               ========================================================= */
            if (fail) {
                mexPrintf("FAIL parameters: x1=%.4f x2=%.4f x3=%.4f x4=%.4f x5=%.4f "
                          "f_p=%.4f T_tr=%.4f f_dd=%.4f\n",
                          x1,
                          x2,
                          x3,
                          x4,
                          x5,
                          f_p,
                          T_tr,
                          f_dd);
                log_fail_params("gr4jB", x1, x2, x3, x4, x5, f_p, T_tr, f_dd);
                break;
            }
            /* =========================================================
               END USER-EDIT BLOCK 4
               ========================================================= */

            double sumWL = 0.0;
            for (int i = 0; i < nvar; ++i) {
                w[i] = 1.0 / (reltol * fabs(ztmp[i]) + abstol);
                const double wl = w[i] * LTE[i];
                sumWL += wl * wl;
            }
            const double wrms = sqrt(sumWL / (double)nvar);

            const int accepted = ((wrms <= 1.0) || (h <= hmin_)) ? 1 : 0;

            if (accepted) {
                for (int i = 0; i < nvar; ++i) {
                    zcur[i] = ztmp[i];
                }
                tcur += h;
                iterCount = 0; // reset on progress
            } else {
                iterCount += 1;
            }

            // new step size
            double fac;
            if (wrms <= 0.0) {
                fac = 5.0;
            } else {
                fac = 0.9 * pow(wrms, -1.0 / order);
                fac = fmin(5.0, fmax(0.2, fac));
            }
            /*
             * Controller recommendation for the next attempted step.
             * Save the unclipped recommendation in hCarry. The current
             * step h may still be shortened below only to land at t2.
             */
            double hNext = h * fac;
            hNext = fmax(hmin_, fmin(hNext, hmax_));
            hCarry = hNext;

            // robust end-of-interval
            const double tleft = t2 - tcur;
            const double tTol = 10.0 * mxGetEps() * fmax(1.0, fmax(fabs(t2), fabs(tcur)));

            if (tleft <= tTol) {
                tcur = t2;
                break;
            }

            // Clip only the current attempted step to the forcing boundary.
            h = fmin(hNext, tleft);
            if (h < tTol) {
                mexPrintf("FAIL: step underflow at t=%.10g, tleft=%.3g\n", tcur, tleft);
                fail = true;
                break;
            }

            /* =========================================================
               USER-EDIT BLOCK 5:
               Replace with parameters of your model
               ========================================================= */
            // maxiter check (consecutive rejections)
            if (iterCount >= (int)maxiter) {
                if (flag == 0) {
                    mexPrintf("WARNING: Max rejection limit reached at t = %.5f\n", tcur);
                    mexPrintf("WARNING: Parameter x1  %8.5f\n", x1);
                    mexPrintf("WARNING: Parameter x2  %8.5f\n", x2);
                    mexPrintf("WARNING: Parameter x3  %8.5f\n", x3);
                    mexPrintf("WARNING: Parameter x4  %8.5f\n", x4);
                    mexPrintf("WARNING: Parameter x5  %8.5f\n", x5);
                    mexPrintf("WARNING: Parameter T_tr %8.5f\n", T_tr);
                    mexPrintf("WARNING: Parameter f_dd %8.5f\n", f_dd);
                    flag = 1;
                }
                mexPrintf("FAIL: maxiter (rejections) at t=%.10g h=%.3g wrms=%.3g iter=%d "
                          "maxiter=%d\n",
                          tcur,
                          h,
                          wrms,
                          iterCount,
                          (int)maxiter);
                fail = true;
                break;
            }
            /* =========================================================
               END USER-EDIT BLOCK 5
               ========================================================= */

        } // while

        if (fail) {
            break;
        }

        // store/stream at print time
        if (mem == 1) {
            for (int j = 0; j < nvar; ++j) {
                Z[(size_t)s + (size_t)nt * (size_t)j] = zcur[j];
            }
            last_stored_row = s;
        } else {
            for (int j = 0; j < nvar; ++j) {
                Z[(size_t)0 + (size_t)1 * (size_t)j] = zcur[j];
            }

            // compute q/J only if s >= ipr
            if (q_out && (s >= ipr)) {
                if (k_out >= n_q) {
                    fail = true;
                    break;
                }

                // by convention: cumulative discharge stored in state m-1
                q_out[k_out] = zcur[m - 1] - prev[m - 1];
                // this is last state variable (= m) but m-1 in C++

                if (J_out) {
                    for (int j = 0; j < d; ++j) {
                        const int idx = (j + 2) * m - 1;
                        J_out[(size_t)k_out + (size_t)n_q * (size_t)j] =
                            zcur[idx] - prev[idx];
                    }
                }
                k_out++;
            }

            for (int i = 0; i < nvar; ++i) {
                prev[i] = zcur[i];
            }
        }
    } // for s

    // Match MATLAB mem==1 failure behavior: repeat the last valid state.
    if (fail && mem == 1) {
        for (int row = last_stored_row + 1; row < nt; ++row) {
            for (int j = 0; j < nvar; ++j) {
                Z[(size_t)row + (size_t)nt * (size_t)j] =
                    Z[(size_t)last_stored_row + (size_t)nt * (size_t)j];
            }
        }
    }

    if (fail_out) {
        *fail_out = fail;
    }

    if (prev) {
        mxFree(prev);
    }
    mxFree(zcur);

    mxFree(LTE);
    mxFree(ztmp);
    mxFree(w);
    mxFree(zdotE);
    mxFree(zE);
    mxFree(zdot);
}

/* =========================================================
   USER-EDIT BLOCK 6:
   Adjust RK2/aug_ode/_odefcn helper functions to match your
   model. For RK2, change block: double x1, ... , double rho
   For _aug_ode change x1, ..., rho  but leave  nvar, m, d
   ========================================================= */
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

/* =========================================================
   END USER-EDIT BLOCK 6
   ========================================================= */

/* =========================================================
   USER-EDIT BLOCK 7:
   Match aug_ode/_odefcn functions to calls in BLOCK 6
   Replace the hydrologic model equations and analytic sensitivities.
   In _odefcn you must fill:
     udot[0:m-1]
     dSdt[0:m*d-1]
   and the final sensitivity propagation pattern
     dSdt = Jx_f*Smat + Jth_f
   can remain unchanged if you redefine Jx_f and Jth_f.
   ========================================================= */
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
/* gr4jB ODE + sensitivities  (MATLAB-consistent version)             */
/* u = [Swe, Sp, Rr, U1(1..n1), U2(1..n2), Qinf]                      */
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
        /* if you don't want mexErrMsgTxt here, replace by return; */
        mexErrMsgTxt("gr4jB_odefcn: m must equal 4 + 2*L for layout [Swe Sp Rr z1 z2 Q].");
    }
    if (d != 8) {
        mexErrMsgTxt(
            "gr4jB_odefcn: d must equal 8 for theta=[x1 x2 x3 x4 x5 f_p T_tr f_dd].");
    }

    /* pointers to UH states */
    double* z1 = &u[iz1];
    double* z2 = &u[iz2];

    /* ----------------------------------------------------------------- */
    /* 4) Smooth positivity on storages (states only)                    */
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

    /* Nonzero state derivatives used directly in the sparse sensitivity
       equations below. The derivatives with respect to the three
       smoothed physical states include their smooth_pos chain rules. */
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
    /* Smat is m x d column-major. Jx_f is never materialized: each      */
    /* equation below evaluates only its actual nonzero dependencies.    */
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

/* =========================================================
   END USER-EDIT BLOCK 7
   ========================================================= */
