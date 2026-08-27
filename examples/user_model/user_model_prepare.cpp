/* ==================================================================
 * USER_MODEL_PREPARE
 *
 * Runtime preparation layer for SAGE user-supplied hydrologic models.
 *
 * This file replaces the MATLAB-side setup formerly done in user_model.m.
 * It converts the generic SAGE MEX call
 *
 *     crr_user_model(par,mdl,data,ode,check)
 *
 * into the solver-core call
 *
 *     sage_user_model::run(tout,z0,dataPrepared,ode)
 *
 * by:
 *   1) transforming parameters from the active SAGE parameter space
 *      to hydrologic space,
 *   2) computing Jth for sensitivity scaling,
 *   3) preparing model-specific constants, routing ordinates, and
 *      auxiliary vectors,
 *   4) defining the physical state dimension m and augmented state
 *      dimension nvar = m*(d+1),
 *   5) constructing the initial augmented state z0 = [u0; vec(S0)],
 *   6) attaching all model-specific fields required by
 *      user_model.cpp to dataPrepared, and
 *   7) defining data.ipr for spin-up-aware q/J extraction.
 *
 * User-edit blocks mark the portions that must be changed when adapting
 * this template to another hydrologic model. The generic SAGE interface,
 * parameter-space transformation, Jth calculation, and augmented-state
 * layout should normally remain unchanged.
 *
 * Written by JA Vrugt, May 2026
 * University of California, Irvine
 * ================================================================== */

#include "mex.h"
#include <cmath>
#include <vector>
#include <string>
#include "user_model.hpp"

/* ---------------------------------------------------------------------
   USER_MODEL_PREPARE
   C++ replacement for the MATLAB-side setup formerly done in user_model.m.

   This file converts the SAGE call

       crr_user_model(par,mdl,data,ode,check)

   into the solver-core call

       sage_user_model::run(tout,Z0,dataPrepared,ode)

   by preparing:
       th, Jth, UH ordinates, dU/dx4, m, nvar, initial Z0,
       constants, parameters, and data.ipr.

   The numerical RK2/model implementation remains in user_model.cpp.
   --------------------------------------------------------------------- */

struct UserModelPrepared {
    double tout = 0.0;
    int d = 0;
    int m = 0;
    int nvar = 0;
    int mem = 1;
    int ipr = 1;

    mxArray* z0 = nullptr;
    mxArray* dataPrepared = nullptr;
    mxArray* Jth = nullptr;
};

static double getScalarField(const mxArray* S, const char* field, const char* where)
{
    const mxArray* A = mxGetField(S, 0, field);
    if (!A || mxGetNumberOfElements(A) != 1 || mxIsComplex(A)) {
        mexErrMsgIdAndTxt("user_model_prepare:badField",
                          "Field '%s' in %s must be a real scalar.",
                          field,
                          where);
    }
    return mxGetScalar(A);
}

static int getIntField(const mxArray* S, const char* field, const char* where)
{
    return static_cast<int>(std::llround(getScalarField(S, field, where)));
}

static const mxArray*
getFieldRequired(const mxArray* S, const char* field, const char* where)
{
    const mxArray* A = mxGetField(S, 0, field);
    if (!A) {
        mexErrMsgIdAndTxt(
            "user_model_prepare:missingField", "Missing field '%s' in %s.", field, where);
    }
    return A;
}

static std::vector<double> vectorFromArray(const mxArray* A, const char* name)
{
    if (!A || mxIsComplex(A) || !(mxIsDouble(A) || mxIsSingle(A))) {
        mexErrMsgIdAndTxt("user_model_prepare:badVector",
                          "%s must be a real single or double vector.",
                          name);
    }

    const mwSize n = mxGetNumberOfElements(A);
    std::vector<double> v((size_t)n);

    if (mxIsDouble(A)) {
        const double* p = mxGetPr(A);
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = p[(size_t)i];
        }
    } else {
        const float* p = static_cast<const float*>(mxGetData(A));
        for (mwSize i = 0; i < n; ++i) {
            v[(size_t)i] = static_cast<double>(p[(size_t)i]);
        }
    }

    return v;
}

static void addScalarField(mxArray* S, const char* name, double value)
{
    mxAddField(S, name);
    mxSetField(S, 0, name, mxCreateDoubleScalar(value));
}

static void addVectorField(mxArray* S, const char* name, const std::vector<double>& v)
{
    mxAddField(S, name);
    mxArray* A = mxCreateDoubleMatrix((mwSize)v.size(), 1, mxREAL);
    double* p = mxGetPr(A);
    for (size_t i = 0; i < v.size(); ++i) {
        p[i] = v[i];
    }
    mxSetField(S, 0, name, A);
}

/* =========================================================
   USER-EDIT BLOCK 1:
   Replace helper functions for a different hydrologic model
   ========================================================= */

/* GR4J cumulative integral function */
static double gr4j_CIF(double t, double x4, double eta)
{
    const double td = t / x4;

    if (t <= 0.0) {
        return 0.0;
    } else if (t < x4) {
        return std::pow(td, eta + 1.0);
    } else if (t < 2.0 * x4) {
        return 1.0 - std::pow(2.0 - td, eta + 1.0);
    } else {
        return 1.0;
    }
}

/* derivative dCIF/dx4 */
static double gr4j_dCIF_dx4(double t, double x4, double eta)
{
    const double p = eta + 1.0;
    const double u = t / x4;

    if (t > 0.0 && t < x4) {
        return -(p / x4) * std::pow(u, p);
    } else if (t >= x4 && t < 2.0 * x4) {
        return -p * std::pow(2.0 - u, p - 1.0) * (t / (x4 * x4));
    } else {
        return 0.0;
    }
}

static void gr4j_UH_ord(double x4,
                        double mts,
                        double eta,
                        int N,
                        std::vector<double>& U,
                        std::vector<double>& dUdx4)
{
    if (N < 1) {
        N = 1;
    }

    std::vector<double> Uraw((size_t)N);
    std::vector<double> dUraw((size_t)N);

    for (int k = 1; k <= N; ++k) {
        const double tk = static_cast<double>(k) * mts;
        const double tkm1 = tk - mts;

        const double Gk = gr4j_CIF(tk, x4, eta);
        const double Gkm1 = gr4j_CIF(tkm1, x4, eta);

        const double dGk = gr4j_dCIF_dx4(tk, x4, eta);
        const double dGkm1 = gr4j_dCIF_dx4(tkm1, x4, eta);

        Uraw[(size_t)(k - 1)] = Gk - Gkm1;
        dUraw[(size_t)(k - 1)] = dGk - dGkm1;
    }

    double maxAbs = 0.0;
    for (double v : Uraw) {
        maxAbs = std::max(maxAbs, std::fabs(v));
    }
    const double eps_U = 1e-12 + 1e-6 * std::max(1.0, maxAbs);

    double s = 0.0;
    double ds = 0.0;

    for (int k = 0; k < N; ++k) {
        const double u = Uraw[(size_t)k];
        const double root = std::sqrt(u * u + eps_U * eps_U);
        const double uPos = 0.5 * (u + root);
        const double dpos = 0.5 * (1.0 + u / root);

        Uraw[(size_t)k] = uPos;
        dUraw[(size_t)k] *= dpos;

        s += Uraw[(size_t)k];
        ds += dUraw[(size_t)k];
    }

    const double sSafe = std::max(s, 1e-30);

    U.assign((size_t)N, 0.0);
    dUdx4.assign((size_t)N, 0.0);

    for (int k = 0; k < N; ++k) {
        U[(size_t)k] = Uraw[(size_t)k] / sSafe;
        dUdx4[(size_t)k] =
            (dUraw[(size_t)k] * sSafe - Uraw[(size_t)k] * ds) / (sSafe * sSafe);
    }
}

/* =========================================================
   END USER-EDIT BLOCK 1:
   ========================================================= */

UserModelPrepared user_model_prepare(const mxArray* par_mx,
                                     const mxArray* mdl,
                                     const mxArray* data,
                                     const mxArray* ode)
{
    UserModelPrepared S;

    const std::vector<double> par = vectorFromArray(par_mx, "par");
    S.d = static_cast<int>(par.size());

    if (S.d < 1) {
        mexErrMsgIdAndTxt("user_model_prepare:emptyPar", "Parameter vector par is empty.");
    }

    const std::vector<double> th_min =
        vectorFromArray(getFieldRequired(mdl, "th_min", "mdl"), "mdl.th_min");
    const std::vector<double> th_max =
        vectorFromArray(getFieldRequired(mdl, "th_max", "mdl"), "mdl.th_max");

    if ((int)th_min.size() != S.d || (int)th_max.size() != S.d) {
        mexErrMsgIdAndTxt(
            "user_model_prepare:badBounds",
            "numel(par), numel(mdl.th_min), and numel(mdl.th_max) must match.");
    }

    const int pspace = getIntField(mdl, "pspace", "mdl");
    S.tout = getScalarField(mdl, "tout", "mdl");
    const std::vector<double> idx =
        vectorFromArray(getFieldRequired(mdl, "idx", "mdl"), "mdl.idx");
    if (idx.size() < 2) {
        mexErrMsgIdAndTxt("user_model_prepare:badIdx",
                          "mdl.idx must have at least two elements.");
    }

    const mxArray* memA = mxGetField(ode, 0, "mem");
    S.mem = memA ? static_cast<int>(std::llround(mxGetScalar(memA))) : 1;
    if (S.mem != 0 && S.mem != 1) {
        mexErrMsgIdAndTxt("user_model_prepare:mem", "ode.mem must be 0 or 1.");
    }
    S.ipr = (S.mem == 0) ? static_cast<int>(std::llround(idx[0])) : 1;

    std::vector<double> th((size_t)S.d, 0.0);
    std::vector<double> Jth((size_t)S.d, 1.0);

    switch (pspace) {
    case 0:
        for (int j = 0; j < S.d; ++j) {
            th[(size_t)j] = par[(size_t)j];
            if (th[(size_t)j] < th_min[(size_t)j] || th[(size_t)j] > th_max[(size_t)j]) {
                mexErrMsgIdAndTxt("user_model_prepare:bounds",
                                  "Hydrologic parameter %d is outside bounds.",
                                  j + 1);
            }
            Jth[(size_t)j] = 1.0;
        }
        break;

    case 1:
        for (int j = 0; j < S.d; ++j) {
            const double nth = par[(size_t)j];
            if (nth < 0.0 || nth > 1.0) {
                mexErrMsgIdAndTxt("user_model_prepare:bounds",
                                  "Normalized parameter %d is outside [0,1].",
                                  j + 1);
            }
            const double dth = th_max[(size_t)j] - th_min[(size_t)j];
            th[(size_t)j] = th_min[(size_t)j] + nth * dth;
            Jth[(size_t)j] = dth;
        }
        break;

    case 2:
        for (int j = 0; j < S.d; ++j) {
            const double varth = par[(size_t)j];
            const double nth = (varth >= 0.0) ? 1.0 / (1.0 + std::exp(-varth))
                                              : std::exp(varth) / (1.0 + std::exp(varth));
            const double dth = th_max[(size_t)j] - th_min[(size_t)j];
            th[(size_t)j] = th_min[(size_t)j] + nth * dth;
            Jth[(size_t)j] = dth * nth * (1.0 - nth);
        }
        break;

    default:
        mexErrMsgIdAndTxt("user_model_prepare:badPspace", "mdl.pspace must be 0, 1, or 2.");
    }

    /* ===========================================
       USER-EDIT BLOCK 2:
       Define model constants and number of states
       =========================================== */
    /* User model constants */
    const double T_sm = 1.0;
    const double eps_m = 1e-6;
    const double eps_s = 1e-12;
    const double kappa = 4.0 / 9.0;
    const double bg = 3.5;
    const double bR = 5.0;
    const double mts = 1.0;
    const double eta = 2.5;
    const double tau = 2.5;
    const double rho = 1e-2;

    if (S.d != 8) {
        mexErrMsgIdAndTxt("user_model_prepare:d",
                          "This GR4JB example user model requires exactly 8 parameters.");
    }

    /* The GR4J unit hydrograph has support [0,2*x4]. The internal state
       dimension may vary between complete model evaluations because SAGE
       only receives discharge and its fixed-size parameter Jacobian. */
    const int N = std::max(1, static_cast<int>(std::ceil((2.0 * th[3]) / mts)));
    std::vector<double> U, dUdx4;
    gr4j_UH_ord(th[3], mts, eta, N, U, dUdx4);

    const int L = static_cast<int>(U.size());
    S.m = 4 + 2 * L;
    S.nvar = S.m * (S.d + 1);
    /* ===========================================
       END USER-EDIT BLOCK 2
       =========================================== */

    /* Initial augmented state */
    S.z0 = mxCreateDoubleMatrix((mwSize)S.nvar, 1, mxREAL);
    double* z0 = mxGetPr(S.z0);

    double y0 = 1e-5;
    const mxArray* y0A = mxGetField(mdl, 0, "y0");
    if (y0A && mxGetNumberOfElements(y0A) >= 1) {
        y0 = mxGetScalar(y0A);
    }

    for (int i = 0; i < S.m; ++i) {
        z0[i] = y0;
    }
    for (int i = S.m; i < S.nvar; ++i) {
        z0[i] = 0.0;
    }

    /* Copy data so we can attach prepared scalar/vector fields */
    S.dataPrepared = mxDuplicateArray(data);

    /* ===============================================================
       USER-EDIT BLOCK 3:
       Define parameter values/input variables used by user_model.cpp
       =============================================================== */
    /* Parameters */
    addScalarField(S.dataPrepared, "x1", th[0]);
    addScalarField(S.dataPrepared, "x2", th[1]);
    addScalarField(S.dataPrepared, "x3", th[2]);
    addScalarField(S.dataPrepared, "x4", th[3]);
    addScalarField(S.dataPrepared, "x5", th[4]);
    addScalarField(S.dataPrepared, "f_p", th[5]);
    addScalarField(S.dataPrepared, "T_tr", th[6]);
    addScalarField(S.dataPrepared, "f_dd", th[7]);

    /* Routing/constant fields */
    addVectorField(S.dataPrepared, "U", U);
    addVectorField(S.dataPrepared, "dUdx4", dUdx4);
    addScalarField(S.dataPrepared, "L", (double)L);
    addScalarField(S.dataPrepared, "eta", eta);
    addScalarField(S.dataPrepared, "tau", tau);
    addScalarField(S.dataPrepared, "kappa", kappa);
    addScalarField(S.dataPrepared, "bg", bg);
    addScalarField(S.dataPrepared, "bR", bR);
    addScalarField(S.dataPrepared, "mts", mts);
    addScalarField(S.dataPrepared, "T_sm", T_sm);
    addScalarField(S.dataPrepared, "eps_m", eps_m);
    addScalarField(S.dataPrepared, "eps_s", eps_s);
    addScalarField(S.dataPrepared, "rho", rho);
    addScalarField(S.dataPrepared, "ipr", (double)S.ipr);
    /* ===============================================================
       END USER-EDIT BLOCK 3
       =============================================================== */

    S.Jth = mxCreateDoubleMatrix((mwSize)S.d, 1, mxREAL);
    double* pJth = mxGetPr(S.Jth);
    for (int j = 0; j < S.d; ++j) {
        pJth[j] = Jth[(size_t)j];
    }

    return S;
}
