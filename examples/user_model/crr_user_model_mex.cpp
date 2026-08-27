/* ==================================================================
 * crr_user_model_mex.cpp
 *
 * Fixed MEX gateway for user-supplied hydrologic models in SAGE.
 *
 * This file provides the single MATLAB/MEX entry point
 *
 *     [q_n,J,Jth,Z,fail] = crr_user_model(par,mdl,data,ode,check)
 *
 * used by the SAGE GUI and pipeline when model = user_model and
 * mcode = 4. It should normally not be edited by users.
 *
 * The gateway connects the generic SAGE calling convention to the
 * modular user-model implementation by:
 *   1) receiving parameters, model metadata, forcing data, and ODE
 *      settings from SAGE,
 *   2) calling user_model_prepare.cpp to transform parameters,
 *      initialize states, construct Jth, and attach model-specific
 *      fields to dataPrepared,
 *   3) calling user_model.cpp to run the hydrologic model
 *      and compute discharge and raw sensitivities,
 *   4) scaling raw sensitivities with Jth so J is returned in the
 *      active SAGE parameter space, and
 *   5) returning q_n, J, Jth, optional state output Z, and fail using the
 *      same convention as built-in SAGE model MEX files.
 *
 * Users implementing a new hydrologic model should normally modify:
 *     make_user_model_info.m
 *     user_model_prepare.cpp
 *     user_model.cpp
 * and leave this gateway unchanged.
 *
 * Written by JA Vrugt, May 2026
 * University of California, Irvine
 * ================================================================== */

#include "mex.h"
#include <cstring>
#include <cmath>
#include "user_model_prepare.hpp"
#include "user_model.hpp"

/* ---------------------------------------------------------------------
   USER_MODEL_GATEWAY
   Single MEX entry point for deployed user_model use.

   MATLAB call:
       [q_n,J,Jth,Z,fail] = crr_user_model(par,mdl,data,ode,check)

   Internally:
       1. user_model_prepare converts par/mdl/data/ode into the old core
          inputs: tout, z0, prepared data, ode
       2. sage_user_model::run executes the RK2 hydrologic model
       3. J is scaled by Jth so sensitivities are returned in the active
          parameter space.
   --------------------------------------------------------------------- */

static mxArray* scaleSensitivity(const mxArray* Jraw, const mxArray* Jth)
{
    if (!Jraw || mxIsEmpty(Jraw)) {
        return Jraw ? mxDuplicateArray(Jraw) : mxCreateDoubleMatrix(0, 0, mxREAL);
    }

    mxArray* J = mxDuplicateArray(Jraw);
    double* pJ = mxGetPr(J);
    const double* pJth = mxGetPr(Jth);
    const mwSize n = mxGetM(J);
    const mwSize d = mxGetN(J);

    if (mxGetNumberOfElements(Jth) != d) {
        mexErrMsgIdAndTxt("crr_user_model:badJth", "numel(Jth) must match size(J,2).");
    }

    for (mwSize j = 0; j < d; ++j) {
        for (mwSize i = 0; i < n; ++i) {
            pJ[i + n * j] *= pJth[j];
        }
    }
    return J;
}

static void getOutputRange(const mxArray* mdl, mwIndex* i0, mwIndex* i1)
{
    const mxArray* idxA = mxGetField(mdl, 0, "idx");
    if (!idxA || mxGetNumberOfElements(idxA) < 2 || mxIsComplex(idxA)) {
        mexErrMsgIdAndTxt("crr_user_model:badIdx",
                          "mdl.idx must contain start and end indices.");
    }
    const double* idx = mxGetPr(idxA);
    const double a = idx[0];
    const double b = idx[1];
    if (!mxIsFinite(a) || !mxIsFinite(b) || a < 1.0 || b <= a || std::floor(a) != a ||
        std::floor(b) != b) {
        mexErrMsgIdAndTxt("crr_user_model:badIdx",
                          "mdl.idx must contain valid integer MATLAB indices [start,end].");
    }
    *i0 = (mwIndex)a - 1;
    *i1 = (mwIndex)b - 1;
}

static mxArray* extractDischargeFromZ(const mxArray* Z, const mxArray* mdl, int m)
{
    mwIndex i0, i1;
    getOutputRange(mdl, &i0, &i1);
    const mwSize nt = mxGetM(Z);
    if (i1 >= nt) {
        mexErrMsgIdAndTxt("crr_user_model:badIdx",
                          "mdl.idx exceeds the stored state trajectory.");
    }
    const mwSize n = (mwSize)(i1 - i0);
    mxArray* q = mxCreateDoubleMatrix(n, 1, mxREAL);
    double* pq = mxGetPr(q);
    const double* pZ = mxGetPr(Z);
    const mwIndex qcol = (mwIndex)(m - 1);
    for (mwIndex k = 0; k < n; ++k) {
        pq[k] = pZ[(i0 + k + 1) + nt * qcol] - pZ[(i0 + k) + nt * qcol];
    }
    return q;
}

static mxArray* extractJacobianFromZ(const mxArray* Z, const mxArray* mdl, int m, int d)
{
    mwIndex i0, i1;
    getOutputRange(mdl, &i0, &i1);
    const mwSize nt = mxGetM(Z);
    if (i1 >= nt) {
        mexErrMsgIdAndTxt("crr_user_model:badIdx",
                          "mdl.idx exceeds the stored state trajectory.");
    }
    const mwSize n = (mwSize)(i1 - i0);
    mxArray* J = mxCreateDoubleMatrix(n, (mwSize)d, mxREAL);
    double* pJ = mxGetPr(J);
    const double* pZ = mxGetPr(Z);
    for (int j = 0; j < d; ++j) {
        const mwIndex col = (mwIndex)((j + 2) * m - 1);
        for (mwIndex k = 0; k < n; ++k) {
            pJ[k + n * (mwIndex)j] = pZ[(i0 + k + 1) + nt * col] - pZ[(i0 + k) + nt * col];
        }
    }
    return J;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 4 || nrhs > 5) {
        mexErrMsgIdAndTxt("crr_user_model:nrhs",
                          "Expected crr_user_model(par,mdl,data,ode,check).");
    }
    if (nlhs < 1 || nlhs > 5) {
        mexErrMsgIdAndTxt("crr_user_model:nlhs",
                          "Expected 1 to 5 outputs: q_n, J, Jth, Z, fail.");
    }

    const mxArray* par = prhs[0];
    const mxArray* mdl = prhs[1];
    const mxArray* data = prhs[2];
    const mxArray* ode = prhs[3];

    UserModelPrepared S = user_model_prepare(par, mdl, data, ode);
    mxArray* coreOut[4] = {nullptr, nullptr, nullptr, nullptr};
    mxArray* Zraw = nullptr;
    mxArray* qraw = nullptr;
    mxArray* Jraw = nullptr;

    if (S.mem == 0) {
        const int coreNlhs = 4;
        sage_user_model::run(coreNlhs, coreOut, S.tout, S.z0, S.dataPrepared, ode);
        Zraw = coreOut[0];
        qraw = coreOut[1];
        Jraw = (coreNlhs >= 3) ? coreOut[2] : nullptr;
    } else {
        sage_user_model::run(4, coreOut, S.tout, S.z0, S.dataPrepared, ode);
        Zraw = coreOut[0];
        qraw = extractDischargeFromZ(Zraw, mdl, S.m);
        if (nlhs >= 2) {
            Jraw = extractJacobianFromZ(Zraw, mdl, S.m, S.d);
        }
    }

    plhs[0] = qraw ? qraw : mxCreateDoubleMatrix(0, 0, mxREAL);

    if (nlhs >= 2) {
        plhs[1] = scaleSensitivity(Jraw, S.Jth);
    }
    if (Jraw) {
        mxDestroyArray(Jraw);
    }

    if (nlhs >= 3) {
        plhs[2] = S.Jth;
    } else {
        mxDestroyArray(S.Jth);
    }

    if (nlhs >= 4) {
        plhs[3] = Zraw ? Zraw : mxCreateDoubleMatrix(0, 0, mxREAL);
    } else if (Zraw) {
        mxDestroyArray(Zraw);
    }

    if (nlhs >= 5) {
        plhs[4] = coreOut[3] ? mxDuplicateArray(coreOut[3]) : mxCreateLogicalScalar(false);
    }
    if (coreOut[3]) {
        mxDestroyArray(coreOut[3]);
    }

    mxDestroyArray(S.z0);
    mxDestroyArray(S.dataPrepared);
}
