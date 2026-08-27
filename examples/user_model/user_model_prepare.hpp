/* DO NOT CHANGE */

#ifndef USER_MODEL_PREPARE_HPP
#define USER_MODEL_PREPARE_HPP

#include "mex.h"

struct UserModelPrepared {
    double tout;
    int d;
    int m;
    int nvar;
    int mem;
    int ipr;

    mxArray* z0;
    mxArray* dataPrepared;
    mxArray* Jth;
};

UserModelPrepared user_model_prepare(const mxArray* par,
                                     const mxArray* mdl,
                                     const mxArray* data,
                                     const mxArray* ode);

#endif
