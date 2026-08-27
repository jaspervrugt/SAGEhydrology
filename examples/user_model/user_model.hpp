/*
 * user_model.hpp
 *
 * Solver interface for a user-supplied SAGE rainfall-runoff model. This
 * header is shared by user_model.cpp and the standalone crr_user_model MEX
 * gateway. Users normally edit the model-specific implementation, not this
 * interface.
 *
 * Written by Jasper A. Vrugt, May 2026.
 * University of California, Irvine.
 */

#pragma once

#include "mex.h"

/*
 * Core solver entry point for the modular crr_user_model MEX.
 *
 * The gateway/user preparation layer supplies the old low-level inputs
 * that the original crr_user_model.cpp mexFunction used to receive:
 *
 *   t_last  : final print time, equal to mdl.tout
 *   z0_mx   : 1 x nvar or nvar x 1 initial augmented state vector
 *   data    : prepared data struct containing P, Ep, T, parameters,
 *             UH ordinates, model constants, and data.ipr
 *   options : ODE settings struct
 *
 * Outputs:
 *   plhs[0] = Z
 *   plhs[1] = q_n    when requested and options.mem == 0
 *   plhs[2] = J_raw  when requested and options.mem == 0
 *   plhs[3] = fail   when requested
 */
namespace sage_user_model {

void run(int nlhs,
         mxArray* plhs[],
         double t_last,
         const mxArray* z0_mx,
         const mxArray* data,
         const mxArray* options);

} // namespace sage_user_model
