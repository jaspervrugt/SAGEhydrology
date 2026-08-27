# SAGE user model

This folder is the independently compiled user-model plug-in for SAGE. It is
not linked into `crr_model_mex`, because a deployed SAGE application cannot
compile or relink user-supplied C++ source code.

## Files users normally edit

- `make_user_model_info.m`: parameter names, descriptions, units, bounds,
  initial states, and parameter-space selection. Running it creates
  `user_model_info.mat`, which the GUI reads.
- `user_model.cpp`: model equations, routing, and analytic sensitivities.
- `user_model_prepare.cpp`: parameter conversion, model constants, routing
  ordinates, dimensions, and augmented initial conditions.

The public solver interface is declared in `user_model.hpp`. The fixed MATLAB
gateway is `crr_user_model_mex.cpp`.

## MATLAB desktop workflow

1. Edit the model-specific files.
2. Run `make_user_model_info`.
3. Configure a supported compiler with `mex -setup C++` if necessary.
4. Run `user_model_compile`.
5. Add this folder to the MATLAB path and select `user_model` in SAGE.

The compiler creates the platform-specific `crr_user_model` MEX binary. SAGE
routes model 8 through the regular MATLAB `crr_model` path, which calls this
independent plug-in.

## Deployed SAGE workflow

A deployed MATLAB application cannot compile the source files. Distribute or
install all of the following together:

- `user_model_info.mat`
- `user_model.m`
- a precompiled `crr_user_model` binary for the target operating system,
  architecture, and compatible MATLAB Runtime release

Do not add the user model to `crr_model_mex`. If the plug-in binary is absent,
SAGE reports that runtime compilation is unavailable instead of attempting to
invoke a compiler.
