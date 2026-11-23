# SAGE user model

This folder is the independently compiled user-model plug-in for SAGE. It is
not linked into `crr_model_mex`, because a deployed SAGE application cannot
compile or relink user-supplied C++ source code.

## Files users normally edit

- `make_user_model_info.m`: parameter names, descriptions, units, bounds,
  initial states, and parameter-space selection. Running it creates
  `user_model_info.mat`, which the GUI reads.
- `user_model.cpp` and `user_model_matlab.m`: matching C++ and MATLAB model equations and analytic sensitivities.
- `user_model_configuration.m`: additional-value names, units, source defaults, editability, limits, and references. `make_user_model_info` stores this catalog in the MAT metadata.
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
routes model 99 through the regular MATLAB `crr_model` path, which calls this
independent plug-in.

## Python workflow

Keep this folder at the root of the source checkout or beside
`SAGEhydrology.exe` in the Windows desktop package. For the macOS desktop
package, keep it beside `SAGEhydrology.app`. Restart the desktop app after
editing the files, then select `user_model`. A saved session with model ID
99 can also run through `Python/run_sage.py`.

The first run after changing the C++ sources needs `g++` or `clang++` to
build the plug-in locally. Use `make_user_model_info.m` in MATLAB to regenerate
the metadata after editing the parameter or state definitions.

## Deployed MATLAB workflow

A deployed MATLAB application cannot compile the source files. Distribute or
install all of the following together:

- `user_model_info.mat`
- `user_model.m`
- `user_model_matlab.m`
- `user_model_configuration.m`
- a precompiled `crr_user_model` binary for the target operating system,
  architecture, and compatible MATLAB Runtime release

Do not add the user model to `crr_model_mex`. If the plug-in binary is absent,
SAGE reports that runtime compilation is unavailable instead of attempting to
invoke a compiler.

## Solver and Model-tab configuration

The supplied GR4JB example supports mcode 1 (MATLAB adaptive RK2), 2 (MATLAB ode45, separately over each forcing interval), 3 (fixed-step MATLAB Euler, 1000 steps per record by default, or ode.int_steps), and 4 (C++ adaptive RK2). The model equations and analytic forward sensitivities are the same; numerical solver tolerances and Euler step size determine numerical agreement.

The template metadata declares its full fixed routing-memory layout: 44 states daily and 204 states hourly or at 15-minute resolution, including cumulative-discharge bookkeeping. ode.y0 overrides the initial vector. Additional values use mdl.config. Both MATLAB and C++ honor these structures. Fixed or inactive coefficients remain read-only.

For a different custom model, edit both implementations and its metadata together. Supply configuration_definition with columns Name, Symbol, Unit, Default, Editable, Lower, Reference in user_model_info.mat, and complete state names/units and y0. Do not retain model_template='gr4jb' when replacing its routing/state equations. The GUI reads the metadata rather than trying to infer executable options from comments or a compiled binary. External custom C++ plug-ins require their own matching MATLAB equations before all four solver modes can be offered.
