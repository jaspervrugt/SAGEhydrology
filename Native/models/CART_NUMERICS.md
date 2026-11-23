# CaRT and CaRT+ numerical options

SAGE and SITE desktop controls now distinguish the model, formulation,
sequential scheme, continuous ODE method, and forcing interval.

| Formulation | Selection | Configuration |
|---|---|---|
| Sequential | Published daily scheme (legacy) | `cart_uniform_input=false` |
| Sequential | Uniform input (Xia, November 2026) | `cart_uniform_input=true` |
| Continuous ODE | DOP853 | `cart_solver="dop853"` |
| Continuous ODE | Adaptive Heun RK2 | `cart_solver="rk2"` |
| Continuous ODE | Bogacki–Shampine RK23 | `cart_solver="rk23"` |
| Continuous ODE | Dormand–Prince RK45 | `cart_solver="rk45"` |
| Continuous ODE | Fixed Heun RK2 | `cart_solver="heun_fixed"` |
| Continuous ODE | Fixed Euler | `cart_solver="euler"` |

The continuous ODE equations remain distinct from the sequential schemes;
selecting another integrator does not turn one formulation into the other.
The existing `rk2` identifier retains its adaptive meaning.

## Desktop controls

In both applications and languages, choose **Model > Formulation** and,
for sequential models, **Model > Sequential scheme**. Continuous methods
are selected on **ODE | Solver**, with relative/absolute tolerances for
adaptive methods and maximum/fixed step in **days**. Minimum adaptive step
is an error threshold: the solver fails if it requires a smaller step;
it does not silently accept an inaccurate trial. A final boundary segment
may be shorter than this minimum.

Python's blank adaptive maximum means automatic; MATLAB uses `1e100` for
an effectively unrestricted maximum. New RK45/fixed selections without an
explicit step default to `1/96` day. Fixed methods divide each forcing or
CaRT snow-event segment into equal steps no larger than the chosen step.
Tolerances do not control fixed-step accuracy. CaRT fixed-step derivatives
include differentiation of the event-delimited segment lengths.

Forcing resolution is selected separately on Period. Daily/hourly/
15-minute forcing is not the integration step. The ODE routing cascade
refinement (`cart_refinement`, default 16 UH1 stages) is also independent
of temporal integration accuracy.

Uniform input retains Xia's pinned `f986af2` defaults: production Lie
splitting, four/two quick-store substeps with Richardson extrapolation,
and one slow substep. These are displayed as fixed implementation details,
not adjustable controls. New GUI selections default to uniform input;
historical sessions without a scheme flag retain the legacy scheme.
The MATLAB uniform scheme and continuous solvers use native C++ through
MEX; they do not require Python.

## Preservation and provenance

New DLLs and MEX gateways use separate names. Existing legacy binaries
and already-running workers are not replaced. Historical Python adaptive
calls without new step controls continue to use their previous libraries.
Saved configurations carry the scheme and solver settings; Python worker
configurations additionally record the chosen DLL name and SHA-256.
Build manifests record the newly compiled source/binary hashes. MATLAB
result identities include the effective numerical settings and revision.
Uniform and legacy networks require separate calibration; switching
formulations is not a continuation of the same training problem.

## Verification and paper scope

`test_cart_numerics.py` checks new native methods, water balance, directional
gradients, fixed-step convergence, 35-column sensitivity mapping, invalid
step settings, and daily/hourly/15-minute subdivision for RK45.
`test_cart_solver_options.py` checks the actual Qt edit/visibility handlers
in each application. Existing scheme-selection tests cover legacy loading
and calibrated parameter mapping.

`make_cart_matlab_fixture.py` writes 48 Python expectations covering both
models, all six ODE methods, both sequential schemes and all three forcing
resolutions. Run `test_cart_numerics_parity` in MATLAB to compare its public
wrappers against these expectations. MEX compilation is verified here;
MATLAB execution is currently blocked by MathWorks services error 5202,
so MATLAB runtime and graphical parity remain unverified.

For the numerics paper, daily calibration plus hourly/subhourly forward
convergence, conservation and timing experiments is sufficient for the
proposed numerical comparison. No extensive hourly training was started.
Euler is primarily a convergence diagnostic; its step must resolve the
fastest model process. At coarse steps, piecewise rate-law switches can
make finite-difference gradient checks sensitive to the perturbation.

## Builds

From the Software directory:

```powershell
python SAGEhydrology/Native/models/build_cart_numerics.py
python SAGEhydrology/Native/models/CaRTPlus/build_uniform.py
python SAGEhydrology/Native/models/test_cart_numerics.py
python SAGEhydrology/Native/models/make_cart_matlab_fixture.py
```

MATLAB: add `SAGEhydrology/Matlab/models` to the path, run
`build_cart_numerics`, then add `Matlab/tests` and run
`test_cart_numerics_parity`. A Windows standalone build helper,
`build_cart_mex_windows.py`, can use an installed MATLAB SDK and MinGW
without starting MATLAB.
