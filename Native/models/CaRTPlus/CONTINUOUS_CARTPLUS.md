# Fully coupled continuous-time CaRT+

`cartplus_continuous.py` is a separate Python research reference. All
catchment states, both routing cascades, quick storage and flux integrals
are advanced together with DOP853. The older `cartplus_ode.py/.cpp` remains
the continuous-store hybrid with discrete routing. Existing GUI selections
and running jobs are unaffected by this separate backend.

```python
from cartplus_continuous import simulate
q, jac, diagnostics = simulate(
    forcing, elevation, area, physical_parameters,
    jacobian=True, sensitivity_coordinates='calibration',
    refinement=16, return_diagnostics=True,
)
```

Inputs are precipitation depth (mm/interval), temperature (degrees C), PET
depth (mm/interval), and a 0/1 training mask. Five elevation-band fractions
must sum to one. The 42 physical parameters follow
`cartplus_inputs.PARAMETER_MAP`. Times and all rate laws use days.
`resolution=1,24,96` selects standard intervals; `interval_days` can instead
be a positive scalar or an arbitrary vector of observation durations.
Discharge q is reported interval depth, after integrating channel-loss
rates. The Jacobian is (intervals,42), with physical-coordinate columns by
default or upstream calibration-coordinate columns when requested. The
latter uses dp/dtheta=p for log parameters, including seven shared inputs.

The state vector is N1..N5,H,S,N,B,A,J,V,W,R followed by m primary routing
stores and 2m symmetric routing stores. Diagnostic wetness B and memories
J,V are excluded from basin water storage. Conditional stage exit fractions
form mixtures of Erlang densities approximating the published routing
kernels; their rates are m/routing_days. Connectivity is instantaneous,
including in quick release; no interval-mean connectivity or discrete
hydrograph convolution remains. Channel loss is an instantaneous rate,
not an operation on an interval depth after integration.

The existing upstream rate regularizations remain explicit: snow/surface
availability time .02 days, one-percent surface transition and Pareto
retention transition. They are fixed constitutive controls, not solver
steps and not scaled by forcing resolution. They can be varied through
`availability_days` and `transition_fraction`. Thus this is a specified
continuous counterpart of the original daily algorithm, not an assertion
of a unique continuous interpretation. No artificial finite-step updates
are hidden inside the ODE. Internal min/max branches have continuous fluxes
and no state resets; phase-threshold coincidences at forcing boundaries
are flagged because two-sided parameter derivatives may not exist.

Default initialization differentiates capacity fractions, the
duration-weighted first-year vegetation target and Wref. Explicit initial
physical states (length 14+3m) are detached. No hidden spin-up is performed.
The independent upstream pilot uses the same initial upstream states and
is differentiated to obtain training mean recharge and Wref. Training
means use total depth divided by marked duration. The recharge floor has
zero derivative on its active branch, while Wref retains its residence-
time dependence. Full state sensitivities and Wref's gradient are exposed
in diagnostics, alongside water-balance residuals and channel/pre-loss
interval depths. Water conservation does not alone establish trajectory
accuracy.

Verification:

```
python -m unittest discover -s SAGEhydrology/Native/models/CaRTPlus -p test_cartplus_continuous.py -v
python SAGEhydrology/benchmarks/cart_numerics/continuous_cartplus.py
```

The driver verifies all 42 Jacobian columns on a synthetic forcing record,
forcing subdivision, tolerance and cascade refinement, and daily
comparisons on two US, two GB and two Australian basins. The original
author's unmodified forward functions run through the existing CPU shim;
this checks algorithm parity, not GPU rounding or GPU runtime. Comparisons
match seven captured initial physical stores, zero snow/cascade prehistory
and stated memory initializations. Each formulation retains its own
training recharge statistic. Native-vs-ODE hydrograph/Jacobian differences
are formulation differences, not errors against observations. All basins
use illustrative common parameters, without recalibration.

CSV/JSON checkpoints record revision, versions, parameters, initial-state
convention, backend hash and measured results. `--resume` reuses completed
comparisons only for an unchanged backend. Timings are single CPU samples,
including the pilot; the Python reference is substantially slower than the
compiled daily model. Measured subdaily validation, regularization
sensitivity, long-record Jacobian comparisons and optimization remain
further work.

### C++ and MATLAB implementations

The coupled C++ and MATLAB solvers use the same DOP853 tableau and integrate
the same stores, cascades, ledgers, and 42 tangent directions. C++ uses
forward dual numbers. MATLAB advances the tangent directions together in
vectorized complex-step lanes, with the real trajectory selecting the
active branches. Both support detached initial states and irregular interval
durations. The derivative of the upstream recharge pilot is included.

Build the Windows DLL from the Software directory with MinGW:

```powershell
g++ -std=c++17 -O3 -shared -static -o SAGEhydrology/Native/models/CaRTPlus/cartplus_continuous_cpu.dll SAGEhydrology/Native/models/CaRTPlus/cartplus_continuous.cpp
```

Import `simulate` from `cartplus_continuous_cpu` with the same inputs
as the Python reference. Set `calibrated_only=True` to integrate the 35
fitted directions and hold the seven shared constants fixed. Build that
DLL with the same command, adding `-DCARTPLUS_CALIBRATED_ONLY` and naming
the output `cartplus_continuous_cpu_35.dll`. The timing builds use
`-march=x86-64-v3` for the measured CPU. The default build command above
also works without that instruction-set requirement.

The C++ binding's optional diagnostics return the interval
balance residual. For MATLAB, add `SAGEhydrology/Matlab/models/cartplus`
to the path and call `[q,J,balance] = cartplus_continuous(data,p)`.
The `data` struct accepts forcing, elevation, area and either resolution
or interval_days. Set `sensitivities=false` to request only discharge and
balance. Pass `seed` to rescale the default physical sensitivities.

Run the matched comparison with a fresh scratch directory:

```powershell
python SAGEhydrology/benchmarks/cart_numerics/coupled_runtime.py .codex-work/cartplus-runtime-new
```

The benchmark uses 16 primary and 32 symmetric routing stores and shared
tolerances, first trial step, inputs, and initial conditions. It checks
all 42 calibration-coordinate derivatives on six 120-day records and
discharge on three 10-year records. The CSV reports elapsed time per
simulation, water balance, and cross-language differences.

Use `--stage cpp-decades` and `--stage cpp-decades-35` with the completed
benchmark directory to time the two C++ Jacobian widths on all six decades.
These stages check the short-record Jacobians against the saved reference
and report the decade discharge agreement and water balance.

For the adaptive Heun (RK2) test, rebuild the DLLs and call the C++
binding with `solver='rk2', jacobian=True, calibrated_only=True`.
The embedded Euler difference controls the error in both states and
sensitivities. DOP853 remains the default. Use `--stage rk2` with the
completed benchmark directory to compare three RK2 tolerance pairs
against tight DOP853 on the US, GB, and Australian decade records.
The comparison writes `coupled_rk2.csv`.

The same native binding also accepts `solver='rk23'` for the
BogackiÃ¢â‚¬â€œShampine 3(2) pair. An RK23 trial outside the CaRT+ storage
acceptance bounds is retried at half the step size. The joint CaRT and
CaRT+ decade comparison is in `benchmarks/cart_numerics/coupled_runge_kutta.py`.

### Experimental native RK2 variants

The Python CPU binding accepts the following opt-in variants. The standard
`solver='rk2'` entry point and its DLLs are retained.

- `rk2_workspace` reuses the Euler and Heun work arrays over forcing
  intervals in both the recharge pilot and the full solve.
- `rk2_exact_vegetation` also evaluates the analytic vegetation relaxation
  at the RK2 stage times. It integrates the same governing equations with
  a different discretization and adaptive step sequence.
- `rk2_analytic_routing` propagates the existing linear cascades with a
  constant inlet in the embedded Euler estimate and a linear inlet in
  the Heun estimate. The downstream stages use the corresponding mean
  releases from the cascade water balance. This changes the numerical
  discretization and may change discharge and sensitivities. It is a
  benchmark option and is not the default solver.

- `rk2_replay` records the accepted state-only Heun steps in the recharge
  pilot and the full solve, then recomputes the stages with tangents on
  that fixed mesh. Discharge and the balance diagnostic come from the
  state-only pass. With diagnostics enabled, `rk2_trials` gives accepted
  and rejected counts in two rows (pilot, full solve). A state-only mesh
  can produce unstable sensitivities. Use this option for benchmarking
  with a coupled reference.

- `rk2_optimized` reuses RK2 arrays and stops the recharge pilot after
  the last training interval. It retains the coupled error controller.
- `rk2_analytic_quick` additionally uses the closed-form power-law
  recession with a coefficient frozen from the two RK2 stages. It
  splits the mean non-recession flow around that update, uses the
  integrated release in the channel-loss calculation, and falls back
  to the standard RK2 trial if the split would make quick storage
  negative. This is an experimental discretization.

Both options report accepted and rejected RK2 trial counts for the
pilot and full solve when diagnostics are enabled.

The variants accept daily, hourly and 15-minute forcing via
`resolution=1`, `24` or `96`, as well as explicit `interval_days`.
Use `jacobian=True, calibrated_only=True` for the 35 calibration
directions. Build the separate experimental DLLs from the Software directory:

```powershell
$p = 'SAGEhydrology/Native/models/CaRTPlus'
g++ -std=c++17 -O3 -march=x86-64-v3 -shared -static -DCARTPLUS_CALIBRATED_ONLY -o "$p/cartplus_continuous_cpu_35_experimental.dll" "$p/cartplus_continuous.cpp"
g++ -std=c++17 -O3 -march=x86-64-v3 -shared -static -o "$p/cartplus_continuous_cpu_experimental.dll" "$p/cartplus_continuous.cpp"
```

### Production RK2 work reuse and pilot window

The standard native RK2 entry point now reuses its two workspaces and
stops the recharge pilot after the last training interval. The full
simulation still covers every requested interval. Python prefers the
versioned production DLLs `cartplus_continuous_cpu_v3.dll` and
`cartplus_continuous_cpu_35_v3.dll` when present. These builds omit
the experimental scalar, constant-cache, and power macros.

The pure Python and MATLAB continuous references and Julia coupled
solver also stop the pilot after the last training interval. Julia
already allocates its integration workspace once per pilot/full pass.
MATLAB now reuses separate pilot/full DOP853 stage buffers. Its reference
integrator uses DOP853, so the buffer change preserves that method.
The Python reference uses SciPy DOP853 and its internal workspace.
