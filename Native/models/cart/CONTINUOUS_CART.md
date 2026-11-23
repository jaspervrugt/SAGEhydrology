# Coupled continuous-time CaRT

`cart_ode.py` implements the continuous counterpart in Section 2.2 of
`papers/cart_numerics/main.tex`. It is a separate research backend; the SAGE
GUI and existing running jobs continue to use their selected implementation.

```python
from cart_ode import simulate
q, jac, diagnostics = simulate(
    forcing, elevation, area, parameters,
    jacobian=True, refinement=16, return_diagnostics=True,
)
```

Forcing has four columns: precipitation depth (mm/interval), temperature
(degrees C), PET depth (mm/interval), and training indicator (0 or 1).
`resolution` is 1, 24, or 96 intervals/day. Depths become rates internally.
Alternatively, `interval_days` accepts a positive scalar or one duration per
row for arbitrary, including irregular, observation intervals; the recharge
statistic uses the total duration of the marked training intervals.
All times and rate parameters use days, independently of forcing resolution.
The eleven parameter coordinates follow `cart_python.py`, including log
coordinates. Returned q is interval depth (mm); jac has shape (intervals, 11)
and differentiates the same interval depths in those parameter coordinates.

All water stores are integrated together: N1..N5, S, R, W, primary cascade,
symmetric cascade. Refinement m adds m and 2m routing stores. Each stage
releases at rate m/X4; conditional exits reproduce a mixture of Erlang
densities converging to the two target unit hydrographs. Finite m adds
dispersion. The finite cascade is continuous time but is an approximation
to the published convolution kernel, not an exact realization of it.

Smooth sensitivities use forward automatic differentiation of the RHS and
are integrated as augmented states. Snow depletion is located analytically
within each constant-forcing interval; saltation updates all sensitivities,
including output/loss accumulators. Rain/snow threshold coincidences and
depletion at forcing boundaries are flagged because a two-sided derivative
need not exist. Clipped bypass is continuous with a kink; finite-difference
verification requires a solver-tolerance and perturbation-size check.

The default initialization is differentiated: S=.7*C, R=.5*X3, W=gbar*tau,
zero snow and routing prehistory. There is no hidden warm-up. Alternatively,
provide a fixed physical initial vector of length 8+3*m; its derivatives are
zero. The independent upstream recharge pilot uses the same initial snow
and soil and is differentiated even with detached physical initial stores.
At the recharge floor the branch derivative is zero. Model bounds are
X2<=0, X4>0 and snowfall multiplier>=0; band fractions must sum to one.

The five accumulators integrate discharge, atmospheric losses, external
loss, effective precipitation, and recharge. Diagnostics expose full water
store endpoints, interval balance residuals, gbar and its gradient, events,
physical-state sensitivities when requested,
and the actual initialization convention. Negative accepted stores beyond
roundoff cause a failure rather than silent clipping. Conservation does not
alone guarantee a sufficiently accurate trajectory.

Run regression tests from the software workspace:

```
python -m unittest discover -s SAGEhydrology/Native/models/cart -p test_cart_ode.py -v
python SAGEhydrology/benchmarks/cart_numerics/continuous_cart.py
```

The benchmark uses the existing `.codex-work/cartplus-upstream` checkout
(revision recorded in JSON), unmodified author's forward functions through
the existing NumPy/Triton CPU compatibility layer, and installed US/AU/GB
forcing. It checks author-port parity before comparisons. Daily comparisons
match detached author-port S,R,W, use zero snow/cascade prehistory, and let
each formulation compute its own training-recharge statistic. Output CSV
and JSON record measured results and versions. Single CPU timing samples
do not measure GPU performance. Synthetic subdaily checks subdivide the
same forcing; they are not validation with measured subdaily rainfall.

The native binding `cart_ode_cpu.simulate` uses the same inputs and
integrates all eleven calibrated sensitivities with `jacobian=True`.
It accepts `solver='dop853'` (default), `'rk2'`, and `'rk23'`. From the
software workspace, build the Windows library with:

```powershell
g++ -std=c++17 -O3 -march=x86-64-v3 -shared -static SAGEhydrology/Native/models/cart/cart_ode.cpp -o SAGEhydrology/Native/models/cart/cart_continuous_cpu_v2.dll
```

Run `python SAGEhydrology/benchmarks/cart_numerics/coupled_runge_kutta.py`
to compare all three solvers on six ten-year CaRT and CaRT+ records.
The output CSV reports time for discharge and Jacobian together,
with differences from tight DOP853 calculated after one year of spin-up.
