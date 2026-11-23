# Xia uniform-input sequential formulations

SAGE Python's new sequential model selections and Project 8 use the
`hourly` branch of <https://github.com/xiaxilin/cartplus>, pinned to
`f986af265633a7a980836c2e0f2124077d632420` (10 November 2026).
This branch changes the daily formulation as well as the hourly formulation.
The implemented upstream defaults are `prod_lie=True`, `nsub=4`, and
`nsub_slow=1`. Arbitrary upstream ablation configurations are not implemented.

Select **Model > Sequential scheme > Uniform input (Xia, Nov 2026)**.
The selector applies to both CaRT and CaRT+. Old reports/sessions without a
scheme flag retain the published/legacy behavior. Project 8 imports explicitly
select the uniform-input scheme. New parameter defaults use Xia's one-hour
minimum X4 (`1/24` day), including daily calibration. Changing resolution preserves every physical parameter bound, including X4;
its one-hour minimum is not divided by the number of intervals per day.
Saved parameter bounds are retained when loading existing sessions.

The flag is `model_config.cart_uniform_input`. New run configurations also
record `model_config.cart_upstream_revision`. Uniform and legacy networks
must be calibrated separately; changing this flag is a model change, not a
way to resume an existing optimization under identical dynamics.

## Implementation

CaRT's existing `cart_uniform` entry point implements the revised snow,
production-spin-up, uniform-input unit hydrographs, quick-store splitting,
and Riccati groundwater update. The release-matching tests cover this path.

CaRT+ now has `cartplus_uniform.hpp` and `cartplus_uniform_catchment.hpp`.
They implement the piecewise exact surface update, proportional snow losses,
mean infiltration wetness, symmetric production/depression spin-up, Lie
production/depression updates, midpoint modulation states, uniform-input unit
hydrographs, groundwater uptake splitting, and the revised quick release
with local exponent and four/two-substep Richardson extrapolation.
Forward sensitivities differentiate all these operations, with spin-up stores
detached exactly as in Xia's reference. The Python reference implementation
is `cartplus_uniform_python.py`; its numerical zone solvers remain available.

`build_uniform.py` creates separate versioned libraries for 42 physical
sensitivities and 35 calibrated sensitivities and writes a source/binary hash
manifest. It does not replace the DLLs used by older or already-running jobs.
Legacy `simulate_us` and the SAGE coupled ODE formulations are unchanged.
SAGE and SITE Python dispatch to these libraries. MATLAB now exposes the
same uniform scheme through isolated MEX gateways built from these sources.
See ../CART_NUMERICS.md for controls, builds and current validation limits.

## Validation

Run from the Software directory:

```powershell
python SAGEhydrology/Native/models/CaRTPlus/build_uniform.py
python SAGEhydrology/Native/models/CaRTPlus/test_xia_uniform.py .codex-work/cartplus-upstream-audit-20261119
```

The second command requires a checkout of the pinned upstream branch. It
executes Xia's unmodified forward kernel bodies through a NumPy CPU shim,
preserving their float32 output rounding. It tests daily, hourly, and
15-minute forcing, including more than four years of daily inputs, cold and
warm periods, full and capped spin-up, and finite-memory filter truncation.
It also compares native/Python discharge and analytic derivatives, checks
the 35-column mapping against the full 42-column Jacobian, verifies all
42 derivatives using five-point differences, and checks the water balance.
This is CPU equivalence testing, not a GPU performance benchmark or a rerun
of the regional calibration experiments.
