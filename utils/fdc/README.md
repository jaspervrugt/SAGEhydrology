# Flow-duration-curve utilities

This directory contains the basin-specific flow-duration-curve observation
operators used by SAGEhydrology:

- `fit_kosugi_fdc.m`: robust single-Kosugi fit;
- `fit_dual_kosugi_fdc.m`: five-parameter dual-Kosugi mixture fit and inverse cache;
- `dual_kosugi_objective_mex.cpp`: compiled objective kernel used by the dual fit;
- `build_dual_kosugi_mex.m`: platform-local MEX build helper.

Run `build_dual_kosugi_mex` on each supported platform to produce the matching
MEX binary. If no compatible binary is available, the dual fitter automatically
falls back to its MATLAB objective implementation.

General FDC loss and performance-metric functions remain in `utils/metrics`.
