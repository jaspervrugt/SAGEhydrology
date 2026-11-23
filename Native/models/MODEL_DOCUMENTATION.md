# CaRT and CaRT+ source documentation

The CaRT and CaRT+ kernels follow the older SAGE model convention:
file identity and attribution, function purpose, input arguments,
output arguments, and formulation-specific notes. Source lines are
limited to 76 characters, including comments and documentation.

## Source locations

- `Native/models/cart`: CaRT and the related CaRTsplit native wrappers.
- `Native/models/CaRTPlus`: sequential and coupled CaRT+ implementations.
- `Matlab/models/cart` and `Matlab/models/cartplus`: MATLAB models.
- `Julia`: standalone coupled model-release sources. These files are
  separate from the SAGE Python package; no Julia package dependency
  was introduced by this documentation update.

The update covers 53 model, helper, build, and validation source files.
It preserves existing upstream license and reference notices.

## Interface conventions

Forcing inputs use `[P, T, PET, training]` in the native model wrappers.
Precipitation and PET are depths per observation interval (mm).
Temperature is in degrees Celsius. Five band elevations are in metres;
their nonnegative area fractions sum to one.

`resolution` is the number of forcing intervals per day: 1, 24, or 96.
Continuous integration steps and explicit interval durations use days.
Native pointer-buffer layout and optional outputs follow the individual
wrapper; MATLAB gateways perform their own array marshalling.

CaRT uses 11 calibration coordinates, including logarithmic entries.
CaRTsplit adds a twelfth routing-partition coordinate. CaRT+ uses 42
physical parameters in the order defined by `cartplus_inputs.PARAMETER_MAP`.
The latter maps named static, law, and shared parameter dictionaries to
the physical kernel order. Selected native builds retain only 35
calibrated derivative directions.

Initialization and differentiation are formulation-specific. Sequential
fixed-spin-up Jacobians and coupled pilot/initialization sensitivities
must not be interchanged. The function notes describe the chosen policy.
Supplying detached initial stores is an explicit derivative convention.

## Validation of this documentation update

- All 53 source files satisfy the 76-character limit.
- Python implementation syntax trees match the originals after excluding
  documentation strings. Every Python function has input/output help.
- C++ token sequences match the originals after excluding comments and
  whitespace. Nine base/alternate syntax-check builds passed.
- MATLAB numerical expressions match the originals after excluding
  documentation, whitespace, and continuation formatting.
- CaRT and CaRT+ discharge/Jacobian regressions passed. CaRT+ checks cover
  daily, hourly, and quarter-hourly forcing and solver cases 1 through 4.
- Julia CaRT and CaRT+ discharge, Jacobians, and balance residuals matched
  the original sources bit for bit on the comparison fixture.
- MATLAB runtime validation was unavailable: startup could not contact
  required MathWorks services (error 5202).

Existing compiled libraries and build provenance were not regenerated.
The source edits change documentation and formatting, not model equations.
