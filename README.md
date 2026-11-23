# SAGEhydrology v2.0.0

Sensitivity-Aware Gradient Estimation (SAGE) trains conceptual hydrologic
models across basins using analytic forward sensitivities and
attribute-conditioned parameter regionalization.

This release adds headless Python alongside updated MATLAB source. Both use
the shared public C++ kernels in `Native/`. The desktop interface is
distributed as a separately licensed compiled application through
[GitHub Releases](https://github.com/jaspervrugt/SAGEhydrology/releases).

## Python

Use Python 3.11 or later and a C++17 compiler available as `g++`. From the
checkout directory, run:

```sh
python -m pip install -r Python/requirements.txt
python Python/run_sage.py --build-native
python Python/examples/synthetic_training.py --output ../sage-example
python Python/check_candidate.py
```

The example trains HYMOD for three iterations on four synthetic basins.
For an existing JSON session, update the data paths and run
`python Python/run_sage.py --session session.json --output ../sage-run`.
Use a fresh output directory for each run. The headless Python session runner
supports HYMOD, HMODEL, SAC-SMA, Xinanjiang, GR4J-A, HBV, CFE-NWM and
the editable `user_model`.
The shared CaRT and CaRT+ model sources and adapters are also included under
`Native/models/`.

## MATLAB

Open `Matlab/examples/demo_SAGE.mlx`. Review the installation path, region,
date range and local forcing, discharge and attribute paths before running.
The examples and project scripts are under `Matlab/examples/` and
`Matlab/projects/`. Keep the sibling `Native/` directory in place. Configure
a supported C++ compiler with `mex -setup C++` if using MEX kernels. For
the shared CaRT numerical gateways, add `Matlab/models` to the MATLAB path
and run `build_cart_numerics`.

Obtain provider datasets separately and follow their licence and citation
terms. SAGE's regional files specify the supported layouts and basin lists.

## Custom model

The top-level [user_model](user_model/) directory contains a working GR4JB
template with MATLAB and C++ equations and analytic sensitivities. The
Windows desktop ZIP includes the same editable folder beside the executable.
Select `user_model` in the desktop model list. Python builds the edited
sources on first use with a local C++ compiler. The plug-in README covers
MATLAB setup and changes to model metadata.

## Citation and licensing

See [CITATION.cff](CITATION.cff) for software citation information.
The public computational source uses the [BSD 3-Clause License](LICENSE).
Keep the adjacent MIT notices for CaRT and CaRT+ and the notices for other
third-party components. See [GUI-LICENSE-NOTICE.md](GUI-LICENSE-NOTICE.md)
for the terms of the compiled application.

[Project blog](https://jaspervrugt.github.io/SAGE-SITE-blog/#welcome)
