# SAGEhydrology

Sensitivity-Aware Gradient Estimation (SAGE) is a framework for scalable,
attribute-conditioned training of conceptual hydrologic models using analytic
forward sensitivities.

This repository contains the public computational source for SAGEhydrology.
The graphical user interface is distributed as a compiled application through
the [GitHub Releases](https://github.com/jaspervrugt/SAGEhydrology/releases)
page; its source code is not part of this repository.

The SAGE GUI can automatically download, extract, organize, and register the
supported regional hydrologic and meteorological datasets. Users therefore do
not need to locate and arrange these data files manually for the standard
regional workflows.

## Repository contents

```text
docs/                    Documentation graphics and model schematics
examples/demo_SAGE.mlx   Illustrated SAGE Live Script
examples/user_model/     Editable user-model plug-in template
flags/                   Regional flag assets
maps/                    Map assets and Natural Earth metadata
models/                  Hydrologic models and analytic sensitivity kernels
regions/                 Regional configuration and basin inventories
src/                     Main SAGE training and postprocessing routines
utils/                   Shared readers, metrics, plotting, and utilities
```

Hydrologic and meteorological datasets, run results, caches, and GUI source
files are intentionally excluded.

## Installation

### Compiled graphical application

Download the appropriate installer from GitHub Releases. Windows builds use
MATLAB Runtime R2026a, which the installer can obtain from MathWorks. macOS
builds are published separately because standalone applications and MEX files
are platform-specific.

After installation, select a supported region and temporal resolution in the
GUI. SAGE identifies missing data and offers the corresponding download and
installation controls. The GUI downloads the source archives, extracts them,
applies the directory and naming conventions expected by SAGE, and prepares
the data for basin selection and quality screening. Manual data installation
remains available for users who already maintain local dataset copies.

### MATLAB source

1. Clone or download this repository.
2. Use the GUI's automatic data installation, or place existing regional
   datasets manually in the layout selected by SAGE.
3. Open `examples/demo_SAGE.mlx`, select the SAGE root directory, and review
   the configuration before running it.
4. Compile platform-specific MEX kernels when required by the selected model
   and execution backend.

### User-defined hydrologic model

The repository includes a clean user-model template in
`examples/user_model/`. Copy this directory to the SAGE workspace root so that
it sits next to `Data/` and `SAGEhydrology/`:

```text
SAGE workspace/
|-- Data/
|-- SAGEhydrology/
`-- user_model/
```

Edit the copied C++ model and metadata files, run `make_user_model_info`, and
compile the platform-specific plug-in with `user_model_compile`. Compiled MEX
files are intentionally excluded from the template because they must match the
user's operating system and MATLAB release. The compiled GUI distribution can
use the same external `user_model/` workspace directory.

The software does not bundle or redistribute CAMELS and other regional
datasets inside the source repository or application installer. Instead, the
GUI automates retrieval from the supported original data sources. Dataset use
remains subject to each provider's availability, citation requirements, and
license.

## Citation

Please cite the relevant SAGE publications listed in `CITATION.cff`. Additional
paper-specific citation information is included in the Live Script.

## Licensing

The computational source in this repository is licensed under the BSD
3-Clause License; see `LICENSE`.

Compiled SAGE GUI applications are separately licensed and are not covered by
the repository's BSD license. MATLAB Runtime and third-party datasets/assets
remain subject to their respective licenses.

## Contact

Jasper A. Vrugt  
University of California, Irvine  
jasper@uci.edu
