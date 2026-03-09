# SAGEhydrology: SAGE + SITE Software Report (with GUI)
**Version:** v1.0  
**Maintainer:** Jasper A. Vrugt (UC Irvine)

<img width="1524" height="780" alt="image" src="https://github.com/user-attachments/assets/bd058ee2-c78d-44a8-9c1a-59cf9c438dfb" />

This report describes the SAGEhydrology software distribution, its graphical user interfaces (GUIs), expected directory structure, and how to install external CAMELS data so that both the GUIs and the Live Script demos work out of the box.

---

## 1. Overview

**SAGEhydrology** is a MATLAB-based software package for large-sample hydrologic modeling, training, and evaluation across CONUS basins. The distribution provides two complementary workflows:

- **SAGE:** large-sample training/benchmarking across many basins (and postprocessing/visualization)
- **SITE:** single-site (basin-by-basin) workflows for model training and inspection

Both workflows include:
- a **GUI** for configuration and execution
- source and helper code in `src/`
- output written to `results/`
- **Live Script demos** for daily and hourly data

---

## 2. Installation and Quick Start

### 2.1 Requirements

#### Core SAGEhydrology / SITE training code (non-GUI)
- **MATLAB:** developed and tested in **R2025b (Version 25.2)**.
- The core modeling and calibration pipeline uses standard MATLAB language features and built-in functions.
- **Optional toolbox:**
  - **Parallel Computing Toolbox** (optional): enables `parpool`, `parfor`, and queue-based streaming used by
    `setup_parpool`, `run_CONUS`, and `run_model`. If unavailable, the code runs in **serial** mode.
- **Compilers / external:**
  - **MEX compiler** (required only for C/C++ MEX acceleration; optional if running MATLAB-only solvers)

#### Graphical User Interfaces (SAGE_ui, SITE_ui)
- **MATLAB:** developed and tested in **R2025b (Version 25.2)**.
- The GUIs rely on the modern UI framework (`uifigure`, `uigridlayout`, etc.) and additionally use:
  - `jsonencode` (SAGE_ui and SITE_ui) for configuration preview
  - `exportgraphics` (SITE_ui) for figure snapshotting into the GUI

**Recommended minimum versions**
- **Core (non-GUI):** **R2019b+**
- **GUI (SAGE_ui / SITE_ui):** **R2020b+**

### 2.2 Install the code

SAGEhydrology is distributed as a zip file. The recommended convention is to choose a **root directory** that will contain both:

- `root/SAGEhydrology/`  (code)
- `root/Data/`           (external CAMELS data)

1. Create (or choose) a folder that will be your `root/`.
2. Place the `SAGEhydrology.zip` file inside `root/`.
3. Extract/unzip it **in place**. After extraction, you should have:

```text
root/
  SAGEhydrology/
    SAGE/
    SITE/
    models/
    utils/
    ...
  Data/              (may be created by the installer, or created later)
    daily/
    hourly/
```

If `Data/` is not created automatically, you can create it manually.

### 2.3 Install the data (required before running SAGE/SITE)

SAGEhydrology does **not** ship with CAMELS data. You must download the CAMELS datasets separately and place them in `Data/`.

Two mechanisms are supported:

#### Mechanism 1: Manual (self-installation)
Follow the download, unzip, and rename procedures in Section **7**.

#### Mechanism 2: GUI-assisted download and installation (SITE_ui)
If you prefer not to manually download and rename folders, the **SITE GUI** provides buttons that call the installer and place files into the canonical layout.

1. In MATLAB, go to:
   - `root/SAGEhydrology/SITE/gui/`
2. Run:
   - `SITE_ui`
3. Go to the **Paths** tab.
4. Under the **Data** path, use the download buttons:
   - **Download Daily Data**
   - **Download Hourly Data**
5. Button visibility and enable/disable logic:
   - The **Daily** button is enabled only if daily data are not detected at  
     `root/Data/daily/v1p2/forcing` and `root/Data/daily/v1p2/streamflow`.
   - The **Hourly** button is enabled only if hourly data are not detected at  
     `root/Data/hourly/forcing` and `root/Data/hourly/streamflow`.
   - The **Hourly** button may only become visible after selecting **Hourly** data under the **Period** tab.

**Notes**
- Hourly downloads are large (~18 GB compressed and much larger after extraction and installation).
- Daily downloads are smaller (~3.3 GB compressed and much larger after extraction and installation).
- The GUI writes progress messages to the Log tab during installation.
- The installer does not redistribute CAMELS; it downloads from the official URLs listed in Section 7.

### 2.4 Run a demo (requires CAMELS data)

**Important:** you must install CAMELS data first (Section **2.3**). Otherwise:

- The GUIs will open, but **Run SAGE** / **Run SITE** will not produce results because forcing and streamflow data are missing.
- The Live Script demos will error when they reach `read_meteo` because the required data are missing.

After data are installed, you can run:

- `root/SAGEhydrology/SAGE/examples/demo_SAGE_daily.mlx`
- `root/SAGEhydrology/SAGE/examples/demo_SAGE_hourly.mlx`
- `root/SAGEhydrology/SITE/examples/demo_SITE_daily.mlx`
- `root/SAGEhydrology/SITE/examples/demo_SITE_hourly.mlx`

Or launch the GUIs:
- `root/SAGEhydrology/SAGE/gui/SAGE_ui.m` (or `SAGE_ui()`)
- `root/SAGEhydrology/SITE/gui/SITE_ui.m` (or `SITE_ui()`)

---

## 3. Canonical Directory Structure

### 3.1 Top-level layout (recommended)

The **root directory** is the folder containing both `Data/` and `SAGEhydrology/`:

```text
root/
  Data/
    daily/
    hourly/
  SAGEhydrology/
    basins/
    docs/
    models/
    SAGE/
    SITE/
    utils/
```

**Important:** The SAGE and SITE software assume this layout by default.

### 3.2 `Data/` layout (external; user-provided)

The `root/Data/` directory contains CAMELS-derived datasets and is **not included** with SAGEhydrology.

```text
Data/
  daily/
    ... (daily forcing, streamflow, attributes, etc.)
  hourly/
    ... (hourly forcing, streamflow, attributes, etc.)
```

### 3.3 `SAGEhydrology/` layout (shipped with the code)

```text
root/SAGEhydrology/
  basins/        % basin selection files (e.g., 516 and 531 basin lists)
  docs/          % model schematic figures (PNG/PDF) used by the GUI
  models/        % model implementations organized by model name
  SAGE/          % large-sample workflow (GUI + src + results + examples)
  SITE/          % single-site workflow (GUI + src + results + examples)
  utils/         % shared utilities used by both SAGE and SITE
```

### 3.4 SAGE workflow layout

```text
root/SAGEhydrology/SAGE/
  examples/
    demo_SAGE_daily.mlx
    demo_SAGE_hourly.mlx
  gui/
    ... (SAGE_ui)
  src/
    ... (core functions)
  results/
    ... (default run output location)
```

### 3.5 SITE workflow layout

```text
root/SAGEhydrology/SITE/
  examples/
    demo_SITE_daily.mlx
    demo_SITE_hourly.mlx
  gui/
    ... (SITE_ui)
  src/
    ... (core functions)
  results/
    ... (default run output location)
```

---

## 4. How the GUI Locates Files

The SAGE and SITE GUIs use a **single root directory** selected by the user. Under the recommended layout:

- code is found at: `root/SAGEhydrology/...`
- data is found at: `root/Data/...`

If the GUI cannot locate required files, runs will fail when forcing or streamflow files are accessed. The expected layout is:

- `root/SAGEhydrology/...`
- `root/Data/...`

---

## 5. Typical User Workflows

### 5.1 GUI workflow
1. Set root directory
2. Select daily or hourly data source
3. Choose model and options
4. Select basins and training/validation settings
5. Click **Run**
6. Review results and plots; outputs are written to:
   - `root/SAGEhydrology/SAGE/results/` (SAGE)
   - `root/SAGEhydrology/SITE/results/` (SITE)

### 5.2 Advanced workflow (Live Scripts)
1. Open the appropriate demo in:
   - `root/SAGEhydrology/SAGE/examples/` or `root/SAGEhydrology/SITE/examples/`
2. Run sections interactively
3. Adjust configuration fields as needed (paths, basins, model selection, solver options)

---

## 6. Outputs and Reproducibility

SAGEhydrology writes results to the `results/` directories by default. For reproducibility, users should retain:
- the configuration used (paths, model choice, basin sets, solver settings)
- model parameters and training history
- performance metrics and summary figures

---

## 7. Getting CAMELS Data (manual procedure + official sources)

CAMELS data must be downloaded separately from official sources. After downloading, place the extracted data under:

- `root/Data/daily/`
- `root/Data/hourly/`

### 7.1 CAMELS-US Daily Daymet/Maurer/NLDAS forcing + daily USGS streamflow

**Official sources**
- CAMELS product page (NCAR/UCAR RAL):  
  https://ral.ucar.edu/solutions/products/camels
- CAMELS-US daily data (Zenodo record / direct download):  
  https://zenodo.org/records/15529996
- Dataset DOI:  
  https://dx.doi.org/10.5065/D6MW2F4D

**What to do**
0. Make sure you have at least 30 GB of free space on your computer.
1. Download the CAMELS-US daily dataset from the Zenodo record above. In particular, download:
   - `basin_timeseries_v1p2_metForcing_obsFlow.zip`
2. Extract/unzip the archive into `root/Data/daily/`.
3. After extraction, clean up and rename folders to match what SAGEhydrology expects:
   - Remove (delete) the folder `basin_dataset_public` (if present).
   - Rename `basin_dataset_public_v1p2` → `v1p2`
   - Inside `v1p2`, rename:
     - `basin_mean_forcing`      → `forcing`
     - `usgs_streamflow`         → `streamflow`
     - `basin_metadata`          → `metadata`
     - `hru_forcing`             → `forcing_hru`
     - `elev_bands_forcing`      → `forcing_elev_bands`
4. Ensure the internal folder structure matches the CAMELS v1.2 layout used by SAGEhydrology.

**Expected daily layout (example)**
- `root/Data/daily/v1p2/forcing/daymet/...`
- `root/Data/daily/v1p2/forcing/maurer/...`
- `root/Data/daily/v1p2/forcing/nldas/...`
- `root/Data/daily/v1p2/streamflow/...`
- `root/Data/daily/v1p2/metadata/...`
- `root/Data/daily/v1p2/shapefiles/...`
- `root/Data/daily/v1p2/forcing_hru/...` (if present)
- `root/Data/daily/v1p2/forcing_elev_bands/...` (if present)

**References**
- Dataset (recommended citation):  
  Newman, A.; Sampson, K.; Clark, M. P.; Bock, A.; Viger, R. J.; Blodgett, D. (2014). *A large-sample watershed-scale hydrometeorological dataset for the contiguous USA.* Boulder, CO: UCAR/NCAR. https://dx.doi.org/10.5065/D6MW2F4D
- Description paper:  
  Newman, A. J., Clark, M. P., Sampson, K., Wood, A., Hay, L. E., Bock, A., Viger, R. J., Blodgett, D., Brekke, L., Arnold, J. R., Hopson, T., & Duan, Q. (2015). *Development of a large-sample watershed-scale hydrometeorological dataset for the contiguous USA: dataset characteristics and assessment of regional variability in hydrologic model performance.* *Hydrol. Earth Syst. Sci.*, 19, 209–223. https://doi.org/10.5194/hess-19-209-2015

---

### 7.2 CAMELS-US Hourly NLDAS forcing + hourly USGS streamflow

SAGEhydrology supports an hourly CAMELS-US dataset distributed by the NeuralHydrology project.

**Official sources**
- Dataset description page (NeuralHydrology):  
  https://neuralhydrology.github.io/post/datasets/camels-us-hourly-nldas-and-streamflow/
- CAMELS-US Hourly data (Zenodo record / direct download):  
  https://zenodo.org/records/4072701
- Related models and predictions archive (Zenodo record):  
  https://zenodo.org/records/4095485

**What to do**
0. Make sure you have at least 100 GB of free space on your computer.
1. Download the hourly CAMELS-US dataset from the Zenodo record above. In particular, download:
   - `nldas_hourly_csv.tar.gz`
   - `usgs-streamflow_csv.tar.gz`
2. Extract both archives into `root/Data/hourly/`.
3. After extraction, clean up and rename folders to match what SAGEhydrology expects:
   - Rename `nldas_hourly` → `forcing`
   - Rename `usgs_streamflow` → `streamflow`
4. Ensure the internal folder structure matches the layout used by SAGEhydrology.

**Expected hourly layout (example)**
- `root/Data/hourly/forcing/01013500_hourly_nldas.csv`
- `root/Data/hourly/forcing/14400000_hourly_nldas.csv`
- `root/Data/hourly/streamflow/01022500-usgs-hourly.csv`
- `root/Data/hourly/streamflow/14400000-usgs-hourly.csv`

**References**
- Dataset (recommended citation):  
  Gauch, M., Kratzert, F., Klotz, D., Nearing, G., Lin, J., & Hochreiter, S. (2020). *Data for "Rainfall-Runoff Prediction at Multiple Timescales with a Single Long Short-Term Memory Network"* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.4072701
- Related software/artifacts (models and predictions):  
  Gauch, M., Kratzert, F., Klotz, D., Nearing, G., Lin, J., & Hochreiter, S. (2020). *Models and Predictions for "Rainfall-Runoff Prediction at Multiple Timescales with a Single Long Short-Term Memory Network"*. Zenodo. https://zenodo.org/records/4095485

---

## 7.3 Notes on disk space and downloads

- The daily and hourly CAMELS-US datasets are large; ensure you have sufficient free disk space before extracting.
- If you store `Data/` on a separate drive, that is fine, provided the GUI/configuration points to the correct root directory layout.

---

## 8. Contact and Issue Reporting

For bugs, feature requests, or questions, please contact:
- Jasper A. Vrugt, University of California Irvine, jasper@uci.edu
- Jonathan M. Frame, University of Alabama, Tuscaloosa, jmframe@ua.edu

## 9. Related papers

**Paper 1**: https://arxiv.org/pdf/2602.06429<br>
**Title**: Reclaiming First Principles: A Differentiable Framework for Conceptual Hydrologic Models<br>
**Authors**: Jasper A. Vrugt, Jonathan M. Frame and E. Bollman<br>
**Addresses**: Dept. of Civil and Environmental Engineering, UC Irvine, USA, E-mail: jasper@uci.edu<br>
           Dept. of Geological Sciences, University of Alabama, Tuscaloosa, USA, E-mail: jmframe@ua.edu<br>
           
**Abstract**: Conceptual hydrologic models remain the cornerstone of rainfall-runoff modeling, yet their calibration is often slow and numerically fragile. Most gradient-based parameter estimation methods rely on finite-difference approximations or automatic differentiation frameworks (e.g., JAX, PyTorch and TensorFlow), which are computationally demanding and introduce truncation errors, solver instabilities, and substantial overhead. These limitations are particularly acute for the ODE systems of conceptual watershed models. Here we introduce a fully analytic and computationally efficient framework for differentiable hydrologic modeling based on exact parameter sensitivities. By augmenting the governing ODE system with sensitivity equations, we jointly evolve the model states and the Jacobian matrix with respect to all parameters. This Jacobian then provides fully analytic gradient vectors for any differentiable loss function. These include classical objective functions such as the sum of absolute and squared residuals, widely used hydrologic performance metrics such as the Nash-Sutcliffe and Kling-Gupta efficiencies, robust loss functions that down-weight extreme events, and hydrograph-based functionals such as flow-duration and recession curves. The analytic sensitivities eliminate the step-size dependence and noise inherent to numerical differentiation, while avoiding the instability of adjoint methods and the overhead of modern machine-learning autodiff toolchains. The resulting gradients are deterministic, physically interpretable, and straightforward to embed in gradient-based optimizers. Overall, this work enables rapid, stable, and transparent gradient-based calibration of conceptual hydrologic models, unlocking the full potential of differentiable modeling without reliance on external, opaque, or CPU-intensive automatic-differentiation libraries. <br>

**Paper 2**: https://egusphere.copernicus.org/preprints/2026/egusphere-2026-693/egusphere-2026-693.pdf <br>
**Title**: CONUS Hydrologic Modeling Using Analytic Gradients <br>
**Authors**: Jasper A. Vrugt and Jonathan Frame<br>
**Addresses**: Dept. of Civil and Environmental Engineering, UC Irvine, USA, E-mail: jasper@uci.edu<br>
           Dept. of Geological Sciences, University of Alabama, Tuscaloosa, USA, E-mail: jmframe@ua.edu<br>
           
Abstract: We introduce SAGE (Sensitivity-Aware Gradient Estimation), a new framework for scalable and physics-consistent training of hydrologic models that leverages analytic forward sensitivities to enable exact and efficient gradient-based learning of model parameters from catchment attributes. Unlike existing approaches that rely on finite-difference approximations, automatic differentiation, or surrogate emulators, SAGE propagates exact derivatives through physically based dynamical systems using analytically derived sensitivity equations. This eliminates the need for repeated model evaluations, substantially reduces computational cost, and preserves the interpretability and structural integrity of process-based hydrologic models. We demonstrate SAGE in a large-sample hydrology experiment using the CAMELS data set, comprising 531 hydrologically valid catchments across the contiguous United States. A feedforward neural network maps static catchment attributes to the parameter space of a conceptual rainfall-runoff model, while exact gradients of the loss function with respect to network weights are computed through analytic sensitivity propagation of the governing ordinary differential equations. Compared to conventional training strategies based on numerical differentiation or automatic differentiation, SAGE achieves machine-precision agreement with reference gradients while reducing computational cost by several orders of magnitude. To assess cross-basin model performance, we further introduce a new integrated distributional skill score based on the empirical cumulative distribution function of Nash-Sutcliffe efficiency (NSE) values across basins. Rather than summarizing performance using a single quantile such as the median NSE, the proposed score quantifies the distance between the observed basin-wise NSE distribution and the ideal degenerate distribution at NSE = 1. This distributional skill score provides a more robust and informative measure of large-sample model skill and enables objective comparison of learning strategies at continental scale. Together, SAGE and the proposed Vrugt-Frame loss score form a unified framework for both training and evaluating physics-based hydrologic models in large-sample settings and offer a new pathway toward continental-scale, attribute-conditioned calibration that is both computationally tractable and physically interpretable.<br>

**Paper 3**: Continental-Scale Hydrologic Model Training at Hourly Resolution with Sensitivity-Aware Gradient Estimation <br>
**Title**: Jasper A. Vrugt, Jonathan Frame and Yifu Gao <br>
**Addresses**: Dept. of Civil and Environmental Engineering, UC Irvine, USA, E-mail: jasper@uci.edu <br>
           Dept. of Geological Sciences, University of Alabama, Tuscaloosa, USA, E-mail: jmframe@ua.edu <br>
           School of Atmospheric Sciences, Nanjing University of Information Science and Technology, Nanjing, Jiangsu, China, E-mail: gaoyifu@nuist.edu.cn <br>

Abstract: We implement SAGE (Sensitivity-Aware Gradient Estimation), a framework for scalable and physics-consistent training of hydrologic models, to the demanding setting of hourly large-sample hydrology. Training process-based rainfall--runoff models at hourly resolution remains computationally prohibitive for many automatic differentiation (AD) workflows due to the length of multi-year trajectories and the associated memory and computational overhead of differentiating through dynamical systems. SAGE circumvents these limitations by propagating \emph{analytic forward sensitivities} through the governing ordinary differential equations, enabling exact gradients of the loss function with respect to both model parameters and attribute-to-parameter mapping weights without repeated model evaluations, numerical differencing, or reverse-mode AD through the full time history. We demonstrate SAGE using approximately six years of hourly forcing and discharge data per basin (including one year of spin-up, three years of training, and two years of validation) for roughly 500 CAMELS catchments across the contiguous United States. On an older desktop CPU (Intel(R) Core(TM) i7-10700T @ 2.00\,GHz; 8 cores), a single SAGE iteration over all 500 basins requires only about 50-60 seconds, establishing that exact, physics-consistent gradient-based learning at continental scale is computationally tractable at hourly resolution on commodity hardware. We benchmark SAGE against single-site training experiments for six conceptual hydrologic models including \texttt{hymod}, \texttt{hmodel}, \texttt{sacsma}, \texttt{Xinanjiang}, \texttt{gr4j}, and \texttt{hbv}, and assess generalization under temporal validation, spatial validation, and spatiotemporal validation using conventional goodness-of-fit metrics and the recently proposed Vrugt-Frame interbasin skill score. Overall, hourly-resolution skill is qualitatively consistent with daily-resolution findings, with basin-wise performance distributions shifted downward by roughly 0.1 in Nash-Sutcliffe efficiency (NSE). Among the tested models, Xinanjiang achieves the strongest predictive performance and learns the most effective attribute-parameter relationships for validation basins, whereas HBV generalizes poorly in comparison. These results position SAGE as a practical pathway toward rapid, continental-scale, attribute-conditioned calibration at hourly resolution, a regime that remains challenging for conventional AD-based differentiable modeling frameworks.
