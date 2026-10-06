# SAGEhydrology v1.0.3

This Windows release retains the public capabilities of v1.0.2 and improves standalone deployment, reproducible projects, custom-model loading, and the presentation of training results. It provides seven built-in models and the manual `user_model` interface. AI-assisted model authoring and private model implementations are excluded.

## Changes from v1.0.2

- **Projects:** the compiled GUI includes publication Projects 1–5 and the new-user-project workflow. “Open source” opens readable MATLAB exports; scripts for Projects 1–2 are included rather than pointing to missing files. Project settings are locked while a SAGE run is active.
- **Custom C++ models:** the manual `user_model` option is available independently of the AI builder. A sibling `user_model` folder supplies parameter metadata and a compatible precompiled MEX. Runtime loading no longer calls `addpath`, and the selected external binary is installed into Runtime's resolved MEX location before evaluation. Users compile replacement C++ models externally, update `user_model_info.mat`, and restart SAGE before using the replacement.
- **Paths and deployment:** Software root can be selected in the executable. Startup recognizes an installation beside Data and SAGEhydrology instead of treating the Runtime cache as the user's software installation. Asset lookup avoids a mandatory call to MATLAB's unavailable internal executable-path API.
- **Installed datasets and maps:** Russia/HydroCIS and daily CAMELS-KR checks resolve their bundled basin inventories in Runtime. Existing 10 m map files under `SAGEhydrology/maps` are detected using the selected software root. The GUI retains the compact 50 m map for responsive interactive rendering.
- **Training preparation and logs:** discharge-statistics progress updates are throttled, and model-specific Kosugi requirements are reset when changing models. Detailed log messages use consistent indentation and wrapping.
- **Run summary and controls:** model and network parameter counts are shown; the two summary columns have improved spacing, and text/control heights avoid clipped letter descenders. Log buttons use improved vertical spacing.
- **Figures:** more consistent ECDF fonts/layouts and tick lengths; shorter hydrograph ticks; readable time-resolution badges; upright variogram panel labels with tighter spacing. Parameter maps are paginated at 20 parameters per figure, with a scale in each figure's top-left panel and its own colorbar.
- **Source examples:** exactly three examples are supplied: `demo_SAGE.mlx`, `verify_jacobians_and_gradients.m`, and `run_SAGE_export_hymod.m`. Project 3–5 exports reside in their project folders. Example paths derive the installation root from their file location.

## Which Windows download should I use?

| Download | Purpose |
| --- | --- |
| `SAGE-v1.0.3-Windows-x64-installer.exe` | Guided installation; obtains MATLAB Runtime R2026a through web delivery when required. |
| `SAGE-v1.0.3-Windows-x64-portable.zip` | Recommended for an existing Software/Data installation. Includes SAGE.exe, public supporting files, maps, projects, and the manual user-model template/kernel. Requires MATLAB Runtime R2026a. |
| `SAGE-v1.0.3-Windows-x64.exe` | The application executable alone, for replacing an existing installation's executable. Requires MATLAB Runtime R2026a and separately maintained Data/user-model files. |

The installer and portable executable deliver the same v1.0.3 GUI and model capabilities. The difference is installation and Runtime setup, not functionality.

Keep `SAGE.exe`, `Data`, `SAGEhydrology`, and `user_model` beside one another in your Software folder. In Paths, select that Software root and your Data/results locations. Existing datasets do not need to be downloaded again. Run results and large regional datasets are not included in these downloads.

Project exports run in MATLAB with the computational source installed; MATLAB Runtime cannot execute arbitrary `.m` scripts. Compiled GUI use is covered by `GUI-LICENSE-NOTICE.md`; the public computational source has the BSD 3-Clause License. Readable GUI source, private models, AI-builder files, results, and development caches are excluded from public assets.

The public repository directories are updated for v1.0.3, including src, models, projects, regions, utils, the manual user-model template, and the three-file examples directory. Current public source is available on the main branch.

## Platform and validation

Windows x64 application and installer packaging passed verification, including nine native MEX inputs. A compiled Runtime test loaded an independently rebuilt external user-model MEX and produced finite discharge and Jacobian values. Source/dependency and GUI startup checks passed. Full training through the final GUI and installation on a clean Windows machine have not been independently verified in this release preparation.

macOS v1.0.3 builds for Apple Silicon are now available in this release.

## macOS downloads available

SAGE v1.0.3 now includes macOS builds for **Apple Silicon (arm64)** alongside the existing Windows downloads.

- `SAGE-v1.0.3-macOS-arm64-installer.zip`: unzip and run the installer application; MATLAB Runtime R2026a is obtained when required.
- `SAGE-v1.0.3-macOS-arm64.dmg`: disk image containing the application; requires the Apple Silicon MATLAB Runtime R2026a.

Download both options from [the v1.0.3 release](https://github.com/jaspervrugt/SAGEhydrology/releases/tag/v1.0.3). These packages target Apple Silicon; they are not Intel Mac builds.

The supplied macOS build report records arm64 launcher, native MEX, ZIP and DMG checks. Downloaded file hashes were verified against that report. Windows downloads are unchanged.
