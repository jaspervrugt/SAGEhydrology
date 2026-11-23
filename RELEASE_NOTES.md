# SAGEhydrology v2.0.0

- Adds headless Python source with JSON session execution, a synthetic
  training example, reader metadata and regression tests.
- Updates the MATLAB training, postprocessing, regional readers and model
  source. The repository now organizes the language sources in `Matlab/`
  and `Python/`.
- Includes the shared C++ backend, MATLAB MEX gateways, Python bindings, and
  CaRT and CaRT+ solver implementations.
- Includes an editable GR4JB example in `user_model/` for custom models in
  MATLAB, headless Python and the desktop application.
- Distributes the desktop GUI as a compiled application. The Python Windows
  bundle uses a standalone runtime. macOS availability will be announced
  after the build is verified on a Mac.

## Python desktop 1.0

The Windows GUI is a separate portable download. Unzip the complete folder and open SAGEhydrology.exe. A macOS download can be added after a separate build and verification on a Mac.

The log-quantile diagnostic handles near-zero negative discharge within three times the solver absolute tolerance. The training objective and saved discharge series retain their original definitions.
