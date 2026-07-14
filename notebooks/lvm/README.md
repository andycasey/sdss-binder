# Local Volume Mapper (LVM) Notebooks

## Getting Started

- **[hilder-demo.ipynb](hilder-demo.ipynb)**: Introduction to working with LVM data. Shows how to read `lvmSFrame` FITS files (both manually with `astropy.io.fits` and using `lvm_tools`), plot calibrated sky-subtracted spectra, make integrated flux maps, and query tiles by coordinates or MJD. By _Thomas Hilder_.

## DR20 DAP / DRP tutorials (`dr20/`)

Tutorials for working with SDSS-V LVM DR20 data products (DRP `SFrame` and DAP files). All by _Sebastian F. Sanchez and the LVM team_.

- **[dr20/lvm_dr20_SFrame_view.ipynb](dr20/lvm_dr20_SFrame_view.ipynb)**: How to open, inspect, visualize, and interactively explore an LVM DRP `SFrame` product, using the DR20 Helix Nebula (NGC 7293) exposure as an example. Covers loading via `sdss_access`/`Tree` with a local fallback, reconstructing narrow-band slices from the RSS, an interactive fiber/spectrum viewer, and exporting selected spectra.
- **[dr20/lvm_dr20_helix_nebula_DAP_tutorial.ipynb](dr20/lvm_dr20_helix_nebula_DAP_tutorial.ipynb)**: How to handle a DR20 DAP file, using the Helix Nebula (NGC 7293) as an example: inspecting the product, generating RGB and individual line maps, separating parametric and non-parametric products, and building integrated line catalogs. Recreates plots from Sánchez et al. 2025.
- **[dr20/lvm_dr20_Mosaic_DAP_tutorial.ipynb](dr20/lvm_dr20_Mosaic_DAP_tutorial.ipynb)**: How to build mosaics from several DR20 DAP files, using the M33 pointings as an example: RGB and line maps, BPT diagrams, kinematics and radial structure, and stellar-continuum products.
- **[dr20/lvm_dr20_drpall_dapall_qa_bitmask.ipynb](dr20/lvm_dr20_drpall_dapall_qa_bitmask.ipynb)**: How to read the DR20 `drpall`, `dapall`, and `qcall` summary products, merge them into a single table, and decode the packed `qa_bitmask` values into human-readable quality flags for science selection.
- **[dr20/lvm_dr20_generate_model_tutorial.ipynb](dr20/lvm_dr20_generate_model_tutorial.ipynb)**: How to reconstruct DAP output-model products (observed spectra, stellar model, non-parametric and parametric gas models, and joint model) by reading configuration parameters directly from the `INFO` extension of a DAP FITS file. A tutorial version of the `gen_output_model.py` workflow.
- **[dr20/lvm_dr20_read_model_tutorial.ipynb](dr20/lvm_dr20_read_model_tutorial.ipynb)**: How to read an LVM DR20 DAP model file (Helix Nebula example) and visualize the reconstructed spectral components: observed spectrum, stellar-continuum model, gas-emission models, joint models, and residuals.
- **[dr20/lvm_dr20_read_model_tutorial_documented_v2.ipynb](dr20/lvm_dr20_read_model_tutorial_documented_v2.ipynb)**: A more heavily documented version of the DR20 model-file reading tutorial, explaining the FITS structure and each model component (stellar continuum, parametric and non-parametric gas, residuals, and combined models).

## `spectracles` (spectrospatial modelling)

Notebooks using the [`spectracles`](https://github.com/thomashilder/spectracles) library for joint spectral-spatial modelling of emission lines in LVM IFU data. `spectracles` models spectral quantities (e.g., line flux, radial velocity, velocity dispersion) as continuous fields on the sky using accelerated Gaussian Processes and JAX. See Hilder+2026 for details.

### Single emission line tutorial (`spectracles_single_line/`)

A self-contained tutorial for fitting a single emission line model to a small LVM dataset (3 tiles covering the Flame Nebula). The line to fit is configured in `config.py` — just change the `LINE` variable to switch between H-alpha, [N II] 6583, or [O III] 5007.

- **[1_fit_model.ipynb](spectracles_single_line/1_fit_model.ipynb)**: Loads data, builds a single-line `spectracles` model with three GP fields (amplitude, velocity, velocity dispersion) plus nuisance parameters (continuum, flux calibration, wavelength calibration), defines a block coordinate descent schedule, and runs the optimisation. Saves the fitted model to disk. By _Thomas Hilder_.
- **[2_plot_results.ipynb](spectracles_single_line/2_plot_results.ipynb)**: Loads the fitted model and produces science maps (amplitude, velocity, velocity dispersion), calibration diagnostics (flux cal, wavelength cal), and fit quality plots (spectra overlays, reduced chi-squared maps, residual histograms, worst-fit spaxels). By _Thomas Hilder_.

### W28 (`spectracles_W28/`)

Spatial-kinematic decomposition of the supernova remnant W28, using two-component emission line models to separate line-of-sight velocity components. Each notebook fits a different line, and results from earlier fits inform later ones (e.g., star masks from [N II] are used in subsequent fits).

- **[1_nii.ipynb](spectracles_W28/1_nii.ipynb)**: Fits [N II] 6583 with a two-component model. This is done first because [N II] is unaffected by stellar line contamination. The per-spaxel continuum offsets from this fit are used to identify and mask star locations for subsequent lines. By _Thomas Hilder_.
- **[2_ha.ipynb](spectracles_W28/2_ha.ipynb)**: Fits H-alpha with star masking derived from the [N II] fit. By _Thomas Hilder_.
- **[3_oiii.ipynb](spectracles_W28/3_oiii.ipynb)**: Fits [O III] 5007. By _Thomas Hilder_.
- **[4_sii.ipynb](spectracles_W28/4_sii.ipynb)**: Fits [S II] 6716. By _Thomas Hilder_.
