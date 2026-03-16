# Local Volume Mapper (LVM) Notebooks

## Getting Started

- **[hilder-demo.ipynb](hilder-demo.ipynb)**: Introduction to working with LVM data. Shows how to read `lvmSFrame` FITS files (both manually with `astropy.io.fits` and using `lvm_tools`), plot calibrated sky-subtracted spectra, make integrated flux maps, and query tiles by coordinates or MJD. By _Thomas Hilder_.

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
