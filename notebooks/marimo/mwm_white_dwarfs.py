import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # A JAX-based spectrum fitter for DA-type White Dwarfs
    For experimentation and understanding.

    Andy Casey (acasey@flatironinstitute.org)
    """)
    return


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _():
    import jax

    jax.config.update("jax_enable_x64", True)

    import jax.numpy as jnp
    import interpax
    import h5py as h5
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    C_KMS = 299792.458  # speed of light [km/s]


    # ---------------------------------------------------------------------------
    # Data loading
    # ---------------------------------------------------------------------------

    def load_grid(filepath="wd_da_spectra.h5"):
        """Load the pre-computed WD DA spectral grid from HDF5.

        Returns dict with keys teff_grid (65,), logg_grid (36,),
        wavelength (4299,), flux (65, 36, 4299).
        """
        with h5.File(filepath, "r") as f:
            # .astype ensures native byte order (HDF5 may store big-endian)
            return dict(
                teff_grid=jnp.asarray(np.asarray(f["Teff"]).astype(np.float64)),
                logg_grid=jnp.asarray(np.asarray(f["logg"]).astype(np.float64)),
                wavelength=jnp.asarray(np.asarray(f["wavelength"]).astype(np.float64)),
                flux=jnp.asarray(np.asarray(f["flux"]).astype(np.float64)),
            )


    # ---------------------------------------------------------------------------
    # Bicubic interpolation
    # ---------------------------------------------------------------------------

    def interpolate_spectrum(teff, logg, teff_grid, logg_grid, flux_grid):
        """Bicubic interpolation of the spectral grid at a single (Teff, logg).

        Uses interpax.interp2d with method='cubic' on the (Teff, logg) axes.
        The flux grid has shape (n_teff, n_logg, n_wave); the trailing wavelength
        dimension is interpolated in batch, returning a (n_wave,) spectrum.
        """
        return interpax.interp2d(
            jnp.atleast_1d(teff),
            jnp.atleast_1d(logg),
            teff_grid,
            logg_grid,
            flux_grid,
            method="cubic",
            extrap=True,
        )[0]  # squeeze out the Nq=1 leading dimension


    # ---------------------------------------------------------------------------
    # Radial-velocity shift
    # ---------------------------------------------------------------------------

    def apply_rv_and_resample(template, model_wavelength, obs_wavelength, v_rad):
        """Doppler-shift *template* and resample onto *obs_wavelength*.

        For each observed pixel the rest-frame wavelength is
        lambda_rest = lambda_obs / (1 + v/c).  The template (defined on
        *model_wavelength*) is evaluated at these rest wavelengths via cubic
        spline interpolation.
        """
        wave_rest = obs_wavelength / (1.0 + v_rad / C_KMS)
        return interpax.interp1d(
            wave_rest,
            model_wavelength,
            template,
            method="cubic",
            extrap=True,
        )


    # ---------------------------------------------------------------------------
    # Fourier continuum
    # ---------------------------------------------------------------------------

    def build_continuum_design_matrix(wavelength, n_coeffs):
        """Fourier basis: {1, cos 2pi*x, sin 2pi*x, cos 4pi*x, sin 4pi*x, ...}

        x = (lambda - lambda_min) / (lambda_max - lambda_min)  in [0, 1] so that
        harmonic *k* completes *k* full cycles over the fitting range.

        Returns an (n_wave, n_coeffs) design matrix.
        """
        x = (wavelength - wavelength[0]) / (wavelength[-1] - wavelength[0])
        basis = [jnp.ones_like(x)]
        k = 1
        while len(basis) < n_coeffs:
            basis.append(jnp.cos(2.0 * jnp.pi * k * x))
            if len(basis) < n_coeffs:
                basis.append(jnp.sin(2.0 * jnp.pi * k * x))
            k += 1
        return jnp.stack(basis[:n_coeffs], axis=-1)


    # ---------------------------------------------------------------------------
    # Forward model
    # ---------------------------------------------------------------------------

    class WhiteDwarfModel:
        """Forward model for fitting white dwarf DA spectra.

        Parameters
        ----------
        grid_path : str
            Path to ``wd_da_spectra.h5``.
        n_continuum : int
            Number of Fourier continuum coefficients.
        obs_wavelength : array-like or None
            Observed wavelength grid [Angstrom].  Defaults to the full model grid.

        The parameter vector layout is::

            params = [teff, logg, v_rad, c_0, c_1, ..., c_{n_continuum-1}]

        where *c_0* is the DC (constant) amplitude and higher indices alternate
        cosine / sine harmonics.

        All returned functions are JIT-compiled and differentiable.
        """

        def __init__(self, grid_path="wd_da_spectra.h5", n_continuum=8,
                     obs_wavelength=None):
            grid = load_grid(grid_path)
            self.teff_grid = grid["teff_grid"]
            self.logg_grid = grid["logg_grid"]
            self.flux_grid = grid["flux"]
            self.model_wavelength = grid["wavelength"]
            self.n_continuum = n_continuum
            self.n_params = 3 + n_continuum

            if obs_wavelength is not None:
                obs_wavelength = np.asarray(obs_wavelength, dtype=np.float64)
                # Clip to wavelengths covered by the model grid
                wmin = float(self.model_wavelength[0])
                wmax = float(self.model_wavelength[-1])
                self.pixel_mask = (obs_wavelength >= wmin) & (obs_wavelength <= wmax)
                self.obs_wavelength = jnp.asarray(obs_wavelength[self.pixel_mask])
            else:
                self.pixel_mask = None
                self.obs_wavelength = self.model_wavelength

            self.design_matrix = build_continuum_design_matrix(
                self.obs_wavelength, n_continuum
            )

            self.param_names = ["teff", "logg", "v_rad"] + [
                f"c_{i}" for i in range(n_continuum)
            ]

            # ---- build JIT-compiled closures once ----
            # All functions are pure closures over the grid arrays so they
            # compose freely with jax.jit / vmap / grad / lax.scan.
            tg = self.teff_grid
            lg = self.logg_grid
            fg = self.flux_grid
            mw = self.model_wavelength
            ow = self.obs_wavelength
            dm = self.design_matrix

            def _forward(params):
                teff, logg, v_rad = params[0], params[1], params[2]
                cont_coeffs = params[3:]
                template = interpolate_spectrum(teff, logg, tg, lg, fg)
                shifted = apply_rv_and_resample(template, mw, ow, v_rad)
                continuum = dm @ cont_coeffs
                return shifted * continuum

            def _weights(flux_err):
                """1/sigma weights, zeroed for masked pixels (sigma >= 1e5)."""
                return jnp.where(flux_err < 1e5, 1.0 / flux_err, 0.0)

            def _loss(params, flux_obs, flux_err):
                model = _forward(params)
                w = _weights(flux_err)
                return jnp.sum(((flux_obs - model) * w) ** 2)

            def _solve_continuum(phys, flux_obs, flux_err):
                teff, logg, v_rad = phys[0], phys[1], phys[2]
                template = interpolate_spectrum(teff, logg, tg, lg, fg)
                shifted = apply_rv_and_resample(template, mw, ow, v_rad)
                w = _weights(flux_err)
                A = (shifted[:, None] * dm) * w[:, None]
                b = flux_obs * w
                return jnp.linalg.lstsq(A, b, rcond=None)[0]

            def _profile_loss(phys, flux_obs, flux_err):
                c = _solve_continuum(phys, flux_obs, flux_err)
                params = jnp.concatenate([phys, c])
                model = _forward(params)
                w = _weights(flux_err)
                return jnp.sum(((flux_obs - model) * w) ** 2)

            def _fit(flux_obs, flux_err, lower, upper,
                     teff_search, logg_search, vrad_search):
                """Vmapped grid search + Levenberg-Marquardt polish."""
                # Stage 1: vectorised grid search with profiled chi-squared
                te_g, lg_g, vr_g = jnp.meshgrid(
                    teff_search, logg_search, vrad_search, indexing="ij")
                grid_phys = jnp.stack(
                    [te_g.ravel(), lg_g.ravel(), vr_g.ravel()], axis=-1)
                chi2_vals = jax.vmap(
                    lambda p: _profile_loss(p, flux_obs, flux_err)
                )(grid_phys)
                phys_best = grid_phys[jnp.argmin(chi2_vals)]
                c_best = _solve_continuum(phys_best, flux_obs, flux_err)
                p0 = jnp.concatenate([phys_best, c_best])

                # Stage 2: Levenberg-Marquardt (Gauss-Newton + adaptive damping)
                w = _weights(flux_err)

                def lm_step(carry, _):
                    params, damping = carry
                    model_flux = _forward(params)
                    r = (flux_obs - model_flux) * w
                    chi2_cur = jnp.sum(r ** 2)
                    J = jax.jacfwd(_forward)(params)
                    Jw = J * w[:, None]
                    JtJ = Jw.T @ Jw
                    Jtr = Jw.T @ r
                    dp = jnp.linalg.solve(
                        JtJ + damping * jnp.diag(jnp.diag(JtJ)), Jtr)
                    p_new = jnp.clip(params + dp, lower, upper)
                    chi2_new = _loss(p_new, flux_obs, flux_err)
                    improved = chi2_new < chi2_cur
                    p_out = jnp.where(improved, p_new, params)
                    d_out = jnp.where(improved, damping * 0.1, damping * 10.0)
                    d_out = jnp.clip(d_out, 1e-12, 1e8)
                    return (p_out, d_out), chi2_cur

                (p_final, _), _ = jax.lax.scan(
                    lm_step, (p0, 1e-3), None, length=30)
                return p_final, p0

            self._forward_jit = jax.jit(_forward)
            self._loss_jit = jax.jit(_loss)
            self._loss_and_grad_jit = jax.jit(jax.value_and_grad(_loss))
            self._hessian_jit = jax.jit(jax.hessian(_loss))
            self._solve_continuum_jit = jax.jit(_solve_continuum)
            self._profile_loss_jit = jax.jit(_profile_loss)
            self._fit_jit = jax.jit(_fit)
            # Expose raw closures for external JIT/vmap/grad composition
            self._forward_raw = _forward
            self._solve_continuum = _solve_continuum
            self._profile_loss = _profile_loss

        # ---- public API ---------------------------------------------------------

        def __call__(self, params):
            """Evaluate the forward model.  Returns model flux array."""
            return self._forward_jit(jnp.asarray(params, dtype=jnp.float64))

        def loss(self, params, flux_obs, flux_err):
            """chi^2 = sum[(data - model)^2 / sigma^2]."""
            return self._loss_jit(
                jnp.asarray(params, dtype=jnp.float64),
                jnp.asarray(flux_obs, dtype=jnp.float64),
                jnp.asarray(flux_err, dtype=jnp.float64),
            )

        def loss_and_grad(self, params, flux_obs, flux_err):
            """Return (chi^2, d(chi^2)/d(params)) via reverse-mode AD."""
            return self._loss_and_grad_jit(
                jnp.asarray(params, dtype=jnp.float64),
                jnp.asarray(flux_obs, dtype=jnp.float64),
                jnp.asarray(flux_err, dtype=jnp.float64),
            )

        def hessian(self, params, flux_obs, flux_err):
            """Full Hessian of chi^2 w.r.t. params (forward-over-reverse AD)."""
            return self._hessian_jit(
                jnp.asarray(params, dtype=jnp.float64),
                jnp.asarray(flux_obs, dtype=jnp.float64),
                jnp.asarray(flux_err, dtype=jnp.float64),
            )

        def initial_params(self, teff=10000.0, logg=8.0, v_rad=0.0):
            """Sensible starting vector (flat continuum at unit amplitude)."""
            cont = jnp.zeros(self.n_continuum).at[0].set(1.0)
            return jnp.concatenate([jnp.array([teff, logg, v_rad]), cont])

        def _apply_mask(self, flux_obs, flux_err):
            """Apply pixel_mask (if set) and return JAX arrays."""
            flux_obs = np.asarray(flux_obs, dtype=np.float64)
            flux_err = np.asarray(flux_err, dtype=np.float64)
            if self.pixel_mask is not None:
                flux_obs = flux_obs[self.pixel_mask]
                flux_err = flux_err[self.pixel_mask]
            return jnp.asarray(flux_obs), jnp.asarray(flux_err)

        def fit(self, flux_obs, flux_err, p0=None, bounds=None):
            """JIT-compiled two-stage fit: vmapped grid search + Levenberg-Marquardt.

            Stage 1: Vectorised grid search over (teff, logg, v_rad) using
            ``jax.vmap``.  At each trial point the continuum coefficients are
            solved analytically, so only 3 parameters are searched.

            Stage 2: Levenberg-Marquardt iterations over all parameters using
            JAX forward-mode Jacobians, compiled via ``jax.lax.scan``.

            The entire pipeline is JIT-compiled on first call.

            Parameters
            ----------
            flux_obs, flux_err : array (n_obs,)
                Full observed arrays — pixels outside the model wavelength range
                are automatically masked out.
            p0 : array (n_params,), optional
                Not used (kept for API compatibility).
            bounds : list of (lo, hi), optional
                Bounds for (teff, logg, v_rad).  Defaults to interior of the
                grid for Teff/logg, +/-500 km/s for v_rad.

            Returns
            -------
            result : SimpleNamespace
                ``.x`` best-fit params, ``.x_de`` grid-search init,
                ``.fun`` final chi-squared, ``.success`` always True.
            """
            from types import SimpleNamespace

            flux_obs_j, flux_err_j = self._apply_mask(flux_obs, flux_err)

            teff_lo, teff_hi = float(self.teff_grid[1]), float(self.teff_grid[-2])
            logg_lo, logg_hi = float(self.logg_grid[1]), float(self.logg_grid[-2])

            if bounds is None:
                phys_bounds = [
                    (teff_lo, teff_hi),
                    (logg_lo, logg_hi),
                    (-500.0, 500.0),
                ]
            else:
                phys_bounds = list(bounds[:3])

            lower = jnp.array(
                [phys_bounds[0][0], phys_bounds[1][0], phys_bounds[2][0]]
                + [-jnp.inf] * self.n_continuum)
            upper = jnp.array(
                [phys_bounds[0][1], phys_bounds[1][1], phys_bounds[2][1]]
                + [jnp.inf] * self.n_continuum)

            teff_search = jnp.linspace(*phys_bounds[0], 25)
            logg_search = jnp.linspace(*phys_bounds[1], 20)
            vrad_search = jnp.linspace(*phys_bounds[2], 11)

            p_final, p_init = self._fit_jit(
                flux_obs_j, flux_err_j, lower, upper,
                teff_search, logg_search, vrad_search,
            )

            chi2 = float(self._loss_jit(p_final, flux_obs_j, flux_err_j))
            return SimpleNamespace(
                x=np.asarray(p_final),
                x_de=np.asarray(p_init),
                fun=chi2,
                success=True,
            )

        def fit_with_uncertainties(self, flux_obs, flux_err, p0=None, bounds=None):
            """Fit + Hessian-based parameter uncertainties.

            Returns (OptimizeResult, covariance, uncertainties).
            For a chi^2 loss the covariance is 2 * H^{-1} at the minimum.
            """
            result = self.fit(flux_obs, flux_err, p0, bounds)
            p_best = jnp.asarray(result.x, dtype=jnp.float64)
            flux_obs_j, flux_err_j = self._apply_mask(flux_obs, flux_err)

            H = self.hessian(p_best, flux_obs_j, flux_err_j)
            cov = 2.0 * jnp.linalg.inv(H)
            unc = jnp.sqrt(jnp.abs(jnp.diag(cov)))
            return result, np.asarray(cov), np.asarray(unc)




    def make_plot(model, flux_obs, flux_err, result, cov, unc, landscape=True, n_teff_grid=30, n_logg_grid=25, n_zoom=40):


        wave = np.asarray(model.obs_wavelength)
        # Use the same masked arrays the fitter saw
        flux_obs_np, flux_err_np = (np.asarray(a) for a in model._apply_mask(flux_obs, flux_err))
        flux_fit = np.asarray(model(jnp.asarray(result.x)))

        # Propagate parameter covariance to model flux uncertainty:
        #   sigma_model^2 = J @ cov @ J^T  (diagonal entries)
        # where J = d(model) / d(params) is the (n_wave, n_params) Jacobian.
        J = np.asarray(jax.jacobian(model._forward_jit)(jnp.asarray(result.x)))
        model_var = np.einsum("ij,jk,ik->i", J, cov, J)
        model_sig = np.sqrt(np.abs(model_var))

        # --- Compute chi^2 landscape over (Teff, logg) ---
        if landscape:
            print("Computing chi^2 landscape ...")
            flux_obs_j, flux_err_j = model._apply_mask(flux_obs, flux_err)
            profile_fn = jax.jit(model._profile_loss)
            # Fix v_rad at the best-fit value and scan Teff, logg
            v_rad_best = float(result.x[2])
            teff_scan = np.linspace(float(model.teff_grid[1]), float(model.teff_grid[-2]), n_teff_grid)
            logg_scan = np.linspace(float(model.logg_grid[1]), float(model.logg_grid[-2]), n_logg_grid)
            chi2_grid = np.empty((n_logg_grid, n_teff_grid))
            for i, lg_val in enumerate(logg_scan):
                for j, te_val in enumerate(teff_scan):
                    chi2_grid[i, j] = float(
                        profile_fn(jnp.array([te_val, lg_val, v_rad_best]), flux_obs_j, flux_err_j)
                    )
            print("  Done.")

            # --- Compute zoomed chi^2 landscape around best fit ---
            teff_final, logg_final = float(result.x[0]), float(result.x[1])
            zoom_nsig = 5
            teff_zoom = np.linspace(teff_final - zoom_nsig * unc[0],
                                    teff_final + zoom_nsig * unc[0], n_zoom)
            logg_zoom = np.linspace(logg_final - zoom_nsig * unc[1],
                                    logg_final + zoom_nsig * unc[1], n_zoom)
            print("Computing zoomed chi^2 landscape ...")
            chi2_zoom = np.empty((n_zoom, n_zoom))
            from tqdm import tqdm
            for i, lg_val in enumerate(tqdm(logg_zoom)):
                for j, te_val in enumerate(teff_zoom):
                    chi2_zoom[i, j] = float(
                        profile_fn(jnp.array([te_val, lg_val, v_rad_best]),
                                flux_obs_j, flux_err_j)
                    )
            print("  Done.")

        # --- Figure ---
        fig = plt.figure(figsize=(16, 8))
        # Left 3/4: spectrum (top) + residuals (bottom) with 3:1 height ratio
        # Right 1/4: two equal panels stacked vertically
        outer_gs = fig.add_gridspec(1, 2, width_ratios=[3, 1], wspace=0.30)
        left_gs = outer_gs[0].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
        right_gs = outer_gs[1].subgridspec(2, 1, height_ratios=[1, 1], hspace=0.35)

        ax1 = fig.add_subplot(left_gs[0])              # spectrum (top-left)
        ax2 = fig.add_subplot(left_gs[1], sharex=ax1)  # residuals (bottom-left)
        ax3 = fig.add_subplot(right_gs[0])             # chi^2 landscape (top-right)
        ax4 = fig.add_subplot(right_gs[1])             # likelihood zoomed (bottom-right)

        # -- top-left panel: data + model --
        ax1.plot(wave, flux_obs_np, "k-", lw=0.4, label="data", drawstyle="steps-mid")
        ax1.fill_between(wave, flux_fit - model_sig, flux_fit + model_sig,
                        color="red", alpha=0.25, label="model $\\pm 1\\sigma$")
        ax1.plot(wave, flux_fit, "r-", lw=0.6, label="model fit")
        ylim = ax1.get_ylim()

        ax1.fill_between(wave, flux_obs_np - flux_err_np, flux_obs_np + flux_err_np,
                        color="grey", alpha=0.4, label="data $\\pm 1\\sigma$", zorder=-10)
        ax1.set_ylim(ylim)
        ax1.set_ylabel("Flux [$\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{\AA}^{-1}$]")
        ax1.set_title(
            f"Teff = {result.x[0]:.0f} K +/- {unc[0]:.0f},  "
            f"log g = {result.x[1]:.2f} +/- {unc[1]:.2f},  "
            f"v_rad = {result.x[2]:.1f} km/s +/- {unc[2]:.1f}"
        )

        # -- bottom-left panel: residuals --
        resid = (flux_obs_np - flux_fit) / flux_err_np
        model_resid_sig = model_sig / flux_err_np  # model uncertainty in sigma units
        ax2.fill_between(wave, -1, 1, color="grey", alpha=0.2)
        ax2.fill_between(wave, -model_resid_sig, model_resid_sig,
                        color="red", alpha=0.3, label="model $\\pm 1\\sigma$")
        ax2.plot(wave, resid, "k-", lw=0.4, drawstyle="steps-mid")
        ax2.axhline(0, color="red", lw=0.5)
        ax2.set_ylabel("Residual ($\\sigma$)")
        ax2.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
        ax2.set_ylim(-3, 3)

        for ax in (ax1, ax2):
            ax.set_xlim(wave[0], wave[-1])
        plt.setp(ax1.get_xticklabels(), visible=False)

        # -- Shared: covariance ellipse parameters for Teff-logg --
        from matplotlib.patches import Ellipse
        cov_teff_logg = cov[:2, :2]
        eigvals, eigvecs = np.linalg.eigh(cov_teff_logg)
        order = eigvals.argsort()[::-1]
        eigvals, eigvecs = eigvals[order], eigvecs[:, order]
        angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))

        if landscape:

            # -- top-right panel: chi^2 landscape --
            chi2_min = chi2_grid.min()
            ln_chi2 = np.log(chi2_grid)
            ln_min = np.log(chi2_min)
            levels = np.linspace(ln_min, ln_min + 5, 30)
            cf = ax3.contourf(teff_scan, logg_scan, ln_chi2,
                            levels=levels, cmap="Blues_r", extend="max")
            cbar = fig.colorbar(cf, ax=ax3, orientation="horizontal", location="top", pad=0.03)
            cbar.set_label("$\\ln\\,\\chi^2$", fontsize=9)
            cbar.ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))

            # Mark DE starting point and final result
            teff_de, logg_de = result.x_de[0], result.x_de[1]

            ax3.plot(teff_de, logg_de, "o", color="none", ms=9, mec="k", mew=1.0, zorder=5)
            ax3.plot(teff_final, logg_final, "o", color="red", ms=9, mec="k", mew=0.8, zorder=5)

            for n_sig in (1, 2, 3):
                width = 2 * n_sig * np.sqrt(eigvals[0])
                height = 2 * n_sig * np.sqrt(eigvals[1])
                ell = Ellipse((teff_final, logg_final), width, height, angle=angle,
                            facecolor="none", edgecolor="red", lw=1.0,
                            ls=("solid" if n_sig == 1 else "dashed" if n_sig == 2 else "dotted"),
                            zorder=4)
                ax3.add_patch(ell)

            ax3.set_xlabel("$T_\\mathrm{eff}$ (K)")
            ax3.set_ylabel("$\\log\\,g$")
            ax3.set_aspect(np.ptp(ax3.get_xlim()) / np.ptp(ax3.get_ylim()))
            ax3.tick_params(axis="x", rotation=45)

            # -- bottom-right panel: likelihood (zoomed) --
            # Use best-fit chi^2 as reference so exp(-0.5 * delta) = 1 at optimum
            chi2_best = float(model.loss(jnp.asarray(result.x), flux_obs_j, flux_err_j))
            dchi2_zoom = chi2_zoom - chi2_best
            likelihood_zoom = np.exp(-0.5 * dchi2_zoom)

            cf4 = ax4.contourf(teff_zoom, logg_zoom, likelihood_zoom,
                                levels=30, cmap="Blues")
            cbar4 = fig.colorbar(cf4, ax=ax4, orientation="horizontal",
                                location="top", pad=0.03)
            cbar4.set_label(r"$\exp(-\frac{1}{2}\,\Delta\chi^2)$", fontsize=9)
            cbar4.ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))

            # Delta-chi^2 contour levels for 2-parameter confidence regions
            delta_chi2_sigma = np.array([2.30, 6.18, 11.83])  # 1σ, 2σ, 3σ
            like_levels = np.exp(-0.5 * delta_chi2_sigma)      # descending: ~0.32, 0.045, 0.003
            cs4 = ax4.contour(teff_zoom, logg_zoom, likelihood_zoom,
                            levels=sorted(like_levels),
                            colors=["red"], linewidths=1.0,
                            linestyles=["dotted", "dashed", "solid"])
            ax4.clabel(cs4, fmt={v: f"${n}\\sigma$"
                                for n, v in zip([3, 2, 1], sorted(like_levels))},
                    fontsize=8)

            ax4.plot(teff_final, logg_final, "o", color="red", ms=9,
                    mec="k", mew=0.8, zorder=5)
            ax4.set_xlabel("$T_\\mathrm{eff}$ (K)")
            ax4.set_ylabel("$\\log\\,g$")
            ax4.set_aspect(np.ptp(ax4.get_xlim()) / np.ptp(ax4.get_ylim()))
            ax4.tick_params(axis="x", rotation=45)
        else:
            ax3.set_visible(False)
            ax4.set_visible(False)
        return fig


    return WhiteDwarfModel, h5, make_plot, np


@app.cell
def _(WhiteDwarfModel, np):
    λ = 10**(1e-4 * np.arange(4648) + 3.5523)

    model = WhiteDwarfModel(
        "data/sandbox/wd_da_spectra.h5", 
        n_continuum=8,
        obs_wavelength=λ
    )
    return (model,)


@app.cell
def _(np):
    # Get a DA-type WD spectrum from the block file
    from astropy.io import fits
    DR20_DIR = "data/release/dr20/spectro/astra/0.8.1"
    snow_white = fits.open(f"{DR20_DIR}/summary/astraAllStarSnowWhite-0.8.1.fits.gz")[1].data

    is_da_wd = (snow_white["classification"] == "DA")
    da_wd_sdss_ids = np.unique(snow_white["sdss_id"][is_da_wd])

    da_wd_sdss_ids = np.sort(da_wd_sdss_ids)
    print(f"{len(da_wd_sdss_ids):,} DA-type white dwarfs")
    return DR20_DIR, da_wd_sdss_ids, snow_white


@app.cell
def _(DR20_DIR, da_wd_sdss_ids, h5, np, snow_white):
    np.random.seed(8)

    da_wd_index = np.random.choice(len(da_wd_sdss_ids))
    sdss_id = da_wd_sdss_ids[da_wd_index]

    snow_white_results = snow_white[(snow_white["sdss_id"] == sdss_id)]

    with h5.File(f"{DR20_DIR}/spectra/block/mwmStarBlock-0.8.1.h5") as fp:
        sdss_ids = fp["boss/meta/sdss_id"][:]

        indices = np.searchsorted(sdss_ids, da_wd_sdss_ids)
        assert len(indices) == len(da_wd_sdss_ids)
        index = indices[da_wd_index]

        assert sdss_id == sdss_ids[index]
        flux = fp["boss/spectra/flux"][index]
        flux_err = 1.0/np.sqrt(fp["boss/spectra/ivar"][index])
        bad = (flux <= 0.0) | ~np.isfinite(flux) | ~np.isfinite(flux_err)
        flux[bad] = 0
        flux_err[bad] = 1e10
    return flux, flux_err, snow_white_results


@app.cell
def _(flux, flux_err, model):
    r, cov, unc = model.fit_with_uncertainties(flux, flux_err)
    return cov, r, unc


@app.cell
def _(model, r, unc):
    print(f"  Converged: {r.success}")
    print(f"{'Parameter':>12s} {'Fitted':>12s} {'+-1sig':>10s}")
    print("-" * 40)
    for i, name in enumerate(model.param_names):
        fmt = " >10.1f" if i < 3 else " >10.6e"
        print(f"{name:>12s} {float(r.x[i]):{fmt}} {float(unc[i]):{fmt}}")
    return


@app.cell
def _(cov, flux, flux_err, make_plot, model, r, unc):
    fig = make_plot(model, flux, flux_err, r, cov, unc, landscape=False)
    fig
    return


@app.cell
def _(mo, snow_white_results):
    from astropy.table import Table
    t = Table(snow_white_results)[["sdss_id", "teff", "e_teff", "logg", "e_logg"]]
    try:
        del t["sdss5_target_flags"]
    except:
        None
    mo.ui.table(t.to_pandas())
    return


@app.cell
def _(cov, flux, flux_err, make_plot, model, r, unc):
    # this takes a few minutes to compute
    fig_full = make_plot(model, flux, flux_err, r, cov, unc, landscape=True)
    fig_full
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
