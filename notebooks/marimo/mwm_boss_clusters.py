import marimo

__generated_with = "0.21.1"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import pandas as pd
    import numpy as np
    import h5py
    from scipy.ndimage import gaussian_filter1d
    from io import StringIO
    from matplotlib.colors import LogNorm, Normalize
    from pathlib import Path
    import fitsio

    base_dir = Path(__file__).resolve().parent
    font_path = str((base_dir.parent / "static" / "GoogleSans-Regular.ttf").resolve())
    mpl.font_manager.fontManager.addfont(font_path)
    font_prop = mpl.font_manager.FontProperties(fname=font_path)

    mpl.rcParams.update(
        {
            "xtick.top": True,
            "ytick.right": True,
            "xtick.bottom": True,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 8,
            "ytick.major.size": 8,
            "xtick.major.width": 1.5,
            "ytick.major.width": 1.5,
            "xtick.minor.visible": True,
            "ytick.minor.visible": True,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
            "xtick.minor.width": 1.5,
            "ytick.minor.width": 1.5,
            "xtick.minor.ndivs": 5,
            "ytick.minor.ndivs": 5,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "mathtext.default": "regular",
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.color": "gainsboro",
            "grid.linewidth": 0.75,
            "font.family": font_prop.get_name(),
        }
    )

    def decode_hdf5_bytes(df):

        for col in df.columns:
            if df[col].apply(lambda x: isinstance(x, bytes)).any():
                df[col] = df[col].apply(lambda x: x.decode("utf-8") if isinstance(x, bytes) else x)
        return df

    return (
        LogNorm,
        Normalize,
        Path,
        decode_hdf5_bytes,
        fitsio,
        gaussian_filter1d,
        h5py,
        mo,
        np,
        pd,
        plt,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <div style="display: flex; align-items: center; justify-content: center; position: relative; height: 80px; width: 100%;">
        <img src="https://www.sdss.org/wp-content/uploads/2022/09/sdss-new-logo-72dpi.png"
             style="height: 72px; position: absolute; right: 0;">

    <div style="text-align: center; max-width: 900px; margin: 0 auto;">
      <h1 style="font-size: 30px; margin-bottom: 5px; font-weight: 500;">
        MWM DR20 BOSS Star Cluster Explorer
      </h1>
      <p style="font-size: 16px; margin-top: 0; font-weight: 500;">
        Hunt & Reffert (2024) members with BOSS spectra
      </p>
    </div>

    </div>
    """)
    return


@app.cell(hide_code=True)
def _(decode_hdf5_bytes, fitsio, h5py, mo, pd):
    with mo.status.spinner("Loading spectra...") as load_spinner:

        with h5py.File("data/sandbox/dr20_boss_clusters.h5", "r") as f:

            hclu = f["/hclu"][()]
            allstar_hmem = f["/allstar_mem"][()]
            flux = f["/spectra/flux"][:]
            ivar = f["/spectra/ivar"][:]
            nmf = f["/spectra/nmf_rectified_model_flux"][:]
            continuum = f["/spectra/continuum"][:]
            wavelength = f["/spectra/wavelength"][()]

        foc = flux/continuum

        hclu = decode_hdf5_bytes(pd.DataFrame(hclu))
        allstar_hmem = decode_hdf5_bytes(pd.DataFrame(allstar_hmem))


        hclu = hclu.sort_values("N_stars_w_BOSS_spectrum", ascending=False).reset_index(drop=True)


        load_spinner.update("Grabbing BOSSNet information...")

    
        bn_path = 'data/release/dr20/spectro/astra/0.8.1/summary/astraAllStarBossNet-0.8.1.fits.gz'
        bn_cols = ['sdss_id', 'telescope', 'teff', 'e_teff', 'logg', 'e_logg', 'fe_h', 
                   'e_fe_h', 'result_flags', 'flag_warn', 'flag_bad', 'bn_v_r', 'e_bn_v_r']

        with fitsio.FITS(bn_path) as f:
        
            bn = pd.DataFrame(f[1].read(columns=bn_cols))

        allstar_hmem = allstar_hmem.merge(bn, on=["sdss_id", "telescope"], how="inner")
        del bn
        
    return allstar_hmem, continuum, flux, foc, hclu, ivar, nmf, wavelength


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <h2 style="text-align: left; font-weight: bold;"> Clusters and Moving Groups with DR20 BOSS Spectra </h2>
    """)
    return


@app.cell(hide_code=True)
def _(hclu, np):
    hclu['Association Name'] = hclu.Name.str.replace('_', ' ')
    hclu['HR24 Age (Myr)'] = np.round(10**hclu.logAge50 / 1e6, 1)
    hclu['HR24 Distance (pc)'] = np.round(hclu.dist50, 1)
    hclu['Number of Members with BOSS Spectrum'] = hclu.N_stars_w_BOSS_spectrum
    hclu["HR24 Association Type"] = hclu["Type"].map({
        "o": "open cluster",
        "m": "moving group",
        "g": "globular cluster",
        "d": "too distant to classify",
        "r": "rejected"
    })


    display_cols = ['Association Name', 'HR24 Association Type', 'Number of Members with BOSS Spectrum', 'HR24 Distance (pc)', 'HR24 Age (Myr)']

    hclu[display_cols]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <h2 style="text-align: left; font-weight: bold;"> Plot and Select </h2>

    Available columns are

    `['sdss_id', 'gaia_dr3_source_id', 'ra', 'dec', 'l', 'b', 'plx', 'e_plx', 'pmra', 'e_pmra', 'pmde', 'e_pmde', 'gaia_v_rad', 'gaia_e_v_rad', 'g_mag', 'bp_mag', 'rp_mag', 'j_mag', 'h_mag', 'k_mag', 'w1_mag', 'w2_mag', 'telescope', 'n_good_visits', 'v_rad', 'e_v_rad', 'std_v_rad', 'median_e_v_rad', 'xcsao_teff', 'xcsao_e_teff', 'xcsao_logg', 'xcsao_e_logg', 'xcsao_fe_h', 'xcsao_e_fe_h', 'xcsao_meanrxc', 'snr', 'zwarning_flags', 'nmf_rchi2', 'nmf_flags', 'HR24_cluster_name', 'HR24_mem_prob', 'RUWE', 'teff', 'e_teff', 'logg', 'e_logg', 'fe_h', 'e_fe_h', 'result_flags', 'flag_warn', 'flag_bad', 'bn_v_r', 'e_bn_v_r']`


    `['HR24_cluster_name', 'HR24_mem_prob', 'RUWE']` are taken from the [Hunt & Reffert (2024) catalog](https://ui.adsabs.harvard.edu/abs/2024yCat..36860042H/abstract).

    `['teff', 'e_teff', 'logg', 'e_logg', 'fe_h', 'e_fe_h', 'result_flags', 'flag_warn', 'flag_bad', 'bn_v_r', 'e_bn_v_r']` come from `astraAllStarBossNet-0.8.1.fits`.

    All other columns come from `mwmAllStar-0.8.1.fits`.

    To view the spectra, make a box selection by clicking and dragging on the plot, or hold `shift` while doing so to make a lasso selection.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    cluster_name = mo.ui.text(label="cluster name:", value="NGC 2516")
    x_col = mo.ui.text_area(value="g_mag - rp_mag", label="x")
    y_col = mo.ui.text_area(value="g_mag + 5*np.log10(plx/100)", label="y")

    cuts = mo.ui.text_area(
        placeholder="e.g.\ng_mag + 5*np.log10(plx/100) > 8.1\nsnr > 10", label="✂️"
    )

    x_label = mo.ui.text(value="G-G_{RP}", label="x label")
    y_label = mo.ui.text(value="M_G", label="y label")

    flip_x = mo.ui.checkbox(label="↔️")
    flip_y = mo.ui.checkbox(label="↕️", value=True)

    log_x = mo.ui.checkbox(label="log x")
    log_y = mo.ui.checkbox(label="log y")

    x_range = mo.ui.text(placeholder="x min, x max", label="x range")
    y_range = mo.ui.text(placeholder="y min, y max", label="y range")

    observatory_options = ["all", "APO", "LCO"]
    observatory_prompt = "observatory:"
    observatory = mo.ui.radio(
        options=observatory_options, value="all", label=observatory_prompt, inline=True
    )

    colorbar = mo.ui.checkbox(label="colorbar&nbsp;🌈")

    mo.vstack(
        [
            mo.hstack([cluster_name, observatory], justify="start", gap=2),
            mo.hstack([x_col, y_col, cuts], justify="start", gap=2),
            mo.hstack([x_label, y_label], justify="start", gap=2),
            mo.hstack(
                [x_range, y_range, log_x, log_y, flip_x, flip_y],
                justify="start",
                gap=2,
            ),
            colorbar,
        ],
        gap=2,
    )
    return (
        cluster_name,
        colorbar,
        cuts,
        flip_x,
        flip_y,
        log_x,
        log_y,
        observatory,
        x_col,
        x_label,
        x_range,
        y_col,
        y_label,
        y_range,
    )


@app.cell(hide_code=True)
def _(colorbar, mo):
    if colorbar.value:

        cb_col = mo.ui.text(value="snr", label="value")
        cb_label = mo.ui.text(value="snr", label="label")
        cb_range = mo.ui.text(placeholder="vmin, vmax", label="range")
        cb_cmap = mo.ui.text(value="turbo", label="colormap")
        cb_flip = mo.ui.checkbox(label="flip plotting order")
        cb_log = mo.ui.checkbox(label="log")

        cb_out = mo.hstack([cb_col, cb_label, cb_range, cb_cmap, cb_log, cb_flip], justify="start", gap=2)

    else:
        cb_out = mo.md('')

    cb_out
    return cb_cmap, cb_col, cb_flip, cb_label, cb_log, cb_range


@app.cell(hide_code=True)
def _(cluster_name):
    cluster_name_val = cluster_name.value.replace(" ", "_").strip()
    return (cluster_name_val,)


@app.cell(hide_code=True)
def _(cluster_name_val, hclu):
    row = hclu.loc[hclu.Name == cluster_name_val].iloc[0]

    age50_val = 10**row.logAge50 / 1e6
    age16_val = 10**row.logAge16 / 1e6
    age84_val = 10**row.logAge84 / 1e6

    age50 = f"{age50_val:,.1f}"
    age50_16 = f"{(age50_val - age16_val):,.1f}"
    age84_50 = f"{(age84_val - age50_val):,.1f}"

    age_str = rf"Age (Myr): ${age50}^{{+{age84_50}}}_{{-{age50_16}}}$"

    dist50_val = row.dist50
    dist16_val = row.dist16
    dist84_val = row.dist84

    dist50 = f"{dist50_val:,.1f}"
    dist50_16 = f"{(dist50_val - dist16_val):,.1f}"
    dist84_50 = f"{(dist84_val - dist50_val):,.1f}"

    dist_str = rf"Distance (pc): ${dist50}^{{+{dist84_50}}}_{{-{dist50_16}}}$"
    return age_str, dist_str


@app.cell(hide_code=True)
def _(
    LogNorm,
    Normalize,
    age_str,
    allstar_hmem,
    cb_cmap,
    cb_col,
    cb_flip,
    cb_label,
    cb_log,
    cb_range,
    cluster_name,
    cluster_name_val,
    colorbar,
    cuts,
    dist_str,
    flip_x,
    flip_y,
    log_x,
    log_y,
    mo,
    np,
    observatory,
    plt,
    x_col,
    x_label,
    x_range,
    y_col,
    y_label,
    y_range,
):
    allstar_hrd = allstar_hmem.loc[allstar_hmem.HR24_cluster_name == cluster_name_val].copy()

    n_stars_hrd = len(np.unique(allstar_hrd['sdss_id']))
    n_spectra_hrd = len(allstar_hrd)

    if observatory.value == "LCO":
        allstar_hrd = allstar_hrd[allstar_hrd["telescope"].str.startswith("lco")]
    elif observatory.value == "APO":
        allstar_hrd = allstar_hrd[allstar_hrd["telescope"].str.startswith("apo")]

    user_cuts = cuts.value.strip()

    if user_cuts:
        namespace_temp = {col: allstar_hrd[col] for col in allstar_hrd.columns}
        namespace_temp["np"] = np

        mask_cuts = np.ones(len(allstar_hrd), dtype=bool)
        for line in user_cuts.split("\n"):
            line = line.strip()
            if line:
                mask_cuts &= eval(line, {"__builtins__": {}}, namespace_temp)

        allstar_hrd = allstar_hrd[mask_cuts]

    namespace = {col: allstar_hrd[col] for col in allstar_hrd.columns}
    namespace["np"] = np

    hrd_fig, hrd_ax = plt.subplots(figsize=(10,6), constrained_layout=True)

    x_vals = eval(x_col.value, {"__builtins__": {}}, namespace)
    y_vals = eval(y_col.value, {"__builtins__": {}}, namespace)

    allstar_hrd = allstar_hrd.assign(x_vals=x_vals, y_vals=y_vals)

    hrd_fontsize = 16

    if colorbar.value:
        cb_vals = eval(cb_col.value, {"__builtins__": {}}, namespace)
        allstar_hrd = allstar_hrd.assign(cb_vals=cb_vals)

        ix_cb_vals = np.argsort(cb_vals)
        if cb_flip.value:
            ix_cb_vals = np.flip(ix_cb_vals)

        allstar_hrd = allstar_hrd.iloc[ix_cb_vals]

        if cb_log.value:
            norm = LogNorm(*tuple(float(x.strip()) for x in cb_range.value.split(","))) if cb_range.value else LogNorm()
        else:
            norm = Normalize(*tuple(float(x.strip()) for x in cb_range.value.split(","))) if cb_range.value else Normalize()

        sc = hrd_ax.scatter(
            allstar_hrd["x_vals"],
            allstar_hrd["y_vals"],
            c=allstar_hrd["cb_vals"],
            s=10,
            edgecolors="gainsboro",
            lw=.5,
            cmap=cb_cmap.value,
            norm=norm,
        )
        cbar = plt.colorbar(sc, ax=hrd_ax, pad=0.01)
        cbar.set_label(f"${cb_label.value}$", fontsize=hrd_fontsize)

    else:
        hrd_ax.scatter(
            allstar_hrd["x_vals"],
            allstar_hrd["y_vals"],
            c="k",
            s=10,
            edgecolors="gainsboro",
            lw=.5,
        )

    if x_range.value:
        xmin, xmax = tuple(float(x.strip()) for x in x_range.value.split(","))
        hrd_ax.set_xlim(xmin, xmax)

    if y_range.value:
        ymin, ymax = tuple(float(x.strip()) for x in y_range.value.split(","))
        hrd_ax.set_ylim(ymin, ymax)

    if flip_x.value:
        hrd_ax.invert_xaxis()
    if flip_y.value:
        hrd_ax.invert_yaxis()

    if log_x.value:
        hrd_ax.semilogx()
    if log_y.value:
        hrd_ax.semilogy()

    if x_label.value:
        hrd_ax.set_xlabel(f"${x_label.value}$", fontsize=hrd_fontsize)
    if y_label.value:
        hrd_ax.set_ylabel(f"${y_label.value}$", fontsize=hrd_fontsize)

    hrd_ax.grid(True, which="both", zorder=-100)

    hrd_tit_str_left = cluster_name.value
    hrd_tit_str_right = dist_str+'    '+age_str
    hrd_ax.set_title(hrd_tit_str_left, loc='left', fontsize=hrd_fontsize)
    hrd_ax.set_title(hrd_tit_str_right, loc='right', fontsize=hrd_fontsize)

    hrd = mo.ui.matplotlib(plt.gca(), debounce=True)

    n_stars_plot_hrd = len(np.unique(allstar_hrd['sdss_id']))


    n_stars_display = mo.stat(
        value=n_stars_hrd,
        label=f"Stars in {cluster_name.value} with DR20 BOSS spectra",
    )

    n_spectra_display = mo.stat(
        value=n_spectra_hrd,
        label=f'Star-level DR20 BOSS Spectra in {cluster_name.value}',
        caption='"Star-level" --> 1 spectrum per SDSS ID per telescope'
    )

    n_stars_plot_display = mo.stat(
        value=n_stars_plot_hrd,
        label=f"Stars shown on scatter plot",
    )

    mo.vstack([
        mo.hstack([n_stars_display, n_spectra_display, n_stars_plot_display], justify='center', gap=10),
        mo.hstack([hrd], justify="center"),  
    ], gap=2)
    return allstar_hrd, hrd


@app.cell(hide_code=True)
def _(allstar_hrd, hrd):
    select_mask = hrd.value.get_mask(
        allstar_hrd["x_vals"], allstar_hrd["y_vals"]
    )

    allstar_select = allstar_hrd[select_mask]
    return (allstar_select,)


@app.cell(hide_code=True)
def _(mo):
    spec_color = mo.ui.text(value="g_mag - rp_mag", label="color spectra by")
    spec_cmap = mo.ui.text(value="turbo", label="colormap")
    smoothing = mo.ui.checkbox(label="smooth spectra with gaussian filter", value=False)
    spec_ranges = mo.ui.checkbox(label="customize axis bounds")
    spec_int = mo.ui.checkbox(label="make interactive")



    pan1_xrange = mo.ui.text(label=r"top", value="3550, 10450")
    pan2_xrange = mo.ui.text(label=r"left", value="6523, 6603")
    pan3_xrange = mo.ui.text(label=r"center", value="8150, 8230")
    pan4_xrange = mo.ui.text(label=r"right", value="8460, 8700")

    pan1_yrange = mo.ui.text(label=r"top", value="0, 1.5")
    pan2_yrange = mo.ui.text(label=r"left", value="0, 1.5")
    pan3_yrange = mo.ui.text(label=r"center", value="0, 1.5")
    pan4_yrange = mo.ui.text(label=r"right", value="0, 1.5")
    return (
        pan1_xrange,
        pan1_yrange,
        pan2_xrange,
        pan2_yrange,
        pan3_xrange,
        pan3_yrange,
        pan4_xrange,
        pan4_yrange,
        smoothing,
        spec_cmap,
        spec_color,
        spec_int,
        spec_ranges,
    )


@app.cell(hide_code=True)
def _(
    allstar_select,
    mo,
    smoothing,
    spec_cmap,
    spec_color,
    spec_int,
    spec_ranges,
):
    if len(allstar_select):
        spec_color_prompt = mo.hstack([spec_color, spec_cmap, smoothing, spec_int, spec_ranges], justify="start", gap=2)
    else:
        spec_color_prompt = mo.md('')

    spec_color_prompt
    return


@app.cell(hide_code=True)
def _(allstar_select, np, spec_color):
    namespace_select = {col: allstar_select[col] for col in allstar_select.columns}
    namespace_select["np"] = np

    spec_color_vals = eval(spec_color.value, {"__builtins__": {}}, namespace_select)

    allstar_so = allstar_select.assign(spec_color_vals=spec_color_vals)

    ix_spec_color_vals = np.argsort(spec_color_vals)
    allstar_so = allstar_so.iloc[ix_spec_color_vals].reset_index(drop=True)
    return (allstar_so,)


@app.cell(hide_code=True)
def _(
    allstar_select,
    mo,
    pan1_xrange,
    pan1_yrange,
    pan2_xrange,
    pan2_yrange,
    pan3_xrange,
    pan3_yrange,
    pan4_xrange,
    pan4_yrange,
    spec_ranges,
):
    if spec_ranges.value and len(allstar_select):
        pan_bounds = mo.hstack(
            [
                mo.vstack(
                    [
                        mo.md(r"$\mathrm{\lambda}$:"),
                        mo.md(r"flux / continuum:"),
                    ],
                    gap=2,
                ),
                mo.vstack([pan1_xrange, pan1_yrange], gap=2),
                mo.vstack([pan2_xrange, pan2_yrange], gap=2),
                mo.vstack([pan3_xrange, pan3_yrange], gap=2),
                mo.vstack([pan4_xrange, pan4_yrange], gap=2),
            ],
            gap=2,
        )
    else:
        pan_bounds = mo.md('')

    pan_bounds
    return


@app.cell(hide_code=True)
def _(
    np,
    pan1_xrange,
    pan1_yrange,
    pan2_xrange,
    pan2_yrange,
    pan3_xrange,
    pan3_yrange,
    pan4_xrange,
    pan4_yrange,
    wavelength,
):
    pan1_xmin, pan1_xmax = tuple(float(x.strip()) for x in pan1_xrange.value.split(","))
    pan2_xmin, pan2_xmax = tuple(float(x.strip()) for x in pan2_xrange.value.split(","))
    pan3_xmin, pan3_xmax = tuple(float(x.strip()) for x in pan3_xrange.value.split(","))
    pan4_xmin, pan4_xmax = tuple(float(x.strip()) for x in pan4_xrange.value.split(","))

    pan1_ymin, pan1_ymax = tuple(float(x.strip()) for x in pan1_yrange.value.split(","))
    pan2_ymin, pan2_ymax = tuple(float(x.strip()) for x in pan2_yrange.value.split(","))
    pan3_ymin, pan3_ymax = tuple(float(x.strip()) for x in pan3_yrange.value.split(","))
    pan4_ymin, pan4_ymax = tuple(float(x.strip()) for x in pan4_yrange.value.split(","))

    ax1_ix = np.where((wavelength > pan1_xmin) & (wavelength < pan1_xmax))
    ax2_ix = np.where((wavelength > pan2_xmin) & (wavelength < pan2_xmax))
    ax3_ix = np.where((wavelength > pan3_xmin) & (wavelength < pan3_xmax))
    ax4_ix = np.where((wavelength > pan4_xmin) & (wavelength < pan4_xmax))

    ax_ixs = [ax1_ix, ax2_ix, ax3_ix, ax4_ix]
    ax_lams = [wavelength[ax_ix] for ax_ix in ax_ixs]
    return (
        ax_ixs,
        ax_lams,
        pan1_xmax,
        pan1_xmin,
        pan1_ymax,
        pan1_ymin,
        pan2_xmax,
        pan2_xmin,
        pan2_ymax,
        pan2_ymin,
        pan3_xmax,
        pan3_xmin,
        pan3_ymax,
        pan3_ymin,
        pan4_xmax,
        pan4_xmin,
        pan4_ymax,
        pan4_ymin,
    )


@app.cell(hide_code=True)
def _(
    allstar_select,
    allstar_so,
    ax_ixs,
    ax_lams,
    foc,
    gaussian_filter1d,
    mo,
    np,
    pan1_xmax,
    pan1_xmin,
    pan1_ymax,
    pan1_ymin,
    pan2_xmax,
    pan2_xmin,
    pan2_ymax,
    pan2_ymin,
    pan3_xmax,
    pan3_xmin,
    pan3_ymax,
    pan3_ymin,
    pan4_xmax,
    pan4_xmin,
    pan4_ymax,
    pan4_ymin,
    plt,
    smoothing,
    spec_cmap,
    spec_int,
):
    if len(allstar_select):

        fig = plt.figure(figsize=(10,6), constrained_layout=True)

        gs = fig.add_gridspec(2, 3)

        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[1, 1])
        ax4 = fig.add_subplot(gs[1, 2])

        axes = [ax1, ax2, ax3, ax4]

        spec_color_cmap = plt.get_cmap(spec_cmap.value)

        color_positions = np.linspace(0, 1, len(allstar_so))

        ix_spec = allstar_so.ix_spectrum.values

        flux_sel = foc[ix_spec]

        for spec_i in reversed(range(len(allstar_so))):

            color = spec_color_cmap(color_positions[spec_i])

            for ax_i, (ax_ix, spec_ax, ax_lam) in enumerate(zip(ax_ixs, axes, ax_lams)):

                spec_ax.plot(
                    ax_lam,
                    gaussian_filter1d(flux_sel[spec_i][ax_ix], sigma=3) if (smoothing.value) else flux_sel[spec_i][ax_ix],
                    lw=0.5,
                    color=color,
                )

        ax1.set_xlim(pan1_xmin, pan1_xmax)
        ax2.set_xlim(pan2_xmin, pan2_xmax)
        ax3.set_xlim(pan3_xmin, pan3_xmax)
        ax4.set_xlim(pan4_xmin, pan4_xmax)

        ax1.set_ylim(pan1_ymin, pan1_ymax)
        ax2.set_ylim(pan2_ymin, pan2_ymax)
        ax3.set_ylim(pan3_ymin, pan3_ymax)
        ax4.set_ylim(pan4_ymin, pan4_ymax)




        n_stars = len(np.unique(allstar_so['sdss_id']))
        n_spectra = len(allstar_so)
        tit_str_left = f"{n_stars:,} stars" if n_stars > 1 else f"{n_stars} star"
        tit_str_right = f"{n_spectra:,} spectra" if n_spectra > 1 else f"{n_spectra} spectrum"
        ax1.set_title(tit_str_left, loc='left', fontsize=16)
        ax1.set_title(tit_str_right, loc='right', fontsize=16)

        fig.supylabel("flux / continuum", fontsize=14)
        ax3.set_xlabel(r"$\lambda\ (\AA)$", fontsize=14)

        for spec_ax in axes:
            spec_ax.grid(True, which="both", alpha=0.4, zorder=-100)
            spec_ax.tick_params(axis="both", labelsize=12)

        if spec_int.value:
            specfig = mo.vstack([mo.md('To zoom, select the 🔍 emoji below the figure and click and drag in any of the panels.'), 
                                 mo.mpl.interactive(fig)], gap=2)
        else:
            specfig = mo.hstack([fig], justify="center")

    else:
        specfig = mo.md("")

    specfig
    return


@app.cell(hide_code=True)
def _(allstar_select, mo):
    if len(allstar_select):
        spec_df_display_check = mo.ui.checkbox(label=r"display `mwmAllStar-0.8.1.fits` information for selected subset")
    else:
        spec_df_display_check = mo.md('')

    spec_df_display_check
    return (spec_df_display_check,)


@app.cell(hide_code=True)
def _():
    allstar_cols = ['sdss_id', 'gaia_dr3_source_id', 'ra', 'dec', 'l', 'b', 'plx', 'e_plx',
           'pmra', 'e_pmra', 'pmde', 'e_pmde', 'gaia_v_rad', 'gaia_e_v_rad',
           'g_mag', 'bp_mag', 'rp_mag', 'j_mag', 'h_mag', 'k_mag', 'w1_mag',
           'w2_mag', 'telescope', 'n_good_visits', 'v_rad', 'e_v_rad', 'std_v_rad',
           'median_e_v_rad', 'xcsao_teff', 'xcsao_e_teff', 'xcsao_logg',
           'xcsao_e_logg', 'xcsao_fe_h', 'xcsao_e_fe_h', 'xcsao_meanrxc', 'snr',
           'zwarning_flags', 'nmf_rchi2', 'nmf_flags', 'ix_spectrum',
           'HR24_cluster_name', 'HR24_mem_prob', 'RUWE', 'teff', 'e_teff', 'logg',
           'e_logg', 'fe_h', 'e_fe_h', 'result_flags', 'flag_warn', 'flag_bad',
           'bn_v_r', 'e_bn_v_r']
    return (allstar_cols,)


@app.cell(hide_code=True)
def _(allstar_cols, allstar_so, mo, spec_df_display_check):
    if spec_df_display_check.value:
        spec_df = allstar_so[allstar_cols]
    else:
        spec_df = mo.md('')
    spec_df
    return


@app.cell(hide_code=True)
def _(cluster_name_val, mo):
    outfilename = mo.ui.text(label="output file name:", value=f"dr20_boss_{cluster_name_val}")
    return (outfilename,)


@app.cell(hide_code=True)
def _(cluster_name, mo, outfilename):
    outfile = f"home/data/{outfilename.value}.h5"

    save_options = [f"all {cluster_name.value} members", "subset selected above"]
    save_subset_option = mo.ui.radio(
        options=save_options, value=f"all {cluster_name.value} members", inline=True
    )

    save_spectra_button = mo.ui.run_button(
        label=f"save spectra to `{outfile}`&nbsp;&nbsp;💾"
    )

    download_md_body = mo.md("After saving, click the file tree icon in the top left of this page and download your file.")
    return download_md_body, outfile, save_spectra_button, save_subset_option


@app.cell(hide_code=True)
def _(
    Path,
    allstar_cols,
    allstar_hmem,
    allstar_select,
    cluster_name,
    cluster_name_val,
    continuum,
    flux,
    h5py,
    ivar,
    mo,
    nmf,
    outfile,
    outfilename,
    save_spectra_button,
    save_subset_option,
    wavelength,
):
    if save_spectra_button.value:

        if outfilename.value:

            with mo.status.spinner(f"Saving {outfile}..."):

                outfile_path = Path(outfile)
                outfile_path.parent.mkdir(parents=True, exist_ok=True)

                with h5py.File(outfile_path, "w") as outfile_f:
                    if save_subset_option.value == f"all {cluster_name.value} members":

                        allstar_allmem = allstar_hmem.loc[allstar_hmem.HR24_cluster_name==cluster_name_val]
                        ix_spec_allmem = allstar_allmem.ix_spectrum.values                                    

                        allstar_rec = allstar_allmem.assign(
                                                        telescope=lambda df: df["telescope"].astype("S6"),
                                                        HR24_cluster_name=lambda df: df["HR24_cluster_name"].astype("S32")
                                                    ).to_records(index=False)

                        outfile_f.create_dataset("/allstar", data=allstar_rec[allstar_cols])
                        outfile_f.create_dataset("/spectra/wavelength", data=wavelength)
                        outfile_f.create_dataset("/spectra/flux", data=flux[ix_spec_allmem], compression="gzip")
                        outfile_f.create_dataset("/spectra/ivar", data=ivar[ix_spec_allmem], compression="gzip")
                        outfile_f.create_dataset("/spectra/nmf_rectified_model_flux", data=nmf[ix_spec_allmem], compression="gzip")
                        outfile_f.create_dataset("/spectra/continuum", data=continuum[ix_spec_allmem], compression="gzip")
                        savemessage = mo.md(f"Done saving to `{outfile}` ✅ ")
                    else:
                        if len(allstar_select):
                            ix_spec_select = allstar_select.ix_spectrum.values
                            allstar_rec = allstar_select.assign(
                                                        telescope=lambda df: df["telescope"].astype("S6"),
                                                        HR24_cluster_name=lambda df: df["HR24_cluster_name"].astype("S32")
                                                    ).to_records(index=False)

                            outfile_f.create_dataset("/allstar", data=allstar_rec[allstar_cols])
                            outfile_f.create_dataset("/spectra/wavelength", data=wavelength)
                            outfile_f.create_dataset("/spectra/flux", data=flux[ix_spec_select], compression="gzip")
                            outfile_f.create_dataset("/spectra/ivar", data=ivar[ix_spec_select], compression="gzip")
                            outfile_f.create_dataset("/spectra/nmf_rectified_model_flux", data=nmf[ix_spec_select], compression="gzip")
                            outfile_f.create_dataset("/spectra/continuum", data=continuum[ix_spec_select], compression="gzip")


                            savemessage = mo.md(f"Done saving to `{outfile}` ✅ ")

                        else:
                            savemessage = mo.md(
                                '<span style="color: red;">Please make a selection on the plot above</span>'
                            )

        else:
            savemessage = mo.md(
                '<span style="color: red;">Please enter an output file name</span>'
            )

    else:
        savemessage = mo.md("")
    return (savemessage,)


@app.cell(hide_code=True)
def _(
    download_md_body,
    mo,
    outfilename,
    save_spectra_button,
    save_subset_option,
    savemessage,
):
    download_md = mo.vstack(
        [
            save_subset_option,
            outfilename,
            save_spectra_button,
            download_md_body,
            savemessage,
        ],
        gap=2,
    )


    mo.accordion({'<h2 style="text-align: left; font-weight: bold;"> Download the spectra </h2>': download_md})
    return


@app.cell(hide_code=True)
def _(mo, outfilename):
    access_md_filename = mo.md(f"""`{outfilename.value}.h5` is structured """)

    access_md_body1 = mo.md(
        r"""
    * `/allstar`
    * `/spectra`
        * `/flux  (n, 4648)  float32`
        * `/continuum  (n, 4648)  float32`
        * `/ivar  (n, 4648)  float32`
        * `/nmf_rectified_model_flux  (n, 4648)  float32`
        * `/wavelength (4648,) float32`

    Access the spectra and `mwmAllStar-0.8.1.fits` information with
    """
    )
    access_md_code = mo.md(
        f"""
    ```python
    with h5py.File(f"/home/jovyan/home/data/{outfilename.value}.h5", "r") as spectra_f:

        allstar_loaded = spectra_f["/allstar"][()]
        flux_loaded = spectra_f["/spectra/flux"][:]
        ivar_loaded = spectra_f["/spectra/ivar"][:]
        nmf_loaded = spectra_f["/spectra/nmf_rectified_model_flux"][:]
        cont_loaded = spectra_f["/spectra/continuum"][:]
        wavelength_loaded = spectra_f["/spectra/wavelength"][()]
    ```
    """
    )

    access_md_body2 = mo.md(
        r"""
    Code snippets can be pasted right into this notebook after saving `dr20_boss_spectra.h5` if you wish 🤠

    `/allstar` comes with 

    `['sdss_id', 'gaia_dr3_source_id', 'ra', 'dec', 'l', 'b', 'plx', 'e_plx', 'pmra', 'e_pmra', 'pmde', 'e_pmde', 'gaia_v_rad', 'gaia_e_v_rad', 'g_mag', 'bp_mag', 'rp_mag', 'j_mag', 'h_mag', 'k_mag', 'w1_mag', 'w2_mag', 'telescope', 'n_good_visits', 'v_rad', 'e_v_rad', 'std_v_rad', 'median_e_v_rad', 'xcsao_teff', 'xcsao_e_teff', 'xcsao_logg', 'xcsao_e_logg', 'xcsao_fe_h', 'xcsao_e_fe_h', 'xcsao_meanrxc', 'snr', 'zwarning_flags', 'nmf_rchi2', 'nmf_flags', 'HR24_cluster_name', 'HR24_mem_prob', 'RUWE', 'teff', 'e_teff', 'logg', 'e_logg', 'fe_h', 'e_fe_h', 'result_flags', 'flag_warn', 'flag_bad', 'bn_v_r', 'e_bn_v_r']`


    Tip 💡:  `/allstar` can immediately be turned into a `pd.DataFrame()`

    ```python
    allstar_df = pd.DataFrame(allstar_loaded)
    ```

    Grab a subset of spectra with

    ```python
    # Examples
    # cond = (allstar_loaded['sdss_id']==111751124)
    # cond = (allstar_loaded['g_mag'] < 15)
    # cond = ((allstar_loaded['g_mag'] < 15) & (allstar_loaded['dec'] > -30))

    cond = (allstar_loaded['telescope']==b'lco25m')
    ix_stars = np.where(cond)[0]

    flux_stars = flux_loaded[ix_stars, :]
    cont_stars = cont_loaded[ix_stars, :]
    normalized_flux = flux_stars / cont_stars
    ```
    """
    )

    access_md = mo.vstack([access_md_filename, access_md_body1, access_md_code, access_md_body2], gap=0)

    acc_out = mo.accordion({'<h2 style="text-align: left; font-weight: bold;"> Access the spectra </h2>': access_md})

    acc_out
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Please send questions/issues/feedback about this notebook to Kayvon Sharifi (ksharifi1@gsu.edu)
    """)
    return


if __name__ == "__main__":
    app.run()
