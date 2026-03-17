import marimo

__generated_with = "0.20.4"
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
    import fitsio
    from io import StringIO
    from matplotlib.colors import LogNorm, Normalize
    from pathlib import Path

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
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "mathtext.default": "regular",
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.color": "gainsboro",
            "grid.linewidth": 0.75,
            "font.family": font_prop.get_name(),
        }
    )
    return (
        LogNorm,
        Normalize,
        Path,
        StringIO,
        base_dir,
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

        <h1 style="margin: 0; text-align: center; width: 100%;">MWM DR20 BOSS Explorer</h1>
    </div>
    """)
    return


@app.cell(hide_code=True)
def _(base_dir):
    block_path = "data/release/dr20/spectro/astra/0.8.1/spectra/block/mwmStarBlock-0.8.1.h5"
    allstar_path = "data/release/dr20/spectro/astra/0.8.1/summary/mwmAllStar-0.8.1.fits.gz"
    block_meta_arr_path = str((base_dir.parent / "static" / "block_meta.npy").resolve())
    return allstar_path, block_meta_arr_path, block_path


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <h2 style="text-align: left; font-weight: bold;"> Upload list of SDSS or <i>Gaia</i> DR3 source IDs </h2>
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    usr_sdss_id_list_upload_button = mo.ui.file(
        kind="button", label="upload SDSS IDs (.txt, .csv)", filetypes=[".txt", ".csv"]
    )
    usr_gaia_id_list_upload_button = mo.ui.file(
        kind="button",
        label="upload Gaia DR3 source IDs (.txt, .csv)",
        filetypes=[".txt", ".csv"],
    )

    id_list_buttons = mo.hstack(
        [
            usr_sdss_id_list_upload_button,
            usr_gaia_id_list_upload_button,
        ],
        justify="start",
        gap=2,
    )

    id_list_buttons
    return usr_gaia_id_list_upload_button, usr_sdss_id_list_upload_button


@app.cell(hide_code=True)
def _(
    StringIO,
    allstar_path,
    fitsio,
    mo,
    np,
    usr_gaia_id_list_upload_button,
    usr_sdss_id_list_upload_button,
):
    sdss_files = usr_sdss_id_list_upload_button.value
    gaia_files = usr_gaia_id_list_upload_button.value

    n_sdss = len(sdss_files)
    n_gaia = len(gaia_files)

    uploaded_cond = (n_sdss > 0) or (n_gaia > 0)
    too_many = (n_sdss > 0) and (n_gaia > 0)

    tipmessage = (
        mo.md(
            "Tip 💡: "
            '<span style="color: red;">Click to clear files</span> above before uploading a new ID list'
        )
        if uploaded_cond
        else mo.md("")
    )

    toomanymessage = (
        mo.md('<span style="color: red;">Too many ID lists uploaded </span>')
        if too_many
        else mo.md("")
    )

    allstar_rows = np.array([], dtype=int)

    if uploaded_cond and not too_many:

        files = sdss_files if n_sdss else gaia_files
        column_name = "sdss_id" if n_sdss else "gaia_dr3_source_id"

        usr_id_list_txt = files[0].contents.decode("utf-8")

        try:
            usr_ids = np.loadtxt(StringIO(usr_id_list_txt), dtype=np.int64)
        except ValueError:
            usr_ids = np.loadtxt(StringIO(usr_id_list_txt), dtype=np.int64, skiprows=1)

        with fitsio.FITS(allstar_path) as allstar_f:
            allstar_ids = allstar_f[1][column_name][:]

        allstar_rows = np.where(np.isin(allstar_ids, usr_ids))[0]

    nospectramessage = (
        mo.md(
            '<span style="color: red;">Sorry, no DR20 BOSS spectra for those IDs</span> 😵‍💫'
        )
        if uploaded_cond and len(allstar_rows) == 0 and not too_many
        else mo.md("")
    )
    matched_cond = len(allstar_rows) > 0

    upload_messages = mo.vstack([tipmessage, toomanymessage, nospectramessage], gap=2)

    upload_messages
    return allstar_rows, matched_cond


@app.cell
def _():
    allstar_cols = ['sdss_id', 'gaia_dr3_source_id', 'ra', 'dec', 'l', 'b', 'plx', 'e_plx', 'pmra', 'e_pmra', 'pmde', 'e_pmde', 'gaia_v_rad', 'gaia_e_v_rad', 'g_mag', 'bp_mag', 'rp_mag', 'j_mag', 'h_mag', 'k_mag', 'w1_mag', 'w2_mag', 'telescope', 'n_good_visits', 'v_rad', 'e_v_rad', 'std_v_rad', 'median_e_v_rad', 'xcsao_teff', 'xcsao_e_teff', 'xcsao_logg', 'xcsao_e_logg', 'xcsao_fe_h', 'xcsao_e_fe_h', 'xcsao_meanrxc', 'snr', 'zwarning_flags', 'nmf_rchi2', 'nmf_flags']


    # Below are all the columns you can grab from mwmAllStar-0.8.1.fits

    # WARNING: The sdss5_target_flags column is two-dimensional and reading it in will interfere with the functionality of this notebook

    # ['sdss_id', 'sdss4_apogee_id', 'gaia_dr2_source_id', 'gaia_dr3_source_id', 'tic_v8_id', 'healpix', 'lead', 'version_id', 'catalogid', 'catalogid21', 'catalogid25', 'catalogid31', 'n_associated', 'n_neighborhood', 'crossmatch_flags', 'sdss5_target_flags', 'sdss4_apogee_target1_flags', 'sdss4_apogee_target2_flags', 'sdss4_apogee2_target1_flags', 'sdss4_apogee2_target2_flags', 'sdss4_apogee2_target3_flags', 'sdss4_apogee_member_flags', 'sdss4_apogee_extra_target_flags', 'sdss5_dr19_apogee_flag', 'ra', 'dec', 'l', 'b', 'plx', 'e_plx', 'pmra', 'e_pmra', 'pmde', 'e_pmde', 'gaia_v_rad', 'gaia_e_v_rad', 'g_mag', 'bp_mag', 'rp_mag', 'j_mag', 'e_j_mag', 'h_mag', 'e_h_mag', 'k_mag', 'e_k_mag', 'ph_qual', 'bl_flg', 'cc_flg', 'w1_mag', 'e_w1_mag', 'w1_flux', 'w1_dflux', 'w1_frac', 'w2_mag', 'e_w2_mag', 'w2_flux', 'w2_dflux', 'w2_frac', 'w1uflags', 'w2uflags', 'w1aflags', 'w2aflags', 'mag4_5', 'd4_5m', 'rms_f4_5', 'sqf_4_5', 'mf4_5', 'csf', 'zgr_teff', 'zgr_e_teff', 'zgr_logg', 'zgr_e_logg', 'zgr_fe_h', 'zgr_e_fe_h', 'zgr_e', 'zgr_e_e', 'zgr_plx', 'zgr_e_plx', 'zgr_teff_confidence', 'zgr_logg_confidence', 'zgr_fe_h_confidence', 'zgr_ln_prior', 'zgr_chi2', 'zgr_quality_flags', 'r_med_geo', 'r_lo_geo', 'r_hi_geo', 'r_med_photogeo', 'r_lo_photogeo', 'r_hi_photogeo', 'bailer_jones_flags', 'ebv', 'e_ebv', 'ebv_flags', 'ebv_zhang_2023', 'e_ebv_zhang_2023', 'ebv_sfd', 'e_ebv_sfd', 'ebv_rjce_glimpse', 'e_ebv_rjce_glimpse', 'ebv_rjce_allwise', 'e_ebv_rjce_allwise', 'ebv_bayestar_2019', 'e_ebv_bayestar_2019', 'ebv_edenhofer_2023', 'e_ebv_edenhofer_2023', 'c_star', 'u_jkc_mag', 'u_jkc_mag_flag', 'b_jkc_mag', 'b_jkc_mag_flag', 'v_jkc_mag', 'v_jkc_mag_flag', 'r_jkc_mag', 'r_jkc_mag_flag', 'i_jkc_mag', 'i_jkc_mag_flag', 'u_sdss_mag', 'u_sdss_mag_flag', 'g_sdss_mag', 'g_sdss_mag_flag', 'r_sdss_mag', 'r_sdss_mag_flag', 'i_sdss_mag', 'i_sdss_mag_flag', 'z_sdss_mag', 'z_sdss_mag_flag', 'y_ps1_mag', 'y_ps1_mag_flag', 'n_boss_visits', 'boss_min_mjd', 'boss_max_mjd', 'n_apogee_visits', 'apogee_min_mjd', 'apogee_max_mjd', 'created', 'modified', 'spectrum_pk', 'source', 'release', 'filetype', 'v_astra', 'run2d', 'telescope', 'min_mjd', 'max_mjd', 'n_visits', 'n_good_visits', 'n_good_rvs', 'v_rad', 'e_v_rad', 'std_v_rad', 'median_e_v_rad', 'xcsao_teff', 'xcsao_e_teff', 'xcsao_logg', 'xcsao_e_logg', 'xcsao_fe_h', 'xcsao_e_fe_h', 'xcsao_meanrxc', 'snr', 'gri_gaia_transform_flags', 'zwarning_flags', 'nmf_rchi2', 'nmf_flags']
    return (allstar_cols,)


@app.cell(hide_code=True)
def _(np):
    def read_h5_datasets(f, datasets, indices):
        dsets = [f[d] for d in datasets]

        n = len(indices)
        outputs = [np.empty((n,) + d.shape[1:], dtype=d.dtype) for d in dsets]

        breaks = np.where(np.diff(indices) != 1)[0] + 1
        runs = np.split(indices, breaks)

        pos = 0
        for r in runs:
            n_run = len(r)
            s = slice(r[0], r[-1] + 1)

            for d, out in zip(dsets, outputs):
                out[pos:pos+n_run] = d[s]

            pos += n_run

        return outputs


    datasets = [
        "boss/spectra/flux",
        "boss/spectra/continuum",
        "boss/spectra/ivar",
        "boss/spectra/nmf_rectified_model_flux",
    ]
    return datasets, read_h5_datasets


@app.cell(hide_code=True)
def _(
    allstar_cols,
    allstar_path,
    allstar_rows,
    block_meta_arr_path,
    block_path,
    datasets,
    fitsio,
    h5py,
    matched_cond,
    mo,
    np,
    pd,
    read_h5_datasets,
):
    mo.stop(not matched_cond)

    with mo.status.spinner("Loading mwmAllStar information...") as load_spinner:

        with fitsio.FITS(allstar_path) as f:
            allstar = f[1].read(columns=allstar_cols, rows=allstar_rows)
            allstar = pd.DataFrame(allstar.astype(allstar.dtype.newbyteorder("="))).replace(np.inf, np.nan)

        block_meta_arr = np.load(block_meta_arr_path)

        load_spinner.update("Mapping mwmAllStar rows to block file indices...")

        allstar_tel = allstar["telescope"].astype(str).values
        block_tel = block_meta_arr["telescope"].astype(str)

        block_keys = np.rec.fromarrays(
            [block_meta_arr["sdss_id"], block_tel],
            names="sdss_id,telescope"
        )

        allstar_keys = np.rec.fromarrays(
            [allstar["sdss_id"].values, allstar_tel],
            names="sdss_id,telescope"
        )

        block_sort = np.argsort(block_keys)
        block_keys_sorted = block_keys[block_sort]

        ix = block_sort[np.searchsorted(block_keys_sorted, allstar_keys)]

        sort_ix = np.argsort(ix)
        ix_sorted = ix[sort_ix]

        allstar = allstar.iloc[sort_ix].reset_index(drop=True)

        allstar["ix_in_blockfile"] = ix_sorted
        allstar["nb_flux_arr_ix"] = np.arange(len(allstar))

        load_spinner.update("Reading in spectra...")

        with h5py.File(block_path, "r", rdcc_nbytes=2*1024**3) as f:
            flux, continuum, ivar, nmf_rectified_model_flux = read_h5_datasets(
                f, datasets, ix_sorted
            )


        flux_over_cont = flux / continuum


    topmessage = mo.md("Uploaded IDs matched to `mwmAllStar-0.8.1.fits` are shown below")

    valuetipmessage = mo.md("Tip 💡: `['sdss_id', 'telescope']` value pairs are unique across rows")

    outmessage = mo.vstack([topmessage, allstar[allstar_cols], valuetipmessage], gap=2)

    outmessage
    return (
        allstar,
        continuum,
        flux,
        flux_over_cont,
        ivar,
        nmf_rectified_model_flux,
    )


@app.cell(hide_code=True)
def _(block_path, h5py):
    with h5py.File(block_path, "r") as wavelength_block_f:
        wavelength = wavelength_block_f["boss/spectra/wavelength"][()]
    return (wavelength,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <h2 style="text-align: left; font-weight: bold;"> Plot and Select </h2>

    Make a scatter plot with any of the `mwmAllStar-0.8.1.fits` columns shown above.

    To view the spectra, make a box selection by clicking and dragging on the plot, or hold `shift` while doing so to make a lasso selection.
    """)
    return


@app.cell(hide_code=True)
def _(matched_cond, mo):
    mo.stop(not matched_cond)

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
            mo.hstack([x_col, y_col, cuts, observatory], justify="start", gap=2),
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
def _(
    LogNorm,
    Normalize,
    allstar,
    cb_cmap,
    cb_col,
    cb_flip,
    cb_label,
    cb_log,
    cb_range,
    colorbar,
    cuts,
    flip_x,
    flip_y,
    log_x,
    log_y,
    matched_cond,
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
    mo.stop(not matched_cond)

    filtered_allstar = allstar.copy()

    if observatory.value == "LCO":
        filtered_allstar = filtered_allstar[filtered_allstar["telescope"].str.startswith("lco")]
    elif observatory.value == "APO":
        filtered_allstar = filtered_allstar[filtered_allstar["telescope"].str.startswith("apo")]

    user_cuts = cuts.value.strip()

    if user_cuts:
        namespace_temp = {col: filtered_allstar[col] for col in filtered_allstar.columns}
        namespace_temp["np"] = np

        mask_cuts = np.ones(len(filtered_allstar), dtype=bool)
        for line in user_cuts.split("\n"):
            line = line.strip()
            if line:
                mask_cuts &= eval(line, {"__builtins__": {}}, namespace_temp)

        filtered_allstar = filtered_allstar[mask_cuts]

    namespace = {col: filtered_allstar[col] for col in filtered_allstar.columns}
    namespace["np"] = np

    hrd_fig, hrd_ax = plt.subplots(figsize=(10,6))


    x_vals = eval(x_col.value, {"__builtins__": {}}, namespace)
    y_vals = eval(y_col.value, {"__builtins__": {}}, namespace)

    filtered_allstar = filtered_allstar.assign(x_vals=x_vals, y_vals=y_vals)

    hrd_fontsize = 16

    if colorbar.value:
        cb_vals = eval(cb_col.value, {"__builtins__": {}}, namespace)
        filtered_allstar = filtered_allstar.assign(cb_vals=cb_vals)

        ix_cb_vals = np.argsort(cb_vals)
        if cb_flip.value:
            ix_cb_vals = np.flip(ix_cb_vals)

        filtered_allstar = filtered_allstar.iloc[ix_cb_vals]

        if cb_log.value:
            norm = LogNorm(*tuple(float(x.strip()) for x in cb_range.value.split(","))) if cb_range.value else LogNorm()
        else:
            norm = Normalize(*tuple(float(x.strip()) for x in cb_range.value.split(","))) if cb_range.value else Normalize()

        sc = hrd_ax.scatter(
            filtered_allstar["x_vals"],
            filtered_allstar["y_vals"],
            c=filtered_allstar["cb_vals"],
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
            filtered_allstar["x_vals"],
            filtered_allstar["y_vals"],
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
    n_stars_hrd = len(np.unique(filtered_allstar['sdss_id']))
    n_spectra_hrd = len(filtered_allstar)
    hrd_tit_str_left = f"{n_stars_hrd:,} stars" if n_stars_hrd > 1 else f"{n_stars_hrd} star"
    hrd_tit_str_right = f"{n_spectra_hrd:,} spectra" if n_spectra_hrd > 1 else f"{n_spectra_hrd} spectrum"
    hrd_ax.set_title(hrd_tit_str_left, loc='left', fontsize=hrd_fontsize)
    hrd_ax.set_title(hrd_tit_str_right, loc='right', fontsize=hrd_fontsize)

    hrd = mo.ui.matplotlib(plt.gca(), debounce=True)

    mo.hstack([hrd], justify="center")
    return filtered_allstar, hrd


@app.cell(hide_code=True)
def _(filtered_allstar, hrd):
    select_mask = hrd.value.get_mask(
        filtered_allstar["x_vals"], filtered_allstar["y_vals"]
    )
    selected_allstar = filtered_allstar[select_mask]
    return (selected_allstar,)


@app.cell(hide_code=True)
def _(mo):
    spec_color = mo.ui.text(value="g_mag - rp_mag", label="color spectra by")
    spec_cmap = mo.ui.text(value="turbo", label="colormap")
    smoothing = mo.ui.checkbox(label="smooth spectra with gaussian filter", value=False)
    spec_ranges = mo.ui.checkbox(label="customize axis bounds")



    pan1_xrange = mo.ui.text(label=r"top", value="3550, 10450")
    pan2_xrange = mo.ui.text(label=r"left", value="4821, 4901")
    pan3_xrange = mo.ui.text(label=r"center", value="6523, 6603")
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
        spec_ranges,
    )


@app.cell(hide_code=True)
def _(mo, selected_allstar, smoothing, spec_cmap, spec_color, spec_ranges):
    if len(selected_allstar):
        spec_color_prompt = mo.hstack([spec_color, spec_cmap, smoothing, spec_ranges], justify="start", gap=2)
    else:
        spec_color_prompt = mo.md('')

    spec_color_prompt
    return


@app.cell(hide_code=True)
def _(np, selected_allstar, spec_color):
    namespace_select = {col: selected_allstar[col] for col in selected_allstar.columns}
    namespace_select["np"] = np

    spec_color_vals = eval(spec_color.value, {"__builtins__": {}}, namespace_select)

    allstar_so = selected_allstar.assign(spec_color_vals=spec_color_vals)

    ix_spec_color_vals = np.argsort(spec_color_vals)
    allstar_so = allstar_so.iloc[ix_spec_color_vals].reset_index(drop=True)
    return (allstar_so,)


@app.cell(hide_code=True)
def _(
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
    if spec_ranges.value:
        pan_bounds = mo.hstack(
            [
                mo.vstack(
                    [
                        mo.md(r"$\mathrm{\lambda}$ ranges:"),
                        mo.md(r"flux/continuum ranges:"),
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
    allstar_so,
    ax_ixs,
    ax_lams,
    flux_over_cont,
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
    selected_allstar,
    smoothing,
    spec_cmap,
):
    if len(selected_allstar):

        fig = plt.figure(figsize=(10, 6), constrained_layout=True)

        gs = fig.add_gridspec(2, 3)

        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[1, 1])
        ax4 = fig.add_subplot(gs[1, 2])

        axes = [ax1, ax2, ax3, ax4]

        spec_color_cmap = plt.get_cmap(spec_cmap.value)

        color_positions = np.linspace(0, 1, len(allstar_so))

        flux_sel = flux_over_cont[allstar_so["nb_flux_arr_ix"]]

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

        specfig = mo.hstack([fig], justify="center")

    else:
        specfig = mo.md("")

    specfig
    return


@app.cell(hide_code=True)
def _(mo, selected_allstar):
    if len(selected_allstar):
        spec_df_display_check = mo.ui.checkbox(label=r"display `mwmAllStar-0.8.1.fits` information for selected subset")
    else:
        spec_df_display_check = mo.md('')

    spec_df_display_check
    return (spec_df_display_check,)


@app.cell(hide_code=True)
def _(allstar_cols, allstar_so, mo, spec_df_display_check):
    if spec_df_display_check.value:
        spec_df = allstar_so[allstar_cols]
    else:
        spec_df = mo.md('')
    spec_df
    return


@app.cell(hide_code=True)
def _(mo):
    outfilename = mo.ui.text(label="output file name:", value="dr20_boss_spectra")
    return (outfilename,)


@app.cell(hide_code=True)
def _(mo, outfilename):
    outfile = f"home/data/{outfilename.value}.h5"

    save_options = ["all uploaded IDs", "subset selected above"]
    save_subset_option = mo.ui.radio(
        options=save_options, value="all uploaded IDs", inline=True
    )

    save_spectra_button = mo.ui.run_button(
        label=f"save spectra to `{outfile}`&nbsp;&nbsp;💾"
    )

    download_md_body = mo.md("After saving, download the spectra from the file tree in your main BinderHub page")
    return download_md_body, outfile, save_spectra_button, save_subset_option


@app.cell(hide_code=True)
def _(
    Path,
    allstar,
    allstar_cols,
    continuum,
    flux,
    h5py,
    ivar,
    matched_cond,
    mo,
    nmf_rectified_model_flux,
    outfile,
    outfilename,
    save_spectra_button,
    save_subset_option,
    selected_allstar,
    wavelength,
):
    if save_spectra_button.value:
        if matched_cond:
            if outfilename.value:
                with mo.status.spinner(f"Saving {outfile}..."):

                    outfile_path = Path(outfile)
                    outfile_path.parent.mkdir(parents=True, exist_ok=True)

                    with h5py.File(outfile_path, "w") as outfile_f:
                        if save_subset_option.value == "all uploaded IDs":
                            allstar_rec = allstar.assign(telescope=lambda df: df["telescope"].astype("S6")).to_records(index=False)

                            outfile_f.create_dataset("/allstar", data=allstar_rec[allstar_cols])
                            outfile_f.create_dataset("/spectra/wavelength", data=wavelength)
                            outfile_f.create_dataset("/spectra/flux", data=flux, compression="gzip")
                            outfile_f.create_dataset("/spectra/ivar", data=ivar, compression="gzip")
                            outfile_f.create_dataset("/spectra/nmf_rectified_model_flux", data=nmf_rectified_model_flux, compression="gzip")
                            outfile_f.create_dataset("/spectra/continuum", data=continuum, compression="gzip")
                            savemessage = mo.md(f"Done saving to `{outfile}` ✅ ")
                        else:
                            if len(selected_allstar):
                                allstar_rec = selected_allstar.assign(telescope=lambda df: df["telescope"].astype("S6")).to_records(index=False)

                                outfile_f.create_dataset("/allstar", data=allstar_rec[allstar_cols])
                                outfile_f.create_dataset("/spectra/wavelength", data=wavelength)
                                outfile_f.create_dataset("/spectra/flux",data=flux[allstar_rec["nb_flux_arr_ix"]], compression="gzip")
                                outfile_f.create_dataset("/spectra/ivar", data=ivar[allstar_rec["nb_flux_arr_ix"]], compression="gzip")
                                outfile_f.create_dataset("/spectra/nmf_rectified_model_flux", data=nmf_rectified_model_flux[allstar_rec["nb_flux_arr_ix"]], compression="gzip")
                                outfile_f.create_dataset("/spectra/continuum", data=continuum[allstar_rec["nb_flux_arr_ix"]], compression="gzip")



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
            savemessage = mo.md(
                '<span style="color: red;">Sorry, no DR20 BOSS spectra for those IDs</span> 😵‍💫'
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
def _(matched_cond, mo, outfilename):
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

    Tip 💡: `allstar`, `wavelength`, `flux`, `ivar`, `nmf_rectified_model_flux`, `continuum` and `flux_over_cont` are already in memory and can be directly accessed in this notebook without saving❗



    By default, `/allstar` comes with 

    `['sdss_id', 'gaia_dr3_source_id', 'ra', 'dec', 'l', 'b', 'plx', 'e_plx', 'pmra', 'e_pmra', 'pmde', 'e_pmde', 'gaia_v_rad', 'gaia_e_v_rad', 'g_mag', 'bp_mag', 'rp_mag', 'j_mag', 'h_mag', 'k_mag', 'w1_mag', 'w2_mag', 'telescope', 'n_good_visits', 'v_rad', 'e_v_rad', 'std_v_rad', 'median_e_v_rad', 'xcsao_teff', 'xcsao_e_teff', 'xcsao_logg', 'xcsao_e_logg', 'xcsao_fe_h', 'xcsao_e_fe_h', 'xcsao_meanrxc', 'snr', 'zwarning_flags', 'nmf_rchi2', 'nmf_flags']`


    The column list may be edited in `allstar_cols` in the code toward the top of the notebook.

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

    if matched_cond:
        acc_out = mo.accordion({'<h2 style="text-align: left; font-weight: bold;"> Access the spectra </h2>': access_md})
    else:
        acc_out = mo.md('')

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
