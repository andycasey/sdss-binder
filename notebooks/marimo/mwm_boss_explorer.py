import marimo

__generated_with = "0.20.2"
app = marimo.App(width="full")


@app.cell
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
    from numpy.lib import recfunctions as rfn
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
    return (
        LogNorm,
        Normalize,
        Path,
        StringIO,
        fitsio,
        gaussian_filter1d,
        h5py,
        mo,
        np,
        pd,
        plt,
        rfn,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <div style="display: flex; align-items: center; justify-content: center; position: relative; height: 80px; width: 100%;">
        <!-- Logo left-aligned -->
        <img src="https://www.sdss.org/wp-content/uploads/2022/09/sdss-new-logo-72dpi.png"
             style="height: 72px; position: absolute; left: 0;">

        <!-- Header centered -->
        <h1 style="margin: 0; text-align: center; width: 100%;">MWM DR20 BOSS Explorer</h1>
    </div>

    Click `Run all stale cells` in the bottom right, then `Toggle app view`
    """
    )
    return


@app.cell
def _():
    # LIKELY NEED TO CHANGE THE BELOW PATHS
    block_path = "/home/jovyan/data/release/dr20/spectro/astra/0.8.1/spectra/block/mwmStarBlock-0.8.1.h5"
    allstar_path = "/home/jovyan/home/data/mwmAllStar-0.8.1.fits"
    block_meta_arr_path = str((base_dir.parent / "static" / "block_meta.npy").resolve())
    return allstar_path, block_meta_arr_path, block_path


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
        gap=1,
    )
    return (
        id_list_buttons,
        usr_gaia_id_list_upload_button,
        usr_sdss_id_list_upload_button,
    )


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
            "Tip 💡: before uploading a new ID list, make sure to "
            '<span style="color: red;">Click to clear files</span> above'
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
    upload_messages = mo.vstack([tipmessage, toomanymessage, nospectramessage])
    return allstar_rows, matched_cond, upload_messages


@app.cell
def _():
    # to see all columns:

    # with fitsio.FITS(allstar_path) as allstar_cols_f:
    #     print(allstar_cols_f[1].get_colnames())

    allstar_cols = [
        "sdss_id",
        "gaia_dr3_source_id",
        "telescope",
        "n_good_visits",
        "v_rad",
        "e_v_rad",
        "std_v_rad",
        "median_e_v_rad",
        "xcsao_teff",
        "xcsao_e_teff",
        "xcsao_logg",
        "xcsao_e_logg",
        "xcsao_fe_h",
        "xcsao_e_fe_h",
        "xcsao_meanrxc",
        "snr",
        "zwarning_flags",
        "nmf_rchi2",
        "nmf_flags",
        "ra",
        "dec",
        "l",
        "b",
        "plx",
        "e_plx",
        "pmra",
        "e_pmra",
        "pmde",
        "e_pmde",
        "gaia_v_rad",
        "gaia_e_v_rad",
        "g_mag",
        "bp_mag",
        "rp_mag",
        "j_mag",
        "h_mag",
        "k_mag",
        "w1_mag",
        "w2_mag",
    ]
    return (allstar_cols,)


@app.cell(hide_code=True)
def _(
    allstar_cols,
    allstar_path,
    allstar_rows,
    block_meta_arr_path,
    block_path,
    fitsio,
    h5py,
    matched_cond,
    mo,
    np,
    pd,
    rfn,
):
    if matched_cond:

        with mo.status.spinner("Loading spectra") as load_spinner:

            with fitsio.FITS(allstar_path) as f:
                allstar = f[1].read(columns=allstar_cols, rows=allstar_rows)

            block_meta_arr = np.load(block_meta_arr_path)

            load_spinner.update("Mapping allstar rows to block indices...")

            allstar_tel = allstar["telescope"].astype(str)
            block_tel = block_meta_arr["telescope"].astype(str)

            block_index_map = {
                (sid, tel): i
                for i, (sid, tel) in enumerate(
                    zip(block_meta_arr["sdss_id"], block_tel)
                )
            }

            try:
                ix_sorted = np.array(
                    [
                        block_index_map[(sid, tel)]
                        for sid, tel in zip(allstar["sdss_id"], allstar_tel)
                    ],
                    dtype=int,
                )
            except KeyError as e:
                missing = e.args[0]
                raise ValueError(
                    f"Some (sdss_id, telescope) pairs in allstar are not in block_meta_arr: {missing}"
                )

            sort_ix = np.argsort(ix_sorted)
            ix_sorted = ix_sorted[sort_ix]
            allstar = allstar[sort_ix]

            allstar = rfn.append_fields(
                allstar, "ix_in_blockfile", ix_sorted, usemask=False
            )

            allstar = rfn.append_fields(
                allstar, "nb_flux_arr_ix", np.arange(len(allstar)), usemask=False
            )

            load_spinner.update("Reading in spectra...")

            def read_h5_in_chunks(f, dataset, indices, chunk_size=500):
                shape = (len(indices),) + f[dataset].shape[1:]
                out = np.empty(shape, dtype=f[dataset].dtype)
                for start in range(0, len(indices), chunk_size):
                    end = start + chunk_size
                    out[start:end] = f[dataset][indices[start:end]]
                return out

            with h5py.File(block_path, "r") as f:
                flux = read_h5_in_chunks(f, "boss/spectra/flux", ix_sorted)
                ivar = read_h5_in_chunks(f, "boss/spectra/ivar", ix_sorted)
                nmf_rectified_model_flux = read_h5_in_chunks(
                    f, "boss/spectra/nmf_rectified_model_flux", ix_sorted
                )
                continuum = read_h5_in_chunks(f, "boss/spectra/continuum", ix_sorted)

            flux_over_cont = flux / continuum

            load_spinner.update("Spectra loaded ✅")

        topmessage = mo.md(
            "Uploaded IDs matched to `mwmAllStar-0.8.1.fits` are shown below"
        )
        outdf = pd.DataFrame(allstar[allstar_cols])
        valuetipmessage = mo.md(
            "Tip 💡: `['sdss_id', 'telescope']` value pairs are unique across rows"
        )
        outmessage = mo.vstack([topmessage, outdf, valuetipmessage], gap=2)
    else:
        outmessage = mo.md("")
    return (
        allstar,
        continuum,
        flux,
        flux_over_cont,
        ivar,
        nmf_rectified_model_flux,
        outmessage,
    )


@app.cell(hide_code=True)
def _(block_path, h5py):
    with h5py.File(block_path, "r") as wavelength_block_f:
        wavelength = wavelength_block_f["boss/spectra/wavelength"][()]
    return (wavelength,)


@app.cell(hide_code=True)
def _(id_list_buttons, mo, outmessage, upload_messages):
    upload_md = mo.vstack([id_list_buttons, upload_messages, outmessage], gap=1)
    return (upload_md,)


@app.cell(hide_code=True)
def _(mo, upload_md):
    mo.accordion(
        {
            '<h2 style="text-align: left; font-weight: bold;"> Upload list of SDSS or <i>Gaia</i> DR3 source IDs </h2>': upload_md,
        }
    )
    return


@app.cell(hide_code=True)
def _(mo):
    figsize = mo.ui.text(value="8,6", label="fig size")
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
    return (
        colorbar,
        cuts,
        figsize,
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

        # cb_params = mo.vstack([mo.hstack([cb_col, cb_label, cb_range], justify='start', gap=2),
        # mo.hstack([cb_cmap, cb_log, cb_flip], justify='start', gap=2)])
        cb_params = mo.hstack(
            [cb_col, cb_label, cb_range, cb_cmap, cb_log, cb_flip],
            justify="start",
            gap=2,
        )
    else:
        cb_params = mo.md("")
    return cb_cmap, cb_col, cb_flip, cb_label, cb_log, cb_params, cb_range


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
    figsize,
    flip_x,
    flip_y,
    log_x,
    log_y,
    matched_cond,
    mo,
    np,
    observatory,
    plt,
    rfn,
    x_col,
    x_label,
    x_range,
    y_col,
    y_label,
    y_range,
):
    if matched_cond:
        filtered_allstar = allstar

        if observatory.value == "LCO":
            obs_mask = np.char.startswith(
                filtered_allstar["telescope"].astype(str), "lco"
            )
            filtered_allstar = filtered_allstar[obs_mask]

        if observatory.value == "APO":
            obs_mask = np.char.startswith(
                filtered_allstar["telescope"].astype(str), "apo"
            )
            filtered_allstar = filtered_allstar[obs_mask]

        user_cuts = cuts.value.strip()

        if user_cuts:
            namespace_temp = {
                name: filtered_allstar[name] for name in filtered_allstar.dtype.names
            }
            namespace_temp["np"] = np

            mask_cuts = np.ones(len(filtered_allstar), dtype=bool)

            for line in user_cuts.split("\n"):
                line = line.strip()
                if line:
                    mask_cuts &= eval(line, {"__builtins__": {}}, namespace_temp)

            filtered_allstar = filtered_allstar[mask_cuts]

        namespace = {
            name: filtered_allstar[name] for name in filtered_allstar.dtype.names
        }
        namespace["np"] = np

        figsizeval = tuple(float(x.strip()) for x in figsize.value.split(","))

        hrd_fig, hrd_ax = plt.subplots(figsize=figsizeval)

        x_vals = eval(x_col.value, {"__builtins__": {}}, namespace)
        y_vals = eval(y_col.value, {"__builtins__": {}}, namespace)

        filtered_allstar = rfn.append_fields(
            filtered_allstar, "x_vals", x_vals, usemask=False
        )
        filtered_allstar = rfn.append_fields(
            filtered_allstar, "y_vals", y_vals, usemask=False
        )

        hrd_fontsize = 14

        if colorbar.value:
            cb_vals = eval(cb_col.value, {"__builtins__": {}}, namespace)

            filtered_allstar = rfn.append_fields(
                filtered_allstar, "cb_vals", cb_vals, usemask=False
            )

            ix_cb_vals = np.argsort(cb_vals)

            if cb_flip.value:
                ix_cb_vals = np.flip(ix_cb_vals)

            filtered_allstar = filtered_allstar[ix_cb_vals]

            if cb_log.value:
                if cb_range.value:
                    vmin, vmax = tuple(
                        float(x.strip()) for x in cb_range.value.split(",")
                    )
                    norm = LogNorm(vmin=vmin, vmax=vmax)
                else:
                    norm = LogNorm()
            else:
                if cb_range.value:
                    vmin, vmax = tuple(
                        float(x.strip()) for x in cb_range.value.split(",")
                    )
                    norm = Normalize(vmin=vmin, vmax=vmax)
                else:
                    norm = Normalize()

            sc = hrd_ax.scatter(
                filtered_allstar["x_vals"],
                filtered_allstar["y_vals"],
                c=filtered_allstar["cb_vals"],
                s=20,
                edgecolors="gainsboro",
                lw=0.5,
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
                s=20,
                edgecolors="gainsboro",
                lw=1,
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

        hrd = mo.ui.matplotlib(plt.gca())

        hrd_widget = mo.hstack([hrd], justify="center")
    return filtered_allstar, hrd, hrd_widget


@app.cell(hide_code=True)
def _(mo):
    spec_color = mo.ui.text(value="g_mag - rp_mag", label="color spectra by")
    spec_cmap = mo.ui.text(value="turbo", label="colormap")
    smoothing = mo.ui.checkbox(label="smooth spectra with gaussian filter", value=False)
    spec_ranges = mo.ui.checkbox(label="customize axis bounds")

    spec_color_prompt = mo.hstack(
        [spec_color, spec_cmap, smoothing, spec_ranges], justify="start", gap=2
    )
    return smoothing, spec_cmap, spec_color, spec_color_prompt, spec_ranges


@app.cell(hide_code=True)
def _(filtered_allstar, hrd, np, rfn, spec_color):
    if ("hrd" in globals()) and (hrd.value):
        mask = hrd.value.get_mask(
            filtered_allstar["x_vals"], filtered_allstar["y_vals"]
        )
        selected_allstar = filtered_allstar[mask]

        namespace_select = {
            name: selected_allstar[name] for name in selected_allstar.dtype.names
        }
        namespace_select["np"] = np
        spec_color_vals = eval(spec_color.value, {"__builtins__": {}}, namespace_select)

        selected_allstar = rfn.append_fields(
            selected_allstar, "spec_color_vals", spec_color_vals, usemask=False
        )

        ix_spec_color_vals = np.argsort(spec_color_vals)

        selected_allstar = selected_allstar[ix_spec_color_vals]
    return (selected_allstar,)


@app.cell(hide_code=True)
def _(mo, spec_ranges):
    pan1_xrange = mo.ui.text(label=r"top", value="3500, 10500")
    pan2_xrange = mo.ui.text(label=r"left", value="4821, 4901")
    pan3_xrange = mo.ui.text(label=r"center", value="6523, 6603")
    pan4_xrange = mo.ui.text(label=r"right", value="8460, 8700")

    pan1_yrange = mo.ui.text(label=r"top", placeholder="flux min, flux max")
    pan2_yrange = mo.ui.text(label=r"left", placeholder="flux min, flux max")
    pan3_yrange = mo.ui.text(label=r"center", placeholder="flux min, flux max")
    pan4_yrange = mo.ui.text(label=r"right", placeholder="flux min, flux max")

    if not spec_ranges.value:
        pan_bounds = mo.md("")
    else:
        pan_bounds = mo.hstack(
            [
                mo.vstack(
                    [
                        mo.md(r"panel $\mathrm{\lambda}$ ranges:"),
                        mo.md(r"panel flux ranges:"),
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
    return (
        pan1_xrange,
        pan1_yrange,
        pan2_xrange,
        pan2_yrange,
        pan3_xrange,
        pan3_yrange,
        pan4_xrange,
        pan4_yrange,
        pan_bounds,
    )


@app.cell(hide_code=True)
def _(np, pan1_xrange, pan2_xrange, pan3_xrange, pan4_xrange, wavelength):
    pan1_xmin, pan1_xmax = tuple(float(x.strip()) for x in pan1_xrange.value.split(","))
    pan2_xmin, pan2_xmax = tuple(float(x.strip()) for x in pan2_xrange.value.split(","))
    pan3_xmin, pan3_xmax = tuple(float(x.strip()) for x in pan3_xrange.value.split(","))
    pan4_xmin, pan4_xmax = tuple(float(x.strip()) for x in pan4_xrange.value.split(","))

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
        pan2_xmax,
        pan2_xmin,
        pan3_xmax,
        pan3_xmin,
        pan4_xmax,
        pan4_xmin,
    )


@app.cell(hide_code=True)
def _(
    ax_ixs,
    ax_lams,
    flux_over_cont,
    gaussian_filter1d,
    mo,
    np,
    pan1_xmax,
    pan1_xmin,
    pan1_yrange,
    pan2_xmax,
    pan2_xmin,
    pan2_yrange,
    pan3_xmax,
    pan3_xmin,
    pan3_yrange,
    pan4_xmax,
    pan4_xmin,
    pan4_yrange,
    plt,
    selected_allstar,
    smoothing,
    spec_cmap,
):
    if "selected_allstar" in globals():
        fig = plt.figure(figsize=(10, 6), constrained_layout=True)

        gs = fig.add_gridspec(2, 3)

        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[1, 1])
        ax4 = fig.add_subplot(gs[1, 2])

        axes = [ax1, ax2, ax3, ax4]

        spec_color_cmap = plt.get_cmap(spec_cmap.value)

        color_positions = np.linspace(0, 1, len(selected_allstar))

        flux_sel = flux_over_cont[selected_allstar["nb_flux_arr_ix"]]

        for spec_i in reversed(range(len(selected_allstar))):

            color = spec_color_cmap(color_positions[spec_i])

            for ax_i, (ax_ix, spec_ax, ax_lam) in enumerate(zip(ax_ixs, axes, ax_lams)):

                spec_ax.plot(
                    ax_lam,
                    (
                        gaussian_filter1d(flux_sel[spec_i][ax_ix], sigma=5)
                        if (smoothing.value and ax_i == 0)
                        else flux_sel[spec_i][ax_ix]
                    ),
                    lw=0.5,
                    color=color,
                )

        ax1.set_xlim(pan1_xmin, pan1_xmax)
        ax2.set_xlim(pan2_xmin, pan2_xmax)
        ax3.set_xlim(pan3_xmin, pan3_xmax)
        ax4.set_xlim(pan4_xmin, pan4_xmax)

        if pan1_yrange.value:
            pan1_ymin, pan1_ymax = tuple(
                float(x.strip()) for x in pan1_yrange.value.split(",")
            )
            ax1.set_ylim(pan1_ymin, pan1_ymax)
        if pan2_yrange.value:
            pan2_ymin, pan2_ymax = tuple(
                float(x.strip()) for x in pan2_yrange.value.split(",")
            )
            ax2.set_ylim(pan2_ymin, pan2_ymax)
        if pan3_yrange.value:
            pan3_ymin, pan3_ymax = tuple(
                float(x.strip()) for x in pan3_yrange.value.split(",")
            )
            ax3.set_ylim(pan3_ymin, pan3_ymax)
        if pan4_yrange.value:
            pan4_ymin, pan4_ymax = tuple(
                float(x.strip()) for x in pan4_yrange.value.split(",")
            )
            ax4.set_ylim(pan4_ymin, pan4_ymax)

        tit_str = f"selected {len(np.unique(selected_allstar['sdss_id']))} stars"
        ax1.set_title(tit_str, loc="left", fontsize=12)

        fig.supylabel("flux / continuum", fontsize=12)
        ax3.set_xlabel(r"$\lambda\ (\AA)$", fontsize=12)

        for spec_ax in axes:
            spec_ax.grid(True, which="both", alpha=0.4, zorder=-100)
            spec_ax.tick_params(axis="both", labelsize=12)

        if len(selected_allstar):
            specfig = mo.hstack([fig], justify="center")
        else:
            specfig = mo.md("")
    return (specfig,)


@app.cell(hide_code=True)
def _(
    cb_params,
    colorbar,
    cuts,
    figsize,
    flip_x,
    flip_y,
    hrd_widget,
    log_x,
    log_y,
    matched_cond,
    mo,
    observatory,
    x_col,
    x_label,
    x_range,
    y_col,
    y_label,
    y_range,
):
    hrd_fig_and_params = (
        mo.vstack(
            [
                mo.hstack([x_col, y_col, cuts, observatory], justify="start", gap=2),
                mo.hstack([figsize, x_label, y_label], justify="start", gap=2),
                mo.hstack(
                    [x_range, y_range, log_x, log_y, flip_x, flip_y],
                    justify="start",
                    gap=2,
                ),
                mo.hstack([colorbar]),
                mo.hstack([cb_params]),
                hrd_widget,
                # # mo.hstack([observatory], justify="start", gap=1)
            ],
            gap=2,
        )
        if matched_cond
        else mo.md("")
    )
    return (hrd_fig_and_params,)


@app.cell(hide_code=True)
def _(mo):
    spec_df_display_check = mo.ui.checkbox(
        label=r"display `mwmAllStar-0.8.1.fits` information for selected subset"
    )
    return (spec_df_display_check,)


@app.cell(hide_code=True)
def _(
    allstar_cols,
    mo,
    pd,
    selected_allstar,
    spec_color_prompt,
    spec_df_display_check,
):
    if "selected_allstar" in globals():
        spec_color_prompt_display = (
            spec_color_prompt if len(selected_allstar) else mo.md("")
        )

        if spec_df_display_check.value:
            spec_df = pd.DataFrame(selected_allstar[allstar_cols])

        else:
            spec_df = mo.md("")
    else:
        spec_df = mo.md("")
    return spec_color_prompt_display, spec_df


@app.cell(hide_code=True)
def _(
    mo,
    pan_bounds,
    selected_allstar,
    spec_color_prompt_display,
    spec_df,
    spec_df_display_check,
    specfig,
):
    if "selected_allstar" in globals():

        spec_color_pan_bounds_display = mo.vstack(
            [spec_color_prompt_display, pan_bounds], gap=2
        )

        specfig_display = mo.hstack([specfig], justify="center")
        spec_df_display = (
            mo.vstack([spec_df_display_check, spec_df], gap=2)
            if len(selected_allstar)
            else mo.md("")
        )

    else:
        spec_color_pan_bounds_display = mo.md("")
        specfig_display = mo.md("")
        spec_df_display = mo.md("")
    return spec_color_pan_bounds_display, spec_df_display, specfig_display


@app.cell(hide_code=True)
def _(
    hrd_fig_and_params,
    matched_cond,
    mo,
    spec_color_pan_bounds_display,
    spec_df_display,
    specfig_display,
):
    visualize_md_body1 = mo.md(
        r"Make a scatter plot with any of the `mwmAllStar-0.8.1.fits` columns shown above. Make a box selection by clicking and dragging on the plot. Hold `shift` while doing so to make a lasso selection."
    )
    visualize_md = (
        mo.vstack(
            [
                visualize_md_body1,
                hrd_fig_and_params,
                spec_color_pan_bounds_display,
                specfig_display,
                spec_df_display,
            ],
            gap=2,
        )
        if matched_cond
        else mo.md("")
    )
    return (visualize_md,)


@app.cell(hide_code=True)
def _(mo, visualize_md):
    mo.accordion(
        {
            '<h2 style="text-align: left; font-weight: bold;"> Visualization and Selection </h2>': visualize_md,
        }
    )
    return


@app.cell(hide_code=True)
def _(mo):
    outfilename = mo.ui.text(label="output file name:", value="dr20_boss_spectra")
    return (outfilename,)


@app.cell(hide_code=True)
def _(mo, outfilename):
    outfile = f"/home/jovyan/home/data/{outfilename.value}.h5"

    save_options = ["all uploaded IDs", "subset selected above"]
    save_subset_option = mo.ui.radio(
        options=save_options, value="all uploaded IDs", inline=True
    )

    save_spectra_button = mo.ui.run_button(
        label=f"save spectra to `{outfile}`&nbsp;&nbsp;💾"
    )

    download_md_body1 = mo.md(
        r"""
    Once saved, the spectra can be downloaded from the file tree in your main binder page:  https://binder.flatironinstitute.org/hub/user/YOUR_GMAIL@gmail.com/lab
    """
    )
    return download_md_body1, outfile, save_spectra_button, save_subset_option


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
    pd,
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
                            allstar_rec = (
                                pd.DataFrame(allstar[allstar_cols])
                                .assign(
                                    telescope=lambda df: df["telescope"].astype("S6")
                                )
                                .to_records(index=False)
                            )
                            outfile_f.create_dataset("/allstar", data=allstar_rec)
                            outfile_f.create_dataset(
                                "/spectra/wavelength", data=wavelength
                            )
                            outfile_f.create_dataset(
                                "/spectra/flux", data=flux, compression="gzip"
                            )
                            outfile_f.create_dataset(
                                "/spectra/ivar", data=ivar, compression="gzip"
                            )
                            outfile_f.create_dataset(
                                "/spectra/nmf_rectified_model_flux",
                                data=nmf_rectified_model_flux,
                                compression="gzip",
                            )
                            outfile_f.create_dataset(
                                "/spectra/continuum", data=continuum, compression="gzip"
                            )
                            savemessage = mo.md(f"Done saving to `{outfile}` ✅ ")
                        else:
                            if "selected_allstar" in globals():
                                allstar_rec = (
                                    pd.DataFrame(selected_allstar)
                                    .assign(
                                        telescope=lambda df: df["telescope"].astype(
                                            "S6"
                                        )
                                    )
                                    .to_records(index=False)
                                )
                                outfile_f.create_dataset(
                                    "/allstar", data=allstar_rec[allstar_cols]
                                )
                                outfile_f.create_dataset(
                                    "/spectra/wavelength",
                                    data=wavelength[allstar_rec["nb_flux_arr_ix"]],
                                )
                                outfile_f.create_dataset(
                                    "/spectra/flux",
                                    data=flux[allstar_rec["nb_flux_arr_ix"]],
                                    compression="gzip",
                                )
                                outfile_f.create_dataset(
                                    "/spectra/ivar",
                                    data=ivar[allstar_rec["nb_flux_arr_ix"]],
                                    compression="gzip",
                                )
                                outfile_f.create_dataset(
                                    "/spectra/nmf_rectified_model_flux",
                                    data=nmf_rectified_model_flux[
                                        allstar_rec["nb_flux_arr_ix"]
                                    ],
                                    compression="gzip",
                                )
                                outfile_f.create_dataset(
                                    "/spectra/continuum",
                                    data=continuum[allstar_rec["nb_flux_arr_ix"]],
                                    compression="gzip",
                                )
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
    download_md_body1,
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
            download_md_body1,
            savemessage,
        ],
        gap=1,
    )
    return (download_md,)


@app.cell(hide_code=True)
def _(download_md, mo):
    mo.accordion(
        {
            '<h2 style="text-align: left; font-weight: bold;"> Download the spectra </h2>': download_md,
        }
    )
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
    access_md_code1 = mo.md(
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
    Code snippets can be pasted right into this notebook after creating `dr20_boss_spectra.h5` if you wish 🤠

    Tip 💡: `allstar`, `wavelength`, `flux`, `ivar`, `nmf_rectified_model_flux`, `continuum` and `flux_over_cont` are already in memory and can be directly accessed in this notebook without saving❗



    `/allstar` by default comes with

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

    access_md = mo.vstack(
        [access_md_filename, access_md_body1, access_md_code1, access_md_body2], gap=0
    )
    return (access_md,)


@app.cell(hide_code=True)
def _(access_md, mo):
    mo.accordion(
        {
            '<h2 style="text-align: left; font-weight: bold;"> Access the spectra </h2>': access_md,
        }
    )
    return


if __name__ == "__main__":
    app.run()
