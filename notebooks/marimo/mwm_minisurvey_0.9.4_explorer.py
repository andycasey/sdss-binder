import marimo

__generated_with = "0.23.3"
app = marimo.App(width="full", app_title="mwm-0.9.4 Explorer")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <div style="display: flex; align-items: center; justify-content: center; position: relative; height: 80px; width: 100%;">
        <img src="https://www.sdss.org/wp-content/uploads/2022/09/sdss-new-logo-72dpi.png"
             style="height: 72px; position: absolute; right: 0;">

        <h1 style="margin: 0; text-align: center; width: 100%;">MWM Mini-Survey 0.9.4 Explorer</h1>
    </div>
    """)
    return


@app.cell
def _():
    import numpy as np
    import pandas as pd
    import h5py
    from io import StringIO
    from scipy.ndimage import gaussian_filter1d
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    from sdss_semaphore.targeting import TargetingFlags

    font_path = "notebooks/static/GoogleSans-Regular.ttf"
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

    block_path = "data/sandbox/minisurvey_block_file.h5"
    return (
        LogNorm,
        Normalize,
        StringIO,
        TargetingFlags,
        block_path,
        gaussian_filter1d,
        h5py,
        np,
        pd,
        plt,
    )


@app.cell
def _(TargetingFlags, block_path, h5py, np, pd):
    block = {}

    with h5py.File(block_path, "r") as block_f:

        block["sdss_id"] = block_f["meta/sdss_id"][:]
        block["telescope"] = block_f["meta/telescope"][:]
        block["summary"] = block_f["summary"][()]
        block["flux"] = block_f["spectra/flux"][:]
        block["continuum"] = block_f["spectra/continuum"][:]
        block["model_flux"] = block_f["spectra/model_flux"][:]
        block["model_continuum"] = block_f["spectra/model_continuum"][:]
        block["wavelength"] = block_f["spectra/wavelength"][()]

    flags = TargetingFlags(block["summary"]["sdss5_target_flags"], sdssc2bv=3)
    block["carton_names"] = flags.get_carton_name()

    summary_colnames = [
        summary_colname
        for summary_colname in block["summary"].dtype.names
        if (len(block["summary"][summary_colname].shape) == 1)
    ]
    summary_df = pd.DataFrame(
        block["summary"].astype(block["summary"].dtype.newbyteorder("="))[
            summary_colnames
        ]
    )
    summary_df["sdss_id"] = block["sdss_id"].astype(
        block["sdss_id"].dtype.newbyteorder("=")
    )
    summary_df["telescope"] = block["telescope"].astype(
        block["telescope"].dtype.newbyteorder("=")
    )
    summary_df["carton_names"] = np.array(
        [", ".join(carts) for carts in block["carton_names"]]
    )
    return block, flags, summary_df


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <h2 style="text-align: left; font-weight: bold;"> Select cartons or upload a list of IDs </h2>
    """)
    return


@app.cell
def _(mo):
    selection_type_options = ["select cartons", "upload a list of IDs"]
    selection_type_radio = mo.ui.radio(
        options=selection_type_options, value="select cartons", inline=True
    )

    selection_type_radio
    return (selection_type_radio,)


@app.cell
def _(flags, mo, selection_type_radio):
    carton_selection_options = sorted(flags.all_carton_names)
    # sorted(set().union(*block["carton_names"]))

    # carton_selection_val = [cart for cart in carton_selection_options if "snc" in cart]
    carton_selection_val = ["mwm_snc_100pc"]

    carton_selection_multiselect = mo.ui.multiselect(
        options=carton_selection_options,
        label="select cartons",
        value=carton_selection_val,
    )

    carton_selection_display = mo.md("")

    if selection_type_radio.value == "select cartons":
        carton_selection_display = mo.hstack(
            [carton_selection_multiselect], justify="start", widths=[5]
        )

    carton_selection_display
    return (carton_selection_multiselect,)


@app.cell
def _(mo, selection_type_radio):
    id_list_options = ["SDSS ID", "Gaia DR3 source ID"]

    id_list_radio = mo.ui.radio(
        options=id_list_options, value="SDSS ID", inline=True, label="choose ID type: "
    )

    id_list_upload_button = mo.ui.file(
        kind="button", label="upload ID list (.txt, .csv)", filetypes=[".txt", ".csv"]
    )

    id_list_display = mo.md("")

    if selection_type_radio.value == "upload a list of IDs":
        id_list_display = mo.hstack(
            [id_list_radio, id_list_upload_button], justify="start", gap=2
        )

    id_list_display
    return id_list_radio, id_list_upload_button


@app.cell
def _(
    StringIO,
    block,
    carton_selection_multiselect,
    id_list_radio,
    id_list_upload_button,
    mo,
    np,
    selection_type_radio,
    summary_df,
):
    no_spectra_callout = mo.md("")

    if selection_type_radio.value == "select cartons":

        selected_cartons = set(carton_selection_multiselect.value)

        block_ix = [
            i
            for i, names in enumerate(block["carton_names"])
            if not selected_cartons.isdisjoint(names)
        ]

        if (len(block_ix) == 0) and (len(carton_selection_multiselect.value)):
            no_spectra_callout = mo.md(
                text=f"carton(s) not in `mwmAllStar-0.9.4.fits` 🤔"
            ).callout(kind="danger")

        elif np.all(np.all(np.isnan(block["flux"][block_ix]), axis=1)) and (
            len(carton_selection_multiselect.value)
        ):
            no_spectra_callout = mo.md(
                text=f"carton(s) in `mwmAllStar-0.9.4.fits` but no corresponding `mwmStar-0.9.4` spectra 🤔"
            ).callout(kind="danger")

        block_ix = np.array(block_ix)

        if not len(block_ix) == 0:

            block_ix = block_ix[~np.all(np.isnan(block["flux"][block_ix]), axis=1)]

    if selection_type_radio.value == "upload a list of IDs":

        upload_file = id_list_upload_button.value

        upload_id_list_str = upload_file[0].contents.decode("utf-8")

        try:

            upload_id = np.loadtxt(StringIO(upload_id_list_str), dtype=np.int64)

        except ValueError:

            upload_id = np.loadtxt(
                StringIO(upload_id_list_str), dtype=np.int64, skiprows=1
            )

        if id_list_radio.value == "SDSS ID":
            block_id = summary_df.sdss_id.values
        else:
            block_id = summary_df.gaia_dr3_source_id.values

        block_ix = np.where(np.isin(block_id, upload_id))[0]

        if len(block_ix) == 0:

            no_spectra_callout = mo.md(
                text=f"{id_list_radio.value}s not in `mwmAllStar-0.9.4.fits` 🤔"
            ).callout(kind="danger")

        elif np.all(np.all(np.isnan(block["flux"][block_ix]), axis=1)):
            no_spectra_callout = mo.md(
                text=f"{id_list_radio.value}s in `mwmAllStar-0.9.4.fits` but no corresponding `mwmStar-0.9.4` spectra 🤔"
            ).callout(kind="danger")

        block_ix = block_ix[~np.all(np.isnan(block["flux"][block_ix]), axis=1)]

    no_spectra_cond = len(block_ix) == 0

    no_spectra_callout
    return block_ix, no_spectra_cond


@app.cell
def _(block_ix, mo, no_spectra_cond, summary_df):
    mo.stop(no_spectra_cond)

    df_match = summary_df.iloc[block_ix].copy()
    df_match["ix_spectrum"] = block_ix
    return (df_match,)


@app.cell
def _(mo, no_spectra_cond, selection_type_radio, summary_df):
    mo.stop(no_spectra_cond)

    if selection_type_radio.value == "select cartons":

        match_message = mo.md(f"""
        Rows in `mwmAllStar-0.9.4` matching to carton(s) are shown below, matched to 
    
        * `astraAllStarASPCAP-0.9.4.fits`,
        * `astraAllStarApogeeNet-0.9.4.fits`,
        * `astraAllStarAstroNNdist-0.9.4.fits`, and 
        * `astraAllStarAstroNN-0.9.4.fits`.
    
        **Astra pipeline column names begin with their pipeline name, e.g. `'teff'` in `astraAllStarASPCAP-0.9.4.fits` is `'aspcap_teff'`**.
    
        """)

    if selection_type_radio.value == "upload a list of IDs":

        match_message = mo.md(f"""
        Rows in `mwmAllStar-0.9.4` matching to uploaded IDs are shown below, matched to 
    
        * `astraAllStarASPCAP-0.9.4.fits`,
        * `astraAllStarApogeeNet-0.9.4.fits`,
        * `astraAllStarAstroNNdist-0.9.4.fits`, and 
        * `astraAllStarAstroNN-0.9.4.fits`.
    
        **Astra pipeline column names begin with their pipeline name, e.g. `'teff'` in `astraAllStarASPCAP-0.9.4.fits` is `'aspcap_teff'`**.
    
        """)

    match_col_sel_val = [
        "sdss_id",
        "telescope",
        "snr",
        "g_mag",
        "aspcap_flag_bad",
        "apogeenet_flag_bad",
        "astroNN_flag_bad",
        "aspcap_teff",
        "apogeenet_teff",
        "astroNN_teff",
        "aspcap_logg",
        "apogeenet_logg",
        "astroNN_logg",
        "aspcap_m_h_atm",
        "aspcap_fe_h",
        "apogeenet_fe_h",
        "astroNN_fe_h",
    ]

    match_col_multiselect = col_selector = mo.ui.multiselect(
        options=summary_df.columns,
        label="select columns to display",
        value=match_col_sel_val,
    )

    mo.vstack([match_message, match_col_multiselect], justify="start")
    return (match_col_multiselect,)


@app.cell
def _(df_match, match_col_multiselect, mo, no_spectra_cond):
    mo.stop(no_spectra_cond)

    df_match[match_col_multiselect.value]
    return


@app.cell
def _(mo, no_spectra_cond):
    mo.stop(no_spectra_cond)
    mo.md(r"""
    <h2 style="text-align: left; font-weight: bold;"> Plot and Select to View Spectra</h2>

    make a scatter plot with any of the columns available below (and above).
    """)
    return


@app.cell
def _(mo, no_spectra_cond, summary_df):
    mo.stop(no_spectra_cond)
    column_md = mo.md(f"`{str(summary_df.columns.to_list())}`")

    all_cols_acc = mo.accordion({"see all available columns": column_md})

    all_cols_acc
    return


@app.cell
def _(mo, no_spectra_cond):
    mo.stop(no_spectra_cond)
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
            mo.hstack(
                [x_range, y_range, log_x, log_y, flip_x, flip_y],
                justify="start",
                gap=2,
            ),
            mo.hstack([x_label, y_label], justify="start", gap=2),
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


@app.cell
def _(colorbar, mo, no_spectra_cond):
    mo.stop(no_spectra_cond)

    if colorbar.value:

        cb_col = mo.ui.text(value="snr", label="value")
        cb_label = mo.ui.text(value="snr", label="label")
        cb_range = mo.ui.text(placeholder="vmin, vmax", label="range")
        cb_cmap = mo.ui.text(value="turbo", label="colormap")
        cb_flip = mo.ui.checkbox(label="flip plotting order")
        cb_log = mo.ui.checkbox(label="log")

        cb_out = mo.hstack(
            [cb_col, cb_label, cb_range, cb_cmap, cb_log, cb_flip],
            justify="start",
            gap=2,
        )

    else:
        cb_out = mo.md("")

    cb_out
    return cb_cmap, cb_col, cb_flip, cb_label, cb_log, cb_range


@app.cell
def _(mo, no_spectra_cond):
    mo.stop(no_spectra_cond)
    mo.md(r"""
    **To view the spectra, make a box selection by clicking and dragging on the plot, or hold `shift` while doing so to make a lasso selection.**
    """)
    return


@app.cell
def _(
    LogNorm,
    Normalize,
    cb_cmap,
    cb_col,
    cb_flip,
    cb_label,
    cb_log,
    cb_range,
    colorbar,
    cuts,
    df_match,
    flip_x,
    flip_y,
    log_x,
    log_y,
    mo,
    no_spectra_cond,
    np,
    observatory,
    plt,
    selection_type_radio,
    x_col,
    x_label,
    x_range,
    y_col,
    y_label,
    y_range,
):
    mo.stop(no_spectra_cond)
    allstar_hrd = df_match.copy().reset_index(drop=True)

    n_spectra_hrd = len(allstar_hrd)

    if observatory.value == "LCO":
        allstar_hrd = allstar_hrd[allstar_hrd["telescope"] == b"lco25m"]
    elif observatory.value == "APO":
        allstar_hrd = allstar_hrd[allstar_hrd["telescope"] == b"apo25m"]

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

    n_spectra_post_cuts = len(allstar_hrd)

    namespace = {col: allstar_hrd[col] for col in allstar_hrd.columns}
    namespace["np"] = np

    hrd_fig, hrd_ax = plt.subplots(figsize=(10, 6), constrained_layout=True)

    x_vals = eval(x_col.value, {"__builtins__": {}}, namespace)
    y_vals = eval(y_col.value, {"__builtins__": {}}, namespace)

    allstar_hrd = allstar_hrd.assign(x_vals=x_vals, y_vals=y_vals)

    hrd_fontsize = 16

    if colorbar.value:
        cb_vals = eval(cb_col.value, {"__builtins__": {}}, namespace)
        allstar_hrd = allstar_hrd.assign(cb_vals=cb_vals)
        allstar_hrd = allstar_hrd.loc[allstar_hrd[cb_vals.name].notna()]
        cb_vals = cb_vals.loc[cb_vals.notna()]

        ix_cb_vals = np.argsort(cb_vals)

        if cb_flip.value:
            ix_cb_vals = np.flip(ix_cb_vals)

        allstar_hrd = allstar_hrd.iloc[ix_cb_vals]

        if cb_log.value:
            norm = (
                LogNorm(*tuple(float(x.strip()) for x in cb_range.value.split(",")))
                if cb_range.value
                else LogNorm()
            )
        else:
            norm = (
                Normalize(*tuple(float(x.strip()) for x in cb_range.value.split(",")))
                if cb_range.value
                else Normalize()
            )

        sc = hrd_ax.scatter(
            allstar_hrd["x_vals"],
            allstar_hrd["y_vals"],
            c=allstar_hrd["cb_vals"],
            s=10,
            edgecolors="gainsboro",
            lw=0.5,
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
            lw=0.5,
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

    hrd = mo.ui.matplotlib(plt.gca(), debounce=True)

    y_lo, y_hi = sorted(hrd_ax.get_ylim())
    x_lo, x_hi = sorted(hrd_ax.get_xlim())

    n_stars_plot_hrd_mask = (
        (allstar_hrd["x_vals"] >= x_lo)
        & (allstar_hrd["x_vals"] <= x_hi)
        & (allstar_hrd["y_vals"] >= y_lo)
        & (allstar_hrd["y_vals"] <= y_hi)
    )

    n_stars_plot_hrd = len(np.unique(allstar_hrd.loc[n_stars_plot_hrd_mask, "sdss_id"]))

    n_spectra_post_cuts_display = mo.md("")

    if selection_type_radio.value == "select cartons":

        n_spectra_display = mo.stat(
            value=n_spectra_hrd,
            label=f"spectra in carton(s)",
        )

        if user_cuts:
            n_spectra_post_cuts_display = mo.stat(
                value=n_spectra_post_cuts, label="spectra in carton(s) after cuts"
            )
    else:

        n_spectra_display = mo.stat(
            value=n_spectra_hrd,
            label=f"spectra matching to uploaded IDs",
        )

        if user_cuts:
            n_spectra_post_cuts_display = mo.stat(
                value=n_spectra_post_cuts, label="spectra after cuts"
            )

    n_spectra_plot_display = mo.stat(
        value=n_stars_plot_hrd,
        label=f"stars shown on scatter plot",
    )

    mo.vstack(
        [
            mo.hstack(
                [
                    n_spectra_display,
                    n_spectra_post_cuts_display,
                    n_spectra_plot_display,
                ],
                justify="center",
                gap=10,
            ),
            mo.hstack([hrd], justify="center"),
        ],
        gap=2,
    )
    return allstar_hrd, hrd


@app.cell
def _(allstar_hrd, hrd, mo, no_spectra_cond):
    mo.stop(no_spectra_cond)
    select_mask = hrd.value.get_mask(allstar_hrd["x_vals"], allstar_hrd["y_vals"])
    selected_allstar = allstar_hrd[select_mask]
    return (selected_allstar,)


@app.cell
def _():
    DEFAULT_RANGES = {
        "pan1_x": (15100, 17000),
        "pan2_x": (15350, 15750),
        "pan3_x": (16000, 16250),
        "pan4_x": (16750, 16900),
        "pan1_y": (0.4, 1.6),
        "pan2_y": (0.4, 1.6),
        "pan3_y": (0.4, 1.6),
        "pan4_y": (0.4, 1.6),
    }
    return (DEFAULT_RANGES,)


@app.cell
def _(mo):
    spec_color = mo.ui.text(value="g_mag - rp_mag", label="color spectra by")
    spec_cmap = mo.ui.text(value="turbo", label="colormap")
    smoothing = mo.ui.checkbox(label="smooth spectra with gaussian filter", value=False)
    spec_ranges = mo.ui.checkbox(label="customize axis bounds")
    spec_int = mo.ui.checkbox(label="interactive zoom")
    spec_continuum_normalize = mo.ui.checkbox(label="continuum normalize spectra")
    return (
        smoothing,
        spec_cmap,
        spec_color,
        spec_continuum_normalize,
        spec_int,
        spec_ranges,
    )


@app.cell
def _(DEFAULT_RANGES, mo, spec_ranges):
    spec_ranges.value

    def _fmt(key):
        lo, hi = DEFAULT_RANGES[key]
        return f"{lo}, {hi}"

    pan1_xrange = mo.ui.text(label=r"top", value=_fmt("pan1_x"))
    pan2_xrange = mo.ui.text(label=r"left", value=_fmt("pan2_x"))
    pan3_xrange = mo.ui.text(label=r"center", value=_fmt("pan3_x"))
    pan4_xrange = mo.ui.text(label=r"right", value=_fmt("pan4_x"))
    pan1_yrange = mo.ui.text(label=r"top", value=_fmt("pan1_y"))
    pan2_yrange = mo.ui.text(label=r"left", value=_fmt("pan2_y"))
    pan3_yrange = mo.ui.text(label=r"center", value=_fmt("pan3_y"))
    pan4_yrange = mo.ui.text(label=r"right", value=_fmt("pan4_y"))
    return (
        pan1_xrange,
        pan1_yrange,
        pan2_xrange,
        pan2_yrange,
        pan3_xrange,
        pan3_yrange,
        pan4_xrange,
        pan4_yrange,
    )


@app.cell
def _(
    mo,
    no_spectra_cond,
    selected_allstar,
    smoothing,
    spec_cmap,
    spec_color,
    spec_continuum_normalize,
    spec_int,
    spec_ranges,
):
    mo.stop(no_spectra_cond)

    if len(selected_allstar):
        spec_color_prompt = mo.hstack(
            [
                spec_color,
                spec_cmap,
                smoothing,
                spec_int,
                spec_ranges,
                spec_continuum_normalize,
            ],
            justify="start",
            gap=2,
        )
    else:
        spec_color_prompt = mo.md("")

    spec_color_prompt
    return


@app.cell
def _(spec_color_vals):
    spec_color_vals.name
    return


@app.cell
def _(np, selected_allstar, spec_color):
    namespace_select = {col: selected_allstar[col] for col in selected_allstar.columns}
    namespace_select["np"] = np

    spec_color_vals = eval(spec_color.value, {"__builtins__": {}}, namespace_select)

    allstar_so = selected_allstar.assign(spec_color_vals=spec_color_vals)

    allstar_so = allstar_so.loc[spec_color_vals.notna()]
    spec_color_vals = spec_color_vals.loc[spec_color_vals.notna()]
    ix_spec_color_vals = np.argsort(spec_color_vals)
    allstar_so = allstar_so.iloc[ix_spec_color_vals].reset_index(drop=True)
    return allstar_so, spec_color_vals


@app.cell
def _(
    mo,
    no_spectra_cond,
    pan1_xrange,
    pan1_yrange,
    pan2_xrange,
    pan2_yrange,
    pan3_xrange,
    pan3_yrange,
    pan4_xrange,
    pan4_yrange,
    spec_continuum_normalize,
    spec_ranges,
):
    mo.stop(no_spectra_cond)
    if spec_ranges.value:
        if spec_continuum_normalize.value:
            ylabelprompt = r"flux / continuum:"
        else:
            ylabelprompt = r"flux / flux(16,000$\AA$)"
        pan_bounds = mo.hstack(
            [
                mo.vstack(
                    [
                        mo.md(r"$\mathrm{\lambda}$:"),
                        mo.md(ylabelprompt),
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
        pan_bounds = mo.md("")

    pan_bounds
    return


@app.cell
def _(mo, no_spectra_cond, smoothing):
    mo.stop(no_spectra_cond)
    smooth_sigma_display = mo.md("")

    if smoothing.value:
        smooth_sigma = mo.ui.number(start=0, stop=20, label="", value=3)
        ss_text = mo.md(r"""
    [sigma](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html)
    """)
        smooth_sigma_display = mo.hstack(
            [ss_text, smooth_sigma], justify="start", gap=2
        )

    smooth_sigma_display
    return (smooth_sigma,)


@app.cell
def _(
    DEFAULT_RANGES,
    block,
    mo,
    no_spectra_cond,
    np,
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
    mo.stop(no_spectra_cond)

    wavelength = block["wavelength"]

    if spec_ranges.value:
        pan1_xmin, pan1_xmax = tuple(
            float(x.strip()) for x in pan1_xrange.value.split(",")
        )
        pan2_xmin, pan2_xmax = tuple(
            float(x.strip()) for x in pan2_xrange.value.split(",")
        )
        pan3_xmin, pan3_xmax = tuple(
            float(x.strip()) for x in pan3_xrange.value.split(",")
        )
        pan4_xmin, pan4_xmax = tuple(
            float(x.strip()) for x in pan4_xrange.value.split(",")
        )
        pan1_ymin, pan1_ymax = tuple(
            float(x.strip()) for x in pan1_yrange.value.split(",")
        )
        pan2_ymin, pan2_ymax = tuple(
            float(x.strip()) for x in pan2_yrange.value.split(",")
        )
        pan3_ymin, pan3_ymax = tuple(
            float(x.strip()) for x in pan3_yrange.value.split(",")
        )
        pan4_ymin, pan4_ymax = tuple(
            float(x.strip()) for x in pan4_yrange.value.split(",")
        )
    else:
        pan1_xmin, pan1_xmax = DEFAULT_RANGES["pan1_x"]
        pan2_xmin, pan2_xmax = DEFAULT_RANGES["pan2_x"]
        pan3_xmin, pan3_xmax = DEFAULT_RANGES["pan3_x"]
        pan4_xmin, pan4_xmax = DEFAULT_RANGES["pan4_x"]
        pan1_ymin, pan1_ymax = DEFAULT_RANGES["pan1_y"]
        pan2_ymin, pan2_ymax = DEFAULT_RANGES["pan2_y"]
        pan3_ymin, pan3_ymax = DEFAULT_RANGES["pan3_y"]
        pan4_ymin, pan4_ymax = DEFAULT_RANGES["pan4_y"]

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


@app.cell
def _(
    allstar_so,
    ax_ixs,
    ax_lams,
    block,
    gaussian_filter1d,
    mo,
    no_spectra_cond,
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
    smooth_sigma,
    smoothing,
    spec_cmap,
    spec_continuum_normalize,
    spec_int,
):
    mo.stop(no_spectra_cond)
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

        if spec_continuum_normalize.value:

            flux_sel = (
                block["flux"][allstar_so["ix_spectrum"]]
                / block["continuum"][allstar_so["ix_spectrum"]]
            )

        else:
            flux_sel = block["flux"][allstar_so["ix_spectrum"]]
            flux_sel = flux_sel / flux_sel[:, 4187][:, None]

        for spec_i in reversed(range(len(allstar_so))):

            color = spec_color_cmap(color_positions[spec_i])

            for ax_i, (ax_ix, spec_ax, ax_lam) in enumerate(zip(ax_ixs, axes, ax_lams)):

                spec_ax.plot(
                    ax_lam,
                    (
                        gaussian_filter1d(
                            flux_sel[spec_i][ax_ix], sigma=smooth_sigma.value
                        )
                        if (smoothing.value)
                        else flux_sel[spec_i][ax_ix]
                    ),
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

        n_stars = len(np.unique(allstar_so["sdss_id"]))
        n_spectra = len(allstar_so)
        tit_str_left = f"{n_stars:,} stars" if n_stars > 1 else f"{n_stars} star"
        tit_str_right = (
            f"{n_spectra:,} spectra" if n_spectra > 1 else f"{n_spectra} spectrum"
        )
        ax1.set_title(tit_str_left, loc="left", fontsize=16)
        ax1.set_title(tit_str_right, loc="right", fontsize=16)
        if spec_continuum_normalize.value:
            spec_ylabel = "flux / continuum"
        else:
            spec_ylabel = r"flux / flux($\lambda=16,000\AA$)"
        fig.supylabel(spec_ylabel, fontsize=14)
        ax3.set_xlabel(r"$\lambda\ (\AA)$", fontsize=14)

        for spec_ax in axes:
            spec_ax.grid(True, which="both", alpha=0.4, zorder=-100)
            spec_ax.tick_params(axis="both", labelsize=12)

        if spec_int.value:
            specfig = mo.vstack(
                [
                    mo.md(
                        "To zoom, select the 🔍 emoji below the figure and click and drag in any of the panels."
                    ),
                    mo.mpl.interactive(fig),
                ],
                gap=2,
            )
        else:
            specfig = mo.hstack([fig], justify="center")
    else:
        specfig = mo.md("")

    specfig
    return


@app.cell
def _(mo, no_spectra_cond, selected_allstar):
    mo.stop(no_spectra_cond)
    if len(selected_allstar):
        spec_df_display_check = mo.ui.checkbox(
            label=r"display table information for selected subset"
        )
    else:
        spec_df_display_check = mo.md("")

    spec_df_display_check
    return (spec_df_display_check,)


@app.cell
def _(mo, no_spectra_cond, spec_df_display_check, summary_df):
    mo.stop(no_spectra_cond)
    if spec_df_display_check.value:
        col_sel_val_ = [
            "sdss_id",
            "telescope",
            "snr",
            "g_mag",
            "aspcap_flag_bad",
            "apogeenet_flag_bad",
            "astroNN_flag_bad",
            "aspcap_teff",
            "apogeenet_teff",
            "astroNN_teff",
            "aspcap_logg",
            "apogeenet_logg",
            "astroNN_logg",
            "aspcap_m_h_atm",
            "aspcap_fe_h",
            "apogeenet_fe_h",
            "astroNN_fe_h",
        ]

        col_selector_ = mo.ui.multiselect(
            options=summary_df.columns,
            label="select columns to display",
            value=col_sel_val_,
        )
    else:
        col_selector_ = mo.md("")
    col_selector_
    return (col_selector_,)


@app.cell
def _(allstar_so, col_selector_, mo, no_spectra_cond, spec_df_display_check):
    mo.stop(no_spectra_cond)
    if spec_df_display_check.value:
        spec_df = allstar_so[col_selector_.value]
    else:
        spec_df = mo.md("")
    spec_df
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    __Please send questions/issues/requests in the #binderhub slack channel or to__ ksharifi1@gsu.edu
    """)
    return


if __name__ == "__main__":
    app.run()
