import marimo

__generated_with = "0.23.3"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import sdss_access.path
    import os
    os.environ['SAS_BASE_DIR'] = '/home/jovyan/data/release'

    # load packages
    from astropy.table import Table
    import numpy as np
    import h5py
    import matplotlib.pylab as plt
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    return Table, go, make_subplots, np, os, sdss_access


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # MWM APOGEE Explorer Notebook

    Please use this notebook to inspect reduced MWM APOGEE spectra and their ASPCAP fits.

    Note: As of 7/2026 this notebook is set to inspect DR19 spectra only
    """)
    return


@app.cell
def _(np):
    elements = ['Al','C_1','C_12_13','Ca','Ce','Co','Cr','Cu','Fe','K','Mg','Mn','N','Na','Nd','Ni','O','P','S','Si','Ti','Ti_2','V']

    def get_apogee_segment_indices():
        segment_log10_wl_start = np.array([4.180476e+00, 4.20051e+00, 4.217064e+00])
        segment_pixels = np.array([3028, 2495, 1991])
        start_indices = np.round((segment_log10_wl_start - 4.179) / 6e-6).astype(int)
        return (start_indices, segment_pixels)

    def get_ferre_pixel_mask():
        mask = np.zeros(8575, dtype=bool)
        for si, p in zip(*get_apogee_segment_indices()):
            mask[si:si+p] = True
        assert mask.sum() == 7514
        return mask

    ferre_mask = get_ferre_pixel_mask()

    mask_dic = {}

    for e in elements:
        ferre_e_mask = np.loadtxt("/home/jovyan/notebooks/marimo/elem_masks/"+e+".mask")
        e_mask = np.zeros(ferre_mask.size, dtype=bool)
        e_mask[ferre_mask] = ferre_e_mask
        mask_dic[e] = e_mask
    return elements, mask_dic


@app.cell
def _(mo):
    sdss_id = mo.ui.number(label="Enter SDSS ID", value=69405638)
    return (sdss_id,)


@app.cell
def _(mo, sdss_id):
    mo.hstack([sdss_id, mo.md(f"SDSS_ID: {sdss_id.value}")])
    return


@app.cell
def _(os, sdss_access, sdss_id):
    sdss_path = sdss_access.path.Path(release='DR19',
                                      preserve_envvars=True)

    path = sdss_path.full('mwmStar', v_astra='0.6.0', sdss_id=sdss_id.value)
    if not os.path.exists(path):
        find_star = False
    else: find_star = True
    return find_star, sdss_path


@app.cell
def _(find_star, mo, sdss_id):
    mo.callout(mo.md(f"No star found with SDSS_ID {sdss_id.value}. Try another star."), kind="danger") if not find_star else mo.md("")
    return


@app.cell
def _(Table, find_star, np, sdss_id, sdss_path):
    if find_star:
        # download obs and fit spectrum
        hdu=3 # first try at APO
        mwmStar = Table.read(sdss_path.full('mwmStar', v_astra='0.6.0', sdss_id=sdss_id.value), hdu=hdu)
        aspcapStar = Table.read(sdss_path.full('astraStarASPCAP', v_astra='0.6.0', sdss_id=sdss_id.value), hdu=hdu)

        if len(mwmStar)==0:
            hdu=4 # then try LCO
            mwmStar = Table.read(sdss_path.full('mwmStar', v_astra='0.6.0', sdss_id=sdss_id.value), hdu=hdu)
            aspcapStar = Table.read(sdss_path.full('astraStarASPCAP', v_astra='0.6.0', sdss_id=sdss_id.value), hdu=hdu)

        # pull out relevant information
        wavelength = aspcapStar["wavelength"][0]
        model = aspcapStar["model_flux"][0]
        continuum = aspcapStar["continuum"][0]
        obs_wl = mwmStar["wavelength"][0]
        obs_flux = mwmStar["flux"][0]
        obs_err = mwmStar["ivar"][0]**(-0.5)
        mask = obs_err < 0.1 * obs_flux

        # continuum is masked around edges of obs wavelength ranges so small regions not normalized
        cont_mask = ~np.ma.getmaskarray(continuum)
    return cont_mask, continuum, model, obs_flux, obs_wl, wavelength


@app.cell
def _(mo):
    obs_spec_check = mo.ui.checkbox(label="Plot Observed Spectrum", value=True)
    fit_spec_check = mo.ui.checkbox(label="Plot Best Fit Spectrum", value=True)

    mo.hstack([obs_spec_check, fit_spec_check])
    return fit_spec_check, obs_spec_check


@app.cell
def _(elements, mo):
    element_checks = mo.ui.dictionary({
        el: mo.ui.checkbox(label=el, value=False) for el in elements
    })

    mo.hstack([
        mo.vstack(list(element_checks.elements.values())[i::4])
        for i in range(4)
    ])
    return (element_checks,)


@app.cell
def _(
    cont_mask,
    continuum,
    element_checks,
    elements,
    find_star,
    fit_spec_check,
    go,
    mask_dic,
    model,
    obs_flux,
    obs_spec_check,
    obs_wl,
    wavelength,
):
    if find_star:

        active_elements = [el for el, checked in element_checks.value.items() if checked]

        _palette = [
            "rgba(251, 201, 123, 128)",
            "rgba(240, 157, 107, 128)",
            "rgba(217, 116, 119, 128)",
            "rgba(177, 83, 136, 128)",
            "rgba(133, 53, 141, 128)",
            "rgba(86, 29, 131, 128)",
            "rgba(66, 42, 131, 128)",
            "rgba(64, 77, 175, 128)",
            "rgba(53, 117, 192, 128)",
            "rgba(50, 158, 186, 128)",
            "rgba(82, 195, 174, 128)",
            "rgba(145, 226, 184, 128)",
        ]

        fig = go.Figure()

        if obs_spec_check.value:
            fig.add_trace(go.Scatter(
                x=obs_wl[cont_mask],
                y=obs_flux[cont_mask]/continuum[cont_mask],
                mode="lines",
                line=dict(color="black", width=1),
                name="Observed",
            ))

        if fit_spec_check.value:
            fig.add_trace(go.Scatter(
                x=wavelength,
                y=model,
                mode="lines",
                line=dict(color="red", width=1, dash="dash"),
                name="Model",
            ))

        for el in active_elements:
            color = _palette[elements.index(el) % len(_palette)]
            el_mask = mask_dic[el]

            in_region = False
            start_wl = None
            for j, (w, flag) in enumerate(zip(wavelength, el_mask)):
                if flag and not in_region:
                    start_wl = w
                    in_region = True
                elif not flag and in_region:
                    fig.add_vrect(
                        x0=start_wl, x1=wavelength[j - 1],
                        fillcolor=color, opacity=.3, layer="below", line_width=0,
                        annotation_text=el, annotation_position="top left",
                        annotation_font_size=11,
                    )
                    in_region = False
            if in_region:
                fig.add_vrect(
                    x0=start_wl, x1=wavelength[-1],
                    fillcolor=color, opacity=.3, layer="below", line_width=0,
                    annotation_text=el, annotation_position="top left",
                    annotation_font_size=11,
                )

        fig.update_layout(
            template="plotly_white",
            xaxis_title="λ [Å]",
            yaxis_title="Normalized Flux",
            height=500,
            margin=dict(l=60, r=20, t=30, b=60),
        )

    fig
    return


@app.cell
def _(elements, mo):
    elements_no_fe = [el for el in elements if el != 'Fe']

    element_select = mo.ui.dropdown(
        options=elements_no_fe,
        value=elements_no_fe[0],
        label="Element",
    )
    element_select
    return (element_select,)


@app.cell
def _(
    continuum,
    element_select,
    find_star,
    go,
    make_subplots,
    mask_dic,
    mo,
    model,
    np,
    obs_flux,
    wavelength,
):
    if find_star:
        el_2 = element_select.value
        el_mask_2 = mask_dic[el_2]

        # find contiguous windows
        windows = []
        in_region_2 = False
        for k, flag_2 in enumerate(el_mask_2):
            if flag_2 and not in_region_2:
                start = k
                in_region_2 = True
            elif not flag_2 and in_region_2:
                windows.append((start, k - 1))
                in_region_2 = False
        if in_region_2:
            windows.append((start, len(el_mask_2) - 1))

        n = len(windows)
        if n == 0:
            mo.md(f"No windows found for **{el_2}**")
        else:
            ncols = 3
            nrows = int(np.ceil(n / ncols))

            fig2 = make_subplots(
                rows=nrows, cols=ncols,
                shared_xaxes=False,
                vertical_spacing=0.14,
                horizontal_spacing=0.08,
                subplot_titles=[f"{el_2} window {i+1}" for i in range(n)],
            )

            for i, (s, e_2) in enumerate(windows):
                row = i // ncols + 1
                col = i % ncols + 1

                # slice to window with small padding
                pad = 5 #max(10, (e_2 - s) // 4)
                sl = slice(max(0, s - pad), min(len(wavelength), e_2 + pad + 1))

                wl_w = wavelength[sl]

                residuals = obs_flux[sl]/continuum[sl] - model[sl]

                showlegend = i == 0

                fig2.add_trace(go.Scatter(
                    x=wavelength[sl], y=obs_flux[sl]/continuum[sl],
                    mode="lines", line=dict(color="black", width=1),
                    name="Observed", legendgroup="obs", showlegend=showlegend,
                ), row=row, col=col)

                fig2.add_trace(go.Scatter(
                    x=wavelength[sl], y=model[sl],
                    mode="lines", line=dict(color="red", width=1, dash="dash"),
                    name="Model", legendgroup="model", showlegend=showlegend,
                ), row=row, col=col)

                fig2.add_trace(go.Scatter(
                    x=wavelength[sl], y=residuals + 0.4,
                    mode="lines", line=dict(color="steelblue", width=1),
                    name="Residuals", legendgroup="resid", showlegend=showlegend,
                ), row=row, col=col)

                # shade the actual element window
                fig2.add_vrect(
                    x0=wavelength[s], x1=wavelength[e_2],
                    fillcolor="rgba(255,200,0,0.15)",
                    opacity=1, layer="below", line_width=0,
                    row=row, col=col,
                )

                # zero line for residuals reference
                fig2.add_hline(y=0.4, line=dict(color="steelblue", width=0.5, dash="dot"),
                              row=row, col=col)

            fig2.update_layout(
                template="plotly_white",
                height=280 * nrows,
                margin=dict(l=50, r=20, t=40, b=40),
                legend=dict(orientation="h", y=1, x=0.5, xanchor="center", yanchor="middle"),
            )

            fig2.update_xaxes(title_text="λ [Å]", title_font_size=11, tickfont_size=10)
            fig2.update_yaxes(title_text="flux", title_font_size=11, tickfont_size=10)

    fig2
    return


if __name__ == "__main__":
    app.run()
