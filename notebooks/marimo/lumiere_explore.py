import marimo

__generated_with = "0.23.3"
app = marimo.App(width="medium")


@app.cell
def _():

    import marimo as mo
    import h5py as h5
    import numpy as np
    import plotly.graph_objects as go
    from sdss_binder.samp_widget import SampClient

    return SampClient, go, h5, mo, np


@app.cell
def _(mo, samp):
    # Cell A — instructions. Reads .value, renders NO widget.
    connected = samp.value.get("registered", False)
    mo.md("# SDSS-V: Lumiere / TOPCAT explorer") if connected else mo.md(
    """# SDSS-V: Lumiere / TOPCAT explorer
    1. Open TOPCAT.
    2. Click **Connect to hub** below. TOPCAT will ask for your approval.
    3. Open a compatible FITS file in TOPCAT.
    4. On first row selection, approve the activation action.

    Compatible FITS files are currently only those produced by Lumiere."""
    )
    return


@app.cell
def _(SampClient, mo):

    client = SampClient(
        subscriptions=["sdss.spectrum.view"],
        client_name="Marimo @ the Flatiron Institute's BinderHub",
        client_description="Plots SDSS-V spectra for the row selected in TOPCAT",
        icon_url="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-32x32.png",
    )
    samp = mo.ui.anywidget(client)
    samp
    return (samp,)


@app.cell(hide_code=True)
def _(LUMIERE_BASE_DIR, go, h5, mo, np, row_for, samp):
    # Cell — load + render

    from plotly.subplots import make_subplots

    WAVELENGTH = 10 ** (4.179 + 6e-6 * np.arange(8575))
    CHIP_NAMES = ("Blue", "Green", "Red")
    CHIP_LIMITS = ((15150, 15800.0), (15860.0, 16425.0), (16480.0, 16875.0))

    message = samp.value.get("message")
    mo.stop(message is None, mo.md("*Waiting for a row selection in TOPCAT…*"))

    _raw_id = message["params"].get("sdss_id")
    spectrum_path = message["params"].get("spectrum_path")
    results_path = message["params"].get("results_path")
    mo.stop(_raw_id is None or spectrum_path is None or results_path is None,
            mo.md("⚠️ Message missing `sdss_id`, `spectrum_path`, or `results_path`."))

    sdss_id = int(_raw_id)
    index = row_for(spectrum_path, sdss_id)          # cached after first hit on this file
    mo.stop(index is None, mo.md(f"⚠️ `sdss_id` **{sdss_id}** not in `{spectrum_path}`."))

    with h5.File(LUMIERE_BASE_DIR + spectrum_path) as s:
        flux = s["spectra/flux"][index]
        ivar = s["spectra/ivar"][index]

    with h5.File(LUMIERE_BASE_DIR + results_path) as fp:
        model_index = np.where(fp["sdss_id"][:] == sdss_id)[0][0]
        model_flux = fp["model_flux"][model_index]
        #continuum = fp["continuum"][index]

    mo.stop(flux.size != WAVELENGTH.size,
            mo.md(f"⚠️ Spectrum has {flux.size} pixels, expected {WAVELENGTH.size}."))
    mo.stop(model_flux.size != flux.size,
            mo.md(f"⚠️ Model has {model_flux.size} pixels, spectrum has {flux.size}."))

    # Mask bad pixels so gaps break the line instead of dropping to zero
    _good = ivar > 0
    _y = np.where(_good, flux, np.nan)
    _sigma = np.where(_good, 1.0 / np.sqrt(np.where(_good, ivar, 1.0)), np.nan)

    fig = make_subplots(rows=3, cols=1, vertical_spacing=0.08)

    for _row, (_lo_w, _hi_w) in enumerate(CHIP_LIMITS, start=1):
        _m = (WAVELENGTH >= _lo_w) & (WAVELENGTH <= _hi_w)
        _xs = WAVELENGTH[_m]
        _first = _row == 1                      # legend entries only once

        fig.add_trace(go.Scattergl(x=_xs, y=_y[_m] - _sigma[_m], mode="lines",
                                   line=dict(width=0), hoverinfo="skip",
                                   showlegend=False), row=_row, col=1)
        fig.add_trace(go.Scattergl(x=_xs, y=_y[_m] + _sigma[_m], mode="lines",
                                   line=dict(width=0), fill="tonexty",
                                   fillcolor="rgba(100,100,100,0.25)", hoverinfo="skip",
                                   showlegend=False), row=_row, col=1)
        fig.add_trace(go.Scattergl(x=_xs, y=_y[_m], mode="lines", name="flux",
                                   legendgroup="flux", showlegend=_first,
                                   line=dict(width=1, color="#000000")), row=_row, col=1)
        fig.add_trace(go.Scattergl(x=_xs, y=model_flux[_m], mode="lines", name="model",
                                   legendgroup="model", showlegend=_first,
                                   line=dict(width=1, color="#d62728")), row=_row, col=1)

        fig.update_xaxes(range=[_lo_w, _hi_w], row=_row, col=1)
        fig.update_yaxes(title_text="Flux", row=_row, col=1)

    fig.update_xaxes(title_text="Wavelength [Å]", row=3, col=1)
    fig.update_layout(
        height=760, margin=dict(l=60, r=20, t=60, b=50),
        hovermode="x unified", template="plotly_white", dragmode="pan",
        uirevision="keep",          # ← add this
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.04, xanchor="right", x=1.0),
    )

    # Shared y-range across chips; spectra have spikes, so clip to something readable
    #_finite = np.concatenate([_y[np.isfinite(_y)], continuum[np.isfinite(continuum)]])
    #if _finite.size:
    #    _lo, _hi = np.percentile(_finite, [0.5, 99.5])
    #    _pad = 0.1 * (_hi - _lo)
    #    fig.update_yaxes(range=[_lo - _pad, _hi + _pad])

    mo.vstack([
        mo.md(f"## SDSS ID {sdss_id}"),
        mo.md(f"*from `{spectrum_path}`, row {index} — via {message.get('sender_name', 'SAMP')}*"),
        fig,
    ])
    return


@app.cell
def _(h5, mo, np):
    # Cell — lookup helpers
    LUMIERE_BASE_DIR = "data/work/acasey/lumiere/"

    @mo.cache
    def sdss_id_index(spectrum_path):
        """Load and cache the sdss_id column for one Lumiere file."""
        with h5.File(LUMIERE_BASE_DIR + spectrum_path) as s:
            ids = s["meta/sdss_id"][:]
        # searchsorted needs sorted keys; only pay for argsort if we have to
        order = None if np.all(ids[:-1] <= ids[1:]) else np.argsort(ids)
        return ids, order


    def row_for(spectrum_path, sdss_id):
        """Row index of sdss_id in that file, or None if absent."""
        ids, order = sdss_id_index(spectrum_path)
        keys = ids if order is None else ids[order]
        i = int(np.searchsorted(keys, sdss_id))
        if i >= keys.size or keys[i] != sdss_id:
            return None
        return i if order is None else int(order[i])

    return LUMIERE_BASE_DIR, row_for


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()