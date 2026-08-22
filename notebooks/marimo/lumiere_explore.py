import marimo

__generated_with = "0.23.3"
app = marimo.App(width="medium")


@app.cell
def _():

    import marimo as mo
    from sdss_binder.samp_widget import SampClient

    client = SampClient(
        subscriptions=["sdss.spectrum.view"],
        client_name="Marimo @ the Flatiron Institute's BinderHub",
        client_description="Plots SDSS-V spectra for the row selected in TOPCAT",
        icon_url="https://www.sdss.org/wp-content/uploads/2022/04/cropped-cropped-site_icon_SDSSV-32x32.png",
    )
    samp = mo.ui.anywidget(client)

    mo.vstack([
        mo.md(
    """# SDSS-V: Lumiere / TOPCAT explorer
    1. Open TOPCAT.
    2. Click the **Connect to hub** button below. TOPCAT will ask for your approval.
    3. Open a compatible FITS file in TOPCAT. 
    4. When you first select a star or row in that file, TOPCAT will ask you to approve or delete an activation action. After approving, any subsequent row selections will show the spectra in Marimo.

    Currently, compatible FITS files are only those that have been produced by Lumiere."""
        ),
        samp
    ])

    return mo, samp


@app.cell(hide_code=True)
def _(mo, samp):
    _msg = samp.value["message"]
    sdss_id = int(_msg["params"].get("sdss_id")) if _msg else None
    mo.md(f"**SDSS ID:** {sdss_id}")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
