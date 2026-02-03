import marimo

__generated_with = "0.10.0"
app = marimo.App(width="full", app_title="SDSS-V Sky Map")


# =============================================================================
# UI CELLS (displayed in order)
# =============================================================================

@app.cell(hide_code=True)
def _(mo):
    return mo.md(
        """
        # SDSS-V Sky Coverage Map

        Interactive all-sky visualization of sources with APOGEE spectra from ApogeeReduction.jl.
        
        - **Projection:** Mollweide (Galactic coordinates)
        - **Interaction:** Pan by dragging, zoom with scroll wheel
        - **Click on a point** to see detailed source information below
        """
    )


@app.cell(hide_code=True)
def _(mo, n_sources, n_total_exposures):
    return mo.md(f"**Loaded:** {n_sources:,} unique sources from {n_total_exposures:,} total exposures")


@app.cell(hide_code=True)
def _(display_note, mo):
    return mo.md(f"_{display_note}_")


@app.cell(hide_code=True)
def _(aladin_widget, mo):
    # Display the Aladin widget
    return mo.ui.anywidget(aladin_widget)


@app.cell(hide_code=True)
def _(mo, selected_source_input):
    return mo.hstack([
        mo.md("**Enter SDSS ID:**"),
        selected_source_input,
        mo.md("_(press Enter to load source details)_")
    ], justify="start", gap=1)


@app.cell(hide_code=True)
def _(mo, n_sources):
    # Note about subsampling (only shown if subsampled)
    _max_pts = 100000
    if n_sources > _max_pts:
        return mo.md(f"""
        > **Note:** The map shows a random subset of {_max_pts:,} sources for browser performance. 
        > If you enter an SDSS ID that isn't in the displayed subset, it won't appear on the map, 
        > but the details will still load below.
        """)


@app.cell(hide_code=True)
def _(mo, selected_sdss_id):
    # Prompt to select a source
    if selected_sdss_id is None:
        return mo.md("_Enter an SDSS ID above to see source details._")


@app.cell(hide_code=True)
def _(mo, n_exposures, selected_sdss_id, source_data):
    # Source summary header
    mo.stop(selected_sdss_id is None)
    
    _ra = source_data.get('gaia_ra')
    _dec = source_data.get('gaia_dec')
    _links = []
    if _ra is not None and _dec is not None:
        _links.append(f"[SIMBAD](https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={_ra}+{_dec}&Radius=2&Radius.unit=arcsec)")
        _links.append(f"[VizieR](https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-c={_ra}+{_dec}&-c.rs=2)")
    _links_str = " · ".join(_links) if _links else ""
    
    return mo.md(f"""
    ---
    
    ## Source: SDSS ID {selected_sdss_id}
    
    **{n_exposures} exposure{'s' if n_exposures != 1 else ''}** · {_links_str}
    """)


@app.cell(hide_code=True)
def _(mo, selected_sdss_id, source_data):
    # Identifiers table
    mo.stop(selected_sdss_id is None)
    
    return mo.md(f"""
    ### Identifiers
    
    | Field | Value |
    |-------|-------|
    | **SDSS ID** | {source_data.get('sdss_id', '—')} |
    | **Catalog ID** | {source_data.get('catalogid', '—')} |
    | **Gaia DR3** | {source_data.get('gaia_dr3_source', '—')} |
    | **TIC v8** | {source_data.get('tic_v8', '—')} |
    | **2MASS** | {source_data.get('twomass_designation', '—')} |
    """)


@app.cell(hide_code=True)
def _(format_value, mo, selected_sdss_id, source_data):
    # Astrometry and photometry
    mo.stop(selected_sdss_id is None)
    
    return mo.md(f"""
    ### Astrometry & Photometry
    
    | Property | Value | Property | Value |
    |----------|-------|----------|-------|
    | **RA** | {format_value(source_data.get('gaia_ra'), 6)}° | **Dec** | {format_value(source_data.get('gaia_dec'), 6)}° |
    | **Parallax** | {format_value(source_data.get('gaia_parallax'), 3)} ± {format_value(source_data.get('gaia_parallax_error'), 3)} mas | **Distance** | {format_value(source_data.get('gaia_distance_gspphot'), 1)} pc |
    | **PM RA** | {format_value(source_data.get('gaia_pmra'), 3)} ± {format_value(source_data.get('gaia_pmra_error'), 3)} mas/yr | **PM Dec** | {format_value(source_data.get('gaia_pmdec'), 3)} ± {format_value(source_data.get('gaia_pmdec_error'), 3)} mas/yr |
    | **G** | {format_value(source_data.get('gaia_phot_g_mean_mag'), 3)} mag | **BP** | {format_value(source_data.get('gaia_phot_bp_mean_mag'), 3)} mag |
    | **RP** | {format_value(source_data.get('gaia_phot_rp_mean_mag'), 3)} mag | **J** | {format_value(source_data.get('twomass_j_m'), 3)} mag |
    | **H** | {format_value(source_data.get('twomass_h_m'), 3)} mag | **K** | {format_value(source_data.get('twomass_k_m'), 3)} mag |
    """)


@app.cell(hide_code=True)
def _(format_value, mo, selected_sdss_id, source_data):
    # Stellar parameters
    mo.stop(selected_sdss_id is None)
    
    return mo.md(f"""
    ### Stellar Parameters (Gaia GSP-Phot)
    
    | Parameter | Value | Parameter | Value |
    |-----------|-------|-----------|-------|
    | **Teff** | {format_value(source_data.get('gaia_teff_gspphot'), 0)} K | **log g** | {format_value(source_data.get('gaia_logg_gspphot'), 2)} |
    | **[M/H]** | {format_value(source_data.get('gaia_mh_gspphot'), 2)} | **RV (Gaia)** | {format_value(source_data.get('gaia_radial_velocity'), 2)} ± {format_value(source_data.get('gaia_radial_velocity_error'), 2)} km/s |
    """)


@app.cell(hide_code=True)
def _(exposure_data, mo, np, selected_sdss_id):
    # Exposure table
    mo.stop(selected_sdss_id is None)
    
    import pandas as _pd
    
    _exp_df_data = {'#': list(range(1, len(exposure_data.get('mjd', [])) + 1))}
    _field_labels = {
        'mjd': 'MJD', 'observatory': 'Obs', 'exposure': 'Exposure',
        'snr': 'SNR', 'v_rad': 'v_rad (km/s)', 'e_v_rad': 'σ_v_rad'
    }
    
    for _field, _label in _field_labels.items():
        if _field in exposure_data:
            _vals = exposure_data[_field]
            if hasattr(_vals, 'tolist'):
                _vals = _vals.tolist()
            if len(_vals) > 0 and isinstance(_vals[0], (float, np.floating)):
                _vals = [round(v, 2) if np.isfinite(v) else None for v in _vals]
            _exp_df_data[_label] = _vals
    
    return mo.vstack([
        mo.md("### Exposures"),
        mo.ui.table(_pd.DataFrame(_exp_df_data), selection=None, pagination=True, page_size=10)
    ])


@app.cell(hide_code=True)
def _(mo):
    return mo.md(
        """
        ---
        
        ### Map Controls
        
        | Action | How |
        |--------|-----|
        | **Pan** | Click and drag |
        | **Zoom** | Mouse scroll wheel |
        | **Select source** | Click on a point |
        | **Toggle layers** | Layers control (top right) |
        | **Full screen** | Fullscreen button |
        """
    )


# =============================================================================
# BACKEND CODE (imports, data loading, helper functions)
# =============================================================================

@app.cell(hide_code=True)
def _():
    import marimo as mo
    import h5py as h5
    import numpy as np
    import os
    from pathlib import Path
    from ipyaladin import Aladin
    from astropy.table import Table
    from astropy.coordinates import SkyCoord, Angle
    import astropy.units as u
    return Aladin, Angle, Path, SkyCoord, Table, h5, mo, np, os, u


@app.cell(hide_code=True)
def _(Path, os):
    # Configuration
    DATA_DIR = Path(os.environ.get("DATA_DIR", os.path.expanduser("~/../sdssv/ceph/work/scifest/0.2.0")))
    EXPOSURES_FILE = DATA_DIR / "exposures.h5"
    
    SOURCE_IDENTIFIER_FIELDS = (
        "sdss_id", "catalogid", "gaia_dr3_source", "tic_v8", "twomass_designation",
    )
    SOURCE_DISPLAY_FIELDS = (
        "gaia_ra", "gaia_dec", "gaia_parallax", "gaia_parallax_error",
        "gaia_pmra", "gaia_pmra_error", "gaia_pmdec", "gaia_pmdec_error",
        "gaia_phot_g_mean_mag", "gaia_phot_bp_mean_mag", "gaia_phot_rp_mean_mag",
        "gaia_radial_velocity", "gaia_radial_velocity_error",
        "gaia_teff_gspphot", "gaia_logg_gspphot", "gaia_mh_gspphot",
        "gaia_distance_gspphot", "twomass_j_m", "twomass_h_m", "twomass_k_m",
    )
    return DATA_DIR, EXPOSURES_FILE, SOURCE_DISPLAY_FIELDS, SOURCE_IDENTIFIER_FIELDS


@app.cell(hide_code=True)
def _():
    # Helper function for formatting values
    def format_value(val, decimals=4):
        if val is None:
            return "—"
        if isinstance(val, float):
            return f"{val:.0f}" if decimals == 0 else f"{val:.{decimals}f}"
        return str(val)
    return (format_value,)


@app.cell(hide_code=True)
def _(EXPOSURES_FILE, h5, mo, np):
    # Load unique sources from exposures file
    mo.stop(not EXPOSURES_FILE.exists(), mo.md(f"❌ **Error:** Exposures file not found at `{EXPOSURES_FILE}`"))
    
    with h5.File(EXPOSURES_FILE, "r", locking=False) as _h5file:
        sdss_id_all = _h5file["sdss_id"][:]
        _ra_all = _h5file["gaia_ra"][:]
        _dec_all = _h5file["gaia_dec"][:]
    
    _unique_sdss_ids, _unique_indices = np.unique(sdss_id_all, return_index=True)
    
    _sdss_id = sdss_id_all[_unique_indices]
    ra = _ra_all[_unique_indices]
    dec = _dec_all[_unique_indices]
    
    _valid_mask = np.isfinite(ra) & np.isfinite(dec)
    sdss_id = _sdss_id[_valid_mask]
    ra = ra[_valid_mask]
    dec = dec[_valid_mask]
    
    _unique_ids, _counts = np.unique(sdss_id_all, return_counts=True)
    exposure_counts = dict(zip(_unique_ids, _counts))
    
    n_sources = len(sdss_id)
    n_total_exposures = len(sdss_id_all)
    return dec, exposure_counts, n_sources, n_total_exposures, ra, sdss_id, sdss_id_all


@app.cell(hide_code=True)
def _(Table, dec, n_sources, np, ra, sdss_id, u):
    # Prepare display data with subsampling if needed
    _MAX_POINTS = 100000
    
    if n_sources > _MAX_POINTS:
        _rng = np.random.default_rng(42)
        _sample_idx = _rng.choice(n_sources, _MAX_POINTS, replace=False)
        display_ra = ra[_sample_idx]
        display_dec = dec[_sample_idx]
        display_sdss_id = sdss_id[_sample_idx]
        display_note = f"Showing {_MAX_POINTS:,} of {n_sources:,} sources (random sample)"
    else:
        display_ra = ra
        display_dec = dec
        display_sdss_id = sdss_id
        display_note = f"Showing all {n_sources:,} sources"
    
    # Create astropy table for ipyaladin
    catalog_table = Table({
        'ra': display_ra * u.deg,
        'dec': display_dec * u.deg,
        'sdss_id': display_sdss_id,
    })
    
    return catalog_table, display_note


@app.cell(hide_code=True)
def _(Aladin, Angle, catalog_table):
    # Create ipyaladin widget
    aladin_widget = Aladin(
        fov=Angle(180, "deg"),
        target="0 0 Galactic",
        survey="P/2MASS/color",
        coo_frame="Galactic",
        projection="MOL",
        show_coo_grid=True,
        show_fullscreen_control=True,
        show_layers_control=True,
        show_zoom_control=True,
        show_settings_control=True,
        height=600,
    )
    
    # Add the catalog
    aladin_widget.add_table(
        catalog_table,
        name="SDSS-V Sources",
        color="#ff6b6b",
        source_size=8,
        shape="circle",
    )
    
    return (aladin_widget,)


@app.cell(hide_code=True)
def _(mo):
    # UI state: text input for source selection
    selected_source_input = mo.ui.text(value="", label="", full_width=False)
    return (selected_source_input,)


@app.cell(hide_code=True)
def _(selected_source_input):
    # Parse selected SDSS ID from input
    selected_sdss_id = None
    if selected_source_input.value:
        try:
            selected_sdss_id = int(selected_source_input.value.strip())
        except (ValueError, AttributeError):
            pass
    return (selected_sdss_id,)


@app.cell(hide_code=True)
def _(
    EXPOSURES_FILE, SOURCE_DISPLAY_FIELDS, SOURCE_IDENTIFIER_FIELDS,
    h5, mo, np, sdss_id_all, selected_sdss_id
):
    # Load source details when selected
    mo.stop(selected_sdss_id is None)
    
    _idx = np.where(sdss_id_all == selected_sdss_id)[0]
    
    if len(_idx) == 0:
        mo.stop(True, mo.md(f"⚠️ **SDSS ID {selected_sdss_id} not found in the dataset.**"))
    
    with h5.File(EXPOSURES_FILE, "r", locking=False) as _src_file:
        source_data = {}
        for _f in SOURCE_IDENTIFIER_FIELDS:
            if _f in _src_file:
                _v = _src_file[_f][_idx[0]]
                if isinstance(_v, (bytes, np.bytes_)):
                    _v = _v.decode('utf-8').strip()
                source_data[_f] = _v
        
        for _f in SOURCE_DISPLAY_FIELDS:
            if _f in _src_file:
                _v = _src_file[_f][_idx[0]]
                if isinstance(_v, (np.floating, float)) and not np.isfinite(_v):
                    _v = None
                source_data[_f] = _v
        
        exposure_data = {}
        for _f in ['mjd', 'snr', 'v_rad', 'e_v_rad', 'observatory', 'exposure']:
            if _f in _src_file:
                _vs = _src_file[_f][_idx]
                if len(_vs) > 0 and isinstance(_vs[0], bytes):
                    _vs = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in _vs]
                exposure_data[_f] = _vs
    
    n_exposures = len(_idx)
    return exposure_data, n_exposures, source_data


if __name__ == "__main__":
    app.run()
