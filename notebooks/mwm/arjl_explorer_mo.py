import marimo

__generated_with = "0.10.0"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
        # SDSS-V Explorer

        Interactive notebook for exploring APOGEE spectra reduced by ApogeeReduction.jl

        **Features:**
        - Search sources by SDSS ID or other identifiers (with autocomplete)
        - View source information (catalog IDs, astrometry, photometry, stellar parameters)
        - Radial velocity plot with phase folding
        - Interactive spectrum visualization (show/hide individual spectra)
        - Exposure data table
        """
    )
    return


@app.cell
def _():
    # Required imports
    import h5py as h5
    import numpy as np
    import pickle
    import os
    from pathlib import Path
    from collections import defaultdict
    import pandas as pd
    import plotly.graph_objects as go
    return Path, defaultdict, go, h5, np, os, pd, pickle


@app.cell
def _(Path, os):
    # Configuration - modify DATA_DIR if needed
    DATA_DIR = Path(os.environ.get("DATA_DIR", os.path.expanduser("~/../sdssv/ceph/work/scifest/0.2.0")))
    EXPOSURES_FILE = DATA_DIR / "exposures.h5"
    SPECTRA_FILE = DATA_DIR / "arMADGICS_out_x_starLines_v0.h5"
    IDENTIFIERS_FILE = DATA_DIR / "identifiers.pkl"
    return DATA_DIR, EXPOSURES_FILE, IDENTIFIERS_FILE, SPECTRA_FILE


@app.cell
def _(DATA_DIR, EXPOSURES_FILE, IDENTIFIERS_FILE, SPECTRA_FILE, mo):
    mo.md(
        f"""
        ## Configuration

        - **Data directory:** `{DATA_DIR}`
        - **Exposures file exists:** {EXPOSURES_FILE.exists()}
        - **Spectra file exists:** {SPECTRA_FILE.exists()}
        - **Identifiers file exists:** {IDENTIFIERS_FILE.exists()}
        """
    )
    return


@app.cell
def _():
    # Field name definitions
    SOURCE_IDENTIFIER_FIELD_NAMES = (
        "sdss_id", "catalogid", "allstar_dr17_synspec_rev1", "gaia_dr3_source",
        "sdss_dr17_specobj", "tic_v8", "twomass_designation",
    )

    SOURCE_FIELD_NAMES = (
        *SOURCE_IDENTIFIER_FIELD_NAMES,
        "version_id", "lead", "allwise", "catwise", "catwise2020", "sdss_dr13_photoobj",
        "gaia_dr2_source", "glimpse", "guvcat", "panstarrs1", "ps1_g18", "skymapper_dr2",
        "supercosmos", "twomass_psc", "tycho2", "unwise",
        # Gaia DR3 Data
        "gaia_source_id", "gaia_ra", "gaia_ra_error", "gaia_dec", "gaia_dec_error",
        "gaia_parallax", "gaia_parallax_error", "gaia_pm", "gaia_pmra", "gaia_pmra_error",
        "gaia_pmdec", "gaia_pmdec_error", "gaia_ruwe", "gaia_duplicated_source",
        "gaia_phot_g_mean_mag", "gaia_phot_bp_mean_mag", "gaia_phot_rp_mean_mag",
        "gaia_phot_bp_rp_excess_factor", "gaia_radial_velocity", "gaia_radial_velocity_error",
        "gaia_rv_nb_transits", "gaia_rv_nb_deblended_transits", "gaia_rv_visibility_periods_used",
        "gaia_rv_expected_sig_to_noise", "gaia_rv_renormalised_gof", "gaia_rv_chisq_pvalue",
        "gaia_rv_time_duration", "gaia_rv_amplitude_robust", "gaia_rv_template_teff",
        "gaia_rv_template_logg", "gaia_rv_template_fe_h", "gaia_rv_atm_param_origin",
        "gaia_vbroad", "gaia_vbroad_error", "gaia_vbroad_nb_transits",
        "gaia_grvs_mag", "gaia_grvs_mag_error", "gaia_grvs_mag_nb_transits",
        "gaia_rvs_spec_sig_to_noise", "gaia_teff_gspphot", "gaia_logg_gspphot",
        "gaia_mh_gspphot", "gaia_distance_gspphot", "gaia_azero_gspphot",
        "gaia_ag_gspphot", "gaia_ebpminrp_gspphot",
        # 2MASS
        "twomass_j_m", "twomass_j_cmsig", "twomass_j_msigcom", "twomass_j_snr",
        "twomass_h_m", "twomass_h_cmsig", "twomass_h_msigcom", "twomass_h_snr",
        "twomass_k_m", "twomass_k_cmsig", "twomass_k_msigcom", "twomass_k_snr",
        "twomass_ph_qual", "twomass_rd_flg", "twomass_bl_flg", "twomass_cc_flg",
        # ToO
        "too", "too_id", "too_program",
        # Targeting
        "sdss5_target_flags",
    )
    return SOURCE_FIELD_NAMES, SOURCE_IDENTIFIER_FIELD_NAMES


@app.cell
def _(IDENTIFIERS_FILE, np, pickle):
    # Helper functions
    def decode_identifier(identifier):
        """Convert identifier to string, decoding bytes if necessary."""
        if isinstance(identifier, bytes):
            return identifier.decode('utf-8').strip()
        elif isinstance(identifier, np.bytes_):
            return identifier.decode('utf-8').strip()
        return str(identifier).strip()

    # Load identifier index
    if IDENTIFIERS_FILE.exists():
        print(f"Loading identifiers from {IDENTIFIERS_FILE}...")
        with open(IDENTIFIERS_FILE, "rb") as f:
            lookup_identifier, identifier_list, trigram_index = pickle.load(f)
        print(f"Loaded {len(lookup_identifier)} identifiers")
    else:
        raise FileNotFoundError(
            f"Identifiers file not found: {IDENTIFIERS_FILE}\n"
            "Please run the web app first to generate the identifiers.pkl file, "
            "or run the Jupyter notebook version which can build it."
        )
    return (
        decode_identifier,
        identifier_list,
        lookup_identifier,
        trigram_index,
    )


@app.cell
def _(identifier_list, np, trigram_index):
    # Autocomplete function
    def autocomplete(query, limit=20):
        """Return identifiers matching the query string using trigram index."""
        if not query or len(query) < 2:
            return []

        q_lower = query.lower()

        if len(q_lower) >= 3:
            trigrams = [q_lower[i:i+3] for i in range(len(q_lower) - 2)]
            first_trigram = trigrams[0]
            if first_trigram not in trigram_index:
                return []

            candidate_ids = trigram_index[first_trigram].copy()
            for trigram in trigrams[1:]:
                if trigram in trigram_index:
                    candidate_ids &= trigram_index[trigram]
                else:
                    return []
                if not candidate_ids:
                    return []

            matches = []
            for ident_id in candidate_ids:
                ident = identifier_list[ident_id]
                if q_lower in ident.lower():
                    matches.append(ident)
                    if len(matches) >= limit:
                        break
        else:
            bigram_key = f"_b_{q_lower}"
            if bigram_key in trigram_index:
                candidate_ids = trigram_index[bigram_key]
                matches = [identifier_list[i] for i in list(candidate_ids)[:limit]]
            else:
                matches = []

        return matches


    def convert_value(val):
        """Convert numpy types for display."""
        if isinstance(val, np.bool_):
            return bool(val)
        elif isinstance(val, (np.integer, np.floating)):
            v = val.item()
            if isinstance(v, float) and not np.isfinite(v):
                return None
            return v
        elif isinstance(val, bytes):
            return val.decode('utf-8')
        elif isinstance(val, np.ndarray):
            return [convert_value(x) for x in val]
        return val


    def format_value(val, decimals=4):
        """Format a value for display."""
        if val is None:
            return 'N/A'
        if isinstance(val, float):
            if not np.isfinite(val):
                return 'N/A'
            if val == int(val):
                return str(int(val))
            return f"{val:.{decimals}f}"
        if isinstance(val, str) and val.strip() == '':
            return 'N/A'
        return str(val)


    def get_wavelength_grid(num_pixels=8575):
        """Generate APOGEE wavelength grid."""
        crval = 4.179  # log10(lambda) starting value
        cdelt = 6e-6   # log10(lambda) step
        return np.array([10 ** (crval + i * cdelt) for i in range(num_pixels)])
    return autocomplete, convert_value, format_value, get_wavelength_grid


@app.cell
def _(
    EXPOSURES_FILE,
    SOURCE_FIELD_NAMES,
    SPECTRA_FILE,
    convert_value,
    h5,
    lookup_identifier,
):
    # Data query functions
    def query_source(identifier):
        """Query data for a given source identifier."""
        indices = None
        for t in (str, int):
            try:
                indices = lookup_identifier[t(identifier)]
                break
            except:
                continue

        if indices is None:
            raise ValueError(f"Identifier '{identifier}' not found")

        with h5.File(EXPOSURES_FILE, "r", locking=False) as fp:
            meta = {}
            for k in fp.keys():
                if k in SOURCE_FIELD_NAMES:
                    value = convert_value(fp[k][indices[0]])
                else:
                    value = convert_value(fp[k][indices])
                meta[k] = value

        return {
            "sdss_id": meta["sdss_id"],
            "num_exposures": len(indices),
            "indices": indices.tolist(),
            "exposures": meta,
        }


    def query_spectrum(index):
        """Query a single spectrum by its index."""
        with h5.File(SPECTRA_FILE, "r", locking=False) as fp:
            spectrum = fp["x_starLines_v0"][index]
        return convert_value(spectrum)
    return query_source, query_spectrum


@app.cell
def _(mo):
    mo.md("## Search for a Source")
    return


@app.cell
def _(autocomplete, mo):
    # Search input with autocomplete
    search_input = mo.ui.text(
        placeholder="Enter SDSS ID or identifier (e.g., 66899483)",
        label="Search",
        full_width=True
    )

    # Get suggestions based on current input
    def get_suggestions(query):
        if len(query) >= 2:
            return autocomplete(query, limit=20)
        return []

    search_input
    return get_suggestions, search_input


@app.cell
def _(get_suggestions, mo, search_input):
    # Show autocomplete suggestions
    suggestions = get_suggestions(search_input.value)

    if suggestions and search_input.value:
        suggestion_buttons = mo.ui.dropdown(
            options={s: s for s in suggestions},
            label="Suggestions (select one):",
            full_width=True
        )
        mo.vstack([
            mo.md(f"**Found {len(suggestions)} matches:**"),
            suggestion_buttons
        ])
    else:
        suggestion_buttons = None
        mo.md("_Type at least 2 characters to see suggestions_")
    return suggestion_buttons, suggestions


@app.cell
def _(mo, search_input, suggestion_buttons):
    # Determine which identifier to use
    selected_identifier = None
    if suggestion_buttons is not None and suggestion_buttons.value:
        selected_identifier = suggestion_buttons.value
    elif search_input.value and len(search_input.value) >= 2:
        selected_identifier = search_input.value

    # Load button
    load_button = mo.ui.run_button(label="Load Source")

    mo.hstack([
        mo.md(f"**Selected:** `{selected_identifier}`" if selected_identifier else "_No identifier selected_"),
        load_button
    ], justify="start", gap=2)
    return load_button, selected_identifier


@app.cell
def _(load_button, mo, query_source, selected_identifier):
    # Load source data when button is clicked
    mo.stop(not load_button.value or not selected_identifier)

    try:
        source_data = query_source(selected_identifier)
        load_error = None
    except Exception as e:
        source_data = None
        load_error = str(e)

    if load_error:
        mo.md(f"**Error:** {load_error}")
    else:
        mo.md(f"**Loaded:** SDSS ID {source_data['sdss_id']} with {source_data['num_exposures']} exposures")
    return load_error, source_data


@app.cell
def _(format_value, load_error, mo, source_data):
    # Source information display
    mo.stop(source_data is None or load_error is not None)

    e = source_data['exposures']

    # Summary card
    summary_html = f"""
    ## Source Summary

    | Property | Value | Property | Value |
    |----------|-------|----------|-------|
    | **SDSS ID** | {source_data['sdss_id']} | **Exposures** | {source_data['num_exposures']} |
    | **RA** | {format_value(e.get('gaia_ra'), 6)}° | **Dec** | {format_value(e.get('gaia_dec'), 6)}° |
    | **G mag** | {format_value(e.get('gaia_phot_g_mean_mag'), 3)} | **Teff** | {format_value(e.get('gaia_teff_gspphot'), 0)} K |
    | **Parallax** | {format_value(e.get('gaia_parallax'), 4)} mas | **Distance** | {format_value(e.get('gaia_distance_gspphot'), 1)} pc |
    """

    # External links
    ra, dec = e.get('gaia_ra'), e.get('gaia_dec')
    if ra is not None and dec is not None:
        links_html = f"""
    ### External Links

    [SIMBAD](https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={ra}+{dec}&Radius=2&Radius.unit=arcsec) |
    [VizieR](https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-c={ra}+{dec}&-c.rs=2) |
    [SDSS IR Spectrum](https://dr17.sdss.org/infrared/spectrum/search?ra={ra}&dec={dec})
    """
    else:
        links_html = ""

    mo.md(summary_html + links_html)
    return dec, e, links_html, ra, summary_html


@app.cell
def _(e, format_value, load_error, mo, source_data):
    # Detailed source information in tabs
    mo.stop(source_data is None or load_error is not None)

    identifiers_md = f"""
    | Field | Value | Field | Value |
    |-------|-------|-------|-------|
    | SDSS ID | {e.get('sdss_id')} | Catalog ID | {e.get('catalogid')} |
    | Gaia DR3 Source | {e.get('gaia_dr3_source')} | TIC v8 | {e.get('tic_v8')} |
    | 2MASS Designation | {e.get('twomass_designation')} | Lead Catalog | {e.get('lead')} |
    """

    astrometry_md = f"""
    | Field | Value | Field | Value |
    |-------|-------|-------|-------|
    | RA (deg) | {format_value(e.get('gaia_ra'), 6)} | Dec (deg) | {format_value(e.get('gaia_dec'), 6)} |
    | Parallax (mas) | {format_value(e.get('gaia_parallax'), 4)} ± {format_value(e.get('gaia_parallax_error'), 4)} | Distance (pc) | {format_value(e.get('gaia_distance_gspphot'), 1)} |
    | PM RA (mas/yr) | {format_value(e.get('gaia_pmra'), 3)} ± {format_value(e.get('gaia_pmra_error'), 3)} | PM Dec (mas/yr) | {format_value(e.get('gaia_pmdec'), 3)} ± {format_value(e.get('gaia_pmdec_error'), 3)} |
    | RUWE | {format_value(e.get('gaia_ruwe'), 3)} | Duplicated | {'Yes' if e.get('gaia_duplicated_source') else 'No'} |
    """

    photometry_md = f"""
    | Band | Magnitude | Band | Magnitude |
    |------|-----------|------|-----------|
    | G | {format_value(e.get('gaia_phot_g_mean_mag'), 3)} | BP | {format_value(e.get('gaia_phot_bp_mean_mag'), 3)} |
    | RP | {format_value(e.get('gaia_phot_rp_mean_mag'), 3)} | BP-RP Excess | {format_value(e.get('gaia_phot_bp_rp_excess_factor'), 3)} |
    | J | {format_value(e.get('twomass_j_m'), 3)} ± {format_value(e.get('twomass_j_msigcom'), 3)} | H | {format_value(e.get('twomass_h_m'), 3)} ± {format_value(e.get('twomass_h_msigcom'), 3)} |
    | K | {format_value(e.get('twomass_k_m'), 3)} ± {format_value(e.get('twomass_k_msigcom'), 3)} | 2MASS Qual | {e.get('twomass_ph_qual')} |
    """

    stellar_params_md = f"""
    | Parameter | Value | Parameter | Value |
    |-----------|-------|-----------|-------|
    | Teff (K) | {format_value(e.get('gaia_teff_gspphot'), 0)} | log g | {format_value(e.get('gaia_logg_gspphot'), 2)} |
    | [M/H] | {format_value(e.get('gaia_mh_gspphot'), 2)} | A_G (mag) | {format_value(e.get('gaia_ag_gspphot'), 3)} |
    | E(BP-RP) | {format_value(e.get('gaia_ebpminrp_gspphot'), 3)} | | |
    """

    gaia_rv_md = f"""
    | Parameter | Value | Parameter | Value |
    |-----------|-------|-----------|-------|
    | RV (km/s) | {format_value(e.get('gaia_radial_velocity'), 2)} ± {format_value(e.get('gaia_radial_velocity_error'), 2)} | RV Transits | {e.get('gaia_rv_nb_transits')} |
    | Template Teff | {format_value(e.get('gaia_rv_template_teff'), 0)} | Template log g | {format_value(e.get('gaia_rv_template_logg'), 2)} |
    | Vbroad (km/s) | {format_value(e.get('gaia_vbroad'), 1)} | | |
    """

    mo.ui.tabs({
        "Identifiers": mo.md(identifiers_md),
        "Astrometry": mo.md(astrometry_md),
        "Photometry": mo.md(photometry_md),
        "Stellar Parameters": mo.md(stellar_params_md),
        "Gaia RV": mo.md(gaia_rv_md),
    })
    return (
        astrometry_md,
        gaia_rv_md,
        identifiers_md,
        photometry_md,
        stellar_params_md,
    )


@app.cell
def _(load_error, mo, source_data):
    mo.stop(source_data is None or load_error is not None)
    mo.md("## Radial Velocity")
    return


@app.cell
def _(load_error, mo, source_data):
    # RV plot controls
    mo.stop(source_data is None or load_error is not None)

    phase_fold_toggle = mo.ui.checkbox(label="Phase Fold", value=False)
    period_input = mo.ui.number(
        label="Period (days)",
        start=0.001,
        step=0.001,
        value=None
    )
    t0_input = mo.ui.number(
        label="T₀ (MJD)",
        value=None
    )

    mo.hstack([phase_fold_toggle, period_input, t0_input], justify="start", gap=2)
    return period_input, phase_fold_toggle, t0_input


@app.cell
def _(
    go,
    load_error,
    mo,
    period_input,
    phase_fold_toggle,
    source_data,
    t0_input,
):
    # RV Plot
    mo.stop(source_data is None or load_error is not None)

    _e = source_data['exposures']
    _mjd = _e.get('mjd_mid_exposure', [])
    _v_rad = _e.get('v_rad', [])
    _e_v_rad = _e.get('e_v_rad', [])
    _snr = _e.get('snr', [])

    # Filter valid data
    _valid_data = []
    for _i in range(len(_mjd)):
        if _mjd[_i] is not None and _v_rad[_i] is not None:
            _valid_data.append({
                'index': _i,
                'mjd': _mjd[_i],
                'v_rad': _v_rad[_i],
                'e_v_rad': _e_v_rad[_i] if _e_v_rad and _i < len(_e_v_rad) and _e_v_rad[_i] is not None else 0,
                'snr': _snr[_i] if _snr and _i < len(_snr) and _snr[_i] is not None else 0
            })

    if _valid_data:
        # Phase folding
        _phase_fold = phase_fold_toggle.value
        _period = period_input.value
        _t0 = t0_input.value

        if _t0 is None or _t0 == 0:
            _t0 = min(d['mjd'] for d in _valid_data)

        if _phase_fold and _period and _period > 0:
            _x_values = []
            for d in _valid_data:
                _phase = ((d['mjd'] - _t0) % _period) / _period
                if _phase < 0:
                    _phase += 1
                _x_values.append(_phase)
            _x_title = 'Phase'
            _title_suffix = f' (P = {_period} days)'
            _x_range = [-0.02, 1.02]
        else:
            _x_values = [d['mjd'] for d in _valid_data]
            _x_title = 'MJD (mid-exposure)'
            _title_suffix = ''
            _x_range = None

        _rv_fig = go.Figure()
        _rv_fig.add_trace(go.Scatter(
            x=_x_values,
            y=[d['v_rad'] for d in _valid_data],
            error_y=dict(
                type='data',
                array=[d['e_v_rad'] for d in _valid_data],
                visible=True,
                color='rgba(100, 100, 100, 0.5)'
            ),
            mode='markers',
            marker=dict(
                size=10,
                color=[d['snr'] for d in _valid_data],
                colorscale='Viridis',
                colorbar=dict(title='SNR'),
                line=dict(color='white', width=1)
            ),
            text=[f"Exp {d['index']+1}<br>MJD: {d['mjd']:.4f}<br>v_rad: {d['v_rad']:.2f} km/s<br>SNR: {d['snr']:.1f}"
                  for d in _valid_data],
            hoverinfo='text'
        ))

        _rv_fig.update_layout(
            title=f"Radial Velocity{_title_suffix} - SDSS ID {source_data['sdss_id']}",
            xaxis_title=_x_title,
            yaxis_title='v_rad (km/s)',
            hovermode='closest',
            height=400
        )

        if _x_range:
            _rv_fig.update_xaxes(range=_x_range)

        _rv_fig
    else:
        mo.md("_No radial velocity data available_")
    return


@app.cell
def _(load_error, mo, source_data):
    mo.stop(source_data is None or load_error is not None)
    mo.md("## Spectrum Visualization")
    return


@app.cell
def _(load_error, mo, source_data):
    # Spectrum controls
    mo.stop(source_data is None or load_error is not None)

    _num_exp = source_data['num_exposures']
    _max_checkboxes = min(_num_exp, 30)

    # Create checkboxes for each exposure
    spectrum_checkboxes = mo.ui.array(
        [mo.ui.checkbox(label=f"{i+1}", value=(i == 0)) for i in range(_max_checkboxes)]
    )

    show_all_button = mo.ui.button(label="Show All (up to 20)")
    clear_all_button = mo.ui.button(label="Clear All")

    mo.vstack([
        mo.hstack([show_all_button, clear_all_button], justify="start", gap=1),
        mo.md("**Toggle individual spectra:**"),
        mo.hstack(spectrum_checkboxes, wrap=True, gap=0.5),
        mo.md(f"_Showing checkboxes for {_max_checkboxes} of {_num_exp} exposures_") if _num_exp > 30 else None
    ])
    return clear_all_button, show_all_button, spectrum_checkboxes


@app.cell
def _(
    clear_all_button,
    get_wavelength_grid,
    go,
    load_error,
    mo,
    query_spectrum,
    show_all_button,
    source_data,
    spectrum_checkboxes,
):
    # Spectrum plot
    mo.stop(source_data is None or load_error is not None)

    # Handle show all / clear all buttons
    _show_all_clicked = show_all_button.value
    _clear_all_clicked = clear_all_button.value

    # Determine which spectra to show
    if _clear_all_clicked:
        _visible_indices = []
    elif _show_all_clicked:
        _visible_indices = list(range(min(20, len(spectrum_checkboxes))))
    else:
        _visible_indices = [i for i, cb in enumerate(spectrum_checkboxes) if cb.value]

    if not _visible_indices:
        mo.md("_No spectra selected. Use checkboxes above to show spectra._")
    else:
        _wavelength = get_wavelength_grid()

        _colors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
        ]

        _spec_fig = go.Figure()

        for _idx, _exp_idx in enumerate(_visible_indices):
            try:
                _file_index = source_data['indices'][_exp_idx]
                _spectrum = query_spectrum(_file_index)
                _spec_fig.add_trace(go.Scatter(
                    x=_wavelength,
                    y=_spectrum,
                    mode='lines',
                    name=f'Exp {_exp_idx + 1}',
                    line=dict(width=1, color=_colors[_idx % len(_colors)]),
                    opacity=0.8
                ))
            except Exception as _e:
                print(f"Error loading exposure {_exp_idx + 1}: {_e}")

        _n_shown = len(_visible_indices)
        _spec_fig.update_layout(
            title=f"SDSS ID {source_data['sdss_id']} - {_n_shown} Exposure{'s' if _n_shown > 1 else ''}",
            xaxis_title='Wavelength (Å)',
            yaxis_title='1 + x_starLines_v0',
            hovermode='closest',
            height=500,
            showlegend=True
        )

        _spec_fig
    return


@app.cell
def _(load_error, mo, source_data):
    mo.stop(source_data is None or load_error is not None)
    mo.md("## Exposure Data Table")
    return


@app.cell
def _(load_error, mo, pd, source_data):
    # Exposure table
    mo.stop(source_data is None or load_error is not None)

    _e = source_data['exposures']

    # Get array fields (exposure-level data)
    _array_fields = [k for k, v in _e.items() if isinstance(v, list)]

    # Primary fields to show first
    _primary_fields = [
        'observatory', 'mjd', 'exposure', 'plate_id', 'config_id',
        'fiber_id', 'snr', 'mjd_mid_exposure', 'v_rel', 'v_rad', 'e_v_rad'
    ]

    # Sort fields
    def _sort_key(field):
        if field in _primary_fields:
            return (0, _primary_fields.index(field))
        return (1, field)

    _sorted_fields = sorted(_array_fields, key=_sort_key)

    # Build DataFrame
    _df_data = {'#': list(range(1, source_data['num_exposures'] + 1))}
    for _field in _sorted_fields:
        _df_data[_field] = _e[_field]

    exposure_df = pd.DataFrame(_df_data)

    mo.vstack([
        mo.md(f"**{len(exposure_df)} exposures, {len(_sorted_fields)} columns**"),
        mo.ui.table(exposure_df, selection=None, pagination=True, page_size=20)
    ])
    return (exposure_df,)


@app.cell
def _(mo):
    mo.md(
        """
        ---

        ## Direct Function Usage

        You can also use the data query functions directly in additional cells:

        ```python
        # Load a source
        data = query_source('66899483')

        # Get autocomplete suggestions
        matches = autocomplete('6689', limit=10)

        # Load a spectrum
        spectrum = query_spectrum(data['indices'][0])
        wavelength = get_wavelength_grid()
        ```
        """
    )
    return


if __name__ == "__main__":
    app.run()
