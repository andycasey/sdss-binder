import marimo

__generated_with = "0.10.0"
app = marimo.App(width="full", app_title="SDSS-V Sky Map")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        # SDSS-V Sky Coverage Map

        Interactive all-sky visualization of sources with APOGEE spectra from ApogeeReduction.jl.
        
        - **Projection:** Mollweide (Galactic coordinates)
        - **Interaction:** Pan by dragging, zoom with scroll wheel
        - **Tooltips:** Hover over points to see SDSS ID
        """
    )
    return


@app.cell(hide_code=True)
def _():
    import h5py as h5
    import numpy as np
    import os
    from pathlib import Path
    return Path, h5, np, os


@app.cell(hide_code=True)
def _(Path, os):
    # Configuration - modify DATA_DIR if needed
    DATA_DIR = Path(os.environ.get("DATA_DIR", os.path.expanduser("~/../sdssv/ceph/work/scifest/0.2.0")))
    EXPOSURES_FILE = DATA_DIR / "exposures.h5"
    return DATA_DIR, EXPOSURES_FILE


@app.cell(hide_code=True)
def _(EXPOSURES_FILE, h5, mo, np):
    # Load unique sources from exposures file
    mo.stop(not EXPOSURES_FILE.exists(), mo.md(f"❌ **Error:** Exposures file not found at `{EXPOSURES_FILE}`"))
    
    with h5.File(EXPOSURES_FILE, "r", locking=False) as fp:
        # Load the columns we need
        sdss_id_all = fp["sdss_id"][:]
        ra_all = fp["gaia_ra"][:]
        dec_all = fp["gaia_dec"][:]
    
    # Get unique sources (first occurrence of each sdss_id)
    _, unique_indices = np.unique(sdss_id_all, return_index=True)
    
    sdss_id = sdss_id_all[unique_indices]
    ra = ra_all[unique_indices]
    dec = dec_all[unique_indices]
    
    # Filter out invalid coordinates
    valid_mask = np.isfinite(ra) & np.isfinite(dec)
    sdss_id = sdss_id[valid_mask]
    ra = ra[valid_mask]
    dec = dec[valid_mask]
    
    n_sources = len(sdss_id)
    n_total_exposures = len(sdss_id_all)
    
    mo.md(f"""
    **Loaded:** {n_sources:,} unique sources from {n_total_exposures:,} total exposures
    """)
    return dec, n_sources, n_total_exposures, ra, sdss_id


@app.cell(hide_code=True)
def _(mo, np, ra, dec, sdss_id, n_sources):
    # Convert coordinates to Galactic using astropy-style transformation
    # For simplicity, we'll pass ICRS coordinates to Aladin and let it handle the frame
    # But we need to prepare the catalog data for JavaScript
    
    # Limit to reasonable number for browser performance
    MAX_POINTS = 100000
    
    if n_sources > MAX_POINTS:
        # Random subsample for display
        rng = np.random.default_rng(42)
        sample_idx = rng.choice(n_sources, MAX_POINTS, replace=False)
        display_ra = ra[sample_idx]
        display_dec = dec[sample_idx]
        display_sdss_id = sdss_id[sample_idx]
        display_note = f"Showing {MAX_POINTS:,} of {n_sources:,} sources (random sample for performance)"
    else:
        display_ra = ra
        display_dec = dec
        display_sdss_id = sdss_id
        display_note = f"Showing all {n_sources:,} sources"
    
    # Create JavaScript array of sources for Aladin catalog
    # Format: [{ra: float, dec: float, sdss_id: int}, ...]
    sources_json = "[\n"
    for i in range(len(display_ra)):
        sources_json += f'  {{"ra": {display_ra[i]:.6f}, "dec": {display_dec[i]:.6f}, "sdss_id": {display_sdss_id[i]}}}'
        if i < len(display_ra) - 1:
            sources_json += ","
        sources_json += "\n"
    sources_json += "]"
    
    mo.md(f"_{display_note}_")
    return sources_json, display_note


@app.cell(hide_code=True)
def _(mo, sources_json):
    # Create the Aladin Lite viewer embedded in HTML
    # Using Aladin Lite v3 with MOL projection and galactic frame
    
    aladin_html = f"""
    <div id="aladin-lite-div" style="width: 100%; height: 700px;"></div>
    
    <link rel="stylesheet" href="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.min.css" />
    <script type="text/javascript" src="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"></script>
    
    <script type="text/javascript">
    (function() {{
        // Wait for Aladin to be ready
        A.init.then(() => {{
            // Create Aladin Lite instance
            let aladin = A.aladin('#aladin-lite-div', {{
                survey: 'P/2MASS/color',  // 2MASS all-sky as background
                projection: 'MOL',         // Mollweide projection
                cooFrame: 'galactic',      // Galactic coordinate frame
                fov: 360,                  // Full sky view
                showCooGrid: true,         // Show coordinate grid
                showCooGridControl: true,
                showSettingsControl: true,
                showShareControl: false,
                showFullscreenControl: true,
                showLayersControl: true,
                showGotoControl: true,
                showZoomControl: true,
                showFrame: true
            }});
            
            // Source data from Python
            let sourceData = {sources_json};
            
            // Create a catalog for our sources
            let catalog = A.catalog({{
                name: 'SDSS-V Sources',
                sourceSize: 8,
                color: '#ff6b6b',
                shape: 'circle',
                onClick: 'showPopup'
            }});
            
            // Add sources to catalog
            let sources = sourceData.map(s => {{
                return A.source(s.ra, s.dec, {{
                    sdss_id: s.sdss_id,
                    popupTitle: 'SDSS-V Source',
                    popupDesc: '<b>SDSS ID:</b> ' + s.sdss_id + '<br><b>RA:</b> ' + s.ra.toFixed(5) + '°<br><b>Dec:</b> ' + s.dec.toFixed(5) + '°'
                }});
            }});
            
            catalog.addSources(sources);
            aladin.addCatalog(catalog);
            
            // Set initial view to galactic center
            aladin.gotoRaDec(266.417, -29.008);  // Galactic center in ICRS
        }});
    }})();
    </script>
    """
    
    mo.Html(aladin_html)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        ---
        
        ### Controls
        
        | Action | How |
        |--------|-----|
        | **Pan** | Click and drag |
        | **Zoom** | Mouse scroll wheel |
        | **View source info** | Click on a point |
        | **Toggle layers** | Use the layers control (top right) |
        | **Change survey** | Click survey selector (bottom) |
        | **Full screen** | Click fullscreen button |
        
        ### About
        
        This map shows the sky positions of all unique sources with APOGEE spectra reduced by 
        [ApogeeReduction.jl](https://github.com/andycasey/ApogeeReduction.jl). The background 
        shows the 2MASS all-sky survey in false color.
        
        Coordinates are displayed in the **Galactic** frame with a **Mollweide** equal-area projection.
        """
    )
    return


if __name__ == "__main__":
    app.run()
