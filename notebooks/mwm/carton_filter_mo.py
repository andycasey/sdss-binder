import marimo

__generated_with = "0.19.7"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # SDSS-V Exposure Filter

    Create subsets of ApogeeReduction.jl-reduced spectra based on SDSS-5 targeting flags (cartons and programs).

    **How to use:**
    1. Select one or more **programs** and/or **cartons** below
    2. Selections within each category are combined with **OR** (any match)
    3. Click **Apply Filter** to see matching exposures
    4. Export the filtered subset to a file
    5. The file will be temporarily saved to the working directory of your BinderHub. You can download it to your computer from there (e.g., by drag-and-drop, or right-click to download).
    """)
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import h5py as h5
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from sdss_semaphore.targeting import TargetingFlags
    # Configuration
    DATA_DIR = Path("data/work/mwm/0.2.0")
    EXPOSURES_FILE = DATA_DIR / "exposures.h5"

    # Load targeting flags and get available cartons/programs
    mo.stop(not EXPOSURES_FILE.exists(), mo.md("❌ **Error:** Exposures file not found. Please check DATA_DIR."))

    with h5.File(EXPOSURES_FILE, "r", locking=False) as fp:
        # Load the targeting flags
        _flags_data = fp["sdss5_target_flags"][:]
        _n_exposures = len(_flags_data)

    # Create TargetingFlags object to get available options
    _tf_sample = TargetingFlags(np.array([[0, 0, 0, 0]], dtype=np.uint64))
    all_programs = sorted([str(p) for p in _tf_sample.all_programs if str(p) not in ('na', 'SKY')])
    all_cartons = sorted([str(c) for c in _tf_sample.all_carton_names])
    mo.md(f"""
    - **Data directory:** `{DATA_DIR}`
    - **Exposures file:** `{EXPOSURES_FILE}`
    - **File exists:** {'✅ Yes' if EXPOSURES_FILE.exists() else '❌ No'}

    **Loaded data:**
    - {_n_exposures:,} total exposures
    - {len(all_programs)} programs available
    - {len(all_cartons)} cartons available
    """)
    return (
        EXPOSURES_FILE,
        TargetingFlags,
        all_cartons,
        all_programs,
        h5,
        mo,
        np,
        pd,
    )


@app.cell(hide_code=True)
def _(mo):
    # Program selection UI
    mo.md("""
    ## Filter Selection

    Select programs and/or cartons to filter exposures. All selections are combined with **OR** logic.
    ### Programs
    """)
    return


@app.cell(hide_code=True)
def _(all_programs, mo):
    # Group programs by prefix for better organization
    _program_groups = {}
    for p in all_programs:
        prefix = p.split('_')[0] if '_' in p else 'other'
        if prefix not in _program_groups:
            _program_groups[prefix] = []
        _program_groups[prefix].append(p)

    program_selector = mo.ui.multiselect(
        options=all_programs,
        label="Select programs (OR logic):",
        full_width=True
    )

    mo.vstack([
        mo.md("_Programs are grouped by survey: `bhm_*` (Black Hole Mapper), `mwm_*` (Milky Way Mapper), `ops_*` (Operations)_"),
        program_selector
    ])
    return (program_selector,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Cartons
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    # Carton selection with search
    carton_search = mo.ui.text(
        placeholder="Type to filter cartons (e.g., 'yso', 'binary', 'rv')...",
        label="Search cartons:",
        full_width=True
    )

    carton_search
    return (carton_search,)


@app.cell(hide_code=True)
def _(all_cartons, carton_search, mo):
    # Filter cartons based on search
    _search_term = carton_search.value.lower().strip() if carton_search.value else ""

    if _search_term:
        filtered_cartons = [c for c in all_cartons if _search_term in c.lower()]
    else:
        filtered_cartons = all_cartons

    carton_selector = mo.ui.multiselect(
        options=filtered_cartons,
        label=f"Select cartons ({len(filtered_cartons)} shown):",
        full_width=True
    )

    mo.vstack([
        mo.md(f"_Showing {len(filtered_cartons)} of {len(all_cartons)} cartons_") if _search_term else None,
        carton_selector
    ])
    return (carton_selector,)


@app.cell(hide_code=True)
def _(carton_selector, mo, program_selector):
    # Selection summary
    _selected_programs = program_selector.value or []
    _selected_cartons = carton_selector.value or []

    _n_programs = len(_selected_programs)
    _n_cartons = len(_selected_cartons)
    _total_selections = _n_programs + _n_cartons

    if _total_selections == 0:
        _summary = "⚠️ **No filters selected** - Select at least one program or carton above."
        _can_filter = False
    else:
        _parts = []
        if _n_programs > 0:
            _parts.append(f"**{_n_programs} program{'s' if _n_programs > 1 else ''}:** {', '.join(_selected_programs)}")
        if _n_cartons > 0:
            _carton_display = ', '.join(_selected_cartons[:5])
            if _n_cartons > 5:
                _carton_display += f" ... and {_n_cartons - 5} more"
            _parts.append(f"**{_n_cartons} carton{'s' if _n_cartons > 1 else ''}:** {_carton_display}")

        _summary = "### Current Selection\n\n" + "\n\n".join(_parts) + "\n\n_All selections combined with OR logic_"
        _can_filter = True

    mo.md(_summary)
    return


@app.cell(hide_code=True)
def _(carton_selector, mo, program_selector):
    # Apply filter button
    _selected_programs = program_selector.value or []
    _selected_cartons = carton_selector.value or []
    _can_filter = len(_selected_programs) > 0 or len(_selected_cartons) > 0

    apply_button = mo.ui.run_button(
        label="🔍 Apply Filter",
        disabled=not _can_filter
    )

    mo.hstack([apply_button], justify="start")
    return (apply_button,)


@app.cell(hide_code=True)
def _(
    EXPOSURES_FILE,
    TargetingFlags,
    apply_button,
    carton_selector,
    h5,
    mo,
    np,
    program_selector,
):
    # Apply the filter
    mo.stop(not apply_button.value)

    _selected_programs = program_selector.value or []
    _selected_cartons = carton_selector.value or []

    mo.stop(
        len(_selected_programs) == 0 and len(_selected_cartons) == 0,
        mo.md("⚠️ No filters selected")
    )

    # Load data and apply filter
    with h5.File(EXPOSURES_FILE, "r", locking=False) as __fp:
        flags_data = __fp["sdss5_target_flags"][:]
        n_total = len(flags_data)

        # Load other useful columns for preview
        preview_columns = {}
        for c in ['sdss_id', 'gaia_ra', 'gaia_dec', 'mjd', 'snr', 'observatory']:
            if c in __fp.keys():
                preview_columns[c] = __fp[c][:]

    # Create TargetingFlags object
    tf = TargetingFlags(flags_data)

    # Build combined mask (OR logic)
    combined_mask = np.zeros(n_total, dtype=bool)

    # Apply program filters
    for prog in _selected_programs:
        prog_mask = tf.in_program(prog)
        combined_mask |= prog_mask

    # Apply carton filters
    for carton in _selected_cartons:
        carton_mask = tf.in_carton_name(carton)
        combined_mask |= carton_mask

    n_matched = np.sum(combined_mask)
    match_indices = np.where(combined_mask)[0]

    mo.md(f"""
    ---
    ## Filter Results

    **Matched:** {n_matched:,} of {n_total:,} exposures ({100*n_matched/n_total:.1f}%)
    """)
    return match_indices, n_matched, preview_columns


@app.cell(hide_code=True)
def _(apply_button, match_indices, mo, n_matched, np, pd, preview_columns):
    # Preview table
    mo.stop(not apply_button.value or n_matched == 0)

    # Build preview DataFrame
    _preview_data = {'index': match_indices[:1000]}  # Limit preview to first 1000

    for _col, _data in preview_columns.items():
        if isinstance(_data[0], bytes):
            _preview_data[_col] = [_data[i].decode('utf-8') if isinstance(_data[i], bytes) else str(_data[i])
                                  for i in match_indices[:1000]]
        elif isinstance(_data[0], np.ndarray):
            # Skip array columns for preview
            continue
        else:
            _preview_data[_col] = [_data[i] for i in match_indices[:1000]]

    preview_df = pd.DataFrame(_preview_data)

    mo.vstack([
        mo.md(f"### Preview (first {min(1000, n_matched):,} rows)"),
        mo.ui.table(preview_df, selection=None, pagination=True, page_size=20)
    ])
    return


@app.cell(hide_code=True)
def _(apply_button, mo, n_matched):
    # Export section
    mo.stop(not apply_button.value or n_matched == 0)

    mo.md("""
    ---
    ## Export Filtered Data
    """)
    return


@app.cell(hide_code=True)
def _(mo, n_matched):
    mo.stop(n_matched == 0)

    output_filename = mo.ui.text(
        value="filtered_exposures",
        label="Output filename (without extension):",
        full_width=False
    )

    output_format = mo.ui.dropdown(
        options={"HDF5 (.h5)": "h5", "CSV (.csv)": "csv", "Parquet (.parquet)": "parquet"},
        value="HDF5 (.h5)",
        label="Output format:"
    )
    mo.hstack([output_filename, output_format], gap=2)
    return output_filename, output_format


@app.cell(hide_code=True)
def _(apply_button, mo, n_matched):
    mo.stop(not apply_button.value or n_matched == 0)

    return


@app.cell(hide_code=True)
def _(mo, n_matched):
    export_button = mo.ui.run_button(
        label="Export selection",
        disabled=n_matched == 0
    )
    export_button
    return (export_button,)


@app.cell(hide_code=True)
def _(
    EXPOSURES_FILE,
    export_button,
    h5,
    match_indices,
    mo,
    n_matched,
    output_filename,
    output_format,
    pd,
):
    mo.stop(not export_button.value or n_matched == 0)

    _filename = output_filename.value or "filtered_exposures"
    _format = output_format.value

    _filename = f"home/{_filename}"

    from tqdm import tqdm
    try:
        # Load all data for export
        with h5.File(EXPOSURES_FILE, "r", locking=False) as _fp:
            _columns_to_export = list(_fp.keys())

            _export_data = {}
            for col in tqdm(_columns_to_export, desc="Collecting"):
                _export_data[col] = _fp[col][match_indices]

        # Export based on format
        if _format == "h5":
            _output_path = f"{_filename}.h5"
            with h5.File(_output_path, "w") as out_fp:
                for col, data in tqdm(_export_data.items(), desc="Writing"):
                    out_fp.create_dataset(col, data=data)

        elif _format == "csv":
            _output_path = f"{_filename}.csv"
            # Convert to DataFrame (flatten arrays if needed)
            _df_data = {}
            for col, data in tqdm(_export_data.items(), desc="Collecting"):
                if len(data.shape) == 1:
                    if isinstance(data[0], bytes):
                        _df_data[col] = [d.decode('utf-8') if isinstance(d, bytes) else str(d) for d in data]
                    else:
                        _df_data[col] = data
            pd.DataFrame(_df_data).to_csv(_output_path, index=False)

        elif _format == "parquet":
            _output_path = f"{_filename}.parquet"
            _df_data = {}
            for col, data in _export_data.items():
                if len(data.shape) == 1:
                    if isinstance(data[0], bytes):
                        _df_data[col] = [d.decode('utf-8') if isinstance(d, bytes) else str(d) for d in data]
                    else:
                        _df_data[col] = data
            pd.DataFrame(_df_data).to_parquet(_output_path, index=False)

        export_result = f"✅ **Success!** Exported {n_matched:,} rows to `{_output_path}`"
    except Exception as e:
        export_result = f"❌ **Error:** {str(e)}"

    mo.md(export_result)
    return


@app.cell(hide_code=True)
def _(apply_button, mo, n_matched):
    # Spectra export section
    mo.stop(not apply_button.value or n_matched == 0)

    mo.md("""
    ---
    ## Export Filtered Spectra

    Export flux and inverse variance (ivar) spectra for the filtered exposures.
    """)
    return


@app.cell(hide_code=True)
def _(mo, n_matched):
    mo.stop(n_matched == 0)

    spectra_filename = mo.ui.text(
        value="filtered_spectra",
        label="Output filename (without extension):",
        full_width=False
    )

    spectra_filename
    return (spectra_filename,)


@app.cell(hide_code=True)
def _(mo, n_matched):
    export_spectra_button = mo.ui.run_button(
        label="Export spectra (HDF5)",
        disabled=n_matched == 0
    )
    export_spectra_button
    return (export_spectra_button,)


@app.cell(hide_code=True)
def _(
    export_spectra_button,
    h5,
    match_indices,
    mo,
    n_matched,
    spectra_filename,
):
    from pathlib import Path as _Path
    mo.stop(not export_spectra_button.value or n_matched == 0)

    _spectra_filename = spectra_filename.value or "filtered_spectra"
    _spectra_filename = f"home/{_spectra_filename}"

    # Define the spectra file paths
    _DATA_DIR = _Path("data/work/mwm/0.2.0")
    _flux_file = _DATA_DIR / "arMADGICS_out_x_starLines_v0.h5"
    _ivar_file = _DATA_DIR / "arMADGICS_out_x_starLines_err_v0.h5"

    from tqdm import tqdm as _tqdm
    try:
        _exported_files = []

        # Export flux file if it exists
        if _flux_file.exists():
            _flux_output = f"{_spectra_filename}_flux.h5"
            with h5.File(_flux_file, "r", locking=False) as _flux_fp:
                _flux_datasets = list(_flux_fp.keys())

                with h5.File(_flux_output, "w") as _out_fp:
                    for _dset in _tqdm(_flux_datasets, desc="Exporting flux"):
                        _data = _flux_fp[_dset][match_indices]
                        _out_fp.create_dataset(_dset, data=_data)

            _exported_files.append(f"`{_flux_output}` (flux)")
        else:
            _exported_files.append(f"⚠️ Flux file not found: `{_flux_file}`")

        # Export ivar file if it exists
        if _ivar_file.exists():
            _ivar_output = f"{_spectra_filename}_ivar.h5"
            with h5.File(_ivar_file, "r", locking=False) as _ivar_fp:
                _ivar_datasets = list(_ivar_fp.keys())

                with h5.File(_ivar_output, "w") as _out_fp:
                    for _dset in _tqdm(_ivar_datasets, desc="Exporting ivar"):
                        _data = _ivar_fp[_dset][match_indices]
                        _out_fp.create_dataset(_dset, data=_data)

            _exported_files.append(f"`{_ivar_output}` (ivar)")
        else:
            _exported_files.append(f"⚠️ Ivar file not found: `{_ivar_file}`")

        _spectra_export_result = f"✅ **Success!** Exported {n_matched:,} spectra to:\n\n" + "\n\n".join([f"- {f}" for f in _exported_files])
    except Exception as _e:
        _spectra_export_result = f"❌ **Error:** {str(_e)}"

    mo.md(_spectra_export_result)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---

    ## Advanced Usage

    You can access the filter results programmatically:

    ```python
    # Get indices of matched rows
    matched_indices = match_indices

    # Get the boolean mask
    mask = combined_mask

    # Access the TargetingFlags object directly
    # Example: count by program
    for prog in tf.all_programs:
        count = np.sum(tf.in_program(prog))
        if count > 0:
            print(f"{prog}: {count}")
    ```
    """)
    return


if __name__ == "__main__":
    app.run()
