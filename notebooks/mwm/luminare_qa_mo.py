import marimo

__generated_with = "0.23.3"
app = marimo.App(width="full", app_title="Luminare QA")


@app.cell(hide_code=True)
def title(mo):
    mo.md(r"""
    # Luminare stellar parameter QA

    Every subsection below is **collapsed and un-computed** until you open it.
    Opening one runs its plotting code in the kernel; the contents list in the
    sidebar jumps to any of them.

    Changing the run or the sample cuts invalidates every section, which
    collapses them all again -- nothing is recomputed until you ask for it.
    """)
    return


@app.cell
def imports():
    import marimo as mo
    import numpy as np
    import pandas as pd
    import h5py as h5
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from pathlib import Path
    import re

    return Path, go, h5, make_subplots, mo, np, pd, re


@app.cell
def nav(mo, re):
    # Every title is declared once, here. The sidebar contents, the section headings
    # and the accordion labels are all generated from this dict, so the contents
    # list cannot drift away from the notebook.
    TOC = {
        "1. Sample and fit diagnostics": [
            "1.1 Sample definition",
            "1.2 Kiel diagram",
            "1.3 Chemical plane",
            "1.4 Fit quality",
            "1.5 Edge-of-box flags",
        ],
        "2. ASPCAP DR17 comparison": [
            "2.1 Effective temperature",
            "2.2 Surface gravity",
            "2.3 Metallicity",
            "2.4 Residual trends",
        ],
        "3. Interferometric Teff": [
            "3.1 Angular diameter sample",
        ],
        "4. Asteroseismic logg": [
            "4.1 Seismic vs spectroscopic logg",
            "4.2 Evolutionary state",
        ],
        "5. Star clusters": [
            "5.1 Metallicity scatter",
            "5.2 Abundance trends",
        ],
        "6. Astrometric logg": [
            "6.1 Gaia parallax logg",
        ],
        "7. Internal precision": [
            "7.1 Reported vs empirical uncertainties",
            "7.2 Single star inspection",
        ],
        "8. Data model": [
            "8.1 Fields",
            "8.2 Statistics",
        ],
    }

    SECTION_TITLES = {int(_t.split(".")[0]): _t for _t in TOC}

    def anchor(title):
        """Reproduce marimo's heading-slug rule so links land on the headings."""
        return re.sub(r"\s+", "-", re.sub(r"[^a-z0-9\s]", "", title.lower()).strip())


    _TOC_CSS = """
    <style>
    .qa-toc a { color: inherit; text-decoration: none; display: block; padding: 1px 0; }
    .qa-toc a:hover { text-decoration: underline; }
    .qa-toc .sec { font-weight: 600; margin-top: 0.6rem; }
    .qa-toc .sub { margin-left: 0.85rem; opacity: 0.72; }
    </style>
    """

    _links = []
    for _t, _subs in TOC.items():
        _links.append(f'<div class="sec"><a href="#{anchor(_t)}">{_t}</a></div>')
        _links += [
            f'<div class="sub"><a href="#{anchor(_s)}">{_s}</a></div>' for _s in _subs
        ]

    # Plain <a href="#..."> on purpose. mo.nav_menu rewrites every href to carry the
    # file= query parameter, and that rewrite resolves against the server root, so
    # behind a path prefix (JupyterHub, Binder) its links navigate out of the
    # notebook entirely. A bare fragment is resolved against the current URL.
    mo.sidebar(
        [
            mo.md("### Luminare QA"),
            mo.Html(_TOC_CSS + '<nav class="qa-toc">' + "".join(_links) + "</nav>"),
        ],
        footer=mo.md("<small>Figures are built when a subsection is opened.</small>"),
        width="290px",
    )
    return SECTION_TITLES, TOC, anchor


@app.cell
def run_config(Path, mo):
    LUMIERE_MODELS = Path("/home/jovyan/data/work/acasey/lumiere/models")

    # <date>/<config>/fit/results.h5
    _paths = sorted(LUMIERE_MODELS.glob("*/*/fit/results.h5"))
    _options = {
        str(p.relative_to(LUMIERE_MODELS).parent.parent): str(p) for p in _paths
    }
    _default = "2026-09-03/th_starlines"

    run = mo.ui.dropdown(
        options=_options,
        value=_default if _default in _options else (
            next(iter(_options), None)
        ),
        label="Luminare run",
    )
    run
    return (run,)


@app.cell
def loader(h5, mo, np, pd):
    def file_order(group):
        """Dataset names in the order the writer created them.

        HDF5 only preserves creation order when the file was written asking for it,
        and these files were not; the datasets are chunked, so they have no single
        byte offset either. The object-header address is the remaining honest
        proxy, and it recovers the writing order exactly: identifiers first, then
        the cross-matches, then the fit output. Plain iteration would give
        alphabetical order instead, which scatters related fields.
        """
        _addressed = [
            (h5.h5o.get_info(_v.id).addr, _k)
            for _k, _v in group.items()
            if isinstance(_v, h5.Dataset)
        ]
        return [_k for _, _k in sorted(_addressed)]


    @mo.cache
    def load_catalog(path):
        """Per-star columns, normalised onto the 2026-09-03 flat schema.

        The new files are flat: one dataset per column, lower-case names, boolean
        `flag_*` columns, and a description (usually a unit) on every dataset.
        Older files kept the labels packed in a `theta` matrix, the cross-match
        columns under `meta/`, and the flags bit-packed into an int32; those are
        translated here so that everything downstream sees a single schema.

        Deliberately reads only the 1-D datasets, so a file that still carries
        `model_flux` and `continuum` costs the same as one that does not.
        """
        _d = {}
        with h5.File(path, "r") as fp:
            # attrs still spell the labels the old way; the columns are lower case
            _labels = [_b.decode().lower() for _b in fp.attrs["label_names"]]

            def _take(_group):
                for _k in file_order(_group):
                    _v = _group[_k]
                    if _v.ndim == 1:
                        _d[_k] = np.char.decode(_v[:]) if _v.dtype.kind == "S" else _v[:]

            _take(fp)
            if "meta" in fp:  # pre-2026-09-03 layout
                _take(fp["meta"])

            for _suffix, _key in (
                ("", "theta"), ("_err", "theta_err"),
                ("_init", "theta_init"), ("_s1", "theta_s1"),
            ):
                if _key in fp:
                    _matrix = fp[_key][:]
                    for _i, _name in enumerate(_labels):
                        _d[_name + _suffix] = _matrix[:, _i]

            # unpack the old int32 flags into the boolean columns the new files have
            if "flags" in _d and not any(_k.startswith("flag_") for _k in _d):
                for _i, _bit in enumerate(fp.attrs["flag_bits"]):
                    _name = _bit.decode().split(": ", 1)[1].split(" ")[0].lower()
                    _d["flag_" + _name] = ((_d["flags"] >> _i) & 1).astype(bool)

        _df = pd.DataFrame(_d)
        _df["rchi2"] = _df["chi2"] / _df["dof"]
        if "log10_vmic" in _df:
            _df["vmic"] = 10 ** _df["log10_vmic"]
        _flags = sorted(_k for _k in _df.columns if _k.startswith("flag_"))
        return _df, _labels, _flags

    return file_order, load_catalog


@app.cell
def schema_loader(Path, file_order, h5, load_catalog, mo, np, pd):
    # Paths inside the files are recorded as they were on the cluster.
    CLUSTER_ROOT = "/mnt/home/acasey/sdssv-work/"
    LOCAL_ROOT = "/home/jovyan/data/work/"


    @mo.cache
    def load_schema(path):
        """Field-level metadata for a run, read from the HDF5 attributes.

        Units and descriptions come from each dataset's own attributes. Only a
        couple of fields carry them today, so most rows come through blank and will
        fill in on their own as the writers start emitting them. Anything the
        results file does not describe is looked up in the coadd it was fitted
        from, whose path is recorded in the file's `coadd` attribute.
        """

        def _walk(fp):
            """Every dataset, in the order the writer created them."""
            _out = {}
            for _group_name in [""] + [
                _k for _k, _v in fp.items() if isinstance(_v, h5.Group)
            ]:
                _group = fp[_group_name] if _group_name else fp
                for _k in file_order(_group):
                    _obj = _group[_k]
                    _full = f"{_group_name}/{_k}" if _group_name else _k
                    _out[_full] = (_obj.shape, _obj.dtype, dict(_obj.attrs))
            return _out

        with h5.File(path, "r") as fp:
            _fields = _walk(fp)
            _coadd = fp.attrs.get("coadd", "")

        _coadd = str(_coadd).replace(CLUSTER_ROOT, LOCAL_ROOT)
        _from_coadd = {}
        if _coadd and Path(_coadd).exists():
            with h5.File(_coadd, "r") as fp:
                _from_coadd = _walk(fp)

        def _text(_attrs, *_keys):
            for _k in _keys:
                _v = _attrs.get(_k)
                if isinstance(_v, bytes):
                    _v = _v.decode()
                if _v not in (None, ""):
                    return str(_v)
            return ""

        _rows = []
        for _name, (_shape, _dtype, _attrs) in _fields.items():
            _fallback = _from_coadd.get(_name, ((), None, {}))[2]
            _desc = _text(_attrs, "description") or _text(_fallback, "description")
            _by = _text(_attrs, "derived_by") or _text(_fallback, "derived_by")
            if _by:
                _desc = (_desc + " " if _desc else "") + f"(derived by `{_by}`)"
            _rows.append(
                {
                    "field": _name,
                    "description": _desc,
                    "unit": _text(_attrs, "unit", "units")
                    or _text(_fallback, "unit", "units"),
                    "dtype": str(_dtype),
                    "shape": " x ".join(f"{_s:,}" for _s in _shape),
                }
            )

        # left in file order: sorting by name scatters related fields
        return pd.DataFrame(_rows), _coadd


    @mo.cache
    def field_statistics(path):
        """Percentiles and the null rate for every numeric column of a run.

        0% and 100% are the minimum and maximum; taking them from the same
        `np.percentile` call as the rest keeps one code path. The null column
        counts everything that is not a finite number -- NaN, +/-inf and the nulls
        a missing cross-match leaves behind -- so it reads as "how much of this
        field is missing", not "how much is usable".
        """
        _cat = load_catalog(path)[0]
        _rows = []
        for _c in _cat.columns:
            _v = _cat[_c].to_numpy()
            if _v.dtype == bool:
                _v = _v.astype(float)
            elif not np.issubdtype(_v.dtype, np.number):
                continue
            _v = _v.astype(float, copy=False)
            _ok = np.isfinite(_v)
            _n = int(_ok.sum())
            _q = (
                np.percentile(_v[_ok], [0, 1, 25, 50, 75, 99, 100])
                if _n
                else [np.nan] * 7
            )
            _rows.append(
                {
                    "field": _c,
                    "% null or non-finite": 100.0 * (_v.size - _n) / max(_v.size, 1),
                    "0%": _q[0],
                    "1%": _q[1],
                    "25%": _q[2],
                    "50%": _q[3],
                    "75%": _q[4],
                    "99%": _q[5],
                    "100%": _q[6],
                }
            )
        return pd.DataFrame(_rows)


    return field_statistics, load_schema


@app.cell
def catalog(load_catalog, mo, run):
    mo.stop(run.value is None, mo.md("**No Luminare `fit/results.h5` files found.**"))

    cat, LABEL_NAMES, FLAG_COLUMNS = load_catalog(run.value)
    mo.md(
        f"`{run.value}` &mdash; **{len(cat):,}** stars, **{len(cat.columns)}** columns, "
        f"labels: `{'`, `'.join(LABEL_NAMES)}`"
    )
    return FLAG_COLUMNS, cat


@app.cell
def controls(FLAG_COLUMNS, mo):
    snr_min = mo.ui.number(0, 10000, value=30, step=5, label="min S/N")
    rchi2_max = mo.ui.number(0.0, 1e4, value=10.0, step=0.5, label="max chi2/dof")
    only_converged = mo.ui.checkbox(True, label="converged only")

    # Reject the fits that ran into the edge of a label box, and the ones that
    # stalled. flag_vsini_edge is left on because it fires for about half the
    # catalogue: it only means vsini sat on the 2 km/s floor, which is what a slow
    # rotator should do. flag_pms and friends are classifications, not faults, so
    # they are never rejected by default either.
    _default_reject = [
        _f for _f in FLAG_COLUMNS
        if (_f.endswith("_edge") and _f != "flag_vsini_edge") or _f == "flag_stalled"
    ]

    reject_flags = mo.ui.multiselect(
        options=FLAG_COLUMNS, value=_default_reject, label="reject flags"
    )

    mo.hstack(
        [snr_min, rchi2_max, only_converged, reject_flags],
        justify="start",
        align="center",
        wrap=True,
        gap=1.5,
    )
    return only_converged, rchi2_max, reject_flags, snr_min


@app.cell
def sample(cat, mo, only_converged, rchi2_max, reject_flags, snr_min):
    _m = (cat["snr"] >= snr_min.value) & (cat["rchi2"] <= rchi2_max.value)
    for _f in reject_flags.value:
        _m &= ~cat[_f].to_numpy().astype(bool)
    if only_converged.value:
        _m &= cat["converged"] == 1

    df = cat[_m.to_numpy()].reset_index(drop=True)

    mo.md(
        f"**{len(df):,}** of {len(cat):,} stars pass the cuts "
        f"({100 * len(df) / max(len(cat), 1):.1f}%)."
    )
    return (df,)


@app.cell
def helpers(SECTION_TITLES, TOC, anchor, go, mo, np):
    # Box and lasso select are gone from the toolbar: selection is parked, so the
    # buttons would draw a marquee that nothing reads. `layout.modebar.remove`
    # rides along in the figure JSON, unlike the `config` route, which a plain
    # figure has no way to pass.
    MODEBAR_REMOVE = ("select2d", "lasso2d")


    def section(number, bodies, *, mode="click"):
        """Render one numbered section from the TOC.

        Each subsection body is a zero-argument callable wrapped in `mo.lazy`, so
        none of them run until the reader asks for them.

        mode="click" (default) puts each subsection in a collapsed accordion:
        nothing runs until it is opened. mode="scroll" renders a subsection when it
        scrolls into view, reserving height for the un-rendered ones so they do not
        all fire on the first paint.
        """
        _title = SECTION_TITLES[number]
        _subs = TOC[_title]
        if len(bodies) != len(_subs):
            raise ValueError(
                f"section {number}: TOC lists {len(_subs)} subsections, "
                f"but {len(bodies)} bodies were given"
            )

        _out = [mo.md(f"## {_title}")]
        for _s, _b in zip(_subs, bodies):
            _lazy = mo.lazy(_b, show_loading_indicator=True)
            # The link target has to sit in the cell's own (light) DOM. marimo
            # renders every plugin -- the accordion included -- inside a shadow
            # root, and browsers do not scroll to fragments inside one.
            _target = mo.Html(
                f'<span id="{anchor(_s)}" '
                f'style="display:block;height:0;scroll-margin-top:3.5rem"></span>'
            )
            if mode == "scroll":
                _item = mo.vstack(
                    [mo.md(f"### {_s}"), _lazy.style({"min-height": "440px"})]
                )
            else:
                _item = mo.accordion({f"### {_s}": _lazy})
            _out += [_target, _item]
        return mo.vstack(_out, gap=0.4)


    def heading(title, level=2):
        """Anchor + heading for content that is not built out of lazy accordions."""
        return mo.vstack(
            [
                mo.Html(
                    f'<span id="{anchor(title)}" '
                    f'style="display:block;height:0;scroll-margin-top:3.5rem"></span>'
                ),
                mo.md(f"{'#' * level} {title}"),
            ],
            gap=0,
        )


    def robust(d):
        """Median offset and MAD-based scatter of a residual array."""
        d = np.asarray(d, float)
        d = d[np.isfinite(d)]
        if d.size == 0:
            return np.nan, np.nan
        _med = np.median(d)
        return _med, 1.4826 * np.median(np.abs(d - _med))


    def density(x, y, *, bins=260, extent=None, ids=None, sparse_max=2,
                colorscale="Magma_r", return_rows=False):
        """Traces showing half a million stars without shipping half a million points.

        The dense locus is binned in the kernel and sent as a log-density heatmap;
        every star landing in a near-empty bin is sent individually instead, so the
        outliers -- which is where QA problems live -- stay hoverable. Arrays leave
        as float32, which plotly serialises as base64 rather than decimal text.

        With `return_rows`, also returns the positions in `x`/`y` of the points in
        the scatter trace, in trace order. That is what makes a box or lasso
        selection resolvable: the frontend strips `customdata` from selected
        points, so the only way back to a star is trace index plus point index.
        """
        x, y = np.asarray(x, float), np.asarray(y, float)
        _ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[_ok], y[_ok]
        _ids = None if ids is None else np.asarray(ids)[_ok]

        _H, _xe, _ye = np.histogram2d(x, y, bins=bins, range=extent)
        _nx, _ny = _H.shape
        _ix = np.clip(np.searchsorted(_xe, x, "right") - 1, 0, _nx - 1)
        _iy = np.clip(np.searchsorted(_ye, y, "right") - 1, 0, _ny - 1)
        _in = (x >= _xe[0]) & (x <= _xe[-1]) & (y >= _ye[0]) & (y <= _ye[-1])
        _sparse = _in & (_H[_ix, _iy] <= sparse_max)

        _hover = "%{x:.5g}, %{y:.5g}<extra></extra>"
        if _ids is not None:
            _hover = "sdss_id %{customdata}<br>%{x:.5g}, %{y:.5g}<extra></extra>"

        _traces = [
            go.Heatmap(
                x=(0.5 * (_xe[:-1] + _xe[1:])).astype("float32"),
                y=(0.5 * (_ye[:-1] + _ye[1:])).astype("float32"),
                z=np.log10(np.where(_H > sparse_max, _H, np.nan)).T.astype("float32"),
                colorscale=colorscale,
                showscale=False,
                hoverinfo="skip",
            ),
            go.Scattergl(
                x=x[_sparse].astype("float32"),
                y=y[_sparse].astype("float32"),
                customdata=None if _ids is None else _ids[_sparse],
                mode="markers",
                marker=dict(size=3, color="#333333", opacity=0.75),
                hovertemplate=_hover,
                showlegend=False,
                name="",
            ),
        ]
        if return_rows:
            # positions in the original arrays, in the same order as the trace
            return _traces, np.flatnonzero(_ok)[_sparse]
        return _traces


    def histogram(values, *, bins=120, extent=None, name="", color="#4c78a8"):
        """Bin in the kernel and ship the bars, not the values."""
        _v = np.asarray(values, float)
        _v = _v[np.isfinite(_v)]
        _n, _e = np.histogram(_v, bins=bins, range=extent)
        return go.Bar(
            x=(0.5 * (_e[:-1] + _e[1:])).astype("float32"),
            y=_n.astype("float32"),
            name=name,
            marker=dict(color=color, line=dict(width=0)),
            hovertemplate="%{x:.4g}: %{y:,.0f} stars<extra></extra>",
        )


    def figure(traces, *, height=450, xlabel="", ylabel="", xrange=None, yrange=None,
               title=None, legend=False):
        """Shared layout: pan by default, and zoom preserved across re-renders."""
        fig = go.Figure(traces)
        fig.update_layout(
            height=height,
            template="plotly_white",
            dragmode="pan",
            hovermode="closest",
            uirevision="keep",
            showlegend=legend or any(getattr(_t, "showlegend", False) for _t in traces),
            margin=dict(l=65, r=20, t=42 if title else 14, b=48),
            title=dict(text=title, font=dict(size=13), x=0.01) if title else None,
            legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="right", x=1.0),
            modebar=dict(remove=MODEBAR_REMOVE),
        )
        fig.update_xaxes(title_text=xlabel, range=xrange)
        fig.update_yaxes(title_text=ylabel, range=yrange)
        return fig


    def view_ranges(element):
        """What the reader is currently looking at, so the kernel can re-bin to it.

        Zoom and pan push {"xaxis": {"range": [...]}, ...} into a `mo.ui.plotly`
        element's value, and marimo keeps that payload verbatim. Returns the ranges,
        the string "auto" when the reader hit autoscale, or None when the view
        carries no range information.
        """
        _d = getattr(element, "_selection_data", None) or {}
        _x, _y = _d.get("xaxis") or {}, _d.get("yaxis") or {}
        if _x.get("autorange") or _y.get("autorange"):
            return "auto"
        if not _x.get("range") or not _y.get("range"):
            return None
        return (
            [float(_v) for _v in _x["range"]],
            [float(_v) for _v in _y["range"]],
        )


    def markdown_table(frame, *, formats=None, na=""):
        """Render a DataFrame as a markdown table.

        Written out by hand because `tabulate` is not installed here, and because
        these tables are read, not sorted -- they live inside a lazy accordion, and
        a `mo.ui.table` there would render but never be interactive.
        """
        _formats = formats or {}

        def _cell(_col, _v):
            if isinstance(_v, (float, np.floating)):
                if not np.isfinite(_v):
                    return na
                return _formats.get(_col, "{:.4g}").format(_v)
            _s = "" if _v is None else str(_v)
            return _s.replace("|", "\\|")

        _cols = list(frame.columns)
        _lines = [
            "| " + " | ".join(_cols) + " |",
            # ":---" rather than "---": left alignment, rather than whatever the
            # renderer would pick for a column that looks numeric
            "|" + "|".join(":---" for _ in _cols) + "|",
        ]
        for _row in frame.itertuples(index=False):
            _lines.append(
                "| " + " | ".join(_cell(_c, _v) for _c, _v in zip(_cols, _row)) + " |"
            )
        return mo.md("\n".join(_lines))


    return (
        MODEBAR_REMOVE,
        density,
        figure,
        histogram,
        markdown_table,
        robust,
        section,
        view_ranges,
    )


@app.cell
def comparison(MODEBAR_REMOVE, density, go, make_subplots, np, robust):
    def one_to_one(x, y, *, xerr=None, yerr=None, ids=None, xlabel="reference",
                   ylabel="Luminare", lim=None, rlim=None, bins=220, height=640,
                   overlay=None):
        """Reference (x) against Luminare (y), with a residual panel and robust stats.

        Above 5000 stars this draws the density map; below it, individual points
        with error bars. The same call therefore serves 30 interferometric targets
        and 400k survey stars.
        """
        x, y = np.asarray(x, float), np.asarray(y, float)
        _ok = np.isfinite(x) & np.isfinite(y)
        _xe = np.asarray(xerr, float)[_ok] if xerr is not None else None
        _ye = np.asarray(yerr, float)[_ok] if yerr is not None else None
        _ids = None if ids is None else np.asarray(ids)[_ok]
        x, y = x[_ok], y[_ok]
        _d = y - x
        _bias, _sig = robust(_d)

        if lim is None:
            if x.size:
                _lo, _hi = np.percentile(np.concatenate([x, y]), [0.5, 99.5])
            else:
                _lo, _hi = 0.0, 1.0
            _pad = 0.05 * (_hi - _lo) or 1.0
            lim = (_lo - _pad, _hi + _pad)
        if rlim is None:
            rlim = (_bias - 5 * _sig, _bias + 5 * _sig) if x.size else (-1.0, 1.0)

        fig = make_subplots(
            rows=2, cols=1, shared_xaxes=True, row_heights=[0.66, 0.34],
            vertical_spacing=0.03,
        )
        if x.size > 5000:
            for _t in density(x, y, bins=bins, extent=[list(lim), list(lim)], ids=_ids):
                fig.add_trace(_t, row=1, col=1)
            for _t in density(x, _d, bins=bins, extent=[list(lim), list(rlim)], ids=_ids):
                fig.add_trace(_t, row=2, col=1)
        else:
            _hover = "%{x:.5g}, %{y:.5g}<extra></extra>"
            if _ids is not None:
                _hover = "sdss_id %{customdata}<br>%{x:.5g}, %{y:.5g}<extra></extra>"
            for _row, _yv in ((1, y), (2, _d)):
                fig.add_trace(
                    go.Scatter(
                        x=x, y=_yv, customdata=_ids, mode="markers",
                        marker=dict(size=7, color="#4c78a8",
                                    line=dict(width=0.5, color="white")),
                        error_x=dict(array=_xe, thickness=0.8, width=0) if _xe is not None else None,
                        error_y=dict(array=_ye, thickness=0.8, width=0) if _ye is not None else None,
                        hovertemplate=_hover, showlegend=False, name="",
                    ),
                    row=_row, col=1,
                )

        fig.add_trace(
            go.Scatter(x=list(lim), y=list(lim), mode="lines",
                       line=dict(color="crimson", width=1), hoverinfo="skip",
                       showlegend=False),
            row=1, col=1,
        )
        for _t in overlay or []:
            fig.add_trace(_t, row=1, col=1)

        fig.add_hline(y=0, row=2, col=1, line=dict(color="crimson", width=1))
        fig.add_hline(y=_bias, row=2, col=1,
                      line=dict(color="#1f77b4", width=1, dash="dash"))

        fig.update_layout(
            height=height, template="plotly_white", dragmode="pan",
            hovermode="closest", uirevision="keep", showlegend=bool(overlay),
            modebar=dict(remove=MODEBAR_REMOVE),
            margin=dict(l=70, r=20, t=42, b=48),
            title=dict(
                text=f"N = {x.size:,}   offset = {_bias:+.4g}   scatter = {_sig:.4g}",
                font=dict(size=13), x=0.01,
            ),
        )
        fig.update_xaxes(range=list(lim), row=1, col=1)
        fig.update_xaxes(range=list(lim), title_text=xlabel, row=2, col=1)
        fig.update_yaxes(range=list(lim), title_text=ylabel, row=1, col=1)
        fig.update_yaxes(range=list(rlim), title_text="residual", row=2, col=1)
        return fig

    return (one_to_one,)


@app.cell
def kiel_state(mo):
    # The Kiel diagram re-bins to whatever is on screen, so it needs somewhere to
    # remember the view. A plain `mo.ui.plotly` cannot live inside a `mo.lazy` body
    # and stay interactive, so the element is built here, at cell level, and section
    # 1.2 renders it.
    KIEL_LIMITS = ([8000.0, 2800.0], [5.6, -0.2])

    get_kiel_view, set_kiel_view = mo.state(None)
    return KIEL_LIMITS, get_kiel_view, set_kiel_view


@app.cell
def kiel_chart(KIEL_LIMITS, density, df, figure, get_kiel_view, mo):
    _view = get_kiel_view() or KIEL_LIMITS
    _xl, _yl = _view

    # 200 bins across the visible range: zooming in buys real resolution instead of
    # magnifying the same coarse cells, and the smaller bins turn more stars into
    # individually drawn, hoverable points.
    _kiel_fig = figure(
        density(df["teff"], df["logg"], bins=200, extent=[sorted(_xl), sorted(_yl)],
                ids=df["sdss_id"].to_numpy()),
        height=560, xlabel="Teff [K]", ylabel="log g", xrange=_xl, yrange=_yl,
    )
    # drag to zoom rather than pan, and a fresh uirevision per view so the explicit
    # ranges are applied instead of plotly preserving the previous ones
    _kiel_fig.update_layout(dragmode="zoom", uirevision=repr(_view))

    kiel_plot = mo.ui.plotly(_kiel_fig)
    None
    return (kiel_plot,)


@app.cell
def kiel_zoom(
    KIEL_LIMITS,
    get_kiel_view,
    kiel_plot,
    set_kiel_view,
    view_ranges,
):
    # Zoom, pan or autoscale the Kiel diagram and the kernel re-bins to the new
    # limits. The guard matters: without it the rebuilt figure's own relayout would
    # feed straight back in here and the two cells would ping-pong.
    _v = view_ranges(kiel_plot)
    if _v == "auto":
        set_kiel_view(None)
    elif _v is not None:
        _old = get_kiel_view() or KIEL_LIMITS
        _moved = any(
            abs(_new - _was) > 0.02 * (abs(_span[1] - _span[0]) or 1.0)
            for _new, _was, _span in zip(
                _v[0] + _v[1], _old[0] + _old[1], [_old[0]] * 2 + [_old[1]] * 2
            )
        )
        if _moved:
            set_kiel_view(_v)
    return


@app.cell(hide_code=True)
def sec_overview(
    FLAG_COLUMNS,
    cat,
    density,
    df,
    figure,
    go,
    histogram,
    kiel_plot,
    mo,
    np,
    only_converged,
    rchi2_max,
    reject_flags,
    run,
    section,
    snr_min,
):
    def _sample_definition():
        _rows = [
            ("run", f"`{run.value}`"),
            ("S/N", f"&ge; {snr_min.value}"),
            ("chi2/dof", f"&le; {rchi2_max.value}"),
            ("converged", "required" if only_converged.value else "not required"),
            ("rejected flags", ", ".join(reject_flags.value) or "none"),
            ("stars passing", f"**{len(df):,}** of {len(cat):,} ({100 * len(df) / len(cat):.1f}%)"),
        ]
        _table = "\n".join(f"| {_k} | {_v} |" for _k, _v in _rows)
        return mo.md(f"| | |\n|---|---|\n{_table}")


    def _chemical_plane():
        return figure(
            density(df["m_h"], df["alpha_m"], bins=260, extent=[[-2.6, 1.1], [-1.1, 1.1]],
                    ids=df["sdss_id"].to_numpy()),
            height=520, xlabel="[M/H]", ylabel="[alpha/M]",
        )


    def _fit_quality():
        _snr = np.log10(np.clip(df["snr"].to_numpy(), 1e-3, None))
        _rc = np.log10(np.clip(df["rchi2"].to_numpy(), 1e-3, None))
        _scatter = figure(
            density(_snr, _rc, bins=240, ids=df["sdss_id"].to_numpy()),
            height=430, xlabel="log10 S/N", ylabel="log10 chi2/dof",
        )
        _hist = figure(
            [
                histogram(cat["rchi2"], bins=150, extent=(0, 20), name="all", color="#c8c8c8"),
                histogram(df["rchi2"], bins=150, extent=(0, 20), name="sample"),
            ],
            height=320, xlabel="chi2/dof", ylabel="stars", legend=True,
        )
        _hist.update_layout(barmode="overlay")
        return mo.vstack([_scatter, _hist])


    def _flag_rates():
        _frac = [float(cat[_f].to_numpy().astype(bool).mean()) for _f in FLAG_COLUMNS]
        return figure(
            [go.Bar(x=[_f.replace("flag_", "") for _f in FLAG_COLUMNS], y=_frac,
                    marker=dict(color="#e45756"),
                    hovertemplate="%{x}: %{y:.1%} of the catalogue<extra></extra>")],
            height=420, ylabel="fraction of all stars",
        )


    section(
        1,
        [_sample_definition, kiel_plot, _chemical_plane, _fit_quality, _flag_rates],
    )
    return


@app.cell(hide_code=True)
def sec_dr17(
    MODEBAR_REMOVE,
    density,
    df,
    make_subplots,
    mo,
    one_to_one,
    robust,
    section,
):
    # The 2026-09-03 files dropped the ASPCAP columns in favour of the Gaia and
    # 2MASS cross-matches, so this section reports that rather than failing.
    DR17_PAIRS = [
        ("dr17_teff", "teff", "ASPCAP Teff [K]", "Luminare Teff [K]"),
        ("dr17_logg", "logg", "ASPCAP log g", "Luminare log g"),
        ("dr17_x_h", "m_h", "ASPCAP [M/H]", "Luminare [M/H]"),
    ]
    HAVE_DR17 = all(_ref in df.columns for _ref, _, _, _ in DR17_PAIRS)


    def _missing():
        return mo.md(
            "*This run does not carry the ASPCAP DR17 columns (`dr17_teff`, "
            "`dr17_logg`, `dr17_x_h`). The 2026-09-03 files replaced them with the "
            "Gaia DR3 and 2MASS cross-matches: `gaia_teff_gspphot`, "
            "`gaia_logg_gspphot`, `gaia_mh_gspphot` and the 2MASS photometry are "
            "all present, if you want this comparison repointed at those.*"
        )


    def _dr17(ref, new, xlabel, ylabel):
        return one_to_one(
            df[ref], df[new], yerr=df.get(new + "_err"), ids=df["sdss_id"].to_numpy(),
            xlabel=xlabel, ylabel=ylabel,
        )


    def _teff():
        return _dr17(*DR17_PAIRS[0]) if HAVE_DR17 else _missing()


    def _logg():
        return _dr17(*DR17_PAIRS[1]) if HAVE_DR17 else _missing()


    def _metallicity():
        return _dr17(*DR17_PAIRS[2]) if HAVE_DR17 else _missing()


    def _residual_trends():
        if not HAVE_DR17:
            return _missing()
        _fig = make_subplots(
            rows=1, cols=3, horizontal_spacing=0.07,
            subplot_titles=[f"delta {_n}" for _, _n, _, _ in DR17_PAIRS],
        )
        for _i, (_ref, _new, _, _) in enumerate(DR17_PAIRS, start=1):
            _d = (df[_new] - df[_ref]).to_numpy()
            _b, _s = robust(_d)
            for _t in density(df["teff"], _d, bins=180,
                              extent=[[2800, 8000], [_b - 5 * _s, _b + 5 * _s]],
                              ids=df["sdss_id"].to_numpy()):
                _fig.add_trace(_t, row=1, col=_i)
            _fig.add_hline(y=0, row=1, col=_i, line=dict(color="crimson", width=1))
            _fig.update_xaxes(range=[8000, 2800], title_text="Luminare Teff [K]",
                              row=1, col=_i)
        _fig.update_layout(
            height=430, template="plotly_white", dragmode="pan", hovermode="closest",
            uirevision="keep", modebar=dict(remove=MODEBAR_REMOVE),
            margin=dict(l=60, r=20, t=40, b=48),
        )
        return _fig


    section(2, [_teff, _logg, _metallicity, _residual_trends])
    return


@app.cell(hide_code=True)
def sec_interferometry(mo, section):
    def _angular_diameters():
        return mo.md(
            """
            **TODO.** Angular-diameter temperatures for the usual samples
            (Karovicova et al. benchmarks, Boyajian et al. dwarfs).

            Join a reference table on `sdss_id` or `gaiaedr3_source_id`, then:

            ```python
            one_to_one(ref["teff"], m["Teff"], xerr=ref["e_teff"], yerr=m["e_Teff"],
                       ids=m["sdss_id"], xlabel="interferometric Teff [K]",
                       ylabel="Luminare Teff [K]")
            ```

            Below 5000 stars `one_to_one` draws points with error bars, so the small
            sample needs no special handling.
            """
        )


    section(3, [_angular_diameters])
    return


@app.cell(hide_code=True)
def sec_seismology(mo, section):
    def _seismic_logg():
        return mo.md(
            """
            **TODO.** Asteroseismic log g (APOKASC-3, K2/GAP, TESS) against the
            spectroscopic log g, as a function of Teff and [M/H].
            """
        )


    def _evolutionary_state():
        return mo.md(
            """
            **TODO.** Split by evolutionary state (RGB / RC / secondary clump) and
            check whether the log g offset differs between them.
            """
        )


    section(4, [_seismic_logg, _evolutionary_state])
    return


@app.cell(hide_code=True)
def sec_clusters(mo, section):
    def _cluster_metallicity():
        return mo.md(
            """
            **TODO.** Open and globular cluster members: the star-to-star [M/H]
            scatter within a cluster is a direct measure of the internal precision.
            """
        )


    def _cluster_trends():
        return mo.md(
            """
            **TODO.** Abundance trends along the giant branch within each cluster,
            which should be flat.
            """
        )


    section(5, [_cluster_metallicity, _cluster_trends])
    return


@app.cell(hide_code=True)
def sec_astrometric(mo, section):
    def _astrometric_logg():
        return mo.md(
            """
            **TODO.** Astrometric log g from the Gaia parallax and a bolometric
            correction, against the spectroscopic log g. `gaiaedr3_parallax` and the
            Bailer-Jones `gaiaedr3_r_med_geo` distances are already columns of `df`.
            """
        )


    section(6, [_astrometric_logg])
    return


@app.cell(hide_code=True)
def sec_repeats(mo, section):
    def _uncertainties():
        return mo.md(
            """
            **TODO.** Reported `theta_err` against the observed visit-to-visit
            scatter for repeat observations, per label and as a function of S/N.
            """
        )


    def _single_star():
        return mo.md(
            """
            **TODO.** Pick a star and show its spectrum and model.

            This one wants a *reactive* selection rather than a lazy body: a UI
            element created inside `mo.lazy` cannot drive other cells. Gate it on a
            `mo.ui.switch` plus `mo.stop` in an ordinary cell instead, which keeps
            the on-demand behaviour without giving up reactivity.
            """
        )


    section(7, [_uncertainties, _single_star])
    return


@app.cell
def sec_schema(
    cat,
    field_statistics,
    load_schema,
    markdown_table,
    mo,
    run,
    section,
):
    def _fields():
        _schema, _coadd = load_schema(run.value)
        _described = int((_schema["description"] != "").sum())
        _note = mo.md(
            f"**{len(_schema)}** datasets in `{run.value}`. "
            f"Units and descriptions are read straight from the HDF5 attributes; "
            f"**{_described}** of them carry one today, and the rest fill in as the "
            f"writers start emitting them. Missing entries fall back to the coadd "
            f"this run was fitted from"
            + (f": `{_coadd}`." if _coadd else ", which is not readable from here.")
        )
        return mo.vstack([_note, markdown_table(_schema)], gap=0.5)


    def _statistics():
        _stats = field_statistics(run.value)
        _note = mo.md(
            f"Every numeric column of the catalogue as loaded, all "
            f"{len(cat):,} rows, before the sample cuts at the top."
        )
        return mo.vstack(
            [
                _note,
                markdown_table(_stats, formats={"% null or non-finite": "{:.2f}"}),
            ],
            gap=0.5,
        )


    section(8, [_fields, _statistics])
    return


@app.cell
def parked_selection(mo):
    # Lasso/box selection and named subsets, parked.
    #
    # Everything below drove the shared selector: pick a plot, box or lasso the
    # individually drawn points, name the group, and see it overlaid on every
    # other figure. Removed from the running notebook, kept here verbatim so it
    # can be brought back by splitting these blocks back into their own cells.
    #
    # Two things it depended on, if you restore it: a UI element is only
    # interactive when a cell-level binding can reach it (hence mo.ui.array for
    # the per-plot buttons), and nothing inside a mo.lazy body qualifies -- which
    # is why there was one shared selector rather than one per figure.
    #
    # --- axis and selectable-plot configuration (lived in the `nav` cell) ---
    #
    # # Axis labels, and which of them read better reversed.
    # AXES = {
    #     "Teff": "Teff [K]",
    #     "logg": "log g",
    #     "m_H": "[M/H]",
    #     "alpha_m": "[alpha/M]",
    #     "C_m": "[C/M]",
    #     "N_m": "[N/M]",
    #     "vmic": "vmic [km/s]",
    #     "vsini": "vsini [km/s]",
    #     "snr": "S/N",
    #     "rchi2": "chi2/dof",
    #     "dr17_teff": "ASPCAP Teff [K]",
    #     "dr17_logg": "ASPCAP log g",
    #     "dr17_x_h": "ASPCAP [M/H]",
    # }
    # REVERSED_AXES = ("Teff", "logg", "dr17_teff", "dr17_logg")
    #
    # # Which subsections can hand their axes to the shared selector. A subsection
    # # whose y is not a plain column (a residual, a histogram) is simply absent.
    # SELECTABLE = {
    #     "1.2 Kiel diagram": ("Teff", "logg"),
    #     "1.3 Chemical plane": ("m_H", "alpha_m"),
    #     "1.4 Fit quality": ("snr", "rchi2"),
    #     "2.1 Effective temperature": ("dr17_teff", "Teff"),
    #     "2.2 Surface gravity": ("dr17_logg", "logg"),
    #     "2.3 Metallicity": ("dr17_x_h", "m_H"),
    # }
    #
    #
    #
    # --- selection helpers (lived in the `helpers` cell) ---
    #
    # def select_buttons(number):
    #     """Buttons that hand a plot's axes to the shared selector at the top.
    #
    #     Returned unbound on purpose: the calling cell must wrap these in
    #     `mo.ui.array` and bind that to a global name, or the clicks go nowhere.
    #     """
    #     _out = []
    #     for _s in TOC[SECTION_TITLES[number]]:
    #         _axes = SELECTABLE.get(_s)
    #         _out.append(
    #             mo.ui.button(
    #                 label="select",
    #                 value=0,
    #                 disabled=_axes is None,
    #                 on_click=lambda _c: (_c or 0) + 1,
    #                 on_change=(
    #                     lambda _c, _t=_s, _a=_axes: (
    #                         set_active((_t, _a[0], _a[1])),
    #                         set_view(None),
    #                     )
    #                 ),
    #             )
    #         )
    #     return _out
    #
    #
    # def selected_rows(element, curve_rows):
    #     """Row positions of the points inside the box or lasso.
    #
    #     The frontend keeps only x/y/curveNumber/pointNumber/pointIndex on a selected
    #     point -- `customdata` is dropped on the way back -- so the reliable route to
    #     a star is the trace it came from plus its index within that trace, mapped
    #     through the row positions that built the trace. `curve_rows` is
    #     {trace index: row positions}; traces not in it (the heatmap) are ignored.
    #     """
    #     _rows = []
    #     for _p in getattr(element, "value", None) or []:
    #         if not isinstance(_p, dict):
    #             continue
    #         _map = curve_rows.get(_p.get("curveNumber"))
    #         if _map is None:
    #             continue
    #         _i = _p.get("pointIndex")
    #         if _i is None:
    #             _i = _p.get("pointNumber")
    #         try:
    #             _i = int(_i)
    #         except (TypeError, ValueError):
    #             continue
    #         if 0 <= _i < len(_map):
    #             _rows.append(int(_map[_i]))
    #     return np.unique(_rows) if _rows else np.empty(0, dtype=np.int64)
    #
    #
    # def view_ranges(element):
    #     """What the reader is currently looking at, so the kernel can re-bin to it.
    #
    #     Zoom and pan push {"xaxis": {"range": [...]}, ...} into the element's value
    #     and marimo keeps that payload verbatim. Returns the ranges, the string
    #     "auto" when the reader hit autoscale, or None when the view is unchanged.
    #     """
    #     _d = getattr(element, "_selection_data", None) or {}
    #     _x, _y = _d.get("xaxis") or {}, _d.get("yaxis") or {}
    #     if _x.get("autorange") or _y.get("autorange"):
    #         return "auto"
    #     if not _x.get("range") or not _y.get("range"):
    #         return None
    #     return (
    #         [float(_v) for _v in _x["range"]],
    #         [float(_v) for _v in _y["range"]],
    #     )
    #
    #
    # def overlays(x, y, ids, subsets):
    #     """Traces marking each visible saved subset.
    #
    #     `x`, `y` and `ids` must be aligned with each other -- in this notebook that
    #     means aligned with `df` -- so the same subset can be shown on any pair of
    #     columns, or on a residual, without re-deriving it.
    #     """
    #     _traces = []
    #     for _name, _s in subsets.items():
    #         if not _s.get("visible", True):
    #             continue
    #         _m = np.isin(ids, _s["ids"])
    #         if not _m.any():
    #             continue
    #         _traces.append(
    #             go.Scattergl(
    #                 x=np.asarray(x, float)[_m],
    #                 y=np.asarray(y, float)[_m],
    #                 customdata=np.asarray(ids)[_m],
    #                 mode="markers",
    #                 marker=dict(
    #                     size=7, color=_s["color"], line=dict(width=0.6, color="white")
    #                 ),
    #                 name=f"{_name} ({int(_m.sum())})",
    #                 hovertemplate="sdss_id %{customdata}<br>%{x:.5g}, %{y:.5g}<extra></extra>",
    #                 showlegend=True,
    #             )
    #         )
    #     return _traces
    #
    # --- `section()` also took a `selectors=` mo.ui.array and, for any subsection
    # listed in SELECTABLE, rendered its button beside the accordion:
    #
    #     if selectors is not None and _s in SELECTABLE:
    #         _item = mo.hstack(
    #             [_item, selectors[_i]], widths=[12, 1], align='center', gap=0.5
    #         )
    #
    # and each section cell bound its own array before calling it:
    #
    #     sec1_select = mo.ui.array(select_buttons(1))
    #     section(1, [...], selectors=sec1_select)
    #
    # Figures added `marks(x, y)` to their trace list to draw the overlays, and
    # `one_to_one` took them through its `overlay=` argument, which is still there.
    # `density(..., return_rows=True)` is also still there: it returns the row
    # positions behind the scatter trace, which is what made a selection
    # resolvable. ---
    #
    #
    # --- cell `subset_store` ---
    #
    # # Named subsets live in reactive state so that a click in one figure can change
    # # what every other figure draws. allow_self_loops lets the control bar rebuild
    # # its own buttons after you hide or delete a subset.
    # get_subsets, set_subsets = mo.state({}, allow_self_loops=True)
    #
    # # The zoomed view of the selector, so the kernel can re-bin to whatever is on
    # # screen. This one must NOT self-loop: the watcher cell writes it.
    # get_view, set_view = mo.state(None)
    #
    # SUBSET_COLORS = (
    #     "#e45756", "#4c78a8", "#f58518", "#54a24b",
    #     "#b279a2", "#9d755d", "#ff9da6", "#17becf",
    # )
    #
    # # Which plot is loaded into the shared selector: (subsection title, x, y).
    # get_active, set_active = mo.state(None)
    #
    # # Whether the "name this subset" prompt is showing.
    # get_naming, set_naming = mo.state(False, allow_self_loops=True)
    #
    # --- cell `marks_cell` ---
    #
    # sample_ids = df["sdss_id"].to_numpy()
    #
    #
    # def marks(x, y):
    #     """Overlay traces for the visible saved subsets.
    #
    #     `x` and `y` are arrays aligned with `df`, so this works on any pair of
    #     columns and on residuals. Every figure that calls it re-runs when a subset
    #     is added, hidden or deleted -- which is the point.
    #     """
    #     return overlays(x, y, sample_ids, get_subsets())
    #
    # --- cell `subset_bar` ---
    #
    # def _toggle(name):
    #     def _cb(_count):
    #         set_subsets(
    #             lambda _s: {
    #                 _k: ({**_v, "visible": not _v["visible"]} if _k == name else _v)
    #                 for _k, _v in _s.items()
    #             }
    #         )
    #
    #     return _cb
    #
    #
    # def _drop(name):
    #     def _cb(_count):
    #         set_subsets(lambda _s: {_k: _v for _k, _v in _s.items() if _k != name})
    #
    #     return _cb
    #
    #
    # def _bump(_count):
    #     return (_count or 0) + 1
    #
    #
    # _subsets = get_subsets()
    # _names = list(_subsets)
    #
    # # A UI element is only interactive if it is bound to a global; buttons built
    # # inline inside a container are rendered but never reach the kernel. mo.ui.array
    # # is how you bind a number of them that is not known ahead of time.
    # subset_buttons = mo.ui.array(
    #     [
    #         mo.ui.button(
    #             label="hide" if _subsets[_n]["visible"] else "show",
    #             value=0, on_click=_bump, on_change=_toggle(_n),
    #         )
    #         for _n in _names
    #     ]
    #     + [
    #         mo.ui.button(
    #             label="delete", kind="danger",
    #             value=0, on_click=_bump, on_change=_drop(_n),
    #         )
    #         for _n in _names
    #     ]
    #     + [
    #         mo.ui.button(
    #             label="delete all", kind="danger",
    #             value=0, on_click=_bump, on_change=lambda _c: set_subsets({}),
    #         )
    #     ]
    # )
    #
    # if not _names:
    #     _bar = mo.md(
    #         "*No saved subsets yet. Lasso or box some points in "
    #         "[8. Selection explorer](#8-selection-explorer) and give them a name.*"
    #     )
    # else:
    #     _rows = []
    #     for _i, _n in enumerate(_names):
    #         _s = _subsets[_n]
    #         _rows.append(
    #             mo.hstack(
    #                 [
    #                     mo.Html(
    #                         f'<span style="display:inline-block;width:0.85rem;'
    #                         f"height:0.85rem;border-radius:2px;background:{_s['color']};"
    #                         f'opacity:{1.0 if _s["visible"] else 0.25}"></span>'
    #                     ),
    #                     mo.md(f"**{_n}** &nbsp;{len(_s['ids']):,} stars"),
    #                     subset_buttons[_i],
    #                     subset_buttons[len(_names) + _i],
    #                 ],
    #                 justify="start", align="center", gap=0.6,
    #             )
    #         )
    #     _rows.append(subset_buttons[-1])
    #     _bar = mo.vstack(_rows, gap=0.3)
    #
    # mo.vstack([mo.md("**Saved subsets**"), _bar], gap=0.3)
    #
    # --- cell `selector_head` ---
    #
    # heading("Selector")
    #
    # --- cell `sel_chart` ---
    #
    # # A UI element created inside `mo.lazy` cannot drive other cells, and would not
    # # even be interactive: lazy output is rendered HTML with no cell-level binding.
    # # So there is exactly one selector, here, and every plot hands its axes to it.
    # _active = get_active()
    # if _active is None:
    #     sel = None
    #     curve_rows = {}
    #     _out = mo.md(
    #         "*Nothing loaded. Press **select** beside any plot below to explore it "
    #         "here, then box or lasso points to save them as a named subset.*"
    #     )
    # else:
    #     _title, _xc, _yc = _active
    #     _xv = df[_xc].to_numpy(float)
    #     _yv = df[_yc].to_numpy(float)
    #
    #     _view = get_view()
    #     if _view is None:
    #         _xl = list(np.percentile(_xv[np.isfinite(_xv)], [0.2, 99.8]))
    #         _yl = list(np.percentile(_yv[np.isfinite(_yv)], [0.2, 99.8]))
    #         if _xc in REVERSED_AXES:
    #             _xl = _xl[::-1]
    #         if _yc in REVERSED_AXES:
    #             _yl = _yl[::-1]
    #     else:
    #         _xl, _yl = _view
    #
    #     # 200 bins across whatever is on screen: zooming in buys real resolution,
    #     # and shrinking the bins turns more stars into individually selectable points.
    #     _traces, _sparse_rows = density(
    #         _xv, _yv, bins=200, extent=[sorted(_xl), sorted(_yl)],
    #         ids=df["sdss_id"].to_numpy(), return_rows=True,
    #     )
    #     # trace 0 is the heatmap, trace 1 the selectable points; subset overlays
    #     # follow, and are already saved, so they are not resolvable.
    #     curve_rows = {1: _sparse_rows}
    #
    #     _fig = figure(
    #         _traces + marks(_xv, _yv),
    #         height=560,
    #         xlabel=AXES.get(_xc, _xc),
    #         ylabel=AXES.get(_yc, _yc),
    #         xrange=_xl,
    #         yrange=_yl,
    #         title=_title,
    #     )
    #     _fig.update_layout(dragmode="select")
    #     sel = mo.ui.plotly(_fig)
    #     _out = sel
    #
    # _out
    #
    # --- cell `sel_zoom` ---
    #
    # # Zoom or pan the selector and the kernel re-bins to the new limits. The guard
    # # matters: without it, the rebuilt figure's own relayout would feed straight
    # # back in here and the two cells would ping-pong.
    # if sel is not None:
    #     _v = view_ranges(sel)
    #     if _v == "auto":
    #         set_view(None)
    #     elif _v is not None:
    #         _old = get_view()
    #         _moved = _old is None or any(
    #             abs(_n - _o) > 0.02 * (abs(_or[1] - _or[0]) or 1.0)
    #             for _n, _o, _or in zip(
    #                 _v[0] + _v[1], _old[0] + _old[1], [_old[0]] * 2 + [_old[1]] * 2
    #             )
    #         )
    #         if _moved:
    #             set_view(_v)
    #
    # --- cell `sel_name` ---
    #
    # subset_name = mo.ui.text(
    #     placeholder="name this group", label="", debounce=True, full_width=False
    # )
    #
    # --- cell `sel_save` ---
    #
    # mo.stop(sel is None, mo.md(""))
    #
    # _ids = sample_ids[selected_rows(sel, curve_rows)]
    #
    #
    # def _bump(_count):
    #     return (_count or 0) + 1
    #
    #
    # def _create(_count):
    #     if len(_ids) == 0:
    #         return
    #     _existing = get_subsets()
    #     _name = (subset_name.value or "").strip() or f"subset {len(_existing) + 1}"
    #     _color = SUBSET_COLORS[len(_existing) % len(SUBSET_COLORS)]
    #     set_subsets(
    #         lambda _s: {**_s, _name: {"ids": _ids, "color": _color, "visible": True}}
    #     )
    #     set_naming(False)
    #
    #
    # # marimo 0.23 has no dialog primitive, so the prompt is inline: the button opens
    # # it, and it collapses back to the button on Create or Cancel.
    # start_button = mo.ui.button(
    #     label=f"create subset ({len(_ids):,})",
    #     kind="success",
    #     disabled=len(_ids) == 0,
    #     value=0,
    #     on_click=_bump,
    #     on_change=lambda _c: set_naming(True),
    # )
    # create_button = mo.ui.button(
    #     label="Create", kind="success", value=0, on_click=_bump, on_change=_create
    # )
    # cancel_button = mo.ui.button(
    #     label="Cancel", value=0, on_click=_bump, on_change=lambda _c: set_naming(False)
    # )
    #
    # mo.hstack(
    #     [subset_name, create_button, cancel_button]
    #     if get_naming()
    #     else [start_button],
    #     justify="end",
    #     align="center",
    #     gap=0.6,
    # )
    #
    # --- cell `sel_result` ---
    #
    # mo.stop(sel is None, mo.md(""))
    #
    # selected = df.iloc[selected_rows(sel, curve_rows)].reset_index(drop=True)
    # mo.stop(len(selected) == 0, mo.md(""))
    #
    #
    # def _crossfilter():
    #     """Where the current selection and the saved subsets sit elsewhere."""
    #     _panels = [("Teff", "logg"), ("m_H", "alpha_m"), ("snr", "rchi2")]
    #     _fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.07)
    #     for _i, (_xc, _yc) in enumerate(_panels, start=1):
    #         _x, _y = df[_xc].to_numpy(float), df[_yc].to_numpy(float)
    #         _xr = list(np.percentile(_x[np.isfinite(_x)], [0.2, 99.8]))
    #         _yr = list(np.percentile(_y[np.isfinite(_y)], [0.2, 99.8]))
    #         if _xc in REVERSED_AXES:
    #             _xr = _xr[::-1]
    #         if _yc in REVERSED_AXES:
    #             _yr = _yr[::-1]
    #
    #         for _t in density(_x, _y, bins=170, extent=[sorted(_xr), sorted(_yr)],
    #                           colorscale="Greys"):
    #             _t.update(showlegend=False)
    #             _fig.add_trace(_t, row=1, col=_i)
    #         for _t in marks(_x, _y):
    #             _t.update(showlegend=_i == 1)
    #             _fig.add_trace(_t, row=1, col=_i)
    #
    #         _m = np.isin(sample_ids, selected["sdss_id"].to_numpy())
    #         _fig.add_trace(
    #             go.Scattergl(
    #                 x=_x[_m], y=_y[_m], customdata=sample_ids[_m], mode="markers",
    #                 marker=dict(size=7, color="#000000", symbol="circle-open"),
    #                 name="current selection",
    #                 hovertemplate="sdss_id %{customdata}<extra></extra>",
    #                 showlegend=_i == 1,
    #             ),
    #             row=1, col=_i,
    #         )
    #         _fig.update_xaxes(range=_xr, title_text=AXES.get(_xc, _xc), row=1, col=_i)
    #         _fig.update_yaxes(range=_yr, title_text=AXES.get(_yc, _yc), row=1, col=_i)
    #
    #     _fig.update_layout(
    #         height=440, template="plotly_white", dragmode="pan", hovermode="closest",
    #         uirevision="keep", showlegend=True,
    #         legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
    #         margin=dict(l=60, r=20, t=40, b=48),
    #     )
    #     return _fig
    #
    #
    # selected_table = mo.ui.table(
    #     selected[
    #         ["sdss_id", "apogee_id", "Teff", "logg", "m_H", "alpha_m", "vsini",
    #          "snr", "rchi2", "flags"]
    #     ],
    #     selection=None,
    #     page_size=10,
    # )
    #
    # mo.vstack(
    #     [mo.md("**Where the selection sits**"), _crossfilter(), selected_table],
    #     gap=0.6,
    # )

    mo.md("""
    *Lasso/box selection and named subsets are parked in this cell, commented out. Split the blocks back into their own cells to restore them.*
    """)
    return


if __name__ == "__main__":
    app.run()
