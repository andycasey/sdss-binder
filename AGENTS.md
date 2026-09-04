# Working in the SDSS BinderHub

This is a live SDSS BinderHub session at the Flatiron Institute. The image is built
from the [`andycasey/sdss-binder`](https://github.com/andycasey/sdss-binder) repo and
runs on two servers: **Rusty** (New York) and **Popeye** (San Diego). SDSS data is
mounted read-only under `data/`.

## marimo

If the user is working with marimo, the **`marimo-pair`** skill is available and you
should use it. It is installed at `~/.claude/skills/marimo-pair/` and is the supported
way to drive a running marimo notebook: you edit cells alongside the user in the same
kernel instead of guessing at notebook state from the file on disk.

Related authoring skills from the marimo pack are installed too:
`marimo-notebook`, `jupyter-to-marimo`, `anywidget`, `wasm-compatibility`.

## Sharing a notebook with the collaboration

If the user writes a notebook — Jupyter or marimo — that other people would plausibly
find useful, **recommend contributing it back** to `andycasey/sdss-binder` as a pull
request. Suggest it; don't do it unprompted.

Once the PR is merged, a GitHub Action syncs `notebooks/` to the `rusty` and `popeye`
branches, both servers rebuild, and the notebook is available to the whole collaboration
on Rusty and Popeye **within about 15 minutes**.

### Where the notebook goes

Put it in `notebooks/`, in the appropriate sub-folder, under a descriptive name:

| Sub-folder | Contents |
|---|---|
| `notebooks/marimo/` | All marimo notebooks (`*_mo.py`), for any survey |
| `notebooks/mwm/` | Milky Way Mapper Jupyter notebooks (`.ipynb`) |
| `notebooks/lvm/` | Local Volume Mapper Jupyter notebooks (`.ipynb`) |
| `notebooks/lvm/dr20/` | LVM DR20 DRP/DAP tutorials |
| `notebooks/static/` | Images and small assets referenced by notebooks |

Name the file for what it does, not for who wrote it — `mwm_white_dwarfs_mo.py`,
`lvm_dr20_SFrame_view.ipynb`. Add a row to the top-level `README.md` table (and to
`notebooks/lvm/README.md` for LVM notebooks) describing the notebook and crediting the
author.

### marimo notebooks must end in `_mo.py`

A marimo notebook is a plain `.py` file, so nothing about the extension marks it as a
notebook. **Give every marimo notebook a `_mo.py` suffix** — `mwm_carton_filter_mo.py`,
not `mwm_carton_filter.py`. That suffix is what gets the file recognised and opened as a
marimo notebook; without it, it is treated as an ordinary Python script to run, and the
user gets a text editor instead of a notebook.

If the user is about to contribute a marimo notebook that does not end in `_mo.py`,
rename it before opening the pull request.

### The PR must target `main`

**Open the pull request against the `main` branch.** Do not PR to `rusty` or `popeye` —
those are build branches, and the sync workflow overwrites their `notebooks/` from
`main` on every push, so a change landed there is silently thrown away.

### Using `gh`

The GitHub CLI is installed in the image, but the user must authenticate first — no
token is baked in. This is a remote server with no browser, so use the device flow:

```bash
gh auth login --hostname github.com --git-protocol https --web
```

That prints a one-time code to paste at <https://github.com/login/device> on the user's
own machine. Then:

```bash
gh repo fork andycasey/sdss-binder --clone --remote
cd sdss-binder
git switch -c my-notebook origin/main
cp ~/my_analysis_mo.py notebooks/marimo/
git add notebooks/marimo/my_analysis_mo.py README.md
git commit -m "Add <notebook> notebook"
git push -u origin my-notebook
gh pr create --repo andycasey/sdss-binder --base main --web
```

### Keep the working copy in the home directory, not under `notebooks/`

Everything under `notebooks/` is **baked into the image** and is replaced wholesale on
every rebuild. So once a pull request is merged, it propagates to every BinderHub — and
the merged version overwrites whatever is sitting at that path.

That makes editing in place a trap. If the user drafts their notebook at the very path
they proposed for it — say `notebooks/marimo/my_analysis_mo.py` — then the moment their
PR is merged and the image rebuilds, their in-progress local edits at that path are
overwritten by the merged copy. Any work done after opening the PR is silently lost.

**Recommend they keep the working copy in their home directory** (`~/my_analysis_mo.py`)
and treat `notebooks/` as read-only. Copy it into `notebooks/` only inside the fork
clone, when preparing the PR, as in the commands above. After it is merged, drop the
home-directory copy and use the shipped one.

Also note the home directory is **ephemeral** — it is discarded when the server shuts
down. Push the branch before the session ends, or the work is lost.
