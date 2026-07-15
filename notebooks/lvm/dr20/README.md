# LVM DR20 DAP / DRP tutorials

Tutorials for working with SDSS-V LVM DR20 data products (DRP `SFrame` and DAP
files), by _Sebastian F. Sanchez and the LVM team_. See the
[LVM notebooks README](../README.md) for a description of each notebook.

## Getting the data

DR20 is not public yet, so the SAS asks for the SDSS collaboration credentials
before it will serve these files.

The notebooks fetch data through [`binder_access.py`](binder_access.py), which
exposes the same `Access` API as `sdss_access` but downloads over HTTPS:

```python
from binder_access import Access
access = Access(release='dr20', verbose=True)
```

The first cell that downloads something will prompt once per kernel for the SDSS
password (the collaboration-wide `sdss5` account). To avoid the prompt, either
create a `~/.netrc` with an entry for `data.sdss5.org`:

```
machine data.sdss5.org
login sdss5
password <the collaboration password>
```

(then `chmod 600 ~/.netrc`), or set `SDSS_USERNAME` / `SDSS_PASSWORD` in the
environment before starting the kernel.

Files are written under `$SAS_BASE_DIR`, which defaults to `~/sas`. That lives in
your home directory, so downloads survive between sessions and you only pay for
them once. **They are large** — the Helix `SFrame` is ~0.5 GB and the DAP model
file is ~0.8 GB — so expect the first run of a notebook to take a few minutes.

## Why not `sdss_access` directly?

`sdss_access.Access` shells out to `rsync` on Linux. That does not work on the
BinderHub image for two reasons:

1. There is no `rsync` binary in the image.
2. `sdss_access` replaces the whole subprocess environment with
   `{"RSYNC_PASSWORD": ...}` when it shells out, so `PATH` is lost and `rsync`
   would not be found even if it were installed somewhere non-standard.

The SAS serves the same files over HTTPS, so `binder_access.py` sidesteps both
problems and needs no setup in the container.

## Dependencies

Each notebook opens with a _Required Python dependencies_ (or _Dependencies_)
section listing what it needs, and some include a `pip install` or `conda
install` line to set that up.

On the BinderHub you can skip that step — everything is already installed in the
image via the repository's [`requirements.txt`](../../../requirements.txt), which
is what the Binder builds from. Those sections are there for running the
notebooks outside the Binder, where you do need to install things yourself.

One package is worth a special mention on the Binder: `ipympl`, below.

## Interactive figures

The notebooks that use `%matplotlib widget` need `ipympl`, which is listed in the
repository's [`requirements.txt`](../../../requirements.txt) so that the
`jupyter-matplotlib` frontend extension is present when the image is built.
Installing `ipympl` from a terminal in a running session does **not** work: the
extension is only registered at server start, so you would get

```
Failed to load model class 'MPLCanvasModel' from module 'jupyter-matplotlib'
```

If you see that error, the image predates the `ipympl` requirement — restart your
server so it picks up the current environment.
