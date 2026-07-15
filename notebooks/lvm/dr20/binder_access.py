"""HTTPS data access for the LVM DR20 tutorials on the SDSS BinderHub.

`sdss_access.Access` shells out to ``rsync`` on Linux, which does not work on the
BinderHub image for two reasons:

1. The image has no ``rsync`` binary.
2. ``sdss_access`` replaces the whole subprocess environment with
   ``{"RSYNC_PASSWORD": ...}`` when it shells out, so ``PATH`` is lost and the
   binary would not be found even if it were installed somewhere non-standard.

The SAS serves the same files over HTTPS with basic auth, so this module provides
an `Access` class with the same API as `sdss_access.Access` (``full``, ``url``,
``remote``, ``add``, ``set_stream``, ``commit``) that downloads with `requests`
instead. Notebooks only need to swap the import:

    from binder_access import Access
    access = Access(release="dr20", verbose=True)

Credentials are taken from ``~/.netrc`` if it has an entry for data.sdss5.org,
otherwise from the ``SDSS_USERNAME``/``SDSS_PASSWORD`` environment variables,
otherwise you are prompted once per session.
"""

from __future__ import annotations

import getpass
import netrc
import os
from pathlib import Path

import requests
from sdss_access.path import Path as SDSSPath
from tqdm.auto import tqdm

__all__ = ["Access", "fetch", "get_credentials"]

MACHINE = "data.sdss5.org"
DEFAULT_USERNAME = "sdss5"

# Cached so that a notebook only ever prompts once per kernel.
_credentials: tuple[str, str] | None = None


def _credentials_from_netrc() -> tuple[str, str] | None:
    try:
        authenticators = netrc.netrc().authenticators(MACHINE)
    except (FileNotFoundError, netrc.NetrcParseError):
        return None
    if not authenticators:
        return None
    login, _, password = authenticators
    if not password:
        return None
    return (login or DEFAULT_USERNAME, password)


def get_credentials(username=None, password=None) -> tuple[str, str]:
    """Return (username, password) for data.sdss5.org, prompting if necessary."""
    global _credentials

    if username and password:
        _credentials = (username, password)
    if _credentials is not None:
        return _credentials

    _credentials = _credentials_from_netrc()
    if _credentials is not None:
        return _credentials

    username = username or os.environ.get("SDSS_USERNAME") or DEFAULT_USERNAME
    password = password or os.environ.get("SDSS_PASSWORD")
    if not password:
        password = getpass.getpass(f"SDSS password for {username}@{MACHINE}: ")
    _credentials = (username, password)
    return _credentials


def _download(url: str, destination: str, auth: tuple[str, str], verbose: bool = True) -> str:
    """Stream `url` to `destination`, writing to a partial file until complete."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.parent / (destination.name + ".part")

    with requests.get(url, auth=auth, stream=True, timeout=60) as response:
        if response.status_code == 401:
            raise PermissionError(
                f"Not authorised for {url}. Check your SDSS credentials; "
                "if you mistyped the password, restart the kernel to be re-prompted."
            )
        response.raise_for_status()
        total = int(response.headers.get("content-length", 0))
        with open(partial, "wb") as file_object, tqdm(
            total=total or None,
            unit="B",
            unit_scale=True,
            desc=destination.name,
            disable=not verbose,
            leave=False,
        ) as progress:
            for chunk in response.iter_content(chunk_size=1 << 20):
                file_object.write(chunk)
                progress.update(len(chunk))

    partial.replace(destination)
    if verbose:
        print(f"Downloaded {destination}")
    return str(destination)


class Access:
    """Drop-in replacement for `sdss_access.Access` that fetches over HTTPS.

    Only the parts of the API the LVM DR20 tutorials use are implemented.
    """

    def __init__(self, release=None, verbose=False, public=False, label=None, **kwargs):
        self.path = SDSSPath(release=release, public=public, verbose=verbose)
        self.release = release
        self.verbose = verbose
        self._queue: list[tuple[str, dict]] = []
        self._auth: tuple[str, str] | None = None

    def __repr__(self):
        return f'<Access(access_mode="https", using="{MACHINE}")>'

    # -- path resolution (delegated straight to sdss_access) ------------------
    def full(self, product, **kwargs):
        return self.path.full(product, **kwargs)

    def url(self, product, **kwargs):
        return self.path.url(product, **kwargs)

    # -- download API --------------------------------------------------------
    def remote(self, username=None, password=None, inquire=None):
        self._auth = get_credentials(username, password)

    def add(self, product, **kwargs):
        self._queue.append((product, kwargs))

    def set_stream(self):
        """No-op. Kept so the sdss_access call sequence still works."""

    def get_paths(self):
        return [self.full(product, **kwargs) for product, kwargs in self._queue]

    def get_urls(self):
        return [self.url(product, **kwargs) for product, kwargs in self._queue]

    def reset(self):
        self._queue = []

    def commit(self, offset=None, limit=None, follow_symlinks=True):
        """Download everything queued with `add`, skipping files already present."""
        auth = self._auth or get_credentials()
        for product, kwargs in self._queue:
            destination = self.full(product, **kwargs)
            if os.path.isfile(destination):
                if self.verbose:
                    print(f"Already downloaded: {destination}")
                continue
            _download(self.url(product, **kwargs), destination, auth, verbose=self.verbose)
        self._queue = []


def fetch(product, release="dr20", verbose=True, **kwargs) -> str:
    """Return the local path to `product`, downloading it if it isn't there yet."""
    access = Access(release=release, verbose=verbose)
    destination = access.full(product, **kwargs)
    if not os.path.isfile(destination):
        access.remote()
        access.add(product, **kwargs)
        access.set_stream()
        access.commit()
    return destination
