#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Opening a DASH tree that lives somewhere else.

A tree is a few hundred megabytes of published data, so it is not distributed
with RDKit. GetDASHTree() takes either a local container or a URL to one,
fetches it once into a cache directory, and opens it from there.

The container has to be in RDKit's own '.dash' format already; nothing here
converts a tree. Code/GraphMol/DASHTree/tools/dash_convert.py does that.
"""

import hashlib
import os
import shutil
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from rdkit.Chem import rdDASHTree

#: fetched when GetDASHTree() is given no source at all
defaultTreeUrl = "http://localhost:8000/default.dash"  # TODO: replace with actual default URL. For example denmarc.ethz.ch/dash/default.dash


def _cacheDir():
    """Where fetched containers are kept. $RDKIT_DASH_DIR overrides it."""
    fromEnv = os.environ.get("RDKIT_DASH_DIR")
    if fromEnv:
        return Path(fromEnv)
    base = os.environ.get("LOCALAPPDATA") or os.environ.get("XDG_CACHE_HOME")
    return Path(base) / "rdkit-dash" if base else Path.home() / ".cache/rdkit-dash"


def _fetch(url, cacheDir):
    """The cached copy of url, downloading it if it is not there yet."""
    # the name is hashed because two hosts can both serve a 'default.dash',
    # and quietly opening the wrong tree is worse than downloading twice
    name = Path(urlparse(url).path).name or "tree.dash"
    dest = Path(cacheDir) / f"{hashlib.sha256(url.encode()).hexdigest()[:8]}-{name}"
    if dest.is_file():
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name + ".part")
    print(f"fetching {url}\n      to {dest}")
    # a plain urllib user agent is rejected by some archives, RDKit's is not
    with urlopen(Request(url, headers={"User-Agent": "RDKit"})) as response:
        with open(partial, "wb") as out:
            shutil.copyfileobj(response, out)
    # only now does it become the cached copy, so an interrupted download
    # cannot leave a truncated container behind to be opened later
    partial.replace(dest)
    return dest


def GetDASHTree(source=None, properties=None, prefetch=False, cacheDir=None):
    """Opens a DASH tree, downloading it first if source is a URL.

    ARGUMENTS:
      - source: path to a '.dash' container, or an http(s) URL of one.
        Defaults to rdkit.Chem.DASHTree.defaultTreeUrl.
      - properties: the property columns to resolve, None meaning every
        column the file has. Naming the ones you want keeps the rest off
        the heap.
      - prefetch: read the mapped arrays through once, worth it before a
        large batch and pointless for a handful of molecules.
      - cacheDir: where downloads are kept, defaulting to $RDKIT_DASH_DIR
        or a per-user cache directory.

    RETURNS: an rdDASHTree.DASHTree
    """
    source = source or defaultTreeUrl
    # a windows path has a single letter scheme, so test for the ones we serve
    if urlparse(str(source)).scheme in ("http", "https"):
        source = _fetch(str(source), cacheDir or _cacheDir())
    return rdDASHTree.DASHTree(str(source), properties, prefetch)
