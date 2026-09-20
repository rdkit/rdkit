#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Opening a DASH tree that lives somewhere else, and reading what it matched.

A tree is a few hundred megabytes of published data, so it is not distributed
with RDKit. GetDASHTree() takes either a local container or a URL to one,
fetches it once into a cache directory, and opens it from there.

The container has to be in RDKit's own '.dash' format already; nothing here
converts a tree. Code/GraphMol/DASHTree/tools/dash_convert.py does that.

The rest turns the nodes a match descended through back into the substructure
they stand for: NodePathToSmarts() writes it as a pattern, and
GetAtomMatchSmarts() does the match and the writing in one call.
"""

import hashlib
import os
import shutil
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from rdkit import Chem
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


#: the SMARTS bond for each DASH bond descriptor. DASH records a conjugated bond
#: as 4 with its order erased, and SMARTS has no conjugation primitive, so that
#: one becomes the wildcard: the pattern then matches everything the tree node
#: matches, never less.
_bondSmarts = {1: "-", 2: "=", 3: "#", 4: "~"}


def _atomSmarts(feature, mapNum=0):
    """One DASH atom-feature class as a SMARTS atom.

    Element, total connectivity (X) and total hydrogen count (H) are exact, and
    both count hydrogens whether or not they are explicit atoms, so the pattern
    reads the same on either kind of molecule. The class's conjugation flag has
    no SMARTS primitive and is left out.
    """
    atomicNum, degree, charge, _conjugated, numHs = feature
    mapped = f":{mapNum}" if mapNum else ""
    return f"[#{atomicNum}X{degree}H{numHs}{charge:+d}{mapped}]"


def _writePattern(tree, nodePath, foldHydrogens):
    """The SMARTS for a node path, and the order its atoms were written in.

    Returns (smarts, order): order[q] is the path position of the q-th atom in
    the string, which is also the q-th atom of Chem.MolFromSmarts(smarts).
    """
    branch, ids = nodePath[0], list(nodePath[1:])
    if not ids:
        raise ValueError("an empty path describes no substructure")
    nodes = [tree.GetNode(branch, i) for i in ids]

    atoms, children = [], []  # per written atom: its SMARTS, its (bond, child)s

    def addAtom(feature, mapNum=0):
        atoms.append(_atomSmarts(feature, mapNum))
        children.append([])
        return len(atoms) - 1

    root = nodes[0].GetFeature()
    if root[0] == 1:
        # a hydrogen is matched through its heavy neighbour: the path spends its
        # first step getting there, and later steps count positions from it
        if len(nodes) < 2:
            raise ValueError("a hydrogen's path has to reach its heavy neighbour")
        hydrogen = addAtom(root, 1)
        heavy = addAtom(nodes[1].GetFeature())
        children[hydrogen].append(("-", heavy))
        positions, rest = [heavy], nodes[2:]
    else:
        positions, rest = [addAtom(root, 1)], nodes[1:]

    for node in rest:
        feature = node.GetFeature()
        if foldHydrogens and feature[0] == 1:
            positions.append(None)  # its neighbour's H count already says so
            continue
        conAtom, bond = node.GetConAtom(), _bondSmarts.get(node.GetConType())
        if conAtom < 0 or bond is None:
            raise ValueError(f"node {node.GetId()} of branch {branch} attaches "
                             f"as ({conAtom}, {node.GetConType()}), which is "
                             f"not a position and a bond descriptor 1-4")
        parent = positions[conAtom]
        if parent is None:
            raise ValueError("a node attaches to a hydrogen that was folded away")
        atom = addAtom(feature)
        children[parent].append((bond, atom))
        positions.append(atom)

    # written depth first, which is not the path order once a branch has a
    # subtree of its own -- so the order is recorded for NodePathToQueryMol
    order = []

    def write(atom):
        order.append(atom)
        out = atoms[atom]
        for bond, child in children[atom][:-1]:
            out += "(" + bond + write(child) + ")"
        if children[atom]:
            bond, child = children[atom][-1]
            out += bond + write(child)
        return out

    return write(0), order


def NodePathToSmarts(tree, nodePath, foldHydrogens=True):
    """The substructure a node path describes, as SMARTS.

    Each node adds one atom bonded to one atom already matched, so a path is a
    tree of atoms and writes out without ring closures. The atom the path is
    about carries atom map 1, so a match of the pattern says which target atom
    it landed on.

    ARGUMENTS:
      - tree: the rdDASHTree.DASHTree the path came from
      - nodePath: what GetAtomNodePath returned, the branch index followed by
        node ids
      - foldHydrogens: leave hydrogen nodes out, since their neighbour's H count
        already requires them; the pattern then also works on molecules without
        explicit hydrogens. A hydrogen the path is *about* stays, as [#1:1],
        and needs them.

    What the pattern cannot say is that an atom or a bond is conjugated: SMARTS
    has no primitive for it, so conjugated bonds are written as '~'. It matches
    everything the node matches and possibly more; it never misses.
    """
    return _writePattern(tree, nodePath, foldHydrogens)[0]


def NodePathToQueryMol(tree, nodePath, foldHydrogens=True):
    """NodePathToSmarts, parsed into a query molecule, with its atoms in path
    order: atom k is path node k -- and so, with foldHydrogens=False, the atom
    that GetMatchedSubstructure reports at position k. Folded hydrogens are
    simply absent.
    """
    smarts, order = _writePattern(tree, nodePath, foldHydrogens)
    query = Chem.MolFromSmarts(smarts)
    newOrder = [0] * len(order)  # newOrder[k] = which parsed atom becomes atom k
    for parsed, written in enumerate(order):
        newOrder[written] = parsed
    return Chem.RenumberAtoms(query, newOrder)


def GetAtomMatchSmarts(tree, mol, atomIdx, params=None, foldHydrogens=True):
    """SMARTS for the substructure the tree matched around one atom.

    Matches the atom the way GetPartialCharges would and writes the nodes it
    descended through as a pattern; see NodePathToSmarts.
    """
    if params is None:
        path = tree.GetAtomNodePath(mol, atomIdx)
    else:
        path = tree.GetAtomNodePath(mol, atomIdx, params)
    if len(path) < 2:
        raise ValueError(f"no node in the tree describes atom {atomIdx}")
    return NodePathToSmarts(tree, path, foldHydrogens)
