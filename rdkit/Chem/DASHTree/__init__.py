#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Opening a DASH tree, and reading what it matched.

A full tree is a few hundred megabytes of published data, so it is not
distributed with RDKit. GetDASHTree() opens a local container, or fetches one
from a URL once into a cache directory and opens it from there. Given nothing,
it uses $RDKIT_DASH_TREE.

The container has to be in RDKit's own '.dash' format: python -m
rdkit.Chem.DASHTree.convert writes one from a legacy tree.

The rest turns the nodes a match descended through back into the substructure
they stand for: NodePathToSmarts() writes it as a pattern, and
GetAtomMatchSmarts() does the match and the writing in one call.
"""

import hashlib
import os
import shutil
import sys
import uuid
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from rdkit.Chem import rdDASHTree


def _cacheDir():
  """The directory for downloaded containers. $RDKIT_DASH_DIR overrides it."""
  fromEnv = os.environ.get("RDKIT_DASH_DIR")
  if fromEnv:
    return Path(fromEnv)
  base = os.environ.get("LOCALAPPDATA") or os.environ.get("XDG_CACHE_HOME")
  return Path(base) / "rdkit-dash" if base else Path.home() / ".cache/rdkit-dash"


def _fetch(url, cacheDir):
  """Returns the cached copy of url, and downloads it first if there is none."""
  # the name starts with a hash of url because two hosts can both serve a
  # 'default.dash'. A second download is better than a wrong tree with no error.
  name = Path(urlparse(url).path).name or "tree.dash"
  dest = Path(cacheDir) / f"{hashlib.sha256(url.encode()).hexdigest()[:8]}-{name}"
  if dest.is_file():
    return dest

  dest.parent.mkdir(parents=True, exist_ok=True)
  # each call uses its own name, so concurrent first-time calls cannot share a download
  partial = dest.with_name(f"{dest.name}.{uuid.uuid4().hex}.part")
  print(f"download {url}\n      to {dest}", file=sys.stderr)
  try:
    # some archives reject the plain urllib user agent but accept RDKit's
    with urlopen(Request(url, headers={"User-Agent": "RDKit"}), timeout=60) as response:
      with open(partial, "wb") as out:
        shutil.copyfileobj(response, out)
    rdDASHTree.DASHTree(str(partial))  # urllib does not notice a connection that closes early
    # only now does it become the cached copy, so an interrupted download
    # cannot leave a truncated container in the cache for later calls
    partial.rename(dest)
  except OSError:
    # a concurrent call made dest first, and windows does not rename onto it
    if not dest.is_file():
      raise
  finally:
    partial.unlink(missing_ok=True)
  return dest


def GetDASHTree(source=None, properties=None, cacheDir=None):
  """Opens a DASH tree, and downloads it first if it is not on disk.

    ARGUMENTS:
      - source: path to a '.dash' container, or an http(s) URL of one.
        Defaults to $RDKIT_DASH_TREE.
      - properties: the property columns to resolve, None meaning every
        column the file has.
      - cacheDir: the directory for downloads. Defaults to $RDKIT_DASH_DIR
        or a per-user cache directory.

    RETURNS: an rdDASHTree.DASHTree
    """
  if source is None:
    source = os.environ.get("RDKIT_DASH_TREE")
  if not source:
    raise ValueError("no DASH tree given: pass a container or set RDKIT_DASH_TREE")
  # a windows path has a single letter scheme, so test for the ones we serve
  if urlparse(str(source)).scheme in ("http", "https"):
    source = _fetch(str(source), cacheDir or _cacheDir())
  return rdDASHTree.DASHTree(str(source), properties)


#: the SMARTS bond for each DASH bond descriptor. DASH records a conjugated bond
#: as 4 with its order erased, and SMARTS has no conjugation primitive. So 4
#: becomes the wildcard, and the pattern matches everything the tree node
#: matches, never less.
_bondSmarts = {1: "-", 2: "=", 3: "#", 4: "~"}


def _atomSmarts(feature, mapNum=0):
  """One DASH atom-feature class as a SMARTS atom.

    Element, total connectivity (X) and total hydrogen count (H) are exact. X
    and H count hydrogens whether or not they are explicit atoms, so the pattern
    reads the same on either kind of molecule. The class's conjugation flag has
    no SMARTS primitive, so the pattern omits it.
    """
  atomicNum, degree, charge, _conjugated, numHs = feature
  mapped = f":{mapNum}" if mapNum else ""
  return f"[#{atomicNum}X{degree}H{numHs}{charge:+d}{mapped}]"


def NodePathToSmarts(tree, nodePath, foldHydrogens=True):
  """The substructure a node path describes, as SMARTS.

    Each node adds one atom bonded to one atom already matched, so a path is a
    tree of atoms and its SMARTS has no ring closures. The atom the path is
    about carries atom map 1, so a match of the pattern says which target atom
    that is.

    ARGUMENTS:
      - tree: the rdDASHTree.DASHTree the path came from
      - nodePath: what GetAtomNodePath returned, the branch index followed by
        node ids
      - foldHydrogens: omit hydrogen nodes, because their neighbour's H count
        already requires them. The pattern then also works on molecules without
        explicit hydrogens. A hydrogen the path is *about* stays, with atom map
        1, and that pattern needs explicit hydrogens.

    The pattern cannot say that an atom or a bond is conjugated, because SMARTS
    has no primitive for it. A conjugated bond becomes '~'. The pattern matches
    everything the node matches, and possibly more.
    """
  if len(nodePath) < 2:
    raise ValueError("an empty path describes no substructure")
  branch, ids = nodePath[0], list(nodePath[1:])
  nodes = [tree.GetNode(branch, i) for i in ids]

  atoms, children = [], []  # per written atom: its SMARTS, its (bond, child)s

  def addAtom(feature, mapNum=0):
    atoms.append(_atomSmarts(feature, mapNum))
    children.append([])
    return len(atoms) - 1

  root = tree.GetAtomFeature(nodes[0].GetAtomFeatureIndex())
  if root[0] == 1:
    # the tree matches a hydrogen through its heavy neighbour. The path uses its
    # first step to reach that neighbour, and later steps count positions from it
    if len(nodes) < 2:
      raise ValueError("a hydrogen's path has to reach its heavy neighbour")
    hydrogen = addAtom(root, 1)
    heavy = addAtom(tree.GetAtomFeature(nodes[1].GetAtomFeatureIndex()))
    children[hydrogen].append(("-", heavy))
    positions, rest = [heavy], nodes[2:]
  else:
    positions, rest = [addAtom(root, 1)], nodes[1:]

  for node in rest:
    feature = tree.GetAtomFeature(node.GetAtomFeatureIndex())
    if foldHydrogens and feature[0] == 1:
      positions.append(None)  # its neighbour's H count already says so
      continue
    conAtom, bond = node.GetConAtom(), _bondSmarts.get(node.GetConType())
    if not 0 <= conAtom < len(positions) or bond is None:
      raise ValueError(f"node {node.GetId()} of branch {branch} attaches "
                       f"as ({conAtom}, {node.GetConType()}), which is "
                       f"not a position and a bond descriptor 1-4")
    parent = positions[conAtom]
    if parent is None:
      raise ValueError("a node attaches to a folded hydrogen")
    atom = addAtom(feature)
    children[parent].append((bond, atom))
    positions.append(atom)

  def write(atom):
    out = atoms[atom]
    for bond, child in children[atom][:-1]:
      out += "(" + bond + write(child) + ")"
    if children[atom]:
      bond, child = children[atom][-1]
      out += bond + write(child)
    return out

  return write(0)


def GetAtomMatchSmarts(tree, mol, atomIdx, params=None, foldHydrogens=True):
  """SMARTS for the substructure the tree matched around one atom.

    Matches the atom the way GetPartialCharges does and writes the nodes it
    descended through as a pattern. See NodePathToSmarts.
    """
  path = tree.GetAtomNodePath(mol, atomIdx, params or rdDASHTree.DASHParams())
  if len(path) < 2:
    raise ValueError(f"no node in the tree describes atom {atomIdx}")
  return NodePathToSmarts(tree, path, foldHydrogens)
