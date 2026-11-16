#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Convert a legacy DASH tree into a .dash container.

A legacy tree is how the DASH-tree python package stores one: a folder of 0.gz,
0.h5, 1.gz, 1.h5, ..., per atom-feature branch a gzipped pickle of the nodes and
an HDF5 table of their values. The .gz files are pickles, so only convert trees
from a source you trust. Needs pandas and pytables, not the DASH-tree package.

The tree must be built on the 122 atom-feature classes of the DASH-tree
package, as the published trees are. The hydrogen class (Z = 1) is checked
against the data, since the descent treats it specially: it is the only branch
whose root's children attach to nothing.
"""

import argparse
import gzip
import os
import pickle
import sys

import numpy as np

from .container import (DEFAULT_ATTENTION_THRESHOLD, DTYPE_F16, DTYPE_F32, DTYPE_F64, DTYPE_I32,
                        FEATURE_DTYPE, MAX_CHILDREN, MAX_CON_ATOM, MAX_CON_TYPE, MAX_NODES,
                        NODE_FLAG_STOP, NP_DTYPE, RECORD_DTYPE, packKey, writeContainer)

# (atomicNum, degree, formalCharge, conjugated, numHs) of the published trees'
# 122 classes, in branch order: serenityff.charge.tree.atom_features
# .AtomFeatures.feature_list. The order is part of those trees' data.
STANDARD_FEATURES = [(5, 1, 0, 0, 2), (5, 3, -1, 1, 0), (5, 3, 0, 0, 0), (5, 3, 0, 0, 2),
                     (5, 4, -1, 0, 0), (5, 4, -1, 0, 2), (35, 1, 0, 0, 0), (6, 1, -1, 0, 0),
                     (6, 1, -1, 1, 0), (6, 1, 0, 0, 1), (6, 1, 0, 0, 2), (6, 1, 0, 0, 3),
                     (6, 1, 0, 1, 1), (6, 2, -1, 1, 1), (6, 2, 0, 0, 0), (6, 2, 0, 0, 1),
                     (6, 2, 0, 0, 2), (6, 2, 0, 1, 0), (6, 2, 0, 1, 1), (6, 3, -1, 0, 0),
                     (6, 3, -1, 1, 0), (6, 3, -1, 0, 1), (6, 3, -1, 1, 1), (6, 3, 0, 0, 0),
                     (6, 3, 0, 0, 1), (6, 3, 0, 0, 2), (6, 3, 0, 1, 0), (6, 3, 0, 1, 1),
                     (6, 3, 0, 1, 2), (6, 3, 1, 0, 0), (6, 3, 1, 1, 0), (6, 4, 0, 0, 0),
                     (6, 4, 0, 0, 1), (6, 4, 0, 0, 2), (6, 4, 0, 0, 3), (17, 1, 0, 0, 0),
                     (9, 1, 0, 0, 0), (1, 1, 0, 0, 0), (53, 1, 0, 0, 0), (53, 3, 0, 0, 0),
                     (53, 3, 0, 0, 1), (53, 4, 0, 0, 0), (7, 1, -1, 0, 0), (7, 1, -1, 1, 0),
                     (7, 1, 0, 0, 0), (7, 1, 0, 1, 0), (7, 1, 0, 1, 1), (7, 1, 0, 1, 2),
                     (7, 1, 1, 0, 3), (7, 2, -1, 0, 0), (7, 2, -1, 0, 1), (7, 2, -1, 1, 0),
                     (7, 2, -1, 1, 1), (7, 2, 0, 0, 0), (7, 2, 0, 0, 1), (7, 2, 0, 1, 0),
                     (7, 2, 0, 1, 1), (7, 2, 1, 0, 0), (7, 2, 1, 0, 1), (7, 2, 1, 0, 2),
                     (7, 2, 1, 1, 0), (7, 2, 1, 1, 1), (7, 3, 0, 0, 0), (7, 3, 0, 0, 1),
                     (7, 3, 0, 0, 2), (7, 3, 0, 1, 0), (7, 3, 0, 1, 1), (7, 3, 0, 1, 2),
                     (7, 3, 1, 0, 0), (7, 3, 1, 0, 1), (7, 3, 1, 1, 0), (7, 3, 1, 1, 1),
                     (7, 4, 1, 0, 0), (7, 4, 1, 0, 1), (7, 4, 1, 0, 2), (7, 4, 1, 0, 3),
                     (8, 1, -1, 0, 0), (8, 1, -1, 1, 0), (8, 1, 0, 0, 0), (8, 1, 0, 0, 1),
                     (8, 1, 0, 1, 0), (8, 1, 0, 1, 1), (8, 2, 0, 0, 0), (8, 2, 0, 0, 1),
                     (8, 2, 0, 0, 2), (8, 2, 0, 1, 0), (8, 2, 0, 1, 1), (8, 2, 0, 1, 2),
                     (8, 2, 1, 0, 0), (8, 2, 1, 0, 1), (8, 2, 1, 1, 0), (8, 3, 1, 0, 0),
                     (8, 3, 1, 0, 1), (15, 1, 0, 0, 1), (15, 2, 0, 0, 0), (15, 2, 0, 0, 1),
                     (15, 2, 0, 1, 0), (15, 3, 0, 0, 0), (15, 4, 0, 0, 0), (15, 4, 0, 0, 1),
                     (15, 4, 0, 1, 0), (15, 4, 0, 1, 1), (15, 4, 1, 0, 0), (15, 5, 0, 0, 0),
                     (15, 5, 0, 0, 1), (16, 1, -1, 0, 0), (16, 1, -1, 1, 0), (16, 1, 0, 0, 0),
                     (16, 1, 0, 1, 0), (16, 2, 0, 0, 0), (16, 2, 0, 1, 0), (16, 2, 0, 0, 1),
                     (16, 2, 1, 1, 0), (16, 3, 0, 0, 0), (16, 3, 0, 1, 0), (16, 3, 1, 1, 0),
                     (16, 3, 1, 0, 0), (16, 4, 0, 0, 0), (16, 4, 0, 0, 1), (16, 4, 0, 1, 0),
                     (16, 4, 1, 0, 0), (16, 4, 1, 1, 0)]
# h5 columns that duplicate the topology and are therefore not stored
REDUNDANT_COLUMNS = ("level", "atom_type", "con_atom", "con_type", "max_attention")


def loadBranch(folder, b):
  with gzip.open(os.path.join(folder, f"{b}.gz"), "rb") as f:
    return pickle.load(f)


def looksLikeHydrogenRoot(tree):
  """Does every child of this branch's root attach to nothing (con_atom = -1)?

  A hydrogen's descent starts at its heavy neighbour, and the node for that
  step has no attachment position; no other root has such children.
  """
  children = [c for c in tree[0][5] if c != 0]
  return bool(children) and all(int(tree[c][2]) == -1 for c in children)


def main(argv=None):
  ap = argparse.ArgumentParser(prog="python -m rdkit.Chem.DASHTree.convert", description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument("tree_folder", help="the legacy tree: the folder holding 0.gz, 0.h5, ...")
  ap.add_argument("out", help="the .dash file to write")
  ap.add_argument("--props", default=None,
                  help="comma separated property columns to keep (default: all of them)")
  args = ap.parse_args(argv)

  try:
    import pandas as pd
    import tables  # noqa: F401, pandas reads HDF5 through it
  except ImportError as e:
    sys.exit(f"converting a legacy tree needs pandas and pytables: {e}")

  folder = args.tree_folder
  if not os.path.isdir(folder):
    sys.exit(f"{folder} is not a directory")
  gz = sorted(int(f[:-3]) for f in os.listdir(folder) if f.endswith(".gz") and f[:-3].isdigit())
  nBranches = len(gz)
  if nBranches == 0:
    sys.exit(f"no <n>.gz files in {folder}: a legacy tree is 0.gz, 0.h5, 1.gz, 1.h5, ... "
             f"directly inside the folder you pass")
  if gz != list(range(nBranches)):
    sys.exit(f"the .gz files in {folder} are not numbered 0..{nBranches - 1} without gaps")
  missingH5 = [b for b in gz if not os.path.exists(os.path.join(folder, f"{b}.h5"))]
  if missingH5:
    sys.exit(f"no .h5 data file for branches {missingH5[:10]}")

  # ---- pass 1: sizes, the hydrogen branch and the property columns ---------
  branchN = np.zeros(nBranches, dtype=np.int64)
  hCandidates = []
  for b in range(nBranches):
    tree = loadBranch(folder, b)
    branchN[b] = len(tree)
    if looksLikeHydrogenRoot(tree):
      hCandidates.append(b)
    del tree
  nNodes = int(branchN.sum())
  if nNodes > MAX_NODES:
    sys.exit(f"{nNodes} nodes exceeds the container's 32-bit node ids")
  print(f"{nBranches} branches, {nNodes:,} nodes")

  if len(STANDARD_FEATURES) != nBranches:
    sys.exit(f"the feature list has {len(STANDARD_FEATURES)} classes but the tree has "
             f"{nBranches} branches")
  hBranch = STANDARD_FEATURES.index((1, 1, 0, 0, 0))
  if hCandidates != [hBranch]:
    sys.exit(f"the data puts hydrogen at branches {hCandidates}: the tree was not built on "
             f"the 122 classes of the published trees")
  featureTable = np.array([(*f, (0, 0, 0)) for f in STANDARD_FEATURES], dtype=FEATURE_DTYPE)
  branchRoot = np.zeros(nBranches, dtype=np.uint32)
  branchRoot[1:] = np.cumsum(branchN)[:-1]

  df0 = pd.read_hdf(os.path.join(folder, "0.h5"), key="df", mode="r")
  if args.props:
    props = [p.strip() for p in args.props.split(",") if p.strip()]
    missing = [p for p in props if p not in df0.columns]
    if missing:
      sys.exit(f"requested columns not in the data: {missing}; available: {list(df0.columns)}")
    redundant = [p for p in props if p in REDUNDANT_COLUMNS]
    if redundant:
      sys.exit(f"these columns duplicate the topology and are not stored: {redundant}")
  else:
    props = [c for c in df0.columns if c not in REDUNDANT_COLUMNS]

  def outDtype(col):
    """Each column keeps the width its source needs: the MBIS charges are
    float16, the DASH-properties columns genuine float64."""
    d = df0[col].dtype
    if np.issubdtype(d, np.integer):
      return DTYPE_I32
    if d == np.float32:
      return DTYPE_F32
    return DTYPE_F16 if d == np.float16 else DTYPE_F64

  propDtypes = [outDtype(c) for c in props]
  print("properties: " + ", ".join(f"{c} ({NP_DTYPE[dt].name})"
                                   for c, dt in zip(props, propDtypes)))

  # ---- pass 2: renumber breadth-first, one branch at a time --------------
  # Processing nodes in id order and appending their children makes each
  # node's children a contiguous run, so no child-index array is needed, and
  # gathers the shallow levels every descent crosses at the start of the branch.
  records = np.zeros(nNodes, dtype=RECORD_DTYPE)
  nodeAttn = np.zeros(nNodes, dtype=np.float32)
  nodeSource = np.zeros(nNodes, dtype=np.uint32)
  propArrays = [np.empty(nNodes, dtype=NP_DTYPE[dt]) for dt in propDtypes]

  for b in range(nBranches):
    tree = loadBranch(folder, b)
    base = int(branchRoot[b])
    count = len(tree)
    cum = np.zeros(count, dtype=np.float64)  # a root's own attention never counts
    root = tree[0]
    records["key"][base] = packKey(int(root[1]), int(root[2]), int(root[3]))
    nodeAttn[base] = np.float32(root[4])

    nextLocal = 1
    for local in range(count):
      src = int(nodeSource[base + local])
      # Every root lists itself as its first child. That edge can never match
      # (its key attaches to nothing, which only a hydrogen's heavy neighbour
      # does, and that is never the hydrogen class), so it is dropped.
      children = [c for c in tree[src][5] if c != src]
      if tree[src][5].count(src) > (src == 0):
        sys.exit(f"branch {b}: node {src} has an unexpected reference to itself")
      if len(children) > MAX_CHILDREN:
        sys.exit(f"node {src} of branch {b} has {len(children)} children, "
                 f"more than the container's {MAX_CHILDREN}")
      records["nChildren"][base + local] = len(children)
      records["firstChild"][base + local] = base + nextLocal

      # the first step out of the hydrogen root only walks to the heavy
      # neighbour, and the descent does not accumulate its attention
      accumulate = not (src == 0 and b == hBranch)
      for c in children:
        cn = tree[c]
        at, ca, ct = int(cn[1]), int(cn[2]), int(cn[3])
        if not 0 <= at < nBranches:
          sys.exit(f"atom_type {at} is not one of the {nBranches} classes")
        if not -1 <= ca <= MAX_CON_ATOM or not -1 <= ct <= MAX_CON_TYPE:
          sys.exit(f"con_atom {ca} or con_type {ct} does not fit the match key")
        attention = np.float32(cn[4])
        records["key"][base + nextLocal] = packKey(at, ca, ct)
        nodeSource[base + nextLocal] = c
        nodeAttn[base + nextLocal] = attention
        cum[nextLocal] = cum[local] + (float(attention) if accumulate else 0.0)
        nextLocal += 1

    if nextLocal != count:
      sys.exit(f"branch {b}: only {nextLocal} of its {count} nodes are reachable from the root")
    records["flags"][base:base + count] = np.where(cum > DEFAULT_ATTENTION_THRESHOLD,
                                                   NODE_FLAG_STOP, 0)
    del tree, cum

    # the property columns, indexed by the new numbering
    d = pd.read_hdf(os.path.join(folder, f"{b}.h5"), key="df", mode="r")
    if len(d) != count:
      sys.exit(f"branch {b}: {count} nodes but {len(d)} data rows")
    rows = nodeSource[base:base + count]
    for arr, col, dt in zip(propArrays, props, propDtypes):
      arr[base:base + count] = d[col].to_numpy()[rows].astype(NP_DTYPE[dt], copy=False)
    del d

  try:
    writeContainer(args.out, branchRoot, featureTable, records, nodeAttn,
                   list(zip(props, propDtypes, propArrays)), DEFAULT_ATTENTION_THRESHOLD)
  except ValueError as e:
    sys.exit(str(e))
  print(f"wrote {args.out}: {os.path.getsize(args.out) / 1e6:.1f} MB")


if __name__ == "__main__":
  main()
