#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Prune a .dash container into a smaller one.

A published DASH tree has 9.4 million nodes and is 150 MB even with nothing but
the charges in it. That is more than a user who wants a quick number needs.
This tool writes a container with a subset of the nodes, renumbered
breadth-first again so it reads like any other.

The tool drops a node together with its whole subtree. An atom whose descent
reaches the cut reads the value of the deepest node that is still there, as it
would at the end of a full path. So a pruned tree never fails on a molecule the
full tree accepts. It answers with a shallower, more general substructure.

For that to hold the descent has to stop where the cut is. It takes the first
child in tree order that matches any candidate. If that child were gone, a
later sibling might match instead and the path would continue in another
subtree. A removed child with a kept sibling after it therefore stays as a
value-less leaf. The descent still takes it and stops there, and the property
walk reads the parent.

Two rules, combined with "and", choose what stays:

  --tolerance EPS   drop a subtree when every value in it lies within EPS of the
                    parent's value (column --value-property, default 'result').
                    On the default MBIS tree, 0.05 e keeps 6.5 % of the nodes
                    and moves charges by 0.014 e RMS.
  --max-level L     keep levels 0..L only

Roots always stay, and so do the children of every Z = 1 root. The descent
matches a hydrogen through its heavy neighbour, and the C++ refuses a hydrogen
whose heavy-neighbour node is missing. The tool copies the stop flags unchanged,
since every ancestor of a kept node stays and the cumulative attention along
its path is the same.
"""

import argparse
import os
import sys

import numpy as np

from .container import MISSING_INT32, Container, writeContainer


def toleranceMask(parent, level, values, eps):
  """True where a node changes the value read by at least eps.

  err(n) = max |value(x) - value(parent(n))| over n and its subtree. A node
  whose parent carries no value stays if its subtree has a value, since the
  path above it may have no value at all. A subtree without any value in this
  column changes nothing that this column reads, so it goes. A populated tree
  has such subtrees where no atom of the data went.
  """
  vals = values.astype(np.float64)
  if not np.issubdtype(values.dtype, np.floating):
    vals[values == MISSING_INT32] = np.nan  # an int32 column's "no value"
  smin = np.where(np.isnan(vals), np.inf, vals)
  smax = np.where(np.isnan(vals), -np.inf, vals)
  for depth in range(int(level.max()), 0, -1):
    nodes = np.nonzero(level == depth)[0]
    np.minimum.at(smin, parent[nodes], smin[nodes])
    np.maximum.at(smax, parent[nodes], smax[nodes])
  parentVal = np.where(parent >= 0, vals[np.maximum(parent, 0)], np.nan)
  err = np.maximum(np.abs(smax - parentVal), np.abs(smin - parentVal))
  err[np.isnan(parentVal)] = np.inf
  err[smin > smax] = 0.0  # only a subtree without values keeps the start values
  return err >= eps


def renumber(keep, parent, roots):
  """Old ids in new order, and per new node its child count and first child.

  The kept nodes are closed under parents, and the source keeps children
  contiguous after their parent. So the old order is already a valid,
  breadth-first numbering.
  """
  old = np.nonzero(keep)[0]
  newOfOld = np.full(keep.size, -1, dtype=np.int64)
  newOfOld[old] = np.arange(old.size)
  kids = old[parent[old] >= 0]
  parentNew = newOfOld[parent[kids]]
  nChildren = np.bincount(parentNew, minlength=old.size)
  firstChild = np.zeros(old.size, dtype=np.int64)
  withKids, first = np.unique(parentNew, return_index=True)
  firstChild[withKids] = newOfOld[kids[first]]
  return old, nChildren, firstChild, newOfOld[roots]


def main(argv=None):
  ap = argparse.ArgumentParser(prog="python -m rdkit.Chem.DASHTree.prune", description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument("inp", help="the .dash file to prune")
  ap.add_argument("out", help="the .dash file to write")
  ap.add_argument("--tolerance", type=float, default=None,
                  help="drop subtrees whose values all lie within this of the parent's")
  ap.add_argument("--value-property", default="result",
                  help="the column that --tolerance compares (default: result)")
  ap.add_argument("--max-level", type=int, default=None, help="keep levels 0 up to this one")
  ap.add_argument("--props", default=None,
                  help="comma-separated property columns to keep (default: all)")
  args = ap.parse_args(argv)
  if (args.tolerance, args.max_level) == (None, None):
    ap.error("give --tolerance, --max-level or both")
  if os.path.exists(args.out) and os.path.samefile(args.inp, args.out):
    ap.error("the output would overwrite the input")

  src = Container(args.inp)
  names = list(src.props)
  if args.props:
    names = [p.strip() for p in args.props.split(",") if p.strip()]
    missing = [p for p in names if p not in src.props]
    if missing:
      sys.exit(f"requested columns not in the file: {missing} (available: {list(src.props)})")
    if len(set(names)) != len(names):
      sys.exit(f"--props names a column twice: {names}")
  print(f"{args.inp}: {src.nBranches} branches, {src.nNodes:,} nodes, "
        f"properties {list(src.props)}")
  parent, level = src.parentAndLevel()
  print(f"  levels 0..{int(level.max())}")

  # ---- a node stays if it passes all the rules -----------------------------
  passes = np.ones(src.nNodes, dtype=bool)
  if args.tolerance is not None:
    if args.value_property not in src.props:
      sys.exit(f"no property '{args.value_property}' in the file (available: {list(src.props)})")
    changed = toleranceMask(parent, level, src.props[args.value_property][1], args.tolerance)
    print(f"  tolerance {args.tolerance} on '{args.value_property}': "
          f"{int((~changed).sum()):,} nodes change the value by less")
    passes &= changed
  if args.max_level is not None:
    passes &= level <= args.max_level
    print(f"  max level {args.max_level}: {int((level > args.max_level).sum()):,} nodes deeper")

  roots = src.branchRoot.astype(np.int64)
  # the reader matches every Z = 1 atom through its heavy neighbour, so those nodes stay
  passes[np.isin(parent, roots[src.features["atomicNum"] == 1])] = True

  # closure: a node stays only if every ancestor stays, so visit the shallow levels first
  keep = np.zeros(src.nNodes, dtype=bool)
  keep[roots] = True
  for depth in range(1, int(level.max()) + 1):
    nodes = np.nonzero(level == depth)[0]
    keep[nodes] = keep[parent[nodes]] & passes[nodes]
  nReal = int(keep.sum())

  # ---- placeholders: a removed child with a kept sibling after it stays as a
  # value-less leaf, so the descent still takes it and stops there ----------
  kids = np.nonzero(parent >= 0)[0]
  kids = kids[keep[parent[kids]]]
  keptKid = keep[kids]
  lastKept = np.full(src.nNodes, -1, dtype=np.int64)
  np.maximum.at(lastKept, parent[kids[keptKid]], kids[keptKid])
  placeholder = np.zeros(src.nNodes, dtype=bool)
  placeholder[kids[~keptKid & (kids < lastKept[parent[kids]])]] = True
  keep |= placeholder
  print(f"  keep {nReal:,} of {src.nNodes:,} nodes "
        f"({100.0 * nReal / src.nNodes:.2f} %) plus {int(placeholder.sum()):,} "
        f"placeholders")

  # ---- renumber and gather ----------------------------------------------
  old, nChildren, firstChild, newRoots = renumber(keep, parent, roots)
  records = src.records[old]
  records["nChildren"] = nChildren
  records["firstChild"] = firstChild
  attn = src.attn[old]
  blank = placeholder[old]
  props = []
  for name in names:
    tag, arr = src.props[name]
    arr = arr[old]
    arr[blank] = np.nan if np.issubdtype(arr.dtype, np.floating) else MISSING_INT32
    props.append((name, tag, arr))

  writeContainer(args.out, newRoots, src.features, records, attn, props, src.threshold)
  print(f"wrote {args.out}: {os.path.getsize(args.out) / 1e6:.1f} MB")


if __name__ == "__main__":
  main()
