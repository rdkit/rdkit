#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Fill a DASH tree with new atomic properties.

The tool keeps the nodes of a tree and adds property columns to it. It matches
each atom of each molecule, as GetMolProperty does, and gives the atom's value
to every node on the atom's path. For each property P, each node then holds:

  P          the mean of the values it got
  P_std      their sample standard deviation
  P_count    the number of values

A node without values has no value in the three columns. A node with one value
has no value in P_std. For each column, the reader uses the deepest node on the
atom's path that has a value, as it does at the cut of a pruned tree. So P_count
gives the number of values in the mean that the reader returns. P_std can come
from a node above that node, which has more values.

The molecules come from SDF files and must have explicit hydrogens. Each record
gives the values of a property as an RDKit atom property list. This is a data
field 'atom.dprop.P' with one value per atom. Spaces separate the values, and
'n/a' marks an atom without a value. Chem.CreateAtomDoublePropertyList writes
such a field.

  python -m rdkit.Chem.DASHTree.populate tree.dash out.dash mols.sdf --props P,Q

The descent decides which nodes an atom reaches. Populate a tree with the
DASHParams that you read it with. For a tree that you read like the DASH-tree
python package, use --legacy-tie-breaking.

From python, collectStatistics() takes any iterable of (mol, {P: values})
pairs, and writeTree() writes the result.
"""

import argparse
import os
import sys

import numpy as np

from rdkit import Chem
from rdkit.Chem import rdDASHTree

from .container import (DTYPE_F16, DTYPE_F32, DTYPE_I32, MISSING_INT32, PROP_NAME_LEN, Container,
                        writeContainer)

# the number of (node, value) pairs to collect before they go into the sums
_BATCH = 1 << 24


def columnNames(names):
  """The three column names of each property. Raises ValueError if a name
  repeats or is too long for the container."""
  columns = [c for p in names for c in (p, f"{p}_std", f"{p}_count")]
  if len(set(columns)) != len(columns):
    raise ValueError(f"column names repeat: {columns}")
  tooLong = [c for c in columns if len(c.encode("utf-8")) >= PROP_NAME_LEN]
  if tooLong:
    raise ValueError(f"column names are longer than {PROP_NAME_LEN - 1} bytes: {tooLong}")
  return columns


def collectStatistics(tree, records, names, params=None):
  """The mean, sample standard deviation and count of each property per node.

    ARGUMENTS:
      - tree: path of the '.dash' file whose nodes collect the values
      - records: an iterable of (mol, values) pairs. values maps a property
        name to one value per atom. Use NaN for an atom without a value. Each
        entry of values can also be a 2D array with one row per conformer.
        Each row then adds one value per atom on the same paths. A record can
        leave out properties.
      - names: the properties to collect
      - params: the rdDASHTree.DASHParams of the descent

    RETURNS: a list of (column name, dtype tag, array) for writeTree(), and the
      number of molecules with an atom outside the tree's classes. These
      molecules add nothing.

    Raises ValueError for a molecule without explicit hydrogens or with a stale
    property cache, and for values that do not fit the molecule.
    """
  columnNames(names)
  src = Container(tree)
  reader = rdDASHTree.DASHTree(str(tree), [])
  params = params or rdDASHTree.DASHParams()
  roots = src.branchRoot.astype(np.int64)
  nNodes = src.nNodes
  count = {p: np.zeros(nNodes, dtype=np.int64) for p in names}
  total = {p: np.zeros(nNodes) for p in names}
  squares = {p: np.zeros(nNodes) for p in names}
  pending = {p: [] for p in names}
  nPending = 0
  failed = 0

  def addPending():
    for p in names:
      if not pending[p]:
        continue
      nodes = np.concatenate([n for n, _ in pending[p]])
      values = np.concatenate([v for _, v in pending[p]])
      pending[p] = []
      has = ~np.isnan(values)
      nodes, values = nodes[has], values[has]
      count[p] += np.bincount(nodes, minlength=nNodes)
      total[p] += np.bincount(nodes, values, minlength=nNodes)
      squares[p] += np.bincount(nodes, values * values, minlength=nNodes)

  for mol, values in records:
    unknown = set(values) - set(names)
    if unknown:
      raise ValueError(f"values for properties that are not collected: {sorted(unknown)}")
    if mol.GetNumAtoms() == 0:
      continue
    if mol.NeedsUpdatePropertyCache() or any(a.GetTotalNumHs() for a in mol.GetAtoms()):
      raise ValueError("DASH needs sanitized molecules with explicit hydrogens (Chem.AddHs)")
    try:
      paths = [reader.GetAtomNodePath(mol, i, params) for i in range(mol.GetNumAtoms())]
    except ValueError:
      failed += 1  # an atom is in none of the tree's classes
      continue
    nodes = np.concatenate(
      [roots[path[0]] + np.asarray(path[1:], dtype=np.int64) for path in paths])
    atomOfNode = np.repeat(np.arange(len(paths)), [len(path) - 1 for path in paths])
    for p, perAtom in values.items():
      perAtom = np.atleast_2d(np.asarray(perAtom, dtype=np.float64))
      if perAtom.ndim != 2 or perAtom.shape[1] != mol.GetNumAtoms():
        raise ValueError(f"property {p} has shape {perAtom.shape} for {mol.GetNumAtoms()} atoms")
      if np.isinf(perAtom).any():
        raise ValueError(f"property {p} has an infinite value")
      pending[p].append((np.tile(nodes, perAtom.shape[0]), perAtom[:, atomOfNode].ravel()))
      nPending += nodes.size * perAtom.shape[0]
    if nPending > _BATCH:
      addPending()
      nPending = 0
  addPending()

  columns = []
  for p in names:
    n = count[p]
    with np.errstate(invalid="ignore", divide="ignore"):
      mean = total[p] / n
      variance = (squares[p] - total[p] * mean) / (n - 1)
    std = np.sqrt(np.maximum(variance, 0.0))
    std[n < 2] = np.nan
    if np.nanmax(std, initial=0.0) > np.finfo(np.float16).max:
      raise ValueError(f"a standard deviation of {p} is too large for a float16 column")
    columns += [(p, DTYPE_F32, mean), (f"{p}_std", DTYPE_F16, std),
                (f"{p}_count", DTYPE_I32, np.where(n > 0, n, MISSING_INT32))]
  return columns, failed


def writeTree(tree, out, columns, keepExisting=True):
  """Writes the nodes of tree with the new columns to out.

    A new column replaces an existing column of the same name. With
    keepExisting False, out has the new columns only.
    """
  if os.path.exists(out) and os.path.samefile(tree, out):
    raise ValueError("the output would overwrite the input")
  src = Container(tree)
  new = {name for name, _, _ in columns}
  kept = [(name, tag, arr) for name, (tag, arr) in src.props.items()
          if keepExisting and name not in new]
  writeContainer(out, src.branchRoot, src.features, src.records, src.attn, kept + columns,
                 src.threshold)


def _readSdf(paths, names, bad):
  """(mol, values) for each record of the SDF files. bad[0] counts the records
  that do not parse."""
  for path in paths:
    for mol in Chem.SDMolSupplier(path, removeHs=False):
      if mol is None:
        bad[0] += 1
        continue
      atoms = list(mol.GetAtoms())
      yield mol, {p: [a.GetDoubleProp(p) if a.HasProp(p) else np.nan for a in atoms] for p in names}


def main(argv=None):
  ap = argparse.ArgumentParser(prog="python -m rdkit.Chem.DASHTree.populate", description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument("inp", help="the .dash file whose nodes collect the values")
  ap.add_argument("out", help="the .dash file to write")
  ap.add_argument("sdf", nargs="+", help="SDF files with atom property lists")
  ap.add_argument("--props", required=True, help="comma-separated property names")
  ap.add_argument("--legacy-tie-breaking", action="store_true",
                  help="match atoms as the DASH-tree python package does")
  ap.add_argument("--drop-existing", action="store_true",
                  help="write the new columns only, not the columns of the input")
  args = ap.parse_args(argv)
  missing = [p for p in [args.inp] + args.sdf if not os.path.isfile(p)]
  if missing:
    ap.error(f"no such file: {missing}")
  if os.path.exists(args.out) and os.path.samefile(args.inp, args.out):
    ap.error("the output would overwrite the input")
  names = [p.strip() for p in args.props.split(",") if p.strip()]
  if not names:
    ap.error("--props names no property")
  try:
    columnNames(names)
  except ValueError as err:
    ap.error(str(err))

  params = rdDASHTree.DASHParams()
  params.legacyTieBreaking = args.legacy_tie_breaking
  bad = [0]
  columns, failed = collectStatistics(args.inp, _readSdf(args.sdf, names, bad), names, params)
  for name, _, arr in columns:
    if name in names:
      print(f"  {name}: {int(np.count_nonzero(~np.isnan(arr))):,} nodes have a value")
  if bad[0] or failed:
    print(
      f"  skipped {bad[0]} records that do not parse and {failed} molecules with an atom "
      f"outside the tree's classes", file=sys.stderr)
  writeTree(args.inp, args.out, columns, keepExisting=not args.drop_existing)
  print(f"wrote {args.out}: {os.path.getsize(args.out) / 1e6:.1f} MB")


if __name__ == "__main__":
  main()
