#!/usr/bin/env python
#  Copyright (C) 2026 Greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Build the histogram table used by the NormalizedDescriptors extension.

This ports descriptastorus' RDKit2DHistogramNormalized
The raw descriptor values over a reference set of molecules are
binned with numpy.histogram (min(1000, number of distinct values)
equal-width bins) exactly as descriptastorus'
data/d_descriptors/make_histdists.py does, giving (left edge, cumulative
fraction) pairs, and a value v is looked up the way descriptastorus does it:

    p = bisect.bisect(bins, (v,))       # number of edges < v
    normalized(v) = bins[p][1] if p < len(bins) else 1.0

see https://github.com/bp-kelley/descriptastorus.git

Table (relative to External/NormalizedDescriptors):

  data/normalized_descriptor_histograms.txt

Subcommands:

  import-hists  seed or update the table from a descriptastorus hists.py
  compute       compute descriptor values for a SMILES file (descriptastorus'
                data/chembl_100k.smi) and write one d_<name>.gz sample file
                per descriptor, in the same format as descriptastorus'
                data/d_descriptors
  histogram     bin each sample file and add it to the table
  missing       list RDKit descriptors that have no histogram yet
  check         compare table lookups against a descriptastorus hists.py

Adding new descriptors (e.g. everything RDKit has that isn't in the table):

  python make_normalized_histograms.py compute \\
      --smiles /path/to/descriptastorus/data/chembl_100k.smi \\
      --out samples $(python make_normalized_histograms.py missing)
  python make_normalized_histograms.py histogram samples/d_*.gz

Rebuilding every histogram with the current RDKit (use --replace):

  python make_normalized_histograms.py compute --smiles ... --out samples \\
      $(python make_normalized_histograms.py names)
  python make_normalized_histograms.py histogram --replace samples/d_*.gz

Requires numpy; and the RDKit python wrappers
"""

import argparse
import ast
import bisect
import gzip
import math
import multiprocessing
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(HERE, "..", "data"))
DEFAULT_HISTS = os.path.join(DATA_DIR, "normalized_descriptor_histograms.txt")

TABLE_FORMAT_VERSION = 1

# descriptastorus make_histdists.py uses at most this many bins
MAX_HIST_BINS = 1000


def log(*args):
  print(*args, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# samples (descriptastorus data/d_descriptors format: gzipped python list)


def read_sample(path):
  with gzip.open(path, "rt") as f:
    txt = f.read()
  try:
    vals = ast.literal_eval(txt)
  except ValueError:
    # older files can contain bare inf / nan
    vals = eval(txt, {"__builtins__": {}}, {"inf": math.inf, "nan": math.nan})
  vals = np.asarray(vals, dtype=float)
  return vals[np.isfinite(vals)]


def sample_name(path):
  base = os.path.basename(path)
  if base.endswith(".gz"):
    base = base[:-3]
  if base.startswith("d_"):
    base = base[2:]
  return base


def write_sample(path, values):
  with gzip.open(path, "wt") as f:
    f.write(repr([float(v) for v in values]))


# ---------------------------------------------------------------------------
# compute


def rdkit_descriptor_functions():
  from rdkit.Chem import Descriptors
  return dict(Descriptors._descList)


_worker_names = None


def _compute_init(names):
  global _worker_names
  _worker_names = names


def _compute_chunk(smiles):
  from rdkit import Chem, RDLogger
  RDLogger.DisableLog("rdApp.*")
  funcs = rdkit_descriptor_functions()
  res = {n: [] for n in _worker_names}
  nfailed = 0
  for smi in smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
      nfailed += 1
      continue
    for n in _worker_names:
      try:
        res[n].append(float(funcs[n](mol)))
      except Exception:
        pass
  return res, nfailed


def read_smiles(path):
  smiles = []
  with open(path) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith("#"):
        continue
      smi = line.split()[0]
      if smi.lower() == "smiles":
        continue
      smiles.append(smi)
  return smiles


def cmd_compute(args):
  funcs = rdkit_descriptor_functions()
  unknown = [n for n in args.descriptors if n not in funcs]
  if unknown:
    log(f"Skipping names RDKit doesn't compute: {' '.join(unknown)}")
    args.descriptors = [n for n in args.descriptors if n in funcs]
  smiles = read_smiles(args.smiles)
  log(f"Computing {len(args.descriptors)} descriptors for {len(smiles)} "
      f"molecules from {args.smiles}")
  chunks = [smiles[i:i + 500] for i in range(0, len(smiles), 500)]
  values = {n: [] for n in args.descriptors}
  nfailed = 0
  with multiprocessing.Pool(args.nprocs or None, _compute_init,
                            (args.descriptors, )) as pool:
    for i, (res, nf) in enumerate(pool.imap(_compute_chunk, chunks)):
      if i % 20 == 19:
        log(f"  {min((i + 1) * 500, len(smiles))} molecules done")
      nfailed += nf
      for n, vs in res.items():
        values[n].extend(vs)
  if nfailed:
    log(f"{nfailed} SMILES failed to parse and were skipped")
  os.makedirs(args.out, exist_ok=True)
  for n, vs in values.items():
    path = os.path.join(args.out, f"d_{n}.gz")
    write_sample(path, vs)
    log(f"Wrote {len(vs)} values to {path}")


def cmd_missing(args):
  have = read_hist_table(args.table)
  for n in rdkit_descriptor_functions():
    if n not in have:
      print(n)


def cmd_names(args):
  for n in sorted(read_hist_table(args.table)):
    print(n)


# ---------------------------------------------------------------------------
# histograms


def make_histogram(vals):
  """descriptastorus data/d_descriptors/make_histdists.py, unchanged:
  (left bin edge, cumulative fraction of values up to and including that
  bin) for min(1000, distinct values) equal-width bins."""
  vals = np.asarray(vals, dtype=float)
  vals = vals[np.isfinite(vals)]
  n = min(MAX_HIST_BINS, len(set(vals)))
  hist, xaxis = np.histogram(vals, bins=n)
  total = hist.sum()
  bins = []
  last = 0.0
  for value, x in zip(hist, xaxis):
    last += value
    bins.append((float(x), float(last / total)))
  return bins


def hist_lookup(bins, v):
  """descriptastorus rdNormalizedDescriptors.histcdf, unchanged."""
  p = bisect.bisect(bins, (v, ))
  if p < len(bins):
    return bins[p][1]
  return 1.0


def read_hist_table(path):
  hists = {}
  if not os.path.exists(path):
    return hists
  with open(path) as f:
    it = (l for l in f if l.strip() and not l.startswith("#"))
    for line in it:
      tag, name, n = line.split()
      assert tag == "histogram", line
      bins = []
      for _ in range(int(n)):
        x, c = next(it).split()
        bins.append((float(x), float(c)))
      hists[name] = bins
  return hists


def write_hist_table(hists, path):
  lines = [
    "# NormalizedDescriptors histogram table, format version "
    f"{TABLE_FORMAT_VERSION}",
    "# Generated by tools/make_normalized_histograms.py; do not edit.",
    "#",
    "# Raw-data histograms (no fitted distribution), as descriptastorus'",
    "# RDKit2DHistogramNormalized. For each descriptor, nbins lines of",
    "# <left edge> <cumulative fraction>, edges ascending. normalized(v) is",
    "# the fraction of the first bin whose edge is >= v, or 1.0 when v is",
    "# greater than every edge (python: bins[bisect(bins, (v,))][1]).",
    "# No interpolation and no clipping.",
    "#",
    "# histogram <name> <nbins>",
    "# <edge> <cumulative fraction>   (nbins lines)",
  ]
  for name in sorted(hists):
    bins = hists[name]
    lines.append(f"histogram {name} {len(bins)}")
    lines.extend(f"{x!r} {c!r}" for x, c in bins)
  with open(path, "w") as f:
    f.write("\n".join(lines))
    f.write("\n")
  log(f"Wrote {len(hists)} histograms to {path}")


def read_hists_py(path):
  """Read a descriptastorus hists.py without importing descriptastorus."""
  ns = {"inf": math.inf, "nan": math.nan}
  with open(path) as f:
    exec(f.read(), ns)
  return ns["hists"]


def cmd_import_hists(args):
  hists = read_hist_table(args.table)
  src = read_hists_py(args.hists_py)
  n = 0
  for name, bins in src.items():
    if name in hists and not args.replace:
      continue
    hists[name] = [(float(x), float(c)) for x, c in bins]
    n += 1
  log(f"Imported {n} of {len(src)} histograms from {args.hists_py}")
  write_hist_table(hists, args.table)


def cmd_histogram(args):
  hists = read_hist_table(args.table)
  for path in args.samples:
    name = sample_name(path)
    if name in hists and not args.replace:
      log(f"Skipping {name}: already has a histogram (use --replace)")
      continue
    vals = read_sample(path)
    if len(vals) == 0:
      log(f"Skipping {name}: no finite values in {path}")
      continue
    hists[name] = make_histogram(vals)
    log(f"{name}: {len(hists[name])} bins from {len(vals)} values")
  write_hist_table(hists, args.table)


def cmd_check_hists(args):
  table = read_hist_table(args.table)
  ref = read_hists_py(args.hists_py)
  rng = np.random.default_rng(0)
  bad = 0
  for name, bins in ref.items():
    if name not in table:
      log(f"{name}: in {args.hists_py} but not the table")
      bad += 1
      continue
    edges = np.array([b[0] for b in bins])
    probes = np.concatenate([
      edges, edges + 1e-9,
      rng.uniform(edges[0] - 1, edges[-1] + 1, 2000)
    ])
    diff = max(
      abs(hist_lookup(table[name], v) - hist_lookup(bins, v)) for v in probes)
    if diff:
      log(f"{name}: lookups differ by up to {diff}")
      bad += 1
  log(f"Checked {len(ref)} histograms, {bad} differ")
  if bad:
    sys.exit(1)


# ---------------------------------------------------------------------------


def main(argv=None):
  p = argparse.ArgumentParser(
    description=__doc__.split("\n\n")[0],
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__.split("\n\n", 1)[1])
  p.add_argument("--table", default=DEFAULT_HISTS,
                 help="histogram table (default: %(default)s)")
  sub = p.add_subparsers(dest="cmd", required=True)

  s = sub.add_parser("import-hists", help="import a descriptastorus hists.py")
  s.add_argument("hists_py")
  s.add_argument("--replace", action="store_true",
                 help="overwrite descriptors already in the table")
  s.set_defaults(func=cmd_import_hists)

  s = sub.add_parser("compute", help="compute descriptor samples")
  s.add_argument("descriptors", nargs="+", help="RDKit descriptor names")
  s.add_argument("--smiles", required=True,
                 help="SMILES file, e.g. descriptastorus/data/chembl_100k.smi")
  s.add_argument("--out", default="samples",
                 help="directory for d_<name>.gz files (default: %(default)s)")
  s.add_argument("--nprocs", type=int, default=0)
  s.set_defaults(func=cmd_compute)

  s = sub.add_parser("histogram", help="add histograms of sample files")
  s.add_argument("samples", nargs="+", help="d_<name>.gz sample files")
  s.add_argument("--replace", action="store_true",
                 help="rebuild histograms already in the table")
  s.set_defaults(func=cmd_histogram)

  s = sub.add_parser("missing",
                     help="list RDKit descriptors without a histogram")
  s.set_defaults(func=cmd_missing)

  s = sub.add_parser("names", help="list descriptors in the table")
  s.set_defaults(func=cmd_names)

  s = sub.add_parser("check",
                     help="check table lookups against a descriptastorus "
                     "hists.py")
  s.add_argument("hists_py")
  s.set_defaults(func=cmd_check_hists)

  args = p.parse_args(argv)
  args.func(args)


if __name__ == "__main__":
  main()
