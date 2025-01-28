# Copyright (C) Brian Kelley 2021
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""
Free Wilson Analysis

A Mathematical Contribution to Structurre-Activity Studies
Spencer M, Free, Jr and James W Wilson
Journal of Medicinal Chemistry
Vol 7, Number 4, July 6, 1964

Basic usage: get a scaffold (or scaffolds) some compounds and their scores, 
then run freewilson:


Read some example data:
```
>>> import os, sys
>>> DATA_PATH = "data"
>>> smilesfile = os.path.join(DATA_PATH, "CHEMBL2321810.smi")
>>> scaffoldfile = os.path.join(DATA_PATH, "CHEMBL2321810_scaffold.mol")
>>> csvfile = os.path.join(DATA_PATH, "CHEMBL2321810_act.csv")
>>> mols = []
>>> for line in open(smilesfile):
...     smiles, name = line.strip().split()
...     m = Chem.MolFromSmiles(smiles)
...     m.SetProp("_Name", name)
...     mols.append(m)
>>> scaffold = Chem.MolFromMolBlock(open(scaffoldfile).read())
>>> data = {k:float(v) for k,v in list(csv.reader(open(csvfile)))[1:]}
>>> scores = [data[m.GetProp("_Name")] for m in mols]

```

And do the decomposition:
```
>>> from freewilson import FWDecompose, FWBuild, predictions_to_csv
>>> decomp = FWDecompose(scaffold, mols, scores)

```

Scores need to be in an appropriate form for regression analysis, i.e. pIC50s as opposed to IC50s.
Enumerations are a lot faster if you know a prediction cutoff or molecular weight cutoff.

Scaffolds can be any scaffold or smarts pattern or list of either.  The system automatically 
enumerates matching cores in the analysis.

To see if the decomposition is useful, check the r squared value

```
>>> print(f"Training R^2 is {decomp.r2:0.2f}")
Training R^2 is 0.81

```
Finally you can build the decomposition into new molecules:

```
>>> for pred in FWBuild(decomp):                            # doctest: +SKIP
...     print(pred.smiles, pred.prediction)

```

Now this builds both well and poorly predicted molecules.  To prevent
this you can use the following filters while building:

   1. pred_filter:  given a prediction, return True to keep the molecule
   2. hvy_filter: given a heavy atom count, return True to keep the molecule
   3. mw_filter: given a molecular weight, return True to keep the molecule
   4. mol_filter: given a sanitized molecule, return True to keep the molecule

Here are some examples (see using Molecular Filters below)

```
>>> preds = FWBuild(decomp, 
...                 pred_filter=lambda x: x > 8, 
...                 mw_filter=lambda mw: 100<mw<550)
>>> predictions_to_csv(sys.stdout, decomp, preds)

```

More Info
---------
You can also get some more information by setting logging to INFO

```
>>> import logging
>>> logging.getLogger().setLevel(logging.INFO)
>>> preds = list(FWBuild(decomp, pred_filter=lambda x: x > 8))

```

You'll see something like this:

```
INFO:root:Enumerating rgroups with no broken cycles...
INFO:root:	Core	1
INFO:root:	R1	73
INFO:root:	R10	2
INFO:root:	R3	445
100%|███████████████████████████████████████████████████████| 64970/64970 [00:05<00:00, 11247.38it/s]
INFO:root:Wrote 3 results out of 64970
	In Training set: 628
	Bad MW: 0
	Bad Pred: 64339
	Bad Filters: 0
	Bad smi: 0
	min mw: 380.429
	max mw: 772.4030000000001
	
	min pred: 2.8927324757327666
	max pred: 8.148182660251193
```

Using FMCS to find a scaffold
-----------------------------

If you don't have a good feel for the dataset, you can generate the scaffold by using
the rdFMCS package.  The following tries to find a decent MCS that covers 80% of the
input structures


```
>>> from rdkit.Chem import rdFMCS
>>> mcs = rdFMCS.FindMCS(mols[0:8], threshold=0.8, atomCompare=rdFMCS.AtomCompare.CompareAny,
...                      completeRingsOnly=True)                           # doctest: +SKIP
>>> decomp = FWDecompose(mcs.queryMol, mols, scores)                       # doctest: +SKIP

```
"""
import csv
import itertools
import logging
import math
import re
import sys
from collections import defaultdict, namedtuple
from typing import Generator, List

from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from tqdm import tqdm

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, molzip
from rdkit.Chem import rdRGroupDecomposition as rgd

logger = logging.getLogger("freewilson")

FreeWilsonPrediction = namedtuple("FreeWilsonPrediction", ['prediction', 'smiles', 'rgroups'])

# match dummy atoms in a smiles string to extract atom maps
dummypat = re.compile(r"\*:([0-9]+)")


# molzip doesn't handle some of the forms that the RGroupDecomposition
#  returns, this solves these issues.
def molzip_smi(smiles):
  """Fix a rgroup smiles for molzip, note that the core MUST come first
    in the smiles string, ala core.rgroup1.rgroup2 ...
    """
  dupes = set()
  sl = []
  for s in smiles.split("."):
    if s.count("*") >= 1:
      if s in dupes:
        continue
      else:
        dupes.add(s)
    sl.append(s)

  smiles = ".".join(sl)
  m = Chem.RWMol(Chem.MolFromSmiles(smiles, sanitize=False))
  frags = Chem.GetMolFrags(m)
  core = frags[0]
  atommaps = {}
  counts = defaultdict(int)
  for idx in core:
    atommap = m.GetAtomWithIdx(idx).GetAtomMapNum()
    if atommap:
      atommaps[atommap] = idx
      counts[atommap] += 1

  next_atommap = max(atommaps) + 1
  add_atommap = []
  for fragment in frags[1:]:
    for idx in fragment:
      atommap = m.GetAtomWithIdx(idx).GetAtomMapNum()
      if atommap:
        count = counts[atommap] = counts[atommap] + 1
        if count > 2:
          m.GetAtomWithIdx(idx).SetAtomMapNum(next_atommap)
          add_atommap.append((atommaps[atommap], next_atommap))
          next_atommap += 1

  for atomidx, atommap in add_atommap:
    atom = m.GetAtomWithIdx(atomidx)
    bonds = list(atom.GetBonds())
    if len(bonds) == 1:
      oatom = bonds[0].GetOtherAtom(atom)
      xatom = Chem.Atom(0)
      idx = m.AddAtom(xatom)
      xatom = m.GetAtomWithIdx(idx)
      xatom.SetAtomMapNum(atommap)
      m.AddBond(oatom.GetIdx(), xatom.GetIdx(), Chem.BondType.SINGLE)
  return Chem.molzip(m)


class RGroup:
  """FreeWilson RGroup
        smiles - isomeric smiles of the rgroup
        rgroup - Rgroup name ['Core', 'R1', 'R2'...]
        count - number of molecules with the rgroup
        coefficient - ridge regression coefficient
        idx - one-hot encoding for the rgroup
    """

  def __init__(self, smiles, rgroup, count, coefficient, idx=None):
    self.smiles = smiles  # smiles for the sidechain (n.b. can be a core as well)
    self.rgroup = rgroup  # rgroup Core, R1, R2,...
    self.count = count  # num molecules with this rgruop
    self.coefficient = coefficient  # ridge coefficient
    self.idx = idx  # descriptor index
    self.dummies = tuple([int(x) for x in sorted(dummypat.findall(smiles))])

    # Assemble some additive properties

    # will fail on some structures
    m = Chem.MolFromSmiles(smiles)
    if m:
      self.mw = Descriptors.MolWt(m)
      self.hvyct = Descriptors.HeavyAtomCount(m)
    else:
      # guess the MW if we can't sanitize
      table = Chem.GetPeriodicTable()
      m = Chem.MolFromSmiles(smiles, sanitize=False)
      if m:
        self.mw = 0.
        self.hvyct = 0
        for atom in m.GetAtoms():
          atomicnum = atom.GetAtomicNum()
          self.mw += table.GetAtomicWeight(atomicnum)
          if atomicnum > 1:
            self.hvyct += 1
        self.hvyct = Descriptors.HeavyAtomCount(m)
      else:
        self.mw = 0  # dunno
        self.hvyct = 0

  def __str__(self):
    return f"RGroup(smiles={repr(self.smiles)}, rgroup={repr(self.rgroup)}, count={repr(self.count)}, coefficient={repr(self.coefficient)}, idx={repr(self.idx)})"

  def __repr__(self):
    return self.__str__()


class FreeWilsonDecomposition:
  """FreeWilson decomposition
      rgroups - dictionary of rgroup to list of RGroups
                 i.e. {'Core': [RGroup(...), ...]
                       'R1': [ RGroup(...), RGroup(...)],
                      }
      rgroup_to_descriptor_idx - one hot encoding of smiles to descriptor index
      fitter - scikit learn compatible fitter used for regression
      N - number of rgroups
      r2 - regression r squared
      descriptors - set of the descriptors for molecules in the training set
                    used to not enumerate existing molecules
      row_decomposition - original rgroup decomposition (With row key 'molecule' is an rdkit molecule)
    """

  def __init__(self, rgroups, rgroup_to_descriptor_idx, fitter, r2, descriptors, row_decomposition,
               num_training, num_reconstructed):
    self.rgroups = rgroups  # dictionary 'Core':[core1, core1], 'R1': [rgroup1, rgroup2], ...
    self.rgroup_to_descriptor_idx = rgroup_to_descriptor_idx  # dictionary {smi:descriptor_idx}
    self.fitter = fitter  # fitter rgroup indices -> prediction
    self.N = len(rgroup_to_descriptor_idx)
    self.r2 = r2
    self.descriptors = set([tuple(d) for d in descriptors])
    self.row_decomposition = row_decomposition
    self.num_training = num_training
    self.num_reconstructed = num_reconstructed


default_decomp_params = rgd.RGroupDecompositionParameters()

# The default decomposition uses the GeneticAlgorirthm
#   fingerprint analysis to break symmetry and make
#   more consistent rgroup decompositions
default_decomp_params.matchingStrategy = rgd.GA

# Also the decomposition is allowed to add new rgroups
default_decomp_params.onlyMatchAtRGroups = False

# use the fingerprint variance method for scoring
default_decomp_params.scoreMethod = rgd.RGroupScore.FingerprintVariance

# we need to keep hydrogens so molzip will work
default_decomp_params.removeHydrogensPostMatch = False


def FWDecompose(scaffolds, mols, scores,
                decomp_params=default_decomp_params) -> FreeWilsonDecomposition:
  """
    Perform a free wilson analysis
        : param scaffolds : scaffold or list of scaffolds to use for the rgroup decomposition
        : param mols : molecules to decompose
        : param scores : list of floating point numbers for the regression (
                             you may need convert these to their logs in some cases)
        : param decomp_params : RgroupDecompositionParams default [
                                    default_decomp_params = rdkit.Chem.rdRGroupDecomposition.RGroupDecompositionParameters()
                                    default_decomp_params.matchingStrategy = rgd.GA
                                    default_decomp_params.onlyMatchAtRGroups = False
                                   ]
                                If you only want to decompose on specific group locations
                                set onlyMatchAtRGroups to True


        >>> from rdkit import Chem
        >>> from freewilson import FWBuild, FWDecompose
        >>> from rdkit.Chem import Descriptors
        >>> scaffold = Chem.MolFromSmiles("c1cccnc1")
        >>> mols = [Chem.MolFromSmiles("c1cccnc1"+"C"*(i+1)) for i in range(100)]
        >>> scores = [Descriptors.MolLogP(m) for m in mols]
        >>> fw = FWDecompose(scaffold, mols, scores)
        >>> for pred in FWBuild(fw):
        ...   print(pred)

    For an easy way to report predictions see 

       >>> from freewilson import FWBuild, predictions_to_csv
       >>> import sys
       >>> predictions_to_csv(sys.stdout, FWBuild(fw))

   
    See FWBuild docs to see how to filter predictions, molecular weight or molecular properties.
    """
  descriptors = []  # list of descriptors, one per matched molecules
  #    descriptors are 1/0 if a sidechain is present
  matched_scores = []  # scores from the matching molecules
  rgroup_idx = {}  # rgroup index into descriptor { smiles: idx }
  rgroups = defaultdict(list)  # final list of rgrups/sidechains

  if len(mols) != len(scores):
    raise ValueError(
      f"The number of molecules must match the number of scores #mols {len(mols)} #scores {len(scores)}"
    )
  # decompose the rgroups
  logger.info(f"Decomposing {len(mols)} molecules...")
  decomposer = rgd.RGroupDecomposition(scaffolds, decomp_params)
  matched = []
  matched_indices = []
  for i, (mol, score) in enumerate(tqdm(zip(mols, scores))):
    if decomposer.Add(mol) >= 0:
      matched_scores.append(score)
      matched.append(mol)
      matched_indices.append(i)

  decomposer.Process()
  logger.info(f"Matched {len(matched_scores)} out of {len(mols)}")
  if not (matched_scores):
    logger.error("No scaffolds matched the input molecules")
    return

  decomposition = decomposer.GetRGroupsAsRows(asSmiles=True)

  logger.info("Get unique rgroups...")
  blocker = rdBase.BlockLogs()
  rgroup_counts = defaultdict(int)
  num_reconstructed = 0
  for num_mols, (row, idx) in enumerate(zip(decomposition, matched_indices)):
    row_smiles = []
    for rgroup, smiles in row.items():
      row_smiles.append(smiles)
      rgroup_counts[smiles] += 1
      if smiles not in rgroup_idx:
        rgroup_idx[smiles] = len(rgroup_idx)
        rgroups[rgroup].append(RGroup(smiles, rgroup, 0, 0))
    row['original_idx'] = idx
    reconstructed = ".".join(row_smiles)
    try:
      blocker = rdBase.BlockLogs()
      mol = molzip_smi(reconstructed)
      num_reconstructed += 1
    except:
      print("failed:", Chem.MolToSmiles(matched[num_mols]), reconstructed)

  logger.info(f"Descriptor size {len(rgroup_idx)}")
  logger.info(f"Reconstructed {num_reconstructed} out of {num_mols}")

  # get the descriptors list, one-hot encoding per rgroup
  if num_reconstructed == 0:
    logging.warning("Could only reconstruct %s out of %s training molecules", num_mols,
                    num_reconstructed)

  for mol, row in zip(matched, decomposition):
    row['molecule'] = mol
    descriptor = [0] * len(rgroup_idx)
    descriptors.append(descriptor)
    for smiles in row.values():
      if smiles in rgroup_idx:
        descriptor[rgroup_idx[smiles]] = 1

  assert len(descriptors) == len(
    matched_scores
  ), f"Number of descriptors({len(descriptors)}) doesn't match number of matcved scores({len(matched_scores)})"

  # Perform the Ridge Regression
  logger.info("Ridge Regressing...")
  lm = Ridge()
  lm.fit(descriptors, matched_scores)
  preds = lm.predict(descriptors)
  r2 = r2_score(matched_scores, preds)
  logger.info(f"R2 {r2}")
  logger.info(f"Intercept = {lm.intercept_:.2f}")

  for sidechains in rgroups.values():
    for rgroup in sidechains:
      rgroup.count = rgroup_counts[rgroup.smiles]
      rgroup.coefficient = lm.coef_[rgroup_idx[rgroup.smiles]]
      rgroup.idx = rgroup_idx[rgroup.smiles]

  return FreeWilsonDecomposition(rgroups, rgroup_idx, lm, r2, descriptors, decomposition, num_mols,
                                 num_reconstructed)


def _enumerate(rgroups, fw, mw_filter=None, hvy_filter=None, pred_filter=None, mol_filter=None):
  N = fw.N
  fitter = fw.fitter
  num_products = 1
  for r in rgroups:
    num_products *= len(r)

  wrote = 0
  in_training_set = 0
  rejected_pred = 0
  rejected_mw = 0
  rejected_hvy = 0
  good_pred = 0
  rejected_filters = 0
  rejected_bad = 0
  min_mw = 10000000
  max_mw = 0
  min_hvy = 10000000
  max_hvy = 0
  max_pred = -1e10
  min_pred = 1e10
  delta = num_products // 10 or 1
  for i, groups in tqdm(enumerate(itertools.product(*rgroups)), total=num_products):
    if i and i % delta == 0:
      logging.debug(
        f"Wrote {wrote} results out of {num_products}\n\t\n\tIn Training Set{in_training_set}\n\tBad MW: {rejected_mw}\n\tBad Pred: {rejected_pred}\n\tBad Filters: {rejected_filters}\n\tBad smi: {rejected_bad}\n\tmin mw: {min_mw}\n\tmax mw: {max_mw}\n\t\n\tmin pred: {min_pred}\n\tmax pred: {max_pred}",
        file=sys.stderr)

    mw = 0
    hvy = 0
    descriptors = [0] * N
    for g in groups:
      descriptors[g.idx] = 1
      mw += g.mw
      hvy += g.hvyct

    if tuple(descriptors) in fw.descriptors:
      in_training_set += 1
      continue

    min_mw = min(min_mw, mw)
    max_mw = max(max_mw, mw)
    if mw_filter and not mw_filter(mw):
      rejected_mw += 1
      continue

    min_hvy = min(min_hvy, hvy)
    max_hvy = max(max_hvy, hvy)
    if hvy_filter and not hvy_filter(hvy):
      rejected_hvy += 1
      continue

    pred = fitter.predict([descriptors])[0]
    max_pred = max(pred, max_pred)
    min_pred = min(pred, min_pred)
    if pred_filter and not pred_filter(pred):
      rejected_pred += 1
      continue
    good_pred += 1
    smiles = set([g.smiles for g in groups])  # remove dupes
    smi = ".".join(set([g.smiles for g in groups]))
    try:
      mol = molzip_smi(smi)
    except:
      rejected_bad += 1
      continue

    rejected = False
    if mol_filter and not mol_filter(mol):
      rejected_filters += 1
      continue

    out_smi = Chem.MolToSmiles(mol)
    yield FreeWilsonPrediction(pred, out_smi, groups)
    wrote += 1
  logging.info(
    f"Wrote {wrote} results out of {num_products}\n\tIn Training set: {in_training_set}\n\tBad MW: {rejected_mw}\n\tBad Pred: {rejected_pred}\n\tBad Filters: {rejected_filters}\n\tBad smi: {rejected_bad}\n\tmin mw: {min_mw}\n\tmax mw: {max_mw}\n\tBad HVY: {rejected_hvy}\n\tBad Pred: {rejected_pred}\n\tBad Filters: {rejected_filters}\n\tBad smi: {rejected_bad}\n\tmin mw: {min_mw}\n\tmax mw: {max_mw}\n\tmin hvy: {min_hvy}\n\tmax hvy: {max_hvy}\n\t\n\tmin pred: {min_pred}\n\tmax pred: {max_pred}"
  )


def FWBuild(fw: FreeWilsonDecomposition, pred_filter=None, mw_filter=None, hvy_filter=None,
            mol_filter=None) -> Generator[FreeWilsonPrediction, None, None]:
  """Enumerate the freewilson decomposition and return their predictions

       :param fw: FreeWilsonDecomposition generated from FWDecompose
       :param pred_filter: return True if the prediction is in a desireable range
                           e.g.  lambda pic50: pic50>=8
       :param mw_filter: return True if the enumerated molecular weight is in a desireable rrange
                           e.g. lambda mw: 150 < mw < 550
       :param hvy_filter: return True if the enumerated heavy couont is in a desireable rrange
                           e.g. lambda hvy: 10 < hvy < 50
       :param mol_filter: return True if the molecule is ok to be enumerated
                           e.g. lambda mol: -3 < Descriptors.MolLogp(mol) < 5
    """
  blocker = rdBase.BlockLogs()
  # check groups for cycles
  cycles = set()
  rgroups_no_cycles = defaultdict(list)
  rgroup_cycles = defaultdict(list)

  # we can handle cycles now?
  for key, rgroup in fw.rgroups.items():
    if key == 'Core':
      rgroups_no_cycles[key] = rgroup
      continue
    no_cycles = rgroups_no_cycles[key]
    for g in rgroup:
      no_cycles.append(g)
      continue

      if len(g.dummies) > 1:
        cycles.add(g.dummies)
        rgroup_cycles[g.dummies].append(g)
      else:
        no_cycles.append(g)

  logging.info("Enumerating rgroups with no broken cycles...")
  for k, v in rgroups_no_cycles.items():
    logging.info(f"\t{k}\t{len(v)}")
  # do the ones that have no cycles first
  rgroups = [rgroup for key, rgroup in sorted(rgroups_no_cycles.items())]
  # core is always first
  for res in _enumerate(rgroups, fw, pred_filter=pred_filter, mw_filter=mw_filter,
                        hvy_filter=hvy_filter, mol_filter=mol_filter):
    yield res

  # iterate on rgroups with cycles
  #  basically only let one set of RGroups show up once.
  indices = set()
  for k in fw.rgroups:
    if k[0] == "R":
      indices.add(int(k[1:]))
  if cycles:
    logging.info("Enumerating rgroups with broken cycles...")

  for rgroup_indices in cycles:
    missing = indices - set(rgroup_indices)
    rgroups = {'Core': fw.rgroups['Core']}
    rgroups["R%s" % ".".join([str(x) for x in rgroup_indices])] = rgroup_cycles[rgroup_indices]
    for m in missing:
      k = "R%s" % m
      rgroups[k] = rgroups_no_cycles[k]

    for k, v in rgroups.items():
      logging.info(f"\t{k}\t{len(v)}")

    for res in _enumerate([rgroup for key, rgroup in sorted(rgroups.items())], fw,
                          pred_filter=pred_filter, mw_filter=mw_filter, hvy_filter=hvy_filter,
                          mol_filter=mol_filter):
      yield res


def _rgroup_sort(r):
  """Sort groups like R1 R2 R10 not R1 R10 R2
    """
  if r[0] == "R":
    return ("R", int(r[1:]))
  return (r, None)


def predictions_to_csv(outstream, decomposition: FreeWilsonDecomposition, predictions):
  """Output predictions in csv format to the output stream

       :param outstream: output stream to write results
       :param decomposition: freewillson decomposition
       :param predictions: list of Predictions to output
    """
  writer = None
  for pred in predictions:
    if not writer:
      rgroups = set()
      for rgroup in decomposition.rgroups:
        rgroups.add(rgroup)
      rgroups = sorted(rgroups, key=_rgroup_sort)

      lookup = {}
      for i, rg in enumerate(rgroups):
        lookup[rg] = i
      writer = csv.writer(outstream)
      header = ['smiles', 'prediction'] + [f"{rg}_smiles" for rg in list(rgroups)]
      writer.writerow(header)
    rg = [""] * len(lookup)
    for s in pred.rgroups:
      rg[lookup[s.rgroup]] = s.smiles

    row = [pred.smiles, repr(pred.prediction)] + rg
    writer.writerow(row)
  return header


def test_freewilson():
  # some simple tests
  from rdkit import Chem
  from rdkit.Chem import Descriptors
  assert dummypat.findall("C[*:1]N.[H][*:2]") == ['1', '2']
  assert dummypat.findall("C[*:1]N.[HH][*:2]") == ['1', '2']
  assert dummypat.findall("C[*:1]N.[2H][*:2]") == ['1', '2']
  assert dummypat.findall("C[*:1]N.[CH2][*:2]") == ['1', '2']

  scaffold = Chem.MolFromSmiles("[*:2]c1cccnc1[*:1]")
  mols = [Chem.MolFromSmiles("N" * (i + 1) + "c1cccnc1" + "C" * (i + 1)) for i in range(10)]
  scores = [Descriptors.MolLogP(m) for m in mols]
  fw = FWDecompose(scaffold, mols, scores)

  for pred in FWBuild(fw, pred_filter=lambda x: -3 < x < 3, mw_filter=lambda mw: 100 < mw < 450,
                      hvy_filter=lambda hvy: 10 < hvy < 50,
                      mol_filter=lambda m: -3 < Descriptors.MolLogP(m) < 3):
    rgroups = set()
    for sidechain in pred.rgroups:
      rgroups.add(sidechain.rgroup)
    rgroups = sorted(rgroups)
    assert list(rgroups) == ['Core', 'R1', 'R2']


def test_rgroups():
  smiles = "[*:2]c1ccccc1[*:1]"
  rg = RGroup(smiles=smiles, rgroup="Core", count=1, coefficient=0.1)
  mol = Chem.MolFromSmiles(smiles)
  assert rg.mw == Descriptors.MolWt(mol)
  assert rg.hvyct == Descriptors.HeavyAtomCount(mol)
  assert rg.dummies == (1, 2)

  # try a non parseable smiles
  smiles = "[*:2]cccccc[*:1]"
  rg = RGroup(smiles=smiles, rgroup="Core", count=1, coefficient=0.1)
  assert rg.mw == 72.06599999999999  # almost but not quite benzene
  assert rg.hvyct == 6
  assert rg.dummies == (1, 2)

  # finally aa non parseble one
  smiles = "Nope"
  rg = RGroup(smiles=smiles, rgroup="Core", count=1, coefficient=0.1)
  assert rg.mw == 0
  assert rg.hvyct == 0
  assert rg.dummies == tuple()
