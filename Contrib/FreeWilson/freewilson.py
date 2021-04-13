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


```
>>> import sys
>>> freewilson impor FWDecompose, FWBuild, predictions_to_csv
>>> decomp = FWDecompose(scaffold, mols, scores)
>>> preds = FWBuild(fdecomp, 
...                    pred_filter=lambda x: x > 8, 
...                    mw_filter=lambda mw: 100<mw<550)
>>> predictions_to_csv(sys.stdout, preds)
```

Scores need to be in an appropriate form for regrerssion analysis, i.e. pIC50s as opposed to IC50s.
Enumerations are a lot faster if you know a prediction cutoff or molecular weight cutoff.

Scaffolds can be any scaffold or smarts pattern or list of either.  The system automatically 
enumerates matching cores in the analysis.

To see if the decomposition is useful, check the r squared value

```
>>> print(f"Training R^2 is {decomp.r2}")
```

You can also get some more information by setting logging to INFO

```
>>> import logging
>>> logging.getLogger().setLevel(logging.INFO)
>>> preds = list(fw.FWBuild(free, pred_filter=lambda x: x > 8))
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
>>> mcs = rdFMCS.FindMCS(mols, threshold=0.8, atomCompare=rdFMCS.AtomCompare.CompareAny,
...                      completeRingsOnly=True)
>>> decomp = FWDecompose(mcs.queryMol, mols, scores)
```

Using Molecule filters
----------------------
Finally, the molecule filter can be used to prune greasy or otherrwise undesirable
molecules:

```
from rdkit.Chem import Descriptors
>>> for pred in fw.FWBuild(free, mol_filter=lambda mol: -3 < Descriptors.MolLogP(mol) < 3):
...   print(pred.smiles, pred.prediction)
"""

from rdkit.Chem import rdRGroupDecomposition as rgd
from rdkit.Chem import molzip, Descriptors
from rdkit import rdBase
from rdkit import Chem
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from collections import defaultdict, namedtuple
from tqdm import tqdm
import itertools
import math
import sys
from typing import List
import logging
import csv
import re
from typing import Generator
logger = logging.getLogger("rdkit.Chem.FreeWilson")

FreeWilsonPrediction = namedtuple("FreeWilsonPrediction", ['prediction', 'smiles', 'rgroups'])

# hydrogens with dummy atoms
hspat = re.compile(r"\[[H0-9]+\]\[\*:([0-9]+)")
# all dummy atoms
dummypat = re.compile(r"\*:([0-9]+)")


class RGroup:
    """FreeWilson RGroup
        smiles - isomeric smiles of the rgroup
        rgroup - Rgroup name ['Core', 'R1', 'R2'...]
        count - number of molecules with the rgroup
        coefficient - ridge regression coefficient
        idx - one-hot encoding for the rgroup
    """
    def __init__(self, smiles, rgroup, count, coefficient, idx=None):
        self.smiles = smiles           # smiles for the sidechain (n.b. can be a core as well)
        self.rgroup  = rgroup          # rgroup Core, R1, R2,...
        self.count = count             # num molecules with this rgruop
        self.coefficient = coefficient # ridge coefficient
        self.idx = idx                # descriptor index
        hs = [int(x) for x in hspat.findall(smiles)]
        assert len(hs) <= 1
        self.hs = hs                   # hydrogen with dummy atoms    
        self.dummies = tuple([int(x) for x in sorted(dummypat.findall(smiles))])

        # Assemble some additive properties 
        try:
            # will fail on some structures
            self.mw = Descriptors.MolWt(Chem.MolFromSmiles(smiles))
        except:
            # guess the MW if we can't sanitize
            table = Chem.GetPeriodicTable()
            try:
                self.mw = 0.
                for atom in Chem.MolFromSmiles(smiles, sanitize=False).GetAtoms():
                    self.mw += table.GetAtomicWeight(atom.GetAtomicNum())
            except:
                self.mw = 0 # dunno
        
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
    """
                      
    def __init__(self, rgroups, rgroup_to_descriptor_idx, fitter, r2, descriptors):
        self.rgroups = rgroups # dictionary 'Core':[core1, core1], 'R1': [rgroup1, rgroup2], ...
        self.rgroup_to_descriptor_idx = rgroup_to_descriptor_idx # dictionary {smi:descriptor_idx}
        self.fitter = fitter # fitter rgroup indices -> prediction
        self.N = len(rgroup_to_descriptor_idx)
        self.r2 = r2
        self.descriptors = set([tuple(d) for d in descriptors])


default_decomp_params = rgd.RGroupDecompositionParameters()

# The default decomposition uses the GeneticAlgorirthm
#   fingerprint analysis to break symmetry and make
#   more consistent rgroup decompositions
default_decomp_params.matchingStrategy = rgd.GA

# Also the decomposition is allowed to add new rgroups
default_decomp_params.onlyMatchAtRGroups = False

def FWDecompose(scaffolds, mols, scores, decomp_params=default_decomp_params) -> FreeWilsonDecomposition:
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
        >>> from rdkit.Chem import Descriptors, FreeWilson
        >>> scaffold = Chem.MolFromSmiles("c1cccnc1")
        >>> mols = [Chem.MolFromSmiles("c1cccnc1"+"C"*(i+1)) for i in range(100)]
        >>> scores = [Descriptors.MolLogP(m) for m in mols]
        >>> fw = FWDecompose(scaffold, mols, scores)
        >>> for pred in FWBuild(fw):
        >>>  ...

    To filter out molecules in the training set, send the freewilson decomposition into
    the builder.

       >>> for pred in FWBuild(fw, decomposition=fw):
       >>>    ...

    For an easy way to report predictions see 

       >>> import sys
       >>> predictions_to_csv(sys.stdout, fw, FWBuild(fw))

   
    See FWBuild docs to see how to filter predictions, molecular weight or molecular properties.
    """
    descriptors = [] # list of descriptors, one per matched molecules
                                   #    descriptors are 1/0 if a sidechain is present
    matched_scores = []            # scores from the matching molecules
    rgroup_idx = {}   # rgroup index into descriptor { smiles: idx }
    rgroups = defaultdict(list) # final list of rgrups/sidechains

    if len(mols) != len(scores):
        raise ValueError(f"The number of molecules must match the number of scores #mols {len(mols)} #scores {len(scores)}")
    # decompose the rgroups
    logger.info(f"Decomposing {len(mols)} molecules...")
    decomposer = rgd.RGroupDecomposition(scaffolds, decomp_params)
    for mol, score in tqdm(zip(mols, scores)):
        if decomposer.Add(mol) >= 0:
            matched_scores.append(score)
    decomposer.Process()
    logger.info(f"Matched {len(matched_scores)} out of {len(mols)}")
    if not(matched_scores):
        logger.error("No scaffolds matched the input molecules")
        return

    decomposition = decomposition = decomposer.GetRGroupsAsRows(asSmiles=True)
    logger.info("Get unique rgroups...")
    rgroup_counts = defaultdict(int)
    for row in decomposition:
        for rgroup,smiles in row.items():
            rgroup_counts[smiles] += 1
            if smiles not in rgroup_idx:
                rgroup_idx[smiles] = len(rgroup_idx)
                rgroups[rgroup].append(RGroup(smiles, rgroup, 0, 0))

    logger.info(f"Descriptor size {len(rgroup_idx)}")
    # get the descriptors list, one-hot encoding per rgroup
    for row in decomposition:
        descriptor = [0] * len(rgroup_idx)
        descriptors.append(descriptor)
        for smiles in row.values():
            if smiles in rgroup_idx:
                descriptor[rgroup_idx[smiles]] = 1

    assert len(descriptors) == len(matched_scores), f"Number of descriptors({len(descriptors)}) doesn't match number of matcved scores({len(matched_scores)})"

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

    return FreeWilsonDecomposition(rgroups, rgroup_idx, lm, r2, descriptors)

def _enumerate(rgroups, fw, 
               mw_filter=None, pred_filter=None, mol_filter=None):
    N = fw.N
    fitter = fw.fitter
    num_products = 1
    for r in rgroups:
        num_products*=len(r)

    wrote = 0
    in_training_set = 0
    rejected_pred = 0
    rejected_mw = 0
    good_pred = 0
    rejected_filters = 0
    rejected_bad = 0
    min_mw = 10000000
    max_mw = 0
    max_pred = -1e10
    min_pred = 1e10
    delta = num_products//10 or 1
    for i,groups in tqdm(enumerate(itertools.product(*rgroups)), total=num_products):
        if i and i % delta == 0:
            logging.debug(f"Wrote {wrote} results out of {num_products}\n\t\n\tIn Training Set{in_training_set}\n\tBad MW: {rejected_mw}\n\tBad Pred: {rejected_pred}\n\tBad Filters: {rejected_filters}\n\tBad smi: {rejected_bad}\n\tmin mw: {min_mw}\n\tmax mw: {max_mw}\n\t\n\tmin pred: {min_pred}\n\tmax pred: {max_pred}", file=sys.stderr)

        mw = 0
        descriptors = [0] * N
        for g in groups:
            descriptors[g.idx] = 1
            mw += g.mw
        if tuple(descriptors) in fw.descriptors:
            in_training_set += 1
            continue
        
        min_mw = min(min_mw, mw)
        max_mw = max(max_mw, mw)
        if mw_filter and not mw_filter(mw):
            rejected_mw += 1
            continue
        pred = fitter.predict([descriptors])[0]
        max_pred = max(pred, max_pred)
        min_pred = min(pred, min_pred)
        if pred_filter and not pred_filter(pred):
            rejected_pred += 1
            continue
        good_pred+=1
        smiles = set([g.smiles for g in groups]) # remove dupes
        smi = ".".join(set([g.smiles for g in groups]))
        try:
            mol = molzip(Chem.MolFromSmiles(smi))
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
    logging.info(f"Wrote {wrote} results out of {num_products}\n\tIn Training set: {in_training_set}\n\tBad MW: {rejected_mw}\n\tBad Pred: {rejected_pred}\n\tBad Filters: {rejected_filters}\n\tBad smi: {rejected_bad}\n\tmin mw: {min_mw}\n\tmax mw: {max_mw}\n\t\n\tmin pred: {min_pred}\n\tmax pred: {max_pred}")

    
def FWBuild(fw: FreeWilsonDecomposition,
            pred_filter=None,
            mw_filter=None,
            mol_filter=None) -> Generator[FreeWilsonPrediction,None,None]:
    """Enumerate the freewilson decomposition and return their predictions

       :param fw: FreeWilsonDecomposition generated from FWDecompose
       :param pred_filter: return True if the prediction is in a desireable range
                           e.g.  lambda pic50: pic50>=8
       :param mw_filter: return True if the enumerated molecular weight is in a desireable rrange
                           e.g. lambda mw: 150 < mw < 550
       :param mol_filter: return True if the molecule is ok to be enumerated
                           e.g. lambda mol: -3 < Descriptors.MolLogp(mol) < 5
    """
    blocker = rdBase.BlockLogs()
    # check groups for cycles
    cycles = set()
    rgroups_no_cycles = defaultdict(list)
    rgroup_cycles = defaultdict(list)
    for key, rgroup in fw.rgroups.items():
        if key == 'Core':
            rgroups_no_cycles[key] = rgroup
            continue
        no_cycles = rgroups_no_cycles[key]
        for g in rgroup:
            if len(g.dummies) > 1:
                cycles.add(g.dummies)
                rgroup_cycles[g.dummies].append(g)
            else:
                no_cycles.append(g)

    logging.info("Enumerating rgroups with no broken cycles...")
    for k,v in rgroups_no_cycles.items():
        logging.info(f"\t{k}\t{len(v)}")
    # do the ones that have no cycles first
    rgroups = [rgroup for key, rgroup in sorted(rgroups_no_cycles.items())]
    # core is always first
    for res in _enumerate(rgroups,
                          fw,
                          pred_filter=pred_filter,
                          mw_filter=mw_filter,
                          mol_filter=mol_filter):
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
        rgroups["R%s"%".".join([str(x) for x in rgroup_indices])] = rgroup_cycles[rgroup_indices]
        for m in missing:
            k = "R%s"%m
            rgroups[k] = rgroups_no_cycles[k]
            
        for k,v in rgroups.items():
            logging.info(f"\t{k}\t{len(v)}")
            
        for res in _enumerate([rgroup for key, rgroup in sorted(rgroups.items())],
                              fw,
                              pred_filter=pred_filter,
                              mw_filter=mw_filter,
                              mol_filter=mol_filter):
            yield res    

def _rgroup_sort(r):
    """Sort groups like R1 R2 R10 not R1 R10 R2
    """
    if r[0] == "R": return ("R", int(r[1:]))
    return (r, None)

def predictions_to_csv(outstream, predictions):
    """Output predictions in csv format to the output stream

       :param outstream: output stream to write results
       :param fw: freewillson decomposition
       :param predictions: list of Predictions to output
    """
    writer = None
    for pred in predictions:
        if not writer:
            rgroups = set()
            for sidechain in pred.rgroups:
                rgroups.add(sidechain.rgroup)
            rgroups = sorted(rgroups, key=_rgroup_sort)

            lookup = {}
            for i, rg in enumerate(rgroups):
                lookup[rg] = i
    
            writer = csv.writer(outstream)
            writer.writerow(['smiles', 'prediction'] + [f"{rg}_smiles" for rg in list(rgroups)])
        rg = [""] * len(lookup)
        for s in pred.rgroups:
            rg[lookup[s.rgroup]] = s.smiles
        row = [pred.smiles, repr(pred.prediction)] + rg
        writer.writerow(row)
            
def test_freewilson():
    # some simple tests
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    assert hspat.findall("C[*:1]N.[H][*:2]") == ['2']
    assert hspat.findall("C[*:1]N.[HH][*:2]") == ['2']
    assert hspat.findall("C[*:1]N.[2H][*:2]") == ['2']
    assert hspat.findall("C[*:1]N.[CH2][*:2]") == []

    assert dummypat.findall("C[*:1]N.[H][*:2]") == ['1', '2']
    assert dummypat.findall("C[*:1]N.[HH][*:2]") == ['1', '2']
    assert dummypat.findall("C[*:1]N.[2H][*:2]") == ['1', '2']
    assert dummypat.findall("C[*:1]N.[CH2][*:2]") == ['1', '2']

    scaffold = Chem.MolFromSmiles("[*:2]c1cccnc1[*:1]")
    mols = [Chem.MolFromSmiles("N"*(i+1) + "c1cccnc1"+"C"*(i+1)) for i in range(10)]
    scores = [Descriptors.MolLogP(m) for m in mols]
    fw = FWDecompose(scaffold, mols, scores)

    for pred in FWBuild(fw,
                        pred_filter=lambda x: -3 < x < 3,
                        mw_filter=lambda mw: 100<mw<450,
                        mol_filter=lambda m: -3 < Descriptors.MolLogP(m) < 3):
        rgroups = set()
        for sidechain in pred.rgroups:
            rgroups.add(sidechain.rgroup)
        rgroups = sorted(rgroups)
        assert list(rgroups) == ['Core', 'R1', 'R2']

