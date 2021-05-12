FreeWilson
==========
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

Scores need to be in an appropriate form for regrerssion analysis, i.e. pIC50s as opposed to IC50s.
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
>>> for pred in FWBuild(decomp):
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

Using Molecule filters
----------------------
Finally, the molecule filter can be used to prune greasy or otherwise undesirable
molecules:

```
>>> from rdkit.Chem import Descriptors
>>> for pred in FWBuild(decomp, pred_filter=lambda x: x > 8,
...                             mol_filter=lambda mol: -3 < Descriptors.MolLogP(mol) < 3):
...   print(pred.smiles, pred.prediction)

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
>>> mcs = rdFMCS.FindMCS(mols, threshold=0.8, atomCompare=rdFMCS.AtomCompare.CompareAny,
...                      completeRingsOnly=True)
>>> decomp = FWDecompose(mcs.queryMol, mols, scores)

```

Note that the MCS returned can generate multiple Cores, this freewilson implementatino
explicitly supports this.

Using Molecule filters
----------------------
Finally, the molecule filter can be used to prune greasy or otherrwise undesirable
molecules:

```
>>> from rdkit.Chem import Descriptors
>>> for pred in FWBuild(decomp, pred_filter=lambda x: x > 8,
...                             mol_filter=lambda mol: -3 < Descriptors.MolLogP(mol) < 3):
...   print(pred.smiles, pred.prediction)
```


Testing
=======

From this directory:

```
PYTHONPATH=`pwd` pytest
```
