from __future__ import print_function
from collections import OrderedDict
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries
from rdkit import Chem
import sys
from collections import namedtuple


import SideChains
print (SideChains.__file__)
from SideChains import Extract
from SideChains import Extract_fixed

# Test Harness
"""
Test Harness
============

The options for the Depomposition are fairly complicated.  To make testing easier, I have made a concrete class that controls the options.

I suggest having this options class in C++ as well.
"""

# this is just like a struct in C++
DecompositionOptionsBase = namedtuple(
    "DecompositionOptionsBase",
    ["cores", 
     "labelledCores", # allow cores to have attachment points labelled
     "requireLabels", # only allow substitutions at labelled points 
     "silent", # be quiet
     "symmetrize", # attempt to symmetrize the R group assignments (experimental)
     "rejectDoubleAttachments", # do not reject molecules where sidechains end up with two attachment points
     "nonLabelledSubstituentHandling" # ['PASS','IGNORE','FAIL']
 ])
        
class DecompositionOptions(DecompositionOptionsBase):
    def __new__(cls, 
                 cores, 
                 labelledCores=False,
                 requireLabels=False,
                 silent=False,
                 symmetrize=False,
                 rejectDoubleAttachments=True,
                 nonLabelledSubstituentHandling="PASS"):
        assert labelledCores in [False, True]
        assert requireLabels in [False, True]
        assert silent in [False, True]
        assert symmetrize in [False, True]
        assert rejectDoubleAttachments in [False, True]
        assert nonLabelledSubstituentHandling in ['PASS','IGNORE','FAIL']
        d = locals().copy()
        del d["cls"]
        self = super(DecompositionOptions, cls).__new__(cls, **d)
        return self


"""
Outputting results
=============
The lists returned by the decomposition object are a bit hard to understand, so let's make a dump function that can print them out in a readable way.

Here is a rundown of the results:

```
decomposition = RunDecomposition(mols, options)
decomposition[mol_idx] is the result for mol[mol_idx]
decomposition[mol_idx][core_idx] is the result for mol[mol_idx] 
                                 and options.cores[core_idx]

side_chain = decomposition[mol_idx][core_idx]
side_chain => [ (core_atom_idx, sidechain), ... ]
  is the list of sidechains originating at core_atom_idx for core core_idx.
```
"""

def dumpResults(decomp_results):
    mols, results = decomp_results
    for index, mol_results in enumerate(results):
        print ("Molecule: %s"%index)
        
        for core_index, sidechains in enumerate(mol_results):
            print("\tcore: %s"%core_index,)
            #print ("\t", sidechains)
            for core_atom_idx, sidechain in sidechains:
                print ("\t\tcore_atom_idx: %s sidechain: %s"%(
                    core_atom_idx, sidechain))

print ("*"*44)
labelledCores = False
print ("Smarts Cores: Original Broken Results")
print ("Dummy atoms are considered part of the scaffold, not the sidechain")

smis = ['CN']
smarts_cores = ['C*']
print("Smarts Cores: ", smarts_cores)
print("Input Smiles:", smis)
print("LabelledCores:", labelledCores)
mols = [Chem.MolFromSmiles(x) for x in smis]
cores = [Chem.MolFromSmarts(core) for core in smarts_cores]

dumpResults(Extract.RunDecomposition(mols,
            DecompositionOptions(cores, labelledCores=labelledCores)))
print()

smis = ['CNN']
smarts_cores = ['C*']
print("Smarts Cores: ", smarts_cores)
print("Input Smiles:", smis)
print("LabelledCores:", labelledCores)
mols = [Chem.MolFromSmiles(x) for x in smis]
cores = [Chem.MolFromSmarts(core) for core in smarts_cores]
dumpResults(Extract.RunDecomposition(mols,
                DecompositionOptions(cores, labelledCores=labelledCores)))
print()


print ("*"*44)
print ("Fixed Results (Expected Results)")

smis = ['CN']
smarts_cores = ['C*']
print("Smarts Cores: ", smarts_cores)
print("Input Smiles:", smis)
print("LabelledCores:", labelledCores)
mols = [Chem.MolFromSmiles(x) for x in smis]
cores = [Chem.MolFromSmarts(core) for core in smarts_cores]

dumpResults(Extract_fixed.RunDecomposition(mols,
            DecompositionOptions(cores, labelledCores=labelledCores)))
print()

smis = ['CNN']
smarts_cores = ['C*']
print("Smarts Cores: ", smarts_cores)
print("Input Smiles:", smis)
print("LabelledCores:", labelledCores)
mols = [Chem.MolFromSmiles(x) for x in smis]
cores = [Chem.MolFromSmarts(core) for core in smarts_cores]
dumpResults(Extract_fixed.RunDecomposition(mols,
                DecompositionOptions(cores, labelledCores=labelledCores)))
print()


print("*"*44)
print ("Smiles Cores")
print ("Broken Results: Atom mappings can be lost if two mappings are attached ")
print("  to the same atom")
print ("")

labelledCores = True

smis = ['C(N)(B)']
smiles_cores = ['C([1*])([2*])']
print("Smiles Cores: ", smiles_cores)
print("Input Smiles:", smis)
print("LabelledCores:", labelledCores)

mols = [Chem.MolFromSmiles(x) for x in smis]
cores = [Chem.MolFromSmiles(smi) for smi in smiles_cores]
dumpResults(Extract.RunDecomposition(mols,
                DecompositionOptions(cores, labelledCores=labelledCores)))

print ("Fixed Results")
mols = [Chem.MolFromSmiles(x) for x in smis]
cores = [Chem.MolFromSmiles(smi) for smi in smiles_cores]
dumpResults(Extract_fixed.RunDecomposition(mols,
                DecompositionOptions(cores, labelledCores=labelledCores)))

smis = ['C(N)(B)']
smiles_cores = ['C([2*])([1*])']
print("Smiles Cores: ", smiles_cores)
print("Input Smiles:", smis)
print("LabelledCores:", labelledCores)
mols = [Chem.MolFromSmiles(x) for x in smis]
cores = [Chem.MolFromSmiles(smi) for smi in smiles_cores]
dumpResults(Extract_fixed.RunDecomposition(mols,
                DecompositionOptions(cores, labelledCores=labelledCores)))

print("*"*44)
print("Stereo substructure searching works in github as of 1/13/2016")
smis = ["S[C@](Br)([F])"]
smiles_cores = ['[1*][C@]([2*])([3*])']
print("Smiles Cores: ", smiles_cores)
print("Input Smiles:", smis)
print("LabelledCores:", labelledCores)

mols = [Chem.AddHs(Chem.MolFromSmiles(x)) for x in smis]
cores = [Chem.MolFromSmiles(smi) for smi in smiles_cores]
dumpResults(Extract_fixed.RunDecomposition(mols,
                DecompositionOptions(cores, labelledCores=labelledCores)))

print()
print("However, reordering to match input order seems to bork")
print(" stereo.  I expect for stero centers, the best we can")
print(" do is rotate... not simply reorder")
smis = ["S[C@]([F])(Br)"]
smiles_cores = ['[1*][C@]([2*])([3*])']
print("Smiles Cores: ", smiles_cores)
print("Input Smiles:", smis)
print("LabelledCores:", labelledCores)
mols = [Chem.AddHs(Chem.MolFromSmiles(x)) for x in smis]
cores = [Chem.MolFromSmiles(smi) for smi in smiles_cores]
dumpResults(Extract_fixed.RunDecomposition(mols,
                DecompositionOptions(cores, labelledCores=labelledCores)))

