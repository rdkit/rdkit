# DASH trees in the RDKit

DASH (Dynamic Attention-based Substructure Hierarchy) assigns partial charges
and other atomic properties by walking a precomputed decision tree: each level
adds one neighbouring atom to the matched substructure, and the node where the
descent stops carries the value. It is the method of the python
[DASH-tree](https://github.com/rinikerlab/DASH-tree) package
([JCIM 2023](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00800),
[J. Chem. Phys. 2024](https://doi.org/10.1063/5.0218154)), reimplemented here
in C++, a few hundred times faster, and reading its tree from a memory-mapped
file instead of 3 GB of heap. It gives the package's results to the last bit,
except where several atoms fit the same node, which it handles differently on
purpose; see the section of that name.

The tree is data, not code. It lives in a `.dash` container, RDKit's own format,
whose layout is defined by the structs in `DASHTreeImpl.h` and mirrored in
`rdkit.Chem.DASHTree.container`. The python package's trees come in a different,
older layout, called the *legacy* format here, which the converter below turns
into a container.

## Quick start

```python
from rdkit import Chem
from rdkit.Chem.DASHTree import GetDASHTree

mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))  # hydrogens carry charge: make them explicit
tree = GetDASHTree()                                       # the tree shipped with the RDKit
charges = tree.GetPartialCharges(mol)                      # one float per atom, summing to the formal charge
```

With no argument `GetDASHTree()` uses `$RDKIT_DASH_TREE` if set and otherwise
the tree in `Data/DASHTree`. That tree is the legacy MBIS
charge tree pruned to 38 MB, see below; its charges stay within 0.044 e of the
full tree's (0.0068 e RMS), and it carries the `result` and `std` columns only.
For the full tree, or for properties other than MBIS charges, convert one and
pass it in, as a path or an http(s) URL that is fetched once into a cache:

```python
tree = GetDASHTree("default.dash", ["result", "std"])      # resolve only the columns you need
```

## Making a tree

Two command-line tools come with the python package, and `--help` explains
every option.

**Convert a legacy tree.** A legacy tree is a folder of `0.gz, 0.h5, 1.gz,
1.h5, ...`, one pair per atom-feature branch, as the DASH-tree package stores
them. In a checkout of its repository the MBIS charge tree is
`serenityff/charge/data/default_dash_tree`; the tree with every property
(AM1BCC, RESP, Mulliken, DFT-D4 C6 and polarizability, dual descriptors,
dipoles) is downloaded by the package into
`serenityff/charge/data/additional_data/dashProps`. The converter needs
`pandas` and `pytables`, not the DASH-tree package. The `.gz` files are
pickles, so only convert trees from a source you trust.

```
python -m rdkit.Chem.DASHTree.convert <tree_folder> default.dash
python -m rdkit.Chem.DASHTree.convert <tree_folder> charges.dash --props result,std
```

The first keeps every column (189 MB for the MBIS tree, a few minutes); the
second is the 151 MB charges-only container to hand around.

**Prune a tree.** The pruner writes a smaller container. Dropping a node
drops its subtree; an atom whose descent reaches the cut reads the deepest node
still there, so a pruned tree answers for every molecule the full one accepts,
with a more general substructure where nodes are missing.

```
python -m rdkit.Chem.DASHTree.prune default.dash small.dash --tolerance 0.02 --props result,std
```

Measured against the full MBIS tree on the 100 molecules in `test_data`
(5019 atoms):

| `--tolerance` | container | RMSD | max error |
|---|---|---|---|
| 0.01 e | 63.7 MB | 0.0036 e | 0.036 e |
| 0.02 e (shipped) | 38.2 MB | 0.0068 e | 0.044 e |
| 0.05 e | 14.3 MB | 0.014 e | 0.058 e |
| 0.1 e | 4.2 MB | 0.024 e | 0.10 e |

`Data/DASHTree/default_pruned.dash` is the 0.02 e output of the first convert
command and this prune command, run on the MBIS tree of the DASH-tree package
(MIT licence), the most accurate of these under GitHub's 50 MB per-file warning
(22 MB as git stores it). Both steps are deterministic, so it is reproducible
byte for byte.

**Train one.** The DASH-tree package builds trees from a trained attention
model and writes them in the legacy format.

## Using a tree from python

```python
from rdkit import Chem
from rdkit.Chem import rdDASHTree
from rdkit.Chem.DASHTree import GetDASHTree, GetAtomMatchSmarts

tree = GetDASHTree("charges.dash", ["result", "std"])
mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))

# charges, and what they were made from
charges = tree.GetPartialCharges(mol)
details = tree.GetPartialChargesDetails(mol)        # 'charges', 'raw', 'std', 'matchDepth'

# how the charges are made to sum to the formal charge
options = rdDASHTree.NormalizationOptions()
options.normalization = rdDASHTree.Normalization.SYMMETRIC   # or NONE, or STD_WEIGHTED (default)
charges = tree.GetPartialCharges(mol, options)

# many molecules at once; 0 means every core
batch = tree.GetPartialChargesBatch(mols, numThreads=0)

# the substructure an atom's descent matched, as SMARTS with the atom mapped as :1
smarts = GetAtomMatchSmarts(tree, mol, 0)           # [#6X4H3+0:1]-[#6X3H0+0](~[#8X1H0+0])~[#7X3H1+0]~...
path = tree.GetAtomNodePath(mol, 0)                 # branch index, then node ids
node = tree.GetNode(path[0], path[-1])              # GetValue('result'), iterable children
```

A tree with more columns works the same way; name the columns you want when you
open it and only those are read:

```python
props = GetDASHTree("props.dash", ["AM1BCC", "AM1BCC_std", "DFTD4:C6"])
c6 = props.GetMolProperty(mol, "DFTD4:C6")          # one value per atom, no normalisation

# any per-atom property with a molecule-wide total normalises like the charges do
options = rdDASHTree.NormalizationOptions()
options.valueProperty, options.stdProperty = "AM1BCC", "AM1BCC_std"
am1bcc = props.GetNormalizedMolProperty(mol, Chem.GetFormalCharge(mol), options)
```

An atom gets the value of the deepest node on its path that has one, for every
property. A node without one holds NaN, or the smallest int32 in an integer
column.

Three things to know. Hydrogens have charges of their own, so molecules need
explicit hydrogens (`Chem.AddHs`). The atom classes read the conjugation flags
that sanitization sets, so an unsanitized molecule is assigned with a warning.
The tree knows a fixed set of atom classes, 122 in the published trees; a
molecule with an atom outside them (silicon, a metal, an unusual valence) is
refused with the atom named.

## Using a tree from C++

```cpp
#include <GraphMol/DASHTree/DASHTree.h>

RDKit::DASH::DASHTree tree("charges.dash", {"result", "std"});
std::vector<double> charges;
tree.getPartialCharges(*mol, charges);
tree.getNormalizedMolProperty(*mol, target, values, options);
tree.getPartialChargesBatch(mols, batch, RDKit::DASH::NormalizationOptions(), 0);
```

The library target is `DASHTree`, gated on `RDK_USE_BOOST_IOSTREAMS` for the
memory mapping. A tree is immutable once opened and safe to query from any
number of threads.

## Atoms that fit the same node

Each level of the descent binds an atom carrying the key of the next node: its
atom class, the matched atom it hangs off, and the bond. Sometimes several
atoms do, like the two CH2 next to the ring NH of 5-aminoisoindoline, or the
ortho carbons of a ring. The descent then follows each of them until the tree
tells them apart, and takes the child that comes first; a node's children are
in falling order of the attention they were built with. So the result does not
depend on the atom order, and symmetry-equivalent atoms get the same value.

The DASH-tree package binds the first such atom in atom order instead, so two
SMILES of one molecule can get different values, and so can pyridine's two
meta carbons. `DASHParams.legacyTieBreaking` does the same, for exact agreement
with the package; over the 100 reference molecules the two rules put
std-weighted charges 0.011 e apart RMS, 0.24 e at most.

```python
options = rdDASHTree.NormalizationOptions()
options.params.legacyTieBreaking = True
charges = tree.GetPartialCharges(mol, options)      # the DASH-tree package's values
```
