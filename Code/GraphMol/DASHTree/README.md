# DASH trees in the RDKit

DASH (Dynamic Attention-based Substructure Hierarchy) assigns partial charges
and other atomic properties with a precomputed decision tree. Each level adds
one neighbour atom to the matched substructure, and the node where the descent
stops carries the value. It is the method of the python
[DASH-tree](https://github.com/rinikerlab/DASH-tree) package
([JCIM 2023](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00800),
[J. Chem. Phys. 2024](https://doi.org/10.1063/5.0218154)), reimplemented here
in C++. This module is a few hundred times faster, and it reads its tree from a
memory-mapped file instead of 3 GB of heap. It gives the package's results to
the last bit, except where several atoms fit the same node. It handles that
case differently on purpose, as
[Atoms that fit the same node](#atoms-that-fit-the-same-node) explains.

The tree is a data file, a `.dash` container in RDKit's own format. The structs
in `DASHTreeImpl.h` define its layout, and `rdkit.Chem.DASHTree.container`
mirrors it. The DASH-tree package stores its trees in an older layout, which
this README calls the *legacy* format. The converter below writes a container
from a legacy tree.

## Quick start

```python
from rdkit import Chem
from rdkit.Chem.DASHTree import GetDASHTree

mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))  # hydrogens carry charge, so make them explicit
tree = GetDASHTree()                                       # the initial_test tree, see below
charges = tree.GetPartialCharges(mol)                      # one float per atom, with the formal charge as their sum
```

With no argument, `GetDASHTree()` uses `$RDKIT_DASH_TREE` if it is set, and
`initial_test` if not. The RDKit repository does not contain the trees.
`GetDASHTree()` knows them by name (`rdkit.Chem.DASHTree.availableTrees`). The
build downloads `initial_test` into `Data/DASHTree`. If a named tree is not
there, `GetDASHTree()` downloads it once into a cache.

`initial_test` is the legacy MBIS charge tree pruned to 38 MB, see below. It
carries the `result` and `std` columns only. For the full tree, or for
properties other than MBIS charges, convert a tree yourself. Then give
`GetDASHTree()` its path, or an http(s) URL that it downloads once into the
cache:

```python
tree = GetDASHTree("default.dash", ["result", "std"])      # resolve only the columns you need
```

## Make a tree

`rdkit.Chem.DASHTree` has two command-line tools, and `--help` explains every
option.

**Convert a legacy tree.** A legacy tree is a folder of `0.gz, 0.h5, 1.gz,
1.h5, ...`, one pair per atom-feature branch, as the DASH-tree package stores
them. In a checkout of its repository, the MBIS charge tree is
`serenityff/charge/data/default_dash_tree`. The package downloads the tree with
every property (AM1BCC, RESP, Mulliken, DFT-D4 C6 and polarizability, dual
descriptors, dipoles) into `serenityff/charge/data/additional_data/dashProps`.
The converter needs `pandas` and `pytables`, not the DASH-tree package. Only
convert trees from a source you trust, because the `.gz` files are pickles.

```
python -m rdkit.Chem.DASHTree.convert <tree_folder> default.dash
python -m rdkit.Chem.DASHTree.convert <tree_folder> charges.dash --props result,std
```

The first command keeps every column (189 MB for the MBIS tree, a few minutes).
The second writes the 151 MB charges-only container to share.

**Prune a tree.** The pruner writes a smaller container. It drops a node
together with its subtree. An atom whose descent reaches the cut reads the
deepest node that is still there. So a pruned tree assigns values to every
molecule the full one accepts, with a more general substructure where nodes are
missing.

```
python -m rdkit.Chem.DASHTree.prune default.dash small.dash --tolerance 0.02 --props result,std
```

This table compares the pruned trees with the full MBIS tree on the 100
molecules in `test_data` (5019 atoms):

| `--tolerance` | container | RMSD | max error |
|---|---|---|---|
| 0.01 e (`mbis_pruned_0.01`) | 63.7 MB | 0.0036 e | 0.036 e |
| 0.02 e (`initial_test`) | 38.2 MB | 0.0068 e | 0.044 e |
| 0.05 e | 14.3 MB | 0.014 e | 0.058 e |
| 0.1 e (`mbis_pruned_0.1`) | 4.2 MB | 0.024 e | 0.10 e |

`initial_test` is the 0.02 e output of the first convert command and this prune
command, run on the MBIS tree of the DASH-tree package (MIT licence). Both steps
are deterministic, so it is reproducible byte for byte.

**Train one.** The DASH-tree package builds trees from a trained attention
model and writes them in the legacy format.

## Use a tree from python

```python
from rdkit import Chem
from rdkit.Chem import rdDASHTree
from rdkit.Chem.DASHTree import GetDASHTree, GetAtomMatchSmarts

tree = GetDASHTree("charges.dash", ["result", "std"])
mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))

# charges, and the data DASH computes them from
charges = tree.GetPartialCharges(mol)
details = tree.GetPartialChargesDetails(mol)        # 'charges', 'raw', 'std', 'matchDepth'

# the rule that makes the charges sum to the formal charge
options = rdDASHTree.NormalizationOptions()
options.normalization = rdDASHTree.Normalization.SYMMETRIC   # or NONE, or STD_WEIGHTED (default)
charges = tree.GetPartialCharges(mol, options)

# many molecules at once, numThreads=0 uses every core
mols = [mol, Chem.AddHs(Chem.MolFromSmiles("c1ccncc1"))]
batch = tree.GetPartialChargesBatch(mols, numThreads=0)

# the substructure an atom's descent matched, as SMARTS with the atom mapped as :1
smarts = GetAtomMatchSmarts(tree, mol, 0)           # [#6X4H3+0:1]-[#6X3H0+0](~[#8X1H0+0])~[#7X3H1+0]~...
path = tree.GetAtomNodePath(mol, 0)                 # branch index, then node ids
node = tree.GetNode(path[0], path[-1])              # GetValue('result'), iterable children
```

A tree with more columns works the same way. When you open it, name the columns
you want. The tree reads only those:

```python
props = GetDASHTree("props.dash", ["AM1BCC", "AM1BCC_std", "DFTD4:C6"])
c6 = props.GetMolProperty(mol, "DFTD4:C6")          # one value per atom, no normalization

# any per-atom property with a molecule-wide total normalizes like the charges do
options = rdDASHTree.NormalizationOptions()
options.valueProperty, options.stdProperty = "AM1BCC", "AM1BCC_std"
am1bcc = props.GetNormalizedMolProperty(mol, Chem.GetFormalCharge(mol), options)
```

An atom gets the value of the deepest node on its path that has one, for every
property. A node without one reads as NaN.

Three things to know about the input:

- Hydrogens have charges of their own, so molecules need explicit hydrogens
  (`Chem.AddHs`). The tree refuses a molecule without them.
- The atom classes read the conjugation flags that sanitization sets, so
  sanitize molecules first. The tree refuses a molecule with a stale property
  cache. It assigns values to an unsanitized molecule with an updated cache,
  but these values can be wrong.
- The tree knows a fixed set of atom classes, 122 in the published trees. It
  refuses a molecule with an atom outside them (silicon, a metal, an unusual
  valence), and the error names the atom.

## Use a tree from C++

```cpp
#include <GraphMol/DASHTree/DASHTree.h>

RDKit::DASH::DASHTree tree("charges.dash", {"result", "std"});
auto charges = tree.getPartialCharges(*mol);
auto values = tree.getNormalizedMolProperty(*mol, target, options);
auto batch = tree.getPartialChargesBatch(mols, RDKit::DASH::NormalizationOptions(), 0);
```

The library target is `DASHTree`. It needs `RDK_USE_BOOST_IOSTREAMS` for the
memory mapping and a little-endian host. An open tree is immutable and safe to
query from any number of threads.

## Atoms that fit the same node

Each level of the descent binds an atom that has the key of the next node: its
atom class, the matched atom it is bonded to, and the bond. Sometimes several
atoms have that key, like the two CH2 next to the ring NH of
5-aminoisoindoline, or the ortho carbons of a ring. The descent then follows
each of them until the tree distinguishes them, and takes the first child
that any of them fits. So the result does not depend on the atom order, and
symmetry-equivalent atoms get the same value. The tree builder created the
children of a node in order of attention, highest first, so the descent takes
the best-attended child that fits.

The DASH-tree package binds the first such atom in atom order instead. So two
SMILES of one molecule can get different values, and so can pyridine's two meta
carbons. `DASHParams.legacyTieBreaking` does the same, for exact agreement with
the package. Over the 100 reference molecules, the std-weighted charges of the
two rules differ by 0.011 e RMS and by 0.24 e at most.

```python
options = rdDASHTree.NormalizationOptions()
options.params.legacyTieBreaking = True
charges = tree.GetPartialCharges(mol, options)      # the DASH-tree package's values
```
