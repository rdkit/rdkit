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
carries the `result` and `std` columns only. For seven atomic properties, use
`dash_props` or `dash_props_pruned_0.02`, see
[Trees with more properties](#trees-with-more-properties). For the full MBIS
tree, convert it yourself. Then give `GetDASHTree()` its path, or an http(s)
URL that it downloads once into the cache:

```python
tree = GetDASHTree("default.dash", ["result", "std"])      # resolve only the columns you need
```

## Make a tree

`rdkit.Chem.DASHTree` has three command-line tools, and `--help` explains every
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

**Populate a tree.** The populator adds property columns to the nodes of a
tree you have. It matches each atom of each molecule and gives the atom's
value to every node on the atom's path. For each property `P`, each node then
stores the mean (`P`), the sample standard deviation (`P_std`) and the number
of values (`P_count`). The molecules come from SDF files with explicit
hydrogens. Each record gives one value per atom as an RDKit atom property list:

```python
mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
for atom, value in zip(mol.GetAtoms(), values):     # your reference values, one per atom
  atom.SetDoubleProp("q", value)
Chem.CreateAtomDoublePropertyList(mol, "q")         # writes the field 'atom.dprop.q'
with Chem.SDWriter("train.sdf") as writer:
  writer.write(mol)
```

```
python -m rdkit.Chem.DASHTree.populate default.dash q.dash train.sdf --props q
```

The descent decides which nodes an atom reaches. Populate a tree with the
tie-breaking rule that you read it with (see
[Atoms that fit the same node](#atoms-that-fit-the-same-node)). For a tree that
you read with `DASHParams.legacyTieBreaking`, add `--legacy-tie-breaking`. In
python, `rdkit.Chem.DASHTree.populate.collectStatistics()` takes any iterable
of `(mol, {name: values})` pairs. Each value array can have one row per
conformer. `writeTree()` writes the result.

**Train one.** The DASH-tree package builds trees from a trained attention
model and writes them in the legacy format.

## Trees with more properties

`dash_props` is the MBIS tree, populated with seven atomic properties of the
DASH-properties data set (ETH Research Collection,
[doi:10.3929/ethz-b-000670546](https://doi.org/10.3929/ethz-b-000670546),
871,906 conformers of 293,822 molecules). The columns are the MBIS charge
(`result`, with `std` and `size`), the `mulliken`, `resp1`, `resp2` and
`am1bcc` charges, the dual descriptor (`dual`) and the `mbis_dipole_strength`.
Each of the six other properties has a `_std` and a `_count` column. The
populator used the default
tie-breaking rule, so read these trees with the default `DASHParams`.
`dash_props_pruned_0.02` is `dash_props` pruned at 0.02 e on `result`, with
the means and standard deviations only.

A tree populated from 90 % of the molecules gives these RMS errors against
every conformer of the other 10 % (29,469 molecules, 1.23 million atoms):

| column | `dash_props` (774 MB) | `dash_props_pruned_0.02` (94 MB) |
|---|---|---|
| `result` (MBIS charge) | 0.0250 e | 0.0254 e |
| `result` after `GetPartialCharges()` | 0.0243 e | 0.0247 e |
| `mulliken` | 0.0374 e | 0.0377 e |
| `resp1` | 0.0676 e | 0.0673 e |
| `resp2` | 0.0699 e | 0.0683 e |
| `am1bcc` | 0.0278 e | 0.0285 e |
| `dual` | 0.0477 | 0.0479 |
| `mbis_dipole_strength` | 0.0106 | 0.0107 |

The conformers of a molecule scatter by 0.012 e RMS around their mean MBIS
charge. No tree can be more accurate than that, because a tree reads the
molecular graph only. The two trees carry the licence of the data set, CC BY-SA
4.0. [github.com/MTLehner/dash-trees](https://github.com/MTLehner/dash-trees)
says how they were made.

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
props = GetDASHTree("dash_props_pruned_0.02", ["am1bcc", "am1bcc_std", "dual"])
dual = props.GetMolProperty(mol, "dual")            # one value per atom, no normalization

# any per-atom property with a molecule-wide total normalizes like the charges do
options = rdDASHTree.NormalizationOptions()
options.valueProperty, options.stdProperty = "am1bcc", "am1bcc_std"
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

The rule decides which nodes collect the values when you populate a tree. So
populate and read a tree with the same rule. With the DASH-properties data,
the default rule gives a 0.3 % smaller MBIS charge error than the legacy rule.
The difference is small but statistically significant. For the six other
properties, the default rule is as accurate or more accurate. If you populate a
tree with one rule and read it with the other, its MBIS charge error is 8 % to
42 % larger.
