# Backwards incompatible changes

Though we really try to maintain backwards compatibility across releases, there
are rare instances when we need to break compatibility. The purpose of this
document is to capture information about backwards incompatible changes that
have been made in the RDKit. It's broken down by release cycle and does *not*
generally include changes in results that arise due to bug fixes; we try to call
those out in the release notes. The release notes, in general, contain a more
comprehensive version of this list.

## Release 2022.09

### Changes in canonicalization

Changes to the way atomic chirality is used in the canonicalization algorithm
mean that canonical atom ranking and canonical SMILES generated with this RDKit
version can be different from those generated with previous versions

### Improvements to GetBestRMS() and CalcRMS()

`GetBestRMS() and CalcRMS()` by default now treat terminal conjugated functional
groups like carboxylate and nitro symmetrically. For example, the group
`C(=[O:1])[O-:2]` can match in either orientation. The SMARTS pattern which is
used to recognize affected groups is:
`[{atomP};$([{atomP}]-[*]=[{atomP}]),$([{atomP}]=[*]-[{atomP}])]~[*]` where
`{atomP}` is `O,N;D1`. The previous behavior can be restored using by setting
the `symmetrizeConjugatedTerminalGroups` argument to false when calling
`GetBestRMS() and CalcRMS()`


## Release 2022.03

### Changes to the handling of dummy atoms and aromaticity

The rules for aromaticity in rings containing dummy atoms have been changed. The
general intention of the new handling is that aromaticity will not be perceived
for rings containing dummy atoms unless it's clear that the dummies should be
aromatic. As an example: the SMILES `C1=C*2=CC=CC=*2C=C1` is perceived to be
aromatic while the SMILES `C1=C*2C=CC=C*2C=C1` (which does not have any double
bonds to the dummy atoms) is not; in previous RDKit releases both of these
structures were aromatic. There's more information about this in the discussion
of PR #4722 (https://github.com/rdkit/rdkit/pull/4722) and Issue #4721
(https://github.com/rdkit/rdkit/issues/4721).

### Tautomer enumeration improvements

The rules for tautomer enumeration in `MolStandardize` have been updated to more
closely match the rules in the original publication. These changes primarily
consist of making the rules more specific; the consequence is that less
tautomers will be generated with this version. The previous rules can still be
accessed via the function `GetV1TautomerEnumerator()` (Python) or
`getV1TautomerEnumerator()` (C++)


## Release 2021.03

### Changes to conformer generation

In the interest of providing higher quality results by default, the conformer
generator now uses symmetry when doing RMS pruning. The pure distance geometry
based conformer generator also now generates trans(oid) conformations for
amides, esters, and related structures.

### Changes to R-Group Decomposition results

There have been a large number of changes in the R-group decomposition code
which we believe improve the quality of the results and the handling of the many
strange edge cases which arise when doing RGD.

### Reactions now produce STEREONONE instead of STEREOANY for unspecified double bonds

Double bonds in the products of reactions which have unspecified stereochemistry
are now set to be `STEREONONE` instead of `STEREOANY`. Using `STEREONONE` is
consistent with the rest of the behavior of the RDKit with respect to double
bonds.

### Moving MolToSVG()

The function `Chem.MolToSVG()` (or `Chem.rdmolops.MolToSVG()`) has been moved to
the `Chem.Draw` package along with the other drawing code.


## Release 2020.09

### Changes in the chemistry model [UPDATED for release 2020.09.2]

We've added additional allowed valences for I (now 1, 3, 5), At (now 1, 3, 5),
Xe (now 0, 2, 4, 6), and Po (now 2, 4, 6). Molecules with atoms in the new
valence states will no longer generate sanitization errors. Note that this has
an impact on the chemistry of molecules containing 3-valent I and at least one
implicit H (present 24 times in ChEMBL 27): previously this was incorrectly
assigned two implicit Hs, now it has no implicit Hs. 

Aromaticity perception of molecules like `Cc1nnc2n1c1ccccc1n1c(C)nnc12` now
correctly recognizes the full outer envelope, i.e. the bonds joining the rings
are now also aromatic. In 2020.03 the three rings were connected by single bonds.

### FindMCS can now return single atom MCSs
FindMCS() may return single atom MCSs, whereas previously it returned an empty
MCS unless there was at least one commond bond across the input structures.
So the MCS between molecules `CC` and `CO` is now `[#6]` rather than being null.

### New return type for tautomer enumeration
Due to improvements in the tautomer enumeration code, the method
`TautomerEnumerator::enumerate` now returns a `TautomerEnumeratorResult`
object instead of a vector of molecules. Note that if you are iterating over
the results of a call to `enumerate()` you shouldn't need to change your code.
If you want to invoke the old (and deprecated, see below) form from C++, call
`TautomerNumerator::enumerate(mol, nullptr)` or explicitly pass a
`boost::dynamic_bitset*` to capture the modified atoms.

### Default precision for coordgen
The default precision setting for coordgen has been changed. The new default
was selected to greatly reduce the number of molecules for which it takes a
very long time to generate coordinates while still producing nice looking
structures. We may continue to tweak this default value if/when problems
with it are reported. If you would like to go back to the previous setting, set 
`CoordgenParams.minimizerPrecision` to `CoordgenParams.sketcherStandardPrecision` 
when you invoke `rdCoordGen.AddCoords()`

### Change in uncharge() behavior for inorganic acids
Uncharger::uncharge() will now neutralize `[Cl,Br,I][O-], [Cl,Br,I](=O)[O-],
[Cl,Br,I](=O)(=O)[O-], [Cl,Br,I](=O)(=O)(=O)[O-], [O-]N=N[O-], [N,P](=O)[O-],
[N+](=O)([O-])[O-], P(=O)([O-])[O-], P(=O)([O-])([O-])[O-], S([O-])[O-],
S(=O)([O-])[O-], S(=O)(=O)([O-])[O-], S(=O)(=O)([O-])OOS(=O)(=O)[O-]`.
Previously not all of these inorganic acid counterions were consistently
neutralized.

## Release 2020.03

### Use of chirality with `@=` in PostgreSQL cartridge
Searches for equal molecules (i.e. `mol1 @= mol2`) in the PostgreSQL cartridge
now use the `do_chiral_sss` option. So if `do_chiral_sss` is false (the
default), the molecules `CC(F)Cl` and `C[C@H](F)Cl` will be considered to be
equal. Previously these molecules were always considered to be different.

### Handling of empty/invalid files in MolSuppliers
Attempting to create a MolSupplier from a filename pointing to an empty file, a
file that does not exist or sometihing that is not a standard file (i.e.
something like a directory) now generates an exception.

### Renaming a cmake option
The cmake option `RDK_OPTIMIZE_NATIVE` has been renamed to `RDK_OPTIMIZE_POPCNT`



## Release 2019.09

### Updated van der Waals radii
The atomic van der Waals radii used by the RDKit were corrected/updated in #2154.
This leads to different results when generating conformations, molecular volumes,
and molecular shapes.

### "R" atoms in mol blocks are now converted into * queries
When reading Mol blocks, atoms with the symbol "R" are now converted into
queries that match any atom when doing a substructure search (analogous to "*"
in SMARTS). The previous behavior was to only match other dummy atoms.

### PandasTools.LoadSDF() uses isomeric SMILES
When loading SDF files using PandasTools.LoadSDF(), we now default to
producing isomeric smiles in pandas tables.  To reproduce the original
behavior, use isomericSmiles=False in the call to the function.

### New fingerprinter defaults changed
The default values for some of the new fingerprint generators have been changed so
that they more closely resemble the original fingerprinting code. In
particular most fingerprinters no longer do count simulation by default and
the RDKit fingerprint now sets two bits per feature by default.


## Release 2019.03
As of this realease (2019.03.1) the RDKit no longer supports Python 2.

## Release 2018.09

### `GetAtomSmiles()` generates isomeric SMILES by default
In earlier releases, if you wanted to include information about stereochemistry
or isotopic labels in the output SMILES it was necessary to set the optional
`isomericSmiles` argument to true. The default value of this is now true. If you
want to go back to the old behavior and get SMILES without stereochemistry
information you can just set `isomericSmiles` to false.

### Changes to ringMatchesRingOnly option in the MCS code
The ringMatchesRingOnly option to the FindMCS() function now applies to
atom-atom matches as well as bond-bond matches.

### The conformation generation code now uses ETKDG by default when called from Python
The Python functions EmbedMolecule() and EmbedMultipleConfs() now use the ETKDG algorithm by default instead of standard distance geometry.

## Release 2018.03

### `MolToSmiles()` generates isomeric SMILES by default
In earlier releases, if you wanted to include information about stereochemistry
or isotopic labels in the output SMILES it was necessary to set the optional
`isomericSmiles` argument to true. The default value of this is now true. If you
want to go back to the old behavior and get SMILES without stereochemistry
information you can just set `isomericSmiles` to false.

### `MolToMolBlock()` generates a 2D conformation when the `includeStereo` flag is set
If you want to capture stereochemistry in Mol blocks, it's necessary to have
coordinates in the output. Earlier versions of the RDKit required you to
remember to generate coordinates yourself. This is now done by default for
molecules that don't have a conformation when `includeStereo` is set.

### The conformation generation code now uses ETKDG by default
Earlier RDKit releases would, by default, generate conformations using standard
distance geometry. The new default is to use Sereina Riniker's ETKDG algorithm,
which is a bit slower but which has been shown to generate much better results.


## License

This document is copyright (C) 2017-2019 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
License. To view a copy of this license, visit
<http://creativecommons.org/licenses/by-sa/4.0/> or send a letter to Creative
Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.

The intent of this license is similar to that of the RDKit itself. In simple
words: “Do whatever you want with it, but please give us some credit.”
