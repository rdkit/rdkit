# Backwards incompatible changes

Though we really try to maintain backwards compatibility across releases, there
are rare instances when we need to break compatibility. The purpose of this
document is to capture information about backwards incompatible changes that
have been made in the RDKit. It's broken down by release cycle and does *not*
generally include changes in results that arise due to bug fixes; we try to call
those out in the release notes. The release notes, in general, contain a more
comprehensive version of this list.

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
