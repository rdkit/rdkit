# SGroups in RDKit
Ricardo Rodriguez-Schmidt

September 2018


## Overview

MDL SD files support SGroups, "Substance groups" in MDL terminology, as a way to store some types of extended intra-molecular information. Some usage examples of these SGroups are the following:

1. Identification of repeating groups in polymers, i.e. identification of monomer units and, in case of copolymers, how these monomers distribute in the polymer chain.
1. Labeling of relevant sections of a molecule.
1. Detailing depiction information for the molecule or other defined SGroups.
1. Storing information relative to parts of the molecule in the form of data fields.

This document describes an intent to make a `SubstanceGroup` data structure available to RDKit molecules, initiating support for parsing, storing, serializing and writing back this information. The goal is to achieve conversion of MDL SD files that contain SGroups using RDkit without loss of information in these, with some limitations. At this point, no intent will be made to allow or support manipulation or depiction of the information contained in SGroups, but this might be added in the future.

## Documentation

As a reference for the design and implementation of the `SubstanceGroup` data class, the following document has been used:

[http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf](http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf)

The document describes the syntax for specifying SGroup information both in version V2000 and V3000 of MDL SD files. Both of these have been implemented for input and output, as well as interconversion between them. Best effort has been made to support as many features of SGroups as possible, as well as following specifications as close as possible. Nevertheless, the following features are not currently included in this implementation:

- V2000 "CRS" lines.
- V3000 "XBHEAD" and "XBCORR" labels.

## Testing

A Test file has been added to the SGroup code to make sure that it works in the expected way. Nonetheless, the sample SGroups which are checked in these tests have been designed to check as many features as possible with as less SGroups as possible. They do not make any sense, and should not be taken as a reference on how to use SGroups.

Besides this test file, the SGroups code has been tested on some sample molecules, which have been parsed and then written back to disk to check data consistency.

**Remarks**:

1. The `SubstanceGroup` class contains elements consisting of pointers to atoms, bonds and other SGroup objects. When making a copy of a molecule that has SubstanceGroups, it is required to update these pointers to objects contained in the new owning molecule. This is done by using the private method `void updateOwningMol(ROMol *other_mol);`.

   This method is to be used **only** after all atoms, bonds, and SubstanceGroups on the new molecule have been created, as the pointers will be updated with the objects that are at the same indexes as those in the original molecule. This is especially important for SubstanceGroups, as these can have a parent relationship pointing to another SubstanceGroup with a higher index (i.e. that has been added to the molecule after the current one), and therefore, might not be initialized yet).

2. This implementation of SubstanceGroups is only meant for parsing and writing from/to SD files and serialization, and not for manipulation. Please take into account that:
    - SubstanceGroups cannot be removed individually (but all of them can be cleared from a molecule).
    - Adding, replacing or removing atoms/bonds from a RWMol will drop the SubstanceGroups.