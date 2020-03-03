# Enhanced Stereochemistry in the RDKit
Greg Landrum (greg.landrum@t5informatics.com)

September 2018

*This is still a DRAFT*

## Overview

We are going to follow, at least for the initial implementation, the enhanced stereo representation used in V3k mol files: groups of atoms with specified stereochemistry with an `ABS`, `AND`, or `OR` marker indicating what is known. The general idea is that `AND` indicates mixtures and `OR` indicates unknown single substances.

Here are some illustrations of what the various combinations mean:

|   What's drawn  | Mixture? | What it means |
|-----------------|-----------------------------|:--------------|
|![img1a](images/enhanced_stereo_and1_and2_base.png) | mixture | ![img1b](images/enhanced_stereo_and1_and2_expand.png) |
|![img2a](images/enhanced_stereo_and1_cis_base.png) | mixture | ![img2b](images/enhanced_stereo_and1_cis_expand.png) |
|![img3a](images/enhanced_stereo_and1_trans_base.png) | mixture | ![img3b](images/enhanced_stereo_and1_trans_expand.png)|
|![img4a](images/enhanced_stereo_or1_or2_base.png) | single | ![img4b](images/enhanced_stereo_and1_and2_expand.png) |
|![img5a](images/enhanced_stereo_or1_cis_base.png) | single | ![img5b](images/enhanced_stereo_and1_cis_expand.png) |
|![img5a](images/enhanced_stereo_or1_trans_base.png) | single | ![img5b](images/enhanced_stereo_and1_trans_expand.png)|
|![img6a](images/enhanced_stereo_abs_and_base.png) | mixture | ![img6b](images/enhanced_stereo_abs_and_expand.png)|
|![img7a](images/enhanced_stereo_abs_or_base.png) | single | ![img7b](images/enhanced_stereo_abs_and_expand.png)|


## Use cases

The initial target is to not lose data on an V3k mol -> RDKit -> V3k mol round trip. Manipulation,
depiction, and searching is a secondary goal.

## Representation

Stored as a vector of `StereoGroup` objects.

A `StereoGroup` contains an enum with the type as well as pointers to the atoms involved. We will need to adjust this when atoms are removed or replaced.

## Enumeration

The Python-only function `rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers` can be used to enumerate
all stereoisomers represented by the StereoGroups. It is possible to enumerate _only_ the Stereo Groups,
ignoring other unspecified centers. 

## Reactions

The reaction runner is aware of stereo groups associated with a Mol. Atoms in a Stereo Group are copied to
all products in which they appear, unless the reaction created or destroyed stereochemistry at that atom.
If an atom from a reactant StereoGroup appears multiple times in the product, all copies of that atom are
put in the same product StereoGroup.

## Searching

This will be handled by searching in MolBundle objects produced by the enumeration code. Substructure searching and canonicalization do not yet support stereo groups.

## Depiction

Something needs to be added to the depiction code to allow the groups to be seen.
