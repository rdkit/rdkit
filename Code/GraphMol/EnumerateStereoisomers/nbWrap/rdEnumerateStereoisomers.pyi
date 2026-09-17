"""
Module containing functions to enumerate stereoisomers of a molecule.
Chiral centers and double bonds will be enumerated if unassigned, or,
if the appropriate option is set, if assigned.  Atropisomers will only
be enumerated if assigned.  There is, as yet, no means of finding
unassigned atropisomers.
"""

from typing import overload

import rdkit.Chem.rdchem


class StereoEnumerationOptions:
    """EnumerateSteroisomers options."""

    def __init__(self) -> None: ...

    @property
    def tryEmbedding(self) -> bool:
        """
        If true, the process attempts to generate a standard RDKit distance geometry
        conformation for the stereoisomer.  If this fails, we assume that the stereoisomer is
        non-physical and don't return it.  NOTE that this is computationally expensive and is
        just a heuristic that could result in stereoisomers being lost.  Default=False
        """

    @tryEmbedding.setter
    def tryEmbedding(self, arg: bool, /) -> None: ...

    @property
    def onlyUnassigned(self) -> bool:
        """
        If true, stereocenters which have a specified stereochemistry will not be
        perturbed unless they are part of a relative stereo group.  Default=True.
        """

    @onlyUnassigned.setter
    def onlyUnassigned(self, arg: bool, /) -> None: ...

    @property
    def onlyStereoGroups(self) -> bool:
        """
        If true, only find stereoisomers that differ at the StereoGroups associated with
        the molecule.  Default=False.
        """

    @onlyStereoGroups.setter
    def onlyStereoGroups(self, arg: bool, /) -> None: ...

    @property
    def unique(self) -> bool:
        """
        If true, only stereoisomers that differ in canonical CXSmiles will be
        returned.  Default=True.
        """

    @unique.setter
    def unique(self, arg: bool, /) -> None: ...

    @property
    def maxIsomers(self) -> int:
        """
        The maximum number of isomers to yield.  If the number of possible isomers
        is greater than maxIsomers, a random subset will be yielded.  If 0, there
        is no maximum.  Since every additional stereocenter doubles the number of
        results (and execution time) it's important to keep an eye on this.
        """

    @maxIsomers.setter
    def maxIsomers(self, arg: int, /) -> None: ...

    @property
    def randomSeed(self) -> int:
        """Seed for random number generator.  Default=-1 means no seed."""

    @randomSeed.setter
    def randomSeed(self, arg: int, /) -> None: ...

class StereoisomerEnumerator:
    """Stereoisomer enumerator."""

    @overload
    def __init__(self, mol: rdkit.Chem.rdchem.Mol, verbose: bool = False) -> None: ...

    @overload
    def __init__(self, mol: rdkit.Chem.rdchem.Mol, options: StereoEnumerationOptions, verbose: bool = False) -> None: ...

    def next(self) -> rdkit.Chem.rdchem.Mol | None:
        """Get next isomer in the sequence, or None if at the end."""

    def GetStereoisomerCount(self) -> int:
        """Get the number of stereoisomers."""
