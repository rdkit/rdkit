"""Module containing functions for creating a Scaffold Network"""

from collections.abc import Sequence
import enum
from typing import overload

import rdkit.Chem.rdchem


class ScaffoldNetworkParams:
    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, bondBreakerSmartsList: Sequence[str]) -> None:
        """
        Constructor taking a list of Reaction SMARTS for the fragmentation reactions
        """

    @property
    def includeGenericScaffolds(self) -> bool:
        """include scaffolds with all atoms replaced by dummies"""

    @includeGenericScaffolds.setter
    def includeGenericScaffolds(self, arg: bool, /) -> None: ...

    @property
    def includeGenericBondScaffolds(self) -> bool:
        """include scaffolds with all bonds replaced by single bonds"""

    @includeGenericBondScaffolds.setter
    def includeGenericBondScaffolds(self, arg: bool, /) -> None: ...

    @property
    def includeScaffoldsWithoutAttachments(self) -> bool:
        """remove attachment points from scaffolds and include the result"""

    @includeScaffoldsWithoutAttachments.setter
    def includeScaffoldsWithoutAttachments(self, arg: bool, /) -> None: ...

    @property
    def includeScaffoldsWithAttachments(self) -> bool:
        """Include the version of the scaffold with attachment points"""

    @includeScaffoldsWithAttachments.setter
    def includeScaffoldsWithAttachments(self, arg: bool, /) -> None: ...

    @property
    def includeNames(self) -> bool:
        """Include molecules names of the input molecules"""

    @includeNames.setter
    def includeNames(self, arg: bool, /) -> None: ...

    @property
    def keepOnlyFirstFragment(self) -> bool:
        """keep only the first fragment from the bond breaking rule"""

    @keepOnlyFirstFragment.setter
    def keepOnlyFirstFragment(self, arg: bool, /) -> None: ...

    @property
    def pruneBeforeFragmenting(self) -> bool:
        """Do a pruning/flattening step before starting fragmenting"""

    @pruneBeforeFragmenting.setter
    def pruneBeforeFragmenting(self, arg: bool, /) -> None: ...

    @property
    def flattenIsotopes(self) -> bool:
        """remove isotopes when flattening"""

    @flattenIsotopes.setter
    def flattenIsotopes(self, arg: bool, /) -> None: ...

    @property
    def flattenChirality(self) -> bool:
        """remove chirality and bond stereo when flattening"""

    @flattenChirality.setter
    def flattenChirality(self, arg: bool, /) -> None: ...

    @property
    def flattenKeepLargest(self) -> bool:
        """keep only the largest fragment when doing flattening"""

    @flattenKeepLargest.setter
    def flattenKeepLargest(self, arg: bool, /) -> None: ...

    @property
    def collectMolCounts(self) -> bool:
        """keep track of the number of molecules each scaffold was found in"""

    @collectMolCounts.setter
    def collectMolCounts(self, arg: bool, /) -> None: ...

class EdgeType(enum.Enum):
    Fragment = 1

    Generic = 2

    GenericBond = 3

    RemoveAttachment = 4

    Initialize = 5

class NetworkEdge:
    @property
    def beginIdx(self) -> int:
        """index of the begin node in node list"""

    @property
    def endIdx(self) -> int:
        """index of the end node in node list"""

    @property
    def type(self) -> EdgeType:
        """type of the edge"""

    def __str__(self) -> str: ...

class ScaffoldNetwork:
    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, pkl: bytes) -> None: ...

    @overload
    def __init__(self, pkl: str) -> None: ...

    @property
    def nodes(self) -> list[str]:
        """the sequence of SMILES defining the nodes"""

    @property
    def counts(self) -> list[int]:
        """
        the number of times each node was encountered while building the network.
        """

    @property
    def molCounts(self) -> list[int]:
        """the number of moleclues each node was found in."""

    @property
    def edges(self) -> list[NetworkEdge]:
        """the sequence of network edges"""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

def CreateScaffoldNetwork(mols: Sequence[rdkit.Chem.rdchem.Mol], params: ScaffoldNetworkParams) -> ScaffoldNetwork:
    """create (and return) a new network from a sequence of molecules"""

def UpdateScaffoldNetwork(mols: Sequence[rdkit.Chem.rdchem.Mol], network: ScaffoldNetwork, params: ScaffoldNetworkParams) -> None:
    """update an existing network by adding molecules"""

def BRICSScaffoldParams() -> ScaffoldNetworkParams:
    """
    Returns parameters for generating scaffolds using BRICS fragmentation rules
    """
