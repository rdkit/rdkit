"""Module containing a C++ implementation of the FMCS algorithm"""

from collections.abc import Sequence
import enum
from typing import overload

import rdkit.Chem.rdchem


class MCSResult:
    """used to return MCS results"""

    @property
    def numAtoms(self) -> int:
        """number of atoms in MCS"""

    @property
    def numBonds(self) -> int:
        """number of bonds in MCS"""

    @property
    def queryMol(self) -> rdkit.Chem.rdchem.Mol | None:
        """query molecule for the MCS"""

    @property
    def smartsString(self) -> str:
        """SMARTS string for the MCS"""

    @property
    def canceled(self) -> bool:
        """if True, the MCS calculation did not finish"""

    @property
    def degenerateSmartsQueryMolDict(self) -> dict:
        """
        Dictionary collecting all degenerate (SMARTS, queryMol) pairs (empty if MCSParameters.StoreAll is False)
        """

class AtomCompare(enum.Enum):
    CompareAny = 0

    CompareElements = 1

    CompareIsotopes = 2

    CompareAnyHeavyAtom = 3

class BondCompare(enum.Enum):
    CompareAny = 0

    CompareOrder = 1

    CompareOrderExact = 2

class RingCompare(enum.Enum):
    IgnoreRingFusion = 0

    PermissiveRingFusion = 1

    StrictRingFusion = 2

@overload
def FindMCS(mols: Sequence[rdkit.Chem.rdchem.Mol], maximizeBonds: bool = True, threshold: float = 1.0, timeout: int = 3600, verbose: bool = False, matchValences: bool = False, ringMatchesRingOnly: bool = False, completeRingsOnly: bool = False, matchChiralTag: bool = False, atomCompare: AtomCompare = AtomCompare.CompareElements, bondCompare: BondCompare = BondCompare.CompareOrder, ringCompare: RingCompare = RingCompare.IgnoreRingFusion, seedSmarts: str = '') -> MCSResult: ...

@overload
def FindMCS(mols: Sequence[rdkit.Chem.rdchem.Mol], parameters: MCSParameters) -> MCSResult:
    """Find the MCS for a set of molecules"""

class MCSParameters:
    """Parameters controlling how the MCS is constructed"""

    def __init__(self) -> None: ...

    @property
    def MaximizeBonds(self) -> bool:
        """
        toggles maximizing the number of bonds (instead of the number of atoms)
        """

    @MaximizeBonds.setter
    def MaximizeBonds(self, arg: bool, /) -> None: ...

    @property
    def Threshold(self) -> float:
        """fraction of the dataset that must contain the MCS"""

    @Threshold.setter
    def Threshold(self, arg: float, /) -> None: ...

    @property
    def Timeout(self) -> int:
        """timeout (in seconds) for the calculation"""

    @Timeout.setter
    def Timeout(self, arg: int, /) -> None: ...

    @property
    def Verbose(self) -> bool:
        """toggles verbose mode"""

    @Verbose.setter
    def Verbose(self, arg: bool, /) -> None: ...

    @property
    def AtomCompareParameters(self) -> MCSAtomCompareParameters:
        """parameters for comparing atoms"""

    @AtomCompareParameters.setter
    def AtomCompareParameters(self, arg: MCSAtomCompareParameters, /) -> None: ...

    @property
    def BondCompareParameters(self) -> MCSBondCompareParameters:
        """parameters for comparing bonds"""

    @BondCompareParameters.setter
    def BondCompareParameters(self, arg: MCSBondCompareParameters, /) -> None: ...

    @property
    def AtomTyper(self) -> object:
        """
        atom typer to be used. Must be one of the
        members of the rdFMCS.AtomCompare class or
        an instance of a user-defined subclass of
        rdFMCS.MCSAtomCompare
        """

    @AtomTyper.setter
    def AtomTyper(self, arg: object, /) -> None: ...

    @property
    def BondTyper(self) -> object:
        """
        bond typer to be used. Must be one of the
        members of the rdFMCS.BondCompare class or
        an instance of a user-defined subclass of
        rdFMCS.MCSBondCompare
        """

    @BondTyper.setter
    def BondTyper(self, arg: object, /) -> None: ...

    @property
    def ProgressCallback(self) -> object:
        """
        progress callback class. Must be a
        user-defined subclass of rdFMCS.MCSProgress
        """

    @ProgressCallback.setter
    def ProgressCallback(self, arg: object, /) -> None: ...

    @property
    def FinalMatchChecker(self) -> object:
        """
        seed final match checker callback class. Must be a
        user-defined subclass of rdFMCS.MCSFinalMatchCheck
        """

    @FinalMatchChecker.setter
    def FinalMatchChecker(self, arg: object, /) -> None: ...

    @property
    def ShouldAcceptMCS(self) -> object:
        """
        MCS acceptance callback class. Must be a
        user-defined subclass of rdFMCS.MCSAcceptance
        """

    @ShouldAcceptMCS.setter
    def ShouldAcceptMCS(self, arg: object, /) -> None: ...

    @property
    def InitialSeed(self) -> str:
        """SMILES string to be used as the seed of the MCS"""

    @InitialSeed.setter
    def InitialSeed(self, arg: str, /) -> None: ...

    @property
    def StoreAll(self) -> bool:
        """toggles storage of degenerate MCSs"""

    @StoreAll.setter
    def StoreAll(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class MCSAtomCompareParameters:
    """Parameters controlling how atom-atom matching is done"""

    def __init__(self) -> None: ...

    @property
    def MatchValences(self) -> bool:
        """include atom valences in the match"""

    @MatchValences.setter
    def MatchValences(self, arg: bool, /) -> None: ...

    @property
    def MatchChiralTag(self) -> bool:
        """include atom chirality in the match"""

    @MatchChiralTag.setter
    def MatchChiralTag(self, arg: bool, /) -> None: ...

    @property
    def MaxDistance(self) -> float:
        """Require atoms to be within this many angstroms in 3D"""

    @MaxDistance.setter
    def MaxDistance(self, arg: float, /) -> None: ...

    @property
    def MatchFormalCharge(self) -> bool:
        """include formal charge in the match"""

    @MatchFormalCharge.setter
    def MatchFormalCharge(self, arg: bool, /) -> None: ...

    @property
    def RingMatchesRingOnly(self) -> bool:
        """ring atoms are only allowed to match other ring atoms"""

    @RingMatchesRingOnly.setter
    def RingMatchesRingOnly(self, arg: bool, /) -> None: ...

    @property
    def CompleteRingsOnly(self) -> bool:
        """results cannot include lone ring atoms"""

    @CompleteRingsOnly.setter
    def CompleteRingsOnly(self, arg: bool, /) -> None: ...

    @property
    def MatchIsotope(self) -> bool:
        """use isotope atom queries in MCSResults"""

    @MatchIsotope.setter
    def MatchIsotope(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class MCSBondCompareParameters:
    """Parameters controlling how bond-bond matching is done"""

    def __init__(self) -> None: ...

    @property
    def RingMatchesRingOnly(self) -> bool:
        """ring bonds are only allowed to match other ring bonds"""

    @RingMatchesRingOnly.setter
    def RingMatchesRingOnly(self, arg: bool, /) -> None: ...

    @property
    def CompleteRingsOnly(self) -> bool:
        """results cannot include partial rings"""

    @CompleteRingsOnly.setter
    def CompleteRingsOnly(self, arg: bool, /) -> None: ...

    @property
    def MatchFusedRings(self) -> bool:
        """
        enforce check on ring fusion, i.e. alpha-methylnaphthalene
        won't match beta-methylnaphtalene, but decalin
        will match cyclodecane unless MatchFusedRingsStrict is True
        """

    @MatchFusedRings.setter
    def MatchFusedRings(self, arg: bool, /) -> None: ...

    @property
    def MatchFusedRingsStrict(self) -> bool:
        """
        only enforced if MatchFusedRings is True; the ring fusion
        must be the same in both query and target, i.e. decalin
        won't match cyclodecane
        """

    @MatchFusedRingsStrict.setter
    def MatchFusedRingsStrict(self, arg: bool, /) -> None: ...

    @property
    def MatchStereo(self) -> bool:
        """include bond stereo in the comparison"""

    @MatchStereo.setter
    def MatchStereo(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class MCSProgressData:
    """Information about the MCS progress"""

    def __init__(self) -> None: ...

    @property
    def numAtoms(self) -> int:
        """number of atoms in MCS"""

    @property
    def numBonds(self) -> int:
        """number of bonds in MCS"""

    @property
    def seedProcessed(self) -> int:
        """number of processed seeds"""

class MCSAtomCompare:
    """
    Base class. Subclass and override
    MCSAtomCompare.__call__() to define custom
    atom compare functions, then set MCSParameters.AtomTyper
    to an instance of the subclass
    """

    def __init__(self) -> None: ...

    def CheckAtomRingMatch(self, parameters: MCSAtomCompareParameters, mol1: rdkit.Chem.rdchem.Mol, atom1: int, mol2: rdkit.Chem.rdchem.Mol, atom2: int) -> bool:
        """Return True if both atoms are, or are not, in a ring"""

    def CheckAtomCharge(self, parameters: MCSAtomCompareParameters, mol1: rdkit.Chem.rdchem.Mol, atom1: int, mol2: rdkit.Chem.rdchem.Mol, atom2: int) -> bool:
        """Return True if both atoms have the same formal charge"""

    def CheckAtomChirality(self, parameters: MCSAtomCompareParameters, mol1: rdkit.Chem.rdchem.Mol, atom1: int, mol2: rdkit.Chem.rdchem.Mol, atom2: int) -> bool:
        """Return True if both atoms have, or have not, a chiral tag"""

    def __call__(self, parameters: MCSAtomCompareParameters, mol1: rdkit.Chem.rdchem.Mol, atom1: int, mol2: rdkit.Chem.rdchem.Mol, atom2: int) -> bool:
        """override to implement custom atom comparison"""

class MCSBondCompare:
    """
    Base class. Subclass and override
    MCSBondCompare.__call__() to define custom
    bond compare functions, then set MCSParameters.BondTyper
    to an instance of the subclass
    """

    def __init__(self) -> None: ...

    def CheckBondStereo(self, parameters: MCSBondCompareParameters, mol1: rdkit.Chem.rdchem.Mol, bond1: int, mol2: rdkit.Chem.rdchem.Mol, bond2: int) -> bool:
        """Return True if both bonds have, or have not, a stereo descriptor"""

    def CheckBondRingMatch(self, parameters: MCSBondCompareParameters, mol1: rdkit.Chem.rdchem.Mol, bond1: int, mol2: rdkit.Chem.rdchem.Mol, bond2: int) -> bool:
        """Return True if both bonds are, or are not, part of a ring"""

    def __call__(self, parameters: MCSBondCompareParameters, mol1: rdkit.Chem.rdchem.Mol, bond1: int, mol2: rdkit.Chem.rdchem.Mol, bond2: int) -> bool:
        """override to implement custom bond comparison"""

class MCSProgress:
    """
    Base class. Subclass and override
    MCSProgress.__call__()
    to define a custom callback function
    """

    def __init__(self) -> None: ...

    def __call__(self, stat: MCSProgressData, parameters: MCSParameters) -> bool:
        """override to implement a custom progress callback"""

class MCSFinalMatchCheck:
    """
    Base class. Subclass and override
    MCSFinalMatchCheck.__call__()
    to define a custom boolean callback function.
    Returning True will cause the growing seed to be accepted,
    False to be rejected
    """

    def __init__(self) -> None: ...

    def __call__(self) -> bool:
        """override to implement a custom seed final match checker callback"""

class MCSAcceptance:
    """
    Base class. Subclass and override
    MCSAcceptance.__call__()
    to define a custom boolean callback function.
    Returning True will cause the MCS candidate to be accepted,
    False to be rejected
    """

    def __init__(self) -> None: ...

    def __call__(self) -> bool:
        """override to implement a custom MCS acceptance callback"""
