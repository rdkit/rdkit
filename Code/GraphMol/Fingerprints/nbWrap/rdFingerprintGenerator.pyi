from collections.abc import Iterable
import enum
from typing import Annotated

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdchem
import rdkit.DataStructs.cDataStructs


class AtomInvariantsGenerator:
    pass

class BondInvariantsGenerator:
    pass

class AdditionalOutput:
    def __init__(self) -> None: ...

    def AllocateAtomToBits(self) -> None:
        """synonym for CollectAtomToBits()"""

    def AllocateBitInfoMap(self) -> None:
        """synonym for CollectBitInfoMap()"""

    def AllocateBitPaths(self) -> None:
        """synonym for CollectBitPaths()"""

    def AllocateAtomCounts(self) -> None:
        """synonym for CollectAtomCounts()"""

    def AllocateAtomsPerBit(self) -> None:
        """synonym for CollectAtomsPerBit()"""

    def CollectAtomToBits(self) -> None:
        """
        toggle collection of information mapping each atom to the bits it is involved in.
        """

    def CollectBitInfoMap(self) -> None:
        """
        toggles collection of information mapping each atom to more detail about the atom environment (not available from all fingerprints)
        """

    def CollectBitPaths(self) -> None:
        """
        toggles collection of information matching each atom to information about the paths it is involved in (not available from all fingerprints).
        """

    def CollectAtomCounts(self) -> None:
        """
        toggles collection of information about the number of bits each atom is involved in
        """

    def CollectAtomsPerBit(self) -> None:
        """
        toggles collection of information about all atoms involved in setting each bit
        """

    def GetAtomToBits(self) -> object: ...

    def GetBitInfoMap(self) -> object: ...

    def GetBitPaths(self) -> object: ...

    def GetAtomCounts(self) -> object: ...

    def GetAtomsPerBit(self) -> object: ...

class FingerprintOptions:
    @property
    def countSimulation(self) -> bool:
        """use count simulation"""

    @countSimulation.setter
    def countSimulation(self, arg: bool, /) -> None: ...

    @property
    def includeChirality(self) -> bool:
        """include chirality in atom invariants (not for all fingerprints)"""

    @includeChirality.setter
    def includeChirality(self, arg: bool, /) -> None: ...

    @property
    def fpSize(self) -> int:
        """size of the fingerprints created"""

    @fpSize.setter
    def fpSize(self, arg: int, /) -> None: ...

    @property
    def numBitsPerFeature(self) -> int:
        """number of bits to set for each feature"""

    @numBitsPerFeature.setter
    def numBitsPerFeature(self, arg: int, /) -> None: ...

    def SetCountBounds(self, bounds: Iterable[int]) -> None:
        """set the bins for the count bounds"""

class FingerprintGenerator32:
    def GetSparseCountFingerprint(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> rdkit.DataStructs.cDataStructs.UIntSparseIntVect:
        """
        Generates a sparse count fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseIntVect containing fingerprint
        """

    def GetSparseFingerprint(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> rdkit.DataStructs.cDataStructs.SparseBitVect:
        """
        Generates a sparse fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseBitVect containing fingerprint
        """

    def GetCountFingerprint(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> rdkit.DataStructs.cDataStructs.UIntSparseIntVect:
        """
        Generates a count fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseIntVect containing fingerprint
        """

    def GetCountFingerprintAsNumPy(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> Annotated[NDArray[numpy.uint32], dict(shape=(None,))]:
        """
        Generates a count fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a numpy array containing the fingerprint
        """

    def GetFingerprint(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
        """
        Generates a fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a ExplicitBitVect containing fingerprint
        """

    def GetFingerprintAsNumPy(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> Annotated[NDArray[numpy.uint8], dict(shape=(None,))]:
        """
        Generates a fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a numpy array containing the fingerprint
        """

    def GetFingerprints(self, mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int = 1) -> tuple[rdkit.DataStructs.cDataStructs.ExplicitBitVect, ...]:
        """
        Generates fingerprints for a sequence of molecules

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - numThreads: number of threads to use

        RETURNS: a tuple of ExplicitBitVects
        """

    def GetCountFingerprints(self, mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int = 1) -> tuple[rdkit.DataStructs.cDataStructs.UIntSparseIntVect, ...]:
        """
        Generates count fingerprints for a sequence of molecules

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - numThreads: number of threads to use

        RETURNS: a tuple of SparseIntVects
        """

    def GetSparseFingerprints(self, mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int = 1) -> tuple[rdkit.DataStructs.cDataStructs.SparseBitVect, ...]:
        """
        Generates sparse fingerprints for a sequence of molecules

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - numThreads: number of threads to use

        RETURNS: a tuple of SparseBitVects
        """

    def GetSparseCountFingerprints(self, mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int = 1) -> tuple[rdkit.DataStructs.cDataStructs.UIntSparseIntVect, ...]:
        """
        Generates sparse count fingerprints for a sequence of molecules

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - numThreads: number of threads to use

        RETURNS: a tuple of SparseIntVects
        """

    def GetInfoString(self) -> str:
        """
        Returns a string containing information about the fingerprint generator

        RETURNS: an information string
        """

    def GetOptions(self) -> FingerprintOptions:
        """return the fingerprint options object"""

    def ToJSON(self) -> str:
        """Serialize a FingerprintGenerator to JSON"""

class FingerprintGenerator64:
    def GetSparseCountFingerprint(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> rdkit.DataStructs.cDataStructs.ULongSparseIntVect:
        """
        Generates a sparse count fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseIntVect containing fingerprint
        """

    def GetSparseFingerprint(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> rdkit.DataStructs.cDataStructs.SparseBitVect:
        """
        Generates a sparse fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseBitVect containing fingerprint
        """

    def GetCountFingerprint(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> rdkit.DataStructs.cDataStructs.UIntSparseIntVect:
        """
        Generates a count fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseIntVect containing fingerprint
        """

    def GetCountFingerprintAsNumPy(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> Annotated[NDArray[numpy.uint32], dict(shape=(None,))]:
        """
        Generates a count fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a numpy array containing the fingerprint
        """

    def GetFingerprint(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
        """
        Generates a fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a ExplicitBitVect containing fingerprint
        """

    def GetFingerprintAsNumPy(self, mol: rdkit.Chem.rdchem.Mol, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, confId: int = -1, customAtomInvariants: Iterable[int] | None = None, customBondInvariants: Iterable[int] | None = None, additionalOutput: AdditionalOutput | None = None) -> Annotated[NDArray[numpy.uint8], dict(shape=(None,))]:
        """
        Generates a fingerprint

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - fromAtoms: only environments starting at or centered on these atoms will be included
            - ignoreAtoms: environments including these atoms will be excluded
            - confId: 3D confirmation to use, only used by AtomPair fingerprint
            - customAtomInvariants: custom atom invariants to be used,
              overrides invariants from the invariant generator
            - customBondInvariants: custom bond invariants to be used,
              overrides invariants from the invariant generator
            - additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a numpy array containing the fingerprint
        """

    def GetFingerprints(self, mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int = 1) -> tuple[rdkit.DataStructs.cDataStructs.ExplicitBitVect, ...]:
        """
        Generates fingerprints for a sequence of molecules

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - numThreads: number of threads to use

        RETURNS: a tuple of ExplicitBitVects
        """

    def GetCountFingerprints(self, mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int = 1) -> tuple[rdkit.DataStructs.cDataStructs.UIntSparseIntVect, ...]:
        """
        Generates count fingerprints for a sequence of molecules

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - numThreads: number of threads to use

        RETURNS: a tuple of SparseIntVects
        """

    def GetSparseFingerprints(self, mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int = 1) -> tuple[rdkit.DataStructs.cDataStructs.SparseBitVect, ...]:
        """
        Generates sparse fingerprints for a sequence of molecules

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - numThreads: number of threads to use

        RETURNS: a tuple of SparseBitVects
        """

    def GetSparseCountFingerprints(self, mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int = 1) -> tuple[rdkit.DataStructs.cDataStructs.ULongSparseIntVect, ...]:
        """
        Generates sparse count fingerprints for a sequence of molecules

        ARGUMENTS:
            - mol: molecule to be fingerprinted
            - numThreads: number of threads to use

        RETURNS: a tuple of SparseIntVects
        """

    def GetInfoString(self) -> str:
        """
        Returns a string containing information about the fingerprint generator

        RETURNS: an information string
        """

    def GetOptions(self) -> FingerprintOptions:
        """return the fingerprint options object"""

    def ToJSON(self) -> str:
        """Serialize a FingerprintGenerator to JSON"""

class FPType(enum.Enum):
    RDKitFP = 2

    MorganFP = 1

    AtomPairFP = 0

    TopologicalTorsionFP = 3

RDKitFP: FPType = FPType.RDKitFP

MorganFP: FPType = FPType.MorganFP

AtomPairFP: FPType = FPType.AtomPairFP

TopologicalTorsionFP: FPType = FPType.TopologicalTorsionFP

def GetSparseCountFPs(molecules: Iterable[rdkit.Chem.rdchem.Mol] | None = None, fpType: FPType = FPType.MorganFP) -> list[rdkit.DataStructs.cDataStructs.ULongSparseIntVect]: ...

def GetSparseFPs(molecules: Iterable[rdkit.Chem.rdchem.Mol] | None = None, fpType: FPType = FPType.MorganFP) -> list[rdkit.DataStructs.cDataStructs.SparseBitVect]: ...

def GetCountFPs(molecules: Iterable[rdkit.Chem.rdchem.Mol] | None = None, fpType: FPType = FPType.MorganFP) -> list[rdkit.DataStructs.cDataStructs.UIntSparseIntVect]: ...

def GetFPs(molecules: Iterable[rdkit.Chem.rdchem.Mol] | None = None, fpType: FPType = FPType.MorganFP) -> list[rdkit.DataStructs.cDataStructs.ExplicitBitVect]: ...

def FingerprintGeneratorFromJSON(jsonString: str) -> FingerprintGenerator64:
    """Deserialize a FingerprintGenerator from a JSON string"""

class AtomPairFingerprintOptions(FingerprintOptions):
    @property
    def use2D(self) -> bool:
        """use 2D distances"""

    @use2D.setter
    def use2D(self, arg: bool, /) -> None: ...

    @property
    def minDistance(self) -> int:
        """minimum distance to be included"""

    @minDistance.setter
    def minDistance(self, arg: int, /) -> None: ...

    @property
    def maxDistance(self) -> int:
        """maximum distance to be included"""

    @maxDistance.setter
    def maxDistance(self, arg: int, /) -> None: ...

def GetAtomPairGenerator(minDistance: int = 1, maxDistance: int = 30, includeChirality: bool = False, use2D: bool = True, countSimulation: bool = True, countBounds: Iterable[int] | None = None, fpSize: int = 2048, atomInvariantsGenerator: AtomInvariantsGenerator | None = None) -> FingerprintGenerator64:
    """
    Get an atom pair fingerprint generator

    ARGUMENTS:
        - minDistance: minimum distance between atoms to be considered in a
          pair, default is 1 bond
        - maxDistance: maximum distance between atoms to be considered in a
          pair, default is maxPathLen-1 bonds
        - includeChirality: if set, chirality will be used in the atom
          invariants, this is ignored if atomInvariantsGenerator is provided
        - use2D: if set, the 2D (topological) distance matrix will be used
        - countSimulation: if set, use count simulation while generating the fingerprint
        - countBounds: boundaries for count simulation, corresponding bit
          will be set if the count is higher than the number provided for that spot
        - fpSize: size of the generated fingerprint, does not affect the sparse versions
        - atomInvariantsGenerator: atom invariants to be used during fingerprint generation

    This generator supports the following AdditionalOutput types:
        - atomToBits: which bits each atom is involved in
        - atomCounts: how many bits each atom sets
        - bitInfoMap: map from bitId to (atomId, radius) pairs

    RETURNS: FingerprintGenerator
    """

def GetAtomPairAtomInvGen(includeChirality: bool = False) -> AtomInvariantsGenerator:
    """
    Get an atom pair atom-invariant generator

    ARGUMENTS:
        - includeChirality: if set, chirality will be taken into account for invariants
    RETURNS: AtomInvariantsGenerator
    """

class MorganFingerprintOptions(FingerprintOptions):
    @property
    def onlyNonzeroInvariants(self) -> bool:
        """use include atoms which have nonzero invariants"""

    @onlyNonzeroInvariants.setter
    def onlyNonzeroInvariants(self, arg: bool, /) -> None: ...

    @property
    def radius(self) -> int:
        """the radius of the fingerprints to generate"""

    @radius.setter
    def radius(self, arg: int, /) -> None: ...

    @property
    def includeRedundantEnvironments(self) -> bool:
        """include redundant environments in the fingerprint"""

    @includeRedundantEnvironments.setter
    def includeRedundantEnvironments(self, arg: bool, /) -> None: ...

def GetMorganGenerator(radius: int = 3, countSimulation: bool = False, includeChirality: bool = False, useBondTypes: bool = True, onlyNonzeroInvariants: bool = False, includeRingMembership: bool = True, countBounds: Iterable[int] | None = None, fpSize: int = 2048, atomInvariantsGenerator: AtomInvariantsGenerator | None = None, bondInvariantsGenerator: BondInvariantsGenerator | None = None, includeRedundantEnvironments: bool = False) -> FingerprintGenerator64:
    """
    Get a morgan fingerprint generator

    ARGUMENTS:
        - radius: the number of iterations to grow the fingerprint
        - countSimulation: if set, use count simulation while generating the fingerprint
        - includeChirality: if set, chirality information will be added to
          the generated fingerprint
        - useBondTypes: if set, bond types will be included as a part of
          the default bond invariants
        - countBounds: boundaries for count simulation, corresponding bit
          will be set if the count is higher than the number provided for that spot
        - fpSize: size of the generated fingerprint, does not affect the sparse versions
        - atomInvariantsGenerator: atom invariants to be used during fingerprint generation

    This generator supports the following AdditionalOutput types:
        - atomToBits: which bits each atom is the center of
        - atomCounts: how many bits each atom sets
        - bitInfoMap: map from bitId to (atomId1, radius) pairs

    RETURNS: FingerprintGenerator
    """

def GetMorganAtomInvGen(includeRingMembership: bool) -> AtomInvariantsGenerator:
    """
    Get a morgan atom invariants generator

    ARGUMENTS:
        - includeRingMembership: if set, whether or not the atom is in a ring
          will be used in the invariant list

    RETURNS: AtomInvariantsGenerator
    """

def GetMorganFeatureAtomInvGen(patterns: Iterable[rdkit.Chem.rdchem.Mol] | None = None) -> AtomInvariantsGenerator:
    """
    Get a morgan feature atom invariants generator

    ARGUMENTS:
        - patterns: if provided should contain the queries used to assign
          atom-types. if not provided, feature definitions adapted from
          reference: Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998) will
          be used for Donor, Acceptor, Aromatic, Halogen, Basic, Acidic.

    RETURNS: AtomInvariantsGenerator
    """

def GetMorganBondInvGen(useBondTypes: bool = True, useChirality: bool = False) -> BondInvariantsGenerator:
    """
    Get a morgan bond invariants generator

    ARGUMENTS:
        - useBondTypes: if set, bond types will be included as a part of the bond invariants
        - useChirality: if set, chirality information will be included as a
          part of the bond invariants

    RETURNS: BondInvariantsGenerator
    """

class RDKitFingerprintOptions(FingerprintOptions):
    @property
    def minPath(self) -> int:
        """minimum path length (in bonds) to be included"""

    @minPath.setter
    def minPath(self, arg: int, /) -> None: ...

    @property
    def maxPath(self) -> int:
        """maximum path length (in bonds) to be included"""

    @maxPath.setter
    def maxPath(self, arg: int, /) -> None: ...

    @property
    def useHs(self) -> bool:
        """use explicit Hs in the paths (if molecule has explicit Hs)"""

    @useHs.setter
    def useHs(self, arg: bool, /) -> None: ...

    @property
    def branchedPaths(self) -> bool:
        """generate branched subgraphs, not just linear ones"""

    @branchedPaths.setter
    def branchedPaths(self, arg: bool, /) -> None: ...

    @property
    def useBondOrder(self) -> bool:
        """include bond orders in the path hashes"""

    @useBondOrder.setter
    def useBondOrder(self, arg: bool, /) -> None: ...

def GetRDKitFPGenerator(minPath: int = 1, maxPath: int = 7, useHs: bool = True, branchedPaths: bool = True, useBondOrder: bool = True, countSimulation: bool = False, countBounds: Iterable[int] | None = None, fpSize: int = 2048, numBitsPerFeature: int = 2, atomInvariantsGenerator: AtomInvariantsGenerator | None = None) -> FingerprintGenerator64:
    """
    Get an RDKit fingerprint generator

    ARGUMENTS:
        - minPath: the minimum path length (in bonds) to be included
        - maxPath: the maximum path length (in bonds) to be included
        - useHs: toggles inclusion of Hs in paths (if the molecule has explicit Hs)
        - branchedPaths: toggles generation of branched subgraphs, not just linear paths
        - useBondOrder: toggles inclusion of bond orders in the path hashes
        - countSimulation: if set, use count simulation while generating the fingerprint
        - countBounds: boundaries for count simulation, corresponding bit
          will be set if the count is higher than the number provided for that spot
        - fpSize: size of the generated fingerprint, does not affect the sparse versions
        - numBitsPerFeature: the number of bits set per path/subgraph found
        - atomInvariantsGenerator: atom invariants to be used during fingerprint generation

    This generator supports the following AdditionalOutput types:
        - atomToBits: which bits each atom is involved in
        - atomCounts: how many bits each atom sets
        - bitPaths: map from bitId to vectors of bond indices for the individual subgraphs

    RETURNS: FingerprintGenerator
    """

def GetRDKitAtomInvGen() -> AtomInvariantsGenerator:
    """
    Get an RDKit atom invariants generator

    RETURNS: AtomInvariantsGenerator
    """

class TopologicalTorsionFingerprintOptions(FingerprintOptions):
    @property
    def torsionAtomCount(self) -> int:
        """number of atoms to be included in the paths"""

    @torsionAtomCount.setter
    def torsionAtomCount(self, arg: int, /) -> None: ...

    @property
    def onlyShortestPaths(self) -> bool:
        """
        whether or not to only include paths which are the shortest path between the start and end atoms
        """

    @onlyShortestPaths.setter
    def onlyShortestPaths(self, arg: bool, /) -> None: ...

def GetTopologicalTorsionGenerator(includeChirality: bool = False, torsionAtomCount: int = 4, countSimulation: bool = True, countBounds: Iterable[int] | None = None, fpSize: int = 2048, atomInvariantsGenerator: AtomInvariantsGenerator | None = None) -> FingerprintGenerator64:
    """
    Get an atom pair fingerprint generator

    ARGUMENTS:
        - includeChirality: includeChirality argument for both the default
          atom invariants generator and the fingerprint arguments
        - torsionAtomCount: the number of atoms to include in the "torsions"
        - countSimulation: if set, use count simulation while generating the fingerprint
        - countBounds: boundaries for count simulation, corresponding bit
          will be set if the count is higher than the number provided for that spot
        - fpSize: size of the generated fingerprint, does not affect the sparse versions
        - atomInvariantsGenerator: atom invariants to be used during fingerprint generation

    This generator supports the following AdditionalOutput types:
        - atomToBits: which bits each atom is involved in
        - atomCounts: how many bits each atom sets
        - bitPaths: map from bitId to vectors of atom indices

    RETURNS: FingerprintGenerator
    """
