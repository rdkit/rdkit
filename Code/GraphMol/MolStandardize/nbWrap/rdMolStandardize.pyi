"""Module containing functions for molecular standardization"""

from collections.abc import Callable, Iterable, Iterator, Sequence
import enum
from typing import Final, overload

import rdkit.Chem.rdchem


class CleanupParameters:
    """Parameters controlling molecular standardization"""

    def __init__(self) -> None: ...

    @property
    def normalizationsFile(self) -> str:
        """file containing the normalization transformations"""

    @normalizationsFile.setter
    def normalizationsFile(self, arg: str, /) -> None: ...

    @property
    def acidbaseFile(self) -> str:
        """file containing the acid and base definitions"""

    @acidbaseFile.setter
    def acidbaseFile(self, arg: str, /) -> None: ...

    @property
    def fragmentFile(self) -> str:
        """file containing the acid and base definitions"""

    @fragmentFile.setter
    def fragmentFile(self, arg: str, /) -> None: ...

    @property
    def tautomerTransformsFile(self) -> str:
        """file containing the tautomer transformations"""

    @tautomerTransformsFile.setter
    def tautomerTransformsFile(self, arg: str, /) -> None: ...

    @property
    def maxRestarts(self) -> int:
        """maximum number of restarts"""

    @maxRestarts.setter
    def maxRestarts(self, arg: int, /) -> None: ...

    @property
    def preferOrganic(self) -> bool:
        """prefer organic fragments to inorganic ones when deciding what to keep"""

    @preferOrganic.setter
    def preferOrganic(self, arg: bool, /) -> None: ...

    @property
    def doCanonical(self) -> bool:
        """
        apply atom-order dependent normalizations (like uncharging) in a canonical order
        """

    @doCanonical.setter
    def doCanonical(self, arg: bool, /) -> None: ...

    @property
    def maxTautomers(self) -> int:
        """maximum number of tautomers to generate (defaults to 1000)"""

    @maxTautomers.setter
    def maxTautomers(self, arg: int, /) -> None: ...

    @property
    def maxTransforms(self) -> int:
        """
        maximum number of transforms to apply during tautomer enumeration (defaults to 1000)
        """

    @maxTransforms.setter
    def maxTransforms(self, arg: int, /) -> None: ...

    @property
    def tautomerRemoveSp3Stereo(self) -> bool:
        """
        remove stereochemistry from sp3 centers involved in tautomerism (defaults to True)
        """

    @tautomerRemoveSp3Stereo.setter
    def tautomerRemoveSp3Stereo(self, arg: bool, /) -> None: ...

    @property
    def tautomerRemoveBondStereo(self) -> bool:
        """
        remove stereochemistry from double bonds involved in tautomerism (defaults to True)
        """

    @tautomerRemoveBondStereo.setter
    def tautomerRemoveBondStereo(self, arg: bool, /) -> None: ...

    @property
    def tautomerRemoveIsotopicHs(self) -> bool:
        """
        remove isotopic Hs from centers involved in tautomerism (defaults to True)
        """

    @tautomerRemoveIsotopicHs.setter
    def tautomerRemoveIsotopicHs(self, arg: bool, /) -> None: ...

    @property
    def tautomerReassignStereo(self) -> bool:
        """
        call AssignStereochemistry on all generated tautomers (defaults to True)
        """

    @tautomerReassignStereo.setter
    def tautomerReassignStereo(self, arg: bool, /) -> None: ...

    @property
    def largestFragmentChooserUseAtomCount(self) -> bool:
        """
        Whether LargestFragmentChooser should use atom count as main criterion before MW (defaults to True)
        """

    @largestFragmentChooserUseAtomCount.setter
    def largestFragmentChooserUseAtomCount(self, arg: bool, /) -> None: ...

    @property
    def largestFragmentChooserCountHeavyAtomsOnly(self) -> bool:
        """
        whether LargestFragmentChooser should only count heavy atoms (defaults to False)
        """

    @largestFragmentChooserCountHeavyAtomsOnly.setter
    def largestFragmentChooserCountHeavyAtomsOnly(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

def UpdateParamsFromJSON(params: CleanupParameters, json: str) -> None:
    """updates the cleanup parameters from the provided JSON string"""

def Cleanup(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None) -> rdkit.Chem.rdchem.Mol:
    """Standardizes a molecule"""

@overload
def CleanupInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None) -> None:
    """Standardizes a molecule in place"""

@overload
def CleanupInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None) -> None:
    """Standardizes multiple molecules in place"""

def StandardizeSmiles(smiles: str) -> str:
    """Convenience function for standardizing a SMILES"""

def TautomerParent(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None, skipStandardize: bool = False) -> rdkit.Chem.rdchem.Mol:
    """
    Returns the tautomer parent of a given molecule. The fragment parent is
    the standardized canonical tautomer of the molecule
    """

@overload
def TautomerParentInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the tautomer parent in place"""

@overload
def TautomerParentInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the tautomer parent in place for multiple molecules"""

def FragmentParent(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None, skipStandardize: bool = False) -> rdkit.Chem.rdchem.Mol:
    """Returns the largest fragment after doing a cleanup"""

@overload
def FragmentParentInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the largest fragment in place"""

@overload
def FragmentParentInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the largest fragment in place for multiple molecules"""

def StereoParent(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None, skipStandardize: bool = False) -> rdkit.Chem.rdchem.Mol:
    """Returns the stereo parent of the molecule"""

@overload
def StereoParentInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the stereo parent in place"""

@overload
def StereoParentInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the stereo parent in place for multiple molecules"""

def IsotopeParent(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None, skipStandardize: bool = False) -> rdkit.Chem.rdchem.Mol:
    """removes all isotopes specifications from the given molecule"""

@overload
def IsotopeParentInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the isotope parent in place"""

@overload
def IsotopeParentInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the isotope parent in place for multiple molecules"""

def ChargeParent(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None, skipStandardize: bool = False) -> rdkit.Chem.rdchem.Mol:
    """Returns the uncharged version of the largest fragment"""

@overload
def ChargeParentInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the charge parent in place"""

@overload
def ChargeParentInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the chargeparent in place for multiple molecules"""

def SuperParent(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None, skipStandardize: bool = False) -> rdkit.Chem.rdchem.Mol:
    """
    Returns the super parent. The super parent is the fragment, charge,
    isotope, stereo, and tautomer parent of the molecule.
    """

@overload
def SuperParentInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the super parent in place"""

@overload
def SuperParentInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None, skipStandardize: bool = False) -> None:
    """Generates the super parent in place for multiple molecules"""

def Normalize(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None) -> rdkit.Chem.rdchem.Mol:
    """
    Applies a series of standard transformations to correct functional
    groups and recombine charges
    """

@overload
def NormalizeInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None) -> None:
    """
    Applies a series of standard transformations to correct functional
    groups and recombine charges, modifies the input molecule
    """

@overload
def NormalizeInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None) -> None:
    """Normalizes multiple molecules in place"""

def Reionize(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None) -> rdkit.Chem.rdchem.Mol:
    """Ensures the strongest acid groups are charged first"""

@overload
def ReionizeInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None) -> None:
    """
    Ensures the strongest acid groups are charged first, modifies the input molecule
    """

@overload
def ReionizeInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None) -> None:
    """Reionizes multiple molecules in place"""

def RemoveFragments(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None) -> rdkit.Chem.rdchem.Mol:
    """Removes fragments from the molecule"""

@overload
def RemoveFragmentsInPlace(mol: rdkit.Chem.rdchem.Mol, params: CleanupParameters | None = None) -> None:
    """Removes fragments from the molecule, modifies the input molecule"""

@overload
def RemoveFragmentsInPlace(mols: Iterable[rdkit.Chem.rdchem.Mol], numThreads: int, params: CleanupParameters | None = None) -> None:
    """Removes fragments from multiple molecules in place"""

def CanonicalTautomer(mol: rdkit.Chem.rdchem.Mol | None, params: CleanupParameters | None = None) -> rdkit.Chem.rdchem.Mol:
    """Returns the canonical tautomer for the molecule"""

def DisconnectOrganometallics(mol: rdkit.Chem.rdchem.Mol, params: MetalDisconnectorOptions | None = None) -> rdkit.Chem.rdchem.Mol:
    """Returns the molecule disconnected using the organometallics rules."""

def DisconnectOrganometallicsInPlace(mol: rdkit.Chem.rdchem.Mol, params: MetalDisconnectorOptions | None = None) -> None:
    """
    Disconnects the molecule using the organometallics rules, modifies the input molecule
    """

class ValidationMethod:
    def __init__(self) -> None: ...

    def validate(self, mol: rdkit.Chem.rdchem.Mol, reportAllFailures: bool = False) -> list[str]: ...

class RDKitValidation(ValidationMethod):
    def __init__(self, allowEmptyMolecules: bool = False) -> None: ...

    @property
    def allowEmptyMolecules(self) -> bool: ...

    @allowEmptyMolecules.setter
    def allowEmptyMolecules(self, arg: bool, /) -> None: ...

class NoAtomValidation(ValidationMethod):
    def __init__(self) -> None: ...

class FragmentValidation(ValidationMethod):
    def __init__(self) -> None: ...

class NeutralValidation(ValidationMethod):
    def __init__(self) -> None: ...

class IsotopeValidation(ValidationMethod):
    def __init__(self, strict: bool = False) -> None: ...

    @property
    def strict(self) -> bool: ...

    @strict.setter
    def strict(self, arg: bool, /) -> None: ...

class MolVSValidation(ValidationMethod):
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, validations: Iterable[ValidationMethod]) -> None: ...

class AllowedAtomsValidation(ValidationMethod):
    def __init__(self, atoms: Iterable[rdkit.Chem.rdchem.Atom]) -> None: ...

class DisallowedAtomsValidation(ValidationMethod):
    def __init__(self, atoms: Iterable[rdkit.Chem.rdchem.Atom]) -> None: ...

class FeaturesValidation(ValidationMethod):
    def __init__(self, allowEnhancedStereo: bool = False, allowAromaticBondType: bool = False, allowDativeBondType: bool = False, allowQueries: bool = False, allowDummies: bool = False, allowAtomAliases: bool = False) -> None: ...

    @property
    def allowEnhancedStereo(self) -> bool: ...

    @allowEnhancedStereo.setter
    def allowEnhancedStereo(self, arg: bool, /) -> None: ...

    @property
    def allowAromaticBondType(self) -> bool: ...

    @allowAromaticBondType.setter
    def allowAromaticBondType(self, arg: bool, /) -> None: ...

    @property
    def allowDativeBondType(self) -> bool: ...

    @allowDativeBondType.setter
    def allowDativeBondType(self, arg: bool, /) -> None: ...

    @property
    def allowQueries(self) -> bool: ...

    @allowQueries.setter
    def allowQueries(self, arg: bool, /) -> None: ...

    @property
    def allowDummies(self) -> bool: ...

    @allowDummies.setter
    def allowDummies(self, arg: bool, /) -> None: ...

    @property
    def allowAtomAliases(self) -> bool: ...

    @allowAtomAliases.setter
    def allowAtomAliases(self, arg: bool, /) -> None: ...

class DisallowedRadicalValidation(ValidationMethod):
    def __init__(self) -> None: ...

class Is2DValidation(ValidationMethod):
    def __init__(self, threshold: float = 0.001) -> None: ...

    @property
    def threshold(self) -> float: ...

    @threshold.setter
    def threshold(self, arg: float, /) -> None: ...

class Layout2DValidation(ValidationMethod):
    def __init__(self, clashLimit: float = 0.15, bondLengthLimit: float = 25.0, allowLongBondsInRings: bool = True, allowAtomBondClashExemption: bool = True, minMedianBondLength: float = False) -> None: ...

    @property
    def clashLimit(self) -> float: ...

    @clashLimit.setter
    def clashLimit(self, arg: float, /) -> None: ...

    @property
    def bondLengthLimit(self) -> float: ...

    @bondLengthLimit.setter
    def bondLengthLimit(self, arg: float, /) -> None: ...

    @property
    def allowLongBondsInRings(self) -> bool: ...

    @allowLongBondsInRings.setter
    def allowLongBondsInRings(self, arg: bool, /) -> None: ...

    @property
    def allowAtomBondClashExemption(self) -> bool: ...

    @allowAtomBondClashExemption.setter
    def allowAtomBondClashExemption(self, arg: bool, /) -> None: ...

    @property
    def minMedianBondLength(self) -> float: ...

    @minMedianBondLength.setter
    def minMedianBondLength(self, arg: float, /) -> None: ...

class StereoValidation(ValidationMethod):
    def __init__(self) -> None: ...

def ValidateSmiles(mol: str) -> list[str]: ...

class ChargeCorrection:
    def __init__(self, name: str, smarts: str, charge: int) -> None: ...

    @property
    def Name(self) -> str: ...

    @Name.setter
    def Name(self, arg: str, /) -> None: ...

    @property
    def Smarts(self) -> str: ...

    @Smarts.setter
    def Smarts(self, arg: str, /) -> None: ...

    @property
    def Charge(self) -> int: ...

    @Charge.setter
    def Charge(self, arg: int, /) -> None: ...

def CHARGE_CORRECTIONS() -> list[ChargeCorrection]: ...

class Reionizer:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, acidbaseFile: str) -> None: ...

    @overload
    def __init__(self, acidbaseFile: str, ccs: Sequence[ChargeCorrection]) -> None: ...

    def reionize(self, mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol: ...

    def reionizeInPlace(self, mol: rdkit.Chem.rdchem.Mol) -> None:
        """modifies the input molecule"""

def ReionizerFromData(paramData: str | bytes, chargeCorrections: Sequence[ChargeCorrection] | None = None) -> Reionizer:
    """
    creates a reionizer from a string containing parameter data and a list of charge corrections
    """

class Uncharger:
    def __init__(self, canonicalOrder: bool = True, force: bool = False, protonationOnly: bool = False) -> None: ...

    def uncharge(self, mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol: ...

    def unchargeInPlace(self, mol: rdkit.Chem.rdchem.Mol) -> None:
        """modifies the input molecule"""

class MetalDisconnectorOptions:
    """Metal Disconnector Options"""

    def __init__(self) -> None: ...

    @property
    def splitGrignards(self) -> bool:
        """Whether to split Grignard-type complexes. Default false."""

    @splitGrignards.setter
    def splitGrignards(self, arg: bool, /) -> None: ...

    @property
    def splitAromaticC(self) -> bool:
        """Whether to split metal-aromatic C bonds.  Default false."""

    @splitAromaticC.setter
    def splitAromaticC(self, arg: bool, /) -> None: ...

    @property
    def adjustCharges(self) -> bool:
        """Whether to adjust charges on ligand atoms.  Default true."""

    @adjustCharges.setter
    def adjustCharges(self, arg: bool, /) -> None: ...

    @property
    def removeHapticDummies(self) -> bool:
        """
        Whether to remove the dummy atoms representing haptic bonds.  Such dummies are bonded to the metal with a bond that has the MolFileBondEndPts prop set.  Default false.
        """

    @removeHapticDummies.setter
    def removeHapticDummies(self, arg: bool, /) -> None: ...

class MetalDisconnector:
    """
    a class to disconnect metals that are defined as covalently bonded to non-metals
    """

    def __init__(self, options: MetalDisconnectorOptions | None = None) -> None: ...

    @property
    def MetalNof(self) -> str:
        """
        SMARTS defining the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine
        """

    @property
    def MetalNon(self) -> str:
        """SMARTS defining the metals to disconnect other inorganic elements"""

    def SetMetalNon(self, mol: rdkit.Chem.rdchem.Mol) -> None:
        """
        Set the query molecule defining the metals to disconnect from other inorganic elements.
        """

    def SetMetalNof(self, mol: rdkit.Chem.rdchem.Mol) -> None:
        """
        Set the query molecule defining the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine.
        """

    def Disconnect(self, mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
        """performs the disconnection"""

    def DisconnectInPlace(self, mol: rdkit.Chem.rdchem.Mol) -> None:
        """performs the disconnection, modifies the input molecule"""

class FragmentRemover:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, fragmentFilename: str = '', leave_last: bool = True, skip_if_all_match: bool = False) -> None: ...

    def remove(self, mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol: ...

    def removeInPlace(self, mol: rdkit.Chem.rdchem.Mol) -> None:
        """modifies the molecule in place"""

def FragmentRemoverFromData(fragmentData: str | bytes, leave_last: bool = True, skip_if_all_match: bool = False) -> FragmentRemover:
    """creates a FragmentRemover from a string containing parameter data"""

class LargestFragmentChooser:
    @overload
    def __init__(self, preferOrganic: bool = False) -> None: ...

    @overload
    def __init__(self, params: CleanupParameters) -> None: ...

    def choose(self, mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol: ...

    def chooseInPlace(self, mol: rdkit.Chem.rdchem.Mol) -> None:
        """modifies the molecule in place"""

class Normalizer:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, normalizeFilename: str, maxRestarts: int) -> None: ...

    def normalize(self, mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol: ...

    def normalizeInPlace(self, mol: rdkit.Chem.rdchem.Mol) -> None:
        """modifies the input molecule"""

def NormalizerFromData(paramData: str | bytes, params: CleanupParameters) -> Normalizer:
    """creates a Normalizer from a string containing normalization SMARTS"""

def NormalizerFromParams(params: CleanupParameters) -> Normalizer:
    """creates a Normalizer from CleanupParameters"""

class TautomerEnumeratorStatus(enum.Enum):
    Completed = 0

    MaxTautomersReached = 1

    MaxTransformsReached = 2

    Canceled = 3

class TautomerEnumeratorCallback:
    """
    Create a derived class from this abstract base class and
    implement the __call__() method.
    The __call__() method is called in the innermost loop of the
    algorithm, and provides a mechanism to monitor or stop
    its progress.

    To have your callback called, pass an instance of your
    derived class to TautomerEnumerator.SetCallback()
    """

    def __init__(self) -> None: ...

    def __call__(self, mol: rdkit.Chem.rdchem.Mol, res: TautomerEnumeratorResult) -> bool:
        """
        This must be implemented in the derived class. Return True if the tautomer enumeration should continue; False if the tautomer enumeration should stop.
        """

class Tautomer:
    """used to hold the aromatic and kekulized versions of each tautomer"""

    @property
    def tautomer(self) -> rdkit.Chem.rdchem.Mol:
        """aromatic version of the tautomer"""

    @property
    def kekulized(self) -> rdkit.Chem.rdchem.Mol:
        """kekulized version of the tautomer"""

class SmilesTautomerMap:
    """maps SMILES strings to the respective Tautomer objects"""

    def keys(self) -> tuple[str, ...]: ...

    def values(self) -> tuple[Tautomer, ...]: ...

    def items(self) -> tuple[tuple[str, Tautomer], ...]: ...

    def __len__(self) -> int: ...

class TautomerEnumeratorResult:
    """used to return tautomer enumeration results"""

    @property
    def tautomers(self) -> list[rdkit.Chem.rdchem.Mol]:
        """tautomers generated by the enumerator"""

    @property
    def smiles(self) -> list[str]:
        """SMILES of tautomers generated by the enumerator"""

    @property
    def smilesTautomerMap(self) -> SmilesTautomerMap:
        """dictionary mapping SMILES strings to the respective Tautomer objects"""

    @property
    def status(self) -> TautomerEnumeratorStatus:
        """
        whether the enumeration completed or not; see TautomerEnumeratorStatus for possible values
        """

    @property
    def modifiedAtoms(self) -> tuple[int, ...]:
        """tuple of atom indices modified by the transforms"""

    @property
    def modifiedBonds(self) -> tuple[int, ...]:
        """tuple of bond indices modified by the transforms"""

    def __call__(self) -> list[rdkit.Chem.rdchem.Mol]:
        """tautomers generated by the enumerator"""

    def __iter__(self) -> Iterator[rdkit.Chem.rdchem.Mol]: ...

    def __getitem__(self, pos: int) -> rdkit.Chem.rdchem.Mol: ...

    def __len__(self) -> int: ...

class TautomerEnumerator:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, params: CleanupParameters) -> None: ...

    @overload
    def __init__(self, other: TautomerEnumerator) -> None: ...

    def Enumerate(self, mol: rdkit.Chem.rdchem.Mol) -> TautomerEnumeratorResult:
        """
        Generates the tautomers for a molecule.

        The enumeration rules are inspired by the publication:
        M. Sitzmann et al., "Tautomerism in Large Databases.", JCAMD 24:521 (2010)
        https://doi.org/10.1007/s10822-010-9346-4

        Note: the definitions used here are that the atoms modified during
        tautomerization are the atoms at the beginning and end of each tautomer
        transform (the H "donor" and H "acceptor" in the transform) and the bonds
        modified during transformation are any bonds whose order is changed during
        the tautomer transform (these are the bonds between the "donor" and the
        "acceptor").
        """

    @overload
    def Canonicalize(self, mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
        """
        Returns the canonical tautomer for a molecule.

        The default scoring scheme is inspired by the publication:
        M. Sitzmann et al., "Tautomerism in Large Databases.", JCAMD 24:521 (2010)
        https://doi.org/10.1007/s10822-010-9346-4

        Note that the canonical tautomer is very likely not the most stable tautomer
        for any given conditions. The default scoring rules are designed to produce
        "reasonable" tautomers, but the primary concern is that the results are
        canonical: you always get the same canonical tautomer for a molecule
        regardless of what the input tautomer or atom ordering were.
        """

    @overload
    def Canonicalize(self, mol: rdkit.Chem.rdchem.Mol, scoreFunc: Callable[[rdkit.Chem.rdchem.Mol], int]) -> rdkit.Chem.rdchem.Mol:
        """
        picks the canonical tautomer from an iterable of molecules using a custom scoring function
        """

    @overload
    def PickCanonical(self, iterable: object) -> rdkit.Chem.rdchem.Mol:
        """picks the canonical tautomer from an iterable of molecules"""

    @overload
    def PickCanonical(self, iterable: object, scoreFunc: Callable[[rdkit.Chem.rdchem.Mol], int]) -> rdkit.Chem.rdchem.Mol:
        """
        returns the canonical tautomer for a molecule using a custom scoring function
        """

    @staticmethod
    def ScoreTautomer(mol: rdkit.Chem.rdchem.Mol) -> int:
        """returns the score for a tautomer using the default scoring scheme."""

    def SetMaxTautomers(self, maxTautomers: int) -> None:
        """set the maximum number of tautomers to be generated."""

    def GetMaxTautomers(self) -> int:
        """returns the maximum number of tautomers to be generated."""

    def SetMaxTransforms(self, maxTransforms: int) -> None:
        """
        set the maximum number of transformations to be applied. This limit is usually hit earlier than the maxTautomers limit and leads to a more linear scaling of CPU time with increasing number of tautomeric centers (see Sitzmann et al.).
        """

    def GetMaxTransforms(self) -> int:
        """returns the maximum number of transformations to be applied."""

    def SetRemoveSp3Stereo(self, removeSp3Stereo: bool) -> None:
        """
        set to True if you wish stereochemistry information to be removed from sp3 atoms involved in tautomerism. This means that S-aminoacids will lose their stereochemistry after going through tautomer enumeration because of the amido-imidol tautomerism. This defaults to True in RDKit, and to False in the workflow described by Sitzmann et al.
        """

    def GetRemoveSp3Stereo(self) -> bool:
        """
        returns whether stereochemistry information will be removed from sp3 atoms involved in tautomerism.
        """

    def SetRemoveBondStereo(self, removeBondStereo: bool) -> None:
        """
        set to True if you wish stereochemistry information to be removed from double bonds involved in tautomerism. This means that enols will lose their E/Z stereochemistry after going through tautomer enumeration because of the keto-enolic tautomerism. This defaults to True in the RDKit and also in the workflow described by Sitzmann et al.
        """

    def GetRemoveBondStereo(self) -> bool:
        """
        returns whether stereochemistry information will be removed from double bonds involved in tautomerism.
        """

    def SetReassignStereo(self, reassignStereo: bool) -> None:
        """
        set to True if you wish AssignStereochemistry to be called on each tautomer generated by the Enumerate() method. This defaults to True.
        """

    def GetReassignStereo(self) -> bool:
        """
        returns whether AssignStereochemistry will be called on each tautomer generated by the Enumerate() method.
        """

    def SetCallback(self, callback: object) -> None:
        """
        Pass an instance of a class derived from
        TautomerEnumeratorCallback, which must implement the
        __call__() method.
        """

    def GetCallback(self) -> object:
        """
        Get the TautomerEnumeratorCallback subclass instance,
        or None if none was set.
        """

    tautomerScoreVersion: Final[str] = ...
    """(arg: object, /) -> str"""

def GetV1TautomerEnumerator() -> TautomerEnumerator:
    """return a TautomerEnumerator using v1 of the enumeration rules"""

def ScoreRings(mol: rdkit.Chem.rdchem.Mol) -> int:
    """
    scores the ring system of the tautomer for canonicalization
    Aromatic rings score 100, all carbon aromatic rings score 250
    """

def ScoreHeteroHs(mol: rdkit.Chem.rdchem.Mol) -> int:
    """
    scores the number of heteroHs of the tautomer for canonicalization
    This gives a negative penalty to hydrogens attached to S,P, Se and Te
    """

class SubstructTerm:
    """
    Sets the score of this particular tautomer substructure, higher scores are more preferable
    Aromatic rings score 100, all carbon aromatic rings score 250
    """

    def __init__(self, name: str, smarts: str, score: int) -> None: ...

    @property
    def name(self) -> str: ...

    @property
    def smarts(self) -> str: ...

    @property
    def score(self) -> int: ...

class SubstructTermVector:
    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, arg: SubstructTermVector) -> None:
        """Copy constructor"""

    @overload
    def __init__(self, arg: Iterable[SubstructTerm], /) -> None:
        """Construct from an iterable object"""

    def __len__(self) -> int: ...

    def __bool__(self) -> bool:
        """Check whether the vector is nonempty"""

    def __repr__(self) -> str: ...

    def __iter__(self) -> Iterator[SubstructTerm]: ...

    @overload
    def __getitem__(self, arg: int, /) -> SubstructTerm: ...

    @overload
    def __getitem__(self, arg: slice, /) -> SubstructTermVector: ...

    def clear(self) -> None:
        """Remove all items from list."""

    def append(self, arg: SubstructTerm, /) -> None:
        """Append ``arg`` to the end of the list."""

    def insert(self, arg0: int, arg1: SubstructTerm, /) -> None:
        """Insert object ``arg1`` before index ``arg0``."""

    def pop(self, index: int = -1) -> SubstructTerm:
        """Remove and return item at ``index`` (default last)."""

    def extend(self, arg: SubstructTermVector, /) -> None:
        """Extend ``self`` by appending elements from ``arg``."""

    @overload
    def __setitem__(self, arg0: int, arg1: SubstructTerm, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: SubstructTermVector, /) -> None: ...

    @overload
    def __delitem__(self, arg: int, /) -> None: ...

    @overload
    def __delitem__(self, arg: slice, /) -> None: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    @overload
    def __contains__(self, arg: SubstructTerm, /) -> bool: ...

    @overload
    def __contains__(self, arg: object, /) -> bool: ...

    def count(self, arg: SubstructTerm, /) -> int:
        """Return number of occurrences of ``arg``."""

    def remove(self, arg: SubstructTerm, /) -> None:
        """Remove first occurrence of ``arg``."""

def ScoreSubstructs(mol: rdkit.Chem.rdchem.Mol, terms: SubstructTermVector | None = None) -> int:
    """scores the tautomer substructures"""

def GetDefaultTautomerScoreSubstructs() -> SubstructTermVector:
    """Return the default tautomer substructure scoring terms"""

class PipelineOptions:
    def __init__(self) -> None: ...

    @property
    def strictParsing(self) -> bool: ...

    @strictParsing.setter
    def strictParsing(self, arg: bool, /) -> None: ...

    @property
    def reportAllFailures(self) -> bool: ...

    @reportAllFailures.setter
    def reportAllFailures(self, arg: bool, /) -> None: ...

    @property
    def allowEmptyMolecules(self) -> bool: ...

    @allowEmptyMolecules.setter
    def allowEmptyMolecules(self, arg: bool, /) -> None: ...

    @property
    def allowEnhancedStereo(self) -> bool: ...

    @allowEnhancedStereo.setter
    def allowEnhancedStereo(self, arg: bool, /) -> None: ...

    @property
    def allowAromaticBondType(self) -> bool: ...

    @allowAromaticBondType.setter
    def allowAromaticBondType(self, arg: bool, /) -> None: ...

    @property
    def allowDativeBondType(self) -> bool: ...

    @allowDativeBondType.setter
    def allowDativeBondType(self, arg: bool, /) -> None: ...

    @property
    def is2DZeroThreshold(self) -> float: ...

    @is2DZeroThreshold.setter
    def is2DZeroThreshold(self, arg: float, /) -> None: ...

    @property
    def atomClashLimit(self) -> float: ...

    @atomClashLimit.setter
    def atomClashLimit(self, arg: float, /) -> None: ...

    @property
    def minMedianBondLength(self) -> float: ...

    @minMedianBondLength.setter
    def minMedianBondLength(self, arg: float, /) -> None: ...

    @property
    def bondLengthLimit(self) -> float: ...

    @bondLengthLimit.setter
    def bondLengthLimit(self, arg: float, /) -> None: ...

    @property
    def allowLongBondsInRings(self) -> bool: ...

    @allowLongBondsInRings.setter
    def allowLongBondsInRings(self, arg: bool, /) -> None: ...

    @property
    def allowAtomBondClashExemption(self) -> bool: ...

    @allowAtomBondClashExemption.setter
    def allowAtomBondClashExemption(self, arg: bool, /) -> None: ...

    @property
    def metalNof(self) -> str: ...

    @metalNof.setter
    def metalNof(self, arg: str, /) -> None: ...

    @property
    def metalNon(self) -> str: ...

    @metalNon.setter
    def metalNon(self, arg: str, /) -> None: ...

    @property
    def normalizerData(self) -> str: ...

    @normalizerData.setter
    def normalizerData(self, arg: str, /) -> None: ...

    @property
    def normalizerMaxRestarts(self) -> int: ...

    @normalizerMaxRestarts.setter
    def normalizerMaxRestarts(self, arg: int, /) -> None: ...

    @property
    def scaledMedianBondLength(self) -> float: ...

    @scaledMedianBondLength.setter
    def scaledMedianBondLength(self, arg: float, /) -> None: ...

    @property
    def outputV2000(self) -> bool: ...

    @outputV2000.setter
    def outputV2000(self, arg: bool, /) -> None: ...

class PipelineStatus(enum.IntFlag):
    __str__ = __repr__

    def __repr__(self, /):
        """Return repr(self)."""

    NO_EVENT = 0

    INPUT_ERROR = 1

    PREPARE_FOR_VALIDATION_ERROR = 2

    FEATURES_VALIDATION_ERROR = 4

    BASIC_VALIDATION_ERROR = 8

    IS2D_VALIDATION_ERROR = 16

    LAYOUT2D_VALIDATION_ERROR = 32

    STEREO_VALIDATION_ERROR = 64

    VALIDATION_ERROR = 124

    PREPARE_FOR_STANDARDIZATION_ERROR = 128

    METAL_STANDARDIZATION_ERROR = 256

    NORMALIZER_STANDARDIZATION_ERROR = 512

    FRAGMENT_STANDARDIZATION_ERROR = 1024

    CHARGE_STANDARDIZATION_ERROR = 2048

    STANDARDIZATION_ERROR = 3840

    OUTPUT_ERROR = 4096

    PIPELINE_ERROR = 8191

    METALS_DISCONNECTED = 8388608

    NORMALIZATION_APPLIED = 16777216

    FRAGMENTS_REMOVED = 33554432

    PROTONATION_CHANGED = 67108864

    STRUCTURE_MODIFICATION = 125829120

class PipelineStage(enum.Enum):
    PARSING_INPUT = 1

    PREPARE_FOR_VALIDATION = 2

    VALIDATION = 3

    PREPARE_FOR_STANDARDIZATION = 4

    STANDARDIZATION = 5

    SERIALIZING_OUTPUT = 9

    COMPLETED = 10

class PipelineLogEntry:
    @property
    def status(self) -> PipelineStatus: ...

    @property
    def detail(self) -> str: ...

class PipelineLog:
    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, arg: PipelineLog) -> None:
        """Copy constructor"""

    @overload
    def __init__(self, arg: Iterable[PipelineLogEntry], /) -> None:
        """Construct from an iterable object"""

    def __len__(self) -> int: ...

    def __bool__(self) -> bool:
        """Check whether the vector is nonempty"""

    def __repr__(self) -> str: ...

    def __iter__(self) -> Iterator[PipelineLogEntry]: ...

    @overload
    def __getitem__(self, arg: int, /) -> PipelineLogEntry: ...

    @overload
    def __getitem__(self, arg: slice, /) -> PipelineLog: ...

    def clear(self) -> None:
        """Remove all items from list."""

    def append(self, arg: PipelineLogEntry, /) -> None:
        """Append ``arg`` to the end of the list."""

    def insert(self, arg0: int, arg1: PipelineLogEntry, /) -> None:
        """Insert object ``arg1`` before index ``arg0``."""

    def pop(self, index: int = -1) -> PipelineLogEntry:
        """Remove and return item at ``index`` (default last)."""

    def extend(self, arg: PipelineLog, /) -> None:
        """Extend ``self`` by appending elements from ``arg``."""

    @overload
    def __setitem__(self, arg0: int, arg1: PipelineLogEntry, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: PipelineLog, /) -> None: ...

    @overload
    def __delitem__(self, arg: int, /) -> None: ...

    @overload
    def __delitem__(self, arg: slice, /) -> None: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    @overload
    def __contains__(self, arg: PipelineLogEntry, /) -> bool: ...

    @overload
    def __contains__(self, arg: object, /) -> bool: ...

    def count(self, arg: PipelineLogEntry, /) -> int:
        """Return number of occurrences of ``arg``."""

    def remove(self, arg: PipelineLogEntry, /) -> None:
        """Remove first occurrence of ``arg``."""

class PipelineResult:
    @property
    def status(self) -> PipelineStatus: ...

    @property
    def stage(self) -> PipelineStage: ...

    @property
    def log(self) -> PipelineLog: ...

    @property
    def inputMolData(self) -> str: ...

    @property
    def outputMolData(self) -> str: ...

    @property
    def parentMolData(self) -> str: ...

class Pipeline:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, options: PipelineOptions) -> None: ...

    def run(self, molData: str) -> PipelineResult: ...
