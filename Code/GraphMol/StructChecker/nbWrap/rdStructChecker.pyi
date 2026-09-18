import enum
from typing import overload

import rdkit.Chem.rdchem


class StructureFlags(enum.IntEnum):
    NO_CHANGE = 0

    BAD_MOLECULE = 1

    ALIAS_CONVERSION_FAILED = 2

    STEREO_ERROR = 4

    STEREO_FORCED_BAD = 8

    ATOM_CLASH = 16

    ATOM_CHECK_FAILED = 32

    SIZE_CHECK_FAILED = 64

    TRANSFORMED = 256

    FRAGMENTS_FOUND = 512

    EITHER_WARNING = 1024

    DUBIOUS_STEREO_REMOVED = 2048

    RECHARGED = 4096

    STEREO_TRANSFORMED = 8192

    TEMPLATE_TRANSFORMED = 16384

    TAUTOMER_TRANSFORMED = 32768

class StructCheckerOptions:
    def __init__(self) -> None: ...

    @property
    def AcidityLimit(self) -> float: ...

    @AcidityLimit.setter
    def AcidityLimit(self, arg: float, /) -> None: ...

    @property
    def RemoveMinorFragments(self) -> bool: ...

    @RemoveMinorFragments.setter
    def RemoveMinorFragments(self, arg: bool, /) -> None: ...

    @property
    def DesiredCharge(self) -> int: ...

    @DesiredCharge.setter
    def DesiredCharge(self, arg: int, /) -> None: ...

    @property
    def CheckCollisions(self) -> bool: ...

    @CheckCollisions.setter
    def CheckCollisions(self, arg: bool, /) -> None: ...

    @property
    def CollisionLimitPercent(self) -> int: ...

    @CollisionLimitPercent.setter
    def CollisionLimitPercent(self, arg: int, /) -> None: ...

    @property
    def MaxMolSize(self) -> int: ...

    @MaxMolSize.setter
    def MaxMolSize(self, arg: int, /) -> None: ...

    @property
    def ConvertSText(self) -> bool: ...

    @ConvertSText.setter
    def ConvertSText(self, arg: bool, /) -> None: ...

    @property
    def StripZeros(self) -> bool: ...

    @StripZeros.setter
    def StripZeros(self, arg: bool, /) -> None: ...

    @property
    def CheckStereo(self) -> bool: ...

    @CheckStereo.setter
    def CheckStereo(self, arg: bool, /) -> None: ...

    @property
    def ConvertAtomTexts(self) -> bool: ...

    @ConvertAtomTexts.setter
    def ConvertAtomTexts(self, arg: bool, /) -> None: ...

    @property
    def GroupsToSGroups(self) -> bool: ...

    @GroupsToSGroups.setter
    def GroupsToSGroups(self, arg: bool, /) -> None: ...

    @property
    def Verbose(self) -> bool: ...

    @Verbose.setter
    def Verbose(self, arg: bool, /) -> None: ...

    def LoadGoodAugmentedAtoms(self, path: str) -> bool:
        """Load the set of good augmented atoms from the specified file path"""

    def LoadAcidicAugmentedAtoms(self, path: str) -> bool:
        """
        Load the set of acidic augmented atoms from the specified file
        path
        """

    def LoadAugmentedAtomTranslations(self, path: str) -> bool:
        """
        Load the set of acidic augmented atoms from the specified file
        path
        """

class StructChecker:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: StructCheckerOptions, /) -> None: ...

    def CheckMolStructure(self, mol: rdkit.Chem.rdchem.Mol) -> int:
        """Check the structure and return a set of structure flags"""

    @staticmethod
    def StructureFlagsToString(flags: int) -> str:
        """Return the structure flags as a human readable string"""

    @staticmethod
    def StringToStructureFlags(str: str) -> int:
        """
        Convert a comma separated string to the appropriate structure
        flags
        """
