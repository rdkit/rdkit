"""Module containing classes and functions for enumerating molecules"""

import enum
from typing import overload

import rdkit.Chem.rdchem


class EnumeratorType(enum.Enum):
    LinkNode = 0

    PositionVariation = 1

    RepeatUnit = 2

class MolEnumeratorParams:
    """Molecular enumerator parameters"""

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, typ: EnumeratorType) -> None: ...

    @property
    def sanitize(self) -> bool:
        """sanitize molecules after enumeration"""

    @sanitize.setter
    def sanitize(self, arg: bool, /) -> None: ...

    @property
    def maxToEnumerate(self) -> int:
        """maximum number of molecules to enumerate"""

    @maxToEnumerate.setter
    def maxToEnumerate(self, arg: int, /) -> None: ...

    @property
    def doRandom(self) -> bool:
        """do random enumeration (not yet implemented)"""

    @doRandom.setter
    def doRandom(self, arg: bool, /) -> None: ...

    @property
    def randomSeed(self) -> int:
        """seed for the random enumeration (not yet implemented)"""

    @randomSeed.setter
    def randomSeed(self, arg: int, /) -> None: ...

    def SetEnumerationOperator(self, typ: EnumeratorType) -> None:
        """set the operator to be used for enumeration"""

@overload
def Enumerate(mol: rdkit.Chem.rdchem.Mol, maxPerOperation: int = 0) -> rdkit.Chem.rdchem.MolBundle:
    """
    do an enumeration and return a MolBundle.
    If maxPerOperation is >0 that will be used as the maximum number of molecules which
    can be returned by any given operation.
    Limitations:
      - the current implementation does not support molecules which include both
        SRUs and LINKNODEs
      - Overlapping SRUs, i.e. where one monomer is contained within another, are
        not supported
    """

@overload
def Enumerate(mol: rdkit.Chem.rdchem.Mol, enumParams: MolEnumeratorParams) -> rdkit.Chem.rdchem.MolBundle:
    """
    do an enumeration for the supplied parameter type and return a MolBundle
    Limitations:
      - Overlapping SRUs, i.e. where one monomer is contained within another, are
        not supported
    """
