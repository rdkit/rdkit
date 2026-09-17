"""Module containing functions for working with groups of molecules"""

import os
from typing import overload

import rdkit.Chem.rdFingerprintGenerator


class SupplierOptions:
    """Supplier Options"""

    def __init__(self) -> None: ...

    @property
    def numThreads(self) -> int:
        """the number of threads to use while working"""

    @numThreads.setter
    def numThreads(self, arg: int, /) -> None: ...

    @property
    def sanitize(self) -> bool: ...

    @sanitize.setter
    def sanitize(self, arg: bool, /) -> None: ...

    @property
    def removeHs(self) -> bool: ...

    @removeHs.setter
    def removeHs(self, arg: bool, /) -> None: ...

    @property
    def strictParsing(self) -> bool: ...

    @strictParsing.setter
    def strictParsing(self, arg: bool, /) -> None: ...

    @property
    def delimiter(self) -> str:
        """used for SMILES files"""

    @delimiter.setter
    def delimiter(self, arg: str, /) -> None: ...

    @property
    def smilesColumn(self) -> int:
        """used for SMILES files"""

    @smilesColumn.setter
    def smilesColumn(self, arg: int, /) -> None: ...

    @property
    def nameColumn(self) -> int:
        """used for SMILES files"""

    @nameColumn.setter
    def nameColumn(self, arg: int, /) -> None: ...

    @property
    def titleLine(self) -> bool:
        """used for SMILES files"""

    @titleLine.setter
    def titleLine(self, arg: bool, /) -> None: ...

    @property
    def nameRecord(self) -> str:
        """used for TDT files"""

    @nameRecord.setter
    def nameRecord(self, arg: str, /) -> None: ...

    @property
    def confId2D(self) -> int:
        """used for TDT files"""

    @confId2D.setter
    def confId2D(self, arg: int, /) -> None: ...

    @property
    def confId3D(self) -> int:
        """used for TDT files"""

    @confId3D.setter
    def confId3D(self, arg: int, /) -> None: ...

@overload
def GetFingerprintsForMolsInFile(filename: str | os.PathLike, generator: rdkit.Chem.rdFingerprintGenerator.FingerprintGenerator32 | None = None, options: SupplierOptions = ...) -> tuple:
    """returns the fingerprints for the molecules in a file (32 bit version)"""

@overload
def GetFingerprintsForMolsInFile(filename: str | os.PathLike, generator: rdkit.Chem.rdFingerprintGenerator.FingerprintGenerator64 | None = None, options: SupplierOptions = ...) -> tuple:
    """returns the fingerprints for the molecules in a file (64 bit version)"""
