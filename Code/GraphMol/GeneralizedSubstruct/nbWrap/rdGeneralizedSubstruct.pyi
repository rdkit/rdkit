"""Module containing functions for generalized substructure searching"""

from typing import overload

import rdkit.Chem.rdchem
import rdkit.Chem.rdmolops
import rdkit.DataStructs.cDataStructs


class ExtendedQueryMol:
    """Extended query molecule for use in generalized substructure searching."""

    @overload
    def __init__(self, text: str, isJSON: bool = False) -> None:
        """
        constructor from either a binary string (from ToBinary()) or a JSON string.
        """

    @overload
    def __init__(self, data: bytes) -> None:
        """constructor from binary data returned by ToBinary()."""

    def InitFromBinary(self, pkl: str) -> None: ...

    def InitFromJSON(self, text: str) -> None: ...

    def ToBinary(self) -> bytes: ...

    def ToJSON(self) -> str: ...

    def PatternFingerprintQuery(self, fingerprintSize: int = 2048) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect: ...

def MolHasSubstructMatch(mol: rdkit.Chem.rdchem.Mol, query: ExtendedQueryMol, params: rdkit.Chem.rdchem.SubstructMatchParameters = ...) -> bool:
    """
    determines whether or not a molecule is a match to a generalized substructure query
    """

def MolGetSubstructMatch(mol: rdkit.Chem.rdchem.Mol, query: ExtendedQueryMol, params: rdkit.Chem.rdchem.SubstructMatchParameters = ...) -> list[int]:
    """
    returns first match (if any) of a molecule to a generalized substructure query
    """

def MolGetSubstructMatches(mol: rdkit.Chem.rdchem.Mol, query: ExtendedQueryMol, params: rdkit.Chem.rdchem.SubstructMatchParameters = ...) -> list[list[int]]:
    """
    returns all matches (if any) of a molecule to a generalized substructure query
    """

def PatternFingerprintTarget(target: rdkit.Chem.rdchem.Mol, fingerprintSize: int = 2048) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """
    Creates a pattern fingerprint for a target molecule that is compatible with an extended query
    """

def CreateExtendedQueryMol(mol: rdkit.Chem.rdchem.Mol, doEnumeration: bool = True, doTautomers: bool = True, adjustQueryProperties: bool = False, adjustQueryParameters: rdkit.Chem.rdmolops.AdjustQueryParameters = ...) -> ExtendedQueryMol:
    """
    Creates an ExtendedQueryMol from the input molecule

    This takes a query molecule and, conceptually, performs the following steps to
    produce an ExtendedQueryMol:

      1. Enumerates features like Link Nodes and SRUs
      2. Converts everything into TautomerQueries
      3. Runs adjustQueryProperties()

    Each step is optional
    """
