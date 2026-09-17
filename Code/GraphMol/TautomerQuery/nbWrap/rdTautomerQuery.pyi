"""
Module for tautomer-aware substructure searching.

Provides the TautomerQuery class which enables substructure searching
that accounts for tautomeric forms of the query molecule.
"""

from typing import overload

import rdkit.Chem.rdchem
import rdkit.DataStructs.cDataStructs


class TautomerQuery:
    """
    The Tautomer Query Class.
    Creates a query that enables structure search accounting for matching of
    Tautomeric forms
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, pickle: bytes) -> None:
        """Construct a TautomerQuery from a pickle bytes."""

    @overload
    def __init__(self, pickle: str) -> None:
        """Construct a TautomerQuery from a pickle string."""

    @overload
    def __init__(self, mol: rdkit.Chem.rdchem.Mol, tautomerTransformFile: str = '') -> None:
        """Construct a TautomerQuery from a molecule."""

    @overload
    def IsSubstructOf(self, target: rdkit.Chem.rdchem.Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
        Check if this tautomer query is a substructure of the target molecule.

        ARGUMENTS:
         - target: the target molecule
         - recursionPossible: (optional) allow recursive queries (default True)
         - useChirality: (optional) use chirality in matching (default False)
         - useQueryQueryMatches: (optional) use query-query matching logic (default False)

        RETURNS: True or False
        """

    @overload
    def IsSubstructOf(self, target: rdkit.Chem.rdchem.Mol, params: rdkit.Chem.rdchem.SubstructMatchParameters) -> bool:
        """
        Check if this tautomer query is a substructure of the target molecule.

        ARGUMENTS:
         - target: the target molecule
         - params: SubstructMatchParameters object

        RETURNS: True or False
        """

    @overload
    def GetSubstructMatch(self, target: rdkit.Chem.rdchem.Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> list[int]:
        """
        Return the first substructure match of this tautomer query in the target.

        ARGUMENTS:
         - target: the target molecule
         - useChirality: (optional) use chirality in matching (default False)
         - useQueryQueryMatches: (optional) use query-query matching logic (default False)

        RETURNS: a tuple of atom indices on match, or empty tuple on no match
        """

    @overload
    def GetSubstructMatch(self, target: rdkit.Chem.rdchem.Mol, params: rdkit.Chem.rdchem.SubstructMatchParameters) -> list[int]:
        """
        Return the first substructure match of this tautomer query in the target.

        ARGUMENTS:
         - target: the target molecule
         - params: SubstructMatchParameters object

        RETURNS: a tuple of atom indices on match, or empty tuple on no match
        """

    @overload
    def GetSubstructMatches(self, target: rdkit.Chem.rdchem.Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> list[list[int]]:
        """
        Return all substructure matches of this tautomer query in the target.

        ARGUMENTS:
         - target: the target molecule
         - uniquify: (optional) only return unique matches (default True)
         - useChirality: (optional) use chirality in matching (default False)
         - useQueryQueryMatches: (optional) use query-query matching logic (default False)
         - maxMatches: (optional) maximum number of matches to return (default 1000)

        RETURNS: a tuple of tuples of atom indices
        """

    @overload
    def GetSubstructMatches(self, target: rdkit.Chem.rdchem.Mol, params: rdkit.Chem.rdchem.SubstructMatchParameters) -> list[list[int]]:
        """
        Return all substructure matches of this tautomer query in the target.

        ARGUMENTS:
         - target: the target molecule
         - params: SubstructMatchParameters object

        RETURNS: a tuple of tuples of atom indices
        """

    @overload
    def GetSubstructMatchesWithTautomers(self, target: rdkit.Chem.rdchem.Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> list[tuple[list[int], rdkit.Chem.rdchem.Mol]]:
        """
        Return all substructure matches with their matching tautomers.

        ARGUMENTS:
         - target: the target molecule
         - uniquify: (optional) only return unique matches (default True)
         - useChirality: (optional) use chirality in matching (default False)
         - useQueryQueryMatches: (optional) use query-query matching logic (default False)
         - maxMatches: (optional) maximum number of matches to return (default 1000)

        RETURNS: a list of (match, tautomer) pairs
        """

    @overload
    def GetSubstructMatchesWithTautomers(self, target: rdkit.Chem.rdchem.Mol, params: rdkit.Chem.rdchem.SubstructMatchParameters) -> list[tuple[list[int], rdkit.Chem.rdchem.Mol]]:
        """
        Return all substructure matches with their matching tautomers.

        ARGUMENTS:
         - target: the target molecule
         - params: SubstructMatchParameters object

        RETURNS: a list of (match, tautomer) pairs
        """

    def PatternFingerprintTemplate(self, fingerprintSize: int = 2048) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
        """
        Return the pattern fingerprint of the template molecule.

        ARGUMENTS:
         - fingerprintSize: (optional) size of the fingerprint (default 2048)

        RETURNS: an ExplicitBitVect fingerprint
        """

    def GetTemplateMolecule(self) -> rdkit.Chem.rdchem.Mol:
        """Return the template molecule used for substructure searching."""

    def GetModifiedAtoms(self) -> list[int]:
        """Return the indices of tautomeric atoms."""

    def GetModifiedBonds(self) -> list[int]:
        """Return the indices of tautomeric bonds."""

    def GetTautomers(self) -> list[rdkit.Chem.rdchem.Mol]:
        """Return the list of tautomers of the query molecule."""

    def ToBinary(self) -> bytes:
        """Return a binary string (pickle) representation of this TautomerQuery."""

    def ToStream(self, fileobj: object) -> None:
        """Serialize this TautomerQuery to a file-like object."""

    def InitFromStream(self, fileobj: object) -> None:
        """Initialize this TautomerQuery from a file-like object."""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

def PatternFingerprintTautomerTarget(target: rdkit.Chem.rdchem.Mol, fingerprintSize: int = 2048) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """
    Return the pattern fingerprint of a target molecule for tautomer searching.

    ARGUMENTS:
     - target: the target molecule
     - fingerprintSize: (optional) size of the fingerprint (default 2048)

    RETURNS: an ExplicitBitVect fingerprint
    """

def TautomerQueryCanSerialize() -> bool:
    """
    Returns True if the TautomerQuery is serializable
    (requires that the RDKit was built with boost::serialization)
    """
