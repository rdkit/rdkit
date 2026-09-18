"""Module containing the core chemistry functionality of the RDKit"""

from collections.abc import Callable, Iterable, Sequence
import enum
from typing import Annotated, Final, overload

import numpy
from numpy.typing import NDArray

import rdkit.DataStructs.cDataStructs
import rdkit.Geometry.rdGeometry


def tossit() -> None: ...

class MolSanitizeException(_BaseSanitException):
    pass

class AtomSanitizeException(MolSanitizeException):
    pass

class AtomValenceException(AtomSanitizeException):
    pass

class AtomKekulizeException(AtomSanitizeException):
    pass

class KekulizeException(MolSanitizeException):
    pass

class PeriodicTable:
    """
    A class which stores information from the Periodic Table.

      It is not possible to create a PeriodicTable object directly from Python,
      use GetPeriodicTable() to get the global table.

      The PeriodicTable object can be queried for a variety of properties:

        - GetAtomicWeight

        - GetAtomicNumber

        - GetElementSymbol

        - GetElementName

        - GetRow

        - GetRvdw (van der Waals radius)

        - GetRCovalent (covalent radius)

        - GetDefaultValence

        - GetValenceList

        - GetNOuterElecs (number of valence electrons)

        - GetMostCommonIsotope

        - GetMostCommonIsotopeMass

        - GetRb0

        - GetAbundanceForIsotope

        - GetMassForIsotope

      When it makes sense, these can be queried using either an atomic number (integer)
      or an atomic symbol (string)
    """

    @overload
    def GetAtomicWeight(self, atomicNumber: int) -> float: ...

    @overload
    def GetAtomicWeight(self, elementSymbol: str) -> float: ...

    def GetAtomicNumber(self, elementSymbol: str) -> int: ...

    def GetElementSymbol(self, atomicNumber: int) -> str: ...

    def GetElementName(self, atomicNumber: int) -> str: ...

    @overload
    def GetRow(self, atomicNumber: int) -> int: ...

    @overload
    def GetRow(self, elementSymbol: str) -> int: ...

    @overload
    def GetRvdw(self, atomicNumber: int) -> float: ...

    @overload
    def GetRvdw(self, elementSymbol: str) -> float: ...

    @overload
    def GetRcovalent(self, atomicNumber: int) -> float: ...

    @overload
    def GetRcovalent(self, elementSymbol: str) -> float: ...

    @overload
    def GetDefaultValence(self, atomicNumber: int) -> int: ...

    @overload
    def GetDefaultValence(self, elementSymbol: str) -> int: ...

    @overload
    def GetValenceList(self, atomicNumber: int) -> list[int]: ...

    @overload
    def GetValenceList(self, elementSymbol: str) -> list[int]: ...

    @overload
    def GetNOuterElecs(self, atomicNumber: int) -> int: ...

    @overload
    def GetNOuterElecs(self, elementSymbol: str) -> int: ...

    @overload
    def GetMostCommonIsotope(self, atomicNumber: int) -> int: ...

    @overload
    def GetMostCommonIsotope(self, elementSymbol: str) -> int: ...

    @overload
    def GetMostCommonIsotopeMass(self, atomicNumber: int) -> float: ...

    @overload
    def GetMostCommonIsotopeMass(self, elementSymbol: str) -> float: ...

    @overload
    def GetRb0(self, atomicNumber: int) -> float: ...

    @overload
    def GetRb0(self, elementSymbol: str) -> float: ...

    @overload
    def GetAbundanceForIsotope(self, atomicNumber: int, isotope: int) -> float: ...

    @overload
    def GetAbundanceForIsotope(self, elementSymbol: str, isotope: int) -> float: ...

    @overload
    def GetMassForIsotope(self, atomicNumber: int, isotope: int) -> float: ...

    @overload
    def GetMassForIsotope(self, elementSymbol: str, isotope: int) -> float: ...

    def GetMaxAtomicNumber(self) -> int: ...

def GetPeriodicTable() -> PeriodicTable:
    """Returns the application's PeriodicTable instance."""

class Atom:
    """
    The class to store Atoms.
    Note that, though it is possible to create one, having an Atom on its own
    (i.e not associated with a molecule) is not particularly useful.
    """

    @overload
    def __init__(self, what: str) -> None: ...

    @overload
    def __init__(self, other: Atom) -> None: ...

    @overload
    def __init__(self, num: int) -> None:
        """Constructor, takes the atomic number"""

    NOATOM: Final[int] = ...
    """marker for unspecified int values"""

    def __copy__(self) -> Atom:
        """Create a copy of the atom"""

    def GetAtomicNum(self) -> int:
        """Returns the atomic number."""

    def SetAtomicNum(self, newNum: int) -> None:
        """Sets the atomic number, takes an integer value as an argument"""

    def GetSymbol(self) -> str:
        """Returns the atomic symbol (a string)"""

    def GetIdx(self) -> int:
        """Returns the atom's index (ordering in the molecule)"""

    def GetDegree(self) -> int:
        """
        Returns the degree of the atom in the molecule.

          The degree of an atom is defined to be its number of
          directly-bonded neighbors.
          The degree is independent of bond orders, but is dependent
            on whether or not Hs are explicit in the graph.
        """

    def GetTotalDegree(self) -> int:
        """
        Returns the degree of the atom in the molecule including Hs.

          The degree of an atom is defined to be its number of
          directly-bonded neighbors.
          The degree is independent of bond orders.
        """

    def GetTotalNumHs(self, includeNeighbors: bool = False) -> int:
        """
        Returns the total number of Hs (explicit and implicit) on the atom.

          ARGUMENTS:

            - includeNeighbors: (optional) toggles inclusion of neighboring H atoms in the sum.
              Defaults to 0.
        """

    def GetNumImplicitHs(self) -> int:
        """Returns the total number of implicit Hs on the atom."""

    def GetExplicitValence(self) -> int:
        """
        DEPRECATED, please use GetValence(Chem.ValenceType,EXPLICIT) instead.
        Returns the explicit valence of the atom.
        """

    def GetImplicitValence(self) -> int:
        """
        DEPRECATED, please use getValence(Chem.ValenceType,IMPLICIT) instead.
        Returns the number of implicit Hs on the atom.
        """

    def GetValence(self, which: ValenceType) -> int:
        """Returns the valence (explicit or implicit) of the atom."""

    def GetTotalValence(self) -> int:
        """Returns the total valence (explicit + implicit) of the atom."""

    def HasValenceViolation(self) -> bool:
        """Returns whether the atom has a valence violation or not."""

    def GetFormalCharge(self) -> int: ...

    def SetFormalCharge(self, what: int) -> None: ...

    def SetNoImplicit(self, what: bool) -> None:
        """
        Sets a marker on the atom that *disallows* implicit Hs.
          This holds even if the atom would otherwise have implicit Hs added.
        """

    def GetNoImplicit(self) -> bool:
        """Returns whether or not the atom is *allowed* to have implicit Hs."""

    def SetNumExplicitHs(self, what: int) -> None: ...

    def GetNumExplicitHs(self) -> int: ...

    def SetIsAromatic(self, what: bool) -> None: ...

    def GetIsAromatic(self) -> bool: ...

    def GetMass(self) -> float: ...

    def SetIsotope(self, what: int) -> None: ...

    def GetIsotope(self) -> int: ...

    def SetNumRadicalElectrons(self, num: int) -> None: ...

    def GetNumRadicalElectrons(self) -> int: ...

    def GetQueryType(self) -> str: ...

    def SetChiralTag(self, what: ChiralType) -> None: ...

    def InvertChirality(self) -> bool: ...

    def GetChiralTag(self) -> ChiralType: ...

    def SetHybridization(self, what: HybridizationType) -> None:
        """
        Sets the hybridization of the atom.
          The argument should be a HybridizationType
        """

    def GetHybridization(self) -> HybridizationType:
        """Returns the atom's hybridization."""

    def HasOwningMol(self) -> bool:
        """Returns whether or not this instance belongs to a molecule."""

    def GetOwningMol(self) -> Mol:
        """Returns the Mol that owns this atom."""

    def GetNeighbors(self) -> _AtomSeqHolder2:
        """Returns a sequence-like object of the atom's neighbors."""

    def GetBonds(self) -> _BondSeqHolder2:
        """Returns a sequence-like object of the atom's bonds."""

    def Match(self, other: Atom) -> bool:
        """
        Returns whether or not this atom matches another Atom.

          Each Atom (or query Atom) has a query function which is
          used for this type of matching.

          ARGUMENTS:
            - other: the other Atom to which to compare
        """

    def IsInRingSize(self, size: int) -> bool:
        """
        Returns whether or not the atom is in a ring of a particular size.

          ARGUMENTS:
            - size: the ring size to look for
        """

    def IsInRing(self) -> bool:
        """Returns whether or not the atom is in a ring"""

    def HasQuery(self) -> bool:
        """Returns whether or not the atom has an associated query"""

    def DescribeQuery(self) -> str:
        """
        returns a text description of the query. Primarily intended for debugging purposes.
        """

    def GetSmarts(self, doKekule: bool = False, allHsExplicit: bool = False, isomericSmiles: bool = True) -> str:
        """returns the SMARTS (or SMILES) string for an Atom"""

    def SetProp(self, key: str, val: str) -> None:
        """
        Sets an atomic property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value (a string).
        """

    @overload
    def GetProp(self, key: str, autoConvert: bool = False) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - autoConvert: if True attempt to convert the property into a python object

          RETURNS: a string

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetProp(self, key: str, autoConvert: bool = False, default: object | None = None) -> object:
        """
        Returns the value of the property.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

                  - autoConvert: if True attempt to convert the property into a python object

                  - default: value to return if the property is not present.

             RETURNS: the property value, or default if the property is not present.
        """

    def SetIntProp(self, key: str, val: int) -> None:
        """
        Sets an atomic property

          ARGUMENTS:
            - key: the name of the property to be set (a int).
            - value: the property value (a int).
        """

    def SetUnsignedProp(self, key: str, val: int) -> None:
        """
        Sets an atomic property

          ARGUMENTS:
            - key: the name of the property to be set (an unsigned integer).
            - value: the property value (a int >= 0).
        """

    @overload
    def GetIntProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (an int).

          RETURNS: an int

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetIntProp(self, key: str, default: object) -> object:
        """
        Returns the value of the property.

             ARGUMENTS:
                  - key: the name of the property to return (an int).

                  - default: value to return if the property is not present.

             RETURNS: an int, or default if the property is not present.
        """

    @overload
    def GetUnsignedProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (an unsigned integer).

          RETURNS: an integer (Python has no unsigned type)

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetUnsignedProp(self, key: str, default: object) -> object:
        """
        Returns the value of the property.

             ARGUMENTS:
                  - key: the name of the property to return (an unsigned integer).

                  - default: value to return if the property is not present.

             RETURNS: an integer, or default if the property is not present.
        """

    def SetDoubleProp(self, key: str, val: float) -> None:
        """
        Sets an atomic property

          ARGUMENTS:
            - key: the name of the property to be set (a double).
            - value: the property value (a double).
        """

    @overload
    def GetDoubleProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a double).

          RETURNS: a double

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetDoubleProp(self, key: str, default: object) -> object:
        """
        Returns the value of the property.

             ARGUMENTS:
                  - key: the name of the property to return (a double).

                  - default: value to return if the property is not present.

             RETURNS: a double, or default if the property is not present.
        """

    def SetBoolProp(self, key: str, val: bool) -> None:
        """
        Sets an atomic property

          ARGUMENTS:
            - key: the name of the property to be set (a bool).
            - value: the property value (a bool).
        """

    @overload
    def GetBoolProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a bool).

          RETURNS: a bool

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetBoolProp(self, key: str, default: object) -> object:
        """
        Returns the value of the property.

             ARGUMENTS:
                  - key: the name of the property to return (a bool).

                  - default: value to return if the property is not present.

             RETURNS: a bool, or default if the property is not present.
        """

    def SetExplicitBitVectProp(self, key: str, val: rdkit.DataStructs.cDataStructs.ExplicitBitVect) -> None:
        """
        Sets an atomic property

          ARGUMENTS:
            - key: the name of the property to be set (an ExplicitBitVect).
            - value: the property value (an ExplicitBitVect).
        """

    def GetExplicitBitVectProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a ExplicitBitVect).

          RETURNS: an ExplicitBitVect 

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    def HasProp(self, key: str) -> int:
        """
        Queries a Atom to see if a particular property has been assigned.

          ARGUMENTS:
            - key: the name of the property to check for (a string).
        """

    def ClearProp(self, key: str) -> None:
        """
        Removes a particular property from an Atom (does nothing if not already set).

          ARGUMENTS:
            - key: the name of the property to be removed.
        """

    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> list[str]:
        """Returns a list of the properties set on the Atom."""

    def GetPropsAsDict(self, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict:
        """
        Returns a dictionary populated with properties.
        When possible, string values will be converted to integers or doubles (trimming if necessary)
         n.b. Some properties are not able to be converted to python types.

          ARGUMENTS:
            - includePrivate: (optional) toggles inclusion of private properties in the result set.
                              Defaults to False.
            - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                              Defaults to False.

            - autoConvertStrings: (optional) toggles automatic conversion of string properties to integers or doubles.
                              Defaults to True.

          RETURNS: a dictionary
        """

    def UpdatePropertyCache(self, strict: bool = True) -> None:
        """
        Regenerates computed properties like implicit valence and ring information.
        """

    def NeedsUpdatePropertyCache(self) -> bool:
        """
        Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.
        """

    def ClearPropertyCache(self) -> None:
        """Clears implicit and explicit valence information."""

    def GetMonomerInfo(self) -> AtomMonomerInfo:
        """Returns the atom's MonomerInfo object, if there is one."""

    def GetPDBResidueInfo(self) -> AtomPDBResidueInfo:
        """Returns the atom's MonomerInfo object, if there is one."""

    def SetMonomerInfo(self, info: AtomMonomerInfo) -> None:
        """Sets the atom's MonomerInfo object."""

    def SetPDBResidueInfo(self, info: AtomMonomerInfo) -> None:
        """Sets the atom's MonomerInfo object."""

    def GetAtomMapNum(self) -> int:
        """Gets the atoms map number, returns 0 if not set"""

    def SetAtomMapNum(self, mapno: int, strict: bool = False) -> None:
        """Sets the atoms map number, a value of 0 clears the atom map"""

class HybridizationType(enum.Enum):
    UNSPECIFIED = 0

    S = 1

    SP = 2

    SP2 = 3

    SP3 = 4

    SP2D = 5

    SP3D = 6

    SP3D2 = 7

    OTHER = 8

class ChiralType(enum.Enum):
    CHI_UNSPECIFIED = 0

    CHI_TETRAHEDRAL_CW = 1

    CHI_TETRAHEDRAL_CCW = 2

    CHI_OTHER = 3

    CHI_TETRAHEDRAL = 4

    CHI_ALLENE = 5

    CHI_SQUAREPLANAR = 6

    CHI_TRIGONALBIPYRAMIDAL = 7

    CHI_OCTAHEDRAL = 8

CHI_UNSPECIFIED: ChiralType = ChiralType.CHI_UNSPECIFIED

CHI_TETRAHEDRAL_CW: ChiralType = ChiralType.CHI_TETRAHEDRAL_CW

CHI_TETRAHEDRAL_CCW: ChiralType = ChiralType.CHI_TETRAHEDRAL_CCW

CHI_OTHER: ChiralType = ChiralType.CHI_OTHER

CHI_TETRAHEDRAL: ChiralType = ChiralType.CHI_TETRAHEDRAL

CHI_ALLENE: ChiralType = ChiralType.CHI_ALLENE

CHI_SQUAREPLANAR: ChiralType = ChiralType.CHI_SQUAREPLANAR

CHI_TRIGONALBIPYRAMIDAL: ChiralType = ChiralType.CHI_TRIGONALBIPYRAMIDAL

CHI_OCTAHEDRAL: ChiralType = ChiralType.CHI_OCTAHEDRAL

class ValenceType(enum.Enum):
    IMPLICIT = 0

    EXPLICIT = 1

IMPLICIT: ValenceType = ValenceType.IMPLICIT

EXPLICIT: ValenceType = ValenceType.EXPLICIT

class CompositeQueryType(enum.Enum):
    COMPOSITE_AND = 0

    COMPOSITE_OR = 1

    COMPOSITE_XOR = 2

COMPOSITE_AND: CompositeQueryType = CompositeQueryType.COMPOSITE_AND

COMPOSITE_OR: CompositeQueryType = CompositeQueryType.COMPOSITE_OR

COMPOSITE_XOR: CompositeQueryType = CompositeQueryType.COMPOSITE_XOR

class QueryAtom(Atom):
    """
    The class to store QueryAtoms.
    These cannot currently be constructed directly from Python
    """

    def ExpandQuery(self, other: QueryAtom, how: CompositeQueryType = CompositeQueryType.COMPOSITE_AND, maintainOrder: bool = True) -> None:
        """combines the query from other with ours"""

    def SetQuery(self, other: QueryAtom) -> None:
        """Replace our query with a copy of the other query"""

def GetAtomRLabel(atom: Atom) -> int:
    """Returns the atom's MDL AtomRLabel (this is an integer from 0 to 99)"""

def SetAtomRLabel(atom: Atom, rlabel: int) -> None:
    """
    Sets the atom's MDL RLabel (this is an integer from 0 to 99).
    Setting to 0 clears the rlabel.
    """

def GetAtomAlias(atom: Atom) -> str:
    """Returns the atom's MDL alias text"""

def SetAtomAlias(atom: Atom, rlabel: str) -> None:
    """
    Sets the atom's MDL alias text.
    Setting to an empty string clears the alias.
    """

def GetAtomValue(atom: Atom) -> str:
    """Returns the atom's MDL alias text"""

def SetAtomValue(atom: Atom, rlabel: str) -> None:
    """
    Sets the atom's MDL alias text.
    Setting to an empty string clears the alias.
    """

def GetSupplementalSmilesLabel(atom: Atom) -> str:
    """
    Gets the supplemental smiles label on an atom, returns an empty string if not present.
    """

def SetSupplementalSmilesLabel(atom: Atom, label: str) -> None:
    """
    Sets a supplemental label on an atom that is written to the smiles string.

    >>> m = Chem.MolFromSmiles("C")
    >>> Chem.SetSupplementalSmilesLabel(m.GetAtomWithIdx(0), '<xxx>')
    >>> Chem.MolToSmiles(m)
    'C<xxx>'
    """

def GetNumPiElectrons(atom: Atom) -> int:
    """Returns the number of electrons an atom is using for pi bonding"""

class Bond:
    """
    The class to store Bonds.
    Note: unlike Atoms, is it currently impossible to construct Bonds from
    Python.
    """

    def HasOwningMol(self) -> bool:
        """Returns whether or not this instance belongs to a molecule."""

    def GetOwningMol(self) -> Mol:
        """Returns the Mol that owns this bond."""

    def GetBondType(self) -> BondType:
        """Returns the type of the bond as a BondType"""

    def SetBondType(self, bT: BondType) -> None:
        """Set the type of the bond as a BondType"""

    def GetBondTypeAsDouble(self) -> float:
        """
        Returns the type of the bond as a double (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)
        """

    def GetBondDir(self) -> BondDir:
        """Returns the type of the bond as a BondDir"""

    def SetBondDir(self, what: BondDir) -> None:
        """Set the type of the bond as a BondDir"""

    def GetStereo(self) -> BondStereo:
        """Returns the stereo configuration of the bond as a BondStereo"""

    def SetStereo(self, what: BondStereo) -> None:
        """Set the stereo configuration of the bond as a BondStereo"""

    def GetStereoAtoms(self) -> list[int]:
        """Returns the indices of the atoms setting this bond's stereochemistry."""

    def SetStereoAtoms(self, bgnIdx: int, endIdx: int) -> None:
        """Set the indices of the atoms setting this bond's stereochemistry."""

    def InvertChirality(self) -> bool: ...

    def GetValenceContrib(self, at: Atom) -> float:
        """
        Returns the contribution of the bond to the valence of an Atom.

          ARGUMENTS:

            - atom: the Atom to consider.
        """

    def GetIsAromatic(self) -> bool: ...

    def SetIsAromatic(self, what: bool) -> None: ...

    def GetIsConjugated(self) -> bool:
        """Returns whether or not the bond is considered to be conjugated."""

    def SetIsConjugated(self, what: bool) -> None: ...

    def GetIdx(self) -> int:
        """Returns the bond's index (ordering in the molecule)"""

    def GetBeginAtomIdx(self) -> int:
        """Returns the index of the bond's first atom."""

    def GetEndAtomIdx(self) -> int:
        """Returns the index of the bond's first atom."""

    def GetOtherAtomIdx(self, thisIdx: int) -> int:
        """
        Given the index of one of the bond's atoms, returns the
        index of the other.
        """

    def GetBeginAtom(self) -> Atom:
        """Returns the bond's first atom."""

    def GetEndAtom(self) -> Atom:
        """Returns the bond's second atom."""

    def GetOtherAtom(self, what: Atom) -> Atom:
        """Given one of the bond's atoms, returns the other one."""

    def Match(self, what: Bond) -> bool:
        """
        Returns whether or not this bond matches another Bond.

          Each Bond (or query Bond) has a query function which is
          used for this type of matching.

          ARGUMENTS:
            - other: the other Bond to which to compare
        """

    def IsInRingSize(self, size: int) -> bool:
        """
        Returns whether or not the bond is in a ring of a particular size.

          ARGUMENTS:
            - size: the ring size to look for
        """

    def IsInRing(self) -> bool:
        """Returns whether or not the bond is in a ring of any size."""

    def HasQuery(self) -> bool:
        """Returns whether or not the bond has an associated query"""

    def DescribeQuery(self) -> str:
        """
        returns a text description of the query. Primarily intended for debugging purposes.
        """

    def GetSmarts(self, allBondsExplicit: bool = False) -> str:
        """returns the SMARTS (or SMILES) string for a Bond"""

    def SetProp(self, key: str, val: str) -> None:
        """
        Sets a bond property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value (a string).
        """

    @overload
    def GetProp(self, key: str, autoConvert: bool = False) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - autoConvert: if True attempt to convert the property into a python object

          RETURNS: a string

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetProp(self, key: str, autoConvert: bool = False, default: object | None = None) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - autoConvert: if True attempt to convert the property into a python object

            - default: value to return if the property is not present.

          RETURNS: the property value, or default if the property is not present.
        """

    def SetIntProp(self, key: str, val: int) -> None:
        """
        Sets a bond property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value (an int).
        """

    def SetUnsignedProp(self, key: str, val: int) -> None:
        """
        Sets a bond property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value (an int >= 0).
        """

    @overload
    def GetIntProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (an int).

          RETURNS: an int

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetIntProp(self, key: str, default: object) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (an int).

            - default: value to return if the property is not present.

          RETURNS: an int, or default if the property is not present.
        """

    @overload
    def GetUnsignedProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (an unsigned integer).

          RETURNS: an int (Python has no unsigned type)

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetUnsignedProp(self, key: str, default: object) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (an unsigned integer).

            - default: value to return if the property is not present.

          RETURNS: an integer, or default if the property is not present.
        """

    def SetDoubleProp(self, key: str, val: float) -> None:
        """
        Sets a bond property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value (a double).
        """

    @overload
    def GetDoubleProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a double).

          RETURNS: a double

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetDoubleProp(self, key: str, default: object) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a double).

            - default: value to return if the property is not present.

          RETURNS: a double, or default if the property is not present.
        """

    def SetBoolProp(self, key: str, val: bool) -> None:
        """
        Sets a bond property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value (a boolean).
        """

    @overload
    def GetBoolProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a boolean).

          RETURNS: a boolean

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetBoolProp(self, key: str, default: object) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a boolean).

            - default: value to return if the property is not present.

          RETURNS: a bool, or default if the property is not present.
        """

    def HasProp(self, key: str) -> int:
        """
        Queries a Bond to see if a particular property has been assigned.

          ARGUMENTS:
            - key: the name of the property to check for (a string).
        """

    def ClearProp(self, key: str) -> None:
        """
        Removes a particular property from an Bond (does nothing if not already set).

          ARGUMENTS:
            - key: the name of the property to be removed.
        """

    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> list[str]:
        """Returns a list of the properties set on the Bond."""

    def GetPropsAsDict(self, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict:
        """
        Returns a dictionary populated with properties.
        When possible, string values will be converted to integers or doubles (trimming if necessary)
         n.b. Some properties are not able to be converted to python types.

          ARGUMENTS:
            - includePrivate: (optional) toggles inclusion of private properties in the result set.
                              Defaults to False.
            - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                              Defaults to False.

            - autoConvertStrings: (optional) toggles automatic conversion of string properties to integers or doubles.
                              Defaults to True.

          RETURNS: a dictionary
        """

class BondType(enum.Enum):
    UNSPECIFIED = 0

    SINGLE = 1

    DOUBLE = 2

    TRIPLE = 3

    QUADRUPLE = 4

    QUINTUPLE = 5

    HEXTUPLE = 6

    ONEANDAHALF = 7

    TWOANDAHALF = 8

    THREEANDAHALF = 9

    FOURANDAHALF = 10

    FIVEANDAHALF = 11

    AROMATIC = 12

    IONIC = 13

    HYDROGEN = 14

    THREECENTER = 15

    DATIVEONE = 16

    DATIVE = 17

    DATIVEL = 18

    DATIVER = 19

    OTHER = 20

    ZERO = 21

class BondDir(enum.Enum):
    NONE = 0

    BEGINWEDGE = 1

    BEGINDASH = 2

    ENDDOWNRIGHT = 3

    ENDUPRIGHT = 4

    EITHERDOUBLE = 5

    UNKNOWN = 6

class BondStereo(enum.Enum):
    STEREONONE = 0

    STEREOANY = 1

    STEREOZ = 2

    STEREOE = 3

    STEREOCIS = 4

    STEREOTRANS = 5

    STEREOATROPCW = 6

    STEREOATROPCCW = 7

class QueryBond(Bond):
    """
    The class to store QueryBonds.
    These cannot currently be constructed directly from Python
    """

    def ExpandQuery(self, other: QueryBond, how: CompositeQueryType = CompositeQueryType.COMPOSITE_AND, maintainOrder: bool = True) -> None:
        """combines the query from other with ours"""

    def SetQuery(self, other: QueryBond) -> None:
        """Replace our query with a copy of the other query"""

class ConformerException(ValueError):
    pass

class PropertyPickleOptions(enum.IntEnum):
    NoProps = 0

    MolProps = 1

    AtomProps = 2

    BondProps = 4

    QueryAtomData = 2

    PrivateProps = 16

    ComputedProps = 32

    AllProps = 65535

    CoordsAsDouble = 65536

    NoConformers = 131072

NoProps: PropertyPickleOptions = PropertyPickleOptions.NoProps

MolProps: PropertyPickleOptions = PropertyPickleOptions.MolProps

AtomProps: PropertyPickleOptions = PropertyPickleOptions.AtomProps

BondProps: PropertyPickleOptions = PropertyPickleOptions.BondProps

PrivateProps: PropertyPickleOptions = PropertyPickleOptions.PrivateProps

ComputedProps: PropertyPickleOptions = PropertyPickleOptions.ComputedProps

AllProps: PropertyPickleOptions = PropertyPickleOptions.AllProps

CoordsAsDouble: PropertyPickleOptions = PropertyPickleOptions.CoordsAsDouble

NoConformers: PropertyPickleOptions = PropertyPickleOptions.NoConformers

def GetDefaultPickleProperties() -> int:
    """Get the current global mol pickler options."""

def SetDefaultPickleProperties(arg1: int) -> None:
    """Set the current global mol pickler options."""

class AtomCoordsMatcher:
    """Allows using atom coordinates as part of substructure matching"""

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, refConfId: int = -1, queryConfId: int = -1, tol: float = 0.0001) -> None:
        """
        constructor taking reference and query conformer IDs and a distance tolerance
        """

    def __call__(self, arg0: Atom, arg1: Atom, /) -> bool: ...

    @property
    def refConfId(self) -> int:
        """reference conformer ID"""

    @refConfId.setter
    def refConfId(self, arg: int, /) -> None: ...

    @property
    def queryConfId(self) -> int:
        """query conformer ID"""

    @queryConfId.setter
    def queryConfId(self, arg: int, /) -> None: ...

    @property
    def tol2(self) -> float:
        """squared distance tolerance"""

    @tol2.setter
    def tol2(self, arg: float, /) -> None: ...

class SubstructMatchParameters:
    """Parameters controlling substructure matching"""

    def __init__(self) -> None:
        """Constructor"""

    @property
    def useChirality(self) -> bool:
        """Use chirality in determining whether or not atoms/bonds match"""

    @useChirality.setter
    def useChirality(self, arg: bool, /) -> None: ...

    @property
    def useEnhancedStereo(self) -> bool:
        """
        take enhanced stereochemistry into account while doing the match. This only has an effect if useChirality is also True.
        """

    @useEnhancedStereo.setter
    def useEnhancedStereo(self, arg: bool, /) -> None: ...

    @property
    def aromaticMatchesConjugated(self) -> bool:
        """aromatic and conjugated bonds match each other"""

    @aromaticMatchesConjugated.setter
    def aromaticMatchesConjugated(self, arg: bool, /) -> None: ...

    @property
    def aromaticMatchesSingleOrDouble(self) -> bool:
        """aromatic and single or double bonds match each other"""

    @aromaticMatchesSingleOrDouble.setter
    def aromaticMatchesSingleOrDouble(self, arg: bool, /) -> None: ...

    @property
    def useGenericMatchers(self) -> bool:
        """
        use generic groups (=homology groups) as a post-filtering step (if any are present in the molecule)
        """

    @useGenericMatchers.setter
    def useGenericMatchers(self, arg: bool, /) -> None: ...

    @property
    def useQueryQueryMatches(self) -> bool:
        """Consider query-query matches, not just simple matches"""

    @useQueryQueryMatches.setter
    def useQueryQueryMatches(self, arg: bool, /) -> None: ...

    @property
    def recursionPossible(self) -> bool:
        """Allow recursive queries"""

    @recursionPossible.setter
    def recursionPossible(self, arg: bool, /) -> None: ...

    @property
    def uniquify(self) -> bool:
        """uniquify (by atom index) match results"""

    @uniquify.setter
    def uniquify(self, arg: bool, /) -> None: ...

    @property
    def maxMatches(self) -> int:
        """maximum number of matches to return"""

    @maxMatches.setter
    def maxMatches(self, arg: int, /) -> None: ...

    @property
    def maxRecursiveMatches(self) -> int:
        """maximum number of recursive matches to find"""

    @maxRecursiveMatches.setter
    def maxRecursiveMatches(self, arg: int, /) -> None: ...

    @property
    def numThreads(self) -> int:
        """
        number of threads to use when multi-threading is possible.0 selects the number of concurrent threads supported by thehardware. negative values are added to the number of concurrentthreads supported by the hardware.
        """

    @numThreads.setter
    def numThreads(self, arg: int, /) -> None: ...

    @property
    def bondProperties(self) -> list[str]:
        """bond properties that must be equivalent in order to match."""

    @bondProperties.setter
    def bondProperties(self, arg: Sequence[str], /) -> None: ...

    @property
    def atomProperties(self) -> list[str]:
        """atom properties that must be equivalent in order to match."""

    @atomProperties.setter
    def atomProperties(self, arg: Sequence[str], /) -> None: ...

    @property
    def specifiedStereoQueryMatchesUnspecified(self) -> bool:
        """
        If set, query atoms and bonds with specified stereochemistry will match atoms and bonds with unspecified stereochemistry.
        """

    @specifiedStereoQueryMatchesUnspecified.setter
    def specifiedStereoQueryMatchesUnspecified(self, arg: bool, /) -> None: ...

    def setExtraFinalCheck(self, func: Callable[[Mol, Sequence[int]], bool]) -> None:
        """
        allows you to provide a function that will be called
                       with the molecule
                   and a vector of atom IDs containing a potential match.
                   The function should return true or false indicating whether or not
                   that match should be accepted.
        """

    def setExtraAtomCheckFunc(self, func: Callable[[Atom, Atom], bool]) -> None:
        """
        allows you to provide a function that will be called
                   for each atom pair that matches during substructure searching,
                   after all other comparisons have passed.
                   The function should return true or false indicating whether or not
                   that atom-match should be accepted.
        """

    @property
    def extraAtomCheckOverridesDefaultCheck(self) -> bool:
        """
        if set, only the extraAtomCheck will be used to determine whether or not atoms match
        """

    @extraAtomCheckOverridesDefaultCheck.setter
    def extraAtomCheckOverridesDefaultCheck(self, arg: bool, /) -> None: ...

    def setExtraBondCheckFunc(self, func: Callable[[Bond, Bond], bool]) -> None:
        """
        allows you to provide a function that will be called
                   for each bond pair that matches during substructure searching,
                   after all other comparisons have passed.
                   The function should return true or false indicating whether or not
                   that bond-match should be accepted.
        """

    @property
    def extraBondCheckOverridesDefaultCheck(self) -> bool:
        """
        if set, only the extraBondCheck will be used to determine whether or not bonds match
        """

    @extraBondCheckOverridesDefaultCheck.setter
    def extraBondCheckOverridesDefaultCheck(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class Mol:
    """
    The Molecule class.

         In addition to the expected Atoms and Bonds, molecules contain:
              - a collection of Atom and Bond bookmarks indexed with integers
                        that can be used to flag and retrieve particular Atoms or Bonds
                        using the {get|set}{Atom|Bond}Bookmark() methods.

              - a set of string-valued properties. These can have arbitrary string
                        labels and can be set and retrieved using the {set|get}Prop() methods
                        Molecular properties can be tagged as being *computed*, in which case
                             they will be automatically cleared under certain circumstances (when the
                             molecule itself is modified, for example).
                        Molecules also have the concept of *private* properties, which are tagged
                             by beginning the property name with an underscore (_).
    """

    @overload
    def __init__(self) -> None:
        """Constructor, takes no arguments"""

    @overload
    def __init__(self, pklString: bytes) -> None: ...

    @overload
    def __init__(self, pklString: bytes, propertyFlags: int) -> None: ...

    @overload
    def __init__(self, pklString: str) -> None:
        """Constructor from a binary string"""

    @overload
    def __init__(self, pklString: str, propertyFlags: int) -> None:
        """Constructor from a binary string with property flags"""

    @overload
    def __init__(self, mol: Mol, quickCopy: bool = False, confId: int = -1) -> None:
        """Constructor from another molecule"""

    def __copy__(self) -> object: ...

    def __deepcopy__(self, memo: dict) -> object: ...

    def GetAtoms(self) -> _AtomSeqHolder1:
        """Returns a sequence-like object of the molecule's atoms."""

    def GetBonds(self) -> _BondSeqHolder1:
        """Returns a sequence-like object of the molecule's bonds."""

    @overload
    def GetNumAtoms(self) -> int:
        """Returns the number of atoms in the molecule."""

    @overload
    def GetNumAtoms(self, onlyExplicit: bool = True) -> int:
        """
        Returns the number of atoms in the molecule. Optionally, only count explicit atoms.
        """

    def GetNumHeavyAtoms(self) -> int:
        """Returns the number of heavy atoms (atomic number >1) in the molecule."""

    def GetAtomWithIdx(self, idx: int) -> Atom:
        """
        Returns a particular Atom.

             ARGUMENTS:
                  - idx: which Atom to return

             NOTE: atom indices start at 0
        """

    def GetNumBonds(self, onlyHeavy: bool = True) -> int:
        """
        Returns the number of Bonds in the molecule.

             ARGUMENTS:
                  - onlyHeavy: (optional) include only bonds to heavy atoms (not Hs)
                                                     defaults to True.
        """

    def GetBondWithIdx(self, idx: int) -> Bond:
        """
        Returns a particular Bond.

             ARGUMENTS:
                  - idx: which Bond to return

             NOTE: bond indices start at 0
        """

    def GetNumConformers(self) -> int:
        """Return the number of conformations on the molecule"""

    def AddConformer(self, conf: Conformer, assignId: bool = False) -> int:
        """Add a conformer to the molecule and return the conformer ID"""

    def GetConformer(self, id: int = -1) -> Conformer:
        """Get the conformer with a specified ID"""

    def GetConformers(self) -> _ROConformerSeq:
        """
        Returns a read-only sequence containing all of the molecule's Conformers.
        """

    def RemoveAllConformers(self) -> None:
        """Remove all the conformations on the molecule"""

    def RemoveConformer(self, id: int) -> None:
        """Remove the conformer with the specified ID"""

    def GetBondBetweenAtoms(self, idx1: int, idx2: int) -> Bond | None:
        """
        Returns the bond between two atoms, if there is one.

             ARGUMENTS:
                  - idx1,idx2: the Atom indices

             Returns:
                  The Bond between the two atoms, if such a bond exists.
                  If there is no Bond between the atoms, None is returned instead.

             NOTE: atom indices start at 0
        """

    def HasQuery(self) -> bool:
        """Returns if any atom or bond in molecule has a query"""

    @overload
    def HasSubstructMatch(self, query: Mol, params: SubstructMatchParameters | None = None) -> bool:
        """
        Queries whether or not the molecule contains a particular substructure.

                  ARGUMENTS:
                  - query: a Molecule

                  - params: parameters controlling the substructure match

                  RETURNS: True or False
        """

    @overload
    def HasSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters | None = None) -> bool: ...

    @overload
    def HasSubstructMatch(self, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
        Queries whether or not the molecule contains a particular substructure.

             ARGUMENTS:
                  - query: a Molecule

                  - recursionPossible: (optional)

                  - useChirality: enables the use of stereochemistry in the matching

                  - useQueryQueryMatches: use query-query matching logic

             RETURNS: True or False
        """

    @overload
    def HasSubstructMatch(self, query: MolBundle, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool: ...

    @overload
    def GetSubstructMatch(self, query: Mol, params: SubstructMatchParameters | None = None) -> list[int]:
        """
        Returns the indices of the molecule's atoms that match a substructure query.

          ARGUMENTS:
            - query: a Molecule

            - params: parameters controlling the substructure match

          RETURNS: a list of integers

          NOTES:
             - only a single match is returned
             - the ordering of the indices corresponds to the atom ordering
                 in the query. For example, the first index is for the atom in
                 this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters | None = None) -> list[int]: ...

    @overload
    def GetSubstructMatch(self, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> list[int]:
        """
        Returns the indices of the molecule's atoms that match a substructure query.

             ARGUMENTS:
                  - query: a Molecule

                  - useChirality: enables the use of stereochemistry in the matching

                  - useQueryQueryMatches: use query-query matching logic

             RETURNS: a list of integers

             NOTES:
                   - only a single match is returned
                   - the ordering of the indices corresponds to the atom ordering
                             in the query. For example, the first index is for the atom in
                             this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatch(self, query: MolBundle, useChirality: bool = False, useQueryQueryMatches: bool = False) -> list[int]: ...

    @overload
    def GetSubstructMatches(self, query: Mol, params: SubstructMatchParameters | None = None) -> list[list[int]]:
        """
        Returns lists of the indices of the molecule's atoms that match a substructure query.

          ARGUMENTS:
            - query: a Molecule.

            - params: parameters controlling the substructure match

          RETURNS: a list of lists of integers

          NOTE:
             - the ordering of the indices corresponds to the atom ordering
                 in the query. For example, the first index is for the atom in
                 this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatches(self, query: MolBundle, params: SubstructMatchParameters | None = None) -> list[list[int]]: ...

    @overload
    def GetSubstructMatches(self, query: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> list[list[int]]:
        """
        Returns lists of the indices of the molecule's
                            atoms that match a substructure query.

                     ARGUMENTS:
                          - query: a Molecule.
                          - uniquify: (optional) determines whether or not the
                          matches are uniquified.
                                                        Defaults to 1.

                          - useChirality: enables the use of stereochemistry in the
                          matching

                          - useQueryQueryMatches: use query-query matching logic

                          - maxMatches: The maximum number of matches that will be
                          returned.
                                                             In high-symmetry cases
                                                             with medium-sized
                                                             molecules, it is very
                                                             easy to end up with a
                                                             combinatorial explosion
                                                             in the number of
                                                             possible matches. This
                                                             argument prevents that
                                                             from having unintended
                                                             consequences

                     RETURNS: a list of lists of integers

                     NOTE:
                           - the ordering of the indices corresponds to the atom
                           ordering
                                     in the query. For example, the first index is
                                     for the atom in this molecule that matches the
                                     first atom in the query.
        """

    @overload
    def GetSubstructMatches(self, query: MolBundle, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> list[list[int]]: ...

    def SetProp(self, key: str, val: str, computed: bool = False) -> None:
        """
        Sets a molecular property

             ARGUMENTS:
                  - key: the name of the property to be set (a string).
                  - value: the property value (a string).
                  - computed: (optional) marks the property as being computed.
                                                Defaults to False.
        """

    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None:
        """
        Sets a double valued molecular property

             ARGUMENTS:
                  - key: the name of the property to be set (a string).
                  - value: the property value as a double.
                  - computed: (optional) marks the property as being computed.
                                                Defaults to 0.
        """

    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None:
        """
        Sets an integer valued molecular property

             ARGUMENTS:
                  - key: the name of the property to be set (an unsigned number).
                  - value: the property value as an integer.
                  - computed: (optional) marks the property as being computed.
                                                Defaults to False.
        """

    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None:
        """
        Sets an unsigned integer valued molecular property

             ARGUMENTS:
                  - key: the name of the property to be set (a string).
                  - value: the property value as an unsigned integer.
                  - computed: (optional) marks the property as being computed.
                                                Defaults to False.
        """

    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None:
        """
        Sets a boolean valued molecular property

             ARGUMENTS:
                  - key: the name of the property to be set (a string).
                  - value: the property value as a bool.
                  - computed: (optional) marks the property as being computed.
                                                Defaults to False.
        """

    def HasProp(self, key: str) -> int:
        """
        Queries a molecule to see if a particular property has been assigned.

             ARGUMENTS:
                  - key: the name of the property to check for (a string).
        """

    def SetName(self, name: str) -> None:
        """
        Sets the molecule name stored as the _Name property.
             ARGUMENTS:
                  - name: the name to set (a string).
        """

    def GetName(self) -> str:
        """
        Returns the molecule name stored as the _Name property.
             NOTE:
                  - If the _Name property has not been set, an empty string will be returned.
        """

    @overload
    def GetProp(self, key: str, autoConvert: bool = False) -> object:
        """
        Returns the value of the property.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

                  - autoConvert: if True attempt to convert the property into a python object

             RETURNS: a string

             NOTE:
                  - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetProp(self, key: str, autoConvert: bool = False, default: object | None = None) -> object:
        """
        Returns the value of the property.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

                  - autoConvert: if True attempt to convert the property into a python object

                  - default: value to return if the property is not present.

             RETURNS: the property value, or default if the property is not present.
        """

    @overload
    def GetDoubleProp(self, key: str) -> object:
        """
        Returns the double value of the property if possible.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

             RETURNS: a double

             NOTE:
                  - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetDoubleProp(self, key: str, default: object) -> object:
        """
        Returns the double value of the property if possible.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

                  - default: value to return if the property is not present.

             RETURNS: a double, or default if the property is not present.
        """

    @overload
    def GetIntProp(self, key: str) -> object:
        """
        Returns the integer value of the property if possible.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

             RETURNS: an integer

             NOTE:
                  - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetIntProp(self, key: str, default: object) -> object:
        """
        Returns the integer value of the property if possible.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

                  - default: value to return if the property is not present.

             RETURNS: an integer, or default if the property is not present.
        """

    @overload
    def GetUnsignedProp(self, key: str) -> object:
        """
        Returns the unsigned int value of the property if possible.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

             RETURNS: an unsigned integer

             NOTE:
                  - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetUnsignedProp(self, key: str, default: object) -> object:
        """
        Returns the unsigned int value of the property if possible.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

                  - default: value to return if the property is not present.

             RETURNS: an unsigned integer, or default if the property is not present.
        """

    @overload
    def GetBoolProp(self, key: str) -> object:
        """
        Returns the Bool value of the property if possible.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

             RETURNS: a bool

             NOTE:
                  - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetBoolProp(self, key: str, default: object) -> object:
        """
        Returns the Bool value of the property if possible.

             ARGUMENTS:
                  - key: the name of the property to return (a string).

                  - default: value to return if the property is not present.

             RETURNS: a bool, or default if the property is not present.
        """

    def ClearProp(self, key: str) -> None:
        """
        Removes a property from the molecule.

             ARGUMENTS:
                  - key: the name of the property to clear (a string).
        """

    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> list[str]:
        """
        Returns a tuple with all property names for this molecule.

             ARGUMENTS:
                  - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                                               Defaults to 0.
                  - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                                               Defaults to 0.

             RETURNS: a tuple of strings
        """

    def GetPropsAsDict(self, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict:
        """
        Returns a dictionary populated with properties.
        When possible, string values will be converted to integers or doubles (trimming if necessary)
         n.b. Some properties are not able to be converted to python types.

          ARGUMENTS:
            - includePrivate: (optional) toggles inclusion of private properties in the result set.
                              Defaults to False.
            - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                              Defaults to False.

            - autoConvertStrings: (optional) toggles automatic conversion of string properties to integers or doubles.
                              Defaults to True.

          RETURNS: a dictionary
        """

    def GetStereoGroups(self) -> list[StereoGroup]:
        """
        Returns a list of StereoGroups defining the relative stereochemistry of the atoms.)
        """

    def GetAromaticAtoms(self) -> _ROQAtomSeq:
        """
        Returns a read-only sequence containing all of the molecule's aromatic Atoms.
        """

    def GetAtomsMatchingQuery(self, qa: QueryAtom) -> _ROQAtomSeq:
        """
        Returns a read-only sequence containing all of the atoms in a molecule that match the query atom.
                          Atom query options are defined in the rdkit.Chem.rdqueries module.
        """

    def ClearComputedProps(self, includeRings: bool = True) -> None:
        """Removes all computed properties from the molecule."""

    def UpdatePropertyCache(self, strict: bool = True) -> None:
        """
        Regenerates computed properties like implicit valence and ring information.
        """

    def NeedsUpdatePropertyCache(self) -> bool:
        """
        Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.
        """

    def ClearPropertyCache(self) -> None:
        """Clears implicit and explicit valence information from all atoms."""

    def Debug(self, useStdout: bool = True) -> None:
        """Prints debugging information about the molecule."""

    @overload
    def ToBinary(self) -> bytes:
        """Returns a binary string representation of the molecule."""

    @overload
    def ToBinary(self, propertyFlags: int) -> bytes:
        """
        Returns a binary string representation of the molecule pickling the specified properties.
        """

    def GetRingInfo(self) -> RingInfo:
        """Returns the number of molecule's RingInfo object."""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

class RWMol(Mol):
    """
    The RW molecule class (read/write)

         This class is a more-performant version of the EditableMolecule class in that
         it is a 'live' molecule and shares the interface from the Mol class.
         All changes are performed without the need to create a copy of the
         molecule using GetMol() (this is still available, however).

         n.b. Eventually this class may become a direct replacement for EditableMol
    """

    @overload
    def __init__(self) -> None:
        """Constructor, takes no arguments"""

    @overload
    def __init__(self, pklString: bytes) -> None: ...

    @overload
    def __init__(self, pklString: bytes, propertyFlags: int) -> None: ...

    @overload
    def __init__(self, pklString: str) -> None:
        """Constructor from a binary string"""

    @overload
    def __init__(self, pklString: str, propertyFlags: int) -> None:
        """Constructor from a binary string with property flags"""

    @overload
    def __init__(self, mol: Mol, quickCopy: bool = False, confId: int = -1) -> None:
        """Constructor from an ROMol"""

    def __copy__(self) -> object: ...

    def __deepcopy__(self, memo: dict) -> object: ...

    def __enter__(self) -> RWMol: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def RemoveAtom(self, idx: int) -> None:
        """Remove the specified atom from the molecule"""

    def RemoveBond(self, idx1: int, idx2: int) -> None:
        """Remove the specified bond from the molecule"""

    def AddBond(self, beginAtomIdx: int, endAtomIdx: int, order: BondType = BondType.UNSPECIFIED) -> int:
        """add a bond, returns the new number of bonds"""

    def AddAtom(self, atom: Atom) -> int:
        """add an atom, returns the index of the newly added atom"""

    def ReplaceAtom(self, index: int, newAtom: Atom, updateLabel: bool = False, preserveProps: bool = False) -> None:
        """
        replaces the specified atom with the provided one
        If updateLabel is True, the new atom becomes the active atom
        If preserveProps is True preserve keep the existing props unless explicit set on the new atom
        """

    def ReplaceBond(self, index: int, newBond: Bond, preserveProps: bool = False, keepSGroups: bool = True) -> None:
        """
        replaces the specified bond with the provided one.
        If preserveProps is True preserve keep the existing props unless explicit set on the new bond. If keepSGroups is False, allSubstance Groups referencing the bond will be dropped.
        """

    def GetMol(self) -> Mol:
        """Returns a Mol (a normal molecule)"""

    def SetStereoGroups(self, stereo_groups: Iterable[StereoGroup]) -> None:
        """Set the stereo groups"""

    def InsertMol(self, mol: Mol) -> None:
        """Insert (add) the given molecule into this one"""

    def BeginBatchEdit(self) -> None:
        """starts batch editing"""

    def RollbackBatchEdit(self) -> None:
        """cancels batch editing"""

    def CommitBatchEdit(self) -> None:
        """finishes batch editing and makes the actual changes"""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

class Conformer:
    """The class to store 2D or 3D conformation of a molecule"""

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, numAtoms: int) -> None:
        """Constructor with the number of atoms specified"""

    @overload
    def __init__(self, other: Conformer) -> None: ...

    def GetNumAtoms(self) -> int:
        """Get the number of atoms in the conformer"""

    def HasOwningMol(self) -> bool:
        """Returns whether or not this instance belongs to a molecule."""

    def GetOwningMol(self) -> Mol:
        """Get the owning molecule"""

    def GetId(self) -> int:
        """Get the ID of the conformer"""

    def SetId(self, id: int) -> None:
        """Set the ID of the conformer"""

    def GetAtomPosition(self, aid: int) -> rdkit.Geometry.rdGeometry.Point3D:
        """Get the position of an atom"""

    def GetPositions(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
        """Get positions of all the atoms"""

    def SetPositions(self, positions: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C', writable=False)]) -> None:
        """
        Set positions of all the atoms given a 2D or 3D numpy array of type double
        """

    @overload
    def SetAtomPosition(self, aid: int, loc: Iterable[float]) -> None: ...

    @overload
    def SetAtomPosition(self, atomId: int, position: rdkit.Geometry.rdGeometry.Point3D) -> None:
        """Set the position of the specified atom"""

    def Set3D(self, v: bool) -> None:
        """Set the 3D flag of the conformer"""

    def Is3D(self) -> bool:
        """returns the 3D flag of the conformer"""

    def SetProp(self, key: str, val: str, computed: bool = False) -> None:
        """
        Sets a molecular property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value (a string).
            - computed: (optional) marks the property as being computed.
                        Defaults to False.
        """

    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None:
        """
        Sets a double valued molecular property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value as a double.
            - computed: (optional) marks the property as being computed.
                        Defaults to 0.
        """

    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None:
        """
        Sets an integer valued molecular property

          ARGUMENTS:
            - key: the name of the property to be set (an unsigned number).
            - value: the property value as an integer.
            - computed: (optional) marks the property as being computed.
                        Defaults to False.
        """

    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None:
        """
        Sets an unsigned integer valued molecular property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value as an unsigned integer.
            - computed: (optional) marks the property as being computed.
                        Defaults to False.
        """

    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None:
        """
        Sets a boolean valued molecular property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value as a bool.
            - computed: (optional) marks the property as being computed.
                        Defaults to False.
        """

    def HasProp(self, key: str) -> int:
        """
        Queries a conformer to see if a particular property has been assigned.

          ARGUMENTS:
            - key: the name of the property to check for (a string).
        """

    @overload
    def GetProp(self, key: str, autoConvert: bool = False) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - autoConvert: if True attempt to convert the property into a python object

          RETURNS: a string

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetProp(self, key: str, autoConvert: bool = False, default: object | None = None) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - autoConvert: if True attempt to convert the property into a python object

            - default: value to return if the property is not present.

          RETURNS: the property value, or default if the property is not present.
        """

    @overload
    def GetDoubleProp(self, key: str) -> object:
        """
        Returns the double value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

          RETURNS: a double

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetDoubleProp(self, key: str, default: object) -> object:
        """
        Returns the double value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - default: value to return if the property is not present.

          RETURNS: a double, or default if the property is not present.
        """

    @overload
    def GetIntProp(self, key: str) -> object:
        """
        Returns the integer value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

          RETURNS: an integer

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetIntProp(self, key: str, default: object) -> object:
        """
        Returns the integer value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - default: value to return if the property is not present.

          RETURNS: an integer, or default if the property is not present.
        """

    @overload
    def GetUnsignedProp(self, key: str) -> object:
        """
        Returns the unsigned int value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

          RETURNS: an unsigned integer

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetUnsignedProp(self, key: str, default: object) -> object:
        """
        Returns the unsigned int value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - default: value to return if the property is not present.

          RETURNS: an unsigned integer, or default if the property is not present.
        """

    @overload
    def GetBoolProp(self, key: str) -> object:
        """
        Returns the Bool value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

          RETURNS: a bool

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetBoolProp(self, key: str, default: object) -> object:
        """
        Returns the Bool value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - default: value to return if the property is not present.

          RETURNS: a bool, or default if the property is not present.
        """

    def ClearProp(self, key: str) -> None:
        """
        Removes a property from the conformer.

          ARGUMENTS:
            - key: the name of the property to clear (a string).
        """

    def ClearComputedProps(self) -> None:
        """Removes all computed properties from the conformer."""

    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> list[str]:
        """
        Returns a tuple with all property names for this conformer.

          ARGUMENTS:
            - includePrivate: (optional) toggles inclusion of private properties in the result set.
                              Defaults to 0.
            - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                              Defaults to 0.

          RETURNS: a tuple of strings
        """

    def GetPropsAsDict(self, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict:
        """
        Returns a dictionary populated with properties.
        When possible, string values will be converted to integers or doubles (trimming if necessary)
         n.b. Some properties are not able to be converted to python types.

          ARGUMENTS:
            - includePrivate: (optional) toggles inclusion of private properties in the result set.
                              Defaults to False.
            - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                              Defaults to False.

            - autoConvertStrings: (optional) toggles automatic conversion of string properties to integers or doubles.
                              Defaults to True.

          RETURNS: a dictionary
        """

class StereoGroupType(enum.Enum):
    STEREO_ABSOLUTE = 0

    STEREO_OR = 1

    STEREO_AND = 2

STEREO_ABSOLUTE: StereoGroupType = StereoGroupType.STEREO_ABSOLUTE

STEREO_OR: StereoGroupType = StereoGroupType.STEREO_OR

STEREO_AND: StereoGroupType = StereoGroupType.STEREO_AND

class StereoGroup:
    """
    A collection of atoms with a defined stereochemical relationship.

    Used to help represent a sample with unknown stereochemistry, or that is a mix
    of diastereomers.
    """

    def GetGroupType(self) -> StereoGroupType:
        """Returns the StereoGroupType."""

    def GetAtoms(self) -> tuple[Atom, ...]:
        """access the atoms in the StereoGroup."""

    def GetBonds(self) -> tuple[Bond, ...]:
        """access the bonds in the StereoGroup."""

    def GetReadId(self) -> int:
        """
        return the StereoGroup's original ID.
        Note that the ID only makes sense for AND/OR groups.
        """

    def GetWriteId(self) -> int:
        """
        return the StereoGroup's ID that will be exported.
        Note that the ID only makes sense for AND/OR groups.
        """

    def SetWriteId(self, id: int) -> None:
        """
        return the StereoGroup's ID that will be exported.
        Note that the ID only makes sense for AND/OR groups.
        """

def CreateStereoGroup(stereoGroupType: StereoGroupType, mol: Mol, atomIds: Iterable[int] = [], bondIds: Iterable[int] = [], readId: int = 0) -> StereoGroup:
    """
    creates a StereoGroup associated with a molecule from a list of atom Ids
    """

def ForwardStereoGroupIds(mol: Mol) -> None:
    """Forward the original Stereo Group IDs when exporting the Mol."""

class EditableMol:
    """
    The EditableMol class.

       This class can be used to add/remove bonds and atoms to
       a molecule.
       In order to use it, you need to first construct an EditableMol
       from a standard Mol:

       >>> m = Chem.MolFromSmiles('CCC')
       >>> em = Chem.EditableMol(m)
       >>> em.AddAtom(Chem.Atom(8))
       >>> em.AddBond(0,3,Chem.BondType.SINGLE)
       >>> m2 = em.GetMol()
       >>> Chem.SanitizeMol(m2)
       >>> Chem.MolToSmiles(m2)
       'CCCO'

       *Note*: It is very, very easy to shoot yourself in the foot with
         this class by constructing an unreasonable molecule.
    """

    def __init__(self, m: Mol) -> None:
        """Construct from a Mol"""

    def RemoveAtom(self, idx: int) -> None:
        """Remove the specified atom from the molecule"""

    def RemoveBond(self, idx1: int, idx2: int) -> None:
        """Remove the specified bond from the molecule"""

    def AddBond(self, beginAtomIdx: int, endAtomIdx: int, order: BondType = BondType.UNSPECIFIED) -> int:
        """add a bond, returns the total number of bonds"""

    def AddAtom(self, atom: Atom) -> int:
        """add an atom, returns the index of the newly added atom"""

    def ReplaceAtom(self, index: int, newAtom: Atom, updateLabel: bool = False, preserveProps: bool = False) -> None:
        """
        replaces the specified atom with the provided one
        If updateLabel is True, the new atom becomes the active atom
        If preserveProps is True preserve keep the existing props unless explicit set on the new atom
        """

    def ReplaceBond(self, index: int, newBond: Bond, preserveProps: bool = False) -> None:
        """
        replaces the specified bond with the provided one.
        If preserveProps is True preserve keep the existing props unless explicit set on the new bond
        """

    def BeginBatchEdit(self) -> None:
        """starts batch editing"""

    def RollbackBatchEdit(self) -> None:
        """cancels batch editing"""

    def CommitBatchEdit(self) -> None:
        """finishes batch editing and makes the actual edits"""

    def GetMol(self) -> Mol:
        """Returns a Mol (a normal molecule)"""

class RingInfo:
    """contains information about a molecule's rings"""

    def IsAtomInRingOfSize(self, idx: int, size: int) -> bool: ...

    def MinAtomRingSize(self, idx: int) -> int: ...

    def AreAtomsInSameRing(self, idx1: int, idx2: int) -> bool: ...

    def AreAtomsInSameRingOfSize(self, idx1: int, idx2: int, size: int) -> bool: ...

    def IsBondInRingOfSize(self, idx: int, size: int) -> bool: ...

    def MinBondRingSize(self, idx: int) -> int: ...

    def AreBondsInSameRing(self, idx1: int, idx2: int) -> bool: ...

    def AreBondsInSameRingOfSize(self, idx1: int, idx2: int, size: int) -> bool: ...

    def NumAtomRings(self, idx: int) -> int: ...

    def NumBondRings(self, idx: int) -> int: ...

    def NumRings(self) -> int: ...

    def IsRingFused(self, ringIdx: int) -> bool: ...

    def AreRingsFused(self, ring1Idx: int, ring2Idx: int) -> bool: ...

    def NumFusedBonds(self, ringIdx: int) -> int: ...

    def AtomRings(self) -> tuple[tuple[int, ...], ...]: ...

    def BondRings(self) -> tuple[tuple[int, ...], ...]: ...

    def AtomMembers(self, idx: int) -> tuple[int, ...]: ...

    def BondMembers(self, idx: int) -> tuple[int, ...]: ...

    def AtomRingSizes(self, idx: int) -> tuple[int, ...]: ...

    def BondRingSizes(self, idx: int) -> tuple[int, ...]: ...

    def NumRingFamilies(self) -> int: ...

    def NumRelevantCycles(self) -> int: ...

    def AtomRingFamilies(self) -> tuple[tuple[int, ...], ...]: ...

    def BondRingFamilies(self) -> tuple[tuple[int, ...], ...]: ...

    def AreRingFamiliesInitialized(self) -> bool: ...

    def AddRing(self, atomIds: Iterable[int], bondIds: Iterable[int]) -> None:
        """Adds a ring to the set. Be very careful with this operation."""

class AtomMonomerInfo:
    """The class to store monomer information attached to Atoms"""

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, type: AtomMonomerType, name: str = '', residueName: str = '', resNum: int = 0, chainId: str = '', monomerClass: str = '') -> None: ...

    def GetName(self) -> str: ...

    def GetMonomerType(self) -> AtomMonomerType: ...

    def GetResidueName(self) -> str: ...

    def GetResidueNumber(self) -> int: ...

    def GetChainId(self) -> str: ...

    def GetMonomerClass(self) -> str: ...

    def SetName(self, nm: str) -> None: ...

    def SetMonomerType(self, typ: AtomMonomerType) -> None: ...

    def SetResidueName(self, val: str) -> None: ...

    def SetResidueNumber(self, val: int) -> None: ...

    def SetChainId(self, val: str) -> None: ...

    def SetMonomerClass(self, val: str) -> None: ...

class AtomMonomerType(enum.Enum):
    UNKNOWN = 0

    PDBRESIDUE = 1

    OTHER = 2

UNKNOWN: AtomMonomerType = AtomMonomerType.UNKNOWN

PDBRESIDUE: AtomMonomerType = AtomMonomerType.PDBRESIDUE

OTHER: AtomMonomerType = AtomMonomerType.OTHER

class AtomPDBResidueInfo(AtomMonomerInfo):
    """The class to store PDB residue information attached to Atoms"""

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, atomName: str, serialNumber: int = 1, altLoc: str = '', residueName: str = '', residueNumber: int = 0, chainId: str = '', insertionCode: str = '', occupancy: float = 1.0, tempFactor: float = 0.0, isHeteroAtom: bool = False, secondaryStructure: int = 0, segmentNumber: int = 0, monomerClass: str = '') -> None: ...

    def GetSerialNumber(self) -> int: ...

    def GetAltLoc(self) -> str: ...

    def GetResidueName(self) -> str: ...

    def GetResidueNumber(self) -> int: ...

    def GetChainId(self) -> str: ...

    def GetInsertionCode(self) -> str: ...

    def GetOccupancy(self) -> float: ...

    def GetTempFactor(self) -> float: ...

    def GetIsHeteroAtom(self) -> bool: ...

    def GetSecondaryStructure(self) -> int: ...

    def GetSegmentNumber(self) -> int: ...

    def GetMonomerClass(self) -> str: ...

    def SetSerialNumber(self, val: int) -> None: ...

    def SetAltLoc(self, val: str) -> None: ...

    def SetResidueName(self, val: str) -> None: ...

    def SetResidueNumber(self, val: int) -> None: ...

    def SetChainId(self, val: str) -> None: ...

    def SetInsertionCode(self, val: str) -> None: ...

    def SetOccupancy(self, val: float) -> None: ...

    def SetTempFactor(self, val: float) -> None: ...

    def SetIsHeteroAtom(self, val: bool) -> None: ...

    def SetSecondaryStructure(self, val: int) -> None: ...

    def SetSegmentNumber(self, val: int) -> None: ...

    def SetMonomerClass(self, val: str) -> None: ...

class ResonanceFlags(enum.IntEnum):
    ALLOW_INCOMPLETE_OCTETS = 1

    ALLOW_CHARGE_SEPARATION = 2

    KEKULE_ALL = 4

    UNCONSTRAINED_CATIONS = 8

    UNCONSTRAINED_ANIONS = 16

ALLOW_INCOMPLETE_OCTETS: ResonanceFlags = ResonanceFlags.ALLOW_INCOMPLETE_OCTETS

ALLOW_CHARGE_SEPARATION: ResonanceFlags = ResonanceFlags.ALLOW_CHARGE_SEPARATION

KEKULE_ALL: ResonanceFlags = ResonanceFlags.KEKULE_ALL

UNCONSTRAINED_CATIONS: ResonanceFlags = ResonanceFlags.UNCONSTRAINED_CATIONS

UNCONSTRAINED_ANIONS: ResonanceFlags = ResonanceFlags.UNCONSTRAINED_ANIONS

class ResonanceMolSupplierCallback:
    """
    Create a derived class from this abstract base class and
              implement the __call__() method.
              The __call__() method is called at each iteration of the
              algorithm, and provides a mechanism to monitor or stop
              its progress.

              To have your callback called, pass an instance of your
              derived class to ResonanceMolSupplier.SetProgressCallback()
    """

    def __init__(self) -> None: ...

    def GetNumConjGrps(self) -> int:
        """Returns the number of individual conjugated groups in the molecule."""

    def GetMaxStructures(self) -> int:
        """Get the number of conjugated groups this molecule has."""

    def GetNumStructures(self, conjGrpIdx: int) -> int:
        """
        Get the number of resonance structures generated so far for the passed conjugated group index.
        """

    def GetNumDiverseStructures(self, conjGrpIdx: int) -> int:
        """
        Get the number of non-degenrate resonance structures generated so far for the passed conjugated group index.
        """

    def __call__(self) -> bool:
        """
        This must be implemented in the derived class. Return True if the resonance structure generation should continue; False if the resonance structure generation should stop.
        """

class ResonanceMolSupplier:
    """
    A class which supplies resonance structures (as mols) from a mol.

         Usage examples:

              1) Lazy evaluation: the resonance structures are not constructed
                    until we ask for them:

                    >>> suppl = ResonanceMolSupplier(mol)
                    >>> for resMol in suppl:
                    ...    resMol.GetNumAtoms()

              2) Lazy evaluation 2:

                    >>> suppl = ResonanceMolSupplier(mol)
                    >>> resMol1 = next(suppl)
                    >>> resMol2 = next(suppl)
                    >>> suppl.reset()
                    >>> resMol3 = next(suppl)
                    # resMol3 and resMol1 are the same:
                    >>> MolToSmiles(resMol3)==MolToSmiles(resMol1)

              3) Random Access:

                    >>> suppl = ResonanceMolSupplier(mol)
                    >>> resMol1 = suppl[0]
                    >>> resMol2 = suppl[1]

                    NOTE: this will generate an IndexError if the supplier doesn't have that many
                    molecules.

              4) Random Access 2: looping over all resonance structures
                    >>> suppl = ResonanceMolSupplier(mol)
                    >>> nResMols = len(suppl)
                    >>> for i in range(nResMols):
                    ...   suppl[i].GetNumAtoms()
    """

    def __init__(self, mol: Mol, flags: int = 0, maxStructs: int = 1000) -> None: ...

    def __iter__(self) -> ResonanceMolSupplier: ...

    def __next__(self) -> Mol | None:
        """
        Returns the next resonance structure in the supplier. Raises _StopIteration_ on end.
        """

    def __getitem__(self, idx: int) -> Mol | None: ...

    def reset(self) -> None:
        """
        Resets our position in the resonance structure supplier to the beginning.
        """

    def __len__(self) -> int: ...

    def atEnd(self) -> bool:
        """
        Returns whether or not we have hit the end of the resonance structure supplier.
        """

    def GetNumConjGrps(self) -> int:
        """Returns the number of individual conjugated groups in the molecule."""

    def GetBondConjGrpIdx(self, bi: int) -> int:
        """
        Given a bond index, it returns the index of the conjugated groupthe bond belongs to, or -1 if it is not conjugated.
        """

    def GetAtomConjGrpIdx(self, ai: int) -> int:
        """
        Given an atom index, it returns the index of the conjugated groupthe atom belongs to, or -1 if it is not conjugated.
        """

    def SetNumThreads(self, numThreads: int) -> None:
        """
        Sets the number of threads to be used to enumerate resonance
        structures (defaults to 1; 0 selects the number of concurrent
        threads supported by the hardware; negative values are added
        to the number of concurrent threads supported by the hardware).
        """

    def SetProgressCallback(self, callback: object | None) -> None:
        """
        Pass an instance of a class derived from
        ResonanceMolSupplierCallback, which must implement the
        __call__() method.
        """

    def GetProgressCallback(self) -> object:
        """
        Get the ResonanceMolSupplierCallback subclass instance,
        or None if none was set.
        """

    def WasCanceled(self) -> bool:
        """Returns True if the resonance structure generation was canceled."""

    def Enumerate(self) -> None:
        """
        Ask ResonanceMolSupplier to enumerate resonance structures(automatically done as soon as any attempt to access them is made).
        """

    def GetIsEnumerated(self) -> bool:
        """Returns true if resonance structure enumeration has already happened."""

    def GetSubstructMatch(self, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> list[int]:
        """
        Returns the indices of the molecule's atoms that match a substructure query,
        taking into account all resonance structures in ResonanceMolSupplier.

          ARGUMENTS:
            - query: a Molecule

            - useChirality: enables the use of stereochemistry in the matching

            - useQueryQueryMatches: use query-query matching logic

          RETURNS: a tuple of integers

          NOTES:
             - only a single match is returned
             - the ordering of the indices corresponds to the atom ordering
                 in the query. For example, the first index is for the atom in
                 this molecule that matches the first atom in the query.
        """

    def GetSubstructMatches(self, query: Mol, uniquify: bool = False, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000, numThreads: int = 1) -> list[list[int]]:
        """
        Returns tuples of the indices of the molecule's atoms that match a substructure query,
        taking into account all resonance structures in ResonanceMolSupplier.

          ARGUMENTS:
            - query: a Molecule.
            - uniquify: (optional) determines whether or not the matches are uniquified.
                        Defaults to 1.

            - useChirality: enables the use of stereochemistry in the matching

            - useQueryQueryMatches: use query-query matching logic

            - maxMatches: The maximum number of matches that will be returned.
                          In high-symmetry cases with medium-sized molecules, it is
                          very easy to end up with a combinatorial explosion in the
                          number of possible matches. This argument prevents that from
                          having unintended consequences

            - numThreads: The number of threads to be used (defaults to 1; 0 selects the
                          number of concurrent threads supported by the hardware; negative
                          values are added to the number of concurrent threads supported
                          by the hardware).

          RETURNS: a tuple of tuples of integers

          NOTE:
             - the ordering of the indices corresponds to the atom ordering
                 in the query. For example, the first index is for the atom in
                 this molecule that matches the first atom in the query.
        """

class MolBundle:
    """A class for storing groups of related molecules."""

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, pklString: bytes) -> None:
        """Constructor from a binary string"""

    def ToBinary(self) -> bytes:
        """Returns a binary string representation of the MolBundle."""

    def __getitem__(self, idx: int) -> Mol: ...

    def __len__(self) -> int: ...

    def AddMol(self, nmol: Mol) -> int: ...

    def GetMol(self, idx: int) -> Mol: ...

    def Size(self) -> int: ...

    @overload
    def HasSubstructMatch(self, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
        Queries whether or not any molecule in the bundle contains a particular substructure.

             ARGUMENTS:
                  - query: a Molecule

                  - recursionPossible: (optional)

                  - useChirality: enables the use of stereochemistry in the matching

                  - useQueryQueryMatches: use query-query matching logic

             RETURNS: True or False
        """

    @overload
    def HasSubstructMatch(self, query: MolBundle, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
        Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.

             ARGUMENTS:
                  - query: a MolBundle

                  - recursionPossible: (optional)

                  - useChirality: enables the use of stereochemistry in the matching

                  - useQueryQueryMatches: use query-query matching logic

             RETURNS: True or False
        """

    @overload
    def HasSubstructMatch(self, query: Mol, params: SubstructMatchParameters | None) -> bool:
        """
        Queries whether or not any molecule in the bundle contains a particular substructure.

             ARGUMENTS:
                  - query: a Molecule

                  - params: parameters controlling the substructure match

             RETURNS: True or False
        """

    @overload
    def HasSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters | None) -> bool:
        """
        Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.

             ARGUMENTS:
                  - query: a MolBundle

                  - params: parameters controlling the substructure match

             RETURNS: True or False
        """

    @overload
    def GetSubstructMatch(self, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> list[int]:
        """
        Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query.

             ARGUMENTS:
                  - query: a Molecule

                  - useChirality: enables the use of stereochemistry in the matching

                  - useQueryQueryMatches: use query-query matching logic

             RETURNS: a tuple of integers

             NOTES:
                   - only a single match is returned
                   - the ordering of the indices corresponds to the atom ordering
                             in the query. For example, the first index is for the atom in
                             this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatch(self, query: MolBundle, useChirality: bool = False, useQueryQueryMatches: bool = False) -> list[int]:
        """
        Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.

             ARGUMENTS:
                  - query: a MolBundle

                  - useChirality: enables the use of stereochemistry in the matching

                  - useQueryQueryMatches: use query-query matching logic

             RETURNS: a tuple of integers

             NOTES:
                   - only a single match is returned
                   - the ordering of the indices corresponds to the atom ordering
                             in the query. For example, the first index is for the atom in
                             this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatch(self, query: Mol, params: SubstructMatchParameters | None) -> list[int]:
        """
        Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query.

             ARGUMENTS:
                  - query: a Molecule

                  - params: parameters controlling the substructure match

             RETURNS: a tuple of integers

             NOTES:
                   - only a single match is returned
                   - the ordering of the indices corresponds to the atom ordering
                             in the query. For example, the first index is for the atom in
                             this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters | None) -> list[int]:
        """
        Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.

             ARGUMENTS:
                  - query: a MolBundle

                  - params: parameters controlling the substructure match

             RETURNS: a tuple of integers

             NOTES:
                   - only a single match is returned
                   - the ordering of the indices corresponds to the atom ordering
                             in the query. For example, the first index is for the atom in
                             this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatches(self, query: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> list[list[int]]:
        """
        Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query.

             ARGUMENTS:
                  - query: a molecule.
                  - uniquify: (optional) determines whether or not the matches are uniquified.
                                                Defaults to 1.

                  - useChirality: enables the use of stereochemistry in the matching

                  - useQueryQueryMatches: use query-query matching logic

                  - maxMatches: The maximum number of matches that will be returned.
                                                     In high-symmetry cases with medium-sized molecules, it is
                                                     very easy to end up with a combinatorial explosion in the
                                                     number of possible matches. This argument prevents that from
                                                     having unintended consequences

             RETURNS: a tuple of tuples of integers

             NOTE:
                   - the ordering of the indices corresponds to the atom ordering
                             in the query. For example, the first index is for the atom in
                             this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatches(self, query: MolBundle, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> list[list[int]]:
        """
        Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.

             ARGUMENTS:
                  - query: a MolBundle.
                  - uniquify: (optional) determines whether or not the matches are uniquified.
                                                Defaults to 1.

                  - useChirality: enables the use of stereochemistry in the matching

                  - useQueryQueryMatches: use query-query matching logic

                  - maxMatches: The maximum number of matches that will be returned.
                                                     In high-symmetry cases with medium-sized molecules, it is
                                                     very easy to end up with a combinatorial explosion in the
                                                     number of possible matches. This argument prevents that from
                                                     having unintended consequences

             RETURNS: a tuple of tuples of integers

             NOTE:
                   - the ordering of the indices corresponds to the atom ordering
                             in the query. For example, the first index is for the atom in
                             this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatches(self, query: Mol, params: SubstructMatchParameters | None) -> list[list[int]]:
        """
        Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query.

             ARGUMENTS:
                  - query: a molecule.
                  - params: parameters controlling the substructure match

             RETURNS: a tuple of tuples of integers

             NOTE:
                   - the ordering of the indices corresponds to the atom ordering
                             in the query. For example, the first index is for the atom in
                             this molecule that matches the first atom in the query.
        """

    @overload
    def GetSubstructMatches(self, query: MolBundle, params: SubstructMatchParameters | None) -> list[list[int]]:
        """
        Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.

             ARGUMENTS:
                  - query: a MolBundle.
                  - params: parameters controlling the substructure match

             RETURNS: a tuple of tuples of integers

             NOTE:
                   - the ordering of the indices corresponds to the atom ordering
                             in the query. For example, the first index is for the atom in
                             this molecule that matches the first atom in the query.
        """

    def __getstate__(self) -> tuple[bytes]: ...

    def __setstate__(self, arg: tuple[bytes], /) -> None: ...

class FixedMolSizeMolBundle(MolBundle):
    """
    A class for storing groups of related molecules.
              Here related means that the molecules have to have the same number of atoms.
    """

    def __init__(self) -> None: ...

def MolBundleCanSerialize() -> bool:
    """
    Returns True if the MolBundle is serializable (requires boost serialization)
    """

class SubstanceGroupCState:
    """CSTATE for a SubstanceGroup"""

    def __init__(self) -> None: ...

    @property
    def bondIdx(self) -> int: ...

    @property
    def vector(self) -> rdkit.Geometry.rdGeometry.Point3D: ...

class SubstanceGroupAttach:
    """AttachPoint for a SubstanceGroup"""

    def __init__(self) -> None: ...

    @property
    def aIdx(self) -> int:
        """attachment index"""

    @property
    def lvIdx(self) -> int:
        """leaving atom or index (0 for implied)"""

    @property
    def id(self) -> str:
        """attachment id"""

class SubstanceGroup:
    """A collection of atoms and bonds with associated properties"""

    def GetOwningMol(self) -> Mol:
        """returns the molecule owning this SubstanceGroup"""

    def GetIndexInMol(self) -> int:
        """
        returns the index of this SubstanceGroup in the owning molecule's list.
        """

    def GetAtoms(self) -> list[int]:
        """returns a list of the indices of the atoms in this SubstanceGroup"""

    def GetParentAtoms(self) -> list[int]:
        """
        returns a list of the indices of the parent atoms in this SubstanceGroup
        """

    def GetBonds(self) -> list[int]:
        """returns a list of the indices of the bonds in this SubstanceGroup"""

    def SetAtoms(self, iterable: Iterable[int] | None) -> None:
        """
        Set the list of the indices of the atoms in this SubstanceGroup.
        Note that this does not update properties, CStates or Attachment Points.
        """

    def SetParentAtoms(self, iterable: Iterable[int] | None) -> None:
        """
        Set the list of the indices of the parent atoms in this SubstanceGroup.
        Note that this does not update properties, CStates or Attachment Points.
        """

    def SetBonds(self, iterable: Iterable[int] | None) -> None:
        """
        Set the list of the indices of the bonds in this SubstanceGroup.
        Note that this does not update properties, CStates or Attachment Points.
        """

    def AddAtomWithIdx(self, idx: int) -> None: ...

    def AddBondWithIdx(self, idx: int) -> None: ...

    def AddParentAtomWithIdx(self, idx: int) -> None: ...

    def AddAtomWithBookmark(self, mark: int) -> None: ...

    def AddParentAtomWithBookmark(self, mark: int) -> None: ...

    def AddCState(self, bondIdx: int, vector: rdkit.Geometry.rdGeometry.Point3D) -> None: ...

    def GetCStates(self) -> tuple[SubstanceGroupCState, ...]: ...

    def AddBondWithBookmark(self, mark: int) -> None: ...

    def AddAttachPoint(self, aIdx: int, lvIdx: int, idStr: str) -> None: ...

    def GetAttachPoints(self) -> tuple[SubstanceGroupAttach, ...]: ...

    def AddBracket(self, pts: Iterable[rdkit.Geometry.rdGeometry.Point3D]) -> None: ...

    def GetBrackets(self) -> tuple[tuple[rdkit.Geometry.rdGeometry.Point3D, rdkit.Geometry.rdGeometry.Point3D, rdkit.Geometry.rdGeometry.Point3D], ...]: ...

    def ClearBrackets(self) -> None:
        """Clear bracket definitions."""

    def ClearCStates(self) -> None:
        """Clear CSTATE entries."""

    def ClearAttachPoints(self) -> None:
        """Clear attachment points."""

    def SetProp(self, key: str, val: str, computed: bool = False) -> None:
        """sets the value of a particular property"""

    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None:
        """sets the value of a particular property"""

    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None:
        """sets the value of a particular property"""

    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None:
        """sets the value of a particular property"""

    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None:
        """sets the value of a particular property"""

    def HasProp(self, key: str) -> bool:
        """returns whether or not a particular property exists"""

    @overload
    def GetProp(self, key: str, autoConvert: bool = False) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - autoConvert: if True attempt to convert the property into a python object

          RETURNS: a string

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    @overload
    def GetProp(self, key: str, autoConvert: bool = False, default: object | None = None) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a string).

            - autoConvert: if True attempt to convert the property into a python object

            - default: value to return if the property is not present.

          RETURNS: the property value, or default if the property is not present.
        """

    @overload
    def GetIntProp(self, key: str) -> int:
        """returns the value of a particular property"""

    @overload
    def GetIntProp(self, key: str, default: object) -> object:
        """returns the value of a particular property, or default if not present"""

    @overload
    def GetUnsignedProp(self, key: str) -> int:
        """returns the value of a particular property"""

    @overload
    def GetUnsignedProp(self, key: str, default: object) -> object:
        """returns the value of a particular property, or default if not present"""

    @overload
    def GetDoubleProp(self, key: str) -> float:
        """returns the value of a particular property"""

    @overload
    def GetDoubleProp(self, key: str, default: object) -> object:
        """returns the value of a particular property, or default if not present"""

    @overload
    def GetBoolProp(self, key: str) -> bool:
        """returns the value of a particular property"""

    @overload
    def GetBoolProp(self, key: str, default: object) -> object:
        """returns the value of a particular property, or default if not present"""

    def GetUnsignedVectProp(self, key: str) -> list[int]:
        """returns the value of a particular property"""

    def GetStringVectProp(self, key: str) -> list[str]:
        """returns the value of a particular property"""

    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> list[str]:
        """Returns a list of the properties set on the SubstanceGroup."""

    def GetPropsAsDict(self, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict:
        """
        Returns a dictionary of the properties set on the SubstanceGroup.
         n.b. some properties cannot be converted to python types.
        """

    def ClearProp(self, key: str) -> None:
        """Removes a particular property (does nothing if not set)."""

def GetMolSubstanceGroups(mol: Mol) -> list[SubstanceGroup]:
    """returns a copy of the molecule's SubstanceGroups (if any)"""

def GetMolSubstanceGroupWithIdx(mol: Mol, idx: int) -> SubstanceGroup:
    """returns a particular SubstanceGroup from the molecule"""

def ClearMolSubstanceGroups(mol: Mol) -> None:
    """removes all SubstanceGroups from a molecule (if any)"""

def CreateMolSubstanceGroup(mol: Mol, type: str) -> SubstanceGroup:
    """
    creates a new SubstanceGroup associated with a molecule, returns the new SubstanceGroup
    """

def CreateMolDataSubstanceGroup(mol: Mol, fieldName: str, value: str) -> SubstanceGroup:
    """
    creates a new DATA SubstanceGroup associated with a molecule, returns the new SubstanceGroup
    """

def AddMolSubstanceGroup(mol: Mol, sgroup: SubstanceGroup) -> SubstanceGroup:
    """
    adds a copy of a SubstanceGroup to a molecule, returns the new SubstanceGroup
    """

class StereoType(enum.Enum):
    Unspecified = 0

    Atom_Tetrahedral = 1

    Atom_SquarePlanar = 2

    Atom_TrigonalBipyramidal = 3

    Atom_Octahedral = 4

    Bond_Double = 5

    Bond_Cumulene_Even = 6

    Bond_Atropisomer = 7

Unspecified: StereoSpecified = StereoSpecified.Unspecified

Atom_Tetrahedral: StereoType = StereoType.Atom_Tetrahedral

Atom_SquarePlanar: StereoType = StereoType.Atom_SquarePlanar

Atom_TrigonalBipyramidal: StereoType = StereoType.Atom_TrigonalBipyramidal

Atom_Octahedral: StereoType = StereoType.Atom_Octahedral

Bond_Double: StereoType = StereoType.Bond_Double

Bond_Cumulene_Even: StereoType = StereoType.Bond_Cumulene_Even

Bond_Atropisomer: StereoType = StereoType.Bond_Atropisomer

class StereoSpecified(enum.IntEnum):
    Unspecified = 0

    Specified = 1

    Unknown = 2

Specified: StereoSpecified = StereoSpecified.Specified

Unknown: StereoSpecified = StereoSpecified.Unknown

class StereoDescriptor(enum.Enum):
    NoValue = 0

    Tet_CW = 1

    Tet_CCW = 2

    Bond_Cis = 3

    Bond_Trans = 4

NoValue: StereoDescriptor = StereoDescriptor.NoValue

Tet_CW: StereoDescriptor = StereoDescriptor.Tet_CW

Tet_CCW: StereoDescriptor = StereoDescriptor.Tet_CCW

Bond_Cis: StereoDescriptor = StereoDescriptor.Bond_Cis

Bond_Trans: StereoDescriptor = StereoDescriptor.Bond_Trans

class StereoInfo:
    """Class describing stereochemistry"""

    NOATOM: Final[int] = ...
    """marker for unspecified int values"""

    @property
    def type(self) -> StereoType:
        """the type of stereo"""

    @type.setter
    def type(self, arg: StereoType, /) -> None: ...

    @property
    def specified(self) -> StereoSpecified:
        """whether or not it is specified"""

    @specified.setter
    def specified(self, arg: StereoSpecified, /) -> None: ...

    @property
    def centeredOn(self) -> int:
        """index of the item the stereo concerns"""

    @centeredOn.setter
    def centeredOn(self, arg: int, /) -> None: ...

    @property
    def descriptor(self) -> StereoDescriptor:
        """stereo descriptor"""

    @descriptor.setter
    def descriptor(self, arg: StereoDescriptor, /) -> None: ...

    @property
    def permutation(self) -> int:
        """permutation index (used for non-tetrahedral chirality)"""

    @permutation.setter
    def permutation(self, arg: int, /) -> None: ...

    @property
    def controllingAtoms(self) -> list[int]:
        """indices of the atoms controlling the stereo"""
