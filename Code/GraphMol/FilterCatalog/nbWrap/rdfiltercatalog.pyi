"""
Module containing FilterCatalog functionality for filtering molecules based on structural patterns.
"""

from collections.abc import Iterable, Sequence
import enum
from typing import TypeAlias, overload

import rdkit.Chem.rdchem
from rdkit.Chem.rdfiltercatalog import (
    FilterMatchOps as FilterMatchOps
)


def IntPair(a: int, b: int) -> tuple[int, int]:
    """Create an integer pair (tuple) for use in MatchTypeVect"""

class FilterMatch:
    """
    Object that holds the result of running FilterMatcherBase::GetMatches

     - filterMatch holds the FilterMatchBase that triggered the match
     - atomPairs holds the [ (query_atom_idx, target_atom_idx) ] pairs for the matches.

    Note that some matches may not have atom pairs (especially matches that
    use FilterMatchOps.Not
    """

    def __init__(self, filter: FilterMatcher, atomPairs: Sequence[tuple[int, int]]) -> None: ...

    @property
    def filterMatch(self) -> FilterMatcher: ...

    @property
    def atomPairs(self) -> list[tuple[int, int]]: ...

class FilterMatcher:
    """
    Base class for matching molecules to filters.

     A FilterMatcherBase supplies the following API
     - IsValid() returns True if the matcher is valid for use, False otherwise
     - HasMatch(mol) returns True if the molecule matches the filter
     - GetMatches(mol) -> [FilterMatch, FilterMatch] returns all the FilterMatch data
           that matches the molecule

    print( FilterMatcherBase ) will print user-friendly information about the filter
    Note that a FilterMatcherBase can be combined from many FilterMatcherBases
    This is why GetMatches can return multiple FilterMatcherBases.
    >>> from rdkit.Chem.FilterCatalog import *
    >>> carbon_matcher = SmartsMatcher('Carbon', '[#6]', 0, 1)
    >>> oxygen_matcher = SmartsMatcher('Oxygen', '[#8]', 0, 1)
    >>> co_matcher = FilterMatchOps.Or(carbon_matcher, oxygen_matcher)
    >>> mol = Chem.MolFromSmiles('C')
    >>> matches = co_matcher.GetMatches(mol)
    >>> len(matches)
    1
    >>> print(matches[0].filterMatch)
    Carbon
    """

    def __init__(self, name: str) -> None: ...

    def IsValid(self) -> bool:
        """Return True if the filter matcher is valid, False otherwise"""

    def HasMatch(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """Returns True if mol matches the filter"""

    def GetMatches(self, mol: rdkit.Chem.rdchem.Mol) -> list[FilterMatch]:
        """Returns the list of matching subfilters mol matches any filter"""

    def GetName(self) -> str: ...

    def __str__(self) -> str: ...

FilterMatcherBase: TypeAlias = FilterMatcher

PythonFilterMatcher: TypeAlias = FilterMatcher

class SmartsMatcher(FilterMatcher):
    """
    Smarts Matcher Filter
     basic constructors:
       SmartsMatcher( name, smarts_pattern, minCount=1, maxCount=UINT_MAX )
       SmartsMatcher( name, molecule, minCount=1, maxCount=UINT_MAX )

      note: If the supplied smarts pattern is not valid, the IsValid() function will
       return False
    >>> from rdkit.Chem.FilterCatalog import *
    >>> minCount, maxCount = 1,2
    >>> carbon_matcher = SmartsMatcher('Carbon', '[#6]', minCount, maxCount)
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CC')))
    True
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))
    False
    >>> carbon_matcher.SetMinCount(2)
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('C')))
    False
    >>> carbon_matcher.SetMaxCount(3)
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))
    True
    """

    @overload
    def __init__(self, name: str) -> None: ...

    @overload
    def __init__(self, rhs: rdkit.Chem.rdchem.Mol) -> None:
        """Construct from a molecule"""

    @overload
    def __init__(self, name: str, mol: rdkit.Chem.rdchem.Mol, minCount: int = 1, maxCount: int = 4294967295) -> None:
        """Construct from a name, molecule, minimum and maximum count"""

    @overload
    def __init__(self, name: str, smarts: str, minCount: int = 1, maxCount: int = 4294967295) -> None:
        """Construct from a name, smarts pattern, minimum and maximum count"""

    def IsValid(self) -> bool:
        """Returns True if the SmartsMatcher is valid"""

    @overload
    def SetPattern(self, pat: rdkit.Chem.rdchem.Mol) -> None:
        """Set the pattern molecule for the SmartsMatcher"""

    @overload
    def SetPattern(self, pat: str) -> None:
        """
        Set the smarts pattern for the Smarts Matcher (warning: MinimumCount is not reset)
        """

    def GetPattern(self) -> rdkit.Chem.rdchem.Mol | None: ...

    def GetMinCount(self) -> int:
        """Get the minimum times pattern must appear for the filter to match"""

    def SetMinCount(self, count: int) -> None:
        """Set the minimum times pattern must appear to match"""

    def GetMaxCount(self) -> int:
        """Get the maximum times pattern can appear for the filter to match"""

    def SetMaxCount(self, count: int) -> None:
        """Set the maximum times pattern can appear for the filter to match"""

class ExclusionList(FilterMatcher):
    def __init__(self) -> None: ...

    def SetExclusionPatterns(self, list: Iterable[FilterMatcher]) -> None:
        """Set a list of FilterMatcherBases that should not appear in a molecule"""

    def AddPattern(self, base: FilterMatcher) -> None:
        """Add a FilterMatcherBase that should not appear in a molecule"""

class FilterHierarchyMatcher(FilterMatcher):
    """
    Hierarchical Filter
     basic constructors:
       FilterHierarchyMatcher( matcher )
       where can be any FilterMatcherBase (SmartsMatcher, etc)
     FilterHierarchyMatcher's have children and can form matching
      trees.  When GetFilterMatches is called, the most specific (
      i.e. lowest node in a branch) is returned.

     n.b. A FilterHierarchicalMatcher of functional groups is returned
      by calling GetFunctionalGroupHierarchy()

    >>> from rdkit.Chem import MolFromSmiles
    >>> from rdkit.Chem.FilterCatalog import *
    >>> functionalGroups = GetFunctionalGroupHierarchy()
    >>> [match.filterMatch.GetName()
    ...     for match in functionalGroups.GetFilterMatches(
    ...         MolFromSmiles('c1ccccc1Cl'))]
    ['Halogen.Aromatic', 'Halogen.NotFluorine.Aromatic']
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, matcher: FilterMatcher) -> None:
        """Construct from a filtermatcher"""

    def SetPattern(self, matcher: FilterMatcher) -> None:
        """
        Set the filtermatcher pattern for this node. An empty node is
        considered a root node and passes along the matches to the children.
        """

    def AddChild(self, hierarchy: FilterHierarchyMatcher) -> FilterHierarchyMatcher:
        """Add a child node to this hierarchy."""

class FilterCatalogEntry:
    """
    FilterCatalogEntry
    A filter catalog entry is an entry in a filter catalog.
    Each filter is named and is used to flag a molecule usually for some
    undesirable property.

    For example, a PAINS (Pan Assay INterference) catalog entry be appear as
    follows:

    >>> from rdkit.Chem.FilterCatalog import *
    >>> params = FilterCatalogParams()
    >>> params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    True
    >>> catalog = FilterCatalog(params)
    >>> mol = Chem.MolFromSmiles('O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2')
    >>> entry = catalog.GetFirstMatch(mol)
    >>> print (entry.GetProp('Scope'))
    PAINS filters (family A)
    >>> print (entry.GetDescription())
    hzone_phenol_A(479)
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, name: str, matcher: FilterMatcher) -> None: ...

    def IsValid(self) -> bool: ...

    def GetDescription(self) -> str:
        """Get the description of the catalog entry"""

    def SetDescription(self, description: str) -> None:
        """Set the description of the catalog entry"""

    def GetFilterMatches(self, mol: rdkit.Chem.rdchem.Mol) -> list[FilterMatch]:
        """Retrieve the list of filters that match the molecule"""

    def HasFilterMatch(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        Returns True if the catalog entry contains filters that match the molecule
        """

    def Serialize(self) -> bytes: ...

    def GetPropList(self) -> list[str]: ...

    def SetProp(self, key: str, val: str) -> None: ...

    def GetProp(self, key: str) -> str: ...

    def ClearProp(self, key: str) -> None: ...

def GetFunctionalGroupHierarchy() -> FilterCatalog:
    """Returns the functional group hierarchy filter catalog"""

def GetFlattenedFunctionalGroupHierarchy(normalized: bool = False) -> dict[str, rdkit.Chem.rdchem.Mol]:
    """
    Returns the flattened functional group hierarchy as a dictionary  of name:ROMOL_SPTR substructure items
    """

class FilterCatalogParams:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, catalogs: FilterCatalogs) -> None:
        """
        Construct from a FilterCatalogs identifier (i.e. FilterCatalogParams.PAINS)
        """

    def AddCatalog(self, catalogs: FilterCatalogs) -> bool: ...

    class FilterCatalogs(enum.IntEnum):
        PAINS_A = 2

        PAINS_B = 4

        PAINS_C = 8

        PAINS = 14

        BRENK = 16

        NIH = 32

        ZINC = 64

        CHEMBL_Glaxo = 128

        CHEMBL_Dundee = 256

        CHEMBL_BMS = 512

        CHEMBL_SureChEMBL = 1024

        CHEMBL_MLSMR = 2048

        CHEMBL_Inpharmatica = 4096

        CHEMBL_LINT = 8192

        CHEMBL = 16256

        ALL = 16382

class FilterCatalog:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, pickle: bytes) -> None: ...

    @overload
    def __init__(self, pickle: str) -> None: ...

    @overload
    def __init__(self, params: FilterCatalogParams) -> None: ...

    @overload
    def __init__(self, catalogs: FilterCatalogParams.FilterCatalogs) -> None: ...

    def Serialize(self) -> bytes: ...

    def AddEntry(self, entry: FilterCatalogEntry) -> None:
        """Add a FilterCatalogEntry to the catalog"""

    def RemoveEntry(self, obj: object) -> bool:
        """Remove the given entry from the catalog"""

    def GetNumEntries(self) -> int:
        """Returns the number of entries in the catalog"""

    def GetEntryWithIdx(self, idx: int) -> FilterCatalogEntry:
        """Return the FilterCatalogEntry at the specified index"""

    def GetEntry(self, idx: int) -> FilterCatalogEntry:
        """Return the FilterCatalogEntry at the specified index"""

    def HasMatch(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """Returns True if the catalog has an entry that matches mol"""

    def GetFirstMatch(self, mol: rdkit.Chem.rdchem.Mol) -> FilterCatalogEntry:
        """Return the first catalog entry that matches mol"""

    def GetMatches(self, mol: rdkit.Chem.rdchem.Mol) -> list[FilterCatalogEntry]:
        """Return all catalog entries that match mol"""

    def GetFilterMatches(self, mol: rdkit.Chem.rdchem.Mol) -> list[FilterMatch]:
        """Return every matching filter from all catalog entries that match mol"""

    def __setstate__(self, arg: object, /) -> None: ...

    def __getstate__(self) -> tuple[bytes, dict]: ...

def FilterCatalogCanSerialize() -> bool:
    """
    Returns True if the FilterCatalog is serializable (requires boost serialization)
    """

def RunFilterCatalog(filterCatalog: FilterCatalog, smiles: Sequence[str], numThreads: int = 1) -> list[list[FilterCatalogEntry]]:
    """
    Run the filter catalog on the input list of smiles strings.
    Use numThreads=0 to use all available processors.
    Returns a vector of vectors. For each input smiles, a vector of
    FilterCatalogEntry objects are returned for each matched filter. If a
    molecule matches no filter, the vector will be empty. If a smiles string
    can't be parsed, a 'Bad smiles' entry is returned.
    """
