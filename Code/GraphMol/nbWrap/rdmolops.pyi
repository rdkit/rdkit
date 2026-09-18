"""Module containing RDKit functionality for manipulating molecules."""

from collections.abc import Iterable, Iterator, Sequence
import enum
from typing import Annotated, Final, overload

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdchem
import rdkit.DataStructs.cDataStructs
import rdkit.Geometry.rdGeometry


class SanitizeFlags(enum.IntEnum):
    SANITIZE_NONE = 0

    SANITIZE_CLEANUP = 1

    SANITIZE_PROPERTIES = 2

    SANITIZE_SYMMRINGS = 4

    SANITIZE_KEKULIZE = 8

    SANITIZE_FINDRADICALS = 16

    SANITIZE_SETAROMATICITY = 32

    SANITIZE_SETCONJUGATION = 64

    SANITIZE_SETHYBRIDIZATION = 128

    SANITIZE_CLEANUPCHIRALITY = 256

    SANITIZE_CLEANUPATROPISOMERS = 2048

    SANITIZE_ADJUSTHS = 512

    SANITIZE_CLEANUP_ORGANOMETALLICS = 1024

    SANITIZE_ALL = 268435455

SANITIZE_NONE: SanitizeFlags = SanitizeFlags.SANITIZE_NONE

SANITIZE_CLEANUP: SanitizeFlags = SanitizeFlags.SANITIZE_CLEANUP

SANITIZE_PROPERTIES: SanitizeFlags = SanitizeFlags.SANITIZE_PROPERTIES

SANITIZE_SYMMRINGS: SanitizeFlags = SanitizeFlags.SANITIZE_SYMMRINGS

SANITIZE_KEKULIZE: SanitizeFlags = SanitizeFlags.SANITIZE_KEKULIZE

SANITIZE_FINDRADICALS: SanitizeFlags = SanitizeFlags.SANITIZE_FINDRADICALS

SANITIZE_SETAROMATICITY: SanitizeFlags = SanitizeFlags.SANITIZE_SETAROMATICITY

SANITIZE_SETCONJUGATION: SanitizeFlags = SanitizeFlags.SANITIZE_SETCONJUGATION

SANITIZE_SETHYBRIDIZATION: SanitizeFlags = SanitizeFlags.SANITIZE_SETHYBRIDIZATION

SANITIZE_CLEANUPCHIRALITY: SanitizeFlags = SanitizeFlags.SANITIZE_CLEANUPCHIRALITY

SANITIZE_CLEANUPATROPISOMERS: SanitizeFlags = SanitizeFlags.SANITIZE_CLEANUPATROPISOMERS

SANITIZE_ADJUSTHS: SanitizeFlags = SanitizeFlags.SANITIZE_ADJUSTHS

SANITIZE_CLEANUP_ORGANOMETALLICS: SanitizeFlags = SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS

SANITIZE_ALL: SanitizeFlags = SanitizeFlags.SANITIZE_ALL

def DetectBondStereoChemistry(mol: rdkit.Chem.rdchem.Mol, conformer: rdkit.Chem.rdchem.Conformer) -> None:
    """
    Assign stereochemistry to bonds based on coordinates and a conformer.
            DEPRECATED

      ARGUMENTS:

        - mol: the molecule to be modified
        - conformer: Conformer providing the coordinates
    """

def DetectBondStereochemistry(mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> None:
    """
    DEPRECATED
        - mol: the molecule to be modified
        - confId: Conformer to use for the coordinates
    """

def SetDoubleBondNeighborDirections(mol: rdkit.Chem.rdchem.Mol, conf: rdkit.Chem.rdchem.Conformer | None = None) -> None:
    """
    Uses the stereo info on double bonds to set the directions of neighboring single bonds

      ARGUMENTS:

        - mol: the molecule to be modified
    """

def SetBondStereoFromDirections(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    Uses the directions of neighboring bonds to set cis/trans stereo on double bonds.

      ARGUMENTS:

        - mol: the molecule to be modified
    """

def SanitizeMol(mol: rdkit.Chem.rdchem.Mol, sanitizeOps: int = 268435455, catchErrors: bool = False) -> SanitizeFlags:
    """
    Kekulize, check valencies, set aromaticity, conjugation and hybridization

        - The molecule is modified in place.

        - If sanitization fails, an exception will be thrown unless catchErrors is set

      ARGUMENTS:

        - mol: the molecule to be modified
        - sanitizeOps: (optional) sanitization operations to be carried out
          these should be constructed by or'ing together the
          operations in rdkit.Chem.SanitizeFlags
        - catchErrors: (optional) if provided, instead of raising an exception
          when sanitization fails (the default behavior), the 
          first operation that failed (as defined in rdkit.Chem.SanitizeFlags)
          is returned. Zero is returned on success.
    """

def GetSSSR(mol: rdkit.Chem.rdchem.Mol, includeDativeBonds: bool = False, includeHydrogenBonds: bool = False) -> list[list[int]]:
    """
    Get the smallest set of simple rings for a molecule.

      ARGUMENTS:

        - mol: the molecule to use.
        - includeDativeBonds: whether or not dative bonds should be included in the ring finding.
        - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring finding.

      RETURNS: a sequence of sequences containing the rings found as atom ids
             The length of this will be equal to NumBonds-NumAtoms+1 for single-fragment molecules.
    """

class SymmetrizeSSSRAlgorithm(enum.Enum):
    DEFAULT = 0

    LEGACY = 1

    RDL = 2

def SetUseLegacyRingFinding(val: bool) -> None:
    """sets usage of the legacy symmetric SSSR code during sanitization"""

def GetUseLegacyRingFinding() -> bool:
    """
    returns whether or not the legacy symmetric SSSR code is being used during sanitization
    """

def GetSymmSSSR(mol: rdkit.Chem.rdchem.Mol, includeDativeBonds: bool = False, includeHydrogenBonds: bool = False, algorithm: SymmetrizeSSSRAlgorithm = SymmetrizeSSSRAlgorithm.DEFAULT, recalcSSSR: bool = True) -> list[list[int]]:
    """
    Get a symmetrized SSSR for a molecule.

      The symmetrized SSSR is at least as large as the SSSR for a molecule.
      In certain highly-symmetric cases (e.g. cubane), the symmetrized SSSR can be
      a bit larger (i.e. the number of symmetrized rings is >= NumBonds-NumAtoms+1).

      ARGUMENTS:

        - mol: the molecule to use.
        - includeDativeBonds: whether or not dative bonds should be included in the ring finding.
        - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring finding.

      RETURNS: a sequence of sequences containing the rings found as atom ids
    """

def SetTerminalAtomCoords(mol: rdkit.Chem.rdchem.Mol, idx: int, otherIdx: int) -> None:
    """
    Sets Cartesian coordinates for a terminal atom.

      Useful for growing an atom off a molecule with sensible 
      coordinates based on the geometry of the neighbor.

      NOTE: this sets the appropriate coordinates in all of the molecule's conformers 
      ARGUMENTS:

        - mol: the molecule the atoms belong to.
        - idx: index of the terminal atom whose coordinates are set.
        - mol: index of the bonded neighbor atom.

      RETURNS: Nothing
    """

def FastFindRings(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    Does a non-SSSR ring finding for a molecule.

      ARGUMENTS:

        - mol: the molecule to use.

      RETURNS: Nothing
    """

def FindRingFamilies(mol: rdkit.Chem.rdchem.Mol, includeDativeBonds: bool = False, includeHydrogenBonds: bool = False) -> None:
    """
    Generate Unique Ring Families.

      ARGUMENTS:

        - mol: the molecule to use.
        - includeDativeBonds: whether or not dative bonds should be included in the ring families finding.
        - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring families .

      RETURNS: Nothing
    """

class AddHsParameters:
    """Parameters controlling H addition."""

    def __init__(self) -> None: ...

    @property
    def explicitOnly(self) -> bool:
        """only add explict Hs"""

    @explicitOnly.setter
    def explicitOnly(self, arg: bool, /) -> None: ...

    @property
    def addCoords(self) -> bool:
        """add coordinates for the Hs"""

    @addCoords.setter
    def addCoords(self, arg: bool, /) -> None: ...

    @property
    def addResidueInfo(self) -> bool:
        """add residue info to the Hs"""

    @addResidueInfo.setter
    def addResidueInfo(self, arg: bool, /) -> None: ...

    @property
    def skipQueries(self) -> bool:
        """do not add Hs to query atoms or atoms with query bonds"""

    @skipQueries.setter
    def skipQueries(self, arg: bool, /) -> None: ...

@overload
def AddHs(mol: rdkit.Chem.rdchem.Mol, params: AddHsParameters, onlyOnAtoms: Iterable[int] | None = None) -> rdkit.Chem.rdchem.Mol:
    """
    Adds hydrogens to the graph of a molecule.

      ARGUMENTS:

        - mol: the molecule to be modified

        - params: AddHsParameters object controlling the addition.

        - onlyOnAtoms: (optional) if this sequence is provided, only these atoms will be
          considered to have Hs added to them

      RETURNS: a new molecule with added Hs

      NOTES:

        - The original molecule is *not* modified.

        - Much of the code assumes that Hs are not included in the molecular
          topology, so be *very* careful with the molecule that comes back from
          this function.\\n
    """

@overload
def AddHs(mol: rdkit.Chem.rdchem.Mol, explicitOnly: bool = False, addCoords: bool = False, onlyOnAtoms: Iterable[int] | None = None, addResidueInfo: bool = False) -> rdkit.Chem.rdchem.Mol:
    """
    Adds hydrogens to the graph of a molecule.

      ARGUMENTS:

        - mol: the molecule to be modified

        - explicitOnly: (optional) if this toggle is set, only explicit Hs will
          be added to the molecule.  Default value is 0 (add implicit and explicit Hs).

        - addCoords: (optional) if this toggle is set, The Hs will have 3D coordinates
          set.  Default value is 0 (no 3D coords).

        - onlyOnAtoms: (optional) if this sequence is provided, only these atoms will be
          considered to have Hs added to them

        - addResidueInfo: (optional) if this is true, add residue info to
          hydrogen atoms (useful for PDB files).

      RETURNS: a new molecule with added Hs

      NOTES:

        - The original molecule is *not* modified.

        - Much of the code assumes that Hs are not included in the molecular
          topology, so be *very* careful with the molecule that comes back from
          this function.
    """

@overload
def RemoveHs(mol: rdkit.Chem.rdchem.Mol, implicitOnly: bool = False, updateExplicitCount: bool = False, sanitize: bool = True) -> rdkit.Chem.rdchem.Mol:
    """
    Removes any hydrogens from the graph of a molecule.

      ARGUMENTS:

        - mol: the molecule to be modified

        - implicitOnly: (optional) if this toggle is set, only implicit Hs will
          be removed from the graph.  Default value is 0 (remove implicit and explicit Hs).

        - updateExplicitCount: (optional) if this toggle is set, the explicit H count on atoms with 
          Hs will be updated. Default value is 0 (do not update explicit H count).

        - sanitize: (optional) if this toggle is set, the molecule will be sanitized after the Hs
          are removed. Default value is 1 (do sanitize).

      RETURNS: a new molecule with the Hs removed

      NOTES:

        - The original molecule is *not* modified.
        - Hydrogens which aren't connected to a heavy atom will not be
          removed.  This prevents molecules like [H][H] from having
          all atoms removed.
        - Labelled hydrogen (e.g. atoms with atomic number=1, but isotope > 1),
          will not be removed.
        - two coordinate Hs, like the central H in C[H-]C, will not be removed
        - Hs connected to dummy atoms will not be removed
        - Hs that are part of the definition of double bond Stereochemistry
          will not be removed
        - Hs that are not connected to anything else will not be removed
    """

@overload
def RemoveHs(mol: rdkit.Chem.rdchem.Mol, params: RemoveHsParameters, sanitize: bool = True) -> rdkit.Chem.rdchem.Mol:
    """
    Returns a copy of the molecule with Hs removed. Which Hs are removed is controlled by the params argument
    """

class RemoveHsParameters:
    """Parameters controlling which Hs are removed."""

    def __init__(self) -> None: ...

    @property
    def removeDegreeZero(self) -> bool:
        """hydrogens that have no bonds"""

    @removeDegreeZero.setter
    def removeDegreeZero(self, arg: bool, /) -> None: ...

    @property
    def removeHigherDegrees(self) -> bool:
        """hydrogens with two (or more) bonds"""

    @removeHigherDegrees.setter
    def removeHigherDegrees(self, arg: bool, /) -> None: ...

    @property
    def removeOnlyHNeighbors(self) -> bool:
        """hydrogens with bonds only to other hydrogens"""

    @removeOnlyHNeighbors.setter
    def removeOnlyHNeighbors(self, arg: bool, /) -> None: ...

    @property
    def removeIsotopes(self) -> bool:
        """hydrogens with non-default isotopes"""

    @removeIsotopes.setter
    def removeIsotopes(self, arg: bool, /) -> None: ...

    @property
    def removeAndTrackIsotopes(self) -> bool:
        """
        hydrogens with non-default isotopes and store them in the _isotopicHs atom property such that AddHs() can add the same isotope at a later stage
        """

    @removeAndTrackIsotopes.setter
    def removeAndTrackIsotopes(self, arg: bool, /) -> None: ...

    @property
    def removeDummyNeighbors(self) -> bool:
        """hydrogens with at least one dummy-atom neighbor"""

    @removeDummyNeighbors.setter
    def removeDummyNeighbors(self, arg: bool, /) -> None: ...

    @property
    def removeDefiningBondStereo(self) -> bool:
        """hydrogens defining bond stereochemistry"""

    @removeDefiningBondStereo.setter
    def removeDefiningBondStereo(self, arg: bool, /) -> None: ...

    @property
    def removeWithWedgedBond(self) -> bool:
        """hydrogens with wedged bonds to them"""

    @removeWithWedgedBond.setter
    def removeWithWedgedBond(self, arg: bool, /) -> None: ...

    @property
    def removeWithQuery(self) -> bool:
        """hydrogens with queries defined"""

    @removeWithQuery.setter
    def removeWithQuery(self, arg: bool, /) -> None: ...

    @property
    def removeMapped(self) -> bool:
        """mapped hydrogens"""

    @removeMapped.setter
    def removeMapped(self, arg: bool, /) -> None: ...

    @property
    def removeInSGroups(self) -> bool:
        """hydrogens involved in SubstanceGroups"""

    @removeInSGroups.setter
    def removeInSGroups(self, arg: bool, /) -> None: ...

    @property
    def removeNonimplicit(self) -> bool:
        """DEPRECATED"""

    @removeNonimplicit.setter
    def removeNonimplicit(self, arg: bool, /) -> None: ...

    @property
    def removeHydrides(self) -> bool:
        """hydrogens with formal charge -1"""

    @removeHydrides.setter
    def removeHydrides(self, arg: bool, /) -> None: ...

    @property
    def removeNontetrahedralNeighbors(self) -> bool:
        """hydrogens with neighbors that have non-tetrahedral stereochemistry"""

    @removeNontetrahedralNeighbors.setter
    def removeNontetrahedralNeighbors(self, arg: bool, /) -> None: ...

    @property
    def showWarnings(self) -> bool:
        """display warning messages for some classes of removed Hs"""

    @showWarnings.setter
    def showWarnings(self, arg: bool, /) -> None: ...

    @property
    def updateExplicitCount(self) -> bool:
        """DEPRECATED"""

    @updateExplicitCount.setter
    def updateExplicitCount(self, arg: bool, /) -> None: ...

def RemoveAllHs(mol: rdkit.Chem.rdchem.Mol, sanitize: bool = True) -> rdkit.Chem.rdchem.Mol:
    """Returns a copy of the molecule with all Hs removed."""

def MergeQueryHs(mol: rdkit.Chem.rdchem.Mol, mergeUnmappedOnly: bool = False, mergeIsotopes: bool = False) -> rdkit.Chem.rdchem.Mol:
    """merges hydrogens into their neighboring atoms as queries"""

def HasQueryHs(mol: rdkit.Chem.rdchem.Mol) -> tuple[bool, bool]:
    """
    Check to see if the molecule has query Hs, this is normally used on query molecules
    such as those returned from MolFromSmarts
    Example: 
          (hasQueryHs, hasUnmergeableQueryHs) = HasQueryHs(mol)

    if hasUnmergeableQueryHs, these query hs cannot be removed by calling
    MergeQueryHs
    """

def DeleteSubstructs(mol: rdkit.Chem.rdchem.Mol, query: rdkit.Chem.rdchem.Mol, onlyFrags: bool = False, useChirality: bool = False) -> rdkit.Chem.rdchem.Mol:
    """
    Removes atoms matching a substructure query from a molecule

      ARGUMENTS:

        - mol: the molecule to be modified

        - query: the molecule to be used as a substructure query

        - onlyFrags: (optional) if this toggle is set, atoms will only be removed if
          the entire fragment in which they are found is matched by the query.
          See below for examples.
          Default value is 0 (remove the atoms whether or not the entire fragment matches)

        - useChirality: (optional) match the substructure query using chirality

      RETURNS: a new molecule with the substructure removed

      NOTES:

        - The original molecule is *not* modified.

      EXAMPLES:

       The following examples substitute SMILES/SMARTS strings for molecules, you'd have
       to actually use molecules:

        - DeleteSubstructs('CCOC','OC') -> 'CC'

        - DeleteSubstructs('CCOC','OC',1) -> 'CCOC'

        - DeleteSubstructs('CCOCCl.Cl','Cl',1) -> 'CCOCCl'

        - DeleteSubstructs('CCOCCl.Cl','Cl') -> 'CCOC'
    """

def MurckoDecompose(mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
    """Do a Murcko decomposition and return the scaffold"""

def CombineMols(mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol, offset: rdkit.Geometry.rdGeometry.Point3D | None = None) -> rdkit.Chem.rdchem.Mol:
    """Combine the atoms from two molecules to produce a third"""

def ReplaceSubstructs(mol: rdkit.Chem.rdchem.Mol, query: rdkit.Chem.rdchem.Mol, replacement: rdkit.Chem.rdchem.Mol, replaceAll: bool = False, replacementConnectionPoint: int = 0, useChirality: bool = False) -> tuple[rdkit.Chem.rdchem.Mol, ...]:
    """
    Replaces atoms matching a substructure query in a molecule

      ARGUMENTS:

        - mol: the molecule to be modified

        - query: the molecule to be used as a substructure query

        - replacement: the molecule to be used as the replacement

        - replaceAll: (optional) if this toggle is set, all substructures matching
          the query will be replaced in a single result, otherwise each result will
          contain a separate replacement.
          Default value is False (return multiple replacements)
        - replacementConnectionPoint: (optional) index of the atom in the replacement that
          the bond should be made to.
        - useChirality: (optional) match the substructure query using chirality

      RETURNS: a tuple of new molecules with the substructures replaced removed

      NOTES:

        - The original molecule is *not* modified.
        - A bond is only formed to the remaining atoms, if any, that were bonded 
          to the first atom in the substructure query. (For finer control over
          substructure replacement, consider using ChemicalReaction.)

      EXAMPLES:

       The following examples substitute SMILES/SMARTS strings for molecules, you'd have
       to actually use molecules:

        - ReplaceSubstructs('CCOC','O[CH3]','NC') -> ('CCNC',)

        - ReplaceSubstructs('COCCOC','O[CH3]','NC') -> ('COCCNC','CNCCOC')

        - ReplaceSubstructs('COCCOC','O[CH3]','NC',True) -> ('CNCCNC',)

        - ReplaceSubstructs('COCCOC','O[CH3]','CN',True,1) -> ('CNCCNC',)

        - ReplaceSubstructs('CCOC','[CH3]O','NC') -> ('CC.CN',)
    """

def GetMostSubstitutedCoreMatch(mol: rdkit.Chem.rdchem.Mol, core: rdkit.Chem.rdchem.Mol, matches: object) -> list[int]:
    """
    Postprocesses the results of a mol.GetSubstructMatches(core) call 
    where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups). 
    It returns the match with the largest number of non-hydrogen matches to 
    the terminal dummy atoms.

      ARGUMENTS:

        - mol: the molecule GetSubstructMatches was run on

        - core: the molecule used as a substructure query

        - matches: the result returned by GetSubstructMatches

      RETURNS: the tuple where terminal dummy atoms in the core match the largest 
               number of non-hydrogen atoms in mol
    """

def SortMatchesByDegreeOfCoreSubstitution(mol: rdkit.Chem.rdchem.Mol, core: rdkit.Chem.rdchem.Mol, matches: object) -> list[list[int]]:
    """
    Postprocesses the results of a mol.GetSubstructMatches(core) call 
    where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups). 
    It returns a copy of matches sorted by decreasing number of non-hydrogen matches 
    to the terminal dummy atoms.

      ARGUMENTS:

        - mol: the molecule GetSubstructMatches was run on

        - core: the molecule used as a substructure query

        - matches: the result returned by GetSubstructMatches

      RETURNS: a copy of matches sorted by decreasing number of non-hydrogen matches 
               to the terminal dummy atoms
    """

def MolAddRecursiveQueries(mol: rdkit.Chem.rdchem.Mol, queries: dict[str, rdkit.Chem.rdchem.Mol], propName: str) -> None:
    """Adds named recursive queries to atoms"""

def ParseMolQueryDefFile(fileobj: object, standardize: bool = True, delimiter: str = '\t', comment: str = '//', nameColumn: int = 0, smartsColumn: int = 1) -> dict[str, rdkit.Chem.rdchem.Mol]:
    """reads query definitions from a simply formatted file"""

def GetDistanceMatrix(mol: rdkit.Chem.rdchem.Mol, useBO: bool = False, useAtomWts: bool = False, force: bool = False, prefix: str = '') -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
    """
    Returns the molecule's topological distance matrix.

      ARGUMENTS:

        - mol: the molecule to use

        - useBO: (optional) toggles use of bond orders in calculating the distance matrix.
          Default value is 0.

        - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
          matrix (to return a "Balaban" distance matrix).
          Default value is 0.

        - force: (optional) forces the calculation to proceed, even if there is a cached value.
          Default value is 0.

        - prefix: (optional, internal use) sets the prefix used in the property cache
          Default value is .

      RETURNS: a Numeric array of floats with the distance matrix
    """

def Get3DDistanceMatrix(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomWts: bool = False, force: bool = False, prefix: str = '') -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
    """
    Returns the molecule's 3D distance matrix.

      ARGUMENTS:

        - mol: the molecule to use

        - confId: (optional) chooses the conformer Id to use
          Default value is -1.

        - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
          matrix (to return a "Balaban" distance matrix).
          Default value is 0.

        - force: (optional) forces the calculation to proceed, even if there is a cached value.
          Default value is 0.

        - prefix: (optional, internal use) sets the prefix used in the property cache
          Default value is .

      RETURNS: a Numeric array of floats with the distance matrix
    """

def GetAdjacencyMatrix(mol: rdkit.Chem.rdchem.Mol, useBO: bool = False, emptyVal: int = 0, force: bool = False, prefix: str = '') -> object:
    """
    Returns the molecule's adjacency matrix.

      ARGUMENTS:

        - mol: the molecule to use

        - useBO: (optional) toggles use of bond orders in calculating the matrix.
          Default value is 0.

        - emptyVal: (optional) sets the elements of the matrix between non-adjacent atoms
          Default value is 0.

        - force: (optional) forces the calculation to proceed, even if there is a cached value.
          Default value is 0.

        - prefix: (optional, internal use) sets the prefix used in the property cache
          Default value is .

      RETURNS: a Numeric array of floats containing the adjacency matrix
    """

def Kekulize(mol: rdkit.Chem.rdchem.Mol, clearAromaticFlags: bool = False, canonical: bool = True) -> None:
    """
    Kekulizes the molecule

      ARGUMENTS:

        - mol: the molecule to use

        - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the
          molecule will be marked non-aromatic following the kekulization.
          Default value is False.

        - canonical: (optional) if true, uses canonical atom ranking so
          that the kekulization result is independent of the atom ordering in the
          molecule.  Set to false to skip the ranking step for better performance
          when deterministic output is not required (e.g. during sanitization).
          Note, this "canonical" order only really makes sense when the molecule's
          chemistry is sane, like after sanitization. If stereochemistry hasn't been
          perceived, the chemistry of the molecule is inconsistent, and
          "canonical" atom ranks are only a technical artifact.

      NOTES:

        - The molecule is modified in place.

        - this does not modify query bonds which have bond type queries (like those
          which come from SMARTS) or rings containing them.

        - even if clearAromaticFlags is False the BondType for all modified
          aromatic bonds will be changed from AROMATIC to SINGLE or DOUBLE
          Kekulization.
    """

def KekulizeIfPossible(mol: rdkit.Chem.rdchem.Mol, clearAromaticFlags: bool = False, canonical: bool = True) -> None:
    """
    Kekulizes the molecule if possible. Otherwise the molecule is not modified

      ARGUMENTS:

        - mol: the molecule to use

        - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the 
          molecule will be marked non-aromatic if the kekulization succeds.
          Default value is False.

        - canonical: (optional) if true  uses canonical atom ranking so
          that the kekulization result is independent of the atom ordering in the
          molecule.  Set to false to skip the ranking step for better performance
          when deterministic output is not required (e.g. during sanitization).
          Note, this "canonical" order only really makes sense when the molecule's
          chemistry is sane, like after sanitization. If stereochemistry hasn't been
          perceived, the chemistry of the molecule is inconsistent, and
          "canonical" atom ranks are only a technical artifact.

    \\n\\
        - canonical: (optional) if True, uses canonical atom ranking so that the\\n\\
          kekulization result is independent of the atom ordering in the molecule.\\n\\
          Default value is False.\\n\\
    \\n\\
      NOTES:\\n\\
    \\n\\
        - The molecule is modified in place.\\n\\
    """

def Cleanup(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    cleans up certain common bad functionalities in the molecule

      ARGUMENTS:

        - mol: the molecule to use

      NOTES:

        - The molecule is modified in place.
    """

def CleanupChirality(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    removes bogus chirality markers (e.g. tetrahedral flags on non-sp3 centers)

      ARGUMENTS:

        - mol: the molecule to use

      NOTES:

        - The molecule is modified in place.
    """

def CleanupAtropisomers(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    removes bogus atropisomeric markers (e.g. those without sp2 begin and end atoms)

      ARGUMENTS:

        - mol: the molecule to use

      NOTES:

        - The molecule is modified in place.
    """

def CleanupOrganometallics(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    cleans up certain common bad functionalities in the organometallic molecule

      Note that this function is experimental and may either change in behavior
      or be replaced with something else in future releases.

            ARGUMENTS :

     - mol : the molecule to use

     NOTES :

     - The molecule is modified in place.
    """

class AromaticityModel(enum.Enum):
    AROMATICITY_DEFAULT = 0

    AROMATICITY_RDKIT = 1

    AROMATICITY_SIMPLE = 2

    AROMATICITY_MDL = 4

    AROMATICITY_MMFF94 = 8

    AROMATICITY_CUSTOM = 268435455

AROMATICITY_DEFAULT: AromaticityModel = AromaticityModel.AROMATICITY_DEFAULT

AROMATICITY_RDKIT: AromaticityModel = AromaticityModel.AROMATICITY_RDKIT

AROMATICITY_SIMPLE: AromaticityModel = AromaticityModel.AROMATICITY_SIMPLE

AROMATICITY_MDL: AromaticityModel = AromaticityModel.AROMATICITY_MDL

AROMATICITY_MMFF94: AromaticityModel = AromaticityModel.AROMATICITY_MMFF94

AROMATICITY_CUSTOM: AromaticityModel = AromaticityModel.AROMATICITY_CUSTOM

def SetAromaticity(mol: rdkit.Chem.rdchem.Mol, model: AromaticityModel = AromaticityModel.AROMATICITY_DEFAULT) -> None:
    """
    does aromaticity perception

      ARGUMENTS:

        - mol: the molecule to use
        - model: the model to use

      NOTES:

        - The molecule is modified in place.
    """

def SetConjugation(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    finds conjugated bonds

      ARGUMENTS:

        - mol: the molecule to use

      NOTES:

        - The molecule is modified in place.
    """

def SetHybridization(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    Assigns hybridization states to atoms

      ARGUMENTS:

        - mol: the molecule to use

      NOTES:

        - The molecule is modified in place.
    """

def AssignRadicals(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    Assigns radical counts to atoms

      ARGUMENTS:

        - mol: the molecule to use

      NOTES:

        - The molecule is modified in place.
    """

def HapticBondsToDative(mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
    """
    One way of showing haptic bonds (such as cyclopentadiene to
    iron in ferrocene) is to use a dummy atom with a dative bond to the
    iron atom with the bond labelled with the atoms involved in the
    organic end of the bond.  Another way is to have explicit dative
    bonds from the atoms of the haptic group to the metal atom.  This
    function converts the former representation to the latter.

    ARGUMENTS:

      - mol: the molecule to use

    RETURNS:
      a modified copy of the molecule
    """

def DativeBondsToHaptic(mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
    """
    Does the reverse of hapticBondsToDative.  If there are multiple
    contiguous atoms attached by dative bonds to an atom (probably a metal
    atom), the dative bonds will be replaced by a dummy atom in their
    centre attached to the (metal) atom by a dative bond, which is
    labelled with ENDPTS of the atoms that had the original dative bonds.

    ARGUMENTS:

      - mol: the molecule to use

    RETURNS:
      a modified copy of the molecule
    """

def FindAllSubgraphsOfLengthN(mol: rdkit.Chem.rdchem.Mol, length: int, useHs: bool = False, rootedAtAtom: int = -1) -> list[list[int]]:
    """
    Finds all subgraphs of a particular length in a molecule

      ARGUMENTS:

        - mol: the molecule to use

        - length: an integer with the target number of bonds for the subgraphs.

        - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
          should be included in the results.
          Defaults to 0.

        - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified
          atom will be returned.

      RETURNS: a tuple of 2-tuples with bond IDs

      NOTES: 

       - Difference between _subgraphs_ and _paths_ :: 

           Subgraphs are potentially branched, whereas paths (in our 
           terminology at least) cannot be.  So, the following graph: 

                C--0--C--1--C--3--C
                      |
                      2
                      |
                      C
      has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)
      but only 2 _paths_ of length 3: (0,1,3),(2,1,3)
    """

def FindAllSubgraphsOfLengthMToN(mol: rdkit.Chem.rdchem.Mol, min: int, max: int, useHs: bool = False, rootedAtAtom: int = -1) -> tuple[list[list[int]], ...]:
    """
    Finds all subgraphs of a particular length in a molecule
      See documentation for FindAllSubgraphsOfLengthN for definitions
    """

def FindUniqueSubgraphsOfLengthN(mol: rdkit.Chem.rdchem.Mol, length: int, useHs: bool = False, useBO: bool = True, rootedAtAtom: int = -1) -> list[list[int]]:
    """
    Finds unique subgraphs of a particular length in a molecule

      ARGUMENTS:

        - mol: the molecule to use

        - length: an integer with the target number of bonds for the subgraphs.

        - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
          should be included in the results.
          Defaults to 0.

        - useBO: (optional) Toggles use of bond orders in distinguishing one subgraph from
          another.
          Defaults to 1.

        - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified
          atom will be returned.

      RETURNS: a tuple of tuples with bond IDs
    """

def FindAllPathsOfLengthN(mol: rdkit.Chem.rdchem.Mol, length: int, useBonds: bool = True, useHs: bool = False, rootedAtAtom: int = -1, onlyShortestPaths: bool = False) -> list[list[int]]:
    """
    Finds all paths of a particular length in a molecule

      ARGUMENTS:

        - mol: the molecule to use

        - length: an integer with the target length for the paths.

        - useBonds: (optional) toggles the use of bond indices in the paths.
          Otherwise atom indices are used.  *Note* this behavior is different
          from that for subgraphs.
          Defaults to 1.

        - rootedAtAtom: (optional) if nonzero, only paths from the specified
          atom will be returned.

        - onlyShortestPaths: (optional) if set then only paths which are <= the shortest
          path between the begin and end atoms will be included in the results

      RETURNS: a tuple of tuples with IDs for the bonds.

      NOTES: 

       - Difference between _subgraphs_ and _paths_ :: 

           Subgraphs are potentially branched, whereas paths (in our 
           terminology at least) cannot be.  So, the following graph: 

                C--0--C--1--C--3--C
                      |
                      2
                      |
                      C

           has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)
           but only 2 _paths_ of length 3: (0,1,3),(2,1,3)
    """

def FindAtomEnvironmentOfRadiusN(mol: rdkit.Chem.rdchem.Mol, radius: int, rootedAtAtom: int, useHs: bool = False, enforceSize: bool = True, atomMap: dict[int, int] | None = None) -> list[int]:
    """
    Find bonds of a particular radius around an atom. 
             Return empty result if there is no bond at the requested radius.

      ARGUMENTS:

        - mol: the molecule to use

        - radius: an integer with the target radius for the environment.

        - rootedAtAtom: the atom to consider

        - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
          should be included in the results.
          Defaults to 0.

        - enforceSize (optional) If set to False, all bonds within the requested radius is 
          collected. Defaults to 1. 

        - atomMap: (optional) If provided, it will measure the minimum distance of the atom 
          from the rooted atom (start with 0 from the rooted atom). The result is a pair of 
          the atom ID and the distance. 

      RETURNS: a vector of bond IDs
    """

def PathToSubmol(mol: rdkit.Chem.rdchem.Mol, path: Sequence[int], useQuery: bool = False, atomMap: dict[int, int] | None = None) -> rdkit.Chem.rdchem.Mol: ...

def GetMolFrags(mol: rdkit.Chem.rdchem.Mol, asMols: bool = False, sanitizeFrags: bool = True, frags: list[int] | None = None, fragsMolAtomMapping: list[tuple[int, ...]] | None = None) -> tuple:
    """
    Finds the disconnected fragments from a molecule.

      For example, for the molecule 'CC(=O)[O-].[NH3+]C' GetMolFrags() returns
      ((0, 1, 2, 3), (4, 5))

      ARGUMENTS:

        - mol: the molecule to use
        - asMols: (optional) if this is provided and true, the fragments
          will be returned as molecules instead of atom ids.
        - sanitizeFrags: (optional) if this is provided and true, the fragments
          molecules will be sanitized before returning them.
        - frags: (optional, defaults to None) if asMols is true and this is provided
           as an empty list, the result will be mol.GetNumAtoms() long on return and
           will contain the fragment assignment for each Atom
        - fragsMolAtomMapping: (optional, defaults to None) if asMols is true and this
          is provided as an empty list, the result will be numFrags long on 
          return, and each entry will contain the indices of the Atoms in that fragment:
          [(0, 1, 2, 3), (4, 5)]

      RETURNS: a tuple of tuples with IDs for the atoms in each fragment
               or a tuple of molecules.
    """

def SplitMolByPDBResidues(mol: rdkit.Chem.rdchem.Mol, whiteList: Sequence[str] | None = None, negateList: bool = False) -> dict[str, rdkit.Chem.rdchem.Mol]:
    """
    Splits a molecule into pieces based on PDB residue information.

              ARGUMENTS:

              - mol: the molecule to use
              - whiteList: only residues in this list will be returned
              - negateList: if set, negates the white list inclusion logic

              RETURNS: a dictionary keyed by residue name with molecules as the values
    """

def SplitMolByPDBChainId(mol: rdkit.Chem.rdchem.Mol, whiteList: Sequence[str] | None = None, negateList: bool = False) -> dict[str, rdkit.Chem.rdchem.Mol]:
    """
    Splits a molecule into pieces based on PDB chain information.

      ARGUMENTS:

        - mol: the molecule to use
        - whiteList: only residues in this list will be returned
        - negateList: if set, negates the white list inclusion logic

      RETURNS: a dictionary keyed by chain id with molecules as the values
    """

def GetFormalCharge(mol: rdkit.Chem.rdchem.Mol) -> int:
    """
    Returns the formal charge for the molecule.

      ARGUMENTS:

        - mol: the molecule to use
    """

def GetShortestPath(mol: rdkit.Chem.rdchem.Mol, aid1: int, aid2: int) -> tuple[int, ...]:
    """
    Find the shortest path between two atoms using the Bellman-Ford algorithm.

      ARGUMENTS:

        - mol: the molecule to use
        - idx1: index of the first atom
        - idx2: index of the second atom
    """

def AssignStereochemistry(mol: rdkit.Chem.rdchem.Mol, cleanIt: bool = False, force: bool = False, flagPossibleStereoCenters: bool = False) -> None:
    """
    Assign stereochemistry tags to atoms and bonds.
      If useLegacyStereoPerception is true, it also does the CIP stereochemistry
      assignment for the molecule's atoms (R/S) and double bonds (Z/E).
      This assignment is based on legacy code which is fast, but is
      known to incorrectly assign CIP labels in some cases.
      instead, to assign CIP labels based on an accurate, though slower,
      implementation of the CIP rules, call CIPLabeler::assignCIPLabels().
      Chiral atoms will have a property '_CIPCode' indicating their chiral code.

      ARGUMENTS:

        - mol: the molecule to use
        - cleanIt: (optional) if provided, any existing values of the property `_CIPCode`
            will be cleared, atoms with a chiral specifier that aren't
          actually chiral (e.g. atoms with duplicate substituents or only 2 substituents,
          etc.) will have their chiral code set to CHI_UNSPECIFIED. Bonds with 
          STEREOCIS/STEREOTRANS specified that have duplicate substituents based upon the CIP 
          atom ranks will be marked STEREONONE. 
        - force: (optional) causes the calculation to be repeated, even if it has already
          been done
        - flagPossibleStereoCenters (optional)   set the _ChiralityPossible property on
          atoms that are possible stereocenters
    """

def ComputeAtomCIPRanks(mol: rdkit.Chem.rdchem.Mol) -> tuple[int, ...]:
    """
    Computes the CIP ranks for the atoms in a molecule.
      The ranks are stored as an atom property '_CIPRank' and returned as a tuple.
    """

def AssignChiralTypesFromBondDirs(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, replaceExistingTags: bool = True) -> None:
    """
    Uses bond directions to assign ChiralTypes to a molecule's atoms.

      ARGUMENTS:

        - mol: the molecule to use
        - confId: (optional) the conformation to use 
        - replaceExistingTags: (optional) replace any existing information about stereochemistry
    """

def AssignStereochemistryFrom3D(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, replaceExistingTags: bool = True) -> None:
    """
    Uses a conformer (should be 3D) to assign ChiralTypes to a molecule's atoms
            and stereo flags to its bonds

      ARGUMENTS:

        - mol: the molecule to use
        - confId: (optional) the conformation to use 
        - replaceExistingTags: (optional) replace any existing information about stereochemistry
    """

def FindPotentialStereoBonds(mol: rdkit.Chem.rdchem.Mol, cleanIt: bool = False) -> None:
    """
    Find bonds than can be cis/trans in a molecule and mark them as 'any'.
             This function finds any double bonds that can potentially be part
             of a cis/trans system. No attempt is made here to mark them cis or trans

      ARGUMENTS:

        - mol: the molecule to use
        - cleanIt: (optional) if this option is set to true, any previous marking of _CIPCode
                   on the bond is cleared - otherwise it is left untouched
    """

def RemoveStereochemistry(mol: rdkit.Chem.rdchem.Mol) -> None:
    """Removes all stereochemistry info from the molecule."""

def AssignAtomChiralTagsFromStructure(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, replaceExistingTags: bool = True) -> None:
    """
    Sets the chiral tags on a molecule's atoms based on
      a 3D conformation.
      NOTE that this does not check to see if atoms are chiral centers (i.e. all
      substituents are different), it merely sets the chiral type flags based on the
      coordinates and atom ordering. Use AssignStereochemistryFrom3D() if you
      want chiral flags only on actual stereocenters.

      ARGUMENTS:

        - mol: the molecule to use
        - confId: the conformer id to use, -1 for the default 
        - replaceExistingTags: if True, existing stereochemistry information will be cleared
        before running the calculation.
    """

def AssignAtomChiralTagsFromMolParity(mol: rdkit.Chem.rdchem.Mol, replaceExistingTags: bool = True) -> None:
    """
    Sets the chiral tags on a molecule's atoms based on
      the molParity atom property.

      ARGUMENTS:

        - mol: the molecule to use
        - replaceExistingTags: if True, existing stereochemistry information will be cleared
        before running the calculation.
    """

def FindMesoCenters(mol: rdkit.Chem.rdchem.Mol, includeIsotopes: bool = True, includeAtomMaps: bool = False) -> tuple[tuple[int, int], ...]:
    """
    returns the meso centers in a molecule (if any).

      ARGUMENTS:

        - mol: the molecule to use
        - includeIsotopes: (optional) toggles whether or not isotopes should be included in the
          calculation.
        - includeAtomMaps: (optional) toggles whether or not atom maps should be included in the
          calculation.
    """

def RDKFingerprint(mol: rdkit.Chem.rdchem.Mol, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, nBitsPerHash: int = 2, useHs: bool = True, tgtDensity: float = 0.0, minSize: int = 128, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: Iterable[int] | None = None, fromAtoms: Iterable[int] | None = None, atomBits: list[list[int]] | None = None, bitInfo: dict[int, list[list[int]]] | None = None) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """
    Returns an RDKit topological fingerprint for a molecule

      Explanation of the algorithm below.

      ARGUMENTS:

        - mol: the molecule to use

        - minPath: (optional) minimum number of bonds to include in the subgraphs
          Defaults to 1.

        - maxPath: (optional) maximum number of bonds to include in the subgraphs
          Defaults to 7.

        - fpSize: (optional) number of bits in the fingerprint
          Defaults to 2048.

        - nBitsPerHash: (optional) number of bits to set per path
          Defaults to 2.

        - useHs: (optional) include paths involving Hs in the fingerprint if the molecule
          has explicit Hs.
          Defaults to True.

        - tgtDensity: (optional) fold the fingerprint until this minimum density has
          been reached
          Defaults to 0.

        - minSize: (optional) the minimum size the fingerprint will be folded to when
          trying to reach tgtDensity
          Defaults to 128.

        - branchedPaths: (optional) if set both branched and unbranched paths will be
          used in the fingerprint.
          Defaults to True.

        - useBondOrder: (optional) if set both bond orders will be used in the path hashes
          Defaults to True.

        - atomInvariants: (optional) a sequence of atom invariants to use in the path hashes
          Defaults to empty.

        - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs 
          starting from these atoms will be used.
          Defaults to empty.

        - atomBits: (optional) an empty list. If provided, the result will contain a list 
          containing the bits each atom sets.
          Defaults to empty.

        - bitInfo: (optional) an empty dict. If provided, the result will contain a dict 
          with bits as keys and corresponding bond paths as values.
          Defaults to empty.

      RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits

      ALGORITHM:

       This algorithm functions by find all subgraphs between minPath and maxPath in
       length.  For each subgraph:

         1) A hash is calculated.

         2) The hash is used to seed a random-number generator

         3) _nBitsPerHash_ random numbers are generated and used to set the corresponding
            bits in the fingerprint
    """

def UnfoldedRDKFingerprintCountBased(mol: rdkit.Chem.rdchem.Mol, minPath: int = 1, maxPath: int = 7, useHs: bool = True, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: Iterable[int] | None = None, fromAtoms: Iterable[int] | None = None, atomBits: list[list[int]] | None = None, bitInfo: dict[int, list[list[int]]] | None = None) -> rdkit.DataStructs.cDataStructs.ULongSparseIntVect:
    """
    Returns an unfolded count-based version of the RDKit fingerprint for a molecule

    ARGUMENTS:

            - mol: the molecule to use

            - minPath: (optional) minimum number of bonds to include in the subgraphs
              Defaults to 1.

            - maxPath: (optional) maximum number of bonds to include in the subgraphs
              Defaults to 7.

            - useHs: (optional) include paths involving Hs in the fingerprint if the molecule
              has explicit Hs.
              Defaults to True.

            - branchedPaths: (optional) if set both branched and unbranched paths will be
              used in the fingerprint.
              Defaults to True.

            - useBondOrder: (optional) if set both bond orders will be used in the path hashes
              Defaults to True.

            - atomInvariants: (optional) a sequence of atom invariants to use in the path hashes
              Defaults to empty.

            - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs 
              starting from these atoms will be used.
              Defaults to empty.

            - atomBits: (optional) an empty list. If provided, the result will contain a list 
              containing the bits each atom sets.
              Defaults to empty.

            - bitInfo: (optional) an empty dict. If provided, the result will contain a dict 
              with bits as keys and corresponding bond paths as values.
              Defaults to empty.
    """

def LayeredFingerprint(mol: rdkit.Chem.rdchem.Mol, layerFlags: int = 4294967295, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, atomCounts: object | None = None, setOnlyBits: rdkit.DataStructs.cDataStructs.ExplicitBitVect | None = None, branchedPaths: bool = True, fromAtoms: Iterable[int] | None = None) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """
    Returns a layered fingerprint for a molecule

      NOTE: This function is experimental. The API or results may change from
        release to release.

      Explanation of the algorithm below.

      ARGUMENTS:

        - mol: the molecule to use

        - layerFlags: (optional) which layers to include in the fingerprint
          See below for definitions. Defaults to all.

        - minPath: (optional) minimum number of bonds to include in the subgraphs
          Defaults to 1.

        - maxPath: (optional) maximum number of bonds to include in the subgraphs
          Defaults to 7.

        - fpSize: (optional) number of bits in the fingerprint
          Defaults to 2048.

        - atomCounts: (optional) 
          if provided, this should be a list at least as long as the number of atoms
          in the molecule. It will be used to provide the count of the number 
          of paths that set bits each atom is involved in.
          NOTE: the list is not zeroed out here.

        - setOnlyBits: (optional) 
          if provided, only bits that are set in this bit vector will be set
          in the result. This is essentially the same as doing:
               res &= setOnlyBits
          but also has an impact on the atomCounts (if being used)

        - branchedPaths: (optional) if set both branched and unbranched paths will be
          used in the fingerprint.
          Defaults to True.

        - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs 
          starting from these atoms will be used.
          Defaults to empty.

      RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits

      Layer definitions:
         - 0x01: pure topology
         - 0x02: bond order
         - 0x04: atom types
         - 0x08: presence of rings
         - 0x10: ring sizes
         - 0x20: aromaticity
    """

LayeredFingerprint_substructLayers: int = 7

@overload
def PatternFingerprint(mol: rdkit.Chem.rdchem.Mol, fpSize: int = 2048, atomCounts: object | None = None, setOnlyBits: rdkit.DataStructs.cDataStructs.ExplicitBitVect | None = None, tautomerFingerprints: bool = False) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect: ...

@overload
def PatternFingerprint(mol: rdkit.Chem.rdchem.MolBundle, fpSize: int = 2048, setOnlyBits: rdkit.DataStructs.cDataStructs.ExplicitBitVect | None = None, tautomerFingerprints: bool = False) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """
    A fingerprint using SMARTS patterns 

      NOTE: This function is experimental. The API or results may change from
        release to release.
    """

class BondWedgingParameters:
    """Parameters controlling how bond wedging is done."""

    @property
    def wedgeTwoBondsIfPossible(self) -> bool:
        """
        If this is enabled then two bonds will be wedged at chiral
          centers subject to the following constraints:
            1. ring bonds will not be wedged
            2. bonds to chiral centers will not be wedged
            3. bonds separated by more than 120 degrees will not be
                wedged
        """

    @wedgeTwoBondsIfPossible.setter
    def wedgeTwoBondsIfPossible(self, arg: bool, /) -> None: ...

def WedgeMolBonds(mol: rdkit.Chem.rdchem.Mol, conformer: rdkit.Chem.rdchem.Conformer, params: BondWedgingParameters | None = None) -> None:
    """
    Set the wedging on single bonds in a molecule.
       The wedging scheme used is that from Mol files.

      ARGUMENTS:

        - molecule: the molecule to update
        - conformer: the conformer to use to determine wedge direction
    """

def ReapplyMolBlockWedging(mol: rdkit.Chem.rdchem.Mol, allBondTypes: bool = True, verify: bool = False) -> None:
    """
    Set the wedging to that which was read from the original
         MolBlock, over-riding anything that was originally there.

              ARGUMENTS:

                - molecule: the molecule to update
                - allBondTypes: reapply the wedging also on bonds other
                  than single and aromatic ones
                - verify: if true, the function will check that the wedges are only
                  applied in sensible places (i.e.single bonds connected to chiral
                  centers or atropisomeric bonds)
    """

def RemoveNonExplicit3DChirality(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    Remove chiral markings that were derived from a 3D mol but were not 
            explicity marked in the mol block. (wedge bond or CFG indication

              ARGUMENTS:

                - molecule: the molecule to update
    """

class StereoGroupAbsOptions(enum.Enum):
    OnlyIncludeWhenOtherGroupsExist = 0

    NeverInclude = 1

    AlwaysInclude = 2

def CanonicalizeStereoGroups(mol: rdkit.Chem.rdchem.Mol, outputAbsoluteGroups: StereoGroupAbsOptions = StereoGroupAbsOptions.OnlyIncludeWhenOtherGroupsExist, maxStereoGroups: int = 12) -> rdkit.Chem.rdchem.Mol:
    """
    Rationalize Enhanced Stereo indications to a canonical form 

              ARGUMENTS:

                - molecule: the molecule to update
                -StereoGroupAbsOptions outputAbsoluteGroups: controls output of abs groups: 
                  one of: OnlyIncludeWhenOtherGroupsExist, NeverInclude, AlwaysInclude 
                 maxStereoGroups: maximm number of OR or AND stereo groups to process (default is 12):
    """

class StereoBondThresholds:
    """
    Constants used to set the thresholds for which single bonds can be made wavy.
    """

    DBL_BOND_NO_STEREO: Final[int] = ...
    """neighboring double bond without stereo info"""

    DBL_BOND_SPECIFIED_STEREO: Final[int] = ...
    """neighboring double bond with stereo specified"""

    CHIRAL_ATOM: Final[int] = ...
    """atom with specified chirality"""

    DIRECTION_SET: Final[int] = ...
    """single bond with the direction already set"""

def AddWavyBondsForStereoAny(mol: rdkit.Chem.rdchem.Mol, clearDoubleBondFlags: bool = True, addWhenImpossible: int = 1000) -> None:
    """
    set wavy bonds around double bonds with STEREOANY stereo
      ARGUMENTS :
        - molecule : the molecule to update\\n -
        - conformer : the conformer to use to determine wedge direction
    """

def WedgeBond(bond: rdkit.Chem.rdchem.Bond, fromAtomIdx: int, conf: rdkit.Chem.rdchem.Conformer) -> None:
    """
    Set the wedging on an individual bond from a molecule.
       The wedging scheme used is that from Mol files.
      ARGUMENTS:
        - bond: the bond to update
        - atom ID: the atom from which to do the wedging
        - conformer: the conformer to use to determine wedge direction
    """

def ReplaceSidechains(mol: rdkit.Chem.rdchem.Mol, coreQuery: rdkit.Chem.rdchem.Mol, useChirality: bool = False) -> rdkit.Chem.rdchem.Mol | None:
    """
    Replaces sidechains in a molecule with dummy atoms for their attachment points.

      ARGUMENTS:

        - mol: the molecule to be modified

        - coreQuery: the molecule to be used as a substructure query for recognizing the core

        - useChirality: (optional) match the substructure query using chirality

      RETURNS: a new molecule with the sidechains removed

      NOTES:

        - The original molecule is *not* modified.

      EXAMPLES:

       The following examples substitute SMILES/SMARTS strings for molecules, you'd have
       to actually use molecules:

        - ReplaceSidechains('CCC1CCC1','C1CCC1') -> '[Xa]C1CCC1'

        - ReplaceSidechains('CCC1CC1','C1CCC1') -> ''

        - ReplaceSidechains('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'
    """

@overload
def ReplaceCore(mol: rdkit.Chem.rdchem.Mol, core: rdkit.Chem.rdchem.Mol, matches: object, replaceDummies: bool = True, labelByIndex: bool = False, requireDummyMatch: bool = False) -> rdkit.Chem.rdchem.Mol | None:
    """
    Removes the core of a molecule and labels the sidechains with dummy atoms based on
    The matches indices given in the matching vector matches.
    Calling:
      ReplaceCore(mol,core,mol.GetSubstructMatch(core))

      ARGUMENTS:

        - mol: the molecule to be modified

        - coreQuery: the molecule to be used as a substructure query for recognizing the core

        - matches: a matching vector of the type returned by mol.GetSubstructMatch(...)

        - replaceDummies: toggles replacement of atoms that match dummies in the query

        - labelByIndex: toggles labeling the attachment point dummy atoms with 
          the index of the core atom they're attached to.

        - requireDummyMatch: if the molecule has side chains that attach at points not
          flagged with a dummy, it will be rejected (None is returned)

      RETURNS: a new molecule with the core removed

      NOTES:

        - The original molecule is *not* modified.
    EXAMPLES:

        >>> from rdkit.Chem import MolToSmiles, MolFromSmiles, ReplaceCore
        >>> mol = MolFromSmiles('C1ONNCC1')
        >>> core = MolFromSmiles('NN')

        >>> MolToSmiles(ReplaceCore(mol, core, mol.GetSubstructMatch(core)))
        '[1*]OCCC[2*]'

        Since NN is symmetric, we should actually get two matches here if we don't
        uniquify the matches.

        >>> [MolToSmiles(ReplaceCore(mol, core, match))
        ...     for match in mol.GetSubstructMatches(core, uniquify=False)]
        ['[1*]OCCC[2*]', '[1*]CCCO[2*]']
    """

@overload
def ReplaceCore(mol: rdkit.Chem.rdchem.Mol, coreQuery: rdkit.Chem.rdchem.Mol, replaceDummies: bool = True, labelByIndex: bool = False, requireDummyMatch: bool = False, useChirality: bool = False) -> rdkit.Chem.rdchem.Mol | None:
    """
    Removes the core of a molecule and labels the sidechains with dummy atoms.

      ARGUMENTS:

        - mol: the molecule to be modified

        - coreQuery: the molecule to be used as a substructure query for recognizing the core

        - replaceDummies: toggles replacement of atoms that match dummies in the query

        - labelByIndex: toggles labeling the attachment point dummy atoms with 
          the index of the core atom they're attached to.

        - requireDummyMatch: if the molecule has side chains that attach at points not
          flagged with a dummy, it will be rejected (None is returned)

        - useChirality: use chirality matching in the coreQuery

      RETURNS: a new molecule with the core removed

      NOTES:

        - The original molecule is *not* modified.

      EXAMPLES:

       >>> from rdkit.Chem import MolToSmiles, MolFromSmiles, MolFromSmarts, ReplaceCore

       Basic usage: remove a core as specified by SMILES (or another molecule).
       To get the atom labels which are stored as an isotope of the matched atom, 
       the output must be written as isomeric smiles.  
       A small confusion is that atom isotopes of 0 aren't shown in smiles strings.

       Here we remove a ring and leave the decoration (r-group) behind.

       >>> MolToSmiles(ReplaceCore(MolFromSmiles('CCCC1CCC1'),MolFromSmiles('C1CCC1')),
       ...             isomericSmiles=True)
       '[1*]CCC'

       The isotope label by default is matched by the first connection found. In order to
       indicate which atom the decoration is attached in the core query, use labelByIndex=True.
       Here the attachment is from the third atom in the smiles string, which is indexed by 3
       in the core, like all good computer scientists expect, atoms indices start at 0.

       >>> MolToSmiles(ReplaceCore(MolFromSmiles('CCN1CCC1'),MolFromSmiles('C1CCN1'),
       ...                         labelByIndex=True),
       ...   isomericSmiles=True)
       '[3*]CC'

       Non-core matches just return None

       >>> ReplaceCore(MolFromSmiles('CCC1CC1'),MolFromSmiles('C1CCC1'))

       The bond between atoms are considered part of the core and are removed as well

       >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CC2C1CCC2'),MolFromSmiles('C1CCC1')),
       ...             isomericSmiles=True)
       '[1*]CCC[2*]'
       >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmiles('N')),
       ...             isomericSmiles=True)
       '[1*]CCCC[2*]'

       When using dummy atoms, cores should be read in as SMARTS.  When read as SMILES
       dummy atoms only match other dummy atoms.
       The replaceDummies flag indicates whether matches to the dummy atoms should be considered as part
       of the core or as part of the decoration (r-group)

       >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmarts('[*]N[*]'),
       ...                         replaceDummies=True),
       ...             isomericSmiles=True)
       '[1*]CC[2*]'
       >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmarts('[*]N[*]'),
       ...                         replaceDummies=False),
       ...             isomericSmiles=True)
       '[1*]CCCC[2*]'


       >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CCC1CN'),MolFromSmarts('C1CCC1[*]'),
       ...                         replaceDummies=False),
       ...             isomericSmiles=True)
       '[1*]CN'
    """

def FragmentOnBRICSBonds(mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
    """Return a new molecule with all BRICS bonds broken"""

def FragmentOnBonds(mol: rdkit.Chem.rdchem.Mol, bondIndices: Iterable[int], addDummies: bool = True, dummyLabels: Sequence[Sequence[int]] | None = None, bondTypes: Sequence[rdkit.Chem.rdchem.BondType] | None = None, cutsPerAtom: object | None = None) -> rdkit.Chem.rdchem.Mol:
    """
    Return a new molecule with all specified bonds broken

      ARGUMENTS:

          - mol            - the molecule to be modified
          - bondIndices    - indices of the bonds to be broken
          - addDummies  - toggles addition of dummy atoms to indicate where 
            bonds were broken
          - dummyLabels - used to provide the labels to be used for the dummies.
            the first element in each pair is the label for the dummy
            that replaces the bond's beginAtom, the second is for the 
            dummy that replaces the bond's endAtom. If not provided, the
            dummies are labeled with atom indices.
          - bondTypes - used to provide the bond type to use between the
            fragments and the dummy atoms. If not provided, defaults to single. 
          - cutsPerAtom - used to return the number of cuts made at each atom. 

      RETURNS:
          a new Mol with the modifications
    """

def FragmentOnSomeBonds(mol: rdkit.Chem.rdchem.Mol, bondIndices: Iterable[int], numToBreak: int = 1, addDummies: bool = True, dummyLabels: Sequence[Sequence[int]] | None = None, bondTypes: Sequence[rdkit.Chem.rdchem.BondType] | None = None, returnCutsPerAtom: bool = False) -> tuple:
    """fragment on some bonds"""

class MolzipLabel(enum.Enum):
    AtomMapNumber = 0

    Isotope = 1

    FragmentOnBonds = 2

    AtomType = 3

class MolzipParams:
    """
    Parameters controlling how to zip molecules together

      OPTIONS:
          label : set the MolzipLabel option [default MolzipLabel.AtomMapNumber]

      MolzipLabel.AtomMapNumber: atom maps are on dummy atoms, zip together the corresponding
         attached atoms, i.e.  zip 'C[*:1]' 'N[*:1]' results in 'CN'

      MolzipLabel.Isotope: isotope labels are on dummy atoms, zip together the corresponding
         attached atoms, i.e.  zip 'C[1*]' 'N[1*]' results in 'CN'

      MolzipLabel.FragmentOnBonds: zip together molecules generated by fragment on bonds.
        Note the atom indices cannot change or be reordered from the output of fragmentOnBonds

      MolzipLabel.AtomTypes: choose the atom types to act as matching dummy atoms.
        i.e.  'C[V]' and 'N[Xe]' with atoms pairs [('V', 'Xe')] results in 'CN'
    """

    def __init__(self) -> None: ...

    @property
    def label(self) -> MolzipLabel:
        """Set the atom labeling system to zip together"""

    @label.setter
    def label(self, arg: MolzipLabel, /) -> None: ...

    @property
    def enforceValenceRules(self) -> bool:
        """
        If true (default) enforce valences after zipping
        Setting this to false allows assembling chemically incorrect fragments.
        """

    @enforceValenceRules.setter
    def enforceValenceRules(self, arg: bool, /) -> None: ...

    @property
    def generateCoordinates(self) -> bool:
        """
        If true will add depiction coordinates to input molecules and
        zipped molecule (for molzipFragments only)
        """

    @generateCoordinates.setter
    def generateCoordinates(self, arg: bool, /) -> None: ...

    @property
    def alignCoordinates(self) -> bool:
        """
        if true and the input fragments have coordinates, the fragments
        will be aligned along connection vectors in the output molecule
        """

    @alignCoordinates.setter
    def alignCoordinates(self, arg: bool, /) -> None: ...

    def setAtomSymbols(self, symbols: Sequence[str] | None) -> None:
        """
        Set the atom symbols used to zip mols together when using AtomType labeling
        """

    def __setattr__(self, name: str, value: object | None) -> None: ...

@overload
def molzip(a: rdkit.Chem.rdchem.Mol, b: rdkit.Chem.rdchem.Mol, params: MolzipParams | None = None) -> rdkit.Chem.rdchem.Mol | None:
    """zip together two molecules using the given matching parameters"""

@overload
def molzip(a: rdkit.Chem.rdchem.Mol, params: MolzipParams | None = None) -> rdkit.Chem.rdchem.Mol | None:
    """
    zip together multiple molecules within a combined molecule using the given matching parameters
    """

@overload
def molzip(row: dict[str, rdkit.Chem.rdchem.Mol], params: MolzipParams | None = None) -> rdkit.Chem.rdchem.Mol | None:
    """
    zip an RGroupRow together to recreate the original molecule.  This correctly handles
    broken cycles that can occur in decompositions.
     example:

      >>> from rdkit import Chem
      >>> from rdkit.Chem import rdRGroupDecomposition as rgd
      >>> core = Chem.MolFromSmiles('CO')
      >>> mols = [Chem.MolFromSmiles('C1NNO1')]
      >>> rgroups, unmatched = rgd.RGroupDecompose(core, mols)
      >>> for rgroup in rgroups:
      ...     mol = rgd.molzip(rgroup)
    """

def molzipFragments(mols: Iterable[rdkit.Chem.rdchem.Mol] | None, params: MolzipParams | None = None) -> rdkit.Chem.rdchem.Mol | None:
    """
    zip together multiple molecules from an R group decomposition 
    using the given matching parameters.  The first molecule in the list
    must be the core
    """

def AddRecursiveQuery(mol: rdkit.Chem.rdchem.Mol, query: rdkit.Chem.rdchem.Mol, atomIdx: int, preserveExistingQuery: bool = True) -> None:
    """
    Adds a recursive query to an atom

      ARGUMENTS:

        - mol: the molecule to be modified

        - query: the molecule to be used as the recursive query (this will be copied)

        - atomIdx: the atom to modify

        - preserveExistingQuery: (optional) if this is set, existing query information on the atom will be preserved

      RETURNS: None
    """

def RenumberAtoms(mol: rdkit.Chem.rdchem.Mol, newOrder: Iterable[int]) -> rdkit.Chem.rdchem.Mol:
    """
    Returns a copy of a molecule with renumbered atoms

      ARGUMENTS:

        - mol: the molecule to be modified

        - newOrder: the new ordering the atoms (should be numAtoms long)
          for example: if newOrder is [3,2,0,1], then atom 3 in the original 
          molecule will be atom 0 in the new one
    """

class AdjustQueryWhichFlags(enum.IntEnum):
    ADJUST_IGNORENONE = 0

    ADJUST_IGNORECHAINS = 1

    ADJUST_IGNORERINGS = 4

    ADJUST_IGNOREDUMMIES = 2

    ADJUST_IGNORENONDUMMIES = 8

    ADJUST_IGNOREMAPPED = 16

    ADJUST_IGNOREALL = 268435455

ADJUST_IGNORENONE: AdjustQueryWhichFlags = AdjustQueryWhichFlags.ADJUST_IGNORENONE

ADJUST_IGNORECHAINS: AdjustQueryWhichFlags = AdjustQueryWhichFlags.ADJUST_IGNORECHAINS

ADJUST_IGNORERINGS: AdjustQueryWhichFlags = AdjustQueryWhichFlags.ADJUST_IGNORERINGS

ADJUST_IGNOREDUMMIES: AdjustQueryWhichFlags = AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES

ADJUST_IGNORENONDUMMIES: AdjustQueryWhichFlags = AdjustQueryWhichFlags.ADJUST_IGNORENONDUMMIES

ADJUST_IGNOREMAPPED: AdjustQueryWhichFlags = AdjustQueryWhichFlags.ADJUST_IGNOREMAPPED

ADJUST_IGNOREALL: AdjustQueryWhichFlags = AdjustQueryWhichFlags.ADJUST_IGNOREALL

class AdjustQueryParameters:
    """
    Parameters controlling which components of the query atoms/bonds are adjusted.

    Note that some of the options here are either directly contradictory or make
      no sense when combined with each other. We generally assume that client code
      is doing something sensible and don't attempt to detect possible conflicts or
      problems.

    A note on the flags controlling which atoms/bonds are modified: 
       These generally limit the set of atoms/bonds to be modified.
       For example:
           - ADJUST_IGNORERINGS atoms/bonds in rings will not be modified.
           - ADJUST_IGNORENONE causes all atoms/bonds to be modified
           - ADJUST_IGNOREALL no atoms/bonds will be modified
       Some of the options obviously make no sense for bonds
    """

    def __init__(self) -> None: ...

    @property
    def adjustDegree(self) -> bool:
        """add degree queries"""

    @adjustDegree.setter
    def adjustDegree(self, arg: bool, /) -> None: ...

    @property
    def adjustDegreeFlags(self) -> int:
        """controls which atoms have their degree queries changed"""

    @adjustDegreeFlags.setter
    def adjustDegreeFlags(self, arg: int, /) -> None: ...

    @property
    def adjustHeavyDegree(self) -> bool:
        """adjust the heavy-atom degree"""

    @adjustHeavyDegree.setter
    def adjustHeavyDegree(self, arg: bool, /) -> None: ...

    @property
    def adjustHeavyDegreeFlags(self) -> int:
        """controls which atoms have their heavy-atom degree queries changed"""

    @adjustHeavyDegreeFlags.setter
    def adjustHeavyDegreeFlags(self, arg: int, /) -> None: ...

    @property
    def adjustRingCount(self) -> bool:
        """add ring-count queries"""

    @adjustRingCount.setter
    def adjustRingCount(self, arg: bool, /) -> None: ...

    @property
    def adjustRingCountFlags(self) -> int:
        """controls which atoms have ring-count queries added"""

    @adjustRingCountFlags.setter
    def adjustRingCountFlags(self, arg: int, /) -> None: ...

    @property
    def makeDummiesQueries(self) -> bool:
        """convert dummy atoms without isotope labels to any-atom queries"""

    @makeDummiesQueries.setter
    def makeDummiesQueries(self, arg: bool, /) -> None: ...

    @property
    def aromatizeIfPossible(self) -> bool:
        """perceive and set aromaticity"""

    @aromatizeIfPossible.setter
    def aromatizeIfPossible(self, arg: bool, /) -> None: ...

    @property
    def makeBondsGeneric(self) -> bool:
        """converts bonds to generic queries (any bonds)"""

    @makeBondsGeneric.setter
    def makeBondsGeneric(self, arg: bool, /) -> None: ...

    @property
    def makeBondsGenericFlags(self) -> int:
        """controls which bonds are converted to generic queries"""

    @makeBondsGenericFlags.setter
    def makeBondsGenericFlags(self, arg: int, /) -> None: ...

    @property
    def makeAtomsGeneric(self) -> bool:
        """convert atoms to generic queries (any atoms)"""

    @makeAtomsGeneric.setter
    def makeAtomsGeneric(self, arg: bool, /) -> None: ...

    @property
    def makeAtomsGenericFlags(self) -> int:
        """controls which atoms are converted to generic queries"""

    @makeAtomsGenericFlags.setter
    def makeAtomsGenericFlags(self, arg: int, /) -> None: ...

    @property
    def adjustRingChain(self) -> bool:
        """add ring-chain queries to atoms"""

    @adjustRingChain.setter
    def adjustRingChain(self, arg: bool, /) -> None: ...

    @property
    def adjustRingChainFlags(self) -> int:
        """controls which atoms have ring-chain queries added"""

    @adjustRingChainFlags.setter
    def adjustRingChainFlags(self, arg: int, /) -> None: ...

    @property
    def useStereoCareForBonds(self) -> bool:
        """
        if this is set sterochemistry information will be removed from double bonds that do not have the stereoCare property set
        """

    @useStereoCareForBonds.setter
    def useStereoCareForBonds(self, arg: bool, /) -> None: ...

    @property
    def adjustConjugatedFiveRings(self) -> bool:
        """set bond queries in conjugated five-rings to SINGLE|DOUBLE|AROMATIC"""

    @adjustConjugatedFiveRings.setter
    def adjustConjugatedFiveRings(self, arg: bool, /) -> None: ...

    @property
    def setMDLFiveRingAromaticity(self) -> bool:
        """
        uses the 5-ring aromaticity behavior of the (former) MDL software as documented in the Chemical Representation Guide
        """

    @setMDLFiveRingAromaticity.setter
    def setMDLFiveRingAromaticity(self, arg: bool, /) -> None: ...

    @property
    def adjustSingleBondsToDegreeOneNeighbors(self) -> bool:
        """
        set single bonds bewteen aromatic or conjugated atoms and degree-one neighbors to SINGLE|AROMATIC
        """

    @adjustSingleBondsToDegreeOneNeighbors.setter
    def adjustSingleBondsToDegreeOneNeighbors(self, arg: bool, /) -> None: ...

    @property
    def adjustSingleBondsBetweenAromaticAtoms(self) -> bool:
        """
        sets non-ring single bonds between two aromatic or conjugated atoms to SINGLE|AROMATIC
        """

    @adjustSingleBondsBetweenAromaticAtoms.setter
    def adjustSingleBondsBetweenAromaticAtoms(self, arg: bool, /) -> None: ...

    @staticmethod
    def NoAdjustments() -> AdjustQueryParameters:
        """
        Returns an AdjustQueryParameters object with all parameters set to false
        """

def AdjustQueryProperties(mol: rdkit.Chem.rdchem.Mol, params: AdjustQueryParameters | None = None) -> rdkit.Chem.rdchem.Mol:
    """
    Returns a new molecule where the query properties of atoms have been modified.
    """

def AdjustQueryPropertiesWithGenericGroups(mol: rdkit.Chem.rdchem.Mol, params: AdjustQueryParameters | None = None) -> rdkit.Chem.rdchem.Mol:
    """
    Returns a new molecule where the query properties of atoms have been modified and generic group queries have been prepared.
    """

def DetectChemistryProblems(mol: rdkit.Chem.rdchem.Mol, sanitizeOps: int = SanitizeFlags.SANITIZE_ALL) -> tuple[rdkit.Chem.rdchem._cppMolSanitizeException, ...]:
    """checks for chemistry problems"""

def SetGenericQueriesFromProperties(mol: rdkit.Chem.rdchem.Mol, useAtomLabels: bool = True, useSGroups: bool = True) -> None:
    """documentation"""

def ConvertGenericQueriesToSubstanceGroups(mol: rdkit.Chem.rdchem.Mol) -> None:
    """documentation"""

def SetAllowNontetrahedralChirality(val: bool) -> None:
    """toggles recognition of non-tetrahedral chirality from 3D structures"""

def GetAllowNontetrahedralChirality() -> bool:
    """
    returns whether or not recognition of non-tetrahedral chirality from 3D structures is enabled
    """

def SetUseLegacyStereoPerception(val: bool) -> None:
    """sets usage of the legacy stereo perception code"""

def GetUseLegacyStereoPerception() -> bool:
    """returns whether or not the legacy stereo perception code is being used"""

def TranslateChiralFlagToStereoGroups(mol: rdkit.Chem.rdchem.Mol, zeroFlagGroupType: rdkit.Chem.rdchem.StereoGroupType = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND) -> None:
    """
    Generate enhanced stereo groups based on the status of the chiral flag property.

                Arguments:
                - mol: molecule to be modified
                - zeroFlagGroupType: how to handle non-grouped stereo centers when the
                chiral flag is set to zero

                If the chiral flag is set to a value of 1 then all specified tetrahedral
                chiral centers which are not already in StereoGroups will be added to an
                ABS StereoGroup.

                If the chiral flag is set to a value of 0 then all specified tetrahedral
                chiral centers will be added to a StereoGroup of the type zeroFlagGroupType

                If there is no chiral flag set (i.e. the property is not present), the
                molecule will not be modified.
    """

def ExpandAttachmentPoints(mol: rdkit.Chem.rdchem.Mol, addAsQueries: bool = True, addCoords: bool = True) -> None:
    """
    attachment points encoded as attachPt properties are added to the graph as dummy atoms

      Arguments:
       - mol: molecule to be modified
       - addAsQueries: if true, the dummy atoms will be added as null queries
            (i.e. they will match any atom in a substructure search)
       - addCoords: if true and the molecule has one or more conformers, 
            positions for the attachment points will be added to the conformer(s)
    """

def CollapseAttachmentPoints(mol: rdkit.Chem.rdchem.Mol, markedOnly: bool = True) -> None:
    """
    dummy atoms in the graph are removed and replaced with attachment point annotations on the attached atoms

      Arguments:
       - mol: molecule to be modified
       - markedOnly: if true, only dummy atoms with the _fromAttachPoint
         property will be collapsed

      In order for a dummy atom to be considered for collapsing it must have:
       - degree 1 with a single or unspecified bond
       - the bond to it can not be wedged
       - either no query or be an AtomNullQuery
    """

def AddStereoAnnotations(mol: rdkit.Chem.rdchem.Mol, absLabel: str = 'abs ({cip})', orLabel: str = 'or{id}', andLabel: str = 'and{id}', cipLabel: str = '({cip})', bondLabel: str = '({cip})') -> None:
    """
    add R/S, relative stereo, and E/Z annotations to atoms and bonds

      Arguments:
       - mol: molecule to modify
       - absLabel: label for atoms in an ABS stereo group
       - orLabel: label for atoms in an OR stereo group
       - andLabel: label for atoms in an AND stereo group
       - cipLabel: label for chiral atoms that aren't in a stereo group.
       - bondLabel: label for CIP stereochemistry on bonds

     If any label is empty, the corresponding annotations will not be added.

     The labels can contain the following placeholders:
       - {id} - the stereo group's index
       - {cip} - the atom or bond's CIP stereochemistry

     Note that CIP labels will only be added if CIP stereochemistry has been
     assigned to the molecule.
    """

def SimplifyEnhancedStereo(mol: rdkit.Chem.rdchem.Mol, removeAffectedStereoGroups: bool = True) -> None:
    """
    Simplifies the stereochemical representation of a molecule where all
    specified stereocenters are in the same StereoGroup

      Arguments:
       - mol: molecule to modify
       - removeAffectedStereoGroups: if set then the affected StereoGroups will be removed

    If all specified stereocenters are in the same AND or OR stereogroup, a
    moleculeNote property will be set on the molecule with the value "AND
    enantiomer" or "OR enantiomer". CIP labels, if present, are removed.
    """

def NeedsHs(mol: rdkit.Chem.rdchem.Mol) -> bool:
    """returns whether or not the molecule needs to have Hs added"""

def CountAtomElec(atom: rdkit.Chem.rdchem.Atom) -> int:
    """
    returns the number of electrons available on an atom to donate for aromaticity
    """

def AtomHasConjugatedBond(atom: rdkit.Chem.rdchem.Atom) -> bool:
    """returns whether or not the atom is involved in a conjugated bond"""

class SubsetMethod(enum.Enum):
    BONDS_BETWEEN_ATOMS = 0

    BONDS = 1

class BoolVector:
    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, arg: BoolVector) -> None:
        """Copy constructor"""

    @overload
    def __init__(self, arg: Iterable[bool], /) -> None:
        """Construct from an iterable object"""

    def __len__(self) -> int: ...

    def __bool__(self) -> bool:
        """Check whether the vector is nonempty"""

    def __repr__(self) -> str: ...

    def __iter__(self) -> Iterator[bool]: ...

    @overload
    def __getitem__(self, arg: int, /) -> bool: ...

    @overload
    def __getitem__(self, arg: slice, /) -> BoolVector: ...

    def clear(self) -> None:
        """Remove all items from list."""

    def append(self, arg: bool, /) -> None:
        """Append ``arg`` to the end of the list."""

    def insert(self, arg0: int, arg1: bool, /) -> None:
        """Insert object ``arg1`` before index ``arg0``."""

    def pop(self, index: int = -1) -> bool:
        """Remove and return item at ``index`` (default last)."""

    def extend(self, arg: BoolVector, /) -> None:
        """Extend ``self`` by appending elements from ``arg``."""

    @overload
    def __setitem__(self, arg0: int, arg1: bool, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: BoolVector, /) -> None: ...

    @overload
    def __delitem__(self, arg: int, /) -> None: ...

    @overload
    def __delitem__(self, arg: slice, /) -> None: ...

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    @overload
    def __contains__(self, arg: bool, /) -> bool: ...

    @overload
    def __contains__(self, arg: object, /) -> bool: ...

    def count(self, arg: bool, /) -> int:
        """Return number of occurrences of ``arg``."""

    def remove(self, arg: bool, /) -> None:
        """Remove first occurrence of ``arg``."""

class SubsetOptions:
    def __init__(self) -> None: ...

    @property
    def sanitize(self) -> bool:
        """Sanitize the resulting subset"""

    @sanitize.setter
    def sanitize(self, arg: bool, /) -> None: ...

    @property
    def clearComputedProps(self) -> bool:
        """clear all computed props on the subsetted molecule"""

    @clearComputedProps.setter
    def clearComputedProps(self, arg: bool, /) -> None: ...

    @property
    def copyAsQuery(self) -> bool:
        """Return the subset as a query"""

    @copyAsQuery.setter
    def copyAsQuery(self, arg: bool, /) -> None: ...

    @property
    def copyCoordinates(self) -> bool:
        """Copy the active coordinates from the molecule"""

    @copyCoordinates.setter
    def copyCoordinates(self, arg: bool, /) -> None: ...

    @property
    def conformerIdx(self) -> int:
        """What conformer idx to use for the coordinates default is -1"""

    @conformerIdx.setter
    def conformerIdx(self, arg: int, /) -> None: ...

    @property
    def method(self) -> SubsetMethod:
        """Subsetting method to use"""

    @method.setter
    def method(self, arg: SubsetMethod, /) -> None: ...

class UIntUIntMap:
    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, arg: UIntUIntMap) -> None:
        """Copy constructor"""

    @overload
    def __init__(self, arg: dict[int, int], /) -> None:
        """Construct from a dictionary"""

    def __len__(self) -> int: ...

    def __bool__(self) -> bool:
        """Check whether the map is nonempty"""

    def __repr__(self) -> str: ...

    @overload
    def __contains__(self, arg: int, /) -> bool: ...

    @overload
    def __contains__(self, arg: object, /) -> bool: ...

    def __iter__(self) -> Iterator[int]: ...

    def __getitem__(self, arg: int, /) -> int: ...

    def __delitem__(self, arg: int, /) -> None: ...

    def clear(self) -> None:
        """Remove all items"""

    def __setitem__(self, arg0: int, arg1: int, /) -> None: ...

    def update(self, arg: UIntUIntMap, /) -> None:
        """Update the map with element from ``arg``"""

    def __eq__(self, arg: object, /) -> bool: ...

    def __ne__(self, arg: object, /) -> bool: ...

    class ItemView:
        def __len__(self) -> int: ...

        def __iter__(self) -> Iterator[tuple[int, int]]: ...

    class KeyView:
        @overload
        def __contains__(self, arg: int, /) -> bool: ...

        @overload
        def __contains__(self, arg: object, /) -> bool: ...

        def __len__(self) -> int: ...

        def __iter__(self) -> Iterator[int]: ...

    class ValueView:
        def __len__(self) -> int: ...

        def __iter__(self) -> Iterator[int]: ...

    def keys(self) -> KeyView:
        """Returns an iterable view of the map's keys."""

    def values(self) -> ValueView:
        """Returns an iterable view of the map's values."""

    def items(self) -> ItemView:
        """Returns an iterable view of the map's items."""

class SubsetInfo:
    def __init__(self) -> None: ...

    @property
    def atomMapping(self) -> UIntUIntMap:
        """mapping from the original atom index to the subset atom index"""

    @atomMapping.setter
    def atomMapping(self, arg: UIntUIntMap, /) -> None: ...

    @property
    def bondMapping(self) -> UIntUIntMap:
        """mapping from the original bond index to the subset bond index"""

    @bondMapping.setter
    def bondMapping(self, arg: UIntUIntMap, /) -> None: ...

@overload
def CopyMolSubset(mol: rdkit.Chem.rdchem.Mol, path: Iterable[int], subsetInfo: SubsetInfo, options: SubsetOptions | None = None) -> rdkit.Chem.rdchem.Mol: ...

@overload
def CopyMolSubset(mol: rdkit.Chem.rdchem.Mol, path: Iterable[int], options: SubsetOptions | None = None) -> rdkit.Chem.rdchem.Mol:
    """copy a subset of a molecule"""

@overload
def CopyMolSubset(mol: rdkit.Chem.rdchem.Mol, atomIndices: Iterable[int], bondIndices: Iterable[int], subsetInfo: SubsetInfo, options: SubsetOptions | None = None) -> rdkit.Chem.rdchem.Mol: ...

@overload
def CopyMolSubset(mol: rdkit.Chem.rdchem.Mol, atomIndices: Iterable[int], bondIndices: Iterable[int], options: SubsetOptions | None = None) -> rdkit.Chem.rdchem.Mol:
    """
    Extract a subgraph from an ROMol. Bonds, atoms, substance groups and 
    stereo groups are only extracted to the subgraph if all participant entities 
    are contained within the given atoms and bonds. 

    ARGUMENTS:
     - mol - starting mol 
     - atoms - indices atoms to extract 
     - bonds - indices bonds to extract 
     - subsetInfo - optional subsetInfo to record the atoms and bonds used 
     - options - optional subset options, note the method is ignored since all the atoms and bonds are sp
    """

def FindPotentialStereo(mol: rdkit.Chem.rdchem.Mol, cleanIt: bool = False, flagPossible: bool = True) -> list[rdkit.Chem.rdchem.StereoInfo]:
    """
    find potential stereo elements in a molecule and returns them as StereoInfo objects
    Note that this function is still somewhat experimental and the API
    and results may change in a future release.
    """

def CleanupStereoGroups(mol: rdkit.Chem.rdchem.Mol) -> None:
    """removes atoms without specified chirality from stereo groups"""
