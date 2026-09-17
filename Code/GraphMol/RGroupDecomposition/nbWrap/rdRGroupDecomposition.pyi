"""TEST!!! Module containing RGroupDecomposition classes and functions."""

from collections.abc import Iterable
import enum
from typing import overload

import rdkit.Chem.rdchem


class RGroupLabels(enum.IntEnum):
    IsotopeLabels = 1

    AtomMapLabels = 2

    AtomIndexLabels = 4

    RelabelDuplicateLabels = 8

    MDLRGroupLabels = 16

    DummyAtomLabels = 32

    AutoDetect = 255

class RGroupMatching(enum.IntEnum):
    Greedy = 1

    GreedyChunks = 2

    Exhaustive = 4

    NoSymmetrization = 8

    GA = 16

class RGroupLabelling(enum.IntEnum):
    AtomMap = 1

    Isotope = 2

    MDLRGroup = 4

class RGroupCoreAlignment(enum.IntEnum):
    NoAlignment = 0

    MCS = 1

class RGroupScore(enum.IntEnum):
    Match = 1

    FingerprintVariance = 4

class RGroupDecompositionParameters:
    """
    RGroupDecompositionParameters controls how the RGroupDecomposition
    sets labelling and matches structures
    OPTIONS:
      - RGroupCoreAlignment: can be one of RGroupCoreAlignment.None_ or
    RGroupCoreAlignment.MCS
                             If set to MCS, cores labels are mapped to
    each other using their
                             Maximum common substructure overlap.
      - RGroupLabels: optionally set where the rgroup labels to use are
    encoded.
                       RGroupLabels.IsotopeLabels - labels are stored
    on isotopes
                       RGroupLabels.AtomMapLabels - labels are stored
    on atommaps
                       RGroupLabels.MDLRGroupLabels - labels are stored
    on MDL R-groups
                       RGroupLabels.DummyAtomLabels - labels are stored
    on dummy atoms
                       RGroupLabels.AtomIndexLabels - use the atom index
    as the label
                       RGroupLabels.RelabelDuplicateLabels - fix any
    duplicate labels
                       RGroupLabels.AutoDetect - auto detect the label
    [default]
         Note: in all cases, any rgroups found on unlabelled atoms will
    be automatically
                labelled.
      - RGroupLabelling: choose where the rlabels are stored on the
    decomposition
                          RGroupLabelling.AtomMap - store rgroups as atom
    maps (for smiles)
                          RGroupLabelling.Isotope - store rgroups on the
    isotope
                          RGroupLabelling.MDLRGroup - store rgroups as mdl
    rgroups (for molblocks)
                         default: AtomMap | MDLRGroup
      - onlyMatchAtRGroups: only allow rgroup decomposition at the
    specified rgroups
      - removeAllHydrogenRGroups: remove all user-defined rgroups that
    only have hydrogens
      - removeAllHydrogenRGroupsAndLabels: remove all user-defined
    rgroups that only have hydrogens, and also remove the corresponding
    labels from the core
      - removeHydrogensPostMatch: remove all hydrogens from the output
    molecules
      - allowNonTerminalRGroups: allow labelled Rgroups of degree 2 or
    more
      - doTautomers: match all tautomers of a core against each
    input structure
      - doEnumeration: expand input cores into enumerated mol bundles
      - allowMultipleRGroupsOnUnlabelled: permit more than one rgroup to
    be attached to an unlabelled core atom
      - allowMultipleCoresInSameMol: permit a core to match more than
    once in the same molecule if the sets of matched atoms are not equal
    (default=False)
    """

    def __init__(self) -> None:
        """Constructor, takes no arguments"""

    @property
    def labels(self) -> int: ...

    @labels.setter
    def labels(self, arg: int, /) -> None: ...

    @property
    def matchingStrategy(self) -> int: ...

    @matchingStrategy.setter
    def matchingStrategy(self, arg: int, /) -> None: ...

    @property
    def scoreMethod(self) -> int: ...

    @scoreMethod.setter
    def scoreMethod(self, arg: int, /) -> None: ...

    @property
    def rgroupLabelling(self) -> int: ...

    @rgroupLabelling.setter
    def rgroupLabelling(self, arg: int, /) -> None: ...

    @property
    def alignment(self) -> int: ...

    @alignment.setter
    def alignment(self, arg: int, /) -> None: ...

    @property
    def chunkSize(self) -> int: ...

    @chunkSize.setter
    def chunkSize(self, arg: int, /) -> None: ...

    @property
    def onlyMatchAtRGroups(self) -> bool: ...

    @onlyMatchAtRGroups.setter
    def onlyMatchAtRGroups(self, arg: bool, /) -> None: ...

    @property
    def removeAllHydrogenRGroups(self) -> bool: ...

    @removeAllHydrogenRGroups.setter
    def removeAllHydrogenRGroups(self, arg: bool, /) -> None: ...

    @property
    def removeHydrogensPostMatch(self) -> bool: ...

    @removeHydrogensPostMatch.setter
    def removeHydrogensPostMatch(self, arg: bool, /) -> None: ...

    @property
    def timeout(self) -> float: ...

    @timeout.setter
    def timeout(self, arg: float, /) -> None: ...

    @property
    def gaPopulationSize(self) -> int: ...

    @gaPopulationSize.setter
    def gaPopulationSize(self, arg: int, /) -> None: ...

    @property
    def gaMaximumOperations(self) -> int: ...

    @gaMaximumOperations.setter
    def gaMaximumOperations(self, arg: int, /) -> None: ...

    @property
    def gaNumberOperationsWithoutImprovement(self) -> int: ...

    @gaNumberOperationsWithoutImprovement.setter
    def gaNumberOperationsWithoutImprovement(self, arg: int, /) -> None: ...

    @property
    def gaRandomSeed(self) -> int: ...

    @gaRandomSeed.setter
    def gaRandomSeed(self, arg: int, /) -> None: ...

    @property
    def gaNumberRuns(self) -> int: ...

    @gaNumberRuns.setter
    def gaNumberRuns(self, arg: int, /) -> None: ...

    @property
    def gaParallelRuns(self) -> bool: ...

    @gaParallelRuns.setter
    def gaParallelRuns(self, arg: bool, /) -> None: ...

    @property
    def allowNonTerminalRGroups(self) -> bool: ...

    @allowNonTerminalRGroups.setter
    def allowNonTerminalRGroups(self, arg: bool, /) -> None: ...

    @property
    def removeAllHydrogenRGroupsAndLabels(self) -> bool: ...

    @removeAllHydrogenRGroupsAndLabels.setter
    def removeAllHydrogenRGroupsAndLabels(self, arg: bool, /) -> None: ...

    @property
    def allowMultipleRGroupsOnUnlabelled(self) -> bool: ...

    @allowMultipleRGroupsOnUnlabelled.setter
    def allowMultipleRGroupsOnUnlabelled(self, arg: bool, /) -> None: ...

    @property
    def allowMultipleCoresInSameMol(self) -> bool: ...

    @allowMultipleCoresInSameMol.setter
    def allowMultipleCoresInSameMol(self, arg: bool, /) -> None: ...

    @property
    def doTautomers(self) -> bool: ...

    @doTautomers.setter
    def doTautomers(self, arg: bool, /) -> None: ...

    @property
    def doEnumeration(self) -> bool: ...

    @doEnumeration.setter
    def doEnumeration(self, arg: bool, /) -> None: ...

    @property
    def substructMatchParams(self) -> rdkit.Chem.rdchem.SubstructMatchParameters: ...

    @property
    def includeTargetMolInResults(self) -> bool: ...

    @includeTargetMolInResults.setter
    def includeTargetMolInResults(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class RGroupDecomposition:
    @overload
    def __init__(self, cores: object) -> None:
        """Construct from a molecule or sequence of molecules"""

    @overload
    def __init__(self, cores: object, params: RGroupDecompositionParameters) -> None:
        """
        Construct from a molecule or sequence of molecules and a parameters object
        """

    def Add(self, mol: rdkit.Chem.rdchem.Mol) -> int: ...

    def GetMatchingCoreIdx(self, mol: rdkit.Chem.rdchem.Mol, matches: list[tuple[tuple[int, int], ...]] | None = None) -> int: ...

    def Process(self) -> bool:
        """
        Process the rgroups (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)
        """

    def ProcessAndScore(self) -> tuple:
        """
        Process the rgroups and returns the score (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)
        """

    def GetRGroupLabels(self) -> list:
        """
        Return the current list of found rgroups.
        Note, Process() should be called first
        """

    def GetRGroupsAsRows(self, asSmiles: bool = False) -> list:
        """
        Return the rgroups as rows (note: can be fed directly into a
        pandas datatable)
        ARGUMENTS:
         - asSmiles: if True return smiles strings, otherwise return
        molecules [default: False]
          Row structure:
             rows[idx] = {rgroup_label: molecule_or_smiles}
        """

    def GetRGroupsAsColumns(self, asSmiles: bool = False) -> dict:
        """
        Return the rgroups as columns (note: can be fed directly into a
        pandas datatable)
        ARGUMENTS:
         - asSmiles: if True return smiles strings, otherwise return
        molecules [default: False]
          Column structure:
             columns[rgroup_label] = [ mols_or_smiles ]
        """

def RGroupDecompose(cores: object, mols: Iterable[rdkit.Chem.rdchem.Mol], asSmiles: bool = False, asRows: bool = True, options: RGroupDecompositionParameters = ...) -> object:
    """
    Decompose a collection of molecules into their Rgroups
    ARGUMENTS:
      - cores: a set of cores from most to least specific.
               See RGroupDecompositionParameters for more details
               on how the cores can be labelled
      - mols: the molecules to be decomposed
      - asSmiles: if True return smiles strings, otherwise return
    molecules [default: False]
      - asRows: return the results as rows (default) otherwise return
    columns
      - options: RGroupDecompositionParameters object that defines
    the parameters for the decomposition.
               See RGroupDecompositionParameters for defaults

    RETURNS: row_or_column_results, unmatched

      Row structure:
         rows[idx] = {rgroup_label: molecule_or_smiles}
      Column structure:
         columns[rgroup_label] = [ mols_or_smiles ]

      unmatched is a vector of indices in the input mols that were not
    matched.
    """

def RelabelMappedDummies(mol: rdkit.Chem.rdchem.Mol, inputLabels: int = 7, outputLabels: int = 4) -> None:
    """
    Relabel dummy atoms bearing an R-group mapping (as
    atom map number, isotope or MDLRGroup label) such that
    they will be displayed by the rendering code as R# rather
    than #*, *:#, #*:#, etc. By default, only the MDLRGroup label
    is retained on output; this may be configured through the
    outputLabels parameter.
    In case there are multiple potential R-group mappings,
    the priority on input is Atom map number > Isotope > MDLRGroup.
    The inputLabels parameter allows to configure which mappings
    are taken into consideration.
    """
