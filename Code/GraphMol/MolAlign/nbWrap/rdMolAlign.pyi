"""Module containing functions to align a molecule to a second molecule"""

from collections.abc import Sequence
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdchem
import rdkit.ForceField.rdForceField


class BestAlignmentParams:
    """Parameters controlling RMSD alignment"""

    def __init__(self) -> None: ...

    @property
    def maxMatches(self) -> int:
        """maximum number of substructure matches to consider"""

    @maxMatches.setter
    def maxMatches(self, arg: int, /) -> None: ...

    @property
    def symmetrizeConjugatedTerminalGroups(self) -> bool:
        """
        if true, conjugated terminal functional groups (like nitro or carboxylate)
        will be considered symmetrically.
        """

    @symmetrizeConjugatedTerminalGroups.setter
    def symmetrizeConjugatedTerminalGroups(self, arg: bool, /) -> None: ...

    @property
    def ignoreHs(self) -> bool:
        """if true, hydrogens will be ignored in the alignment"""

    @ignoreHs.setter
    def ignoreHs(self, arg: bool, /) -> None: ...

    @property
    def numThreads(self) -> int:
        """number of threads to use"""

    @numThreads.setter
    def numThreads(self, arg: int, /) -> None: ...

    @property
    def map(self) -> tuple:
        """the atom-atom mapping(s) used in the alignment"""

    @map.setter
    def map(self, arg: Sequence[Sequence[Sequence[int]]], /) -> None: ...

    @property
    def weights(self) -> tuple:
        """the weights used in the alignment"""

    @weights.setter
    def weights(self, arg: Sequence[float], /) -> None: ...

class O3A:
    """Open3DALIGN object"""

    def Align(self) -> float:
        """aligns probe molecule onto reference molecule"""

    def Trans(self) -> tuple[float, Annotated[NDArray[numpy.float64], dict(shape=(None, None))]]:
        """
        returns the transformation which aligns probe molecule onto reference molecule
        """

    def Score(self) -> float:
        """returns the O3AScore of the alignment"""

    def Matches(self) -> list[list[int]]:
        """returns the AtomMap as found by Open3DALIGN"""

    def Weights(self) -> list[float]:
        """returns the weight vector as found by Open3DALIGN"""

def GetAlignmentTransform(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, prbCid: int = -1, refCid: int = -1, atomMap: Sequence[Sequence[int]] | None = None, weights: Sequence[float] | None = None, reflect: bool = False, maxIters: int = 50) -> tuple[float, Annotated[NDArray[numpy.float64], dict(shape=(None, None))]]:
    """
    Compute the transformation required to align a molecule

    The 3D transformation required to align the specied conformation in the probe molecule
    to a specified conformation in the reference molecule is computed so that the root mean
    squared distance between a specified set of atoms is minimized

    ARGUMENTS
     - prbMol    molecule that is to be aligned
     - refMol    molecule used as the reference for the alignment
     - prbCid    ID of the conformation in the probe to be used
                      for the alignment (defaults to first conformation)
     - refCid    ID of the conformation in the ref molecule to which
                      the alignment is computed (defaults to first conformation)
     - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                      used to compute the alignments. If this mapping is
                      not specified an attempt is made to generate one by
                      substructure matching
     - weights   Optionally specify weights for each of the atom pairs
     - reflect   if true reflect the conformation of the probe molecule
     - maxIters  maximum number of iterations used in minimizing the RMSD

    RETURNS
    a tuple of (RMSD value, transform matrix)
    """

@overload
def GetBestAlignmentTransform(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, prbCid: int = -1, refCid: int = -1, map: Sequence[Sequence[Sequence[int]]] | None = None, maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: Sequence[float] | None = None, reflect: bool = False, maxIters: int = 50, numThreads: int = 1) -> tuple[float, Annotated[NDArray[numpy.float64], dict(shape=(None, None))], list[tuple[int, int]]]:
    """
    Compute the optimal RMS, transformation and atom map for aligning
    two molecules, taking symmetry into account. Molecule coordinates
    are left unaltered.

    This function will attempt to align all permutations of matching atom
    orders in both molecules, for some molecules it will lead to 'combinatorial
    explosion' especially if hydrogens are present.
    Use 'GetAlignmentTransform' to align molecules without changing the atom order.

    ARGUMENTS
     - prbMol      molecule that is to be aligned
     - refMol      molecule used as the reference for the alignment
     - prbCid      ID of the conformation in the probe to be used
                   for the alignment (defaults to first conformation)
     - refCid      ID of the conformation in the ref molecule to which
                   the alignment is computed (defaults to first conformation)
     - map:        (optional) a list of lists of (probeAtomId, refAtomId)
                   tuples with the atom-atom mappings of the two
                   molecules. If not provided, these will be generated
                   using a substructure search.
     - maxMatches  (optional) if atomMap is empty, this will be the max number of
                   matches found in a SubstructMatch().
     - symmetrizeConjugatedTerminalGroups (optional) if set, conjugated
                   terminal functional groups (like nitro or carboxylate)
                   will be considered symmetrically.
     - weights     Optionally specify weights for each of the atom pairs
     - reflect     if true reflect the conformation of the probe molecule
     - maxIters    maximum number of iterations used in minimizing the RMSD
     - numThreads  (optional) number of threads to use

    RETURNS
    a tuple of (RMSD value, best transform matrix, best atom map)
    """

@overload
def GetBestAlignmentTransform(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, params: BestAlignmentParams, prbCid: int = -1, refCid: int = -1, reflect: bool = False, maxIters: int = 50) -> tuple[float, Annotated[NDArray[numpy.float64], dict(shape=(None, None))], list[tuple[int, int]]]:
    """
    Compute the optimal RMS, transformation and atom map for aligning
    two molecules using a BestAlignmentParams object.

    RETURNS
    a tuple of (RMSD value, best transform matrix, best atom map)
    """

def AlignMol(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, prbCid: int = -1, refCid: int = -1, atomMap: Sequence[Sequence[int]] | None = None, weights: Sequence[float] | None = None, reflect: bool = False, maxIters: int = 50) -> float:
    """
    Optimally (minimum RMSD) align a molecule to another molecule

    The 3D transformation required to align the specied conformation in the probe molecule
    to a specified conformation in the reference molecule is computed so that the root mean
    squared distance between a specified set of atoms is minimized.
    This transform is then applied to the specified conformation in the probe molecule

    ARGUMENTS
     - prbMol    molecule that is to be aligned
     - refMol    molecule used as the reference for the alignment
     - prbCid    ID of the conformation in the probe to be used
                      for the alignment (defaults to first conformation)
     - refCid    ID of the conformation in the ref molecule to which
                      the alignment is computed (defaults to first conformation)
     - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                      used to compute the alignments. If this mapping is
                      not specified an attempt is made to generate one by
                      substructure matching
     - weights   Optionally specify weights for each of the atom pairs
     - reflect   if true reflect the conformation of the probe molecule
     - maxIters  maximum number of iterations used in minimizing the RMSD

    RETURNS
    RMSD value
    """

@overload
def GetBestRMS(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, prbId: int = -1, refId: int = -1, map: Sequence[Sequence[Sequence[int]]] | None = None, maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: Sequence[float] | None = None, numThreads: int = 1) -> float:
    """
    Returns the optimal RMS for aligning two molecules, taking
    symmetry into account. As a side-effect, the probe molecule is
    left in the aligned state.

    Note:
    This function will attempt to align all permutations of matching atom
    orders in both molecules, for some molecules it will lead to
    'combinatorial explosion' especially if hydrogens are present.
    Use 'rdkit.Chem.AllChem.AlignMol' to align molecules without changing
    the atom order.

    ARGUMENTS
     - prbMol:      the molecule to be aligned to the reference
     - refMol:      the reference molecule
     - prbId:       (optional) probe conformation to use
     - refId:       (optional) reference conformation to use
     - map:         (optional) a list of lists of (probeAtomId,refAtomId)
                   tuples with the atom-atom mappings of the two
                   molecules. If not provided, these will be generated
                   using a substructure search.
     - maxMatches:  (optional) if map isn't specified, this will be
                   the max number of matches found in a SubstructMatch()
     - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated
                   terminal functional groups (like nitro or carboxylate)
                   will be considered symmetrically
     - weights:     (optional) weights for mapping
     - numThreads:  (optional) number of threads to use

    RETURNS
    The best RMSD found
    """

@overload
def GetBestRMS(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, params: BestAlignmentParams, prbId: int = -1, refId: int = -1) -> float:
    """
    Returns the optimal RMS for aligning two molecules using a BestAlignmentParams object.

    RETURNS
    The best RMSD found
    """

@overload
def GetAllConformerBestRMS(mol: rdkit.Chem.rdchem.Mol, numThreads: int = 1, map: Sequence[Sequence[Sequence[int]]] | None = None, maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: Sequence[float] | None = None) -> tuple[float, ...]:
    """
    Returns the symmetric distance matrix between the conformers of a molecule.
    getBestRMS() is used to calculate the inter-conformer distances

    ARGUMENTS
     - mol:       the molecule to be considered
     - numThreads:  (optional) number of threads to use
     - map:         (optional) a list of lists of (probeAtomId,refAtomId)
                   tuples with the atom-atom mappings of the two
                   molecules. If not provided, these will be generated
                   using a substructure search.
     - maxMatches:  (optional) if map isn't specified, this will be
                   the max number of matches found in a SubstructMatch()
     - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated
                   terminal functional groups (like nitro or carboxylate)
                   will be considered symmetrically
     - weights:     (optional) weights for mapping

    RETURNS
    A tuple with the best RMSDS. The ordering is [(1,0),(2,0),(2,1),(3,0),... etc]
    """

@overload
def GetAllConformerBestRMS(mol: rdkit.Chem.rdchem.Mol, params: BestAlignmentParams) -> tuple[float, ...]:
    """
    Returns the symmetric distance matrix between the conformers of a molecule
    using a BestAlignmentParams object.

    RETURNS
    A tuple with the best RMSDS. The ordering is [(1,0),(2,0),(2,1),(3,0),... etc]
    """

def CalcRMS(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, prbId: int = -1, refId: int = -1, map: Sequence[Sequence[Sequence[int]]] | None = None, maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: Sequence[float] | None = None) -> float:
    """
    Returns the RMS between two molecules, taking symmetry into account.
    In contrast to getBestRMS, the RMS is computed 'in place', i.e.
    probe molecules are not aligned to the reference ahead of the
    RMS calculation. This is useful, for example, to compute
    the RMSD between docking poses and the co-crystallized ligand.

    Note:
    This function will attempt to match all permutations of matching atom
    orders in both molecules, for some molecules it will lead to
    'combinatorial explosion' especially if hydrogens are present.

    ARGUMENTS
     - prbMol:      the molecule to be aligned to the reference
     - refMol:      the reference molecule
     - prbCId:      (optional) probe conformation to use
     - refCId:      (optional) reference conformation to use
     - map:         (optional) a list of lists of (probeAtomId, refAtomId)
                   tuples with the atom-atom mappings of the two
                   molecules. If not provided, these will be generated
                   using a substructure search.
     - maxMatches:  (optional) if map isn't specified, this will be
                   the max number of matches found in a SubstructMatch()
     - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated
                   terminal functional groups (like nitro or carboxylate)
                   will be considered symmetrically
     - weights:     (optional) weights for mapping

    RETURNS
    The best RMSD found
    """

def AlignMolConformers(mol: rdkit.Chem.rdchem.Mol, atomIds: Sequence[int] | None = None, confIds: Sequence[int] | None = None, weights: Sequence[float] | None = None, reflect: bool = False, maxIters: int = 50, RMSlist: object | None = None) -> None:
    """
    Align conformations in a molecule to each other

    The first conformation in the molecule is used as the reference

    ARGUMENTS
     - mol          molecule of interest
     - atomIds      List of atom ids to use a points for alignment - defaults to all atoms
     - confIds      Ids of conformations to align - defaults to all conformers
     - weights      Optionally specify weights for each of the atom pairs
     - reflect      if true reflect the conformation of the probe molecule
     - maxIters     maximum number of iterations used in minimizing the RMSD
     - RMSlist      if provided, fills in the RMS values between the reference
                    conformation and the other aligned conformations
    """

def RandomTransform(mol: rdkit.Chem.rdchem.Mol, cid: int = -1, seed: int = -1) -> None:
    """
    Perform a random transformation on a molecule

    ARGUMENTS
     - mol    molecule that is to be transformed
     - cid    ID of the conformation in the mol to be transformed
              (defaults to first conformation)
     - seed   seed used to initialize the random generator
              (defaults to -1, that is no seeding)
    """

def GetO3A(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, prbMMFFMolProperties: rdkit.ForceField.rdForceField.MMFFMolProperties | None = None, refMMFFMolProperties: rdkit.ForceField.rdForceField.MMFFMolProperties | None = None, prbCid: int = -1, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: Sequence[Sequence[int]] | None = None, constraintWeights: Sequence[float] | None = None) -> O3A:
    """
    Get an O3A object with atomMap and weights vectors to overlay
    the probe molecule onto the reference molecule based on
    MMFF atom types and charges

    ARGUMENTS
     - prbMol                   molecule that is to be aligned
     - refMol                   molecule used as the reference for the alignment
     - prbMMFFMolProperties   MMFFMolProperties object for the probe molecule as returned
                                by MMFFGetMoleculeProperties()
     - refMMFFMolProperties   MMFFMolProperties object for the reference molecule as returned
                                by MMFFGetMoleculeProperties()
     - prbCid                   ID of the conformation in the probe to be used
                                for the alignment (defaults to first conformation)
     - refCid                   ID of the conformation in the ref molecule to which
                                the alignment is computed (defaults to first conformation)
     - reflect                  if true reflect the conformation of the probe molecule
                                (defaults to false)
     - maxIters                 maximum number of iterations used in minimizing the RMSD
                                (defaults to 50)
     - options                  least 2 significant bits encode accuracy
                                (0: maximum, 3: minimum; defaults to 0)
                                bit 3 triggers local optimization of the alignment
                                (no computation of the cost matrix; defaults: off)
     - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                                which shall be used for the alignment (defaults to [])
     - constraintWeights        optionally specify weights for each of the constraints
                                (weights default to 100.0)

    RETURNS
    The O3A object
    """

def GetCrippenO3A(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, prbCrippenContribs: Sequence[Sequence[float]] | None = None, refCrippenContribs: Sequence[Sequence[float]] | None = None, prbCid: int = -1, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: Sequence[Sequence[int]] | None = None, constraintWeights: Sequence[float] | None = None) -> O3A:
    """
    Get an O3A object with atomMap and weights vectors to overlay
    the probe molecule onto the reference molecule based on
    Crippen logP atom contributions

    ARGUMENTS
     - prbMol                   molecule that is to be aligned
     - refMol                   molecule used as the reference for the alignment
     - prbCrippenContribs       Crippen atom contributions for the probe molecule
                                as a list of (logp, mr) tuples, as returned
                                by _CalcCrippenContribs()
     - refCrippenContribs       Crippen atom contributions for the reference molecule
                                as a list of (logp, mr) tuples, as returned
                                by _CalcCrippenContribs()
     - prbCid                   ID of the conformation in the probe to be used
                                for the alignment (defaults to first conformation)
     - refCid                   ID of the conformation in the ref molecule to which
                                the alignment is computed (defaults to first conformation)
     - reflect                  if true reflect the conformation of the probe molecule
                                (defaults to false)
     - maxIters                 maximum number of iterations used in minimizing the RMSD
                                (defaults to 50)
     - options                  least 2 significant bits encode accuracy
                                (0: maximum, 3: minimum; defaults to 0)
                                bit 3 triggers local optimization of the alignment
                                (no computation of the cost matrix; defaults: off)
     - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                                which shall be used for the alignment (defaults to [])
     - constraintWeights        optionally specify weights for each of the constraints
                                (weights default to 100.0)

    RETURNS
    The O3A object
    """

def GetO3AForProbeConfs(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, numThreads: int = 1, prbMMFFMolProperties: rdkit.ForceField.rdForceField.MMFFMolProperties | None = None, refMMFFMolProperties: rdkit.ForceField.rdForceField.MMFFMolProperties | None = None, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: Sequence[Sequence[int]] | None = None, constraintWeights: Sequence[float] | None = None) -> tuple[O3A, ...]:
    """
    Get a vector of O3A objects for the overlay of all
    the probe molecule's conformations onto the reference molecule based on
    MMFF atom types and charges

    ARGUMENTS
     - prbMol                   molecule that is to be aligned
     - refMol                   molecule used as the reference for the alignment
     - numThreads :             the number of threads to use, only has an effect if
                                the RDKit was built with thread support (defaults to 1)
                                If set to zero, the max supported by the system will be used.
     - prbMMFFMolProperties   MMFFMolProperties object for the probe molecule as returned
                                by MMFFGetMoleculeProperties()
     - refMMFFMolProperties   MMFFMolProperties object for the reference molecule as returned
                                by MMFFGetMoleculeProperties()
     - refCid                   ID of the conformation in the ref molecule to which
                                the alignment is computed (defaults to first conformation)
     - reflect                  if true reflect the conformation of the probe molecule
                                (defaults to false)
     - maxIters                 maximum number of iterations used in minimizing the RMSD
                                (defaults to 50)
     - options                  least 2 significant bits encode accuracy
                                (0: maximum, 3: minimum; defaults to 0)
                                bit 3 triggers local optimization of the alignment
                                (no computation of the cost matrix; defaults: off)
     - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                                which shall be used for the alignment (defaults to [])
     - constraintWeights        optionally specify weights for each of the constraints
                                (weights default to 100.0)

    RETURNS
    A vector of O3A objects
    """

def GetCrippenO3AForProbeConfs(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, numThreads: int = 1, prbCrippenContribs: Sequence[Sequence[float]] | None = None, refCrippenContribs: Sequence[Sequence[float]] | None = None, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: Sequence[Sequence[int]] | None = None, constraintWeights: Sequence[float] | None = None) -> tuple[O3A, ...]:
    """
    Get a vector of O3A objects for the overlay of all
    the probe molecule's conformations onto the reference molecule based on
    Crippen logP atom contributions

    ARGUMENTS
     - prbMol                   molecule that is to be aligned
     - refMol                   molecule used as the reference for the alignment
     - numThreads :             the number of threads to use, only has an effect if
                                the RDKit was built with thread support (defaults to 1)
     - prbCrippenContribs       Crippen atom contributions for the probe molecule
                                as a list of (logp, mr) tuples, as returned
                                by _CalcCrippenContribs()
     - refCrippenContribs       Crippen atom contributions for the reference molecule
                                as a list of (logp, mr) tuples, as returned
                                by _CalcCrippenContribs()
     - refCid                   ID of the conformation in the ref molecule to which
                                the alignment is computed (defaults to first conformation)
     - reflect                  if true reflect the conformation of the probe molecule
                                (defaults to false)
     - maxIters                 maximum number of iterations used in minimizing the RMSD
                                (defaults to 50)
     - options                  least 2 significant bits encode accuracy
                                (0: maximum, 3: minimum; defaults to 0)
                                bit 3 triggers local optimization of the alignment
                                (no computation of the cost matrix; defaults: off)
     - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                                which shall be used for the alignment (defaults to [])
     - constraintWeights        optionally specify weights for each of the constraints
                                (weights default to 100.0)

    RETURNS
    A vector of O3A objects
    """

def GetAllConformerBestRMSToRef(prbMol: rdkit.Chem.rdchem.Mol, refMol: rdkit.Chem.rdchem.Mol, params: BestAlignmentParams | None = None) -> tuple[float, ...]:
    """
    Get the RMSD matrix between all conformers of refMol\\n\\
    and all the conformers of prbMol.
    getBestRMS() is used to calculate the inter-conformer distances
      This function will attempt to align all permutations of matching atom
      orders in both molecules, for some molecules it will lead to 'combinatorial' 
      explosion' especially if hydrogens are present.

    ARGUMENTS
     - prbMol        the probe molecule\\n\\
     - refMol        the reference molecule\\n\\
     - params     parameters for the matching\\n\\

    RETURNS
     Vector with the RMSD values stored in the order:
     [(0, 0), (0, 1), (0, 2), (1, 0), (2, 1), ...]
     where the first idx is a conformerID of the refMol and the second is the
     confid of the prbMol
    """
