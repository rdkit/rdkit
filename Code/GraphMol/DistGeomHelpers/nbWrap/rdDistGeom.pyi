"""
Module containing functions to compute atomic coordinates in 3D using
distance geometry
"""

from collections.abc import Mapping
import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdchem
import rdkit.Geometry.rdGeometry


@overload
def GetExperimentalTorsions(mol: rdkit.Chem.rdchem.Mol, useExpTorsionAnglePrefs: bool = True, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, useBasicKnowledge: bool = True, ETversion: int = 2, printExpTorsionAngles: bool = False) -> list: ...

@overload
def GetExperimentalTorsions(mol: rdkit.Chem.rdchem.Mol, embedParams: EmbedParameters) -> list:
    """
    returns information about the bonds corresponding to experimental torsions
    """

@overload
def EmbedMolecule(mol: rdkit.Chem.rdchem.Mol, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, coordMap: Mapping[int, rdkit.Geometry.rdGeometry.Point3D] = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, ETversion: int = 2, useMacrocycle14config: bool = True) -> int:
    """
    Use distance geometry to obtain initial
    coordinates for a molecule

    ARGUMENTS:

       - mol : the molecule of interest
       - maxAttempts : maximum number of embedding attempts to use for a single conformation
       - randomSeed : provide a seed for the random number generator
                      so that the same coordinates can be obtained
                      for a molecule on multiple runs. If -1, the
                      RNG will not be seeded.
       - clearConfs : clear all existing conformations on the molecule
       - useRandomCoords : Start the embedding from random coordinates instead of
                           using eigenvalues of the distance matrix.
       - boxSizeMult :  Determines the size of the box that is used for
                        random coordinates. If this is a positive number, the
                        side length will equal the largest element of the distance
                        matrix times boxSizeMult. If this is a negative number,
                        the side length will equal -boxSizeMult (i.e. independent
                        of the elements of the distance matrix).
       - randNegEig : If the embedding yields a negative eigenvalue,
                      pick coordinates that correspond
                      to this component at random
       - numZeroFail : fail embedding if we have at least this many zero eigenvalues
       - coordMap : a dictionary mapping atom IDs->coordinates. Use this to
                    require some atoms to have fixed coordinates in the resulting
                    conformation.
       - forceTol : tolerance to be used during the force-field minimization with
                    the distance geometry force field.
       - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
                    of the bounds matrix fails.
       - enforceChirality : enforce the correct chirality if chiral centers are present.
       - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
       - useBasicKnowledge : impose basic knowledge such as flat rings
       - printExpTorsionAngles : print the output from the experimental torsion angles
       - useMacrocycleTorsions : use additional torsion profiles for macrocycles
       - ETversion : version of the standard torsion definitions to use. NOTE for both
                     ETKDGv2 and ETKDGv3 this should be 2 since ETKDGv3 uses the ETKDGv2
                     definitions for standard torsions
       - useMacrocycle14config : This forces amides and esters to be trans in macrocycles.
                                 This does not affect chain amides / esters!\\n\\

    RETURNS:

       ID of the new conformation added to the molecule or -1 if the embedding fails.
    """

@overload
def EmbedMolecule(mol: rdkit.Chem.rdchem.Mol, params: EmbedParameters) -> int:
    """
    Use distance geometry to obtain initial
    coordinates for a molecule

    ARGUMENTS:

       - mol : the molecule of interest
       - params : an EmbedParameters object

    RETURNS:

       ID of the new conformation added to the molecule or -1 if the embedding fails.
    """

@overload
def EmbedMultipleConfs(mol: rdkit.Chem.rdchem.Mol, numConfs: int = 10, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, pruneRmsThresh: float = -1.0, coordMap: Mapping[int, rdkit.Geometry.rdGeometry.Point3D] = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, numThreads: int = 1, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, ETversion: int = 2, useMacrocycle14config: bool = True) -> list[int]:
    """
    Use distance geometry to obtain multiple sets of
    coordinates for a molecule

    ARGUMENTS:

      - mol : the molecule of interest
      - numConfs : the number of conformers to generate
      - maxAttempts : maximum number of embedding attempts to use for a single conformation
      - randomSeed : provide a seed for the random number generator
                     so that the same coordinates can be obtained
                     for a molecule on multiple runs. If -1, the
                     RNG will not be seeded.
      - clearConfs : clear all existing conformations on the molecule
      - useRandomCoords : Start the embedding from random coordinates instead of
                          using eigenvalues of the distance matrix.
      - boxSizeMult    Determines the size of the box that is used for
                       random coordinates. If this is a positive number, the
                       side length will equal the largest element of the distance
                       matrix times boxSizeMult. If this is a negative number,
                       the side length will equal -boxSizeMult (i.e. independent
                       of the elements of the distance matrix).
      - randNegEig : If the embedding yields a negative eigenvalue,
                     pick coordinates that correspond
                     to this component at random
      - numZeroFail : fail embedding if we have at least this many zero eigenvalues
      - pruneRmsThresh : Retain only the conformations out of 'numConfs'
                        after embedding that are at least
                        this far apart from each other.
                        RMSD is computed on the heavy atoms.
                        Pruning is greedy; i.e. the first embedded conformation
                        is retained and from then on only those that are at
                        least pruneRmsThresh away from all retained conformations
                        are kept. The pruning is done after embedding and
                        bounds violation minimization. No pruning by default.
      - coordMap : a dictionary mapping atom IDs->coordinates. Use this to
                   require some atoms to have fixed coordinates in the resulting
                   conformation.
      - forceTol : tolerance to be used during the force-field minimization with
                   the distance geometry force field.
      - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
                   of the bounds matrix fails.
      - enforceChirality : enforce the correct chirality if chiral centers are present.
      - numThreads : number of threads to use while embedding. This only has an effect if the RDKit
                   was built with multi-thread support.
                  If set to zero, the max supported by the system will be used.
      - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
      - useBasicKnowledge : impose basic knowledge such as flat rings
      - printExpTorsionAngles : print the output from the experimental torsion angles

    RETURNS:

       Iterator which yields new conformation IDs
    """

@overload
def EmbedMultipleConfs(mol: rdkit.Chem.rdchem.Mol, numConfs: int, params: EmbedParameters) -> list[int]:
    """
    Use distance geometry to obtain multiple sets of
    coordinates for a molecule

    ARGUMENTS:

      - mol : the molecule of interest
      - numConfs : the number of conformers to generate
      - params : an EmbedParameters object

    RETURNS:

       Iterator which yields new conformation IDs
    """

class EmbedFailureCauses(enum.IntEnum):
    INITIAL_COORDS = 0

    FIRST_MINIMIZATION = 1

    CHECK_TETRAHEDRAL_CENTERS = 2

    CHECK_CHIRAL_CENTERS = 3

    MINIMIZE_FOURTH_DIMENSION = 4

    ETK_MINIMIZATION = 5

    FINAL_CHIRAL_BOUNDS = 6

    FINAL_CENTER_IN_VOLUME = 7

    LINEAR_DOUBLE_BOND = 8

    BAD_DOUBLE_BOND_STEREO = 9

    CHECK_CHIRAL_CENTERS2 = 10

    EXCEEDED_TIMEOUT = 11

    MINIMIZATION = 12

    KTERM_VIOLATION = 13

    CLASH = 14

INITIAL_COORDS: EmbedFailureCauses = EmbedFailureCauses.INITIAL_COORDS

FIRST_MINIMIZATION: EmbedFailureCauses = EmbedFailureCauses.FIRST_MINIMIZATION

CHECK_TETRAHEDRAL_CENTERS: EmbedFailureCauses = EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS

CHECK_CHIRAL_CENTERS: EmbedFailureCauses = EmbedFailureCauses.CHECK_CHIRAL_CENTERS

MINIMIZE_FOURTH_DIMENSION: EmbedFailureCauses = EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION

ETK_MINIMIZATION: EmbedFailureCauses = EmbedFailureCauses.ETK_MINIMIZATION

FINAL_CHIRAL_BOUNDS: EmbedFailureCauses = EmbedFailureCauses.FINAL_CHIRAL_BOUNDS

FINAL_CENTER_IN_VOLUME: EmbedFailureCauses = EmbedFailureCauses.FINAL_CENTER_IN_VOLUME

LINEAR_DOUBLE_BOND: EmbedFailureCauses = EmbedFailureCauses.LINEAR_DOUBLE_BOND

BAD_DOUBLE_BOND_STEREO: EmbedFailureCauses = EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO

CHECK_CHIRAL_CENTERS2: EmbedFailureCauses = EmbedFailureCauses.CHECK_CHIRAL_CENTERS2

EXCEEDED_TIMEOUT: EmbedFailureCauses = EmbedFailureCauses.EXCEEDED_TIMEOUT

MINIMIZATION: EmbedFailureCauses = EmbedFailureCauses.MINIMIZATION

KTERM_VIOLATION: EmbedFailureCauses = EmbedFailureCauses.KTERM_VIOLATION

CLASH: EmbedFailureCauses = EmbedFailureCauses.CLASH

def OrderedEmbedFailureCauses(legacyImplementation: bool = True) -> list: ...

class EmbedParameters:
    """Parameters controlling embedding"""

    def __init__(self) -> None: ...

    @property
    def maxIterations(self) -> int:
        """
        maximum number of embedding attempts to use for a
        single conformation
        """

    @maxIterations.setter
    def maxIterations(self, arg: int, /) -> None: ...

    @property
    def numThreads(self) -> int:
        """number of threads to use when embedding multiple conformations"""

    @numThreads.setter
    def numThreads(self, arg: int, /) -> None: ...

    @property
    def timeout(self) -> int:
        """
        maximum time in seconds to generate a conformer for a
        single molecule fragment. If set to 0, no timeout is set
        """

    @timeout.setter
    def timeout(self, arg: int, /) -> None: ...

    @property
    def randomSeed(self) -> int:
        """seed for the random number generator"""

    @randomSeed.setter
    def randomSeed(self, arg: int, /) -> None: ...

    @property
    def clearConfs(self) -> bool:
        """clear all existing conformations on the molecule"""

    @clearConfs.setter
    def clearConfs(self, arg: bool, /) -> None: ...

    @property
    def useRandomCoords(self) -> bool:
        """
        start the embedding from random coordinates instead of
        using eigenvalues of the distance matrix
        """

    @useRandomCoords.setter
    def useRandomCoords(self, arg: bool, /) -> None: ...

    @property
    def boxSizeMult(self) -> float:
        """determines the size of the box used for random coordinates"""

    @boxSizeMult.setter
    def boxSizeMult(self, arg: float, /) -> None: ...

    @property
    def randNegEig(self) -> bool:
        """
        if the embedding yields a negative eigenvalue, pick
        coordinates that correspond to this component at random
        """

    @randNegEig.setter
    def randNegEig(self, arg: bool, /) -> None: ...

    @property
    def numZeroFail(self) -> int:
        """fail embedding if we have at least this many zero eigenvalues"""

    @numZeroFail.setter
    def numZeroFail(self, arg: int, /) -> None: ...

    @property
    def optimizerForceTol(self) -> float:
        """
        the tolerance to be used during the distance-geometry
        force field minimization
        """

    @optimizerForceTol.setter
    def optimizerForceTol(self, arg: float, /) -> None: ...

    @property
    def basinThresh(self) -> float:
        """set the basin threshold for the DGeom force field."""

    @basinThresh.setter
    def basinThresh(self, arg: float, /) -> None: ...

    @property
    def ignoreSmoothingFailures(self) -> bool:
        """
        try and embed the molecule if if triangle smoothing of
        the bounds matrix fails
        """

    @ignoreSmoothingFailures.setter
    def ignoreSmoothingFailures(self, arg: bool, /) -> None: ...

    @property
    def enforceChirality(self) -> bool:
        """enforce correct chirilaty if chiral centers are present"""

    @enforceChirality.setter
    def enforceChirality(self, arg: bool, /) -> None: ...

    @property
    def useExpTorsionAnglePrefs(self) -> bool:
        """impose experimental torsion angle preferences"""

    @useExpTorsionAnglePrefs.setter
    def useExpTorsionAnglePrefs(self, arg: bool, /) -> None: ...

    @property
    def useBasicKnowledge(self) -> bool:
        """impose basic-knowledge constraints such as flat rings"""

    @useBasicKnowledge.setter
    def useBasicKnowledge(self, arg: bool, /) -> None: ...

    @property
    def ETversion(self) -> int:
        """version of the experimental torsion-angle preferences"""

    @ETversion.setter
    def ETversion(self, arg: int, /) -> None: ...

    @property
    def verbose(self) -> bool:
        """be verbose about configuration"""

    @verbose.setter
    def verbose(self, arg: bool, /) -> None: ...

    @property
    def pruneRmsThresh(self) -> float:
        """
        used to filter multiple conformations: keep only
        conformations that are at least this far apart from each other
        """

    @pruneRmsThresh.setter
    def pruneRmsThresh(self, arg: float, /) -> None: ...

    @property
    def onlyHeavyAtomsForRMS(self) -> bool:
        """Only consider heavy atoms when doing RMS filtering"""

    @onlyHeavyAtomsForRMS.setter
    def onlyHeavyAtomsForRMS(self, arg: bool, /) -> None: ...

    @property
    def embedFragmentsSeparately(self) -> bool:
        """split the molecule into fragments and embed them separately"""

    @embedFragmentsSeparately.setter
    def embedFragmentsSeparately(self, arg: bool, /) -> None: ...

    @property
    def useSmallRingTorsions(self) -> bool:
        """impose small ring torsion angle preferences"""

    @useSmallRingTorsions.setter
    def useSmallRingTorsions(self, arg: bool, /) -> None: ...

    @property
    def useMacrocycleTorsions(self) -> bool:
        """impose macrocycle torsion angle preferences"""

    @useMacrocycleTorsions.setter
    def useMacrocycleTorsions(self, arg: bool, /) -> None: ...

    @property
    def useMacrocycle14config(self) -> bool:
        """
        This forces amides and esters to be trans in macrocycles. This does not affect chain amides / esters!
        """

    @useMacrocycle14config.setter
    def useMacrocycle14config(self, arg: bool, /) -> None: ...

    @property
    def useLegacyImplementation(self) -> bool:
        """whether to use the combined minimization approach"""

    @useLegacyImplementation.setter
    def useLegacyImplementation(self, arg: bool, /) -> None: ...

    @property
    def boundsMatForceScaling(self) -> float:
        """
        scale the weights of the atom pair distance restraints relative to
        the other types of restraints
        """

    @boundsMatForceScaling.setter
    def boundsMatForceScaling(self, arg: float, /) -> None: ...

    @property
    def useSymmetryForPruning(self) -> bool:
        """
        use molecule symmetry when doing the RMSD pruning. Note that this
        option automatically also sets onlyHeavyAtomsForRMS to true.
        """

    @useSymmetryForPruning.setter
    def useSymmetryForPruning(self, arg: bool, /) -> None: ...

    def SetBoundsMat(self, boundsMatArg: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C')]) -> None:
        """
        set the distance-bounds matrix to be used (no triangle smoothing
        will be done on this) from a Numpy array
        """

    def SetCPCI(self, CPCIdict: dict[tuple[int, int], float]) -> None:
        """
        set the customised pairwise Columb-like interaction to atom pairs.
        used during structural minimisation stage
        """

    @property
    def forceTransAmides(self) -> bool:
        """
        This forces chain amides and esters to be trans. This does not affect amides / esters in macrocycles!
        """

    @forceTransAmides.setter
    def forceTransAmides(self, arg: bool, /) -> None: ...

    @property
    def trackFailures(self) -> bool:
        """keep track of which checks during the embedding process fail"""

    @trackFailures.setter
    def trackFailures(self, arg: bool, /) -> None: ...

    def GetFailureCounts(self) -> tuple:
        """returns the counts of each failure type"""

    @property
    def enableSequentialRandomSeeds(self) -> bool:
        """
        handle random number seeds so that conformer generation can be restarted
        """

    @enableSequentialRandomSeeds.setter
    def enableSequentialRandomSeeds(self, arg: bool, /) -> None: ...

    @property
    def symmetrizeConjugatedTerminalGroupsForPruning(self) -> bool:
        """symmetrize terminal conjugated groups for RMSD pruning"""

    @symmetrizeConjugatedTerminalGroupsForPruning.setter
    def symmetrizeConjugatedTerminalGroupsForPruning(self, arg: bool, /) -> None: ...

    def SetCoordMap(self, arg: Mapping[int, rdkit.Geometry.rdGeometry.Point3D], /) -> None:
        """sets the coordmap to be used"""

    def __setattr__(self, name: str, value: object | None) -> None: ...

def ETKDG() -> EmbedParameters:
    """Returns an EmbedParameters object for the ETKDG method - version 1."""

def ETKDGv2() -> EmbedParameters:
    """Returns an EmbedParameters object for the ETKDG method - version 2."""

def srETKDGv3() -> EmbedParameters:
    """
    Returns an EmbedParameters object for the ETKDG method -
    version 3 (small rings).
    """

def ETKDGv3() -> EmbedParameters:
    """
    Returns an EmbedParameters object for the ETKDG method -
    version 3 (macrocycles).
    """

def ETDG() -> EmbedParameters:
    """Returns an EmbedParameters object for the ETDG method."""

def ETDGv2() -> EmbedParameters:
    """Returns an EmbedParameters object for the ETDG method - version 2."""

def KDG() -> EmbedParameters:
    """Returns an EmbedParameters object for the KDG method."""

def DG() -> EmbedParameters:
    """Returns an EmbedParameters object for plain distance geometry."""

@overload
def GetMoleculeBoundsMatrix(mol: rdkit.Chem.rdchem.Mol, set15bounds: bool = True, scaleVDW: bool = False, doTriangleSmoothing: bool = True, useMacrocycle14config: bool = False, forceTransAmides: bool = True, set14bounds: bool = True, set13bounds: bool = True) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
    """
    Returns the distance bounds matrix for a molecule

    ARGUMENTS:

       - mol : the molecule of interest
       - set15bounds : set bounds for 1-5 atom distances based on
                       topology (otherwise stop at 1-4s)
       - scaleVDW : scale down the sum of VDW radii when setting the
                    lower bounds for atoms less than 5 bonds apart
       - doTriangleSmoothing : run triangle smoothing on the bounds
                    matrix before returning it
       - useMacrocycle14config : use 1-4 distance bound heuristics for macrocycles
       - forceTransAmides : constrain amide bonds to be trans
       - set14bounds : set bounds for 1-4 atom distances based on
                       topology
       - set13bounds : set bounds for 1-3 atom distances based on
                       topology

    RETURNS:

       the bounds matrix as a Numeric array with lower bounds in
       the lower triangle and upper bounds in the upper triangle
    """

@overload
def GetMoleculeBoundsMatrix(mol: rdkit.Chem.rdchem.Mol, embedParameters: EmbedParameters, doTriangleSmoothing: bool = True, scaleVDW: bool = False, set15bounds: bool = True, set14bounds: bool = True, set13bounds: bool = True) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
    """
    Returns the distance bounds matrix for a molecule

    ARGUMENTS:

       - mol : the molecule of interest
       - embedParameters : an EmbedParameters object
       - doTriangleSmoothing : run triangle smoothing on the bounds
       matrix before returning it
       - scaleVDW : scale down the sum of VDW radii when setting the
       lower bounds for atoms less than 5 bonds apart
       - set15bounds : set bounds for 1-5 atom distances based on
                       topology (otherwise stop at 1-4s)
       - set14bounds : set bounds for 1-4 atom distances based on
                       topology
       - set13bounds : set bounds for 1-3 atom distances based on
                       topology

    RETURNS:

       the bounds matrix as a Numeric array with lower bounds in
       the lower triangle and upper bounds in the upper triangle
    """

def EmbedParametersToJSON(embedParameters: EmbedParameters) -> str:
    """
    Returns json string containing embedParameters attributes

    ARGUMENTS:

      - embedParameters : the Params object you want serialized

    RETURNS:

      The Params object as json string
    """
