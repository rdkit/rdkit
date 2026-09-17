"""
Module containing implementation of Gaussian-based shape overlay and scoring.
NOTE: This functionality is experimental and the API and/or results may change in future releases.
"""

from collections.abc import Iterable, Sequence
import enum
from typing import overload

import rdkit.Chem.rdchem
import rdkit.Geometry.rdGeometry


class StartMode(enum.Enum):
    ROTATE_0 = 0

    ROTATE_180 = 1

    ROTATE_180_WIGGLE = 2

    ROTATE_45 = 3

    ROTATE_0_FRAGMENT = 4

    ROTATE_180_FRAGMENT = 5

    ROTATE_45_FRAGMENT = 6

    A_LA_PUBCHEM = 7

ROTATE_0: StartMode = StartMode.ROTATE_0

ROTATE_180: StartMode = StartMode.ROTATE_180

ROTATE_180_WIGGLE: StartMode = StartMode.ROTATE_180_WIGGLE

ROTATE_45: StartMode = StartMode.ROTATE_45

ROTATE_0_FRAGMENT: StartMode = StartMode.ROTATE_0_FRAGMENT

ROTATE_180_FRAGMENT: StartMode = StartMode.ROTATE_180_FRAGMENT

ROTATE_45_FRAGMENT: StartMode = StartMode.ROTATE_45_FRAGMENT

A_LA_PUBCHEM: StartMode = StartMode.A_LA_PUBCHEM

class OptimMode(enum.Enum):
    SHAPE_ONLY = 0

    SHAPE_PLUS_COLOR_SCORE = 1

    SHAPE_PLUS_COLOR = 2

SHAPE_ONLY: OptimMode = OptimMode.SHAPE_ONLY

SHAPE_PLUS_COLOR_SCORE: OptimMode = OptimMode.SHAPE_PLUS_COLOR_SCORE

SHAPE_PLUS_COLOR: OptimMode = OptimMode.SHAPE_PLUS_COLOR

class ShapeInputOptions:
    """ShapeInputOptions - options for setting up ShapeInput objects."""

    def __init__(self) -> None: ...

    @property
    def useColors(self) -> bool:
        """Whether to use color features in overlay.  Default=True."""

    @useColors.setter
    def useColors(self, arg: bool, /) -> None: ...

    @property
    def allCarbonRadii(self) -> bool:
        """
        Whether to use the same radius, appropriate for Carbon, for all atoms.  There is a
        slight accuracy penalty but significant speed gain if used.  Default=True.
        """

    @allCarbonRadii.setter
    def allCarbonRadii(self, arg: bool, /) -> None: ...

    @property
    def atomSubset(self) -> tuple:
        """
        If not empty, use just these atoms in the molecule to form the ShapeInput object.
        """

    @atomSubset.setter
    def atomSubset(self, arg: Iterable[int] | None) -> None: ...

    @property
    def customFeatures(self) -> tuple[tuple[tuple[int, rdkit.Geometry.rdGeometry.Point3D, float, list[int]], ...], ...]:
        """
        Custom features for the shape.  Requires a list of lists of tuples of
        int (the feature type), Point3D (the coordinates), float (the radius)
        and optionally a list of indices of the atoms that the feature was derived from.
        """

    @customFeatures.setter
    def customFeatures(self, arg: object, /) -> None: ...

    @property
    def atomRadii(self) -> tuple:
        """
        Non-standard radii to use for the atoms specified by their indices
        in the molecule.  Not all atoms need have a radius specified.
        A list of tuples of [int, float].
        """

    @atomRadii.setter
    def atomRadii(self, arg: Sequence[tuple[int, float]], /) -> None: ...

    @property
    def shapePruneThreshold(self) -> float:
        """
        If there is more than 1 conformer for the input molecule, prune the shapes so that none of them are more similar to each other than the threshold.  Default -1.0 means no pruning.
        """

    @shapePruneThreshold.setter
    def shapePruneThreshold(self, arg: float, /) -> None: ...

    @property
    def sortShapes(self) -> bool:
        """
        If True (the default), the shapes are sorted into descending order of total volume.
        """

    @sortShapes.setter
    def sortShapes(self, arg: bool, /) -> None: ...

    @property
    def includeDummies(self) -> bool:
        """Whether to include dummy atoms in the shape or not.  Default=True."""

    @includeDummies.setter
    def includeDummies(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class ShapeOverlayOptions:
    """
    ShapeOverlayOptions - options for controlling the shape overlay process.
    """

    def __init__(self) -> None: ...

    @property
    def startMode(self) -> StartMode:
        """
        Start modes for optimisation.  Default is A_LA_PUBCHEM - as used by the
        PubChem code - either ROTATE_180_WIGGLE or ROTATE_45 depending on the shape
        of the two molecules.  ROTATE_180_WIGGLE means 180 rotations about
        the x, y and z axes, then a small
        rotation about each axis from that point, using the best scoring one of
        those. ROTATE_180 uses 180 degree rotations for 4 start points,
        ROTATE_45 uses 45 degree rotations for 9 start points and ROTATE_0
        leaves the relative orientations of the 2 molecules as passed in before
        optimisation.  There are also ROTATE_0_FRAGMENT, ROTATE_45_FRAGMENT
        and ROTATE_180_FRAGMENT that as well as the above move the fit
        molecule to the ends of each of the principal axes and then does
        the appropriate rotations.  This is useful when the fit molecule is
        a lot smaller than the reference molecule, but requires a large number
        of optimisations so is relatively slow.
        """

    @startMode.setter
    def startMode(self, arg: StartMode, /) -> None: ...

    @property
    def optimMode(self) -> OptimMode:
        """
        Optimisation mode, controlling what parameters are used
        to drive the overlay.  Default=SHAPE_PLUS_COLOR_SCORE which
        optimises using just the overlap of shape, but uses the
        color to decide which is the best overlay.  Other options
        are SHAPE_ONLY and SHAPE_AND_COLOR with the latter using
        the overlap of color features as well.
        """

    @optimMode.setter
    def optimMode(self, arg: OptimMode, /) -> None: ...

    @property
    def simAlpha(self) -> float:
        """
        When doing a Tversky similarity, the alpha value.  If alpha and
        beta are both the default 1.0, it's a Tanimoto similarity.  A
        high alpha and low beta emphasize the fit volume in the
        similarity and vice versa. Tversky is O / (A * (R - O) + B * (F
        - O) + O) where O is the overlap volume, R is the reference's
        volume and F is the fit's volume.  This is different from that
        used by OpenEye (O / (A * R + B * F)).
        """

    @simAlpha.setter
    def simAlpha(self, arg: float, /) -> None: ...

    @property
    def simBeta(self) -> float:
        """When doing a Tversky similarity, the beta value."""

    @simBeta.setter
    def simBeta(self, arg: float, /) -> None: ...

    @property
    def optParam(self) -> float:
        """
        If using colors, the relative weights of the shape and color scores,
        as a fraction of 1.  Default=0.5.
        """

    @optParam.setter
    def optParam(self, arg: float, /) -> None: ...

    @property
    def nSteps(self) -> int:
        """Maximum number of steps for the shape overlay process. Default=100."""

    @nSteps.setter
    def nSteps(self, arg: int, /) -> None: ...

    @property
    def normalize(self) -> bool:
        """
        Whether to normalize the shapes before overlay by putting them into their
        canonical orientation (centred on the origin, aligned along its
        principal axes.  Default=True.
        """

    @normalize.setter
    def normalize(self, arg: bool, /) -> None: ...

    @property
    def useDistCutoff(self) -> bool:
        """
        Whether to use distance cutoff when calculating the shape volumes.  If used,
        there will be a small penalty in accuracy but a significant increase in speed.
        Default=True.
        """

    @useDistCutoff.setter
    def useDistCutoff(self, arg: bool, /) -> None: ...

    @property
    def distCutoff(self) -> float:
        """
        If using a distance cutoff, this is the value used.  Default=4.5 of whatever
        units the coordinates are in.
        """

    @distCutoff.setter
    def distCutoff(self, arg: float, /) -> None: ...

    @property
    def shapeConvergenceCriterion(self) -> float:
        """
        Optimisation stops when the shape Tversky score changes by less
        than this amount after an optimisation step.  A larger number is
        faster but gives less precise overlays.  Default=0.001.
        """

    @shapeConvergenceCriterion.setter
    def shapeConvergenceCriterion(self, arg: float, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class ShapeInput:
    """ShapeInput object"""

    def __init__(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1, shapeOpt: ShapeInputOptions = ..., overlayOpts: ShapeOverlayOptions = ...) -> None: ...

    @property
    def GetSmiles(self) -> str:
        """Get the SMILES string for the molecule that the shape relates to."""

    def setActiveShape(self, arg: int, /) -> None:
        """Set the active shape, the one that will be used for overlays etc."""

    def getActiveShape(self) -> int:
        """Return the number of the active shape."""

    @property
    def NumAtoms(self) -> int:
        """Get the number of atoms defining the shape."""

    @property
    def NumFeatures(self) -> int:
        """Get the number of features in the shape."""

    @property
    def NumShapes(self) -> int:
        """
        Get the number of shapes.  There will be a shape for each conformation of the input molecule, unless shape pruning was carried out in which case there may be fewer.
        """

    @property
    def ShapeVolume(self) -> float:
        """Get the shape's volume due to the atoms."""

    @property
    def ColorVolume(self) -> float:
        """Get the volume of the shape's color features."""

    def NormalizeCoords(self) -> None:
        """
        Align the principal axes to the cartesian axes and centre on the origin. Doesn't require that the shape was created from a molecule.  Creates the necessary transformation if not already done.
        """

    def ShapeToMol(self, includeColors: bool = False, withBonds: bool = True) -> rdkit.Chem.rdchem.Mol:
        """
        Return a molecule with coordinates of the current active shape.  If includeColors is True, (default is False) the color features will be added as xenon atoms.  If withBonds is True (the default) a molecule with bonds will be created, if not then just atoms at the appropriate positions will be produced.
        """

    def BestSimilarity(self, fitShape: ShapeInput, threshold: float = -1.0, overlayOpts: ShapeOverlayOptions | None = None) -> tuple[tuple[float, float, float], int, int, list[float]]:
        """
        Find the best similarity score between all shapes in this shape and the other one. Stops as soon as it gets something above the threshold. The score runs between 0.0 and 1.0, so the default threshold of -1.0 means no threshold. Fills in the shape numbers of the two that were responsible if there is something above the threshold, and the transformation that did it. Returns a tuple of the similarity scores ((-1.0, -1.0, -1.0) if there was nothing above the threshold), the number of the shape for this object and the shape number of the fitShape that gave the best similarity and the transformation matrix (as a list of 16 floats) that will reproduce the best overlay.  The shapes won't necessarily be left in the state that gave the best similarity.  Note that the shape numbers are not necessarily the same as the original molecule conformation numbers.
        """

    def MaxPossibleSimilarity(self, fitShape: ShapeInput, overlayOpts: ShapeOverlayOptions | None = None) -> float:
        """
        Get the maximum possible similarity score between all shapes in this shape and all shapes in the fitShape.  The maximum similarity is when one shape is entirely inside the other.  This returns the similarity in that case, which is the upper bound on what is achievable between these 2 shapes.
        """

    def __setattr__(self, name: str, value: object | None) -> None: ...

@overload
def AlignMol(ref: rdkit.Chem.rdchem.Mol, fit: rdkit.Chem.rdchem.Mol, refOpts: ShapeInputOptions | None = None, fitOpts: ShapeInputOptions | None = None, overlayOpts: ShapeOverlayOptions | None = None, refConfId: int = -1, fitConfId: int = -1) -> tuple:
    """
    Aligns a fit molecule onto a reference molecule.  The fit is modified.

    Parameters
    ----------
    ref: RDKit.ROMol
        Reference molecule
    fit: RDKit.ROMol
        Fit molecule that will be overlaid
    refOpts: ShapeInputOptions, optional
        Options for building the ref shape
    fitOpts: ShapeInputOptions, optional
        Options for building the fit shape
    overlayOpts: ShapeOverlayOptions, optional
        Options for controlling the overlay
    refConfId : int, optional
        Reference conformer ID (default is -1)
    fitConfId : int, optional
        fit conformer ID (default is -1)

    Returns
    -------
    3-tuple of floats
        The results are (combo_score, shape_score, color_score).  The color_score is
        0.0 if color features not used, in which case combo_score and shape_score will
        be the same.
    """

@overload
def AlignMol(refShape: ShapeInput, fit: rdkit.Chem.rdchem.Mol, fitOpts: ShapeInputOptions | None = None, overlayOpts: ShapeOverlayOptions | None = None, fitConfId: int = -1) -> tuple:
    """
    Aligns a fit molecule onto a reference shape.  The fit is modified.

    Parameters
    ----------
    refShape: ShapeInput
        Reference shape
    fit: RDKit.ROMol
        Fit molecule that will be overlaid
    fitOpts: ShapeInputOptions, optional
        Options for building the fit shape
    overlayOpts: ShapeOverlayOptions, optional
        Options for controlling the overlay
    fitConfId : int, optional
        Fit conformer ID (default is -1)

    Returns
    -------
    3-tuple of floats
        The results are (combo_score, shape_score, color_score).  The color_score is
        0.0 if color features not used, in which case combo_score and shape_score will
        be the same.
    """

def AlignShapes(refShape: ShapeInput, fitShape: ShapeInput, overlayOpts: ShapeOverlayOptions | None = None) -> tuple:
    """
    Aligns a fit shape to a reference shape. The fit is modified.

    Parameters
    ----------
    refShape : ShapeInput
        Reference shape
    fitShape : ShapeInput
        fit shape
    overlayOpts: ShapeOverlayOptions, optional
        Options for controlling the overlay

    Returns
    -------
     4-tuple of float, float, list of floats
        The results are (combo_score, shape_score, color_score, matrix)
        The matrix is a 16-float list giving the transformation matrix that
        overlays the fit onto the reference.
    """

@overload
def ScoreMol(ref: rdkit.Chem.rdchem.Mol, fit: rdkit.Chem.rdchem.Mol, refOpts: ShapeInputOptions | None = None, fitOpts: ShapeInputOptions | None = None, overlayOpts: ShapeOverlayOptions | None = None, refConfId: int = -1, fitConfId: int = -1) -> tuple:
    """
    Calculates the scores between a reference molecule and a fit
    molecule without overlay.

    Parameters
    ----------
    ref: RDKit.ROMol
        Reference molecule
    fit: RDKit.ROMol
        Fit molecule that will be scored
    refOpts: ShapeInputOptions, optional
        Options for building the ref shape
    fitOpts: ShapeInputOptions, optional
        Options for building the fit shape
    overlayOpts: ShapeOverlayOptions, optional
        Options for controlling the volume calculation
    refConfId : int, optional
        Reference conformer ID (default is -1)
    fitConfId : int, optional
        fit conformer ID (default is -1)

    Returns
    -------
    3-tuple of floats
        The results are (combo_score, shape_score, color_score).  The color_score is
        0.0 if color features not used, in which case combo_score and shape_score will
        be the same.
    """

@overload
def ScoreMol(refShape: ShapeInput, fit: rdkit.Chem.rdchem.Mol, fitOpts: ShapeInputOptions | None = None, overlayOpts: ShapeOverlayOptions | None = None, fitConfId: int = -1) -> tuple:
    """
    Calculates the scores between a reference shape and a fit molecule
    without overlay.

    Parameters
    ----------
    refShape: ShapeInput
        Reference shape
    fit: RDKit.ROMol
        Fit molecule that will be scored
    fitOpts: ShapeInputOptions, optional
        Options for building the fit shape
    overlayOpts: ShapeOverlayOptions, optional
        Options for controlling the volume calculation
    fitConfId : int, optional
        fit conformer ID (default is -1)

    Returns
    -------
    3-tuple of floats
        The results are (combo_score, shape_score, color_score).  The color_score is
        0.0 if color features not used, in which case combo_score and shape_score will
        be the same.
    """

def ScoreShape(refShape: ShapeInput, fitShape: ShapeInput, overlayOpts: ShapeOverlayOptions | None = None) -> tuple:
    """
    Calculates the scores between a reference shape and a fit shape without
    overlay.

    Parameters
    ----------
    refShape: ShapeInput
        Reference shape
    fitShape: ShapeInput
        Fit shape
    fitOpts: ShapeInputOptions, optional
        Options for building the fit shape
    overlayOpts: ShapeOverlayOptions, optional
        Options for controlling the volume calculation

    Returns
    -------
    3-tuple of floats
        The results are (combo_score, shape_score, color_score).  The color_score is
        0.0 if color features not used, in which case combo_score and shape_score will
        be the same.
    """

def ScoreMoleculeAllConformers(ref: rdkit.Chem.rdchem.Mol, fit: rdkit.Chem.rdchem.Mol, refOpts: ShapeInputOptions | None = None, fitOpts: ShapeInputOptions | None = None, overlayOpts: ShapeOverlayOptions | None = None) -> tuple[list[list[float]], int, int, list[float]]:
    """
    Calculate the scores for the alignment of all conformers
     of the fit molecule onto the reference.  The molecules themselves are not
     altered.

    Parameters
    ----------
    ref: RDKit.ROMol
        Reference molecule
    fit: RDKit.ROMol
        Fit molecule that will be scored
    refOpts: ShapeInputOptions, optional
        Options for building the ref shape
    fitOpts: ShapeInputOptions, optional
        Options for building the fit shape
    overlayOpts: ShapeOverlayOptions, optional
        Options for controlling the volume calculation

    Returns
    -------
    A complex tuple containing:
        A tuple of tuples containing the scores from aligning the fit conformations
        onto the reference conformations.  scores[0][1] is the score of aligning
        fit conformation 1 onto ref conformation 0.
        The ID of the ref conformer from the best-scoring alignment
        The ID of the fit conformer from the best-scoring alignment
        The transformation that gives the best-scoring alignment for those
        conformers as a 16-float tuple.
    """
