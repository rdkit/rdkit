"""Module containing PubChem shape alignment functionality."""

from collections.abc import Iterable, Sequence
from typing import overload

import rdkit.Chem.rdchem


class ShapeInputOptions:
    """Shape Input Options"""

    def __init__(self) -> None: ...

    @property
    def useColors(self) -> bool:
        """
        Whether to use colors (pharmacophore features) in the score.  Default=True.
        """

    @useColors.setter
    def useColors(self, arg: bool, /) -> None: ...

    @property
    def includeDummies(self) -> bool:
        """Whether to use dummy atoms in the alignment. Default=False."""

    @includeDummies.setter
    def includeDummies(self, arg: bool, /) -> None: ...

    @property
    def dummyRadius(self) -> float:
        """
        If using dummy atoms in the alignment, what radius to use for them.
          Default=2.16 (the radius of Xe).
        """

    @dummyRadius.setter
    def dummyRadius(self, arg: float, /) -> None: ...

    @property
    def atomSubset(self) -> tuple:
        """
        If not empty, use just these atoms in the molecule to form the ShapeInput object.
        """

    @atomSubset.setter
    def atomSubset(self, arg: Iterable[int] | None) -> None: ...

    @property
    def notColorAtoms(self) -> tuple:
        """
        Any atoms mentioned here by index should not be used in a color feature.
        """

    @notColorAtoms.setter
    def notColorAtoms(self, arg: Iterable[int] | None) -> None: ...

    @property
    def atomRadii(self) -> tuple:
        """
        Non-standard radii to use for the atoms specified by their indices in the molecule.  A list of tuples of [int, float].
        """

    @atomRadii.setter
    def atomRadii(self, arg: Sequence[tuple[int, float]], /) -> None: ...

    @property
    def customFeatures(self) -> tuple:
        """Custom features for the shape."""

    @customFeatures.setter
    def customFeatures(self, arg: object, /) -> None: ...

    @property
    def normalize(self) -> bool:
        """
        Whether to normalise the shape by putting into
        its inertial frame.  Default=True.
        """

    @normalize.setter
    def normalize(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class ShapeInput:
    @property
    def coord(self) -> list[float]: ...

    @coord.setter
    def coord(self, arg: Sequence[float], /) -> None: ...

    @property
    def alpha_vector(self) -> list[float]: ...

    @alpha_vector.setter
    def alpha_vector(self, arg: Sequence[float], /) -> None: ...

    @property
    def atom_type_vector(self) -> list[int]: ...

    @atom_type_vector.setter
    def atom_type_vector(self, arg: Sequence[int], /) -> None: ...

    @property
    def volumeAtomIndexVector(self) -> list[int]: ...

    @volumeAtomIndexVector.setter
    def volumeAtomIndexVector(self, arg: Sequence[int], /) -> None: ...

    @property
    def shift(self) -> list:
        """Translation of centre of shape coordinates to origin."""

    @shift.setter
    def shift(self, arg: Iterable[float] | None) -> None: ...

    @property
    def inertialRot(self) -> list:
        """
        Rotation applied to put the shape into its principal axes frame of reference.
        """

    @property
    def sov(self) -> float: ...

    @sov.setter
    def sov(self, arg: float, /) -> None: ...

    @property
    def sof(self) -> float: ...

    @sof.setter
    def sof(self, arg: float, /) -> None: ...

@overload
def AlignMol(ref: rdkit.Chem.rdchem.Mol, probe: rdkit.Chem.rdchem.Mol, refConfId: int = -1, probeConfId: int = -1, useColors: bool = True, opt_param: float = 1.0, max_preiters: int = 10, max_postiters: int = 30) -> tuple:
    """
    Aligns a probe molecule to a reference molecule. The probe is modified.

    Parameters
    ----------
    ref : RDKit.ROMol
        Reference molecule
    probe : RDKit.ROMol
        Probe molecule
    refConfId : int, optional
        Reference conformer ID (default is -1)
    probeConfId : int, optional
        Probe conformer ID (default is -1)
    useColors : bool, optional
        Whether or not to use colors in the scoring (default is True)
    opt_param : float, optional
        Balance of shape and color for optimization.
        0 is only color, 0.5 is equal weight, and 1.0 is only shape.
        Default is 1.0.
    max_preiters : int, optional
        In the two phase optimization, the maximum iterations done on all poses.
    max_postiters : int, optional
        In the two phase optimization, the maximum iterations during the second phase on
        only the best poses from the first phase


    Returns
    -------
     2-tuple of doubles
        The results are (shape_score, color_score)
        The color_score is zero if opt_param is 1.0.
    """

@overload
def AlignMol(ref: rdkit.Chem.rdchem.Mol, probe: rdkit.Chem.rdchem.Mol, refShapeOpts: ShapeInputOptions, probeShapeOpts: ShapeInputOptions, refConfId: int = -1, probeConfId: int = -1, opt_param: float = 1.0, max_preiters: int = 10, max_postiters: int = 30) -> tuple:
    """
    Aligns a probe molecule to a reference molecule. The probe is modified.

    Parameters
    ----------
    ref : RDKit.ROMol
        Reference molecule
    probe : RDKit.ROMol
        Probe molecule
    refShapeOpts : ShapeInputOptions
        Options for constructing the shape for the reference molecule
    probeShapeOpts : ShapeInputOptions
        Options for constructing the shape for the probe molecule
    refConfId : int, optional
        Reference conformer ID (default is -1)
    probeConfId : int, optional
        Probe conformer ID (default is -1)
    opt_param : float, optional
        Balance of shape and color for optimization.
        0 is only color, 0.5 is equal weight, and 1.0 is only shape.
        Default is 1.0.
    max_preiters : int, optional
        In the two phase optimization, the maximum iterations done on all poses.
    max_postiters : int, optional
        In the two phase optimization, the maximum iterations during the second phase on
        only the best poses from the first phase


    Returns
    -------
     2-tuple of doubles
        The results are (shape_score, color_score)
        The color_score is zero if opt_param is 1.0.
    """

@overload
def AlignMol(refShape: ShapeInput, probe: rdkit.Chem.rdchem.Mol, probeConfId: int = -1, useColors: bool = True, opt_param: float = 1.0, max_preiters: int = 10, max_postiters: int = 30, applyRefShift: bool = False) -> tuple:
    """
    Aligns a probe molecule to a reference shape. The probe is modified.
    Assumes the shapes are both centred on the origin.

    Parameters
    ----------
    refShape : ShapeInput
        Reference shape
    probe : RDKit.ROMol
        Probe molecule
    probeConfId : int, optional
        Probe conformer ID (default is -1)
    useColors : bool, optional
        Whether or not to use colors in the scoring (default is True)
    opt_param : float, optional
        Balance of shape and color for optimization.
        0 is only color, 0.5 is equal weight, and 1.0 is only shape.
        Default is 1.0.
    max_preiters : int, optional
        In the two phase optimization, the maximum iterations done on all poses.
    max_postiters : int, optional
        In the two phase optimization, the maximum iterations during the second phase on
        only the best poses from the first phase
    applyRefShift : bool, optional
        If True, apply the reference shape's shift translation to the final
        coordinates.


    Returns
    -------
     2-tuple of doubles
        The results are (shape_score, color_score)
        The color_score is zero if opt_param is 1.0.
    """

def AlignShapes(refShape: ShapeInput, probeShape: ShapeInput, opt_param: float = 1.0, max_preiters: int = 10, max_postiters: int = 30) -> tuple:
    """
    Aligns a probe shape to a reference shape. The probe is modified.

    Parameters
    ----------
    refShape : ShapeInput
        Reference shape
    probeShape : ShapeInput
        Probe shape
    opt_param : float, optional
        Balance of shape and color for optimization.
        0 is only color, 0.5 is equal weight, and 1.0 is only shape
    max_preiters : int, optional
        In the two phase optimization, the maximum iterations done on all poses.
    max_postiters : int, optional
        In the two phase optimization, the maximum iterations during the second phase on
        only the best poses from the first phase


    Returns
    -------
     3-tuple of double, double, list of doubles
        The results are (shape_score, color_score, matrix)
        The matrix is a 12-float list giving the transformation matrix that
        overlays the probe onto the reference.
    """

def TransformConformer(finalTrans: Iterable[float], finalRot: Iterable[float], matrix: Iterable[float], probeShape: ShapeInput, probeConformer: rdkit.Chem.rdchem.Conformer) -> None:
    """
    Assuming that probeShape has been overlaid onto refShape to give
    the supplied transformation matrix, applies that transformation to the
     given conformer.

    Parameters
    ----------
    finalTrans : list[float * 3]
        The final translation to apply to conformer.
    matrix: list[float * 12]
        The transformation matrix
    probeShape : ShapeInput
        Probe shape
    probeConformer : Conformer
        Probe conformer
    """

def PrepareConformer(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, opts: ShapeInputOptions | None = None) -> ShapeInput:
    """
    Generates a ShapeInput object for a molecule

    Parameters
    ----------
    mol : RDKit.ROMol
        Reference molecule
    confId : int, optional
        Conformer ID to use (default is -1)
    opts : ShapeInputOptions, optional
        Options for Shapeinput

    Returns
    -------
     a ShapeInput for the molecule
    """

@overload
def ScoreMol(mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol, mol1ShapeOpts: ShapeInputOptions, mol2ShapeOpts: ShapeInputOptions, mol1ConfId: int = -1, mol2ConfId: int = -1) -> tuple:
    """
    Calculate the scores between a shape and a molecule without moving them.

    Parameters
    ----------
    mol1 : RDKit.ROMol
        First molecule
    mol2 : RDKit.ROMol
        Second molecule
    mol1ShapeOptions:
        Options for constructing the shape for molecule 1
    mol2ShapeOptions:
        Options for constructing the shape for molecule 2
    mol1ConfId : int, optional
        First molecule conformer ID (default is -1)
    mol2ConfId : int, optional
        Second conformer ID (default is -1)


    Returns
    -------
     2-tuple of doubles
        The results are (shape_score, color_score)
        The color_score is zero if useColors is False for either of the
    shape options
    """

@overload
def ScoreMol(shape: ShapeInput, mol: rdkit.Chem.rdchem.Mol, molShapeOpts: ShapeInputOptions, molConfId: int = -1) -> tuple:
    """
    Calculate the scores between 2 molecules without moving them.

    Parameters
    ----------
    shape : ShapeInput
        Shape
    mol : RDKit.ROMol
        Molecule
    molShapeOptions:
        Options for constructing the shape for molecule
    molConfId : int, optional
        Molecule conformer ID (default is -1)

    Returns
    -------
     2-tuple of doubles
        The results are (shape_score, color_score)
        The color_score is zero if shape.useColors is False
    """

def ScoreShape(shape1: ShapeInput, shape2: ShapeInput, useColors: bool = False) -> tuple:
    """
    Calculate the scores between 2 shapes without moving them.

    Parameters
    ----------
    shape1 : ShapeInput
        Shape
    shape2 : ShapeInput
        Shape
    useColors : bool
        Whether to use colors for the score or not.
    Returns
    -------
     2-tuple of doubles
        The results are (shape_score, color_score)
        The color_score is zero if useColors is False
    """
