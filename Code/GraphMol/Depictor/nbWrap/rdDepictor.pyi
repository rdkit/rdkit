"""
Module containing the functionality to compute 2D coordinates for a molecule
"""

from collections.abc import Iterable, Sequence
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdchem
import rdkit.Geometry.rdGeometry


class DepictException(ValueError):
    pass

def IsCoordGenSupportAvailable() -> bool:
    """Returns whether RDKit was built with CoordGen support."""

def SetPreferCoordGen(val: bool) -> None:
    """
    Sets whether or not the CoordGen library should be preferred to the RDKit depiction library.
    """

def SetRingSystemTemplates(templatePath: str) -> None:
    """
    Loads the ring system templates from the specified file to be
    used in 2D coordinate generation. Each template must be a single
    line in the file represented using CXSMILES, and the structure should
    be a single ring system. Throws a DepictException if any templates
    are invalid.
    """

def AddRingSystemTemplates(templatePath: str) -> None:
    """
    Adds the ring system templates from the specified file to be
    used in 2D coordinate generation. If there are duplicates, the most
    recently added template will be used. Each template must be a single
    line in the file represented using CXSMILES, and the structure should
    be a single ring system. Throws a DepictException if any templates
    are invalid.
    """

def LoadDefaultRingSystemTemplates() -> None:
    """
    Loads the default ring system templates and removes existing ones, if present.
    """

def GetPreferCoordGen() -> bool:
    """
    Return whether or not the CoordGen library is used for coordinate generation in the RDKit depiction library.
    """

class UsingCoordGen:
    """
    Context manager to temporarily set CoordGen library preference in RDKit depiction.
    """

    def __init__(self, temp_state: bool) -> None:
        """Constructor"""

    def __enter__(self) -> None: ...

    def __exit__(self, exc_type: object | None = None, exc_val: object | None = None, traceback: object | None = None) -> None: ...

class ConstrainedDepictionParams:
    """Parameters controlling constrained depiction"""

    def __init__(self) -> None: ...

    @property
    def acceptFailure(self) -> bool:
        """
        if False (default), a DepictException is thrown if the molecule
        does not have a substructure match to the reference;
        if True, an unconstrained depiction will be generated
        """

    @acceptFailure.setter
    def acceptFailure(self, arg: bool, /) -> None: ...

    @property
    def forceRDKit(self) -> bool:
        """
        if True, use RDKit to generate coordinates even if preferCoordGen
        is set to True; defaults to False
        """

    @forceRDKit.setter
    def forceRDKit(self, arg: bool, /) -> None: ...

    @property
    def allowRGroups(self) -> bool:
        """
        if True, terminal dummy atoms in the reference are ignored
        if they match an implicit hydrogen in the molecule or if they are
        attached top a query atom; defaults to False
        """

    @allowRGroups.setter
    def allowRGroups(self, arg: bool, /) -> None: ...

    @property
    def alignOnly(self) -> bool:
        """
        if False (default), a part of the molecule is hard-constrained
        to have the same coordinates as the reference, and the rest of
        the molecule is built around it; if True, coordinates
        from conformation existingConfId are preserved (if they exist)
        or generated without constraints (if they do not exist), then
        the conformation is rigid-body aligned to the reference
        """

    @alignOnly.setter
    def alignOnly(self, arg: bool, /) -> None: ...

    @property
    def adjustMolBlockWedging(self) -> bool:
        """
        if True (default), existing wedging information
        will be updated or cleared as required; if False,
        existing molblock wedging information will
        always be preserved
        """

    @adjustMolBlockWedging.setter
    def adjustMolBlockWedging(self, arg: bool, /) -> None: ...

    @property
    def existingConfId(self) -> int:
        """
        conformation id whose 2D coordinates should be
        rigid-body aligned to the reference (if alignOnly is True),
        or used to determine whether existing molblock wedging information
        can be preserved following the constrained depiction (if
        adjustMolBlockWedging is True
        """

    @existingConfId.setter
    def existingConfId(self, arg: int, /) -> None: ...

    @property
    def useRingTemplates(self) -> bool:
        """use templates to generate coordinates of complex ring systems"""

    @useRingTemplates.setter
    def useRingTemplates(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

def Compute2DCoords(mol: rdkit.Chem.rdchem.Mol, canonOrient: bool = True, clearConfs: bool = True, coordMap: dict[int, rdkit.Geometry.rdGeometry.Point2D] = {}, nFlipsPerSample: int = 0, nSample: int = 0, sampleSeed: int = 0, permuteDeg4Nodes: int = 0, bondLength: float = -1.0, forceRDKit: bool = False, useRingTemplates: bool = False) -> int:
    """
    Compute 2D coordinates for a molecule.
      The resulting coordinates are stored on each atom of the molecule

      ARGUMENTS:

         mol - the molecule of interest
         canonOrient - orient the molecule in a canonical way
         clearConfs - if true, all existing conformations on the molecule
                 will be cleared
         coordMap - a dictionary mapping atom Ids -> Point2D objects
                    with starting coordinates for atoms that should
                    have their positions locked.
         nFlipsPerSample - number of rotatable bonds that are
                    flipped at random at a time.
         nSample - Number of random samplings of rotatable bonds.
         sampleSeed - seed for the random sampling process.
         permuteDeg4Nodes - allow permutation of bonds at a degree 4
                     node during the sampling process
         bondLength - change the default bond length for depiction
         forceRDKit - use RDKit to generate coordinates even if
                      preferCoordGen is set to true
         useRingTemplates - use templates to generate coordinates of complex
                      ring systems

      RETURNS:

         ID of the conformation added to the molecule
    """

def Compute2DCoordsMimicDistmat(mol: rdkit.Chem.rdchem.Mol, distMat: Annotated[NDArray[numpy.float64], dict(shape=(None,))], canonOrient: bool = False, clearConfs: bool = True, weightDistMat: float = 0.5, nFlipsPerSample: int = 3, nSample: int = 100, sampleSeed: int = 100, permuteDeg4Nodes: bool = True, bondLength: float = -1.0, forceRDKit: bool = False) -> int:
    """
    Compute 2D coordinates for a molecule such
      that the inter-atom distances mimic those in a user-provided
      distance matrix.
      The resulting coordinates are stored on each atom of the molecule

      ARGUMENTS:

         mol - the molecule of interest
         distMat - distance matrix that we want the 2D structure to mimic
         canonOrient - orient the molecule in a canonical way
         clearConfs - if true, all existing conformations on the molecule
                 will be cleared
         weightDistMat - weight assigned in the cost function to mimicking
                         the distance matrix.
                         This must be between (0.0,1.0). (1.0-weightDistMat)
                         is then the weight assigned to improving
                         the density of the 2D structure i.e. try to
                         make it spread out
         nFlipsPerSample - number of rotatable bonds that are
                    flipped at random at a time.
         nSample - Number of random samplings of rotatable bonds.
         sampleSeed - seed for the random sampling process.
         permuteDeg4Nodes - allow permutation of bonds at a degree 4
                     node during the sampling process
         bondLength - change the default bond length for depiction
         forceRDKit - use RDKit to generate coordinates even if
                      preferCoordGen is set to true

      RETURNS:

         ID of the conformation added to the molecule
    """

@overload
def GenerateDepictionMatching2DStructure(mol: rdkit.Chem.rdchem.Mol, reference: rdkit.Chem.rdchem.Mol, confId: int = -1, refPatt: rdkit.Chem.rdchem.Mol | None = None, params: ConstrainedDepictionParams | None = None) -> tuple:
    """
    Generate a depiction for a molecule where a piece of the
      molecule is constrained to have the same coordinates as a reference.

      The constraint can be hard (default) or soft.

      Hard (default, ConstrainedDepictionParams.alignOnly=False):
      Existing molecule coordinates, if present, are discarded;
      new coordinates are generated constraining a piece of the molecule
      to have the same coordinates as the reference, while the rest of
      the molecule is built around it.
      If ConstrainedDepictionParams.adjustMolBlockWedging is False
      (default), existing molblock wedging information is always preserved.
      If ConstrainedDepictionParams.adjustMolBlockWedging is True,
      existing molblock wedging information is preserved in case it
      only involves the invariant core and the core conformation has not
      changed, while it is cleared in case the wedging is also outside
      the invariant core, or core coordinates were changed.
      If ConstrainedDepictionParams.acceptFailure is set to True and no
      substructure match is found, coordinates will be recomputed from
      scratch, hence molblock wedging information will be cleared.

      Soft (ConstrainedDepictionParams.alignOnly=True):
      Existing coordinates in the conformation identified by
      ConstrainedDepictionParams.existingConfId are preserved if present,
      otherwise unconstrained new coordinates are generated.
      Subsequently, coodinates undergo a rigid-body alignment to the reference.
      If ConstrainedDepictionParams.adjustMolBlockWedging is False
      (default), existing molblock wedging information is always preserved.
      If ConstrainedDepictionParams.adjustMolBlockWedging is True,
      existing molblock wedging information is inverted in case the rigid-body
      alignment involved a flip around the Z axis.

      This is useful, for example, for generating depictions of SAR data
      sets such that the cores of the molecules are all oriented the same way.
      ARGUMENTS:

      mol -    the molecule to be aligned, this will come back
               with a single conformer.
      reference -    a molecule with the reference atoms to align to;
                     this should have a depiction.
      confId -       (optional) the id of the reference conformation to use
      refPatt -      (optional) a query molecule to be used to generate
                     the atom mapping between the molecule and the reference
      params - (optional) a ConstrainedDepictionParams instance

      RETURNS: a tuple of (refIdx, molIdx) tuples corresponding to the atom
               indices in mol constrained to have the same coordinates as atom
               indices in reference.
    """

@overload
def GenerateDepictionMatching2DStructure(mol: rdkit.Chem.rdchem.Mol, reference: rdkit.Chem.rdchem.Mol, confId: int = -1, refPatt: rdkit.Chem.rdchem.Mol | None = None, acceptFailure: bool = False, forceRDKit: bool = False, allowRGroups: bool = False) -> tuple:
    """
    Generate a depiction for a molecule where a piece of the
      molecule is constrained to have the same coordinates as a reference.

      This is useful, for example, for generating depictions of SAR data
      sets such that the cores of the molecules are all oriented the same way.
      ARGUMENTS:

      mol -    the molecule to be aligned, this will come back
               with a single conformer.
      reference -    a molecule with the reference atoms to align to;
                     this should have a depiction.
      confId -       the id of the reference conformation to use
      refPatt -      a query molecule to be used to generate
                     the atom mapping between the molecule and the reference
      acceptFailure - if True, standard depictions will be generated
                      for molecules that don't have a substructure match to the
                      reference; if False, throws a DepictException.
      forceRDKit -    (optional) use RDKit to generate coordinates even if
                      preferCoordGen is set to true
      allowRGroups -  (optional) if True, terminal dummy atoms in the
                      reference are ignored if they match an implicit
                      hydrogen in the molecule, and a constrained
                      depiction is still attempted

      RETURNS: a tuple of (refIdx, molIdx) tuples corresponding to the atom
               indices in mol constrained to have the same coordinates as atom
               indices in reference.
    """

@overload
def GenerateDepictionMatching2DStructure(mol: rdkit.Chem.rdchem.Mol, reference: rdkit.Chem.rdchem.Mol, atomMap: Iterable[Sequence[int]], confId: int = -1, params: ConstrainedDepictionParams | None = None) -> None:
    """
    Generate a depiction for a molecule where a piece of the
      molecule is constrained to have the same coordinates as a reference.

      This is useful for, for example, generating depictions of SAR data
      sets so that the cores of the molecules are all oriented the same way.
      ARGUMENTS:

      mol -    the molecule to be aligned, this will come back
               with a single conformer.
      reference -    a molecule with the reference atoms to align to;
                     this should have a depiction.
      atomMap -      a sequence of (queryAtomIdx, molAtomIdx) pairs that will
                     be used to generate the atom mapping between the molecule
                     and the reference. Note that this sequence can be shorter
                     than the number of atoms in the reference.
      confId -       (optional) the id of the reference conformation to use
      params -       (optional) an instance of ConstrainedDepictionParams
    """

@overload
def GenerateDepictionMatching2DStructure(mol: rdkit.Chem.rdchem.Mol, reference: rdkit.Chem.rdchem.Mol, atomMap: Iterable[Sequence[int]], confId: int, forceRDKit: bool) -> None:
    """
    Generate a depiction for a molecule where a piece of the
      molecule is constrained to have the same coordinates as a reference.

      This is useful for, for example, generating depictions of SAR data
      sets so that the cores of the molecules are all oriented the same way.
      ARGUMENTS:

      mol -    the molecule to be aligned, this will come back
               with a single conformer.
      reference -    a molecule with the reference atoms to align to;
                     this should have a depiction.
      atomMap -      a sequence of (queryAtomIdx, molAtomIdx) pairs that will
                     be used to generate the atom mapping between the molecule
                     and the reference. Note that this sequence can be shorter
                     than the number of atoms in the reference.
      confId -       the id of the reference conformation to use
      forceRDKit -   use RDKit to generate coordinates even if
                     preferCoordGen is set to true
    """

def GenerateDepictionMatching3DStructure(mol: rdkit.Chem.rdchem.Mol, reference: rdkit.Chem.rdchem.Mol, confId: int = -1, refPatt: rdkit.Chem.rdchem.Mol | None = None, acceptFailure: bool = False, forceRDKit: bool = False) -> None:
    """
    Generate a depiction for a molecule where a piece of the molecule
      is constrained to have coordinates similar to those of a 3D reference
      structure.
      ARGUMENTS:

      mol -    the molecule to be aligned, this will come back
               with a single conformer containing the 2D coordinates.
      reference -    a molecule with the reference atoms to align to.
                     By default this should be the same as mol, but with
                     3D coordinates
      confId -       (optional) the id of the reference conformation to use
      referencePattern -  (optional) a query molecule to map a subset of
                          the reference onto the mol, so that only some of the
                          atoms are aligned.
      acceptFailure - (optional) if True, standard depictions will be generated
                      for molecules that don't match the reference or the
                      referencePattern; if False, throws a DepictException.
      forceRDKit -    (optional) use RDKit to generate coordinates even if
                      preferCoordGen is set to true
    """

def StraightenDepiction(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, minimizeRotation: bool = False) -> None:
    """
    Rotate the 2D depiction such that the majority of bonds have a
      30-degree angle with the X axis.
      ARGUMENTS:

      mol              - the molecule to be rotated.
      confId           - (optional) the id of the reference conformation to use.
      minimizeRotation - (optional) if False (the default), the molecule
                         is rotated such that the majority of bonds have an angle
                         with the X axis of 30 or 90 degrees. If True, the minimum
                         rotation is applied such that the majority of bonds have
                         an angle with the X axis of 0, 30, 60, or 90 degrees,
                         with the goal of altering the initial orientation as
                         little as possible .
    """

def NormalizeDepiction(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, canonicalize: int = 1, scaleFactor: float = -1.0) -> float:
    """
    Normalizes the 2D depiction.
    If canonicalize is != 0, the depiction is subjected to a canonical
    transformation such that its main axis is aligned along the X axis
    (canonicalize >0, the default) or the Y axis (canonicalize <0).
    If canonicalize is 0, no canonicalization takes place.
    If scaleFactor is <0.0 (the default) the depiction is scaled such
    that bond lengths conform to RDKit standards. The applied scaling
    factor is returned.

    ARGUMENTS:

    mol          - the molecule to be normalized
    confId       - (optional) the id of the reference conformation to use
    canonicalize - (optional) if != 0, a canonical transformation is
                   applied: if >0 (the default), the main molecule axis is
                   aligned to the X axis, if <0 to the Y axis.
                   If 0, no canonical transformation is applied.
    scaleFactor  - (optional) if >0.0, the scaling factor to apply. The default
                   (-1.0) means that the depiction is automatically scaled
                   such that bond lengths are the standard RDKit ones.

    RETURNS: the applied scaling factor.
    """
