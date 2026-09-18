"""Module containing functions to compute molecular descriptors"""

from collections.abc import Iterable, Sequence
import enum
from typing import Final, overload

import rdkit.Chem.rdchem
import rdkit.DataStructs.cDataStructs
import rdkit.Geometry.rdGeometry


class AtomPairsParameters:
    version: Final[str] = ...
    """(arg: object, /) -> str"""

    numTypeBits: Final[int] = ...
    """(arg: object, /) -> int"""

    numPiBits: Final[int] = ...
    """(arg: object, /) -> int"""

    numBranchBits: Final[int] = ...
    """(arg: object, /) -> int"""

    numChiralBits: Final[int] = ...
    """(arg: object, /) -> int"""

    codeSize: Final[int] = ...
    """(arg: object, /) -> int"""

    atomTypes: Final[list[int]] = ...
    """(arg: object, /) -> list[int]"""

    numPathBits: Final[int] = ...
    """(arg: object, /) -> int"""

    numAtomPairFingerprintBits: Final[int] = ...
    """(arg: object, /) -> int"""

    def __setattr__(self, name: str, value: object | None) -> None: ...

def GetAtomPairAtomCode(atom: rdkit.Chem.rdchem.Atom, branchSubtract: int = 0, includeChirality: bool = False) -> int:
    """Returns the atom code (hash) for an atom"""

def GetAtomPairCode(atom1Code: int, atom2Code: int, distance: int, includeChirality: bool = False) -> int:
    """
    Returns the atom-pair code (hash) for a pair of atoms separated by a
    certain number of bonds
    """

def GetAtomPairFingerprint(mol: rdkit.Chem.rdchem.Mol, minLength: int = 1, maxLength: int = 30, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, atomInvariants: Iterable[int] | None = None, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> rdkit.DataStructs.cDataStructs.IntSparseIntVect:
    """
    Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect
    """

def GetHashedAtomPairFingerprint(mol: rdkit.Chem.rdchem.Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, atomInvariants: Iterable[int] | None = None, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> rdkit.DataStructs.cDataStructs.IntSparseIntVect:
    """
    Returns the hashed atom-pair fingerprint for a molecule as an IntSparseIntVect
    """

def GetHashedAtomPairFingerprintAsBitVect(mol: rdkit.Chem.rdchem.Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, atomInvariants: Iterable[int] | None = None, nBitsPerEntry: int = 4, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """Returns the atom-pair fingerprint for a molecule as an ExplicitBitVect"""

def GetTopologicalTorsionFingerprint(mol: rdkit.Chem.rdchem.Mol, targetSize: int = 4, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, atomInvariants: Iterable[int] | None = None, includeChirality: bool = False) -> rdkit.DataStructs.cDataStructs.LongSparseIntVect:
    """
    Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect
    """

def GetHashedTopologicalTorsionFingerprint(mol: rdkit.Chem.rdchem.Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, atomInvariants: Iterable[int] | None = None, includeChirality: bool = False) -> rdkit.DataStructs.cDataStructs.LongSparseIntVect:
    """
    Returns the hashed topological-torsion fingerprint for a molecule as a LongIntSparseIntVect
    """

def GetHashedTopologicalTorsionFingerprintAsBitVect(mol: rdkit.Chem.rdchem.Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: Iterable[int] | None = None, ignoreAtoms: Iterable[int] | None = None, atomInvariants: Iterable[int] | None = None, nBitsPerEntry: int = 4, includeChirality: bool = False) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """
    Returns the topological-torsion fingerprint for a molecule as an ExplicitBitVect
    """

def GetMorganFingerprint(mol: rdkit.Chem.rdchem.Mol, radius: int, invariants: Sequence[int] | None = None, fromAtoms: Sequence[int] | None = None, useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, useCounts: bool = True, bitInfo: dict[int, tuple[tuple[int, int], ...]] | None = None, includeRedundantEnvironments: bool = False) -> rdkit.DataStructs.cDataStructs.UIntSparseIntVect:
    """Returns a Morgan fingerprint for a molecule"""

def GetHashedMorganFingerprint(mol: rdkit.Chem.rdchem.Mol, radius: int, nBits: int = 2048, invariants: Sequence[int] | None = None, fromAtoms: Sequence[int] | None = None, useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, bitInfo: dict[int, tuple[tuple[int, int], ...]] | None = None, includeRedundantEnvironments: bool = False) -> rdkit.DataStructs.cDataStructs.UIntSparseIntVect:
    """Returns a hashed Morgan fingerprint for a molecule"""

def GetMorganFingerprintAsBitVect(mol: rdkit.Chem.rdchem.Mol, radius: int, nBits: int = 2048, invariants: Sequence[int] | None = None, fromAtoms: Iterable[int] | None = None, useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, bitInfo: dict[int, tuple[tuple[int, int], ...]] | None = None, includeRedundantEnvironments: bool = False) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """Returns a Morgan fingerprint for a molecule as a bit vector"""

def GetConnectivityInvariants(mol: rdkit.Chem.rdchem.Mol, includeRingMembership: bool = True) -> list[int]:
    """Returns connectivity invariants (ECFP-like) for a molecule."""

def GetFeatureInvariants(mol: rdkit.Chem.rdchem.Mol) -> list[int]:
    """Returns feature invariants (FCFP-like) for a molecule."""

def GetUSR(mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> list[float]:
    """Returns a USR descriptor for one conformer of a molecule"""

def GetUSRDistributions(coords: Sequence[rdkit.Geometry.rdGeometry.Point3D], points: list[rdkit.Geometry.rdGeometry.Point3D] | None = None) -> list[list[float]]:
    """Returns the four USR distance distributions for a set of coordinates"""

def GetUSRDistributionsFromPoints(coords: Sequence[rdkit.Geometry.rdGeometry.Point3D], points: Sequence[rdkit.Geometry.rdGeometry.Point3D]) -> list[list[float]]:
    """
    Returns the USR distance distributions for a set of coordinates and points
    """

def GetUSRFromDistributions(distances: Sequence[Sequence[float]]) -> list[float]:
    """Returns the USR descriptor from a set of distance distributions"""

def GetUSRScore(descriptor1: Sequence[float], descriptor2: Sequence[float], weights: Sequence[float] = []) -> float:
    """Returns the USR score for two USR or USRCAT descriptors"""

def GetUSRCAT(mol: rdkit.Chem.rdchem.Mol, atomSelections: Sequence[Sequence[int]] | None = None, confId: int = -1) -> list[float]:
    """Returns a USRCAT descriptor for one conformer of a molecule"""

def CalcCrippenDescriptors(mol: rdkit.Chem.rdchem.Mol, includeHs: bool = True, force: bool = False) -> tuple[float, float]:
    """returns a 2-tuple with the Wildman-Crippen logp,mr values"""

def CalcLabuteASA(mol: rdkit.Chem.rdchem.Mol, includeHs: bool = True, force: bool = False) -> float:
    """returns the Labute ASA value for a molecule"""

def CalcTPSA(mol: rdkit.Chem.rdchem.Mol, force: bool = False, includeSandP: bool = False) -> float:
    """returns the TPSA value for a molecule"""

def CalcExactMolWt(mol: rdkit.Chem.rdchem.Mol, onlyHeavy: bool = False) -> float:
    """returns the molecule's exact molecular weight"""

def CalcMolFormula(mol: rdkit.Chem.rdchem.Mol, separateIsotopes: bool = False, abbreviateHIsotopes: bool = True) -> str:
    """returns the molecule's formula"""

def CalcNumLipinskiHBD(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of Lipinski H-bond donors for a molecule"""

def CalcNumLipinskiHBA(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of Lipinski H-bond acceptors for a molecule"""

def CalcNumHBD(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of H-bond donors for a molecule"""

def CalcNumHBA(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of H-bond acceptors for a molecule"""

class NumRotatableBondsOptions(enum.Enum):
    """
    Options for generating rotatable bonds
    NonStrict - standard loose definitions
    Strict - stricter definition excluding amides, esters, etc
    StrictLinkages - adds rotors between rotatable bonds
    Default - Current RDKit default
    """

    NonStrict = 0

    Strict = 1

    StrictLinkages = 2

    Default = -1

@overload
def CalcNumRotatableBonds(mol: rdkit.Chem.rdchem.Mol, strict: bool) -> int: ...

@overload
def CalcNumRotatableBonds(mol: rdkit.Chem.rdchem.Mol, strict: NumRotatableBondsOptions = NumRotatableBondsOptions.Default) -> int:
    """
    returns the number of rotatable bonds for a molecule.
    strict = NumRotatableBondsOptions.NonStrict - Simple rotatable bond definition.
    strict = NumRotatableBondsOptions.Strict - (default) does not count things like
             amide or ester bonds
    strict = NumRotatableBondsOptions.StrictLinkages - handles linkages between ring
       systems.
       - Single bonds between aliphatic ring Cs are always rotatable. This
         means that the central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now
         considered rotatable; it was not before
       - Heteroatoms in the linked rings no longer affect whether or not
         the linking bond is rotatable
       - the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is now
          considered non-rotatable
    """

def CalcNumRings(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of rings for a molecule"""

def CalcNumAromaticRings(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of aromatic rings for a molecule"""

def CalcNumSaturatedRings(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of saturated rings for a molecule"""

def CalcNumHeterocycles(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of heterocycles for a molecule"""

def CalcNumAromaticHeterocycles(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of aromatic heterocycles for a molecule"""

def CalcNumAromaticCarbocycles(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of aromatic carbocycles for a molecule"""

def CalcNumSaturatedHeterocycles(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of saturated heterocycles for a molecule"""

def CalcNumSaturatedCarbocycles(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of saturated carbocycles for a molecule"""

def CalcNumAliphaticRings(mol: rdkit.Chem.rdchem.Mol) -> int:
    """
    returns the number of aliphatic (containing at least one non-aromatic bond) rings for a molecule
    """

def CalcNumAliphaticHeterocycles(mol: rdkit.Chem.rdchem.Mol) -> int:
    """
    returns the number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule
    """

def CalcNumAliphaticCarbocycles(mol: rdkit.Chem.rdchem.Mol) -> int:
    """
    returns the number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule
    """

def CalcNumHeavyAtoms(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of heavy atoms for a molecule"""

def CalcNumAtoms(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the total number of atoms for a molecule"""

def CalcNumHeteroatoms(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of heteroatoms for a molecule"""

def CalcNumAmideBonds(mol: rdkit.Chem.rdchem.Mol) -> int:
    """returns the number of amide bonds in a molecule"""

def CalcFractionCSP3(mol: rdkit.Chem.rdchem.Mol) -> float:
    """returns the fraction of C atoms that are SP3 hybridized"""

def CalcChiNv(mol: rdkit.Chem.rdchem.Mol, n: int, force: bool = False) -> float:
    """
    From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    """

def CalcChi0v(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    """

def CalcChi1v(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    """

def CalcChi2v(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    """

def CalcChi3v(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    """

def CalcChi4v(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    """

def CalcChiNn(mol: rdkit.Chem.rdchem.Mol, n: int, force: bool = False) -> float:
    """
    Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.
    """

def CalcChi0n(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.
    """

def CalcChi1n(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.
    """

def CalcChi2n(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.
    """

def CalcChi3n(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.
    """

def CalcChi4n(mol: rdkit.Chem.rdchem.Mol, force: bool = False) -> float:
    """
    Similar to ChiXv, but uses uses nVal instead of valence. This makes a big difference after we get out of the first row.
    """

def CalcHallKierAlpha(mol: rdkit.Chem.rdchem.Mol, atomContribs: list[float] | None = None) -> float:
    """
    From equation (58) of Rev. Comp. Chem. vol 2, 367-422, (1991).
    NOTE: Because hybridization is used to calculate this, results may
    differ from other implementations which have different conventions for
    assigning hybridization
    """

def CalcKappa1(mol: rdkit.Chem.rdchem.Mol) -> float:
    """
    From equations (58) and (59) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    NOTE: Because hybridization is used to calculate this, results may
    differ from other implementations which have different conventions for
    assigning hybridization
    """

def CalcKappa2(mol: rdkit.Chem.rdchem.Mol) -> float:
    """
    From equations (58) and (60) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    NOTE: Because hybridization is used to calculate this, results may
    differ from other implementations which have different conventions for
    assigning hybridization
    """

def CalcKappa3(mol: rdkit.Chem.rdchem.Mol) -> float:
    """
    From equations (58), (61) and (62) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    NOTE: Because hybridization is used to calculate this, results may
    differ from other implementations which have different conventions for
    assigning hybridization
    """

def CalcPhi(mol: rdkit.Chem.rdchem.Mol) -> float:
    """
    From Quantitative Structure-Activity Relationships 8, 221-224 (1989). NOTE: Because hybridization is used to calculate this, results may differ from other implementations which have different conventions for assigning hybridization
    """

def GetMACCSKeysFingerprint(mol: rdkit.Chem.rdchem.Mol) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """Returns the MACCS keys for a molecule as an ExplicitBitVect"""

def GetAtomFeatures(mol: rdkit.Chem.rdchem.Mol, atomid: int, addchiral: bool = False) -> list[float]:
    """Returns the Atom Features vector"""

def CalcNumSpiroAtoms(mol: rdkit.Chem.rdchem.Mol, atoms: list[int] | None = None) -> int:
    """
    Returns the number of spiro atoms (atoms shared between rings that share exactly one atom)
    """

def CalcNumBridgeheadAtoms(mol: rdkit.Chem.rdchem.Mol, atoms: list[int] | None = None) -> int:
    """
    Returns the number of bridgehead atoms (atoms shared between rings that share at least two bonds)
    """

def CalcNumAtomStereoCenters(mol: rdkit.Chem.rdchem.Mol) -> int:
    """
    Returns the total number of atomic stereocenters (specified and unspecified)
    """

def CalcNumUnspecifiedAtomStereoCenters(mol: rdkit.Chem.rdchem.Mol) -> int:
    """Returns the number of unspecified atomic stereocenters"""

class Properties:
    """
    Property computation and registry system.  To compute all registered properties:
    mol = Chem.MolFromSmiles('c1ccccc1')
    properties = rdMolDescriptors.Properties()
    for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):
      print(name, value)

    To compute a subset
    properties = rdMolDescriptors.Properties(['exactmw', 'lipinskiHBA'])
    for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):
      print(name, value)
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, propNames: Sequence[str]) -> None: ...

    def GetPropertyNames(self) -> list[str]:
        """Return the property names computed by this instance"""

    def ComputeProperties(self, mol: rdkit.Chem.rdchem.Mol, annotateMol: bool = False) -> list[float]:
        """
        Return a list of computed properties, if annotateMol==True, annotate the molecule with the computed properties.
        """

    def AnnotateProperties(self, mol: rdkit.Chem.rdchem.Mol) -> None:
        """
        Annotate the molecule with the computed properties.  These properties will be available as SDData or from mol.GetProp(prop)
        """

    @staticmethod
    def GetAvailableProperties() -> list[str]:
        """Return all available property names that can be computed"""

    @staticmethod
    def GetProperty(propName: str) -> PythonPropertyFunctor:
        """Return the named property if it exists"""

    @staticmethod
    def RegisterProperty(propertyFunctor: PythonPropertyFunctor) -> int:
        """Register a new property object (not thread safe)"""

class PythonPropertyFunctor:
    def __init__(self, name: str, version: str) -> None: ...

    def __call__(self, mol: rdkit.Chem.rdchem.Mol) -> float:
        """Compute the property for the specified molecule"""

    def GetName(self) -> str:
        """Return the name of the property to calculate"""

    def GetVersion(self) -> str:
        """Return the version of the calculated property"""

class PropertyRangeQuery:
    """Property Range Query for a molecule.  Match(mol) -> true if in range"""

    def Match(self, what: rdkit.Chem.rdchem.Mol) -> bool: ...

def MakePropertyRangeQuery(name: str, min: float, max: float) -> PropertyRangeQuery:
    """
    Generates a Range property for the specified property, between min and max
    query = MakePropertyRangeQuery('exactmw', 0, 500)
    query.Match( mol )
    """

class DoubleCubicLatticeVolume:
    """Class for the Double Cubic Lattice Volume method"""

    @overload
    def __init__(self, mol: rdkit.Chem.rdchem.Mol, isProtein: bool = False, includeLigand: bool = True, probeRadius: float = 1.4, confId: int = -1) -> None: ...

    @overload
    def __init__(self, mol: rdkit.Chem.rdchem.Mol, radii: Iterable[float], isProtein: bool = False, includeLigand: bool = True, probeRadius: float = 1.4, confId: int = -1) -> None:
        """
        ARGUMENTS:
           - mol: molecule or protein under consideration
           - radii: radii for atoms of input mol (get using GetPeriodicTable or provide custom list)
           - isProtein: flag to indicate if the input is a protein (default=False, free ligand).
           - includeLigand: flag to include or exclude a bound ligand when input is a protein (default=True)
           - probeRadius: radius of the solvent probe (default=1.2)
           - confId: conformer ID to consider (default=-1)
        """

    def GetSurfaceArea(self) -> float:
        """Get the Surface Area of the Molecule or Protein"""

    def GetAtomSurfaceArea(self, atom_idx: int) -> float:
        """Get the surface area of atom with atom_idx"""

    def GetPolarSurfaceArea(self, includeSandP: bool = False, includeHs: bool = False) -> float:
        """Get the Polar Surface Area of the Molecule or Protein"""

    def GetPartialSurfaceArea(self, atomIndices: Iterable[int] | None) -> float:
        """
        Get the Partial Surface Area of the Molecule or Protein for specified subset of atoms
        """

    def GetSurfacePoints(self, allPoints: bool = False) -> dict[int, list[rdkit.Geometry.rdGeometry.Point3D]]:
        """
        Get the set of points representing the surface. If allPoints is True, returns all surface points; otherwise, returns the standard surface points.
        """

    def GetVolume(self) -> float:
        """Get the Total Volume of the Molecule or Protein"""

    def GetVDWVolume(self) -> float:
        """Get the van der Waals Volume of the Molecule or Protein"""

    def GetAtomVolume(self, atomIdx: int, solventRadius: float) -> float:
        """Get the volume atom of atom_idx with volume for specified Probe Radius"""

    def GetPolarVolume(self, includeSandP: bool = False, includeHs: bool = False) -> float:
        """Get the Polar Volume of the Molecule or Protein"""

    def GetPartialVolume(self, atomIdx: Iterable[int] | None) -> float:
        """
        Get the Partial Volume of the Molecule or Protein for specified subset of atoms
        """

    def GetCompactness(self) -> float:
        """Get the Compactness of the Protein"""

    def GetPackingDensity(self) -> float:
        """Get the PackingDensity of the Protein"""

def CalcCoulombMat(mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> list[list[float]]:
    """Returns severals Coulomb randomized matrices"""

def CalcEEMcharges(mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> list[float]:
    """Returns EEM atomic partial charges"""

def CalcWHIM(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, thresh: float = 0.001, CustomAtomProperty: str = '') -> list[float]:
    """Returns the WHIM descriptors vector"""

def CalcGETAWAY(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, precision: float = 2, CustomAtomProperty: str = '') -> list[float]:
    """Returns the GETAWAY descriptors vector"""

def CalcRDF(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, CustomAtomProperty: str = '') -> list[float]:
    """Returns radial distribution fonction descriptors (RDF)"""

def CalcMORSE(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, CustomAtomProperty: str = '') -> list[float]:
    """
    Returns Molecule Representation of Structures based on Electron diffraction descriptors
    """

def CalcAUTOCORR3D(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, CustomAtomProperty: str = '') -> list[float]:
    """Returns 3D Autocorrelation descriptors vector"""

def CalcPBF(mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> float:
    """
    Returns the PBF (plane of best fit) descriptor (https://doi.org/10.1021/ci300293f)
    """

def CalcNPR1(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float: ...

def CalcNPR2(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float: ...

def CalcPMI1(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float: ...

def CalcPMI2(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float: ...

def CalcPMI3(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float: ...

def CalcRadiusOfGyration(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float: ...

def CalcInertialShapeFactor(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float: ...

def CalcEccentricity(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float: ...

def CalcAsphericity(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float: ...

def CalcSpherocityIndex(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, force: bool = True) -> float: ...

def CalcAUTOCORR2D(mol: rdkit.Chem.rdchem.Mol, CustomAtomProperty: str = '') -> list[float]:
    """Returns 2D Autocorrelation descriptors vector"""

@overload
def BCUT2D(mol: rdkit.Chem.rdchem.Mol) -> list[float]:
    """
    Implements BCUT descriptors From J. Chem. Inf. Comput. Sci., Vol. 39, No. 1, 1999
    Diagonal elements are (currently) atomic mass, gasteiger charge,
    crippen logP and crippen MRReturns the 2D BCUT2D descriptors vector as described in
    returns [mass eigen value high, mass eigen value low,
             gasteiger charge eigenvalue high, gasteiger charge low,
             crippen lowgp  eigenvalue high, crippen lowgp  low,
             crippen mr eigenvalue high, crippen mr low]
    """

@overload
def BCUT2D(mol: rdkit.Chem.rdchem.Mol, atom_propname: str) -> tuple[float, float]:
    """
    Returns a 2D BCUT (eigen value high, eigen value low) given the
    molecule and the specified atom prop name
    atom_propname must exist on each atom and be convertible to a float
    """

@overload
def BCUT2D(mol: rdkit.Chem.rdchem.Mol, atom_props: Iterable[float]) -> tuple[float, float]:
    """
    Returns a 2D BCUT (eigen value hi, eigenvalue low) given the molecule
    and the specified atom props
     atom_props must be an iterable of floats equal in
    size to the number of atoms in mol
    """

def CalcOxidationNumbers(mol: rdkit.Chem.rdchem.Mol) -> None:
    """
    Adds the oxidation number/state to the atoms of a molecule as property OxidationNumber on each atom.  Use Pauling electronegativities.  This is experimental code, still under development.
    """
