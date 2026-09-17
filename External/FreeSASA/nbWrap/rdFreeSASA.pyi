"""Module containing rdFreeSASA classes and functions."""

from collections.abc import Sequence
import enum
from typing import overload

import rdkit.Chem.rdchem


class SASAAlgorithm(enum.Enum):
    LeeRichards = 0

    ShrakeRupley = 1

LeeRichards: SASAAlgorithm = SASAAlgorithm.LeeRichards

ShrakeRupley: SASAAlgorithm = SASAAlgorithm.ShrakeRupley

class SASAClassifier(enum.Enum):
    Protor = 0

    NACCESS = 1

    OONS = 2

Protor: SASAClassifier = SASAClassifier.Protor

NACCESS: SASAClassifier = SASAClassifier.NACCESS

OONS: SASAClassifier = SASAClassifier.OONS

class SASAClass(enum.Enum):
    Unclassified = 0

    APolar = 1

    Polar = 2

Unclassified: SASAClass = SASAClass.Unclassified

APolar: SASAClass = SASAClass.APolar

Polar: SASAClass = SASAClass.Polar

class SASAOpts:
    @overload
    def __init__(self) -> None:
        """Constructor takes no arguments"""

    @overload
    def __init__(self, alg: SASAAlgorithm, cls: SASAClassifier) -> None: ...

    @overload
    def __init__(self, alg: SASAAlgorithm, cls: SASAClassifier, pr: float) -> None: ...

    @property
    def algorithm(self) -> SASAAlgorithm: ...

    @algorithm.setter
    def algorithm(self, arg: SASAAlgorithm, /) -> None: ...

    @property
    def classifier(self) -> SASAClassifier: ...

    @classifier.setter
    def classifier(self, arg: SASAClassifier, /) -> None: ...

    @property
    def probeRadius(self) -> float: ...

    @probeRadius.setter
    def probeRadius(self, arg: float, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

def classifyAtoms(mol: rdkit.Chem.rdchem.Mol, options: SASAOpts = ...) -> list:
    """
    Classify the atoms in the molecule returning their radii if possible.
    ARGUMENTS:
       - mol: molecule to classify
       - options: FreeSASA options class specifying the classification method.
                   Current classifiers are Protor, NACCESS and OONS
                   classification is stored as atom property 'SASAClass' for the integer value
                    and 'SASAClassName' for the string name of the class, Polar, APolar...

    RETURNS:
      list of radii where radii[atom.GetIdx()] is the radii of the atom.
      If classification fails, NONE is returned
    """

def CalcSASA(mol: rdkit.Chem.rdchem.Mol, radii: Sequence[float], confIdx: int = -1, query: rdkit.Chem.rdchem.Atom | None = None, opts: SASAOpts = ...) -> float:
    """
    Compute the Solvent Accessible Surface Area using the FreeSASA library
    ARGUMENTS:
      - mol: The molecule to compute.
      - radii:  A list of atom raddii where radii[atom.GetIdx()] is the radius of the atom
                These can be passed in or calculated with classifyAtoms for some proteins
      - confIdx: Specify the conformer to use for the 3D geometry  [default -1]
      - query: Pass along a query atom to compute the SASA for a subset of atoms.
               precanned query atoms can be made with MakeFreeSasaPolarAtomQuery and
               MakeFreeSasaAPolarAtomQuery for classified polar and apolar atoms respectively.
      - opts: a SASAOpts class specifying the algorithm to use

    RETURNS:
    The computed solvent accessible surface area.
    """

def MakeFreeSasaAPolarAtomQuery() -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns an APolar atom query for use with CalcSASA.  An apolar atom has the SASAClass
    and SASAClassName set to the APOLAR class.  (see classifyAtoms)
    """

def MakeFreeSasaPolarAtomQuery() -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a polar atom query for use with CalcSASA.  An polar atom has the SASAClass
    and SASAClassName set to the POLAR class.  (see classifyAtoms)
    """
