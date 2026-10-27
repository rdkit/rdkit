# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Calculation of topological/topochemical descriptors.

The calculations are done in C++ (rdkit.Chem.rdMolDescriptors); the
functions here provide the historical python API.

"""

import numpy

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ptable = Chem.GetPeriodicTable()

# WARNING: this data should probably go somewhere else...
hallKierAlphas = {
  'Br': [None, None, 0.48],
  'C': [-0.22, -0.13, 0.0],
  'Cl': [None, None, 0.29],
  'F': [None, None, -0.07],
  'H': [0.0, 0.0, 0.0],
  'I': [None, None, 0.73],
  'N': [-0.29, -0.2, -0.04],
  'O': [None, -0.2, -0.04],
  'P': [None, 0.3, 0.43],
  'S': [None, 0.22, 0.35]
}


def _asDistMat(dMat):
  """ *Internal Use Only*

  a C-contiguous float64 copy (if needed) of a user-supplied distance matrix

  """
  return numpy.ascontiguousarray(dMat, dtype=numpy.float64)


def Ipc(mol, avg=False, dMat=None, forceDMat=False):
  """This returns the information content of the coefficients of the characteristic
    polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule.

    'avg = True' returns the information content divided by the total population.

    From Eq 6 of D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67, 4517-4533 (1977)

  """
  if forceDMat or dMat is None:
    return rdMolDescriptors.CalcIpc(mol, avg=bool(avg))
  return rdMolDescriptors.CalcIpc(mol, avg=bool(avg), dMat=_asDistMat(dMat))


Ipc.version = rdMolDescriptors._CalcIpc_version


def AvgIpc(mol, dMat=None, forceDMat=False):
  """This returns the average information content of the coefficients of the characteristic
    polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule.

    From Eq 7 of D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67, 4517-4533 (1977)

  """
  return Ipc(mol, avg=True, dMat=dMat, forceDMat=forceDMat)


AvgIpc.version = rdMolDescriptors._CalcAvgIpc_version

HallKierAlpha = lambda x: rdMolDescriptors.CalcHallKierAlpha(x)
HallKierAlpha.version = rdMolDescriptors._CalcHallKierAlpha_version
Kappa1 = lambda x: rdMolDescriptors.CalcKappa1(x)
Kappa1.version = rdMolDescriptors._CalcKappa1_version
Kappa2 = lambda x: rdMolDescriptors.CalcKappa2(x)
Kappa2.version = rdMolDescriptors._CalcKappa2_version
Kappa3 = lambda x: rdMolDescriptors.CalcKappa3(x)
Kappa3.version = rdMolDescriptors._CalcKappa3_version


def Chi0(mol):
  """ From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  return rdMolDescriptors.CalcChi0(mol)


Chi0.version = rdMolDescriptors._CalcChi0_version


def Chi1(mol):
  """ From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  return rdMolDescriptors.CalcChi1(mol)


Chi1.version = rdMolDescriptors._CalcChi1_version

Chi0v = lambda x: rdMolDescriptors.CalcChi0v(x)
Chi0v.version = rdMolDescriptors._CalcChi0v_version
Chi1v = lambda x: rdMolDescriptors.CalcChi1v(x)
Chi1v.version = rdMolDescriptors._CalcChi1v_version
Chi2v = lambda x: rdMolDescriptors.CalcChi2v(x)
Chi2v.version = rdMolDescriptors._CalcChi2v_version
Chi3v = lambda x: rdMolDescriptors.CalcChi3v(x)
Chi3v.version = rdMolDescriptors._CalcChi3v_version
Chi4v = lambda x: rdMolDescriptors.CalcChi4v(x)
Chi4v.version = rdMolDescriptors._CalcChi4v_version
ChiNv_ = lambda x, y: rdMolDescriptors.CalcChiNv(x, y)
ChiNv_.version = rdMolDescriptors._CalcChiNv_version

Chi0n = lambda x: rdMolDescriptors.CalcChi0n(x)
Chi0n.version = rdMolDescriptors._CalcChi0n_version
Chi1n = lambda x: rdMolDescriptors.CalcChi1n(x)
Chi1n.version = rdMolDescriptors._CalcChi1n_version
Chi2n = lambda x: rdMolDescriptors.CalcChi2n(x)
Chi2n.version = rdMolDescriptors._CalcChi2n_version
Chi3n = lambda x: rdMolDescriptors.CalcChi3n(x)
Chi3n.version = rdMolDescriptors._CalcChi3n_version
Chi4n = lambda x: rdMolDescriptors.CalcChi4n(x)
Chi4n.version = rdMolDescriptors._CalcChi4n_version
ChiNn_ = lambda x, y: rdMolDescriptors.CalcChiNn(x, y)
ChiNn_.version = rdMolDescriptors._CalcChiNn_version


def BalabanJ(mol, dMat=None, forceDMat=False):
  """ Calculate Balaban's J value for a molecule

  **Arguments**

    - mol: a molecule

    - dMat: (optional) a distance/adjacency matrix for the molecule, if this
      is not provide, one will be calculated

    - forceDMat: (optional) if this is set, the distance/adjacency matrix
      will be recalculated regardless of whether or not _dMat_ is provided
      or the molecule already has one

  **Returns**

    - a float containing the J value

  We follow the notation of Balaban's paper:
    Chem. Phys. Lett. vol 89, 399-404, (1982)

  """
  if forceDMat or dMat is None:
    return rdMolDescriptors.CalcBalabanJ(mol)
  return rdMolDescriptors.CalcBalabanJ(mol, dMat=_asDistMat(dMat))


BalabanJ.version = rdMolDescriptors._CalcBalabanJ_version


def BertzCT(mol, cutoff=100, dMat=None, forceDMat=1):
  """ A topological index meant to quantify "complexity" of molecules.

     Consists of a sum of two terms, one representing the complexity
     of the bonding, the other representing the complexity of the
     distribution of heteroatoms.

     From S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981)

     "cutoff" is an integer value used to limit the computational
     expense.  A cutoff value tells the program to consider vertices
     topologically identical if their distance vectors (sets of
     distances to all other vertices) are equal out to the "cutoff"th
     nearest-neighbor.

     **NOTE**  The original implementation had the following comment:
         > this implementation treats aromatic rings as the
         > corresponding Kekule structure with alternating bonds,
         > for purposes of counting "connections".
       Upon further thought, this is the WRONG thing to do.  It
        results in the possibility of a molecule giving two different
        CT values depending on the kekulization.  For example, in the
        old implementation, these two SMILES:
           CC2=CN=C1C3=C(C(C)=C(C=N3)C)C=CC1=C2C
           CC3=CN=C2C1=NC=C(C)C(C)=C1C=CC2=C3C
        which correspond to differentk kekule forms, yield different
        values.
       The new implementation uses consistent (aromatic) bond orders
        for aromatic bonds.

       THIS MEANS THAT THIS IMPLEMENTATION IS NOT BACKWARDS COMPATIBLE.

       Any molecule containing aromatic rings will yield different
       values with this implementation.  The new behavior is the correct
       one, so we're going to live with the breakage.

     "dMat" is an optional distance matrix used to assign the symmetry
     classes. It is only used if "forceDMat" is not set; by default the
     bond-order weighted (Balaban) distance matrix is used.

  """
  if forceDMat:
    return rdMolDescriptors.CalcBertzCT(mol, int(cutoff))
  if dMat is None:
    # without forceDMat the symmetry classes have always come from the
    # topological distance matrix
    dMat = Chem.GetDistanceMatrix(mol, useBO=False, useAtomWts=False, force=True)
  return rdMolDescriptors.CalcBertzCT(mol, int(cutoff), dMat=_asDistMat(dMat))


BertzCT.version = rdMolDescriptors._CalcBertzCT_version
# Recent Revisions:
#  1.0.0 -> 2.0.0:
#    - force distance matrix updates properly (Fixed as part of Issue 125)
#    - handle single-atom fragments (Issue 136)
