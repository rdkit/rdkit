#
# Copyright (C) 2001-2008 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Exposes functionality for MOE-like approximate molecular surface area
descriptors.

  The MOE-like VSA descriptors are also calculated here

"""

import bisect

import numpy

from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors, rdPartialCharges

ptable = Chem.GetPeriodicTable()
bondScaleFacts = [.1, 0, .2, .3]  # aromatic,single,double,triple


def _LabuteHelper(mol, includeHs=1, force=0):
  """ *Internal Use Only*
    helper function for LabuteASA calculation
    returns an array of atomic contributions to the ASA

  **Note:** Changes here affect the version numbers of all ASA descriptors

  """
  if not force:
    try:
      res = mol._labuteContribs
    except AttributeError:
      pass
    else:
      if res:
        return res
  tpl = rdMolDescriptors._CalcLabuteASAContribs(mol, includeHs)
  ats, hs = tpl
  Vi = [hs] + list(ats)
  mol._labuteContribs = Vi
  return Vi


def _pyLabuteHelper(mol, includeHs=1, force=0):
  """ *Internal Use Only*
    helper function for LabuteASA calculation
    returns an array of atomic contributions to the ASA

  **Note:** Changes here affect the version numbers of all ASA descriptors

  """
  import math
  if not force:
    try:
      res = mol._labuteContribs
    except AttributeError:
      pass
    else:
      if res.all():
        return res

  nAts = mol.GetNumAtoms()
  Vi = numpy.zeros(nAts + 1, 'd')
  rads = numpy.zeros(nAts + 1, 'd')

  # 0 contains the H information
  rads[0] = ptable.GetRb0(1)
  for i in range(nAts):
    rads[i + 1] = ptable.GetRb0(mol.GetAtomWithIdx(i).GetAtomicNum())

  # start with explicit bonds
  for bond in mol.GetBonds():
    idx1 = bond.GetBeginAtomIdx() + 1
    idx2 = bond.GetEndAtomIdx() + 1
    Ri = rads[idx1]
    Rj = rads[idx2]

    if not bond.GetIsAromatic():
      bij = Ri + Rj - bondScaleFacts[bond.GetBondType()]
    else:
      bij = Ri + Rj - bondScaleFacts[0]
    dij = min(max(abs(Ri - Rj), bij), Ri + Rj)
    Vi[idx1] += Rj * Rj - (Ri - dij)**2 / dij
    Vi[idx2] += Ri * Ri - (Rj - dij)**2 / dij

  # add in hydrogens
  if includeHs:
    j = 0
    Rj = rads[j]
    for i in range(1, nAts + 1):
      Ri = rads[i]
      bij = Ri + Rj
      dij = min(max(abs(Ri - Rj), bij), Ri + Rj)
      Vi[i] += Rj * Rj - (Ri - dij)**2 / dij
      Vi[j] += Ri * Ri - (Rj - dij)**2 / dij

  for i in range(nAts + 1):
    Ri = rads[i]
    Vi[i] = 4 * math.pi * Ri**2 - math.pi * Ri * Vi[i]

  mol._labuteContribs = Vi
  return Vi


# def SMR_VSA(mol,bins=[0.11,0.26,0.35,0.39,0.44,0.485,0.56]):
# original default bins from assuming Labute values are logs
# mrBins=[1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63]
mrBins = [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0]


def pySMR_VSA_(mol, bins=None, force=1):
  """ *Internal Use Only*
  """
  if not force:
    try:
      res = mol._smrVSA
    except AttributeError:
      pass
    else:
      if res.all():
        return res

  if bins is None:
    bins = mrBins
  Crippen._Init()
  propContribs = Crippen._GetAtomContribs(mol, force=force)
  volContribs = _LabuteHelper(mol)

  ans = numpy.zeros(len(bins) + 1, 'd')
  for i in range(len(propContribs)):
    prop = propContribs[i]
    vol = volContribs[i + 1]
    if prop is not None:
      bin_ = bisect.bisect_right(bins, prop[1])
      ans[bin_] += vol

  mol._smrVSA = ans
  return ans


SMR_VSA_ = rdMolDescriptors.SMR_VSA_

#
# Original bins (from Labute paper) are:
#  [-0.4,-0.2,0,0.1,0.15,0.2,0.25,0.3,0.4]
#
logpBins = [-0.4, -0.2, 0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]


def pySlogP_VSA_(mol, bins=None, force=1):
  """ *Internal Use Only*
  """
  if not force:
    try:
      res = mol._slogpVSA
    except AttributeError:
      pass
    else:
      if res.all():
        return res

  if bins is None:
    bins = logpBins
  Crippen._Init()
  propContribs = Crippen._GetAtomContribs(mol, force=force)
  volContribs = _LabuteHelper(mol)

  ans = numpy.zeros(len(bins) + 1, 'd')
  for i in range(len(propContribs)):
    prop = propContribs[i]
    vol = volContribs[i + 1]
    if prop is not None:
      bin_ = bisect.bisect_right(bins, prop[0])
      ans[bin_] += vol

  mol._slogpVSA = ans
  return ans


SlogP_VSA_ = rdMolDescriptors.SlogP_VSA_

chgBins = [-.3, -.25, -.20, -.15, -.10, -.05, 0, .05, .10, .15, .20, .25, .30]


def pyPEOE_VSA_(mol, bins=None, force=1):
  """ *Internal Use Only*
  """
  if not force:
    try:
      res = mol._peoeVSA
    except AttributeError:
      pass
    else:
      if res.all():
        return res
  if bins is None:
    bins = chgBins
  Crippen._Init()
  # print('\ts:',repr(mol.GetMol()))
  # print('\t\t:',len(mol.GetAtoms()))
  rdPartialCharges.ComputeGasteigerCharges(mol)

  # propContribs = [float(x.GetProp('_GasteigerCharge'))  for x in mol.GetAtoms()]
  propContribs = []
  for at in mol.GetAtoms():
    p = at.GetProp('_GasteigerCharge')
    try:
      v = float(p)
    except ValueError:
      v = 0.0
    propContribs.append(v)
  volContribs = _LabuteHelper(mol)

  ans = numpy.zeros(len(bins) + 1, 'd')
  for i in range(len(propContribs)):
    prop = propContribs[i]
    vol = volContribs[i + 1]
    if prop is not None:
      bin_ = bisect.bisect_right(bins, prop)
      ans[bin_] += vol

  mol._peoeVSA = ans
  return ans


PEOE_VSA_ = rdMolDescriptors.PEOE_VSA_


# -------------------------------------------------
# install the various VSA descriptors in the namespace
def _InstallDescriptors():
  for i in range(len(mrBins)):
    fn = lambda x, y=i: SMR_VSA_(x, force=0)[y]
    if i > 0:
      fn.__doc__ = "MOE MR VSA Descriptor %d (% 4.2f <= x < % 4.2f)" % (i + 1, mrBins[i - 1],
                                                                        mrBins[i])
    else:
      fn.__doc__ = "MOE MR VSA Descriptor %d (-inf < x < % 4.2f)" % (i + 1, mrBins[i])
    name = "SMR_VSA%d" % (i + 1)
    fn.version = "1.0.1"
    globals()[name] = fn
  i += 1
  fn = lambda x, y=i: SMR_VSA_(x, force=0)[y]
  fn.__doc__ = "MOE MR VSA Descriptor %d (% 4.2f <= x < inf)" % (i + 1, mrBins[i - 1])
  fn.version = "1.0.1"
  name = "SMR_VSA%d" % (i + 1)
  globals()[name] = fn

  for i in range(len(logpBins)):
    fn = lambda x, y=i: SlogP_VSA_(x, force=0)[y]
    if i > 0:
      fn.__doc__ = "MOE logP VSA Descriptor %d (% 4.2f <= x < % 4.2f)" % (i + 1, logpBins[i - 1],
                                                                          logpBins[i])
    else:
      fn.__doc__ = "MOE logP VSA Descriptor %d (-inf < x < % 4.2f)" % (i + 1, logpBins[i])
    name = "SlogP_VSA%d" % (i + 1)
    fn.version = "1.0.1"
    globals()[name] = fn
  i += 1
  fn = lambda x, y=i: SlogP_VSA_(x, force=0)[y]
  fn.__doc__ = "MOE logP VSA Descriptor %d (% 4.2f <= x < inf)" % (i + 1, logpBins[i - 1])
  fn.version = "1.0.1"
  name = "SlogP_VSA%d" % (i + 1)
  globals()[name] = fn

  for i in range(len(chgBins)):
    fn = lambda x, y=i: PEOE_VSA_(x, force=0)[y]
    if i > 0:
      fn.__doc__ = "MOE Charge VSA Descriptor %d (% 4.2f <= x < % 4.2f)" % (i + 1, chgBins[i - 1],
                                                                            chgBins[i])
    else:
      fn.__doc__ = "MOE Charge VSA Descriptor %d (-inf < x < % 4.2f)" % (i + 1, chgBins[i])
    name = "PEOE_VSA%d" % (i + 1)
    fn.version = "1.0.1"
    globals()[name] = fn
  i += 1
  fn = lambda x, y=i: PEOE_VSA_(x, force=0)[y]
  fn.version = "1.0.1"
  fn.__doc__ = "MOE Charge VSA Descriptor %d (% 4.2f <= x < inf)" % (i + 1, chgBins[i - 1])
  name = "PEOE_VSA%d" % (i + 1)
  globals()[name] = fn
  fn = None
  # Change log for the MOE-type descriptors:
  #  version 1.0.1: optimizations, values unaffected


_InstallDescriptors()


def pyLabuteASA(mol, includeHs=1):
  """ calculates Labute's Approximate Surface Area (ASA from MOE)

    Definition from P. Labute's article in the Journal of the Chemical Computing Group
    and J. Mol. Graph. Mod.  _18_ 464-477 (2000)

  """
  Vi = _LabuteHelper(mol, includeHs=includeHs)
  return sum(Vi)


pyLabuteASA.version = "1.0.1"
# Change log for LabuteASA:
#  version 1.0.1: optimizations, values unaffected
LabuteASA = lambda *x, **y: rdMolDescriptors.CalcLabuteASA(*x, **y)
LabuteASA.version = rdMolDescriptors._CalcLabuteASA_version


def _pyTPSAContribs(mol, verbose=False):
  """ DEPRECATED: this has been reimplmented in C++
  calculates atomic contributions to a molecules TPSA

   Algorithm described in:
    P. Ertl, B. Rohde, P. Selzer
     Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-based
     Contributions and Its Application to the Prediction of Drug Transport
     Properties, J.Med.Chem. 43, 3714-3717, 2000

   Implementation based on the Daylight contrib program tpsa.c

   NOTE: The JMC paper describing the TPSA algorithm includes
   contributions from sulfur and phosphorus, however according to
   Peter Ertl (personal communication, 2010) the correlation of TPSA
   with various ADME properties is better if only contributions from
   oxygen and nitrogen are used. This matches the daylight contrib
   implementation.

  """
  res = [0] * mol.GetNumAtoms()
  for i in range(mol.GetNumAtoms()):
    atom = mol.GetAtomWithIdx(i)
    atNum = atom.GetAtomicNum()
    if atNum in [7, 8]:
      nHs = atom.GetTotalNumHs()
      chg = atom.GetFormalCharge()
      in3Ring = atom.IsInRingSize(3)

      bonds = atom.GetBonds()
      numNeighbors = atom.GetDegree()
      nSing = 0
      nDoub = 0
      nTrip = 0
      nArom = 0
      for bond in bonds:
        otherAt = bond.GetOtherAtom(atom)
        if otherAt.GetAtomicNum() != 1:
          if bond.GetIsAromatic():
            nArom += 1
          else:
            order = bond.GetBondType()
            if order == Chem.BondType.SINGLE:
              nSing += 1
            elif order == Chem.BondType.DOUBLE:
              nDoub += 1
            elif order == Chem.BondType.TRIPLE:
              nTrip += 1
        else:
          numNeighbors -= 1
          nHs += 1
      tmp = -1
      if atNum == 7:
        if numNeighbors == 1:
          if nHs == 0 and nTrip == 1 and chg == 0:
            tmp = 23.79
          elif nHs == 1 and nDoub == 1 and chg == 0:
            tmp = 23.85
          elif nHs == 2 and nSing == 1 and chg == 0:
            tmp = 26.02
          elif nHs == 2 and nDoub == 1 and chg == 1:
            tmp = 25.59
          elif nHs == 3 and nSing == 1 and chg == 1:
            tmp = 27.64
        elif numNeighbors == 2:
          if nHs == 0 and nSing == 1 and nDoub == 1 and chg == 0:
            tmp = 12.36
          elif nHs == 0 and nTrip == 1 and nDoub == 1 and chg == 0:
            tmp = 13.60
          elif nHs == 1 and nSing == 2 and chg == 0:
            if not in3Ring:
              tmp = 12.03
            else:
              tmp = 21.94
          elif nHs == 0 and nTrip == 1 and nSing == 1 and chg == 1:
            tmp = 4.36
          elif nHs == 1 and nDoub == 1 and nSing == 1 and chg == 1:
            tmp = 13.97
          elif nHs == 2 and nSing == 2 and chg == 1:
            tmp = 16.61
          elif nHs == 0 and nArom == 2 and chg == 0:
            tmp = 12.89
          elif nHs == 1 and nArom == 2 and chg == 0:
            tmp = 15.79
          elif nHs == 1 and nArom == 2 and chg == 1:
            tmp = 14.14
        elif numNeighbors == 3:
          if nHs == 0 and nSing == 3 and chg == 0:
            if not in3Ring:
              tmp = 3.24
            else:
              tmp = 3.01
          elif nHs == 0 and nSing == 1 and nDoub == 2 and chg == 0:
            tmp = 11.68
          elif nHs == 0 and nSing == 2 and nDoub == 1 and chg == 1:
            tmp = 3.01
          elif nHs == 1 and nSing == 3 and chg == 1:
            tmp = 4.44
          elif nHs == 0 and nArom == 3 and chg == 0:
            tmp = 4.41
          elif nHs == 0 and nSing == 1 and nArom == 2 and chg == 0:
            tmp = 4.93
          elif nHs == 0 and nDoub == 1 and nArom == 2 and chg == 0:
            tmp = 8.39
          elif nHs == 0 and nArom == 3 and chg == 1:
            tmp = 4.10
          elif nHs == 0 and nSing == 1 and nArom == 2 and chg == 1:
            tmp = 3.88
        elif numNeighbors == 4:
          if nHs == 0 and nSing == 4 and chg == 1:
            tmp = 0.00
        if tmp < 0.0:
          tmp = 30.5 - numNeighbors * 8.2 + nHs * 1.5
          if tmp < 0.0:
            tmp = 0.0
      elif atNum == 8:
        if numNeighbors == 1:
          if nHs == 0 and nDoub == 1 and chg == 0:
            tmp = 17.07
          elif nHs == 1 and nSing == 1 and chg == 0:
            tmp = 20.23
          elif nHs == 0 and nSing == 1 and chg == -1:
            tmp = 23.06
        elif numNeighbors == 2:
          if nHs == 0 and nSing == 2 and chg == 0:
            if not in3Ring:
              tmp = 9.23
            else:
              tmp = 12.53
          elif nHs == 0 and nArom == 2 and chg == 0:
            tmp = 13.14

        if tmp < 0.0:
          tmp = 28.5 - numNeighbors * 8.6 + nHs * 1.5
          if tmp < 0.0:
            tmp = 0.0
      if verbose:
        print('\t', atom.GetIdx(), atom.GetSymbol(), atNum, nHs, nSing, nDoub, nTrip, nArom, chg,
              tmp)

      res[atom.GetIdx()] = tmp
  return res


def _pyTPSA(mol, verbose=False):
  """ DEPRECATED: this has been reimplmented in C++
   calculates the polar surface area of a molecule based upon fragments

   Algorithm in:
    P. Ertl, B. Rohde, P. Selzer
     Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-based
     Contributions and Its Application to the Prediction of Drug Transport
     Properties, J.Med.Chem. 43, 3714-3717, 2000

   Implementation based on the Daylight contrib program tpsa.c
  """
  contribs = _pyTPSAContribs(mol, verbose=verbose)
  res = 0.0
  for contrib in contribs:
    res += contrib
  return res


_pyTPSA.version = "1.0.1"

TPSA = lambda *x, **y: rdMolDescriptors.CalcTPSA(*x, **y)
TPSA.version = rdMolDescriptors._CalcTPSA_version

if __name__ == '__main__':  # pragma: nocover
  smis = ['C', 'CC', 'CCC', 'CCCC', 'CO', 'CCO', 'COC']
  smis = ['C(=O)O', 'c1ccccc1']
  for smi in smis:
    m = Chem.MolFromSmiles(smi)
    # print(smi, LabuteASA(m))
    print('-----------\n', smi)
    # print('M:',['% 4.2f'%x for x in SMR_VSA_(m)])
    # print('L:',['% 4.2f'%x for x in SlogP_VSA_(m)])
    print('P:', ['% 4.2f' % x for x in PEOE_VSA_(m)])
    print('P:', ['% 4.2f' % x for x in PEOE_VSA_(m)])
    print()
