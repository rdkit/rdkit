# $Id$
#
# Copyright (C)2003-2010 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Hybrid EState-VSA descriptors (like the MOE VSA descriptors)

"""

import bisect

import numpy

from rdkit.Chem.EState.EState import EStateIndices as EStateIndices_
from rdkit.Chem.MolSurf import _LabuteHelper as VSAContribs_

"""

These default VSA bins were chosen using the PP3K solubility data
set.  An arbitrary number of bins were selected and the
boundaries were selected to give an approximately equal number of
atoms per bin

"""

vsaBins = [4.78, 5.00, 5.410, 5.740, 6.00, 6.07, 6.45, 7.00, 11.0]


def VSA_EState_(mol, bins=None, force=1):
  """ *Internal Use Only*
  """
  if not force and hasattr(mol, '_vsaEState'):
    return mol._vsaEState

  if bins is None:
    bins = vsaBins
  propContribs = EStateIndices_(mol, force=force)
  volContribs = VSAContribs_(mol)

  ans = numpy.zeros(len(bins) + 1, dtype=numpy.float64)
  for i, prop in enumerate(propContribs):
    if prop is not None:
      nbin = bisect.bisect_right(bins, volContribs[i + 1])
      ans[nbin] += prop
  mol._vsaEState = ans
  return ans


"""

These default EState bins were chosen using the PP3K solubility data
set.  An arbitrary number of bins (10) were selected and the
boundaries were selected to give an approximately equal number of
atoms per bin

"""
estateBins = [-0.390, 0.290, 0.717, 1.165, 1.540, 1.807, 2.05, 4.69, 9.17, 15.0]


def EState_VSA_(mol, bins=None, force=1):
  """ *Internal Use Only*
  """
  if not force and hasattr(mol, '_eStateVSA'):
    return mol._eStateVSA

  if bins is None:
    bins = estateBins
  propContribs = EStateIndices_(mol, force=force)
  volContribs = VSAContribs_(mol)

  ans = numpy.zeros(len(bins) + 1, dtype=numpy.float64)
  for i, prop in enumerate(propContribs):
    if prop is not None:
      nbin = bisect.bisect_right(bins, prop)
      ans[nbin] += volContribs[i + 1]
  mol._eStateVSA = ans
  return ans


def _descriptorDocstring(name, nbin, bins):
  """ Create a docstring for the descriptor name """
  if nbin == 0:
    interval = "-inf < x <  {0:.2f}".format(bins[nbin])
  elif nbin < len(bins):
    interval = " {0:.2f} <= x <  {1:.2f}".format(bins[nbin - 1], bins[nbin])
  else:
    interval = " {0:.2f} <= x < inf".format(bins[nbin - 1])
  return '{0} Descriptor {1} ({2})'.format(name, nbin + 1, interval)


def _descriptor_VSA_EState(nbin):

  def VSA_EState_bin(mol):
    return VSA_EState_(mol, force=False)[nbin]

  name = "VSA_EState{0}".format(nbin + 1)
  fn = VSA_EState_bin
  fn.__doc__ = _descriptorDocstring('VSA EState', nbin, vsaBins)
  fn.version = '1.0.0'
  return name, fn


def _descriptor_EState_VSA(nbin):

  def EState_VSA_bin(mol):
    return EState_VSA_(mol, force=False)[nbin]

  name = "EState_VSA{0}".format(nbin + 1)
  fn = EState_VSA_bin
  fn.__name__ = name
  if hasattr(fn, '__qualname__'):
    fn.__qualname__ = name
  fn.__doc__ = _descriptorDocstring('EState VSA', nbin, estateBins)
  fn.version = '1.0.1'
  return name, fn


def _InstallDescriptors():
  for nbin in range(len(vsaBins) + 1):
    name, fn = _descriptor_VSA_EState(nbin)
    globals()[name] = fn

  for nbin in range(len(estateBins) + 1):
    name, fn = _descriptor_EState_VSA(nbin)
    globals()[name] = fn


# Change log for EState_VSA descriptors:
#  version 1.0.1: optimizations, values unaffected
_InstallDescriptors()
