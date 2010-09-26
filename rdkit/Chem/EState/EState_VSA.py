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
import numpy
from rdkit.Chem.EState.EState import EStateIndices as EStateIndices_
from rdkit.Chem.MolSurf import _LabuteHelper as VSAContribs_
import bisect

"""

These default VSA bins were chosen using the PP3K solubility data
set.  An arbitrary number of bins were selected and the
boundaries were selected to give an approximately equal number of
atoms per bin

"""
vsaBins=[4.78,5.00,5.410,5.740,6.00,6.07,6.45,7.00,11.0]
def VSA_EState_(mol,bins=None,force=1):
  """ *Internal Use Only*
  """
  if not force and hasattr(mol,'_vsaEState'):
    return mol._vsaEState

  if bins is None: bins = estateBins
  propContribs = EStateIndices_(mol,force=force)
  volContribs = VSAContribs_(mol)

  ans = numpy.zeros(len(bins)+1,numpy.float)
  for i,prop in enumerate(propContribs):
    if prop is not None:
      bin = bisect.bisect_right(bins,volContribs[i+1])
      ans[bin] += prop
  mol._vsaEState=ans
  return ans    


"""

These default EState bins were chosen using the PP3K solubility data
set.  An arbitrary number of bins (10) were selected and the
boundaries were selected to give an approximately equal number of
atoms per bin

"""
estateBins=[-0.390,0.290,0.717,1.165,1.540,1.807,2.05,4.69,9.17,15.0]
def EState_VSA_(mol,bins=None,force=1):
  """ *Internal Use Only*
  """
  if not force and hasattr(mol,'_eStateVSA'):
    return mol._eStateVSA

  if bins is None: bins = estateBins
  propContribs = EStateIndices_(mol,force=force)
  volContribs = VSAContribs_(mol)

  ans = numpy.zeros(len(bins)+1,numpy.float)
  for i,prop in enumerate(propContribs):
    if prop is not None:
      bin = bisect.bisect_right(bins,prop)
      ans[bin] += volContribs[i+1]
  mol._eStateVSA=ans
  return ans    
def _InstallDescriptors():
  for i in range(len(vsaBins)):
    fn = lambda x,y=i:VSA_EState_(x,force=0)[y]
    if i > 0:
      fn.__doc__="VSA EState Descriptor %d (% 4.2f <= x < % 4.2f)"%(i+1,vsaBins[i-1],vsaBins[i])
    else:
      fn.__doc__="VSA EState Descriptor %d (-inf < x < % 4.2f)"%(i+1,vsaBins[i])
    name="VSA_EState%d"%(i+1)
    fn.version="1.0.0"
    globals()[name]=fn
  i+=1
  fn = lambda x,y=i:VSA_EState_(x,force=0)[y]
  fn.__doc__="VSA EState Descriptor %d (% 4.2f <= x < inf)"%(i+1,vsaBins[i-1])
  name="VSA_EState%d"%(i+1)
  fn.version="1.0.0"
  globals()[name]=fn
  fn=None

  for i in range(len(estateBins)):
    fn = lambda x,y=i:EState_VSA_(x,force=0)[y]
    if i > 0:
      fn.__doc__="EState VSA Descriptor %d (% 4.2f <= x < % 4.2f)"%(i+1,estateBins[i-1],estateBins[i])
    else:
      fn.__doc__="EState VSA Descriptor %d (-inf < x < % 4.2f)"%(i+1,estateBins[i])
    name="EState_VSA%d"%(i+1)
    fn.version="1.0.1"
    globals()[name]=fn
  i+=1
  fn = lambda x,y=i:EState_VSA_(x,force=0)[y]
  fn.__doc__="EState VSA Descriptor %d (% 4.2f <= x < inf)"%(i+1,estateBins[i-1])
  name="EState_VSA%d"%(i+1)
  fn.version="1.0.1"
  globals()[name]=fn
  fn=None
# Change log for EState_VSA descriptors:
#  version 1.0.1: optimizations, values unaffected
_InstallDescriptors()
