#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
from rdkit import RDConfig
from rdkit import six
import sys,os
from rdkit import Chem
from rdkit.VLib.Filter import FilterNode

class DupeFilter(FilterNode):
  """ canonical-smiles based duplicate filter

  Assumptions:

    - inputs are molecules


  Sample Usage:
    >>> from rdkit.VLib.NodeLib.SDSupply import SDSupplyNode
    >>> fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib',\
                             'test_data','NCI_aids.10.sdf')
    >>> suppl = SDSupplyNode(fileN)
    >>> filt = DupeFilter()
    >>> filt.AddParent(suppl)
    >>> ms = [x for x in filt]
    >>> len(ms)
    10
    >>> ms[0].GetProp("_Name")
    '48'
    >>> ms[1].GetProp("_Name")
    '78'
    >>> filt.reset()
    >>> filt.next().GetProp("_Name")
    '48'


  """
  def __init__(self,**kwargs):
    FilterNode.__init__(self,func=self.filter,**kwargs)
    self._smisSeen = []
    
  def reset(self):
    FilterNode.reset(self)
    self._smisSeen = []

  def filter(self,cmpd):
    smi = Chem.MolToSmiles(cmpd)
    if smi not in self._smisSeen:
      self._smisSeen.append(smi)
      return 1
    else:
      return 0
  
#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)

  
