#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
import sys,os.path
from rdkit import six
from rdkit import RDConfig
from rdkit.VLib.Supply import SupplyNode
from rdkit import Chem

class SDSupplyNode(SupplyNode):
  """ SD supplier

  Sample Usage:
    >>> fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib',\
                             'test_data','NCI_aids.10.sdf')
    >>> suppl = SDSupplyNode(fileN)
    >>> ms = [x for x in suppl]
    >>> len(ms)
    10
    >>> ms[0].GetProp("_Name")
    '48'
    >>> ms[1].GetProp("_Name")
    '78'
    >>> suppl.reset()
    >>> suppl.next().GetProp("_Name")
    '48'
    >>> suppl.next().GetProp("_Name")
    '78'

  
  """
  def __init__(self,fileName,**kwargs):
    SupplyNode.__init__(self,**kwargs)
    self._fileName = fileName
    self._supplier = Chem.SDMolSupplier(self._fileName)

  def reset(self):
    SupplyNode.reset(self)
    self._supplier.reset()
  def next(self):
    """

    """
    return next(self._supplier)

if six.PY3:
    SDSupplyNode.__next__ = SDSupplyNode.next

  
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


  
