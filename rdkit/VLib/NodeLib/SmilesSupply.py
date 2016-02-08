#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
import sys,os.path
from rdkit import RDConfig
from rdkit.VLib.Supply import SupplyNode
from rdkit import Chem
from rdkit import six


class SmilesSupplyNode(SupplyNode):
  """ Smiles supplier

  Sample Usage:
    >>> fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib',\
                             'test_data','pgp_20.txt')
    >>> suppl = SmilesSupplyNode(fileN,delim="\\t",smilesColumn=2,nameColumn=1,titleLine=1)
    >>> ms = [x for x in suppl]
    >>> len(ms)
    20
    >>> ms[0].GetProp("_Name")
    'ALDOSTERONE'
    >>> ms[0].GetProp("ID")
    'RD-PGP-0001'
    >>> ms[1].GetProp("_Name")
    'AMIODARONE'
    >>> ms[3].GetProp("ID")
    'RD-PGP-0004'
    >>> suppl.reset()
    >>> suppl.next().GetProp("_Name")
    'ALDOSTERONE'
    >>> suppl.next().GetProp("_Name")
    'AMIODARONE'
    >>> suppl.reset()
  
  """
  def __init__(self,fileName,delim="\t",nameColumn=1,smilesColumn=0,titleLine=0,
               **kwargs):
    SupplyNode.__init__(self,**kwargs)
    self._fileName = fileName
    self._supplier = Chem.SmilesMolSupplier(self._fileName,delimiter=delim,
                                            smilesColumn=smilesColumn,
                                            nameColumn=nameColumn,
                                            titleLine=titleLine)

  def reset(self):
    SupplyNode.reset(self)
    self._supplier.reset()
  def next(self):
    """

    """
    r = None
    while not r:
      r = next(self._supplier)
    return r

if six.PY3:
    SmilesSupplyNode.__next__ = SmilesSupplyNode.next
  
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

  
