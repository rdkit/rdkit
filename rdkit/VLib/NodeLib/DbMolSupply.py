#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
import os.path
import sys

from rdkit import Chem, RDConfig
from rdkit.Chem.Suppliers import DbMolSupplier
from rdkit.VLib.Supply import SupplyNode


class DbMolSupplyNode(SupplyNode):
  """ Supplies molecules from a db result set:

  Sample Usage:
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> dbName = os.path.join(RDConfig.RDCodeDir,'Chem','Fingerprints',\
                             'test_data','data.gdb')
    >>> conn = DbConnect(dbName,'simple_mols')
    >>> dataset = conn.GetData()
    >>> suppl = DbMolSupplyNode(dataset)
    >>> ms = [x for x in suppl]
    >>> len(ms)
    12
    >>> ms[0].GetProp("ID")
    'ether-1'
    >>> ms[10].GetProp("ID")
    'acid-4'
    >>> suppl.reset()
    >>> suppl.next().GetProp("ID")
    'ether-1'
    >>> suppl.next().GetProp("ID")
    'acid-1'
    >>> suppl.reset()
  
  """

  def __init__(self, dbResults, **kwargs):
    SupplyNode.__init__(self, **kwargs)
    self._dbResults = dbResults
    self._supplier = DbMolSupplier.RandomAccessDbMolSupplier(self._dbResults, **kwargs)

  def reset(self):
    SupplyNode.reset(self)
    self._supplier.Reset()

  def next(self):
    """

    """
    return self._supplier.next()


def GetNode(dbName, tableName):
  from rdkit.Dbase.DbConnection import DbConnect
  conn = DbConnect(dbName, tableName)
  return DbMolSupplyNode(conn.GetData())


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest
  import sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed, tried = _test()
  sys.exit(failed)
