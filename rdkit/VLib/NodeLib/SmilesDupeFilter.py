#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
from rdkit import Chem
from rdkit.VLib.Filter import FilterNode


class DupeFilter(FilterNode):
  """ canonical-smiles based duplicate filter

  Assumptions:

    - inputs are molecules


  Sample Usage:
    >>> import os
    >>> from rdkit import RDConfig
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

  def __init__(self, **kwargs):
    FilterNode.__init__(self, func=self.filter, **kwargs)
    self._smisSeen = set()

  def reset(self):
    FilterNode.reset(self)
    self._smisSeen = set()

  def filter(self, cmpd):
    smi = Chem.MolToSmiles(cmpd)
    if smi not in self._smisSeen:
      self._smisSeen.add(smi)
      return 1
    else:
      return 0


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
