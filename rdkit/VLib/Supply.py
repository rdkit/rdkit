#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
from rdkit.VLib.Node import VLibNode


class SupplyNode(VLibNode):
  """ base class for nodes which supply things

    Assumptions:
      1) no parents

    Usage Example:
    
      >>> supplier = SupplyNode(contents=[1,2,3])
      >>> supplier.next()
      1
      >>> supplier.next()
      2
      >>> supplier.next()
      3
      >>> supplier.next()
      Traceback (most recent call last):
          ...
      StopIteration
      >>> supplier.reset()
      >>> supplier.next()
      1
      >>> [x for x in supplier]
      [1, 2, 3]


    """

  def __init__(self, contents=None, **kwargs):
    VLibNode.__init__(self, **kwargs)
    if contents is not None:
      self._contents = contents
    else:
      self._contents = []
    self._pos = 0

  def reset(self):
    VLibNode.reset(self)
    self._pos = 0

  def next(self):
    if self._pos == len(self._contents):
      raise StopIteration

    res = self._contents[self._pos]
    self._pos += 1
    return res

  def AddParent(self, parent, notify=1):
    raise ValueError('SupplyNodes do not have parents')


SupplyNode.__next__ = SupplyNode.next


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
