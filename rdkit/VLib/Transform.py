#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
from rdkit.VLib.Node import VLibNode


class TransformNode(VLibNode):
  """ base class for nodes which filter their input

    Assumptions:

      - transform function takes a number of arguments equal to the
        number of inputs we have.  We return whatever it returns

      - inputs (parents) can be stepped through in lockstep

    Usage Example:

      >>> from rdkit.VLib.Supply import SupplyNode
      >>> def func(a,b):
      ...   return a+b
      >>> tform = TransformNode(func)
      >>> suppl1 = SupplyNode(contents=[1,2,3,3])
      >>> suppl2 = SupplyNode(contents=[1,2,3,1])
      >>> tform.AddParent(suppl1)
      >>> tform.AddParent(suppl2)
      >>> v = [x for x in tform]
      >>> v
      [2, 4, 6, 4]
      >>> tform.reset()
      >>> v = [x for x in tform]
      >>> v
      [2, 4, 6, 4]

    If we don't provide a function, just return the inputs:

      >>> tform = TransformNode()
      >>> suppl1 = SupplyNode(contents=[1,2,3,3])
      >>> suppl2 = SupplyNode(contents=[1,2,3,1])
      >>> tform.AddParent(suppl1)
      >>> tform.AddParent(suppl2)
      >>> v = [x for x in tform]
      >>> v
      [(1, 1), (2, 2), (3, 3), (3, 1)]

    """

  def __init__(self, func=None, **kwargs):
    VLibNode.__init__(self, **kwargs)
    self._func = func

  def next(self):
    parent = self.GetParents()[0]
    args = []
    try:
      for parent in self.GetParents():
        args.append(parent.next())
    except StopIteration:
      raise StopIteration
    args = tuple(args)
    if self._func is not None:
      res = self._func(*args)
    else:
      res = args
    return res


TransformNode.__next__ = TransformNode.next


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
