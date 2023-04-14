#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#

from rdkit.VLib.Node import VLibNode


class FilterNode(VLibNode):
  """ base class for nodes which filter their input

    Assumptions:

      - filter function takes a number of arguments equal to the
        number of inputs we have.  It returns a bool

      - inputs (parents) can be stepped through in lockstep

      - we return a tuple if there's more than one input

    Usage Example:

      >>> from rdkit.VLib.Supply import SupplyNode
      >>> def func(a,b):
      ...   return a+b < 5
      >>> filt = FilterNode(func=func)
      >>> suppl1 = SupplyNode(contents=[1,2,3,3])
      >>> suppl2 = SupplyNode(contents=[1,2,3,1])
      >>> filt.AddParent(suppl1)
      >>> filt.AddParent(suppl2)
      >>> v = [x for x in filt]
      >>> v
      [(1, 1), (2, 2), (3, 1)]
      >>> filt.reset()
      >>> v = [x for x in filt]
      >>> v
      [(1, 1), (2, 2), (3, 1)]
      >>> filt.Destroy()

      Negation is also possible:

      >>> filt = FilterNode(func=func,negate=1)
      >>> suppl1 = SupplyNode(contents=[1,2,3,3])
      >>> suppl2 = SupplyNode(contents=[1,2,3,1])
      >>> filt.AddParent(suppl1)
      >>> filt.AddParent(suppl2)
      >>> v = [x for x in filt]
      >>> v
      [(3, 3)]
      >>> filt.Destroy()

      With no function, just return the inputs:

      >>> filt = FilterNode()
      >>> suppl1 = SupplyNode(contents=[1,2,3,3])
      >>> filt.AddParent(suppl1)
      >>> v = [x for x in filt]
      >>> v
      [1, 2, 3, 3]
      >>> filt.Destroy()

    """

  def __init__(self, func=None, negate=0, **kwargs):
    VLibNode.__init__(self, **kwargs)
    self._func = func
    self._negate = negate

  def SetNegate(self, state):
    self._negate = state

  def Negate(self):
    return self._negate

  def next(self):
    parents = self.GetParents()
    while 1:
      args = []
      try:
        for parent in parents:
          args.append(next(parent))
      except StopIteration:
        raise StopIteration
      args = tuple(args)
      if self._func is not None:
        r = self._func(*args)
        if self._negate:
          r = not r
          # sys.stderr.write('\t\tNEGATE -> %d\n'%(r))
        if r:
          res = args
          break
      else:
        res = args
        break
    if len(parents) == 1:
      res = res[0]
    return res


FilterNode.__next__ = FilterNode.next


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
