#  $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#     All Rights Reserved
#
import sys
from rdkit import six

from rdkit.VLib.Node import VLibNode

class OutputNode(VLibNode):
  """ base class for nodes which dump output

  Assumptions:

    - destination supports a write() method

    - strFunc, if provided, returns a string representation of
      the input

    - inputs (parents) can be stepped through in lockstep


  Usage Example:
    >>> from rdkit.VLib.Supply import SupplyNode
    >>> supplier = SupplyNode(contents=[1,2,3])
    >>> from rdkit.six import StringIO
    >>> sio = StringIO()
    >>> node = OutputNode(dest=sio,strFunc=lambda x:'%s '%(str(x)))
    >>> node.AddParent(supplier)
    >>> node.next()
    1
    >>> sio.getvalue()
    '1 '
    >>> node.next()
    2
    >>> sio.getvalue()
    '1 2 '

  """
  def __init__(self,dest=None,strFunc=None,**kwargs):
    VLibNode.__init__(self,**kwargs)
    self._dest = dest
    self._func=strFunc
  def next(self):
    parents = self.GetParents()
    args = []
    for parent in parents:
      try:
        args.append(parent.next())
      except StopIteration:
        raise StopIteration
    if len(args)>1:
      args = tuple(args)
    else:
      args = args[0]
    if self._func is not None:
      outp = self._func(args)
    else:
      outp = str(args)
    if self._dest:
      self._dest.write(outp)
    return args

if six.PY3:
    OutputNode.__next__ = OutputNode.next

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


  
