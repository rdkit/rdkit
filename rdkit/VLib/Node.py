#  $Id$
#
#  Copyright (C) 2003-2006 Rational Discovery LLC
#     All Rights Reserved
#
import sys
from rdkit import six

class VLibNode(object):
  """ base class for all virtual library nodes,
  defines minimal required interface

  """
  def __init__(self,*args,**kwargs):
    self._children = []
    self._parents = []

  #------------------------------------
  #
  #  Iteration
  #
  def __iter__(self):
    """ part of the iterator interface """
    self.reset()
    return self
  def next(self):
    """ part of the iterator interface

      raises StopIteration on failure
    """
    pass
  def reset(self):
    """ resets our iteration state

    """
    for parent in self.GetParents():
      parent.reset()


  #------------------------------------
  #
  #  Library graph operations
  #  Probably most of these won't need to be reimplemented in
  #  child classes
  #
  def AddChild(self,child,notify=1):
    """

    >>> p1 = VLibNode()
    >>> p2 = VLibNode()
    >>> c1 = VLibNode()
    >>> p1.AddChild(c1)
    >>> len(c1.GetParents())
    1
    >>> len(p1.GetChildren())
    1
    >>> p2.AddChild(c1,notify=0)
    >>> len(c1.GetParents())
    1
    >>> len(p2.GetChildren())
    1
    >>> c1.AddParent(p2,notify=0)
    >>> len(c1.GetParents())
    2
    >>> len(p2.GetChildren())
    1

    """
    self._children.append(child)
    if notify:
      child.AddParent(self,notify=0)
  def RemoveChild(self,child,notify=1):
    """
    >>> p1 = VLibNode()
    >>> c1 = VLibNode()
    >>> p1.AddChild(c1)
    >>> len(c1.GetParents())
    1
    >>> len(p1.GetChildren())
    1
    >>> p1.RemoveChild(c1)
    >>> len(c1.GetParents())
    0
    >>> len(p1.GetChildren())
    0
    """
    self._children.remove(child)
    if notify:
      child.RemoveParent(self,notify=0)
  def GetChildren(self):
    return tuple(self._children)
  
  def AddParent(self,parent,notify=1):
    """
    >>> p1 = VLibNode()
    >>> p2 = VLibNode()
    >>> c1 = VLibNode()
    >>> c1.AddParent(p1)
    >>> len(c1.GetParents())
    1
    >>> len(p1.GetChildren())
    1
    >>> c1.AddParent(p2,notify=0)
    >>> len(c1.GetParents())
    2
    >>> len(p2.GetChildren())
    0
    >>> p2.AddChild(c1,notify=0)
    >>> len(c1.GetParents())
    2
    >>> len(p2.GetChildren())
    1
    """
    self._parents.append(parent)
    if notify:
      parent.AddChild(self,notify=0)
  def RemoveParent(self,parent,notify=1):
    """
    >>> p1 = VLibNode()
    >>> c1 = VLibNode()
    >>> p1.AddChild(c1)
    >>> len(c1.GetParents())
    1
    >>> len(p1.GetChildren())
    1
    >>> c1.RemoveParent(p1)
    >>> len(c1.GetParents())
    0
    >>> len(p1.GetChildren())
    0
    """
    self._parents.remove(parent)
    if notify:
      parent.RemoveChild(self,notify=0)
  def GetParents(self):
    return tuple(self._parents)

  def Destroy(self,notify=1,propagateDown=0,propagateUp=0):
    """
    >>> p1 = VLibNode()
    >>> p2 = VLibNode()
    >>> c1 = VLibNode()
    >>> c2 = VLibNode()
    >>> p1.AddChild(c1)
    >>> p2.AddChild(c1)
    >>> p2.AddChild(c2)
    >>> len(c1.GetParents())
    2
    >>> len(c2.GetParents())
    1
    >>> len(p1.GetChildren())
    1
    >>> len(p2.GetChildren())
    2
    >>> c1.Destroy(propagateUp=1)
    >>> len(p2.GetChildren())
    0
    >>> len(c1.GetParents())
    0
    >>> len(c2.GetParents())
    0
    
    """
    #sys.stderr.write('DESTROY: %s\n'%(str(self)))
    if hasattr(self,'_destroyed'): return
    self._destroyed=1

    if notify:
      for o in self.GetChildren():
        o.RemoveParent(self,notify=0)
        if propagateDown:
          o.Destroy(notify=1,propagateDown=1,propagateUp=propagateUp)
      for o in self.GetParents():
        #sys.stderr.write('\tparent: %s\n'%(str(o)))
        o.RemoveChild(self,notify=0)
        if propagateUp:
          o.Destroy(notify=1,propagateDown=propagateDown,propagateUp=1)
    self._children = []
    self._parents = []
    
if six.PY3:
    VLibNode.__next__ = VLibNode.next


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


