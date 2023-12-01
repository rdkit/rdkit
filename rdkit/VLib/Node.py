#  $Id$
#
#  Copyright (C) 2003-2006 Rational Discovery LLC
#     All Rights Reserved
#
import sys


class VLibNode(object):
  """ base class for all virtual library nodes,
    defines minimal required interface

    """

  def __init__(self, *args, **kwargs):
    self._children = []
    self._parents = []

  # ------------------------------------
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

  # ------------------------------------
  #
  #  Library graph operations
  #  Probably most of these won't need to be reimplemented in
  #  child classes
  #
  def AddChild(self, child, notify=1):
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
      child.AddParent(self, notify=0)

  def RemoveChild(self, child, notify=1):
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
      child.RemoveParent(self, notify=0)

  def GetChildren(self):
    return tuple(self._children)

  def AddParent(self, parent, notify=True):
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
      parent.AddChild(self, notify=False)

  def RemoveParent(self, parent, notify=True):
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
      parent.RemoveChild(self, notify=False)

  def GetParents(self):
    return tuple(self._parents)

  def Destroy(self, notify=True, propagateDown=False, propagateUp=False):
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
        >>> c1.Destroy(propagateUp=True)
        >>> len(p2.GetChildren())
        0
        >>> len(c1.GetParents())
        0
        >>> len(c2.GetParents())
        0

        """
    if hasattr(self, '_destroyed'):
      return
    self._destroyed = True

    if notify:
      for o in self.GetChildren():
        o.RemoveParent(self, notify=False)
        if propagateDown:
          o.Destroy(notify=True, propagateDown=True, propagateUp=propagateUp)
      for o in self.GetParents():
        o.RemoveChild(self, notify=False)
        if propagateUp:
          o.Destroy(notify=True, propagateDown=propagateDown, propagateUp=True)
    self._children = []
    self._parents = []


VLibNode.__next__ = VLibNode.next


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
