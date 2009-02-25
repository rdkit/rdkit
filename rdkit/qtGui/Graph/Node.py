#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#


class GraphNode:
  def __init__(self,type='Node'):
    self._type = type
    self._linksTo = []
  def addLinkTo(self,other):
    if other not in self._linksTo:
      self._linksTo.append(other)
  def delLinkTo(self,other):
    if other in self._linksTo:
      self._linksTo.remove(other)
  def destroy(self):
    for other in self._linksTo:
      other.delLinkTo(self)
    self._linksTo = []    
  def getChildren(self):
    return tuple(self._linksTo)
  def getType(self):
    return self._type


