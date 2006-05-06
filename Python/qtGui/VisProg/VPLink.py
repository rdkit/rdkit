#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
from qtcanvas import *
from qtGui.qtUtils import _objsSoFar

class VPLink:
  """ this is the link between objects on the canvas
      it also takes care of drawing the connecting line

  """
  def __init__(self,canvas,obj1,obj2,obj1Slot=None,obj2Slot=None,id=None):
    """

    Arguments:

      - canvas: where we belong

      - obj1,obj2: the objects to which this links (not the LinkPts, 
        the objects themselves)

      - obj1Slot,obj2Slot: (optional) the LinkPts to which this connects
        (determines endpoints)

      - id: (optional) an id to identify this guy on the canvas  

    """
    global _objsSoFar
    if id is None:
      id = _objsSoFar
      _objsSoFar += 1
    self._canvas = canvas  
    self._selectable=1
    self._moveable=0
    self._linkable=0
    self._id = id
    self._obj1 = obj1
    self._obj2 = obj2
    self._children=[]
    self._visible = 0
    
    if obj1Slot:
      self.slot1 = obj1Slot
    else:
      self.slot1 = obj1
    if obj2Slot:
      self.slot2 = obj2Slot
    else:
      self.slot2 = obj2

    self.generatePoints()
    self.show()
    self.setZ(-10)
    self.setPen(QPen(QColor(0,0,0),2))
    self.setBrush(QBrush(QColor(0,0,0)))
  def obj1(self):
    return self._obj1
  def obj2(self):
    return self._obj2
  def getOtherObj(self,obj):
    """ given one object to which this is connected, returns the other
        returns None on failure

    """    
    if obj == self._obj1:
      return self._obj2
    elif obj == self._obj2:
      return self._obj1
    else:
      return None

  def destroy(self):
    self._obj1.removeLink(self)
    self._obj2.removeLink(self)
    self._obj1 = None
    self._obj2 = None
    self._children = []

  def generatePoints(self):
    if not len(self._children):
      for i in range(3):
        self._children.append(VPLinkSegment(self._canvas,self))

    p1x,p1y = self.slot1.getPos()
    p2x,p2y = self.slot2.getPos()

    dx = p2x-p1x
    offset = dx/5

    
    self._children[0].setPoints(p1x,p1y,p1x+offset,p1y)
    self._children[1].setPoints(p1x+offset,p1y,p2x-offset,p2y)
    self._children[2].setPoints(p2x-offset,p2y,p2x,p2y)


  def pen(self):
    return self._pen
  def brush(self):
    return self._brush
  
  def setPen(self,pen):
    self._pen = pen
    for child in self._children:
      child.setPen(self._pen)

  def setBrush(self,brush):
    self._brush = brush
    for child in self._children:
      child.setBrush(self._brush)

  def visible(self):
    return self._visible

  def setVisible(self,state):
    self._visible = state
    for child in self._children:
      child.setVisible(state)
    if self._obj1._gParent is not None and self._obj2._gParent is not None:
      if state:
        self._obj1._gParent.addLinkTo(self._obj2._gParent)
      else:
        try:
          self._obj1._gParent.delLinkTo(self._obj2._gParent)
        except ValueError:
          self._obj2._gParent.delLinkTo(self._obj1._gParent)
        
  def setZ(self,z):
    self._z = z
    for child in self._children:
      child.setZ(self._z)
  def show(self):
    self._visible=1
    for child in self._children:
      child.show()

      
class VPLinkSegment(QCanvasLine):
  def __init__(self,canvas,link,id=None):
    """

    Arguments:

      - canvas: where we belong

      - id: (optional) an id to identify this guy on the canvas  

    """
    QCanvasLine.__init__(self,canvas)
    global _objsSoFar
    if id is None:
      id = _objsSoFar
      _objsSoFar += 1
    self._link = link
    self._selectable=-1
    self._moveable=0
    self._linkable=0
    self._id = id
      
  def parent(self):
    return self._link
