#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
from qtcanvas import *
import math
from VPLink import VPLink
from qtGui.qtUtils import _objsSoFar


class VPItemMixin:
  """ mixin class for things that go on the visual programming canvas

  """
  def __init__(self,id=None,graphParent=None):
    """

    Arguments:

      - id: a unique id for the item on the canvas

    """
    self._gParent = graphParent
    
    self._selectable=1
    self._moveable=1
    self._linkable=0
    
    self._vpLinks=[]
    self._vpLinkPts=[]
    self._inputIndices=[]
    self._outputIndices=[]
    self._pos = [0,0]
    self._text = None
    self._label = ''
    self._icon = 0
    self.setZ(10)
    global _objsSoFar
    if id is None:
      id = _objsSoFar
      _objsSoFar += 1
    self._id = id
    self._propDlg = None

  def _addLabel(self):
    self._font = QFont("Helvetica",12)
    self._text = QCanvasText(self._canvas)
    self._text.setText(self._label)
    self._text.setTextFlags(Qt.AlignHCenter|Qt.AlignBottom)
    self._text.setFont(self._font)
    self._text.moveBy(self._width/2,self._height-10)
    self._text.setZ(self.z()+2)
    self._text.show()

  def setLabelText(self,what):
    self._label = what
    self._text.setText(what)

  def destroy(self):
    """ destroy this object and any of its links

    """
    try:
      for link in self._vpLinks[:]:
        link.destroy()
    except:
      pass
    self._vpLinks = None
    self._vpLinkPts = None
    self._propDlg = None

    if self._gParent is not None:
      self._gParent.destroy()
      self._gParent = None

    if self._icon is not None:
      self._icon.destroy()
      self._icon = None
      
  def createLink(self,other,thisID=-1,otherID=-1):
    """ adds a new link

    Arguments:

     - other: the other object to link to

     - thisID: (Optional) numeric ID of the LinkPt to use here

     - otherID: (Optional) numeric ID of the LinkPt to use on the
       other object

    """
    # FIX: need link validation (i.e. don't link to ourself,
    #   don't duplicate links, don't link output to output, etc.)
    if thisID >= 0:
      thisSlot = self.getLinkPts()[thisID]
    if otherID >= 0:
      otherSlot = other.getLinkPts()[otherID]

    alreadyThere = 0
    for link in self.getLinks():
      if link.getOtherObj(self) == other and link.visible():
        alreadyThere = 1
        break

    if not alreadyThere:
      link = VPLink(self.canvas(),self,other,thisSlot,otherSlot)
      self.addLink(link)
      other.addLink(link)
    else:
      link.setVisible(1)
    print self,self._gParent,other._gParent
    if self._gParent is not None and other._gParent is not None:
      print 'linky'
      self._gParent.addLinkTo(other._gParent)
    return link

  def removeLink(self,link):
    """ removes a link from the list

      Arguments:

        - link: a VPLink

    """
    self._vpLinks.remove(link)
    
  def addLink(self,link):
    """ appends a link to the list

      Arguments:

        - link: a VPLink

    """
    self._vpLinks.append(link)

  def getLinkTo(self,other):
    """ returns the link we have to object other

      returns None on failure
    """
    for link in self._vpLinks:
      objs = link.obj1(),link.obj2()
      if other in objs:
        return link
    return None
  def getLinks(self):
    """ returns our list of links """
    return self._vpLinks

  def getLinkPts(self):
    """ returns our list of linkPts """
    return self._vpLinkPts

  def getInputs(self):
    res = map(lambda x,y=self._vpLinks,z=self:y[x].getOtherObj(z),
              self._inputIndices)
    return res

  def getOutputs(self):
    res = map(lambda x,y=self._vpLinks,z=self: y[x].getOtherObj(z),
              self._outputIndices)
    return res
  
  def setPos(self,x,y):
    """ sets our current position """
    self._pos=[x,y]
  def getPos(self):
    """ returns our position as a list  """
    return self._pos

  def _updateLinks(self):
    """ PRIVATE
       updates the positions of the end points of our links

    """
    for link in self.getLinks():
      link.generatePoints()

  def move(self,x,y):
    """  moves this object to a point and updates its links  """
    self._parentClass.move(self,x,y)
    self._updateLinks()

  def moveBy(self,dx,dy):
    """  moves this object by an offset and updates its links  """
    self._parentClass.moveBy(self,dx,dy)
    attach = self.getPos()
    attach[0] += dx
    attach[1] += dy
    for pt in self._vpLinkPts:
      pt.moveBy(dx,dy)
    self._updateLinks()
    if self._text: self._text.moveBy(dx,dy)
    if self._icon: self._icon.moveBy(dx,dy)
  def showProperties(self):
    print 'show properties'

  def setVisible(self,state):
    self._parentClass.setVisible(self,state)
    for pt in self._vpLinkPts:
      pt.setVisible(state)
    for link in self.getLinks():
      if state:
        # only show links if both end points are visible
        #  (this is relevant when undoing things on the canvas)
        if link.obj1().isVisible() and link.obj2().isVisible():
          link.setVisible(state)
      else:
        link.setVisible(state)
    try:
      self._text.setVisible(state)
    except AttributeError:
      pass
    try:
      self._icon.setVisible(state)
    except AttributeError:
      pass
    
class VPLinkPoint(QCanvasRectangle):
  """ these are the points used to establish where links connect to nodes on the
    canvas

  """
  def __init__(self,canvas,pos,parent,parentIdx,height=10,width=10,id=None):
    """

      Arguments:

        - canvas: the canvas on which it's drawn

        - pos: drawing position

        - parent: who owns us (i.e. which object)

        - parentIdx:  our index in the parents list of linkPts

        - height,width: (Optional) dimensions of our box

        - id: (Optional) unique id on the canvas

    """    
    QCanvasRectangle.__init__(self,0,0,width,height,canvas)
    self.setBrush(QBrush(QColor(180,180,180)))
    self.setPen(QPen(QColor(0,0,0),1))
    self.setZ(15)
    self._pos = pos
    self._parent = parent
    self._parentIdx = parentIdx
    self.move(pos[0],pos[1])
    # move so that the damn thing is centered
    self.moveBy(-width/2,-height/2)
    self.show()

    self._selectable=1
    self._moveable=0
    self._linkable=1

  def getPos(self):
    """ returns the location of our *center*  """
    return self.x()+self.width()/2,self.y()+self.height()/2

  def getParentInfo(self):
    """ returns a (parent,parentIdx) 2-tuple """
    return self._parent,self._parentIdx

  def destroy(self):
    self._parent = None

  
