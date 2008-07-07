#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" mixin providing a visual programming canvas

"""
import RDConfig
from qt import *
from qtcanvas import *
import math
import VPItems,VPLib
from qtGui import qtUtils,GuiUndo

REQUIRED_MIXINS = []
MODULES_ALTERED = ['qtGui.GuiBase']
METHODS_DEFINED = {
  '__canvasUndo':'qtGui.GuiBase.GuiBase.vpcUndo',
  '__canvasRedo':'qtGui.GuiBase.GuiBase.vpcRedo',
  '__saveCanvas':'qtGui.GuiBase.GuiBase.vpcSaveCanvas',
  }
VARS_TO_SAVE=['self.items']

class MoveUndoObj(GuiUndo.UndoObj):
  def __init__(self,objs,delta):
    self._name = 'move'
    self._objs = tuple(objs)
    self._delta = delta
  def undo(self):
    for obj in self._objs:
      try:
        obj.moveBy(-self._delta[0],-self._delta[1])
      except AttributeError:
        pass
  def redo(self):
    for obj in self._objs:
      try:
        obj.moveBy(self._delta[0],self._delta[1])
      except AttributeError:
        pass
    
  def accepts(self,other):
    if isinstance(other,MoveUndoObj) and other._objs == self._objs:
      return 1
  def accept(self,other):
    self._delta = (self._delta[0]+other._delta[0],
                   self._delta[1]+other._delta[1])
    


class CanvasView(QCanvasView):
  """  A canvas for visual programming applications


  Notes:
    - member variables that things on the canvas might have and
      what they mean:

      - _selectable: can be selected
      - _moveable: can be dragged about the canvas
      - _linkable: can be linked to other objects (for LinkPts)

    - selections are handled using the
      _active list.  Multiple selections are possible.


  """
  def __init__(self,parent,name,canvas=None,size=(400,400),undoLevels=10):
    QCanvasView.__init__(self,parent,name,Qt.WSubWindow)
    if canvas is None:
      canvas = QCanvas(size[0],size[1])

    self.setCanvas(canvas)
    self._canvas = canvas
    self.items = []
    self._active = []

    # NOTE: According to the docs, this will slow stuff down
    self._scale = 1.0
    self._minScale = 0.1
    self._tformMat = QWMatrix()
    self._tformMat.setMatrix(self._scale,0.0,0.0,self._scale,0.0,0.0)
    self.setWorldMatrix(self._tformMat)
    self._startPos = [0,0]
    
    # we want keyboard events
    self.setFocusPolicy(QWidget.StrongFocus)

    # for highlighting selections
    self.selectPen=QPen(QColor(255,255,0),2,Qt.DashLine)
    self.selectStyle = Qt.Dense5Pattern

    self._undoStack = GuiUndo.UndoStack(undoLevels)
    
  def linkObjects(self,linker1,linker2):
    parent1,id1 = linker1.getParentInfo()
    parent2,id2 = linker2.getParentInfo()
    parent1.createLink(parent2,id1,id2)

  def unlinkObjects(self,linker1,linker2):
    parent1,id1 = linker1.getParentInfo()
    parent2,id2 = linker2.getParentInfo()
    for link in parent1.getLinks():
      if link.getOtherObj(parent1) == parent2:
        link.setVisible(0)
        link.destroy()

  def update(self):
    self._canvas.update()
    
  #-----------------
  #
  #  Undo support
  #
  #-----------------
  def __registerAddNodeUndo(self,node):
    """

      FIX:
      at the moment it's way too much trouble to actually undo
      adding the thing... we'll just cheat by hiding it.  hee hee.

    """
    undoAdd = node.setVisible
    undoArgs = (0,)
    redoAdd = node.setVisible
    redoArgs = (1,)
    obj = GuiUndo.UndoObj('add object',undoAdd,undoArgs,redoAdd,redoArgs)
    self._undoStack.push(obj)

  def __registerMoveUndo(self,moved,delta):
    obj =  MoveUndoObj(moved,delta)
    self._undoStack.push(obj)

  def __registerDeleteUndo(self,which):
    args = (tuple(which),)
    undoFunc = lambda x:map(lambda y:y.setVisible(1),x)
    redoFunc = lambda x:map(lambda y:y.setVisible(0),x)
    obj = GuiUndo.UndoObj('delete objects',
                          undoFunc,args,
                          redoFunc,args)
    self._undoStack.push(obj)

  def __registerLinkUndo(self,linker1,linker2):
    args = (linker1,linker2)
    obj = GuiUndo.UndoObj('link',
                          self.unlinkObjects,args,
                          self.linkObjects,args)
    self._undoStack.push(obj)

  #-----------------
  #
  #  Internal Event management
  #
  #-----------------
  def __selectObject(self,obj):
    if obj:
      obj._origPen = obj.pen()
      obj._origStyle = obj.brush().style()
      obj.setPen(self.selectPen)
      b = obj.brush()
      obj.setBrush(QBrush(b.color(),self.selectStyle))
  def __deselectObjects(self,objs):
    if objs:
      for obj in objs:
        try:
          obj.setPen(obj._origPen)
          b = obj.brush()
          obj.setBrush(QBrush(b.color(),obj._origStyle))
        except:
          pass

  def __leftMouseEvent(self,evt):
    clickPt = self.inverseWorldMatrix().map(evt.pos())
    items = self._canvas.collisions(clickPt)

    extendSel = evt.state() & (Qt.ShiftButton | Qt.ControlButton)
    
    if not extendSel:
      self.__deselectObjects(self._active)
      self._active = []
    go = 0
    # it's so rewarding to write UI code.
    if len(items):
      for item in items:
        try:
          go = item._selectable
        except AttributeError:
          go = 0
        if go:
          if go < 0:
            item = item.parent()
          if item not in self._active:
            self._active.append(item)
            self.__selectObject(item)
          else:
            if not extendSel:
              self.__deselectObjects([item])
              self._active.remove(item)
          break
    if (go and len(items)) or extendSel: self._startPos = clickPt

  def __leftDoubleMouseEvent(self,evt):
    clickPt = self.inverseWorldMatrix().map(evt.pos())
    items = self._canvas.collisions(clickPt)

    if len(items):
      for item in items:
        try:
          go = item._selectable
        except:
          go = 0
        if go < 0:
          item = item.parent()
        if go:
          try:
            # FIX: this doesn't want to stay this way
            item.showProperties()
          except AttributeError:
            pass
          break
          
  def __rightMouseEvent(self,evt):
    clickPt = self.inverseWorldMatrix().map(evt.pos())
    items = self._canvas.collisions(clickPt)
    if len(items):
      go = 0
      for item in items:
        try:
          go = item._selectable
        except AttributeError:
          go = 0
        if go < 0:
          item = item.parent()
        if go:
          self.__itemPopup(clickPt,item)
          break
    else:
      self.__bodyPopup(clickPt)
    

  def __itemPopup(self,pos,item):
    # FIX: this doesn't want to stay this way
    if hasattr(item,'showProperties'): item.showProperties()


  def __bodyPopup(self,pos):
    menu = QPopupMenu(self)
    linkId = menu.insertItem(self.trUtf8('&Link'),self.__linkActive)
    menu.setItemEnabled(linkId,0)
    if len(self._active) == 2:
      linkable = 0
      try:
        linkable = self._active[0]._linkable
        linkable = linkable & self._active[1]._linkable
      except AttributeError:
        pass
      if linkable:
        menu.setItemEnabled(linkId,1)
    menu.exec_loop(self.mapToGlobal(self.worldMatrix().map(pos)))

  def __linkActive(self):
    linker1 = self._active[0]
    linker2 = self._active[1]
    self.linkObjects(linker1,linker2)
    self.__deselectObjects(self._active)
    self._active = []
    self.__registerLinkUndo(linker1,linker2)

  def __delActive(self):
    if 0:
      for item in self._active:
        try:
          item.setVisible(0)
          item.destroy()
        except AttributeError:
          print 'cannot destroy: ',item
        try:
          self.items.remove(item)
        except:
          pass
    else:
      for item in self._active:
        item.setVisible(0)
    self.__registerDeleteUndo(self._active)    
    self.__deselectObjects(self._active)
    self._active = []

  def __moveObjects(self,which,delta):
    goAhead = 0
    for thing in which:
      try:
        if thing._moveable:
          thing.moveBy(delta[0],delta[1])
          goAhead = 1
      except AttributeError:
        thing.moveBy(delta[0],delta[1])
        goAhead = 1
    return goAhead


  def __moveActive(self,delta):
    self.__registerMoveUndo(self._active,delta)
    return self.__moveObjects(self._active,delta)

  def __cycleActive(self,direction):
    if len(self._active) == 1:
      idx = self.items.index(self._active[0])
      nItems = len(self.items)
      if direction > 0:
        idx += 1
        if idx == nItems:
          idx = 0
      elif direction < 0:
        if idx == 0:
          idx = nItems
        idx -= 1  
      item = self.items[idx]
      self.__deselectObjects(self._active)
      self.__selectObject(item)
      self._active = [item]
    elif len(self._active)  == 0:
      item = self.items[0]
      self.__deselectObjects(self._active)
      self.__selectObject(item)
      self._active = [item]

  #-----------------
  #
  #  Slots called by Qt (Event handlers)
  #
  #-----------------
  def resizeEvent(self,evt):
    # delegate, then resize the canvas
    QCanvasView.resizeEvent(self,evt)
    sz = evt.size()
    self._canvas.resize(sz.width(),sz.height())

  def contentsWheelEvent(self,evt):
    """

      if shift is down we scale the canvas.
      otherwise we pass on this event

    """  
    # FIX:  When zooming out, the canvas often ends up being
    #  smaller than the window... this looks DOOFY.
    delta = evt.delta()
    isShift = evt.state() & Qt.ShiftButton
    if isShift:
      scaleOff = .05*float(delta)/120.
      if self._scale + scaleOff < self._minScale:
        self._scale = self._minScale
      else:
        self._scale += scaleOff
      self._tformMat.setMatrix(self._scale,0.0,0.0,self._scale,0.0,0.0)
      self.setWorldMatrix(self._tformMat)
      self.parent().statusBar().message('Scale: %3.0f%%'%(self._scale*100))
    else:
      # pass the event up the chain
      evt.ignore()

  def contentsMouseDoubleClickEvent(self,evt):
    but = evt.button()
    if but == Qt.LeftButton:
      self.__leftDoubleMouseEvent(evt)
    
  def contentsMousePressEvent(self,evt):
    but = evt.button()
    if but == Qt.LeftButton:
      self.__leftMouseEvent(evt)
    elif but == Qt.RightButton:
      self.__rightMouseEvent(evt)
    self._canvas.update()

  def contentsMouseMoveEvent(self,evt):
    clickPt = self.inverseWorldMatrix().map(evt.pos())
    delta = clickPt.x()-self._startPos.x(),\
            clickPt.y()-self._startPos.y()
    goAhead = self.__moveActive(delta)
    self._canvas.update()
    if goAhead:  
      self._startPos = clickPt

  def keyPressEvent(self,evt):
    # NOTE: as of version 3.1, there's some kind of nastiness with
    #  things not getting properly updated after an event.
    #  The maintainers of PyQt appear to have no interest in fixing this.
    code = evt.key()
    if (evt.state() & Qt.ShiftButton):
      scale = 5
    else:
      scale = 1
    if code == Qt.Key_Delete:
      self.__delActive()
      evt.accept()
    elif code == Qt.Key_Up:
      self.__moveActive((0,-1*scale))
      evt.accept()
    elif code == Qt.Key_Down:
      self.__moveActive((0,1*scale))
      evt.accept()
    elif code == Qt.Key_Left:
      self.__moveActive((-1*scale,0))
      evt.accept()
    elif code == Qt.Key_Right:
      self.__moveActive((1*scale,0))
      evt.accept()
    elif code == Qt.Key_Tab:
      self.__cycleActive(1)
      evt.accept()
    elif code == Qt.Key_BackTab:
      self.__cycleActive(-1)
      evt.accept()
    else:
      evt.ignore()
    self._canvas.update()

def __canvasUndo(self):
  #self._vpcView._undoStack.dbg()
  try:
    self._vpcView._undoStack.undo()
  except GuiUndo.UndoError:
    print 'no undo'
  self._vpcView.update()

def __canvasRedo(self):
  #self._vpcView._undoStack.dbg()
  try:
    self._vpcView._undoStack.redo()
  except GuiUndo.UndoError:
    print 'no redo'

  self._vpcView.update()

def __saveCanvas(self):
  fName = str(QFileDialog.getSaveFileName('','Pickles (*.pkl)'))
  if fName:
    import cPickle
    try:
      outF = open(fName,'wb+')
    except IOError:
      qtUtils.error('cannot open save file %s for writing'%(fName))
      return
    
    #cPickle.dump(self._vpcView.items,outF)
    print 'nothing was saved'
    outF.close()

def LocalInit(self):
  """ mixin initialization code """
  self._saveMenu = QPopupMenu(self)
  self._fileMenu.insertItem(self.trUtf8("&Save"),self._saveMenu)

  self._editMenu.insertItem(self.trUtf8("&Undo"),self.vpcUndo,Qt.CTRL+Qt.Key_Z)
  self._editMenu.insertItem(self.trUtf8("&Redo"),self.vpcRedo,Qt.CTRL+Qt.SHIFT+Qt.Key_Z)

  self._vpcView = CanvasView(self,'_vpView',size=(400,400))
  self.setCentralWidget(self._vpcView)
  self._vpcView.show()
  qApp.setMainWidget(self)

  self._saveMenu.insertItem(self.trUtf8("&Canvas"),self.vpcSaveCanvas)
