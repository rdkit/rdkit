# $Id$
#
#  Copyright (C) 2003  Rational Discovery LLC
#    All Rights Reserved
#
""" point-sprinkler canvas

"""
import RDConfig
from qt import *
from qtcanvas import *
from qtGui import qtUtils,GuiUndo
from sping import utils
import PointCanvasItems

REQUIRED_MIXINS = []
MODULES_ALTERED = ['qtGui.GuiBase']
METHODS_DEFINED = {
  '__canvasUndo':'qtGui.GuiBase.GuiBase.psUndo',
  '__canvasRedo':'qtGui.GuiBase.GuiBase.psRedo',
  '__savePoints':'qtGui.GuiBase.GuiBase.psSavePoints',
  '__loadPoints':'qtGui.GuiBase.GuiBase.psLoadPoints',
  '__exportPoints':'qtGui.GuiBase.GuiBase.psExportPoints',
  }
VARS_TO_SAVE=[]

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
  """  A canvas for the point sprinkler


  Notes:
    - member variables that things on the canvas might have and
      what they mean:

      - _selectable: can be selected
      - _moveable: can be dragged about the canvas

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
    self._startPos = QPoint()
    
    # we want keyboard events
    self.setFocusPolicy(QWidget.StrongFocus)

    # for highlighting selections
    self.selectPen=QPen(QColor(150,150,150),1,Qt.DashLine)
    self.selectStyle = Qt.Dense5Pattern
    self.selectBrush=QBrush(QColor(255,255,0),self.selectStyle)
    self._undoStack = GuiUndo.UndoStack(undoLevels)

    self._addMode = 1
    self._nPts = 0

    self._emptyClick=None
    self._selRect = QCanvasRectangle(self._canvas)
    self._selRect.setPen(QPen(QColor(150,150,150),1,Qt.DashLine))

    self._defaultAct = 0
    self._defaultSymbol = 0
    
  def update(self):
    self._canvas.update()

  def addPoints(self,loc=None,rad=10,inPts=None,label=None,act=None,
                undoable=1,shape=None):
    if act is None:
      act = self._defaultAct
    if not inPts:
      self._nPts += 1
      if label is None:
        label = 'point-%d'%(self._nPts)

      if shape is None:
        shape = PointCanvasItems.shapes[self._defaultSymbol%len(PointCanvasItems.shapes)]
      itm = shape(self._canvas,rad,label=label,act=act)
      if hasattr(loc,'x'):
        itm.move(loc.x(),loc.y())
      else:
        itm.move(loc[0],loc[1])
      if undoable:
        self.__registerAddNodeUndo(itm)
      inPts = [itm]
      newPts=1
    else:
      newPts=0
    for itm in inPts:
      self.items.append(itm)
      itm.setVisible(1)
      if not newPts:
        self.__selectObject(itm)
        self._active.append(itm)
    self.update()    
  def delPoints(self,inPts,enableUndo=1):
    for item in inPts:
      item.setVisible(0)
      try:
        self.items.remove(item)
      except:
        pass
    if enableUndo:
      self.__registerDeleteUndo(inPts)    
    self.__deselectObjects(inPts)

  def getPoints(self):
    points = []
    for item in self.items:
      if isinstance(item,PointCanvasItems.PointCanvasMixin) and item.isVisible():
        points.append(item)
    return points

  def getNormalizedPoints(self):
    sz = self._canvas.width(),self._canvas.height()
    points = []
    for item in self.items:
      if isinstance(item,PointCanvasItems.PointCanvasMixin) and item.isVisible():
        loc = item.x()/sz[0],item.y()/sz[1]
        points.append((item.label(),loc,item.activity()))
    return points

  def labelPoints(self,val,activePts=None):
    if activePts is None:
      for obj in self.items:
        if obj.isSelected():
          obj.setActivity(val)
    else:
      for obj in activePts:
        obj.setActivity(val)
    self._defaultAct = val
    self.update()
  #-----------------
  #
  #  Undo support
  #
  #-----------------
  def __registerAddNodeUndo(self,node):
    """

    """
    undoFunc = lambda :self.delPoints((node,),enableUndo=0)
    redoFunc = lambda :self.addPoints(inPts=(node,))
    obj = GuiUndo.UndoObj('add object',undoFunc,(),redoFunc,())
    self._undoStack.push(obj)

  def __registerMoveUndo(self,moved,delta):
    obj =  MoveUndoObj(moved,delta)
    self._undoStack.push(obj)

  def __registerDeleteUndo(self,which):
    args = (tuple(which),)
    undoFunc = lambda x:self.addPoints(inPts=tuple(which))
    redoFunc = lambda x:self.delPoints(tuple(which),enableUndo=0)
    obj = GuiUndo.UndoObj('delete objects',
                          undoFunc,args,
                          redoFunc,args)
    self._undoStack.push(obj)
  #-----------------
  #
  #  Internal Event management
  #
  #-----------------
  def __selectObject(self,obj,labelIt=1):
    if obj:
      obj.setSelected(1)
      obj._origPen = obj.pen()
      obj._origBrush = obj.brush()
      obj.setPen(self.selectPen)
      b = obj.brush()
      obj._origStyle = b.style()
      b.setStyle(self.selectStyle)
      obj.setBrush(b)
      if labelIt and hasattr(obj,'id'):
        widg = self.topLevelWidget()
        if widg and widg.statusBar():
          widg.statusBar().message(str(obj.label()))

  def __deselectObjects(self,objs):
    if objs:
      for obj in objs:
        obj.setSelected(0)
        try:
          obj.setPen(obj._origPen)
          b = obj.brush()
          b.setStyle(obj._origStyle)
          obj.setBrush(b)
        except:
          pass

  def __leftMouseEvent(self,evt):
    clickPt = self.inverseWorldMatrix().map(evt.pos())
    if not (evt.state() & Qt.AltButton):
      items = self._canvas.collisions(clickPt)
      extendSel = evt.state() & (Qt.ShiftButton | Qt.ControlButton)
      if not extendSel:
        self.__deselectObjects(self._active)
        self._active = []
      go = 0
      # it's so rewarding to write UI code.
      if len(items):
        for item in items:
          go = 0
          if hasattr(item,'_selectable'):
            go = item._selectable
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
      #if (go and len(items)) or extendSel: self._startPos = clickPt
      self._startPos = clickPt
      if not len(items):
        self._emptyClick=1
      else:
        self._emptyClick=0
    else:
      self.addPoints(loc=clickPt)
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
    pass


  def __bodyPopup(self,pos):
    pass

  def __delActive(self):
    self.delPoints(self._active)
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
    if not len(self.items):
      return
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

  def _updateSelRect(self,loc):
    pt = [0,0]
    sz = [0,0]
    startP = self._startPos.x(),self._startPos.y()
    if startP[0] < loc[0]:
      pt[0] = startP[0]
      sz[0] = loc[0]-startP[0]
    else:
      pt[0] = loc[0]
      sz[0] = startP[0]-loc[0]
    if startP[1] < loc[1]:
      pt[1] = startP[1]
      sz[1] = loc[1]-startP[1]
    else:
      pt[1] = loc[1]
      sz[1] = startP[1]-loc[1]
    self._selRect.move(pt[0],pt[1])
    self._selRect.setSize(sz[0],sz[1])
    
  def _selectRegion(self,region=None,extendSel=0):
    """ if provided, region should be (x,y,width,height)
    """
    if not extendSel:
      self.__deselectObjects(self._active)
      self._active=[]
      
    if region is None:
      region = (self._selRect.x(),self._selRect.y(),
                self._selRect.width(),self._selRect.height())
    rect = QRect(region[0],region[1],region[2],region[3])
    items = self._canvas.collisions(rect)
    for item in items:
      go = 0
      if hasattr(item,'_selectable'):
        go = item._selectable
      if go:
        if go < 0:
          item = item.parent()
        if item not in self._active:
          self._active.append(item)
          self.__selectObject(item,labelIt=0)
        elif not extendSel:
          self.__deselectObjects([item])
          self._active.remove(item)
    
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
    if not (evt.state() & Qt.AltButton):
      clickPt = self.inverseWorldMatrix().map(evt.pos())
      delta = clickPt.x()-self._startPos.x(),\
              clickPt.y()-self._startPos.y()
      if not self._emptyClick:
        goAhead = self.__moveActive(delta)
      else:
        self._updateSelRect((clickPt.x(),clickPt.y()))
        self._selRect.setVisible(1)
        goAhead = 0
      self._canvas.update()
      if goAhead:  
        self._startPos = clickPt

  def contentsMouseReleaseEvent(self,evt):
    extendSel = evt.state() & (Qt.ShiftButton | Qt.ControlButton)
    if self._selRect.isVisible():
      self._selectRegion(extendSel=extendSel)
      self._selRect.setVisible(0)
      self._canvas.update()
  

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
    elif (evt.state() & Qt.ControlButton) and code==Qt.Key_A:
      region = (0,0,self.canvas().size().width(),self.canvas().size().height())
      self._selectRegion(region=region,extendSel=1)
    else:
      evt.ignore()
      
    self._canvas.update()

def __canvasUndo(self):
  #self._vpcView._undoStack.dbg()
  try:
    self._psView._undoStack.undo()
  except GuiUndo.UndoError:
    print 'no undo'
  self._psView.update()

def __canvasRedo(self):
  #self._vpcView._undoStack.dbg()
  try:
    self._psView._undoStack.redo()
  except GuiUndo.UndoError:
    print 'no redo'

  self._psView.update()

def __savePoints(self,id=None,fName=''):
  if not fName:
    fName = str(QFileDialog.getSaveFileName('','Text files (*.txt *.csv);;All files (*.*)'))
  if fName:
    pts = self._psView.getNormalizedPoints()
    if len(pts):
      outF = open(fName,'w+')
      ext = fName.split('.')[-1]
      if ext in ['csv','CSV']:
        delim = ','
      else:
        delim = '\t'
      for pt in pts:
        locTxt = ['% .4f'%x for x in pt[1]]
        txt = delim.join([str(pt[0])]+locTxt+[str(pt[2])])
        outF.write('%s\n'%(txt))
      outF.close()
      qtUtils.information("%d points saved."%(len(pts)))
    else:
      qtUtils.information("No points to save.")

def __loadPoints(self,id=None,fName=''):
  if not fName:
    fName = str(QFileDialog.getOpenFileName('','Text files (*.txt *.csv);;All files (*.*)'))
  if fName:
    view = self._psView
    width = view._canvas.width()
    height = view._canvas.height()
    inF = open(fName,'r')
    ext = fName.split('.')[-1]
    if ext in ['csv','CSV']:
      delim = ','
    else:
      delim = '\t'
    nPts = 0
    for line in inF.readlines():
      line = line.strip()
      if len(line):
        splitL = line.split(delim)
        if len(splitL)>=3:
          label = splitL[0]
          xP = float(splitL[1])*width
          yP = float(splitL[2])*height
          if len(splitL)>=4:
            act = int(splitL[3])
          else:
            act = 0
          view.addPoints(loc=(xP,yP),label=label,undoable=0,act=act)


          nPts += 1
      else:
        # blank lines trigger series switches.
        view._defaultSymbol+=1  
        
      inF.close()
    view._defaultSymbol+=1  
    #qtUtils.information("%d points loaded."%(nPts))

def __exportPoints(self,id=None,fName='',canvas=None):
  if not canvas:
    sz = self._psView.canvas().size()
    sz = sz.width(),sz.height()
    canvas,fName = qtUtils.GetPiddleCanvas(sz)
  if canvas:
    for pt in self._psView.getPoints():
      pt.spingRender(canvas)
  canvas.drawRect(0,0,sz[0],sz[1])
  canvas.save()

def LocalInit(self):
  """ mixin initialization code """
  self._fileMenu.insertItem(self.trUtf8("&Load"),self.psLoadPoints)
  self._fileMenu.insertItem(self.trUtf8("&Save"),self.psSavePoints)

  self._fileMenu.insertSeparator()
  self._fileMenu.insertItem(self.trUtf8("&Export"),self.psExportPoints)

  self._editMenu.insertItem(self.trUtf8("&Undo"),self.psUndo,Qt.CTRL+Qt.Key_Z)
  self._editMenu.insertItem(self.trUtf8("&Redo"),self.psRedo,Qt.CTRL+Qt.SHIFT+Qt.Key_Z)

  self._psView = CanvasView(self,'_psView',size=(400,400))

  self._psLabelMenu = QPopupMenu(self)
  self._editMenu.insertSeparator()
  self._editMenu.insertItem(self.trUtf8("&Label Selection"),self._psLabelMenu)
  self._psLabelFn0 = lambda x,y=self._psView:y.labelPoints(0)
  self._psLabelMenu.insertItem("&0",self._psLabelFn0)
  self._psLabelFn1 = lambda x,y=self._psView:y.labelPoints(1)
  self._psLabelMenu.insertItem("&1",self._psLabelFn1)
  self._psLabelFn2 = lambda x,y=self._psView:y.labelPoints(2)
  self._psLabelMenu.insertItem("&2",self._psLabelFn2)

  
  
  widg = self.centralWidget()
  widg.hide()
  self.setCentralWidget(self._psView)
  self._psView.show()
  qApp.setMainWidget(self)

