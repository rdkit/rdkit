#
#  Copyright (C) 2003  Rational Discovery LLC
#    All Rights Reserved
#
""" defines classes required for interacting with hierarchical catalogs

#DOC entire module

"""    
import RDConfig
from qt import *
from qtGui.PiddleWindowImpl import PiddleWindow
from qtGui.HierarchyBrowser import BrowserInputImpl as BrowserInput
from qtGui.DbQueryWidgetImpl import DbQueryWidget
from Chem import FragmentCatalog
from DataStructs import HierarchyVis
import cPickle,os,types

class CatalogWindow(PiddleWindow):
  """ a window class for interacting with hierarchical catalogs

  import attributes:

    - catalog
    - bits
    
  
  """
  def __init__(self,parent=None,name='HierarchyVis',**kwargs):
    if kwargs.has_key('size'):
      sz = kwargs['size']
      if type(sz) not in [types.TupleType(),types.ListType()]:
        size = [sz.width(),sz.height()]
      else:
        size = tuple(sz)
      del kwargs['size']
    else:
      size = (600,600)
    PiddleWindow.__init__(self,parent,name,**kwargs)

    self._catalog = None
    self._bits = None
    self._oBits = None
    self._highlights = []
    self._dir = '.'
    self._fileN = ''
    self._initMenubar()

    self.view().addLeftHandler(self.leftCanvasClick)
    self.view().addLeftDblHandler(self.leftCanvasDblClick)
    self.view().addRightHandler(self.rightCanvasClick)

    self.resizeCanvas(size)
    self.setCaption(name)

    self._initHighlightControls()
    self._loadDlg = None

    
  def _initMenubar(self):
    mb = self.menuBar()

    self._actionMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&Actions'),self._actionMenu,-1,1)
    id = self._actionMenu.insertItem(self.trUtf8('&Redraw'),self.draw,
                                     QKeySequence(self.trUtf8("Ctrl+R")))
    self._redrawId = id
    self._actionMenu.setItemEnabled(self._redrawId,0)
    id = self._actionMenu.insertItem(self.trUtf8('&Params'),self.loadCatalog)
    self._paramsId = id
    self._actionMenu.setItemEnabled(self._paramsId,0)

    id = self._actionMenu.insertItem(self.trUtf8('&Hide Highlighted'),self.hideHighlights)
    self._hideId = id
    self._actionMenu.setItemEnabled(self._hideId,0)

    self._optionsMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&Options'),self._optionsMenu,-1,2)
    self._colorToggle = QCheckBox(self.trUtf8('&Color Nodes'),self._optionsMenu)
    self._optionsMenu.insertItem(self._colorToggle)
    self._colorToggle.setChecked(0)
    self._byClassToggle = QCheckBox(self.trUtf8('Color by C&lass'),self._optionsMenu)
    self._optionsMenu.insertItem(self._byClassToggle)
    self._byClassToggle.setChecked(1)
    self._scaleGainsToggle = QCheckBox(self.trUtf8('&Scale Gains'),self._optionsMenu)
    self._optionsMenu.insertItem(self._scaleGainsToggle)
    self._scaleGainsToggle.setChecked(1)

    self.fileMenu.insertItem(self.trUtf8('&Open'),self.fileOpen,
                             QKeySequence( self.trUtf8("Ctrl+O", "File|Open")),
                             -1,0)



  def _initHighlightControls(self):
    lout = QHBoxLayout(None,2,2,"hbox_layout_1")

    self.ptsLine = QLineEdit(self.centralWidget())
    lout.addWidget(self.ptsLine)
    self.connect(self.ptsLine,SIGNAL('returnPressed()'),self.highlightPointsSlot)

    QToolTip.add(self.ptsLine,
                 self.trUtf8("Enter a comma-separated list of bit IDs to be highlighted."))

    self.ptsButton = QPushButton(self.centralWidget(),"ptsbutt")
    self.ptsButton.setText(self.trUtf8("Highlight Points"))
    lout.addWidget(self.ptsButton)
    self.connect(self.ptsButton,SIGNAL('clicked()'),self.highlightPointsSlot)
    QToolTip.add(self.ptsButton,
                 self.trUtf8("Highlight the specified bits in the hierarchy."))

    self._mainLayout.addLayout(lout)
    self._highlightLayout = lout
    self.updateGeometry()
    self._highlightedPoints = []
    
    
  def catalog(self):
    return self._catalog
  def setCatalog(self,val):
    self._catalog = val
  def bits(self):
    return self._bits
  def setBits(self,val,reset=0):
    if reset:
      self._oBits = val[:]
    self._bits = val

  def _findEntryAtPoint(self,loc,catalog):
    """ INTERNAL USE ONLY
    Finds the node at a given location

      - nodeLocs: a dictionary mapping node id -> location information
        ((x,y),rad)

    """
    if hasattr(catalog,'_drawLocs'):
      nodeLocs = catalog._drawLocs
    else:
      nodeLocs = {}
    for entry in nodeLocs.keys():
      nLoc,rad = nodeLocs[entry]
      p = [nLoc[0]-rad,nLoc[1]-rad,nLoc[0]+rad,nLoc[1]+rad]
      if p is not None:
        if (loc[0]>=p[0] and loc[0]<=p[2]) and \
           (loc[1]>=p[1] and loc[1]<=p[3]):
          return entry
    return None

  def draw(self,id=-1,minLevel=-1,maxLevel=-1,bitList=None):
    """ causes the hierarchy to be drawn in a given window

    **Notes**

      - the destination window should have the following methods:

         - _canvas()_ : which returns a Piddle/Sping canvas on which
           we can draw

         - _size()_ : returns a 2-tuple with the size of the canvas

         - _view()_ : returns an object (used for book-keeping, at the
           moment this should be a _PiddleCanvasView_)

         - _show()_ : return value ignored

      - drawing is done using _HierarchyVis.DrawHierarchy()_

    """
    if minLevel < 0:
      minLevel = self.minLevel
    if maxLevel < 0:
      maxLevel = self.maxLevel
    if not bitList:
      bitList = self.bits()
      if len(bitList)>self.numBits:
        bitList = bitList[:self.numBits]
      
    catalog = self.catalog()
    if not hasattr(catalog,'adjList') or not catalog.adjList:
      adjList,levels = FragmentCatalog.BuildAdjacencyList(catalog,bitList)
      catalog.adjList = adjList
      catalog.levelList = levels
    self.view().reset()
    #dims = self.view().size()
    #canvWidth = dims.width()
    #canvHeight = dims.height()
    canvWidth = self.view().visibleWidth()
    canvHeight = self.view().visibleHeight()
    minSize = HierarchyVis.GetMinCanvasSize(catalog.adjList,catalog.levelList)
    canvWidth = max(canvWidth,minSize[0])
    canvHeight = max(canvHeight,minSize[1])
    self.resizeCanvas((canvWidth,canvHeight))
    canv = self.canvas()
    canv.clear()

    drawColors = {}
    if self._colorToggle.isChecked():
      colorKlass = HierarchyVis.piddle.Color
      if self._scaleGainsToggle.isChecked():
        minGain = 1.0
        maxGain = 0.0
        for entry in bitList:
          g = entry.gain
          minGain = min(g,minGain)
          maxGain = max(g,maxGain)
        dGain = maxGain-minGain
      else:
        minGain = 0.0
        dGain = 1.0
      for entry in bitList:
        # FIX: this should not be specific to 2-class problems:
        if self._byClassToggle.isChecked() and entry.nPerClass[0] > entry.nPerClass[1]:
          c2 = colorKlass(.9,.1,.1)
          c1 = colorKlass(.9,.8,.8)
        else:
          c2 = colorKlass(.1,.1,.9)
          c1 = colorKlass(.8,.8,.9)
        drawColors[entry.id] = c1 + (c2-c1)*(entry.gain-minGain)/dGain
    catalog._drawLocs = HierarchyVis.DrawHierarchy(catalog.adjList,
                                                   catalog.levelList,
                                                   canv,minLevel=minLevel,
                                                   maxLevel=maxLevel,
                                                   entryColors=drawColors)
    canv.flush()
    self.view()._catalog = catalog
    self.show()
    self._highlights = []
    
  def _doHighlights(self,ids):
    catalog = self.catalog()
    visOpts = HierarchyVis.visOpts
    canv = self.canvas()
    for id in ids:
      try:
        pos,rad = catalog._drawLocs[id]
      except KeyError:
        pass
      else:
        canv.drawEllipse(pos[0]-rad,pos[1]-rad,pos[0]+rad,pos[1]+rad,
                         visOpts.highlightColor,visOpts.highlightWidth)
    self._highlights = ids[:]
    self.view().update()
    self.view().canvas().update()

  def highlightPoints(self,which):
    catalog = self.catalog()
    if catalog is None or not catalog or \
       not hasattr(catalog,'_drawLocs'):
      return
    if type(which) not in (types.ListType,types.TupleType):
      which = (which,)
    self._doHighlights(which)

  def loadCatalog(self,id=-1,catalog=None):
    if catalog is None:
      if self._loadDlg is None:
        dlg = QDialog(self,'Catalog Loader',1,Qt.WType_Dialog)
        dlg.setCaption('Catalog Loader')
        widg = BrowserInput.insertBrowserInputWidget(dlg)
        lout = QHBoxLayout()
        okButton = QPushButton('Ok',dlg)
        dlg.connect(okButton,SIGNAL("clicked()"),dlg.accept)
        cancelButton = QPushButton('Cancel',dlg)
        dlg.connect(cancelButton,SIGNAL("clicked()"),dlg.reject)
        lout.addWidget(okButton)
        lout.addStretch()
        lout.addWidget(cancelButton)
        dlg.layout().addLayout(lout)
        dlg.widg = widg
        self._loadDlg = dlg
      else:
        dlg = self._loadDlg
        widg = dlg.widg
      
      res = dlg.exec_loop()
      if res==QDialog.Accepted:
        self._actionMenu.setItemEnabled(self._paramsId,1)
        self._actionMenu.setItemEnabled(self._redrawId,1)
        self._actionMenu.setItemEnabled(self._hideId,1)
        cat = widg.catalog()
        cat.adjList=None
        self.setCatalog(cat)
        self.setBits(widg.bits(),reset=1)
        self.numBits = widg.numBitsSpin.value()
        self.minLevel = widg.minLevelSpin.value()
        self.maxLevel = widg.maxLevelSpin.value()
        bits = self.bits()[:self.numBits]
        self.draw(minLevel=self.minLevel,maxLevel=self.maxLevel,bitList=bits)
      else:
        print 'nope'

  def linkDown(self,which):
    catalog = self.catalog()
    if catalog is None or not catalog or \
       not hasattr(catalog,'_drawLocs'):
      return
    if type(which) not in (types.ListType,types.TupleType):
      stack = [which]
    else:
      stack = list(which)
    highlights = []
    while stack:
      # this is super easy because we can use use the catalog's
      # own down entries functionality
      top = stack.pop()
      highlights.append(top)
      for down in catalog.GetEntryDownIds(catalog.GetBitEntryId(top)):
        down = catalog.GetEntryBitId(down)
        if down not in stack and \
           down not in highlights:
          stack.append(down)
    if highlights:
      self.draw()
      self.highlightPoints(highlights)
      
  def linkUp(self,which):
    catalog = self.catalog()
    if catalog is None or not catalog or \
       not hasattr(catalog,'_drawLocs'):
      return
    if type(which) not in (types.ListType,types.TupleType):
      stack = [which]
    else:
      stack = list(which)
    highlights = []
    while stack:
      # FIX: this won't skip a level if an intermediate bit is missing.
      top = stack.pop()
      highlights.append(top)
      orderHere = catalog.GetBitOrder(top)
      upOrder = orderHere-1
      if catalog.levelList.has_key(upOrder):
        ups = catalog.levelList[upOrder]
        for up in ups:
          if top in catalog.adjList[up] and \
             up not in stack and\
             up not in highlights:
            stack.append(up)
    if highlights:
      self.draw()
      self.highlightPoints(highlights)
      
  def hideHighlights(self,id=-1,highlights=None):
    if highlights is None:
      highlights = self._highlights
    if not highlights:
      return
    if len(self._bits) > self.numBits:
      bits = self._bits[:self.numBits]
    else:
      bits = self._bits[:]
    for highlight in highlights:
      for bit in bits:
        if bit.id == highlight:
          bits.remove(bit)
          break
    self.setBits(bits)
    catalog = self.catalog()
    adjList,levels = FragmentCatalog.BuildAdjacencyList(catalog,bits)
    catalog.adjList = adjList
    catalog.levelList = levels
    self.draw(bitList=bits)    

      
  #
  # Callbacks and slots and stuff
  #
  def fileClose(self):
    """ callback for File->Close

    """
    self.close()


  def leftCanvasClick(self,view,evt):
    if self.catalog() is None:
      return
    clickPt = view.inverseWorldMatrix().map(evt.pos())
    loc = clickPt.x(),clickPt.y()

    entry = self._findEntryAtPoint(loc,self.catalog())
    if entry is not None:
      descr = self.catalog().GetBitDescription(entry)
      for bit in self.bits():
        if bit.id==entry:
          break
      text = "Bit %d: %s | Gain: %.4f | %s"%(entry,descr,bit.gain,str(bit.nPerClass))
      self.statusBar().message(text)
    else:
      self.statusBar().clear()
      
  def leftCanvasDblClick(self,view,evt):
    """ #DOC

    """
    pass
      
  def rightCanvasClick(self,view,evt):
    """ callback for right mouse clicks
      (method of _GuiBase_)    

    """
    if self.catalog() is None:
      return
    clickPt = view.inverseWorldMatrix().map(evt.pos())
    loc = clickPt.x(),clickPt.y()
    entry = self._findEntryAtPoint(loc,self.catalog())
    if entry is not None:
      if evt.state() & Qt.ShiftButton:
        self.linkDown([entry])
      else:
        self.linkUp([entry])


  def fileOpen(self):
    self.loadCatalog()

  def fileSave(self):
    pass

  def fileSaveAs(self):
    pass

  def highlightPointsSlot(self):
    text = str(self.ptsLine.text())
    if text:
      ptList=[int(x.strip()) for x in text.split(',')]
      self.highlightPoints(ptList)

if __name__ == '__main__':
  from qtGui import Gui
  
  app,widg = Gui.Launcher(CatalogWindow,None)
  app.exec_loop()
  widg.destroy(1)

