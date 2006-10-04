## Automatically adapted for numpy.oldnumeric Sep 23, 2006 by alter_code1.py

#  $Id$
#
#  Copyright (C) 2003-2005  Rational Discovery LLC
#    All Rights Reserved
#
""" defines classes required for interacting with Dec Trees

#DOC entire module

"""    
import RDConfig
from qt import *
from qtGui.PiddleWindowImpl import PiddleWindow
from qtGui import qtUtils
from ML.DecTree import TreeVis
import cPickle,os,copy,types


class FloatItem(QListViewItem):
  def compare(self,other,col,ascending=0):
    v1 = str(self.text(col))
    v2 = str(other.text(col))
    try:
      v1 = float(v1)
    except ValueError:
      pass
    else:
      v2 = float(v2)
    return cmp(v1,v2)
    
      
      
    
class TreeWindow(PiddleWindow):
  """ a window class for interacting with decision trees

  import attributes:

    - tree
    
  
  """
  def __init__(self,parent=None,name='TreeVis',**kwargs):
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

    self._tree = None
    self._dir = '.'
    self._fileN = ''
    self._initMenubar()

    self.view().addLeftHandler(self.leftCanvasClick)
    self.view().addLeftDblHandler(self.leftCanvasDblClick)
    self.view().addRightHandler(self.rightCanvasClick)

    self.resizeCanvas(size)
    self.setCaption(name)

    self._initHighlightControls()

    self._exampleListCtrls = []
    
  def _initMenubar(self):
    mb = self.menuBar()
    self._actionMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&Actions'),self._actionMenu,-1,1)
    self._actionMenu.insertItem(self.trUtf8('&Redraw'),self.draw)
    self._optionsMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&Options'),self._optionsMenu,-1,2)
    self._scaleToggle = QCheckBox(self.trUtf8('&Scale Leaves'),
                                  self._optionsMenu)
    self._optionsMenu.insertItem(self._scaleToggle)
    self._scaleToggle.setChecked(0)
    self._pureToggle = QCheckBox(self.trUtf8('Show &Purity'),
                                 self._optionsMenu)
    self._optionsMenu.insertItem(self._pureToggle)
    self._pureToggle.setChecked(1)

  def _initHighlightControls(self):
    lout = QHBoxLayout(None,2,2,"hbox_layout_1")

    self.ptsLine = QLineEdit(self.centralWidget())
    lout.addWidget(self.ptsLine)
    self.connect(self.ptsLine,SIGNAL('returnPressed()'),self.highlightPointsSlot)

    self.ptsButton = QPushButton(self.centralWidget(),"ptsbutt")
    self.ptsButton.setText(self.trUtf8("Highlight Points"))
    lout.addWidget(self.ptsButton)
    self.connect(self.ptsButton,SIGNAL('pressed()'),self.highlightPointsSlot)

    self._mainLayout.addLayout(lout)
    self.updateGeometry()

    self._highlightedPoints = []

  def tree(self):
    return self._tree
  def setTree(self,tree):
    self._tree = tree

  def scaleLeaves(self):
    return self._scaleToggle.isChecked()
  def setScaleLeaves(self,val):
    self._scaleToggle.setChecked(val)
    
  def showPurity(self):
    return self._pureToggle.isChecked()
  def setShowPurity(self,val):
    self._pureToggle.setChecked(val)
    
  def _findNodeAtPoint(self,node,loc):
    """ INTERNAL USE ONLY
    Finds the node at a given location

    """
    try:
      p = node._bBox
    except AttributeError:
      p = None

    if p is not None:
      if (loc[0]>=p[0] and loc[0]<=p[2]) and \
         (loc[1]>=p[1] and loc[1]<=p[3]):
        return node
      else:
        for child in node.GetChildren():
          r = self._findNodeAtPoint(child,loc)
          if r is not None:
            return r
    return None

  def loadPickledTree(self,fileN=None):
    """ loads a pickled tree and displays it

    **Arguments**

      - fileN: (optional) the filename from which to load the
        tree.  If this is not provided, a _QFileDialog_ will be
        launched to prompt for the name

    """
    if not fileN:
      fileN = str(QFileDialog.getOpenFileName(self._dir,'Pickled files (*.pkl);;All files (*.*)'))
    if not fileN:
      return None
    else:
      try:
        inF = open(fileN,'rb')
      except IOError:
        qtUtils.error('could not open file %s for reading'%(fileN))
        return None
      QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
      try:
        tree = cPickle.load(inF)
        splitName = os.path.split(fileN)
        nameBase = splitName[-1]
        self._dir = os.sep.join(splitName[0:-1])
        self.setTree(tree)
      except:
        QApplication.restoreOverrideCursor()
        qtUtils.error('problems encountered loading a tree from pickle file: %s'%fileN,exc_info=True)
        return None
      QApplication.restoreOverrideCursor()
      self.draw()
      return tree

  def saveTree(self,fileN=None):
    if not fileN:
      fileN = str(QFileDialog.getSaveFileName(self._dir,'Pickled files (*.pkl);;All files (*.*)'))
    if not fileN:
      return None
    else:
      try:
        outF = open(fileN,'wb+')
      except IOError:
        qtUtils.error('could not open file %s for writing'%(fileN))
      self._fileN = fileN
      cPickle.dump(self.tree(),outF)
      outF.close()
    

  def draw(self,idx=None,nRes=2):
    """ causes this model (Tree) to be drawn in a given window

    **Arguments**

      - nRes: (optional) the number of possible result codes for the model

    **Notes**

      - the destination window should have the following methods:

         - _canvas()_ : which returns a Piddle/Sping canvas on which
           we can draw

         - _size()_ : returns a 2-tuple with the size of the canvas

         - _view()_ : returns an object (used for book-keeping, at the
           moment this should be a _PiddleCanvasView_)

         - _show()_ : return value ignored

      - drawing is done using _TreeVis.DrawTree()_

    """
    tree = self._tree

    self.view().reset()
    dims = self.view().size()
    canvWidth = dims.width()
    canvHeight = dims.height()

    if not hasattr(tree,'totNChildren') or not hasattr(tree,'nLevelsBelow'):
      TreeVis.CalcTreeNodeSizes(tree)
    TreeVis.visOpts = TreeVis.VisOpts()
    totWidth = tree.totNChildren * (TreeVis.visOpts.circRad+TreeVis.visOpts.horizOffset)
    if totWidth > canvWidth:
      while totWidth > canvWidth:
        TreeVis.visOpts.circRad *= .75
        TreeVis.visOpts.horizOffset *= .75
        totWidth = tree.totNChildren * (TreeVis.visOpts.circRad+TreeVis.visOpts.horizOffset)
      canvWidth = totWidth

    totHeight = (tree.nLevelsBelow+1) * TreeVis.visOpts.vertOffset + TreeVis.visOpts.circRad
    if totHeight > canvHeight:
      while totHeight > canvHeight:
        TreeVis.visOpts.vertOffset *= .75
        totHeight = (tree.nLevelsBelow+1)* TreeVis.visOpts.vertOffset + TreeVis.visOpts.circRad
      canvHeight = totHeight
      
    self.resizeCanvas((canvWidth,canvHeight))
    canv = self.canvas()
    canv.clear()

    TreeVis.DrawTree(tree,canv,nRes=nRes,scaleLeaves=self.scaleLeaves(),
                     showPurity=self.showPurity())
    canv.flush()
    self.view()._currentTree = tree
    self.show()
  def _doHighlights(self,pts):
    visOpts = TreeVis.visOpts
    canv = self.canvas()
    for pt in pts:
      try:
        pos = pt._bBox
      except AttributeError:
        pass
      else:
        canv.drawEllipse(pos[0],pos[1],pos[2],pos[3],
                         visOpts.highlightColor,visOpts.highlightWidth)
    self.view().update()
    self.view().canvas().update()

  def _pointFinder(self,tree,which,col=0,tagAll=0):
    if tree is None or not tree:
      return []
    highlights = []
    if tree.GetTerminal() or tagAll:
      pts = tree.GetExamples()
      for pt in pts:
        v = str(pt[col])
        if v in which:
          highlights.append(tree)
          break
    if not tree.GetTerminal():
      for child in tree.GetChildren():
        highlights += self._pointFinder(child,which,col=col,tagAll=tagAll)
    return highlights
  def highlightPoints(self,which):
    if self._tree is None or not self._tree:
      return
    if type(which) not in (types.ListType,types.TupleType):
      which = (which,)
    highlights=self._pointFinder(self._tree,which,col=0)
    if highlights:
      self._doHighlights(highlights)
        

  def launchExamplesWin(self,node,newWin=0):
    """ opens a window to display examples

    **Arguments**

      - node: the node to use

      - newWin: (optional) if nonzero, a new window will be opened,
        otherwise the current window will be reused

    """
    examples = node.GetExamples()
    if examples:
      if not newWin and len(self._exampleListCtrls):
        listCtrl = self._exampleListCtrls[-1]
      else:
        # FIX: add icon
        listCtrl = QListView()
        listCtrl.setCaption('Examples')
        self._exampleListCtrls.append(listCtrl)

      listCtrl.show()
      listCtrl.clear()
      while listCtrl.columns():
        listCtrl.removeColumn(0)
      nCols = len(examples[0])
      for i in range(nCols):
        # FIX: add reasonable column headers or tool tips
        listCtrl.addColumn('%d'%(i))
      for example in examples:
        itm = FloatItem(listCtrl)
        for i in range(nCols):
          if isinstance(example[i], types.FloatType):
            itm.setText(i,'%4.3g'%example[i])
          else:
            itm.setText(i,str(example[i]))
      self.statusBar().message('%d Examples'%(len(examples)))
    else:
      self.statusBar().message('No Examples')


  #
  # Callbacks and slots and stuff
  #
  def leftCanvasClick(self,view,evt):
    if self._tree is None:
      return
    clickPt = view.inverseWorldMatrix().map(evt.pos())
    loc = clickPt.x(),clickPt.y()

    node = self._findNodeAtPoint(self._tree,loc)
    if node is not None:
      if not node.GetTerminal():
        txt = 'Descriptor: %s'%(node.GetName())
        try:
          qb = node.GetQuantBounds()
        except AttributeError:
          pass
        else:
          txt += '  Bounds: [ ' + ' '.join(['%6.3g'%x for x in qb]) + ']'
      else:    
        txt = 'Terminal: %s'%(node.GetName())
      self.statusBar().message(txt)
    else:
      self.statusBar().clear()
      

  def leftCanvasDblClick(self,view,evt):
    """ #DOC

    """
    pass
      
  def rightCanvasClick(self,view,evt):
    """ callback for right mouse clicks in canvases containing 
      trees
      (method of _GuiBase_)    

    launches a (possibly new) examples window using the information from
    the clicked node (if there is one)


    """
    if self._tree is None:
      return
    clickPt = view.inverseWorldMatrix().map(evt.pos())
    loc = clickPt.x(),clickPt.y()

    node = self._findNodeAtPoint(self._tree,loc)
    if node is not None:
      newWin = evt.state()&Qt.ShiftButton
      self.launchExamplesWin(node,newWin=newWin)

  def fileOpen(self):
    pass

  def fileSave(self):
    self.saveTree(self._fileN)

  def fileSaveAs(self):
    self.saveTree()

  def highlightPointsSlot(self):
    ptList=[x.strip() for x in str(self.ptsLine.text()).split(',')]
    self.highlightPoints(ptList)



if __name__ == '__main__':
  from qtGui import Gui
  
  app,widg = Gui.Launcher(TreeWindow,None)

  widg.loadPickledTree(os.path.join(RDConfig.RDDataDir,'Ferromag','ferromag_tree.pkl'))
  app.exec_loop()
  widg.destroy(1)

