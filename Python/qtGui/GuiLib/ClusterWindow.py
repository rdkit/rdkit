#  $Id$
#
#  Copyright (C) 2002-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" defines classes required for interacting with Clusters

#DOC entire module

"""    
import RDConfig
from qt import *
from qtGui.PiddleWindowImpl import PiddleWindow
from qtGui import qtUtils
from ML.Cluster import ClusterVis,Clusters

import cPickle,os,copy,types,re

class ClusterWindow(PiddleWindow):
  """ a window class for interacting with cluster trees

  import attributes:

    - cluster
    
  
  """
  def __init__(self,parent=None,name='ClusterVis',**kwargs):
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

    self._cluster = None
    self._dists = None
    self._children = []
    self._clickSlack = 3
    self._dir = '.'
    self._fileN = ''
    self._tooClose=2
    self._initMenubar()

    self.view().addLeftHandler(self.leftCanvasClick)
    self.view().addLeftDblHandler(self.leftCanvasDblClick)
    self.view().addRightHandler(self.rightCanvasClick)

    self.resize(size[0]+50,size[1]+50)
    self.resizeCanvas(size)
    self.setCaption(name)

    self._initHighlightControls()

    self._clusterClickCallbacks = []
    self._pointClickCallbacks = []

  def _initMenubar(self):
    mb = self.menuBar()
    self._actionMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&Actions'),self._actionMenu,-1,1)
    self._actionMenu.insertItem(self.trUtf8('&Redraw'),self.draw,Qt.CTRL+Qt.Key_L)
    self._divPickMenuId = self._actionMenu.insertItem(self.trUtf8('&Diversity Pick'),
                                                      self.diversityPick)
    self._actionMenu.setItemEnabled(self._divPickMenuId,False)

    self._optionsMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&Options'),self._optionsMenu,-1,2)
    self._logToggle = QCheckBox(self.trUtf8('&Log Scale'),self._optionsMenu)
    self._optionsMenu.insertItem(self._logToggle)
    self._logToggle.setChecked(0)
    self._labelToggle = QCheckBox(self.trUtf8('L&abel Points'),self._optionsMenu)
    self._optionsMenu.insertItem(self._labelToggle)
    self._labelToggle.setChecked(0)

    self._optionsMenu.insertItem(self.trUtf8('&Too Close'),self.updateTooClose)

  def _initHighlightControls(self):
    lout = QHBoxLayout(None,2,2,"hbox_layout_1")

    self.ptsLine = QLineEdit(self.centralWidget())
    lout.addWidget(self.ptsLine)
    self.connect(self.ptsLine,SIGNAL('returnPressed()'),self.highlightPointsSlot)
    QToolTip.add(self.ptsLine,self.trUtf8('Enter a list of point names (comma-delimited) to be highlighted in the cluster tree.\nHit <Enter> or press the "Highlight Points" button to highlight them.'))

    self.ptsButton = QPushButton(self.centralWidget(),"ptsbutt")
    self.ptsButton.setText(self.trUtf8("Highlight Points"))
    lout.addWidget(self.ptsButton)
    self.connect(self.ptsButton,SIGNAL('pressed()'),self.highlightPointsSlot)
    QToolTip.add(self.ptsButton,self.trUtf8("Highlight the named points."))

    self.batchPtsButton = QPushButton(self.centralWidget(),"batchptsbutt")
    self.batchPtsButton.setText(self.trUtf8("Highlight Points From Clipboard"))
    lout.addWidget(self.batchPtsButton)
    self.connect(self.batchPtsButton,SIGNAL('pressed()'),self.highlightBatchPointsSlot)
    QToolTip.add(self.batchPtsButton,
                 self.trUtf8("""Grabs text currently on the clipboard and highlights the named points.
Clipboard text can be comma, space or tab delimted and may be spread across multiple lines.
Copied blocks from Excel can also be used (e.g. copy a column of point names and
  then highlight them with this button)."""))

    self.clustsLine = QLineEdit(self.centralWidget())
    lout.addWidget(self.clustsLine)
    self.connect(self.clustsLine,SIGNAL('returnPressed()'),self.highlightClustersSlot)
    QToolTip.add(self.clustsLine,self.trUtf8('Enter a list of cluster IDs (comma-delimited) to be highlighted in the cluster tree.\nHit <Enter> or press the "Highlight Clusters" button to highlight them.'))

    self.clustsButton = QPushButton(self.centralWidget(),"clustsbutt")
    self.clustsButton.setText(self.trUtf8("Highlight Clusters"))
    lout.addWidget(self.clustsButton)
    self.connect(self.clustsButton,SIGNAL('pressed()'),self.highlightClustersSlot)
    QToolTip.add(self.clustsButton,self.trUtf8("Highlight the selected clusters."))

    self._mainLayout.addLayout(lout)
    self.updateGeometry()

    self._highlightedPoints = []

  def addClusterClickCallback(self,cb):
    self._clusterClickCallbacks.append(cb)
  def clearClusterClickCallbacks(self):
    self._clusterClickCallbacks = []
  def addPointClickCallback(self,cb):
    self._pointClickCallbacks.append(cb)
  def clearPointClickCallbacks(self):
    self._pointClickCallbacks = []

  def cluster(self):
    return self._cluster
  def setCluster(self,cluster):
    self._cluster = cluster
    if hasattr(cluster,'_ptLabels'):
      self._clusterLabels = cluster._ptLabels
    else:
      self._clusterLabels = []

    if hasattr(cluster,'_ptValues'):
      self._clusterActs = cluster._ptValues
    else:
      self._clusterActs = []

  def distances(self):
    return self._dists
  def setDistances(self,dists):
    self._dists = dists
    if dists:
      self._actionMenu.setItemEnabled(self._divPickMenuId,True)
    else:
      self._actionMenu.setItemEnabled(self._divPickMenuId,False)
      
  def slack(self):
    return self._clickSlack
  def setSlack(self,slack):
    self._clickSlack = slack

  def logScale(self):
    return self._logToggle.isChecked()
  def setLogScale(self,val):
    self._logToggle.setChecked(val)
  def showLabels(self):
    return self._labelToggle.isChecked()
  def setShowLabels(self,val):
    self._labelToggle.setChecked(val)
    
  def _findClusterAtPoint(self,cluster,loc):
    """ INTERNAL USE ONLY
    Finds the cluster at a given location

    """
    if hasattr(cluster,'_drawPos'):
      p = cluster._drawPos
    else:
      p = None

    if p is not None:
      if abs(p[0]-loc[0]) < self._clickSlack and \
         abs(p[1]-loc[1]) < self._clickSlack:
        return cluster
      else:
        for child in cluster.GetChildren():
          r = self._findClusterAtPoint(child,loc)
          if r is not None:
            return r
    return None

  def diversityPick(self,idx=0):
    if not (self._cluster and self._dists):
      return

    if not hasattr(self,'_divPickDlg') or not self._divPickDlg:
      dlg = QDialog()
      self._divPickDlg = dlg
      dlg.setCaption('Diversity Picker')
      lout = QVBoxLayout(dlg)
      dlg.sb = QSpinBox(dlg)
      dlg.sb.setMinValue(0)
      dlg.sb.setMaxValue(len(self._cluster.GetPoints()))
      dlg.sb.setValue(5)
      lout.addWidget(dlg.sb)

      hout = QHBoxLayout()

      b1 = QPushButton(dlg)
      b1.setText('Pick')
      l1 = lambda y=dlg,z=self:z.pickCenters(y.sb.value())
      dlg.b1 = b1
      dlg.l1 = l1
      dlg.connect(b1,SIGNAL('clicked()'),l1)
      hout.addWidget(b1)

      b2 = QPushButton(dlg)
      b2.setText('Cancel')
      l2 = lambda y=dlg:y.hide()
      dlg.b2 = b2
      dlg.l2 = l2
      dlg.connect(b2,SIGNAL('clicked()'),l2)

      hout.addWidget(b2)

      lout.addLayout(hout)

    self._divPickDlg.show()

  def pickCenters(self,nToPick,method=None):
    import SimDivFilters
    if nToPick <= 0 or not (self._cluster and self._dists):
      return
    if method is None:
      method = SimDivFilters.ClusterMethod.WARD
    picker = SimDivFilters.HierarchicalClusterPicker(method)
    sz = len(self._cluster.GetPoints())
    picks = picker.Pick(self._dists,sz,nToPick)
    if self._clusterLabels:
      picks = [self._clusterLabels[x] for x in picks]
    else:
      picks = list(picks)
    print 'Picks:'
    for pick in picks:
      print '\t',pick
    self.draw()
    self.highlightPoints(picks)

  def launchExamplesWin(self,view,cluster,newWin=1):
    if not newWin:
      cluster.Print()
    else:
      pts = cluster.GetPoints()
      for i,pt in enumerate(pts):
        print i+1,'\t',pt.GetName()

  def loadPickledClusterTree(self,fileN=None):
    """ loads a pickled cluster and displays it

    **Arguments**

      - idx: (optional) not used

      - fileN: (optional) the filename from which to load the
        cluster tree.  If this is not provided, a _QFileDialog_ will be
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
        clust = cPickle.load(inF)
        if not isinstance(clust,Clusters.Cluster):
          clust=None
          raise ValueError,'invalid cluster pickle'
        splitName = os.path.split(fileN)
        nameBase = splitName[-1]
        self._dir = os.sep.join(splitName[0:-1])
        self.setCluster(clust)
      except:
        QApplication.restoreOverrideCursor()
        qtUtils.error("Problems encounted constructing cluster from pickle file %s\n"%fileN,
                      exc_info=True)
        return None
      QApplication.restoreOverrideCursor()
      self.draw()
      return clust

  def saveClusterTree(self,fileN=None):
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
      cPickle.dump(self.cluster(),outF)
      outF.close()
    

  def draw(self,idx=0):
    self.view().reset()
    clust = self.cluster()
    if clust is None or not clust:
      return
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    self.resizeCanvas((self.view().size().width()-10,self.view().size().height()-10))
    canv = self.canvas()
    ptColors = []
    if self._clusterActs:
      ptColors = ClusterVis.VisOpts.terminalColors
    try:
      ClusterVis.DrawClusterTree(clust,canv,canv.size,ptColors=ptColors,
                                 logScale=self.logScale(),
                                 tooClose=self._tooClose,
                                 showIndices=self.showLabels())
    except:
      QApplication.restoreOverrideCursor()
      qtUtils.error('Problems encountered drawing cluster tree.',exc_info=True)
      return
    self.view().update()
    QApplication.restoreOverrideCursor()

  def _doHighlights(self,pts):
    VisOpts = ClusterVis.VisOpts
    canv = self.canvas()
    for pt in pts:
      if hasattr(pt,'_drawPos'):
        pos = pt._drawPos
        canv.drawEllipse(pos[0]-VisOpts.highlightRad/2,pos[1]-VisOpts.highlightRad/2,
                         pos[0]+VisOpts.highlightRad/2,pos[1]+VisOpts.highlightRad/2,
                         VisOpts.highlightColor,0,VisOpts.highlightColor)
    self.view().update()
    self.view().canvas().update()

  def highlightClusters(self,which):
    clust = self.cluster()
    if clust is None or not clust:
      return
    pts = clust.GetPoints()
    if type(which) not in (types.ListType,types.TupleType):
      which = (which,)
    highlights = []
    for val in which:
      try:
        tmp = int(val)
      except ValueError:
        pass
      else:
        pt = clust.FindSubtree(tmp)
        if pt:
          highlights.append(pt)
    if highlights:
      self._doHighlights(highlights)
  def _pointFinder(self,which,methodName):
    missing = which[:]
    clust = self.cluster()
    if clust is None or not clust:
      return
    pts = clust.GetPoints()
    highlights = []
    for pt in pts:
      fn = getattr(pt,methodName)
      v = str(fn())
      if v in which:
        highlights.append(pt)
        try:
          missing.remove(v)
        except:
          qtUtils.logger.debug("could not remove point '%s'"%v,exc_info=True)
        
    return highlights,missing
  def highlightPoints(self,which):
    if type(which) not in (types.ListType,types.TupleType):
      which = (which,)
    highlights,missing=self._pointFinder(which,"GetName")
    if highlights:
      self._doHighlights(highlights)
    if missing:
      qtUtils.logger.info("Could not find points: %s"%(', '.join([str(x) for x in missing])))
        
  def updateTooClose(self,idx=0):
    val,ok = QInputDialog.getInteger(self.trUtf8("Enter new value for Too Close"),
                                  self.trUtf8("Too Close"),
                                  self._tooClose,-1,100)
    if ok:
      self._tooClose = val

  def tooClose(self):
    return self._tooClose
  def setTooClose(self,val):
    self._tooClose = val

  #
  # Callbacks and slots and stuff
  #
  def leftCanvasClick(self,view,evt):
    if self._cluster is None:
      return
    clickPt = view.inverseWorldMatrix().map(evt.pos())
    loc = clickPt.x(),clickPt.y()

    activeCluster = self._findClusterAtPoint(self._cluster,loc)
    if activeCluster:
      isTerminal = len(activeCluster.GetChildren())==0
      if isTerminal:
        for cb in self._pointClickCallbacks:
          cb(self,activeCluster)
      else:
        for cb in self._clusterClickCallbacks:
          cb(self,activeCluster)
      
      msg = '%s:  Idx: %d; Metric: %f'%(
        activeCluster.GetName(),
        activeCluster.GetIndex(),
        activeCluster.GetMetric())
      if activeCluster.GetData() not in ["",None]:
        msg += '; Data: %s'%(str(activeCluster.GetData()))
      self.statusBar().message(msg)
    else:
      self.statusBar().clear()
      

  def leftCanvasDblClick(self,view,evt):
    """ #DOC

    """
    if self._cluster is None:
      return
    
    clickPt = view.inverseWorldMatrix().map(evt.pos())
    loc = clickPt.x(),clickPt.y()

    activeCluster = self._findClusterAtPoint(self._cluster,loc)
    if activeCluster and len(activeCluster)>1:
      w = ClusterWindow(size=self.size())
      newCluster = copy.deepcopy(activeCluster)
      newCluster._ptLabels = self._clusterLabels
      newCluster._ptValues = self._clusterActs
      w.setCluster(newCluster)
      w.setLogScale(self.logScale())
      w.setTooClose(self.tooClose())
      w.setShowLabels(self.showLabels())
      w.setCaption('ClusterZoom: %s'%(newCluster.GetName()))
      w._dir = self._dir
      w.show()
      w.draw()
      self._children.append(w)
      
  def rightCanvasClick(self,view,evt):
    """ callback for right mouse clicks in canvases containing cluster
      trees
      (method of _GuiBase_)    

    launches a (possibly new) examples window using the information from
    the clicked node (if there is one)


    """
    if self._cluster is None:
      return
    clickPt = view.inverseWorldMatrix().map(evt.pos())
    loc = clickPt.x(),clickPt.y()

    activeCluster = self._findClusterAtPoint(self._cluster,loc)
    if activeCluster is not None:
      newWin = evt.state()&Qt.ShiftButton
      self.launchExamplesWin(view,activeCluster,newWin=newWin)

  def fileOpen(self):
    pass

  def fileSave(self):
    self.saveClusterTree(self._fileN)

  def fileSaveAs(self):
    self.saveClusterTree()

  def highlightPointsSlot(self):
    ptList=[x.strip() for x in str(self.ptsLine.text()).split(',')]
    self.highlightPoints(ptList)

  def highlightBatchPointsSlot(self):
    clip = qApp.clipboard()
    txt = clip.text()
    if txt.isEmpty():
      qtUtils.warning("no text found on the clipboard")
    txt = str(txt)
    ptList = []
    for line in txt.split('\n'):
      for entry in re.split(r'[\ \t,]+',line):
        entry = entry.strip()
        if entry:
          ptList.append(entry)
    self.highlightPoints(ptList)

  def highlightClustersSlot(self):
    ptList=[x.strip() for x in str(self.clustsLine.text()).split(',')]
    self.highlightClusters(ptList)



if __name__ == '__main__':
  from qtGui import Gui
  
  app,widg = Gui.Launcher(ClusterWindow,None)

  widg.loadPickledClusterTree(os.path.join(RDConfig.RDDataDir,'Ferromag','clusts.pkl'))
  app.exec_loop()
  widg.destroy(1)

