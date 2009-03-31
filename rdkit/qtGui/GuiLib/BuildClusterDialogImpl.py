# $Id$
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for BuildClusterDialogs

#DOC this entire file

"""
from rdkit import RDConfig
from rdkit.ML.Cluster import Murtagh
from qt import *
from qttable import *
from rdkit.qtGui.GuiLib.forms.BuildClusterDialog import ClusterBuilderDialog as _Form
from rdkit.qtGui.GuiLib.FingerprintDbWidget import insertFingerprintDbWidget
from rdkit.qtGui import GuiTable,qtUtils
from rdkit.qtGui.GuiLib.ClusterWindow import ClusterWindow
from rdkit.ML.Cluster import Murtagh,Resemblance,Standardize
from rdkit.Chem.Fingerprints import ClusterMols
from rdkit import DataStructs
import cPickle,types


#FIX: we don't necessarily want to have to do this:
from Numeric import *

class BuildClusterDialog(_Form):
  """ Defines the class which is to be used to build cluster trees

    The base widget is defined in forms.BuildClusterDialog
  
  """
  def __init__(self,parent=None,initDir=''):
    _Form.__init__(self)
    self._dbDir = initDir
    self.mainTab.setTabEnabled(self.dataPage,0)
    self.mainTab.setTabEnabled(self.colsPage,0)

    # add the dbQuery widget to the db page
    #fn = lambda x=self.mainTab,y=self.dataPage:x.setTabEnabled(y,1)
    self.dbQueryWidget = insertFingerprintDbWidget(self.dbPage,
                                                 clickCallback=self.dbFileClick)
    self.dbQueryWidget.show()

    self.data_table = GuiTable.insertTable(self.dataPage,GuiTable.GuiTable)    
    self._initMethods()

    # it's stupid that we need to do this, but I can't figure out a work-around
    #  at the moment and this does work
    self.paramsTypeFingerprints.setOn(1)
    self.paramsTypeFingerprints.setOn(0)
    self.paramsTypeDescriptors.setOn(1)

    
    self._clusterTree = None
    
    
  def _initMethods(self):
    """  #DOC

    """
    self.params_algorithmCombo.clear()
    self._algorithms = {}
    algorithms = Murtagh.methods
    for name,constant,desc in algorithms:
      self._algorithms[name] = (constant,desc)
      self.params_algorithmCombo.insertItem(name)

    self.params_distanceCombo.clear()
    self._metrics = {}
    metrics = Resemblance.methods
    for name,func,desc in metrics:
      self._metrics[name] = (func,desc)
      self.params_distanceCombo.insertItem(name)

    self.params_standardCombo.clear()
    self._standards = {}
    methods = Standardize.methods
    for name,func,desc in methods:
      self._standards[name] = (func,desc)
      self.params_standardCombo.insertItem(name)


    self._fpMetrics = {}
    for name,func,desc in DataStructs.similarityFunctions:
      self._fpMetrics[name] = (func,desc)
      self.params_fpMetricCombo.insertItem(name)
      
    

  def updateTable(self):
    """ updates our table of input data

    **Notes**

      - if _self.data_table_ already exists, it will be overwritten
        here.

    """
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    GuiTable.loadTable(self.data_table,
                       self.dbQueryWidget.getData(),
                       self._colNames)
    QApplication.restoreOverrideCursor()
                       
      
  def refreshDataContents(self):
    """  fills our table of input data using the database information
    provided 

    **Notes**

      - if _self.data_table_ already exists, it will be overwritten
        here.

    """
    tmp = self.dbQueryWidget.getColumnNamesAndTypes()
    self._colNames = [x[0] for x in tmp]
    self._colTypes = [x[1] for x in tmp]
    self.updateTable()
    self.mainTab.setTabEnabled(self.dataPage,1)

  def clusterFromDataTable(self,show=1,labelTree=1):
    """ sets up the parameters and builds the cluster tree using
    descriptor-based clustering
      #DOC

    **Returns**

      the cluster tree
      
    **Notes**

    """
    rawD = self.data_table.contentsAsList()

    labels = self.params_labelCheck.isChecked()
    acts = self.params_activityCheck.isChecked()
    methodId,methodDesc = self._algorithms[str(self.params_algorithmCombo.currentText())]
    stdFunc,stdDesc = self._standards[str(self.params_standardCombo.currentText())]
    distFunc,distDesc = self._metrics[str(self.params_distanceCombo.currentText())]

    if labels or acts:
      if labels:
        startP = 1
        labels = [x[0] for x in rawD]
      else:
        startP = 0

      if acts:
        endP = -1
        acts = [x[-1] for x in rawD]
        bounds = ()
        if self.params_quantizeCheck.isChecked():
          txt = str(self.params_activityQuant.text())
          if txt:
            try:
              bounds = eval(txt)
            except:
              qtUtils.warning('Activity quantization bounds cannot be converted to a meaningful value.')
            else:
              if type(bounds) not in (types.ListType,types.TupleType):
                bounds = (bounds,)
        if bounds:
          for i in range(len(acts)):
            act = acts[i]
            placed = 0
            bound = 0
            while not placed and bound < len(bounds):
              if act < bounds[bound]:
                acts[i] = bound
                placed = 1
              else:
                bound += 1
            if not placed:
              acts[i] = bound
      else:
        endP = len(rawD[0])+1
      try:
        data = array([x[startP:endP] for x in rawD])
      except ValueError:
        qtUtils.error('Unable to copy data table.\nThe problem is probably due to one or more of the following:\n  1) Attempting to normally cluster fingerprint data.\n  2) The "Activities Present?" box is checked for a dataset where no activities are present.')
        return
        
    else:
      data = array(rawD)
      #data = rawD


    nPts = len(data)
    
    # FIX:  the clustering code should maybe be in the mixin class?
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    try:
      qtUtils.logger.info('standardizing')
      data = stdFunc(data)
      qtUtils.logger.info('distance')
      dists = distFunc(data)
      data = None
      qtUtils.logger.info('clustering')
      res = Murtagh.ClusterData(dists,nPts,methodId,isDistData=1)
      tree = res[0]

      if labelTree:
        pts = tree.GetPoints()
      if labels:
        tree._ptLabels = labels
        if labelTree:
          for pt in pts:
            idx = pt.GetIndex()-1
            pt.SetName(labels[idx])
      if acts:
        tree._ptValues = acts
        if labelTree:
          for pt in pts:
            idx = pt.GetIndex()-1
            try:
              pt.SetData(int(acts[idx]))
            except ValueError:
              pass

      self._clusterTree = tree
    except:
      qtUtils.logger.info('problems encountered while clustering',exc_info=True)
      self._clusterTree=None  
    QApplication.restoreOverrideCursor()

    if show:
      w = ClusterWindow()
      w.show()
      w.setCluster(self._clusterTree)
      w.setDistances(dists)
      w.draw()
    
    return self._clusterTree

  def clusterFromFingerprints(self,show=1,labelTree=1):
    """ sets up the parameters and builds the cluster tree
      using fingerprint-based clustering

    **Returns**

      the cluster tree
      
    """
    haveFps = 1
    if not self.dbQueryWidget.getFpColumn():
      haveFps=0
      qtUtils.information("No fingerprint column specified.\nClustering cannot be carried out without fingerprints.")
      #ans = QMessageBox.warning(self,
      #                          self.trUtf8("No Fingerprint Column"),
      #                          self.trUtf8("No fingerprint column specified.\nThe clustering may be very time consuming.\nAre you certain that you wish to continue?"),
      #                          1,2)
      #if ans!=1:
      return None

    bounds = None
    if self.params_quantizeCheck.isChecked():
      txt = str(self.params_activityQuant.text())
      if txt:
        try:
          bounds = eval(txt)
        except:
          qtUtils.warning('Activity quantization bounds cannot be converted to a meaningful value.')
        else:
          if type(bounds) not in (types.ListType,types.TupleType):
            bounds = (bounds,)

    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    # we're going to ignore the data currently in the table and just pull it from the
    #  db again:
    colNames = [x.upper() for x in self.dbQueryWidget.getColumnNames()]
    fpName = self.dbQueryWidget.getFpColumn().upper()
    pruneNames = colNames[:]
    try:
      fpIdx = colNames.index(fpName)
    except:
      QApplication.restoreOverrideCursor()
      qtUtils.error('could not find fingerprint column (%s) in database with columns:\n%s'%(fpName,str(colNames)),
                    exc_info=True)
      return None
    while pruneNames.count(fpName):
      pruneNames.remove(fpName)
    toTake = []
    haveLabels=0
    haveActs=0
    if self.params_labelCheck.isChecked():
      if not len(pruneNames):
        haveLabels=False
        nameCol=None
      else:
        haveLabels=True
        nameCol = pruneNames[0]
        toTake.append(colNames.index(nameCol))
        while pruneNames.count(nameCol):
          pruneNames.remove(nameCol)
    if self.params_activityCheck.isChecked():
      if not len(pruneNames):
        haveActs=False
        actCol=None
      else:
        haveActs=True
        actCol = pruneNames[-1]
        toTake.append(colNames.index(actCol))
        while pruneNames.count(actCol):
          pruneNames.remove(actCol)
    toTake.append(fpIdx)
    data = self.dbQueryWidget.getData()
    pts = []
    usePickle=False
    for pt in data:
      pt = [pt[x] for x in toTake]
      pos = 0
      label = None
      act = None
      if haveLabels:
        label = pt[pos]
        pos += 1
      if haveActs:
        act = pt[pos]
        pos += 1
      fpTxt = str(pt[pos])
      if not usePickle:
        conv = DataStructs.ExplicitBitVect
      else:
        conv = cPickle.loads
      try:
        fp = conv(fpTxt)
      except cPickle.UnpicklingError:
        fp = None
        usePickle = False
        try:
          fp = conv(fpTxt)
        except:
          fp = None
          qtUtils.logger.debug('problems creating a cluster w/o pickle',exc_info=True)
      except:
        fp = None
        qtUtils.logger.debug('problems creating a cluster with pickle',exc_info=True)
        
      if not fp:
        QApplication.restoreOverrideCursor()
        qtUtils.error('Could not construct fingerprint for point %s.\n\nAborting clustering.\n'%(label))
        return None

      if bounds and act is not None:
        placed=0
        bound=0
        while not placed and bound < len(bounds):
          if act < bounds[bound]:
            act = bound
            placed = 1
          else:
            bound += 1
        if not placed:
          act = bound
      pts.append((label,fp,act))
    metricNm = str(self.params_fpMetricCombo.currentText())
    metricFn,desc = self._fpMetrics[metricNm]
    methodId,methodDesc = self._algorithms[str(self.params_algorithmCombo.currentText())]

    try:
      clusterTree,dists = ClusterMols.ClusterPoints(pts,metricFn,methodId,
                                                    haveLabels=haveLabels,
                                                    haveActs=haveActs,
                                                    returnDistances=True)
    except:
      clusterTree=None
      QApplication.restoreOverrideCursor()
      qtUtils.error('problems encountered clustering data',exc_info=True)
    else:
      QApplication.restoreOverrideCursor()

    if show and clusterTree:
      self._clusterTree = clusterTree
      w = ClusterWindow()
      w.setCluster(clusterTree)
      w.setDistances(dists)
      w.show()
      w.draw()
      
    return clusterTree
  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def refreshClick(self):
    """ callback

    """
    self.refreshDataContents()
    self.buttonBuild.setEnabled(1)

  def buildClick(self):
    """ callback

    """
    if self.paramsTypeDescriptors.isOn():
      self.clusterFromDataTable()
    else:
      self.clusterFromFingerprints()

  def accept(self):
    """ close dialog with the ok button

    """
    _Form.accept(self)

  def dbFileClick(self):
    self.refreshButton.setEnabled(1)

  def quantizeChecked(self):
    if self.params_quantizeCheck.isChecked():
      self.params_activityQuant.setEnabled(1)
    else:
      self.params_activityQuant.setEnabled(0)

  def descriptorsChecked(self,v):
    if v:
      self.params_standardCombo.setEnabled(1)
      self.params_distanceCombo.setEnabled(1)
    else:
      self.params_standardCombo.setEnabled(0)
      self.params_distanceCombo.setEnabled(0)
  def fingerprintsChecked(self,v):
    if v:
      self.params_fpMetricCombo.setEnabled(1)
    else:
      self.params_fpMetricCombo.setEnabled(0)

      
    
if __name__ == '__main__':
  from rdkit import RDLogger
  RDLogger.EnableLog('rdApp.debug')
  from rkdit.qtGui import Gui

  app,widg = Gui.Launcher(BuildClusterDialog)
  app.exec_loop()
  widg.destroy(1)
