# $Id: BuildCompositeDialogImpl.py 5175 2006-05-06 17:58:57Z glandrum $
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for BuildCompositeDialogs
#DOC: multiple composite stuff

"""    
import RDConfig
from ML import BuildComposite
#from ML.SVM import SVMClassificationModel as SVM
from qt import *
from qttable import *
from qtGui.GuiLib.forms.BuildCompositeDialog import CompositeBuilderDialog as _Form
from qtGui.DbQueryWidgetImpl import insertQueryWidget,DbQueryWidget
from qtGui.GuiLib import CompositeUtils
from qtGui.GuiTextViewer import GuiTextViewer
from qtGui import qtUtils

_nameCol = 0
_useCol = 1
_roleCol = 2
_boundsCol = 3
_numDescCols = 4

_roles=['Label','Info','Descriptor','Activity']
_roleItems = QStringList()
for thing in _roles: _roleItems.append(thing)

_nBounds=['-1','0','1','2','3']
_boundItems = QStringList()
for thing in _nBounds: _boundItems.append(thing)

class BuildCompositeDialog(_Form):
  """ Defines the class which is to be used to build composites

    The base widget is defined in forms.BuildCompositeDialog
  
  """
  def __init__(self,parent=None,initDir=''):
    _Form.__init__(self)
    self._dbDir = initDir
    self._colNames = None
    self._colTypes = None
    self._composites = None
    self.parent = parent
    self.mainTab.setTabEnabled(self.descriptorPage,0)

    # add the dbQuery widget to the db page
    fn = lambda x=self.mainTab,y=self.descriptorPage:x.setTabEnabled(y,1)
    self.queryWidget = insertQueryWidget(self.dbPage,
                                         clickCallback=fn)
    self.queryWidget.show()
    
  def updateTable(self):
    """ updates our table with descriptor info

    **Notes**

      - if _self.descs_table_ already exists, it will be overwritten
        here.

    """
    if not self._colNames:
      return
    nDescs = len(self._colNames)

    tbl = self.descs_table
    tbl.setNumCols(_numDescCols)
    tbl.setNumRows(nDescs)
    for i in range(nDescs):
      tbl.setText(i,_nameCol,self._colNames[i])
      chk = QCheckTableItem(tbl,None)
      tbl.setItem(i,_useCol,chk)
      chk.setChecked(1)
      role = QComboTableItem(tbl,_roleItems)
      tbl.setItem(i,_roleCol,role)
      bound = QComboTableItem(tbl,_boundItems)
      tbl.setItem(i,_boundsCol,bound)
      if self._colTypes[i] in ['integer','float']:
        role.setCurrentItem(2)
        bound.setCurrentItem(2)
      else:
        role.setCurrentItem(1)
        bound.setCurrentItem(0)
    # adjust the roles
    tbl.item(0,_roleCol).setCurrentItem(0)
    tbl.item(nDescs-1,_roleCol).setCurrentItem(3)
    # now the nBounds
    tbl.item(0,_boundsCol).setCurrentItem(1)
    tbl.item(nDescs-1,_boundsCol).setCurrentItem(1)
    for i in range(_numDescCols):
      tbl.adjustColumn(i)
    
  def refreshDescriptorContents(self):
    """  fills our table of input data using the database information
    provided 

    **Notes**

      - if _self.descs_table_ already exists, it will be overwritten
        here.

    """
    conn = self.queryWidget.getConn()
    tmp = conn.GetColumnNamesAndTypes(join=str(self.queryWidget.sqlJoin()),
                                      what=str(self.queryWidget.sqlWhat()))
    self._colNames = [x[0] for x in tmp]
    self._colTypes = [x[1] for x in tmp]
    self.updateTable()

  def _setRunDetails(self,details=None):
    """ INTERNAL USE ONLY:  builds a _CompositeRun.CompositeRun_
    object using the information in the dialog

    **Arguments**

     - details: (optional) a _CompositeRun.CompositeRun_ to be
       modified.  If this is not provided, one will be constructed
       using _ML.BuildComposite.SetDefaults()_

    **Returns**

      the details object

    """
    if details is None:
      details = BuildComposite.SetDefaults()
    details.dbName = self.queryWidget.dbName()
    details.dbUser = self.queryWidget.user()
    details.dbPassword = self.queryWidget.password()
    details.tableName = self.queryWidget.tableName()
    details.dbWhere = self.queryWidget.sqlWhere() 
    details.dbWhat =  self.queryWidget.sqlWhat()
    details.dbJoin =  self.queryWidget.sqlJoin()
    details.nModels = self.params_count.value()
    details.threshold = float(str(self.params_threshold.text()))
    details.splitRun = self.params_split.isChecked()
    details.splitFrac = float(str(self.params_splitFrac.text()))
    details.shuffleActivities = self.params_dataShuffle.isChecked()
    details.randomActivities = self.params_dataRandomize.isChecked()
    details.note = str(self.params_note.text())
    details.persistTblName = str(self.params_persisttable.text())
    details.detailedRes=1
    if self.params_multiplebuilds.isChecked():
      details.nRuns = int(str(self.params_multiplebuildcount.text()))

    if self.params_lock.isChecked():
      details.lockRandom=1
      v1 = int(str(self.params_lockV1.text()))
      v2 = int(str(self.params_lockV2.text()))
      details.randomSeed=(v1,v2)

    if self.params_autoBounds.isChecked():
      details.qBounds = []
      tbl = self.descs_table
      for i in range(tbl.numRows()):
        if tbl.item(i,_useCol).isChecked():
          itm = tbl.item(i,_boundsCol)
          nBnds = int(str(itm.currentText()))
        else:
          nBnds = -1
        details.qBounds.append(nBnds)
          
      details.qBoundCount = str(details.qBounds)
    else:
      details.qBounds = []
      details.qBoundCount = ''

    if self.params_filter.isChecked():
      filtVal = float(self.params_filterVal.value())
      filtFrac = float(str(self.params_filterFrac.text()))
      if self.params_filterModels.isChecked():
        details.modelFilterFrac = filtFrac
        details.modelFilterVal = filtVal
        details.filterFrac = 0.0
        details.filterVal = 0
      else:
        details.modelFilterFrac = 0.0
        details.modelFilterVal = 0
        details.filterFrac = filtFrac
        details.filterVal = filtVal
        

    # ------------------------
    # Tree options
    if self.params_treeRadio.isChecked():
      details.useTrees=1
      details.useKNN=0
      details.useSVM=0
      details.limitDepth = self.params_depth.value()
      details.lessGreedy = self.params_greedy.isChecked()
      details.pruneIt = self.params_prune.isChecked()
      details.recycleVars = self.params_recycle.isChecked()
      details.randomDescriptors = self.params_randomDescriptors.value()
    # ------------------------
    # KNN options
    elif self.params_knnRadio.isChecked():
      details.useTrees=0
      details.useKNN=1
      details.useNaiveBayes=0
      details.useSVM=0
      details.knnNeighs = self.params_knnKVal.value()
      details.knnDistFunc= self.params_knnMetric.currentText()
    # ------------------------
    # Naive Bayes options
    elif self.params_naiveBayesRadio.isChecked():
      details.useTrees=0
      details.useKNN=0
      details.useNaiveBayes=1
      details.useSVM=0
      if self.params_naiveBayesM.text():
        details.mEstimateVal = float(str(self.params_naiveBayesM.text()))
    # ------------------------
    # SVM options
    elif self.params_svmRadio.isChecked():
      details.useTrees=0
      details.useKNN=0
      details.useNaiveBayes=0
      details.useSVM=1
      val = str(self.params_svmKernel.currentText())
      if val not in SVM.kernels.keys():
        qtUtils.warning('Bad SVM kernel value.')
      else:
        details.svmKernel=SVM.kernels[val]
      val = str(self.params_svmType.currentText())
      if val not in SVM.machineTypes.keys():
        qtUtils.warning('Bad SVM type value.')
      else:
        details.svmType=SVM.machineTypes[val]
        
      val = str(self.params_svmWeights.text())
      if val:
        # FIX: dangerous eval
        vs = eval(val)
        nPts = len(vs)
        details.svmWeights = zip(range(nPts),vs)

      details.svmDegree = self.params_svmDegree.value()
      if self.params_svmEps.text():
        details.svmEps = float(str(self.params_svmEps.text()))

    txt = str(self.params_actBounds.text())
    if txt:
      try:
        qBounds = eval(txt)
      except:
        qtUtils.warning('Bad activity bounds value.',exc_info=True)
      else:
        if type(qBounds) not in [type([]),type(())]:
          newVal = '(%s,)'%(txt)
          try:
            qBounds = eval(newVal)
          except:
            qtUtils.warning('Bad activity bounds value.',exc_info=True)
            qBounds = None
          else:
            txt = newVal

        if type(qBounds) in [type([]),type(())]:
          details.activityBounds = qBounds
          details.activityBoundsVals = txt
    return details

  def buildIt(self):
    """ sets up the parameters and builds the composite

    **Returns**

      the composite
      
    **Notes**

      - Parameters are set using _self._setRunDetails()_

      - composite building is doing using our parent's
        _cbBuildCompositeFromDetails()_ method.

      - the fresh composites are stored in _self._composites_

    """
    self._runDetails = self._setRunDetails()
    cmd = BuildComposite.GetCommandLine(self._runDetails)
    qtUtils.logger.info('Command line:\n%s\n\n'%cmd)
    self._runDetails.cmd=cmd
    if self.parent:
      self._composites = self.parent.cbBuildCompositeFromDetails(self._runDetails)

    return self._composites
  
  def screenIt(self):
    """ screens our current composite against the current data

    **Notes**

      - Parameters are set using _self._setRunDetails()_

      - _self._composites[-1]_ is screened

      - composite screening and output window construction is handled
        by _CompositeUtils.ScreenCompositeFromDetails_

    """
    self._setRunDetails(details=self._runDetails)
    if self.parent:
      screenWin = self.parent._ciScreenWin
    else:
      screenWin = GuiTextViewer()
      self._screenWin = screenWin
    if len(self._composites):  
      CompositeUtils.ScreenCompositeFromDetails(self._runDetails,self._composites[-1],screenWin)

  def inspectIt(self):
    """ allows inspection of the composite... currently does nothing

    """
    pass
  
  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def toggleMult(self):
    """ callback

    """
    state = self.params_multiplebuilds.isChecked()
    self.params_multiplebuildcount.setEnabled(state)


  def toggleFilter(self):
    """ callback

    """
    state = self.params_filter.isChecked()
    self.params_filterVal.setEnabled(state)
    self.params_filterFrac.setEnabled(state)
    self.params_filterModels.setEnabled(state)
    
  def toggleLock(self):
    """ callback

    """
    state = self.params_lock.isChecked()
    self.params_lockV1.setEnabled(state)
    self.params_lockV2.setEnabled(state)

  def toggleSplit(self):
    """ callback

    """
    state = self.params_split.isChecked()
    self.params_splitFrac.setEnabled(state)

  def refreshClick(self):
    """ callback

    """
    self.refreshDescriptorContents()
    self.buttonBuild.setEnabled(1)

  def buildClick(self):
    """ callback

    """
    self.buildIt()
    if len(self._composites):
      self.buttonScreen.setEnabled(1)
      self.buttonInspect.setEnabled(1)

  def screenClick(self):
    """ callback

    """
    self.screenIt()
    
  def inspectClick(self):
    """ callback

    """
    self.inspectIt()
    
  def accept(self):
    """ close dialog with the ok button

    """
    _Form.accept(self)

  def _setTreeState(self,state):
    self.params_depth.setEnabled(state)
    self.params_autoBounds.setEnabled(state)
    self.params_greedy.setEnabled(state)
    self.params_recycle.setEnabled(state)
    self.params_prune.setEnabled(state)

  def _setKNNState(self,state):
    self.params_knnKVal.setEnabled(state)
    self.params_knnMetric.setEnabled(state)

  def _setNaiveBayesState(self,state):
    self.params_naiveBayesM.setEnabled(state)

  def _setSVMState(self,state):
    self.params_svmKernel.setEnabled(state)
    self.params_svmType.setEnabled(state)
    self.params_svmWeights.setEnabled(state)
    self.params_svmDegree.setEnabled(state)
    self.params_svmEps.setEnabled(state)
    self.params_svmNu.setEnabled(state)
    self.params_svmCoeff.setEnabled(state)
    self.params_svmGamma.setEnabled(state)
    self.params_svmCost.setEnabled(state)

  def selectTrees(self):
    state = self.params_treeRadio.isChecked()
    self._setTreeState(state)
    self._setKNNState(not state)
    self._setSVMState(not state)
    self._setNaiveBayesState(not state)
    
  def selectKNN(self):
    state = self.params_knnRadio.isChecked()
    self._setKNNState(state)
    self._setTreeState(not state)
    self._setSVMState(not state)
    self._setNaiveBayesState(not state)

  def selectNaiveBayes(self):
    state = self.params_naiveBayesRadio.isChecked()
    self._setNaiveBayesState(state)
    self._setSVMState(not state)
    self._setTreeState(not state)
    self._setKNNState(not state)
    
  def selectSVM(self):
    state = self.params_svmRadio.isChecked()
    self._setSVMState(state)
    self._setTreeState(not state)
    self._setKNNState(not state)
    self._setNaiveBayesState(not state)
    


if __name__ == '__main__':
  from qtGui import Gui

  app,widg = Gui.Launcher(BuildCompositeDialog)
  app.exec_loop()
  widg.destroy(1)
