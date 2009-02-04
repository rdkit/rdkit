# $Id$
#
#  Copyright (C) 2002-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for Similarity Searching 

"""    
from qt import *
import re
import os,os.path
from qtGui.GuiLib.forms.SimilaritySearch import SimilaritySearchWidget as _Form
from qtGui.GuiLib.FingerprintDbWidget import insertFingerprintDbWidget
from qtGui.GuiLib.SimilarityParamsImpl import SimilarityParamsWidget
from qtGui.GuiLib.MolEditDisplayWidget2 import MolEditDisplayWidget
from qtGui.GuiLib.MolCanvas import MolCanvasView
from qtGui.GuiLib import MolTable 
from qtGui import qtUtils

import Chem
from Chem.Fingerprints import MolSimilarity,FingerprintMols
from Chem.Suppliers.DbMolSupplier import ForwardDbMolSupplier as Supplier


class SimilaritySearchWidget(_Form):
  def __init__(self,*args,**kwargs):
    """

    """
    _Form.__init__(self,*args,**kwargs)
    self._initQueryWidget()
    self._initParamsWidget()
    self._initDbWidget()
    self._resultsCounter = 0

    self._drawTarget = MolCanvasView(caption="Similarity Results")
    self._drawTarget.initCanvas((300,300))
    self._drawTarget.show()
    self._drawTarget.hide()

    
  def _initQueryWidget(self):
    layout = self.queryPage.layout()
    if not layout:
      layout = QVBoxLayout(self.queryPage,2,2,'queryWidgetLayout')
    self.queryWidget = MolEditDisplayWidget(self.queryPage)
    self.queryWidget.insertMolUpdateCallback(self.enableSearch)
    layout.insertWidget(0,self.queryWidget)

  def _initParamsWidget(self):
    layout = self.paramsPage.layout()
    if not layout:
      layout = QVBoxLayout(self.paramsPage,2,2,'paramsWidgetLayout')
    self.simParams = SimilarityParamsWidget(self.paramsPage)
    layout.insertWidget(0,self.simParams)

  def _initDbWidget(self):
    # here's the base layout widget:
    self.dbQueryWidget = insertFingerprintDbWidget(self.databasePage,
                                                   clickCallback=self.enableSearch)


  def enableSearch(self,*args):
    """ checks to see if the search button should be enabled
    """
    if self.queryWidget.getMol() and \
       self.dbQueryWidget.dbName():
      self.searchButton.setEnabled(1)
    else:
      self.searchButton.setEnabled(0)
    return 0
  
  def processQueryFile(self,fileN=None):
    if fileN is None:
      fileN = str(QFileDialog.getOpenFileName(self._dir,
                                              "CDXML Files (*.cdxml);;Mol Files (*.mol);;All files (*.*)"))
    if fileN:
      self._dir,fileN = os.path.split(fileN)
      self._queryFileName = fileN
      try:
        self._queryData = open(os.path.join(self._dir,fileN),'r').read()
      except IOError:
        qtUtils.error('problems encountered reading from file %s'%(fileN))
        self._queryData = ''
      self.queryFile.setText(fileN)
      self._queryFormat = fileN.split('.')[-1]
      self.fileQueryRadio.setChecked(1)
      self.updateQuery()  
  

  def completeScreenDetails(self,details):
    """ fills the contents of the details structure passed in

    **Arguments**:

      - details: a _FingerprinterDetails_ instance

    """
    details.dbName = self.dbQueryWidget.dbName()
    details.tableName = self.dbQueryWidget.tableName()
    details.idName = self.dbQueryWidget.getColumnNames()[0]
    details.fpColName = self.dbQueryWidget.getFpColumn()

    # the Similarity Parameters instance konws how to do
    #  its own updating.
    self.simParams.updateFingerprinterDetails(details)

    
  def checkInput(self):
    if self.simParams.fragmentRadio.isChecked():
      if not self.dbQueryWidget.getFpColumn():
        ans = QMessageBox.warning(self,
                                  self.trUtf8("No Fingerprint Column"),
                                  self.trUtf8("No fingerprint column specified.\nThe similarity search may be very time consuming.\nAre you certain that you wish to continue?"),
                                  1,2)
        if ans!=1:
          return 0
      else:
        tmp = self.dbQueryWidget.getColumnNamesAndTypes()
    return 1     
        
  def doSearch(self):
    res = []
    molData = self.queryWidget.getMol()
    fmt = self.queryWidget.getFormat()
    smi = None
    if fmt == 'chemical/daylight-smiles':
      try:
        mol = Chem.MolFromSmiles(molData)
      except:
        mol = None
    elif fmt== 'chemical/mdl-molfile':
      # assume it's a mol file:
      try:
        mol = Chem.MolFromMolBlock(molData)
      except:
        qtUtils.logger.info('could not construct molecule',exc_info=True)
    if not mol:  
      qtUtils.warning("Could not process molecule, aborting search.\n")
      return None

    details = FingerprintMols.FingerprinterDetails()
    details.probeMol = mol
    self.completeScreenDetails(details)

    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    res = None
    if self.dbQueryWidget.getFpColumn():
      try:
        res = MolSimilarity.ScreenFromDetails(details,mol=details.probeMol)
      except TypeError:
        QApplication.restoreOverrideCursor()
        qtUtils.error('Bad fingerprint column type.\nPlease be sure that the fingerprint column contains pickled fingerprints.')
      except:
        QApplication.restoreOverrideCursor()
        qtUtils.warning('Problems encountered while searching',exc_info=True)
    else:
      try:
        fps = FingerprintMols.FingerprintsFromDetails(details)
        res = MolSimilarity.ScreenFingerprints(details,fps,mol=details.probeMol)
      except:
        QApplication.restoreOverrideCursor()
        qtUtils.warning('Errors encountered while searching',exc_info=True)
    QApplication.restoreOverrideCursor()

    return res
  
  def showResults(self,res):
    """ presents the results of the similarity search

    **Arguments**

      - res: a list of 2-tuples (id,score)
      
    To present the results reasonably, we'll need to query the db again to
    get some SMILES strings.

    """
    if not res or not len(res):
      return

    scores = {}
    for id,score in res:
      scores[id] = score
    textIds = ["'%s'"%(str(x)) for x in scores.keys()]
    
    self._resultsCounter += 1
    widg = QWidget(self.tabWidget,"resultsPage_%d"%(self._resultsCounter))
    self.tabWidget.insertTab(widg,"Results Set &%d"%(self._resultsCounter))

    tbl = MolTable.insertMolTable(widg,molCanvas=self._drawTarget)
    dbWidget = self.dbQueryWidget
    conn = dbWidget.getConn()
    names = dbWidget.getColumnNames()
    where = 'where %s in (%s)'%(names[0],','.join(textIds))
    d = conn.GetData(table=dbWidget.tableName(),
                     where=where,
                     fields=dbWidget.sqlWhat(),
                     join=dbWidget.sqlJoin(),
                     randomAccess=0)
    suppl = Supplier(d)

    tbl.loadFromMolSupplier(suppl,kekulize=True)
    # add the scores:
    scoreCol = 1
    tbl.insertColumns(scoreCol,1)
    tbl.horizontalHeader().setLabel(scoreCol,'SimScore')
    row = 0
    for mol in suppl:
      # FIX: maybe this assumption isn't so hot
      id = mol._fieldsFromDb[0]
      tbl.setText(row,scoreCol,'%0.4f'%(scores[id]))
      row += 1
    tbl.adjustColumn(1)
        
    
  
  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def patternFileButtonPressed(self):
    qtUtils.logger.debug('pattern File')

  def searchButtonClicked(self):
    if self.checkInput():
      res = self.doSearch()
      self.showResults(res)

if __name__ == '__main__':
  from qtGui import Gui

  app,widg = Gui.Launcher(SimilaritySearchWidget)
  app.exec_loop()
  widg.destroy(1)
