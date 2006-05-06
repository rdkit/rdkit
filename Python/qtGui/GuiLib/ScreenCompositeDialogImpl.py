#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for ScreenCompositeDialogs

"""    
import RDConfig
from qt import *
from qttable import *
from qtGui.GuiLib.forms.ScreenCompositeDialog import ScreenCompositeDialog as _Form
from qtGui.DbQueryWidgetImpl import insertQueryWidget
from qtGui.GuiLib import CompositeUtils
from qtGui.GuiTextViewer import GuiTextViewer
from qtGui import GuiTable,qtUtils
from ML import ScreenComposite
from ML.Data import DataUtils
from Dbase.DbConnection import DbConnect

class ScreenCompositeDialog(_Form):
  """ class used to handle composite screening
  
    The base widget is defined in forms.ScreenCompositeDialog

  """
  def __init__(self,parent=None,initDir=''):
    _Form.__init__(self,parent)
    self._dbDir = initDir
    self._colNames = None
    self._composite = None
    self._screenResWin = None
    self._runDetails = None
    # add the dbQuery widget to the db page
    fn = lambda y=self.buttonScreen,z=self.buttonScreenAll:(y.setEnabled(1),z.setEnabled(1))
    self.queryWidget = insertQueryWidget(self.dataPage,
                                         clickCallback=fn)
    self.queryWidget.show()
    self.tabWidget.setTabEnabled(self.filePage,0)

    self._errorTable = None
    self._rejTable = None

    self._composWin = None

  def setCompositeWindow(self,win):
    """  #DOC
      win should be a _CompositeWindow_

    """
    
    self._composWin = win
  def compositeWindow(self):
    """  #DOC  """
    return self._composWin
    
  def setComposite(self,composite,name=None):
    """ setter for composite attribute """
    self._composite = composite
    if name is not None:
      self.options_modelName.setText(name)
  def composite(self):
    """ getter for composite attribute """
    return self._composite
  
  def runDetails(self):
    """  #DOC   """
    return self._runDetails
  
  def setRunDetails(self):
    """ builds a _CompositeRun.CompositeRun_ object using the
    information in the dialog

    **Returns**

      the details object

    **Notes**
    
     - if our _runDetails attribute has not been set, a
       _CompositeRun.CompositeRun_ will be constructed using
       _ML.ScreenComposite.SetDefaults()_

    """

    if self._runDetails is None:
      self._runDetails = ScreenComposite.SetDefaults()
    details = self._runDetails
    widg = self.queryWidget
    details.dbName = widg.dbName()
    details.tableName = widg.tableName()
    details.dbUser = widg.user()
    details.dbPassword = widg.password()
    details.dbWhat = widg.sqlWhat()
    details.dbWhere = widg.sqlWhere()
    details.dbJoin = widg.sqlJoin()
    details.threshold = float(str(self.options_thresholdEdit.text()))
    details.errorEstimate = self.options_oobCheck.isChecked()
    if self.options_sourceTraining.isChecked():
      details.doTraining = 1
      details.doHoldout = 0
    elif self.options_sourceHoldout.isChecked():
      details.doHoldout = 1
      details.doTraining = 0
    else:
      details.doHoldout=0
      details.doTraining=0

    info = details.GetDataSetInfo()
    pickleCol=-1
    col = 0
    for cName,cType in info:
      if cType=='binary': 
        pickleCol=col
      col+=1
    if pickleCol>0: details.pickleCol=pickleCol-1  
    return details
  
  def loadTable(self,colNames,table,votes,dataSet):
    """ loads up a results table
#DOC
    **Arguments**

      - colNames: a sequence with the names of the columns

      - table: the table to load

      - votes: a sequence of 4-tuples:

         1) ignored (this is likely the actual value)

         2) predicted value

         3) confidence

         4) example (a sequence with the actual data point)

    **Notes**

      - the resulting table will contain _len(colNames)_ + 2 columns,
        where the extra 2 columns are the prediction and accuracy
        

    """
    model = self.composite()
    try:
      seed = model._randomSeed
    except AttributeError:
      pass
    else:
      DataUtils.InitRandomNumbers(seed)
    trainIdx,testIdx = ScreenComposite.PrepareDataFromDetails(model,
                                                              self._runDetails,
                                                              dataSet)
    nms = colNames[:]
    nms.extend(['Pred','Conf'])
    data = [dataSet[testIdx[x[3]]]+[x[1],x[2]] for x in votes]
    GuiTable.loadTable(table,data,nms)

      
  def updateVotePages(self,badVotes,noVotes):
    """ updates our error and rejected tables from vote information
#DOC
    **Arguments**

      - badVotes/noVotes: sequences of 4-tuples:

         1) actual value

         2) predicted value

         3) confidence

         4) example (a sequence with the actual data point)

      
    **Notes**

      - if our _errorTable_ or _rejTable_ already exist, they will be
        destroyed here.  Otherwise they will be created (and the
        corresponding pages inserted into our notebook

      - if _badVotes_ / _noVotes_ is empty, the corresponding page
        will not be created.

    """
    data = self._runDetails.GetDataSet()
    tmp = data.GetVarNames()
    self._colNames = [x.strip() for x in tmp]

    # FIX: remove bad/rejected page if there are no examples
    if badVotes:
      if self._errorTable is None:
        pg = QWidget(self)
        tbl = GuiTable.insertTable(pg,GuiTable.GuiTable)
        self._errorTable = tbl,pg
        self.tabWidget.addTab(pg,'Errors')
      else:
        tbl,pg = self._errorTable
      self.loadTable(self._colNames,tbl,badVotes,data)
      self.tabWidget.showPage(pg)
      
    if noVotes:
      if self._rejTable is None:
        pg = QWidget(self)
        tbl = GuiTable.insertTable(pg,GuiTable.GuiTable)
        self._rejTable = tbl,pg
        self.tabWidget.addTab(pg,'Rejected')
      else:
        tbl,pg = self._rejTable
      newPg = QWidget(self)
      self.loadTable(self._colNames,tbl,noVotes)
      self.tabWidget.showPage(pg)

    
  def screenIt(self,screenAll=0):
    """ screens our current composite using our details

    **Notes**

      - we call out to our _self.setRunDetails() to get screening
        parameters 

      - actual screening is done by
        _CompositeUtils.ScreenCompositeFromDetails()_
      
    """
    self.setRunDetails()
    if self._screenResWin is None:
      self._screenResWin = GuiTextViewer()
    if self.options_errorCheck.isChecked():
      badVotes = []
      #self._runDetails.errorAnalysis=1
    else:
      badVotes = None
    if self.options_rejectionCheck.isChecked():
      noVotes = []
      #self._runDetails.errorAnalysis=1
    else:
      noVotes = None
    if not screenAll:
      CompositeUtils.ScreenCompositeFromDetails(self._runDetails,self.composite(),
                                                self._screenResWin,
                                                badVotes=badVotes,noVotes=noVotes)
      self.updateVotePages(badVotes,noVotes)
    else:
      cWin = self.compositeWindow()
      composites = cWin.allComposites()
      CompositeUtils.ScreenCompositeFromDetails(self._runDetails,composites,
                                                self._screenResWin)



  def refreshIt(self):
    """ updates our current composite from our parent

      #DOC
    
    """
    composWin = self.compositeWindow()
    if composWin is not None:
      nm = composWin.tabWidget.label(composWin.tabWidget.currentPageIndex())
      try:
        compos = composWin.activeComposite()
      except AttributeError:
        qtUtils.warning('no composite there')
      else:
        self.setComposite(compos)
        self.options_modelName.setText(nm)
        
    
    

  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def screenClick(self):
    """ callback """
    self.screenIt()

  def screenAllClick(self):
    """ callback """
    self.screenIt(screenAll=1)
    
  def refreshClick(self):
    """ callback """
    self.refreshIt()
    


