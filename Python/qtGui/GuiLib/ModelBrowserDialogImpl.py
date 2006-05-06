# $Id$
#
#  Copyright (C) 2003-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for Model Browsers

"""    
import RDConfig
from qt import *
from qttable import *
from qtGui.GuiLib.forms.ModelBrowserDialog import ModelBrowserDialog as _Form
from qtGui.DbConnWidgetImpl import insertConnWidget
from qtGui import GuiTable,qtUtils
import cPickle

# column numbers in the model table
_idCol      = 0
_useCol     = 1
_noteCol    = 2
_holdoutCol = 3
_overallCol = 4
_tableCol = 5
_commandCol = 6
_colNames=['ID','Use','Note','Holdout Error','Overall Error','Table','Command']

def _stringify(inVal):
  """ #DOC
  """
  try:
    val = float(inVal)
  except:
    val = str(inVal)
  else:
    val = '%.4f'%(val)
  return val

class ModelBrowserDialog(_Form):
  """ class used to handle selecting lists of models from a database
  
    The base widget is defined in forms.ModelBrowserDialog

  """
  def __init__(self,parent=None,initDir=''):
    _Form.__init__(self,parent)
    self._dbDir = initDir
    self._composites = []

    # add the dbQuery widget to the db page
    fn = lambda z=self.buttonRefresh:z.setEnabled(1)
    self.queryWidget = insertConnWidget(self.dbPage,
                                      clickCallback=fn)
    self._sqlWhat = 'Note,holdout_error,overall_error,tablename,cmd'
    self.queryWidget.show()
    self.tabWidget.setTabEnabled(self.modelsPage,0)

    self._errorTable = None
    self._rejTable = None

    # FIX: initialize this
    if parent and hasattr(parent,'ciAddComposite'):
      self.ciAddComposite = parent.ciAddComposite
    else:
      self.ciAddComposite = None

    self._initModelsTable()
  def _initModelsTable(self):
    """ INTERNAL USE ONLY
    """
    self.modelsTable = GuiTable.insertTable(self.modelsPage,GuiTable.GuiTable)
    self.modelsTable.setNumCols(len(_colNames))
    hdr = self.modelsTable.horizontalHeader()
    for i in range(len(_colNames)):
      hdr.setLabel(i,_colNames[i])
    

  def refreshIt(self):
    widget = self.queryWidget
    self._conn = widget.getConn()
    models = []
    nModels = 0
    try:
      models = self._conn.GetData(fields=self._sqlWhat)
      nModels = len(models)
    except:
      qtUtils.error('Database query did not find any models.\nPlease check that the database and table specifications are correct.',exc_info=True)
    if nModels:
      self._origNumRows = nModels
      self.modelsTable.setNumRows(0)
      self.modelsTable.setNumRows(nModels)
      for i in range(nModels):
        model = models[i]
        self.modelsTable.setText(i,_idCol,str(i+1))
        self.modelsTable.setItem(i,_useCol,QCheckTableItem(self.modelsTable,None))
        self.modelsTable.setText(i,_noteCol,str(model[0]))
        self.modelsTable.setText(i,_holdoutCol,_stringify(model[1]))
        self.modelsTable.setText(i,_overallCol,_stringify(model[2]))
        self.modelsTable.setText(i,_tableCol,str(model[3]))
        self.modelsTable.setText(i,_commandCol,str(model[4]))
                                  
      for i in range(_commandCol+1):
        self.modelsTable.adjustColumn(i)
        
      self.tabWidget.setTabEnabled(self.modelsPage,1)
      if self.ciAddComposite is not None:
        self.buttonTransfer.setEnabled(1)

  def transferIt(self):
    indicesToKeep=[]
    notesToKeep={}
    for i in range(self.modelsTable.numRows()):
      if self.modelsTable.item(i,_useCol).isChecked():
        indicesToKeep.append(int(str(self.modelsTable.text(i,_idCol))))
        notesToKeep[str(self.modelsTable.text(i,_noteCol))] = 1
    qtUtils.logger.debug('keep: %s'%(str(indicesToKeep)))

    what = 'model'
    curs = self._conn.GetCursor()
    # FIX: implementation specific
    tbl = self._conn.tableName
    curs.execute('select %s from %s'%(what,tbl))
    nDone = 0
    dlg = QProgressDialog('Extracting Models','',self._origNumRows)
    for i in range(self._origNumRows):
      thing = curs.fetchone()
      if i+1 in indicesToKeep:
        mdl = cPickle.loads(str(thing[0]))
        self.ciAddComposite(mdl,'Composite %d'%(i+1))
        nDone += 1
        if nDone == len(indicesToKeep):
          break
      dlg.setProgress(i)  

  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def transferClick(self):
    """ callback """
    self.transferIt()
    
  def refreshClick(self):
    """ callback """
    self.refreshIt()
    


