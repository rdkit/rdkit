# $Id: MolBrowserImpl.py 4684 2005-05-25 21:50:46Z glandrum $
#
#  Copyright (C) 2003-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for MolBrowsers

"""    
import RDConfig
from qt import *
from qtGui.GuiLib.forms.MolBrowser import MolBrowser as _Form
from qtGui.DbQueryWidgetImpl import insertQueryWidget,DbQueryWidget
from qtGui.DbConnWidgetImpl import insertConnWidget
from qtGui import GuiTable
from qtGui.GuiLib import MolTable
from qtGui.GuiLib.MolCanvas import MolCanvasView

from qtGui import qtUtils
import Chem
from Chem.Suppliers.DbMolSupplier import RandomAccessDbMolSupplier

class MolBrowser(_Form):
  """ Defines the class which is to be used to browse molecules

    The base widget is defined in forms.MolBrowser
  
  """
  def __init__(self,parent=None,initDir=''):
    _Form.__init__(self)
    self._dir = '.'

    self.tabWidget.setTabEnabled(self.dataPage,0)
    self.tabWidget.setTabEnabled(self.outputPage,0)

    # add the dbQuery widget to the db page
    self.queryWidget = insertQueryWidget(self.inputDbFrame,
                                         clickCallback=self.inputUpdate)
    self.queryWidget.show()
    self._drawTarget = None

    self.outputWidget = insertConnWidget(self.outputDbFrame,
                                         clickCallback=self.outputUpdate)
    self.outputWidget.allowTableEdits(1)
    self.outputWidget.show()
    self._drawTarget = MolCanvasView()
    self._drawTarget.initCanvas((300,300))
    self._drawTarget.setIcon(QPixmap(qtUtils.logoImageData))
    self._drawTarget.setCaption('Molecule')
    self._drawTarget.show()
    
    # add a MolTable:
    self.inDataTable = MolTable.insertMolTable(self.dataPage)

    # initialization:
    self.inputLoadButton.setEnabled(0)

    # FIX: get the file parsing stuff working:
    self.inputFileButton.setEnabled(0)
    self.inputFilenameBox.setEnabled(0)
    self.inputFilenameButton.setEnabled(0)
    self.outputFileButton.setEnabled(0)
    self.outputFilenameBox.setEnabled(0)
    self.outputFilenameButton.setEnabled(0)

    self.pdfButton = QPushButton(self.dataPage,"pdfButton")
    self.pdfButton.setText("Write PDF")
    self.dataPage.layout().addWidget(self.pdfButton)
    self.connect(self.pdfButton,SIGNAL("clicked()"),self.inDataTable.contentsToPDF)


  def loadFromDb(self,resultSet=None):
    """

    If provided, resultSet should be a model of DbResultSet
    
    """
    if resultSet is None:
      resultSet = self.queryWidget.getData(randomAccess=1)
    supplier = RandomAccessDbMolSupplier(resultSet)

    # we can now populate our MolTable:
    self.inDataTable.loadFromMolSupplier(supplier,drawTarget=self._drawTarget,
                                         processFunc=lambda x:Chem.Kekulize(x))
    self.tabWidget.setTabEnabled(self.dataPage,1)
    self.tabWidget.setTabEnabled(self.outputPage,1)
    self.tabWidget.setCurrentPage(self.tabWidget.indexOf(self.dataPage))

  def saveToDb(self,delim='\t',force=0):
    """

    """
    widg = self.outputWidget
    conn = widg.getConn()
    tblNames = [str(x).strip().upper() for x in conn.GetTableNames()]
    tgtName = str(widg.tableName()).upper()
    if not force and tgtName in tblNames:
      res = QMessageBox.warning(self,
                                "Table Exists",
                                "A table named:\n%s\nalready exists.\nContinuing will overrwrite it."%(str(widg.tableName())),
                                1,2)
      if res!=1:
        return
    from Dbase import DbUtils
    from cStringIO import StringIO
    # FIX: this is, how you say, hacky:
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    sio = StringIO(self.inDataTable.contentsToDelimText(delim=delim,skipCols=[],
                                                        replaceMolName=1))
    try:
      DbUtils.TextFileToDatabase(widg.dbName(),widg.tableName(),sio,
                                 user=widg.user(),password=widg.password(),
                                 delim=delim,maxColLabelLen=50)
    except:
      QApplication.restoreOverrideCursor()
      qtUtils.error('problems saving table to database',exc_info=True)
    else:
      QApplication.restoreOverrideCursor()
    
  #
  #  Signals and Slots
  #
  def inputUpdate(self):
    # we can now, in theory move on to the next stage, so enable the load button:
    self.inputLoadButton.setEnabled(1)

  def outputUpdate(self):
    # we can now, in theory move on to the next stage, so enable the load button:
    self.outputSaveButton.setEnabled(1)

  def inputLoadClicked(self):
    # load some molecules:
    if self.inputDbButton.isChecked():
      self.loadFromDb()

  def outputSaveClicked(self):
    # save some molecules:
    if self.outputDbButton.isChecked():
      self.saveToDb()
  
    
if __name__ == '__main__':
  from qtGui import Gui

  app,widg = Gui.Launcher(MolBrowser,'MolBrowser',Qt.WDestructiveClose)
  app.exec_loop()
  widg.destroy(1)
