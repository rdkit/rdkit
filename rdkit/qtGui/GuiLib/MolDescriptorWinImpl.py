# $Id$
#
#  Copyright (C) 2002-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for MolDescriptorWins

"""    
from rdkit import RDConfig
from qt import *
from qttable import *
from rdkit.qtGui.GuiLib.forms.MolDescriptorWin import MolDescriptorWin as _Form
from rdkit.qtGui.DbQueryWidgetImpl import insertQueryWidget,DbQueryWidget
from rdkit.qtGui import GuiTable
from rdkit.qtGui.GuiLib import MolTable
from rdkit.qtGui.GuiLib.DescTable import DescTable
from rdkit.qtGui import qtUtils
from rdkit import Chem
from rdkit.Chem import AvailDescriptors
from rdkit.Chem.Suppliers.DbMolSupplier import RandomAccessDbMolSupplier
from rdkit.ML.Descriptors import MoleculeDescriptors

import re,types,os,sys


class MolDescriptorWin(_Form):
  """ Defines the class which is to be used to select/display descriptors

    The base widget is defined in forms.MolDescriptorWin
  
  """
  def __init__(self,parent=None,initDir=''):
    _Form.__init__(self)
    self._dbDir = initDir
    self._dir = '.'

    self.parent = parent
    self.mainTab.setTabEnabled(self.descPage,0)
    self.mainTab.setTabEnabled(self.dataPage,0)

    # add the dbQuery widget to the db page
    fn = lambda x=self.buttonRead:x.setEnabled(1)
    self.queryWidget = insertQueryWidget(self.sourcePage,
                                         clickCallback=fn)
    self.queryWidget.show()
    self.initSimpleTable()
    self.data_table=None
    self.desc_table=None
    self._descCalc=None

    self._dbMenu = QPopupMenu(self)
    self.menuBar().insertItem(self.trUtf8('&Database'),self._dbMenu,-1,2)
    self._dbMenu.insertItem(self.trUtf8('&Save Table'),self.saveActiveTable)

    self._drawTarget = None

  def setDrawTarget(self,tgt):
    self._drawTarget = tgt
  def drawTarget(self):
    return self._drawTarget

  def saveActiveTable(self):
    """ dumps the active table (selected page) to a database
    using that table's _tableToDb()_ method

    """
    # FIX: this needs to be generalized
    pg = self.mainTab.currentPage()
    if pg == self.dataPage:
      self.data_table.tableToDb()
    elif pg == self.descPage:
      self.desc_table.tableToDb()
    else:
      qtUtils.warning('cannot save on that page')
      
  def initSimpleTable(self):
    """ initializes the table of simple descriptors (a _DescTable_
    instance) 

    """
    descs = AvailDescriptors.descList
    nDescs = len(descs)

    layout = QHBoxLayout(self.simplePage,1,1,"simplepagelayout")
    tbl = DescTable(self.simplePage,"simpleDescriptors")
    layout.addWidget(tbl)
    self.simple_table = tbl
    tbl.setNumCols(3)
    tbl.setNumRows(nDescs)
    
    hdr = tbl.horizontalHeader()
    hdr.setLabel(0,self.trUtf8("Name"))
    hdr.setLabel(1,"")
    hdr.setLabel(2,self.trUtf8("Description"))
    
    for i in range(nDescs):
      name,fn = descs[i]
      tbl.setText(i,0,name)
      if fn.__doc__:
        doc = fn.__doc__.split('\n\n')[0].strip()
        doc = re.sub('\ *\n\ *',' ',doc)
        tbl.setText(i,2,doc)
        tbl.adjustRow(i)
      itm = QCheckTableItem(tbl,"")
      tbl.setItem(i,1,itm)
    tbl.adjustColumn(0)
    tbl.adjustColumn(1)
    tbl.adjustColumn(2)


  def loadInfoFromDb(self):
    """ uses the specified parameters to load a bunch of compounds
    from a database

    **Notes**:

      - The data are inserted into _self.data_table_, which will be
        created as a _MolTable_ if necessary.  Any existing data there
        will be mercilessly blown out.

      - If it was not already, the page containing the data table will
        be activated

     **Returns**

       the number of rows (compounds) added

    """
    w = self.queryWidget
    conn = w.getConn()
    data = conn.GetData(fields=w.sqlWhat(),where=w.sqlWhere(),join=w.sqlJoin())
    mols = RandomAccessDbMolSupplier(data)
    if not len(mols):
      return 0
    if self.data_table is None:
      self.data_table = MolTable.insertMolTable(self.dataPage)
    else:
      self.data_table.setNumRows(0)
      self.data_table.setNumCols(0)
    self.data_table.loadFromMolSupplier(mols,drawTarget=self.drawTarget(),
                                        includeCheckboxes=0)
    self.setUpdatesEnabled(1)
    self.mainTab.setTabEnabled(self.dataPage,1)
    return len(mols)

  def buildDescriptorCalc(self):
    descsToCalculate = self.simple_table.getSelected()
    descNames = [AvailDescriptors.descList[x][0] for x in descsToCalculate]
    self._descCalc = MoleculeDescriptors.MolecularDescriptorCalculator(descNames)
    return self._descCalc
  
  def genDescriptorVals(self):
    """ generates descriptor values for all the compounds in
    _self.data_table_
    
    **Notes**:

      - any descriptor calculations which have problems will result in
        a 666 being inserted into the table
        
      - The data are inserted into _self.desc_table_, which will be
        created as a _MolTable_ if necessary.  Any existing data there
        will be mercilessly blown out.

      - If it was not already, the page containing the descriptor
        table will be activated


    """
    self.buildDescriptorCalc()
    descNames = self._descCalc.GetDescriptorNames()
    nDescs = len(descNames)
    if self.desc_table is None:
      self.desc_table = MolTable.insertMolTable(self.descPage)
    else:
      self.desc_table.setNumRows(0)
      self.desc_table.setNumCols(0)
      self.desc_table._includeCheckboxes=0
      self.desc_table._molCol=0
      self.desc_table._checkboxCol=0
    tbl = self.desc_table
    nRows = self.data_table.numRows()
    nCols = self.data_table.numCols()+nDescs-self.data_table._includeCheckboxes
    tbl.setNumRows(nRows)
    tbl.setNumCols(nCols)
    hdr = tbl.horizontalHeader()

    # are there mols in that table?
    molCol = self.data_table._molCol
    if isinstance(self.data_table.item(0,molCol),MolTable.MolTableItem):
      gotMols = molCol+1
      hdr.setLabel(0,'Molecule')
    else:
      gotMols = 0
    if not gotMols:
      qtUtils.error('no molecules in table')
      return
    
    tbl.setUpdatesEnabled(0)
    pos = 1
    oHdr = self.data_table.horizontalHeader()
    for col in range(gotMols,self.data_table.numCols()):
      hdr.setLabel(pos,oHdr.label(col))
      pos += 1
    for col in range(nDescs):
      hdr.setLabel(pos,descNames[col])
      pos += 1
      
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    dlg = QProgressDialog('calculating','Cancel',nRows)
    dlg.setLabelText('Calculating Descriptors')
    dlg.setCancelButtonText('Cancel')
    for row in xrange(nRows):
      pos = 0
      oldItm = self.data_table.item(row,molCol)
      if hasattr(oldItm,'smiles') and hasattr(oldItm,'mol'):
        smi = oldItm.smiles()
        newI = MolTable.MolTableItem(tbl,QTableItem.Never,smi)
        mol = oldItm.mol()
        newI.setMol(mol)
        newI.setDrawTarget(self.drawTarget())
        tbl.setItem(row,0,newI)
        pos += 1
        for col in range(gotMols,self.data_table.numCols()):
          tbl.setText(row,pos,self.data_table.text(row,col))
          pos += 1

        try:
          descVals = self._descCalc.CalcDescriptors(mol)
        except:
          qtUtils.logger.info('Problems encountered with molecule: %s'%(smi),exc_print=True)
          descVals = [666]*len(self._descCalc.GetDescriptorNames())
        for val in descVals:
          tbl.setText(row,pos,str(val))  
          pos += 1
      dlg.setProgress(row)
      qApp.processEvents()
      if dlg.wasCancelled():
        break
    tbl.setUpdatesEnabled(1)
      
    self.mainTab.setTabEnabled(self.descPage,1)
    QApplication.restoreOverrideCursor()
        
  def saveCalculator(self,fileName):
    """ saves our calculator by pickling it into a file

    **Arguments**

      - fileName: the filename to use

    """
    calc = self.buildDescriptorCalc()
    if not calc:
      return
    else:
      calc.SaveState(fileName)
    

  #-------------------
  #
  #  Slots and callbacks
  # 
  #-------------------
  def readClick(self):
    """ callback for clicks on the read button

    """
    nMols = self.loadInfoFromDb()
    if nMols:
      self.buttonCalculate.setEnabled(1)
      
  def calcClick(self):
    """ callback for clicks on the calculate button

    """
    self.genDescriptorVals()

  def fileSaveAs(self):
    """ callback """
    # FIX: this directory nonsense sucks
    fileN = str(QFileDialog.getSaveFileName(self._dir,'Descriptor Calculator files (*.dsc);;All files (*.*)'))
    if fileN:
      self._dir = os.sep.join(os.path.split(fileN)[0:-1])
      self.saveCalculator(fileN)

  def fileClose(self):
    """ callback for selection of the menu item File->Close

    """
    self.hide()

  def editCopy(self,delim='\t'):
    """ callback for selection of the menu item Edit->Copy

      calls the table's _contentsToDelimText()_ method, if appropriate,
      and puts the results on the clipboard.

    """
    activePage = self.mainTab.currentPage()
    if activePage == self.descPage:
      tbl = self.desc_table
    elif activePage == self.dataPage:
      tbl = self.data_table
    else:
      tbl = None

    if tbl is not None:
      txt = tbl.contentsToDelimText(delim=delim,skipCols=[0])
      clip = qApp.clipboard()
      clip.setText(txt)
    else:
      qtUtils.warning("No table")

if __name__ == '__main__':
  from rdkit.qtGui import Gui

  app,widg = Gui.Launcher(MolDescriptorWin,'DescriptorWindow',Qt.WDestructiveClose)
  app.exec_loop()
  widg.destroy(1)
