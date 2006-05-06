#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for SupplierWindows (for editing supplier
    node properties)

  NOTE: this has been put together with the greatest of speed, so
   there ain't much in the way of docs.

"""    
from qt import *
import re
import os,os.path
from forms.LoaderWidget import LoaderWidget as _Form
from qtGui.DbQueryWidgetImpl import insertQueryWidget,DbQueryWidget
from qtGui.GuiLib.MolTable import insertMolTable

sdFileDataExpr = re.compile(r'>\s+<(.*)>\s+(.+)\n\n')

class LoaderWidget(_Form):
  def __init__(self,suppObj=None,initDir=''):
    """

    """
    _Form.__init__(self)

    self.guts = suppObj
    self._dir=initDir
    self._fileName = ''
    self._dbName = ''

    self.queryWidget = insertQueryWidget(self.queryPage,clickCallback=self.refreshContents)
    self.queryWidget.show()

    self.supplierTable = insertMolTable(self.contentsPage,'Contents')
    self.supplierTab.setTabEnabled(self.contentsPage,0)
    self.supplierTable.setAllowColInsertion(1)
    self.supplierTable.setAllowRowInsertion(1)
    self.supplierTable.setAllowColRename(1)
    self.supplierTable.setAllowRowRename(1)

    
  def getLen(self):
    res = 0
    if self.info_molButton.isOn():
      res = 1
    elif self.info_fileButton.isOn() and self._fileName != '':
      fName = os.path.join(self._dir,self._fileName)
      if self._fmt == 'sdf':
        res = int(os.popen(r"egrep -c '\$\$\$\$' %s"%(fName),'r').read())
      elif self._fmt == 'smi':
        res = int(os.popen(r"egrep -vc '^#' %s"%(fName),'r').read())
      elif self._fmt == 'tdt':
        res = int(os.popen(r"egrep -c '^\$SMI' %s"%(fName),'r').read())
      elif self._fmt == 'txt':
        res = int(os.popen(r"egrep -vc '^#' %s"%(fName),'r').read())
    elif self.info_dbButton.isOn() and self.queryWidget.dbName():
      res = len(self.getDbData())
      
    return res

  def getDbData(self,fieldNames=[]):
    conn = self.queryWidget.getConn()
    what = self.queryWidget.sqlWhat()
    where = self.queryWidget.sqlWhere()
    join = self.queryWidget.sqlJoin()
    fieldNames += conn.GetColumnNames(what=what,join=join)
    print '\t',repr(fieldNames)
    return conn.GetData(fields=what,where=where,join=join)

  def updateTable(self):
    if self.info_dbButton.isOn() and self.queryWidget.dbName():
      fieldNames = []
      data = self.getDbData(fieldNames=fieldNames)
      nData = len(data)
      if nData:
        nVals = len(data[0])
        if nVals != len(fieldNames):
          print 'mismatch: ',nVals,len(fieldNames)
          return
    elif self.info_fileButton.isOn() and self._fileName != '':
      fName = os.path.join(self._dir,self._fileName)
      if self._fmt == 'sdf':
        raw = sdFileDataExpr.findall(open(fName,'r').read())
        fieldNames = []
        for name,val in raw:
          if name in fieldNames:
            break
          else:
            fieldNames.append(name)
        nVals = len(fieldNames)    
        nData = len(raw) / nVals
        data = [None]*nData
        for row in range(nData):
          data[row]=['']*nVals
          for col in range(nVals):
            data[row][col] = raw[row*nVals+col][1]
    elif self.info_molButton.isOn():
      smiles = str(self.mol_smilesBox.text())
      nData = 1
      nVals = 4
      fieldNames = ['Smiles','Desc-1','Desc-2','Desc-3']
      data = [[smiles,0,0,0]]
    if nData:
      self.supplierTab.setTabEnabled(self.contentsPage,1)

      tbl = self.supplierTable
      tbl.setNumCols(nVals)
      tbl.setNumRows(nData)
      hdr = tbl.horizontalHeader()
      for col in range(nVals):
        #tbl.setText(0,col,fieldNames[col])
        hdr.setLabel(col,fieldNames[col])
      for row in range(nData):
        for col in range(nVals):
          tbl.setText(row,col,str(data[row][col]))
    else:      
      self.supplierTab.setTabEnabled(self.contentsPage,0)
      
    
    
  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def fileClick(self):
    if self.file_sdfButton.isOn():
      fmt = 'sdf'
      fmtStr = 'SD Files (*.sdf);;All files (*.*)'
    elif self.file_smiButton.isOn():
      fmt = 'smi'
      fmtStr = 'Smiles Files (*.smi);;All files (*.*)'
    elif self.file_tdtButton.isOn():
      fmt = 'tdt'
      fmtStr = 'TDT Files (*.tdt);;All files (*.*)'
    elif self.file_txtButton.isOn():
      fmt = 'txt'
      fmtStr = 'TDT Files (*.txt);;All files (*.*)'
    fileN = str(QFileDialog.getOpenFileName(self._dir,fmtStr))
    if fileN:
      self.info_fileButton.setChecked(1)
      self._dir,fileN = os.path.split(fileN)
      self._fileName = fileN
      self._fmt = fmt
      self.file_nameBox.setText(fileN)
    self.refreshContents()

  def refreshContents(self):
    l = self.getLen()
    self.info_countBox.setMaxValue(l)
    self.info_countBox.setValue(l)
    if l:
      self.info_countBox.setEnabled(1)
      self.updateTable()
    else:
      self.info_countBox.setEnabled(0)
      self.supplierTab.setTabEnabled(self.contentsPage,0)

      
  def accept(self):
    """ accept changes made to the reaction """
    _Form.accept(self)



if __name__ == '__main__':
  from qtGui import Gui

  app,widg = Gui.Launcher(LoaderWidget)
  app.exec_loop()
  widg.destroy(1)
