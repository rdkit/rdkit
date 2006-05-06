#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for SinkWindows

  NOTE: this has been put together with the greatest of speed, so
   there ain't much in the way of docs.

"""    
from qt import *
import re
import os,os.path
from qtGui.VisProg.forms.SinkWindow import SinkWindow as _Form
import VPUtils

from Dbase import DbInfo,DbConnection


class SinkWindow(_Form):
  def __init__(self,obj=None,initDir=''):
    """

    """
    _Form.__init__(self)

    self.guts = obj
    self._dir=initDir
    self._fileName = ''
    self._dbName = ''

  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def dbFileClick(self):
    fileN = str(QFileDialog.getOpenFileName(self._dir,"Interbase files (*.gdb);;All files (*.*)"))
    if fileN:
      self.info_dbButton.setChecked(1)
      self._dir,fileN = os.path.split(fileN)
      self._dbName = fileN
      self.db_nameBox.setText(fileN)
      self.db_tableCombo.setEnabled(1)
      user = str(self.db_userField.text())
      passwd = str(self.db_passwordField.text())
      tblNames = DbInfo.GetTableNames(os.path.join(self._dir,self._dbName),
                                       user=user,
                                       password=passwd)
      while self.db_tableCombo.count() != 0:
        self.db_tableCombo.removeItem(self.db_tableCombo.count()-1,)
      for name in tblNames:
        self.db_tableCombo.insertItem(name)
    elif len(str(self.db_nameBox.text())) == 0:
      self.db_tableCombo.setEnabled(0)
    self.refreshContents()
    

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

      
  def accept(self):
    """ accept changes made to the reaction """
    txt = self.nodeTextBox.text()
    self.guts.setLabelText(txt)
    # FIX: uh, bit too specific here
    self.guts._canvas.update()
    _Form.accept(self)
