#
#  Copyright (C) 2003  Rational Discovery LLC
#    All Rights Reserved
#
""" Implementation of the BrowserInput widget

"""    
import RDConfig
from qt import *
from qtGui import qtUtils
from qtGui.HierarchyBrowser.forms.BrowserInput import BrowserInput as _Form
from Chem import FragmentCatalog
import os,cPickle

def insertBrowserInputWidget(parent,*args,**kwargs):
  """ constructs a BrowserInputWidget and inserts it into a parent

  """
  layout = parent.layout()
  if not layout:
    layout = QVBoxLayout(parent,2,2,'browserinputwidgetlayout')
  obj = BrowserInput(parent,*args,**kwargs)
  layout.insertWidget(0,obj)
  layout.insertStretch(1)
  return obj


class BrowserInput(_Form):
  """
  **Form:**  _BrowserInput_

  """
  def __init__(self,parent=None,initDir='.',updateCallback=None):
    _Form.__init__(self,parent)
    self._catDir = initDir
    self._gainsDir = initDir

    self._updateCallbacks = []
    if updateCallback:
      self.insertUpdateCallback(updateCallback)

    self._catalog = None
    self._bits = None

  def insertUpdateCallback(self,fn):
    self._updateCallbacks.append(fn)

  def catalog(self):
    return self._catalog
  def setCatalog(self,val):
    self._catalog = val
  def bits(self):
    return self._bits
  def setBits(self,val):
    self._bits = val
    
  def updateFields(self,readCatalog=0,readGains=0):
    if readCatalog:
      catName = str(self.catalogEdit.text())
      if catName:
        fileN = os.path.join(self._catDir,catName)
        try:
          inF = open(fileN,'rb')
        except IOError:
          qtUtils.error('Could not open file %s for reading.'%(fileN))
        else:
          QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
          # FIX: some kind of informational dialog should really be displayed
          #  here.
          try:
            cat = cPickle.load(inF)
          except:
            import traceback
            traceback.print_exc()
            qtUtils.error('Problems encountered loading catalog.\n')
          else:
            if not hasattr(cat,'GetCatalogParams'):
              qtUtils.error('Object loaded does not look like a catalog.\n')
            else:
              self.setCatalog(cat)
          QApplication.restoreOverrideCursor()
    if readGains:
      gainsName = str(self.gainsEdit.text())
      if gainsName:
        fileN = os.path.join(self._gainsDir,gainsName)
        try:
          inF = open(fileN,'rb')
        except IOError:
          qtUtils.error('Could not open file %s for reading.'%(fileN))
        else:
          QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
          try:
            bits = FragmentCatalog.ProcessGainsFile(fileN)
          except:
            import traceback
            traceback.print_exc()
            qtUtils.error('Problems encountered loading gains information.\n')
          else:
            self.setBits(bits)
          QApplication.restoreOverrideCursor()
    for callback in self._updateCallbacks:
      try:
        callback()
      except:
        import traceback
        traceback.print_exc()
    if self.catalog() and self.bits():
      # FIX: this is too specific to FragCatalogs:
      self.numBitsSpin.setEnabled(1)
      self.numBitsSpin.setMaxValue(len(self.bits()))
      params = self.catalog().GetCatalogParams()
      self.minLevelSpin.setEnabled(1)
      self.minLevelSpin.setMinValue(params.GetLowerFragLength())
      self.minLevelSpin.setMaxValue(params.GetUpperFragLength())
      self.minLevelSpin.setValue(params.GetLowerFragLength())
      self.maxLevelSpin.setEnabled(1)
      self.maxLevelSpin.setMinValue(params.GetLowerFragLength())
      self.maxLevelSpin.setMaxValue(params.GetUpperFragLength())
      self.maxLevelSpin.setValue(params.GetUpperFragLength())
    else:
      self.numBitsSpin.setEnabled(0)
      self.minLevelSpin.setEnabled(0)
      self.maxLevelSpin.setEnabled(0)
      

  #
  # Slots and Signals
  #
  def catalogBrowseClicked(self,fileN=None):
    if fileN is None:
      fileN = str(QFileDialog.getOpenFileName(self._catDir,"Pickle Files (*.pkl);;All files (*.*)"))
    if fileN:
      self._catDir,fileN = os.path.split(fileN)
      self.catalogEdit.setText(fileN)
      self.updateFields(readCatalog=1)
  def catalogReturnPressed(self):
    self.updateFields(readCatalog=1)


  def gainsBrowseClicked(self,fileN=None):
    if fileN is None:
      fileN = str(QFileDialog.getOpenFileName(self._gainsDir,"Text Files (*.txt *.csv);;All files (*.*)"))
    if fileN:
      self._gainsDir,fileN = os.path.split(fileN)
      self.gainsEdit.setText(fileN)
      self.updateFields(readGains=1)

  def gainsReturnPressed(self):
    self.updateFields(readGains=1)
  
if __name__ == '__main__':
  # build the app and widget
  import os
  from qtGui import Gui
  app,widg = Gui.Launcher(BrowserInput,None)
  app.exec_loop()
  widg.destroy(1)

