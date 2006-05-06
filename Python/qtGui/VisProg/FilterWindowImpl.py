#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for FilterWindows

  NOTE: this has been put together with the greatest of speed, so
   there ain't much in the way of docs.

"""    
from qt import *
import os,os.path,sys
from qtGui.VisProg.forms.FilterWindow import FilterWindow as _Form
import VPUtils

class FilterWindow(_Form):
  def __init__(self,filterObj=None,initDir=''):
    """

    """
    _Form.__init__(self)

    self.guts = filterObj
    self._dir=initDir
    self._fileName = ''
    self._methodName = ''
    self._code = ''
    # fix: too close to the bone
    self.info_numInputsBox.setValue(self.guts._nInputs)
    self.info_numOutputsBox.setValue(self.guts._nOutputs)
    
  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def fileClick(self):
    if sys.platform == 'win32':
      sysLibs="Library files (*.dll *.pyd)"
    else:
      sysLibs="Library files (*.so *.pyo *.pyd)"
      
    types = "Python files (*.py *.pyc);;%s;;All files (*.*)"%sysLibs
    fileN = str(QFileDialog.getOpenFileName(self._dir,types))
    if fileN:
      self._dir,fileN = os.path.split(fileN)
      self._dbName = fileN
      self.info_nameBox.setText(fileN)
      self.info_numInputsBox.setEnabled(1)
      self.info_numOutputsBox.setEnabled(1)
    else:
      pass
    self.refreshContents()
    

  def refreshContents(self):
    pass

      
  def accept(self):
    """ accept changes made to the reaction """
    txt = self.nodeTextBox.text()
    self.guts.setLabelText(txt)
    # FIX: uh, bit too specific here
    self.guts._canvas.update()
    _Form.accept(self)
