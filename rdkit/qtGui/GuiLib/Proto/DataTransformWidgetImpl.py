#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for DataTransformWidgets

"""    
import RDConfig
from qt import *
from qttable import *
from forms.DataTransformWidget import DataTransformWidget as _Form
from ML.Data import Transforms

class DataTransformWidget(_Form):
  """ 

  
  """
  def __init__(self,parent=None,initCount=3,dataSet=None):
    _Form.__init__(self)
    self._dataSet = dataSet

    self._initAvailTransforms()
    self._initGrid(initCount)

  def _initAvailTransforms(self):
    """ INTERNAL USE ONLY

    """
    # FIX: this needs to be massively refactored, the custom stuff shouldn't
    #  go here
    self._availTransforms = list(Transforms.GetAvailTransforms())
    self._availTransforms.append(('Custom',None,'uses the formula defined to the right'))

  def _initGrid(self,count):
    """ INTERNAL USE ONLY

    """
    tbl = self.choiceTable
    tbl.setNumRows(0)
    for i in range(count):
      self.addRow()
    #tbl.horizontalHeader().setStretchEnabled(1,3)
    tbl.setColumnStretchable(3,1)
    
  def availTransforms(self):
    """ getter for _availTransforms_ attribute

    **Returns**

      a tuple of tuples

    """
    return tuple(self._availTransforms)
  
    

  def addRow(self,name=None):
    tbl = self.choiceTable
    nRows = tbl.numRows()
    # create the new row
    tbl.insertRows(nRows,1)

    if name is None:
      name = 'Operation-%d'%(nRows+1)
      
    tbl.setText(nRows,0,name)

    names = ['Columns','Rows']
    if self._dataSet:
      names += list(self._dataSet.GetVarNames())
    stringL = QStringList()
    for name in names:
      stringL.append(name)
    itm = QComboTableItem(tbl,stringL)
    tbl.setItem(nRows,1,itm)

    names = [x[0] for x in self.availTransforms()]
    stringL = QStringList()
    for name in names:
      stringL.append(name)
    itm = QComboTableItem(tbl,stringL)
    tbl.setItem(nRows,2,itm)
    #tbl.setColumnWidth(3,200)
    


  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def addClick(self):
    self.addRow()

  def helpClick(self):
    pass

  def loadClick(self):
    pass

  def saveClick(self):
    pass
    
    
  def accept(self):
    """ close dialog with the ok button

    """
    _Form.accept(self)

if __name__ == '__main__':
  from qtGui import Gui

  app,widg = Gui.Launcher(DataTransformWidget)
  app.exec_loop()
  widg.destroy(1)
