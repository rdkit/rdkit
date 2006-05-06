# $Id: ModelTable.py 5064 2006-03-09 02:03:00Z glandrum $
#
#  Copyright (C) 2002-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" defines classes required to put models in QTables

"""    
import RDConfig
from qt import *
from qttable import *
from qtGui import GuiTable,qtUtils
from ML.DecTree import TreeVis,Tree

class ModelTableItem(QTableItem):
  """ class for containing models in ModelTables

  """
  def __init__(self,*args,**kwargs):
    QTableItem.__init__(self,*args)
    self._model = None

  def setModel(self,model):
    """ setter for model attribute """
    self._model = model
  def model(self):
    """ getter for model attribute """
    return self._model

  def draw(self,where,nRes=2):
    """ causes this model (Tree) to be drawn in a given window

    **Arguments**

      - where: the window in which to draw

      - nRes: (optional) the number of possible result codes for the model

    **Notes**

      - the destination window should have the following methods:

         - _canvas()_ : which returns a Piddle/Sping canvas on which
           we can draw

         - _size()_ : returns a 2-tuple with the size of the canvas

         - _view()_ : returns an object (used for book-keeping, at the
           moment this should be a _PiddleCanvasView_)

         - _show()_ : return value ignored

      - drawing is done using _TreeVis.DrawTree()_

    """
    if isinstance(self._model,Tree.TreeNode):
      where.setTree(self._model)
      where.draw(nRes=nRes)
    else:
      qtUtils.logger.info('Cannot draw "%s" yet.'%(repr(self._model)))

def insertModelTable(where,*args,**kwargs):
  """ inserts a _ModelTable_ into a given window

    The heavy lifting is done using _GuiTable.insertTable()_
    
  """
  
  return GuiTable.insertTable(where,ModelTable,*args,**kwargs)
  
class ModelTable(GuiTable.GuiTable):
  """ a class able to intelligently hold models in a _GuiTable_

    The table has three columns:

      1) Model: contains _ModelTableItems_

      2) Cont: contains ints

      3) Errors: contains floats

  """
  def __init__(self,*args,**kwargs):
    GuiTable.GuiTable.__init__(self,*args,**kwargs)
    self.setNumCols(3)
    hdr = self.horizontalHeader()
    hdr.setLabel(0,'Model')
    hdr.setLabel(1,'Count')
    hdr.setLabel(2,'Errors')
    self._target = None
    self._composite = None
    self.adjustSize()
    
  def destroy(self,win=1,subwins=1):
    """ wipes out this table after disconnecting from all signals

    """
    # FIX: required?
    self._target = None
    self.disconnect(self,SIGNAL("currentChanged(int,int)"),self.currentChanged)
    GuiTable.GuiTable.destroy(self,win,subwins)
    
  def target(self):
    """ getter for target attribute """
    # because we'll maybe need to be able to resize it, our "canvas"
    #  is actually a qtGui.PiddleWindow
    return self._target
  def setTarget(self,target):
    """ setter for target attribute """
    # because we'll maybe need to be able to resize it, our "canvas"
    #  is actually a qtGui.PiddleWindow
    self._target = target
      
  def setModel(self,row,dataTuple):
    """ adds a model to a row

    **Arguments**

      - row: the row to be set

      - dataTuple a 3-tuple:

         1) model

         2) count

         3) err
      

    """
    model,count,err = dataTuple
    try:
      model._varNames = self.composite()._varNames
    except:
      pass
    model._gridName = 'Model %d'%(row)
    itm = ModelTableItem(self,QTableItem.Never,model._gridName)
    itm.setModel(model)
    self.setItem(row,0,itm)
    self.setText(row,1,str(count))
    self.setText(row,2,str(err))
    
  def composite(self):
    """ getter for our composite attribute, the setter is _load()_ """
    return self._composite

  def _setComposite(self,composite):
    """ INTERNAL USE ONLY

    """
    self._composite = composite
  
  def load(self,composite):
    """ loads a composite """
    nModels = len(composite)
    self.setNumRows(nModels)
    self._setComposite(composite)

    for i in range(nModels):
      self.setModel(i,composite.GetDataTuple(i))

  #
  # Signals/slots
  #
  def currentChanged(self,row,col):
    """ callback for selection changes in the widget

    """
    if col == 0:
      if self.target():
        p = self.topLevelWidget()
        self.item(row,col).draw(self.target())
    
