#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" defines the DescTable class

"""    
from rdkit import RDConfig
from qt import *
from qttable import *

class DescTable(QTable):
  """ a table class for selecting descriptors

  """
  checkCol=1
  def __init__(self,*args,**kwargs):
    QTable.__init__(self,*args)
    QObject.connect(self,SIGNAL('selectionChanged()'),self.selectionChanged)
    QObject.connect(self,SIGNAL('currentChanged(int,int)'),self.currentChanged)
    self.setSelectionMode(self.Single)
    self._lastSel=None

  def currentChanged(self,row,col):
    pass
    
  def selectionChanged(self):
    """ handles selections

      currently this ignores stuff that's not in the check column

    """
    # FIX: this does not handle selections particularly properly
    #  I remain less than 100% convinced that this is completely
    #  correct, but it's a lot better than it was.
    sel = self.selection(self.currentSelection())
    anchor = sel.anchorRow(),sel.anchorCol()
    leftCol = sel.leftCol()
    rightCol = sel.rightCol()
    if leftCol==rightCol==self.checkCol:
      state = self.item(anchor[0],anchor[1]).isChecked()
      # they just selected the entire col, this is an inversion
      top,bot = sel.topRow(),sel.bottomRow()
      if self.isColumnSelected(self.checkCol,1):
        state = not state
        self.item(anchor[0],anchor[1]).setChecked(state)
        
      if top != bot:
        for row in range(top,bot+1):
          if row != anchor[0]:
            self.item(row,anchor[1]).setChecked(state)
    self._lastSel = anchor

    
  def getSelected(self):
    """ returns the indices of checked items

    """
    res = []
    for i in xrange(self.numRows()):
      if self.item(i,self.checkCol).isChecked():
        res.append(i)
    return res


  
