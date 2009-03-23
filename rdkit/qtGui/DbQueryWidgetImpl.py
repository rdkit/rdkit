#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for DbQueryWidgets

"""    
from rdkit import RDConfig
from qt import *
from rdkit.qtGui.forms.DbQueryWidget import DbQueryWidget as _Form
from rdkit.qtGui.forms.ListSelectorDialog import ListSelectorDialog
from rdkit.qtGui import DbWidget
import os,os.path
from rdkit.Dbase import DbConnection


def insertQueryWidget(parent,*args,**kwargs):
  """ constructs a _DbQueryWidget_ and inserts it into a parent

   This uses _DbWidget.insertDbWidget_ to do its work
   
  """
  return DbWidget.insertDbWidget(parent,DbQueryWidget,*args,**kwargs)

class DbQueryWidget(DbWidget.DbWidgetMixin,_Form):
  """  Widget for specifying parameters to set up a database query

   The widget includes fields for:

     - the database

     - username

     - password

     - name of the table

     - SQL what

     - SQL where

     - SQL join

   **Form:**  _DbQueryWidget_

  """
  def __init__(self,parent=None,initDir='',clickCallback=None):
    _Form.__init__(self,parent)
    DbWidget.DbWidgetMixin.__init__(self,parent,initDir,clickCallback)

  def sqlWhat(self):
    """ return the _sqlWhat_ value """
    return str(self.sqlWhatBox.text())
  def setWhat(self,val):
    """ return the _sqlWhat_ value """
    self.sqlWhatBox.setText(str(val))
    
  def sqlJoin(self):
    """ return the _sqlJoin_ value """
    return str(self.sqlJoinBox.text())
  def setJoin(self,val):
    self.sqlJoinBox.setText(str(val))
    
  def sqlWhere(self):
    """ return the _sqlWhere_ value """
    return str(self.sqlWhereBox.text())
  def setWhere(self,val):
    self.sqlWhereBox.setText(str(val))


  def pickColumnNames(self,startingSelection='',allowMultiple=1):
    """ launches a dialog allowing the user to specify which columns to load,
        returns the list of column names.

    """
    dlg = ListSelectorDialog()
    lst = dlg.listBox
    names = self.getColumnNames()
    for name in names:
      lst.insertItem(name)
      if startingSelection.find(name)>-1:
        lst.setSelected(lst.count()-1,1)
    if not allowMultiple:
      lst.setSelectionMode(QListBox.Single)
    ok = dlg.exec_loop()
    res = []
    if ok:
      selNames = []
      for i in range(lst.count()):
        if lst.isSelected(i):
          selNames.append(str(lst.text(i)))
      if len(selNames):    
        res = selNames
    return res
    
  #
  # Slots and Signals
  #
  def dbChooseClick(self):
    """ callback for clicking on the "Choose" button
    (for filling the "what" field of the SQL query)

    """
    origWhat = self.sqlWhat()
    self.setWhat('*')
    selNames = self.pickColumnNames(startingSelection=origWhat)
    if selNames:
      self.setWhat(','.join(selNames))
    else:
      self.setWhat(origWhat)

if __name__ == '__main__':
  # build the app and widget
  import os
  from rdkit.qtGui import Gui
  app,widg = Gui.Launcher(DbQueryWidget,None)

  app.exec_loop()
  widg.destroy(1)

