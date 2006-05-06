#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" contains a mixin class for the various Database widgets

"""    
import RDConfig
from qt import *
import os,os.path
from Dbase import DbConnection
from qtGui.qtUtils import logger

def insertDbWidget(parent,klass,*args,**kwargs):
  """ constructs a widget and inserts it into a parent widget

    The widget is inserted into a QVBoxLayout along with a
    vertical spacer.
  

  """
  layout = parent.layout()
  if not layout:
    layout = QVBoxLayout(parent,2,2,'dbwidgetlayout')
  obj = klass(parent,*args,**kwargs)
  layout.insertWidget(0,obj)
  layout.insertStretch(1)
  return obj

class DbWidgetMixin:
  """ mixin class providing common non-GUI functionality
    for DbWidgets

  """
  def __init__(self,parent=None,initDir='',clickCallback=None):
    self._dbDir = initDir
    self.parent = parent
    self._allowTableEdits = 0
    self._fileClickCallbacks = []
    if clickCallback:
      self.insertFileClickCallback(clickCallback)

  def insertFileClickCallback(self,fn):
    self._fileClickCallbacks.append(fn)

  def allowTableEdits(self,allow=1):
    """ setter for allowTableEdits attribute """
    self._allowTableEdits = allow
    
  def dbName(self,fullyQualified=1):
    """ returns the name of the selected database

    **Arguments**

      - fullyQualified: (optional) if this flag is nonzero, the name
        will be returned with an attached directory name

    **Returns**

      a string
      
    """
    base = str(self.nameBox.text())
    if fullyQualified:
      res = os.path.join(self._dbDir,base)
    else:
      res = base
    return res
  def tableName(self):
    """ returns the name of the selected table from the widget

    """
    if self.tableCombo.isEnabled():
      table = str(self.tableCombo.currentText())
    else:
      table = ''
    return table
  def user(self):
    """ returns the user name from the widget """
    return str(self.userField.text())    
  def password(self):
    """ returns the password from the widget  """
    return str(self.passwordField.text())    
  def sqlWhat(self):
    """ return the _sqlWhat_ value """
    return '*'
  def sqlJoin(self):
    """ return the _sqlJoin_ value """
    return ''
  def sqlWhere(self):
    """ return the _sqlWhere_ value """
    return ''

  
  def getConn(self):
    """ builds a _DbConnection.DbConnect_ instance and returns it

    """
    user = self.user()
    passwd = self.password()
    table = self.tableName()
    conn = DbConnection.DbConnect(self.dbName(),
                                  tableName = table,
                                  user=user,password=passwd)
    return conn


  def getData(self,**kwargs):
    """ returns a data set from the db based on the current settings

    """
    if not kwargs.has_key('fields'):
      kwargs['fields'] = self.sqlWhat()
    if not kwargs.has_key('where'):
      kwargs['where'] = self.sqlWhere()
    if not kwargs.has_key('join'):
      kwargs['join'] = self.sqlJoin()
    return self.getConn().GetData(**kwargs)


  def getDataCount(self,**kwargs):
    """ returns the numbea data set from the db based on the current settings

    """
    if not kwargs.has_key('fields'):
      kwargs['fields'] = self.sqlWhat()
    if not kwargs.has_key('where'):
      kwargs['where'] = self.sqlWhere()
    if not kwargs.has_key('join'):
      kwargs['join'] = self.sqlJoin()
    return self.getConn().GetDataCount(**kwargs)

  def getColumnNames(self):
    """ returns the column names for a query based on the current settings

    """
    return self.getConn().GetColumnNames(what=self.sqlWhat(),
                                         where=self.sqlWhere(),
                                         join=self.sqlJoin())

  def getColumnNamesAndTypes(self):
    """ returns the column names and types for a query based on the current settings

    """
    return self.getConn().GetColumnNamesAndTypes(what=self.sqlWhat(),
                                                 where=self.sqlWhere(),
                                                 join=self.sqlJoin())

  def getDbFilename(self):
    if RDConfig.usePgSQL:
      return None
    fileN = str(QFileDialog.getOpenFileName(self._dbDir,"Interbase files (*.gdb);;All files (*.*)"))
    if fileN:
      self._dbDir,fileN = os.path.split(fileN)
      self._dbName = fileN
    return fileN
    
  def getDbName(self):
    if not RDConfig.usePgSQL:
      return None
    from Dbase import DbInfo
    names = QStringList.fromStrList(DbInfo.GetDbNames())
    res,ok = QInputDialog.getItem('Select Database','Select Database',names,0,False)

    if ok:
      return res
    else:
      return ''

  def enableQueries(self):
    self.tableCombo.setEnabled(1)
    self.tableCombo.setEditable(self._allowTableEdits)

    for name in ['sqlWhatBox','sqlChooseButton','sqlWhereBox','sqlJoinBox']:
      if hasattr(self,name):
        getattr(self,name).setEnabled(1)


  #    
  # Slots
  #
  def dbFileClick(self,fileN=None):
    """ callback for clicking on the file button

    """
    if fileN is None:
      if RDConfig.usePgSQL:
        fileN=self.getDbName()
      else:
        fileN=self.getDbFilename()
       
       
       
    if fileN:
      self.nameBox.setText(fileN)
      conn = self.getConn()
      tblNames = conn.GetTableNames(includeViews=1)
      while self.tableCombo.count() != 0:
        self.tableCombo.removeItem(self.tableCombo.count()-1,)
      for name in tblNames:
        name = name.strip()
        self.tableCombo.insertItem(name)

      self.enableQueries()
    elif len(str(self.nameBox.text())) == 0:
      self.tableCombo.setEnabled(0)

    for callback in self._fileClickCallbacks:
      try:
        callback()
      except:
        logger.debug('callback failed',exc_info=True)

        

