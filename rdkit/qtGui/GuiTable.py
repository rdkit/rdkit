# $Id$
#
#  Copyright (C) 2002-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" defines classes to make QTables more useful for us

"""    
from rdkit import RDConfig
from qt import *
from qttable import *
from rdkit.qtGui.DbConnDialog import DbConnDialog
from rdkit.qtGui.qtUtils import logger
try:
  from cStringIO import StringIO
except ImportError:
  from StringIO import StringIO
from rdkit.Dbase import DbUtils
import types


class GuiTableHeaderEventFilter(QObject):
  """ event filter for clicks in _GuiTable_ headers


  This ignores everything except right and left button clicks which
  emit the signals _headerRightClick_ and _headerLeftClick_ respectively.

  Because a _GuiTable_ handles right clicks on headers itself (by
  launching a context menu), these are filtered out.  Left clicks are
  allowed to pass.
  

  """
  HORIZONTAL=0
  VERTICAL=1
  def __init__(self,table,whichHdr,*args,**kwargs):
    QObject.__init__(self,*args)
    self._table = table
    self._whichHdr = whichHdr
    
  def eventFilter(self,obj,evt):
    """

      if an event is handled, return 1 immediately on finish,
      otherwise let it flow through, we'll return 0 at the bottom

    """
    if evt.type()==QEvent.MouseButtonPress:
      b = evt.button()
      if b==Qt.RightButton:
        self._table.emit(PYSIGNAL("headerRightClick"),(evt,self._whichHdr))
        return 1
      elif b==Qt.LeftButton:
        self._table.emit(PYSIGNAL("headerLeftClick"),(evt,self._whichHdr))

        
    return 0

def insertTable(parent,klass,*args,**kwargs):
  """ inserts a table into a _QVBoxLayout_ in a parent

  **Arguments**

    - parent: the widget into which we insert

    - klass: the class (**not** an instance) which is to be inserted

    - if the keyword argument 'addInPos' is provided, the table
      will be inserted into a particular location

  """
  # DOC: finish documenting these
  if kwargs.has_key('addInPos'):
    addInPos = kwargs['addInPos']
    del kwargs['addInPos']
  else:
    addInPos = -1
  if kwargs.has_key('clearExisting'):
    clearExisting = kwargs['clearExisting']
    del kwargs['clearExisting']
  else:
    clearExisting = 1

  tbl = klass(parent,*args,**kwargs)
  layout = parent.layout()
  if not layout:
    layout = QVBoxLayout(parent,1,1,'tablelayout')
  else:
    if not layout.isEmpty():
      if clearExisting and layout.children() and len(layout.children()):
        layout.deleteAllItems()
  if addInPos < 0:
    layout.addWidget(tbl)
  else:
    layout.insertWidget(addInPos,tbl)

  parent.adjustSize()
  parent._tbl = tbl
  return tbl


_acceptableTypes = (types.IntType,types.FloatType,
                    types.StringType,types.UnicodeType,
                    types.LongType,types.ComplexType)
def loadTable(tbl,data,colNames=None,adjust=1,acceptableTypes=None):
  """ slaps the data passed in into the table
  #DOC

  **Arguments**

  **Notes**

    - the contents of the table will be wiped out

  """
  if acceptableTypes is None:
    acceptableTypes = _acceptableTypes
  nRows = len(data)
  if colNames:
    nCols = len(colNames)
  else:
    nCols = len(data[0])
    colNames = map(str,range(nCols))
  tbl.setNumRows(nRows)
  tbl.setNumCols(nCols)
  
  # update the header
  hdr = tbl.horizontalHeader()
  for col in range(nCols):
    hdr.setLabel(col,colNames[col])

  # load the points
  for row in range(nRows):
    pt = data[row]
    for col in range(nCols):
      if type(pt[col]) in acceptableTypes:
        txt = str(pt[col])
      else:
        txt = "<Undisplayable>"
      tbl.setText(row,col,txt)
      
  if adjust:
    for col in range(nCols):
      tbl.adjustColumn(col)

class GuiTable(QTable):
  """ an augmented QTable class

  """
  def __init__(self,*args,**kwargs):
    QTable.__init__(self,*args)

    self._sortOrders = None
    self._allowColInserts = 0
    self._allowRowInserts = 0
    self._allowColRename = 1
    self._allowRowRename = 0
    self._allowColDelete = 1
    self._allowRowDelete = 1
    
    self._hHdr = self.horizontalHeader()
    self._vHdr = self.verticalHeader()
    QObject.connect(self,SIGNAL("currentChanged(int,int)"),self.currentChanged)
    QObject.connect(self,SIGNAL("selectionChanged()"),self.selectionChanged)
    self._hFilt = GuiTableHeaderEventFilter(self,GuiTableHeaderEventFilter.HORIZONTAL)
    self.horizontalHeader().installEventFilter(self._hFilt)
    self._vFilt = GuiTableHeaderEventFilter(self,GuiTableHeaderEventFilter.VERTICAL)
    self.verticalHeader().installEventFilter(self._vFilt)
    QObject.connect(self,PYSIGNAL("headerRightClick"),self.headerRightClick)
    QObject.connect(self,PYSIGNAL("headerLeftClick"),self.headerLeftClick)
    self.setColumnMovingEnabled(1)
    self.setRowMovingEnabled(1)

    QObject.connect(self,SIGNAL("contextMenuRequested(int,int,const QPoint &)"),
                    self.contextMenuRequested)


  def _limitCols(self,skipCols):
    """ #DOC

    """
    cols = range(self.numCols())
    for col in skipCols:
      try:
        cols.remove(col)
      except:
        pass
    return cols
    
  def setAllowColInsertion(self,val):
    self._allowColInserts = val
  def allowColInsertion(self):
    return self._allowColInserts
  def setAllowRowInsertion(self,val):
    self._allowRowInserts = val
  def allowRowInsertion(self):
    return self._allowRowInserts

  def setAllowColRename(self,val):
    self._allowColRename = val
  def allowColRename(self):
    return self._allowColRename
  def setAllowRowRename(self,val):
    self._allowRowRename = val
  def allowRowRename(self):
    return self._allowRowRename
    
  def setAllowColDelete(self,val):
    self._allowColDelete = val
  def allowColDelete(self):
    return self._allowColDelete
  def setAllowRowDelete(self,val):
    self._allowRowDelete = val
  def allowRowDelete(self):
    return self._allowRowDelete
    
  def contentsAsList(self,skipCols=[]):
    """ #DOC

     #EFF: this stuff ain't so super fast because of the type
      conversion code

    """
    cols = self._limitCols(skipCols)
    nCols = len(cols)
    nRows = self.numRows()
    res = [None]*(nRows)
    for row in xrange(nRows):
      tmp = ['']*nCols
      for i in range(nCols):
        col = cols[i]
        itm = self.item(row,col)
        if itm:
          try:
            v = float(str(itm.text()))
          except:
            v = str(itm.text())
          tmp[i] = v
        else:
          tmp[i] = ''
      res[row]=tmp
    return res
    
  def contentsToDelimText(self,delim=',',skipCols=[]):
    """ converts the contents of the table to a delimited text string
       suitable for saving or putting on the clipboard
       
    **Arguments**

      - delim: (optional) the delimiter string to use

      - skipCols: (optional) a sequence containing columns which
        should not be included in the output.

    **Returns**

       a string

    """
    cols = self._limitCols(skipCols)
    nCols = len(cols)
    hdr = self.horizontalHeader()
    tmp = ['']*nCols
    for i in range(nCols):
      tmp[i] = str(hdr.label(cols[i]))
    res = [delim.join(tmp)]
    d = self.contentsAsList(skipCols)

    #EFF: this stuff ain't so super speedy
    for pt in d:
      res.append(delim.join(map(str,pt)))
    return '\n'.join(res)

  def tableToDb(self,**kwargs):
    """ dumps the table into a database
    
    **Notes**

      - the target db details are collected using a _DbConnDialog_

      - Actual insertion of the data into the db is done by converting
        our contents into delimited text (using
        _self.contentsToDelimText()_) and then inserting using
        _DbUtils.TextFileToDatabase()_, which is the worst sort of
        cheating. :-)

    """
    dlg = DbConnDialog(self)
    dlg.dbWidget().allowTableEdits(1)
    res = dlg.exec_loop()
    if res == QDialog.Accepted:
      widget = dlg.dbWidget()
      conn = widget.getConn()
      tableNames = [x.strip() for x in conn.GetTableNames()]
      dlg = None
      name = widget.tableName().strip()
      confirm = 1
      if name in tableNames:
        confirm = QMessageBox.warning(self,
                                      self.trUtf8("Save Table"),
                                      self.trUtf8("Table %s already exists, should it be replaced?"%(name)),
                                      1,2)
      if confirm==1:
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        sio = StringIO()
        kwargs['delim']='\t'
        sio.write(self.contentsToDelimText(**kwargs))
        sio.seek(0)
        DbUtils.TextFileToDatabase(widget.dbName(),widget.tableName(),sio,
                                     delim='\t',user=widget.user(),password=widget.password())
        QApplication.restoreOverrideCursor()
                                     

  def numColumns(self):
    """ returns the number of columns in the table

    **Notes***

     - it's really dumb that we need to have this, but the qt guys,
       in their infinite wisdom, decided to name just about everything else
       using "Column" instead of "Col" but to do the reverse here.  Yay
       for consistency!

    """
    return self.numCols()

  def _sortHelper(self,i1,i2):
    """ Internal use only

      used to sort column contents
      
    """
    k1 = i1.key()
    k2 = i2.key()
    try:
      v1 = float(str(k1))
      v2 = float(str(k2))
    except:
      return k1.localeAwareCompare(k2)
    else:
      return cmp(v1,v2)
  
  def sortColumn(self,col,ascending=0,wholeRows=1):
    """ causes the entire table to be sorted based on the contents of
    an individual column

    **Arguments**

      - col: index of the column on which to sort

      - ascending: (optional) used to set the sort order

      - wholeRows: (optional) toggles motion of entire rows instead of
        just the contents of the column

    **Notes**

      - at the moment whole row sorts are always done (the _wholeRows_
        argument is ignored).

    """
    if ascending is None:
      ascending = self._sortOrders[col]
      self._sortOrders[col] = not ascending
      
    # FIX: this only supports whole row sorts
    items = []
    nRows = self.numRows()
    rows = xrange(nRows)
    for row in rows:
      itm = self.item(row,col)
      if itm:
        items.append(itm)
    nFilled = len(items)
    items.sort(self._sortHelper)
    if ascending: items.reverse()
    self.setUpdatesEnabled(0)
    for i in rows:
      if i < nFilled:
        if items[i].row()!=i:
          self.swapRows(items[i].row(),i)
    self.setUpdatesEnabled(1)
    self.horizontalHeader().setSortIndicator(col,ascending)
    self.redraw()

  def redraw(self):
    self.repaintContents( self.contentsX(),self.contentsY(),
                          self.visibleWidth(),self.visibleHeight(), 0)
    
  def moveColToEnd(self,col):
    """ moves the specified column to the end
    """
    self.setUpdatesEnabled(0)
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    try:
      destPos = self.numCols()
      for i in range(col,destPos-1):
        self.swapColumns(i,i+1,1)
    except:
      logger.debug('problems encountered moving colum',exc_info=True)
    QApplication.restoreOverrideCursor()
    self.setUpdatesEnabled(1)
    self.redraw()
  def insertColumn(self,col):
    col += 1
    self.insertColumns(col,1)
    self.horizontalHeader().setLabel(col,'Col %d'%(col+1))
    self.redraw()
  def insertRow(self,row):
    row += 1
    self.insertRows(row,1)
    self.verticalHeader().setLabel(row,'Row %d'%(row+1))
    self.redraw()
  def renameColumn(self,col,name=None):
    """ rename a table column
    
    """
    hdr = self.horizontalHeader()
    lab = str(hdr.label(col))
    if name is None:
      txt,ok = QInputDialog.getText(self.trUtf8('Enter Column Label'),
                                    self.trUtf8('Enter Column Label'),
                                    QLineEdit.Normal,
                                    str(lab))
      if ok:
        name = str(txt)
      else:
        name = lab
    hdr.setLabel(col,str(txt))
  def renameRow(self,row,name=None):
    """ rename a table row
    
    """
    hdr = self.verticalHeader()
    lab = str(hdr.label(row))
    if name is None:
      txt,ok = QInputDialog.getText(self.trUtf8('Enter Row Label'),
                                    self.trUtf8('Enter Row Label'),
                                    QLineEdit.Normal,
                                    str(lab))
      if ok:
        name = str(txt)
      else:
        name = lab
    hdr.setLabel(row,str(name))
  def delColumn(self,cols=None):
    """ removes a column (or columns) from the table

    **Arguments**

     - cols: either a sequence or integer with the column(s) to be
       deleted

    """
    if cols is None:
      cols = []
      for i in range(self.numColumns()):
        if self.isColumnSelected(i,1):
          cols.append(i)
    if type(cols) == type(1):
      cols = [cols]
    cols.sort()
    cols.reverse()
    self.setUpdatesEnabled(0)
    for col in cols:
      self.removeColumn(col)
    self.setUpdatesEnabled(1)
    self.redraw()

  def delRow(self,rows=None):
    """ removes a row or rows

    **Arguments**

     - rows: either a sequence or integer with the row(s) to be
       deleted

    """
    if rows is None:
      rows = []
      for i in range(self.numRows()):
        if self.isRowSelected(i,1):
          rows.append(i)
    if type(rows) == type(1):
      rows = [rows]
    rows.sort()
    rows.reverse()
    self.setUpdatesEnabled(0)
    for row in rows:
      self.removeRow(row)
    self.setUpdatesEnabled(1)
    self.redraw()


  def horizHeaderRightClick(self,evt):
    if self._sortOrders is None or len(self._sortOrders)!=self.numCols():
      self._sortOrders = [0]*self.numCols()
    # FIX:  I'm probably reading the docs wrong, but it doesn't seem
    #  to me like this bit of trickery should be required.  Still, as of
    #  qt version 3.0.5 it was
    pos = evt.pos()
    col = self.columnAt(pos.x()+self.contentsX())
    if col >= 0:
      self._clickCol=col
      menu = QPopupMenu(self)
      fn1 = lambda a,x=self,y=col:x.sortColumn(y,None,1)
      menu.insertItem(self.trUtf8('&Sort'),fn1)
      menu.insertSeparator()
      hdr = self.horizontalHeader()
      if not hasattr(self.horizontalHeader(),'isMovingEnabled') or \
         self.horizontalHeader().isMovingEnabled():
        fn2 = lambda a,x=self,y=col:x.moveColToEnd(y)
        menu.insertItem(self.trUtf8('&Move To End'),fn2)
        
      if self.allowColRename():
        fn3 = lambda a,x=self,y=col:x.renameColumn(y)
        menu.insertItem(self.trUtf8('&Rename'),fn3)
      if self.allowColInsertion():
        fn4 = lambda a,x=self,y=col:x.insertColumn(y)
        menu.insertItem(self.trUtf8('&Insert'),fn4)
      if self.allowColDelete():
        fn5 = lambda a,x=self,y=col:x.delColumn([y])
        menu.insertItem(self.trUtf8('&Delete This Column'),fn5)
        fn6 = lambda a,x=self,y=col:x.delColumn()
        menu.insertItem(self.trUtf8('Delete &Selected Columns'),fn6)

      menu.exec_loop(evt.globalPos())
    
  def vertHeaderRightClick(self,evt):
    # FIX:  I'm probably reading the docs wrong, but it doesn't seem
    #  to me like this bit of trickery should be required.  Still, as of
    #  qt version 3.0.5 it was
    pos = evt.pos()
    row = self.rowAt(pos.y()+self.contentsY())
    if row >= 0:
      self._clickRow=row
      menu = QPopupMenu(self)
      if self.allowRowRename():
        fn3 = lambda a,x=self,y=row:x.renameRow(y)
        menu.insertItem(self.trUtf8('&Rename'),fn3)
      if self.allowColInsertion():
        fn4 = lambda a,x=self,y=row:x.insertRow(y)
        menu.insertItem(self.trUtf8('&Insert'),fn4)
      if self.allowRowDelete():
        fn5 = lambda a,x=self,y=row:x.delRow([y])
        menu.insertItem(self.trUtf8('&Delete This Row'),fn5)
        fn6 = lambda a,x=self,y=row:x.delRow()
        menu.insertItem(self.trUtf8('Delete &Selected Rows'),fn6)

      menu.exec_loop(evt.globalPos())
  #
  #  Signals 'n stuff
  #
  def headerRightClick(self,evt,whichOne):
    """ callback for right clicks on column headers

    """
    if whichOne == GuiTableHeaderEventFilter.HORIZONTAL:
      self.horizHeaderRightClick(evt)
    elif whichOne == GuiTableHeaderEventFilter.VERTICAL:
      self.vertHeaderRightClick(evt)
    else:
      logger.error('invalid header selected',exc_info=True)
      
  def headerLeftClick(self,evt,whichOne):
    """ callback for left clicks on column headers

      Currently does nothing.

    """
    pass


  def currentChanged(self,row,col):
    pass

  def selectionChanged(self):
    pass

  def contextMenuRequested(self,row,col,pos):
    """ to be overridden in child classes """
    pass



    
