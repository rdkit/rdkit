# $Id: PubmedSearchImpl.py 4684 2005-05-25 21:50:46Z glandrum $
#
#  Copyright (C) 2003-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for Pubmed searching

"""    
import RDConfig
from qt import *
from qttable import *
from qtGui.GuiLib.forms.PubmedSearchWidget import PubmedSearchWidget as _Form
from qtGui import qtUtils
from qtGui.forms import CountDialog
from Dbase.Pubmed import Searches,QueryParams,PubmedUtils
from qtGui import GuiTable
from qtGui.GuiLib.PubmedRecordImpl import PubmedRecord
import os

class PubmedSummaryTable(GuiTable.GuiTable):
  cols=['Keep','PubMedId','Title','Authors','Source','PubYear','Volume',
         'Pages','RecordStatus']
  def __init__(self,*args,**kwargs):
    if kwargs.has_key('cols'):
      self.cols = kwargs['cols']
      del kwargs['cols']
    else:
      self.cols = self.cols[:]
    GuiTable.GuiTable.__init__(self,*args,**kwargs)
    self.setNumCols(len(self.cols))
    header = self.horizontalHeader()
    for i in range(len(self.cols)):
      header.setLabel(i,self.cols[i])
  def addRow(self,summary):
    nRows = self.numRows()
    row = nRows
    self.insertRows(nRows,1)
    nCols = len(self.cols)
    for col in range(nCols):
      if col != 0:
        colName = self.cols[col]
        val = getattr(summary,colName)
        self.setText(row,col,val)
      else:
        self.setItem(row,col,QCheckTableItem(self,''))
  def cleanup(self):
    for col in range(self.numCols()):
      self.adjustColumn(col)

  def findNamedCol(self,label):
    hdr = self.horizontalHeader()
    res = -1
    for idx in range(hdr.count()):
      heading = str(hdr.label(idx))
      if heading==label:
        res = idx
        break
    return res  
    
  def getChecked(self):
    res = []
    checkCol = self.findNamedCol('Keep')
    idCol =  self.findNamedCol('PubMedId')
    if checkCol >=0 and idCol >=0:
      for row in range(self.numRows()):
        itm = self.item(row,checkCol)
        if itm and itm.isChecked():
          id = str(self.text(row,idCol))
          if id:
            res.append(id)
    else:
      if idCol < 0:
        qtUtils.warning('PubMedId Column has been deleted')
    return res
      
    
      
class PubmedSearchWidget(_Form):
  """ DOC
  

  """
  def __init__(self,*args,**kwargs):
    _Form.__init__(self,*args,**kwargs)
    self._dir = '.'
    self._initQuery()
    self._resPages = []
    self._records = {}
    self.connect(self.tabWidget,SIGNAL('currentChanged(QWidget *)'),self.tabPageChanged)
    
  def _insertQueryRow(self,row,enableIt=0):
    grid = self._queryGridLayout
    if row >= grid.numRows():
      grid.expand(row+1,3)
    page = self.queryPage
    # FIX: update the size setting on the combo boxes
    if row != 0:
      mod = QComboBox(page)
      mod.setEditable(0)
      mod.insertItem("AND")
      mod.insertItem("OR")
      #mod.setMaximumSize(QSize(50,field.maximumSize().height()))
      grid.addWidget(mod,row,0)
      mod.show()
      mod.setEnabled(enableIt)
      self.connect(mod,SIGNAL("activated(const QString &)"),
                   self.searchFieldsUpdated)
    else:
      mod = None
    field = QComboBox(page)
    field.setEditable(0)
    field.show()
    field.setEnabled(enableIt)
    for fieldName in QueryParams.searchableFieldsOrder:
      field.insertItem(fieldName)
    self.connect(field,SIGNAL("activated(const QString &)"),
                 self.searchFieldsUpdated)
    grid.addWidget(field,row,1)
    box = QLineEdit(page)
    box.show()
    box.setEnabled(enableIt)
    self.connect(box,SIGNAL("textChanged(const QString &)"),
                 self.searchFieldsUpdated)
    grid.addWidget(box,row,2)
    tt = QueryParams.searchableFields[QueryParams.searchableFieldsOrder[0]][1]
    QToolTip.add(box,tt)
    QToolTip.add(field,tt)
    self._queryFields.append((mod,field,box))
    
  def _initQuery(self,nBoxes=4):
    self.queryLine.setText('penzotti je[au] AND grootenhuis pd[au]')
    self._queryFields = []
    page = self.queryPage
    grid = QGridLayout(nBoxes,3)
    self._queryGridLayout = grid
    for i in range(nBoxes):
      self._insertQueryRow(i,enableIt=(i==0))

    page.layout().insertItem(0,QSpacerItem(0,10,QSizePolicy.Minimum,QSizePolicy.Expanding))
    page.layout().insertLayout(0,grid)
      
  def searchFieldsUpdated(self,arg):
    query = ""
    for i in range(len(self._queryFields)):
      mod,field,box = self._queryFields[i]
      txt = str(box.text()).strip()
      fieldName = str(field.text(field.currentItem()))
      if txt:
        if mod:
          query += ' %s '%(str(mod.text(mod.currentItem())))
        modifier = QueryParams.searchableFields[fieldName][0]
        query += '%s[%s]'%(txt,modifier)

      tt = QueryParams.searchableFields[fieldName][1]
      QToolTip.remove(box)
      QToolTip.add(box,tt)
      QToolTip.remove(field)
      QToolTip.add(field,tt)

      # enable/disable other boxes, as appropriate
      if not txt:
        for j in range(i+1,len(self._queryFields)):
          m2,f2,b2 = self._queryFields[j]
          m2.setEnabled(0)
          f2.setEnabled(0)
          b2.setEnabled(0)
        break
      else:
        if i+1<len(self._queryFields):
          m2,f2,b2 = self._queryFields[i+1]
          m2.setEnabled(1)
          f2.setEnabled(1)
          b2.setEnabled(1)
        else:
          self._insertQueryRow(i+1,enableIt=1)
    self.queryLine.setText(query)
          
      
  def doSearch(self,query=None):
    if not query:
      query = QueryParams.details()
      query['term'] = str(self.queryLine.text())
    if not query:
      qtUtils.information('No query provided\n Search aborted.')
      return []

    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    try:
      nHits = Searches.GetNumHits(query)
    except:
      qtUtils.logger.info('problems encountered in retrieving hits',exc_info=True)
      nHits = 0
    QApplication.restoreOverrideCursor()

    if not nHits:
      qtUtils.information("Search returned no results.")
    else:
      ans = qtUtils.infoPrompt("Search returned %d hits.\nRetrieve Summaries?"%(nHits))
      if ans == QMessageBox.Yes:
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        query['rettype']='uilist'
        query['retmax']=nHits
        try:
          ids = Searches.GetSearchIds(query)
        except:
          qtUtils.error('Problems retrieving results ids.')
          ids = []
        if len(ids) != nHits:
          qtUtils.error('Retrieved %d results, %d were expected.'%(len(ids),nHits))
          ids = []
        QApplication.restoreOverrideCursor()
        if not ids:
          return
        pageTitle = 'Search %d results'%(len(self._resPages)+1)
        widg = QWidget(self.tabWidget,pageTitle)
        layout = QVBoxLayout(widg)
        tbl = GuiTable.insertTable(widg,PubmedSummaryTable)

        # display the query on the results page:
        box = QLineEdit(widg)
        box.setReadOnly(1)
        box.setText(self.queryLine.text())
        widg.layout().addWidget(box)

        self.tabWidget.addTab(widg,pageTitle)
        widg._resTable = tbl
        self._resPages.append(widg)
        self.fillSummaryTable(ids,tbl)
        self.tabWidget.setCurrentPage(len(self._resPages))
  def fillSummaryTable(self,ids,tbl,query=None,summaries=None):
    if not tbl or not ids:
      return
    tbl.setUpdatesEnabled(0)
    tbl.setNumRows(0)
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    if summaries is None:
      try:
        # FIX: break this into pieces and use a progress dialog
        summaries = Searches.GetSummaries(ids)
      except:
        qtUtils.logger.info('problems encountered in retrieving summaries',exc_info=True)
        summaries = None
    if summaries:
      row = 0
      for summary in summaries:
        tbl.addRow(summary)
    QApplication.restoreOverrideCursor()
    tbl.cleanup()
    tbl.setUpdatesEnabled(1)
    # the checkboxes don't seem to do the right thing if they start out selected,
    # this hackery solves that:
    tbl.setCurrentCell(0,1)
    tbl.setCurrentCell(0,0)
    
  def getRecords(self,ids=None,showThem=0):
    if not ids:
      widg = self.tabWidget.currentPage()
      if hasattr(widg,'_resTable'):
        ids = widg._resTable.getChecked()
    for id in ids[:]:
      if self._records.has_key(id):
        self._records[id].show()
        ids.remove(id)
    if ids:
      recs = Searches.GetRecords(ids)
    else:
      recs = []
    for rec in recs:
      widg = PubmedRecord(record=rec)
      widg.setCaption('Pubmed Record: %s'%rec.PubMedId)
      if showThem:
        widg.show()
      self._records[rec.PubMedId] = widg
    return recs
  
  def exportRecords(self,ids=None,fileName=None):
    if not fileName:
      fileName = str(QFileDialog.getSaveFileName(self._dir,
                                                 'Text files (*.txt);;All files (*.*)'))
      if fileName:
        self._dir = os.path.dirname(fileName)
    if not fileName:
      return
    try:
      outF = open(fileName,'w+')
    except:
      qtUtils.error('Could not open file %s for writing.\n'%fileName)
      return
      
    # ok, we're set... retrieve the records:
    records = self.getRecords(ids=ids,showThem=0)
    # ... and write them out:
    res = PubmedUtils.RecordsToPubmedText(records)
    outF.write(res)
    outF.close()
    return res
    
  
  def getRelatedRecords(self,ids=None,query=None):
    if not ids:
      widg = self.tabWidget.currentPage()
      if hasattr(widg,'_resTable'):
        ids = widg._resTable.getChecked()

    if ids:
      QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
      try:
        links = Searches.GetLinks(ids)
      except:
        qtUtils.logger.info('problems encountered in retrieving links',exc_info=True)
        links = []
      QApplication.restoreOverrideCursor()
      nRes = len(links)
      dlg = CountDialog.CountDialog()
      dlg.countBox.setMaxValue(nRes)
      default = min(10,nRes)
      dlg.countBox.setValue(default)
      dlg.labelText.setText("Search returned %d related records.\nSelect number to be retrieved."%(nRes))
      dlg.setCaption("Related Record Count")
      dlg.resize(dlg.minimumSizeHint())
      res = dlg.exec_loop()
      if res == QDialog.Accepted:
        nToPull = dlg.countBox.value()
      else:
        nToPull = 0
      if nToPull:
        pulls = links[:nToPull]
        idsToPull = [x[0] for x in pulls]
        scores = [x[1] for x in pulls]
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        try:
          summs = Searches.GetSummaries(idsToPull)
        except:
          qtUtils.logger.info('problems encountered in retrieving summaries',exc_info=True)
          summs = []
        QApplication.restoreOverrideCursor()
        pageTitle = 'Related results %d'%(len(self._resPages)+1)
        widg = QWidget(self.tabWidget,pageTitle)
        layout = QVBoxLayout(widg)
        tbl = GuiTable.insertTable(widg,PubmedSummaryTable)
        box = QLineEdit(widg)
        box.setReadOnly(1)
        box.setText(','.join(ids))
        widg.layout().addWidget(box)
        self.tabWidget.addTab(widg,pageTitle)
        widg._resTable = tbl
        self._resPages.append(widg)
        self.fillSummaryTable(idsToPull,tbl)
        tbl.insertColumns(1,1)
        tbl.cols.insert(1,'Score')
        tbl.horizontalHeader().setLabel(1,'Score')
        for i in range(len(summs)):
          tbl.setText(i,1,str(scores[i]))
        
        self.tabWidget.setCurrentPage(len(self._resPages))
        
        




  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def searchButtonClicked(self):
    self.doSearch()

  def recordsButtonClicked(self):
    self.getRecords()

  def relatedButtonClicked(self):
    self.getRelatedRecords()

  def exportButtonClicked(self):
    self.exportRecords()

  def tabPageChanged(self,widg):
    if hasattr(widg,'_resTable'):
      self.recordsButton.setEnabled(1)
      self.relatedButton.setEnabled(1)
      self.exportButton.setEnabled(1)
    else:
      self.recordsButton.setEnabled(0)
      self.relatedButton.setEnabled(0)
      self.exportButton.setEnabled(0)

if __name__ == '__main__':
  from qtGui import Gui

  app,widg = Gui.Launcher(PubmedSearchWidget,None,'PubmedSearch')
  app.exec_loop()

