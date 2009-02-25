# $Id$
#
#  Copyright (C) 2005-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for CompositeWindows

"""    
import RDConfig
from qt import *
from qtcanvas import *
from qtGui.GuiLib.forms.CompositeWindow import CompositeWindow as _Form
from qtGui import GuiTable
from qtGui.GuiTextViewer import GuiTextViewer
from qtGui.GuiLib import ModelTable,CompositeUtils
from qtGui.GuiLib.ScreenCompositeDialogImpl import ScreenCompositeDialog
from qtGui import qtUtils,DbQueryDialog
from ML import AnalyzeComposite
from ML.DecTree import TreeVis
import os

class CompositeWindow(_Form):
  """ a class used for interacting with composites

    The base widget is defined in forms.CompositeWindow

  """
  def __init__(self,*args,**kwargs):
    _Form.__init__(self,*args,**kwargs)
    self._count = 1
    self._pages = []
    self._dir = '.'
    self._analysisDepth=4
    #self._queryDialog = DbQueryDialog.DbQueryDialog(self)
    self._screenWin = ScreenCompositeDialog()
    self._screenWin.hide()
    self._screenWin.setCompositeWindow(self)
    
    self._screenResultsWin = GuiTextViewer()
    self._screenResultsWin.hide()

    self._updateMenuBar()
    
  def _updateMenuBar(self):
    """ adds some extra stuff to the menu bar

    """
    self._optionsMenu = QPopupMenu(self)
    self.menuBar().insertItem(self.trUtf8('&Options'),self._optionsMenu,-1,1)
    #self._scaleLeavesCheck = QCheckBox(self.trUtf8('&Scale Leaves'),self._optionsMenu)
    #self._optionsMenu.insertItem(self._scaleLeavesCheck)
    #self._scaleLeavesCheck.setChecked(0)
    self._optionsMenu.setEnabled(0)

    self._actionsMenu = QPopupMenu(self)
    self.menuBar().insertItem(self.trUtf8('&Actions'),self._actionsMenu,-1,2)
    self._actionsMenu.insertItem(self.trUtf8('&Analyze'),self.addAnalysisPage)
    self._actionsMenu.insertItem(self.trUtf8('A&nalyze All'),self.addAllAnalysisPage)
    #self._screenMenu = QPopupMenu(self._actionsMenu)
    #self._screenMenu.insertItem(self.trUtf8('From &Db'),self.screenCompositeFromDetails)
    #self._actionsMenu.insertItem(self.trUtf8('&Screen'),self._screenMenu)
    self._actionsMenu.insertItem(self.trUtf8('&Screen'),self._screenWin.show)

    self._actionsMenu.insertSeparator()
    self._actionsMenu.insertItem(self.trUtf8('&Delete Active Page'),self.delPage)

  def addAnalysisPage(self,idx=0,page=None,composites=None):
    """ inserts an analysis page into our notebook

    **Arguments**

      - idx: (optional) not used

      - page: (optional) the page containing the composite to be
        screened.  If this is not provided, the active page will be
        used

      - composites (optional) composites to be analyzed. If not provided,
        the composite associated with the active notebook page will
        be used. 
        
    **Notes**

      - the heavy lifting for the analysis is done by
        _AnalyzeComposite.ProcessIt_
        
    """
    if page is None:
      page = self.tabWidget.currentPage()
    if not page:
      return
    if not composites:
      if hasattr(page,'_composite'):
        composites = [page._composite]
      else:
        qtUtils.warning('no composite')
        return
      newPageTitle = '%s - Analysis'%(page.name())
    else:
      newPageTitle = 'Multiple Composite Analysis'
      
    res = AnalyzeComposite.ProcessIt(composites,nToConsider=self._analysisDepth,verbose=-1)

    widg = QWidget(self.tabWidget,newPageTitle)
    tbl = GuiTable.insertTable(widg,GuiTable.GuiTable)
    self.tabWidget.addTab(widg,newPageTitle)
    nCols = len(res[0])
    tbl.setNumCols(nCols)
    tbl.setNumRows(len(res))
    tbl.setReadOnly(1)

    hdr = tbl.horizontalHeader()
    hdr.setLabel(0,'Descriptor')
    for i in range(1,nCols):
      hdr.setLabel(i,'Level %d'%(i))
    hdr.setLabel(i,'Sum')

    hdr = tbl.verticalHeader()
    for i in range(len(res)):
      for j in range(nCols):
        tbl.setText(i,j,str(res[i][j]))

    for i in range(nCols):
      tbl.adjustColumn(i)
    self.tabWidget.showPage(widg)

  def addAllAnalysisPage(self,idx=0):
    """ inserts an analysis page into our notebook for all loaded composites

    **Arguments**

      - idx: (optional) not used


    """
    composites = []
    nPages = self.tabWidget.count()
    for i in range(nPages):
      page = self.tabWidget.page(i)
      if hasattr(page,'_composite'):
        composites.append(page._composite)
    if composites:
      self.addAnalysisPage(composites=composites)
    else:
      qtUtils.warning('No Composites found')
                    
  def scaleDrawing(self):
    """ getter for our scaleDrawing property

    """
    return self._scaleLeavesCheck.isChecked()

  def saveComposite(self,fileName):
    """ saves our composite by pickling it into a file

    **Arguments**

      - fileName: the filename to use

    """
    try:
      compos = self.tabWidget.currentPage()._composite
      if not compos:
        return
      else:
        compos.Pickle(fileName)
    except AttributeError:
      pass

  def activeComposite(self):
    """ #DOC

    """
    try:
      compos = self.tabWidget.currentPage()._composite
      if not compos:
        return
      else:
        return compos
    except AttributeError:
      return

  def allComposites(self):
    """ #DOC

    """
    res = []
    for i in range(self.tabWidget.count()):
      page = self.tabWidget.page(i)
      try:
        compos = page._composite
      except AttributeError:
        compos = None
      if compos:
        res.append(compos)
    return res

  def addComposite(self,composite,drawTarget=None,name=None):
    """ adds a composite page to our notebook

    **Arguments**

      - composite: the composite to add

      - drawTarget: (optional) the window to which the composite should draw

      - name: (optional) the name to use for the page
      

    """
    composNum = self._count
    self._count += 1
    if not name:
      name = 'Composite %d'%(composNum)

    widg = QWidget(self.tabWidget,name)
    widg._composite = composite
    tbl = ModelTable.insertModelTable(widg)
    self.tabWidget.addTab(widg,name)
    self.tabWidget.showPage(widg)
    widg._tbl = tbl
    tbl.load(composite)
    if drawTarget:
      tbl.setTarget(drawTarget)
    self._screenWin.setComposite(composite,name)
    
  def delPage(self,idx=0,page=None):
    """ removes a page from the notebook

    **Arguments**

      - idx: (optional) not used

      - page: (optional) the page to delete, if not provided the
        current page is removed.

    **Notes**

      - the first page in the notebook (if there is one) is activated
        after page removal 

    """
    if page is None:
      page = self.tabWidget.currentPage()
    if not page:
      return
    # FIX: this is, how you say, a hack
    try:
      page._tbl.hide()
    except AttributeError:
      return
    self.tabWidget.removePage(page)
    if self.tabWidget.count():
      self.tabWidget.setCurrentPage(0)

  def screenCompositeFromDetails(self,idx=0,composite=None,details=None):
    """ screens the active composite

    **Arguments**

      - idx: (optional) not used

      - composite: (optional) the composite to use.  If not provided,
        the active window's composite will be used

      - details: (optional) a _CompositeRun.CompositeRun_ structure
        containing screening details.  If not provided,
        _ML.BuildCompsite.SetDefaults()_ will be used

     **Notes**

      - our _queryDialog_ is used to prompt for the name of the
        database containing the screening data.

      - the actual screening is done using
        _CompositeUtils.ScreenCompositeFromDetails()_

      - HTML screening results will be shown in _self._screenWin_. If
        this already exists, it will be cleared.

    """
    if composite is None:
      widg = self.tabWidget.currentPage()
      try:
        composite = widg._tbl.composite()
      except AttributeError:
        pass
      if composite is None:
        qtUtils.error('No composite loaded')
        return
    if details is None:
      details = BuildComposite.SetDefaults()
    dlg = self._queryDialog
    res = dlg.exec_loop()
    if res == QDialog.Accepted:
      details.dbName = dlg.widget.dbName()
      details.tableName = dlg.widget.tableName()
      details.dbUser = dlg.widget.user()
      details.dbPassword = dlg.widget.password()
      details.dbWhat = dlg.widget.sqlWhat()
      details.dbWhere = dlg.widget.sqlWhere()
      details.dbJoin = dlg.widget.sqlJoin()
    self._screenWin = GuiTextViewer()
    CompositeUtils.ScreenCompositeFromDetails(details,composite,self._screenWin)

    # reset the visualization options
    #  FIX: this maybe violates encapsulation a bit too much
    tbl = widg._tbl
    for i in range(tbl.numRows()):
      itm = tbl.item(i,0)
      TreeVis.ResetTree(itm.model())
    
  #
  # Slots and signals and stuff
  #
  def fileSaveAs(self):
    """ callback """
    # FIX: this directory nonsense sucks
    fileN = str(QFileDialog.getSaveFileName(self._dir,'Pickle files (*.pkl);;All files (*.*)'))
    if fileN:
      self._dir = os.sep.join(os.path.split(fileN)[0:-1])
      self.saveComposite(fileN)

  def fileClose(self):
    """ callback """
    self.hide()
        
    
