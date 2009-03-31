#
#  Copyright (C) 2003  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Allows similarity searching with molecules

"""
from rdkit import RDConfig
from qt import *
from rdkit.qtGui.GuiLib.SimilaritySearchWidgetImpl import SimilaritySearchWidget
from rdkit.qtGui.qtUtils import logoImageData
import os

REQUIRED_MIXINS = ['MolBrowser']
MODULES_ALTERED = ['rdkit.qtGui.GuiBase']
METHODS_DEFINED = {
  '__LaunchSimilarityWin':'rdkit.qtGui.GuiBase.GuiBase.msLaunchSimilarityWin',
  }


class MolSimilarityWin(QMainWindow):
  def __init__(self,*args,**kwargs):
    QMainWindow.__init__(self,*args,**kwargs)
    self.setCaption('SimilaritySearch')
    self.searchWidget = SimilaritySearchWidget(self)
    self.connect(self.searchWidget.tabWidget,SIGNAL('currentChanged(QWidget*)'),
                 self.pageChanged)
    self.setCentralWidget(self.searchWidget)
    self.updateGeometry()
    self._roPages=(self.searchWidget.queryPage,self.searchWidget.paramsPage,
                self.searchWidget.databasePage)
    self._initMenubar()
    image0 = QPixmap(logoImageData)
    self.setIcon(image0)

  def _initMenubar(self):
    mb = self.menuBar()
    self.fileMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&File'),self.fileMenu)

    self.editMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&Edit'),self.editMenu)
    self._pageDelIdx = self.editMenu.insertItem(self.trUtf8('&Delete Active Page'),
                                           self.delPage)


    self.dbMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&Database'),self.dbMenu)
    self._saveTableIdx = self.dbMenu.insertItem(self.trUtf8('&Save Active Page'),
                                                  self.pageToDb)


    mb.insertSeparator()
    self.helpMenu = QPopupMenu(self)
    mb.insertItem(self.trUtf8('&Help'),self.helpMenu)
    self.editMenu.setItemEnabled(self._pageDelIdx,0)


  def pageChanged(self,page):
    if page in self._roPages:
      self.editMenu.setItemEnabled(self._pageDelIdx,0)
      self.dbMenu.setItemEnabled(self._saveTableIdx,0)
    else:
      self.editMenu.setItemEnabled(self._pageDelIdx,1)
      self.dbMenu.setItemEnabled(self._saveTableIdx,1)

    
  def delPage(self,idx=0,page=None):
    tabWidget = self.searchWidget.tabWidget
    confirm = 1
    if page is None:
      page = tabWidget.currentPage()
      if page in self._roPages:
        page = None
      else:
        idx = tabWidget.currentPageIndex()
        label = tabWidget.label(idx)
        confirm = QMessageBox.warning(self,
                                      self.trUtf8("Delete Page?"),
                                      self.trUtf8("Are you sure you want to delete the page %s?"%(label)),
                                      1,2)
    if page and confirm:
      tabWidget.removePage(page)
      tabWidget.setCurrentPage(0)
      
  def pageToDb(self,idx=0,page=None):
    tabWidget = self.searchWidget.tabWidget
    confirm = 1
    if page is None:
      page = tabWidget.currentPage()
      if page in self._roPages:
        page = None
    if page and hasattr(page,'_tbl'):
      page._tbl.tableToDb(skipCols=[],replaceMolName=1)
    else:
      print 'no _tbl'


def __LaunchSimilarityWin(self):
  if self._msSimilarityWin is None:
    win = MolSimilarityWin()
    self._msSimilarityWin = win
  else:
    win = self._msSimilarityWin
  win.show()
  
def LocalInit(self):
  self._msSimilarityWin = None
  self._mbMolMenu.insertItem(self.trUtf8('&Similarity Search'),self.msLaunchSimilarityWin)


                                
                                
