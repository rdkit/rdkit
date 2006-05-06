#! /usr/bin/env python
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
#  SVN Info:
#   $Author: glandrum $
#   $Date: 2006-01-12 08:20:06 -0800 (Thu, 12 Jan 2006) $
#   $LastChangedRevision: 4899 $
#
import os
env = os.environ
if not env.has_key('RDBASE'):
  env['RDBASE']=os.getcwd()

from qt import *
from RDOutline import RDOutlineImpl
from ToDoTree import ToDoTree
from ToDoListView import ToDoListView
from TextViewImpl import TextView
from HelpData import helpData

__version = '0.9.0'
__date = '$Date: 2006-01-12 08:20:06 -0800 (Thu, 12 Jan 2006) $'
helpData = helpData%locals()

        
class RDToDoImpl(RDOutlineImpl):
  def __init__(self,fileName=None,mainWin=0):
    RDOutlineImpl.__init__(self,fileName,populate=0,mainWin=mainWin)
    self._initListWin()
    self._initFonts()
    self._dispDoneId = self._viewMenu.insertItem(self.trUtf8('Show &Done?'),
                                                 self.toggleShowDone)
    self._viewMenu.setItemChecked(self._dispDoneId,1)
    self._dispListId = self._viewMenu.insertItem(self.trUtf8('Show &Win?'),
                                                 self.toggleListView)
    self._viewMenu.setItemChecked(self._dispListId,not self.listView.isHidden())

    self._smallIconsId = self._optionMenu.insertItem(self.trUtf8('&Small Icons?'),
                                                 self.toggleSmallIcons)
    self._viewMenu.setItemChecked(self._smallIconsId,0)

    self.editMenu.insertItem(self.trUtf8('&Duplicate Item'),
                             self.editDuplicateItem,
                             Qt.SHIFT+Qt.CTRL+Qt.Key_D,-1,7)


    self._populate(os.path.join(self._dirName,self._fileName))
    self.treeView.setNeedsSave(0);
    self.initIcon()

  def _initTreeWin(self):
    self.splitter = QSplitter(self)
    self.treeView = ToDoTree(self.splitter,'ToDo Tree',owner=self)
    self.setCentralWidget(self.splitter)
    tree = self.treeView
    tree.show()
    tree.setRootIsDecorated(1)

  def _initListWin(self):
    self.listView = ToDoListView(self.treeView,self.splitter,
                                 'ToDo List')
    self.listView.setFocusPolicy(QWidget.ClickFocus)
    self.listView.addColumn('Item')
    self.listView.addColumn('Pri')
    self.listView.addColumn('Due')
    self.listView.hide()

  def _initFonts(self,baseFont=None):
    fnt = RDOutlineImpl._initFonts(self,baseFont=baseFont)
    try:
      self.listView.setFont(fnt)
    except AttributeError:
      pass
    
  def initIcon(self):
    from icons.icons import pngs
    pm = QPixmap()
    pm.loadFromData(pngs[7])
    self.setIcon(pm)

  def toggleSmallIcons(self):
    isChecked = not self._optionMenu.isItemChecked(self._smallIconsId)
    if not isChecked:
      print 'nope'
    else:
      print 'yep'
    self._viewMenu.setItemChecked(self._smallIconsId,isChecked)
    
  def toggleListView(self):
    isChecked = not self._viewMenu.isItemChecked(self._dispListId)
    if not isChecked:
      self.listView.hide()
    else:
      self.listView.show()
    self._viewMenu.setItemChecked(self._dispListId,isChecked)

  def toggleShowDone(self):
    isChecked = not self._viewMenu.isItemChecked(self._dispDoneId)
    if not isChecked:
      self.treeView._invisibleThings['percentDone']=[100,'100']
    else:
      try:
        del self.treeView._invisibleThings['percentDone']
      except:
        pass
    self.treeView.updateVisibility()
    self._viewMenu.setItemChecked(self._dispDoneId,isChecked)
      
  def helpContents(self):
    self.helpView = TextView(helpData)
    self.helpView.show()

  def editDuplicateItem(self):
    item = self.treeView.selectedItem()
    if item and not item.attribs().get('duplicate',0):
      from ToDoTreeItem import ToDoTreeDupeItem
      parent = item.parent() or self.treeView
      heads = self.treeView.headings()
      things = ['']*len(heads)
      for i,head in enumerate(heads):
        things[i] = item.text(i)
      dupe = ToDoTreeDupeItem(parent,*things)
      dupe.attribs()['text']=item.attribs()['text']
      dupe.attribs()['duplicate']=item.GUID()
      dupe.setGUID()
      dupe.moveItem(item)
      dupe.setDragEnabled(True)
      dupe._updateStatus(item.getPercentDone())
      self.treeView.setSelected(dupe,True)
      self.treeView.setNeedsSave(True)
      self.treeView.refresh(start=parent)
    
if __name__=='__main__':
  import sys,getopt,os
  fileN = os.path.join(os.getcwd(),'ToDo.opml')
  if len(sys.argv) > 1:
    args,extras = getopt.getopt(sys.argv[1:],'')
    if len(extras):
      fileN = extras[0]
  app = QApplication(sys.argv)
  win = RDToDoImpl(fileName=fileN,mainWin=0)
  app.setMainWidget(win)
  win.show()
  app.exec_loop()


