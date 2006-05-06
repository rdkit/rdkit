#! /usr/bin/env python
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
#  SVN Info:
#   $Author: glandrum $
#   $Date$
#   $LastChangedRevision$
#
from qt import *
from forms.mainform import mainForm as _Form
from TextViewImpl import TextView
from HelpData import helpData
from xml.dom import minidom
import sys,os.path
from OutlineTree import OutlineTree

__version = '0.1.2'

defaultFontFamily="Arial"
defaultFontSize=10

class RDOutlineImpl(_Form):
  def __init__(self,fileName=None,populate=1,mainWin=0):
    _Form.__init__(self)
    self._initTreeWin()
    font = self._initFonts()
    if mainWin:
      qApp.setFont(font)
    self._dirName = ''
    self._fileName = ''
    if populate:
      self._populate(fileName)
    if fileName is not None:
      dir,path=os.path.split(fileName)
      self._dirName = dir
      self._fileName = path
    self._viewMenu = QPopupMenu(self)
    self._viewMenu.insertItem(self.trUtf8('&Font'),self.viewSetFont)
    self._viewMenu.insertItem(self.trUtf8('&Expand All'),self.viewExpandAll,
                              Qt.ALT+Qt.CTRL+Qt.Key_E)
    self._viewMenu.insertItem(self.trUtf8('&Collapse All'),self.viewCollapseAll,
                              Qt.SHIFT+Qt.CTRL+Qt.ALT+Qt.Key_E)
    self.menuBar().insertItem(self.trUtf8('&View'),self._viewMenu,-1,2)

    self._optionMenu = QPopupMenu(self)
    self.menuBar().insertItem(self.trUtf8('&Options'),self._optionMenu,-1,3)

    self.editUndoAction.setEnabled(0)
    self.editRedoAction.setEnabled(0)
    self.editFindAction.setEnabled(0)
    self.helpIndexAction.setEnabled(0)
    self.helpAboutAction.setEnabled(0)
    self.filePrintAction.setEnabled(0)

    self.fileMenu.insertItem(self.trUtf8('&Export'),self.fileExport,
                             0,-1,4)
    self.editMenu.insertItem(self.trUtf8('Copy as &Text'),self.editCopyText,
                             Qt.SHIFT+Qt.CTRL+Qt.Key_C,-1,5)

    self._childWidgets = []
    self.treeView.setNeedsSave(0)
    self.customFont=None
    
  def _initFonts(self,baseFont=None):
    if baseFont is None:
      baseFont=QFont()
      baseFont.setStyleHint(QFont.SansSerif)
      baseFont.setFamily(defaultFontFamily)
      baseFont.setPointSize(defaultFontSize)
    self.treeView.setFont(baseFont)
    return baseFont


  def _initTreeWin(self):
    self.treeView = OutlineTree(self,'Outline View')
    self.setCentralWidget(self)
    tree = self.treeView
    tree.show()
    tree.setRootIsDecorated(1)
    
  def _populate(self,fileN):
    if fileN is None:
      return
    try:
      inD = open(fileN,'r').read()
    except:
      print 'can not open demo file'
      self.blankWindow()
    else:
      self.treeView.addOPML(inD,new=1)

  def blankWindow(self):
    self.treeView._initColSettings()
    self.treeView._colInit()
    doc = minidom.Document()
    node = doc.createElement('outline')
    node.setAttribute('text','New Outline')
    newItm = self.treeView.addDomNode(node,select=1)

  def confirmAction(self):
    if self.treeView.needsSave():
      response = QMessageBox.information(self,sys.argv[0],
                                         "You haven't saved your work\nSave Changes?",
                                         "&Save","&Discard","&Cancel",0,2)
    else:
      response = 1
    return response  

  #
  #  signals and slots and handlers and stuff
  #
  def viewSetFont(self):
    if self.customFont is not None:
      fnt = self.customFont
    else:
      fnt = self.treeView.font()
    newFont,ok = QFontDialog.getFont(fnt,self)
    if ok:
      self.customFont = newFont
      self._initFonts(self.customFont)


  def fileExport(self):
    print 'export'
  def fileExit(self):
    response = self.confirmAction()
    if response == 2:
      return
    elif response == 0:
      self.fileSave()
    qApp.quit()

  def fileSaveAs(self):
    fileN = str(QFileDialog.getSaveFileName(self._dirName,
                                            "OPML files (*.opml);;All files (*.*)"))
    if fileN:
      self._dirName,self._fileName = os.path.split(fileN)
      self.fileSave()

      
  def fileSave(self):
    if self._fileName == '':
      self.fileSaveAs()
    else:
      outF = open(os.path.join(self._dirName,self._fileName),'w+')
      txt = self.treeView.toxml()
      txt = txt.replace('><','>\n<')
      outF.write(txt)
      outF.close()
      self.treeView.setNeedsSave(0)


  def _fileNew(self):
    response = self.confirmAction()
    if response == 2:
      return
    elif response == 0:
      self.fileSave()

    self.treeView.empty()
    self._fileName = ''
    self.blankWindow()

  def fileNew(self):
    child = self.__class__()
    child.show()
    child.blankWindow()
    self._childWidgets.append(child)
  def fileOpen(self):
    response = self.confirmAction()
    if response == 2:
      return
    elif response == 0:
      self.fileSave()

    fileN = str(QFileDialog.getOpenFileName(self._dirName,
                                            "OPML files (*.opml);;All files (*.*)"))
    if fileN:
      self._dirName,self._fileName = os.path.split(fileN)
      inD = open(fileN,'r').read()
      self.treeView.empty()
      self.treeView.addOPML(inD,new=1)
      

  def editCopy(self):
    item = self.treeView.selectedItem()
    if item:
      txt = self.treeView.toxml(item=item,clearIDs=1)
      clip = qApp.clipboard()
      clip.setText(txt)
    
  def editCopyText(self):
    item = self.treeView.selectedItem()
    if item:
      txt = item.toText()
      clip = qApp.clipboard()
      clip.setText(txt)
    
  def editCut(self):
    item = self.treeView.selectedItem()
    if item:
      txt = self.treeView.toxml(item=item)
      clip = qApp.clipboard()
      clip.setText(txt)
      parent = item.parent() or self.treeView
      parent.takeItem(item)
      if hasattr(item,'_duplicates') and item._duplicates:
        for dupe in item._duplicates.values():
          # dupe is weakref:
          dupe = dupe()
          parent = dupe.parent() or self.treeView
          parent.takeItem(dupe)
        item._duplicates=None
          
          
  def editPaste(self):
    txt = QString()
    clip = qApp.clipboard()
    res = QTextDrag.decode(clip.data(),txt)
    if res:
      txt = str(txt)
      if self.treeView.selectedItem():
        self.treeView.addOPML(txt,where=self.treeView.selectedItem())
      else:
        self.treeView.addOPML(txt)

  def helpContents(self):
    self.helpView = TextView(helpData)
    self.helpView.show()

  def viewExpandAll(self):
    self.treeView.collapseAll(state=1)

  def viewCollapseAll(self):
    self.treeView.collapseAll(state=0)


  def closeEvent(self,evt):
    quit = 1
    if self.treeView.needsSave():
      response = QMessageBox.information(self,sys.argv[0],
                                         "You haven't saved your work\nSave Changes?",
                                         "&Save","&Discard","&Cancel",0,2)
      if response == 0:
        self.fileSave()
      elif response == 2:
        quit = 0
    if quit:
      evt.accept()
    else:
      evt.ignore()

    
    
if __name__=='__main__':
  import sys,getopt
  fileN = 'tests/foo.opml'
  if len(sys.argv) > 1:
    args,extras = getopt.getopt(sys.argv[1:],'')
    if len(extras):
      fileN = extras[0]
  app = QApplication(sys.argv)
  win = RDOutlineImpl(fileName=fileN,mainWin=1)
  app.setMainWidget(win)
  win.show()
  app.exec_loop()

