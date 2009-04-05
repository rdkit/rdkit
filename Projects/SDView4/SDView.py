#!/usr/bin/env python
# $Id$
#
#  Copyright (c) 2008-2009, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Greg Landrum
#
import sys,cStringIO,math
from PyQt4 import QtCore, QtGui, QtSvg
from rdkit.sping.Qt.pidQt4 import QtCanvas
from rdkit.Chem.Draw.MolDrawing import MolDrawing
from rdkit.Chem.Draw.MolDrawing import registerCanvas
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.RDLogger import logger
logger = logger()

class SDUtils(object):
  @staticmethod
  def GetSDPropNames(fName=None,data=None,blockSize=1024,nBlocks=2):
    if not (fName or data): raise ValueError,"provide at least one value"
    import re
    nRead=0
    nmMatcher=re.compile(r'^> +<(\w+)>',re.M)
    res=set()
    if data:
      res.update((x.group(1) for x in nmMatcher.finditer(data)))
    else:
      f = file(fName,'r')
      lines = f.readlines(blockSize)
      while lines and nBlocks:
        res.update((x.group(1) for x in nmMatcher.finditer(''.join(lines))))
        nBlocks-=1
        lines = f.readlines(blockSize)
    return res
      
class MolView(QtGui.QGraphicsView):
  def __init__(self,parent=None,scene=None,size=(300,300)):
    QtGui.QGraphicsView.__init__(self,scene,parent)
    scene = QtGui.QGraphicsScene(0,0,size[0],size[1])
    self.setRenderHints(QtGui.QPainter.Antialiasing|QtGui.QPainter.TextAntialiasing)
    self.setScene(scene)
    self._scenes=[scene]
    registerCanvas('sping')
    self.molDrawer=MolDrawing()
    self.molDrawer.atomLabelMinFontSize=6
    self.mols=[]

  def clearMols(self):
    if self.mols:
      self.molDrawer.canvas.clear()
      self.scene().clear()
    self.mols=[]

  def addMol(self,mol):
    if not hasattr(self.molDrawer,'canvas') or self.molDrawer.canvas is None:
      self.molDrawer.canvas=QtCanvas(self.scene())
    self.molDrawer.AddMol(mol)
    self.mols.append(mol)

  def setMol(self,mol):
    self.clearMols()
    self.molDrawer.canvas=QtCanvas(self.scene())
    self.addMol(mol)

  def setMols(self,mols):
    self.clearMols()
    nMols = len(mols)
    nCols = int(math.ceil(math.sqrt(nMols)))
    nRows = int(math.ceil(1.*nMols/nCols))

    sz = self.size()
    self.scene().setSceneRect(0,0,sz.width(),sz.height())
    
    w = 0.9*self.scene().width()/nCols
    h = 0.9*self.scene().height()/nRows
    self.scene().clear()
    for i,mol in enumerate(mols):
      rowIdx=i//nCols
      colIdx=i%nCols
      self.molDrawer.canvas=QtCanvas(self.scene(),size=(w,h))
      self.molDrawer.scaleAndCenter(mol,mol.GetConformer(),canvasSize=(w,h))
      self.molDrawer.AddMol(mol,drawingTrans=(w*colIdx+w/2,-h*rowIdx+h/2),molTrans=self.molDrawer.molTrans,
                            centerIt=False)
      self.mols.append(mol)
      

class MolTableView(QtGui.QTableView):
  def __init__(self,parent=None):
    QtGui.QTableView.__init__(self,parent)

    
class MolTableModel(QtCore.QAbstractTableModel):
  _supplier=None
  _propNames=None
  _molCache=None
  def flags(self,index):
    return QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled
  def rowCount(self,parent):
    return len(self._supplier)

  def columnCount(self,parent):
    return len(self._propNames)
  def data(self,index,role):
    if not index.isValid() or role!=QtCore.Qt.DisplayRole:
      return QtCore.QVariant()
    row,col = index.row(),index.column()
    pn = self._propNames[col]
    mol = self.mol(index,role)
    if not mol:
      res=''
    elif mol.HasProp(pn):
      res=mol.GetProp(pn)
    else:
      res='NA'
    return QtCore.QVariant(res)
  def mol(self,index,role):
    if not index.isValid() or role!=QtCore.Qt.DisplayRole:
      return None
    row = index.row()
    if self._molCache is not None:
      mol = self._molCache.get(row,None)
      if not mol:
        mol = self._supplier[row]
        AllChem.Compute2DCoords(mol)
        Chem.Kekulize(mol)
        self._molCache[row]=mol
    else:
      mol = self._supplier[row]
    return mol

  def headerData(self,section,orientation,role):
    if orientation==QtCore.Qt.Horizontal:
      if role==QtCore.Qt.DisplayRole:
        return QtCore.QVariant(self._propNames[section])
    return QtCore.QVariant()
  def setSupplier(self,supplier,propNames=None):
    self._supplier=supplier
    if propNames:
      self._propNames=propNames
    else:
      for m in supplier:
        if m:
          self._propNames=['_Name']+list(m.GetPropNames())
      supplier.reset()
  def sort(self,col,order):
    print 'SORT!',col,order
    
class MainWindow(QtGui.QMainWindow):
  _applicationName="SDView"
  _vendorName="rdkit.org"
  def __init__(self):
    QtGui.QMainWindow.__init__(self)
    self.curFile = QtCore.QString()

    
    #self.textEdit = QtGui.QTextEdit()
    #self.setCentralWidget(self.textEdit)

    self.createActions()
    self.createMenus()
    self.createStatusBar()

    self.readSettings()

    fName= sys.argv[1]
    suppl=Chem.SDMolSupplier(fName)
    pns = ['_Name']+list(SDUtils.GetSDPropNames(fName=fName))
    self.sdModel=MolTableModel()
    self.sdModel.setSupplier(suppl,propNames=pns)
    self.sdModel._molCache={}

    self.tblView = MolTableView()
    
    self.tblView.setModel(self.sdModel)
    self.setCentralWidget(self.tblView)
    self.connect(self.tblView.selectionModel(),
                 QtCore.SIGNAL('selectionChanged(QItemSelection,QItemSelection)'),
                 self.changeCurrent)
    self.molView=MolView()
    self.molView.show()


  def changeCurrent(self,current,prev):
    idxL = self.tblView.selectedIndexes()
    if not len(idxL):
      return
    mols = []
    for i,idx in enumerate(idxL):
      mol = self.sdModel.mol(idx,QtCore.Qt.DisplayRole)
      mols.append(mol)
    if len(mols)==1:
      self.molView.setMol(mols[0])
    else:
      self.molView.setMols(mols)
    
  def closeEvent(self, event):
    if self.maybeSave():
      self.writeSettings()
      self.molDrawer=None
      self.molView=None
      event.accept()
    else:
      event.ignore()

  def newFile(self):
    if self.maybeSave():
      pass

  def open(self):
    if self.maybeSave():
      fileName = QtGui.QFileDialog.getOpenFileName(self)
      if not fileName.isEmpty():
        pass

  def save(self):
    if self.curFile.isEmpty():
      return self.saveAs()
    else:
      return True

  def saveAs(self):
    fileName = QtGui.QFileDialog.getSaveFileName(self)
    if fileName.isEmpty():
      return False

    return True

  def about(self):
    QtGui.QMessageBox.about(self, self.tr("About Application"),
        self.tr("The <b>Application</b> example demonstrates how to "
                "write modern GUI applications using Qt, with a menu bar, "
                "toolbars, and a status bar."))

  def documentWasModified(self):
    #self.setWindowModified(self.textEdit.document().isModified())
    pass

  def createActions(self):
    self.newAct = QtGui.QAction(self.tr("&New"), self)
    self.newAct.setShortcut(self.tr("Ctrl+N"))
    self.newAct.setStatusTip(self.tr("Create a new file"))
    self.connect(self.newAct, QtCore.SIGNAL("triggered()"), self.newFile)

    self.openAct = QtGui.QAction(self.tr("&Open..."), self)
    self.openAct.setShortcut(self.tr("Ctrl+O"))
    self.openAct.setStatusTip(self.tr("Open an existing file"))
    self.connect(self.openAct, QtCore.SIGNAL("triggered()"), self.open)

    self.saveAct = QtGui.QAction(self.tr("&Save"), self)
    self.saveAct.setShortcut(self.tr("Ctrl+S"))
    self.saveAct.setStatusTip(self.tr("Save the document to disk"))
    self.connect(self.saveAct, QtCore.SIGNAL("triggered()"), self.save)

    self.saveAsAct = QtGui.QAction(self.tr("Save &As..."), self)
    self.saveAsAct.setStatusTip(self.tr("Save the document under a new name"))
    self.connect(self.saveAsAct, QtCore.SIGNAL("triggered()"), self.saveAs)

    self.exitAct = QtGui.QAction(self.tr("E&xit"), self)
    self.exitAct.setShortcut(self.tr("Ctrl+Q"))
    self.exitAct.setStatusTip(self.tr("Exit the application"))
    self.connect(self.exitAct, QtCore.SIGNAL("triggered()"), self, QtCore.SLOT("close()"))

    self.cutAct = QtGui.QAction(self.tr("Cu&t"), self)
    self.cutAct.setShortcut(self.tr("Ctrl+X"))
    self.cutAct.setStatusTip(self.tr("Cut the current selection's contents to the "
                                     "clipboard"))
    #self.connect(self.cutAct, QtCore.SIGNAL("triggered()"), self.textEdit, QtCore.SLOT("cut()"))

    self.copyAct = QtGui.QAction(self.tr("&Copy"), self)
    self.copyAct.setShortcut(self.tr("Ctrl+C"))
    self.copyAct.setStatusTip(self.tr("Copy the current selection's contents to the "
                                      "clipboard"))
    #self.connect(self.copyAct, QtCore.SIGNAL("triggered()"), self.textEdit, QtCore.SLOT("copy()"))

    self.pasteAct = QtGui.QAction(self.tr("&Paste"), self)
    self.pasteAct.setShortcut(self.tr("Ctrl+V"))
    self.pasteAct.setStatusTip(self.tr("Paste the clipboard's contents into the current "
                                       "selection"))
    #self.connect(self.pasteAct, QtCore.SIGNAL("triggered()"), self.textEdit, QtCore.SLOT("paste()"))
    
    self.aboutAct = QtGui.QAction(self.tr("&About"), self)
    self.aboutAct.setStatusTip(self.tr("Show the application's About box"))
    self.connect(self.aboutAct, QtCore.SIGNAL("triggered()"), self.about)

    self.aboutQtAct = QtGui.QAction(self.tr("About &Qt"), self)
    self.aboutQtAct.setStatusTip(self.tr("Show the Qt library's About box"))
    self.connect(self.aboutQtAct, QtCore.SIGNAL("triggered()"), QtGui.qApp, QtCore.SLOT("aboutQt()"))

    self.cutAct.setEnabled(False)
    self.copyAct.setEnabled(False)
    self.pasteAct.setEnabled(False)
    #self.connect(self.textEdit, QtCore.SIGNAL("copyAvailable(bool)"),
    #             self.cutAct, QtCore.SLOT("setEnabled(bool)"))
    #self.connect(self.textEdit, QtCore.SIGNAL("copyAvailable(bool)"),
    #             self.copyAct, QtCore.SLOT("setEnabled(bool)"))

  def createMenus(self):
    self.fileMenu = self.menuBar().addMenu(self.tr("&File"))
    self.fileMenu.addAction(self.newAct)
    self.fileMenu.addAction(self.openAct)
    self.fileMenu.addAction(self.saveAct)
    self.fileMenu.addAction(self.saveAsAct)
    self.fileMenu.addSeparator();
    self.fileMenu.addAction(self.exitAct)

    self.editMenu = self.menuBar().addMenu(self.tr("&Edit"))
    self.editMenu.addAction(self.cutAct)
    self.editMenu.addAction(self.copyAct)
    self.editMenu.addAction(self.pasteAct)

    self.menuBar().addSeparator()

    self.helpMenu = self.menuBar().addMenu(self.tr("&Help"))
    self.helpMenu.addAction(self.aboutAct)
    self.helpMenu.addAction(self.aboutQtAct)

  def createStatusBar(self):
    self.statusBar().showMessage(self.tr("Ready"))

  def readSettings(self):
    settings = QtCore.QSettings(self._vendorName,self._applicationName)
    pos = settings.value("pos", QtCore.QVariant(QtCore.QPoint(200, 200))).toPoint()
    size = settings.value("size", QtCore.QVariant(QtCore.QSize(400, 400))).toSize()
    self.resize(size)
    self.move(pos)

  def writeSettings(self):
    settings = QtCore.QSettings(self._vendorName,self._applicationName)
    settings.setValue("pos", QtCore.QVariant(self.pos()))
    settings.setValue("size", QtCore.QVariant(self.size()))

  def maybeSave(self):
    if False:
      ret = QtGui.QMessageBox.warning(self, self.tr("Application"),
                  self.tr("The document has been modified.\n"
                          "Do you want to save your changes?"),
                  QtGui.QMessageBox.Yes | QtGui.QMessageBox.Default,
                  QtGui.QMessageBox.No,
                  QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Escape)
      if ret == QtGui.QMessageBox.Yes:
        return self.save()
      elif ret == QtGui.QMessageBox.Cancel:
        return False
    return True

  def loadFile(self, fileName):
    file = QtCore.QFile(fileName)
    if not file.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text):
      QtGui.QMessageBox.warning(self, self.tr("Application"),
                  self.tr("Cannot read file %1:\n%2.").arg(fileName).arg(file.errorString()))
      return

    #inf = QtCore.QTextStream(file)
    #QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
    #self.textEdit.setPlainText(inf.readAll())
    #QbtGui.QApplication.restoreOverrideCursor()

    #self.setCurrentFile(fileName)
    #self.statusBar().showMessage(self.tr("File loaded"), 2000)

  def saveFile(self, fileName):
    file = QtCore.QFile(fileName)
    if not file.open(QtCore.QFile.WriteOnly | QtCore.QFile.Text):
      QtGui.QMessageBox.warning(self, self.tr("Application"),
                  self.tr("Cannot write file %1:\n%2.").arg(fileName).arg(file.errorString()))
      return False

    #outf = QtCore.QTextStream(file)
    #QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
    #outf << self.textEdit.toPlainText()
    #QtGui.QApplication.restoreOverrideCursor()

    #self.setCurrentFile(fileName);
    #self.statusBar().showMessage(self.tr("File saved"), 2000)
    return True


if __name__ == "__main__":
  app = QtGui.QApplication(sys.argv)
  mainWin = MainWindow()
  mainWin.show()
  sys.exit(app.exec_())
