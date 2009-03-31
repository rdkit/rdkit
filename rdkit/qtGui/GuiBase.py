# $Id$
#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" base class for the main window in the qt MixinGui

"""
from rdkit import RDConfig
import qt
from rdkit.qtGui import GuiShell,GuiTextViewer
try:
  from rdkit.qtGui import rpcClient
except:
  rpcClient=None
from rdkit.qtGui import qtUtils

class GuiBase(qt.QMainWindow):
  """ Base class for the qt MixinGui

    This opens a QMainWindow with menu and status bars.

    optionally, an rpcServer will be started as well

  """
  def __init__(self,parent=None,name="GuiBase",startServer=False,includeLogo=True):
    apply(qt.QMainWindow.__init__,(self,parent,name,qt.Qt.WDestructiveClose))
    self.setCaption(self.trUtf8(name))
    self._initMenus()
    self.statusBar()
    if startServer and rpcClient:
      client = rpcClient.rpcClient(qt.qApp,self)
      self.servThread = rpcClient.ServThread(client)
      self.servThread.start()
    #self.resize(400,400)
    self._shell = GuiShell.PyShell(locals={'GUI':self})
    self._helpWin = GuiTextViewer.GuiTextViewer()
    self._aboutWin = GuiTextViewer.GuiTextViewer()
    if includeLogo:
      self._addLogo()
    self.adjustSize()
    image0 = qt.QPixmap(qtUtils.logoImageData)
    self.setIcon(image0)


  def _addLogo(self):
    if not self.centralWidget():
      self.setCentralWidget(qt.QWidget(self,'central_widget'))

    # code to build the logo pixmap
    self._displayLabel = qt.QLabel(self.centralWidget(),'picture')
    self._displayLabel.setScaledContents(0)
    self._logoPM = qt.QPixmap()
    self._logoPM.load(RDConfig.RDDocsDir+'/Images/logo.png')
    self._displayLabel.setPixmap(self._logoPM)

    # now insert it nicely
    layout = qt.QHBoxLayout(self.centralWidget(),1,1,'labellayout')
    spacer = qt.QSpacerItem(2,0,qt.QSizePolicy.Expanding,qt.QSizePolicy.Minimum)
    layout.addItem(spacer)
    layout.addWidget(self._displayLabel)
    spacer = qt.QSpacerItem(2,0,qt.QSizePolicy.Expanding,qt.QSizePolicy.Minimum)
    layout.addItem(spacer)

  def aboutBox(self):
    """ displays our about box (in the _aboutWin_)

    """
    self._aboutWin.setSource(RDConfig.RDDocsDir+'/about.qt.html')
    self._aboutWin.show()
    
  def launchHelp(self):
    """ displays our _helpWin_)

    """
    self._helpWin.show()
    
  def launchShell(self):
    """ displays our _shell_ window

    """
    self._shell.show()
    
  def _initMenus(self):
    """ INTERNAL USE ONLY

    """
    self._fileMenu = qt.QPopupMenu(self)
    self._editMenu = qt.QPopupMenu(self)
    self._viewMenu = qt.QPopupMenu(self)
    self._helpMenu = qt.QPopupMenu(self)
 
    self._fileMenu._id=self.menuBar().insertItem(self.trUtf8('&File'),self._fileMenu)
    self._editMenu._id=self.menuBar().insertItem(self.trUtf8('&Edit'),self._editMenu)
    self._viewMenu._id=self.menuBar().insertItem(self.trUtf8('&View'),self._viewMenu)

  def finalizeMenus(self,includeShell=True):
    """ finalizes our menu bar

      This adds:

        - File->PyShell

        - File->Quit

        - the Help menu

        - Help->About

    """
    if includeShell:
      self._fileMenu._pyShellId=self._fileMenu.insertItem(self.trUtf8('&PyShell'),self.launchShell)
    self._fileMenu.insertItem(self.trUtf8('&Quit'),qt.qApp,qt.SLOT("quit()"),qt.Qt.CTRL+qt.Qt.Key_Q )
    self.menuBar().insertItem(self.trUtf8('&Help'),self._helpMenu)

    self._helpMenu.insertItem(self.trUtf8('&About'),self.aboutBox)

  def closeEvent(self,evt):
    """ callback for when we receive a close event

    """
    self._shell.close(1)
    self._helpWin.close(1)
    # FIX: this ain't right
    qt.qApp.quit()

  def customEvent(self,evt):
    """ event handler for custom events

    """
    if rpcClient and evt.type() == rpcClient._rpcEventID:
      rpcClient.HandleRPCEvent(self,evt)

    

if __name__ == '__main__':
  import sys
  qApp = qt.QApplication(sys.argv)
  print qt.qApp
  widg = GuiBase(startServer=0)
  widg.finalizeMenus()
  widg.show()
  qApp.setMainWidget(widg)
  qApp.exec_loop()


