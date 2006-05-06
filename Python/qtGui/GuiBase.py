#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" base class for the main window in the qt MixinGui

"""
import RDConfig
from qt import *
from qtGui import GuiShell,GuiTextViewer
from qtGui import rpcClient
from qtGui import qtUtils

class GuiBase(QMainWindow):
  """ Base class for the qt MixinGui

    This opens a QMainWindow with menu and status bars.

    optionally, an rpcServer will be started as well

  """
  def __init__(self,parent=None,name="GuiBase",startServer=False,includeLogo=True):
    apply(QMainWindow.__init__,(self,parent,name,Qt.WDestructiveClose))
    self.setCaption(self.trUtf8(name))
    self._initMenus()
    self.statusBar()
    if startServer:
      client = rpcClient.rpcClient(qApp,self)
      self.servThread = rpcClient.ServThread(client)
      self.servThread.start()
    #self.resize(400,400)
    self._shell = GuiShell.PyShell(locals={'GUI':self})
    self._helpWin = GuiTextViewer.GuiTextViewer()
    self._aboutWin = GuiTextViewer.GuiTextViewer()
    if includeLogo:
      self._addLogo()
    self.adjustSize()
    image0 = QPixmap(qtUtils.logoImageData)
    self.setIcon(image0)


  def _addLogo(self):
    if not self.centralWidget():
      self.setCentralWidget(QWidget(self,'central_widget'))

    # code to build the logo pixmap
    self._displayLabel = QLabel(self.centralWidget(),'picture')
    self._displayLabel.setScaledContents(0)
    self._logoPM = QPixmap()
    self._logoPM.load(RDConfig.RDDocsDir+'/logo.jpg')
    self._displayLabel.setPixmap(self._logoPM)

    # now insert it nicely
    layout = QHBoxLayout(self.centralWidget(),1,1,'labellayout')
    spacer = QSpacerItem(2,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout.addItem(spacer)
    layout.addWidget(self._displayLabel)
    spacer = QSpacerItem(2,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
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
    self._fileMenu = QPopupMenu(self)
    self._editMenu = QPopupMenu(self)
    self._viewMenu = QPopupMenu(self)
    self._helpMenu = QPopupMenu(self)
 
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
    self._fileMenu.insertItem(self.trUtf8('&Quit'),qApp,SLOT("quit()"),Qt.CTRL+Qt.Key_Q )
    self.menuBar().insertItem(self.trUtf8('&Help'),self._helpMenu)

    self._helpMenu.insertItem(self.trUtf8('&About'),self.aboutBox)

  def closeEvent(self,evt):
    """ callback for when we receive a close event

    """
    self._shell.close(1)
    self._helpWin.close(1)
    # FIX: this ain't right
    qApp.quit()

  def customEvent(self,evt):
    """ event handler for custom events

    """
    if evt.type() == rpcClient._rpcEventID:
      rpcClient.HandleRPCEvent(self,evt)

    

if __name__ == '__main__':
  import sys
  a = QApplication(sys.argv)
  widg = GuiBase(startServer=0)
  widg.finalizeMenus()
  widg.show()
  a.setMainWidget(widg)
  a.exec_loop()


