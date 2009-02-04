#
#  Copyright (C) 2003  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Allows browsing and manipulating sets of molecules

"""
import RDConfig
from qt import *
from qtGui.GuiLib.MolBrowserImpl import MolBrowser
import os

REQUIRED_MIXINS = ['PiddleCanvas']
MODULES_ALTERED = ['qtGui.GuiBase']
METHODS_DEFINED = {
  '__LaunchBrowserWin':'qtGui.GuiBase.GuiBase.mbLaunchBrowserWin',
  }

def __LaunchBrowserWin(self):
  if self._mbBrowserWin is None:
    win = MolBrowser()
    self._mbBrowserWin = win
  else:
    win = self._mbBrowserWin
  win.show()
  
def LocalInit(self):
  self._mbBrowserWin = None

  self._mbMolMenu = QPopupMenu(self)
  self.menuBar().insertItem(self.trUtf8('&Molecules'),self._mbMolMenu)
  self._mbMolMenu.insertItem(self.trUtf8('&Browser'),self.mbLaunchBrowserWin)


                                
                                
